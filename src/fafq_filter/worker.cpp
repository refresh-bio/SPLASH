#include "worker.h"
#include <fstream>

// ******************************************************************************************
bool CWorker::load_fastx(const string& in_fn)
{
	refresh::stream_in_buffered* f_in;

	if (in_fn.empty())
		f_in = new refresh::stream_in_stdin();
	else
	{
		f_in = new refresh::stream_in_file(in_fn);

		if (!(dynamic_cast<refresh::stream_in_file*>(f_in)->is_open()))
		{
			cerr << "Cannot open file " + in_fn + "\n";
			delete f_in;

			return false;
		}
	}

	refresh::stream_decompression sd_in(f_in);
	bool is_fasta = false;
	bool read_allowed = true;
	string line;

	stats_n_reads = 0;
	stats_n_bases = 0;
	stats_n_allowed_reads = 0;
	stats_n_allowed_bases = 0;

	while (!sd_in.eof())
	{
		if (sd_in.getline(line) < 0)
			break;

		if (line.empty())
			continue;

		if (line.front() == '>')
		{
			is_fasta = true;
			read_allowed = !params.check_header || regex_match(line, params.allowed_headers);
			++stats_n_reads;
		}
		else if (line.front() == '@')
		{
			is_fasta = false;
			read_allowed = !params.check_header || regex_match(line, params.allowed_headers);
			++stats_n_reads;
		}
		else
		{
			cerr << "Strange input line in: " + in_fn + " : " + line + "\n";
			delete f_in;

			return false;
		}

		if (sd_in.getline(line) < 0)
		{
			cerr << "Truncated file: " + in_fn + "\n";
			delete f_in;

			return false;
		}

		if (read_allowed)
		{
			look_for_anchors(line);
			++stats_n_allowed_reads;
			stats_n_allowed_bases += line.size();
			stats_n_bases += line.size();
		}
		else
			stats_n_bases += line.size();

		if (!is_fasta)
		{
			if (sd_in.getline(line) < 0)
			{
				cerr << "Truncated file: " + in_fn + "\n";
				delete f_in;

				return false;
			}
			if (sd_in.getline(line) < 0)
			{
				cerr << "Truncated file: " + in_fn + "\n";
				delete f_in;

				return false;
			}
		}
	}

	delete f_in;

	return true;
}

// ******************************************************************************************
bool CWorker::process_file(const uint32_t sample_id, const string& in_fn, const string &out_fn)
{
	anchor_stats.clear();

	if (!load_fastx(in_fn))
		return false;

	if (!store_satc(sample_id, out_fn))
		return false;

	return true;
}

// ******************************************************************************************
bool CWorker::process_stdin(const uint32_t sample_id, const string &out_fn)
{
	anchor_stats.clear();

	if (!load_fastx(""))
		return false;

	if (!store_satc(sample_id, out_fn))
		return false;

	return true;
}

// ******************************************************************************************
bool CWorker::store_satc(const uint32_t sample_id, const string& out_fn)
{
	size_t max_cnt = 1;

	vector<kmer_t> sorted_anchors;
	sorted_anchors.reserve(anchor_stats.size());

	size_t stats_n_anchors_reported = 0;
	size_t stats_n_targets_reported = 0;

	for (const auto& anchor_info : anchor_stats)
	{
		sorted_anchors.push_back(anchor_info.first);

		for (const auto& tc : anchor_info.second)
			if (max_cnt < tc.second)
				max_cnt = tc.second;
	}

	sort(sorted_anchors.begin(), sorted_anchors.end());

	buffered_binary_writer satc_writer(out_fn);

	Header header;

	header.sample_id_size_bytes = no_bytes(sample_id);
	header.barcode_size_bytes = 0;
	header.anchor_size_bytes = (params.anchor_len + 3) / 4;
	header.target_size_bytes = (params.target_len + 3) / 4;
	header.counter_size_bytes = no_bytes(max_cnt);
	header.barcode_len_symbols = 0;
	header.anchor_len_symbols = params.anchor_len;
	header.target_len_symbols = params.target_len;
	header.gap_len_symbols = 0;
	header.ordering = Header::ordering_t::SBATC;

	header.serialize(satc_writer);

	Record rec;
	rec.sample_id = sample_id;
	rec.barcode = 0;

	auto cmp = [](auto& x, auto& y) {
		if (x.second != y.second)
			return x.second > y.second;
		return x.first < y.first;
		};

	for (const auto& anchor : sorted_anchors) {
		const auto& anchor_info = anchor_stats[anchor];

		vector<pair<kmer_t, size_t>> target_cnt(anchor_info.begin(), anchor_info.end());

		if (target_cnt.size() < params.n_top_targets)
			sort(target_cnt.begin(), target_cnt.end(), cmp);
		else
		{
			partial_sort(target_cnt.begin(), target_cnt.begin() + params.n_top_targets, target_cnt.end(), cmp);
			target_cnt.resize(params.n_top_targets);
		}

		rec.anchor = anchor;

		for (const auto& tc : target_cnt)
		{
			rec.count = tc.second;
			rec.target = tc.first;

			rec.serialize(satc_writer, header);
		}

		if (!target_cnt.empty())
		{
			++stats_n_anchors_reported;
			stats_n_targets_reported += target_cnt.size();
		}
	}

	if (!params.stats_json_fn.empty())
	{
		ofstream jf(params.stats_json_fn);
		if (!jf)
		{
			cerr << "Cannot create stats file: " << params.stats_json_fn << endl;
			exit(1);
		}

		jf << "{" << endl;
		jf << "\t\"n_reads\": \"" << stats_n_reads << "\"," << endl;
		jf << "\t\"n_bases\": \"" << stats_n_bases << "\"," << endl;
		jf << "\t\"n_allowed_reads\": \"" << stats_n_allowed_reads << "\"," << endl;
		jf << "\t\"n_allowed_bases\": \"" << stats_n_allowed_bases << "\"," << endl;
		jf << "\t\"n_anchors_reported\": \"" << stats_n_anchors_reported << "\"," << endl;
		jf << "\t\"n_targets_reported\": \"" << stats_n_targets_reported << "\"" << endl;
		jf << "}" << endl;
	}


	//dtor will close

	return true;
}

// ******************************************************************************************
void CWorker::prepare_kmers(const string& str)
{
	rec_kmers.clear();

	if (str.size() < params.anchor_len + params.target_len)
		return;

	kmer_t anchor = 0;
	int no_valid_symbols = 0;

	for (int i = 0; i < str.size() - params.target_len; ++i)
	{
		kmer_t code = char2bits[(uint8_t)str[i]];
		if (code > 3)
		{
			anchor = 0;
			no_valid_symbols = 0;
			if (i >= params.anchor_len)
				rec_kmers.emplace_back((kmer_t)~0ull);
		}
		else
		{
			anchor <<= 2;
			anchor += code;
			++no_valid_symbols;

			if (i + 1 >= params.anchor_len)
			{
				if (no_valid_symbols >= params.anchor_len)
					rec_kmers.emplace_back(anchor & anchor_mask);
				else
					rec_kmers.emplace_back((kmer_t)~0ull);
			}
		}
	}
}

// ******************************************************************************************
kmer_t CWorker::get_target(const string& str, size_t pos)
{
	kmer_t target = 0;

	for (size_t i = pos; i < pos + params.target_len; ++i)
	{
		kmer_t code = char2bits[(uint8_t)str[i]];

		if (code > 3)
			return empty_mask;

		target <<= 2;
		target += code;
	}

	return target;
}

// ******************************************************************************************
void CWorker::look_for_anchors(const string& str)
{
	prepare_kmers(str);

	for (size_t i = 0; i < rec_kmers.size(); ++i)
		if(rec_kmers[i] != empty_mask && accepted_anchors->IsAccepted(rec_kmers[i]))
			anchor_stats[rec_kmers[i]][get_target(str, i + params.anchor_len)]++;
}
