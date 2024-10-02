#include "io.h"
#include "log.h"

#include <fstream>
#include <algorithm>
#include <filesystem>
#include <stdexcept>

using namespace std;
namespace fs = std::filesystem;

ReadLoader::ReadLoader(
	const std::vector<std::string>& fastqFiles,
	input_format_t inputFormat,
	const std::string& anchorFile,
	int homopolymerThreshold,
	int numThreads,
	int readsBufferGb,
	size_t anchorsBatchSize,
	bool keepTemp,
	string tempPath)
	: polyFilter(homopolymerThreshold), anchorFile(anchorFile), anchorsBatchSize(anchorsBatchSize), tempPath(tempPath) {

	std::vector<string> outFastqFiles;

	// check if file exist
	for (const auto& file: fastqFiles) {
		if (!fs::exists(file)) {
			LOG_NORMAL << "Warning: skipping unexisting file " << file << std::endl;
		}
		else {
			outFastqFiles.push_back(file);
		}
	}

	selector.set_input_names(outFastqFiles);
	selector.set_input_format(inputFormat);

	selector.set_max_memory_usage((size_t)readsBufferGb << 30);		// I/O buffer
	selector.set_max_read_len(1 << 20);				// optional 
	selector.set_no_threads(numThreads);					
	selector.set_keep_temps(keepTemp);
	
	fs::path tmp_dir(tempPath);
	bool ok = false;
	if (!fs::exists(tmp_dir)) {
		ok = fs::create_directory(tmp_dir);
	}
	else {
		ok = fs::is_directory(tmp_dir);
	}
	
	if (!ok) {
		throw runtime_error("Unable to create/open temporary directory: " + tempPath);
	}

	selector.set_output_dir(tmp_dir.string());
}

ReadLoader::~ReadLoader() {
	if (!selector.get_keep_temps()) {
		fs::remove_all(fs::path(tempPath));
	}
}

size_t ReadLoader::extractKmers(
	int numFollowers,
	int followerLen,
	bool allAnchors,
	std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
	int& anchorLen) {

	std::unordered_map<kmer_t, bool> anchors;
	
	// load anchors from file
	ifstream ifs;
	const int N = 16 << 20; // 16MB buffer
	char* buf = new char[N];
	ifs.rdbuf()->pubsetbuf(buf, N);

	ifs.open(anchorFile);

	if (!ifs) {
		throw(std::runtime_error("Unable to open anchor file: " + anchorFile));
	}

	string line;
	anchorLen = 0;

	LOG_NORMAL << "Loading anchors..." << endl;

	// read an optional header and a first anchor to establish length
	size_t line_number = 1;
	
	for (; anchorLen == 0 && getline(ifs, line); ++line_number) {
		if (line.length() > 0) {
			stringstream ss;
			string token;
			ss << line;
			ss >> token;
			if (token != "anchor") {
				kmer_t kmer;
				anchorLen = KmerHelper::from_string(line.c_str(), token.length(), kmer);
				if (anchorLen != token.length() || anchorLen > 31) {
					delete[] buf;
					throw std::runtime_error("Incorrect anchor at line " + to_string(line_number) + " (at most 31-mers are supported)");
				}
				anchors[kmer] = false;
			}
		}
	}

	// read remaining anchors
	for (; getline(ifs, line); ++line_number) {
		if (line.length() == 0) {
			continue;
		}

		if (line_number % 10000 == 0) {
			LOG_NORMAL << "\r" << line_number << "..." << std::flush;
		}

		kmer_t kmer;
		int curLen = KmerHelper::from_string(line.c_str(), line.length(), kmer);
		
		if (curLen != anchorLen) {
			delete[] buf;
			throw std::runtime_error("Incorrect anchor at line " + to_string(line_number));
		}
		anchors[kmer] = false;
	}
	LOG_NORMAL << "\r" << anchors.size() << " anchors loaded\n" << std::flush;

	delete[] buf;
	
	// run read selector
	LOG_NORMAL << "Selecting reads in batches..." << endl;

	size_t n_passed = 0;
	size_t n_batches = anchors.size() / anchorsBatchSize + 1;
	
	auto begin = anchors.begin();

	for (int i = 0; i < n_batches; ++i) {
		LOG_NORMAL << "\r" << i + 1 << "/" << n_batches << std::flush;
		
		// activate batch of anchors
		auto it = begin;
		for (int j = 0; (j < anchorsBatchSize) && (it != anchors.end()); ++j, ++it) {
			it->second = true;
		}

		selector.set_dict(anchors);
		auto tmp = selector.process_anchor_followers(ReadSelector::direction_t::forward, anchorLen, followerLen, numFollowers, allAnchors);
		n_passed += translate(tmp, numFollowers, targets);

		// deactivate batch
		auto end = it;
		for (it = begin; it != end; ++it) {
			it->second = false;
		}

		begin = end;
	}
	LOG_NORMAL << "\rAnchor batches processed successfully\n" << std::flush;

	return n_passed;
}


size_t ReadLoader::extractKmers(
	const std::unordered_set<kmer_t>& queries,
	int queryLen,
	int numFollowers,
	int followerLen,
	bool reverse,
	bool allAnchors,
	std::unordered_map<kmer_t, std::vector<kmer_t>>& targets) {

	ReadSelector::direction_t dir = reverse ? ReadSelector::direction_t::reverse : ReadSelector::direction_t::forward;
	
	std::unordered_map<kmer_t, bool> query;
	for (auto a : queries) {
		query[a] = false;
	}

	size_t n_passed = 0;
	size_t n_batches = query.size() / anchorsBatchSize + 1;

	auto begin = query.begin();

	for (int i = 0; i < n_batches; ++i) {
		// activate batch of anchors
		auto it = begin;
		for (int j = 0; (j < anchorsBatchSize) && (it != query.end()); ++j, ++it) {
			it->second = true;
		}

		selector.set_dict(query);
		auto tmp = selector.process_anchor_followers(dir, queryLen, followerLen, numFollowers, allAnchors);
		n_passed += translate(tmp, numFollowers, targets);

		// deactivate batch
		auto end = it;
		for (it = begin; it != end; ++it) {
			it->second = false;
		}

		begin = end;
	}

	return n_passed;
}


size_t ReadLoader::extractPredecessorWithGap(
	const std::unordered_set<kmer_t>& queries,
	int queryLen,
	int gapLen,
	int predecessorLen,
	std::unordered_map<kmer_t, std::vector<kmer_t>>& targets) {

	
	std::unordered_map<kmer_t, bool> query;
	for (auto a : queries) {
		query[a] = false;
	}

	size_t n_passed = 0;
	size_t n_batches = query.size() / anchorsBatchSize + 1;

	auto begin = query.begin();

	for (int i = 0; i < n_batches; ++i) {
		// activate batch of anchors
		auto it = begin;
		for (int j = 0; (j < anchorsBatchSize) && (it != query.end()); ++j, ++it) {
			it->second = true;
		}

		selector.set_dict(query);
		auto tmp = selector.process_anchor_extender(ReadSelector::direction_t::forward, predecessorLen, queryLen, gapLen, true);
		n_passed += translate(tmp, 1, targets);

		// deactivate batch
		auto end = it;
		for (it = begin; it != end; ++it) {
			it->second = false;
		}

		begin = end;
	}

	return n_passed;



}

size_t ReadLoader::translate(
	ReadSelector::dict_t& in,
	int numKmers,
	std::unordered_map<kmer_t, std::vector<kmer_t>>& out) const
{
	size_t n_passed = 0;
	size_t n_total = 0;

	for (auto& entry : in) {

		std::vector<kmer_t>& out_v = out[entry.first];

		for (const vector<uint64_t>& sub_v : entry.second) {

			n_passed += sub_v.size() / numKmers;
			out_v.insert(out_v.end(), sub_v.begin(), sub_v.end());

			/*
			n_total += sub_v.size() / numKmers;
			for (auto it = sub_v.begin(); it < sub_v.end(); it += numKmers) {
				bool ok = none_of(it, it + numKmers, [this](kmer_t kmer)->bool {
					return kmer == EMPTY_KMER || polyFilter.IsPolyACGT(kmer, kmerLen);
					});

				if (ok) {
					out_v.insert(out_v.end(), it, it + numKmers);
					++n_passed;
				}
			}
			*/

		}

		// force memory free
		std::vector<std::vector<kmer_t>>().swap(entry.second);

		if (out_v.size() == 0) {
			out.erase(entry.first);
		}
	}


	LOG_DEBUG << "Finished" << endl
		<< out.size() << " anchors, "
		<< n_passed << " reads (" << n_total << " before filtring)" << endl;

	return n_passed;
}



Output::Output(const std::string& table, const std::string& fasta, bool noSubcompactors, bool cumulatedStats)
	: tableFile(table, std::ios_base::out), fastaFile(fasta, std::ios_base::out), writerQueue(1024, 1), noSubcompactors(noSubcompactors), cumulatedStats(cumulatedStats) {
	if (!tableFile) {
		throw std::runtime_error("Unable to open output file: " + table);
	}

	tableBuffer.resize(16 << 20);
	fastaBuffer.resize(16 << 20);

	std::string header{ 
		"anchor\t"
		"compactor\t"
		"id\t"
		"parent_id\t"
		"support\t"
		"exact_support\t"
		"extender_specificity\t"
		"extender_shift\t"
		"total_length\t"
		"num_extended\t"
		"expected_read_count"
	};

	if (cumulatedStats) {
		header += "\t"
			"cumulated_id\t"
			"cumulated_exact_support\t"
			"cumulated_extender_specificity";
	}

	header += "\n";

	
	tableFile.write(header.c_str(), header.length());

	worker = std::thread([this]() {
		task_t task;
		while (this->writerQueue.pop(task)) {
			this->internal_save(task.begin, task.end);
		}
	});

}

void Output::internal_save(std::deque<Compactor>::iterator begin, std::deque<Compactor>::iterator end) {

	
	char* p_table = tableBuffer.data();
	char* p_fasta = fastaBuffer.data();

	for (auto it = begin; it != end; ++it) {
		Compactor& c = *it;

		// omit subcompactors if disabled
		if (noSubcompactors && c.num_children > 0) {
			continue;
		}

		if (c.total_support >= 0) { // omit anchors which have total support set to -1

			// set identifier
			c.id = numCompactors++;
			
			size_t bytes_needed = (c.num_kmers_total + 1) * c.k  + (1 << 16); // +1 due to anchor column

			// no more space in buffer - write to file
			if (p_table + bytes_needed >= tableBuffer.data() + tableBuffer.size()) {
				tableFile.write(tableBuffer.data(), p_table - tableBuffer.data());
				p_table = tableBuffer.data();
			}

			// leave empty place for an anchor
			char* anchor_begin = p_table;
			p_table += c.ancestor->k;
			*p_table = '\t';
			++p_table;

			// extract compactor sequence 
			char* seq_begin = p_table;
			p_table = c.to_string(p_table, false);
			char* seq_end = p_table;

			// copy anchor
			copy_n(seq_begin, c.ancestor->k, anchor_begin);

			// fill rest of columns
			int n_extensions = ((c.num_kmers_total - 1) / c.num_kmers) - 1;
			int total_len = seq_end - seq_begin;
			p_table += sprintf(p_table, "\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%f",
				c.id,
				c.get_parent_id(),
				c.total_support,
				c.exact_support,
				c.extender_specificity,
				c.extender_shift,
				total_len,
				n_extensions,
				c.calculate_expected_read_count());

			if (cumulatedStats) {
				*p_table = '\t';
				++p_table;
				p_table += c.print_cumulated_id(p_table);
				*p_table = '\t';
				++p_table;
				p_table += c.print_cumulated_exact_support(p_table);
				*p_table = '\t';
				++p_table;
				p_table += c.print_cumulated_extender_specificity(p_table);
			}

			*p_table = '\n';
			++p_table;


			if (fastaFile) {
				// no more space in buffer - write to file
				if (p_fasta + bytes_needed >= fastaBuffer.data() + fastaBuffer.size()) {
					fastaFile.write(fastaBuffer.data(), p_fasta - fastaBuffer.data());
					p_fasta = fastaBuffer.data();
				}

				p_fasta += sprintf(p_fasta, ">compactor-%d\n", c.id);
				p_fasta = copy(seq_begin, seq_end, p_fasta);
				*p_fasta = '\n';
				++p_fasta;
			}
		}
	}

	// write to file remaining parts
	tableFile.write(tableBuffer.data(), p_table - tableBuffer.data());
	tableFile.flush();

	if (fastaFile) {
		fastaFile.write(fastaBuffer.data(), p_fasta - fastaBuffer.data());
		fastaFile.flush();
	}



}