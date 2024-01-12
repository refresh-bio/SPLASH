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
	const std::string& anchorFile,
	int homopolymerThreshold,
	int numThreads,
	int readsBufferGb,
	size_t anchorsBatchSize,
	bool keepTemp)
	: polyFilter(homopolymerThreshold), anchorFile(anchorFile), anchorsBatchSize(anchorsBatchSize) {

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

	selector.set_max_memory_usage((size_t)readsBufferGb << 30);		// I/O buffer
	selector.set_max_read_len(1 << 20);				// optional 
	selector.set_no_threads(numThreads);					
	selector.set_keep_temps(keepTemp);
	
	fs::path tmp_dir(TEMP_PATH);
	bool ok = false;
	if (!fs::exists(tmp_dir)) {
		ok = fs::create_directory(tmp_dir);
	}
	else {
		ok = fs::is_directory(tmp_dir);
	}
	
	if (!ok) {
		throw runtime_error("Unable to create/open temporary directory: " + TEMP_PATH);
	}

	selector.set_output_dir(tmp_dir.string());
}

ReadLoader::~ReadLoader() {
	if (!selector.get_keep_temps()) {
		fs::remove_all(fs::path(TEMP_PATH));
	}
}

size_t ReadLoader::extractKmers(
	int numFollowers,
	int followerLen,
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
		auto tmp = selector.process_anchor_followers(ReadSelector::direction_t::forward, anchorLen, followerLen, numFollowers);
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
	std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
	bool reverse) {

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
		auto tmp = selector.process_anchor_followers(dir, queryLen, followerLen, numFollowers);
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

	for (const auto& entry : in) {

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

		if (out_v.size() == 0) {
			out.erase(entry.first);
		}
	}


	LOG_DEBUG << "Finished" << endl
		<< out.size() << " anchors, "
		<< n_passed << " reads (" << n_total << " before filtring)" << endl;

	return n_passed;
}



Output::Output(const std::string& table, const std::string& fasta)
	: tableFile(table, std::ios_base::out), fastaFile(fasta, std::ios_base::out), writerQueue(1024, 1) {
	if (!tableFile) {
		throw std::runtime_error("Unable to open output file: " + table);
	}

	tableBuffer.resize(16 << 20);
	fastaBuffer.resize(16 << 20);

	std::string header{ "anchor\tcompactor\tsupport\texact_support\textender_specificity\tnum_extended\n" };
	tableFile.write(header.c_str(), header.length());

	worker = std::thread([this]() {
		task_t task;
		while (this->writerQueue.pop(task)) {
			this->internal_save(task.first, task.last);
		}
	});

}

void Output::internal_save(std::deque<Compactor>::iterator first, std::deque<Compactor>::iterator last) {

	
	char* p_table = tableBuffer.data();
	char* p_fasta = fastaBuffer.data();

	int n_compactors = 0;

	for (auto it = first; it != last; ++it) {
		const Compactor& c = *it;
		if (c.total_support >= 0) { // omit anchors which have total support set to -1

			size_t bytes_needed = (c.num_kmers_total + 1) * c.k  + 1024; // +1 due to anchor column

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
			p_table = c.to_string(p_table);
			char* seq_end = p_table;

			// copy anchor
			copy_n(seq_begin, c.ancestor->k, anchor_begin);

			// fill rest of columns
			int n_extensions = ((c.num_kmers_total - 1) / c.num_kmers) - 1;
			p_table += sprintf(p_table, "\t%d\t%d\t%f\t%d\n", c.total_support, c.exact_support, c.extender_specificity, n_extensions);


			if (fastaFile) {
				// no more space in buffer - write to file
				if (p_fasta + bytes_needed >= fastaBuffer.data() + fastaBuffer.size()) {
					fastaFile.write(fastaBuffer.data(), p_fasta - fastaBuffer.data());
					p_fasta = fastaBuffer.data();
				}

				p_fasta += sprintf(p_fasta, ">compactor-%d\n", n_compactors);
				p_fasta = copy(seq_begin, seq_end, p_fasta);
				*p_fasta = '\n';
				++p_fasta;
			}

			++n_compactors;
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