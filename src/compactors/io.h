#pragma once
#include "kmer.h"
#include "../common/poly_ACGT_filter.h"
#include "../../libs/refresh/parallel-queues.h"
#include "read_select.h"
#include "engine.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <thread>


class ReadLoader : public IKmerProvider {

	
	const string TEMP_PATH{ "./tmp-compactors" };
	
	const PolyACGTFilter polyFilter;

	ReadSelector selector;

	std::string anchorFile;

	size_t anchorsBatchSize;

public: 

	ReadLoader(
		const std::vector<std::string>& fastqFiles,
		input_format_t inputFormat,
		const std::string& anchorFile,
		int homopolymerThreshold,
		int numThreads,
		int readsBufferGb,
		size_t anchorsBatchSize,
		bool keepTemp);

	~ReadLoader();

	size_t extractKmers(
		int numFollowers,
		int followerLen,
		bool allAnchors,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
		int& anchorLen) override;

	size_t extractKmers(
		const std::unordered_set<kmer_t>& queries,
		int queryLen,
		int numFollowers,
		int followerLen,
		bool reverse,
		bool allAnchors,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets	
	) override;

	size_t extractPredecessorWithGap(
		const std::unordered_set<kmer_t>& queries,
		int queryLen,
		int gapLen,
		int predecessorLen,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets) override;

	virtual bool supportsExtension() const { return true; }

protected:

	size_t translate(
		ReadSelector::dict_t& in,
		int numKmers,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& out
	) const;
	
};


class Output : public ICompactorWriter {

	struct task_t {
		std::deque<Compactor>::iterator first;
		std::deque<Compactor>::iterator last;
	};

	std::ofstream tableFile;
	std::ofstream fastaFile;

	std::vector<char> tableBuffer;
	std::vector<char> fastaBuffer;

	std::thread worker;

	refresh::parallel_queue<task_t> writerQueue;


public:
	Output(const std::string& table, const std::string& fasta);

	void save(std::deque<Compactor>::iterator first, std::deque<Compactor>::iterator last, bool lastPortion) override {
		writerQueue.push(task_t{ first, last });

		//internal_save(first, last);

		if (lastPortion) {
			writerQueue.mark_completed();
			worker.join();
		}
	}

	~Output() {
		if (!writerQueue.check_completed()) {
			writerQueue.mark_completed();
			worker.join();
		}
	}

protected:
	void internal_save(std::deque<Compactor>::iterator first, std::deque<Compactor>::iterator last);

};
