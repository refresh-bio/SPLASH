#pragma once
#include "kmer.h"
#include "../common/filters/poly_ACGT_filter.h"
#include <refresh/parallel_queues/lib/parallel-queues.h>
#include "read_select.h"
#include "engine.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <thread>


class ReadLoader : public IKmerProvider {
	
	const PolyACGTFilter polyFilter;

	ReadSelector selector;

	std::string anchorFile;

	size_t anchorsBatchSize;
	
	std::string tempPath;

public: 

	ReadLoader(
		const std::vector<std::string>& fastqFiles,
		input_format_t inputFormat,
		const std::string& anchorFile,
		int homopolymerThreshold,
		int numThreads,
		int readsBufferGb,
		size_t anchorsBatchSize,
		bool keepTemp,
		std::string tempPath = "./tmp-compactors");

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
		std::deque<Compactor>::iterator begin;
		std::deque<Compactor>::iterator end;
	};

	std::ofstream tableFile;
	std::ofstream fastaFile;

	std::vector<char> tableBuffer;
	std::vector<char> fastaBuffer;

	std::thread worker;

	refresh::parallel_queue<task_t> writerQueue;

	bool noSubcompactors;
	bool cumulatedStats;

	int numCompactors{ 0 };

public:

	Output(const std::string& table, const std::string& fasta, bool noSubcompactors, bool cumulatedStats);

	void save(
		std::deque<Compactor>::iterator begin,
		std::deque<Compactor>::iterator current,
		std::deque<Compactor>::iterator end, 
		bool lastPortion) override {
		
		if (lastPortion) {
			
			if (noSubcompactors) {
				writerQueue.push(task_t{ begin, end });
			}
			else {
				writerQueue.push(task_t{ current, end });
			}

			writerQueue.mark_completed();
			worker.join();
		}
		else {
			if (!noSubcompactors) {
				writerQueue.push(task_t{ current, end });
			}
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
