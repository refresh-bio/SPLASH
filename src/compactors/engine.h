#pragma once
#include "distributions.h"
#include "compactor.h"
#include "kmer.h"
#include "array.h"
#include "log.h"

#include "chunked_vector.h"
#include "../common/filters/poly_ACGT_filter.h"
#include "../common/edit_distance.h"
#include "../common/filters/artifacts_filter.h"
#include "../common/filters/illumina_adapters_static.h"
#include "../common/types/common_types.h"
#include <refresh/parallel_queues/lib/parallel-queues.h>

#include <thread>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <chrono>
#include <atomic>


class IKmerProvider {
public:
	virtual size_t extractKmers(
		int numFollowers,
		int followerLen,
		bool allAnchors,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
		int& anchorLen) = 0;

	virtual size_t extractKmers(
		const std::unordered_set<kmer_t>& queries,
		int queryLen,
		int numFollowers,
		int followerLen,
		bool reverse,
		bool allAnchors,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets) = 0;

	virtual size_t extractPredecessorWithGap(
		const std::unordered_set<kmer_t>& queries,
		int queryLen,
		int gapLen,
		int predecessorLen,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets) = 0;

	virtual bool supportsExtension() const { return false; }

	virtual ~IKmerProvider() {}
};


class ICompactorWriter {
public:

	virtual void save(
		std::deque<Compactor>::iterator begin, 
		std::deque<Compactor>::iterator current, 
		std::deque<Compactor>::iterator end, 
		bool lastPortion) {}

	virtual ~ICompactorWriter() {}
};


class Engine {

public:
	struct Params {
		int anchorLen{ 27 }; // this will be overwriten by the contents of the anchor table
		int kmerLen{ 27 };
		int numKmers{ 2 };
		bool allAnchors{ false };

		double epsilon{ 0.05 };
		double beta{ 5 };
		int lowerBound{ 10 };
		int maxMismatch{ 4 };
		int polyThreshold{ 6 };

		bool useEditDistance{ false };
		
		bool useRecursion{ false };
		int maxLen{ 2000 };
		double minExtenderSpecificity{ 0.9 };
		int numExtenders{ 1 };
		int extendersShift{ 1 };

		int maxAnchorCompactors{ 1000 };
		int maxChildCompactors{ 20 };
		bool extendAll{ false };
		bool newAcceptanceRule{ false };
	};

protected:

	static const int MAX_READ_LENGTH{ 1 << 20 };

	const NChooseK_log<MAX_READ_LENGTH> n_choose_k;
	const Poisson<MAX_READ_LENGTH> poisson;
	const Binomial<MAX_READ_LENGTH> binomial; 
	const PolyACGTFilter polyFilter;
	ArtifactsFilter artifactsFilter;

	CEditDistanceOneWord editDistance;

	Params params;

	std::shared_ptr<IKmerProvider> kmerProvider;
	std::shared_ptr<ICompactorWriter> compactorWriter;

	std::thread loader;
	std::vector<std::thread> workers;

	std::deque<Compactor> compactors;
	chunked_vector<kmer_t> compactorKmers{ 1000000 };

	

public:

	const Params& getParams() const { return params; }

	int getCompactorLen() const { return params.anchorLen + params.numKmers * params.kmerLen; }


	Engine(
		std::shared_ptr<IKmerProvider> kmerReader,
		std::shared_ptr<ICompactorWriter> writer,
		const Params& params)  
		: 
		kmerProvider(kmerReader),
		compactorWriter(writer),
		binomial(params.epsilon),
		poisson(params.epsilon), 
		polyFilter(params.polyThreshold),
		params(params)
		{
			artifactsFilter.Add(12, IlluminaAdaptersStatic::Get12Mers());
		}

	bool operator()();

	const std::deque<Compactor>& getOutputCompactors() const { return compactors; }

protected:

	bool extendCompactors(
		Compactor& parent,
		Array<kmer_t>& hits, 
		int maxCount,
		std::deque<Compactor>& compactors,
		chunked_vector<kmer_t>& compactorKmers);

};

