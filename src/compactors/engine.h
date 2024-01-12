#pragma once
#include "distributions.h"
#include "compactor.h"
#include "kmer.h"
#include "array.h"
#include "log.h"

#include "chunked_vector.h"
#include "../common/poly_ACGT_filter.h"
#include "../common/edit_distance.h"
#include "../common/artifacts_filter.h"
#include "../common/illumina_adapters_static.h"

#include "../../libs/refresh/parallel-queues.h"

#include <thread>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <chrono>


class IKmerProvider {
public:
	virtual size_t extractKmers(
		int numFollowers,
		int followerLen,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
		int& anchorLen) = 0;

	virtual size_t extractKmers(
		const std::unordered_set<kmer_t>& queries,
		int queryLen,
		int numFollowers,
		int followerLen,
		std::unordered_map<kmer_t, std::vector<kmer_t>>& targets,
		bool reverse) = 0;

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
	virtual void save(std::deque<Compactor>::iterator first, std::deque<Compactor>::iterator last, bool lastPortion) {}

	virtual ~ICompactorWriter() {}
};


class Engine {

public:
	struct Params {
		int anchorLen{ 27 }; // this will be overwriten by the contents of the anchor table
		int kmerLen{ 27 };
		int numKmers{ 2 };

		double epsilon{ 0.05 };
		double beta{ 5 };
		int lowerBound{ 10 };
		int maxMismatch{ 4 };
		int polyThreshold{ 6 };

		bool useEditDistance{ false };
		
		bool useRecursion{ false };
		int maxLen{ 2000 };
		double minExtenderSpecificity{ 0.9 };
		int maxAnchorCompactors{ 1000 };
		int maxChildCompactors{ 20 };
		bool extendAll{ false };
		bool newAcceptanceRule{ false };
	};

protected:

	static const int MaxN{ 128 };
	static const int MaxPoisson{ 1 << 20 };

	const NChooseK<uint64_t, MaxN> n_choose_k;
	const Binomial<MaxN> binomial;
	const Poisson<MaxPoisson> poisson;
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

