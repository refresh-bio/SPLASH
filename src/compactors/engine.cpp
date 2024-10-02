#include "engine.h"
#include "array.h"
#include "kmer.h"

#include <string>
#include <algorithm>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <fstream>
#include <string>
#include <array>
#include <ios>
#include <future>

using namespace std;

bool Engine::operator()() {

	// extract kmers
	LOG_NORMAL << "Querying for anchors..." << endl;
	int numKmers = params.numKmers;

	auto t = std::chrono::high_resolution_clock::now();
	std::unordered_map<kmer_t, std::vector<kmer_t>> anchors2kmers;
	size_t n_hits = kmerProvider->extractKmers(numKmers, params.kmerLen, params.allAnchors, anchors2kmers, params.anchorLen);  // overwrite anchor length
	auto dt = std::chrono::high_resolution_clock::now() - t;

	LOG_NORMAL << n_hits << " reads found for " << anchors2kmers.size() << " anchors in " << chrono::duration<double>(dt).count() << " s" << endl << endl;
	
	if (n_hits == 0) {
		LOG_NORMAL << "Check if your reads have at least " << (params.anchorLen + params.numKmers * params.kmerLen) 
			<< " bases (anchor_len + num_kmers * kmer_len)." << endl;
		return false;
	}

	LOG_NORMAL << "Generating compactors for anchors..." << endl;
	t = std::chrono::high_resolution_clock::now();

	int n_anchors = (int)anchors2kmers.size();
	kmer_t* anchors = compactorKmers.resize_for_additional(n_anchors);

	transform(anchors2kmers.begin(), anchors2kmers.end(), anchors,
		[](const auto& pair)->kmer_t { return pair.first; });
	std::sort(anchors, anchors + n_anchors);

	// add anchor compactors
	for (int i = 0; i < n_anchors; ++i) {
		compactors.emplace_back(&anchors[i], params.anchorLen, params.maxAnchorCompactors);
	}

	// extend compactors
	for (int i = 0; i < n_anchors; ++i) {
		
		if (i % 100 == 0) {
			LOG_NORMAL << "\r" << i + 1 << "/" << anchors2kmers.size() << "..." << std::flush;
		}

		Array<kmer_t> hits(std::move(anchors2kmers[anchors[i]]), numKmers);
		extendCompactors(compactors[i], hits, std::numeric_limits<int>::max(), compactors, compactorKmers);
	}
	dt = std::chrono::high_resolution_clock::now() - t;

	int n_compactors = compactors.size() - n_anchors;

	LOG_NORMAL << "\r" <<  n_compactors << " compactors generated from " << anchors2kmers.size() << " anchors in "
		<< chrono::duration<double>(dt).count() << " s" << endl << endl;

	int i_start = n_anchors;
	int extenderLen = params.anchorLen;

	if (params.useRecursion && kmerProvider->supportsExtension()) {

		LOG_NORMAL << "Extending compactors..." << endl;
		t = std::chrono::high_resolution_clock::now();

		size_t non_extended = compactors.size();
		for (int iter = 1; i_start != compactors.size(); ++iter) {
			
			// verify max compactor length
			if ((compactors.back().num_kmers_total + numKmers) * params.kmerLen > params.maxLen) {
				break;
			}

			unordered_set<kmer_t> forwardQuery;

			int n_to_extend = (int)compactors.size() - i_start;
			LOG_NORMAL << "  Iteration " << iter << ": " << n_to_extend << " compactors" << endl;
			auto t = std::chrono::high_resolution_clock::now();
			// check all offsets
			for (int j = 0; j < params.numExtenders; ++j) {
				int shift = j * params.extendersShift;

				unordered_set<kmer_t> query;

				LOG_NORMAL << "    Shift " << shift << ": ";

				// iterate over compactors added to extract pontential extenders
				for (int i = i_start; i < compactors.size(); ++i) {
					Compactor& c = compactors[i];

					// consider only compactors without extender
					if (c.extender_specificity < params.minExtenderSpecificity) {
						kmer_t extender = c.get_extender(extenderLen, shift);
						query.insert(extender);
					}

					//cout << std::hex << extender << endl;
				}

				// make a query to kmer provider
				anchors2kmers.clear();

				if (!kmerProvider->extractPredecessorWithGap(query, extenderLen, numKmers * params.kmerLen - (extenderLen + shift), extenderLen, anchors2kmers)) {
					// if kmer provider is unable to extract kmers for given query
					break;
				}

				// iterate over compactors again to establish extendors specificity
				int n_extendables = 0;
				size_t n_queries_before = forwardQuery.size();
				for (int i = i_start; i < compactors.size(); ++i) {
					Compactor& c = compactors[i];

					// consider only compactors without extender
					if (c.extender_specificity < params.minExtenderSpecificity) {

						kmer_t extender = c.get_extender(extenderLen, shift);
						kmer_t anchor = c.parent->get_extender(extenderLen, c.parent->extender_shift);

						const std::vector<kmer_t>& predecessors = anchors2kmers[extender];

						// get only leftmost predecessor (at original anchor position)
						int agree = 0;
						int total = 0;
						for (int i = 0; i < predecessors.size(); ++i) {
							//if (!polyFilter.IsPolyACGT(predecessors[i], extenderLen)) {
							++total;
							if (predecessors[i] == anchor) {
								++agree;
							}
							//}
						}
						c.extender_specificity = (total == 0) ? -1.0 : (double)agree / total;
						c.extender_shift = shift;

						// use extendor as a new anchor
						if (c.extender_specificity >= params.minExtenderSpecificity) {
							++n_extendables;
							forwardQuery.insert(extender);
						}
					}
				}

				LOG_NORMAL << (forwardQuery.size() - n_queries_before) << " unique seeds (" << n_extendables << " occurences)" << endl;

				if (n_extendables == n_to_extend) {
					break;
				}
			}

			LOG_NORMAL << "    Trying to extend (" << forwardQuery.size() << " total extenders) ";

			// after calculating specificity one can save the portion of compactors 
			compactorWriter->save(compactors.begin(), compactors.begin() + i_start, compactors.end(), false);

			// Try to extend
			int n_last = compactors.size();
			if (forwardQuery.size()) {
				anchors2kmers.clear();
				kmerProvider->extractKmers(forwardQuery, extenderLen, numKmers, params.kmerLen, false, params.allAnchors, anchors2kmers);

				for (int i = i_start; i < n_last; ++i) {
					Compactor& c = compactors[i];
					kmer_t extender = c.get_extender(extenderLen, c.extender_shift);

					if (c.extender_specificity >= params.minExtenderSpecificity) {
						
						if (params.extendAll) {
							Array<kmer_t> hits(anchors2kmers[extender], numKmers);
							extendCompactors(c, hits, params.maxChildCompactors, compactors, compactorKmers);
						}
						else {
							Array<kmer_t> hits(std::move(anchors2kmers[extender]), numKmers);
							extendCompactors(c, hits, params.maxChildCompactors, compactors, compactorKmers);
						}
					}
				}
			}
			dt = std::chrono::high_resolution_clock::now() - t;

			i_start = n_last;

			LOG_NORMAL << compactors.size() - i_start << " extensions (" << chrono::duration<double>(dt).count() << " s)" << endl;
		}

		dt = std::chrono::high_resolution_clock::now() - t;
		LOG_NORMAL << "Extension finished: " << (compactors.size() - non_extended) << " additional compactors generated in "
			<< chrono::duration<double>(dt).count() << " s" << endl << endl;

	}

	// after calculating specificity one can save the portion of compactors 
	compactorWriter->save(compactors.begin(), compactors.begin() + i_start, compactors.end(), true);
	
	return true;
}


bool Engine::extendCompactors(
	Compactor& parent,
	Array<kmer_t>& hits,
	int maxCount,
	deque<Compactor>& compactors,
	chunked_vector<kmer_t>& compactorKmers) {

	struct KmerInfo {
		int n{ 0 };
		int d{ 0 };
		kmer_t closest{ 0 };
	};

	/*
	static std::array<std::array<double, MaxN>, MaxN> expected_factors = [this]() {
		static std::array<std::array<double, MaxN>, MaxN> out;

		for (int n = 0; n < MaxN; ++n) {
			for (int k = 0; k <= n; ++k) {
				out[n][k] = params.beta * (double)this->n_choose_k(n, k) * poisson.pdf(k);
			}
		}

		return out;
	}();
	*/

	if (hits.getHeight() == 0) { 
		return true;
	}

	int numKmers = params.numKmers;
	int kmerLen = params.kmerLen;

	LOG_DEBUG << "ANCHOR " << KmerHelper::to_string(parent.kmers[0], params.kmerLen) << ", reads: " << hits.getHeight() << endl;

	// generate active sets
	vector<unordered_map<kmer_t, KmerInfo>> kmerInfos(numKmers);
	bool kill = false;

	for (int j = 0; j < numKmers; ++j) {
		LOG_DEBUG << "Active sets, column " << j << ": " << endl;
		
		// local vector of kmers sorted wrt abundance
		unordered_map<kmer_t, KmerInfo>& infos = kmerInfos[j];

		// fill infos structure
		for (int i = 0; i < hits.getHeight(); ++i) {
			/*
			// this will also catch empty kmer as a span of Ts
			if (polyFilter.IsPolyACGT(hits[i][j], kmerLen)) {
				hits[i][j] = EMPTY_KMER;
			} else*/

			if (artifactsFilter.ContainsArtifact(hits[i][j], kmerLen)) {
				//cout << "BOO: " << KmerHelper::to_string(hits[i][j], kmerLen) << endl;
				hits[i][j] = EMPTY_KMER;
			}
			
			if (hits[i][j] != EMPTY_KMER) {
				++infos[hits[i][j]].n;
			}
		}

		// break immediately when no proper kmers found in a current column 
		if (infos.size() == 0) {
			kill = true;
			break;
		}

		vector<kmer_t> V(infos.size());
	
		// fill vector of unique kmers and sort wrt abundance
		transform(infos.begin(), infos.end(), V.begin(), [](const pair<kmer_t, KmerInfo>& k2i) {
			return k2i.first;
			});

		std::sort(V.begin(), V.end(), [&infos](const kmer_t a, const kmer_t b) {
			int na = infos[a].n;
			int nb = infos[b].n;
			return (na == nb) ? (a < b) : (na > nb);
			});

		// break immediately when no kmers 
		if (infos[V[0]].n < params.lowerBound) {
			kill = true;
			break;
		}

		vector<pair<kmer_t, int>> S;
		S.reserve(V.size());

		// add the most abundant kmer and make it pointing to itself
		S.emplace_back(V[0], infos[V[0]].n);
		infos[V[0]].closest = V[0];

		LOG_DEBUG << KmerHelper::to_string(V[0], params.kmerLen) << "(n = " << infos[V[0]].n << ")... ACCEPT" << endl;

		// iterate over remaining kmers
		int n_elems = (int)V.size();
		for (int i = 1; i < n_elems; ++i) {
			kmer_t v = V[i];
			auto& v_info = infos[v];
			int min_d = std::numeric_limits<int>::max();
			kmer_t min_kmer = 0;
			double min_factor = 0;

			LOG_DEBUG << KmerHelper::to_string(v, params.kmerLen) << " (n=" << v_info.n << ")... ";

			//editDistance.Prepare(v, kmerLen);

			// iterate over active set elements
			bool add = true;
			for (auto& r : S) {
				//int d = editDistance.Calculate(r);
				int d = KmerHelper::calculateHamming(r.first, v);

				if (d < min_d) {
					min_d = d;
					min_kmer = r.first;
				}

				if (!params.newAcceptanceRule) {
					// OLD method
					double psn = poisson.pdf(d);
					double n_c_k = (double)n_choose_k(kmerLen, d);

					//double exp_count = r.second * expected_factors[kmerLen][d];
					double exp_count = r.second * params.beta * psn * n_c_k;

					if (exp_count >= v_info.n) {
						add = false;
						LOG_DEBUG << "similar found: [" << KmerHelper::to_string(r.first, params.kmerLen) << ", d = " << d
							<< ", n * k_choose_d * poisson * beta = " << r.second << " * " << n_c_k << " * " << psn << " * " << params.beta << " = " << exp_count << "], ";
					}
				}
				else {
					// new method
					//double k_choose_d = (double)n_choose_k(kmerLen, d);
					double bin = binomial.pdf(kmerLen, d);

					double p_val = 1.0 - poissonCdf(r.second * bin, (double)v_info.n); // get probability of observing n or more reads given that they were produced by sequencing error from 
		
					double p_val_corr = params.beta * p_val;

					if (p_val_corr > 0.01) {
						add = false;
						LOG_DEBUG << "similar found: [" << KmerHelper::to_string(r.first, params.kmerLen) 
							<< ", d = " << d
							<< ", pval(lambda=" << r.second << " * " << bin << ", n=" << v_info.n << ") = " << p_val
							<<  ", P_CORR = " << params.beta << " * pval = " << p_val_corr << "], ";
					}
				}
			} 

			if (add) {
				
				if (v_info.n >= params.lowerBound) {
					S.emplace_back(v, v_info.n);
					v_info.closest = v; // point a kmer to itself

					LOG_DEBUG << "ACCEPT" << endl;
				}
				else {
					v_info.closest = min_kmer;
					v_info.d = min_d;

					LOG_DEBUG << "n < lower_bound: REJECT" << endl;
				}
			}
			else {
				v_info.closest = min_kmer;
				v_info.d = min_d;

				LOG_DEBUG << "REJECT" << endl;
			}
		}
	}

	if (kill) {
		LOG_DEBUG << endl;
		return false;
	}

	// generate compactor candidates
	KmersEqual equalComparer(numKmers);
	KmersHasher hasher(numKmers);

	struct Counters {
		int n_exact{ 0 };
		int n_total{ 0 };
	};

	std::unordered_map<kmer_t*, Counters, decltype(hasher), decltype(equalComparer)> counters(8, hasher, equalComparer);
	
	for (int i = 0; i < hits.getHeight(); ++i) {
		int d_sum = 0;
		for (int j = 0; j < numKmers; ++j) {
			kmer_t& kmer = hits[i][j];

			if (kmer == EMPTY_KMER) {
				d_sum = params.maxMismatch + 1;
				break;
			}

			auto& infos = kmerInfos[j];
			
			// replace original kmers with representatives
			d_sum += infos[kmer].d;
			kmer = infos[kmer].closest;
		}

		if (d_sum <= params.maxMismatch) {
			if (d_sum == 0) {
				++counters[hits[i]].n_exact;
			}
			
			++counters[hits[i]].n_total;
		}

	}

	//LOG_DEBUG << "# candidates: " << counters.size();

	if (counters.size() == 0) {
		LOG_DEBUG << endl;
		return false;
	}
	
	KmersLess lessComparer(numKmers);

	vector<pair<kmer_t*, Counters>> candidates(counters.begin(), counters.end());
	std::sort(candidates.begin(), candidates.end(), [&lessComparer](const pair<kmer_t*, Counters>& a, const pair<kmer_t*, Counters>& b) {
		return (a.second.n_total == b.second.n_total) ? lessComparer(a.first, b.first) : (a.second.n_total > b.second.n_total);
	});

	// repeat active set generation procedure to extract final compactors
	vector<pair<kmer_t*, Counters>> output;
	output.reserve(candidates.size());
	output.emplace_back(candidates[0]); // first candidate accepted in spite of containing homopolymers

	int n_candidates = (int)candidates.size();

	int max_in_output = std::min((int)parent.get_descendants_left(), maxCount);

	for (int i = 1; i < n_candidates; ++i) {
		// iterate over active set elements
		auto& v = candidates[i];
		bool add = true;

		for (auto& r : output) {
			int d = KmerHelper::calculateHamming(r.first, v.first, numKmers);
			double factor = n_choose_k(numKmers * kmerLen, d) * poisson.pdf(d) * params.beta;
			//double factor = expected_factors[kmerLen * numKmers][d];
			if (r.second.n_total * factor >= v.second.n_total) {
				add = false;
			}
		}

		if (add) {
			// this also catches empty kmers as a span of Ts
			bool is_poly = any_of(v.first, v.first + numKmers, [this, kmerLen](kmer_t kmer) {
				return this->polyFilter.IsPolyACGT(kmer, kmerLen); 
			});
			
			// other compactors cannot contain homopolymers
			if (!is_poly) {
				output.emplace_back(v);

				if ((int)output.size() >= max_in_output)
					break;
			}
		}
	}

	// iterate over final compactors
	int n_max_extensions = std::min((int)parent.get_descendants_left(), maxCount);
	int n_compactors = std::min(n_max_extensions, (int)output.size());

	kmer_t* raw = compactorKmers.resize_for_additional(n_compactors * numKmers);
	
	for (int i = 0; i < n_compactors; ++i) {
		auto& o = output[i];
		std::copy_n(o.first, numKmers, raw);
		compactors.emplace_back(parent, raw, (int8_t)kmerLen, (int8_t)numKmers, o.second.n_total, o.second.n_exact);
		raw += numKmers;
	}

	//LOG_DEBUG << "# output: " << n_compactors << endl;

	return true;
}
