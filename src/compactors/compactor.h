#pragma once
#include <cstdlib>
#include <deque>
#include <stdexcept>
#include <ostream>

#include "kmer.h"

#include <limits>

class Compactor {
	
public:

	const int8_t k{ -1 };
	const int8_t num_kmers{ -1 };
	const int16_t num_kmers_total{ -1 };
	const int32_t total_support{ -1 };
	const int32_t exact_support{ -1 };
	
	int16_t descendants_countdown{ -1 };
	int16_t num_children{ 0 };
	
	// Indicates if compactor has children
	Compactor* parent{ nullptr };
	Compactor* ancestor{ nullptr };
	const kmer_t* kmers{ nullptr };
	double extender_specificity{ -1.0 };

	Compactor(kmer_t* kmer, int8_t k, int16_t descendants_countdown)
		: 
		k{ k }, 
		num_kmers{ 1 }, 
		num_kmers_total{ 1 }, 
		descendants_countdown{ descendants_countdown },
		kmers{ kmer },
		ancestor{ this }
	{ }

	Compactor(Compactor& parent, kmer_t *kmers, int8_t k, int8_t num_kmers, int32_t total_support, int32_t exact_support) 
		: 
		k{ k }, 
		num_kmers { num_kmers }, 
		num_kmers_total{ (int16_t)(parent.num_kmers_total + num_kmers) }, 
		total_support(total_support),
		exact_support(exact_support),
		parent{ &parent },
		ancestor{ parent.ancestor },
		kmers{ kmers }
	{
		if (parent.num_children) {
			// decrase counter only when second children
			--(ancestor->descendants_countdown);
		}
		
		++parent.num_children;
	}

	Compactor(const Compactor& rhs) = delete;
	Compactor(Compactor&& rhs) = default;

	kmer_t get_extender(int len) {
		// fast path when follower and extender are of same length
		if (len == k) {
			return kmers[num_kmers - 1];
		}
		
		kmer_t out = 0;
		int collected = 0;

		Compactor* compactor = this;
		int kmer_index = this->num_kmers - 1;

		do {
		
			// if we need to switch to parent compactor
			if (kmer_index < 0) {
				compactor = compactor->parent;
				kmer_index = compactor->num_kmers - 1;
			}

			kmer_t aux = compactor->kmers[kmer_index];
			out |= (aux << collected * 2);
			
			collected += compactor->k;
			--kmer_index;

		} while (collected < len);

		// trail left-most bits
		int to_trail = 64 - 2 * len;

		out <<= to_trail;
		out >>= to_trail;

		return out;
	}

	int16_t get_descendants_left() const {
		return ancestor->descendants_countdown;
	}
	
	char* to_string(char* dst) const {
		
		if (parent) {
			dst = parent->to_string(dst);
		} 

		for (int i = 0; i < num_kmers; ++i) {
			KmerHelper::to_string(kmers[i], dst, k);
			dst += k;
		}
		*dst = 0;
		return dst;
	}

};
