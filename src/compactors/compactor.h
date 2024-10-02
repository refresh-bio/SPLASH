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
	int32_t extender_shift{ 0 };
	int32_t id{ -1 };
	

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

	kmer_t get_extender(int len, int shift_from_end = 0) {
		// fast path when follower and extender are of same length
		int kmer_index = this->num_kmers - 1 - shift_from_end / k;
		int kmer_offset = shift_from_end % k;

		if ((len == k) && (kmer_offset == 0)) {
			return kmers[kmer_index];
		}
		
		kmer_t out = 0;
		int collected = 0;
		Compactor* compactor = this;
		
		do {
		
			// if we need to switch to parent compactor
			if (kmer_index < 0) {
				compactor = compactor->parent;
				kmer_index = compactor->num_kmers - 1;
			}

			kmer_t aux = compactor->kmers[kmer_index] >> (kmer_offset * 2);
			out |= (aux << collected * 2);
			
			collected += compactor->k - kmer_offset;
			--kmer_index;
			kmer_offset = 0;

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

	int32_t get_parent_id() const {
		return parent == nullptr ? -1 : parent->id;
	}

	double calculate_expected_read_count() const {
		return (parent == nullptr || parent->id == -1)
			? exact_support
			: std::min(parent->calculate_expected_read_count(), parent->extender_specificity * exact_support);
	}

	int print_cumulated_id(char* buf) const {
		int out = 0;

		// add parent part and separator
		if (parent != nullptr && parent->id != -1) {
			out += parent->print_cumulated_id(buf);
			buf += out;
			*buf = ',';
			++buf;
			++out;
		}
		
		out += sprintf(buf, "%d", id);

		return out;
	}

	int print_cumulated_exact_support(char* buf) const {
		int out = 0;

		// add parent part and separator
		if (parent != nullptr && parent->id != -1) {
			out += parent->print_cumulated_exact_support(buf);
			buf += out;
			*buf = ',';
			++buf;
			++out;
		}

		out += sprintf(buf, "%d", exact_support);

		return out;
	}

	int print_cumulated_extender_specificity(char* buf) const {
		int out = 0;

		// add parent part and separator
		if (parent != nullptr && parent->id != -1) {
			out += parent->print_cumulated_extender_specificity(buf);
			buf += out;
			*buf = ',';
			++buf;
			++out;
		}

		out += sprintf(buf, "%lf", extender_specificity);

		return out;
	}


	char* to_string(char* dst, bool finishAtExtender) const {
		
		if (parent) {
			dst = parent->to_string(dst, true);
		} 

		// stop printing at extender
		if (finishAtExtender) {

			// process all kmers beside last one
			int last_kmer = this->num_kmers - 1 - extender_shift / k;
			for (int i = 0; i < last_kmer; ++i) {
				KmerHelper::to_string(kmers[i], dst, k);
				dst += k;
			}

			// process last one
			int offset = extender_shift % k;
			kmer_t kmer = kmers[last_kmer];
			kmer >>= offset * 2;

			KmerHelper::to_string(kmer, dst, k - offset);
			dst += k - offset;
		} 
		else {
			for (int i = 0; i < num_kmers; ++i) {
				KmerHelper::to_string(kmers[i], dst, k);
				dst += k;
			}
		}

		*dst = 0;
		return dst;
	}

	

};
