#ifndef _WALK_KMERS
#define _WALK_KMERS
#include "../common/types/satc_data.h"
#include <cinttypes>
#include <string>
#include <string_view>

class WalkKmers {
	std::string_view seq;
	uint32_t kmer_len;
	bool canonical;
	uint32_t rc_shift{};
	uint64_t pos{};
	uint64_t cur_kmer{};
	uint64_t mask{};
	static uint8_t map_symb(uint8_t symb) {
		static auto mapping = []() {
			static uint8_t map[256];
			std::fill_n(map, 256, 255);
			map['A'] = map['a'] = 0;
			map['C'] = map['c'] = 1;
			map['G'] = map['g'] = 2;
			map['T'] = map['t'] = 3;
			return map;
		}();
		return mapping[symb];
	}

	bool start_kmer_prefix() {
		cur_kmer = 0;
		uint32_t in_kmer = 0;
		while (in_kmer < kmer_len - 1) {
			if (pos == seq.length())
				return false;
			auto symb = map_symb(seq[pos++]);
			if (symb == 255) {
				in_kmer = 0;
				cur_kmer = 0;
			}
			else {
				cur_kmer <<= 2;
				cur_kmer |= symb;
				++in_kmer;
			}
		}
		return true;
	}

	uint64_t make_mask() {
		if (kmer_len == 32)
			return -1;
		return (1ull << (2 * kmer_len)) - 1;
	}

public:
	WalkKmers(std::string_view seq,
		uint32_t kmer_len,
		bool canonical) :
		seq(seq),
		kmer_len(kmer_len),
		canonical(canonical),
		rc_shift(get_rev_compl_shift(kmer_len)),
		mask(make_mask()) {
		assert(kmer_len);
		assert(kmer_len <= 32);
		if (seq.size() < kmer_len) {
			pos = seq.size();
		}
		start_kmer_prefix();
	}

	bool Next(uint64_t& kmer) {
		while (pos < seq.size()) {
			auto symb = map_symb(seq[pos++]);
			if (symb == 255) {
				if (!start_kmer_prefix())
					return false;
				continue;
			}
			cur_kmer <<= 2;
			cur_kmer += symb;
			cur_kmer &= mask;
			if (canonical)
			{
				auto rc = get_rev_compl(cur_kmer, rc_shift);
				kmer = cur_kmer < rc ? cur_kmer : rc;
			}
			else
				kmer = cur_kmer;
			return true;
		}
		return false;
	}
};

//same as WalKmers but also walks invalid (with non ACGT) k-mers
class WalkKmersWithInvalid {
	std::string_view seq;
	uint32_t kmer_len;
	bool canonical;
	uint32_t rc_shift{};
	uint64_t pos{};
	uint64_t cur_kmer{};
	uint64_t mask{};
	uint32_t skip_next_n{};
	static uint8_t map_symb(uint8_t symb) {
		static auto mapping = []() {
			static uint8_t map[256];
			std::fill_n(map, 256, 255);
			map['A'] = map['a'] = 0;
			map['C'] = map['c'] = 1;
			map['G'] = map['g'] = 2;
			map['T'] = map['t'] = 3;
			return map;
		}();
		return mapping[symb];
	}

	bool start_kmer() {
		cur_kmer = 0;
		uint32_t in_kmer = 0;
		while (in_kmer < kmer_len - 1) {
			if (pos == seq.length())
				return false;
			auto symb = map_symb(seq[pos++]);
			if (symb == 255) {
				skip_next_n = in_kmer + 1;
				in_kmer = 0;
				cur_kmer = 0;
				return false;
			}
			else {
				cur_kmer <<= 2;
				cur_kmer |= symb;
				++in_kmer;
			}
		}
		return true;
	}

	bool start_kmer_prefix() {
		cur_kmer = 0;
		uint32_t in_kmer = 0;
		while (in_kmer < kmer_len - 1) {
			if (pos == seq.length())
				return false;
			auto symb = map_symb(seq[pos++]);
			if (symb == 255) {
				skip_next_n = in_kmer;
				in_kmer = 0;
				cur_kmer = 0;
			}
			else {
				cur_kmer <<= 2;
				cur_kmer |= symb;
				++in_kmer;
			}
		}
		return true;
	}

	uint64_t make_mask() {
		if (kmer_len == 32)
			return -1;
		return (1ull << (2 * kmer_len)) - 1;
	}

public:
	WalkKmersWithInvalid(std::string_view seq,
		uint32_t kmer_len,
		bool canonical) :
		seq(seq),
		kmer_len(kmer_len),
		canonical(canonical),
		rc_shift(get_rev_compl_shift(kmer_len)),
		mask(make_mask()) {
		assert(kmer_len);
		assert(kmer_len <= 32);
		if (seq.size() < kmer_len) {
			pos = seq.size();
		}
		//start_kmer_prefix();
		start_kmer();
	}

	bool Next(uint64_t& kmer, bool& is_valid) {
		is_valid = true;
		while (pos < seq.size()) {
			if (skip_next_n) {
				--skip_next_n;
				if (skip_next_n == 0)
					start_kmer();
				is_valid = false;
				return true;
			}
			auto symb = map_symb(seq[pos++]);
			if (symb == 255) {
				skip_next_n = kmer_len - 1;
				is_valid = false;
				return true;
			}
			cur_kmer <<= 2;
			cur_kmer += symb;
			cur_kmer &= mask;
			if (canonical)
			{
				auto rc = get_rev_compl(cur_kmer, rc_shift);
				kmer = cur_kmer < rc ? cur_kmer : rc;
			}
			else
				kmer = cur_kmer;
			return true;
		}
		return false;
	}
};

#endif // !_WALK_KMERS
