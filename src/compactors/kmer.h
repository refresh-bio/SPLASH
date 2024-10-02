#pragma once
#include <cstdint>
#include <functional>
#include <string>

#define EMPTY_KMER ~0ULL

using kmer_t = uint64_t;

class KmersHasher {
	int count;
	std::hash<kmer_t> hasher;
public:
	KmersHasher(int count) : count(count) {}

	size_t operator()(kmer_t* kmers) const {
		size_t hs = hasher(kmers[0]);
		for (int i = 1; i < count; ++i) {
			hs ^= hasher(kmers[i]);
		}
		return hs;
	}
};

class KmersEqual {
	int count;

public:
	KmersEqual(int count) : count(count) {}

	bool operator()(const kmer_t* va, const kmer_t* vb) const {
		bool ok = true;
		for (int i = 0; i < count; ++i) {
			ok &= (va[i] == vb[i]);
		}
		return ok;
	}
};

class KmersLess {
	int count;

public:
	KmersLess(int count) : count(count) {}

	bool operator()(const kmer_t* va, const kmer_t* vb) const {
		for (int i = 0; i < count; ++i, ++va, ++vb) {
			if (*va == *vb) {
				continue;
			} else {
				return *va < *vb;
			} 
		}
		return false;
	}
};


class KmerHelper {
public:
	static const int HOMOPOLYMER_LEN = 5;
	static const kmer_t HOMOPOLYMER_MASK = (kmer_t)1ULL << 63;

	template <class T>
	static T base2num(T b) {
		return b == 'A' ? 0 : (b == 'C' ? 1 : (b == 'G' ? 2 : (b == 'T' ? 3 : -1)));

	}

	template <class T>
	static T num2base(T i) {
		return i == 0 ? 'A' : (i == 1 ? 'C' : (i == 2 ? 'G' : (i == 3 ? 'T' : 'N')));

	}

	static int from_string(const char* p, int k, kmer_t& kmer) {
		kmer = 0;
		// bool ok = true;
		const char* beg = p;
		const char* end = p + k;

		for (; p < end; ++p) {
			char symbol = base2num(*p);

			if (symbol < 0) {
				break;
			}

			kmer = (kmer << 2) + symbol;
	
		}

		return p - beg;
	}

	/*
	static bool from_string(const char* p, int k, kmer_t& kmer) {
		kmer = 0;
		bool ok = true;
		const char* end = p + k;

		char prevSymbol = -1;
		int runCounter = 1;
		bool homopolymer = false;

		for (; p < end; ++p) {
			char symbol = base2num(*p);

			if (symbol < 0) {
				ok = false;
				break;
			}

			if (symbol == prevSymbol) {
				++runCounter;
				if (runCounter >= HOMOPOLYMER_LEN) {
					homopolymer = true;
				}
			}
			else {
				runCounter = 1;
			}

			kmer = (kmer << 2) + symbol;
			prevSymbol = symbol;
		}

		if (homopolymer) {
			kmer |= HOMOPOLYMER_MASK;
		}

		return ok;
	}
	*/

	static void to_string(kmer_t kmer, char* p, int k) {
		do {
			p[--k] = (char)num2base(kmer & 0x03);
			kmer >>= 2;
		} while (k);
	}

	static std::string to_string(kmer_t kmer, int k) {
		std::string s;
		s.resize(k);
		to_string(kmer, const_cast<char*>(s.data()), k);
		return s;
	}

	static int calculateHamming(kmer_t a, kmer_t b) {
		
		// generate mask with even bits
		static kmer_t MASK_EVEN_BITS = []() {
			kmer_t v = 0;
			for (int i = 0; i < 32; ++i) {
				v <<= 2;
				v |= 1ULL;
			}
			return v;
		}();

		//a &= ~HOMOPOLYMER_MASK;
		//b &= ~HOMOPOLYMER_MASK;

//		kmer_t diffs_even = (a ^ b) & MASK_EVEN_BITS;
//		kmer_t diffs_odd = ((a >> 1) ^ (b >> 1)) & MASK_EVEN_BITS;
//		kmer_t diffs = diffs_even | diffs_odd;

		kmer_t x = a ^ b;
		x |= x >> 1;
		kmer_t diffs = x & MASK_EVEN_BITS;

		#ifdef _WIN32
		int cnt = __popcnt64(diffs);
		#else
		int cnt = __builtin_popcountll(diffs);
		#endif
		return cnt;
	}


	static int calculateHamming(kmer_t* a, kmer_t* b, int num) {
		int d = 0;
		for (int i = 0; i < num; ++i) {
			d += calculateHamming(a[i], b[i]);
		}
		return d;
	}
};


