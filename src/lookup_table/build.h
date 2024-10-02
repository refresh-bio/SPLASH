#ifndef _BUILD_H
#define _BUILD_H

#include <string>
#include <cinttypes>
#include <thread>
#include "kmer_index.h"

namespace build_mode {
	struct Params
	{
		uint64_t poly_ACGT_len{};
		std::string input;
		uint32_t category_3_threshold = 1;
		size_t n_threads = std::thread::hardware_concurrency() > 8 ? 8 : std::thread::hardware_concurrency();
		std::string tmp_dir = ".";
		std::string output;
		kmer_index_variant_t kmer_index_variant = kmer_index_variant_t::sbwt;

		//for kmer_index_variant_t::sbwt
		//instead of computing sbwt use precomputed one
		std::string precomputed_sbwt;

		void Print(std::ostream& oss) const;
		
		static void Usage(char* prog_name);
	};

	Params read_params(int argc, char** argv);

	void run(const Params& params);
}

#endif // !_BUILD_H

