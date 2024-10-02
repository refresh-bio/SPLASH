#include "build.h"
#include "lookup.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <variant>

#include "../common/filters/poly_ACGT_filter.h"
#include "../common/kmc_api/kmc_file.h"
#include "../common/version.h"
#include "transcriptome_kmer_finder.h"
#include "readers.h"
#include "walk_kmers.h"
#include "seq_id_mapper.h"
#include "utils.h"

namespace build_mode {
	void Params::Print(std::ostream& oss) const
	{
		oss << "input                          : " << input << "\n";
		oss << "output                         : " << output << "\n";

		oss << "n_threads                      : " << n_threads << "\n";
		oss << "tmp_dir                        : " << tmp_dir << "\n";

		oss << "poly_ACGT_len                  : " << poly_ACGT_len << "\n";
		oss << "category_3_threshold           : " << category_3_threshold << "\n";

		//oss << "variant                        : " << to_string(kmer_index_variant) << "\n";

		oss << "precomputed_sbwt               : " << precomputed_sbwt << "\n";
	}

	void Params::Usage(char* prog_name) {
		std::cerr << "lookup_table (build lookup table)\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " build [options] <input> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input>              - path to the file containing fasta list and its precomputed kmc databases, the line format is (<fasta_path>;<all_kmers_kmc_db>;<unique_to_this_fasta_kmc_db>)\n"
			<< "    <output>             - path to the output file\n"
			;
		std::cerr
			<< "Options:\n"
			<< "    --n_threads <int>  - number of threads\n"
			<< "    --tmp_dir <string> - path to directory for temporary files\n"

			<< "    --poly_ACGT_len <int> - all k-mers containing polyACGT of this length will be filtered out (0 means no filtering)\n"
			<< "    --category_3_threshold <int> - accept k-mer in category 3 if its present in a given file <=category_3_threshold times (default: 1)\n"
			//<< "    --variant <sbwt> - k-mer index type (sbwt is current recommendation)\n"

			;

		std::cerr
			//<< "Options specific for variant=sbwt:\n"
			<< "    --precomputed_sbwt <string> - path to precomputed sbwt index (if set lookup_table will use it instead of building own sbwt - must be build for exactly the same set of input files!)\n"
			;

	}

	Params read_params(int argc, char** argv)
	{
		Params res;
		if (argc == 2) {
			Params::Usage(argv[0]);
			exit(0);
		}
		int i = 2;
		for (; i < argc; ++i)
		{
			if (argv[i][0] != '-')
				break;

			std::string param = argv[i];

			if (param == "--n_threads") {
				std::string tmp = argv[++i];
				res.n_threads = std::stoull(tmp);
			}
			else if(param == "--tmp_dir") {
				res.tmp_dir = argv[++i];
			}
			else if (param == "--poly_ACGT_len") {
				std::string tmp = argv[++i];
				res.poly_ACGT_len = std::stoull(tmp);
			}
			else if (param == "--category_3_threshold") {
				std::string tmp = argv[++i];
				res.category_3_threshold = std::stoull(tmp);
			}
			else if (param == "--precomputed_sbwt")
			{
				res.precomputed_sbwt = argv[++i];
			}
			else if(param == "--variant")
			{
				res.kmer_index_variant = kmer_index_variant_from_string(argv[++i]);
			}
		}

		if (i >= argc) {
			std::cerr << "Error: input missing\n";
			exit(1);
		}
		res.input = argv[i++];

		if (i >= argc) {
			std::cerr << "Error: output missing\n";
			exit(1);
		}
		res.output = argv[i++];

		return res;
	}

	std::vector<std::string> extract_fasta_paths(const std::vector<std::string>& lines) {
		std::vector<std::string> res;
		res.reserve(lines.size());

		for (const auto& line : lines) {
			std::istringstream iss(line);
			res.emplace_back();
			if (!std::getline(iss, res.back(), ';')) {
				std::cerr << "Error: cannot parse line:" << line << "\n";
			}
		}
		return res;
	}

	void split_input_line(const std::string& line, std::string& fasta_path, std::string& all_kmers_kmc_db_path, std::string& uniq_kmers_kmc_db_path, bool& is_transcriptome) {
		std::istringstream iss(line);
		if (!std::getline(iss, fasta_path, ';')) {
			std::cerr << "Error: cannot parse line:" << line << "\n";
		}
		if (!std::getline(iss, all_kmers_kmc_db_path, ';')) {
			std::cerr << "Error: cannot parse line:" << line << "\n";
		}
		if (!std::getline(iss, uniq_kmers_kmc_db_path, ';')) {
			std::cerr << "Error: cannot parse line:" << line << "\n";
		}
		std::string tmp;
		if (!std::getline(iss, tmp, ';')) {
			std::cerr << "Error: cannot parse line:" << line << "\n";
		}
		else {
			if (tmp == "True")
				is_transcriptome = true;
			else if (tmp == "False")
				is_transcriptome = false;
			else
				std::cerr << "Error: expected 'True' or 'False' (info if this is transcriptome) but get :" << tmp <<" in line " << line << "\n";
		}
	}

	std::vector<std::string> get_input_lines(const std::string& file_list) {
		std::ifstream input(file_list);
		if (!input) {
			std::cerr << "Error: cannot open file " << file_list << "\n";
			exit(1);
		}

		std::vector<std::string> res;
		std::string line;

		while (std::getline(input, line)) {
			if (line == "") {
				continue;
			}
			res.emplace_back(line);
		}
		return res;
	}

	namespace parallel {
		template<typename BUILDER_T>
		class CategoryDeterminer
		{
			refresh::thread_control& tc;
			CKMCFile* kmc_db_all;
			CKMCFile* kmc_db_unique;
			CKmerAPI kmer_api;
			std::vector<uint64_t> kmer_in_vec;
			refresh::parallel_queue<std::vector<task_desc>>& build_queue;
			uint32_t category_3_threshold;
			BUILDER_T& builder;

			uint32_t query_kmer(CKMCFile& kmc_db, uint64_t kmer) {
				kmer_in_vec[0] = kmer;
				kmer_api.from_long(kmer_in_vec);
				uint32_t res = 0;
				if (!kmc_db.CheckKmer(kmer_api, res))
					return 0;
				return res;
			};

			void determine_category(task_desc& pack) {
				pack.cnt = query_kmer(*kmc_db_unique, pack.kmer);
				pack.cat = -1; //non of 1, 2, 3
				if (pack.cnt) {
					if (pack.cnt > 1) {
						//mkokot_TODO: to wonder, discuss, change in documentation
						//this is the case that some k-mer ocurss multiple times in a single file, but it is pretty rare in this file, so maybe better to store each its ocurrence
						//if (pack.cnt <= category_3_threshold) {
						//  pack.cat = 3;
						//}
						//else {
						pack.cat = 2;

						//}
					}
					else {
						pack.cat = 1;
					}
				}
				else {
					uint32_t counter_all = query_kmer(*kmc_db_all, pack.kmer);
					if (counter_all <= category_3_threshold) {
						pack.cat = 3;
					}
				}
			}

		public:
			CategoryDeterminer(refresh::thread_control& tc, CKMCFile* kmc_db_all, CKMCFile* kmc_db_unique, uint32_t kmer_len,
			                   refresh::parallel_queue<std::vector<task_desc>>& build_queue,
			                   uint32_t category_3_threshold, BUILDER_T& builder) :
				tc(tc),
				kmc_db_all(kmc_db_all),
				kmc_db_unique(kmc_db_unique),
				kmer_api(kmer_len),
				kmer_in_vec(1),
				build_queue(build_queue),
				category_3_threshold(category_3_threshold),
				builder(builder)
			{
				std::vector<task_desc> pack;

				while (build_queue.pop(pack)) {
					tc.execute([&] {
						for (auto& x : pack)
							determine_category(x);

						builder.AddPack(std::move(pack));
					});
				}
			}
		};

		class ParallelBuildHelper {
			const size_t pack_size = 100'000;
			std::vector<task_desc> data;

			std::vector<std::thread> threads;

			refresh::parallel_queue<std::vector<task_desc>> build_queue;

			void add_pack_if_needed() {
				if (data.size() == pack_size) {
					build_queue.push(std::move(data));
					data.reserve(pack_size);
				}
			}

			void add_last_pack() {
				if (data.size()) {
					build_queue.push(std::move(data));
				}
			}
		public:
			template<typename BUILDER_T>
			ParallelBuildHelper(CKMCFile* kmc_db_all, CKMCFile* kmc_db_unique, refresh::thread_control* tc, uint32_t kmer_len,
			                    uint32_t category_3_threshold, BUILDER_T* builder) :
				build_queue(tc->max_running() + 1, 1, "build-queue")
			{
				threads.reserve(tc->max_running());
				for (size_t tid = 0; tid < tc->max_running(); ++tid) {
					threads.emplace_back([kmc_db_all, kmc_db_unique, tc, kmer_len, category_3_threshold, builder, this] {
						CategoryDeterminer<BUILDER_T> category_determiner(*tc, kmc_db_all, kmc_db_unique, kmer_len, build_queue, category_3_threshold, *builder);
					});
				}
			}

			void Add(uint64_t kmer, uint32_t file_id, uint32_t header_id) {
				data.emplace_back(kmer, file_id, header_id);
				add_pack_if_needed();
			}

			void Finish() {

				add_last_pack();
				build_queue.mark_completed();

				for (auto& th : threads)
					th.join();
			}
		};

	};
	//mkokot_TODO: consider
	//maybe it would all work better if we send whole sequence
	//it would allow querying SBWT in streaming mode which should work faster
	//besides It may be worth to just replace global kmc database with computed earlier 
	//SBWT... (if possible)
	
	//I am leaving here a BUILDER_T althoug currently it is the only one variant
	//previously I had also hash table variant
	template<typename BUILDER_T>
	void run_impl(const Params& params) {
		static_assert(
			std::is_same_v<BUILDER_T, KmerIndex_SBWTBuilder>,
			"BUILDER_T must be KmerIndex_SBWTBuilder"
			);

		refresh::thread_control tc(params.n_threads);

		PolyACGTFilter poly_ACGT_filter(params.poly_ACGT_len);

		std::vector<std::string> input_lines = get_input_lines(params.input);

		uint32_t n_bits_for_file_id = lookup_table_utils::bits_required_to_represent(input_lines.size() - 1);

		std::cerr << "Info: n_bits_for_file_id = " << n_bits_for_file_id << "\n";

		SeqIdMapper file_mapper;
		HeaderIdMapper header_mapper;

		std::unique_ptr<BUILDER_T> builder;

		std::string fasta_path;
		std::string all_kmers_kmc_db_path;
		std::string uniq_kmers_kmc_db_path;
		bool is_transcriptome{};

		size_t kmer_len = 0; //0 means we don't know yet

		for (const auto& line : input_lines) {

			split_input_line(line, fasta_path, all_kmers_kmc_db_path, uniq_kmers_kmc_db_path, is_transcriptome);

			std::cerr << "Processing " << fasta_path;
			if (is_transcriptome)
				std::cerr << " (transcriptome)";
			std::cerr << "\t" << all_kmers_kmc_db_path << "\t" << uniq_kmers_kmc_db_path << "\n";

			auto file_id = file_mapper.Add(fasta_path);

			std::cerr << "Opening KMC database " << uniq_kmers_kmc_db_path << "\n";
			CKMCFile kmc_file_uniq;
			if (!kmc_file_uniq.OpenForRA(uniq_kmers_kmc_db_path)) {
				std::cerr << "Error: cannot open KMC database " << uniq_kmers_kmc_db_path << "\n";
				exit(1);
			}

			if (kmer_len == 0) {
				kmer_len = kmc_file_uniq.KmerLength();
				if constexpr (std::is_same_v<BUILDER_T, KmerIndex_SBWTBuilder>) {
					if (params.precomputed_sbwt.empty())
						builder = std::make_unique<KmerIndex_SBWTBuilder>(extract_fasta_paths(input_lines), kmer_len,
							tc, params.tmp_dir);
					else
						builder = std::make_unique<KmerIndex_SBWTBuilder>(params.precomputed_sbwt, tc, params.tmp_dir);
				}
				else
				{
					//cannot static_assert(false): https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
					static_assert(!sizeof(BUILDER_T), "This constexpr if else if must be extended!");
				}
			}

			if (kmer_len != kmc_file_uniq.KmerLength()) {
				std::cerr << "k-mer len in " << uniq_kmers_kmc_db_path << " is " << kmc_file_uniq.KmerLength() << ", but the first read k-mer len was" << kmer_len << "\n";
				exit(1);
			}

			CKmerAPI kmer_api(kmer_len);
			std::vector<uint64_t> kmer_in_vec(1);

			auto query_kmer = [&kmer_api, &kmer_in_vec](CKMCFile& kmc_db, uint64_t kmer) -> uint32_t {
				kmer_in_vec[0] = kmer;
				kmer_api.from_long(kmer_in_vec);
				uint32_t res = 0;
				if (!kmc_db.CheckKmer(kmer_api, res))
					return 0;
				return res;
			};

			if (is_transcriptome) {
				//mkokot_TODO: consider
				//This part is not parallelized, may consider if querying k-mer may be parallelized like in non-transcriptome case
				//alos maybe some parts of TranscriptomeKmerFinder may be parallelized

				TranscriptomeKmerFinder transcriptome_finder(kmer_len, fasta_path);

				const size_t pack_size = 100'000;
				std::vector<task_desc> pack;
				pack.reserve(pack_size);

				auto add_kmer = [&](int cat, uint64_t kmer, uint32_t header_id, uint32_t file_id = std::numeric_limits<uint32_t>::max(), uint32_t cnt = 0) {
					pack.emplace_back(kmer, file_id, header_id, cnt, cat);
					if (pack.size() == pack.capacity())
					{
						builder->AddPack(std::move(pack));
						pack.clear();
						pack.reserve(pack_size);
					}
				};

				auto add_last_pack = [&] {
					if (!pack.empty())
					{
						builder->AddPack(std::move(pack));
						pack.clear();
					}
				};

				auto handle_unique_kmer = [&](uint64_t kmer, uint32_t id)
					{
						uint32_t counter_uniq = query_kmer(kmc_file_uniq, kmer);
						if (counter_uniq)
							add_kmer(1, kmer, id);
						else
							add_kmer(3, kmer, id);
					};

				auto handle_kmer = [&](uint64_t kmer)
					{
						uint32_t counter_uniq = query_kmer(kmc_file_uniq, kmer);
						//if is unique to this file
						if (counter_uniq)
							add_kmer(2, kmer, std::numeric_limits<uint32_t>::max(), file_id, counter_uniq);
					};

				auto register_label = [&](const std::string& gene_name) {
					return header_mapper.Add(gene_name, file_id);
				};

				transcriptome_finder.Process(handle_unique_kmer,
					handle_kmer,
					register_label);

				add_last_pack();
			}
			else {
				MultiFastaReader reader(fasta_path);

				std::cerr << "Opening KMC database " << all_kmers_kmc_db_path << "\n";
				CKMCFile kmc_file_all;
				if (!kmc_file_all.OpenForRA(all_kmers_kmc_db_path)) {
					std::cerr << "Error: cannot open KMC database " << all_kmers_kmc_db_path << "\n";
					exit(1);
				}

				if (kmer_len != kmc_file_all.KmerLength()) {
					std::cerr << "k-mer len in " << all_kmers_kmc_db_path << " is " << kmc_file_all.KmerLength() << ", but the first read k-mer len was" << kmer_len << "\n";
					exit(1);
				}

				std::string header, seq;

				parallel::ParallelBuildHelper parallel_build_helper(&kmc_file_all, &kmc_file_uniq, 
					&tc, kmer_len, params.category_3_threshold, builder.get());

				while (reader.NextSeq(header, seq)) {
					std::cerr << "*";
					WalkKmersWithInvalid kmer_walker(seq, kmer_len, true);
					uint64_t kmer;
					auto header_id = header_mapper.Add(header, file_id);

					bool is_valid;
					for (size_t i = 0; kmer_walker.Next(kmer, is_valid); ++i) {
						if (!is_valid)
							continue;

						if (poly_ACGT_filter.IsPolyACGT(kmer, kmer_len))
							continue;

						parallel_build_helper.Add(kmer, file_id, header_id);
					}
				}
				parallel_build_helper.Finish();
			}
			std::cerr << "\n";
		}
		builder->Build(file_mapper.Size(), header_mapper.Size());

		std::string out_name = params.output;
		std::ofstream out(out_name);
		if (!out) {
			std::cerr << "Error: cannot open file " << out_name << "\n";
			exit(1);
		}

		std::cerr << "Tot. k-mers in category 1: " << builder->GetInCat1() << "\n";
		std::cerr << "Tot. k-mers in category 2: " << builder->GetInCat2() << "\n";
		std::cerr << "Tot. k-mers in category 3: " << builder->GetInCat3() << "\n";

		//mkokot_TODO: In a hurry I'm forcing this, but maybe shoud be controlled via parameters
		file_mapper.TransformNames([](std::string& name) {
			auto p = name.find_last_of("/\\");
			if (p != std::string::npos)
				name = name.substr(p + 1);
		});

		header_mapper.Compact();

		Lookup::Serialize(out, params.poly_ACGT_len, file_mapper,
			header_mapper, BUILDER_T::GetIndexVariant(), 
			[&](std::ostream& out, LookupIndex& lookup_index)
			{
				builder->Serialize(out, lookup_index);
			});
	}

	void run(const Params& params) {
		switch (params.kmer_index_variant)
		{
		case kmer_index_variant_t::sbwt:
			run_impl<KmerIndex_SBWTBuilder>(params);
			break;
		default:
			std::cerr << "Error: switch needs to be extended, contact developers showing this message (" << __FILE__ << " : " << __LINE__ << ")";
			exit(1);
		}
	}
}
