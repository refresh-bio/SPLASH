#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "../common/version.h"
#include "../common/kmc_api/kmer_api.h"
#include "../common/kmc_api/kmc_file.h"
#include <refresh/hash_tables/lib/murmur_hash.h>
#include "../common/types/satc_data.h"
#include "../common/filters/poly_ACGT_filter.h"
#include "../common/filters/artifacts_filter.h"
#include "../common/hamming_filter.h"
#include "../common/filters/illumina_adapters_static.h"
#include "../common/target_count.h"

struct SampleDesc
{
	std::string input_kmc_db_path;
	uint64_t sample_id;
};

struct Params
{
	SampleDesc input_sample;
	uint64_t anchor_len{};
	uint64_t target_len{};
	uint64_t n_bins{};
	uint64_t anchor_sample_counts_threshold{}; //keep only anchors with counts > anchor_sample_counts_threshold
	uint64_t poly_ACGT_len{};
	uint64_t min_hamming_threshold{};
	std::string artifacts;
	bool dont_filter_illumina_adapters = false;
	std::string out_base;
	void Print(std::ostream& oss) const
	{
		oss << "anchor len                     : " << anchor_len << "\n";
		oss << "target len                     : " << target_len << "\n";
		oss << "n bins                         : " << n_bins << "\n";
		oss << "outbase                        : " << out_base << "\n";
		oss << "anchor_sample_counts_threshold : " << anchor_sample_counts_threshold << "\n";
		oss << "poly_ACGT_len                  : " << poly_ACGT_len << "\n";
		oss << "artifacts                      : " << artifacts << "\n";
		oss << "dont_filter_illumina_adapters  : " << std::boolalpha << dont_filter_illumina_adapters << "\n";
		oss << "min_hamming_threshold          : " << min_hamming_threshold << "\n";
		oss << "input sample                   : " << input_sample.input_kmc_db_path << " " << input_sample.sample_id << "\n";
	}
	static void Usage(char* prog_name) {
		std::cerr << "satc (sample anchor target count)\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <outbase> <input_path> <input_id>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <outbase>    - base name of output (will be extended with \".{bin_id}.bin\" \n"
			<< "    <input_path> - path to sorted (kmc1 format) kmc database\n"
			<< "    <input_id>   - id stored with each record\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_len <int> - anchor len\n"
			<< "    --target_len <int> - target len\n"
			<< "    --n_bins <int> - number of output bins\n"
			<< "    --poly_ACGT_len <int> - all anchors containing polyACGT of this length will be filtered out (0 means no filtering)\n"
			<< "    --artifacts <string> - path to artifacts, each anchor containing artifact will be filtered out\n"
			<< "    --dont_filter_illumina_adapters - if used anchors containing Illumina adapters will not be filtered out\n"
			<< "    --anchor_sample_counts_threshold <int> - keep only anchors with counts > anchor_sample_counts_threshold\n"
			<< "    --min_hamming_threshold <int> - keep only anchors with a pair of targets that differ by >= min_hamming_threshold\n";
	}
};

struct Stats {
	uint64_t tot_poly_filtered_out{};
	uint64_t tot_artifacts_filtered_out{};
	uint64_t tot_hamming_distance_filtered_out{};
	uint64_t tot_cnt_threshold_filtered_out{};
	uint64_t tot_unique_anchors{};
	uint64_t tot_out_recs{};
	uint64_t tot_in_recs{};

	void print(std::ostream& oss) const {
		oss << "# poly filtered anchors	               : " << tot_poly_filtered_out << "\n";
		oss << "# artifacts filtered anchors	       : " << tot_artifacts_filtered_out << "\n";
		oss << "# hamming distance filtered anchors    : " << tot_hamming_distance_filtered_out << "\n";
		oss << "# filtered cnt threshold anchors       : " << tot_cnt_threshold_filtered_out << "\n";
		oss << "# unique anchors                       : " << tot_unique_anchors << "\n";
		oss << "# out recs                             : " << tot_out_recs << "\n";
		oss << "# in recs                              : " << tot_in_recs << "\n";
	}
};

Params read_params(int argc, char** argv)
{
	Params res;
	if (argc == 1) {
		Params::Usage(argv[0]);
		exit(0);
	}
	int i = 1;
	for (; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];

		if (param == "--n_bins") {
			std::string tmp = argv[++i];
			res.n_bins = std::stoull(tmp);
		}
		else if (param == "--anchor_len") {
			std::string tmp = argv[++i];
			res.anchor_len = std::stoull(tmp);
		}
		else if (param == "--target_len") {
			std::string tmp = argv[++i];
			res.target_len = std::stoull(tmp);
		}
		else if (param == "--anchor_sample_counts_threshold") {
			std::string tmp = argv[++i];
			res.anchor_sample_counts_threshold = std::stoull(tmp);
		}
		else if (param == "--poly_ACGT_len") {
			std::string tmp = argv[++i];
			res.poly_ACGT_len = std::stoull(tmp);
		}
		else if (param == "--min_hamming_threshold") {
			std::string tmp = argv[++i];
			res.min_hamming_threshold = std::stoull(tmp);
		}
		else if (param == "--artifacts")
			res.artifacts = argv[++i];
		else if (param == "--dont_filter_illumina_adapters")
			res.dont_filter_illumina_adapters = true;
	}
	if (i >= argc) {
		std::cerr << "Error: outbase missing\n";
		exit(1);
	}

	res.out_base = argv[i++];
	if (i >= argc) {
		std::cerr << "Error: input sample must be specified\n";
		exit(1);
	}
	res.input_sample.input_kmc_db_path = argv[i++];
	if (i >= argc) {
		std::cerr << "Error: sample id missing\n";
		exit(1);
	}
	res.input_sample.sample_id = atoi(argv[i++]);

	if (res.n_bins == 0) {
		std::cerr << "Warning: number of bins was not specified, using default (64)\n";
		res.n_bins = 64;
	}
	if (res.anchor_len == 0) {
		std::cerr << "Error: anchor len (--anchor_len) must be specified\n";
		exit(1);
	}
	if (res.target_len == 0) {
		std::cerr << "Error: target len (--target_len) must be specified\n";
		exit(1);
	}
	return res;
}

/* Output file format
* All integers stored in LSB (little-endian)
* Header:
	 sample_id size in bytes   - 1 byte
	 barcode size in bytes    - 1 byte
	 anchor size in bytes      - 1 byte
	 target size in bytes      - 1 byte
	 counter size in bytes     - 1 byte

	 anchor len in symbols     - 1 byte
	 target len in symbols     - 1 byte
	 gap len in symbols        - 1 byte
 * List of records
	[sample_id][barcode][anchor][target][count]
*/


void split(CKmerAPI& to_split, uint64_t& anchor, uint64_t& target, uint8_t anchor_len_symbols, uint8_t target_len_symbols) {
	anchor = to_split.subkmer(0, anchor_len_symbols);
	target = to_split.subkmer(to_split.get_len() - target_len_symbols, target_len_symbols);
}

void sort_merge_and_store(uint64_t anchor,
						  std::vector<TargetCount>& targets_in_current_anchor,
						  const Header& header,
						  Record& rec,
						  std::vector<buffered_binary_writer>& bins,
						  uint64_t anchor_sample_counts_threshold,
						  const PolyACGTFilter& poly_ACGT_filter,
						  const ArtifactsFilter& artifacts_filter,
						  const HammingFilter& hamming_filter,
						  Stats& stats) {
	rec.anchor = anchor;
	if (targets_in_current_anchor.empty())
	{
		std::cerr << "Error: targets list for current anchors is empty, this should not happen\n";
		exit(1);
	}

	++stats.tot_unique_anchors;

	uint64_t tot_count{};
	for (auto& x : targets_in_current_anchor)
		tot_count += x.count;

	if (tot_count <= anchor_sample_counts_threshold) {
		++stats.tot_cnt_threshold_filtered_out;
		return;
	}

	if (poly_ACGT_filter.IsPolyACGT(anchor, header.anchor_len_symbols)) {
		++stats.tot_poly_filtered_out;
		return;
	}

	if (artifacts_filter.ContainsArtifact(anchor, header.anchor_len_symbols)) {
		++stats.tot_artifacts_filtered_out;
		return;
	}

	if (!hamming_filter.ContainsDistantPair(targets_in_current_anchor)) {
		++stats.tot_hamming_distance_filtered_out;
		return;
	}

	uint64_t n_bins = bins.size();
	uint64_t bin_id = refresh::MurMur64Hash{}(anchor) % n_bins;
	auto& bin = bins[bin_id];

	std::sort(targets_in_current_anchor.begin(), targets_in_current_anchor.end(), [](const auto& e1, const auto& e2) {return e1.target < e2.target; }); //mkokot_TODO: use parallel sort? raduls?
	rec.target = targets_in_current_anchor[0].target;
	rec.count = targets_in_current_anchor[0].count;

	for (size_t i = 1; i < targets_in_current_anchor.size(); ++i) {
		if (rec.target == targets_in_current_anchor[i].target)
			rec.count += targets_in_current_anchor[i].count;
		else {
			rec.serialize(bin, header);
			++stats.tot_out_recs;

			rec.count = targets_in_current_anchor[i].count;
			rec.target = targets_in_current_anchor[i].target;
		}
	}

	rec.serialize(bin, header);

	++stats.tot_out_recs;
}

void process_kmc_db(const std::string& path,
					uint64_t sample_id,
					const Header& header,
					std::vector<buffered_binary_writer>& bins,
					uint64_t anchor_sample_counts_threshold,
					const PolyACGTFilter& poly_ACGT_filter,
					const ArtifactsFilter& artifacts_filter,
					const HammingFilter& hamming_filter,
					Stats& stats)
{
	CKMCFile kmc_db;
	if (!kmc_db.OpenForListing(path)) {
		std::cerr << "Error: cannot open kmc db: " << path << "\n";
		exit(1);
	}

	auto kmer_len = kmc_db.KmerLength();
	if (kmc_db.IsKMC2()) {
		std::cerr << "Error: kmc_db must be sorted (use kmc_tools transform <input_db> sort <output_db>\n";
		exit(1);
	}

	CKmerAPI kmer(kmer_len);
	uint32_t count;

	uint64_t tot_kmers = kmc_db.KmerCount();

	if (!kmc_db.ReadNextKmer(kmer, count)) {
		std::cerr << "Warning: no k-mers in " << path << "\n";
		return;
	}
	++stats.tot_in_recs;

	uint64_t processed_kmers = 1;

	uint64_t prev_anchor;

	uint64_t target;
	std::vector<TargetCount> targets_in_current_anchor;

	split(kmer, prev_anchor, target, header.anchor_len_symbols, header.target_len_symbols);

	uint64_t anchor = prev_anchor;

	targets_in_current_anchor.emplace_back(target, count);
	Record rec;
	rec.barcode = 0;
	rec.sample_id = sample_id;
	while (kmc_db.ReadNextKmer(kmer, count)) {
		++stats.tot_in_recs;
		split(kmer, anchor, target, header.anchor_len_symbols, header.target_len_symbols);
		if (prev_anchor == anchor)
			targets_in_current_anchor.emplace_back(target, count);
		else {
			sort_merge_and_store(prev_anchor,
				targets_in_current_anchor,
				header,
				rec,
				bins,
				anchor_sample_counts_threshold,
				poly_ACGT_filter,
				artifacts_filter,
				hamming_filter,
				stats);

			targets_in_current_anchor.clear();
			targets_in_current_anchor.emplace_back(target, count);
			prev_anchor = anchor;
		}

		++processed_kmers;
		if (processed_kmers % 100'000'000 == 0)
			std::cerr << "\r" << processed_kmers << "/" << tot_kmers << "\n";
	}
	std::cerr << "\n";

	sort_merge_and_store(prev_anchor,
		targets_in_current_anchor,
		header,
		rec,
		bins,
		anchor_sample_counts_threshold,
		poly_ACGT_filter,
		artifacts_filter,
		hamming_filter,
		stats);
	targets_in_current_anchor.clear();
}

uint8_t verify_kmc_dbs_and_get_gap_len(SampleDesc& input_sample, uint64_t anchor_len, uint64_t target_len) {

	CKMCFile kmc_file;
	if (!kmc_file.OpenForListing(input_sample.input_kmc_db_path)) {
		std::cerr << "Erorr: cannot open kmc db: " << input_sample.input_kmc_db_path << "\n";
		exit(1);
	}
	uint32_t kmer_len = kmc_file.KmerLength();
	kmc_file.Close();

	if (anchor_len + target_len > kmer_len) {
		std::cerr << "Error: anchor_len + target_len > kmer_len\n";
		exit(1);
	}
	return kmer_len - anchor_len - target_len;
}

//mkokot_TODO: copied from txc/utils.h -> refactor this!
constexpr uint32_t no_bytes(uint32_t x)
{
	if (x < 256)
		return 1;
	if (x < 256 * 256)
		return 2;
	if (x < 256 * 256 * 256)
		return 3;
	return 4;
}

int main(int argc, char** argv)
{
	std::cerr << "Welcome to satc (sample anchor target count)\n";

	auto params = read_params(argc, argv);

	params.Print(std::cerr);

	Header header;

	header.sample_id_size_bytes = no_bytes(params.input_sample.sample_id);
	header.barcode_size_bytes = 0;
	header.counter_size_bytes = 2;
	header.barcode_len_symbols = 0;
	header.anchor_len_symbols = params.anchor_len;
	header.target_len_symbols = params.target_len;
	header.anchor_size_bytes = (header.anchor_len_symbols + 3) / 4;
	header.target_size_bytes = (header.target_len_symbols + 3) / 4;

	header.ordering = Header::ordering_t::SBATC;

	header.gap_len_symbols = verify_kmc_dbs_and_get_gap_len(params.input_sample, params.anchor_len, params.target_len);

	header.print(std::cerr);

	std::vector<buffered_binary_writer> out_files;
	for (size_t i = 0; i < params.n_bins; ++i) {
		auto fname = params.out_base + "." + std::to_string(i) + ".bin";
		out_files.emplace_back(fname);
		auto& outfile = out_files.back();
		if (!outfile) {
			std::cerr << "Error: cannot open file " << fname << "\n";
		}
		header.serialize(outfile);
	}

	Stats stats;

	PolyACGTFilter poly_ACGT_filter(params.poly_ACGT_len);
	ArtifactsFilter artifacts_filter(params.artifacts);
	HammingFilter hamming_filter(params.min_hamming_threshold);

	if (!params.dont_filter_illumina_adapters)
		artifacts_filter.Add(12, IlluminaAdaptersStatic::Get12Mers());

	process_kmc_db(params.input_sample.input_kmc_db_path,
		params.input_sample.sample_id,
		header,
		out_files,
		params.anchor_sample_counts_threshold,
		poly_ACGT_filter,
		artifacts_filter,
		hamming_filter,
		stats);

	stats.print(std::cerr);
}
