#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <memory>
#include <chrono>
#include "../common/version.h"
#include "pvals.h"
#include "extra_stats.h"
#include "../common/accepted_anchors.h"
#include "anchor.h"
#include <set>

class Timer {
	using time_type = decltype(std::chrono::high_resolution_clock::now());
	time_type start;
	time_type get_time() const {
		return std::chrono::high_resolution_clock::now();
	}
public:
	Timer() {
		Start();
	}

	void Start() {
		start = get_time();
	}
	double GetElapsed() const {
		return std::chrono::duration<double>(get_time() - start).count();
	}
};

struct Params
{
	uint64_t anchor_count_threshold{}; //remove all anchors for which the sum of counters across all samples and targets is <= anchor_count_threshold
	uint64_t anchor_unique_targets_threshold{}; // remove all anchors having <= anchor_unique_targets_threshold unique targets
	uint64_t anchor_samples_threshold{}; //remove all anchors having <= anchor_samples_threshold unique samples

	uint32_t n_most_freq_targets{}; //in stats mode print also n_most_freq_targets and their counts in output file

	uint32_t n_most_freq_targets_for_stats{}; // if 0 - keep all, else keep only <n_most_freq_targets_for_stats>

	double opt_train_fraction = 0.25;

	int opt_num_inits = 10;

	int opt_num_iters = 50;

	size_t num_rand_cf = 50;

	size_t num_splits = 1;

	bool without_alt_max = false;

	bool with_effect_size_cts = false;

	bool with_pval_asymp_opt = false;

	bool compute_also_old_base_pvals = false;

	std::string outpath;

	std::string anchor_list;

	std::string cjs_out;

	double max_pval_opt_for_Cjs = 0.1;

	std::string sample_names;

	std::vector<std::string> bins;

	std::string dump_sample_anchor_target_count_txt;
	std::string dump_sample_anchor_target_count_binary;

	RecFmt format = RecFmt::SATC; //only for JustMergeAndDump

	std::string cell_type_samplesheet;

	void print(std::ostream& oss) {
		oss << "Parameters:\n";
		oss << "\tanchor_count_threshold                  : " << anchor_count_threshold << "\n";
		oss << "\tanchor_unique_targets_threshold         : " << anchor_unique_targets_threshold << "\n";
		oss << "\tanchor_samples_threshold                : " << anchor_samples_threshold << "\n";
		oss << "\twithout_alt_max                         : " << std::boolalpha << without_alt_max << "\n";
		oss << "\twith_effect_size_cts                    : " << std::boolalpha << with_effect_size_cts << "\n";
		oss << "\twith_pval_asymp_opt                     : " << std::boolalpha << with_pval_asymp_opt << "\n";
		oss << "\tcompute_also_old_base_pvals             : " << std::boolalpha << compute_also_old_base_pvals << "\n";
		oss << "\toutpath                                 : " << outpath << "\n";
		oss << "\tanchor_list                             : " << anchor_list << "\n";
		oss << "\tcjs_out                                 : " << cjs_out << "\n";
		oss << "\tmax_pval_opt_for_Cjs                    : " << max_pval_opt_for_Cjs << "\n";
		oss << "\tsample_names	                          : " << sample_names << "\n";
		oss << "\tn_most_freq_targets                     : " << n_most_freq_targets << "\n";
		oss << "\tn_most_freq_targets_for_stats           : " << n_most_freq_targets_for_stats << "\n";
		oss << "\topt_train_fraction                      : " << opt_train_fraction << "\n";
		oss << "\topt_num_inits                           : " << opt_num_inits << "\n";
		oss << "\topt_num_iters                           : " << opt_num_iters << "\n";
		oss << "\tnum_rand_cf                             : " << num_rand_cf << "\n";
		oss << "\tnum_splits                              : " << num_splits << "\n";
		oss << "\tcell_type_samplesheet                   : " << cell_type_samplesheet << "\n";
		oss << "\t dump_sample_anchor_target_count_txt    : " << dump_sample_anchor_target_count_txt << "\n";
		oss << "\t dump_sample_anchor_target_count_binary : " << dump_sample_anchor_target_count_binary << "\n";
		oss << "\t format                                 : " << RecFmtConv::to_string(format) << "\n";
		oss << "\tinput bins:\n";
		for (const auto& bin : bins)
			oss << "\t\t" << bin << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_merge\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <path> <outpath> <list_of_bins_to_merge>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <outpath>                           - output path\n"
			<< "    <file_with_list_of_bins_to_merge>   - file with list of paths of bins to be processed\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_count_threshold <int>                    - filter out all anchors for which the total count <= anchor_count_threshold\n"
			<< "    --anchor_unique_targets_threshold <int>           - filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold\n"
			<< "    --n_most_freq_targets <int>                       - output also n_most_freq_targets most frequent targets and their counts (default: 0)\n"
			<< "    --n_most_freq_targets_for_stats <int>             - use at most n_most_freq_targets_for_stats for each contingency table, 0 means use all (default: 0)\n"
			<< "    --opt_train_fraction <double> (0.0;1.0)           - use this fraction to create train X from contingency table\n"
			<< "    --opt_num_inits <int>                             - the number of altMaximize runs (default: 10)\n"
			<< "    --opt_num_iters <int>                             - the number of iteration in altMaximize (default: 50)\n"
			<< "    --num_rand_cf <int>                               - the number of rand cf (default: 50)\n"
			<< "    --num_splits <int>                                - the number of contingency table splits (default: 1)\n"
			<< "    --anchor_samples_threshold <int>                  - filter out all anchors for which the number of unique samples is <= anchor_samples_threshold\n"
			<< "    --anchor_list <string>                            - path to text file containing anchors separated by whitespaces, only anchors from this file will be processed\n"
			<< "    --cjs_out <string>                                - path to output text file where Cjs will be stored\n"
			<< "    --max_pval_opt_for_Cjs <double>                   - dump only Cjs for anchors that have pval_opt <= max_pval_opt_for_Cjs\n"
			<< "    --dump_sample_anchor_target_count_txt <string>    - dump merged anchors in textual representation\n"
			<< "    --dump_sample_anchor_target_count_binary <string> - dump merged anchors in textual satc format\n"
			<< "    --without_alt_max                                 - disable alt max computation\n"
			<< "    --with_effect_size_cts                            - compute effect_size_cts\n"
			<< "    --with_pval_asymp_opt                             - compute pval_asymp_opt\n"
			<< "    --compute_also_old_base_pvals                     - compute old base pvals\n"
			<< "    --sample_names <path>                             - path for decode sample id, each line should contain <sample_name> <sample_id>\n"
			<< "    --format <string>                                 - output format when txt dump, available options: satc, splash (default: satc)\n";
	}
};


struct Stats
{
	uint64_t tot_filtered_out_anchors{};

	uint64_t tot_writen_anchors{};
	uint64_t tot_writen_records{};
	std::pair<uint64_t, uint64_t> max_contignency_matrix_size{};

	void print(std::ostream& oss) {
		oss << "tot writen anchors                                : " << tot_writen_anchors << "\n";
		oss << "tot writen records                                : " << tot_writen_records << "\n";
		oss << "n samples in max contignency matrix               : " << max_contignency_matrix_size.first << "\n";
		oss << "n targets in max contignency matrix               : " << max_contignency_matrix_size.second << "\n";
		oss << "tot_filtered_out_anchors				          : " << tot_filtered_out_anchors << "\n";
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
	for (; i < argc; ++i) {
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];

		if (param == "--anchor_count_threshold") {
			std::string tmp = argv[++i];
			res.anchor_count_threshold = std::stoull(tmp);
		}
		if (param == "--anchor_unique_targets_threshold") {
			std::string tmp = argv[++i];
			res.anchor_unique_targets_threshold = std::stoull(tmp);
		}
		if (param == "--anchor_samples_threshold") {
			std::string tmp = argv[++i];
			res.anchor_samples_threshold = std::stoull(tmp);
		}
		if (param == "--n_most_freq_targets") {
			std::string tmp = argv[++i];
			res.n_most_freq_targets = std::stoull(tmp);
		}
		if (param == "--n_most_freq_targets_for_stats") {
			std::string tmp = argv[++i];
			res.n_most_freq_targets_for_stats = std::stoull(tmp);
		}
		if (param == "--opt_train_fraction") {
			std::string tmp = argv[++i];
			double val = std::stod(tmp);
			if (val <= 0.0 || val >= 1.0) {
				std::cerr << "Error: opt_train_fraction must be in range (0.0;1.0)\n";
				exit(1);
			}
			res.opt_train_fraction = val;
		}
		if (param == "--opt_num_inits") {
			std::string tmp = argv[++i];
			res.opt_num_inits = std::stoull(tmp);
		}
		if (param == "--opt_num_iters") {
			std::string tmp = argv[++i];
			res.opt_num_iters = std::stoull(tmp);
		}
		if (param == "--num_rand_cf") {
			std::string tmp = argv[++i];
			res.num_rand_cf = std::stoull(tmp);
		}
		if (param == "--num_splits") {
			std::string tmp = argv[++i];
			res.num_splits = std::stoull(tmp);
		}
		if (param == "--without_alt_max") {
			res.without_alt_max = true;
		}
		if (param == "--with_effect_size_cts") {
			res.with_effect_size_cts = true;
		}
		if (param == "--with_pval_asymp_opt") {
			res.with_pval_asymp_opt = true;
		}
		if (param == "--compute_also_old_base_pvals") {
			res.compute_also_old_base_pvals = true;
		}
		if (param == "--anchor_list") {
			res.anchor_list = argv[++i];
		}
		if (param == "--cjs_out") {
			res.cjs_out = argv[++i];
		}
		if (param == "--max_pval_opt_for_Cjs") {
			std::string tmp = argv[++i];
			res.max_pval_opt_for_Cjs = std::stod(tmp);
		}
		if (param == "--sample_names") {
			res.sample_names = argv[++i];
		}
		if (param == "--cell_type_samplesheet") {
			res.cell_type_samplesheet = argv[++i];
		}
		if (param == "--dump_sample_anchor_target_count_txt") {
			res.dump_sample_anchor_target_count_txt = argv[++i];
		}
		if (param == "--dump_sample_anchor_target_count_binary") {
			res.dump_sample_anchor_target_count_binary = argv[++i];
		}

		if (param == "--format")
			res.format = RecFmtConv::from_string(argv[++i]);
	}
	if (i >= argc) {
		std::cerr << "Error: outpath missing\n";
		exit(1);
	}

	res.outpath = argv[i++];

	if (i >= argc) {
		std::cerr << "Error: path to file with input list missing\n";
		exit(1);
	}
	std::string input_list_path = argv[i++];
	if (i < argc) {
		std::cerr << "Error: following arguments are unexpected:\n";
		for (int j = i; j < argc; ++j)
			std::cerr << argv[j] << "\n";
		exit(1);
	}
	std::ifstream in(input_list_path);
	if (!in) {
		std::cerr << "Error: cannot open file " << input_list_path << "\n";
		exit(1);
	}
	std::string line;
	while (std::getline(in, line)) {
		res.bins.emplace_back(line);
	}

	if (!res.bins.size()) {
		std::cerr << "Error: at leas one input bin must be specified\n";
		exit(1);
	}

	return res;
}




class CachedRecord {
	bool is_cached = false;
	Record _rec{};
	buffered_binary_reader& in;
public:
	CachedRecord(buffered_binary_reader& in) :
		in(in) {

	}
	bool Peek(const Header& header, Record& rec) {
		if (is_cached) {
			rec = _rec;
			return true;
		}

		if (!_rec.load(in, header))
			return false;

		is_cached = true;
		rec = _rec;
		return true;
	}
	void Skip() {
		is_cached = false;
	}
};

class Bin {
	buffered_binary_reader in;
	Header header;
	Non10SingleSampleAnchor current_anchor;
	bool is_loaded = false;
	CachedRecord cached_rec;
	bool load_anchor() {
		current_anchor.data.clear();
		Record rec;

		if (!cached_rec.Peek(header, rec))
			return false;
		cached_rec.Skip();

		current_anchor.anchor = rec.anchor;
		current_anchor.sample_id = rec.sample_id;
		current_anchor.data.emplace_back(rec.target, rec.count);
		while (cached_rec.Peek(header, rec)) {
			assert(rec.sample_id == current_anchor.sample_id);
			if (rec.anchor == current_anchor.anchor) {
				current_anchor.data.emplace_back(rec.target, rec.count);
				cached_rec.Skip();
			}
			else
				break;
		}
		is_loaded = true;
		return true;
	}
public:
	explicit Bin(const std::string& path) :
		in(path),
		cached_rec(in)
	{
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		header.load(in);
		if (!load_anchor()) {
			std::cerr << "Warning: no anchors in " << path << "\n";
		}
	}

	const Header& get_header() const {
		return header;
	}

	bool PeekAnchor(uint64_t& anchor)
	{
		if (!is_loaded) {
			if (!load_anchor())
				return false;
		}
		anchor = current_anchor.anchor;
		return true;
	}

	void GetAnchor(Non10SingleSampleAnchor& anchor) {
		if (!is_loaded) {
			std::cerr << "Error: wrnog GetAnchor function call\n";
			exit(1);
		}
		anchor = std::move(current_anchor);
		is_loaded = false;
	}

	void Skip() {
		if (!is_loaded) {
			std::cerr << "Error: wrnog Skip function call\n";
			exit(1);
		}
		is_loaded = false;
	}
};


void write_out_header(std::ofstream& out, bool without_alt_max, bool with_effect_size_cts, bool with_pval_asymp_opt, bool compute_also_old_base_pvals, bool is_removing_least_freq_targets_enabled, uint32_t n_most_freq_targets) {
	//mkokot_TODO: czy doimplementowac zakomentowane?
	out
		<< "anchor" << "\t";
	if (compute_also_old_base_pvals) {
		out << "pval_base_old" << "\t";
		out << "effect_size_bin_old" << "\t";
	}
	if (!without_alt_max) {
		out << "pval_opt" << "\t";

		if (with_effect_size_cts)
			out << "effect_size_cts" << "\t";
	}

	if (!without_alt_max) {
		out
			<< "effect_size_bin" << "\t";

		if (with_pval_asymp_opt)
			out << "pval_asymp_opt" << "\t";
	}
	out
		//<< "effect_size_cts_SVD" << "\t"					// TODO: ?
		//<< "pval_rand_init_EM" << "\t"					// TODO: ?
		<< "pval_base" << "\t"
		//<< "effect_size_base" << "\t"						// TODO: ?
		<< "M" << "\t"									// X.sum()
		<< "anch_uniqTargs" << "\t";							// (X.sum(axis=1)>0).sum()
	if (is_removing_least_freq_targets_enabled)
		out
			<< "M_before_filter" << "\t"
			<< "anch_uniqTargs_before_filter" << "\t";
	out
		<< "number_nonzero_samples" << "\t"					// (X.sum(axis=0)>0).sum()
		<< "target_entropy" << "\t"						// scipy.stats.entropy(X.sum(axis=1),base=2)
		//<< "entropy_difference" << "\t"					// TODO: ?
		//<< "mean_target_levenshtein_distance" << "\t"
/* << "mean_target_hamming_distance" << "\t"
		<< "avg_edit_distance" << "\t"
		<< "avg_no_homopolymer_targets" << "\t"*/
		<< "avg_no_homopolymer_targets\t"
		<< "avg_hamming_distance_max_target\t"
		<< "avg_hamming_distance_all_pairs\t"
		<< "avg_edit_distance_max_target\t"
		<< "avg_edit_distance_all_pairs\t";

	for (size_t i = 1; i <= n_most_freq_targets; ++i) {
		out << "most_freq_target_" << i << "\t";
		out << "cnt_most_freq_target_" << i << "\t";
	}
	out
		<< "\n";
}

void write_out_rec(
	std::ofstream& out,
	const AnchorStats& anchor_stats,
	uint64_t anchor,
	size_t anchor_len_symbols,
	size_t target_len_symbols,
	uint64_t tot_cnt,
	uint64_t n_unique_targets,
	uint64_t tot_cnt_before_filter,
	uint64_t n_uniqe_targets_before_filter,
	uint64_t n_unique_samples,
	bool without_alt_max,
	bool with_effect_size_cts,
	bool with_pval_asymp_opt,
	bool compute_also_old_base_pvals,
	bool is_removing_least_freq_targets_enabled,
	uint32_t n_most_freq_targets) {
	out
		<< kmer_to_string(anchor, anchor_len_symbols) << "\t";

	if (compute_also_old_base_pvals) {
		out
			<< anchor_stats.pval_base_old << "\t";
		out
			<< anchor_stats.effect_size_bin_old << "\t";
	}
	if (!without_alt_max) {
		out
			<< anchor_stats.pval_opt << "\t";

		if (with_effect_size_cts)
			out << anchor_stats.effect_size_cts << "\t";
	}

	if (!without_alt_max) {
		out << anchor_stats.effect_size_bin << "\t";
		if (with_pval_asymp_opt)
			out << anchor_stats.pval_asymp_opt << "\t";
	}
	out
		//<< "effect_size_cts_SVD" << "\t"					// TODO: ?
		//<< "pval_rand_init_EM" << "\t"					// TODO: ?
		<< anchor_stats.pval_base << "\t"
		//<< "effect_size_base" << "\t"						// TODO: ?
		<< tot_cnt << "\t"									// X.sum()
		<< n_unique_targets << "\t";							// (X.sum(axis=1)>0).sum()
	if (is_removing_least_freq_targets_enabled) {
		out
			<< tot_cnt_before_filter << "\t"
			<< n_uniqe_targets_before_filter << "\t";
	}
	out
		<< n_unique_samples << "\t"							// (X.sum(axis=0)>0).sum()
		//		<< anchor_stats.entropy << "\t"						// scipy.stats.entropy(X.sum(axis=1),base=2)
				//<< "entropy_difference" << "\t"					// TODO: ?
				//<< "mean_target_levenshtein_distance" << "\t"
		<< anchor_stats.entropy << "\t"
		<< anchor_stats.avg_no_homopolymer_targets << "\t"
		<< anchor_stats.avg_hamming_distance_max_target << "\t"
		<< anchor_stats.avg_hamming_distance_all_pairs << "\t"
		<< anchor_stats.avg_edit_distance_max_target << "\t"
		<< anchor_stats.avg_edit_distance_all_pairs << "\t";
	for (size_t i = 0; i < anchor_stats.most_freq_targets.size(); ++i) {
		out << kmer_to_string(anchor_stats.most_freq_targets[i].kmer, target_len_symbols) << "\t";
		out << anchor_stats.most_freq_targets[i].counter << "\t";
	}
	//if there were less targets than n_most_freq_targets for current anchor
	for (size_t i = anchor_stats.most_freq_targets.size(); i < n_most_freq_targets; ++i) {
		out << "-\t";
		out << "0\t";
	}
	out
		<< "\n";
}

bool anchor_filtered_out(size_t tot_cnt, size_t n_unique_targets, size_t n_unique_samples, const Params& params) {
	return
		tot_cnt <= params.anchor_count_threshold ||
		n_unique_targets <= params.anchor_unique_targets_threshold ||
		n_unique_samples <= params.anchor_samples_threshold;
}

class IAnchorProcessor {
public:
	virtual void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniqe_targets,
		size_t tot_cnt,
		size_t n_uniqe_targets_before_filter,
		size_t tot_cnt_before_filter,
		const std::unordered_set<uint64_t>& unique_samples) = 0;

	virtual ~IAnchorProcessor() = default;
};

class StatsWriter : public IAnchorProcessor {
	std::ofstream out;
	bool without_alt_max;
	bool with_effect_size_cts;
	bool with_pval_asymp_opt;
	bool compute_also_old_base_pvals;
	bool is_removing_least_freq_targets_enabled;
	uint32_t n_most_freq_targets;
	double opt_train_fraction;
	int opt_num_inits;
	int opt_num_iters;
	size_t num_rand_cf;
	size_t num_splits;
	CExtraStats extra_stats;
	AnchorStats anchor_stats;
	CjWriter cj_writer;
	double max_pval_opt_for_Cjs;
	
	const size_t io_buffer_size = 1 << 20;
	char* io_buffer;

public:
	StatsWriter(
		const std::string& outpath,
		bool without_alt_max,
		bool with_effect_size_cts,
		bool with_pval_asymp_opt,
		bool compute_also_old_base_pvals,
		bool is_removing_least_freq_targets_enabled,
		uint32_t n_most_freq_targets,
		double opt_train_fraction,
		int opt_num_inits,
		int opt_num_iters,
		size_t num_rand_cf,
		size_t num_splits,
		const std::string& cjs_out,
		size_t anchor_len_symbols,
		size_t barcode_len_symbols,
		const std::string& sample_names,
		double max_pval_opt_for_Cjs) :
		out(outpath, std::ios_base::binary),
//		out(outpath),
		without_alt_max(without_alt_max),
		with_effect_size_cts(with_effect_size_cts),
		with_pval_asymp_opt(with_pval_asymp_opt),
		compute_also_old_base_pvals(compute_also_old_base_pvals),
		is_removing_least_freq_targets_enabled(is_removing_least_freq_targets_enabled),
		n_most_freq_targets(n_most_freq_targets),
		opt_train_fraction(opt_train_fraction),
		opt_num_inits(opt_num_inits),
		opt_num_iters(opt_num_iters),
		num_rand_cf(num_rand_cf),
		num_splits(num_splits),
		cj_writer(cjs_out, anchor_len_symbols, barcode_len_symbols, sample_names),
		max_pval_opt_for_Cjs(max_pval_opt_for_Cjs) {
		if (!out) {
			std::cerr << "Error: cannot open file " << outpath << "\n";
			exit(1);
		}

		io_buffer = new char[io_buffer_size];
		out.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

		write_out_header(
			out,
			without_alt_max,
			with_effect_size_cts,
			with_pval_asymp_opt,
			compute_also_old_base_pvals,
			is_removing_least_freq_targets_enabled,
			n_most_freq_targets);
	}

	~StatsWriter()
	{
		out.close();
		delete[] io_buffer;
	}

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniqe_targets,
		size_t tot_cnt,
		size_t n_uniqe_targets_before_filter,
		size_t tot_cnt_before_filter,
		const std::unordered_set<uint64_t>& unique_samples
	)  override {

		auto _anchor = anchor.anchor;

		extra_stats.Compute(anchor, target_len_symbols, 5, 200, n_uniqe_targets, unique_samples, anchor_stats);

		compute_stats(
			std::move(anchor),
			anchor_len_symbols,
			n_uniqe_targets,
			unique_samples,
			anchor_stats,
			without_alt_max,
			with_effect_size_cts,
			with_pval_asymp_opt,
			compute_also_old_base_pvals,
			n_most_freq_targets,
			opt_train_fraction,
			opt_num_inits,
			opt_num_iters,
			num_rand_cf,
			num_splits,
			cj_writer,
			max_pval_opt_for_Cjs);

		write_out_rec(
			out,
			anchor_stats,
			_anchor,
			anchor_len_symbols,
			target_len_symbols,
			tot_cnt,
			n_uniqe_targets,
			tot_cnt_before_filter,
			n_uniqe_targets_before_filter,
			unique_samples.size(),
			without_alt_max,
			with_effect_size_cts,
			with_pval_asymp_opt,
			compute_also_old_base_pvals,
			is_removing_least_freq_targets_enabled,
			n_most_freq_targets);
	}
};

class SatcWriter : public IAnchorProcessor {
	buffered_binary_writer out;
	Header out_header;
public:
	SatcWriter(const std::string& outpath,
		uint8_t sample_id_size_bytes,
		uint8_t barcode_size_bytes,
		uint8_t anchor_size_bytes,
		uint8_t target_size_bytes,
		uint8_t counter_size_bytes,
		uint8_t barcode_len_symbols,
		uint8_t anchor_len_symbols,
		uint8_t target_len_symbols,
		uint8_t gap_len_symbols
	) :out(outpath)
	{

		if (!out) {
			std::cerr << "Error: cannot open file " << outpath << "\n";
			exit(1);
		}
		out_header.sample_id_size_bytes = sample_id_size_bytes;
		out_header.barcode_size_bytes = barcode_size_bytes;
		out_header.anchor_size_bytes = anchor_size_bytes;
		out_header.target_size_bytes = target_size_bytes;
		out_header.counter_size_bytes = counter_size_bytes;
		out_header.barcode_len_symbols = barcode_len_symbols;
		out_header.anchor_len_symbols = anchor_len_symbols;
		out_header.target_len_symbols = target_len_symbols;
		out_header.gap_len_symbols = gap_len_symbols;

		out_header.serialize(out);
	}

	void ProcessAnchor(const Anchor& anchor) {
		Record rec;
		rec.anchor = anchor.anchor;

		for (const auto& x : anchor.data) {
			rec.barcode = x.barcode;
			rec.count = x.count;
			rec.sample_id = x.sample_id;
			rec.target = x.target;

			rec.serialize(out, out_header);
		}
	}

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniqe_targets,
		size_t tot_cnt,
		size_t n_uniqe_targets_before_filter,
		size_t tot_cnt_before_filter,
		const std::unordered_set<uint64_t>& unique_samples
	)  override {
		ProcessAnchor(anchor);

		anchor.data.clear();
		anchor.data.shrink_to_fit();
	}
};

class SatcDumpWriter : public IAnchorProcessor {
	std::ofstream out;
	Header header;
	RecFmt format;
	SampleNameDecoder sample_name_decoder;
public:
	SatcDumpWriter(const std::string& outpath,
		RecFmt format,
		uint8_t sample_id_size_bytes,
		uint8_t barcode_size_bytes,
		uint8_t anchor_size_bytes,
		uint8_t target_size_bytes,
		uint8_t counter_size_bytes,
		uint8_t barcode_len_symbols,
		uint8_t anchor_len_symbols,
		uint8_t target_len_symbols,
		uint8_t gap_len_symbols,
		const std::string& sample_names
	) :
	out(outpath),
	format(format),
	sample_name_decoder(sample_names) {
		if (!out) {
			std::cerr << "Error: cannot open file " << outpath << "\n";
			exit(1);
		}

		header.sample_id_size_bytes = sample_id_size_bytes;
		header.barcode_size_bytes = barcode_size_bytes;
		header.anchor_size_bytes = anchor_size_bytes;
		header.target_size_bytes = target_size_bytes;
		header.counter_size_bytes = counter_size_bytes;
		header.barcode_len_symbols = barcode_len_symbols;
		header.anchor_len_symbols = anchor_len_symbols;
		header.target_len_symbols = target_len_symbols;
		header.gap_len_symbols = gap_len_symbols;
	}

	void ProcessAnchor(const Anchor& anchor)
	{
		Record rec;
		rec.anchor = anchor.anchor;

		for (const auto& x : anchor.data) {
			rec.barcode = x.barcode;
			rec.count = x.count;
			rec.sample_id = x.sample_id;
			rec.target = x.target;

			rec.print(out, header, format, sample_name_decoder);
		}
	}

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniqe_targets,
		size_t tot_cnt,
		size_t n_uniqe_targets_before_filter,
		size_t tot_cnt_before_filter,
		const std::unordered_set<uint64_t>& unique_samples
	) override {
		ProcessAnchor(anchor);

		anchor.data.clear();
		anchor.data.shrink_to_fit();
	}
};

class CombinedAnchorWriter : public IAnchorProcessor {
	std::unique_ptr<SatcDumpWriter> satc_dump_writer;
	std::unique_ptr<SatcWriter> satc_writer;
	std::unique_ptr<StatsWriter> stats_writer;

public:
	CombinedAnchorWriter(
		std::unique_ptr<SatcDumpWriter>&& satc_dump_writer,
		std::unique_ptr<SatcWriter>&& satc_writer,
		std::unique_ptr<StatsWriter>&& stats_writer)
		:
		satc_dump_writer(std::move(satc_dump_writer)),
		satc_writer(std::move(satc_writer)),
		stats_writer(std::move(stats_writer))
	{

	}

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniqe_targets,
		size_t tot_cnt,
		size_t n_uniqe_targets_before_filter,
		size_t tot_cnt_before_filter,
		const std::unordered_set<uint64_t>& unique_samples
	)  override {
		if (satc_dump_writer)
			satc_dump_writer->ProcessAnchor(anchor);
		if (satc_writer)
			satc_writer->ProcessAnchor(anchor);

		//must be last because it destroys anchor content
		if (stats_writer)
			stats_writer->ProcessAnchor(std::move(anchor),
				anchor_len_symbols,
				target_len_symbols,
				n_uniqe_targets,
				tot_cnt,
				n_uniqe_targets_before_filter,
				tot_cnt_before_filter,
				unique_samples
			);
	}
};

bool get_top_anchors(std::vector<std::unique_ptr<Bin>>& bins, AcceptedAnchors& anchor_filter, std::vector<uint64_t>& top_anchors) {
	for (size_t bin_id = 0; bin_id < bins.size(); ++bin_id) {
		uint64_t anchor;
		if (!bins[bin_id]->PeekAnchor(anchor)) {
			bins[bin_id--] = std::move(bins.back());
			bins.pop_back();
		}
		else {
			if (anchor_filter.IsAccepted(anchor))
				top_anchors.push_back(anchor);
			else { //skip this anchor
				while (true) {
					bins[bin_id]->Skip();
					uint64_t new_anchor;
					if (!bins[bin_id]->PeekAnchor(new_anchor)) {
						bins[bin_id--] = std::move(bins.back());
						bins.pop_back();
						break;
					}
					else if (new_anchor != anchor) {
						--bin_id;
						break;
					}
				}
			}
		}
	}
	return !bins.empty();
}

std::unique_ptr<IAnchorProcessor> get_anchor_processor(
	const std::string& outpath,
	const std::string& dump_sample_anchor_target_count_txt,
	const std::string dump_sample_anchor_target_count_binary,
	RecFmt format,
	uint8_t sample_id_size_bytes,
	uint8_t barcode_size_bytes,
	uint8_t anchor_size_bytes,
	uint8_t target_size_bytes,
	uint8_t counter_size_bytes,
	uint8_t barcode_len_symbols,
	uint8_t anchor_len_symbols,
	uint8_t target_len_symbols,
	uint8_t gap_len_symbols,
	const std::string& sample_names,
	bool without_alt_max,
	bool with_effect_size_cts,
	bool with_pval_asymp_opt,
	bool compute_also_old_base_pvals,
	bool is_removing_least_freq_targets_enabled,
	uint32_t n_most_freq_targets,
	double opt_train_fraction,
	int opt_num_inits,
	int opt_num_iters,
	size_t num_rand_cf,
	size_t num_splits,
	const std::string& cjs_out,
	double max_pval_opt_for_Cjs) {

	std::unique_ptr<SatcDumpWriter> satc_dump_writer;

	if (dump_sample_anchor_target_count_txt != "") {
		satc_dump_writer = std::make_unique<SatcDumpWriter>(
			dump_sample_anchor_target_count_txt,
			format,
			sample_id_size_bytes,
			barcode_size_bytes,
			anchor_size_bytes,
			target_size_bytes,
			counter_size_bytes,
			barcode_len_symbols,
			anchor_len_symbols,
			target_len_symbols,
			gap_len_symbols,
			sample_names
			);
	}

	std::unique_ptr<SatcWriter> satc_writer;

	if (dump_sample_anchor_target_count_binary != "") {
		satc_writer = std::make_unique<SatcWriter>(
			dump_sample_anchor_target_count_binary,
			sample_id_size_bytes,
			barcode_size_bytes,
			anchor_size_bytes,
			target_size_bytes,
			counter_size_bytes,
			barcode_len_symbols,
			anchor_len_symbols,
			target_len_symbols,
			gap_len_symbols
			);
	}

	std::unique_ptr<StatsWriter> stats_writer = std::make_unique<StatsWriter>(
		outpath,
		without_alt_max,
		with_effect_size_cts,
		with_pval_asymp_opt,
		compute_also_old_base_pvals,
		is_removing_least_freq_targets_enabled,
		n_most_freq_targets,
		opt_train_fraction,
		opt_num_inits,
		opt_num_iters,
		num_rand_cf,
		num_splits,
		cjs_out,
		anchor_len_symbols,
		barcode_len_symbols,
		sample_names,
		max_pval_opt_for_Cjs);

	return std::make_unique<CombinedAnchorWriter>(
		std::move(satc_dump_writer),
		std::move(satc_writer),
		std::move(stats_writer));
}

void run(const Params& params) {
	AcceptedAnchors anchor_filter(params.anchor_list);

	std::vector<std::unique_ptr<Bin>> bins;
	for (auto& path : params.bins)
		bins.emplace_back(std::make_unique<Bin>(path));

	Stats stats;

	const auto bin0_header = bins[0]->get_header();
	uint8_t max_sample_id_size_bytes = bin0_header.sample_id_size_bytes;

	for (size_t i = 1; i < bins.size(); ++i) {
		auto& bin_i_header = bins[i]->get_header();
		if (bin0_header.anchor_len_symbols != bin_i_header.anchor_len_symbols) {
			std::cerr << "Error: bins have different anchor lengths\n";
			exit(1);
		}

		if (bin0_header.target_len_symbols != bin_i_header.target_len_symbols) {
			std::cerr << "Error: bins have different target lengths\n";
			exit(1);
		}

		if (bin0_header.gap_len_symbols != bin_i_header.gap_len_symbols) {
			std::cerr << "Error: bins have different gap (lookahead) lengths\n";
			exit(1);
		}

		if (bin_i_header.sample_id_size_bytes > max_sample_id_size_bytes)
			max_sample_id_size_bytes = bin_i_header.sample_id_size_bytes;
	}

	std::unique_ptr<IAnchorProcessor> anchor_processor = get_anchor_processor(
		params.outpath,
		params.dump_sample_anchor_target_count_txt,
		params.dump_sample_anchor_target_count_binary,
		params.format,
		max_sample_id_size_bytes,
		bin0_header.barcode_size_bytes,
		bin0_header.anchor_size_bytes,
		bin0_header.target_size_bytes,
		bin0_header.counter_size_bytes,
		bin0_header.barcode_len_symbols,
		bin0_header.anchor_len_symbols,
		bin0_header.target_len_symbols,
		bin0_header.gap_len_symbols,
		params.sample_names,
		params.without_alt_max,
		params.with_effect_size_cts,
		params.with_pval_asymp_opt,
		params.compute_also_old_base_pvals,
		params.n_most_freq_targets_for_stats != 0,
		params.n_most_freq_targets,
		params.opt_train_fraction,
		params.opt_num_inits,
		params.opt_num_iters,
		params.num_rand_cf,
		params.num_splits,
		params.cjs_out,
		params.max_pval_opt_for_Cjs
	);

	std::vector<uint64_t> top_anchors;

	while (true) {
		top_anchors.clear();

		if (!get_top_anchors(bins, anchor_filter, top_anchors))
			break;

		uint64_t min_anchor = *std::min_element(top_anchors.begin(), top_anchors.end());

		//merge all anchors that are min		
		std::vector<Non10SingleSampleAnchor> anchors;
		for (size_t i = 0; i < top_anchors.size(); ++i) {
			if (top_anchors[i] == min_anchor) {
				Non10SingleSampleAnchor anch;
				bins[i]->GetAnchor(anch);
				anch.data.shrink_to_fit();
				anchors.emplace_back(std::move(anch));
			}
		}
		uint64_t n_unique_targets{};
		uint64_t tot_cnt{};

		uint64_t n_unique_targets_kept{};
		uint64_t tot_cnt_kept{};

		//Anchor merged = merge_keep_target_order(anchors, n_unique_targets, tot_cnt);
		Anchor merged = merge_keep_target_order_binary_heap(anchors, params.n_most_freq_targets_for_stats, n_unique_targets, tot_cnt, n_unique_targets_kept, tot_cnt_kept);

		anchors.clear();
		anchors.shrink_to_fit();

		//consider better filtering
		std::unordered_set<uint64_t> unique_sample_ids;
		for (const auto& x : merged.data)
			unique_sample_ids.insert(pack_smaple_id_target(x.sample_id, 0)); //mkokot_TODO: for compatibility, not the most elegant solution

		if (anchor_filtered_out(tot_cnt_kept, n_unique_targets_kept, unique_sample_ids.size(), params)) {
			++stats.tot_filtered_out_anchors;
			continue;
		}

		++stats.tot_writen_anchors;
		stats.tot_writen_records += merged.data.size();

		anchor_processor->ProcessAnchor(std::move(merged), bin0_header.anchor_len_symbols, bin0_header.target_len_symbols, n_unique_targets_kept, tot_cnt_kept, n_unique_targets, tot_cnt, unique_sample_ids);

		if (stats.max_contignency_matrix_size.first * stats.max_contignency_matrix_size.second < unique_sample_ids.size() * n_unique_targets_kept) {
			stats.max_contignency_matrix_size.first = unique_sample_ids.size();
			stats.max_contignency_matrix_size.second = n_unique_targets_kept;
		}
	}
	stats.print(std::cerr);
}

int main(int argc, char** argv)
{
#ifdef _WIN32
	_setmaxstdio(2045);
#endif
	auto start_time = std::chrono::high_resolution_clock::now();

	auto params = read_params(argc, argv);
	params.print(std::cerr);

	run(params);

	std::chrono::duration<double> dur = (std::chrono::high_resolution_clock::now() - start_time);
	std::cerr << "Time: " << dur.count() << "s\n";
}