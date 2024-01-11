#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <memory>
#include <chrono>
#include <cblas.h>
#include "../common/version.h"
#include "pvals.h"
#include "extra_stats.h"
#include "ob_utils.h"
#include "../common/accepted_anchors.h"

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

//merge all samples from single bin
enum class RunMode { CalcStats, JustMerge, JustMergeAndDump };
std::string to_string(RunMode run_mode) {
	switch (run_mode)
	{
	case RunMode::CalcStats:
		return "calc_stats";
	case RunMode::JustMerge:
		return "just_merge";
	case RunMode::JustMergeAndDump:
		return "just_merge_and_dump";
	default:
		std::cerr << "Error: unsupported RunMode\n";
		exit(1);
		break;
	}
}

RunMode run_mode_from_string(const std::string& run_mode) {
	if (run_mode == "calc_stats")
		return RunMode::CalcStats;
	else if (run_mode == "just_merge")
		return RunMode::JustMerge;
	else if (run_mode == "just_merge_and_dump")
		return RunMode::JustMergeAndDump;
	else {
		std::cerr << "Error: unknown run mode " << run_mode << "\n";
		exit(1);
	}
}

struct Params
{
	uint64_t anchor_count_threshold{}; //remove all anchors for which the sum of counters across all samples and targets is <= anchor_count_threshold
	uint64_t anchor_unique_targets_threshold{}; // remove all anchors having <= anchor_unique_targets_threshold unique targets
	uint64_t anchor_samples_threshold{}; //remove all anchors having <= anchor_samples_threshold unique samples

	uint32_t n_most_freq_targets{}; //in stats mode print also n_most_freq_targets and their counts in output file

	double train_fraction = 0.25;

	int generate_alt_max_cf_no_tires = 10;

	int altMaximize_iters = 50;

	bool without_SVD = false;

	bool with_effect_size_cts = false;

	std::string outpath;

	std::string anchor_list;

	std::string cjs_out;

	double max_pval_rand_init_alt_max_for_Cjs = 0.1;

	std::string sample_names;

	std::vector<std::string> bins;

	RunMode run_mode = RunMode::CalcStats;
	RecFmt format = RecFmt::SATC; //only for JustMergeAndDump

	void print(std::ostream& oss) {
		oss << "Parameters:\n";
		oss << "\tanchor_count_threshold             : " << anchor_count_threshold << "\n";
		oss << "\tanchor_unique_targets_threshold    : " << anchor_unique_targets_threshold << "\n";
		oss << "\tanchor_samples_threshold           : " << anchor_samples_threshold << "\n";
		oss << "\twithout_SVD                        : " << std::boolalpha << without_SVD << "\n";
		oss << "\twith_effect_size_cts               : " << std::boolalpha << with_effect_size_cts << "\n";
		oss << "\toutpath                            : " << outpath << "\n";
		oss << "\tanchor_list                        : " << anchor_list << "\n";
		oss << "\tcjs_out                            : " << cjs_out << "\n";
		oss << "\tmax_pval_rand_init_alt_max_for_Cjs : " << max_pval_rand_init_alt_max_for_Cjs << "\n";
		oss << "\tsample_names	                     : " << sample_names << "\n";
		oss << "\tn_most_freq_targets                : " << n_most_freq_targets << "\n";
		oss << "\ttrain_fraction                     : " << train_fraction << "\n";
		oss << "\tgenerate_alt_max_cf_no_tires       : " << generate_alt_max_cf_no_tires << "\n";
		oss << "\taltMaximize_iters                  : " << altMaximize_iters << "\n";

		oss << "\trun_mode                           : " << to_string(run_mode) << "\n";
		if (run_mode == RunMode::JustMergeAndDump)
			oss << "\t\t format					       : " << RecFmtConv::to_string(format) << "\n";
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
			<< "    <outpath>                 - output path\n"
			<< "    <list_of_bins_to_merge>   - list of paths of bins to be processed\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_count_threshold <int>                - anchor_count_threshold\n" //mkokot_TODO: fix param desc
			<< "    --anchor_unique_targets_threshold <int>       - anchor_unique_targets_threshold\n" //mkokot_TODO: fix param desc
			<< "    --n_most_freq_targets <int>                   - in calc_stats mode output also n_most_freq_targets most frequent targets and their counts (default: 0)\n"
			<< "    --train_fraction <double> (0.0;1.0)           - in calc_stats mode use this fraction to create train X from contingency table\n"
			<< "    --generate_alt_max_cf_no_tires <int>          - in calc_stats mode the number of altMaximize runs (default: 10)\n"
			<< "    --altMaximize_iters <int>                     - in calc_stats mode the number of iteration in altMaximize (default: 50)\n"
			<< "    --anchor_samples_threshold <int>              - anchor_samples_threshold\n" //mkokot_TODO: fix param desc
			<< "    --anchor_list <string>                        - path to text file containing anchors separated by whitespaces, only anchors from this file will be processed\n"
			<< "    --cjs_out <string>                            - path to output text file where Cjs will be stored\n"
			<< "    --max_pval_rand_init_alt_max_for_Cjs <double> - dump only Cjs for anchors that have pval_rand_init_alt_max <= max_pval_rand_init_alt_max_for_Cjs\n"
			<< "    --run_mode <string>                           - what to do with merged anchors, available options: calc_stats, just_merge, just_merge_and_dump (default: calc_stats)\n"
			<< "    --without_SVD                                 - disable SVD computation\n"
			<< "    --with_effect_size_cts                        - compute with_effect_size_cts\n"
			<< "    --sample_names <path>                         - path for decode sample id, each line should contain <sample_name> <sample_id>\n"
			<< "    --format <string>                             - output format when 'just_merge_and_dump' is used as run_mode, available options: satc, splash (default: satc)\n";
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
		if (param == "--train_fraction") {
			std::string tmp = argv[++i];
			double val = std::stod(tmp);
			if (val <= 0.0 || val >= 1.0) {
				std::cerr << "Error: train_fraction must be in range (0.0;1.0)\n";
				exit(1);
			}
			res.train_fraction = val;
		}
		if (param == "--generate_alt_max_cf_no_tires") {
			std::string tmp = argv[++i];
			res.generate_alt_max_cf_no_tires = std::stoull(tmp);
		}
		if (param == "--altMaximize_iters") {
			std::string tmp = argv[++i];
			res.altMaximize_iters = std::stoull(tmp);
		}
		if (param == "--without_SVD") {
			res.without_SVD = true;
		}
		if (param == "--with_effect_size_cts") {
			res.with_effect_size_cts = true;
		}
		if (param == "--anchor_list") {
			res.anchor_list = argv[++i];
		}
		if (param == "--cjs_out") {
			res.cjs_out = argv[++i];
		}
		if (param == "--max_pval_rand_init_alt_max_for_Cjs") {
			std::string tmp = argv[++i];
			res.max_pval_rand_init_alt_max_for_Cjs = std::stod(tmp);
		}
		if (param == "--sample_names") {
			res.sample_names = argv[++i];
		}
		if (param == "--run_mode") {
			res.run_mode = run_mode_from_string(argv[++i]);
		}
		if (param == "--format")
			res.format = RecFmtConv::from_string(argv[++i]);
	}
	if (i >= argc) {
		std::cerr << "Error: outpath missing\n";
		exit(1);
	}

	res.outpath = argv[i++];

	while (i < argc)
		res.bins.emplace_back(argv[i++]);

	return res;
}

Anchor merge_keep_target_order(const std::vector<Anchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt) {
	assert(to_merge.size());
	n_unique_targets = 0;
	tot_cnt = 0;
	Anchor res;
	res.anchor = to_merge[0].anchor;
	std::vector<uint64_t> read_pos(to_merge.size());

	size_t records_left = 0;
	for (const auto& x : to_merge)
		records_left += x.data.size();

	while (records_left) {
		uint64_t min_target = std::numeric_limits<uint64_t>::max();

		//find min target
		for (size_t i = 0; i < to_merge.size(); ++i) {
			if (read_pos[i] < to_merge[i].data.size()) {
				if (to_merge[i].data[read_pos[i]].target < min_target) {
					min_target = to_merge[i].data[read_pos[i]].target;
				}
			}
		}

		++n_unique_targets;

		//store all record having min target in the result
		for (size_t i = 0; i < to_merge.size(); ++i) {
			if (read_pos[i] < to_merge[i].data.size()) {
				if (to_merge[i].data[read_pos[i]].target == min_target) {
					res.data.push_back(to_merge[i].data[read_pos[i]]);
					tot_cnt += res.data.back().count;

					++read_pos[i];
					--records_left;
				}
			}
		}
	}
	res.data.shrink_to_fit();
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
	Anchor current_anchor;
	bool is_loaded = false;
	CachedRecord cached_rec;
	bool load_anchor() {
		current_anchor.data.clear();
		Record rec;

		if (!cached_rec.Peek(header, rec))
			return false;
		cached_rec.Skip();

		current_anchor.anchor = rec.anchor;
		current_anchor.data.emplace_back(rec.barcode, rec.target, rec.sample_id, rec.count);
		while (cached_rec.Peek(header, rec)) {
			if (rec.anchor == current_anchor.anchor) {
				current_anchor.data.emplace_back(rec.barcode, rec.target, rec.sample_id, rec.count);
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

	void GetAnchor(Anchor& anchor) {
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


void write_out_header(std::ofstream& out, bool without_SVD, bool with_effect_size_cts, uint32_t n_most_freq_targets) {
	//mkokot_TODO: czy doimplementowac zakomentowane?
	out
		<< "anchor" << "\t";
	if (!without_SVD)
		out << "pval_SVD_corrAnalysis" << "\t";				// TODO: ?
	out
		<< "pval_rand_init_alt_max" << "\t";
	if (with_effect_size_cts)
		out << "effect_size_cts" << "\t";
	out
		<< "effect_size_bin" << "\t"
		//<< "effect_size_cts_SVD" << "\t"					// TODO: ?
		//<< "pval_rand_init_EM" << "\t"					// TODO: ?
		<< "pval_base" << "\t"
		//<< "effect_size_base" << "\t"						// TODO: ?
		<< "M" << "\t"									// X.sum()
		<< "anch_uniqTargs" << "\t"							// (X.sum(axis=1)>0).sum()
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

void write_out_rec(std::ofstream& out, const AnchorStats& anchor_stats, uint64_t anchor, size_t anchor_len_symbols, size_t target_len_symbols, uint64_t tot_cnt, uint64_t n_unique_targets, uint64_t n_unique_samples, bool without_SVD, bool with_effect_size_cts, uint32_t n_most_freq_targets) {
	out
		<< kmer_to_string(anchor, anchor_len_symbols) << "\t";
	if (!without_SVD)
		out
		<< anchor_stats.pval_SVD_corrAnalysis << "\t";

	out
		<< anchor_stats.pval_rand_init_alt_max << "\t";
	if (with_effect_size_cts)
		out << anchor_stats.effect_size_cts << "\t";
	out
		<< anchor_stats.effect_size_bin << "\t"
		//<< "effect_size_cts_SVD" << "\t"					// TODO: ?
		//<< "pval_rand_init_EM" << "\t"					// TODO: ?
		<< anchor_stats.pval_base << "\t"
		//<< "effect_size_base" << "\t"						// TODO: ?
		<< tot_cnt << "\t"									// X.sum()
		<< n_unique_targets << "\t"							// (X.sum(axis=1)>0).sum()
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
		size_t n_uniq_targets,
		size_t tot_cnt,
		const std::unordered_set<uint64_t>& unique_samples) = 0;

	virtual ~IAnchorProcessor() = default;
};

class StatsWriter : public IAnchorProcessor {
	std::ofstream out;
	bool without_SVD;
	bool with_effect_size_cts;
	uint32_t n_most_freq_targets;
	double train_fraction;
	int generate_alt_max_cf_no_tires;
	int altMaximize_iters;
	CExtraStats extra_stats;
	AnchorStats anchor_stats;
	CjWriter cj_writer;
	double max_pval_rand_init_alt_max_for_Cjs;
public:
	StatsWriter(
		const std::string& outpath,
		bool without_SVD,
		bool with_effect_size_cts,
		uint32_t n_most_freq_targets,
		double train_fraction,
		int generate_alt_max_cf_no_tires,
		int altMaximize_iters,
		const std::string& cjs_out,
		size_t anchor_len_symbols,
		size_t barcode_len_symbols,
		const std::string& sample_names,
		double max_pval_rand_init_alt_max_for_Cjs) :
		out(outpath),
		without_SVD(without_SVD),
		with_effect_size_cts(with_effect_size_cts),
		n_most_freq_targets(n_most_freq_targets),
		train_fraction(train_fraction),
		generate_alt_max_cf_no_tires(generate_alt_max_cf_no_tires),
		altMaximize_iters(altMaximize_iters),
		cj_writer(cjs_out, anchor_len_symbols, barcode_len_symbols, sample_names),
		max_pval_rand_init_alt_max_for_Cjs(max_pval_rand_init_alt_max_for_Cjs) {
		if (!out) {
			std::cerr << "Error: cannot open file " << outpath << "\n";
			exit(1);
		}

		write_out_header(out, without_SVD, with_effect_size_cts, n_most_freq_targets);
	}

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniq_targets,
		size_t tot_cnt,
		const std::unordered_set<uint64_t>& unique_samples
	)  override {

		auto _anchor = anchor.anchor;

		extra_stats.Compute(anchor, target_len_symbols, 5, 200, n_uniq_targets, unique_samples, anchor_stats);

		compute_stats(std::move(anchor), anchor_len_symbols, n_uniq_targets, unique_samples, anchor_stats, without_SVD, with_effect_size_cts, n_most_freq_targets, train_fraction, generate_alt_max_cf_no_tires, altMaximize_iters, cj_writer, max_pval_rand_init_alt_max_for_Cjs);

		write_out_rec(out, anchor_stats, _anchor, anchor_len_symbols, target_len_symbols, tot_cnt, n_uniq_targets, unique_samples.size(), without_SVD, with_effect_size_cts, n_most_freq_targets);
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

	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniq_targets,
		size_t tot_cnt,
		const std::unordered_set<uint64_t>& unique_samples
	)  override {
		Record rec;
		rec.anchor = anchor.anchor;

		for (const auto& x : anchor.data) {
			rec.barcode = x.barcode;
			rec.count = x.count;
			rec.sample_id = x.sample_id;
			rec.target = x.target;

			rec.serialize(out, out_header);
		}

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
	void ProcessAnchor(
		Anchor&& anchor,
		size_t anchor_len_symbols,
		size_t target_len_symbols,
		size_t n_uniq_targets,
		size_t tot_cnt,
		const std::unordered_set<uint64_t>& unique_samples
	) override {
		Record rec;
		rec.anchor = anchor.anchor;

		for (const auto& x : anchor.data) {
			rec.barcode = x.barcode;
			rec.count = x.count;
			rec.sample_id = x.sample_id;
			rec.target = x.target;

			rec.print(out, header, format, sample_name_decoder);
		}

		anchor.data.clear();
		anchor.data.shrink_to_fit();
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
	RunMode run_mode,
	const std::string& outpath,
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
	bool without_SVD,
	bool with_effect_size_cts,
	uint32_t n_most_freq_targets,
	double train_fraction,
	int generate_alt_max_cf_no_tires,
	int altMaximize_iters,
	const std::string& cjs_out,
	double max_pval_rand_init_alt_max_for_Cjs) {
	switch (run_mode)
	{
	case RunMode::JustMerge:
		return std::make_unique<SatcWriter>(
			outpath,
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
	case RunMode::JustMergeAndDump:
		return std::make_unique<SatcDumpWriter>(
			outpath,
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
	case RunMode::CalcStats:
		return std::make_unique<StatsWriter>(
			outpath,
			without_SVD,
			with_effect_size_cts,
			n_most_freq_targets,
			train_fraction,
			generate_alt_max_cf_no_tires,
			altMaximize_iters,
			cjs_out,
			anchor_len_symbols,
			barcode_len_symbols,
			sample_names,
			max_pval_rand_init_alt_max_for_Cjs);
		break;
	default:
		std::cerr << "Error: this is not supported, mode " << to_string(run_mode) << "\n";
		exit(1);
	}
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

	std::unique_ptr<IAnchorProcessor> anchor_processor = get_anchor_processor(params.run_mode,
		params.outpath,
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
		params.without_SVD,
		params.with_effect_size_cts,
		params.n_most_freq_targets,
		params.train_fraction,
		params.generate_alt_max_cf_no_tires,
		params.altMaximize_iters,
		params.cjs_out,
		params.max_pval_rand_init_alt_max_for_Cjs
	);

	std::vector<uint64_t> top_anchors;

	while (true) {
		top_anchors.clear();

		if (!get_top_anchors(bins, anchor_filter, top_anchors))
			break;

		uint64_t min_anchor = *std::min_element(top_anchors.begin(), top_anchors.end());

		//merge all anchors that are min		
		std::vector<Anchor> anchors;
		for (size_t i = 0; i < top_anchors.size(); ++i) {
			if (top_anchors[i] == min_anchor) {
				Anchor anch;
				bins[i]->GetAnchor(anch);
				anch.data.shrink_to_fit();
				anchors.emplace_back(std::move(anch));
			}
		}
		uint64_t n_unique_targets{};
		uint64_t tot_cnt{};
		Anchor merged = merge_keep_target_order(anchors, n_unique_targets, tot_cnt);

		//consider better filtering
		std::unordered_set<uint64_t> unique_sample_ids;
		for (const auto& x : merged.data)
			unique_sample_ids.insert(pack_smaple_id_target(x.sample_id, 0)); //mkokot_TODO: for compatibility, not the most elegant solution

		if (anchor_filtered_out(tot_cnt, n_unique_targets, unique_sample_ids.size(), params)) {
			++stats.tot_filtered_out_anchors;
			continue;
		}

		++stats.tot_writen_anchors;
		stats.tot_writen_records += merged.data.size();

		anchor_processor->ProcessAnchor(std::move(merged), bin0_header.anchor_len_symbols, bin0_header.target_len_symbols, n_unique_targets, tot_cnt, unique_sample_ids);

		if (stats.max_contignency_matrix_size.first * stats.max_contignency_matrix_size.second < unique_sample_ids.size() * n_unique_targets) {
			stats.max_contignency_matrix_size.first = unique_sample_ids.size();
			stats.max_contignency_matrix_size.second = n_unique_targets;
		}
	}
	stats.print(std::cerr);
}

int main(int argc, char** argv)
{
#ifdef _WIN32
	_setmaxstdio(2045);
#endif
#ifdef USE_OPEN_BLAS
	openblas_set_num_threads(1);
#endif
	auto start_time = std::chrono::high_resolution_clock::now();

	auto params = read_params(argc, argv);
	params.print(std::cerr);

	run(params);

	std::chrono::duration<double> dur = (std::chrono::high_resolution_clock::now() - start_time);
	std::cerr << "Time: " << dur.count() << "s\n";
}