#ifndef _PVALS_H
#define _PVALS_H
#include <cstdio>
#include <unordered_set>
#include "../common/types/satc_data.h"
#include "../common/cbc_to_cell_type.h"
#include "matrix.h"
#include "non_10X_supervised.h"
#include "matrix.h"
#include "anchor.h"

#include <refresh/compression/lib/gz_wrapper.h>
#include <refresh/conversions/lib/conversions.h>

struct KmerAndCounter {
	uint64_t kmer;
	uint64_t counter;
};

struct AnchorStats {
	double pval_base;
	double pval_base_old;
	double pval_opt;
	double pval_sample_spectral_sum;
	double effect_size_cts;
	double effect_size_bin;
	double effect_size_bin_old;
	double pval_asymp_opt;
	double entropy;
	double avg_no_homopolymer_targets;
	double avg_hamming_distance_max_target;
	double avg_hamming_distance_all_pairs;
	double avg_edit_distance_max_target;
	double avg_edit_distance_all_pairs;

	std::vector<KmerAndCounter> most_freq_targets;
	std::vector<std::pair<float, float>> sequence_entropy_targets;
	std::pair<float, float> sequence_entropy_anchor;

	std::vector<double> helmert_decomposition_pvals;
	std::vector<double> helmert_decomposition_effect_size_cts;
	std::vector<double> helmert_decomposition_effect_size_bin;
	std::vector<double> Cjs_pvals; //non-10X supervised
	std::vector<double> Cjs_effect_size_cts; //non-10X supervised
	std::vector<double> Cjs_effect_size_bin; //non-10X supervised
	uint64_t Cjs_num_lt_0; //non-10X supervised
	uint64_t Cjs_num_rest; //non-10X supervised

	std::vector<uint32_t> cell_types_ids;  //unique ids of cell types for a given Xtrain contignency table

	AnchorStats(double pval_base, double pval_base_old, double pval_opt, double pval_sample_spectral_sum,
		double effect_size_cts, double effect_size_bin, double effect_size_bin_old, double pval_asymp_opt, double entropy,
		double avg_no_homopolymer_targets,
		double avg_hamming_distance_max_target, double avg_hamming_distance_all_pairs,
		double avg_edit_distance_max_target, double avg_edit_distance_all_pairs
	) :
		pval_base(pval_base),
		pval_base_old(pval_base_old),
		pval_opt(pval_opt),
		pval_sample_spectral_sum(pval_sample_spectral_sum),
		effect_size_cts(effect_size_cts),
		effect_size_bin(effect_size_bin),
		effect_size_bin_old(effect_size_bin_old),
		pval_asymp_opt(pval_asymp_opt),
		entropy(entropy),
		avg_no_homopolymer_targets(avg_no_homopolymer_targets),
		avg_hamming_distance_max_target(avg_hamming_distance_max_target),
		avg_hamming_distance_all_pairs(avg_hamming_distance_all_pairs),
		avg_edit_distance_max_target(avg_edit_distance_max_target),
		avg_edit_distance_all_pairs(avg_edit_distance_all_pairs)
	{}

	AnchorStats() :
		pval_base(1),
		pval_base_old(1),
		pval_opt(1),
		pval_sample_spectral_sum(1),
		effect_size_cts(0),
		effect_size_bin(0),
		effect_size_bin_old(0),
		pval_asymp_opt(1),
		entropy(0),
		avg_no_homopolymer_targets(0),
		avg_hamming_distance_max_target(0),
		avg_hamming_distance_all_pairs(0),
		avg_edit_distance_max_target(0),
		avg_edit_distance_all_pairs(0)
	{}

	void clear()
	{
		pval_base = 1;
		pval_base_old = 1;
		pval_opt = 1;
		pval_sample_spectral_sum = 1;
		
		clear_extra();
	}

	void clear_extra()
	{
		entropy = 0;
		avg_no_homopolymer_targets = 0;
		avg_hamming_distance_max_target = 0;
		avg_hamming_distance_all_pairs = 0;
		avg_edit_distance_max_target = 0;
		avg_edit_distance_all_pairs = 0;
	}
};


class Calculator_S
{
//	const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major> &X;
	const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> &X;
	refresh::matrix_1d<double> f;
	refresh::matrix_1d<double> c;

	refresh::matrix_1d<double> X_row_sum;
	refresh::matrix_1d<double> X_col_sum;
	refresh::matrix_1d<double> X_col_sum_sqrt;
	double M;

	// Already precomputed
	refresh::matrix_1d<double> pre_fT_Xtild;
	refresh::matrix_1d<double> pre_Xtild_c;
	std::pair<double, bool> pre_S;

	void prepare_X_precalc()
	{
		X.get_row_col_sums(X_row_sum, X_col_sum);

		M = X_row_sum.sum();

		X_col_sum_sqrt = sqrt(X_col_sum);
	}

public:
	Calculator_S() = delete;

//	Calculator_S(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X) : X(X)
	Calculator_S(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X) : X(X)
	{
		prepare_X_precalc();

		pre_S.second = false;
	}

//	Calculator_S(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X,
	Calculator_S(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X,
		const refresh::matrix_1d<double> &f, const refresh::matrix_1d<double> &c) : X(X), f(f), c(c)
	{
		prepare_X_precalc();
	}

	void set_f(const refresh::matrix_1d<double>& f_)
	{
		f = f_;

		pre_S.second = false;
		pre_fT_Xtild.clear();
	}

	void clear_f()
	{
		f.clear();

		pre_S.second = false;
		pre_fT_Xtild.clear();
	}

	void move_f(refresh::matrix_1d<double>&& f_)
	{
		f = std::move(f_);

		pre_S.second = false;
		pre_fT_Xtild.clear();
	}

	void set_c(const refresh::matrix_1d<double>& c_)
	{
		c = c_;

		pre_S.second = false;
		pre_Xtild_c.clear();
	}

	void clear_c()
	{
		c.clear();

		pre_S.second = false;
		pre_Xtild_c.clear();
	}

	void move_c(refresh::matrix_1d<double>&& c_)
	{
		c = std::move(c_);

		pre_S.second = false;
		pre_Xtild_c.clear();
	}

	refresh::matrix_1d<double> get_c()
	{
		return c;
	}

	refresh::matrix_1d<double> get_f()
	{
		return f;
	}

	refresh::matrix_1d<double>& get_X_row_sum()
	{
		return X_row_sum;
	}

	refresh::matrix_1d<double>& get_X_col_sum()
	{
		return X_col_sum;
	}

	double get_M()
	{
		return M;
	}

	refresh::matrix_1d<double> mult_fT_Xtild()
	{
//		if (pre_fT_Xtild.empty())
		{
			// Calculate "base" values
			pre_fT_Xtild = (-dot_product(f, X_row_sum) / M) * X_col_sum_sqrt;

			// Add non-zero X values
			for (auto p = X.begin(); p != X.end(); ++p)
				pre_fT_Xtild(p->first.col) += f(p->first.row) * p->second / X_col_sum_sqrt(p->first.col);
		}

		return pre_fT_Xtild;
	}

	refresh::matrix_1d<double> mult_Xtild_c()
	{
//		if (pre_Xtild_c.empty())
		{
			// Calculate "base" values
			pre_Xtild_c = X_row_sum * -(dot_product(X_col_sum_sqrt, c) / M);

			// Add non-zero X values
			for (auto p = X.begin(); p != X.end(); ++p)
				pre_Xtild_c(p->first.row) += c(p->first.col) * p->second / X_col_sum_sqrt(p->first.col);
		}

		return pre_Xtild_c;
	}

	double mult_fT_Xtild_c()
	{
/*		if (pre_S.second)
			return pre_S.first;

		if(!pre_fT_Xtild.empty())*/
			pre_S.first = dot_product(mult_fT_Xtild(), c);
/*		else
			pre_S.first = dot_product(f, mult_Xtild_c());		

		pre_S.second = true;*/

		return pre_S.first;
	}
};

class CjWriter {
	bool _10X_or_visium;
	uint32_t anchor_len;
	uint32_t barcode_len;
	bool enabled = false;
	SampleNameDecoder sample_name_decoder;

	FILE* out = nullptr;

	const size_t buffer_size = 16 << 20;
	const size_t max_line_len = 128;

	char* buffer;
	char* compressed_buffer;
	size_t compressed_buffer_size;
	size_t in_buf_pos;

	refresh::gz_in_memory giz{ 6, false };

	void store_buffer()
	{
		size_t compressed_size = giz.compress(buffer, in_buf_pos, compressed_buffer, compressed_buffer_size);

		fwrite(compressed_buffer, 1, compressed_size, out);

		in_buf_pos = 0;
	}

	void add_to_buffer(const char* p)
	{
		strcpy(buffer + in_buf_pos, p);
		in_buf_pos += strlen(p);
	}

public:
	CjWriter(
		const std::string& path,
		bool without_header,
		bool _10X_or_visium,
		uint32_t anchor_len,
		uint32_t barcode_len,
		const std::string& sample_names) :
		_10X_or_visium(_10X_or_visium),
		anchor_len(anchor_len),
		barcode_len(barcode_len),
		sample_name_decoder(sample_names) {

		if (path == "")
			return;

		out = fopen(path.c_str(), "wb");

		if (!out) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		setvbuf(out, nullptr, _IOFBF, 16 << 20);

		enabled = true;

		buffer = new char[buffer_size];
		compressed_buffer_size = buffer_size + giz.get_overhead(buffer_size);
		compressed_buffer = new char[compressed_buffer_size];
		in_buf_pos = 0;

		if (!without_header)
		{
			add_to_buffer("anchor\t");
			add_to_buffer("sample\t");
			if (_10X_or_visium)
				add_to_buffer("barcode\t");
			add_to_buffer("Cj\n");
		}
	}

	~CjWriter()
	{
		if (!out)
			return;

		if (in_buf_pos)
			store_buffer();

		fclose(out);

		delete[] buffer;
		delete[] compressed_buffer;
	}

	operator bool() const {
		return enabled;
	}

	void write(uint64_t anchor, uint64_t sample_id, uint64_t barcode, double Cj) {
		if (in_buf_pos + max_line_len >= buffer_size)
			store_buffer();

		in_buf_pos += refresh::kmer_to_pchar(anchor, buffer + in_buf_pos, anchor_len, false, '\t');

		in_buf_pos += sample_name_decoder.store_sample_id(buffer + in_buf_pos, sample_id);
		buffer[in_buf_pos++] = '\t';

		if (_10X_or_visium)
			in_buf_pos += refresh::kmer_to_pchar(barcode, buffer + in_buf_pos, barcode_len, false, '\t');

		in_buf_pos += refresh::real_to_pchar(Cj, buffer + in_buf_pos, 8, '\n');
	}
};


void compute_stats(
	Anchor&& anchor,
	size_t anchor_len_symbols,
	size_t target_len_symbols,
	size_t n_uniq_targets,
	const std::unordered_set<uint64_t>& unique_samples,
	AnchorStats &anchor_stats,
	bool without_alt_max,
	bool without_sample_spectral_sum,
	bool with_effect_size_cts,
	bool with_pval_asymp_opt,
	bool compute_also_old_base_pvals,
	uint32_t n_most_freq_targets,
	double opt_train_fraction,
	int opt_num_inits,
	int opt_num_iters,
	size_t num_rand_cf,
	size_t num_splits,
	CjWriter& cj_writer,
	double max_pval_opt_for_Cjs,
	CBCToCellType* cbc_to_cell_type,
	Non10XSupervised* non_10X_supervised);

void dump_fXc(const refresh::matrix_1d<double>& f, const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& c);
std::pair<double, double> sample_spectral_sum(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& table);


#endif
