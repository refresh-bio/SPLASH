#ifndef _PVALS_H
#define _PVALS_H
#include <unordered_set>
#include "../common/satc_data.h"
#include "../common/cbc_to_cell_type.h"
#include "non_10X_supervised.h"
#include "matrix.h"
#include "anchor.h"



struct KmerAndCounter {
	uint64_t kmer;
	uint64_t counter;
};

struct AnchorStats {
	double pval_base;
	double pval_base_old;
	double pval_opt;
	double effect_size_cts;
	double effect_size_bin;
	double effect_size_bin_old;
	double pval_asymp_base;
	double entropy;
	double avg_no_homopolymer_targets;
	double avg_hamming_distance_max_target;
	double avg_hamming_distance_all_pairs;
	double avg_edit_distance_max_target;
	double avg_edit_distance_all_pairs;

	std::vector<KmerAndCounter> most_freq_targets;

	std::vector<double> helmert_decomposition_pvals;
	std::vector<double> helmert_decomposition_effect_size_cts;
	std::vector<double> helmert_decomposition_effect_size_bin;
	std::vector<double> Cjs_pvals; //non-10X supervised
	std::vector<double> Cjs_effect_size_cts; //non-10X supervised
	std::vector<double> Cjs_effect_size_bin; //non-10X supervised
	uint64_t Cjs_num_lt_0; //non-10X supervised
	uint64_t Cjs_num_rest; //non-10X supervised

	std::vector<uint32_t> cell_types_ids;  //unique ids of cell types for a given Xtrain contignency table

	AnchorStats(double pval_base, double pval_base_old, double pval_opt,
		double effect_size_cts, double effect_size_bin, double effect_size_bin_old, double pval_asymp_base, double entropy,
		double avg_no_homopolymer_targets,
		double avg_hamming_distance_max_target, double avg_hamming_distance_all_pairs,
		double avg_edit_distance_max_target, double avg_edit_distance_all_pairs
	) :
		pval_base(pval_base),
		pval_base_old(pval_base_old),
		pval_opt(pval_opt),
		effect_size_cts(effect_size_cts),
		effect_size_bin(effect_size_bin),
		effect_size_bin_old(effect_size_bin_old),
		pval_asymp_base(pval_asymp_base),
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
		effect_size_cts(0),
		effect_size_bin(0),
		effect_size_bin_old(0),
		pval_asymp_base(1),
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
	bool is_10X;
	uint32_t anchor_len;
	uint32_t barcode_len;
	bool enabled = false;
	SampleNameDecoder sample_name_decoder;
	std::ofstream out;
public:
	CjWriter(
		const std::string& path,
		bool is_10X,
		uint32_t anchor_len,
		uint32_t barcode_len,
		const std::string& sample_names) :
		is_10X(is_10X),
		anchor_len(anchor_len),
		barcode_len(barcode_len),
		sample_name_decoder(sample_names) {

		if (path == "")
			return;

		out.open(path);
		if (!out) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		enabled = true;
		out << "anchor\t";
		out << "sample\t";
		if (is_10X)
			out << "barcode\t";
		out << "Cj\n";
	}

	operator bool() const {
		return enabled;
	}

	void write(uint64_t anchor, uint64_t sample_id, uint64_t barcode, double Cj) {
		out << kmer_to_string(anchor, anchor_len) << "\t";
		//out << sample_id << "\t";
		sample_name_decoder.store_sample_id(out, sample_id);
		out << "\t";
		if (is_10X)
			out << kmer_to_string(barcode, barcode_len) << "\t";
		out << Cj << "\n";
	}
};


void compute_stats(
	Anchor&& anchor,
	size_t anchor_len_symbols,
	size_t n_uniq_targets,
	const std::unordered_set<uint64_t>& unique_samples,
	AnchorStats &anchor_stats,
	bool without_alt_max,
	bool with_effect_size_cts,
	bool compute_also_old_base_pvals,
	uint32_t n_most_freq_targets,
	double opt_train_fraction,
	int opt_num_inits,
	int opt_num_iters,
	size_t num_rand_cf,
	CjWriter& cj_writer,
	double max_pval_opt_for_Cjs,
	CBCToCellType* cbc_to_cell_type,
	Non10XSupervised* non_10X_supervised);

void dump_fXc(const refresh::matrix_1d<double>& f, const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& c);

#endif
