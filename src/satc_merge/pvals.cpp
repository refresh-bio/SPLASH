#include <random>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>

#include "pvals.h"
#include "get_train_mtx.h"

class SampleToColMapper {
	using IndexType = size_t;
	//std::vector<IndexType> mapping;
	std::unordered_map<uint64_t, IndexType> mapping;
	std::vector<uint64_t> decode_mapping;
public:
	SampleToColMapper(const std::unordered_set<uint64_t>& unique_samples) {
		IndexType idx = 0;

		//mkokot_TODO: restore this version
		//for (auto sample_id : unique_samples)
		//	mapping[sample_id] = idx++;

		//mkokot_TODO: remove this version
		{
			std::vector<uint64_t> sorted_samples(unique_samples.begin(), unique_samples.end());
			std::sort(sorted_samples.begin(), sorted_samples.end());

			decode_mapping.resize(unique_samples.size());

			for (auto sample_id : sorted_samples) {
				decode_mapping[idx] = sample_id;
				mapping[sample_id] = idx++;
			}
		}

		//for vector based indexing
		//size_t max_sample_id{};
		//for (auto sample_id : unique_samples)
		//	if (sample_id > max_sample_id)
		//		max_sample_id = sample_id;

		//mapping.resize(max_sample_id + 1, std::numeric_limits<IndexType>::max());

		//std::vector<uint64_t> sorted_samples(unique_samples.begin(), unique_samples.end());
		//std::sort(sorted_samples.begin(), sorted_samples.end());
		//
		//IndexType idx = 0;
		//for (auto sample : sorted_samples)
		//	mapping[sample] = idx++;
	}
	IndexType map(uint64_t sample_id) const {
		//return mapping[sample_id]; //cannot use because of `const`
		return mapping.find(sample_id)->second;
	}

	uint64_t decode(IndexType index) const {
		return decode_mapping[index];
	}
};

bool all_values_the_same(const refresh::matrix_1d<double>& vec) {
	for (size_t i = 1; i < vec.size(); ++i) {
		if (vec(i) != vec(0))
			return false;
	}
	return true;
}

//refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major> remove_zero_cols_rows(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& mtx, std::vector<uint64_t>& non_zero_cols, std::vector<uint64_t>& non_zero_rows)
refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> 
remove_zero_cols_rows(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& mtx, std::vector<uint64_t>& non_zero_cols, std::vector<uint64_t>& non_zero_rows)
{
	auto sums_in_cols = mtx.get_col_sums();
	auto sums_in_rows = mtx.get_row_sums();

	non_zero_cols.clear();
	for (size_t col = 0; col < sums_in_cols.size(); ++col)
		if (sums_in_cols(col))
			non_zero_cols.push_back(col);

	non_zero_rows.clear();
	for (size_t row = 0; row < sums_in_rows.size(); ++row)
		if (sums_in_rows(row))
			non_zero_rows.push_back(row);

	return mtx.compact(non_zero_rows, non_zero_cols);
}

refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>
remove_zero_rows(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& mtx, std::vector<uint64_t>& non_zero_rows)
{
	auto sums_in_rows = mtx.get_row_sums();

	non_zero_rows.clear();
	for (size_t row = 0; row < sums_in_rows.size(); ++row)
		if (sums_in_rows(row))
			non_zero_rows.push_back(row);

	return mtx.compact(non_zero_rows);
}

void normalize_vec_0_1(refresh::matrix_1d<double>& vec) 
{
	auto _min = vec.min_coeff();
	auto _max = vec.max_coeff();
	if (_min == _max) {
		vec.set_to_zero();
		return;
	}

	double dif = _max - _min;

	for (size_t i = 0; i < vec.size(); ++i)
		vec(i) = (vec(i) - _min) / dif;
}

double effectSize_cts(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {
	auto f_abs_sum = (X * abs(cOpt)).sum();

	if (f_abs_sum == 0)
		return 0;

	return fabs(refresh::dot_product(fOpt * X, cOpt) / f_abs_sum);
}

//double effectSize_bin(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {
double effectSize_bin(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {
	auto pos_vec = cOpt > 0;
	auto neg_vec = cOpt < 0;

	if (std::none_of(pos_vec.begin(), pos_vec.end(), [](auto x) {return x; }) || std::none_of(neg_vec.begin(), neg_vec.end(), [](auto x) {return x; }))
		return 0;

#if 0		// Dump for tests
	cerr << "X = np.zeros((" << X.rows() << "," << X.cols() << "))\n";

	for (auto e : X)
		cerr << "X[" << e.first.row << "][" << e.first.col << "] = " << e.second << endl;

	cerr << "fOpt = np.zeros((1," << fOpt.size() << "))\n";
	for (size_t i = 0; i < fOpt.size(); ++i)
		if (fOpt(i) != 0)
			cerr << "fOpt[0][" << i << "] = " << fOpt(i) << "\n";

	cerr << "cOpt = np.zeros((" << cOpt.size() << ", 1))\n";
	for (size_t i = 0; i < cOpt.size(); ++i)
		if (cOpt(i) != 0)
			cerr << "cOpt[" << i << "][0] = " << cOpt(i) << "\n";
#endif
	refresh::matrix_1d<double> cPos(pos_vec.begin(), pos_vec.end());
	refresh::matrix_1d<double> cNeg(neg_vec.begin(), neg_vec.end());

	return fabs(refresh::dot_product(fOpt * X, cPos) / (X * cPos).sum() - refresh::dot_product(fOpt * X, cNeg) / (X * cNeg).sum());
}

double computeAsympSPLASH(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt)
{
	if (cOpt.all_items_same())
		return 1.0;

	double min_f;
	double max_f;

	std::tie(min_f, max_f) = fOpt.min_max_coeff();

//	dump_fXc(fOpt, X, cOpt);
//	if (X.rows() > 50000000)

/*	if (min_f == max_f)		// all items same
		return 1.0;*/

	if (!(max_f <= 1 && min_f >= 0)) 
	{
		std::cerr << "Error: if !(fOpt.maxCoeff() <= 1 && fOpt.minCoeff() >= 0) not supported\n";		// !!! TODO: maybe normalization would be necessary
		exit(1);
	}

	double S, M;

	{
		Calculator_S calc_S(X, fOpt, cOpt);

		S = calc_S.mult_fT_Xtild_c();
		M = calc_S.get_M();
	}

	double muhat = (fOpt * X).sum() / M;

	double varF = dot_product((fOpt - muhat).pow2(), X.get_row_sums()) / M;

	double totalVar = varF * (pow(cOpt.norm(), 2) - pow(dot_product(cOpt, X.get_col_sums().sqrt()), 2) / M);

	if (totalVar <= 0)
		return 1.0;

	double normalizedTestStat = S / sqrt(totalVar);

	// cdf = std::erfc(-x / std::sqrt(2)) / 2;

	double pval = 2 * erfc(abs(normalizedTestStat) / sqrt(2)) / 2;

//	cout << pval << endl;

	return pval;
}

double testPvalOld(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {

//	dump_fXc(fOpt, X, cOpt);

//	if (cOpt.all_items_same() || fOpt.all_items_same())
	if (cOpt.all_items_same())
		return 1.0;

	double min_f;
	double max_f;

	std::tie(min_f, max_f) = fOpt.min_max_coeff();

	if (min_f == max_f)		// all items same
		return 1.0;

	if (!(max_f <= 1 && min_f >= 0)) {
		std::cerr << "Error: if !(fOpt.maxCoeff() <= 1 && fOpt.minCoeff() >= 0) not supported\n";
		exit(1);
	}

	Calculator_S calc_S(X, fOpt, cOpt);

	double S = calc_S.mult_fT_Xtild_c();

	double gamma = 0;
	auto& X_col_sum = calc_S.get_X_col_sum();

	for (size_t i = 0; i < cOpt.size(); ++i)
		gamma += cOpt(i) * sqrt(X_col_sum(i));

	double cOpt_norm = cOpt.norm();

	gamma /= cOpt_norm;
	gamma /= sqrt(calc_S.get_M());

	gamma *= gamma;

	double xi = 0;
	if (gamma > 0)
		xi = 1.0 / (1.0 + 1 / sqrt(gamma));

	double pval = 2 * exp(-2.0 * (1 - xi) * S * S / (cOpt_norm * cOpt_norm));

	if (gamma > 0)
		pval *= 2;

	return pval < 1 ? pval : 1;
}

//double testPval(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {
double testPval(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {

//	if (cOpt.all_items_same() || fOpt.all_items_same())
	if (cOpt.all_items_same())
		return 1.0;

	double min_f;
	double max_f;

	std::tie(min_f, max_f) = fOpt.min_max_coeff();
	
	if (min_f == max_f)		// all items same
		return 1.0;

	if (!(max_f <= 1 && min_f >= 0)) {
		std::cerr << "Error: if !(fOpt.maxCoeff() <= 1 && fOpt.minCoeff() >= 0) not supported\n";
		exit(1);
	}

	Calculator_S calc_S(X, fOpt, cOpt);

	double S = calc_S.mult_fT_Xtild_c();
	
	double denom = pow2(cOpt).sum();

	double denom_part2 = 0;

	auto& X_col_sum = calc_S.get_X_col_sum();
	for (size_t i = 0; i < cOpt.size(); ++i)
		denom_part2 += cOpt(i) * sqrt(X_col_sum(i));

	denom -= denom_part2 * denom_part2 / calc_S.get_M();
	
	auto pval = 2 * exp(-2 * S * S / denom);

	return pval < 1 ? pval : 1;
}

//double altMaximize(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X, refresh::matrix_1d<double>& c, refresh::matrix_1d<double>& f, int n_iters)
double altMaximize(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, refresh::matrix_1d<double>& c, refresh::matrix_1d<double>& f, int n_iters)
{
	if (c.all_items_same())
		c(0) = -c(0);

	Calculator_S calc_S(X, f, c);

	double S = 0;
	double Sold = 0;

	for (int i = 0; i < n_iters; ++i)
	{
		f.clear();

		auto f1 = sign_01(calc_S.mult_Xtild_c());
//		auto f2 = 1.0 - f1;

		calc_S.set_f(f1);
		double v1 = abs(calc_S.mult_fT_Xtild_c());

		f1.negate_sign_01();
//		calc_S.set_f(f2);
		calc_S.set_f(f1);
		double v2 = abs(calc_S.mult_fT_Xtild_c());

		if (v1 >= v2)
		{
			f1.negate_sign_01();
			f = std::move(f1);
		}
		else
		{
//			f = move(f2);
			f = std::move(f1);
		}
		
		calc_S.set_f(f);

		c = calc_S.mult_fT_Xtild();

		auto norm_c = c.norm();
		if (norm_c > 0)
			c /= norm_c;

		calc_S.set_c(c);

		S = calc_S.mult_fT_Xtild_c();

		if (S == Sold)
			return abs(S);

		if (i == n_iters - 1 && S != 0 && abs((S - Sold) / S) < 0.001)
			return abs(S);

		Sold = S;
	}

	c.set_to_zero();
	f.set_to_zero();

	return 0;
}

//void generate_alt_max_cf(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt, int no_tries, int opt_num_iters)
void generate_alt_max_cf(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt, int no_tries, int opt_num_iters)
{
	std::default_random_engine eng;
	std::uniform_int_distribution<int> dist_0_1(0, 1);
	double minus_1_and_1[] = { -1, 1 };

	std::vector<uint64_t> non_zero_cols;
	std::vector<uint64_t> non_zero_rows;

	auto X = remove_zero_cols_rows(X_org, non_zero_cols, non_zero_rows);

	if (X.cols() < 2 || X.rows() < 2)
	{
		cOpt.resize(X_org.cols(), 0);
		fOpt.resize(X_org.rows(), 0);

		return;
	}

	refresh::matrix_1d<double> fMax, cMax;
	double Sbase = 0;

	refresh::matrix_1d<double> c(X.cols());
	refresh::matrix_1d<double> f;

	for (int i = 0; i < no_tries; ++i)
	{
		std::generate(c.data(), c.data() + c.size(), [&] {
			return minus_1_and_1[dist_0_1(eng)];
			});

		double S = altMaximize(X, c, f, opt_num_iters);

		if (S > Sbase)
		{
			fMax = f;
			cMax = c;
			Sbase = S;
		}
	}

	fOpt.resize(X_org.rows());

	std::generate(fOpt.data(), fOpt.data() + fOpt.size(), [&] {
		return dist_0_1(eng);
		});

	for (size_t row = 0; row < fMax.size(); ++row) 
		fOpt(non_zero_rows[row]) = fMax(row);

	cOpt.resize(X_org.cols(), 0);

	for (size_t col = 0; col < cMax.size(); ++col) 
		cOpt(non_zero_cols[col]) = cMax(col);
}

void generate_f_from_c(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X_org, const refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt)
{
	std::default_random_engine eng;
	std::uniform_int_distribution<int> dist_0_1(0, 1);

	std::vector<uint64_t> non_zero_rows;

	auto X = remove_zero_rows(X_org,  non_zero_rows);

	if (X.cols() < 2 || X.rows() < 2)
	{
		fOpt.resize(X_org.rows(), 0);

		return;
	}

	Calculator_S calc_S(X);

	calc_S.set_c(cOpt);

/*	auto fMax = sign(calc_S.mult_Xtild_c());

	auto f1 = (fMax + 1.0) / 2.0;
	auto f2 = (1.0 - fMax) / 2.0;*/

	auto f1 = sign_01(calc_S.mult_Xtild_c());
	auto f2 = 1.0 - f1;

	calc_S.set_f(f1);
	double v1 = abs(calc_S.mult_fT_Xtild_c());
	calc_S.set_f(f2);
	double v2 = abs(calc_S.mult_fT_Xtild_c());

	refresh::matrix_1d<double> fMax;

	if (v1 >= v2)
		fMax = std::move(f1);
	else
		fMax = std::move(f2);

	fOpt.resize(X_org.rows());

	std::generate(fOpt.data(), fOpt.data() + fOpt.size(), [&] {
		return dist_0_1(eng);
		});

	for (size_t row = 0; row < fMax.size(); ++row)
		fOpt(non_zero_rows[row]) = fMax(row);
}


void get_most_freq_targets(
//	const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& anch_contingency_table,
	const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& anch_contingency_table,
	const std::vector<uint64_t> &targets,
	uint32_t n_most_freq_targets,
	std::vector<KmerAndCounter>& out) {

	auto row_sums = anch_contingency_table.get_row_sums();

	std::vector<std::pair<uint32_t, uint32_t>> data; //sum_counts, row_id
	for (size_t i = 0; i < row_sums.size(); ++i)
		data.emplace_back(row_sums(i), i);

	auto begin = data.begin();
	auto mid = n_most_freq_targets > data.size() ? data.end() : data.begin() + n_most_freq_targets;
	auto end = data.end();


	std::partial_sort(begin, mid, end, [](const auto& e1, const auto& e2) {
		return e1.first > e2.first;
	});

	out.clear();

	for (auto it = begin; it != mid; ++it) {
		KmerAndCounter item;
		item.kmer = targets[it->second];
		item.counter = it->first;
		out.emplace_back(item);
	}
}


//double calc_pval_base(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& sp_anch_contingency_table, size_t num_rand_cf, std::default_random_engine& eng) {
double calc_pval_base(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& sp_anch_contingency_table, size_t num_rand_cf, std::default_random_engine& eng) {
	std::uniform_int_distribution<int> dist_0_1(0, 1);
	double minus_1_and_1[] = { -1, 1 };
	double zero_and_1[] = { 0, 1 };
	double min_pval = std::numeric_limits<double>::max();

	refresh::matrix_1d<double> randCs_row(sp_anch_contingency_table.cols());
	refresh::matrix_1d<double> randFs_row(sp_anch_contingency_table.rows());
	for (size_t k = 0; k < num_rand_cf; ++k) {
		std::generate(randCs_row.data(), randCs_row.data() + randCs_row.size(), [&] { return minus_1_and_1[dist_0_1(eng)]; });
		std::generate(randFs_row.data(), randFs_row.data() + randFs_row.size(), [&] { return zero_and_1[dist_0_1(eng)];  });

		auto pval = testPval(sp_anch_contingency_table, randCs_row, randFs_row);
		if (pval < min_pval)
			min_pval = pval;
	}

	auto pval_base = min_pval * num_rand_cf;
	pval_base = pval_base < 1 ? pval_base : 1;

	return pval_base;
}

//double calc_pval_base(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& sp_anch_contingency_table, size_t num_rand_cf, std::default_random_engine& eng) {
// Extended version, returning also old_pval_base and fOptR and cOptR
double calc_pval_base_ext(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& sp_anch_contingency_table, size_t num_rand_cf, std::default_random_engine& eng,
	double &pval_base_old, refresh::matrix_1d<double> &fOptR, refresh::matrix_1d<double> &cOptR) {
	std::uniform_int_distribution<int> dist_0_1(0, 1);
	double minus_1_and_1[] = { -1, 1 };
	double zero_and_1[] = { 0, 1 };
	double min_pval = std::numeric_limits<double>::max();
	double min_pval_old = std::numeric_limits<double>::max();

	refresh::matrix_1d<double> randCs_row(sp_anch_contingency_table.cols());
	refresh::matrix_1d<double> randFs_row(sp_anch_contingency_table.rows());
	for (size_t k = 0; k < num_rand_cf; ++k) {
		std::generate(randCs_row.data(), randCs_row.data() + randCs_row.size(), [&] { return minus_1_and_1[dist_0_1(eng)]; });
		std::generate(randFs_row.data(), randFs_row.data() + randFs_row.size(), [&] { return zero_and_1[dist_0_1(eng)];  });

		auto pval = testPval(sp_anch_contingency_table, randCs_row, randFs_row);
		if (pval < min_pval)
			min_pval = pval;

		auto pval_old = testPvalOld(sp_anch_contingency_table, randCs_row, randFs_row);
		if (pval_old < min_pval_old)
		{
			min_pval_old = pval_old;
			fOptR = randFs_row;
			cOptR = randCs_row;
		}
	}

	auto pval_base = min_pval * num_rand_cf;
	pval_base = pval_base < 1 ? pval_base : 1;

	pval_base_old = min_pval_old * num_rand_cf;
	pval_base_old = pval_base_old < 1 ? pval_base_old : 1;

	return pval_base;
}

void compute_stats(
	Anchor&& anchor,
	size_t anchor_len_symbols,
	size_t n_uniq_targets,
	const std::unordered_set<uint64_t>& unique_samples,
	AnchorStats& anchor_stats,
	bool without_alt_max,
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
	double max_pval_opt_for_Cjs) {

	//std::cerr << "compute stats for anchor: " << kmer_to_string(anchor.anchor, anchor_len_symbols) << "\n";

	//targets are rows 
	//samples are cols

	size_t n_unique_samples = unique_samples.size();
	//lets just load everything, but consider building Xtrain and Xtest during construction
	//mkokot_TODO: anch_contingency_table is build in row major manner while the representation is col major -> validate performance
//	refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major> sp_anch_contingency_table(n_uniq_targets, n_unique_samples);
	refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> sp_anch_contingency_table(n_uniq_targets, n_unique_samples);

	std::default_random_engine eng;
	//std::default_random_engine eng(std::random_device{}()); 

	SampleToColMapper mapper(unique_samples);

	//assumes targets are sorted and for each target samples are sorted

	auto target = anchor.data[0].target;

	std::vector<uint64_t> targets;

	if (n_most_freq_targets) {
		targets.resize(n_uniq_targets); //targets[row_id] = target associated with row_id
		targets[0] = target;
	}

	sp_anch_contingency_table.reserve(anchor.data.size());

	size_t row_id = 0;
	for (const auto& e : anchor.data) {
		if (target != e.target) { //new target -> new row
			++row_id;
			target = e.target;

			if (n_most_freq_targets) {
				targets[row_id] = target;
			}
		}

		auto sample_id = pack_smaple_id_target(e.sample_id, e.barcode);

		if (unique_samples.find(sample_id) == unique_samples.end()) {
			std::cerr << "cannot find sample in unique set\n";
		}

		auto col_id = mapper.map(sample_id);

//		sp_anch_contingency_table(row_id, col_id) = e.count;
		sp_anch_contingency_table.insert_unsafe(row_id, col_id, e.count);
	}

	sp_anch_contingency_table.fix();

	anchor.data.clear();
	anchor.data.shrink_to_fit();

	if (n_most_freq_targets)
		get_most_freq_targets(sp_anch_contingency_table, targets, n_most_freq_targets, anchor_stats.most_freq_targets);

	auto Xtrain = get_train_mtx_2(sp_anch_contingency_table, opt_train_fraction, eng);

	if (compute_also_old_base_pvals)
	{
		refresh::matrix_1d<double> fOptR;
		refresh::matrix_1d<double> cOptR;
		
		//mkokot_TODO: I am calculating this after splitting because the same generator is used and I want to keep new results as similar to the pre-memory optimization as possible, but in general this could be computed before get_train_mtx_2
		anchor_stats.pval_base = calc_pval_base_ext(sp_anch_contingency_table, num_rand_cf, eng, anchor_stats.pval_base_old, fOptR, cOptR);

		anchor_stats.effect_size_bin_old = effectSize_bin(sp_anch_contingency_table, cOptR, fOptR);
	}
	else
		anchor_stats.pval_base = calc_pval_base(sp_anch_contingency_table, num_rand_cf, eng); //mkokot_TODO: I am calculating this after splitting because the same generator is used and I want to keep new results as similar to the pre-memory optimization as possible, but in general this could be computed before get_train_mtx_2

	sp_anch_contingency_table -= Xtrain;

	auto Xtest = std::move(sp_anch_contingency_table);

	refresh::matrix_1d<double> cOpt(Xtrain.cols());
	refresh::matrix_1d<double> fOpt(Xtrain.rows());

	// AltMax
	if (!without_alt_max)
	{
		generate_alt_max_cf(Xtrain, cOpt, fOpt, opt_num_inits, opt_num_iters);

		anchor_stats.pval_opt = testPval(Xtest, cOpt, fOpt);

		if (with_pval_asymp_opt)
			anchor_stats.pval_asymp_opt = computeAsympSPLASH(Xtest, cOpt, fOpt);

		anchor_stats.effect_size_bin = effectSize_bin(Xtest, cOpt, fOpt);

		if (with_effect_size_cts)
			anchor_stats.effect_size_cts = effectSize_cts(Xtest, cOpt, fOpt);

		if (num_splits > 1) {
			//I create new cOpt because in the cOpt I will have the optimal one, which may be needed to print it
			refresh::matrix_1d<double> new_cOpt(Xtrain.cols());

			for (size_t split_no = 1; split_no < num_splits; ++split_no)
			{
				//reconstruct contingency table
				sp_anch_contingency_table = std::move(Xtest += Xtrain);

				//make a new split
				Xtrain = get_train_mtx_2(sp_anch_contingency_table, opt_train_fraction, eng);
				sp_anch_contingency_table -= Xtrain;
				auto Xtest = std::move(sp_anch_contingency_table);

				generate_alt_max_cf(Xtrain, new_cOpt, fOpt, opt_num_inits, opt_num_iters);

				if (with_pval_asymp_opt) {
					auto new_pval_asymp_opt = computeAsympSPLASH(Xtest, new_cOpt, fOpt);
					if (new_pval_asymp_opt < anchor_stats.pval_asymp_opt)
						anchor_stats.pval_asymp_opt = new_pval_asymp_opt;
				}

				auto new_pval_opt = testPval(Xtest, new_cOpt, fOpt);
				if (new_pval_opt < anchor_stats.pval_opt)
				{
					anchor_stats.pval_opt = new_pval_opt;
					anchor_stats.effect_size_bin = effectSize_bin(Xtest, new_cOpt, fOpt);

					if (with_effect_size_cts)
						anchor_stats.effect_size_cts = effectSize_cts(Xtest, new_cOpt, fOpt);

					if (cj_writer)
						cOpt = new_cOpt;
				}
			}

			anchor_stats.pval_opt *= num_splits;
			if (with_pval_asymp_opt) {
				anchor_stats.pval_asymp_opt *= num_splits;
			}
		}
		if (cj_writer && anchor_stats.pval_opt <= max_pval_opt_for_Cjs)
		{
			for (size_t col_id = 0; col_id < cOpt.size(); ++col_id)
			{
				uint64_t packed_sample_id_barcode = mapper.decode(col_id);
				uint64_t sample_id, barcode;
				unpack_sample_id_target(packed_sample_id_barcode, sample_id, barcode);
				cj_writer.write(anchor.anchor, sample_id, barcode, cOpt(col_id));
			}
		}
	}
}


void dump_fXc(const refresh::matrix_1d<double>& f, const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& c)
{
	std::ofstream ofs("fXc_dump.py");

	ofs << "import numpy as np\n";
	ofs << "def prep_fXc():\n";
	ofs << "    fOpt = np.zeros(" << f.size() << ")\n";
	for (size_t i = 0; i < f.size(); ++i)
		ofs << "    fOpt[" << i << "] = " << f(i) << endl;

	ofs << endl << "    cOpt = np.zeros(" << c.size() << ")\n";
	for (size_t i = 0; i < c.size(); ++i)
		ofs << "    cOpt[" << i << "] = " << c(i) << endl;

	ofs << endl << "    X = np.zeros((" << X.rows() << ", " << X.cols() << "))\n";
	for (auto p = X.begin(); p != X.end(); ++p)
		ofs << "    X[" << p->first.row << "][" << p->first.col << "] = " << p->second << endl;

	ofs << "    return [fOpt, X, cOpt]\n";
}
