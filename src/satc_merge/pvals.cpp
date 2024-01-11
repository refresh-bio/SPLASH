#include <random>
#include <unordered_map>
#include <algorithm>
#include <numeric>

#include "pvals.h"
#include "get_train_mtx.h"
#include "ob_utils.h"

//mkokot_TODO: it would be better to just have a vector for mapping, but in the current packed representation for 10X it may not work
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

//void get_spectral_cf_svd(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt) {
void get_spectral_cf_svd(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt) {

	std::vector<uint64_t> non_zero_cols;
	std::vector<uint64_t> non_zero_rows;

	auto X = remove_zero_cols_rows(X_org, non_zero_cols, non_zero_rows);

	if (X.cols() < 2 || X.rows() < 2)
	{
		cOpt.resize(X_org.cols(), 0);
		fOpt.resize(X_org.rows(), 0);

		return;
	}

	refresh::matrix_1d<double> X_col_sum = X.get_col_sums();
	refresh::matrix_1d<double> X_row_sum = X.get_row_sums();

	double M_rcpt = 1.0 / X_row_sum.sum();

	refresh::matrix_dense<double, refresh::matrix_col_major> svdmat(X.rows(), X.cols());

	for (size_t i = 0; i < X.cols(); ++i)
	{
		double mult = X_col_sum(i) * M_rcpt / X_col_sum(i);
		for (size_t j = 0; j < X.rows(); ++j)
			svdmat(j, i) = -X_row_sum(j) * mult;
	}

	for (auto p = X.begin(); p != X.end(); ++p)
		svdmat(p->first.row, p->first.col) += p->second / X_col_sum(p->first.col);

	refresh::matrix_1d<double> fRaw;
	refresh::matrix_1d<double> cRaw;

	ob_svd(svdmat, fRaw, cRaw);

	refresh::matrix_1d<double> cGuess = std::move(cRaw);
	normalize_vec_0_1(fRaw);
	refresh::matrix_1d<double> fGuess = std::move(fRaw);

	assert(fGuess.size() == non_zero_rows.size());

	refresh::matrix_1d<double> fElong(X_org.rows(), 0.5);

	for (size_t row = 0; row < fGuess.size(); ++row) //fElong[relevantTargs] = fGuess
		fElong(non_zero_rows[row]) = fGuess(row);

	fOpt.resize(fElong.size());
	for (size_t i = 0; i < fElong.size(); ++i)
		fOpt(i) = fElong(i);

	refresh::matrix_1d<double> cElong(X_org.cols(), 0);

	assert(cGuess.size() == non_zero_cols.size());
	for (size_t col = 0; col < cGuess.size(); ++col) // cElong[np.arange(tblShape[1])[relevantSamples]] = cGuess ### fancy indexing
		cElong(non_zero_cols[col]) = cGuess(col);

	cOpt = move(cElong);
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

//double testPval(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {
double testPval(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, const refresh::matrix_1d<double>& cOpt, const refresh::matrix_1d<double>& fOpt) {

	if (cOpt.all_items_same() || fOpt.all_items_same())
		return 1.0;
	
	if (!(fOpt.max_coeff() <= 1 && fOpt.min_coeff() >= 0)) {
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

	X.print("X");

	c.print("c");

	for (int i = 0; i < n_iters; ++i)
	{
		f = sign(calc_S.mult_Xtild_c());

		f.print("f");
		auto f1 = (f + 1.0) / 2.0;
		auto f2 = (1.0 - f) / 2.0;
		f1.print("f1");
		f2.print("f2");

		calc_S.set_f(f1);
		double v1 = abs(calc_S.mult_fT_Xtild_c());
		calc_S.set_f(f2);
		double v2 = abs(calc_S.mult_fT_Xtild_c());

//		std::cout << "v1: " << v1 << "   v2: " << v2 << std::endl;

		if (v1 >= v2)
			f = f1;
		else
			f = f2;
		
		calc_S.set_f(f);

		f.print("f");

		c = calc_S.mult_fT_Xtild();

		c.print("c");

		auto norm_c = c.norm();
		if (norm_c > 0)
			c /= norm_c;

		c.print("c after norm");

		calc_S.set_c(c);

		S = calc_S.mult_fT_Xtild_c();

//		std::cout << S << std::endl;

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

//void generate_alt_max_cf(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt, int no_tries, int altMaximize_iters)
void generate_alt_max_cf(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X_org, refresh::matrix_1d<double>& cOpt, refresh::matrix_1d<double>& fOpt, int no_tries, int altMaximize_iters)
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

		double S = altMaximize(X, c, f, altMaximize_iters);

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
	auto mid = data.begin() + n_most_freq_targets;
	auto end = data.end();

	if (mid > end)
		mid = end;

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

void compute_stats(
	Anchor&& anchor,
	size_t anchor_len_symbols,
	size_t n_uniq_targets,
	const std::unordered_set<uint64_t>& unique_samples,
	AnchorStats& anchor_stats,
	bool without_SVD,
	bool with_effect_size_cts,
	uint32_t n_most_freq_targets,
	double train_fraction,
	int generate_alt_max_cf_no_tires,
	int altMaximize_iters,
	CjWriter& cj_writer,
	double max_pval_rand_init_alt_max_for_Cjs) {

	//std::cerr << "compute stats for anchor: " << kmer_to_string(anchor.anchor, anchor_len_symbols) << "\n";

	//mkokot_TODO: make this a parameters
	size_t num_rand_cf = 50;         //--num_rand_cf  in compute_spectral_pvals.py

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

	sp_anch_contingency_table.reserve((size_t) (n_uniq_targets * 1.2));

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

	sp_anch_contingency_table.shrink_to_fit();
	sp_anch_contingency_table.fix();

	anchor.data.clear();
	anchor.data.shrink_to_fit();

	if (n_most_freq_targets)
		get_most_freq_targets(sp_anch_contingency_table, targets, n_most_freq_targets, anchor_stats.most_freq_targets);

	auto Xtrain = get_train_mtx_2(sp_anch_contingency_table, train_fraction, eng);

	anchor_stats.pval_base = calc_pval_base(sp_anch_contingency_table, num_rand_cf, eng); //mkokot_TODO: I am calculating this after splitting because the same generator is used and I want to keep new results as similar to the pre-memory optimization as possible, but in general this could be computed before get_train_mtx_2

	sp_anch_contingency_table -= Xtrain;

	auto Xtest = std::move(sp_anch_contingency_table);

	refresh::matrix_1d<double> cOpt, fOpt;

	// SVD
	if (!without_SVD)
	{
		if (n_uniq_targets < 1000)
		{
			get_spectral_cf_svd(Xtrain, cOpt, fOpt);

			anchor_stats.pval_SVD_corrAnalysis = testPval(Xtest, cOpt, fOpt);
		}
		else
			anchor_stats.pval_SVD_corrAnalysis = 1.0;
	}
	// AltMax

	generate_alt_max_cf(Xtrain, cOpt, fOpt, generate_alt_max_cf_no_tires, altMaximize_iters);

	anchor_stats.pval_rand_init_alt_max = testPval(Xtest, cOpt, fOpt);

	if (with_effect_size_cts)
		anchor_stats.effect_size_cts = effectSize_cts(Xtest, cOpt, fOpt);

	anchor_stats.effect_size_bin = effectSize_bin(Xtest, cOpt, fOpt);

	if (cj_writer && anchor_stats.pval_rand_init_alt_max <= max_pval_rand_init_alt_max_for_Cjs)
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
