#include "extra_stats.h"
#include "../../libs/refresh/deterministic_random.h"
#include <array>
#include <random>
#include <algorithm>
#include <utility>
#include <random>
#include <set>

// *********************************************************************************************
void CEditDistanceOneWord::Prepare(const uint64_t seq1, uint32_t _len)
{
	uint64_t x = seq1;
	len = _len;

	masks[0] = masks[1] = masks[2] = masks[3] = 0;

	for (uint32_t i = 0; i < len; ++i)
	{
		masks[x & 3ull] += 1ull << i;
		x >>= 2;
	}
}

// *********************************************************************************************
// Based on Hyyro's bit-par algorithm "A Bit-Vector Algorithm for Computing Levenshtein and Damerau Edit Distances"
// https://www.semanticscholar.org/paper/A-Bit-Vector-Algorithm-for-Computing-Levenshtein-Hyyr%C3%B6/813e26d8920d17c2afac6bf5a15c537b067a128a/figure/1
uint32_t CEditDistanceOneWord::Calculate(const uint64_t seq2)
{
	uint64_t VN = 0;
	uint64_t VP = ~0ull >> (64 - len);
	uint64_t test_bit = 1ull << (len - 1);

	uint64_t Rp, HPp, HNp, VPp, VNp;
	uint32_t ed = len;

	uint64_t x = seq2;

	for (uint32_t i = 0; i < len; ++i)
	{
		uint64_t symbol = x & 3ull;
		x >>= 2;

		Rp = (((masks[symbol] & VP) + VP) ^ VP) | masks[symbol] | VN;
		HPp = VN | ~(Rp | VP);
		HNp = VP & Rp;
//		VPp = (HNp << 1) | ~(Rp | (HPp << 1));
//		VNp = (HPp << 1) & Rp;
		VPp = (HNp << 1) | ~(Rp | ((HPp << 1) | 1));
		VNp = ((HPp << 1) | 1) & Rp;

		if (HPp & test_bit)
			ed++;
		else if (HNp & test_bit)
			ed--;

		VP = VPp;
		VN = VNp;
	}

	return ed;
}

// *********************************************************************************************
uint32_t CEditDistanceOneWord::Calculate(const uint64_t seq1, const uint64_t seq2, uint32_t _len)
{
	Prepare(seq1, _len);

	return Calculate(seq2);
}

// *********************************************************************************************
//
// *********************************************************************************************


// *********************************************************************************************
void CExtraStats::condense_targets(const Anchor& anchor)
{
	target_counter.clear();

	auto& data = anchor.data;
	size_t size = data.size();

	target_counter.emplace_back(data.front().target, data.front().count);

	for (size_t i = 1; i < size; ++i)
		if (target_counter.back().first == data[i].target)
			target_counter.back().second += data[i].count;
		else
			target_counter.emplace_back(data[i].target, data[i].count);
}

// *********************************************************************************************
void CExtraStats::determine_monte_carlo_trials(size_t no_trials)
{
	monte_carlo_trials.clear();

	std::vector<uint64_t> cum_counts;
	cum_counts.reserve(target_counter.size());

	uint64_t x = 0;

	for (auto& t : target_counter)
	{
		x += t.second;
		cum_counts.emplace_back(x);
	}

	auto tot_count = cum_counts.back();

	std::mt19937_64 mt;
	det_uniform_int_distribution<> dist1(0, tot_count - 1);
	det_uniform_int_distribution<> dist2(0, tot_count - 2);

//	std::multiset<std::pair<uint64_t, uint64_t>> unq_s;		// !!! TODO: zmienic na unordered_multiset, albo od razu wyniki wpisywac do monte_carlo_trials

	monte_carlo_trials.clear();

//	while (unq_s.size() < no_trials)
	while (monte_carlo_trials.size() < no_trials)
	{
		uint64_t sum1 = dist1(mt);
		uint64_t sum2 = dist2(mt);

		if (sum2 >= sum1)
			++sum2;

		uint64_t id1 = distance(cum_counts.begin(), std::upper_bound(cum_counts.begin(), cum_counts.end(), sum1));
		uint64_t id2 = distance(cum_counts.begin(), std::upper_bound(cum_counts.begin(), cum_counts.end(), sum2));

//		if(id1 != id2)
//			unq_s.emplace(id1, id2);
		monte_carlo_trials.emplace_back(id1, id2);
	}

//	monte_carlo_trials.assign(unq_s.begin(), unq_s.end());
}

// *********************************************************************************************
void CExtraStats::compute_hamming_distance()
{
	auto size = target_counter.size();

	if (size < 2)
		return;

	std::array<std::array<uint64_t, 4>, 32> pos_hist{};
	std::array<std::array<uint64_t, 4>, 32> pos_hist_weighted{};

	uint64_t tot_count = 0;
	uint64_t tot_hamming_distance = 0;

	for (size_t i = 0; i < size; ++i)
	{
		uint64_t x = target_counter[i].first;
		uint64_t count = target_counter[i].second;
		tot_count += count;

		for (size_t j = 0; j < target_len; ++j)
		{
			pos_hist[j][x & 3]++;
			pos_hist_weighted[j][x & 3] += count;
			x >>= 2;
		}
	}

	for (size_t j = 0; j < target_len; ++j)
		for (int k = 0; k < 4; ++k)
			tot_hamming_distance += pos_hist_weighted[j][k] * (tot_count - pos_hist_weighted[j][k]);

	stats.avg_hamming_distance_all_pairs = (double)tot_hamming_distance / (tot_count * (tot_count - 1));
}

// *********************************************************************************************
void CExtraStats::compute_entropy()
{
	size_t size = target_counter.size();
	double entropy = 0;

	if (size > 1)
	{
		uint64_t tot_count = 0;

		for (auto& x : target_counter)
			tot_count += x.second;

		double den = tot_count;

		for (auto& x : target_counter)
		{
			double p = x.second / den;
			entropy += -p * log2(p);
		}
	}

	stats.entropy = entropy;
}

// *********************************************************************************************
void CExtraStats::compute_homopolymers(const size_t hp_len_thr)
{
	uint64_t total = 0;
	size_t no_hp = 0;

	for (const auto& rec : target_counter)
	{
		uint64_t x = rec.first;

		size_t hp_len = 0;
		uint64_t p_symbol = 4;		

		total += rec.second;

		for (size_t i = 0; i < target_len; ++i)
		{
			uint64_t symbol = x & 3ull;
			if (symbol == p_symbol)
			{
				if (++hp_len == hp_len_thr)
				{
					no_hp += rec.second;
					break;
				}
			}
			else
			{
				p_symbol = symbol;
				hp_len = 1;
			}

			x >>= 2;
		}
	}

	stats.avg_no_homopolymer_targets = (double)no_hp / total;
}

// *********************************************************************************************
void CExtraStats::compute_distances_max_target()
{
	uint64_t tot_count = 0;
	uint64_t hamming_dist = 0;
	uint64_t edit_dist = 0;

	uint64_t max_target = target_counter[max_target_id].first;

	edow.Prepare(max_target, target_len);

	for (size_t i = 0; i < target_counter.size(); ++i)
	{
		uint64_t count = target_counter[i].second;
		tot_count += count;

		if (i != max_target_id)
		{
			size_t hd = hamming(max_target, target_counter[i].first);
			hamming_dist += hd * count;

			size_t ed = edow.Calculate(target_counter[i].first);
			edit_dist += ed * count;
		}
	}

	stats.avg_hamming_distance_max_target = (double)hamming_dist / tot_count;
	stats.avg_edit_distance_max_target = (double)edit_dist / tot_count;
}

// *********************************************************************************************
void CExtraStats::compute_distances_all_pairs()
{
	uint64_t tot_count = 0;
	uint64_t hamming_dist = 0;
	uint64_t edit_dist = 0;

	auto size = target_counter.size();

	for (size_t i = 0; i < size; ++i)
	{
		uint64_t target_i = target_counter[i].first;
		edow.Prepare(target_i, target_len);
		uint64_t target_i_count = target_counter[i].second;

		tot_count += target_counter[i].second * (target_counter[i].second - 1);

		for (size_t j = i + 1; j < size; ++j)
		{
			uint64_t count = target_counter[j].second;
			tot_count += 2 * count * target_i_count;

			size_t hd = hamming(target_i, target_counter[j].first);
			hamming_dist += 2 * hd * count * target_i_count;

			size_t ed = edow.Calculate(target_counter[j].first);
			edit_dist += 2 * ed * count * target_i_count;
		}
	}

	stats.avg_hamming_distance_all_pairs = (double)hamming_dist / tot_count;
	stats.avg_edit_distance_all_pairs = (double)edit_dist / tot_count;
}

// *********************************************************************************************
void CExtraStats::compute_edit_distance_monte_carlo()
{
//	uint64_t tot_count = 0;
	uint64_t edit_dist = 0;

	auto size = monte_carlo_trials.size();

	for (auto& trial : monte_carlo_trials)
		if(target_counter[trial.first].first != target_counter[trial.second].first)
			edit_dist += edow.Calculate(target_counter[trial.first].first, target_counter[trial.second].first, target_len);

	if (size)
		stats.avg_edit_distance_all_pairs = (double)edit_dist / size;
}

// *********************************************************************************************
size_t CExtraStats::hamming(uint64_t x, uint64_t y)
{
	size_t r = 0;

	for (size_t i = 0; i < target_len; ++i)
	{
		r += (x & 3ull) != (y & 3ull);
		x >>= 2;
		y >>= 2;
	}

	return r;
}

// *********************************************************************************************
size_t CExtraStats::edit_distance_dp(uint64_t x, uint64_t y)
{
	std::array<std::array<int, 33>, 33> dp;

	for (size_t i = 0; i <= target_len; ++i)
		dp[i][0] = dp[0][i] = i;

	for (size_t i = 1; i <= target_len; ++i)
	{
		uint64_t x_s = x & 3ull;
		x >>= 2;

		uint64_t y2 = y;

		for (size_t j = 1; j <= target_len; ++j)
		{
			uint64_t y_s = y2 & 3ull;
			y2 >>= 2;

			int vl = dp[i][j - 1] + 1;
			int vu = dp[i - 1][j] + 1;
			int vd = dp[i - 1][j - 1] + (x_s != y_s);

			dp[i][j] = std::min({ vl, vu, vd });
		}
	}

	return dp[target_len][target_len];
}

// *********************************************************************************************
void CExtraStats::Compute(const Anchor& _anchor, size_t anchor_len_symbols, size_t target_len_symbols, size_t min_hp_len, size_t all2all_max_no_targets,
	size_t n_uniq_targets, const std::unordered_set<uint64_t>& unique_samples, AnchorStats& anchor_stats)
{
	stats = anchor_stats;

	stats.clear_extra();

	if (_anchor.data.empty())
		return;

	condense_targets(_anchor);

	size_t size = target_counter.size();
	need_monte_carlo = size > all2all_max_no_targets;

	anchor_len = anchor_len_symbols;
	target_len = target_len_symbols;

	anchor = _anchor.anchor;

	max_target_id = std::distance(target_counter.begin(), std::max_element(target_counter.begin(), target_counter.end(), [](const auto& x, const auto& y) {return x.second < y.second; }));

	compute_entropy();
	compute_homopolymers(min_hp_len);
	compute_distances_max_target();

	if (need_monte_carlo)
	{
		compute_hamming_distance();
		determine_monte_carlo_trials(all2all_max_no_targets * (all2all_max_no_targets - 1) / 2);
		compute_edit_distance_monte_carlo();
	}
	else
		compute_distances_all_pairs();

	anchor_stats = stats;
}
