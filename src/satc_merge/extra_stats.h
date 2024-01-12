#pragma once

#include "pvals.h"

class CEditDistanceOneWord
{
	uint32_t len;
	uint64_t masks[4];

public:
	CEditDistanceOneWord() = default;

	void Prepare(const uint64_t seq1, uint32_t _len);
	uint32_t Calculate(const uint64_t seq2);
	uint32_t Calculate(const uint64_t seq1, const uint64_t seq2, uint32_t _len);
};

class CExtraStats
{
	AnchorStats stats;
	CEditDistanceOneWord edow;
	size_t max_target_id;
	bool need_monte_carlo;
	size_t target_len;

	std::vector<std::pair<uint64_t, uint64_t>> target_counter;
	std::vector<std::pair<uint64_t, uint64_t>> monte_carlo_trials;

	void condense_targets(const Anchor& anchor);

	void compute_entropy();														// verified
	void compute_homopolymers(const size_t hp_len_thr);							// verified (initially)
	void compute_distances_max_target();										// verified
	void compute_distances_all_pairs();											// verified (initially)
	void compute_edit_distance_monte_carlo();									// verified (initially)
	void compute_hamming_distance();											// verified

	void determine_monte_carlo_trials(size_t no_trials);

	size_t hamming(uint64_t x, uint64_t y);
	size_t edit_distance_dp(uint64_t x, uint64_t y);

public:
	CExtraStats() = default;

	void Compute(const Anchor& anchor, size_t target_len_symbols, size_t min_hp_len, size_t all2all_max_no_targets,
		size_t n_uniq_targets, const std::unordered_set<uint64_t>& unique_samples, AnchorStats& anchor_stats);
};


