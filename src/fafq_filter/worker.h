#pragma once
#include <unordered_map>
#include <map>
#include <iostream>
#include <string>

#include "params.h"
#include <refresh/compression/lib/file_wrapper.h>
#include "../common/types/satc_data.h"
#include "../common/accepted_anchors.h"

class CWorker
{
	size_t stats_n_reads;
	size_t stats_n_bases;
	size_t stats_n_allowed_reads;
	size_t stats_n_allowed_bases;

	CParams params;
	array<uint8_t, 256> char2bits;

	shared_ptr<AcceptedAnchors> accepted_anchors;
	unordered_map<kmer_t, unordered_map<kmer_t, size_t>> anchor_stats;
//	map<kmer_t, unordered_map<kmer_t, size_t>> anchor_stats;
	kmer_t anchor_mask;
	kmer_t target_mask;
	const kmer_t empty_mask = (kmer_t) ~0ull;
	vector<kmer_t> rec_kmers;

	void init()
	{
		fill(char2bits.begin(), char2bits.end(), 4);
		char2bits['A'] = char2bits['a'] = 0;
		char2bits['C'] = char2bits['c'] = 1;
		char2bits['G'] = char2bits['g'] = 2;
		char2bits['T'] = char2bits['t'] = 3;
	}

	uint8_t no_bytes(size_t x)
	{
		uint8_t nb = 1;
		x /= 256;

		while (x)
		{
			++nb;
			x /= 256;
		}

		return nb;
	}

	void look_for_anchors(const string& str);
	void prepare_kmers(const string& str);
	kmer_t get_target(const string& str, size_t pos);

	bool load_fastx(const string& in_fn);

	bool store_satc(const uint32_t sample_id, const string& out_fn);

public:
	CWorker(const CParams& params, const shared_ptr<AcceptedAnchors> accepted_anchors) :
		params(params),
		accepted_anchors(accepted_anchors)
	{
		init();
		anchor_mask = ((kmer_t)~0ull) >> (64 - 2 * params.anchor_len);
		target_mask = ((kmer_t)~0ull) >> (64 - 2 * params.target_len);
	}

	bool process_file(const uint32_t sample_id, const string& in_fn, const string& out_fn);
	bool process_stdin(const uint32_t sample_id, const string& out_fn);
};