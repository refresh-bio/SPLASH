#pragma once
#include <cinttypes>
#include <vector>
#include <string>
#include <regex>

using namespace std;

using kmer_t = uint64_t;

struct CParams
{
	vector<string> input_fn;
	vector<string> output_fn;
	vector<uint32_t> sample_ids;
	string dict_fn;
	string sample_ids_fn;
	string stats_json_fn;
	regex allowed_headers{ ".*" };
	bool check_header = false;
	uint32_t n_top_targets = 1;
	uint32_t anchor_len = 0;
	uint32_t target_len = 0;
	uint32_t n_threads = 1;
	uint32_t verbosity = 0;
};

