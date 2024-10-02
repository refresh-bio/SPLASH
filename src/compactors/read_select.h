#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <cinttypes>
#include <functional>
#include <utility>

#include <zlib-ng/zlib.h>
#include <refresh/hash_tables/lib/hash_map.h>
#include "../common/types/common_types.h"

using namespace std;

class ReadSelector
{
public:
	using dict_t = vector<pair<uint64_t, vector<vector<uint64_t>>>>;
	//	static const size_t max_n_followers = 16;
	enum class direction_t { forward, reverse};
	const uint64_t empty_kmer = ~0ull;

private:
	using map_t = refresh::hash_map_lp<uint64_t, uint64_t>;

	enum class process_mode_t {anchor_followers, extender_anchor};

	uint32_t anchor_len;		
	uint32_t follower_len;		
	uint32_t extender_len;		
	uint32_t gap_len;		
	uint32_t n_followers;
	dict_t dict;
	process_mode_t process_mode;
	bool find_all_mode = false;
	input_format_t input_format = input_format_t::fastq;

	dict_t dict_fwd;
	dict_t dict_rev;
	
	const uint64_t empty_anchor = ~0ull;
	const uint64_t excluded_val = ~0ull;
	map_t dict_map;

	const uint64_t no_bloom_funcs = 2;
	const double bloom_fill_factor = 0.01;
	vector<uint64_t> anchor_bloom;
	uint64_t anchor_bloom_size;
	uint64_t anchor_bloom_mask;
	uint64_t anchor_bloom_no_bits_in_main_part;

	vector<string> file_names;
	vector<pair<uint64_t, uint64_t>> file_size_id;
	string output_dir;

	size_t max_mem;
	size_t mem_in_use;
	size_t no_threads;
	size_t max_read_len;
	bool need_preprocessing;
	direction_t direction;
	bool keep_temps{ false };

	uint8_t base_to_code[256];
	char uint8_to_bases[256][3];
	uint8_t uint8_to_codes[256][3];
	uint64_t uint8_to_codes64[256][3];
	bool code_with_N[256];
	bool code_with_EOR[256];
	uint32_t no_EORs[256];

	vector<pair<uint8_t*, int64_t>> prefetched_files;

	void fill_base_coding_tables();

	string dna_file_name(const string& file_name);
	size_t dna_compress(char* dna_raw, uint8_t* dna_packed);

	bool open_input_dna_file(const string& fn_in, FILE*& dna_in);
	bool open_input_FASTQ_file(const string& fn_in, gzFile& gz_in);
	bool open_output_file(const string& fn_out, FILE*& out);

	bool preprocess_fastq(const string& file_name, size_t file_id);
	bool prefetch_dna_files();
	void release_prefetch_files();

	void determine_read_kmers(vector<uint64_t>& read_kmers, uint32_t k, uint8_t* ptr, const size_t len);
	void determine_clean_read_kmers(vector<uint64_t>& read_kmers, uint32_t k, uint8_t* ptr, const size_t len);
	void determine_read_kmer(uint64_t &read_kmer, uint32_t k, uint8_t* ptr, const size_t pos);
	void determine_clean_read_kmer(uint64_t& read_kmer, uint32_t k, uint8_t* ptr, const size_t pos);

	bool process_file(const string& fn, const size_t file_id);
	void process_read_anchor_followers(vector<uint64_t>& read_kmers, vector<uint64_t>& aux_kmers, uint8_t* ptr, size_t len, bool contains_Ns, size_t file_id);
	void process_read_anchor_extender(vector<uint64_t>& read_kmers, uint8_t* ptr, size_t len, bool contains_Ns, size_t file_id);

	uint64_t MurMur64Hash(uint64_t h) const noexcept
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;

		return h;
	}

	uint64_t MurMur64TruncHash(uint64_t h) const noexcept
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
/*		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;*/

		return h;
	}

	void init_bloom(const unordered_map<uint64_t, bool>& anchors);
	
	bool check_in_bloom(uint64_t x)
	{
//		auto h = MurMur64Hash(x);
		auto h = MurMur64TruncHash(x);

		uint64_t h_main = h & anchor_bloom_mask;
		h >>= anchor_bloom_no_bits_in_main_part;

		uint64_t query = 0;

		for (uint64_t i = 0; i < no_bloom_funcs; ++i)
		{
			query |= 1ull << (h & 63ull);
			h >>= 8;
		}

		return (anchor_bloom[h_main] & query) == query;
	}

	void process_general();

public:
	ReadSelector() :
		anchor_len(27),
		extender_len(27),
		follower_len(27),
		gap_len(0),
		n_followers(3),
		max_mem(8ull << 30),
		mem_in_use(0),
		no_threads(1),
		max_read_len(1 << 20),
		need_preprocessing(true),
		direction(direction_t::forward),
		process_mode(process_mode_t::anchor_followers)
	{
		fill_base_coding_tables();
	}

	~ReadSelector()
	{}

	void set_dict(const unordered_map<uint64_t, bool>& _dict);
	void set_input_names(const vector<string>& _file_names);
	void set_output_dir(const string& _output_dir);
	void set_max_memory_usage(size_t _max_mem);
	void set_no_threads(size_t _no_threads);
	void set_max_read_len(size_t _max_read_len);
	void set_input_format(const input_format_t _input_format);
	void set_keep_temps(bool v) { keep_temps = v; }
	bool get_keep_temps() const { return keep_temps; }


	dict_t& process_anchor_followers(direction_t _direction, uint32_t _anchor_len, uint32_t _follower_len, uint32_t _n_followers, bool _find_all_mode = false);
	dict_t& process_anchor_extender(direction_t _direction, uint32_t _anchor_len, uint32_t _extender_len, uint32_t _gap_len, bool _find_all_mode = true);

	pair<reference_wrapper<ReadSelector::dict_t>, reference_wrapper<ReadSelector::dict_t>> process_both_ways();
};

// EOF
