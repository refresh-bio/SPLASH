#include "read_select.h"
#include <iostream>
#include <algorithm>
#include <thread>
#include <atomic>
#include <filesystem>
#include "../../src/common/kmer.h"

#include <cstring>

using namespace std::filesystem;

// ************************************************************************************
void ReadSelector::set_dict(const unordered_map<uint64_t, bool>& _dict)
{
	dict.clear();
	dict.shrink_to_fit();

	map_t empty_map(empty_anchor, 3*_dict.size(), 0.2);

	uint64_t no_active_anchors = 0;

	dict_map = empty_map;
	for (const auto x : _dict)
	{
		if (x.second)
		{
			dict.emplace_back(x.first, vector<vector<uint64_t>>());
			dict_map.insert_fast(make_pair(x.first, no_active_anchors++));
		}
		else
			dict_map.insert_fast(make_pair(x.first, excluded_val));
	}

	init_bloom(_dict);
}

// ************************************************************************************
void ReadSelector::set_input_format(const input_format_t _input_format)
{
	input_format = _input_format;
}

// ************************************************************************************
void ReadSelector::set_input_names(const vector<string>& _file_names)
{
	file_names = _file_names;

	for (size_t i = 0; i < file_names.size(); ++i)
	{
		std::filesystem::path fp(file_names[i]);
		file_size_id.emplace_back((uint64_t)file_size(fp), (uint64_t)i);
	}

	sort(file_size_id.begin(), file_size_id.end(), greater<pair<uint64_t, uint64_t>>());

	need_preprocessing = true;

	release_prefetch_files();

	mem_in_use = 0;
}

// ************************************************************************************
void ReadSelector::set_output_dir(const string& _output_dir)
{
	output_dir = _output_dir;
}

// ************************************************************************************
void ReadSelector::set_max_memory_usage(size_t _max_mem)
{
	max_mem = _max_mem;
}

// ************************************************************************************
void ReadSelector::set_no_threads(size_t _no_threads)
{
	no_threads = _no_threads;
}

// ************************************************************************************
void ReadSelector::set_max_read_len(size_t _max_read_len)
{
	max_read_len = _max_read_len;
}

// ************************************************************************************
void ReadSelector::init_bloom(const unordered_map<uint64_t, bool>& anchors)
{
	uint64_t bloom_size_in_bits = (uint64_t) ((2.0 * anchors.size()) / bloom_fill_factor);

	while (bloom_size_in_bits & (bloom_size_in_bits - 1))
		bloom_size_in_bits &= (bloom_size_in_bits - 1);

	bloom_size_in_bits <<= 1;
	bloom_size_in_bits = max(bloom_size_in_bits, (uint64_t)1024ull);

	anchor_bloom_size = bloom_size_in_bits / 64;
	anchor_bloom_mask = anchor_bloom_size - 1ull;

	anchor_bloom_no_bits_in_main_part = 0;
	for (auto x = anchor_bloom_size; x; x >>= 1)
		++anchor_bloom_no_bits_in_main_part;
	--anchor_bloom_no_bits_in_main_part;

	anchor_bloom.resize(anchor_bloom_size, 0ull);

	for (const auto x : anchors)
	{
//		auto h = MurMur64Hash(x.first);
		auto h = MurMur64TruncHash(x.first);

		uint64_t h_main = h & anchor_bloom_mask;
		h >>= anchor_bloom_no_bits_in_main_part;

		uint64_t query = 0;

		for (uint64_t i = 0; i < no_bloom_funcs; ++i)
		{
			query |= 1ull << (h & 63ull);
			h >>= 8;
		}

		anchor_bloom[h_main] |= query;
	}
}

// ************************************************************************************
void ReadSelector::fill_base_coding_tables()
{
	fill_n(base_to_code, 256, 5);
	base_to_code['A'] = 0;
	base_to_code['C'] = 1;
	base_to_code['G'] = 2;
	base_to_code['T'] = 3;
	base_to_code['N'] = 4;

	for (int i = 216; i < 256; ++i)
		uint8_to_bases[i][0] = uint8_to_bases[i][1] = uint8_to_bases[i][2] = 0;

	fill_n(code_with_N, 256, false);
	fill_n(code_with_EOR, 256, false);

	fill_n(no_EORs, 256, 0);

	for (int i = 0; i < 216; ++i)
	{
		uint8_to_bases[i][0] = "ACGTN\0"[i / 36];
		uint8_to_bases[i][1] = "ACGTN\0"[(i / 6) % 6];
		uint8_to_bases[i][2] = "ACGTN\0"[i % 6];
	}

	for (int i = 216; i < 256; ++i)
		uint8_to_codes[i][0] = uint8_to_codes[i][1] = uint8_to_codes[i][2] = 5;

	for (int i = 0; i < 216; ++i)
	{
		uint8_to_codes[i][0] = i / 36;
		uint8_to_codes[i][1] = (i / 6) % 6;
		uint8_to_codes[i][2] = i % 6;

		if (i / 36 == 4 || (i / 6) % 6 == 4 || i % 6 == 4)
			code_with_N[i] = true;

		if (i % 6 == 5)
			code_with_EOR[i] = true;

		if (i % 6 == 5)
			no_EORs[i]++;
		if((i / 6) % 6 == 5)
			no_EORs[i]++;
		if (i / 36 == 5)
			no_EORs[i]++;
	}

	for (int i = 0; i < 256; ++i)
	{
		uint8_to_codes64[i][0] = uint8_to_codes[i][0];
		uint8_to_codes64[i][1] = uint8_to_codes[i][1];
		uint8_to_codes64[i][2] = uint8_to_codes[i][2];
	}
}

// ************************************************************************************
string ReadSelector::dna_file_name(const string& file_name)
{
	if (output_dir.empty())
		return file_name + ".dna";
	else
		return output_dir + "/" + path(file_name).filename().string() + ".dna";
}

// ************************************************************************************
size_t ReadSelector::dna_compress(char* dna_raw, uint8_t* dna_packed)
{
	size_t len = strlen(dna_raw);
	
	while (len && (dna_raw[len - 1] == '\n' || dna_raw[len - 1] == '\r'))
		--len;

	dna_raw[len] = 0;

	// Add terminators for easier compression
	dna_raw[len+1] = 0;
	dna_raw[len+2] = 0;

	// Rounding up to multiplicity of 3 - just for easier compression
	len = 3 * ((len + 3) / 3);
	size_t packed_len = 0;

	for (size_t i = 0; i < len; i += 3, ++packed_len)
		dna_packed[packed_len] = base_to_code[dna_raw[i]] * 36 + base_to_code[dna_raw[i + 1]] * 6 + base_to_code[dna_raw[i + 2]];

//	strcpy((char*)dna_packed, dna_raw);

	return packed_len;
//	return strlen(dna_raw) + 1;
}

// ************************************************************************************
void ReadSelector::process_general()
{
	for (auto& ad : dict)
	{
		ad.second.clear();
		ad.second.resize(file_names.size());
	}

	dict_fwd.clear();
	dict_fwd.shrink_to_fit();
	dict_rev.clear();
	dict_fwd.shrink_to_fit();

	vector<thread> threads;
	atomic<int> file_id{ 0 };

	if (need_preprocessing)
	{
		for (size_t i = 0; i < no_threads; ++i)
			threads.emplace_back([&] {
			while (true)
			{
				int c_id = file_id.fetch_add(1);
				if (c_id < file_names.size())
				{
					int f_id = file_size_id[c_id].second;

					// if keep temps - try to load a temporary file
					bool temp_loaded = false;
					if (keep_temps) {
						temp_loaded = process_file(file_names[f_id], f_id);
					}

					if (!temp_loaded) {
						preprocess_fastq(file_names[f_id], f_id);
					}
				}
				else
					break;
			}
				});

		for (auto& th : threads)
			th.join();

		threads.clear();

		need_preprocessing = false;

		prefetch_dna_files();
	}
	else
	{
		for (size_t i = 0; i < no_threads; ++i)
			threads.emplace_back([&] {
			while (true)
			{
				int c_id = file_id.fetch_add(1);
				if (c_id < file_names.size())
				{
					int f_id = file_size_id[c_id].second;
					process_file(file_names[f_id], f_id);
				}
				else
					break;
			}
				});

		for (auto& th : threads)
			th.join();

		threads.clear();
	}
}

// ************************************************************************************
ReadSelector::dict_t& ReadSelector::process_anchor_followers(direction_t _direction, uint32_t _anchor_len, uint32_t _follower_len, uint32_t _n_followers, bool _find_all_mode)
{
	direction = _direction;
	anchor_len = _anchor_len;
	follower_len = _follower_len;
	n_followers = _n_followers;
	process_mode = process_mode_t::anchor_followers;

	find_all_mode = _find_all_mode;

	process_general();

	return dict;
}

// ************************************************************************************
ReadSelector::dict_t& ReadSelector::process_anchor_extender(direction_t _direction, uint32_t _anchor_len, uint32_t _extender_len, uint32_t _gap_len, bool _find_all_mode)
{
	direction = _direction;
	extender_len = _extender_len;
	anchor_len = _anchor_len;
	gap_len = _gap_len;
	process_mode = process_mode_t::extender_anchor;

	find_all_mode = _find_all_mode;

	process_general();

	return dict;
}

// ************************************************************************************
bool ReadSelector::prefetch_dna_files()
{
	prefetched_files.resize(file_names.size(), make_pair(nullptr, 0));

	file_size_id.clear();

	for (size_t i = 0; i < file_names.size(); ++i)
	{
		std::filesystem::path fp(dna_file_name(file_names[i]));
		prefetched_files[i].second = (int64_t) file_size(fp);

		file_size_id.emplace_back((uint64_t)prefetched_files[i].second, i);
	}

	mem_in_use = 0;

	while (true)
	{
		int64_t largest_size = 0;
		int largest_id = -1;

		for (size_t i = 0; i < prefetched_files.size(); ++i)
			if (prefetched_files[i].first == nullptr && prefetched_files[i].second > largest_size)
			{
				largest_id = i;
				largest_size = prefetched_files[i].second;
			}

		if (largest_id < 0)		// all files prefetched
			break;

		if (mem_in_use + largest_size > max_mem)
		{
			prefetched_files[largest_id].second = -prefetched_files[largest_id].second;		// To avoid prefetching this file
			continue;
		}

		FILE* f;
		if (!open_input_dna_file(dna_file_name(file_names[largest_id]), f))
			return false;

		prefetched_files[largest_id].first = new uint8_t[prefetched_files[largest_id].second];
		fread(prefetched_files[largest_id].first, 1, prefetched_files[largest_id].second, f);
		fclose(f);

		mem_in_use += prefetched_files[largest_id].second;
	}

	// Convert file sizes of non-prefetched files to correct values
	for (auto& x : prefetched_files)
		if (x.second < 0)
			x.second = -x.second;

	// Process non-prefetched files first
	for (size_t i = 0; i < prefetched_files.size(); ++i)
		if (prefetched_files[i].first == nullptr)
			file_size_id[i].first += 1ull << 62;

	sort(file_size_id.begin(), file_size_id.end(), greater<pair<uint64_t, uint64_t>>());

	return true;
}

// ************************************************************************************
void ReadSelector::release_prefetch_files()
{
	for (auto& x : prefetched_files)
		if (x.first != nullptr)
			delete[] x.first;

	prefetched_files.clear();
}

// ************************************************************************************
void ReadSelector::process_read_anchor_followers(vector<uint64_t>& read_kmers, vector<uint64_t>& aux_kmers, uint8_t* ptr, size_t len, bool contains_Ns, size_t file_id)
{
	if (contains_Ns)
		determine_read_kmers(read_kmers, anchor_len, ptr, len);
	else
		determine_clean_read_kmers(read_kmers, anchor_len, ptr, len);

	if (direction == direction_t::reverse)
		reverse(read_kmers.begin(), read_kmers.end());

	size_t n_kmers = read_kmers.size();

	if (n_kmers <= follower_len)
		return;

	size_t i;
	size_t max_i = n_kmers - follower_len;
	bool aux_constructed = false;

	if (find_all_mode)
	{
		for (i = 0; i < max_i; ++i)
			if (check_in_bloom(read_kmers[i]) && dict_map.check(read_kmers[i]))
			{
				auto p = dict_map.find(read_kmers[i]);

				if (p->second == excluded_val)
					continue;

				if (read_kmers[i + follower_len] == empty_kmer)
					continue;

				if (!aux_constructed)
				{
					if (anchor_len != follower_len)
					{
						if (contains_Ns)
							determine_read_kmers(aux_kmers, follower_len, ptr, len);
						else
							determine_clean_read_kmers(aux_kmers, follower_len, ptr, len);
						if (direction == direction_t::reverse)
							reverse(aux_kmers.begin(), aux_kmers.end());

						n_kmers = aux_kmers.size();
					}
					else
						aux_kmers = read_kmers;

					aux_constructed = true;
				}

				size_t j_end = i + anchor_len + follower_len * n_followers;

				auto& ext_vec = dict[p->second].second[file_id];

				size_t n_extended = 0;

				for (size_t j = i + anchor_len; j < n_kmers && j < j_end; j += follower_len)
					if (aux_kmers[j] != empty_kmer)
					{
						ext_vec.emplace_back(aux_kmers[j]);
						++n_extended;
					}
					else
						break;

				for (; n_extended < n_followers; ++n_extended)
					ext_vec.emplace_back(empty_kmer);
			}
	}
	else
	{
		for (i = 0; i < max_i; ++i)
		{
			if (check_in_bloom(read_kmers[i]) && dict_map.check(read_kmers[i]))
				break;

			if (++i >= max_i)
				return;

			if (check_in_bloom(read_kmers[i]) && dict_map.check(read_kmers[i]))
				break;
		}

		if (i >= max_i)
			return;

		auto p = dict_map.find(read_kmers[i]);

		if (p->second == excluded_val)
			return;

		if (read_kmers[i + follower_len] == empty_kmer)
			return;

		if (anchor_len != follower_len)
		{
			if (contains_Ns)
				determine_read_kmers(read_kmers, follower_len, ptr, len);
			else
				determine_clean_read_kmers(read_kmers, follower_len, ptr, len);
			if (direction == direction_t::reverse)
				reverse(read_kmers.begin(), read_kmers.end());

			n_kmers = read_kmers.size();
		}

		size_t j_end = i + anchor_len + follower_len * n_followers;

		auto& ext_vec = dict[p->second].second[file_id];

		size_t n_extended = 0;

		for (size_t j = i + anchor_len; j < n_kmers && j < j_end; j += follower_len)
			if (read_kmers[j] != empty_kmer)
			{
				ext_vec.emplace_back(read_kmers[j]);
				++n_extended;
			}
			else
				break;

		for (; n_extended < n_followers; ++n_extended)
			ext_vec.emplace_back(empty_kmer);
	}
}

// ************************************************************************************
void ReadSelector::process_read_anchor_extender(vector<uint64_t>& read_kmers, uint8_t* ptr, size_t len, bool contains_Ns, size_t file_id)
{
	if (contains_Ns)
		determine_read_kmers(read_kmers, extender_len, ptr, len);
	else
		determine_clean_read_kmers(read_kmers, extender_len, ptr, len);

	if (direction == direction_t::reverse)
		reverse(read_kmers.begin(), read_kmers.end());

	size_t n_kmers = read_kmers.size();

	size_t read_len = n_kmers + extender_len - 1;
	size_t prefix_len = gap_len + anchor_len;

	if (n_kmers <= prefix_len)
		return;

	read_kmers.erase(read_kmers.begin(), read_kmers.begin() + prefix_len);
	n_kmers -= prefix_len;

	size_t i;

	for (i = 0; i < n_kmers; ++i)
	{
		if (check_in_bloom(read_kmers[i]) && dict_map.check(read_kmers[i]))
		{
			auto p = dict_map.find(read_kmers[i]);

			if (p->second == excluded_val)
			{
				if (find_all_mode)
					continue;
				else
					return;
			}

			uint64_t anchor_kmer;

			if (contains_Ns)
				determine_read_kmer(anchor_kmer, anchor_len, ptr, i);
			else
				determine_clean_read_kmer(anchor_kmer, anchor_len, ptr, i);

			if (anchor_kmer == empty_anchor)
			{
				if (find_all_mode)
					continue;
				else
					return;
			}

			auto& ext_vec = dict[p->second].second[file_id];
			ext_vec.emplace_back(anchor_kmer);

			if (!find_all_mode)
				return;
		}
	}
}

// ************************************************************************************
void ReadSelector::determine_read_kmers(vector<uint64_t>& read_kmers, uint32_t k, uint8_t* ptr, const size_t len)
{
	read_kmers.clear();
	read_kmers.reserve(3 * len - k + 3);

	CKmer kmer(k, kmer_mode_t::direct);

	uint8_t* str = uint8_to_codes[*ptr++];

	int m_i = 0;

	for (int i = 0; i < 3 * (int)len; ++i)
	{
		if (m_i == 3)
		{
			str = uint8_to_codes[*ptr++];
			m_i = 0;
		}

		uint8_t c = str[m_i++];

		if (c < 4)
			kmer.insert_direct(c);
		else
			kmer.Reset();

		if (i + 1 >= k)
		{
			if (kmer.is_full())
				read_kmers.emplace_back(kmer.data_aligned_dir());
			else
				read_kmers.emplace_back(empty_kmer);
		}
	}
}

// ************************************************************************************
// Variant for read without Ns
void ReadSelector::determine_clean_read_kmers(vector<uint64_t>& read_kmers, uint32_t k, uint8_t* ptr, const size_t len)
{
	CKmer kmer(k, kmer_mode_t::direct);

	size_t len_in_bases = 3 * len - no_EORs[ptr[len-1]];

	uint64_t* str64;

	read_kmers.clear();
	read_kmers.reserve(len_in_bases - k + 1);

	int m_i;
	int i;

	uint64_t kmer_mask = ~0ull >> (64 - 2 * k);
	uint64_t loc_kmer = 0;

	int max_i = (((int)k - 1) / 3) * 3;

	for (i = 0; i < max_i; i += 3)
	{
		str64 = uint8_to_codes64[*ptr++];

		loc_kmer <<= 2;
		loc_kmer += str64[0];

		loc_kmer <<= 2;
		loc_kmer += str64[1];

		loc_kmer <<= 2;
		loc_kmer += str64[2];
	}

	str64 = uint8_to_codes64[*ptr++];
	m_i = 0;

	for (; i < (int)k - 1; ++i)
	{
		loc_kmer <<= 2;
		loc_kmer += str64[m_i++];
	}

	for (; i < (int)len_in_bases && m_i < 3; ++i)
	{
		loc_kmer <<= 2;
		loc_kmer += str64[m_i++];
		loc_kmer &= kmer_mask;

		read_kmers.emplace_back(loc_kmer);
	}

	for (; i < (int)len_in_bases; ++i)
	{
		str64 = uint8_to_codes64[*ptr++];

		loc_kmer <<= 2;
		loc_kmer += str64[0];
		loc_kmer &= kmer_mask;

		read_kmers.emplace_back(loc_kmer);
		
		if (++i >= (int)len_in_bases)
			break;

		loc_kmer <<= 2;
		loc_kmer += str64[1];
		loc_kmer &= kmer_mask;

		read_kmers.emplace_back(loc_kmer);

		if (++i >= (int)len_in_bases)
			break;

		loc_kmer <<= 2;
		loc_kmer += str64[2];
		loc_kmer &= kmer_mask;

		read_kmers.emplace_back(loc_kmer);
	}
}

// ************************************************************************************
void ReadSelector::determine_read_kmer(uint64_t& read_kmer, uint32_t k, uint8_t* ptr, const size_t pos)
{
	CKmer kmer(k, kmer_mode_t::direct);

	ptr += pos / 3;

	uint8_t* str = uint8_to_codes[*ptr++];

	int m_i = pos % 3;

	for (int i = 0; i < (int)k; ++i)
	{
		if (m_i == 3)
		{
			str = uint8_to_codes[*ptr++];
			m_i = 0;
		}

		uint8_t c = str[m_i++];

		if (c < 4)
			kmer.insert_direct(c);
		else
		{
			read_kmer = empty_kmer;
			return;
		}
	}

	read_kmer = kmer.data_aligned_dir();
}

// ************************************************************************************
void ReadSelector::determine_clean_read_kmer(uint64_t& read_kmer, uint32_t k, uint8_t* ptr, const size_t pos)
{
	CKmer kmer(k, kmer_mode_t::direct);

	ptr += pos / 3;

	uint8_t* str = uint8_to_codes[*ptr++];

	int m_i = pos % 3;

	for (int i = 0; i < (int)k; ++i)
	{
		if (m_i == 3)
		{
			str = uint8_to_codes[*ptr++];
			m_i = 0;
		}

		uint8_t c = str[m_i++];

		kmer.insert_direct(c);
	}

	read_kmer = kmer.data_aligned_dir();
}

// ************************************************************************************
bool ReadSelector::process_file(const string& fn, const size_t file_id)
{
	uint8_t *packed_dna = new uint8_t[max_read_len + 3];
	size_t read_pos = 0;

	vector<uint64_t> read_kmers;
	vector<uint64_t> aux_kmers;
	bool contains_Ns = false;

	if (!prefetched_files.empty() && prefetched_files[file_id].first != nullptr)
	{
		size_t pos;
		size_t max_pos = prefetched_files[file_id].second;
		auto file_data = prefetched_files[file_id].first;

		for (pos = 0; pos < max_pos; ++pos)
		{
			packed_dna[read_pos] = file_data[pos];
			contains_Ns |= code_with_N[packed_dna[read_pos]];
//			if (packed_dna[read_pos] % 6 == 5)		// end-of-read marker
			if (code_with_EOR[packed_dna[read_pos]])		// end-of-read marker
			{
				if(process_mode == process_mode_t::anchor_followers)
					process_read_anchor_followers(read_kmers, aux_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				else
					process_read_anchor_extender(read_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				read_pos = 0;
				contains_Ns = false;
			}
			else
				++read_pos;

			if (++pos >= max_pos)
				break;

			packed_dna[read_pos] = file_data[pos];
			contains_Ns |= code_with_N[packed_dna[read_pos]];
//			if (packed_dna[read_pos] % 6 == 5)		// end-of-read marker
			if (code_with_EOR[packed_dna[read_pos]])		// end-of-read marker
			{
				if (process_mode == process_mode_t::anchor_followers)
					process_read_anchor_followers(read_kmers, aux_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				else
					process_read_anchor_extender(read_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				read_pos = 0;
				contains_Ns = false;
			}
			else
				++read_pos;
		}
	}
	else
	{
		FILE* f;
		if (!open_input_dna_file(dna_file_name(fn), f))
		{
			//cerr << "Cannot open " << dna_file_name(fn) << " file\n";
			delete[] packed_dna;
			return false;
		}

		while (!feof(f))
		{
			int c = getc(f);
			if (c == EOF)
				break;

			packed_dna[read_pos] = (uint8_t) c;
			contains_Ns |= code_with_N[c];
//			if (packed_dna[read_pos] % 6 == 5)		// end-of-read marker
			if (code_with_EOR[packed_dna[read_pos]])		// end-of-read marker
			{
				if (process_mode == process_mode_t::anchor_followers)
					process_read_anchor_followers(read_kmers, aux_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				else
					process_read_anchor_extender(read_kmers, packed_dna, read_pos + 1, contains_Ns, file_id);
				read_pos = 0;
				contains_Ns = false;
			}
			else
				++read_pos;
		}

		fclose(f);
	}

	delete[] packed_dna;

	return true;
}

// ***************************************************************************
bool ReadSelector::open_input_FASTQ_file(const string& fn_in, gzFile& gz_in)
{
	gz_in = gzopen(fn_in.c_str(), "r");

	if (!gz_in)
	{
		cerr << "Cannot open file: " << fn_in << endl;
		return false;
	}

	gzbuffer(gz_in, 32 << 20);

	return true;
}

// ***************************************************************************
bool ReadSelector::open_input_dna_file(const string& fn_in, FILE* &dna_in)
{
	dna_in = fopen(fn_in.c_str(), "rb");

	if (!dna_in)
	{
		//cerr << "Cannot open file: " << fn_in << endl;
		return false;
	}

	setvbuf(dna_in, nullptr, _IOFBF, 16 << 20);

	return true;
}

// ***************************************************************************
bool ReadSelector::open_output_file(const string& fn_out, FILE*& out)
{
	out = fopen(fn_out.c_str(), "wb");

	if (!out)
	{
		cerr << "Cannot create file: " << fn_out << endl;
		return false;
	}

	setvbuf(out, nullptr, _IOFBF, 16 << 20);

	return true;
}

// ************************************************************************************
bool ReadSelector::preprocess_fastq(const string& file_name, size_t file_id)
{
	gzFile gz_in;
	FILE* dna_out;

	vector<uint64_t> read_kmers;
	vector<uint64_t> aux_kmers;

	if (!open_input_FASTQ_file(file_name, gz_in))
		return false;
	if (!open_output_file(dna_file_name(file_name), dna_out))
	{
		gzclose(gz_in);
		return false;
	}

	char* tmp = new char[max_read_len + 4];
	char* dna = new char[max_read_len + 4];
	uint8_t* packed = new uint8_t[max_read_len / 3 + 2];

	bool contains_Ns = false;
	size_t n_reads = 0;

	while (!gzeof(gz_in))
	{
		gzgets(gz_in, tmp, max_read_len);	// id
		gzgets(gz_in, dna, max_read_len);	// bases

		if (input_format == input_format_t::fastq)
		{
			gzgets(gz_in, tmp, max_read_len);	// +
			gzgets(gz_in, tmp, max_read_len);	// qualities
		}

//		if (strlen(dna) < k_len * n_followers)
		if (process_mode == process_mode_t::anchor_followers)
		{
			if (strlen(dna) < anchor_len + follower_len * n_followers)
				continue;
		}
		else
		{
			if (strlen(dna) < anchor_len + extender_len + gap_len)
				continue;
		}

		auto packed_len = dna_compress(dna, packed);
		fwrite(packed, 1, packed_len, dna_out);

		for (size_t i = 0; i < packed_len; ++i)
			contains_Ns |= code_with_N[packed[i]];

		if(process_mode == process_mode_t::anchor_followers)
			process_read_anchor_followers(read_kmers, aux_kmers, packed, packed_len, contains_Ns, file_id);
		else
			process_read_anchor_extender(read_kmers, packed, packed_len, contains_Ns, file_id);
		if (!contains_Ns) {
			++n_reads;
		}
		contains_Ns = false;
	}

	delete[] tmp;
	delete[] dna;
	delete[] packed;

	gzclose(gz_in);
	fclose(dna_out);
	
	return true;
}

// EOF
