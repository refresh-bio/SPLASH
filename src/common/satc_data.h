#pragma once
#include <fstream>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <unordered_map>

#include "../../libs/refresh/zstd_file.h"

#define USE_ZSTD_FOR_TEMPS

template<typename T>
void append_int_msb(std::vector<uint8_t>& v, T x, int n_bytes)
{
	switch (n_bytes)
	{
	case 8: v.emplace_back((x >> 56) & 0xff);
	case 7: v.emplace_back((x >> 48) & 0xff);
	case 6: v.emplace_back((x >> 40) & 0xff);
	case 5: v.emplace_back((x >> 32) & 0xff);
	case 4: v.emplace_back((x >> 24) & 0xff);
	case 3: v.emplace_back((x >> 16) & 0xff);
	case 2: v.emplace_back((x >> 8) & 0xff);
	case 1: v.emplace_back(x & 0xff);
	}
}

//for binary streaming writing
class buffered_binary_writer {
#ifdef USE_ZSTD_FOR_TEMPS
	refresh::zstd_file out{ 9 };
#else
	std::ofstream out;
#endif
	std::vector<uint8_t> buff;

	struct {
		std::vector<uint8_t> prev_rec;
		std::vector<uint8_t> diff_rec;
		void compute_diffed(const std::vector<uint8_t>& rec) {
			diff_rec.clear();
			size_t i = 0;
			for (; i < prev_rec.size() && prev_rec[i] == rec[i]; ++i)
				;
			diff_rec.push_back(i);
			for (; i < rec.size(); ++i)
				diff_rec.push_back(rec[i]);
		}
	} compr_serialization;

	void flush() {
		out.write(reinterpret_cast<char*>(buff.data()), buff.size());
		buff.clear();
	}

	void assure_space(size_t size) {
		if (buff.size() + size > buff.capacity()) {
			flush();
			if (size > buff.capacity())
				buff.reserve(size);
		}
	}

	void write(const std::vector<uint8_t>& vec) {
		assure_space(vec.size());
		for (auto x : vec)
			buff.push_back(x);
	}

public:
	buffered_binary_writer(buffered_binary_writer&&) = default;
	buffered_binary_writer& operator=(buffered_binary_writer&&) = default;
	buffered_binary_writer(const std::string& path, size_t buff_size = 1ull << 25) {
#ifdef USE_ZSTD_FOR_TEMPS
		out.open_writing(path);
#else
		out.rdbuf()->pubsetbuf(0, 0);
		out.open(path, std::ios::binary);
#endif
		buff.reserve(buff_size);
	}
	operator bool() const {
#ifdef USE_ZSTD_FOR_TEMPS
		return out.is_opened_for_writing();
#else
		return out.operator bool();
#endif
	}
	void write(const uint8_t* ptr, size_t size) {
		assure_space(size);
		buff.insert(buff.end(), ptr, ptr + size);
	}

	template<typename T>
	void write_little_endian(const T& data)
	{
		assure_space(sizeof(data));

		for (uint32_t b = 0; b < sizeof(data); ++b)
			buff.push_back(data >> (8 * b));
	}

	void write_rec_compr(std::vector<uint8_t>& rec) {
		compr_serialization.compute_diffed(rec);
		write(compr_serialization.diff_rec);
		std::swap(compr_serialization.prev_rec, rec);
	}

	~buffered_binary_writer() {
		if (buff.size())
			flush();
	}
	void close() {
		if (buff.size())
			flush();
		out.close();
	}
};

//for binary streaming reading
class buffered_binary_reader {
#ifdef USE_ZSTD_FOR_TEMPS
	refresh::zstd_file in;
#else
	std::ifstream in;
#endif
	std::vector<uint8_t> buff;
	uint8_t* ptr{};
	size_t in_buff{};

	std::vector<uint8_t> prev_rec;

	void load()
	{
		//copy tail
		auto tail_size = buff.data() + in_buff - ptr;
		std::copy_n(buff.data() + buff.size() - tail_size, tail_size, buff.data());
		//read
#ifdef USE_ZSTD_FOR_TEMPS
		auto readed = in.read(reinterpret_cast<char*>(buff.data() + tail_size), buff.size() - tail_size);
		in_buff = tail_size + readed;
#else
		in.read(reinterpret_cast<char*>(buff.data() + tail_size), buff.size() - tail_size);
		in_buff = tail_size + in.gcount();
#endif
		ptr = buff.data();
	}
	bool assure_data_in_buffer(size_t size) {
		assert(size <= buff.size());
		if (ptr + size > buff.data() + in_buff) {
			load();
			if (ptr + size > buff.data() + in_buff)
				return false;
		}
		return true;
	}

	bool read(std::vector<uint8_t>& vec, size_t to_read) {
		if (!assure_data_in_buffer(to_read))
			return false;

		for (uint32_t b = 0; b < to_read; ++b)
			vec.push_back(*ptr++);

		return true;
	}

	size_t calc_buff_size(size_t fsize, size_t max_size) {
		return fsize < max_size ? fsize : max_size;
	}

	size_t get_file_size(const std::string& path) {
		FILE* f = fopen(path.c_str(), "rb");
		if (!f) {
			std::cerr << "Error: cannot open file " << path << " to get file size\n";
			exit(1);
		}
		fseek(f, 0, SEEK_END);
		size_t fsize = ftell(f);
		fclose(f);

		return fsize;
	}
	buffered_binary_reader(const std::string& path, size_t file_size, size_t max_buff_size) :
#ifdef USE_ZSTD_FOR_TEMPS
	in(9, calc_buff_size(file_size, 1 << 20), 1 << 13),
#endif
	 buff(calc_buff_size(file_size, max_buff_size)) {
#ifdef USE_ZSTD_FOR_TEMPS
		in.open_reading(path);
#else
		in.rdbuf()->pubsetbuf(0, 0);
		in.open(path, std::ios::binary);
#endif
		ptr = buff.data();
	}
public:
	buffered_binary_reader(const std::string& path, size_t buff_size = 1ull << 13) :
		buffered_binary_reader(path, get_file_size(path), buff_size) { //delegate ctor, to get file size only once

	}

	operator bool() const {
#ifdef USE_ZSTD_FOR_TEMPS
		return in.is_opened_for_reading();
#else
		return in.operator bool();
#endif
	}

	bool read(uint8_t*& _ptr, size_t size) {
		if (!assure_data_in_buffer(size))
			return false;

		_ptr = ptr;
		ptr += size;
		return true;
	}

	template<typename T>
	bool read_little_endian(T& data)
	{
		if (!assure_data_in_buffer(sizeof(data)))
			return false;

		data = T{};

		for (uint32_t b = 0; b < sizeof(data); ++b)
			data += (uint64_t)*ptr++ << (8 * b);

		return true;
	}

	bool read_rec_compr(std::vector<uint8_t>& vec, size_t rec_len) {
		vec.clear();
		uint8_t n_common{};
		if (!read_little_endian(n_common)) //endianness does not matter for 1 byte
			return false;

		for (size_t i = 0; i < n_common; ++i)
			vec.push_back(prev_rec[i]);

		if (!read(vec, rec_len - n_common))
		{
			std::cerr << "Error: only part of the record was stored in the file\n";
			exit(1);
		}

		prev_rec = vec;
		return true;
	}

	void close() {
		in.close();
	}
};

template<typename T>
bool LoadBigEndian(std::vector<uint8_t>::iterator& p, T& data, uint8_t data_size_in_bytes)
{
	data = T{};

	for (int i = 0; i < data_size_in_bytes; ++i, ++p)
	{
		data <<= 8;
		data += (T)*p;
	}

	return true;
}

inline uint64_t str_kmer_to_uint64_t(const std::string& kmer) {
	uint64_t res{};
	static auto mapping = []() {
		static uint8_t map[256];
		std::fill_n(map, 256, 0);
		map['A'] = map['a'] = 0;
		map['C'] = map['c'] = 1;
		map['G'] = map['g'] = 2;
		map['T'] = map['t'] = 3;
		return map;
	}();
	for (auto& s : kmer) {
		res <<= 2;
		res += mapping[static_cast<uint8_t>(s)];
	}
	return res;
}

inline std::string kmer_to_string(uint64_t kmer, uint8_t len) {
	std::string str_kmer;
	for (int i = len - 1; i >= 0; --i) {
		auto s = (kmer >> (2 * i)) & 3;
		str_kmer.push_back("ACGT"[s]);
	}
	return str_kmer;
}

struct Header {
	uint8_t sample_id_size_bytes;
	uint8_t barcode_size_bytes;
	uint8_t anchor_size_bytes;
	uint8_t target_size_bytes;
	uint8_t counter_size_bytes;
	uint8_t barcode_len_symbols;
	uint8_t anchor_len_symbols;
	uint8_t target_len_symbols;
	uint8_t gap_len_symbols;

	uint32_t rec_len;

	void serialize(buffered_binary_writer& out) {
		out.write_little_endian(sample_id_size_bytes);
		out.write_little_endian(barcode_size_bytes);
		out.write_little_endian(anchor_size_bytes);
		out.write_little_endian(target_size_bytes);
		out.write_little_endian(counter_size_bytes);
		out.write_little_endian(barcode_len_symbols);
		out.write_little_endian(anchor_len_symbols);
		out.write_little_endian(target_len_symbols);
		out.write_little_endian(gap_len_symbols);
	}

	void load(buffered_binary_reader& in) {
		in.read_little_endian(sample_id_size_bytes);
		in.read_little_endian(barcode_size_bytes);
		in.read_little_endian(anchor_size_bytes);
		in.read_little_endian(target_size_bytes);
		in.read_little_endian(counter_size_bytes);
		in.read_little_endian(barcode_len_symbols);
		in.read_little_endian(anchor_len_symbols);
		in.read_little_endian(target_len_symbols);
		in.read_little_endian(gap_len_symbols);

		rec_len = sample_id_size_bytes + barcode_size_bytes + anchor_size_bytes + target_size_bytes + counter_size_bytes;
	}

	void print(std::ostream& oss) {
		oss << "Header: \n";
		oss << "\tsample_id_size_bytes     : " << (uint64_t)sample_id_size_bytes << "\n";
		oss << "\tbarcode_size_bytes       : " << (uint64_t)barcode_size_bytes << "\n";
		oss << "\tanchor_size_bytes        : " << (uint64_t)anchor_size_bytes << "\n";
		oss << "\ttarget_size_bytes        : " << (uint64_t)target_size_bytes << "\n";
		oss << "\tcounter_size_bytes       : " << (uint64_t)counter_size_bytes << "\n";
		oss << "\tbarcode_len_symbols      : " << (uint64_t)barcode_len_symbols << "\n";
		oss << "\tanchor_len_symbols       : " << (uint64_t)anchor_len_symbols << "\n";
		oss << "\ttarget_len_symbols       : " << (uint64_t)target_len_symbols << "\n";
		oss << "\tgap_len_symbols          : " << (uint64_t)gap_len_symbols << "\n";
	}
};

enum class RecFmt { SATC, NOMAD };
struct RecFmtConv {

	inline static RecFmt from_string(const std::string& str) {
		if (str == "satc")
			return RecFmt::SATC;
		else if (str == "nomad")
			return RecFmt::NOMAD;
		else {
			std::cerr << "Error: cannot convert \"" << str << "\" to RecordPrintFormat\n";
			exit(1);
		}
	}
	inline static std::string to_string(RecFmt type) {
		switch (type)
		{
		case RecFmt::SATC:
			return "satc";
		case RecFmt::NOMAD:
			return "nomad";
		default:
			std::cerr << "Error: unsupported, swich in RecordPrintFormat should be extended";
			exit(1);
		}
	}
};

class SampleNameDecoder {
	std::vector<std::string> decoded;
	bool can_decode = false;
public:
	SampleNameDecoder(const std::string& path) {
		if(path == "") {
			return;
		}
		std::ifstream in(path);
		if(!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		can_decode = true;
		std::string name;
		size_t id;
		while(in >> name >> id)
		{
			if(id >= decoded.size()) {
				decoded.resize(id + 1);
			}
			decoded[id] = name;
		}
	}

	void store_sample_id(std::ostream& oss, size_t id) const {
		if(can_decode) {
			assert(id < decoded.size());
			oss << decoded[id];
		} else {
			oss << id;
		}
	}
};

class SampleNameToId {
	std::unordered_map<std::string, uint32_t> m;
public:
	SampleNameToId(const std::string& path) {
		std::ifstream in(path);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		std::string name;
		size_t id;
		while (in >> name >> id)
			m[name] = id;
	}

	size_t get_n_samples() const {
		return m.size();
	}

	bool get_sample_id(const std::string& sample_name, uint32_t& sample_id) const {
		auto it = m.find(sample_name);
		if (it == m.end())
			return false;
		sample_id = it->second;
		return true;
	}
};

inline uint64_t pack_smaple_id_target(uint64_t sample_id, uint64_t barcode) {
	return (sample_id << 32) + barcode;
}

inline void unpack_sample_id_target(uint64_t packed, uint64_t& sample_id, uint64_t& barcode) {
	sample_id = packed >> 32;
	barcode = packed & ((1ull << 32) - 1);
}
inline std::pair<uint64_t, uint64_t> unpack_sample_id_target(uint64_t packed) {
	std::pair<uint64_t, uint64_t> res;
	unpack_sample_id_target(packed, res.first, res.second);
	return res;
}


struct Record
{
	uint64_t sample_id;
	uint64_t barcode;
	uint64_t anchor;
	uint64_t target;
	uint64_t count;

	std::vector<uint8_t> compr_record;

private:
	void serialize_msb(std::vector<uint8_t>& out, const Header& header) {
		append_int_msb(out, sample_id, header.sample_id_size_bytes);
		append_int_msb(out, barcode, header.barcode_size_bytes);
		append_int_msb(out, anchor, header.anchor_size_bytes);
		append_int_msb(out, target, header.target_size_bytes);
		append_int_msb(out, count, header.counter_size_bytes);
	}

	void load_msb(std::vector<uint8_t>& in, const Header& header) {
		auto p = in.begin();
		LoadBigEndian(p, sample_id, header.sample_id_size_bytes);
		LoadBigEndian(p, barcode, header.barcode_size_bytes);
		LoadBigEndian(p, anchor, header.anchor_size_bytes);
		LoadBigEndian(p, target, header.target_size_bytes);
		LoadBigEndian(p, count, header.counter_size_bytes);
	}
public:
	void serialize(buffered_binary_writer& out, const Header& header) {
		compr_record.clear();
		serialize_msb(compr_record, header);
		out.write_rec_compr(compr_record);
	};

	bool load(buffered_binary_reader& in, const Header& header) {
		if (!in.read_rec_compr(compr_record, header.rec_len))
			return false;

		load_msb(compr_record, header);

		return true;
	}

	void print(std::ostream& oss, const Header& header, RecFmt format = RecFmt::SATC, const SampleNameDecoder& sample_name_decoder = SampleNameDecoder("")) {
		if (format == RecFmt::SATC) {
			sample_name_decoder.store_sample_id(oss, sample_id);
			oss << "\t";
			if (header.barcode_size_bytes)
				oss << kmer_to_string(barcode, header.barcode_len_symbols) << "\t";
			oss << kmer_to_string(anchor, header.anchor_len_symbols) << "\t";
			oss << kmer_to_string(target, header.target_len_symbols) << "\t";
			oss << count << "\n";
		}
		else if (format == RecFmt::NOMAD) {
			oss << count << " ";
			oss << kmer_to_string(anchor, header.anchor_len_symbols);
			oss << kmer_to_string(target, header.target_len_symbols) << " ";
			sample_name_decoder.store_sample_id(oss, sample_id);
			if (header.barcode_size_bytes)
				oss << "_" << kmer_to_string(barcode, header.barcode_len_symbols);

			oss << "\n";
		}
		else {
			std::cerr << "Error: unsupported print format\n";
			exit(1);
		}
	}
};
