#pragma once
#include <fstream>
#include <cinttypes>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <unordered_map>

#include <refresh/compression/lib/zstd_wrapper.h>
#include <refresh/conversions/lib/conversions.h>

#ifdef _WIN32
#define _bswap64(x) _byteswap_uint64(x)
#else
#define _bswap64(x) __builtin_bswap64(x)
#endif

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

template<typename T>
void load_int_msb(std::vector<uint8_t>::iterator &p, T &x, int n_bytes)
{
	x = (T)0;

	switch (n_bytes)
	{
	case 8: x += ((T) (*p++)) << 56;
	case 7: x += ((T)(*p++)) << 48;
	case 6: x += ((T)(*p++)) << 40;
	case 5: x += ((T)(*p++)) << 32;
	case 4: x += ((T)(*p++)) << 24;
	case 3: x += ((T)(*p++)) << 16;
	case 2: x += ((T)(*p++)) << 8;
	case 1: x += ((T)(*p++));
	}
}

//for binary streaming writing
class buffered_binary_writer {
	refresh::zstd_file out{ 3 };
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
	buffered_binary_writer& operator=(buffered_binary_writer&& rhs) {
		if (&rhs == this)
			return *this;
		close();

		out = std::move(rhs.out);
		buff = std::move(rhs.buff);
		compr_serialization = std::move(rhs.compr_serialization);

		return *this;
	}
	buffered_binary_writer(const std::string& path, size_t buff_size = 1ull << 25) {
		out.open_writing(path);
		buff.reserve(buff_size);
	}
	operator bool() const {
		return out.is_opened_for_writing();
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
	refresh::zstd_file in;
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

		auto readed = in.read(reinterpret_cast<char*>(buff.data() + tail_size), buff.size() - tail_size);
		in_buff = tail_size + readed;
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

		vec.insert(vec.end(), ptr, ptr + to_read);
		ptr += to_read;

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
	buffered_binary_reader(const std::string& path, size_t file_size, size_t max_io_buff_size, size_t max_buff_size) :
		in(3, calc_buff_size(file_size, max_io_buff_size), 1 << 13),
		buff(calc_buff_size(file_size, max_buff_size)) {
		in.open_reading(path);
		ptr = buff.data();
	}
public:
	buffered_binary_reader(const std::string& path, size_t max_io_buff_size = 1ull << 20, size_t buff_size = 1ull << 13) :
//	buffered_binary_reader(const std::string& path, size_t max_io_buff_size, size_t buff_size) :
		buffered_binary_reader(path, get_file_size(path), max_io_buff_size, buff_size) { //delegate ctor, to get file size only once

	}

	operator bool() const {
		return in.is_opened_for_reading();
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

		vec.insert(vec.end(), prev_rec.begin(), prev_rec.begin() + n_common);

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

inline uint32_t get_rev_compl_shift(uint32_t len)
{
	return (len + 31) / 32 * 32 - len;
}

//instead of x >> 2*p there is (x>>p) >> p, because
//for p=32 we would have x >> 64 which is UB
inline uint64_t shr_2p(uint64_t x, uint64_t p) {
	return (x >> p) >> p;
}

inline uint64_t get_rev_compl(const uint64_t kmer, const uint32_t shift) {
	uint64_t res;

	const uint64_t mask4 = 0b0000111100001111000011110000111100001111000011110000111100001111ull;
	const uint64_t mask2 = 0b0011001100110011001100110011001100110011001100110011001100110011ull;

	uint64_t sf4 = ((kmer & mask4) << 4) + ((kmer >> 4) & mask4);
	uint64_t sf2 = ((sf4 & mask2) << 2) + ((sf4 >> 2) & mask2);

	res = shr_2p((~_bswap64(sf2)), shift);

	return res;
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

inline void kmer_to_string(uint64_t kmer, uint8_t len, std::string& str_kmer) {
	int out_i{};
	for (int i = len - 1; i >= 0; --i) {
		auto s = (kmer >> (2 * i)) & 3;
		str_kmer[out_i++] = "ACGT"[s];
	}
}

struct Header {
	uint8_t sample_id_size_bytes{};
	uint8_t barcode_size_bytes{};
	uint8_t anchor_size_bytes{};
	uint8_t target_size_bytes{};
	uint8_t counter_size_bytes{};
	uint8_t barcode_len_symbols{};
	uint8_t anchor_len_symbols{};
	uint8_t target_len_symbols{};
	uint8_t gap_len_symbols{};

	uint32_t rec_len{};

	enum class ordering_t {SBATC = 0, ATSBC = 1, TASBC = 2};

	inline static const ordering_t default_ordering = ordering_t::SBATC;

	static std::string to_string(ordering_t ordering)
	{
		switch (ordering)
		{
		case ordering_t::ATSBC:
			return "ATSBC";
		case ordering_t::SBATC:
			return "SBATC";
		case ordering_t::TASBC:
			return "TASBC";
		default:
			assert(false);
		}
	}
	static ordering_t ordering_from_string(const std::string& ordering_str)
	{
		if (ordering_str == "ATSBC")
			return ordering_t::ATSBC;
		if (ordering_str == "SBATC")
			return ordering_t::SBATC;
		if (ordering_str == "TASBC")
			return ordering_t::TASBC;

		std::cerr << "Error: unknown ordering: " << ordering_str << "\n";
		exit(1);
	}

	ordering_t ordering = ordering_t::SBATC;

	bool operator==(const Header& rhs) const {
		return sample_id_size_bytes == rhs.sample_id_size_bytes &&
			barcode_size_bytes == rhs.barcode_size_bytes &&
			anchor_size_bytes == rhs.anchor_size_bytes &&
			target_size_bytes == rhs.target_size_bytes &&
			counter_size_bytes == rhs.counter_size_bytes &&
			barcode_len_symbols == rhs.barcode_len_symbols &&
			anchor_len_symbols == rhs.anchor_len_symbols &&
			target_len_symbols == rhs.target_len_symbols &&
			gap_len_symbols == rhs.gap_len_symbols &&
			rec_len == rhs.rec_len &&
			ordering == rhs.ordering;
	}
	bool operator!=(const Header& rhs) const
	{
		return !operator==(rhs);
	}
	void serialize(buffered_binary_writer& out) const {
		uint8_t marker = 128;
		marker += (uint8_t)ordering;

		out.write_little_endian(marker);						// marker for new format

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
		uint8_t marker{};
		in.read_little_endian(marker);

		if (marker < 128)
			sample_id_size_bytes = marker;
		else
		{
			ordering = (ordering_t)(marker - 128);
			in.read_little_endian(sample_id_size_bytes);
		}
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

	void print(std::ostream& oss) const {
		const std::string ordering_name[] = { "SBATC", "ATSBC", "TASBC" };

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
		oss << "\tordering                 : " << ordering_name[(uint8_t)ordering] << "\n";
	}
};

enum class RecFmt { SATC, SPLASH };
struct RecFmtConv {

	inline static RecFmt from_string(const std::string& str) {
		if (str == "satc")
			return RecFmt::SATC;
		else if (str == "splash")
			return RecFmt::SPLASH;
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
		case RecFmt::SPLASH:
			return "splash";
		default:
			std::cerr << "Error: unsupported, switch in RecordPrintFormat should be extended";
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

	size_t store_sample_id(char *ptr, size_t id) const {
		if(can_decode) {
			assert(id < decoded.size());

			strcpy(ptr, decoded[id].c_str());
			return decoded[id].size();
		} else {
			return refresh::int_to_pchar(id, ptr, 0) - 1;
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
	uint64_t sample_id{};
	uint64_t barcode{};
	uint64_t anchor{};
	uint64_t target{};
	uint64_t count{};

	std::vector<uint8_t> compr_record;

private:
	void serialize_msb(std::vector<uint8_t>& out, const Header& header) {
	
		switch (header.ordering)
		{
		case Header::ordering_t::SBATC:
			append_int_msb(out, sample_id, header.sample_id_size_bytes);
			append_int_msb(out, barcode, header.barcode_size_bytes);
			append_int_msb(out, anchor, header.anchor_size_bytes);
			append_int_msb(out, target, header.target_size_bytes);
			append_int_msb(out, count, header.counter_size_bytes);
			break;
		case Header::ordering_t::ATSBC:
			append_int_msb(out, anchor, header.anchor_size_bytes);
			append_int_msb(out, target, header.target_size_bytes);
			append_int_msb(out, sample_id, header.sample_id_size_bytes);
			append_int_msb(out, barcode, header.barcode_size_bytes);
			append_int_msb(out, count, header.counter_size_bytes);
			break;
		case Header::ordering_t::TASBC:
			append_int_msb(out, target, header.target_size_bytes);
			append_int_msb(out, anchor, header.anchor_size_bytes);
			append_int_msb(out, sample_id, header.sample_id_size_bytes);
			append_int_msb(out, barcode, header.barcode_size_bytes);
			append_int_msb(out, count, header.counter_size_bytes);
			break;
		}
	}

	void load_msb(std::vector<uint8_t>& in, const Header& header) {
		auto p = in.begin();

		switch (header.ordering)
		{
		case Header::ordering_t::SBATC:
			LoadBigEndian(p, sample_id, header.sample_id_size_bytes);
			LoadBigEndian(p, barcode, header.barcode_size_bytes);
			LoadBigEndian(p, anchor, header.anchor_size_bytes);
			LoadBigEndian(p, target, header.target_size_bytes);
			LoadBigEndian(p, count, header.counter_size_bytes);
			break;
		case Header::ordering_t::ATSBC:
			LoadBigEndian(p, anchor, header.anchor_size_bytes);
			LoadBigEndian(p, target, header.target_size_bytes);
			LoadBigEndian(p, sample_id, header.sample_id_size_bytes);
			LoadBigEndian(p, barcode, header.barcode_size_bytes);
			LoadBigEndian(p, count, header.counter_size_bytes);
			break;
		case Header::ordering_t::TASBC:
			LoadBigEndian(p, target, header.target_size_bytes);
			LoadBigEndian(p, anchor, header.anchor_size_bytes);
			LoadBigEndian(p, sample_id, header.sample_id_size_bytes);
			LoadBigEndian(p, barcode, header.barcode_size_bytes);
			LoadBigEndian(p, count, header.counter_size_bytes);
			break;
		}
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
		else if (format == RecFmt::SPLASH) {
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
