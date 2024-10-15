#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <istream>
#include <ostream>
#include <cinttypes>

template<typename T>
void write_little_endian(std::ostream& out, const T& data) {
	for (uint32_t b = 0; b < sizeof(data); ++b) {
		uint8_t x = data >> (8 * b);
		out.write((const char*)&x, 1);
	}
}

template<typename T>
void read_little_endian(std::istream& in, T& data) {
	data = T{};
	uint8_t x;
	for (uint32_t b = 0; b < sizeof(data); ++b) {
		in.read((char*)&x, 1);
		data |= ((T)x) << (8 * b);
	}
}
inline void write_string_no_len(std::ostream& out, const std::string& str) {
	out.write(str.c_str(), str.size());
}

inline void write_string_with_len(std::ostream& out, const std::string& str) {
	size_t len = str.length();
	write_little_endian(out, len);
	write_string_no_len(out, str);
}


inline void read_string(std::istream& in, std::string& data, size_t len) {
	data.resize(len, ' ');
	in.read(const_cast<char*>(data.c_str()), len);
}

inline void read_string(std::istream& in, std::string& data) {
	size_t len;
	read_little_endian(in, len);
	read_string(in, data, len);
}

inline bool read_string_verify_len(std::istream& in, std::string& data, size_t the_len_should_be) {
	size_t len;
	read_little_endian(in, len);
	if (len != the_len_should_be)
		return false;
	read_string(in, data, len);
	return true;
}


#endif // !SERIALIZATION_H
