#ifndef _SEQ_ID_MAPPER
#define _SEQ_ID_MAPPER
#include <vector>
#include <cassert>
#include <string>
#include <numeric>
//#include <iostream> //mkokot_TODO: remove
#include "serialization.h"
#include <limits>
#include <fstream>
#include <functional>



class SeqIdMapper {
	std::vector<std::string> M;
public:
	SeqIdMapper() = default;

	uint32_t Add(const std::string& desc) {
		assert(M.size() < std::numeric_limits<uint32_t>::max());
		M.push_back(desc);
		return M.size() - 1;
	}

	void Dump(std::ostream& out) const {
		for (size_t i = 0; i < M.size(); ++i)
			out << i << "\t" << M[i] << "\n";
	}

	const std::string& Decode(uint32_t id) const {
		return M[id];
	}

	size_t Size() const {
		return M.size();
	}

	void Serialize(std::ostream& out) const {
		write_little_endian(out, M.size());
		for (const auto& x : M)
			write_string_with_len(out, x);
	}

	void Load(std::istream& in) {
		size_t size;
		read_little_endian(in, size);
		M.resize(size);
		for (size_t i = 0; i < M.size(); ++i)
			read_string(in, M[i]);
	}

	void TransformNames(const std::function<void(std::string& /*name*/)>& transform) {
		for (auto& x : M)
			transform(x);
	}
};

class HeaderIdMapper {
	struct SeqAndId {
		std::string header;
		uint32_t file_id;
		SeqAndId(const std::string& header, uint32_t file_id) :
			header(header),
			file_id(file_id)
		{

		}
		SeqAndId() = default;
	};
	std::vector<SeqAndId> M;
public:
	HeaderIdMapper() = default;

	uint32_t Add(const std::string& header, uint32_t file_id) {
		assert(M.size() < std::numeric_limits<uint32_t>::max());
		M.emplace_back(header, file_id);
		return M.size() - 1;
	}

	void Dump(std::ostream& out, const SeqIdMapper& file_mapper) const {
		for (size_t i = 0; i < M.size(); ++i)
			out << i << "\t" << file_mapper.Decode(M[i].file_id) << "\t" << M[i].header << "\n";
	}

	const SeqAndId& Decode(uint32_t id) const {
		return M[id];
	}

	size_t Size() const {
		return M.size();
	}

	void Compact() {
		if (M.empty())
			return;

		auto substr_to_first_whitespace = [](const std::string& str) {
			return str.substr(0, str.find_first_of(" \t"));
		};
		size_t org_size{};
		size_t after_compact_size{};

		for (auto& e : M) {
			org_size += e.header.length();
			e.header = substr_to_first_whitespace(e.header);
			after_compact_size += e.header.length();
		}

		//std::cerr << "org_size: " << org_size << "\n";
		//std::cerr << "after_compact_size: " << after_compact_size << "\n";

	}

	void Serialize(std::ostream& out) const {
		write_little_endian(out, M.size());
		for (const auto& x : M) {
			write_string_with_len(out, x.header);
			write_little_endian(out, x.file_id);
		}
	}

	void Load(std::istream& in) {
		size_t size;
		read_little_endian(in, size);
		M.resize(size);
		for (size_t i = 0; i < M.size(); ++i) {
			read_string(in, M[i].header);
			read_little_endian(in, M[i].file_id);
		}
	}
};

#endif // !_SEQ_ID_MAPPER

