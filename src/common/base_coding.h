#include <algorithm>
#include <string>
#include <algorithm>
#include <cinttypes>
#include <cstring>
#include <vector>

// *********************************************************************************************
// Two bases into single byte
class BaseCoding4
{
	uint8_t base_code[128] = {
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 0-15
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 16-31
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 32-47
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 48-63
		7, 0, 7, 1, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 4, 7,				// 64-79 
		7, 7, 7, 7, 3, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 80-95
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,				// 96-111
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7				// 112-127
	};

	const char *code_base = "ACGTN...........";

	public:
// *********************************************************************************************
	uint64_t encode_bases_2b(const char* begin, const char* end)
	{
		uint64_t code = 0;
		int len = end - begin;
		int shift = 2 * len - 2;

		for (auto p = begin; p != end; ++p, shift -= 2)
		{
			uint64_t c = base_code[(int) *p];
			if (c > 3)
				return ~0ull;

			code += c << shift;
		}

		return code;
	}
	
// *********************************************************************************************
	uint64_t encode_bases_2b(const std::string &str)
	{
		return encode_bases_2b(str.c_str(), str.c_str() + str.size());
	}
	
// *********************************************************************************************
	std::string decode_bases_2b(uint64_t code, uint32_t len)
	{
		std::string str;

		for (int i = 0; i < (int) len; ++i)
		{
			uint8_t c = code & 0x3;

			str.push_back(code_base[c]);
			code >>= 2;
		}

		std::reverse(str.begin(), str.end());

		return str;
	}
};

// *********************************************************************************************
// 3 bases into single byte
class BaseCoding3
{
	uint8_t base_to_code[256];
	char uint8_to_bases[256][3];
	uint8_t uint8_to_codes[256][3];

	void fill_base_coding_tables()
	{
		std::fill_n(base_to_code, 256, 5);
		base_to_code['A'] = 0;
		base_to_code['C'] = 1;
		base_to_code['G'] = 2;
		base_to_code['T'] = 3;
		base_to_code['N'] = 4;

		for (int i = 216; i < 256; ++i)
			uint8_to_bases[i][0] = uint8_to_bases[i][1] = uint8_to_bases[i][2] = 0;

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
		}
	}


public:
	BaseCoding3()
	{
		fill_base_coding_tables();
	}

	// ************************************************************************************
	// Overwrites 2 bytes after the string!
	size_t encode_bases(char* dna_raw, uint8_t* dna_packed)
	{
		size_t len = strlen(dna_raw);

		while (len && (dna_raw[len - 1] == '\n' || dna_raw[len - 1] == '\r'))
			--len;

		dna_raw[len] = 0;

		// Add terminators for easier compression
		dna_raw[len + 1] = 0;
		dna_raw[len + 2] = 0;

		// Rounding up to multiplicity of 3 - just for easier compression
		len = 3 * ((len + 3) / 3);
		size_t packed_len = 0;

		for (size_t i = 0; i < len; i += 3, ++packed_len)
			dna_packed[packed_len] =
				base_to_code[static_cast<uint8_t>(dna_raw[i])] * 36 +
				base_to_code[static_cast<uint8_t>(dna_raw[i + 1])] * 6 +
				base_to_code[static_cast<uint8_t>(dna_raw[i + 2])];

		return packed_len;
	}

	// ************************************************************************************
	size_t encode_bases(char* dna_raw, size_t len, uint8_t* dna_packed)
	{
		size_t packed_len = 0;
		size_t i;
		++len;						// terminating 0 is also compressed

		for (i = 0; i + 2 < len; i += 3, ++packed_len)
			dna_packed[packed_len] =
				base_to_code[static_cast<uint8_t>(dna_raw[i])] * 36 +
				base_to_code[static_cast<uint8_t>(dna_raw[i + 1])] * 6 +
				base_to_code[static_cast<uint8_t>(dna_raw[i + 2])];

		if(i + 1 == len)
			dna_packed[packed_len++] =
				base_to_code[static_cast<uint8_t>(dna_raw[i])] * 36 +
				base_to_code[0] * 6 +
				base_to_code[0];
		else if(i + 2 == len)
			dna_packed[packed_len++] =
				base_to_code[static_cast<uint8_t>(dna_raw[i])] * 36 +
				base_to_code[static_cast<uint8_t>(dna_raw[i + 1])] * 6 +
				base_to_code[0];
/*		else
			dna_packed[packed_len++] = base_to_code[0] * 36 + base_to_code[0] * 6 + base_to_code[0];*/

		return packed_len;
	}

	// ************************************************************************************
	void decode_bases(uint8_t* dna_packed, std::vector<uint8_t>& dna_raw)
	{
		dna_raw.clear();

		for(; *dna_packed % 6 != 5; ++dna_packed)
		{
			auto& x = uint8_to_bases[*dna_packed];
			dna_raw.emplace_back(x[0]);
			dna_raw.emplace_back(x[1]);
			dna_raw.emplace_back(x[2]);
		}

		auto& x = uint8_to_bases[*dna_packed];

		dna_raw.emplace_back(x[0]);
		if (x[0] == 0)
			return;
		dna_raw.emplace_back(x[1]);
		if (x[1] == 0)
			return;
		dna_raw.emplace_back(x[2]);
	}
};

// EOF