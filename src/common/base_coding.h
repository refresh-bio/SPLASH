#include <algorithm>
#include <string>

class BaseCoding
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

// EOF