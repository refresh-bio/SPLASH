#ifndef _UTILS_H
#define _UTILS_H
#include <cinttypes>
#include <fstream>
#include <iostream>
#include <refresh/serialization/lib/serialization.h>

namespace lookup_table_utils {
	inline uint32_t bits_required_to_represent(uint32_t val)
	{
		uint32_t res{};
		while (val)
		{
			++res;
			val >>= 1;
		}
		return res;
	}

	//mkokot_TODO: this may be not necessary because probably what is important is serialized by sdsl
	struct ConfigForSBWT {
		size_t bits_for_file_id;
		size_t bits_for_header_id;
		size_t max_cnt;
		size_t bits_for_counter;
		size_t bits_per_elem_in_cat_1;
		size_t bits_per_elem_in_cat_2;
		size_t bits_per_elem_in_cat_3;

		void Serialize(std::ostream& out) const {
			refresh::serialization::serialize_little_endian(bits_for_file_id, out);
			refresh::serialization::serialize_little_endian(bits_for_header_id, out);
			refresh::serialization::serialize_little_endian(max_cnt, out);
			refresh::serialization::serialize_little_endian(bits_for_counter, out);
			refresh::serialization::serialize_little_endian(bits_per_elem_in_cat_1, out);
			refresh::serialization::serialize_little_endian(bits_per_elem_in_cat_2, out);
			refresh::serialization::serialize_little_endian(bits_per_elem_in_cat_3, out);
		}

		void Load(std::istream& in) {
			refresh::serialization::load_little_endian(bits_for_file_id, in);
			refresh::serialization::load_little_endian(bits_for_header_id, in);
			refresh::serialization::load_little_endian(max_cnt, in);
			refresh::serialization::load_little_endian(bits_for_counter, in);
			refresh::serialization::load_little_endian(bits_per_elem_in_cat_1, in);
			refresh::serialization::load_little_endian(bits_per_elem_in_cat_2, in);
			refresh::serialization::load_little_endian(bits_per_elem_in_cat_3, in);
		}

		void Print() const {
			std::cerr << "bits_for_file_id:" << bits_for_file_id << "\n";
			std::cerr << "bits_for_header_id:" << bits_for_header_id << "\n";
			std::cerr << "max_cnt:" << max_cnt << "\n";
			std::cerr << "bits_for_counter:" << bits_for_counter << "\n";
			std::cerr << "bits_per_elem_in_cat_1:" << bits_per_elem_in_cat_1 << "\n";
			std::cerr << "bits_per_elem_in_cat_2:" << bits_per_elem_in_cat_2 << "\n";
			std::cerr << "bits_per_elem_in_cat_3:" << bits_per_elem_in_cat_3 << "\n";
		}
	};
}

#endif // !_UTILS_H