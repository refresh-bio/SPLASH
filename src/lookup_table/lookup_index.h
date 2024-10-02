#ifndef LOOKUP_INDEX_
#define LOOKUP_INDEX_
struct LookupIndex {
	struct elem {
		size_t start{}, end{};
        auto size() const {
            return end - start;
        }
		void Serialize(std::ostream& out) const {
			write_little_endian(out, start);
			write_little_endian(out, end);
		}
        void Load(std::istream& in) {
            read_little_endian(in, start);
            read_little_endian(in, end);
        }
	};
	elem index{};
	elem file_mapper{};
	elem header_mapper{};
	elem sbwt{};
    elem cfg{};
    elem bv_cat1{};
	elem bv_cat2{};
	elem bv_cat3{};
	elem dense_compr_int_vec_cat_1{};
	elem dense_compr_int_vec_cat_2{};
	elem dense_compr_int_vec_cat_3{};
	elem bv_cat3_delim{};

	void Serialize(std::ostream& out) const {
        index.Serialize(out);
		file_mapper.Serialize(out);
		header_mapper.Serialize(out);
		sbwt.Serialize(out);
        cfg.Serialize(out);
        bv_cat1.Serialize(out);
        bv_cat2.Serialize(out);
        bv_cat3.Serialize(out);
        dense_compr_int_vec_cat_1.Serialize(out);
        dense_compr_int_vec_cat_2.Serialize(out);
        dense_compr_int_vec_cat_3.Serialize(out);
        bv_cat3_delim.Serialize(out);
	}

    void Load(std::istream& in) {
        index.Load(in);
		file_mapper.Load(in);
		header_mapper.Load(in);
		sbwt.Load(in);
        cfg.Load(in);
        bv_cat1.Load(in);
        bv_cat2.Load(in);
        bv_cat3.Load(in);
        dense_compr_int_vec_cat_1.Load(in);
        dense_compr_int_vec_cat_2.Load(in);
        dense_compr_int_vec_cat_3.Load(in);
        bv_cat3_delim.Load(in);
    }

};

#endif