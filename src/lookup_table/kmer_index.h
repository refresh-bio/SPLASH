#ifndef KMER_INDEX_H
#define KMER_INDEX_H

#include "reporters.h"
#include <iostream>

//at some point we had two variant of kmer index, so I am leaving this code 
//if we want to add new variant at some point
enum class kmer_index_variant_t { sbwt };

inline kmer_index_variant_t kmer_index_variant_from_string(const std::string& str) {
	if (str == "sbwt")
		return kmer_index_variant_t::sbwt;
	std::cerr << "Unknown kmer index variant: " << str << std::endl;
	exit(1);
}

inline std::string to_string(kmer_index_variant_t kmer_index_variant) {
	switch (kmer_index_variant) {
	case kmer_index_variant_t::sbwt:
		return "sbwt";
	}
	std::cerr << "Error:  switch needs to be extended, contact developers " << __FILE__ << ":" << __LINE__ << "\n";
	exit(1);
}

class IKmerIndex
{
public:
	virtual void query(const std::string& seq, uint32_t kmer_skip, QueryResult& query_result) = 0;
	virtual ~IKmerIndex() = default;
};


#endif // !KMER_INDEX_H
