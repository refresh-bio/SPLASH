#ifndef _LOOKUP_H
#define _LOOKUP_H
#include "seq_id_mapper.h"
#include "../common/filters/poly_ACGT_filter.h"
#include <string>
#include "utils.h"
#include "kmer_index.h"
#include "kmer_index_sbwt.h"
#include "lookup_index.h"

class Lookup {
	static const std::string MARKER;
	static const uint32_t VERSION_MAJOR;
	static const uint32_t VERSION_MINOR;

	static const uint32_t MIN_SUPPORTED_VERSION_MAJOR;
	static const uint32_t MIN_SUPPORTED_VERSION_MINOR;
	
	SeqIdMapper file_mapper;
	HeaderIdMapper header_mapper;
	PolyACGTFilter poly_ACGT_filter;

	std::unique_ptr<IKmerIndex> kmer_index;

	static std::ifstream open_ifstream(const std::string& path);
	static kmer_index_variant_t read_kmer_index_variant_impl(std::istream& in);

	void TruncatePaths();

public:
	static std::string GetVersion() {
		return std::to_string(VERSION_MAJOR) + "."s + std::to_string(VERSION_MINOR);
	}
	static kmer_index_variant_t ReadKmerIndexVariant(const std::string& path);
	Lookup(const std::string& path, uint32_t n_threads);

	const IKmerIndex* GetKmerIndex() const
	{
		return kmer_index.get();
	}

	void DumpMapping(const std::string& path);

	const SeqIdMapper& GetFileMapper() const { return file_mapper; }

	const HeaderIdMapper& GetHeaderMappers() const { return header_mapper; }

	void query_seq(const std::string& seq, uint32_t kmer_skip, QueryResult& query_result);

	uint64_t GetPolyACGTLen() const { return poly_ACGT_filter.GetLen(); }

	static void Serialize(std::ostream& out,
		uint64_t poly_ACGT_len,
		const SeqIdMapper& file_mapper,
		const HeaderIdMapper& header_mapper,
		kmer_index_variant_t kmer_index_variant,
		const std::function<void(std::ostream&, LookupIndex&)>& serialize_kmer_index);
};
#endif // !_LOOKUP_H
