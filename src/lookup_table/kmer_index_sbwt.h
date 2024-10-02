#ifndef _KMER_INDEX_SBWT_H
#define _KMER_INDEX_SBWT_H

#ifdef _MSC_VER
#pragma warning(disable:4996)
#pragma warning(disable:4146)
#endif // _MSC_VER

#include "kmer_index.h"
#include "utils.h"
#include "seq_id_mapper.h"
#include "transcriptome_kmer_finder.h"
#include "../common/filters/poly_ACGT_filter.h"
#include <refresh/parallel_queues/lib/parallel-queues.h>
#include "../common/thread-control.h"
#include <refresh/archive/lib/archive.h>
#include "task_desc.h"
#include "sbwt/variants.hh"
#include "dense_compr_int_vector.h"
#include "../common/binary_heap_merge.h"
#include "lookup_index.h"

class KmerIndex_SBWT : public IKmerIndex
{
	PolyACGTFilter& poly_ACGT_filter;
	HeaderIdMapper& header_mapper;

	//mkokot_TODO: rrr_matrix_sbwt_t variant? maybe configure this
	sbwt::plain_matrix_sbwt_t sbwt;

	//rand support or maybe change to different representation?
	sdsl::bit_vector bv_cat1;
	sdsl::bit_vector bv_cat2;
	sdsl::bit_vector bv_cat3;

	sdsl::rank_support_v5<> rs_bv_cat1;
	sdsl::rank_support_v5<> rs_bv_cat2;
	sdsl::rank_support_v5<> rs_bv_cat3;

	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_1;
	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_2;
	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_3;

	sdsl::bit_vector bv_cat3_delim;

	sdsl::select_support_mcl<> ss_bv_cat3_delim;

	lookup_table_utils::ConfigForSBWT cfg;
public:
	KmerIndex_SBWT(PolyACGTFilter& poly_ACGT_filter, HeaderIdMapper& header_mapper, const std::string& path, const LookupIndex& lookup_index, std::vector<std::function<void()>>& load_tasks) :
		poly_ACGT_filter(poly_ACGT_filter),
		header_mapper(header_mapper) {

		auto load_report_size = [](const std::string& path, auto& to_load, const std::string& name, const LookupIndex::elem& index_elem)
			{

				std::ifstream in(path, std::ios::binary);
				if (!in) {
					std::cerr << "Error: cannot open file " << path << " to read " << name << "\n";
					exit(1);
				}

				in.seekg(index_elem.start);
				to_load.load(in);

				assert((size_t)in.tellg() == index_elem.end);
				std::cerr << (name + " loaded. Size: " + std::to_string(index_elem.size()) + " bytes\n");
			};


		load_tasks.emplace_back([&, path]{
			std::ifstream in(path, std::ios::binary);
			if (!in) {
				std::cerr << "Error: cannot open file " << path << " to read sbwt\n";
				exit(1);
			}
			in.seekg(lookup_index.sbwt.start);

			auto sbwt_variant = sbwt::load_string(in);
			std::cerr << "sbwt_variant: " + sbwt_variant + "\n";
			if (sbwt_variant != "plain-matrix") {
				std::cerr << "Error: for now only plain-matrix is supported matrix type\n"; //mkokot_TODO: add support for the rest of variants
				exit(1);
			}
			sbwt.load(in);
			assert((size_t)in.tellg() == lookup_index.sbwt.end);
		});

		load_tasks.emplace_back([&, path]{
			std::ifstream in(path, std::ios::binary);
			if (!in) {
				std::cerr << "Error: cannot open file " << path << " to read sbwt\n";
				exit(1);
			}
			in.seekg(lookup_index.cfg.start);
			cfg.Load(in);
			assert((size_t)in.tellg() == lookup_index.cfg.end);
		});

#define LOAD_WRAPPER(TO_LOAD) load_report_size(path, TO_LOAD, #TO_LOAD, lookup_index.TO_LOAD)

		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(bv_cat1);
			rs_bv_cat1 = sdsl::rank_support_v5<>(&bv_cat1);
		});

		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(bv_cat2);
			rs_bv_cat2 = sdsl::rank_support_v5<>(&bv_cat2);
		});

		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(bv_cat3);
			rs_bv_cat3 = sdsl::rank_support_v5<>(&bv_cat3);
		});

		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(dense_compr_int_vec_cat_1);
		});
		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(dense_compr_int_vec_cat_2);
		});
		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(dense_compr_int_vec_cat_3);
		});
		load_tasks.emplace_back([&, path]{
			LOAD_WRAPPER(bv_cat3_delim);
			ss_bv_cat3_delim = sdsl::select_support_mcl<>(&bv_cat3_delim);
		});

#undef LOAD_WRAPPER
	}

	void query(const std::string& seq, uint32_t kmer_skip, QueryResult& query_result) override;

	void print_config() const {
		cfg.Print();
	}
};




//mkokot_TODO: move this somewhere else
//thread safe
class archive_store
{
	refresh::archive_buffered_io archive;
	size_t single_pack_size_bytes;
	std::atomic<size_t> n_streams{};
	std::atomic<size_t> n_elems{};

	template<typename REC_T>
	size_t get_serializd_rec_size(const REC_T& rec)
	{
		std::vector<uint8_t> tmp;
		refresh::serialization::serialize_little_endian(rec, tmp);
		return tmp.size();
	}
public:
	archive_store(const std::string& name, size_t single_pack_size_bytes = 1ull << 21) :
		archive(false),
		single_pack_size_bytes(single_pack_size_bytes) {
		if (!archive.open(name)) {
			std::cerr << "Error: cannt open archive " << name << "\n";
			exit(1);
		}
	}

	template<typename REC_T>
	void Add(const std::vector<REC_T>& data) {
		if (data.empty())
			return;

		auto id = archive.register_stream(std::to_string(n_streams++));
		n_elems += data.size();

		auto serialized_rec_size_bytes = get_serializd_rec_size(data[0]);

		std::vector<uint8_t> serialized;
		serialized.reserve(single_pack_size_bytes / serialized_rec_size_bytes * serialized_rec_size_bytes);

		for (const auto& rec : data) {
			refresh::serialization::serialize_little_endian(rec, serialized);
			if (serialized.size() == serialized.capacity())
			{
				archive.add_part(id, serialized, 0);
				serialized.clear();
			}
		}
		if (!serialized.empty())
			archive.add_part(id, serialized, 0);
	}

	void Close()
	{
		archive.close();
	}

};

class archive_stream_rec_reader {
	refresh::archive_buffered_io& archive;
	int id;
	std::vector<uint8_t> data;
	size_t pos{};
public:
	archive_stream_rec_reader(refresh::archive_buffered_io& archive, const std::string& stream_name) :
		archive(archive),
		id(archive.get_stream_id(stream_name))
	{

	}

	template<typename REC_T>
	bool NextRec(REC_T& rec)
	{
		assert(pos <= data.size());
		if (pos == data.size())
		{
			uint64_t meta;
			if (!archive.get_part(id, data, meta))
				return false;
			pos = 0;
		}
		refresh::serialization::load_little_endian(rec, data, pos);
		return true;
	}
};


template<typename REC_T>
class archive_stream_merge {
	std::string archive_path;
	refresh::archive_buffered_io archive;
	std::vector<archive_stream_rec_reader> rec_readers;
	BinaryHeapMergeStreams<REC_T> heap;

	template<typename CALLBACK_T>
	bool do_with_elem_if_exists(size_t id, const CALLBACK_T& callback)
	{
		REC_T rec;
		if (rec_readers[id].NextRec(rec))
		{
			callback(rec);
			return true;
		}
		return false;
	}


	std::vector<archive_stream_rec_reader> make_rec_readers(const std::string& name)
	{
		if (!archive.open(name)) {
			std::cerr << "Error: cannot open archive file " << name << " for reading\n";
			exit(1);
		}
		auto n_streams = archive.get_no_streams();
		std::vector<archive_stream_rec_reader> res;
		res.reserve(n_streams);
		for (size_t i = 0; i < n_streams; ++i)
			res.emplace_back(archive, std::to_string(i));

		//std::cerr << "number of streams in "s + name + ": " << res.size() + "\n";
		return res;
	}
public:
	archive_stream_merge(const std::string& name) :
		archive_path(name),
		archive(true),
		rec_readers(make_rec_readers(name)),
		heap(rec_readers.size(), [this](size_t id, const auto& callback) { return do_with_elem_if_exists(id, callback); })
	{

	}

	bool NextRec(REC_T& rec)
	{
		if (heap.Empty())
			return false;

		heap.ProcessElem([this](size_t id, const auto& callback) { return do_with_elem_if_exists(id, callback); }, [&](const REC_T& _rec, size_t id)
			{
				rec = _rec;
			});
		return true;
	}

	void CloseAndRemove()
	{
		archive.close();
		std::filesystem::remove(archive_path);
	}
};


//It does not always build sbwt, sometimes just load existing
//in fact it builds sbwt based kmer index
class KmerIndex_SBWTBuilder {

	//mkokot_TODO: rrr_matrix_sbwt_t variant? maybe configure this
	struct sbwt_wrapper
	{
		std::string variant;
		sbwt::plain_matrix_sbwt_t sbwt;

		sbwt_wrapper(const std::string& variant, sbwt::plain_matrix_sbwt_t&& sbwt):
			variant(variant),
			sbwt(std::move(sbwt)) {

		}
	};

	sbwt_wrapper sbwt;

	sdsl::bit_vector bv_cat1;
	sdsl::bit_vector bv_cat2;
	sdsl::bit_vector bv_cat3;

	refresh::thread_control& tc;

	lookup_table_utils::ConfigForSBWT cfg;

	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_1;
	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_2;
	dense_compr_int_vector<uint32_t> dense_compr_int_vec_cat_3;

	sdsl::bit_vector bv_cat3_delim;

	//queue of packs to compute sbwt index
	refresh::parallel_queue<std::vector<task_desc>> sbwt_queue;

	//queue of packs with computed sbwt, the thread
	//will set appropriate bitvectors and add data to vectors
	//used later to create sdsl-lite int_vectors
	refresh::parallel_queue<std::vector<task_desc>> fill_bvs_queue;

	struct sort_and_store_task_data
	{
		std::vector<std::pair<int64_t, uint32_t>> cat_1_data; //idx, header_id
		std::vector<std::tuple<int64_t, uint32_t, uint32_t>> cat_2_data; //idx, file_id, cnt
		std::vector<std::pair<int64_t, uint32_t>> cat_3_data; //in case of cat 3 we may have more entries per single k-mer, but we dont know this during iterations, will need to merge after . idx, header_i

		sort_and_store_task_data() = default;
		sort_and_store_task_data(
			std::vector<std::pair<int64_t, uint32_t>>&& cat_1_data,
			std::vector<std::tuple<int64_t, uint32_t, uint32_t>>&& cat_2_data,
			std::vector<std::pair<int64_t, uint32_t>> cat_3_data
			):
			cat_1_data(std::move(cat_1_data)),
			cat_2_data(std::move(cat_2_data)),
			cat_3_data(std::move(cat_3_data)) {
		}

	};

	std::string archive_file_name;

	//queue of data to be sorted and stored in temp files because keeping it all in memory is quite consuming
	refresh::parallel_queue<sort_and_store_task_data> sort_and_store_queue;

	//mkokot_TODO: consider if this is needed!
	struct
	{
		std::atomic<size_t> cat_1_elems{};
		std::atomic<size_t> cat_2_elems{}; //before global unique!!!
		std::atomic<size_t> cat_3_elems{};
	} temp_data_stats;
	//
	//refresh::archive_buffered_io archive;

	archive_store archive_store_cat_1;
	archive_store archive_store_cat_2;
	archive_store archive_store_cat_3;

	//mkokot_TODO: remove?
	//std::vector<std::pair<int64_t, uint32_t>> cat_1_data; //idx, header_id
	//std::vector<std::tuple<int64_t, uint32_t, uint32_t>> cat_2_data; //idx, file_id, cnt
	//std::vector<std::pair<int64_t, uint32_t>> cat_3_data; //in case of cat 3 we may have more entries per single k-mer, but we dont know this during iterations, will need to merge after . idx, header_i

	//std::atomic<size_t> call_no{}; //for progress
	std::vector<std::thread> threads;

	size_t max_cnt = 0;

	template<typename DATA_T>
	void put_in_queue(refresh::parallel_queue<DATA_T>& queue, DATA_T&& pack) {
		//this prevents hanging on push by allowing to start other thread
		//at some point this other thread will need to be consumer

		tc.hang_if_needed([&queue, &pack] {
			queue.push(std::move(pack));
		});
	}

	int64_t get_index(size_t kmer_len, uint32_t rc_shift, uint64_t canonical_kmer,
	                  std::string& str_kmer_data);

	sbwt_wrapper load_sbwt(const std::string& path);

	sbwt_wrapper build_sbwt(
		const std::vector<std::string>& input_paths,
		size_t kmer_len,
		size_t n_threads,
		const std::string& tmp_dir);

	void start_threads();

	//common ctor
	KmerIndex_SBWTBuilder(sbwt_wrapper&& sbwt, refresh::thread_control& tc, const std::string& tmp_dir) :
		sbwt(std::move(sbwt)),
		bv_cat1(sbwt.sbwt.number_of_subsets()),
		bv_cat2(sbwt.sbwt.number_of_subsets()),
		bv_cat3(sbwt.sbwt.number_of_subsets()),
		tc(tc),
		sbwt_queue( tc.max_running() + 1, 1, "sbwt_queue"),
		fill_bvs_queue(tc.max_running() + 1, tc.max_running(), "fill_bvs_queue"),
		archive_file_name((std::filesystem::path(tmp_dir) / "streams").string()),
		sort_and_store_queue(1, 1, "sort_and_store_queue"),
		archive_store_cat_1(archive_file_name + "_cat1.bin"),
		archive_store_cat_2(archive_file_name + "_cat2.bin"),
		archive_store_cat_3(archive_file_name + "_cat3.bin")
		{
			start_threads();
		}
public:
	static kmer_index_variant_t GetIndexVariant() { return kmer_index_variant_t::sbwt; }
	KmerIndex_SBWTBuilder(
		const std::string& path,
		refresh::thread_control& tc,
		const std::string& tmp_dir) :
		KmerIndex_SBWTBuilder(load_sbwt(path), tc, tmp_dir) {

	}

	KmerIndex_SBWTBuilder(
		const std::vector<std::string>& input_paths,
		size_t kmer_len,
		refresh::thread_control& tc,
		const std::string& tmp_dir):
		KmerIndex_SBWTBuilder(build_sbwt(input_paths, kmer_len, tc.max_running(), tmp_dir), tc, tmp_dir) {

	}

	void Build(size_t file_mapper_size,
		size_t header_mapper_size);

	
	void AddPack(std::vector<task_desc>&& pack) {
		put_in_queue(sbwt_queue, std::move(pack));
	}

	size_t GetInCat1() const { return this->temp_data_stats.cat_1_elems; }
	size_t GetInCat2() const { return this->temp_data_stats.cat_2_elems; }
	size_t GetInCat3() const { return this->temp_data_stats.cat_3_elems; }

	//size_t GetInCat1() const { return cat_1_data.size(); }
	//size_t GetInCat2() const { return cat_2_data.size(); }
	//size_t GetInCat3() const { return cat_3_data.size(); }

	void Serialize(std::ostream& out, LookupIndex& lookup_index) const;
};
#endif // !_KMER_INDEX_SBWT_H
