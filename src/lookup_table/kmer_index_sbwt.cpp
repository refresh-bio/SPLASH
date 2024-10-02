#include "kmer_index_sbwt.h"
#include "walk_kmers.h"

void KmerIndex_SBWT::query(const std::string& seq, uint32_t kmer_skip, QueryResult& query_result)
{
	query_result.Clear();
	uint32_t kmer_len = sbwt.get_k();
	if (seq.length() < kmer_len) {
		return;
	}

	//mkokot_TODO: how will this behave for N?
	std::vector<int64_t> idxs = sbwt.streaming_search(seq);

	//fix reverse complements
	{
		auto get_rc = [](const std::string& seq) {
			auto rev_compl_c = [](char c) {
				switch (c) {
				case 'A': return 'T';
				case 'C': return 'G';
				case 'G': return 'C';
				case 'T': return 'A';
				default: return c;
				}
				};
			std::string res(seq.size(), ' ');
			for (size_t i = 0; i < seq.size(); ++i) {
				res[seq.size() - i - 1] = rev_compl_c(seq[i]);
			}
			return res;
			};
		auto seq_rc = get_rc(seq);
		std::vector<int64_t> rc_idxs = sbwt.streaming_search(seq_rc);
		std::reverse(rc_idxs.begin(), rc_idxs.end());
		for (size_t i = 0; i < idxs.size(); ++i) {
			if (idxs[i] == -1)
				idxs[i] = rc_idxs[i];
			else if (rc_idxs[i] != -1) {
				//in this case I should use the idx of canonical k-mer
				char* kmer = (char*)seq.data() + i;
				char* rev_compl = seq_rc.data() + seq_rc.size() - kmer_len - i;
				//mkokot_TODO: rewrite this to work not on strings?
				auto is_kmer_canonical = [](char* kmer, char* rev_compl, size_t k) {
					for (size_t i = 0; i < k; ++i) {
						if (kmer[i] < rev_compl[i])
							return true;
						if (kmer[i] > rev_compl[i])
							return false;
					}
					return true;
					};
				if (!is_kmer_canonical(kmer, rev_compl, kmer_len))
					idxs[i] = rc_idxs[i];
			}
		}
	}

	WalkKmersWithInvalid walk_kmers(seq, kmer_len, true);
	uint64_t kmer{};
	bool is_valid;

	std::vector<kmer_desc> V;

	size_t bits_for_counter = lookup_table_utils::bits_required_to_represent(cfg.max_cnt);

	for (size_t i = 0; i < idxs.size(); ++i) {
		V.clear();
		if (!walk_kmers.Next(kmer, is_valid)) {
			std::cerr << "Error: please contact authors, where: " << __FILE__ << "(" << __LINE__ << ")";
			exit(1);
		}

		if (i % (kmer_skip + 1))
			continue;

		if (!is_valid) {
			query_result.ReportIllegal(); //k-mer contains N or other non-ACGT symbol
		}
		else if (poly_ACGT_filter.IsPolyACGT(kmer, kmer_len)) {
			query_result.ReportPolyACGT();
		}
		else {
			auto idx = idxs[i];
			if (idx == -1) {
				query_result.ReportUnknown();
			}
			else {
				if (bv_cat1[idx]) {
					auto rank = rs_bv_cat1.rank(idx);
					auto header_id = dense_compr_int_vec_cat_1[rank];
					V.emplace_back(1, header_mapper.Decode(header_id).file_id, header_id);

				}
				else if (bv_cat2[idx]) {
					auto rank = rs_bv_cat2.rank(idx);
					auto file_id_and_cnt = dense_compr_int_vec_cat_2[rank];
					uint32_t file_id = file_id_and_cnt >> bits_for_counter;
					uint32_t cnt = file_id_and_cnt & ((1ul << bits_for_counter) - 1);
					V.emplace_back(2, file_id, cnt);
				}
				else if (bv_cat3[idx]) {
					auto rank = rs_bv_cat3.rank(idx);
					size_t start = ss_bv_cat3_delim.select(rank + 1);
					size_t end = ss_bv_cat3_delim.select(rank + 2);

					while (start < end) {
						auto header_id = dense_compr_int_vec_cat_3[start++];
						//mkokot_TODO: I think here I don't need to decode and store file_id, because it may be do later from header_id
						//although it is used for sorting, so need to be reconsidered
						V.emplace_back(3, header_mapper.Decode(header_id).file_id, header_id);
					}
				}
				else {
					query_result.ReportCategory4();
				}

				if (V.size()) {
					auto category = V[0].category;
					for (const auto& x : V)
						if (x.category != category) {
							std::cerr << "Error: for a k-mer " << kmer_to_string(kmer, kmer_len) << " there were more than one category, which should not happen and means there is a bug\n";
							exit(1);
						}

					if (category == 1) {
						assert(V.size() == 1);
						query_result.ReportCategory1(V[0].meta);
					}
					else if (category == 2) {
						assert(V.size() == 1);
						query_result.ReportCategory2(V[0].file_id, V[0].meta);
					}
					else if (category == 3) {
						std::sort(V.begin(), V.end(), [](const kmer_desc& lhs, const kmer_desc& rhs) { return std::make_pair(lhs.file_id, lhs.meta) < std::make_pair(rhs.file_id, rhs.meta); });
						query_result.ReportCategory3(V);
					}
					else {
						std::cerr << "Error: unknown k-mer (" << kmer_to_string(kmer, kmer_len) << ") category: " << category << "\n";
						exit(1);
					}
				}
			}
		}
	}
}

int64_t KmerIndex_SBWTBuilder::get_index(size_t kmer_len, uint32_t rc_shift,
	uint64_t canonical_kmer, std::string& str_kmer_data)
{

	//auto p = call_no++;
	//if (p % 100'000'000 == 0)
	//	std::cerr << p << "\n";

	kmer_to_string(canonical_kmer, kmer_len, str_kmer_data);

	auto r = sbwt.sbwt.search(str_kmer_data);
	if (r != -1)
		return r;

	auto rev_compl = get_rev_compl(canonical_kmer, rc_shift);
	kmer_to_string(rev_compl, kmer_len, str_kmer_data);
	r = sbwt.sbwt.search(str_kmer_data);
	if (r != -1)
		return r;

	std::cerr << "Error: cannot find k-mer (nor its rev compl) " << str_kmer_data << " in the sbwt\n";
	exit(1);
}

KmerIndex_SBWTBuilder::sbwt_wrapper KmerIndex_SBWTBuilder::load_sbwt(const std::string& path)
{
	std::ifstream in(path, std::ios::binary);
	std::string sbwt_variant = sbwt::load_string(in);
	std::cerr << "variant: " << sbwt_variant << "\n";
	if (sbwt_variant != "plain-matrix") {
		std::cerr << "Error: for now only plain-matrix is supported matrix type\n"; //mkokot_TODO: add support for the rest of variants
		exit(1);
	}
	sbwt::plain_matrix_sbwt_t x;
	x.load(in);
	return { sbwt_variant, std::move(x) };
}

KmerIndex_SBWTBuilder::sbwt_wrapper KmerIndex_SBWTBuilder::build_sbwt(
	const std::vector<std::string>& input_paths,
	size_t kmer_len, size_t n_threads,
	const std::string& tmp_dir)
{
	//mkokot_TODO: support other variants?
	//mkokot_TODO: forcing, but maybe should be possible to set this
	//sbwt_variant = "plain-matrix";
	sbwt::plain_matrix_sbwt_t::BuildConfig config;
	config.input_files = input_paths;
	config.k = kmer_len;

	config.build_streaming_support = true;
	config.n_threads = n_threads;
	//		config.min_abundance = min_abundance;
	//		config.max_abundance = max_abundance;
	config.ram_gigas = 12;
	config.temp_dir = tmp_dir;
	//		config.precalc_k = 0;

	sbwt::get_temp_file_manager().set_dir(config.temp_dir);

	sbwt::plain_matrix_sbwt_t sbwt(config);

	auto precalc_length = 8; //mkokot_TODO: make parameter?
	sbwt.do_kmer_prefix_precalc(precalc_length);

	return { "plain-matrix", std::move(sbwt) };
}

void KmerIndex_SBWTBuilder::start_threads()
{
	//threads computing sbwt indexes
	for (size_t tid = 0; tid < tc.max_running(); ++tid)
		threads.emplace_back([&] {
			size_t kmer_len = sbwt.sbwt.get_k();
			uint32_t rc_shift = get_rev_compl_shift(kmer_len);

			std::vector<task_desc> tasks;
			std::string str_kmer_data(kmer_len, ' ');
			while (sbwt_queue.pop(tasks)) {
				tc.execute([&] {
					for (auto& x : tasks) {
						x.idx = get_index(kmer_len, rc_shift, x.kmer, str_kmer_data);
					}

					put_in_queue(fill_bvs_queue, std::move(tasks));
				});
			}

			fill_bvs_queue.mark_completed();
		});

	//thread that fills bitvectors
	threads.emplace_back([this] {
		std::vector<task_desc> tasks;

		//size_t max_vec_size_bytes = 1ull << 30;
		size_t max_vec_size_bytes = 1ull << 27;

		auto calc_capacity = [max_vec_size_bytes](const auto& vec)
		{
			return max_vec_size_bytes / sizeof(vec[0]);
		};

		std::vector<std::pair<int64_t, uint32_t>> cat_1_data; //idx, header_id
		std::vector<std::tuple<int64_t, uint32_t, uint32_t>> cat_2_data; //idx, file_id, cnt
		std::vector<std::pair<int64_t, uint32_t>> cat_3_data; //in case of cat 3 we may have more entries per single k-mer, but we dont know this during iterations, will need to merge after . idx, header_i

		cat_1_data.reserve(calc_capacity(cat_1_data));
		cat_2_data.reserve(calc_capacity(cat_2_data));
		cat_3_data.reserve(calc_capacity(cat_3_data));

		auto send_if_needed_cat_1 = [&cat_1_data, &calc_capacity, this] {
			if (cat_1_data.size() < cat_1_data.capacity())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data(std::move(cat_1_data), {}, {} ));
			cat_1_data.reserve(calc_capacity(cat_1_data));
		};

		auto send_if_needed_cat_2 = [&cat_2_data, &calc_capacity, this] {
			if (cat_2_data.size() < cat_2_data.capacity())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data({}, std::move(cat_2_data), {}));
			cat_2_data.reserve(calc_capacity(cat_2_data));
		};

		auto send_if_needed_cat_3 = [&cat_3_data, &calc_capacity, this] {
			if (cat_3_data.size() < cat_3_data.capacity())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data({}, {}, std::move(cat_3_data)));
			cat_3_data.reserve(calc_capacity(cat_3_data));
		};

		auto send_last_cat_1 = [&cat_1_data, this] {
			if (cat_1_data.empty())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data(std::move(cat_1_data), {}, {}));
		};

		auto send_last_cat_2 = [&cat_2_data, this] {
			if (cat_2_data.empty())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data({}, std::move(cat_2_data), {}));
		};

		auto send_last_cat_3 = [&cat_3_data, this] {
			if (cat_3_data.empty())
				return;

			put_in_queue(sort_and_store_queue, sort_and_store_task_data({}, {}, std::move(cat_3_data)));
		};

		while (fill_bvs_queue.pop(tasks)) {
			tc.execute([&] {
				for (const auto& x : tasks) {
					if (x.cat == 1) {
						bv_cat1[x.idx] = 1;
						cat_1_data.emplace_back(x.idx, x.header_id);
						send_if_needed_cat_1();
					}
					else if (x.cat == 2) {
						if (x.cnt > max_cnt)
							max_cnt = x.cnt;
						bv_cat2[x.idx] = 1;
						cat_2_data.emplace_back(x.idx, x.file_id, x.cnt);
						send_if_needed_cat_2();
					}
					else if (x.cat == 3) {
						bv_cat3[x.idx] = 1;
						cat_3_data.emplace_back(x.idx, x.header_id);
						send_if_needed_cat_3();
					}
				}
			});
		}

		send_last_cat_1();
		send_last_cat_2();
		send_last_cat_3();

		sort_and_store_queue.mark_completed();
	});


	//threads computing sorting and storing to temp
	for (size_t tid = 0; tid < tc.max_running(); ++tid)
		threads.emplace_back([&] {

			sort_and_store_task_data data;
			while (sort_and_store_queue.pop(data)) {
				tc.execute([&] {
					//mkokot_INFO: I tried replacing std::sort with block quick sort but the gain was neglible (at least when using 64 threads)
					//so to decrease the number of dependant GPL3 code I use std::sort
					if (!data.cat_1_data.empty()) {
						temp_data_stats.cat_1_elems += data.cat_1_data.size();
						std::sort(data.cat_1_data.begin(), data.cat_1_data.end());

						archive_store_cat_1.Add(data.cat_1_data);
					}
					if (!data.cat_2_data.empty()) {
						std::sort(data.cat_2_data.begin(), data.cat_2_data.end());

						//I deduplicate cat 2 because I store this in parts
						//there are still some duplicates possible that needs to be handled on merge!!!!
						data.cat_2_data.erase(std::unique(data.cat_2_data.begin(), data.cat_2_data.end()), data.cat_2_data.end());

						temp_data_stats.cat_2_elems += data.cat_2_data.size();

						archive_store_cat_2.Add(data.cat_2_data);
					}
					if (!data.cat_3_data.empty()) {
						temp_data_stats.cat_3_elems += data.cat_3_data.size();

						std::sort(data.cat_3_data.begin(), data.cat_3_data.end(), std::less<>{});

						archive_store_cat_3.Add(data.cat_3_data);
					}
				});
			}
		});
}



void KmerIndex_SBWTBuilder::Build(size_t file_mapper_size,
                                  size_t header_mapper_size)
{
	//Finish
	sbwt_queue.mark_completed();

	std::cerr << "SBWT builder waits for threads completition\n";
	for (auto& th : threads)
		th.join();

	std::cerr << "Threads completed\n";

	//close all archives
	archive_store_cat_1.Close();
	archive_store_cat_2.Close();
	archive_store_cat_3.Close();


	//Build sdsl int_vectors
	std::cerr << "INFO: Sorting categories data\n";

	//now we may create for each category int_vector with required bits to store values

	cfg.bits_for_file_id = lookup_table_utils::bits_required_to_represent(file_mapper_size);
	cfg.bits_for_header_id = lookup_table_utils::bits_required_to_represent(header_mapper_size);
	cfg.max_cnt = max_cnt;

	std::cerr << "INFO: max cnt: " << cfg.max_cnt << "\n";

	cfg.bits_for_counter = lookup_table_utils::bits_required_to_represent(cfg.max_cnt);

	cfg.bits_per_elem_in_cat_1 = cfg.bits_for_header_id; //header id
	cfg.bits_per_elem_in_cat_2 = cfg.bits_for_file_id + cfg.bits_for_counter; //file_id + cnt
	cfg.bits_per_elem_in_cat_3 = cfg.bits_for_header_id;

	std::cerr << "INFO: create int_vectors\n";

	//std::cerr << "temp_data_stats.cat_1_elems: " << temp_data_stats.cat_1_elems << "\n";
	//std::cerr << "temp_data_stats.cat_2_elems: " << temp_data_stats.cat_2_elems << "\n";
	//std::cerr << "temp_data_stats.cat_3_elems: " << temp_data_stats.cat_3_elems << "\n";

	auto start_time = std::chrono::high_resolution_clock::now();
	std::vector<std::thread> merging_threads;

	merging_threads.emplace_back([this]{
		tc.execute([this]{

			auto start_time = std::chrono::high_resolution_clock::now();

			sdsl::int_vector<> int_vec_cat_1(this->temp_data_stats.cat_1_elems, 0, cfg.bits_per_elem_in_cat_1);
			using cat_1_rec_t = std::pair<int64_t, uint32_t>;
			archive_stream_merge<cat_1_rec_t> cat_1_merge(archive_file_name + "_cat1.bin");
			cat_1_rec_t cat1_rec;
			size_t i = 0;
			while(cat_1_merge.NextRec(cat1_rec))
			{
				int_vec_cat_1[i++] = cat1_rec.second;
			}

			auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
			std::cerr << "int_vec_cat_1 created. Time: "s + std::to_string(time) + ", size: " + std::to_string(sdsl::size_in_bytes(int_vec_cat_1)) + "\n"s;

			dense_compr_int_vec_cat_1 = dense_compr_int_vector<uint32_t>(int_vec_cat_1);

			cat_1_merge.CloseAndRemove();
		});
	});

	merging_threads.emplace_back([this]{
		tc.execute([this]{
			//this is more complicated in streaming mode since there are still duplicates possible, but I think I should be able to create a little larger int_vector and at the end resize it
			
			auto start_time = std::chrono::high_resolution_clock::now();

			sdsl::int_vector<> int_vec_cat_2(this->temp_data_stats.cat_2_elems, 0, cfg.bits_per_elem_in_cat_2);
			using cat_2_rec_t = std::tuple<int64_t, uint32_t, uint32_t>;
			archive_stream_merge<cat_2_rec_t> cat_2_merge(archive_file_name + "_cat2.bin");
			cat_2_rec_t cat2_rec;
			size_t i = 0;
			int64_t prev_idx = -1;
			//first rec because we also need to deduplicate
			if (cat_2_merge.NextRec(cat2_rec))
			{
				prev_idx = std::get<0>(cat2_rec);
				auto file_id = std::get<1>(cat2_rec);
				auto cnt = std::get<2>(cat2_rec);
				if (cnt > cfg.max_cnt) //I think this is never true in the current implementation
					cnt = cfg.max_cnt;
				int_vec_cat_2[i++] = (file_id << cfg.bits_for_counter) + cnt;
			}

			//remaining records, now we can check prev
			while (cat_2_merge.NextRec(cat2_rec))
			{
				//check if that was the same k-mer, and if so skip adding it second time
				int idx = std::get<0>(cat2_rec);
				if(idx == prev_idx)
					continue;
				prev_idx = idx;

				auto file_id = std::get<1>(cat2_rec);
				auto cnt = std::get<2>(cat2_rec);
				if (cnt > cfg.max_cnt) //I think this is never true in the current implementation
					cnt = cfg.max_cnt;
				int_vec_cat_2[i++] = (file_id << cfg.bits_for_counter) + cnt;
			}

			//now we may resize

			//std::cerr << "temp_data_stats.cat_2_elems global before deduplication resize: "s + std::to_string(this->temp_data_stats.cat_2_elems) + "\n"s;

			int_vec_cat_2.resize(i);
			this->temp_data_stats.cat_2_elems = i;

			//std::cerr << "int_vec_cat_2.size() after resize: "s + std::to_string(int_vec_cat_2.size()) + "\n"s;

			dense_compr_int_vec_cat_2 = dense_compr_int_vector<uint32_t>(int_vec_cat_2);

			auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
			std::cerr << "int_vec_cat_2 created. Time: "s + std::to_string(time) + ", size: " + std::to_string(sdsl::size_in_bytes(int_vec_cat_2)) + "\n"s;

			cat_2_merge.CloseAndRemove();
		});
	});
	
	merging_threads.emplace_back([this]{
		tc.execute([this]{
			
			auto start_time = std::chrono::high_resolution_clock::now();

			sdsl::int_vector<> int_vec_cat_3(this->temp_data_stats.cat_3_elems, 0, cfg.bits_per_elem_in_cat_3);
			//+1 for guard
			bv_cat3_delim = sdsl::bit_vector(this->temp_data_stats.cat_3_elems + 1);

			int64_t cur_idx = -1;

			using cat_3_rec_t = std::pair<int64_t, uint32_t>;
			archive_stream_merge<cat_3_rec_t> cat_3_merge(archive_file_name + "_cat3.bin");
			cat_3_rec_t cat3_rec;
			size_t i = 0;

			while (cat_3_merge.NextRec(cat3_rec))
			{
				assert(cat3_rec.first != -1);

				int_vec_cat_3[i++] = cat3_rec.second;
					if (cat3_rec.first != cur_idx) {
					bv_cat3_delim[i-1] = 1;
					cur_idx = cat3_rec.first;
				}
			}

			//guard
			bv_cat3_delim[this->temp_data_stats.cat_3_elems] = 1;

			dense_compr_int_vec_cat_3 = dense_compr_int_vector<uint32_t>(int_vec_cat_3);

			auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
			std::cerr << "int_vec_cat_3 created. Time: "s + std::to_string(time) + ", size: " + std::to_string(sdsl::size_in_bytes(int_vec_cat_3)) + "\n"s;

			start_time = std::chrono::high_resolution_clock::now();

			cat_3_merge.CloseAndRemove();
		});
	});
	
	for (auto& th: merging_threads)
		th.join();

	auto time =std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
	std::cerr << "Total merge time: " << time << "\n";
}

void KmerIndex_SBWTBuilder::Serialize(std::ostream& out, LookupIndex& lookup_index) const
{
	std::cerr << "Serialize SBWT\n";
	lookup_index.sbwt.start = out.tellp();
	sbwt::serialize_string(sbwt.variant, out);
	sbwt.sbwt.serialize(out);
	lookup_index.sbwt.end = out.tellp();

	std::cerr << "Done. Size: " << lookup_index.sbwt.size() << " bytes\n";

	std::cerr << "Serialize config which is:\n";
	cfg.Print();
	lookup_index.cfg.start = out.tellp();
	cfg.Serialize(out);
	lookup_index.cfg.end = out.tellp();
	std::cerr << "Done. Size: " << lookup_index.cfg.size() << " bytes\n";

	auto serialize_report_size = [&out](const auto& to_serialize, const std::string& name, LookupIndex::elem& index_elem)
	{
		std::cerr << "Serialize " << name << "\n";
		index_elem.start = out.tellp();
		to_serialize.serialize(out);
		index_elem.end = out.tellp();
		std::cerr << "Done. Size: " << index_elem.size() << " bytes\n";
	};

#define SERIALIZE_WRAPPER(TO_SERIALIZE) serialize_report_size(TO_SERIALIZE, #TO_SERIALIZE, lookup_index.TO_SERIALIZE)

	SERIALIZE_WRAPPER(bv_cat1);
	SERIALIZE_WRAPPER(bv_cat2);
	SERIALIZE_WRAPPER(bv_cat3);
	SERIALIZE_WRAPPER(dense_compr_int_vec_cat_1);
	SERIALIZE_WRAPPER(dense_compr_int_vec_cat_2);
	SERIALIZE_WRAPPER(dense_compr_int_vec_cat_3);
	SERIALIZE_WRAPPER(bv_cat3_delim);

#undef SERIALIZE_WRAPPER
}
