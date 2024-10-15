#include "lookup.h"
#include <chrono>

const std::string Lookup::MARKER = "KMER_LOOKUP";
const uint32_t Lookup::VERSION_MAJOR = 2;
const uint32_t Lookup::VERSION_MINOR = 2;

const uint32_t Lookup::MIN_SUPPORTED_VERSION_MAJOR = 2;
const uint32_t Lookup::MIN_SUPPORTED_VERSION_MINOR = 2;

//OK so we have version number of index which is the same as version of lookup table that produced index (VERSION_MAJOR.VERSION_MINOR)
//although is not exactly the software version because it is the same as splash version
//we have minimum required version than the current version of lookup table may read (MIN_SUPPORTED_VERSION_MAJOR.MIN_SUPPORTED_VERSION_MINOR)
//we will never support newer version but gently exit

std::ifstream Lookup::open_ifstream(const std::string& path)
{
	std::ifstream in(path, std::ios::binary);
	if (!in) {
		//maybe this is previous version and we need to append ".ht" extension
		in.open(path + ".ht", std::ios::binary);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
	}
	return in;
}

kmer_index_variant_t Lookup::read_kmer_index_variant_impl(std::istream& in)
{
	uint32_t ver_maj;
	uint32_t ver_min;
	std::string marker;
	bool makred_readed = read_string_verify_len(in, marker, MARKER.length());
	if (!makred_readed || marker != MARKER) {
		std::cerr << "Error: wrong marker at the beginning of index file\n";
		exit(1);
	}
	read_little_endian(in, ver_maj);
	read_little_endian(in, ver_min);

	auto ver_to_str = [](const auto& ver) {
		return std::to_string(ver.first) + "." + std::to_string(ver.second);
		};

	const auto soft_version = std::make_pair(VERSION_MAJOR, VERSION_MINOR);
	const auto index_version = std::make_pair(ver_maj, ver_min);

	if (index_version > soft_version)
	{
		std::cerr << "Error: you are using lookup table ver. " << ver_to_str(soft_version) << ", but trying to read newer index ver. " << ver_to_str(index_version) << ". Update your lookup_table version\n";
		exit(1);
	}

	const auto min_supported_index_version = std::make_pair(MIN_SUPPORTED_VERSION_MAJOR, MIN_SUPPORTED_VERSION_MINOR);

	//now lets check if the index is not too old
	if (index_version < min_supported_index_version)
	{
		std::cerr << "Error: you are trying to read lookup table index version " << ver_to_str(index_version)
			<< ", but this version of lookup_table (" << ver_to_str(soft_version)
			<< ") requires index version at least " << ver_to_str(min_supported_index_version) << "\n";
		exit(1);
	}
	kmer_index_variant_t kmer_index_variant;

	std::string index_variant;
	read_string(in, index_variant);
	kmer_index_variant = kmer_index_variant_from_string(index_variant);

	return kmer_index_variant;
}
kmer_index_variant_t Lookup::ReadKmerIndexVariant(const std::string& path)
{
	auto in = open_ifstream(path);
	return read_kmer_index_variant_impl(in);
}

Lookup::Lookup(const std::string& path, uint32_t n_threads) {
	auto in = open_ifstream(path);
	kmer_index_variant_t kmer_index_variant = read_kmer_index_variant_impl(in);

	uint64_t poly_ACGT_len;
	read_little_endian(in, poly_ACGT_len);
	std::cerr << "poly_ACGT_len: " << poly_ACGT_len << "\n";
	poly_ACGT_filter = PolyACGTFilter(poly_ACGT_len);

	LookupIndex lookup_index;

	lookup_index.Load(in);

	std::vector<std::function<void()>> load_tasks;
	load_tasks.emplace_back([&] {
		std::cerr << "Loading file mapper...\n";
		auto in = open_ifstream(path);
		in.seekg(lookup_index.file_mapper.start);
		file_mapper.Load(in);
		assert((size_t)in.tellg() == lookup_index.file_mapper.end);
		std::cerr << "Loading file mapper done\n";
	});

	load_tasks.emplace_back([&] {
		std::cerr << "Loading header mapper...\n";
		auto in = open_ifstream(path);
		in.seekg(lookup_index.header_mapper.start);
		header_mapper.Load(in);
		assert((size_t)in.tellg() == lookup_index.header_mapper.end);
		std::cerr << "Loading header mapper done\n";
	});

	switch (kmer_index_variant)
	{
	case kmer_index_variant_t::sbwt:
		kmer_index = std::make_unique<KmerIndex_SBWT>(poly_ACGT_filter, header_mapper, path, lookup_index, load_tasks);
		break;
	default:
		std::cerr << "Error: not implemented (" << __FILE__ << ":" << __LINE__ << ")";
		exit(1);
	}

	std::atomic<int> task_id{};
	std::vector<std::thread> threads(n_threads);
	assert(n_threads > 0);
	for (size_t i = 0 ; i < n_threads ; ++i) {
		threads[i] = std::thread([&] {
			while(true) {
				auto id = task_id++;
				if(id >= (int)load_tasks.size())
					break;
				load_tasks[id]();
			}
		});
	}

	//wait for threads to end
	for (auto& th: threads)
		th.join();

	if (auto ptr = reinterpret_cast<KmerIndex_SBWT*>(kmer_index.get()))
		ptr->print_config();

	//go to end marker pos
	in.seekg(-(sizeof(size_t) + MARKER.size()), ios_base::end);
	std::string marker;
	read_string(in, marker);
	if (marker != MARKER) {
		std::cerr << "Error: wrong marker at the end of index file\n";
		exit(1);
	}
}

void Lookup::Serialize(std::ostream& out,
	uint64_t poly_ACGT_len,
	const SeqIdMapper& file_mapper,
	const HeaderIdMapper& header_mapper,
	kmer_index_variant_t kmer_index_variant,
	const std::function<void(std::ostream&, LookupIndex&)>& serialize_kmer_index)
{
	write_string_with_len(out, MARKER);
	write_little_endian(out, VERSION_MAJOR);
	write_little_endian(out, VERSION_MINOR);

	write_string_with_len(out, to_string(kmer_index_variant));

	write_little_endian(out, poly_ACGT_len);

	LookupIndex lookup_index;

	lookup_index.index.start = out.tellp();
	lookup_index.Serialize(out); //serialize empty to make space
	lookup_index.index.end = out.tellp();

	lookup_index.file_mapper.start = out.tellp();
	file_mapper.Serialize(out);
	lookup_index.file_mapper.end = out.tellp();

	lookup_index.header_mapper.start = out.tellp();
	header_mapper.Serialize(out);
	lookup_index.header_mapper.end = out.tellp();

	serialize_kmer_index(out, lookup_index);

	auto tmp = out.tellp();
	out.seekp(lookup_index.index.start);
	lookup_index.Serialize(out); //serialize again, now its filled
	assert(lookup_index.index.end == (size_t)out.tellp());
	out.seekp(tmp);

	write_string_with_len(out, MARKER);
}

void Lookup::DumpMapping(const std::string& path) {
	std::ofstream out(path);
	if (!out) {
		std::cerr << "Error: cannot open output file " << path << "\n";
		exit(1);
	}

	out << "---- MAPPING OF FILES ----\n";
	file_mapper.Dump(out);
	out << " ---- MAPPING OF HAEDERS ----\n";
	header_mapper.Dump(out, file_mapper);
}

void Lookup::TruncatePaths() {
	file_mapper.TransformNames([](std::string& name) {
		auto p = name.find_last_of("/\\");
		if (p != std::string::npos)
			name = name.substr(p + 1);
	});
}

void Lookup::query_seq(const std::string& seq, uint32_t kmer_skip, QueryResult& query_result)
{
	assert(kmer_index.get());
	kmer_index->query(seq, kmer_skip, query_result);
}

