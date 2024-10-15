#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <iostream>
#include <thread>
#include <thread>
#include <atomic>

#include "../common/accepted_anchors.h"
#include "../common/version.h"

struct Params {

	inline const static std::string default_output_ext = ".filtered.satc";

	std::vector<std::string> input;
	std::vector<std::string> output;

	std::string anchor_list;

	uint32_t n = 0; // keep n most freq targets per each anchor
	uint32_t t = (std::min)(8u, std::thread::hardware_concurrency());
	uint32_t L = 0;
	bool l = false;

	void Print(std::ostream& oss) const
	{
		oss << "input         :\n";
		for (const auto& x : input)
			oss << "\t" << x << "\n";
		oss << "output        :\n";
		for (const auto& x : output)
			oss << "\t" << x << "\n";
		oss << "anchor_list   : " << anchor_list << "\n";
		oss << "n             : " << n << "\n";
		oss << "t             : " << t << "\n";
		oss << "L             : " << L << "\n";
		oss << "l             : " << std::boolalpha << l << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_filter\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options]\n";

		Params defaults;
		std::cerr
			<< "Options:\n"
			<< "    -i <path> [REQUIRED] - path to input SATC file (or txt with list of paths, if -l used)" << "\n"
			<< "    -o <path>            - path to output SATC file (or txt with list of paths, if -l used) (default: input path with added extension " << default_output_ext << ")" << "\n"
			<< "    -d <path> [REQUIRED] - path to file with anchors dictionary (k-mer per line)" << "\n"
			<< "    -n <int>             - number of most freq targets to keep for each anchor in each sample (0 means keep all), default: " << defaults.n << "\n"
			<< "    -t <int>             - number of threads, default: " << defaults.t << "\n"
			<< "    -L <int>             - keep only samples with number of anchors > L, default: " << defaults.L << "\n"
			<< "    -l                   - if set <input> and <output> are treated as lists of file paths\n";
	}
};


Params read_params(int argc, char** argv)
{
	Params res;
	if (argc == 1) {
		Params::Usage(argv[0]);
		exit(0);
	}
	int i = 1;
	for (; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];

		if (param == "-i" && i + 1 < argc)
			res.input.emplace_back(argv[++i]);
		else if (param == "-o" && i + 1 < argc)
			res.output.emplace_back(argv[++i]);
		else if (param == "-d" && i + 1 < argc)
			res.anchor_list = argv[++i];
		else if (param == "-n" && i + 1 < argc)
			res.n = std::stoul(argv[++i]);
		else if (param == "-t" && i + 1 < argc)
			res.t = std::stoul(argv[++i]);
		else if (param == "-L" && i + 1 < argc)
			res.L = std::stoul(argv[++i]);
		else if (param == "-l")
			res.l = true;
		else {
			std::cerr << "Error: unknown option " << param << "\n";
			exit(1);
		}
	}

	if (res.t == 0) {
		res.t = Params{}.t; //get default
	}

	if (res.anchor_list.empty()) {
		std::cerr << "Error: anchor_list missing (-d)\n";
		exit(1);
	}

	if (res.input.empty()) {
		std::cerr << "Error: input missing (-i)\n";
		exit(1);
	}

	auto read_file_list = [](const std::string& path, std::vector<std::string>& out) {
		std::ifstream in(path);
		if (!in)
		{
			std::cerr << "Error: cannot open file " << path << " for reading\n";
			exit(1);
		}
		out.clear(); //path may be the the first element in out, so clear must be after reading from path...
		std::string line;
		while (std::getline(in, line))
		{
			if (line.empty())
				continue;
			out.emplace_back(std::move(line));
		}
	};

	//input is a file with list of paths
	if (res.l) {
		read_file_list(res.input.front(), res.input);

		//if output was defined it should be list of output files
		if (!res.output.empty())
			read_file_list(res.output.front(), res.output);
	}

	//if output is empty by default just add extension to input files
	if (res.output.empty()) {
		for (const auto& in : res.input)
			res.output.emplace_back(in + Params::default_output_ext);
	}

	if (res.input.size() != res.output.size()) {
		std::cerr << "Error: number of input and output paths is not the same\n";
		exit(1);
	}

	return res;
}

struct SampleAnchorTargetCount {
	uint32_t sample_id;
	uint64_t anchor;
	uint64_t target;
	uint32_t count;
	SampleAnchorTargetCount(
		uint32_t sample_id,
		uint64_t anchor,
		uint64_t target,
		uint32_t count) :
		sample_id(sample_id),
		anchor(anchor),
		target(target),
		count(count)
	{}
	bool operator<(const SampleAnchorTargetCount& rhs) const {
		//count from largest
		return std::make_tuple(sample_id, anchor, rhs.count, target) <
			std::make_tuple(rhs.sample_id, rhs.anchor, count, rhs.target);
	}
};
//for 10X
struct SampleBarcodeAnchorTargetCount {
	uint32_t sample_id;
	uint64_t barcode;
	uint64_t anchor;
	uint64_t target;
	uint32_t count;
	SampleBarcodeAnchorTargetCount(
		uint32_t sample_id,
		uint64_t barcode,
		uint64_t anchor,
		uint64_t target,
		uint32_t count) :
		sample_id(sample_id),
		barcode(barcode),
		anchor(anchor),
		target(target),
		count(count)
	{}
	bool operator<(const SampleBarcodeAnchorTargetCount& rhs) const {
		//count from largest
		return std::make_tuple(sample_id, barcode, anchor, rhs.count, target) <
			std::make_tuple(rhs.sample_id, rhs.barcode, rhs.anchor, count, rhs.target);
	}
};

//mkokot_TODO: this code is very similar to the process_single function, consider refactoring...
void process_single_10X(const AcceptedAnchors& accepted_anchors,
	const Header& header,
	buffered_binary_reader& reader,
	const std::string& output,
	uint64_t n_most_freq_targets,
	uint64_t anchors_per_sample_threshold) {

	Record rec;
	std::vector<SampleBarcodeAnchorTargetCount> data;
	while (rec.load(reader, header)) {
		if (accepted_anchors.IsAccepted(rec.anchor))
			data.emplace_back(rec.sample_id, rec.barcode, rec.anchor, rec.target, rec.count);
	}

	std::sort(data.begin(), data.end());

	buffered_binary_writer out(output);

	if (!out) {
		std::cerr << "Error: cannot open output file " << output << "\n";
		exit(1);
	}
	Header out_header;
	out_header = header;
	out_header.ordering = Header::ordering_t::SBATC;

	out_header.serialize(out);
	Record out_rec;

	if (n_most_freq_targets == 0)
		n_most_freq_targets = std::numeric_limits<decltype(n_most_freq_targets)>::max();

	for (uint64_t i = 0; i < data.size();) {

		out_rec.sample_id = data[i].sample_id;
		out_rec.barcode = data[i].barcode;

		auto anchors_start_idx = i;
		uint64_t num_anchors = 0;

		while (i < data.size() && data[i].sample_id == out_rec.sample_id && data[i].barcode == out_rec.barcode) {
			++num_anchors;
			auto anchor = data[i].anchor;

			while (i < data.size() && data[i].sample_id == out_rec.sample_id && data[i].barcode == out_rec.barcode && anchor == data[i].anchor)
				++i;
		}

		if (num_anchors > anchors_per_sample_threshold) {
			auto anchor_end_idx = i;

			for (uint64_t j = anchors_start_idx; j < anchor_end_idx; ) {
				out_rec.anchor = data[j].anchor;
				uint64_t n_targets = 1;
				while (j < anchor_end_idx && out_rec.anchor == data[j].anchor) {
					if (n_targets <= n_most_freq_targets) {
						out_rec.target = data[j].target;
						out_rec.count = data[j].count;
						out_rec.serialize(out, out_header);
					}

					++n_targets;
					++j;
				}
			}
		}
	}
}

void process_single(
	const AcceptedAnchors& accepted_anchors,
	const std::string& input,
	const std::string& output,
	uint64_t n_most_freq_targets,
	uint64_t anchors_per_sample_threshold)
{

	std::cerr << "processing " + input + " -> " + output + "\n";

	Header header;
	buffered_binary_reader reader(input);
	header.load(reader);

	if (header.barcode_len_symbols != 0) {
		process_single_10X(accepted_anchors, header, reader, output, n_most_freq_targets, anchors_per_sample_threshold);
		return;
	}

	Record rec;
	std::vector<SampleAnchorTargetCount> data;
	while (rec.load(reader, header)) {
		if (accepted_anchors.IsAccepted(rec.anchor))
			data.emplace_back(rec.sample_id, rec.anchor, rec.target, rec.count);
	}

	std::sort(data.begin(), data.end());

	buffered_binary_writer out(output);

	if (!out) {
		std::cerr << "Error: cannot open output file " << output << "\n";
		exit(1);
	}
	Header out_header;
	out_header = header;
	out_header.ordering = Header::ordering_t::SBATC;

	out_header.serialize(out);
	Record out_rec;

	if (n_most_freq_targets == 0)
		n_most_freq_targets = std::numeric_limits<decltype(n_most_freq_targets)>::max();

	for (uint64_t i = 0; i < data.size();) {

		out_rec.sample_id = data[i].sample_id;

		auto anchors_start_idx = i;
		uint64_t num_anchors = 0;

		while (i < data.size() && data[i].sample_id == out_rec.sample_id) {
			++num_anchors;
			auto anchor = data[i].anchor;

			while (i < data.size() && data[i].sample_id == out_rec.sample_id && anchor == data[i].anchor)
				++i;
		}

		if (num_anchors > anchors_per_sample_threshold) {
			auto anchor_end_idx = i;

			for (uint64_t j = anchors_start_idx; j < anchor_end_idx; ) {
				out_rec.anchor = data[j].anchor;
				uint64_t n_targets = 1;
				while (j < anchor_end_idx && out_rec.anchor == data[j].anchor) {
					if (n_targets <= n_most_freq_targets) {
						out_rec.target = data[j].target;
						out_rec.count = data[j].count;
						out_rec.serialize(out, out_header);
					}

					++n_targets;
					++j;
				}
			}
		}

	}
}

void process(const Params& params) {

	const AcceptedAnchors accepted_anchors(params.anchor_list);

	std::vector<std::thread> threads;
	auto n_threads = params.t;

	if (params.input.size() < n_threads)
		n_threads = params.input.size();

	assert(n_threads > 0);
	threads.reserve(n_threads - 1); //main thread will also work

	std::atomic<size_t> task_id{};

	auto task = [&]
	{
		while (true)
		{
			const auto my_task_id = task_id++;
			if (my_task_id >= params.input.size())
				break;
			process_single(accepted_anchors, params.input[my_task_id], params.output[my_task_id], params.n, params.L);
		}
	};

	for (uint64_t tid = 0; tid < n_threads - 1; ++tid)
		threads.emplace_back(task);
	task();

	for (auto& th : threads)
		th.join();
}

int main(int argc, char** argv) {
	auto params = read_params(argc, argv);
	params.Print(std::cerr);

	process(params);

	return 0;
}
