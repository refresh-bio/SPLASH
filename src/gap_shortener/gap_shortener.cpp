#include <iostream>
#include <thread>
#include <string>
#include "../common/version.h"
#include "reader.h"
#include "pseudoreads_maker.h"
#include <algorithm>

struct Params
{
	std::string input_path;
	std::string output_base;
	uint64_t anchor_len{};
	uint64_t gap_len{};
	uint64_t new_gap_len = std::numeric_limits<uint64_t>::max();
	uint64_t target_len{};
	uint64_t n_output_files = 1;
	uint64_t n_threads = 2;
	int compression_level = 3;

	void Print(std::ostream& oss) const
	{
		oss << "input path                     : " << input_path << "\n";
		oss << "output base                    : " << output_base << "\n";
		oss << "anchor len                     : " << anchor_len << "\n";
		oss << "gap len                        : " << gap_len << "\n";
		oss << "new gap len                    : " << new_gap_len << "\n";
		oss << "target len                     : " << target_len << "\n";
		oss << "n output files                 : " << n_output_files << "\n";
		oss << "n threads                      : " << n_threads << "\n";
		oss << "compression level              : " << compression_level << "\n";
	}
	static void Usage(char* prog_name) {
		std::cerr << "gap_shortener - converts input file such that (a+g+t)-mers are the same in original file like (a+new_gap +t)-mers in output of this program\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <input_path> <output_base>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input_path>        - path input fastq/fasta file\n"
			<< "    <output_base>       - base name of output (will be extended with \"_<i>.fa.gz\" where i is in [1, n_output_files]\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_len <int>        - anchor len\n"
			<< "    --gap_len <int>           - gap len\n"
			<< "    --new_gap_len <int>       - new gap len\n"
			<< "    --target_len <int>        - target len\n"
			<< "    --n_output_files <int>    - the number of output files (default: 1)\n"
			<< "    --n_threads <int>         - the number of computing threads (at least two will be created - one for reading and one for shortening gaps) (default: 2)\n"
			<< "    --compression_level <int> - output compression level (default: 3)\n";
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

		if (param == "--n_output_files") {
			std::string tmp = argv[++i];
			res.n_output_files = std::stoull(tmp);
		}
		else if (param == "--anchor_len") {
			std::string tmp = argv[++i];
			res.anchor_len = std::stoull(tmp);
		}
		else if (param == "--gap_len") {
			std::string tmp = argv[++i];
			res.gap_len = std::stoull(tmp);
		}
		else if (param == "--new_gap_len") {
			std::string tmp = argv[++i];
			res.new_gap_len = std::stoull(tmp);
		}
		else if (param == "--target_len") {
			std::string tmp = argv[++i];
			res.target_len = std::stoull(tmp);
		}
		else if (param == "--n_threads") {
			std::string tmp = argv[++i];
			res.n_threads = std::stoull(tmp);
		}
		else if (param == "--compression_level") {
			std::string tmp = argv[++i];
			res.compression_level = std::stoi(tmp);
		}
	}
	if (i >= argc) {
		std::cerr << "Error: input_file missing\n";
		exit(1);
	}

	res.input_path = argv[i++];
	if (i >= argc) {
		std::cerr << "Error: output_base missing\n";
		exit(1);
	}
	res.output_base = argv[i++];


	if (res.anchor_len == 0) {
		std::cerr << "Error: anchor len (--anchor_len) must be specified\n";
		exit(1);
	}
	if (res.target_len == 0) {
		std::cerr << "Error: target len (--target_len) must be specified\n";
		exit(1);
	}
	if (res.gap_len == 0) {
		std::cerr << "Error: gap len (--gap_len) must be specified\n";
		exit(1);
	}
	if (res.new_gap_len == std::numeric_limits<uint64_t>::max()) {
		std::cerr << "Error: gap len (--new_gap_len) must be specified\n";
		exit(1);
	}
	if (res.new_gap_len >= res.gap_len) {
		std::cerr << "Error: new gap len must be lower than gap len\n";
		exit(1);
	}

	return res;
}

int main(int argc, char** argv)
{
	std::cerr << "Welcome to gap_shortener\n";

	auto params = read_params(argc, argv);

	params.Print(std::cerr);

	
	size_t n_pseudoreads_makers = std::max(params.n_threads - 1, (uint64_t)1);

	refresh::parallel_queue<read_pack_t> reads_q(n_pseudoreads_makers * 20, 1, "reads_q");

	refresh::parallel_queue<pseudoreads_pack_t> pseudoreads_pack_q(n_pseudoreads_makers * 2, n_pseudoreads_makers, "pseudoreads_pack_q");

	std::vector<std::thread> threads;

	threads.emplace_back([&] {
		Reader reader(params.input_path, reads_q);
	});

	for (size_t i = 0; i < n_pseudoreads_makers; ++i) {
		threads.emplace_back([&] {
			PseudoreadsMaker pseudoreads_maker(reads_q, pseudoreads_pack_q, params.anchor_len, params.gap_len, params.new_gap_len, params.target_len, params.compression_level);
		});
	}

	threads.emplace_back([&] {
		pseudoreads_pack_t pack;

		std::vector<FILE*> files;
		std::vector<size_t> writen_per_file(params.n_output_files);
		std::vector<std::string> file_names;
		for (size_t i = 1; i <= params.n_output_files; ++i) {
			std::string fname = params.output_base + "_" + std::to_string(i) + ".fa.gz";
			file_names.push_back(fname);
			files.emplace_back(fopen(fname.c_str(), "wb"));
			if (!files.back()) {
				std::cerr << "Error: cannot open output file: " << fname << "\n";
				exit(1);
			}
			setvbuf(files.back(), nullptr, _IOFBF, 16 << 20);
		}
		size_t file_id{};
		while (pseudoreads_pack_q.pop(pack)) {
			auto writen = fwrite(pack.data(), 1, pack.size(), files[file_id]);
			if (writen != pack.size()) {
				std::cerr << "Error: cannot write " << pack.size() << " bytes to output file " << file_names[file_id] << " (only " << writen << " bytes writen)\n";
				exit(1);
			}
			writen_per_file[file_id] += writen;
			file_id = (file_id + 1) % files.size();
		}
		for (size_t i = 0; i < files.size(); ++i) {
			fclose(files[i]);
			//if there were not data writen to a file it is not a correct gz file, lest just create correct gz file
			if (!writen_per_file[i]) {
				std::cerr << "Warning: there is no data writen to " << file_names[i] << ". Creating empty gz file\n";
				gzFile gz_file = gzopen(file_names[i].c_str(), "wb");
				if (!gz_file) {
					std::cerr << "Error: cannot create empty gz file\n";
					exit(1);
				}
				gzclose(gz_file);
			}
		}
	});

	for (auto& x : threads)
		x.join();

}