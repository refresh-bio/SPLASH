#include <iostream>
#include <unordered_set>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cstring>
#include <tuple>

using namespace std;

#include "../common/kmer.h"
#include "zlib.h"
#include "../common/version.h"
#include "anchor_read_collector.h"

enum class running_mode_t {per_input_file, anchor_split};
enum class read_mode_t {suffix, prefix, whole_read};

running_mode_t running_mode;
bool all_anchors;
read_mode_t read_mode = read_mode_t::suffix;
bool single_end_mode;

vector<tuple<string, string, string, string>> FQ_fn;
string K_fn;
string output_dir;

unordered_set<uint64_t> K_set;
uint32_t k_len;

CAnchorReads anchor_reads(1 << 25);

bool usage();
bool parse_args(int argc, char** argv);
bool load_K_file();
bool open_input_file(const string& fn_in, gzFile &gz_in);
bool open_output_file(const string& fn_out, FILE *&out);
void process_reads(const string& fn_in1, const string& fn_in2, const string& fn_out1, const string& fn_out2);

// ***************************************************************************
bool usage()
{
	cerr << "read-selector\n";
	SPLASH_VER_PRINT(cerr);
	cerr << "Usage:" << endl
		<< "read_selector [options] <k> <K_file_name> <PE_FASTQ_file_names|SE_FASTQ_file_names> [output_dir]"
		<< endl << endl
		<< "Options:\n"
		<< "   -prefix-mode        - report read prefixes (instead of suffixes by default) or reads\n"
		<< "   -whole-read-mode     - report whole read with marked anchor position\n"
		<< "   -all-anchors        - if provided all anchors in each read will be determined (only first by default)\n"
		<< "   -single-end-mode    - if provided input data are single-end"
		<< "Parameters:\n"
		<< "   k                   - anchor len\n"
		<< "   K_file_name         - file containing dictionary of anchors (in the first column of each line)\n"
		<< "   PE_FASTQ_file_names - description of tasks - one task per line (for default, paired-end mode)\n"
		<< "      Each task is described by 2 (only if output_dir is provided) or 4 file names (space-separated).\n"
		<< "      The first 2 of them are input file names. The (optional) remaining 2 are output file names.\n"
		<< "      Sample lines:\n"
		<< "         in/data1.fastq.gz in/data2.fastq.gz out/read1.fa out/read2.fa\n"
		<< "         in/data1.fastq.gz in/data2.fastq.gz\n"
		<< "   SE_FASTQ_file_names - description of tasks - one task per line (for single-end mode)\n"
		<< "      Each task is described by 1 (only if output_dir is provided) or 2 file names (space-separated).\n"
		<< "      The first of them are input file names. The (optional) remaining are output file names.\n"
		<< "      Sample lines:\n"
		<< "         in/data.fastq.gz out/read.fa\n"
		<< "         in/data1.fastq.gz\n"
		<< "   output_dir          - if provided the reads (only from 1st file of each task) will be placed in separate files according to anchor.\n";
		
	return false;
}

// ***************************************************************************
bool parse_args(int argc, char** argv)
{
	if (argc < 4)
		return usage();

	int i_arg;

	all_anchors = false;
	read_mode = read_mode_t::suffix;
	single_end_mode = false;

	for (i_arg = 1; i_arg < argc; ++i_arg)
	{
		if (argv[i_arg][0] != '-')
			break;

		if (argv[i_arg] == "-all-anchors"s)
			all_anchors = true;
		else if (argv[i_arg] == "-prefix-mode"s)
			read_mode = read_mode_t::prefix;
		else if (argv[i_arg] == "-whole-read-mode"s)
			read_mode = read_mode_t::whole_read;
		else if (argv[i_arg] == "-single-end-mode"s)
			single_end_mode = true;
	}

	if (i_arg + 3 > argc)
		return false;

	k_len = atoi(argv[i_arg++]);

	K_fn = argv[i_arg++];

	if (i_arg >= argc)
		return usage();

	ifstream ifs(argv[i_arg++]);

	if (ifs.bad())
	{
		cerr << "Cannot open file: " << argv[3] << endl;
		return false;
	}

	string fn_in1, fn_in2, fn_out1, fn_out2;
	string line;

	if (i_arg < argc)
	{
		output_dir = argv[i_arg];
		running_mode = running_mode_t::anchor_split;
	}
	else
		running_mode = running_mode_t::per_input_file;

	while (!ifs.eof())
	{
		getline(ifs, line);

		stringstream ss(line);

		if(single_end_mode)
			ss >> fn_in1 >> fn_out1;
		else
			ss >> fn_in1 >> fn_in2 >> fn_out1 >> fn_out2;

		if (line.empty())
			continue;

		if (single_end_mode)
		{
			if (running_mode == running_mode_t::per_input_file)
			{
				if (fn_in1.empty() || fn_out1.empty())
				{
					cerr << "Wrong specification in line: " << line << endl;
					cerr << "You should give 2 file names (1 for input and 1 for output), space-separated, in each line\n";
					continue;
				}
			}
			else
			{
				if (fn_in1.empty())
				{
					cerr << "Wrong specification in line: " << line << endl;
					cerr << "You should give 1 input file name in each line\n";
					continue;
				}
			}
		}
		else
		{
			if (running_mode == running_mode_t::per_input_file)
			{
				if (fn_in1.empty() || fn_in2.empty() || fn_out1.empty() || fn_out2.empty())
				{
					cerr << "Wrong specification in line: " << line << endl;
					cerr << "You should give 4 file names (2 for input and 2 for output), space-separated, in each line\n";
					continue;
				}
			}
			else
			{
				if (fn_in1.empty() || fn_in2.empty())
				{
					cerr << "Wrong specification in line: " << line << endl;
					cerr << "You should give 2 input file names, space-separated, in each line\n";
					continue;
				}
			}
		}

		if(single_end_mode)
			FQ_fn.emplace_back(fn_in1, "", fn_out1, "");
		else
			FQ_fn.emplace_back(fn_in1, fn_in2, fn_out1, fn_out2);
	}

	return true;
}

// ***************************************************************************
bool load_K_file()
{
	ifstream ifs(K_fn);

	if (ifs.bad())
	{
		cerr << "Cannot open " << K_fn << endl;
		return false;
	}

	CKmer kmer(k_len, kmer_mode_t::direct);

	string kmer_str;

	while (!ifs.eof())
	{
		getline(ifs, kmer_str);
			
		if (kmer_str.empty())
			continue;

		if (kmer_str.length() < k_len)
		{
			cerr << "Too short k-mer: " << kmer_str << endl;
			continue;
		}

		kmer.Reset();

		for (uint32_t i = 0; i < k_len; ++i)
		{
			auto c = dna_code(kmer_str[i]);

			if(c < 4)
				kmer.insert(c);
			else
			{
				cerr << "Wrong base in: " << kmer_str << endl;
				continue;
			}
		}

		if(kmer.is_full())
			K_set.insert(kmer.data());
	}

	return true;
}

// ***************************************************************************
bool open_input_file(const string& fn_in, gzFile& gz_in)
{
	gz_in = gzopen(fn_in.c_str(), "r");

	if (!gz_in)
	{
		cerr << "Cannot open file: " << fn_in << endl;
		return false;
	}

	gzbuffer(gz_in, 64 << 20);

	return true;
}

// ***************************************************************************
bool open_output_file(const string& fn_out, FILE*& out)
{
	out = fopen(fn_out.c_str(), "wb");

	if (!out)
	{
		cerr << "Cannot create file: " << fn_out << endl;
		return false;
	}

	setvbuf(out, nullptr, _IOFBF, 16 << 20);

	return true;
}

// ***************************************************************************
void process_reads(const string& fn_in1, const string& fn_in2, const string& fn_out1, const string& fn_out2)
{
	cerr << "Processing " << fn_in1 << endl;

	gzFile fq1 = nullptr, fq2 = nullptr;
	FILE* out1 = nullptr;
	FILE* out2 = nullptr;

	if (!open_input_file(fn_in1, fq1))
		return;
	if (!single_end_mode)
	{
		if (!open_input_file(fn_in2, fq2))
		{
			gzclose(fq1);
			return;
		}
	}

	if (running_mode == running_mode_t::per_input_file)
	{
		if (!open_output_file(fn_out1, out1))
		{
			gzclose(fq1);
			if(!single_end_mode)
				gzclose(fq2);
			return;
		}

		if (!single_end_mode)
		{
			if (!open_output_file(fn_out2, out2))
			{
				gzclose(fq1);
				gzclose(fq2);
				fclose(out1);
				return;
			}
		}
	}

	const int MAX_LEN = 1 << 20;

	char* tmp = new char[MAX_LEN];
	char* id1 = new char[MAX_LEN];
	char* id2 = new char[MAX_LEN];
	char* dna1 = new char[MAX_LEN];
	char* dna2 = new char[MAX_LEN];

	CKmer kmer(k_len, kmer_mode_t::direct);

	while (!gzeof(fq1) && !gzeof(fq2))
	{
		gzgets(fq1, id1, MAX_LEN);
		gzgets(fq1, dna1, MAX_LEN);
		gzgets(fq1, tmp, MAX_LEN);
		gzgets(fq1, tmp, MAX_LEN);

		if (!single_end_mode && running_mode == running_mode_t::per_input_file)
		{
			gzgets(fq2, id2, MAX_LEN);
			gzgets(fq2, dna2, MAX_LEN);
			gzgets(fq2, tmp, MAX_LEN);
			gzgets(fq2, tmp, MAX_LEN);
		}

		auto len1 = strlen(dna1);

		if (len1 < k_len)
			continue;

		kmer.Reset();

		uint32_t i;
		for (i = 0; i < len1; ++i)
		{
			auto c = dna_code(dna1[i]);

			if (c < 4)
				kmer.insert(c);
			else
				kmer.Reset();

			if(kmer.is_full())
				if (K_set.count(kmer.data()) > 0)
				{
					if(running_mode == running_mode_t::per_input_file)
					{
						::fputs(id1, out1);
						if (read_mode == read_mode_t::prefix)
						{
							fwrite(dna1, 1, i + 1, out1);
							putc('\n', out1);
						}
						else if(read_mode == read_mode_t::suffix)
							fputs(dna1 + i + 1 - k_len, out1);
						else
						{
							fprintf(out1, "%d\t", i + 1 - k_len);
							fputs(dna1, out1);
						}
//						fputs(dna1, out1);
//						fputs("+\n", out1);
//						fputs(dna1, out1);

						if (!single_end_mode)
						{
							::fputs(id2, out2);
							::fputs(dna2, out2);
						}
					}
					else
					{
						if (read_mode == read_mode_t::prefix)
							anchor_reads.Add(kmer.data_aligned(), string(dna1, dna1 + i + 1) + "\n"s, fn_in1);
						else if (read_mode == read_mode_t::suffix)
							anchor_reads.Add(kmer.data_aligned(), string(dna1 + i + 1 - k_len), fn_in1);
						else
							anchor_reads.Add(kmer.data_aligned(), to_string(i + 1 - k_len) + "\t"s + string(dna1), fn_in1);
					}

					if(!all_anchors)
						break;
				}
		}
	}

	delete[] tmp;
	delete[] id1;
	delete[] id2;
	delete[] dna1;
	delete[] dna2;

	gzclose(fq1);

	if(!single_end_mode)
		gzclose(fq2);

	if (running_mode == running_mode_t::per_input_file)
	{
		fclose(out1);
		if(!single_end_mode)
			fclose(out2);
	}
}

// ***************************************************************************
int main(int argc, char **argv)
{
	if (!parse_args(argc, argv))
		return 0;

	if (!load_K_file())
		return 0;

	if (running_mode == running_mode_t::anchor_split)
		anchor_reads.SetParams(k_len, output_dir);

	for (const auto& pe_files : FQ_fn)
		process_reads(get<0>(pe_files), get<1>(pe_files), get<2>(pe_files), get<3>(pe_files));

	return 0;
}

