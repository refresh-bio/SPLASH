// -d ../Danny/all_aws_deduplicated_anchors_for_fasta.txt -i fl -o flo -l --anchor_len 27 --target_len 27 -v 1 -t 8 -n 10 
 
#include "app.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <array>
#include <cinttypes>
#include <thread>
#include <atomic>

#include "worker.h"
#include "../common/version.h"
// ******************************************************************************************
bool CApplication::load_strings(vector<string>& vec, const string& fn)
{
	ifstream ifs(fn);

	if (!ifs)
	{
		cerr << "Error: Cannot open: " << fn << endl;
		return false;
	}

	istream_iterator<string> ie;
	vec.assign(istream_iterator<string>(ifs), ie);

	return true;
}

// ******************************************************************************************
bool CApplication::prepare_anchor_dict(const vector<string>& vec)
{
	vector<kmer_t> anchors;
	anchors.reserve(vec.size());

	for (const auto& s : vec)
	{
		if (s.size() != params.anchor_len)
		{
			cerr << "Error: Wrong anchor length: " << s << endl;
			return false;
		}

		kmer_t anchor = 0;

		for (auto c : s)
		{
			kmer_t x = char2bits[(uint8_t) c];
			if (x > 3)
			{
				cerr << "Error: Anchor contains strange symbols: " << s << endl;
				return false;
			}

			anchor <<= 2;
			anchor += x;
		}

		anchors.emplace_back(anchor);
	}

	accepted_anchors = make_shared<AcceptedAnchors>(anchors);

	return true;
}

// ******************************************************************************************
bool CApplication::parse_args(int argc, char** argv)
{
	string dict_name;
	const string default_output_ext = ".filtered.satc";

	bool input_file_list = false;

	for (int i = 1; i < argc; ++i)
	{
		string par = argv[i];

		if (par == "-i" && i + 1 < argc)
		{
			params.input_fn.clear();
			params.input_fn.emplace_back(argv[++i]);
		}
		else if (par == "-o" && i + 1 < argc)
		{
			params.output_fn.clear();
			params.output_fn.emplace_back(argv[++i]);
		}
		else if (par == "--sample_id" && i + 1 < argc)
		{
			params.sample_ids.clear();
			params.sample_ids_fn.clear();
			params.sample_ids.emplace_back(atoi(argv[++i]));
		}
		else if (par == "--sample_ids" && i + 1 < argc)
		{
			params.sample_ids.clear();
			params.sample_ids_fn = argv[++i];
		}
		else if (par == "-l")
		{
			input_file_list = true;
		}
		else if (par == "-d" && i + 1 < argc)
		{
			dict_name = argv[++i];
		}
		else if ((par == "-n" && i + 1 < argc) || (par == "--top_target" && i + 1 < argc))
		{
			params.n_top_targets = atoi(argv[++i]);
			if (params.n_top_targets < 1 || params.n_top_targets > 1000)
			{
				cerr << "n must be in range [1, 1000]" << endl;
				return false;
			}
		}
		else if (par == "-t" && i + 1 < argc)
		{
			params.n_threads = atoi(argv[++i]);

			if (params.n_threads == 0)
			{
				params.n_threads = thread::hardware_concurrency();
				if (params.n_threads == 0)
					params.n_threads = 1;
			}
		}
		else if (par == "--anchor_len" && i + 1 < argc)
		{
			params.anchor_len = atoi(argv[++i]);
			if (params.anchor_len < 10 || params.anchor_len > 31)
			{
				cerr << "Error: anchor_len must be in range [10, 31]" << endl;
				return false;
			}
		}
		else if (par == "--target_len" && i + 1 < argc)
		{
			params.target_len = atoi(argv[++i]);
			if (params.target_len < 10 || params.target_len > 31)
			{
				cerr << "Error: target_len must be in range [10, 31]" << endl;
				return false;
			}
		}
		else if (par == "--allowed_headers" && i + 1 < argc)
		{
			params.allowed_headers = regex(argv[++i]);
			params.check_header = true;
		}
		else if (par == "--stats_json" && i + 1 < argc)
		{
			params.stats_json_fn = argv[++i];
		}
		else if (par == "-v" && i + 1 < argc)
		{
			params.verbosity = atoi(argv[++i]);
		}
		else
		{
			cerr << "Unknown option: " << par << endl;
			return false;
		}
	}

	if (params.anchor_len == 0)
	{
		cerr << "Error: anchor len must be provided" << endl;
		return false;
	}

	if (params.target_len == 0)
	{
		cerr << "Error: target len must be provided" << endl;
		return false;
	}

	// Prepare dictionary
	vector<string> ad;

	if (dict_name.empty())
	{
		cerr << "Error: dictionary of archors must be provided" << endl;
		return false;
	}

	if (!load_strings(ad, dict_name))
		return false;

	if (!prepare_anchor_dict(ad))
		return false;

	if (ad.empty())
	{
		cerr << "Error: Empty anchor dictionary" << endl;
		return false;
	}

	// Load file lists if provided
	if (input_file_list)
	{
		if (params.input_fn.empty())
		{
			cerr << "Error: A single file file FAST(A/Q) files must be provided with '-i' switch" << endl;
			return false;
		}

		string fl_fn = params.input_fn.front();
		params.input_fn.clear();

		if (!load_strings(params.input_fn, fl_fn))
			return false;

		if (params.output_fn.empty())
		{
			params.output_fn = params.input_fn;
			for_each(params.output_fn.begin(), params.output_fn.end(), [&default_output_ext](auto& s) {s += default_output_ext; });
		}
		else
		{
			fl_fn = params.output_fn.front();
			params.output_fn.clear();

			if (!load_strings(params.output_fn, fl_fn))
				return false;
		}

		if (params.input_fn.size() != params.output_fn.size())
		{
			cerr << "Error: Inconsistent no. of input and output file names" << endl;
			return false;
		}

		if (params.sample_ids_fn.empty())
		{
			params.sample_ids.resize(params.input_fn.size());
			iota(params.sample_ids.begin(), params.sample_ids.end(), 0);
		}
		else
		{
			params.sample_ids.clear();

			vector<string> s_ids;

			if (!load_strings(s_ids, params.sample_ids_fn))
				return false;

			params.sample_ids.resize(s_ids.size());
			transform(s_ids.begin(), s_ids.end(), params.sample_ids.begin(), [](const string& s) {return stoi(s); });
		}

		if (params.input_fn.size() != params.sample_ids.size())
		{
			cerr << "Error: Inconsistent no. of input file names and sample ids" << endl;
			return false;
		}
	}
	else
	{
		if (params.output_fn.empty())
		{
			if (!params.input_fn.empty())
				params.output_fn.emplace_back(params.input_fn.front() + default_output_ext);
			else
				params.output_fn.emplace_back("fafq" + default_output_ext);
		}

		if (params.sample_ids.empty())
			params.sample_ids.emplace_back(0);
	}
	
	return true;
}

// ******************************************************************************************
void CApplication::usage()
{
	cerr
		<< "Filter of anchor-target pairs from fasta/fastq files\n";
	SPLASH_VER_PRINT(cerr);
	cerr << "Usage:" << endl;
	cerr << "  fafq_filter [options]" << endl;
	cerr << "options:" << endl;
	cerr << "  -i <file_name>            - input FAST(A/Q) file or text file with FAST(A/Q) list (if not given, stdin will be used)" << endl;
	cerr << "  -o <file_name>            - output SATC file name or text file with output SATC files (optional)" << endl;
	cerr << "  -d <file_name>            - file with anchors dictionary (k-mer per line)" << endl;
	cerr << "  -l                        - treat input/output files as txt with FAST(A/Q) lists" << endl;
	cerr << "  -n <int>                  - no. of top targets to keep" << endl;
	cerr << "  --top_target <int>" << endl;
	cerr << "  --sample_id <int>         - sample id" << endl;
	cerr << "  --sample_ids <file_name>  - file with sample ids" << endl;
	cerr << "  --anchor_len <int>        - anchor length" << endl;
	cerr << "  --target_len <int>        - target length" << endl;
	cerr << "  --allowed_headers <regex> - regex for allowed headers (if not given, allow all)" << endl;
	cerr << "  --stats_json <file_name>  - name of output JSON file with stats" << endl;
	cerr << "  -t <int>                  - no. of threads (0 is for autodetect)" << endl;
	cerr << "  -v <int>                  - verbosity level (default: 0)" << endl;
}

// ******************************************************************************************
bool CApplication::run()
{
	atomic<size_t> job_id = 0;

	vector<thread> app_workers;
	app_workers.reserve(params.n_threads);

	if(params.input_fn.empty())
	{ 
		CWorker worker(params, accepted_anchors);
		if (!worker.process_stdin(params.sample_ids.front(), params.output_fn.front()))
			exit(1);
	}
	else
		for (int i = 0; i < params.n_threads; ++i)
		{
			app_workers.emplace_back([&] {
				CWorker worker(params, accepted_anchors);

				while (true)
				{
					size_t local_id = job_id.fetch_add(1);
					if (local_id >= params.input_fn.size())
						break;

					if (params.verbosity > 0)
						cerr << "Processing: " + params.input_fn[local_id] + "\n";

					if (!worker.process_file(params.sample_ids[local_id], params.input_fn[local_id], params.output_fn[local_id]))
						exit(1);
				}
				});
		}
	
	for (auto& t : app_workers)
		t.join();

/*	if (!params.sample_to_id_fn.empty())
	{
		ofstream ofs(params.sample_to_id_fn);

		if (!ofs)
		{
			cerr << "Cannot create sample to id mapping file: " << params.sample_to_id_fn << endl;
			return false;
		}

		for (size_t i = 0; i < params.input_fn.size(); ++i)
			ofs << params.input_fn[i] << "\t" << i << endl;
	}*/

	return true;
}
