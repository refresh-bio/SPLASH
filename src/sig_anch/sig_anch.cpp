// --fdr_threshold 0.05 --samplesheet test_samplesheet.csv --outfile_scores output/anchors_pvals.tsv --outfile_all_anchors_pvals output/all_anchors_pvals.tsv --outfile_Cjs output/anchors_Cjs_random_opt.tsv --infile_bins input_bins.txt

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <regex>
#include <memory>
#include <iterator>

#include "../common/version.h"
#include "../common/csv.h"
#include "../common/types/base_coding.h"

//#define ENABLE_CJ_CODE
//#define ENABLE_READ_ALL_MODE

using namespace std;

double fdr_threshold = 0.05;
string samplesheet;
string infile_bins;
string outfile_scores;
string outfile_all_anchors_pvals;
#ifdef ENABLE_CJ_CODE
string outfile_Cjs;
#endif
string col_name;

uint32_t anchor_len = 0;

vector<string> bin_files;

vector<string> csv_header_in;
vector<string> csv_header_out;

BaseCoding4 base_coding;

struct rec_data_t
{
	uint64_t anchor;
	uint64_t file_id;
	double col_to_correct;
#ifdef ENABLE_CJ_CODE
	double pval_samplesheet;
	double pval_aggregated;
#endif
	enum class field_id {col_to_correct, pval_samplesheet, pval_aggregated};

#ifdef ENABLE_CJ_CODE
	rec_data_t(uint64_t anchor, uint64_t file_id, double col_to_correct, double pval_samplesheet = 0, double pval_aggregated = 0) :
		anchor(anchor), file_id(file_id), col_to_correct(col_to_correct), pval_samplesheet(pval_samplesheet), pval_aggregated(pval_aggregated)
#else
	rec_data_t(uint64_t anchor, uint64_t file_id, double col_to_correct) :
		anchor(anchor), file_id(file_id), col_to_correct(col_to_correct)
#endif
	{}

	bool operator<(const rec_data_t& x)
	{
		return anchor < x.anchor;
	}

#ifdef ENABLE_CJ_CODE
	double get(field_id id) const
	{
		if (id == field_id::col_to_correct)
			return col_to_correct;
		else if (id == field_id::pval_samplesheet)
			return pval_samplesheet;
		else
			return pval_aggregated;
	}

	void set(field_id id, double val) 
	{
		if (id == field_id::col_to_correct)
			col_to_correct = val;
		else if (id == field_id::pval_samplesheet)
			pval_samplesheet = val;
		else
			pval_aggregated = val;
	}
#endif
};

vector<rec_data_t> csv_records;

int col_to_correct_id = -1;
#ifdef ENABLE_CJ_CODE
int pval_samplesheet_id = -1;
int pval_aggregated_id = -1;
#endif
int col_to_correct_corrected_id = -1;
#ifdef ENABLE_CJ_CODE
int pval_samplesheet_corrected_id = -1;
int pval_aggregated_corrected_id = -1;
#endif

#ifdef ENABLE_READ_ALL_MODE
refresh::csv_file csv_data("\t", '\t');
#endif

#ifdef ENABLE_CJ_CODE
bool use_sheet_cjs;
#endif

bool parse_params(int argc, char** argv);
bool usage();

#ifdef ENABLE_READ_ALL_MODE
bool load_tsv(const string& fn);
vector<pair<double, bool>> fdr_correction(const vector<double>& pvals, double alpha);
vector<double> fdr_correction_simple(const vector<double>& pvals, double alpha);
void correct_pvals();
#else
bool load_tsv_stream(const string& fn, uint64_t file_id);
#ifdef ENABLE_CJ_CODE
void fdr_correction_simple_stream(vector<rec_data_t>& records, rec_data_t::field_id field, double alpha);
#else
void fdr_correction_simple_stream(vector<rec_data_t>& records, double alpha);
#endif
void fdr_correction_stream(vector<double>& pvals, double alpha);
void correct_pvals_stream();
#endif

// ************************************************************************************
bool usage()
{
	cerr << "sig_anch\n";
	SPLASH_VER_PRINT(cerr);
	cerr << "Usage: sig_anch [options]" << endl
		<< "Options:" << endl
		<< "  --fdr_threshold <value>" << endl
#ifdef ENABLE_CJ_CODE
		<< "  --samplesheet <path>" << endl
#endif
		<< "  --outfile_scores <path>" << endl
		<< "  --infile_bins <path>" << endl
		<< "  --outfile_all_anchors_pvals <path>" << endl
		<< "  --outfile_Cjs <path>" << endl
		<< "  --col_name <string>" << endl;

	return false;
}

// ************************************************************************************
bool parse_params(int argc, char** argv)
{
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == "--fdr_threshold"s && i + 1 < argc)
			fdr_threshold = atof(argv[++i]);
		else if (argv[i] == "--samplesheet"s && i + 1 < argc)
			samplesheet = argv[++i];
		else if (argv[i] == "--infile_bins"s && i + 1 < argc)
			infile_bins = argv[++i];
		else if (argv[i] == "--outfile_scores"s && i + 1 < argc)
			outfile_scores = argv[++i];
		else if (argv[i] == "--outfile_all_anchors_pvals"s && i + 1 < argc)
			outfile_all_anchors_pvals = argv[++i];
#ifdef ENABLE_CJ_CODE
		else if (argv[i] == "--outfile_Cjs"s && i + 1 < argc)
			outfile_Cjs = argv[++i];
#endif
		else if (argv[i] == "--col_name"s && i + 1 < argc)
			col_name = argv[++i];
	}

	if (!(fdr_threshold <= 1 && fdr_threshold >= 0 &&
		!samplesheet.empty() &&
		!infile_bins.empty() &&
		!outfile_scores.empty() &&
		!col_name.empty() &&
		!outfile_all_anchors_pvals.empty()))
#ifdef ENABLE_CJ_CODE
		//		!outfile_Cjs.empty()))
#endif
		return usage();

	ifstream ifs(infile_bins);

	if (ifs.bad())
	{
		cerr << "Cannot open file: " << infile_bins << endl;
		return false;
	}

	bin_files.assign(istream_iterator<string>(ifs), istream_iterator<string>());

	if (bin_files.empty())
	{
		cerr << "No input files in: " << infile_bins << endl;
		return false;
	}

	return true;
}

#ifdef ENABLE_READ_ALL_MODE
// ************************************************************************************
bool load_tsv(const string& fn)
{
	refresh::csv_file cf("\t", '\t');

	if (!cf.load(fn))
		return false;

	return csv_data.merge(cf);
}
#else
// ************************************************************************************
bool load_tsv_stream(const string& fn, uint64_t file_id)
{
	refresh::csv_istream cfs(fn, "\t");

	if (!cfs.is_open())
		return false;

	if (csv_header_in.empty())
	{
		csv_header_in = cfs.get_header();

		col_to_correct_id = cfs.col_id(col_name);
		if (col_to_correct_id < 0)
		{
			cerr << "No column " << col_name << " in " << fn << endl;
			return false;
		}

#ifdef ENABLE_CJ_CODE
		if (use_sheet_cjs)
		{
			pval_samplesheet_id = cfs.col_id("pval_samplesheet");
			pval_aggregated_id = cfs.col_id("pval_aggregated");

			if (pval_samplesheet_id < 0)
			{
				cerr << "No column pval_samplesheet_id in " << fn << endl;
				return false;
			}
			if (pval_aggregated_id < 0)
			{
				cerr << "No column pval_aggregated_id in " << fn << endl;
				return false;
			}
		}
#endif
	}
	else
		if (csv_header_in != cfs.get_header())
			return false;

	vector<string> record; 
	
	while (!cfs.eof())
	{
		if (cfs.get_record(record))
		{
			if (!anchor_len)
				anchor_len = record[0].size();

#ifdef ENABLE_CJ_CODE
			if (use_sheet_cjs)
				csv_records.emplace_back(base_coding.encode_bases_2b(record[0]), file_id, cfs.to_double(record[col_to_correct_id]), cfs.to_double(record[pval_samplesheet_id]), cfs.to_double(record[pval_aggregated_id]));
			else
#endif
				csv_records.emplace_back(base_coding.encode_bases_2b(record[0]), file_id, cfs.to_double(record[col_to_correct_id]));
		}
		else
		{
			if (cfs.eof())
				break;
			else
			{
				cerr << "Problem with parsing " << fn << " file" << endl;
				return false;
			}
		}
	}

	return true;
}
#endif

#ifdef ENABLE_READ_ALL_MODE
// ************************************************************************************
vector<pair<double, bool>> fdr_correction(const vector<double>& pvals, double alpha)
{
	vector<pair<double, int>> pv_sorted;
	vector<pair<double, bool>> pv_corr;

	int size = (int)pvals.size();

	pv_corr.resize(size);

	if (size == 0)
		return pv_corr;

	pv_sorted.resize(size);

	for (int i = 0; i < size; ++i)
		pv_sorted[i] = make_pair(pvals[i], i);

	stable_sort(pv_sorted.begin(), pv_sorted.end());

	double cm = 0;
	for (int i = 1; i <= size; ++i)
		cm += 1.0 / i;

	//    cout << "cm: " << cm << endl;

		// Calculate pvals_corrected_raw
	for (int i = 0; i < size; ++i)
	{
		pv_corr[i].first = pv_sorted[i].first * (size * cm) / (i + 1);
		pv_corr[i].second = pv_corr[i].first <= alpha;
	}

	/*    cout << "pvals_corrected_raw: ";
		for (auto x : pv_corr)
			cout << "[" << x.first << "," << x.second << "]  ";
		cout << endl;*/

		// Calculate pvals_corrected
	double min_pv = 1.0;
	for (int i = size - 1; i >= 0; --i)
	{
		min_pv = min(min_pv, pv_corr[i].first);
		pv_corr[i].first = min_pv;
	}

	/*    cout << "pvals_corrected: ";
		for (auto x : pv_corr)
			cout << "[" << x.first << "," << x.second << "]  ";
		cout << endl;*/

	vector<pair<double, bool>> pv_original_order;

	pv_original_order.resize(size);

	for (int i = 0; i < size; ++i)
		pv_original_order[pv_sorted[i].second] = pv_corr[i];

	return pv_original_order;
}

// ************************************************************************************
vector<double> fdr_correction_simple(const vector<double>& pvals, double alpha)
{
	auto col = fdr_correction(pvals, alpha);

	vector<double> res;
	res.reserve(col.size());

	for (const auto& x : col)
		res.emplace_back(x.first);

	return res;
}

// ************************************************************************************
void correct_pvals()
{
	if (csv_data.empty())
	{
		csv_data.save(outfile_scores);
		return;
	}

	auto &out_pvals = csv_data;

	if (csv_data.col_id(col_name) < 0)
	{
		cerr << "Column " << col_name << " does not exist in the input file\n";
		exit(0);
	}

	auto col_to_correct = csv_data.copy_col_double(col_name);
	auto col_corrected = fdr_correction_simple(col_to_correct, fdr_threshold);
	out_pvals.insert_col(col_name + "_corrected", col_corrected);

	vector<double> pval_samplesheet_corrected;
	vector<double> pval_aggregated_corrected;

#ifdef ENABLE_CJ_CODE
	if (use_sheet_cjs)
	{
		auto pval_samplesheet = csv_data.copy_col_double("pval_samplesheet");
		pval_samplesheet_corrected = fdr_correction_simple(pval_samplesheet, fdr_threshold);
		out_pvals.insert_col("pval_samplesheet_corrected", pval_samplesheet_corrected);

		auto pval_aggregated = csv_data.copy_col_double("pval_aggregated");
		pval_aggregated_corrected = fdr_correction_simple(pval_aggregated, fdr_threshold);
		out_pvals.insert_col("pval_aggregated_corrected", pval_aggregated_corrected);
	}
#endif

	auto nr = out_pvals.no_rows();

#ifdef ENABLE_CJ_CODE
	if (use_sheet_cjs)
	{
		for (size_t i = 0; i < nr; ++i)
			if (pval_aggregated_corrected[i] >= fdr_threshold)
				out_pvals.mark(i);
	}
	else
	{
#endif
		for (size_t i = 0; i < nr; ++i)
			if (col_corrected[i] >= fdr_threshold)
				out_pvals.mark(i);
#ifdef ENABLE_CJ_CODE
}
#endif

	out_pvals.sort();
	out_pvals.save(outfile_all_anchors_pvals);

	out_pvals.remove_marked();
//	out_pvals.remove_col("optHash"); //mkokot: we don't have this column

	auto header = out_pvals.header();

#ifdef ENABLE_CJ_CODE
	vector<size_t> sel_ids_cj;
	vector<size_t> sel_ids_stats;
	regex re_cj("anchor|cj_rand_opt_.*");
	regex re_stats("cj_rand_opt_.*");

	for (size_t i = 0; i < header.size(); ++i)
	{
		if (regex_match(header[i], re_cj))
			sel_ids_cj.emplace_back(i);
		if (regex_match(header[i], re_stats))
			sel_ids_stats.emplace_back(i);
	}

	auto csv_cj = out_pvals.filter(sel_ids_cj.begin(), sel_ids_cj.end());
	auto header_cj = csv_cj.header();

	regex re_repl("cj_rand_opt_");
	for (size_t i = 0; i < header_cj.size(); ++i)
	{
		string new_name = regex_replace(header_cj[i], re_repl, "");
		csv_cj.rename_col(i, new_name);
	}

	sort(sel_ids_stats.rbegin(), sel_ids_stats.rend());
	for (auto x : sel_ids_stats)
		out_pvals.remove_col(x);

#ifndef NO_SORT_CSV
	csv_cj.sort();
	out_pvals.sort();
#endif

	if(!outfile_Cjs.empty())
		csv_cj.save(outfile_Cjs);
#endif
	out_pvals.save(outfile_scores);
}
#else

// ************************************************************************************
void fdr_correction_stream(vector<double>& pvals, double alpha)
{
	vector<pair<double, int>> pv_sorted;

	int size = (int)pvals.size();

	if (size == 0)
		return;

	pv_sorted.resize(size);

	for (int i = 0; i < size; ++i)
		pv_sorted[i] = make_pair(pvals[i], i);

	stable_sort(pv_sorted.begin(), pv_sorted.end());

	double cm = 0;
	for (int i = 1; i <= size; ++i)
		cm += 1.0 / i;

	// Calculate pvals_corrected
	double min_pv = 1.0;
	for (int i = size - 1; i >= 0; --i)
	{
		double pv = pv_sorted[i].first * (size * cm) / (i + 1);
		min_pv = min(min_pv, pv);

		pvals[pv_sorted[i].second] = min_pv;
	}
}

// ************************************************************************************
#ifdef ENABLE_CJ_CODE
void fdr_correction_simple_stream(vector<rec_data_t>& records, rec_data_t::field_id field, double alpha)
#else
void fdr_correction_simple_stream(vector<rec_data_t>& records, double alpha)
#endif
{
	vector<double> to_correct;

	to_correct.reserve(records.size());

	for (const auto& rec : records)
#ifdef ENABLE_CJ_CODE
		to_correct.emplace_back(rec.get(field));
#else
		to_correct.emplace_back(rec.col_to_correct);
#endif

	fdr_correction_stream(to_correct, alpha);

	for (size_t i = 0; i < records.size(); ++i)
#ifdef ENABLE_CJ_CODE
		records[i].set(field, to_correct[i]);
#else
		records[i].col_to_correct = to_correct[i];
#endif
}

// ************************************************************************************
void correct_pvals_stream()
{
	cerr << "Correcting p-values\n";

	csv_records.shrink_to_fit();
	sort(csv_records.begin(), csv_records.end());

	refresh::csv_ostream ofs_scores(outfile_scores, '\t');
	refresh::csv_ostream ofs_all_pvals(outfile_all_anchors_pvals, '\t');
#ifdef ENABLE_CJ_CODE
	refresh::csv_ostream ofs_Cjs(outfile_Cjs, '\t');
#endif

	if (!ofs_scores.is_open())
	{
		cerr << "Cannot open " << outfile_scores << " file\n";
		return;
	}

	if (!ofs_all_pvals.is_open())
	{
		cerr << "Cannot open " << outfile_all_anchors_pvals << " file\n";
		return;
	}

#ifdef ENABLE_CJ_CODE
	if (!outfile_Cjs.empty() && !ofs_Cjs.is_open())
	{
		cerr << "Cannot open " << outfile_Cjs << " file\n";
		return;
	}
#endif

	csv_header_out = csv_header_in;
	csv_header_out.emplace_back(col_name + "_corrected");
	col_to_correct_corrected_id = (int) csv_header_out.size() - 1;

#ifdef ENABLE_CJ_CODE
	if (pval_samplesheet_id >= 0)
	{
		csv_header_out.emplace_back("pval_samplesheet_corrected");
		pval_samplesheet_corrected_id = (int) csv_header_out.size() - 1;
	}
	
	if (pval_aggregated_id >= 0)
	{
		csv_header_out.emplace_back("pval_aggregated_corrected");
		pval_aggregated_corrected_id = (int) csv_header_out.size() - 1;
	}
#endif

	ofs_scores.set_header(csv_header_out);
	ofs_all_pvals.set_header(csv_header_out);

#ifdef ENABLE_CJ_CODE
	vector<string> csv_header_cj;

	csv_header_cj.emplace_back(csv_header_out.front());

	vector<size_t> sel_ids_stats;
	if (use_sheet_cjs)
	{
		regex re_stats("cj_rand_opt_.*");
		regex re_repl("cj_rand_opt_");

		for (size_t i = 0; i < csv_header_out.size(); ++i)
		{
			if (regex_match(csv_header_out[i], re_stats))
			{
				sel_ids_stats.emplace_back(i);
				string new_name = regex_replace(csv_header_out[i], re_repl, "");
				csv_header_cj.emplace_back(new_name);
			}
		}

		ofs_Cjs.set_header(csv_header_cj);
	}
#endif

#ifdef ENABLE_CJ_CODE
	fdr_correction_simple_stream(csv_records, rec_data_t::field_id::col_to_correct, fdr_threshold);
#else
	fdr_correction_simple_stream(csv_records, fdr_threshold);
#endif

#ifdef ENABLE_CJ_CODE
	if (use_sheet_cjs)
	{
		fdr_correction_simple_stream(csv_records, rec_data_t::field_id::pval_samplesheet, fdr_threshold);
		fdr_correction_simple_stream(csv_records, rec_data_t::field_id::pval_aggregated, fdr_threshold);
	}
#endif

	vector<shared_ptr<refresh::csv_istream>> ifs_bins(bin_files.size(), nullptr);
	
	cerr << "Saving output files\n";

	for (size_t i = 0; i < bin_files.size(); ++i)
	{
		ifs_bins[i] = make_shared<refresh::csv_istream>(bin_files[i], '\t');
		if (!ifs_bins[i])
		{
			cerr << "Cannot open " << bin_files[i] << " file\n";
			return;
		}
	}

	vector<string> record;
	vector<string> record_cj;

	for (const auto& rec : csv_records)
	{
		ifs_bins[rec.file_id]->get_record(record);

		record.emplace_back(ofs_all_pvals.to_string(rec.col_to_correct));
#ifdef ENABLE_CJ_CODE
		if (pval_samplesheet_corrected_id >= 0)
			record.emplace_back(ofs_all_pvals.to_string(rec.pval_samplesheet));
		if (pval_aggregated_corrected_id >= 0)
			record.emplace_back(ofs_all_pvals.to_string(rec.pval_aggregated));
#endif

		ofs_all_pvals.add_record(record);

#ifdef ENABLE_CJ_CODE
		if ((use_sheet_cjs && rec.pval_aggregated < fdr_threshold) || (!use_sheet_cjs && rec.col_to_correct < fdr_threshold))
			ofs_scores.add_record(record);

		if (use_sheet_cjs)
		{
			record_cj.clear();
			for (auto id : sel_ids_stats)
				record_cj.emplace_back(record[id]);

			ofs_Cjs.add_record(record_cj);
		}
#else
		if (rec.col_to_correct < fdr_threshold)
			ofs_scores.add_record(record);
#endif
	}
}
#endif

// ************************************************************************************
int main(int argc, char **argv)
{
#ifdef _WIN32
	_setmaxstdio(2045);
#endif
	if (!parse_params(argc, argv))
		return 0;

	cerr << "Loading input files\n";

#ifdef ENABLE_CJ_CODE
	refresh::csv_file csv_sample_sheet(",", ',');
	if (!csv_sample_sheet.load(samplesheet))
	{
		cerr << "Cannot open " << samplesheet << endl;
		return 0;
	}

	use_sheet_cjs = csv_sample_sheet.no_cols() == 2;
#endif

#ifdef ENABLE_READ_ALL_MODE
	// Read-all mode

	for (auto& x : bin_files)
		if (!load_tsv(x))
			return 0;

	correct_pvals();
#else
	// Streaming mode
	for(size_t i = 0; i < bin_files.size(); ++i)
		if (!load_tsv_stream(bin_files[i], i))
			return 0;

	correct_pvals_stream();
#endif

	return 0;
}

