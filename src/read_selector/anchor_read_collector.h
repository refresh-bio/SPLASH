#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cinttypes>

using namespace std;

class CAnchorReads
{
	unordered_map<uint64_t, vector<pair<string, string>>> output_data;
	string output_dir;
	bool first_write = true;
	size_t max_reads_in_buffer;
	size_t no_reads_in_buffer;
	int k_len;

	string strip(string str)
	{
		while (!str.empty() && (str.back() == '\n' || str.back() == '\r'))
			str.pop_back();

		return str;
	}

	void SaveDict()
	{
		for (const auto& anchor_data : output_data)
		{
//			string fn = output_dir + anchor_data.second.front().first.substr(0, k_len) + ".intermediary";
			string fn = output_dir + kmer_to_string(anchor_data.first, k_len) + ".intermediary";
			FILE* out;

			out = fopen(fn.c_str(), first_write ? "wb" : "ab");
			if (!out)
			{
				cerr << "Cannot open output file: " << fn << endl;
				continue;
			}

			setvbuf(out, nullptr, _IOFBF, 16 << 20);

			for (const auto& read : anchor_data.second)
			{
				fputs(read.first.c_str(), out);
				putc('\t', out);
				fputs(read.second.c_str(), out);
				putc('\n', out);
			}

			fclose(out);
		}

		output_data.clear();
		no_reads_in_buffer = 0;

		first_write = false;
	}

	string kmer_to_string(uint64_t kmer, int len)
	{
		string str;

		for (int i = 0; i < len; ++i)
		{
			auto c = kmer & 3;
			str.push_back("ACGT"[c]);
			kmer >>= 2;
		}

		reverse(str.begin(), str.end());

		return str;
	}

public:
	CAnchorReads(const size_t _max_reads_in_buffer = 1 << 25) :
		max_reads_in_buffer(_max_reads_in_buffer),
		no_reads_in_buffer(0),
		k_len(0)
	{}
	~CAnchorReads()
	{
		SaveDict();
	}

	void SetParams(const int _k_len, const string& _output_dir)
	{
		k_len = _k_len;
		output_dir = _output_dir;
		first_write = true;
	}

	void Add(const uint64_t anchor, const string& read, const string &fq_name)
	{
		output_data[anchor].emplace_back(strip(read), strip(fq_name));

		if (++no_reads_in_buffer >= max_reads_in_buffer)
			SaveDict();
	}
};

