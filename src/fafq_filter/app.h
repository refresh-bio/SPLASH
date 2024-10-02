#pragma once

#include <array>
#include <algorithm>
#include <memory>
#include "params.h"
#include "../common/accepted_anchors.h"

class CApplication
{
	array<uint8_t, 256> char2bits;
	
	shared_ptr<AcceptedAnchors> accepted_anchors;

	CParams params;

	void init()
	{
		fill(char2bits.begin(), char2bits.end(), 4);
		char2bits['A'] = char2bits['a'] = 0;
		char2bits['C'] = char2bits['c'] = 1;
		char2bits['G'] = char2bits['g'] = 2;
		char2bits['T'] = char2bits['t'] = 3;
	}

public:
	CApplication()
	{
		init();
	}

	bool parse_args(int argc, char** argv);
	void usage();
	bool load_strings(vector<string>& vec, const string& fn);
	bool prepare_anchor_dict(const vector<string>& vec);
	bool run();
};