#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include "../common/types/common_types.h"
#include "engine.h"

class Console {

	const std::string PARAM_INPUT_FORMAT {"--input_format"};

	const std::string PARAM_NUM_KMERS {"--num_kmers"};
	const std::string PARAM_KMER_LEN{ "--kmer_len" };
	const std::string SWITCH_ALL_ANCHORS{ "--all_anchors"};
	
	const std::string PARAM_EPSILON { "--epsilon" };
	const std::string PARAM_BETA{ "--beta" };
	const std::string PARAM_LOWER_BOUND{ "--lower_bound" };
	const std::string PARAM_MAX_MISMATCH{ "--max_mismatch" };

	const std::string PARAM_OUT_FASTA{ "--out_fasta" };
	const std::string PARAM_POLY_THRSEHOLD{ "--poly_threshold" };

	const std::string SWITCH_EDIT_DISTANCE{ "--use_edit_distance" };
	const std::string SWITCH_NO_EXTENSION{ "--no_extension" };
	const std::string PARAM_MAX_LENGTH{ "--max_length" };
	const std::string PARAM_MIN_EXTENDER_SPECIFICITY{ "--min_extender_specificity" };
	const std::string PARAM_NUM_EXTENDERS{ "--num_extenders" };
	const std::string PARAM_EXTENDERS_SHIFT{ "--extenders_shift"};
	
	const std::string PARAM_MAX_ANCHOR_COMPACTORS{ "--max_anchor_compactors" };
	const std::string PARAM_MAX_CHILD_COMPACTORS{ "--max_child_compactors" };
	const std::string SWITCH_EXTEND_ALL{ "--extend_all" };
	const std::string SWITCH_NEW_ACCEPTANCE_RULE{ "--new_acceptance_rule" };
	
	const std::string PARAM_NUM_THREADS{ "--num_threads" };
	const std::string SWITCH_VERBOSE{ "-v" };
	const std::string PARAM_READS_BUFFER{ "--reads_buffer_gb" };
	const std::string PARAM_ANCHORS_BATCH{ "--anchors_batch" };
	const std::string SWITCH_KEEP_TEMP{ "--keep_temp" };
	const std::string SWITCH_NO_SUBCOMPACTORS{ "--no_subcompactors"};
	const std::string SWITCH_CUMULATED_STATS{ "--cumulated_stats" };
	const std::string SWITCH_INDEPENDENT_OUTPUTS{ "--independent_outputs" };

	const std::string PARAM_LOG{ "--log" };

	Engine::Params params;

public:

	bool verbose{ false };
	int numThreads{ 4 };
	int readsBufferGb{ 24 };
	size_t anchorsBatchSize{ 500000 };

	std::string anchorsTsv;
	std::vector<std::string> sampleFastqs;
	input_format_t inputFormat = input_format_t::fastq;
	std::string outputFasta;
	std::string outputTsv;
	bool keepTemp{ false };
	bool noSubcompactors{ false };
	bool cumulatedStats{ false };
	bool independentOutputs{ false };

	std::string logFile;

	const Engine::Params& getParams() const { return params; }

	bool parse(int argc, char** argv);

	void printUsage();

	bool findSwitch(std::vector<std::string>& params, const std::string& name) {
		auto it = find(params.begin(), params.end(), name); // verbose mode
		if (it != params.end()) {
			params.erase(it);
			return true;
		}

		return false;
	}

	template <typename T>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& v) {
		auto prevToEnd = std::prev(params.end());
		auto it = find(params.begin(), prevToEnd, name); // verbose mode
		if (it != prevToEnd) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				params.erase(it, it + 2);
				return true;
			}
		}

		return false;
	}



};
