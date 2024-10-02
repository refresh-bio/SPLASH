#ifndef QUERY_COMMON_H
#define QUERY_COMMON_H
#include <string>
#include <vector>
#include <thread>
#include <cinttypes>
#include <algorithm>

namespace query_common {
	enum class InputFmt { FASTA, EXTENDORS, COMPACTORS };

	enum class ReportFmt { PLAIN, CONCISE, IDS, EMPTY };

	enum class StatsFmt { EMPTY, WITH_STATS };

	enum class OutputFmt { TXT, EXTENDORS, COMPACTORS };

	std::string to_string(InputFmt input_fmt);
	InputFmt input_format_from_string(const std::string& fmt);

	std::string to_string(ReportFmt report_fmt);
	ReportFmt report_format_from_string(const std::string& fmt);

	std::string to_string(StatsFmt stats_format);
	StatsFmt stats_format_from_string(const std::string& stats_format);

	std::string to_string(OutputFmt output_format);
	OutputFmt output_format_from_string(const std::string& output_format);

	struct SingleInputConfig {
		std::string input;
		std::string output;

		InputFmt input_fmt = InputFmt::FASTA;
		ReportFmt report_fmt = ReportFmt::CONCISE;
		StatsFmt stats_fmt = StatsFmt::EMPTY;
		OutputFmt output_fmt = OutputFmt::TXT;
		uint32_t kmer_skip = 0;

		void Print(std::ostream& oss) const;
		SingleInputConfig(const std::vector<std::string>& desc);
	};

	struct QueryCfg {
		std::vector<SingleInputConfig> input_cfg;
		std::string lookup;

		uint32_t n_threads = std::min(8u, std::thread::hardware_concurrency());

		void Print(std::ostream& oss) const;
	};

	void run(const QueryCfg& query_cfg);
}

#endif // !QUERY_COMMON_H

