#pragma once
#include <iostream>

enum class input_format_t { unknown, fasta, fastq, bam, cram };

inline input_format_t input_format_from_string(const std::string& str) {
	if (str == "fa" || str == "fasta" || str == "FASTA")
		return input_format_t::fasta;
	else if (str == "fq" || str == "fastq" || str == "FASTQ")
		return input_format_t::fastq;
	else
		return input_format_t::unknown;
}

inline std::string to_string(input_format_t input_format) {
	switch (input_format) {
		case input_format_t::fasta:
			return "fasta";
		case input_format_t::fastq:
			return "fastq";
		case input_format_t::bam:
			return "bam";
		case input_format_t::cram:
			return "cram";
		case input_format_t::unknown:
			return "unknown";
	}
}
