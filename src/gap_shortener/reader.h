#pragma once
#include <vector>
#include <zlib-ng/zlib.h>
#include "defs.h"

class Reader
{
	std::vector<read_t> reads;

	uint64_t current_reads_bytes{}, current_header_bytes{};

	refresh::parallel_queue<read_pack_t>& reads_queue;

	read_t to_read_t(const std::string& str, bool& hasN);

	std::string currentLine;
	read_t last_read;
	std::string last_read_header;

	bool is_fastq;

	uint64_t total_bytes{};
	uint64_t total_bases{};
	uint64_t total_symb_header{};

	void addRead();
	
	void porcessFastaOrMultiFasta(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile);
	void processFastq(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile);
public:
	explicit Reader(const std::string& path, refresh::parallel_queue<read_pack_t>& reads_queue);
	void GetStats(uint64_t& total_bytes, uint64_t& total_bases, uint64_t& total_symb_header)
	{
		total_bytes = this->total_bytes;
		total_bases = this->total_bases;
		total_symb_header = this->total_symb_header;
	}
};