#pragma once
#include "defs.h"
#include <iostream>
#include <zlib-ng/zlib.h>

class PseudoreadsMaker {
	refresh::parallel_queue<read_pack_t>& reads_queue;
	refresh::parallel_queue<pseudoreads_pack_t>& out_queue;
	size_t anchor_len;
	size_t gap_len;
	size_t new_gap_len;
	size_t target_len;

	size_t max_out_read_len;

	pseudoreads_pack_t pseudoreads;
	int compression_level;

	z_stream stream;

	pseudoreads_pack_t compress();
	void send() {
		out_queue.push(compress());
		pseudoreads.clear();
	}
	void make_pseudoreads(const read_t& read);
public:
	PseudoreadsMaker(refresh::parallel_queue<read_pack_t>& reads_queue, 
		refresh::parallel_queue<pseudoreads_pack_t>& out_queue,
		size_t anchor_len,
		size_t gap_len,
		size_t new_gap_len,
		size_t target_len,
		int compression_level) :

		reads_queue(reads_queue),
		out_queue(out_queue),
		anchor_len(anchor_len),
		gap_len(gap_len),
		new_gap_len(new_gap_len),
		target_len(target_len),
		pseudoreads(pseudoreads_pack_size),
		compression_level(compression_level)
	{
		max_out_read_len = anchor_len + 2 * new_gap_len + target_len;

		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
		
		if (deflateInit2(&stream, compression_level, Z_DEFLATED, 31, 9, Z_DEFAULT_STRATEGY) != Z_OK) {
			std::cerr << "Error: deflateInit2 failed in " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}

		read_pack_t read_pack;
		while (reads_queue.pop(read_pack))
			for (const auto& read : read_pack)
				make_pseudoreads(read);
		
		if (pseudoreads.size())
			send();

		int status = deflateEnd(&stream);
		if (status != Z_OK) {
			std::cerr << "Error: deflateEnd error in " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}

		out_queue.mark_completed();
	}
};
