#include "pseudoreads_maker.h"

pseudoreads_pack_t PseudoreadsMaker::compress() {
	//mkokot_TODO: make this field of class and use deflateReset (or something) to not allocate buffer each time
	
	stream.avail_in = (uInt)pseudoreads.size();
	stream.next_in = pseudoreads.data();

	size_t max_compressed_size = deflateBound(&stream, pseudoreads.size());

	pseudoreads_pack_t compressed(max_compressed_size);

	stream.avail_out = (uInt)compressed.capacity();
	stream.next_out = compressed.data();

	int status = deflate(&stream, Z_FINISH);
	if (status != Z_STREAM_END) {
		std::cerr << "Error: deflate error in " << __FILE__ << "(" << __LINE__ << ")\n";
		exit(1);
	}
	compressed.resize(stream.total_out);
	status = deflateReset(&stream);
	if (status != Z_OK) {
		std::cerr << "Error: deflateReset error in " << __FILE__ << "(" << __LINE__ << ")\n";
		exit(1);
	}
	return compressed;
}

void PseudoreadsMaker::make_pseudoreads(const read_t& read) {
	if (anchor_len + target_len + gap_len > read.length())
		return;
	size_t last_start_pos = read.length() - anchor_len - target_len - gap_len;
	size_t jump = new_gap_len + 1;

	for (size_t pos = 0; pos <= last_start_pos; pos += jump) {

		//+3 is ">\n" and "\n" after read
		if (pseudoreads.size() + max_out_read_len + 3 > pseudoreads.capacity())
			send();

		if (pseudoreads.size() + max_out_read_len + 3 > pseudoreads.capacity()) {
			std::cerr << "Error: internal buffer to small (" << __FILE__ << "(" << __LINE__ << ")" << ")\n";
			exit(1);
		}

		pseudoreads.push_back('>');
		pseudoreads.push_back('\n');
		std::copy(read.begin() + pos, read.begin() + pos + anchor_len + new_gap_len, pseudoreads.end());
		pseudoreads.resize(pseudoreads.size() + anchor_len + new_gap_len);
		
		size_t end_size = target_len + new_gap_len;
		if (pos + anchor_len + gap_len + end_size > read.size())
			end_size = read.size() - (pos + anchor_len + gap_len);

		std::copy(read.begin() + pos + anchor_len + gap_len, read.begin() + pos + anchor_len + gap_len + end_size, pseudoreads.end());
		pseudoreads.resize(pseudoreads.size() + end_size);
		pseudoreads.push_back('\n');
	}
}