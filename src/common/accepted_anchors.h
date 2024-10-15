#pragma once
#include "../common/types/satc_data.h"
#include <unordered_set>
#include <sstream>

#include <refresh/hash_tables/lib/hash_set.h>
#include <refresh/hash_tables/lib/bloom_set.h>
#include <refresh/hash_tables/lib/murmur_hash.h>

class ReadAnchorsFromPlainOrDSV
{
	uint32_t anchor_len = 0;
	char determine_delim(const std::string& header) {
		if (header.find('\t') != std::string::npos)
			return '\t';
		if (header.find(';') != std::string::npos)
			return ';';
		if (header.find(',') != std::string::npos)
			return ',';
		return '\t';
	}

	size_t find_anchor_col_id(const std::string& header, char delim, const std::string& path) {
		std::istringstream iss(header);
		size_t anchor_col_id = 0;
		std::string col_name;
		while (std::getline(iss, col_name, delim)) {
			if (col_name == "anchor")
				return anchor_col_id;
			++anchor_col_id;
		}
		std::cerr << "Error: cannot fine 'anchor' column in ile: " << path << "\n";
		exit(1);
	}
	void find_anchor(const std::string& line, char delim, std::string& anchor, size_t anchor_col_id, const std::string& path)
	{
		std::istringstream iss(line);
		size_t id = 0;
		while (std::getline(iss, anchor, delim)) {
			if (id == anchor_col_id)
				return;
			++id;
		}
		std::cerr << "Error: cannot read anchor from line " << line << " in a file: " << path << "\n";
		exit(1);
	}
public:
	template<typename ANCHOR_CALLBACK_T>
	ReadAnchorsFromPlainOrDSV(std::ifstream& in, const std::string& path, ANCHOR_CALLBACK_T callback)
	{
		std::string line;

		//empty file?
		if (!std::getline(in, line))
			return;

		//read from tsv file
		if (line.find("anchor") != std::string::npos) {

			char delim = determine_delim(line);
			size_t anchor_col_id = find_anchor_col_id(line, delim, path);
			std::string anchor;
			while (std::getline(in, line)) {
				find_anchor(line, delim, anchor, anchor_col_id, path);
				if (anchor_len == 0)
					anchor_len = anchor.length();
				else if (anchor_len != anchor.length()) {
					std::cerr << "Error: inconsistent anchor length in file " << path << "\n";
					exit(1);
				}

				callback(str_kmer_to_uint64_t(anchor));
			}
		}
		else { //plain text whitespace separated anchors

			//first just reaad all the anchors from first line
			std::istringstream iss(line);

			std::string anchor;
			while (iss >> anchor)
				callback(str_kmer_to_uint64_t(anchor));

			//and then the remaining from the file
			while (in >> anchor) {
				if (anchor_len == 0)
					anchor_len = anchor.length();
				else if (anchor_len != anchor.length()) {
					std::cerr << "Error: inconsistent anchor length in file " << path << "\n";
					exit(1);
				}

				callback(str_kmer_to_uint64_t(anchor));
			}
		}
	}
	uint32_t GetAnchorLen() const {
		return anchor_len;
	}
};

#if 0
class AcceptedAnchors {
	std::unordered_set<uint64_t> accepted_anchors;
	bool use_filter = false;

public:
	//if anchors empty all anchors are accepted
	AcceptedAnchors(const std::vector<uint64_t>& anchors) {
		if (anchors.size() == 0)
			return;

		use_filter = true;

		for (auto anchor : anchors)
			accepted_anchors.emplace(anchor);
	}
	//if path == "" all anchors accepted
	//the input path may be just a pure list of anchors as txt whitespace separated
	//or it may be a TSV file with header and column named 'anchor'
	AcceptedAnchors(const std::string& path) {
		if (path == "")
			return;

		use_filter = true;

		std::ifstream in(path);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		ReadAnchorsFromPlainOrDSV(in, path, [this](uint64_t anchor)
			{
				accepted_anchors.insert(anchor);
			});
	}
	bool IsAccepted(uint64_t anchor) const {
		if (!use_filter)
			return true;
		return accepted_anchors.find(anchor) != accepted_anchors.end();
	}
};
#else
class AcceptedAnchors {
	refresh::bloom_set<uint64_t, refresh::MurMur64Hash, 2> bf_accepted_anchors;
	refresh::hash_set_lp<uint64_t, std::equal_to<uint64_t>, refresh::MurMur64Hash> accepted_anchors;

	bool use_filter = false;

	void insert(const std::vector<uint64_t>& anchors)
	{
		bf_accepted_anchors.resize(anchors.size());
		accepted_anchors.reserve(anchors.size());

		for (const auto anchor : anchors)
		{
			bf_accepted_anchors.insert(anchor);
			accepted_anchors.insert(anchor);
		}
	}

public:
	//if anchors empty all anchors are accepted
	AcceptedAnchors(const std::vector<uint64_t>& anchors) : accepted_anchors(~0ull, 16, 0.5) {
		if (anchors.size() == 0)
			return;

		use_filter = true;

		insert(anchors);
	}
	//if path == "" all anchors accepted
	//the input path may be just a pure list of anchors as txt whitespace separated
	//or it may be a TSV file with header and column named 'anchor'
	AcceptedAnchors(const std::string& path) : accepted_anchors(~0ull, 16, 0.5) {
		if (path == "")
			return;

		use_filter = true;

		std::ifstream in(path);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		std::vector<uint64_t> anchors;

		ReadAnchorsFromPlainOrDSV(in, path, [this, &anchors](uint64_t anchor)
			{
				anchors.emplace_back(anchor);
			});

		insert(anchors);
	}
	bool IsAccepted(uint64_t anchor) const {
		if (!use_filter)
			return true;
		
		return bf_accepted_anchors.check(anchor) && accepted_anchors.check(anchor);
	}
};
#endif