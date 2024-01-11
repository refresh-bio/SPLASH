#pragma once
#include "../common/satc_data.h"
#include <unordered_set>
class AcceptedAnchors {
	std::unordered_set<uint64_t> accepted_anchors;
	bool use_filter = false;
public:
	AcceptedAnchors(const std::string& path) {
		if (path == "")
			return;
		use_filter = true;

		std::ifstream in(path);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		std::string anchor;
		while (in >> anchor)
			accepted_anchors.insert(str_kmer_to_uint64_t(anchor));

	}
	bool IsAccepted(uint64_t anchor) const {
		if (!use_filter)
			return true;
		return accepted_anchors.find(anchor) != accepted_anchors.end();
	}
};
