#ifndef _ARTIFACTS_FILTER_H
#define _ARTIFACTS_FILTER_H
#include "../common/satc_data.h"
#include <unordered_set>
#include <map>

class ArtifactsFilter
{
	std::map<uint32_t, std::unordered_set<uint64_t>> artifacts; //index is artifact len
public:
	ArtifactsFilter(const std::string& path) //empty path means no filtering
	{
		if (path == "")
			return;

		std::ifstream in(path);
		
		if (!in)
		{
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		std::string artifact;

		while (in >> artifact)
		{
			auto len = artifact.length();
			assert(len <= 32);
			artifacts[len].insert(str_kmer_to_uint64_t(artifact));
		}
	}
	ArtifactsFilter() :
		ArtifactsFilter("")
	{

	}

	void Add(uint32_t len, const std::vector<uint64_t>& new_artifacts) {
		auto& art = artifacts[len];
		for (auto x : new_artifacts)
			art.insert(x);
	}

	bool ContainsArtifact(uint64_t anchor, uint32_t len) const
	{
		for (const auto& x : artifacts)
		{
			uint32_t art_len = x.first;

			//skip too long artifacts
			if (art_len > len)
				continue;

			const auto& art = x.second;
			
			uint64_t mask = ((1ull) << (2 * art_len)) - 1;

			auto no_parts = len - art_len + 1;
			uint64_t check = anchor;
			for (uint32_t i = 0; i < no_parts - 1; ++i)
			{
				if (art.count(check & mask))
					return true;
				check >>= 2;
			}
			//last one
			if (art.count(check & mask))
				return true;
		}
		return false;
	}
};


#endif 
