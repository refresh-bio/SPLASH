#ifndef _ARTIFACTS_FILTER_H
#define _ARTIFACTS_FILTER_H
#include "../../common/types/satc_data.h"
#include <unordered_set>
#include <map>

#include <refresh/hash_tables/lib/hash_set.h>
#include <refresh/hash_tables/lib/bloom_set.h>
#include <refresh/hash_tables/lib/murmur_hash.h>

class ArtifactsFilter
{
	std::vector<refresh::hash_set_lp<uint64_t, std::equal_to<uint64_t>, refresh::MurMur64Hash>> artifacts; //index is artifact len
	uint32_t min_art_len = 1000;
	refresh::bloom_set<uint64_t, refresh::MurMur64Hash, 2> bf_artifacts;

	void rebuild_bloom()
	{
		size_t req_size = 0;
		for (auto& x : artifacts)
			req_size += x.size();

		bf_artifacts.resize(req_size, 0.01);

		for (auto& x : artifacts)
			for (auto kmer : x)
				bf_artifacts.insert(kmer);
	}

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
			if (len >= artifacts.size())
				artifacts.resize(len + 1, refresh::hash_set_lp<uint64_t, std::equal_to<uint64_t>, refresh::MurMur64Hash>(~0ull, 16, 0.5));

			if (len < min_art_len)
				min_art_len = len;

			artifacts[len].insert(str_kmer_to_uint64_t(artifact));
		}

		rebuild_bloom();
	}
	ArtifactsFilter() :
		ArtifactsFilter("")
	{

	}

	void Add(uint32_t len, const std::vector<uint64_t>& new_artifacts) {
		if (len >= artifacts.size())
			artifacts.resize(len + 1, refresh::hash_set_lp<uint64_t, std::equal_to<uint64_t>, refresh::MurMur64Hash>(~0ull, 16, 0.5));

		if (len < min_art_len)
			min_art_len = len;

		auto& art = artifacts[len];

		for (auto x : new_artifacts)
			art.insert(x);

		rebuild_bloom();
	}

	bool ContainsArtifact(uint64_t anchor, uint32_t len) const
	{
		for(uint32_t art_len = min_art_len; art_len < (uint32_t) artifacts.size(); ++art_len)
		{
			//skip too long artifacts
			if (art_len > len)
//				continue;
				break;

			const auto& art = artifacts[art_len];
			
			uint64_t mask = ((1ull) << (2 * art_len)) - 1;

			auto no_parts = len - art_len + 1;
			uint64_t check = anchor;
			for (uint32_t i = 0; i < no_parts - 1; ++i)
			{
				uint64_t key = check & mask;
				size_t h_val = refresh::MurMur64Hash{}(key);
				if (bf_artifacts.check(key, h_val) && art.check(key, h_val))
					return true;
				check >>= 2;
			}

			//last one
//			if (art.check(check & mask))
			uint64_t key = check & mask;
			size_t h_val = refresh::MurMur64Hash{}(key);
			if (bf_artifacts.check(key, h_val) && art.check(key, h_val))
				return true;
		}
		return false;
	}
};


#endif 
