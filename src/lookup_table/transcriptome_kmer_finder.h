#ifndef _TRANSCRIPTOME_KMER_FINDER_H
#define _TRANSCRIPTOME_KMER_FINDER_H

#include <string>
#include <cinttypes>
#include <cassert>
#include <chrono>
#include <map>
#include "readers.h"
#include "walk_kmers.h"

namespace build_mode {
	class TranscriptomeKmerFinder {
		size_t kmer_len;
		std::string path;
		MultiFastaReader reader;
		
		//handles format like in GRCh38_latest_rna.fna
		std::string get_gene_name_parentheses(const std::string& header) {
			if (std::count(header.begin(), header.end(), '(') == 1) {
				auto x = header.find_last_of("(");
				assert(x != std::string::npos);
				auto y = header.find_last_of(")");
				assert(y != std::string::npos);
				++x;
				--y;
				return header.substr(x, y - x + 1);
			}
			auto modified = header;
			auto p = modified.find("transcript variant");
			if (p != std::string::npos) {
				modified = header.substr(0, p);
				modified = modified.substr(0, modified.find_last_of(','));
			}

			auto x = modified.find_last_of("(");
			if (x == std::string::npos) {
				std::cerr << "Error: cannot extract gene name from " << header << "\n";
				exit(1);
			}
			auto y = modified.find_last_of(")");
			if (y == std::string::npos) {
				std::cerr << "Error: cannot extract gene name from " << header << "\n";
				exit(1);
			}
			++x;
			--y;
			return modified.substr(x, y - x + 1);
		}

		//handles format like in CAT_liftoff_genes.fa
		std::string get_gene_name_dot_dash(const std::string& header) {

			//get whats before first whitespace
			std::string part = header.substr(0, header.find_first_of(" \t"));

			auto n_dots = std::count(part.begin(), part.end(), '.');

			if (n_dots == 1) {
				auto p = part.find_first_of('.');
				return part.substr(0, p);
			}
			else if (n_dots == 0) {
				auto p = part.find_first_of('-');
				if (p != std::string::npos)
					return part.substr(0, p);
				return part;
			}
			else {
				std::cerr << "Error: don't know how to extract gene name from: " << header << "\n";
				exit(1);
			}
		}

		//mkokot_TODO: make it more generic
		std::string get_gene_name(const std::string& header) {
			if (std::count(header.begin(), header.end(), '(') > 0) {
				return get_gene_name_parentheses(header);
			}
			return get_gene_name_dot_dash(header);
		}

		class GeneIdMaker
		{
			std::map<std::string, uint32_t> _m;

		public:
			template<typename REGISTER_LABEL_T>
			uint32_t get_gene_id(const std::string& gene, const REGISTER_LABEL_T& register_label)
			{
				auto it = _m.find(gene);
				if (it != _m.end())
					return it->second;
				return _m.emplace(gene, register_label(gene)).first->second;
			}

			auto size() const
			{
				return _m.size();
			}
		};

	public:
		TranscriptomeKmerFinder(size_t kmer_len, const std::string& path) : kmer_len(kmer_len), path(path), reader(path) {

		}

		//REPORT_UNIQUE_KMER_T must fulfill std::function<void(uint64_t /*kmer*/, uint32_t /*ids*/)> - this is if we have k-mer that is in only one header (meaning single variant of single gene) - id is variant_id in such a case or in a single gene but in multiple variants, id is gene_id in such a case
		//REPORT_KMER_T must fulfill std::function<void(uint64_t /*kmer*/)> - just tell this k-mer exists here, but it does not uniquely identify any gene
		//REGISTER_LABEL_T must fullfill std::function<uint32_t(const std::string& /*gene_name*/)>
		template<typename REPORT_UNIQUE_KMER_T, typename REPORT_KMER_T, typename REGISTER_LABEL_T>
		void Process(
			const REPORT_UNIQUE_KMER_T& report_unique_kmer,
			const REPORT_KMER_T& report_kmer,
			const REGISTER_LABEL_T& register_label) {

			struct kmer_gene_variant
			{
				uint64_t kmer;
				uint32_t gene_id;
				uint32_t variant_id;
				kmer_gene_variant(uint64_t kmer,uint32_t gene_id, uint32_t variant_id) :
					kmer(kmer),
					gene_id(gene_id),
					variant_id(variant_id)
				{

				}
			};

			GeneIdMaker gene_id_maker;

			size_t num_headers = 0;
			std::string header, seq;

			std::vector<kmer_gene_variant> data;

			while (reader.NextSeq(header, seq)) {
				++num_headers;
				std::string gene_name = get_gene_name(header);
				auto gene_id = gene_id_maker.get_gene_id(gene_name, register_label);
				auto variant_id = register_label(header);

				WalkKmers kmer_walker(seq, kmer_len, true);
				uint64_t kmer;
				while (kmer_walker.Next(kmer))
					data.emplace_back(kmer, gene_id, variant_id);
			}

			std::cerr << "Number of genes in " << path << ": " << gene_id_maker.size() << "\n";
			std::cerr << "Number of FASTA records in " << path << ": " << num_headers << "\n";
			std::cerr << "Number of k-mers: " << data.size() << "\n";

			std::cerr << "Sorting...";
			auto sort_start_time = std::chrono::high_resolution_clock::now();
			std::sort(data.begin(), data.end(), [](const auto& lhs, const auto& rhs) { return std::make_tuple(lhs.kmer, lhs.gene_id, lhs.variant_id) < std::make_tuple(rhs.kmer, rhs.gene_id, rhs.variant_id); });
			std::cerr << "\nDone\n";
			std::cerr << "Time: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - sort_start_time).count() << " s\n";

			if (data.empty())
				return;

			auto cur_kmer = data[0].kmer;
			uint32_t cur_gene_id = data[0].gene_id;
			uint32_t cur_variant_id = data[0].variant_id;
			uint32_t n_genes = 1;
			uint32_t n_variants = 1;

			uint64_t tot_unique_kmers = 1;

			for (size_t i = 1 ; i < data.size() ; ++i)
			{
				if (data[i].kmer != cur_kmer)
				{
					if (n_variants == 1)
						report_unique_kmer(cur_kmer, cur_variant_id);
					else if (n_genes == 1)
						report_unique_kmer(cur_kmer, cur_gene_id);
					else
						report_kmer(cur_kmer);

					cur_kmer = data[i].kmer;
					cur_gene_id = data[i].gene_id;
					cur_variant_id = data[i].variant_id;
					n_genes = 1;
					n_variants = 1;

					++tot_unique_kmers;
				}
				else if (data[i].gene_id != cur_gene_id)
				{
					assert(data[i].variant_id != cur_variant_id);
					++n_genes;
					++n_variants;
					cur_variant_id = data[i].variant_id;
					cur_gene_id = data[i].gene_id;
				}
				else if (data[i].variant_id != cur_variant_id)
				{
					++n_variants;
					cur_variant_id = data[i].variant_id;
				}
			}

			//last one
			if (n_variants == 1)
				report_unique_kmer(cur_kmer, cur_variant_id);
			else if (n_genes == 1)
				report_unique_kmer(cur_kmer, cur_gene_id);
			else
				report_kmer(cur_kmer);

			std::cerr << "Number of unique k-mers: " << tot_unique_kmers << "\n";
		}
	};
}

#endif // !_TRANSCRIPTOME_KMER_FINDER_H

