#ifndef _CBC_TO_CELL_TYPE
#define _CBC_TO_CELL_TYPE

#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include "../common/types/satc_data.h"

class CBCToCellType {
	std::unordered_map<std::string, uint32_t> cell_name_to_cell_id;
	std::unordered_map<uint64_t, uint32_t> sample_id_cbc_to_cell_id;
	static void parse_line(const std::string& line, std::string& cell_type, std::string& cbc, std::string& sample_name) {
		auto trim = [](std::string& str) {
			str.erase(0, str.find_first_not_of(" \n\r\t"));
			str.erase(str.find_last_not_of(" \n\r\t") + 1);
		};
		std::istringstream iss(line);
		if (!std::getline(iss, cell_type, ',')) {
			std::cerr << "Error: cannot read cell type name from line " << line << "\n";
			exit(1);
		}

		if (!std::getline(iss, cbc, ',')) {
			std::cerr << "Error: cannot read cbc from line " << line << "\n";
			exit(1);
		}

		if (!std::getline(iss, sample_name)) {
			std::cerr << "Error: cannot read sample name from line " << line << "\n";
			exit(1);
		}

		trim(cell_type);
		trim(cbc);
		trim(sample_name);
	}
public:
	static std::vector<uint32_t> get_barcodes_for_sample_name(const std::string& cbc_to_cell_type_samplesheet, const std::string& sample_name) {
		std::ifstream in(cbc_to_cell_type_samplesheet);
		if (!in) {
			std::cerr << "Error: cannot open file " << cbc_to_cell_type_samplesheet << "\n";
			exit(1);
		}
		std::string line;
		std::vector<uint32_t> res;
		while (std::getline(in, line)) {
			std::string cell_name;
			std::string str_cbc;
			std::string _sample_name;
			parse_line(line, cell_name, str_cbc, _sample_name);
			if (_sample_name == sample_name)
				res.push_back(str_kmer_to_uint64_t(str_cbc));
		}
		res.shrink_to_fit();
		return res;
	}
	CBCToCellType(const std::string& cbc_to_cell_type_samplesheet, const SampleNameToId& sample_name_to_id) {
		std::ifstream in(cbc_to_cell_type_samplesheet);
		if (!in) {
			std::cerr << "Error: cannot open file " << cbc_to_cell_type_samplesheet << "\n";
			exit(1);
		}
		std::string line;
		uint32_t cur_id{};
		while (std::getline(in, line)) {
			std::string cell_name;
			std::string str_cbc;
			std::string sample_name;
			parse_line(line, cell_name, str_cbc, sample_name);
			uint32_t sample_id;
			if (!sample_name_to_id.get_sample_id(sample_name, sample_id)) {
				std::cerr << "Error: cannot find sample id for sample name " << sample_name << "\n";
				exit(1);
			}

			uint32_t cbc = static_cast<uint32_t>(str_kmer_to_uint64_t(str_cbc));

			auto sample_id_cbc = pack_smaple_id_target(sample_id, cbc);

			auto [it, new_elem] = cell_name_to_cell_id.emplace(cell_name, cur_id);
			if (new_elem)
				++cur_id;
			uint32_t cell_type_id = it->second;

			if (auto [it2, new_elem] = sample_id_cbc_to_cell_id.emplace(sample_id_cbc, cell_type_id); !new_elem) {
				std::cerr << "Error: (barcode, sample) pair " << str_cbc << "," << sample_name << " maps to more than one cell type\n";
				exit(1);
			}
		}
	}

	bool has_sample_id_barcode(uint64_t sample_id_barcode) const {
		return sample_id_cbc_to_cell_id.find(sample_id_barcode) != sample_id_cbc_to_cell_id.end();
	}

	//mkokot_TODO: remove this methid if not needed
	//may not exist
	bool get_cell_type_id(uint64_t sample_id_barcode, uint32_t& cell_type_id) const {
		auto it = sample_id_cbc_to_cell_id.find(sample_id_barcode);
		if (it != sample_id_cbc_to_cell_id.end()) {
			cell_type_id = it->second;
			return true;
		}
		return false;
	}

	//must exist
	uint32_t get_cell_type_id(uint64_t sample_id_barcode) const {
		return sample_id_cbc_to_cell_id.find(sample_id_barcode)->second;
	}

	bool get_cell_type_name_from_sample_id_barcode(uint32_t sample_id_barcode, std::string& cell_type_name) const {
		//mkokot_TODO: this is inefficient and should be only used for testing, if needed in production add appropriate std::unordered_map<uint32_t, std::string>
		uint32_t cell_id_to_find;
		if (!get_cell_type_id(sample_id_barcode, cell_id_to_find))
			return false;

		if (!get_cell_type_name_from_cell_type_id(cell_id_to_find, cell_type_name)) {
			std::cerr << "Error: this should never happen: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
		return true;
	}

	bool get_cell_type_name_from_cell_type_id(uint32_t cell_type_id, std::string& cell_type_name) const {
		for (const auto& [cell_name, cell_id] : cell_name_to_cell_id)
			if (cell_id == cell_type_id) {
				cell_type_name = cell_name;
				return true;
			}
		return false;
	}

	size_t get_n_cell_types() const {
		return cell_name_to_cell_id.size();
	}
};


#endif // !_CBC_TO_CELL_TYPE
