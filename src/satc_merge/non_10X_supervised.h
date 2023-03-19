#ifndef _NON_10X_SUPERVISED_H
#define _NON_10X_SUPERVISED_H

#include <string>
#include <sstream>
#include "../common/satc_data.h"

class Non10XSupervised {	
	std::vector<std::vector<double>> Cjs;//first index is Cj_variant, second is sample_id
	void parse_line(const std::string& line, std::string& sample_name, std::vector<double>& Cjs) {
		Cjs.clear();
		auto trim = [](std::string& str) {
			str.erase(0, str.find_first_not_of(" \n\r\t"));
			str.erase(str.find_last_not_of(" \n\r\t") + 1);
		};
		std::istringstream iss(line);
		if (!std::getline(iss, sample_name, ',')) {
			std::cerr << "Error: cannot read sample name from line " << line << "\n";
			exit(1);
		}
		trim(sample_name);
		double Cj;
		char comma;
		while (iss >> Cj) {
			Cjs.push_back(Cj);
			iss >> comma;
		}
	}
public:
	Non10XSupervised(const std::string& Cj_samplesheet, const SampleNameToId& sample_name_to_id) {
		std::ifstream in(Cj_samplesheet);
		if (!in) {
			std::cerr << "Error: cannot open file " << Cj_samplesheet << "\n";
			exit(1);
		}
		std::string line;
		uint32_t sample_id;
		std::string sample_name;
		std::vector<double> _Cjs;
		
		//first line
		size_t tot_lines = 0;
		std::vector<int> was_seen; //per each sample id 1 if this sample was already read
		if (std::getline(in, line)) {
			++tot_lines;
			parse_line(line, sample_name, _Cjs);
			Cjs.resize(_Cjs.size(), std::vector<double>(sample_name_to_id.get_n_samples()));
			was_seen.resize(sample_name_to_id.get_n_samples());
			
			if (!sample_name_to_id.get_sample_id(sample_name, sample_id)) {
				std::cerr << "Error: cannot find sample id for sample name " << sample_name << "\n";
				exit(1);
			}
			was_seen[sample_id] = true;
			for (size_t i = 0; i < _Cjs.size(); ++i)
				Cjs[i][sample_id] = _Cjs[i];
		}
		while (std::getline(in, line)) {
			++tot_lines;
			parse_line(line, sample_name, _Cjs);
			if (_Cjs.size() != Cjs.size()) {
				std::cerr << "Error: different number of Cjs values in line: " << line << "\n";
				exit(1);
			}
			if (!sample_name_to_id.get_sample_id(sample_name, sample_id)) {
				std::cerr << "Error: cannot find sample id for sample name " << sample_name << "\n";
				exit(1);
			}
			if (was_seen[sample_id]) {
				std::cerr << "Error: sample " << sample_name << " occurs more than once in file " << Cj_samplesheet << "\n";
				exit(1);
			}
			was_seen[sample_id] = true;
			for (size_t i = 0; i < _Cjs.size(); ++i)
				Cjs[i][sample_id] = _Cjs[i];
		}
		if (tot_lines != sample_name_to_id.get_n_samples()) {
			std::cerr << "Error: there are samples for which Cjs are not defined in file: " << Cj_samplesheet << "\n";
			exit(1);
		}
	}
	const std::vector<double>& get(size_t Cj_idx) const {
		return Cjs[Cj_idx];
	}

	size_t get_n_Cjs() const {
		return Cjs.size();
	}
};

#endif // !_NON_10X_SUPERVISED_H

