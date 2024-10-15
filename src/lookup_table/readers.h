#ifndef _READERS_H
#define _READERS_H
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <sstream>

class MultiFastaReader {
	std::ifstream in;
	std::string current_header;

	std::string line; //member to avoid reallocations
public:
	MultiFastaReader(const std::string& path) :
		in(path) {
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		if (!std::getline(in, current_header)) {
			std::cerr << "Error: cannot read first line from " << path << "\n";
			exit(1);
		}
		if (current_header[0] != '>') {
			std::cerr << "Error: wrong symbol, expected '>', but have " << current_header[0] << "\n";
			exit(1);
		}
	}

	bool NextSeq(std::string& header, std::string& seq) {
		header.clear();
		seq.clear();
		if (current_header == "")
			return false;
		header = current_header.substr(1); //without '>;
		current_header.clear();

		while (std::getline(in, line)) {
			if (line[0] == '>') //new header
			{
				current_header = line;
				break;
			}
			seq += line;
		}
		return true;
	}
};

class ISeqReader {
public:
	virtual bool NextSeq(std::string& seq) = 0;
	virtual void Init() { }
	virtual ~ISeqReader() = default;
};

class SeqReaderMultiFasta : public ISeqReader {
	MultiFastaReader reader;
	std::string header;
public:
	SeqReaderMultiFasta(const std::string& path) :
		reader(path) { }

	bool NextSeq(std::string& seq) override {
		return reader.NextSeq(header, seq);
	}
};

class ISeqReaderExtendorsEventListener {
public:
	virtual void NotifyNewLine(const std::string& line) = 0;
	virtual void NotifyLineEnd() = 0;
	virtual ~ISeqReaderExtendorsEventListener() = default;
};

class SeqReaderExtendors : public ISeqReader {
	std::vector<std::string> col_names;
	std::string anchor;
	std::vector<std::string> targets;
	size_t cur_target_id{};
	bool finished = false;

	std::vector<size_t> target_cols_ids;
	std::ifstream in;
	std::string line;

	ISeqReaderExtendorsEventListener* seq_reader_extendors_event_listener{};
	bool is_col_with_target(const std::string& col_name) {
		const static std::regex pattern("most_freq_target_([1-9][0-9]*)");
		return std::regex_match(col_name, pattern);
	}

	bool init_extendor() {
		cur_target_id = 0;

		if (seq_reader_extendors_event_listener)
			seq_reader_extendors_event_listener->NotifyLineEnd();

		if (!(std::getline(in, line)) || line == "")
			return false;

		if (seq_reader_extendors_event_listener)
			seq_reader_extendors_event_listener->NotifyNewLine(line);

		std::istringstream iss(line);

		if (!(iss >> anchor)) {
			std::cerr << "Error: cannot read anchor from line " << line;
			exit(1);
		}
		size_t id = 1;
		size_t i = 0;
		std::string val;
		while (iss >> val && i < target_cols_ids.size()) {
			if (id == target_cols_ids[i]) {
				targets[i++] = val;
			}
			++id;
		}

		if (i != target_cols_ids.size()) {
			std::cerr << "Error: cannot read all targets from line " << line << "\n";
			exit(1);
		}

		return true;
	}

public:
	SeqReaderExtendors(const std::string& path) : in(path) {
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		if (!std::getline(in, line)) {
			std::cerr << "Error: cannot read header from file " << path << "\n";
			exit(1);
		}
		std::istringstream iss(line);

		std::string name;
		while (iss >> name)
			col_names.emplace_back(std::move(name));

		if (col_names[0] != "anchor") {
			std::cerr << "Error: first columns should be \"anchor\" but is " << col_names[0] << "\n";
			exit(1);
		}

		for (size_t i = 0; i < col_names.size(); ++i) {
			if (is_col_with_target(col_names[i]))
				target_cols_ids.push_back(i);
		}

		if (target_cols_ids.empty()) {
			std::cerr << "Error: this file does not contain targets (non of columns start with pattern \"most_freq_target_\"\n";
			exit(1);
		}

		targets.resize(target_cols_ids.size());
	}

	void SetListener(ISeqReaderExtendorsEventListener* seq_reader_extendors_event_listener) {
		this->seq_reader_extendors_event_listener = seq_reader_extendors_event_listener;
	}

	//cannot be done in ctor, because there may be not listener yet
	void Init() override {
		finished = !init_extendor();
	}

	const std::vector<std::string>& GetColNames() const {
		return col_names;
	}

	const size_t GetNTargets() const {
		return target_cols_ids.size();
	}
	const std::string& GetLastLine() const {
		return line;
	}
	bool NextSeq(std::string& seq) override {
		while (true) {

			if (finished)
				return false;

			if (cur_target_id < targets.size()) {
				const std::string& target = targets[cur_target_id++];
				if (target == "-")
					continue;
				seq = anchor + target;
				return true;
			}

			finished = !init_extendor();
		}
	}
};

class SeqReaderCompactors : public ISeqReader {
	std::ifstream in;

	std::string header;
	std::string current_line;

	size_t compactor_pos = std::numeric_limits<size_t>::max();

public:
	SeqReaderCompactors(const std::string& path) :in(path) {
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		if (!std::getline(in, header)) {
			std::cerr << "Error: cannot read header from file " << path << "\n";
			exit(1);
		}

		std::istringstream iss(header);

		size_t i{};
		std::string header_part;
		while (iss >> header_part) {
			if (header_part == "compactor") {
				compactor_pos = i;
				break;
			}
			++i;
		}

		if (compactor_pos == std::numeric_limits<size_t>::max()) {
			std::cerr << "Error: cannot find column with name \"compactor\"\n";
			exit(1);
		}
	}
	bool NextSeq(std::string& seq) override {

		if (!std::getline(in, current_line))
			return false;

		if (current_line == "")
			return false;

		std::istringstream iss(current_line);

		size_t i{};
		while (iss >> seq) {
			if (i == compactor_pos) {
				return true;
			}
			++i;
		}

		std::cerr << "Error: cannot read compactor from line " << current_line << "\n";
		exit(1);
	}

	const std::string& GetHeader() const {
		return header;
	}

	const std::string& GetLastLine() const {
		return current_line;
	}
};

#endif // !_READERS_H

