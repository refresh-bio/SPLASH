#include "reader.h"

#include <iostream>

using namespace std;

void Reader::addRead()
{
	total_bases += currentLine.length();
	current_reads_bytes += currentLine.size();

	reads.emplace_back(currentLine);

	if (current_reads_bytes >= reads_pack_size)
	{
		current_reads_bytes = 0;
		reads_queue.push(std::move(reads));
	}
}

void Reader::porcessFastaOrMultiFasta(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile)
{
	enum class FastaReadState { Header, EOLsAfterHeader, Read, EOLsAfterOrInsideRead };
	FastaReadState fastaReadState = FastaReadState::Header;
	while (readed)
	{
		for (uint32_t pos = 0; pos < readed; ++pos)
		{
			auto symb = buff[pos];
			bool is_eol = symb == '\n' || symb == '\r';
			switch (fastaReadState)
			{
			case FastaReadState::Header:
				if (is_eol)
					fastaReadState = FastaReadState::EOLsAfterHeader;
				break;
			case FastaReadState::EOLsAfterHeader:
				if (!is_eol)
				{
					currentLine.push_back(symb);
					fastaReadState = FastaReadState::Read;
				}
				break;
			case FastaReadState::Read:
				if (is_eol)
					fastaReadState = FastaReadState::EOLsAfterOrInsideRead;
				else
					currentLine.push_back(symb);
				break;
			case FastaReadState::EOLsAfterOrInsideRead:
				if (!is_eol)
				{
					if (symb == '>')
					{
						currentLine.shrink_to_fit();
						addRead();
						currentLine.clear();
						fastaReadState = FastaReadState::Header;
					}
					else
						fastaReadState = FastaReadState::Read;
					currentLine.push_back(symb);
				}
				break;
			default:
				break;
			}
		}
		readed = gzread(gzfile, buff.data(), buf_size);
		total_bytes += readed;
	}
	currentLine.shrink_to_fit();
	addRead();
	currentLine.clear();
}

void Reader::processFastq(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile)
{
	enum class WhereInRead { read_header, read, qual_header, qual };
	WhereInRead whereInRead = WhereInRead::read_header;
	bool was_eol = false;
	uint32_t record_lines = 4; //4 lines form fastq record

	while (readed)
	{
		uint32_t pos = 0;
		while (pos < readed)
		{
			if (buff[pos] == '\n' || buff[pos] == '\r') // EOL reached
			{
				if (was_eol) //we are skipping windows EOL
					++pos;
				else
				{
					if (whereInRead == WhereInRead::read) {
						currentLine.shrink_to_fit();
						addRead();
						currentLine.clear();
					}

					whereInRead = (WhereInRead)(((int)whereInRead + 1) % record_lines);
					++pos;
				}
				was_eol = true;
			}
			else {
				if (whereInRead == WhereInRead::read)
					currentLine.push_back(buff[pos]);
				++pos;
				was_eol = false;
			}
		}
		readed = gzread(gzfile, buff.data(), buf_size);
		total_bytes += readed;
	}
}


Reader::Reader(const std::string& path, refresh::parallel_queue<read_pack_t>& reads_queue) :
	reads_queue(reads_queue)
{
	auto gzfile = gzopen(path.c_str(), "rb");
	if (!gzfile)
	{
		cerr << "Error: cannot open file: " << path << "\n";
		exit(1);
	}

	const uint32_t buf_size = 1ul << 25;
	std::vector<uint8_t> buff(buf_size);

	uint64_t readed = gzread(gzfile, buff.data(), buf_size);
	total_bytes += readed;
	if (!readed)
	{
		int code;
		auto errmsg = gzerror(gzfile, &code);
		if (code < 0)
		{
			std::cerr << "zblib error: " << errmsg << "\n";
			exit(1);
		}
		std::cerr << "Error: file " << path << " is empty\n";
		exit(1);
	}
	if (buff[0] != '@' && buff[0] != '>')
	{
		std::cerr << "Error: unknown file format\n";
		exit(1);
	}
	is_fastq = buff[0] == '@';

	//FASTA or multi-fasta
	if (!is_fastq)
		porcessFastaOrMultiFasta(buff, readed, buf_size, gzfile);
	else
		processFastq(buff, readed, buf_size, gzfile);

	int code;
	auto errmsg = gzerror(gzfile, &code);
	if (code < 0)
	{
		std::cerr << "zblib error: " << errmsg << "\n";
		exit(1);
	}

	if (!currentLine.empty())
	{
		std::cerr << "Error: something went wrong during input reading\n";
		exit(1);
	}
	gzclose(gzfile);


	if (reads.size())
		reads_queue.push(std::move(reads));
	
	reads_queue.mark_completed();
}