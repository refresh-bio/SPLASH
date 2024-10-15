#ifndef _REPORTERS_H
#define _REPORTERS_H
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "seq_id_mapper.h"

struct kmer_desc {
	uint32_t category;
	uint32_t file_id;
	uint32_t meta;
	kmer_desc(uint32_t category, uint32_t file_id, uint32_t meta) :
		category(category),
		file_id(file_id),
		meta(meta)
	{
	}

	bool operator==(const kmer_desc& rhs) const
	{
		return category == rhs.category && file_id == rhs.file_id && meta == rhs.meta;
	}
};

struct SeqStats {
private:
public:
	size_t in_cat1{};
	size_t in_cat2{};
	size_t in_cat3{};
	size_t in_cat4{};
	size_t in_poly{};
	size_t in_illegal{};
	size_t in_unknown{};

	void to_string(std::string& res) const {
		res.clear();
		res.push_back('1'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_cat1);
		res.push_back(','); res.push_back(' '); res.push_back('2'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_cat2);
		res.push_back(','); res.push_back(' '); res.push_back('3'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_cat3);
		res.push_back(','); res.push_back(' '); res.push_back('4'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_cat4);
		res.push_back(','); res.push_back(' '); res.push_back('U'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_unknown);
		res.push_back(','); res.push_back(' '); res.push_back('P'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_poly);
		res.push_back(','); res.push_back(' '); res.push_back('N'); res.push_back(':'); res.push_back(' '); res += std::to_string(in_illegal);
	}
};

class IStatsReportMaker {
public:
	virtual void ReportStats(const SeqStats& stats) = 0;

	virtual ~IStatsReportMaker() = default;

	virtual std::string ToString() = 0;

	virtual std::string GetEmpty() = 0;

	virtual bool IsEnabled() = 0;
};

class EmptyStatsReportMaker : public IStatsReportMaker {
public:
	void ReportStats(const SeqStats& stats) override { }

	std::string ToString() override { return ""; }

	std::string GetEmpty() override { return ""; }

	bool IsEnabled() override { return false; }
};

class SimpleStatsReportMaker : public IStatsReportMaker {
	std::string out;
public:
	void ReportStats(const SeqStats& stats) override {
		stats.to_string(out);
	}

	std::string ToString() override {
		return out;
	}

	std::string GetEmpty() override { return "-"; }

	bool IsEnabled() override { return true; }
};

class ICategoryReportMaker {
public:
	virtual void ReportCategory1(uint32_t header_id) = 0;
	virtual void ReportCategory2(uint32_t file_id, uint32_t cnt) = 0;
	virtual void ReportCategory3(const std::vector<kmer_desc>& V) = 0;
	virtual void ReportCategory4() = 0;

	virtual void ReportIllegal() = 0;
	virtual void ReportPolyACGT() = 0;

	virtual void ReportUnknown() = 0;
	virtual void MakeSeparator() = 0; //mkokot_TODO: remove

	virtual std::string ToString() = 0;

	virtual std::string GetEmpty() = 0;

	virtual bool IsEnabled() = 0;

	virtual ~ICategoryReportMaker() = default;
};

//to make it all work in parallel I need to break some dependencies in the code, so the query will just return a common query result that will be foramtted after
//the data here is compacted
//first we have category, and then it depends on the category what is next
//In general this provides som of the methods of ICategoryReportMaker
class QueryResult
{
	std::vector<uint32_t> data;
public:
	void Clear()
	{
		data.clear();
	}

	void ReportCategory1(uint32_t header_id)
	{
		data.push_back(1); //category
		data.push_back(header_id);
	}
	void ReportCategory2(uint32_t file_id, uint32_t cnt)
	{
		data.push_back(2); //category
		data.push_back(file_id);
		data.push_back(cnt);
	}
	void ReportCategory3(const std::vector<kmer_desc>& V)
	{
		data.push_back(3);
		assert(V.size() <= std::numeric_limits<uint32_t>::max());
		uint32_t casted_size = static_cast<uint32_t>(V.size());
		data.push_back(casted_size);

		for (const auto& elem : V)
			data.push_back(elem.meta);
	}
	void ReportCategory4()
	{
		data.push_back(4);
	}

	void ReportIllegal()
	{
		data.push_back((uint32_t)'N');
	}
	void ReportPolyACGT()
	{
		data.push_back((uint32_t)'P');
	}

	void ReportUnknown()
	{
		data.push_back((uint32_t)'U');
	}

	void Dump(ICategoryReportMaker& reporter, IStatsReportMaker& stats_reporter, const HeaderIdMapper& header_mapper)
	{
		uint64_t i = 0;
		SeqStats stats;
		while(i < data.size())
		{
			auto cat = data[i++];
			if(cat == 1)
			{
				auto header_id = data[i++];
				reporter.ReportCategory1(header_id);
				++stats.in_cat1;
			}
			else if (cat == 2)
			{
				auto file_id = data[i++];
				auto cnt = data[i++];
				reporter.ReportCategory2(file_id, cnt);
				++stats.in_cat2;
			}
			else if(cat == 3)
			{
				auto n = data[i++];
				std::vector<kmer_desc> V;
				V.reserve(n);
				for (uint32_t j = 0; j < n; ++j) {
					uint32_t header_id = data[i++];
					V.emplace_back(3, header_mapper.Decode(header_id).file_id, header_id);
				}

				reporter.ReportCategory3(V);
				++stats.in_cat3;
			}
			else if(cat == 4)
			{
				reporter.ReportCategory4();
				++stats.in_cat4;
			}
			else if(cat == (uint32_t)'N')
			{
				reporter.ReportIllegal();
				++stats.in_illegal;
			}
			else if(cat == (uint32_t)'P')
			{
				reporter.ReportPolyACGT();
				++stats.in_poly;
			}
			else if(cat == (uint32_t)'U')
			{
				reporter.ReportUnknown();
				++stats.in_unknown;
			}

			if (i < data.size())
				reporter.MakeSeparator();
		}
		stats_reporter.ReportStats(stats);
	}
};

class EmptyReportMaker : public ICategoryReportMaker {
public:
	EmptyReportMaker() {

	}
	void ReportCategory1(uint32_t header_id) override { }

	void ReportCategory2(uint32_t file_id, uint32_t cnt) override { }

	void ReportCategory3(const std::vector<kmer_desc>& V) override { }

	void ReportCategory4() override { }

	void ReportIllegal() override { }
	void ReportPolyACGT() override { }

	void ReportUnknown() override { }
	void MakeSeparator() override { }

	std::string ToString() override { return ""; }

	std::string GetEmpty() override { return ""; }

	bool IsEnabled() override { return false; }
};

class PlainReportMaker : public ICategoryReportMaker {
	const SeqIdMapper& file_mapper;
	const HeaderIdMapper& header_mapper;
	std::string out;
public:

	PlainReportMaker(const SeqIdMapper& file_mapper, const HeaderIdMapper& header_mapper) :
		file_mapper(file_mapper),
		header_mapper(header_mapper) {

	}

public:
	void ReportCategory1(uint32_t header_id) override {
		const auto& tmp = header_mapper.Decode(header_id);
		//category
		out.push_back('1'); out.push_back(':'); out.push_back(' ');
		out.push_back('(');
		out.push_back('"'); out += file_mapper.Decode(tmp.file_id); out.push_back('"');
		out.push_back(','); out.push_back(' ');
		out.push_back('"'); out += tmp.header; out.push_back('"');
		out.push_back(')');
	}

	void ReportCategory2(uint32_t file_id, uint32_t cnt) override {
		//category
		out.push_back('2'); out.push_back(':'); out.push_back(' ');
		out.push_back('(');
		out.push_back('"'); out += file_mapper.Decode(file_id);
		out.push_back('"'); out.push_back(','); out.push_back(' ');
		out += std::to_string(cnt);
		out.push_back(')');
	}

	void ReportCategory3(const std::vector<kmer_desc>& V) override {
		//category
		out.push_back('3'); out.push_back(':'); out.push_back(' ');
		out.push_back('[');

		uint32_t cur_file_id = V[0].file_id;
		out.push_back('"'); out += file_mapper.Decode(cur_file_id); out.push_back('"');
		out.push_back(':'); out.push_back(' '); out.push_back('['); out.push_back('"'); out += header_mapper.Decode(V[0].meta).header; out.push_back('"');
		for (size_t i = 1; i < V.size(); ++i) {
			auto file_id = V[i].file_id;
			if (file_id != cur_file_id) {
				cur_file_id = file_id;
				out.push_back(']'); out.push_back(','); out.push_back(' '); out.push_back('"'); out += file_mapper.Decode(cur_file_id); out.push_back('"'); out.push_back(':'); out.push_back(' '); out.push_back('['); out.push_back('"'); out += header_mapper.Decode(V[i].meta).header; out.push_back('"');
			}
			else {
				out.push_back(','); out.push_back(' '); out.push_back('"'); out += header_mapper.Decode(V[i].meta).header; out.push_back('"');
			}
		}
		out.push_back(']');
		out.push_back(']');
	}

	void ReportCategory4() override {
		out.push_back('4');
	}

	void ReportIllegal() override {
		out.push_back('N');
	}
	void ReportPolyACGT() override {
		out.push_back('P');
	}

	void ReportUnknown() override {
		out.push_back('U');
	}
	void MakeSeparator() override {
		out.push_back(',');
		out.push_back(' ');
	}

	std::string ToString() override {
		auto res = out;
		out.clear();
		return res;
	}

	std::string GetEmpty() override { return "-"; }

	bool IsEnabled() override { return true; }
};

class ConciseReportMaker : public ICategoryReportMaker {
	const SeqIdMapper& file_mapper;
	const HeaderIdMapper& header_mapper;
	std::ostringstream out; //mkokot_TODO: !!! change to std::string similar like in PlainReportMaker

	struct
	{
		int cat;
		uint32_t pos = 0; //how many calls for a single query anything was reported
		uint32_t start_pos = 0; //of the current range
		uint32_t n_kmers = 0;
		uint32_t header_id; // cat 1
		uint32_t file_id; // cat 2
		uint32_t cnt; //cat 2
		std::vector<kmer_desc> desc; //cat 3
		bool separator_needed = false;

		void reset()
		{
			cat = -1;
			start_pos = pos;
			n_kmers = 0;
		}
	} current_state;


	void flush()
	{
		auto cat = current_state.cat;
		if (cat == -1)
		{
			current_state.reset();
			return;
		}

		out << cat << ":" << "<" << current_state.start_pos << ":" << current_state.n_kmers << ">";

		if (cat == 1)
		{
			//only header
			const auto& header = header_mapper.Decode(current_state.header_id).header;
			out << ":\"" << header << "\"";
		}
		else if (cat == 2)
		{
			//only file
			const auto& file = file_mapper.Decode(current_state.file_id);
			out << ":\"" << file << "\"";
		}
		else if (cat == 3)
		{
			//only headers
			out << ":[";
			//for (const auto& elem : current_state.desc)
			for (size_t i = 0 ; i < current_state.desc.size() ; ++i)
			{
				out << "\"";
				out << header_mapper.Decode(current_state.desc[i].meta).header;
				out << "\"";
				if (i != current_state.desc.size() - 1)
					out << ", ";
			}
			out << "]";
		}
		else if (cat == 4)
		{
			//for 4 we don't report anything additional
		}

		current_state.reset();
		current_state.separator_needed = true;
	}
public:

	ConciseReportMaker(const SeqIdMapper& file_mapper, const HeaderIdMapper& header_mapper) :
		file_mapper(file_mapper),
		header_mapper(header_mapper) {
		current_state.reset();
	}

public:
	void ReportCategory1(uint32_t header_id) override {
		if (current_state.cat != 1 || current_state.header_id != header_id)
		{
			flush();
			current_state.cat = 1;
			current_state.header_id = header_id;
		}
		++current_state.n_kmers;
		++current_state.pos;
	}

	void ReportCategory2(uint32_t file_id, uint32_t cnt) override {
		//mkokot_TODO: i ignore cnt in this reporting mechanism, not sure if that is fine
		if (current_state.cat != 2 || current_state.file_id != file_id)
		{
			flush();
			current_state.cat = 2;
			current_state.file_id = file_id;
		}
		++current_state.n_kmers;
		++current_state.pos;
	}

	void ReportCategory3(const std::vector<kmer_desc>& V) override {

		if (current_state.cat != 3 || V != current_state.desc)
		{
			flush();
			current_state.cat = 3;
			current_state.desc = V;
		}
		++current_state.n_kmers;
		++current_state.pos;
	}

	void ReportCategory4() override {
		if (current_state.cat != 4)
		{
			flush();
			current_state.cat = 4;
		}
		++current_state.n_kmers;
		++current_state.pos;
		//out << "4";
	}

	void ReportIllegal() override {
		if (current_state.cat != -1)
			flush();
		++current_state.pos;
	}
	void ReportPolyACGT() override {
		if (current_state.cat != -1)
			flush();
		++current_state.pos;
	}

	void ReportUnknown() override {
		if (current_state.cat != -1)
			flush();
		++current_state.pos;
	}
	void MakeSeparator() override {
		if (current_state.separator_needed)
		{
			out << ", ";
			current_state.separator_needed = false;
		}
	}

	std::string ToString() override {
		if (current_state.cat != -1)
			flush();
		auto res = out.str();
		out.clear();
		out.str("");
		current_state.pos = 0;
		current_state.reset();
		current_state.separator_needed = false;
		return res;
	}

	std::string GetEmpty() override
	{
		//mkokot_TODO: reimplement ConciseReportMaker
		return "-";
	}

	bool IsEnabled() override { return true; }
};


class IdsReportMaker : public ICategoryReportMaker {
	const HeaderIdMapper& header_mapper;
	std::ostringstream out; //mkokot_TODO: !!! change to std::string similar like in PlainReportMaker
public:

	IdsReportMaker(const HeaderIdMapper& header_mapper) :
		header_mapper(header_mapper) {

	}
	void ReportCategory1(uint32_t header_id) override {
		out << "1: "; //category
		out << "(";
		out << header_mapper.Decode(header_id).file_id;
		out << ", ";
		out << header_id;
		out << ")";
	}

	void ReportCategory2(uint32_t file_id, uint32_t cnt) override {
		out << "2: "; //category
		out << "(";
		out << file_id;
		out << ", ";
		out << cnt;
		out << ")";
	}

	void ReportCategory3(const std::vector<kmer_desc>& V) override {
		out << "3: "; //category
		out << "[";

		uint32_t cur_file_id = V[0].file_id;
		out << cur_file_id;
		out << ": [" << V[0].meta;
		for (size_t i = 1; i < V.size(); ++i) {
			auto file_id = V[i].file_id;
			if (file_id != cur_file_id) {
				cur_file_id = file_id;
				out << "], " << cur_file_id << ": [" << V[i].meta;
			}
			else {
				out << ", " << V[i].meta;
			}
		}
		out << "]";
		out << "]";
	}

	void ReportCategory4() override {
		out << "4";
	}

	void ReportIllegal() override {
		out << "N";
	}
	void ReportPolyACGT() override {
		out << "P";
	}

	void ReportUnknown() override {
		out << "U";
	}
	void MakeSeparator() override {
		out << ", ";
	}

	std::string ToString() override {
		auto res = out.str();
		out.clear();
		out.str("");
		return res;
	}

	std::string GetEmpty() override { return "-"; }

	bool IsEnabled() override { return true; }
};

#endif // !_REPORTERS_H
