#include <iostream>
#include <chrono>
#include "query_common.h"
#include "lookup.h"
#include "readers.h"

namespace query_common {

	//quite ugly because similar class is in satc_data, but it was extended for compression
	//so here is simplified variant
	class simple_buffered_binary_writer {
		std::ofstream out;
		std::vector<uint8_t> buff;

		void flush() {
			out.write(reinterpret_cast<char*>(buff.data()), buff.size());
			buff.clear();
		}

		void assure_space(size_t size) {
			if (buff.size() + size > buff.capacity()) {
				flush();
				if (size > buff.capacity())
					buff.reserve(size);
			}
		}

		void write(const std::vector<uint8_t>& vec) {
			assure_space(vec.size());
			for (auto x : vec)
				buff.push_back(x);
		}

	public:
		simple_buffered_binary_writer(simple_buffered_binary_writer&&) = default;
		simple_buffered_binary_writer& operator=(simple_buffered_binary_writer&& rhs) noexcept {
			if (&rhs == this)
				return *this;
			close();

			out = std::move(rhs.out);
			buff = std::move(rhs.buff);

			return *this;
		}
		simple_buffered_binary_writer(const std::string& path, size_t buff_size = 1ull << 25) {
			out.rdbuf()->pubsetbuf(nullptr, 0);
			out.open(path, std::ios::binary);
			buff.reserve(buff_size);
		}
		operator bool() const {
			return out.operator bool();
		}
		void write(const uint8_t* ptr, size_t size) {
			assure_space(size);
			buff.insert(buff.end(), ptr, ptr + size);
		}

		void write(const char* ptr, size_t size) {
			write(reinterpret_cast<const uint8_t*>(ptr), size);
		}

		void write(const char* ptr) {
			write(ptr, strlen(ptr));
		}

		void write(const std::string& str) {
			write(str.c_str(), str.size());
		}

		void write(char c)
		{
			assure_space(1);
			buff.push_back(static_cast<uint8_t>(c));
		}

		~simple_buffered_binary_writer() {
			if (!buff.empty())
				flush();
		}
		void close() {
			if (!buff.empty())
				flush();
			out.close();
		}
	};

	std::string to_string(InputFmt input_fmt) {
		switch (input_fmt)
		{
		case InputFmt::FASTA:
			return "fasta";
		case InputFmt::EXTENDORS:
			return "extendors";
		case InputFmt::COMPACTORS:
			return "compactors";
		default:
			std::cerr << "Error: unsupported input format\n";
			exit(1);
		}
	}

	InputFmt input_format_from_string(const std::string& fmt) {
		if (fmt == "fasta")
			return InputFmt::FASTA;
		if (fmt == "extendors")
			return InputFmt::EXTENDORS;
		if (fmt == "compactors")
			return InputFmt::COMPACTORS;

		std::cerr << "Error: unknown input format: " << fmt << "\n";
		exit(1);
	}

	std::string to_string(ReportFmt report_fmt) {
		switch (report_fmt)
		{
		case ReportFmt::PLAIN:
			return "plain";
		case ReportFmt::CONCISE:
			return "concise";
		case ReportFmt::IDS:
			return "ids";
		case ReportFmt::EMPTY:
			return "empty";
		default:
			std::cerr << "Error: unsupported output format\n";
			exit(1);
		}
	}

	ReportFmt report_format_from_string(const std::string& fmt) {
		if (fmt == "plain")
			return ReportFmt::PLAIN;
		if (fmt == "concise")
			return ReportFmt::CONCISE;
		if (fmt == "ids")
			return ReportFmt::IDS;
		if (fmt == "empty")
			return ReportFmt::EMPTY;

		std::cerr << "Error: unknown report format: " << fmt << "\n";
		exit(1);
	}

	std::string to_string(StatsFmt stats_format) {
		switch (stats_format)
		{
		case StatsFmt::EMPTY:
			return "empty";
		case StatsFmt::WITH_STATS:
			return "with_stats";
		default:
			std::cerr << "Error: unsupported stats format\n";
			exit(1);
		}
	}
	StatsFmt stats_format_from_string(const std::string& stats_format) {
		if (stats_format == "empty")
			return StatsFmt::EMPTY;
		if (stats_format == "with_stats")
			return StatsFmt::WITH_STATS;

		std::cerr << "Error: unknown stats format: " << stats_format << "\n";
		exit(1);
	}

	std::string to_string(OutputFmt output_format) {
		switch (output_format)
		{
		case OutputFmt::TXT:
			return "txt";
		case OutputFmt::EXTENDORS:
			return "extendors";
		case OutputFmt::COMPACTORS:
			return "compactors";
		default:
			std::cerr << "Error: unsupported output format\n";
			exit(1);
		}
	}
	OutputFmt output_format_from_string(const std::string& output_format) {
		if (output_format == "txt")
			return OutputFmt::TXT;
		if (output_format == "extendors")
			return OutputFmt::EXTENDORS;
		if (output_format == "compactors")
			return OutputFmt::COMPACTORS;

		std::cerr << "Error: unknown output format: " << output_format << "\n";
		exit(1);
	}

	void SingleInputConfig::Print(std::ostream& oss) const {
		oss << "input                          : " << input << "\n";
		oss << "input_fmt                      : " << to_string(input_fmt) << "\n";
		oss << "report_fmt                     : " << to_string(report_fmt) << "\n";
		oss << "stats_fmt                      : " << to_string(stats_fmt) << "\n";
		oss << "output_fmt                     : " << to_string(output_fmt) << "\n";
		oss << "output                         : " << output << "\n";
		oss << "kmer_skip                      : " << kmer_skip << "\n";
	}

	SingleInputConfig::SingleInputConfig(const std::vector<std::string>& desc) {
		size_t i = 0;
		for (; i < desc.size(); ++i)
		{
			if (desc[i][0] != '-')
				break;

			std::string param = desc[i];

			if (param == "--input_fmt") {
				input_fmt = input_format_from_string(desc[++i]);
			}
			if (param == "--report_fmt") {
				report_fmt = report_format_from_string(desc[++i]);
			}
			if (param == "--stats_fmt") {
				stats_fmt = stats_format_from_string(desc[++i]);
			}
			if (param == "--output_fmt") {
				output_fmt = output_format_from_string(desc[++i]);
			}
			if (param == "--kmer_skip") {
				kmer_skip = std::stoull(desc[++i]);
			}
		}

		if (i >= desc.size()) {
			std::cerr << "Error: input missing\n";
			exit(1);
		}
		input = desc[i++];

		if (i >= desc.size()) {
			std::cerr << "Error: output missing\n";
			exit(1);
		}
		output = desc[i++];

		if (output_fmt == OutputFmt::EXTENDORS && input_fmt != InputFmt::EXTENDORS) {
			std::cerr << "Error: --output_fmt extendors may be used only with --input_format extendors\n";
			exit(1);
		}

		if (output_fmt == OutputFmt::COMPACTORS && input_fmt != InputFmt::COMPACTORS) {
			std::cerr << "Error: --output_fmt compactors may be used only with --input_format compactors\n";
			exit(1);
		}
	}

	void QueryCfg::Print(std::ostream& oss) const {
		oss << "lookup                         : " << lookup << "\n";
		oss << "n_threads                      : " << n_threads << "\n";

		if (input_cfg.size() == 1) {
			input_cfg.front().Print(oss);
		}
		else {
			oss << "Inputs:\n";
			for (const auto& icf : input_cfg) {
				oss << "-------------------------------------------------\n";
				icf.Print(oss);
			}
		}
	}

	std::unique_ptr<ISeqReader> make_reader(InputFmt input_fmt, const std::string& path) {
		switch (input_fmt)
		{
		case InputFmt::FASTA:
			return std::make_unique<SeqReaderMultiFasta>(path);
		case InputFmt::EXTENDORS:
			return std::make_unique<SeqReaderExtendors>(path);
		case InputFmt::COMPACTORS:
			return std::make_unique<SeqReaderCompactors>(path);
		default:
			std::cerr << "Error: unsupported input format, please contact authors showing this message: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
	}
	//mkokot_TODO: workaround by introducing template - probably to be removed
	template<typename LOOKUP_T>
	std::unique_ptr<ICategoryReportMaker> make_report_maker(ReportFmt report_fmt, const LOOKUP_T& lookup) {
		switch (report_fmt)
		{
		case ReportFmt::PLAIN:
			return std::make_unique<PlainReportMaker>(lookup.GetFileMapper(), lookup.GetHeaderMappers());
		case ReportFmt::CONCISE:
			return std::make_unique<ConciseReportMaker>(lookup.GetFileMapper(), lookup.GetHeaderMappers());
		case ReportFmt::IDS:
			return std::make_unique<IdsReportMaker>(lookup.GetHeaderMappers());
		case ReportFmt::EMPTY:
			return std::make_unique<EmptyReportMaker>();
		default:
			std::cerr << "Error: unsupported report format, please contact authors showing this message: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
	}

	std::unique_ptr<IStatsReportMaker> make_stats_maker(StatsFmt stats_fmt) {
		switch (stats_fmt)
		{
		case StatsFmt::EMPTY:
			return std::make_unique<EmptyStatsReportMaker>();
		case StatsFmt::WITH_STATS:
			return std::make_unique<SimpleStatsReportMaker>();
		default:
			std::cerr << "Error: unsupported stats format, please contact authors showing this message: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
	}


	struct QueryDesc
	{
		std::string line; //for extendors and compactors
		std::vector<std::string> queries;
		std::string result;
	};

	struct QueryTaskData
	{
		uint64_t priority{};
		std::vector<QueryDesc> task_data;
	};

	using QueryTaskQueue = refresh::parallel_queue<QueryTaskData, refresh::limited_stl_queue<QueryTaskData>>;
	using QueryResultQueue = refresh::parallel_priority_queue<std::vector<QueryDesc>>;


	class IOutputLineFormatter
	{
	public:
		virtual void BuildString(std::string& res, ICategoryReportMaker& report_maker, IStatsReportMaker& stats_report_maker) = 0;
		virtual void Complete(std::string& res, ICategoryReportMaker& report_maker, IStatsReportMaker& stats_report_maker) = 0;
		virtual ~IOutputLineFormatter() = default;
	};

	class IOutputFormatter {
	public:
		virtual void Wait() = 0;
		virtual std::unique_ptr<IOutputLineFormatter> GetLineFormatter() = 0;
		virtual void RegisterQuery(const std::string& query) = 0;
		virtual ~IOutputFormatter() = default;
	};

	class BasicOutputLineFormatter : public IOutputLineFormatter {
		const std::string report_stats_separator;
	public:
		BasicOutputLineFormatter(const std::string& report_stats_separator):
			report_stats_separator(report_stats_separator)
		{

		}
		void BuildString(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override {

			auto cat_str = report_maker.ToString();
			auto stats_str = stats_report_maker.ToString();
			res += cat_str;
			res += report_stats_separator;
			res += stats_str;
		}

		void Complete(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override
		{
			res.push_back('\n');
		}
	};

	class BasicOutputFormatter : public IOutputFormatter {
		simple_buffered_binary_writer& out;
		ICategoryReportMaker& category_reporter;
		IStatsReportMaker& stats_reporter;
		const std::string report_stats_separator;

		QueryTaskQueue& query_task_queue;
		QueryResultQueue& query_result_queue;
		std::thread worker_thread;

		QueryTaskData query_task_data;
		uint32_t max_pack_size = 100; //mkokot_TODO: adjust
	public:
		BasicOutputFormatter(
			simple_buffered_binary_writer& out,
			ICategoryReportMaker& category_reporter,
			IStatsReportMaker& stats_reporter,
			const std::string& report_stats_separator,
			QueryTaskQueue& query_task_queue,
			QueryResultQueue& query_result_queue) :
			out(out),
			category_reporter(category_reporter),
			stats_reporter(stats_reporter),
			report_stats_separator(report_stats_separator),
			query_task_queue(query_task_queue),
			query_result_queue(query_result_queue)
		{
			worker_thread = std::thread([this]
				{
					std::vector<QueryDesc> query_result;
					while (this->query_result_queue.pop(query_result))
					{
						for (auto& elem : query_result)
						{
							this->out.write(elem.result);
						}
					}
				});
		}

		std::unique_ptr<IOutputLineFormatter> GetLineFormatter() override
		{
			return  std::make_unique<BasicOutputLineFormatter>(report_stats_separator);
		}
		void RegisterQuery(const std::string& query) override
		{
			query_task_data.task_data.emplace_back();
			//query_task_data.task_data.back().line not needed here
			query_task_data.task_data.back().queries.push_back(query);
			if (query_task_data.task_data.size() >= max_pack_size) {
				auto cur_prior = query_task_data.priority;
				query_task_queue.push(std::move(query_task_data));
				query_task_data.priority = cur_prior + 1;
				query_task_data.task_data.clear(); //just in case
			}
		}
		void Wait() override
		{
			if (!query_task_data.task_data.empty())
				query_task_queue.push(std::move(query_task_data));
			query_task_queue.mark_completed();

			worker_thread.join();
		}
	};

	class ExtendorsOutputLineFormatter : public IOutputLineFormatter
	{
		size_t newline_each_n_writes;
		size_t n_build_calls{};
	public:
		ExtendorsOutputLineFormatter(size_t newline_each_n_writes):
			newline_each_n_writes(newline_each_n_writes)
		{

		}
		void BuildString(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override {

			++n_build_calls;

			auto cat_str = report_maker.ToString();
			auto stats_str = stats_report_maker.ToString();
			if (!cat_str.empty()) {
				res.push_back('\t');
				res += cat_str;
			}
			if (!stats_str.empty()) {
				res.push_back('\t');
				res += stats_str;
			}
		}

		void Complete(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override
		{
			if (n_build_calls < newline_each_n_writes)
			{
				auto cat_str = report_maker.GetEmpty();
				auto stats_str = stats_report_maker.GetEmpty();
				for (auto i = n_build_calls; i < newline_each_n_writes; ++i)
				{
					if (cat_str != "") {
						res.push_back('\t');
						res += cat_str;
					}

					if (stats_str != "") {
						res.push_back('\t');
						res += stats_str;
					}
				}
			}
			res.push_back('\n');
			n_build_calls = 0;
		}
	};

	class ExtendorsOutputFormatter : public IOutputFormatter, public ISeqReaderExtendorsEventListener {
		simple_buffered_binary_writer& out;
		ICategoryReportMaker& category_reporter;
		IStatsReportMaker& stats_reporter;
		SeqReaderExtendors& seq_reader_extendors;
		size_t newline_each_n_writes;
		size_t to_next_newline;

		QueryTaskQueue& query_task_queue;
		QueryResultQueue& query_result_queue;
		std::thread worker_thread;

		QueryTaskData query_task_data;
		uint32_t max_pack_size = 100; //mkokot_TODO: adjust
	public:
		ExtendorsOutputFormatter(simple_buffered_binary_writer& out,
			ICategoryReportMaker& category_reporter,
			IStatsReportMaker& stats_reporter,
			SeqReaderExtendors& seq_reader_extendors,
			QueryTaskQueue& query_task_queue,
			QueryResultQueue& query_result_queue
			) :
			out(out),
			category_reporter(category_reporter),
			stats_reporter(stats_reporter),
			seq_reader_extendors(seq_reader_extendors),
			newline_each_n_writes(seq_reader_extendors.GetNTargets()),
			to_next_newline(newline_each_n_writes),
			query_task_queue(query_task_queue),
			query_result_queue(query_result_queue)
		{
			const auto& col_names = seq_reader_extendors.GetColNames();
			for (size_t i = 0; i < col_names.size(); ++i) {
				out.write(col_names[i]);
				if (i + 1 != col_names.size())
					out.write('\t');
			}

			size_t n_targets = seq_reader_extendors.GetNTargets();

			if (category_reporter.IsEnabled() || stats_reporter.IsEnabled()) {
				for (size_t i = 0; i < n_targets; ++i) {

					std::string col_name_query_res = std::string("most_freq_extendor_") + std::to_string(i + 1) + "_query_res";
					std::string col_name_query_stats = std::string("most_freq_extendor_") + std::to_string(i + 1) + "_query_stats";

					if (category_reporter.IsEnabled())
					{
						out.write('\t');
						out.write(col_name_query_res);
					}
					if (stats_reporter.IsEnabled())
					{
						out.write('\t');
						out.write(col_name_query_stats);
					}
				}
			}

			seq_reader_extendors.SetListener(this);

			out.write('\n');

			worker_thread = std::thread([this]
				{
					std::vector<QueryDesc> query_result;
					while (this->query_result_queue.pop(query_result))
					{
						for (auto& elem : query_result)
						{
							this->out.write(elem.line);
							this->out.write(elem.result);
						}
					}
				});
		}

		void NotifyNewLine(const std::string& line) override
		{
			query_task_data.task_data.emplace_back();
			query_task_data.task_data.back().line = line;
		}
		void NotifyLineEnd() override
		{
			if (query_task_data.task_data.size() >= max_pack_size) {
				auto cur_prior = query_task_data.priority;
				query_task_queue.push(std::move(query_task_data));
				query_task_data.priority = cur_prior + 1;
				query_task_data.task_data.clear(); //just in case
			}
		}
		std::unique_ptr<IOutputLineFormatter> GetLineFormatter() override
		{
			return std::make_unique<ExtendorsOutputLineFormatter>(newline_each_n_writes);
		}

		void RegisterQuery(const std::string& query) override
		{
			query_task_data.task_data.back().queries.push_back(query);
		}

		void Wait() override
		{
			if (!query_task_data.task_data.empty())
				query_task_queue.push(std::move(query_task_data));
			query_task_queue.mark_completed();

			worker_thread.join();
		}
	};
	class CompactorsOutputLineFormatter : public IOutputLineFormatter
	{
	public:
		void BuildString(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override {

			if (report_maker.IsEnabled()) {
				res.push_back('\t');
				res += report_maker.ToString();
			}

			if (stats_report_maker.IsEnabled())
			{
				res.push_back('\t');
				res += stats_report_maker.ToString();
			}
		}

		void Complete(std::string& res,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_report_maker) override
		{
			res.push_back('\n');
		}
	};

	class CompactorsOutputFormatter : public IOutputFormatter {
		simple_buffered_binary_writer& out;
		ICategoryReportMaker& category_reporter;
		IStatsReportMaker& stats_reporter;
		SeqReaderCompactors& seq_reader_compactors;

		QueryTaskQueue& query_task_queue;
		QueryResultQueue& query_result_queue;
		std::thread worker_thread;

		QueryTaskData query_task_data;
		uint32_t max_pack_size = 100; //mkokot_TODO: adjust
	public:
		CompactorsOutputFormatter(simple_buffered_binary_writer& out,
			ICategoryReportMaker& category_reporter,
			IStatsReportMaker& stats_reporter,
			SeqReaderCompactors& seq_reader_compactors,
			QueryTaskQueue& query_task_queue,
			QueryResultQueue& query_result_queue
			) :
			out(out),
			category_reporter(category_reporter),
			stats_reporter(stats_reporter),
			seq_reader_compactors(seq_reader_compactors),
			query_task_queue(query_task_queue),
			query_result_queue(query_result_queue) {

			out.write(seq_reader_compactors.GetHeader());
			if (category_reporter.IsEnabled())
				out.write("\tquery_res");
			if (stats_reporter.IsEnabled())
				out.write("\tquery_stats");
			out.write('\n');

			worker_thread = std::thread([this]
				{
					std::vector<QueryDesc> query_result;
					while (this->query_result_queue.pop(query_result))
					{
						for (auto& elem : query_result)
						{
							this->out.write(elem.line);
							this->out.write(elem.result);
						}
					}
				});
		}

		std::unique_ptr<IOutputLineFormatter> GetLineFormatter() override
		{
			return std::make_unique<CompactorsOutputLineFormatter>();
		}
		void RegisterQuery(const std::string& query) override
		{
			query_task_data.task_data.emplace_back();
			query_task_data.task_data.back().line = seq_reader_compactors.GetLastLine();
			query_task_data.task_data.back().queries.push_back(query);

			if (query_task_data.task_data.size() >= max_pack_size) {
				auto cur_prior = query_task_data.priority;
				query_task_queue.push(std::move(query_task_data));
				query_task_data.priority = cur_prior + 1;
				query_task_data.task_data.clear(); //just in case
			}
		}
		void Wait() override
		{
			if (!query_task_data.task_data.empty())
				query_task_queue.push(std::move(query_task_data));
			query_task_queue.mark_completed();

			worker_thread.join();
		}
	};

	std::unique_ptr<IOutputFormatter> make_output_formatter(
			simple_buffered_binary_writer& out,
			ICategoryReportMaker& report_maker,
			IStatsReportMaker& stats_maker,
			OutputFmt output_fmt,
			const std::string& report_stats_separator,
			ISeqReader* seq_reader,
			QueryTaskQueue& query_task_queue,
			QueryResultQueue& query_result_queue) {

		switch (output_fmt)
		{
		case OutputFmt::TXT:
			return std::make_unique<BasicOutputFormatter>(
				out,
				report_maker,
				stats_maker,
				report_stats_separator,
				query_task_queue,
				query_result_queue);
		case OutputFmt::EXTENDORS: {
			//mkokot_TODO: this is terribly wrong desing if I need dynamic_cast at this early stage of a project
			SeqReaderExtendors* seq_reader_extendors = dynamic_cast<SeqReaderExtendors*>(seq_reader);
			assert(seq_reader_extendors);
			return std::make_unique<ExtendorsOutputFormatter>(
				out,
				report_maker,
				stats_maker,
				*seq_reader_extendors,
				query_task_queue,
				query_result_queue
				);
		}
		case OutputFmt::COMPACTORS: {
			SeqReaderCompactors* seq_reader_compactors = dynamic_cast<SeqReaderCompactors*>(seq_reader);
			assert(seq_reader_compactors);
			return std::make_unique<CompactorsOutputFormatter>(
				out,
				report_maker,
				stats_maker,
				*seq_reader_compactors,
				query_task_queue,
				query_result_queue
			);
		}
		default:
			std::cerr << "Error: unsupported output format, please contact authors showing this message: " << __FILE__ << "(" << __LINE__ << ")\n";
			exit(1);
		}
	}

	std::string make_report_stats_separator(ReportFmt report_fmt, StatsFmt stats_fmt) {
		std::string res = "";
		if (report_fmt != ReportFmt::EMPTY && stats_fmt != StatsFmt::EMPTY)
			res = "\t";
		return res;
	}

	class QueryExecutor
	{
		QueryTaskQueue& tasks_queue;
		QueryResultQueue& results_queue;
		Lookup& lookup;
		uint32_t kmer_skip;
		std::unique_ptr<ICategoryReportMaker> report_maker;
		std::unique_ptr<IStatsReportMaker> stats_maker;
		std::unique_ptr<IOutputLineFormatter> output_line_formatter;
	public:
		QueryExecutor(QueryTaskQueue& tasks_queue,
			QueryResultQueue& results_queue,
			Lookup& lookup,
			uint32_t kmer_skip,
			ReportFmt report_fmt,
			StatsFmt stats_fmt,
			std::unique_ptr<IOutputLineFormatter>&& output_line_formatter)
			:
			tasks_queue(tasks_queue),
			results_queue(results_queue),
			lookup(lookup),
			kmer_skip(kmer_skip),
			report_maker(make_report_maker(report_fmt, lookup)),
			stats_maker(make_stats_maker(stats_fmt)),
			output_line_formatter(std::move(output_line_formatter))
		{

		}
		void Process()
		{
			QueryTaskData task;
			QueryResult query_result;
			while (tasks_queue.pop(task))
			{
				auto priority = task.priority;
				for (auto& elem: task.task_data)
				{
					for (auto& query : elem.queries)
					{
						//SBWT requires capital letters
						for (auto& c : query)
							c = std::toupper(c);

						lookup.query_seq(query, kmer_skip, query_result);

						query_result.Dump(*report_maker, *stats_maker, lookup.GetHeaderMappers());

						output_line_formatter->BuildString(elem.result, *report_maker, *stats_maker);
					}
					output_line_formatter->Complete(elem.result, *report_maker, *stats_maker);
				}
				results_queue.push(priority, std::move(task.task_data));
			}

			results_queue.mark_completed();
		}
	};

	void run(const QueryCfg& query_cfg) {
		auto create_lookup_table_object_start = std::chrono::high_resolution_clock::now();
		Lookup lookup(query_cfg.lookup, query_cfg.n_threads);

		std::cerr << "Construct lookup table object time: "
			<< std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - create_lookup_table_object_start).count()
			<< " s\n";

		//lookup.DumpMapping("mapping.log"); //mkokot_TODO: keep this? remove?

		auto start_time_all = std::chrono::high_resolution_clock::now();

		for (const auto& cfg : query_cfg.input_cfg) {

			std::cerr << "start querying " << cfg.input << "...";
			auto start_time = std::chrono::high_resolution_clock::now();

			auto seq_reader = make_reader(cfg.input_fmt, cfg.input);

			simple_buffered_binary_writer out(cfg.output);

			if (!out) {
				std::cerr << "Error: cannot open output file " << cfg.output << "\n";
				exit(1);
			}

			QueryTaskQueue query_task_queue(query_cfg.n_threads * 2, 1, "query_task_queue");
			QueryResultQueue query_result_queue(query_cfg.n_threads * 2, query_cfg.n_threads, "query_result_queue");

			//mkokot_TODO: this is a little strange that I need to create also global object for formatter...
			//I think per threads should be fine, but require some adjustments
			auto report_maker = make_report_maker(cfg.report_fmt, lookup);
			auto report_stats_separator = make_report_stats_separator(cfg.report_fmt, cfg.stats_fmt);
			auto stats_maker = make_stats_maker(cfg.stats_fmt);

			//mkokot_TODO: actualy output_formater is not very nice name now, since it make much more
			auto output_formatter = make_output_formatter(
				out,
				*report_maker,
				*stats_maker,
				cfg.output_fmt,
				report_stats_separator,
				seq_reader.get(),
				query_task_queue,
				query_result_queue);

			//!!! it must be after output formater is created, because there is some listener to be registered
			seq_reader->Init();

			//mkokot_TODO: apply thread_control here?

			//start threads
			std::vector<std::thread> threads(query_cfg.n_threads);
			for(int i = 0 ; i < query_cfg.n_threads; ++i)
			{
				threads[i] = std::thread([&]
					{
						QueryExecutor executor(
							query_task_queue, 
							query_result_queue, 
							lookup,
							cfg.kmer_skip,
							cfg.report_fmt,
							cfg.stats_fmt,
							output_formatter->GetLineFormatter());

						executor.Process();
					});
			}

			using time_diff_t = decltype((std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now()) + (std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now()));
			time_diff_t next_seq_time{};

			auto next_seq_start = std::chrono::high_resolution_clock::now();

			std::string seq;
			while (seq_reader->NextSeq(seq)) {

				output_formatter->RegisterQuery(seq);

				auto next_seq_duration = std::chrono::high_resolution_clock::now() - next_seq_start;
				next_seq_time += next_seq_duration;

				next_seq_start = std::chrono::high_resolution_clock::now();
			}
			//must be before thread joins
			output_formatter->Wait();

			for (auto& th : threads)
				th.join();

			auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
			std::cerr << "\nDone (time: " << time << " s)\n";

			auto next_seq_time_s = std::chrono::duration<double>(next_seq_time).count();
			std::cerr << "Next seq time: " << next_seq_time_s << " s\n";
		}

		auto time_all = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time_all).count();
		std::cerr << "All queries time: " << time_all << " s)\n";
	}
}
