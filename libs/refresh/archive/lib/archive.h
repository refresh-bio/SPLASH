#ifndef _ARCHIVE_H
#define _ARCHIVE_H

#include <cstdio>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>
#include <thread>
#include <mutex>
#include <cinttypes>
#include "buffered_io.h"
#include "unbuffered_io.h"

namespace refresh
{
	const uint32_t REFRESH_BUILD_ARCHIVE = 2;

template<class T_INPUT, class T_OUTPUT>
class archive_basic
{
public:
	struct part_t {
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		std::string stream_name;
		size_t cur_id;
		size_t raw_size;
		size_t packed_size;
		size_t packed_data_size;
		std::vector<part_t> parts;
	} stream_t;

private:
	bool input_mode;
	T_INPUT f_in;
	T_OUTPUT f_out;
	size_t io_buffer_size;

	size_t f_offset;

	std::map<int, std::vector<std::pair<std::vector<uint8_t>, uint64_t>>> m_buffer;

	std::vector<stream_t> v_streams;
	std::unordered_map<std::string, size_t> rm_streams;
	std::mutex mtx;

	bool serialize()
	{
		size_t footer_size = 0;

		// Store stram part offsets
		footer_size += write(v_streams.size());

		for (auto& stream : v_streams)
		{
			size_t p = footer_size;

			footer_size += write(stream.stream_name);
			footer_size += write(stream.parts.size());
			footer_size += write(stream.raw_size);

			for (auto& part : stream.parts)
			{
				footer_size += write(part.offset);
				footer_size += write(part.size);
			}

			stream.packed_size += footer_size - p;
		}

		write_fixed(footer_size);

		return true;
	}

	// *******************************************************************************************
	bool deserialize()
	{
		size_t footer_size;
		size_t file_size_ = f_in.file_size();

		f_in.seek(file_size_ - 8ull);
		read_fixed(footer_size);

		f_in.seek(file_size_ - (size_t)(8 + footer_size));

		// Read stream part offsets
		size_t n_streams;
		read(n_streams);

		v_streams.resize(n_streams, stream_t());

		for (size_t i = 0; i < n_streams; ++i)
		{
			auto& stream_second = v_streams[i];

			read(stream_second.stream_name);
			read(stream_second.cur_id);
			read(stream_second.raw_size);

			stream_second.parts.resize(stream_second.cur_id);
			for (size_t j = 0; j < stream_second.cur_id; ++j)
			{
				read(stream_second.parts[j].offset);
				read(stream_second.parts[j].size);
			}

			stream_second.cur_id = 0;

			rm_streams[stream_second.stream_name] = i;
		}

		f_in.seek(0);

		return true;
	}

	// *******************************************************************************************
	template<typename T>
	size_t write_fixed(const T x)
	{
		f_out.write_uint(static_cast<uint64_t>(x), 8);

		return 8;
	}

	// *******************************************************************************************
	template<typename T>
	size_t write(const T _x)
	{
		int no_bytes = 0;
		uint64_t x = static_cast<uint64_t>(_x);

		for (size_t tmp = x; tmp; tmp >>= 8)
			++no_bytes;

		f_out.put(no_bytes);

		for (int i = no_bytes; i; --i)
			f_out.put((x >> ((i - 1) * 8)) & 0xff);

		return no_bytes + 1;
	}

	// *******************************************************************************************
	size_t write(const std::string& s)
	{
//		f_out.write(s);
		f_out.write((uint8_t*)s.c_str(), s.size());
		f_out.put(0);

		return s.size() + 1;
	}

	// *******************************************************************************************
	template<typename T>
	size_t read_fixed(T& x)
	{
		x = static_cast<T>(f_in.read_uint(8));

		return 8;
	}

	// *******************************************************************************************
	size_t read(std::string& s)
	{
		s.clear();

		while (true)
		{
			int c = f_in.get();
			if (c == EOF)
				return 0;

			if (c == 0)
				return s.size() + 1;

			s.push_back((char)c);
		}

		return 0;
	}

	// *******************************************************************************************
	template<typename T>
	size_t read(T& x)
	{
		int no_bytes = f_in.get();

		x = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			x <<= 8;
			x += static_cast<T>(f_in.get());
		}

		return no_bytes + 1;
	}

	// *******************************************************************************************
	bool add_part_impl(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata)
	{
		v_streams[stream_id].parts.push_back(part_t(f_offset, data_size));

		f_offset += write(metadata);
		f_out.write(data, data_size);

		f_offset += data_size;

		v_streams[stream_id].packed_size += f_offset - v_streams[stream_id].parts.back().offset;
		v_streams[stream_id].packed_data_size += data_size;

		return true;
	}

	// *******************************************************************************************
	bool get_part_impl(const int stream_id, uint8_t* data, uint64_t& metadata)
	{
		auto& p = v_streams[stream_id];

		f_in.seek(p.parts[p.cur_id].offset);

		if (p.parts[p.cur_id].size != 0)
			read(metadata);
		else
		{
			metadata = 0;
			p.cur_id++;
			return true;
		}

		f_in.read(data, p.parts[p.cur_id].size);

		p.cur_id++;

		//if (r != p.parts[p.cur_id-1].size)
			//return false;

		//return r == p.parts[p.cur_id-1].size;
		return true;
	}

	// *******************************************************************************************
	size_t get_part_size_impl(const int stream_id, size_t& size)
	{
		auto& p = v_streams[stream_id];

		if (p.cur_id >= p.parts.size())
			return false;
		size = p.parts[p.cur_id].size;
		return true;
	}

	// *******************************************************************************************
	bool flush_out_buffers_imp()
	{
		for (auto& x : m_buffer)
			for (auto& y : x.second)
				add_part(x.first, y.first, y.second);

		m_buffer.clear();

		return true;
	}

public:
	// *******************************************************************************************
	archive_basic(const bool _input_mode, const size_t _io_buffer_size = 64 << 20)
	{
		input_mode = _input_mode;
		io_buffer_size = _io_buffer_size;
	}

	// *******************************************************************************************
	~archive_basic()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name) 
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (f_in.opened())
			f_in.close();
		if (f_out.opened())
			f_out.close();

		if (input_mode)
			f_in.open(file_name, io_buffer_size);
		else
			f_out.open(file_name, io_buffer_size);

		if (!f_in.opened() && !f_out.opened())
			return false;

		if (input_mode)
			deserialize();

		f_offset = 0;

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (!f_in.opened() && !f_out.opened())
			return false;

		if (input_mode)
			f_in.close();
		else
		{
			flush_out_buffers_imp();
			serialize();
			f_out.close();
		}

		return true;
	}

	// *******************************************************************************************
	int register_stream(const std::string& stream_name)
	{
		std::lock_guard<std::mutex> lck(mtx);

		// Before adding new stream check if stream_name is already registered
		auto p = rm_streams.find(stream_name);
		if (p != rm_streams.end())
			return (int)p->second;

		int id = (int)v_streams.size();

		v_streams.emplace_back(stream_t());

		v_streams[id].cur_id = 0;
		v_streams[id].stream_name = stream_name;
		v_streams[id].raw_size = 0;
		v_streams[id].packed_size = 0;
		v_streams[id].packed_data_size = 0;

		rm_streams[stream_name] = id;

		return id;
	}

	// *******************************************************************************************
	int get_stream_id(const std::string& stream_name)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto p = rm_streams.find(stream_name);
		if (p != rm_streams.end())
			return (int)p->second;

		return -1;
	}

	// *******************************************************************************************
	size_t get_stream_packed_size(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
			return 0;

		return v_streams[stream_id].packed_size;
	}

	// *******************************************************************************************
	size_t get_stream_packed_data_size(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
			return 0;

		return v_streams[stream_id].packed_data_size;
	}

	// *******************************************************************************************
	bool add_part(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_impl(stream_id, v_data.data(), v_data.size(), metadata);
	}

	// *******************************************************************************************
	bool add_part(const int stream_id, const uint8_t* data, const size_t data_size, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return add_part_impl(stream_id, data, data_size, metadata);
	}

	// *******************************************************************************************
	int add_part_prepare(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		v_streams[stream_id].parts.push_back(part_t(0, 0));

		return static_cast<int>(v_streams[stream_id].parts.size()) - 1;
	}

	// *******************************************************************************************
	bool add_part_complete(const int stream_id, const int part_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		v_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

		f_offset += write(metadata);
		f_out.write(v_data.data(), v_data.size());

		f_offset += v_data.size();

		v_streams[stream_id].packed_size += f_offset - v_streams[stream_id].parts[part_id].offset;
		v_streams[stream_id].packed_data_size += v_data.size();

		return true;
	}

	// *******************************************************************************************
	bool add_part_buffered(const int stream_id, const std::vector<uint8_t>& v_data, const uint64_t metadata = 0)
	{
		std::lock_guard<std::mutex> lck(mtx);

		m_buffer[stream_id].emplace_back(v_data, metadata);

		return true;
	}

	// *******************************************************************************************
	bool flush_out_buffers()
	{
		std::lock_guard<std::mutex> lck(mtx);

		return flush_out_buffers_imp();
	}

	// *******************************************************************************************
	bool get_part_size(const int stream_id, size_t& size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return get_part_size_impl(stream_id, size);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		size_t size;
		if (!get_part_size_impl(stream_id, size))
			return false;

		v_data.resize(size);

		return get_part_impl(stream_id, v_data.data(), metadata);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, uint8_t* data, size_t& data_size, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (!get_part_size_impl(stream_id, data_size))
			return false;

		return get_part_impl(stream_id, data, metadata);
	}

	// *******************************************************************************************
	bool get_part(const int stream_id, const int part_id, std::vector<uint8_t>& v_data, uint64_t& metadata)
	{
		std::lock_guard<std::mutex> lck(mtx);

		auto& p = v_streams[stream_id];

		if ((size_t)part_id >= p.parts.size())
			return false;

		v_data.resize(p.parts[part_id].size);

		f_in.seek(p.parts[part_id].offset);

		if (p.parts[part_id].size != 0)
			read(metadata);
		else
		{
			metadata = 0;
			return true;
		}

		f_in.read(v_data.data(), p.parts[part_id].size);

		return true;
	}

	// *******************************************************************************************
	void set_raw_size(const int stream_id, const size_t raw_size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		v_streams[stream_id].raw_size = raw_size;
	}

	// *******************************************************************************************
	size_t get_raw_size(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		return v_streams[stream_id].raw_size;
	}

	// *******************************************************************************************
	size_t get_no_streams()
	{
		std::lock_guard<std::mutex> lck(mtx);

		return v_streams.size();
	}

	// *******************************************************************************************
	size_t get_no_parts(const int stream_id)
	{
		std::lock_guard<std::mutex> lck(mtx);

		if (stream_id < 0 || (size_t)stream_id >= v_streams.size())
			return 0;

		return v_streams[stream_id].parts.size();
	}

	// *******************************************************************************************
	void list_streams(std::vector<stream_t>& _v_streams)
	{
		std::lock_guard<std::mutex> lck(mtx);

		_v_streams = v_streams;
	}
};

using archive_buffered_io = archive_basic<buffered_in_file, buffered_out_file>;
using archive_unbuffered_io = archive_basic<unbuffered_in_file, unbuffered_out_file>;

}
// EOF
#endif