#ifndef _ZSTD_FILE_H
#define _ZSTD_FILE_H

#include <cstdint>
#include <cstdio>
#include <string>

#include "../../libs/zstd/lib/zstd.h"

namespace refresh
{
	class zstd_file
	{
		enum class working_mode_t {none, reading, writing};

		int compression_level;
		working_mode_t working_mode;
		FILE* fio;

		size_t io_buffer_size;
		size_t zstd_buffer_size;

		char *mem_buf[2];

		ZSTD_CStream *zstd_cstream;
		ZSTD_DStream *zstd_dstream;

		ZSTD_inBuffer zstd_read_buffer;
		bool read_eof;

		void flush()
		{
			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			zstd_in_buffer.src = mem_buf[0];
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = 0;

			zstd_out_buffer.dst = mem_buf[1];
			zstd_out_buffer.pos = 0;
			zstd_out_buffer.size = zstd_buffer_size;

			for(bool need_work = true; need_work;)
			{
				zstd_out_buffer.pos = 0;
				need_work = ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_end) != 0;
				fwrite(zstd_out_buffer.dst, 1, zstd_out_buffer.pos, fio);
			}
		}

		void close_reading()
		{
			if (!fio)
				return;
			fclose(fio);
			fio = nullptr;

			ZSTD_freeDStream(zstd_dstream);
			zstd_dstream = nullptr;

			working_mode = working_mode_t::none;
		}

		void close_writing()
		{
			if (!fio)
				return;

			flush();

			fclose(fio);
			fio = nullptr;

			ZSTD_freeCStream(zstd_cstream);
			zstd_cstream = nullptr;

			working_mode = working_mode_t::none;
		}

		void allocate_buffers()
		{
			release_buffers();

			mem_buf[0] = new char[zstd_buffer_size];
			mem_buf[1] = new char[zstd_buffer_size];
		}

		void release_buffers()
		{
			if (mem_buf[0])
			{
				delete[] mem_buf[0];
				mem_buf[0] = nullptr;
			}

			if (mem_buf[1])
			{
				delete[] mem_buf[1];
				mem_buf[1] = nullptr;
			}
		}

	public:
		zstd_file(int compression_level = 9, size_t io_buffer_size = 1 << 20, size_t zstd_buffer_size = 1 << 20) : 
			compression_level(compression_level),
			working_mode(working_mode_t::none), 
			fio(nullptr), 
			io_buffer_size(io_buffer_size),
			zstd_buffer_size(zstd_buffer_size)
		{
			mem_buf[0] = nullptr;
			mem_buf[1] = nullptr;

			zstd_cstream = nullptr;
			zstd_dstream = nullptr;
		}

		zstd_file(zstd_file&& rhs)
		{
			working_mode = rhs.working_mode;
			rhs.working_mode = working_mode_t::none;

			fio = rhs.fio;
			rhs.fio = nullptr;

			compression_level = rhs.compression_level;

			io_buffer_size = rhs.io_buffer_size;
			zstd_buffer_size = rhs.zstd_buffer_size;

			mem_buf[0] = rhs.mem_buf[0];
			rhs.mem_buf[0] = nullptr;

			mem_buf[1] = rhs.mem_buf[1];
			rhs.mem_buf[1] = nullptr;

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			zstd_dstream = rhs.zstd_dstream;
			rhs.zstd_dstream = nullptr;

			zstd_read_buffer = rhs.zstd_read_buffer;
			rhs.zstd_read_buffer.src = nullptr;

			read_eof = rhs.read_eof;
		};

		zstd_file& operator=(zstd_file&& rhs)
		{
			if (fio)
				close();

			working_mode = rhs.working_mode;
			rhs.working_mode = working_mode_t::none;

			fio = rhs.fio;
			rhs.fio = nullptr;

			compression_level = rhs.compression_level;

			io_buffer_size = rhs.io_buffer_size;
			zstd_buffer_size = rhs.zstd_buffer_size;

			if (mem_buf[0])
				delete[] mem_buf[0];

			mem_buf[0] = rhs.mem_buf[0];
			rhs.mem_buf[0] = nullptr;

			if (mem_buf[1])
				delete[] mem_buf[1];

			mem_buf[1] = rhs.mem_buf[1];
			rhs.mem_buf[1] = nullptr;

			if (zstd_cstream)
				ZSTD_freeCStream(zstd_cstream);

			zstd_cstream = rhs.zstd_cstream;
			rhs.zstd_cstream = nullptr;

			if (zstd_dstream)
				ZSTD_freeDStream(zstd_dstream);

			zstd_dstream = rhs.zstd_dstream;
			rhs.zstd_dstream = nullptr;

			zstd_read_buffer = rhs.zstd_read_buffer;
			rhs.zstd_read_buffer.src = nullptr;

			read_eof = rhs.read_eof;

			return *this;
		};

		zstd_file(const zstd_file& rhs) = delete;

		~zstd_file()
		{
			close();

			release_buffers();
		}

		void close()
		{
			if (working_mode == working_mode_t::reading)
				close_reading();
			else if (working_mode == working_mode_t::writing)
				close_writing();

			release_buffers();
		}

		void set_compression_level(size_t _compression_level)
		{
			if(working_mode == working_mode_t::none)
				compression_level = _compression_level;
		}

		void set_io_buffer_size(size_t _io_buffer_size)
		{
			if (working_mode == working_mode_t::none)
				io_buffer_size = _io_buffer_size;
		}

		void set_zstd_buffer_size(size_t _zstd_buffer_size)
		{
			if (working_mode == working_mode_t::none)
			{
				zstd_buffer_size = _zstd_buffer_size;

				allocate_buffers();
			}
		}

		bool eof() const
		{
			return read_eof;
		}

		bool is_opened() const
		{
			return working_mode != working_mode_t::none;
		}

		bool is_opened_for_reading() const
		{
			return working_mode == working_mode_t::reading;
		}
		
		bool is_opened_for_writing() const
		{
			return working_mode == working_mode_t::writing;
		}

		bool open_reading(const std::string& file_name)
		{
			if (working_mode != working_mode_t::none)
				return false;

			fio = fopen(file_name.c_str(), "rb");
			if (!fio)
				return false;

			setvbuf(fio, nullptr, _IOFBF, io_buffer_size);

			working_mode = working_mode_t::reading;

			allocate_buffers();

			zstd_dstream = ZSTD_createDStream();

			zstd_read_buffer.src = mem_buf[0];
			zstd_read_buffer.pos = 0;
			zstd_read_buffer.size = 0;

			read_eof = false;

			ZSTD_initDStream(zstd_dstream);

			return true;
		}

		bool open_writing(const std::string& file_name)
		{
			if (working_mode != working_mode_t::none)
				return false;

			fio = fopen(file_name.c_str(), "wb");
			if (!fio)
				return false;

			setvbuf(fio, nullptr, _IOFBF, io_buffer_size);

			working_mode = working_mode_t::writing;

			allocate_buffers();

			zstd_cstream = ZSTD_createCStream();

			ZSTD_initCStream(zstd_cstream, compression_level);

			return true;
		}

		bool put(char c)
		{
			return write(&c, 1);
		}

		bool write(char* p, size_t size)
		{
			if (working_mode != working_mode_t::writing)
				return false;

			ZSTD_inBuffer zstd_in_buffer;
			ZSTD_outBuffer zstd_out_buffer;

			zstd_in_buffer.src = p;
			zstd_in_buffer.pos = 0;
			zstd_in_buffer.size = size;

			zstd_out_buffer.dst = mem_buf[1];
			zstd_out_buffer.pos = 0;
			zstd_out_buffer.size = zstd_buffer_size;

			while (zstd_in_buffer.pos < zstd_in_buffer.size)
			{
				ZSTD_compressStream2(zstd_cstream, &zstd_out_buffer, &zstd_in_buffer, ZSTD_e_continue);

				fwrite(zstd_out_buffer.dst, 1, zstd_out_buffer.pos, fio);
				zstd_out_buffer.pos = 0;
			}

			return true;
		}

		size_t get(char& c)
		{
			return read(&c, 1);
		}

		size_t read(char* p, size_t size)
		{
			if (working_mode != working_mode_t::reading || read_eof)
				return 0;

			size_t readed = 0;

			ZSTD_outBuffer zstd_out_buffer;

			while (readed < size && !read_eof)
			{
				if (zstd_read_buffer.pos == zstd_read_buffer.size)
				{
					zstd_read_buffer.pos = 0;
					zstd_read_buffer.size = fread((void*)zstd_read_buffer.src, 1, zstd_buffer_size, fio);
				}

				zstd_out_buffer.dst = p + readed;
				zstd_out_buffer.pos = 0;
				zstd_out_buffer.size = size - readed;

				read_eof = ZSTD_decompressStream(zstd_dstream, &zstd_out_buffer, &zstd_read_buffer) == 0;

				readed += zstd_out_buffer.pos;
			}

			return readed;
		}
	};
}

#endif
