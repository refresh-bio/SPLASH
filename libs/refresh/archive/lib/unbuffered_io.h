#ifndef _UNBUFFERED_IO_H
#define _UNBUFFERED_IO_H

#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <cinttypes>

#include <iostream>

#ifndef _WIN32
#define REFRESH_fseek	fseek
#define REFRESH_ftell	ftell
#else
#define REFRESH_fseek	_fseeki64
#define REFRESH_ftell	_ftelli64
#include <fcntl.h>
#include <io.h>
#endif

namespace refresh
{
// *******************************************************************************************
// Unbuffered input file - just a wrapper for FILE* 
class unbuffered_in_file
{
	FILE* f;
	
	size_t file_size_;

	bool use_stdin;
	
public:
	// *******************************************************************************************
	unbuffered_in_file() : f(nullptr), file_size_(0), use_stdin(false)
	{};

	// *******************************************************************************************
	~unbuffered_in_file()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name, const size_t _BUFFER_SIZE = 128 << 20)
	{
		if (f)
			return false;

		bool use_stdin = file_name.empty();

		if (use_stdin)
		{
			f = stdin;
#ifdef _WIN32
			_setmode(_fileno(f), _O_BINARY);
#endif
		}
		else
		{
			if (!(f = fopen(file_name.c_str(), "rb")))
				return false;
		}
		REFRESH_fseek(f, 0, SEEK_END);
		file_size_ = REFRESH_ftell(f);
		REFRESH_fseek(f, 0, SEEK_SET);

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		if (f)
		{
			if(!use_stdin)
				fclose(f);
			f = nullptr;
		}		
		return true;
	}

	// *******************************************************************************************
	bool opened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	int get()
	{
		return getc(f);
	}

	// *******************************************************************************************
	bool un_get()
	{
		if (REFRESH_ftell(f))
		{
			REFRESH_fseek(f, -1, SEEK_CUR);
			return true;
		}

		return false;
	}

	// *******************************************************************************************
	uint64_t read_uint(const int no_bytes)
	{
		uint64_t x = 0;
		uint64_t shift = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			uint64_t c = get();
			x += c << shift;
			shift += 8;
		}

		return x;
	}

	// *******************************************************************************************
	void read(uint8_t* ptr, size_t size)
	{
		fread(ptr, 1, size, f);		
	}

	// *******************************************************************************************
	bool eof() const
	{
		return feof(f);
	}

	// *******************************************************************************************
	bool seek(const size_t requested_pos)
	{
		REFRESH_fseek(f, requested_pos, SEEK_SET);		
		return true;
	}

	// *******************************************************************************************
	size_t file_size() const
	{
		if (f)
			return file_size_;
		else
			return 0;
	}

	// *******************************************************************************************
	size_t get_pos() const
	{
		return REFRESH_ftell(f);
	}
};

// *******************************************************************************************
// Unbuffered output file - just a wrapper for FILE* 
class unbuffered_out_file
{	
	FILE* f;	
	bool success;
	bool use_stdout;

public:
	// *******************************************************************************************
	unbuffered_out_file() : f(nullptr), success(false), use_stdout(false)
	{};

	// *******************************************************************************************
	~unbuffered_out_file()
	{
		close();
	}

	// *******************************************************************************************
	bool open(const std::string& file_name, const size_t _BUFFER_SIZE = 8 << 20)
	{
		if (f)
			return false;

		use_stdout = file_name.empty();

		if (use_stdout)
		{
			f = stdout;
#ifdef _WIN32
			_setmode(_fileno(f), _O_BINARY);
#endif
		}
		else
		{
			if (!(f = fopen(file_name.c_str(), "wb")))
				return false;
		}

		success = true;

		return true;
	}

	// *******************************************************************************************
	bool close()
	{
		if (f)
		{
			fflush(f);

			if(!use_stdout)
				fclose(f);
			f = nullptr;
		}		
		return success;
	}

	// *******************************************************************************************
	bool opened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	void put(const uint8_t c)
	{
		success = putc(c, f) == EOF;		
	}

	// *******************************************************************************************
	void write(const uint8_t* p, size_t n)
	{
		success &= fwrite(p, 1, n, f) == n;
	}

	// *******************************************************************************************
	void write_uint(const uint64_t _x, const int no_bytes)
	{
		uint64_t x = _x;

		for (int i = 0; i < no_bytes; ++i)
		{
			put(x & 0xff);
			x >>= 8;
		}
	}
};
}
// EOF
#endif