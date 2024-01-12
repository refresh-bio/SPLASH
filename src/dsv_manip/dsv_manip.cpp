#include <iostream>
#include <string>
#include <limits>
#include <sstream>
#include "../common/version.h"
#include "sort_mode.h"
#include "limit_mode.h"


enum class Mode { sort, limit };
std::string to_string(Mode mode) {
	switch (mode)
	{
	case Mode::sort:
		return "sort";
	case Mode::limit:
		return "limit";
	default:
		std::cerr << "Error: unhandled mode " << (int)mode << "\n";
		exit(1);
		break;
	}
}

Mode mode_from_string(const std::string& mode) {
	if (mode == "sort") {
		return Mode::sort;
	}
	else if (mode == "limit") {
		return Mode::limit;
	}

	std::cerr << "Error: unknown mode " << mode << "\n";
	exit(1);
}

Mode get_mode_or_print_help(int argc, char** argv) {
	if (argc == 1) {
		std::cerr << "dsv_manip\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "use " << argv[0] << " <mode>\n";
		std::cerr << "Available modes:\n";
		std::cerr << "\tsort - sort (asc/desc) file by a given column\n";
		std::cerr << "\tlimit - limit number of rows to first highests/lowests values of a given column\n";
		exit(1);
	}
	return mode_from_string(argv[1]);
}

int main_sort(int argc, char** argv) {
	auto params = sort_mode::read_params(argc, argv);

	params.Print(std::cerr);

	sort_mode::run(params);
	return 0;
}

int main_limit(int argc, char** argv) {
	auto params = limit_mode::read_params(argc, argv);

	params.Print(std::cerr);

	limit_mode::run(params);
	return 0;
}

int main(int argc, char** argv) {
	std::cerr << "Welcome to dsv_manip (simple manipulation of dsv or delimiter-separated values file format\n";
	auto mode = get_mode_or_print_help(argc, argv);

	switch (mode)
	{
	case Mode::sort:
		return main_sort(argc, argv);
	case Mode::limit:
		return main_limit(argc, argv);
	default:
		std::cerr << "Error: should never be here";
		return 1;
	}
}
