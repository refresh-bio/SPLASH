#include <iostream>
#include <string>
#include <limits>
#include <sstream>
#include "../common/version.h"

#include "build.h"
#include "query.h"
#include "query_many.h"
#include "lookup.h"


enum class Mode { build, query, query_many };

std::string to_string(Mode mode) {
	switch (mode)
	{
	case Mode::build:
		return "build";
	case Mode::query:
		return "query";
	case Mode::query_many:
		return "query_many";
	default:
		std::cerr << "Error: unhandled mode " << (int) mode << "\n";
		exit(1);
		break;
	}
}

Mode mode_from_string(const std::string& mode) {
	if (mode == "build") {
		return Mode::build;
	}
	else if (mode == "query") {
		return Mode::query;
	}
	else if (mode == "query_many") {
		return Mode::query_many;
	}
	std::cerr << "Error: unknown mode " << mode << "\n";
	exit(1);
}

Mode get_mode_or_print_help(int argc, char** argv) {
	if (argc == 1) {
		std::cerr << "lookup_table " << Lookup::GetVersion() << "\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "use " << argv[0] << " <mode>\n";
		std::cerr << "Available modes:\n";
		std::cerr << "\tbuild - build lookup table\n";
		std::cerr << "\tquery - query lookup table\n";
		std::cerr << "\tquery_many - query single lookup table with many queries (to avoid loading index each run)\n";
		exit(1);
	}
	return mode_from_string(argv[1]);
}

int main_build(int argc, char** argv) {
	auto params = build_mode::read_params(argc, argv);

	params.Print(std::cerr);

	build_mode::run(params);
	return 0;
}

int main_query(int argc, char** argv) {
	auto params = query_mode::read_params(argc, argv);

	params.Print(std::cerr);

	query_mode::run(params);
	return 0;
}

int main_query_many(int argc, char** argv) {
	auto params = query_many_mode::read_params(argc, argv);

	params.Print(std::cerr);

	query_many_mode::run(params);
	return 0;
}

int main(int argc, char** argv) {
	std::cerr << "Welcome to lookup_table\n";
	auto mode = get_mode_or_print_help(argc, argv);

	switch (mode)
	{
	case Mode::build:
		return main_build(argc, argv);
	case Mode::query:
		return main_query(argc, argv);
	case Mode::query_many:
		return main_query_many(argc, argv);
	default:
		std::cerr << "Error: should never be here";
		return 1;
	}
}
