#include "app.h"

using namespace std;

// ******************************************************************************************
int main(int argc, char** argv)
{
	CApplication app;

	if (!app.parse_args(argc, argv))
	{
		app.usage();
		return -1;
	}

	app.run();

	return 0;
}