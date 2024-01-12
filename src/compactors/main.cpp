#include "kmer.h"
#include "engine.h"

#include "io.h"
#include "console.h"
#include "../common/version.h"

#include <iostream>
#include <string>
#include <iomanip>

using namespace std;


size_t verify(const string& interFile, const std::vector<string>& compactorsFiles, int compactorLen) {

	ifstream ifs;
	const int N = 16 << 20; // 16MB buffer
	char* buf = new char[N];
	ifs.rdbuf()->pubsetbuf(buf, N);

	ifs.open(interFile);
	unordered_map<string, int> hits;
	string line;

	while (getline(ifs, line)) {
		if ((int)line.length() < compactorLen) {
			continue;
		}
			
		++hits[line.substr(0, compactorLen)];
	}

	ifs.close();

	size_t bad = 0;
	for (const auto& file : compactorsFiles) {
		ifs.open(file);
		string line;

		vector<string> compactors;
		while (getline(ifs, line)) {
			compactors.push_back(line);
		}

		std::vector<int> counts(compactors.size());
		
		for (int i = 0; i < (int)compactors.size(); ++i) {
			auto it = hits.find(compactors[i]);

			if (it == hits.end()) {
				++bad;
			}
			else {
				counts[i] += it->second;
			}
		}
	}

	delete[] buf;
	
	return bad;
}

int main(int argc, char** argv) {

	try {
		cerr << "\nSPLASH: compactors\n";
		SPLASH_VER_PRINT(cerr);
		
		Log::getInstance(Log::LEVEL_NORMAL).enable();

		Console console;
		if (!console.parse(argc, argv)) {
			console.printUsage();
			return 0;
		}

		ofstream logStream;
		if (console.logFile.length()) {
			logStream.open(console.logFile);
			if (logStream) {
				Log::getInstance(Log::LEVEL_DEBUG).setOutStream(logStream);
				Log::getInstance(Log::LEVEL_DEBUG).enable();
				LOG_DEBUG << "Debug mode on" << endl;
			}
		}

		std::shared_ptr<IKmerProvider> loader;
		
		loader = make_shared<ReadLoader>(
			console.sampleFastqs,
			console.anchorsTsv,
			console.getParams().polyThreshold,
			console.numThreads,
			console.readsBufferGb,
			console.anchorsBatchSize,
			console.keepTemp);

		std::shared_ptr<ICompactorWriter> writer = make_shared<Output>(console.outputTsv, console.outputFasta);
		Engine engine(loader, writer, console.getParams());
		engine();
	}
	catch (std::runtime_error& err) {
		cerr << err.what();
	}

	return 0;
}
