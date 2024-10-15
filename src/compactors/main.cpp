#include "kmer.h"
#include "engine.h"

#include "io.h"
#include "console.h"
#include "../common/version.h"
#include "../../libs/refresh/parallel_queues/lib/parallel-queues.h"

#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>
#include <thread>
#include <filesystem>
#include <numeric>
#include <mutex>

using namespace std;

/*
void testExtenders(){

	string sequence{ "AAAAAAAGCCTCCGGAACTTCTTCAAGTCCGGAGGCTTCAATCAAAACGGCGGCTACCTACTCTCCCACTGTGACGCAGTA" };
	int k = 9;
	int L = sequence.size();
	int n = L / k;


	std::vector<kmer_t> kmers(n);

	for (int i = 0; i < n; ++i) {
		KmerHelper::from_string(sequence.c_str() + i * k, k, kmers[i]);
	}
	
	Compactor c0(&kmers[0], k, 0);
	Compactor c1(c0, &kmers[1], k, n - 1, 0, 0);

	char buffer[1024];

	c0.to_string(buffer, true);
	cout << buffer << endl;

	c1.to_string(buffer, true);
	cout << buffer << endl;

	// get extenders
	for (int len = 1; len < 30; ++len) {
		kmer_t e = c1.get_extender(len);
		cout << string(L - len, '-');
		string s = KmerHelper::to_string(e, len);
		cout << s << endl;
	}

	c1.to_string(buffer, true);
	cout << buffer << endl;

	for (int len : {5, 9, 11}) {
		// get extenders
		for (int shift = 0; shift < 30; ++shift) {
			kmer_t e = c1.get_extender(len, shift);
			cout << string(L - len - shift, '-');
			string s = KmerHelper::to_string(e, len);
			cout << s << string(shift, '-') << endl;
			c1.extender_shift = shift;
			c1.to_string(buffer, true);
			cout << buffer << endl;
		}

//		c1.to_string(buffer);
//		cout << buffer << endl;
	}

	return;
}



void testNChooseK() {

	NChooseK<uint64_t, 10000> n_choose_k;
	NChooseK_log<10000> n_choose_k_log;

	std::vector<std::pair<int, int>> nks{{10, 6}, { 20, 11 }, {30, 27}, 
		{ 100, 95 }, { 1000,995 }, { 1000,776 }, { 10000,9997 }};
	
	int n_times = 10000;

	for (const auto& nk : nks) {

		double old;
		double nw;
		
		auto start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < n_times; ++i) {
			old = (double)n_choose_k(nk.first, nk.second);
		}
		std::chrono::duration<double> dt_old = std::chrono::high_resolution_clock::now() - start;
		
		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < n_times; ++i) {
			nw = n_choose_k_log(nk.first, nk.second);
		}
		std::chrono::duration<double> dt_new = std::chrono::high_resolution_clock::now() - start;

		cout << nk.first << " choose " << nk.second << endl
			<< "  Old: " << old << ", t = " << dt_old.count() << endl
			<< "  New: " << nw << ", t = " << dt_new.count() << endl;
	}
}
*/


int main(int argc, char** argv) {

	//testNChooseK();
	//testExtenders();
	
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

		if (console.independentOutputs) {
			LOG_NORMAL << "Running compactors independently on " << console.sampleFastqs.size() << " files" << endl;
			
			refresh::parallel_queue<int> queue(console.sampleFastqs.size(), 1, "fastq-queue");
			
			for (int i = 0; i < console.sampleFastqs.size(); ++i) {
				queue.push(std::move(i));
			}
			queue.mark_completed();

			mutex mtx;

			Log::getInstance(Log::LEVEL_NORMAL).disable();

			vector<thread> threads(console.numThreads);
			for (int tid = 0; tid < threads.size(); ++tid) {
				threads[tid] = thread([&console, &queue, &mtx, tid]() {
					int i;
					while (queue.pop(i)) {

						string inFastq{ console.sampleFastqs[i] };

						string inFastqName = std::filesystem::path(inFastq).filename().string();
						string outTsvName = std::filesystem::path(console.outputTsv).filename().string();
						string outTsv = std::filesystem::path(console.outputTsv).replace_filename(inFastqName + '-' + outTsvName).string();
						
						string outFasta;
						if (console.outputFasta.length()) {
							string outFastaName = std::filesystem::path(console.outputFasta).filename().string();
							outFasta = std::filesystem::path(console.outputFasta).replace_filename(inFastqName + '-' + outFastaName).string();
						}

						std::shared_ptr<ICompactorWriter> writer = make_shared<Output>(outTsv, outFasta, console.noSubcompactors, console.cumulatedStats);

						mtx.lock();
						cerr << "Started " << i + 1 << ": " << inFastqName << endl;
						mtx.unlock();

						std::shared_ptr<IKmerProvider> loader;

						loader = make_shared<ReadLoader>(
							vector<string>(1, inFastq),
							console.inputFormat,
							console.anchorsTsv,
							console.getParams().polyThreshold,
							1,
							console.readsBufferGb,
							console.anchorsBatchSize,
							console.keepTemp,
							"./tmp-compactors-" + to_string(tid));

					
						Engine engine(loader, writer, console.getParams());
						engine();
						auto n_compactors = engine.getOutputCompactors().back().id + 1;
						
						mtx.lock();
						cerr << "Finished " << i + 1 << ": " << inFastqName << " with " << n_compactors << " final compactors" << endl;
						mtx.unlock();
					}
					
				});
			}

			for (auto& t : threads) {
				t.join();
			}

			Log::getInstance(Log::LEVEL_NORMAL).enable();

		}
		else {

			std::shared_ptr<IKmerProvider> loader;

			loader = make_shared<ReadLoader>(
				console.sampleFastqs,
				console.inputFormat,
				console.anchorsTsv,
				console.getParams().polyThreshold,
				console.numThreads,
				console.readsBufferGb,
				console.anchorsBatchSize,
				console.keepTemp);

			std::shared_ptr<ICompactorWriter> writer = make_shared<Output>(console.outputTsv, console.outputFasta, console.noSubcompactors, console.cumulatedStats);
			Engine engine(loader, writer, console.getParams());
			engine();

		}
	}
	catch (std::runtime_error& err) {
		cerr << err.what();
	}

	return 0;
}
