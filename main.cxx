/**
 * main.cxx
 *
 * A Series of Randmoness Tests for Binary Sequence Validation
 * by snovvcrash
 * 04.2017
 */

/*
	$ LD_LIBRARY_PATH=/home/username/gsl/lib/:${LD_LIBRARY_PATH}
	$ export LD_LIBRARY_PATH
	$ g++ -I/home/username/gsl/include/ -Wall -c -std=c++11 main.cxx
	$ g++ -L/home/username/gsl/lib/ main.o -lgsl -lgslcblas -lm
	$ ./a.out
	OR
	$ make GSL_INCLUDE='-I/home/username/gsl/include/' GSL_LIBRARY='-L/home/username/gsl/lib/ main.o -lgsl -lgslcblas -lm'
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sys/types.h> // S_ISREG
#include <sys/stat.h>  // struct stat
#include "randomness_test.hxx"

using std::cout;
using std::endl;
using std::cerr;

#define ERROR_IFILE_OPEN      ( -1)
#define ERROR_OPTION_NUMBER   ( -2)
#define ERROR_IS_REGULAR_FILE ( -3)

int is_regular_file(char const* path);
int prepare_input_file(char const* ifilename, std::ifstream& ifile);

int main(int argc, char* argv[]) {
	char* ifilename;

	if (argc == 2) ifilename = argv[1];
	else {
		cerr << "main: Invalid number of options" << endl;
		return ERROR_OPTION_NUMBER;
	}

	std::ifstream ifile;
	if (int errcode = prepare_input_file(ifilename, ifile))
		return errcode;

	std::vector<bool> seq;
	char c;
	while (ifile.get(c))
		if (c == '0')
			seq.push_back(0);
		else if (c == '1')
			seq.push_back(1);
	ifile.close();

	// M = block size
	// N = number of blocks
	// m = pattern size
	// n = length of sequence

	statistics::tests test(seq);

	test.monobit_frequency();
	test.block_frequency(512);     // M >= 20 and M > 0.01*n and N < 100
	test.runs();
	test.longest_run_of_ones(512); // M = 8 or 128 or 512 or 10000
	test.spectral();
	test.serial(2);                // 0 < m < (int)(log2(n))-2
	test.approximate_entropy(2);   // 0 < m < (int)(log2(n))-5
	test.cumulative_sums();

	test.print_verdict();

	return 0;
}

int is_regular_file(char const* path) {
	struct stat s;
	stat(path, &s);
	return S_ISREG(s.st_mode);
}

int prepare_input_file(char const* ifilename, std::ifstream& ifile) {
	if (!is_regular_file(ifilename)) {
		cerr << "main: No such file or it is not a regular file" << endl;
		return ERROR_IS_REGULAR_FILE;
	}

	ifile.open(ifilename, std::ios::binary);

	if (!ifile.is_open()) {
		cerr << "main: " << strerror(errno) << endl;
		return ERROR_IFILE_OPEN;
	}

	return 0;
}
