# randomness-test
A series of randomness tests for binary sequence validation based on a "Statistical Test Suite for Random and Pseudorandom Number Generators for Cryptographic Applications (Special Publication 800-22, Revision 1a)" by NIST.

Note: this implementation uses *GNU Scientific Library (GSL)*, so the latter must be installed and configured along with the *LD_LIBRARY_PATH* variable.

## How-to-use:

### Build and run

* Compile
  ```
  $ make GSL_INCLUDE='-I/home/username/gsl/include/' GSL_LIBRARY='-L/home/username/gsl/lib/ main.o -lgsl -lgslcblas -lm'
  ```
  where `/home/username/gsl/include/` and `/home/username/gsl/lib/` are paths to your GSL's *include* and *lib* folders respectively
  
* Execute
  ```
  $ ./randomness_test somefile.txt
  ```
  
### Testing

This project includes a *test-files* folder which contains some binary sequences to be tested. To do that one can use *run-test.sh* bash scrip the following way:

  ```
  $ ./run-test.sh ./randomness_test test-files/TrueRandom
  ```
