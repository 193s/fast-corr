GPP ?= g++-14
all: ./basic_tests ./benchmark ./sample
LIB=lib/spearman-algos.hpp lib/kendall-algos.hpp

./basic_tests: $(LIB) test/basic_tests.cpp
	$(GPP) test/basic_tests.cpp -O3 -o basic_tests

./benchmark: $(LIB) test/benchmark.cpp
	$(GPP) test/benchmark.cpp -O3 -o benchmark

./sample: $(LIB) sample.cpp
	$(GPP) sample.cpp -O3 -o sample

clean:
	rm -f basic_tests benchmark sample
