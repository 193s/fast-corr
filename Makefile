GPP ?= g++-14
all: ./basic_tests.out ./benchmark.out ./sample.out
LIB=lib/spearman-algos.hpp lib/kendall-algos.hpp

./basic_tests.out: $(LIB) test/basic_tests.cpp
	$(GPP) test/basic_tests.cpp -O3 -o basic_tests.out

./benchmark.out: $(LIB) test/benchmark.cpp
	$(GPP) test/benchmark.cpp -O3 -o benchmark.out

./sample.out: $(LIB) sample.cpp
	$(GPP) sample.cpp -O3 -o sample.out

clean:
	rm -f *.out
