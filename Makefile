GPP ?= g++-14
all: ./basic_tests ./benchmark ./sample
LIB=lib/rbst_base.hpp lib/lazy_reversible_rbst.hpp lib/fast_corr_base.hpp lib/spearman_algos.hpp lib/kendall_algos.hpp

./basic_tests: $(LIB) test/basic_tests.cpp
	$(GPP) test/basic_tests.cpp -O3 -o basic_tests

./benchmark: $(LIB) test/benchmark.cpp
	$(GPP) test/benchmark.cpp -O3 -o benchmark

./sample: $(LIB) sample.cpp
	$(GPP) sample.cpp -O3 -o sample

clean:
	rm -f basic_tests benchmark sample
