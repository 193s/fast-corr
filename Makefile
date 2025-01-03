GPP ?= g++-14
all: ./basic_tests ./benchmark ./sample
LIB=lib/rbst_base.hpp lib/lazy_rbst.hpp lib/fast_corr_base.hpp lib/spearman_algos.hpp lib/kendall_algos.hpp lib/pearson_algos.hpp

./basic_tests: $(LIB) test/basic_tests.cpp
	$(GPP) test/basic_tests.cpp -O3 -Wall -o basic_tests

./benchmark: $(LIB) test/basic_tests.cpp test/benchmark.cpp
	$(GPP) test/benchmark.cpp -O3 -Wall -o benchmark

./sample: $(LIB) sample.cpp
	$(GPP) sample.cpp -O3 -Wall -o sample

compile_test: $(LIB)
	$(GPP) test/compile_test1.cpp -Wall -Ofast -o compile_test
	./compile_test
	$(GPP) test/compile_test2.cpp -Wall -Ofast -o compile_test
	./compile_test
	$(GPP) test/compile_test3.cpp -Wall -Ofast -o compile_test
	./compile_test
	rm compile_test

clean:
	rm -f basic_tests benchmark sample
