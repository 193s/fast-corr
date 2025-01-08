C_COMPILER ?= g++-14 --std=c++11
all: ./basic_tests ./benchmark ./sample
LIB=lib/rbst_base.hpp lib/lazy_rbst.hpp lib/fast_corr_base.hpp lib/spearman_algos.hpp lib/kendall_algos.hpp lib/pearson_algos.hpp

./basic_tests: $(LIB) test/basic_tests.cpp test/test_base.hpp
	$(C_COMPILER) test/basic_tests.cpp -O3 -Wall -o basic_tests

./benchmark: $(LIB) test/benchmark.cpp test/test_base.hpp
	$(C_COMPILER) test/benchmark.cpp -O3 -Wall -o benchmark

./sample: $(LIB) sample.cpp
	$(C_COMPILER) sample.cpp -O3 -Wall -o sample

compile_test: $(LIB)
	$(C_COMPILER) test/compile_test1.cpp -Wall -Ofast -o compile_test
	./compile_test
	$(C_COMPILER) test/compile_test2.cpp -Wall -Ofast -o compile_test
	./compile_test
	$(C_COMPILER) test/compile_test3.cpp -Wall -Ofast -o compile_test
	./compile_test
	rm compile_test

clean:
	rm -f basic_tests benchmark sample
