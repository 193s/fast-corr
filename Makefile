GPP ?= g++-14
all: ./basic_tests.out ./benchmark.out ./sample.out

./basic_tests.out:
	$(GPP) test/basic_tests.cpp -O3 -o basic_tests.out

./benchmark.out:
	$(GPP) test/benchmark.cpp -O3 -o benchmark.out

./sample.out:
	$(GPP) sample.cpp -O3 -o sample.out

clean:
	rm -f *.out
