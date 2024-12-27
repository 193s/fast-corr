all: ./basic_tests.out ./benchmark.out ./sample.out

./basic_tests.out:
	/usr/local/bin/g++-14 test/basic_tests.cpp -O3 -o basic_tests.out

./benchmark.out:
	/usr/local/bin/g++-14 test/benchmark.cpp -O3 -o benchmark.out

./sample.out:
	/usr/local/bin/g++-14 sample.cpp -O3 -o sample.out

clean:
	rm -f *.out
