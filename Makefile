all: ./a.out
a.out:
	/usr/local/bin/g++-14 test.cpp -O3
clean:
	rm a.out
run:
	./a.out < in.txt
