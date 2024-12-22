# fast-nonparametric-corr
`/usr/local/bin/g++-14 test.cpp -O3`

testing:

- `./a.out r <<< "20000 1000"` : testing on randomized sequence without duplicate values, T=20000, N=1000
- `./a.out d <<< "20000 1000"` : testing on randomized sequence with duplicate values, T=20000, N=1000
