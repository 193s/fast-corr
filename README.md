# fast-nonparametric-corr
[![C/C++ CI](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml)  
Currently works only on g++ with -O1/O2/O3 options.  

## testing:
`g++-14 test/basic_tests.cpp -O3`  

## benchmark:
`g++-14 test/benchmark.cpp -O3`  
<!--
- `./a.out r <<< "20000 1000"` : testing on randomized sequence without duplicate values, T=20000, N=1000
- `./a.out d <<< "20000 1000"` : testing on randomized sequence with duplicate values, T=20000, N=1000
-->

## running sample code:
`g++-14 sample.cpp -O3`  
```c++
#include <iostream>
#include <vector>
#include <limits>
#include "lib/spearman-algos.hpp"
#include "lib/kendall-algos.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  auto x = OnlineKendall<double>();
  x.push_back(0);
  x.push_back(1);
  x.push_back(1);
  x.push_back(2);
  x.push_back(2);
  x.push_back(1);
  x.push_back(321);
  cout << x.kendall_tau() << "\n";
  auto y = OnlineSpearman<double>();
  y.push_back(0);
  y.push_back(1);
  y.push_back(1);
  y.push_back(2);
  y.push_back(2);
  y.push_back(1);
  y.push_back(321);
  cout << y.spearman_r() << "\n";
  return 0;
}
```
