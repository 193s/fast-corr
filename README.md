# fast-nonparametric-corr
[![C/C++ CI](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml)  
Currently works only on g++ with -O1/O2/O3 options.  

## APIs
| Class | time complexity of each operation |
| ---- | ---- |
| OnlineSpearman | O(logN) |
| OnlineSpearmanLinear | O(N) |
| OnlineKendall | O(logN) |

### OnlineSpearman\<T\>
- `void push_front(T x_val)`
- `void push_back(T x_val)`
- `void pop_front()`
- `void pop_back()`
- `double spearman_r()`

OnlineSpearmanLinear\<T\> may work faster on smaller `N`s.  
`T` can be `int`, `double`, ... or any other type with comparison operators defined.

### OnlineKendall\<T\>
- `void push_front(T x_val)`
- `void push_back(T x_val)`
- `void pop_front()`
- `void pop_back()`
- `double kendall_tau()`  


## Installation
`git clone https://github.com/193s/fast-nonparametric-corr/`
### testing:
`g++-14 test/basic_tests.cpp -O3`  

### benchmark:
`g++-14 test/benchmark.cpp -O3`  
<!--
- `./a.out r <<< "20000 1000"` : testing on randomized sequence without duplicate values, T=20000, N=1000
- `./a.out d <<< "20000 1000"` : testing on randomized sequence with duplicate values, T=20000, N=1000
-->

### running sample code:
`g++-14 sample.cpp -O3`  
```c++
#include <iostream>
#include "lib/spearman-algos.hpp"
#include "lib/kendall-algos.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  auto sp = OnlineSpearman<double>({0, 1, 1, 2, 2, 1});
  sp.push_back(321);
  cout<<"spearman([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << sp.spearman_r() << "\n";
  sp.pop_front();
  cout<<"spearman([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << sp.spearman_r() << "\n";

  auto kd = OnlineKendall<double>({0, 1, 1, 2, 2, 1});
  kd.push_back(321);
  cout<<"kendall([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << kd.kendall_tau() << "\n";
  kd.pop_front();
  cout<<"kendall([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << kd.kendall_tau() << "\n";
  return 0;
}
```

result:
```
spearman([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = 0.767193
spearman([1,1,2,2,1,321], [1,2,3,4,5,6]) = 0.617213

kendall([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = 0.688033
kendall([1,1,2,2,1,321], [1,2,3,4,5,6]) = 0.544949
```

verify this with Python's `scipy.stats`:
```python
from scipy import stats
>>> '%.6f' % stats.spearmanr([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]).statistic
'0.767193'
>>> '%.6f' % stats.spearmanr([1,1,2,2,1,321], [1,2,3,4,5,6]).statistic
'0.617213'

>>> '%.6f' % stats.kendalltau([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]).statistic
'0.688033'
>>> '%.6f' % stats.kendalltau([1,1,2,2,1,321], [1,2,3,4,5,6]).statistic
'0.544949'
```
