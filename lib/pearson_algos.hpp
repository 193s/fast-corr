#pragma once
#include <vector>

namespace FastCorr {
  /**
   * straightforward pearson implementation: O(N) time complexity
   */
  double straightforward_pearson_r(const std::vector<double> &X, const std::vector<double> &Y) {
    assert(X.size() == Y.size());
    int n = X.size();
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double sum_X2 = 0, sum_Y2 = 0;

    for (int i=0; i<n; i++) {
        sum_X += X[i];
        sum_Y += Y[i];
        sum_XY += X[i]*Y[i];
        sum_X2 += X[i]*X[i];
        sum_Y2 += Y[i]*Y[i];
    }
    return (double)(n*sum_XY - sum_X*sum_Y)
            / sqrt((n*sum_X2 - sum_X*sum_X) * (n*sum_Y2 - sum_Y*sum_Y));
  }

  // TODO
  namespace OnlineCorr {
    template< class TX, class TY >
    class Pearson : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual double r() const = 0;
        virtual size_t size() const = 0;
    };
  }
}
