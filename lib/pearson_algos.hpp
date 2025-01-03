#pragma once
#include <cassert>
#include <vector>
#include <cmath>
#include "fast_corr_base.hpp"

namespace FastCorr::OfflineCorr {
  /**
   * O(N) straightforward offline algorithm (pearson-r)
   */
  double pearson_r(const std::vector<double> &X, const std::vector<double> &Y) {
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
}
namespace FastCorr::OnlineCorr {
  /**
   * straightforward pearson implementation: O(N) time complexity
   */
  //template< class TX, class TY >
  class Pearson : public Base<double, double> {
    public:
      int N = 0;
      double x_mean = 0, y_mean = 0;
      double cov = 0; //  cov  = sum[i=1..n] (x_i-x_mean)*(y_i-y_mean)
      double var_x = 0; // var_x = sum[i=1..n] (x_i - x_mean)^2
      double var_y = 0; // var_y = sum[i=1..n] (y_i - y_mean)^2
      void add(const double &x_val, const double &y_val) override {
        double x_new_mean = x_mean + (x_val - x_mean) / (N+1); // = (x_val + N*x_mean) / (N+1)
        double y_new_mean = y_mean + (y_val - y_mean) / (N+1);
        cov += (x_val - x_mean) * (y_val - y_new_mean);
        var_x += (x_val - x_mean) * (x_val - x_new_mean);
        var_y += (y_val - y_mean) * (y_val - y_new_mean);
        x_mean = x_new_mean, y_mean = y_new_mean;
        N += 1;
      }
      /**
       * note: this does not check if the pair (x_val, y_val) exists in the set
       */
      void remove(const double &x_val, const double &y_val) override {
        double x_new_mean = x_mean + (x_mean - x_val) / (N-1); // = (N*x_mean - x_val) / (N-1)
        double y_new_mean = y_mean + (y_mean - y_val) / (N-1);
        cov -= (x_val - x_mean) * (y_val - y_new_mean);
        var_x -= (x_val - x_mean) * (x_val - x_new_mean);
        var_y -= (y_val - y_mean) * (y_val - y_new_mean);
        x_mean = x_new_mean, y_mean = y_new_mean;
        N -= 1;
      }
      double pearson_r() const {
        if (N <= 1) return NAN;
        return cov / (sqrt(var_x)*sqrt(var_y));
      }
      virtual double r() const override { return pearson_r(); }
      virtual size_t size() const override { return N; }
  };
}
