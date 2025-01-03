#pragma once
#include "rbst_base.hpp"

namespace FastCorr {
  /**
   * @brief Randomized binary search tree with lazy propergation
   */
  template <typename T, typename E>
  struct LazyRBSTNode {
    typename RBSTBase<LazyRBSTNode>::Ptr l, r;
    T key, sum;
    E lazy;
    int cnt;

    LazyRBSTNode(const T &t = T(), const E &e = E())
        : l(), r(), key(t), sum(t), lazy(e), cnt(1) {}
  };

  template <typename T, typename E, T (*f)(T, T), T (*g)(T, E), E (*h)(E, E),
            T (*ts)(T)>
  struct LazyRBST : RBSTBase<LazyRBSTNode<T, E>> {
    using Node = LazyRBSTNode<T, E>;
    using base = RBSTBase<LazyRBSTNode<T, E>>;
    using base::merge;
    using base::split;
    using typename base::Ptr;

    LazyRBST() = default;

    T fold(Ptr &t, int a, int b) {
      auto x = split(t, a);
      auto y = split(x.second, b - a);
      auto ret = sum(y.first);
      t = merge(x.first, merge(y.first, y.second));
      return ret;
    }

    void apply(Ptr &t, int a, int b, const E &e) {
      auto x = split(t, a);
      auto y = split(x.second, b - a);
      propagate(y.first, e);
      t = merge(x.first, merge(y.first, y.second));
    }

   protected:
    inline T sum(const Ptr t) const { return t ? t->sum : T(); }

    Ptr update(Ptr t) override {
      push(t);
      t->cnt = 1;
      t->sum = t->key;
      if (t->l) t->cnt += t->l->cnt, t->sum = f(t->l->sum, t->sum);
      if (t->r) t->cnt += t->r->cnt, t->sum = f(t->sum, t->r->sum);
      return t;
    }

    void push(Ptr t) override {
      if (t->lazy != E()) {
        if (t->l) propagate(t->l, t->lazy);
        if (t->r) propagate(t->r, t->lazy);
        t->lazy = E();
      }
    }

    void propagate(Ptr t, const E &x) {
      t->lazy = h(t->lazy, x);
      t->key = g(t->key, x);
      t->sum = g(t->sum, x);
    }
  };
}
