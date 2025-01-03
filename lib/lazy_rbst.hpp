#pragma once
#include "rbst_base.hpp"

namespace FastCorr {
  /**
   * Internal node class defined for Spearman algorithm
   */
  struct NodeVal {
    sp_d1_type d1;
    sp_d2_type d2;
    int cnt;
    NodeVal() { d1 = d2 = 0, cnt = 1; }
    NodeVal(sp_d1_type x1, sp_d2_type x2, int c) { d1 = x1, d2 = x2, cnt = c; }
  };
  inline NodeVal f(const NodeVal &x, const NodeVal &y) {
    return NodeVal(x.d1 + y.d1, x.d2 + y.d2, x.cnt + y.cnt);
  }
  inline void g(NodeVal &x, int a) { // this modifies the content of x
    if (a == 0) return;
    int size = x.cnt;
    x.d2 += (sp_d2_type)a*(a*size + 2*x.d1); // d2 += a*a*size + 2*a*d1
    x.d1 += (sp_d1_type)a*size;
  }
  typedef int E;
  typedef NodeVal T;

  /**
   * @brief Randomized binary search tree with lazy propergation, customized for spearman calculation
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

  //template <typename T, typename E, T (*f)(T, T), T (*g)(T, E), E (*h)(E, E)>
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
      t->lazy += x;
      g(t->key, x);
      g(t->sum, x);
    }
  };
}
