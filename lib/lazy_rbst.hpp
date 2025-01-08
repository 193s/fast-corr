#pragma once
#include "rbst_base.hpp"

namespace FastCorr {
  /**
   * Internal node class defined for Spearman algorithm
   */
  template <typename E>
  struct D2LazyRBSTNode {
    typename RBSTBase<D2LazyRBSTNode>::Ptr l, r;
    sp_d1_type key_d1, sum_d1;
    sp_d2_type key_d2, sum_d2;
    E lazy;
    int cnt;
    //D2LazyRBSTNode(const T &t = T(), const E &e = E())
    D2LazyRBSTNode(sp_d1_type d1 = sp_d1_type(), sp_d2_type d2 = sp_d2_type())
        : l(), r(), key_d1(d1), sum_d1(d1), key_d2(d2), sum_d2(d2), lazy(0), cnt(1) {}
  };

  typedef int E;
  /**
   * @brief Randomized binary search tree with lazy propagation, customized for spearman calculation
   */
  struct D2LazyRBST : RBSTBase<D2LazyRBSTNode<E>> {
    using Node = D2LazyRBSTNode<E>;
    using base = RBSTBase<D2LazyRBSTNode<E>>;
    using base::merge;
    using base::split;
    using typename base::Ptr;

    D2LazyRBST() = default;

    void apply(Ptr &t, int a, int b, const E &e) {
      auto x = split(t, a);
      auto y = split(x.second, b - a);
      propagate(y.first, e);
      t = merge(x.first, merge(y.first, y.second));
    }
    void propagate(Ptr t, const E &a) {
      t->lazy += a;
      // add(t->key, x, 1):
      t->key_d2 += (sp_d2_type)a*((sp_d2_type)2*t->key_d1 + a);
      t->key_d1 += a;
      // add(t->sum, x, t->cnt):
      t->sum_d2 += (sp_d2_type)a*((sp_d2_type)2*t->sum_d1 + a*(t->cnt)); // d2 += a^2*n + 2*a*d1
      t->sum_d1 += (sp_d2_type)a*(t->cnt); // d1 += a*n
    }
    inline sp_d2_type sum_d2(const Ptr t) const noexcept { return t ? t->sum_d2 : sp_d2_type(); }
    protected:
      Ptr update(Ptr t) override {
        push(t);
        t->cnt = 1;
        t->sum_d1 = t->key_d1;
        t->sum_d2 = t->key_d2;
        if (t->l) t->cnt += t->l->cnt, t->sum_d1 += t->l->sum_d1, t->sum_d2 += t->l->sum_d2;
        if (t->r) t->cnt += t->r->cnt, t->sum_d1 += t->r->sum_d1, t->sum_d2 += t->r->sum_d2;
        return t;
      }
      void push(Ptr t) override {
        if (t->lazy != E()) {
          if (t->l) propagate(t->l, t->lazy);
          if (t->r) propagate(t->r, t->lazy);
          t->lazy = E();
        }
      }
  };
}
