#pragma once
#include <cstdint>

namespace FastCorr {
  /**
   * @brief RBSTBase<Node>: A base class for randomized (self-balancing) binary search trees
   */
  template <typename Node>
  struct RBSTBase {
    using Ptr = Node *;
    template <typename... Args>
    inline Ptr my_new(Args... args) {
      return new Node(args...);
    }
    inline void my_del(Ptr t) { delete t; }
    inline Ptr make_tree() const { return nullptr; }

    // for avoiding memory leak, activate below
    /*
    using Ptr = shared_ptr<Node>;
    template <typename... Args>
    inline Ptr my_new(Args... args) {
      return make_shared<Node>(args...);
    }
    inline void my_del(Ptr t) {}
    Ptr make_tree() {return Ptr();}
    */

    inline int count(const Ptr t) const noexcept { return t ? t->cnt : 0; }
    inline int size (const Ptr t) const noexcept { return count(t); }

    Ptr merge(Ptr l, Ptr r) {
      if (!l) return r;
      if (!r) return l;
      if (int((rng() * (l->cnt + r->cnt)) >> 32) < l->cnt) {
      //if ((int)(rng() % (uint64_t)(l->cnt + r->cnt)) < l->cnt) {
        push(l);
        l->r = merge(l->r, r);
        return update(l);
      }
      else {
        push(r);
        r->l = merge(l, r->l);
        return update(r);
      }
    }

    std::pair<Ptr, Ptr> split(Ptr t, int k) {
      if (!t) return {nullptr, nullptr};
      push(t);
      if (k <= count(t->l)) {
        auto s = split(t->l, k);
        t->l = s.second;
        return {s.first, update(t)};
      }
      else {
        auto s = split(t->r, k - count(t->l) - 1);
        t->r = s.first;
        return {update(t), s.second};
      }
    }

    // this is equivalent to split(t, 1) but a bit faster
    std::pair<Ptr, Ptr> split_by_first_element(Ptr t) {
      if (FAST_CORR_UNLIKELY(!t)) return {nullptr, nullptr};
      push(t);
      if (FAST_CORR_LIKELY((t->l) != nullptr)) {
        auto s = split_by_first_element(t->l);
        t->l = s.second;
        return {s.first, update(t)};
      }
      else { // t->l == NULL
        if (t->r) {
          auto s = t->r;
          t->r = nullptr;
          return {update(t), s};
        }
        else {
          return {t, nullptr};
        }
      }
    }

    std::pair<Ptr, std::pair<Ptr, Ptr> > split3(Ptr t, int k1, int k2) {
      // split (a|b)|c or a|(b|c) : len(a) = k1, len(b) = k2-k1, len(c) = n-k2
      if (k1 <= count(t)-k2) { // choose the order of two splits to minimize time cost
        // len(a) <= len(c): split by k2 first
        auto s = split(t, k2);
        auto s2 = split(s.first, k1);
        return {s2.first, {s2.second, s.second}};
      }
      else {
        auto s = split(t, k1);
        auto s2 = split(s.second, k2-k1);
        return {s.first, s2};
      }
    }

    Ptr merge3(Ptr t1, Ptr t2, Ptr t3) {
      return merge(merge(t1, t2), t3);
    }
    Ptr merge4(Ptr t1, Ptr t2, Ptr t3, Ptr t4) {
      return merge(merge(t1, t2), merge(t3, t4));
    }

    Ptr build(int l, int r, const std::vector<decltype(Node::key_d1)> &v) {
      if (l + 1 == r) return my_new(v[l]);
      int m = (l + r) >> 1;
      Ptr pm = my_new(v[m]);
      if (l < m) pm->l = build(l, m, v);
      if (m + 1 < r) pm->r = build(m + 1, r, v);
      return update(pm);
    }

    Ptr build(const std::vector<decltype(Node::key_d1)> &v) {
      return build(0, (int)v.size(), v);
    }

    template <typename... Args>
    void insert(Ptr &t, int k, const Args &... args) {
      auto x = split(t, k);
      t = merge(merge(x.first, my_new(args...)), x.second);
    }

    void erase(Ptr &t, int k) {
      auto x = split(t, k);
      auto y = split(x.second, 1);
      my_del(y.first);
      t = merge(x.first, y.second);
    }

   protected:
    static uint64_t rng() {
      static uint64_t x_ = 88172645463325252ULL;
      return x_ ^= x_ << 7, x_ ^= x_ >> 9, x_ & 0xFFFFFFFFull;
    }

    virtual void push(Ptr) = 0;
    virtual Ptr update(Ptr) = 0;
  };
}
