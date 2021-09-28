// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include <string>

class RabinKarp {
    uint64_t kp_;
    uint64_t init_;
    uint64_t rm_;

public:
    uint64_t sigma_;
    uint64_t hash_;
    uint64_t size_;
    const std::basic_string<uint64_t> &ws_;

    RabinKarp(const std::basic_string<uint64_t> &s, uint64_t init, uint64_t size, uint64_t range, uint64_t sigma = 257);

    uint64_t hash() const { return hash_; };
    void next();
};

RabinKarp::RabinKarp(const std::basic_string<uint64_t> &s, uint64_t init, uint64_t size, uint64_t range, uint64_t sigma)
    : sigma_(sigma), size_(size), ws_(s), hash_(0), init_(init), rm_(1), kp_(range) {
    for (auto i = init; i < init + size_; ++i) {
        uint64_t next = ws_[i];
        next = next % kp_;
        hash_ = (sigma_ * hash_ + next) % kp_; // sigma or  little prime
    }

    for (int i = 0; i < size_ - 1; ++i)
        rm_ = (rm_ * sigma_) % kp_;
}

void RabinKarp::next() {
    uint64_t next = ws_[init_];
    next = next % kp_;
    hash_ = (hash_ + kp_ - rm_ * next % kp_) % kp_;
    init_++;
    next = ws_[init_ + size_ - 1];
    next = next % kp_;
    hash_ = (hash_ * sigma_ + next) % kp_;
}
