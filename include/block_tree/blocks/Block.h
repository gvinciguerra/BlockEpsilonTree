// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include "LABlock.hpp"
#include <string>
#include <vector>
#include <unordered_map>

using BlockLAType = LABlock<uint32_t>;
using DataPtr = const uint32_t *;
static uint8_t max_bpc = 14;

double ef_cost(size_t n, size_t u) {
    auto log_n = sdsl::bits::hi(n) + 1;
    auto log_u = sdsl::bits::hi(u) + 1;
    if (log_n == log_u)
        --log_n; // to ensure log_u-log_n > 0
    return n * (log_u - log_n) + n + (1ULL << log_n) + 0.2 * n;
//        return size_t(n * (2 + std::ceil(std::log2(u / double(n)))));
}

class Block {
public:
    Block *parent_;
    size_t start_;
    size_t end_;

    const std::basic_string<uint64_t> &source_;

    Block *first_block_;
    Block *second_block_;
    uint32_t offset_;
    bool left_: 1;
    bool right_: 1;
    uint32_t pointing_to_me_;
    uint8_t level_index_;
    uint32_t first_occurrence_level_index_;

    uint32_t la_block_id_;
    uint32_t level_start_pos_;

    std::vector<Block *> children_;

    size_t length() const { return end_ - start_ + 1; };
    size_t actual_length() const { return std::min<size_t>(end_ + 1, source_.size()) - start_; };
    std::basic_string<uint64_t> represented_string() const { return source_.substr(start_, length()); };

    virtual size_t count_leaf_elements() { return 0; };

    virtual std::vector<Block *> &children(int, int) { return children_; };
    virtual void clean_unnecessary_expansions() {};

    void replace_child(Block *old_child, Block *new_child) {
        for (auto &c : children_) {
            if (c == old_child) {
                c = new_child;
                return;
            }
        }
    }

    virtual bool is_leaf() const { return true; };
    virtual bool is_la_leaf() const { return false; };

    virtual size_t la_visit(DataPtr data, uint8_t log_n, uint8_t log_u, uint8_t depth, uint8_t leaf_approx_cost) = 0;

    virtual bool get_la(DataPtr data, BlockLAType &la) {
        // TODO: replace with exponential search from bpc=0 (skipping 1)
        for (uint8_t bpc = 0; bpc < max_bpc; bpc += 1 + (bpc == 0))
            if (BlockLAType::make(data + start_, data + std::min<size_t>(end_ + 1, source_.size()), bpc, la))
                return true;
        return false;
    }

    virtual std::pair<uint8_t, long double> compute_depth(uint8_t depth) = 0;

    virtual uint32_t number_la_blocks(uint32_t id) {
        if (is_la_leaf()) {
            la_block_id_ = id;
            return id + 1;
        }
        if (is_leaf())
            return id;

        for (auto *&c : children_)
            id = c->number_la_blocks(id);
        return id;
    }

    Block(Block *parent, size_t start_index, size_t end_index, const std::basic_string<uint64_t> &source)
        : parent_(parent),
          start_(start_index),
          end_(end_index),
          source_(source),
          offset_(0),
          left_(false),
          right_(false),
          first_block_(this),
          second_block_(nullptr),
          pointing_to_me_(0),
          level_index_(0),
          first_occurrence_level_index_(0),
          la_block_id_(0),
          level_start_pos_(0) {}

    virtual ~Block() = default;
};