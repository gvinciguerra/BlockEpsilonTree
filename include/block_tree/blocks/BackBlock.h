// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include "Block.h"

class BackBlock final : public Block {
public:

    BackBlock(Block *parent,
              size_t start_index,
              size_t end_index,
              const std::basic_string<uint64_t> &source,
              Block *first_block,
              Block *second_block,
              uint32_t offset)
        : Block(parent, start_index, end_index, source) {
        first_block_ = first_block;
        if (second_block) {
            if (second_block->start_ == start_index && second_block->end_ == end_index) second_block_ = this;
            else second_block_ = second_block;
        }
        offset_ = offset;
        if (first_block_)
            first_block_->pointing_to_me_++;
        if (second_block_)
            second_block_->pointing_to_me_++;
    }

    ~BackBlock() final {
        if (first_block_ && first_block_->pointing_to_me_ > 0)
            first_block_->pointing_to_me_--;
        if (second_block_ && second_block_->pointing_to_me_ > 0)
            second_block_->pointing_to_me_--;
    }

    size_t la_visit(DataPtr, uint8_t log_n, uint8_t log_u, uint8_t depth, uint8_t) final {
        return 2 * log_n + log_u;
    }

    virtual std::pair<uint8_t, long double> compute_depth(uint8_t depth) {
        auto left_pruned = first_block_->parent_->is_la_leaf();
        auto right_pruned = (second_block_ && second_block_->parent_->is_la_leaf()) || (left_pruned && !second_block_);

        auto left_length = length() - offset_;
        auto right_length = offset_;
        if (left_pruned && right_pruned)
            return {depth, depth * actual_length()};
        else if (left_pruned) {
            auto right_cost = second_block_ ? second_block_->compute_depth(depth).second : 0;
            return {depth, depth * left_length + right_cost * right_length / double(length())};
        } else if (right_pruned) {
            auto left_cost = first_block_->compute_depth(depth).second;
            return {depth, depth * right_length + left_cost * left_length / double(length())};
        } else {
            auto left_cost = first_block_->compute_depth(depth).second;
            auto right_cost = second_block_ ? second_block_->compute_depth(depth).second : 0;
            return {depth, left_cost * left_length / double(length()) + right_cost * right_length / double(length())};
        }
    }
};


