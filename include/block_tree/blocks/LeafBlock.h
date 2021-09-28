// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include "Block.h"

class LeafBlock final : public Block {
public:

    LeafBlock(Block *parent, size_t start_index, size_t end_index, const std::basic_string<uint64_t> &source)
        : Block(parent, start_index, end_index, source) {
    }

    ~LeafBlock() final = default;

    size_t la_visit(DataPtr, uint8_t, uint8_t log_u, uint8_t depth, uint8_t leaf_approx_cost) final {
        return leaf_approx_cost * actual_length();
    };

    virtual std::pair<uint8_t, long double> compute_depth(uint8_t depth) {
        return {depth, depth * actual_length()};
    }

    size_t count_leaf_elements() final { return actual_length(); }
};

