// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include "Block.h"
#include "LeafBlock.h"
#include "BackBlock.h"

class InternalBlock final : public Block {
    bool use_la;
    BlockLAType cached_la;

    void collapse_subtree(Block *block, Block *parent) {
        for (const auto &c : block->children_) {
            c->parent_ = parent;
            collapse_subtree(c, parent);
        }
    }

public:

    InternalBlock(Block *, size_t, size_t, const std::basic_string<uint64_t> &);
    ~InternalBlock() final;

    std::vector<Block *> &children(int, int) final;
    void clean_unnecessary_expansions() final;

    bool is_leaf() const final { return false; }
    bool is_la_leaf() const final { return use_la; };

    size_t la_visit(DataPtr data, uint8_t log_n, uint8_t log_u, uint8_t depth, uint8_t leaf_approx_cost) final {
        auto tree_bit_cost = children_.size(); // 1 bit for each child
        tree_bit_cost += log_u * (children_.size() - 1); // for the samples needed to implement rank
        for (auto rit = children_.rbegin(); rit != children_.rend(); ++rit)
            tree_bit_cost += (*rit)->la_visit(data, log_n, log_u, depth + 1, leaf_approx_cost);

        auto pruned_bit_cost = std::numeric_limits<size_t>::max();
        if (get_la(data, cached_la)) {
            auto back_block_cost = 2 * log_n + log_u; // An LA block is encoded as a back block at the same level
            pruned_bit_cost = cached_la.size_in_bytes() * 8 + back_block_cost;
        }

        if (pruned_bit_cost < tree_bit_cost) {
            use_la = true;
            collapse_subtree(this, this);
            this->parent_ = this;
            return pruned_bit_cost;
        }

        use_la = false;
        return tree_bit_cost;
    };

    std::pair<uint8_t, long double> compute_depth(uint8_t depth) final {
        if (is_la_leaf())
            return {depth, depth * actual_length()};

        size_t max_depth = depth;
        long double depth_sum = 0;
        for (auto *&c : children_) {
            auto[subtree_max_depth, subtree_depth_sum] = c->compute_depth(depth + 1);
            max_depth = (uint8_t) std::max<size_t>(max_depth, subtree_max_depth);
            depth_sum += subtree_depth_sum;
        }
        return {max_depth, depth_sum};
    }


    size_t count_leaf_elements() final {
        size_t sum = 0;
        for (auto *&c : children_)
            sum += c->count_leaf_elements();
        return sum;
    }

    bool get_la(DataPtr data, BlockLAType &la) final {
        if (use_la) {
            la = cached_la;
            return true;
        }
        return Block::get_la(data, la);
    }
};

InternalBlock::InternalBlock(Block *parent,
                             size_t start_index,
                             size_t end_index,
                             const std::basic_string<uint64_t> &source)
    : Block(parent, start_index, end_index, source),
      use_la(false),
      cached_la() {
}

InternalBlock::~InternalBlock() {
    for (auto rit = children_.rbegin(); rit != children_.rend(); ++rit)
        delete *rit;
}

std::vector<Block *> &InternalBlock::children(int leaf_length, int r) {
    if (children_.empty()) {
        auto next_length = length() / r;
        if (next_length <= leaf_length) {
            for (auto i = 0; i < r; ++i) {
                auto init = start_ + i * next_length;
                auto end = start_ + (i + 1) * next_length - 1;
                if (init < source_.size()) {
                    Block *child = new LeafBlock(this, init, end, source_);
                    children_.push_back(child);
                }
            }
        } else {
            for (auto i = 0; i < r; ++i) {
                auto init = start_ + i * next_length;
                auto end = start_ + (i + 1) * next_length - 1;
                if (init < source_.size()) {
                    Block *child = new InternalBlock(this, init, end, source_);
                    children_.push_back(child);
                }
            }
        }
    }
    return children_;
}

void InternalBlock::clean_unnecessary_expansions() {
    for (auto rit = children_.rbegin(); rit != children_.rend(); ++rit) {
        Block *b = (*rit);
        b->clean_unnecessary_expansions();
    }

    bool all_children_leaves = true;
    for (Block *child : children_)
        all_children_leaves = all_children_leaves && child->is_leaf();

    if (all_children_leaves && pointing_to_me_ == 0 && first_block_->start_ < start_ && second_block_ != this) {
        auto *bb = new BackBlock(parent_, start_, end_, source_, first_block_, second_block_, offset_);
        bb->level_index_ = level_index_;
        bb->first_occurrence_level_index_ = first_occurrence_level_index_;
        bb->left_ = true;
        bb->right_ = true;
        parent_->replace_child(this, bb);
        delete this;
    } else { //To avoid dangling references
        first_block_ = this;
        second_block_ = nullptr;
    }
}