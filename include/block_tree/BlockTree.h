// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include "RabinKarp.h"
#include "HashString.h"
#include "blocks/Block.h"
#include "blocks/BackBlock.h"
#include "blocks/InternalBlock.h"
#include "blocks/LeafBlock.h"

#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>

class BlockTree {
    using block_map = std::unordered_map<HashString, std::vector<Block *>>;
    using block_pairs_map = std::unordered_map<HashString, std::vector<std::pair<Block *, Block *>>>;

    void block_scan(std::vector<Block *> &, int, block_map &);
public:
    uint8_t r_;
    size_t leaf_length_;
    const std::basic_string<uint64_t> &input_;
    Block *root_block_;

    BlockTree(const std::basic_string<uint64_t> &, uint8_t, size_t);
    ~BlockTree();

    void process_back_pointers(size_t starting_block_size);
    void clean_unnecessary_expansions();

    void process_level(std::vector<Block *> &);

    void forward_window_block_scan(std::vector<Block *> &level,
                                   uint32_t window_size,
                                   uint32_t N,
                                   block_map &hashtable);
    void forward_pair_window_block_scan(
        std::vector<Block *> &level,
        uint32_t pair_window_size,
        uint32_t N,
        block_pairs_map &pair_hashtable);

    std::vector<Block *> next_level(std::vector<Block *> &level, bool skip_la_leaf) const;
    // Returns a vector of levels of nodes of the tree where
    // each level is represented by a vector of its nodes (left-to-right).
    //
    // A simple levelwise (left-to-right) traversal of the tree would be:
    //     for (std::vector<Block*> level : bt->levelwise_iterator()) {
    //         for (Block* b : level) {
    //             ...
    std::vector<std::vector<Block *>> levelwise_iterator();
};

BlockTree::BlockTree(const std::basic_string<uint64_t> &input, uint8_t r, size_t leaf_length)
    : r_(r),
      input_(input),
      leaf_length_(leaf_length) {

    if (input_.size() <= leaf_length_ || input_.size() < r)
        root_block_ = new LeafBlock(nullptr, 0, input_.size() - 1, input_);
    else {
        int number_of_leaves =
            (input_.size() % leaf_length_ == 0) ? input_.size() / leaf_length_ : input_.size() / leaf_length_ + 1;
        int height = 0;

        auto nl = number_of_leaves - 1;
        auto block_length = leaf_length_;
        while (nl) {
            height++;
            block_length *= r_;
            nl /= r_;
        }

        root_block_ = new InternalBlock(nullptr, 0, block_length - 1, input_);
    }

}

BlockTree::~BlockTree() {
    delete root_block_;
}

std::vector<std::vector<Block *>> BlockTree::levelwise_iterator() {
    std::vector<std::vector<Block *>> result = {{root_block_}};
    while (!dynamic_cast<LeafBlock *>(result.back()[0])) {
        std::vector<Block *> next_level;
        for (Block *b : result.back())
            for (Block *child : b->children(leaf_length_, r_))
                next_level.push_back(child);
        result.push_back(next_level);
    }

    return result;
}

void BlockTree::clean_unnecessary_expansions() {
    root_block_->clean_unnecessary_expansions();
    for (std::vector<Block *> level : levelwise_iterator()) {
        for (int i = 0; i < level.size(); ++i) {
            level[i]->level_index_ = i;
            level[i]->first_occurrence_level_index_ = level[i]->first_block_->level_index_;
        }
    }
}

std::vector<Block *> BlockTree::next_level(std::vector<Block *> &level, bool skip_la_leaf = false) const {
    std::vector<Block *> next_level;
    for (auto &b : level) {
        if (skip_la_leaf && b->is_la_leaf())
            continue;
        for (Block *child : b->children(leaf_length_, r_)) { // Do it in order
            child->level_index_ = next_level.size();
            child->first_occurrence_level_index_ = next_level.size();
            next_level.push_back(child);
        }
    }
    return next_level;
}

void BlockTree::forward_pair_window_block_scan(
    std::vector<Block *> &level,
    uint32_t pair_window_size,
    uint32_t N,
    block_pairs_map &pair_hashtable) {
    for (auto it = level.begin(); it != level.end();) {
        Block *b = (*it);
        b->right_ = true;
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_ + offset, pair_window_size, N); // position is always 0 here
        for (; it != level.end() && ((*it) == b || (*(it - 1))->end_ == (*it)->start_ - 1); it++) {
            Block *current = *(it);
            bool last_block = ((it + 1) == level.end() || current->end_ != (*(it + 1))->start_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < pair_window_size)
                    break;
                HashString hS(rk.hash(),
                              input_,
                              current->start_ + offset,
                              current->start_ + offset + pair_window_size - 1);
                auto result = pair_hashtable.find(hS);
                if (result != pair_hashtable.end()) {
                    // Here, It could be that the scanning should have finished with the penultimate, but it never should enter this ''if''
                    // when We're on the penultimate block and the window exceeds the last block because if that is a first occurrence should have been occured before in a pair of blocks
                    // maybe use a condition more like rk's condition below could work fine too
                    // Same logic: for when passing a window of size 2l + 2 over 2 block of length l
                    for (std::pair<Block *, Block *> p: result->second) {
                        if (current->start_ + offset < p.first->start_) {
                            p.first->left_ = true;
                            p.second->right_ = true;
                        }
                    }
                    pair_hashtable.erase(hS);
                }
                if (current->start_ + offset + pair_window_size < input_.size())
                    rk.next();
            }
        }
        (*(it - 1))->left_ = true;
    }
}

void BlockTree::forward_window_block_scan(std::vector<Block *> &level,
                                          uint32_t window_size,
                                          uint32_t N,
                                          block_map &hashtable) {
    int i = 0;
    for (auto it = level.begin(); it != level.end();) {
        Block *b = (*it);
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_ + offset, window_size, N);
        for (; it != level.end() && ((*it) == b || (*(it - 1))->end_ == (*it)->start_ - 1); it++, i++) {
            Block *current = *(it);
            bool last_block = ((it + 1) == level.end() || current->end_ != (*(it + 1))->start_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < window_size)
                    break;
                HashString hS(rk.hash(), input_, current->start_ + offset, current->start_ + offset + window_size - 1);
                auto result = hashtable.find(hS);
                if (result != hashtable.end()) {
                    std::vector<Block *> blocks = result->second;
                    for (Block *b : blocks) {
                        b->first_occurrence_level_index_ = i;
                        b->first_block_ = current;
                        b->offset_ = offset;
                        if (offset + window_size > b->first_block_->length()) b->second_block_ = (*(it + 1));
                        else b->second_block_ = nullptr;
                    }
                    hashtable.erase(hS);
                }
                if (current->start_ + offset + window_size < input_.size()) rk.next();
            }
        }
    }
}

void BlockTree::block_scan(std::vector<Block *> &level, int N, block_map &hashtable) {
    for (Block *b : level) {
        RabinKarp rk(input_, b->start_, b->length(), N);
        HashString hS(rk.hash(), input_, b->start_, b->end_);

        auto result = hashtable.find(hS);
        if (result == hashtable.end())
            hashtable[hS] = {b};
        else
            hashtable[hS].push_back(b);
    }
}

void BlockTree::process_level(std::vector<Block *> &level) {
    auto N = 6700417; //Large prime
    auto level_length = level.front()->length();

    // Block scan
    block_map hashtable;
    block_scan(level, N, hashtable);

    // Pairs of blocks scan
    block_pairs_map pair_hashtable;
    for (auto it = level.begin(); it != level.end();) {
        for (++it; (it != level.end() && (*(it - 1))->end_ == (*it)->start_ - 1); ++it) {
            Block *current = (*(it - 1));
            Block *next = (*it);
            RabinKarp rk(input_, current->start_, current->length() + next->length(), N);
            HashString hS(rk.hash(),
                          input_,
                          current->start_,
                          current->start_ + current->length() + next->length() - 1);

            auto result = pair_hashtable.find(hS);
            if (result == pair_hashtable.end())
                pair_hashtable[hS] = {{current, next}};
            else
                pair_hashtable[hS].push_back({current, next});
        }
    }

    // Window block scan
    //Establishes first occurrences of blocks
    forward_window_block_scan(level, level_length, N, hashtable);

    // Window Pair of blocks scans
    if (level.size() > 1)
        forward_pair_window_block_scan(level, level_length * 2, N, pair_hashtable);

    // BackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        Block *b = level[i];
        if (b->left_ && b->right_ && b->first_occurrence_level_index_ < b->level_index_) {
            // This doesn't have the bug of the dangling reference fixed with first_occurrence_level_index, because it
            // shouldn't happen that  block points back to a BackBlock
            auto *bb = new BackBlock(b->parent_, b->start_, b->end_, input_,
                                     level[b->first_occurrence_level_index_],
                                     b->second_block_ ? level[b->first_occurrence_level_index_ + 1] : nullptr,
                                     b->offset_);
            bb->level_index_ = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            bb->left_ = true;
            bb->right_ = true;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }

}

void BlockTree::process_back_pointers(size_t starting_block_size = std::numeric_limits<size_t>::max()) {
    std::vector<Block *> current_level = {root_block_};

    while (current_level.front()->length() > starting_block_size)
        current_level = next_level(current_level, false);

    std::stack<Block *> none_blocks;
    while (!(current_level = next_level(current_level, false)).empty()) {
        if (current_level[0]->length() < r_ || current_level[0]->length() <= leaf_length_) break;
        while (!current_level.empty() && current_level.back()->end_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        process_level(current_level);
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}