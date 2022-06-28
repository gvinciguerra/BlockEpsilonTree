// This file is part of <https://github.com/gvinciguerra/BlockEpsilonTree>.
// Copyright (c) 2021 Giorgio Vinciguerra.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// NOTE: some methods implemented here are patent pending.

#pragma once

#include "LABlock.hpp"
#include "block_tree/BlockTree.h"
#include "block_tree/blocks/Block.h"

#include <sdsl/bits.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_scan.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

//#define BLOCK_EPS_DEBUG

#ifdef BLOCK_EPS_DEBUG
#define IF_DEBUG(X) { X; }
#define DEBUG_OUT(X) { std::cout << X; }
#define DEBUG_OUTLN(X) { std::cout << X << std::endl; }
#else
#define IF_DEBUG(X)
#define DEBUG_OUT(X)
#define DEBUG_OUTLN(X)
#endif

class BlockEpsilonTree {
    using rank_type = sdsl::rank_support_v5<>;
    uint8_t branching_factor = 0;
    size_t starting_block_size = 0;
    uint8_t number_of_levels = 0;
    uint8_t ptr_shift = 0;
    uint64_t leaves_shift = 0;

    /** Elias-Fano representation without efficient rank support. */
    using elias_fano = sdsl::sd_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::select_support_scan<0>>;

    std::vector<sdsl::bit_vector> bitvectors;          ///< 1 if internal block, 0 if either left pointer or an LA-block
    std::vector<rank_type> bitvectors_rank;            ///< the rank structure associated to the corresponding level bv
    std::vector<sdsl::int_vector<>> pointers;          ///< one pointer for each 0 in the corresponding level bv
    std::vector<sdsl::int_vector<>> difference_values; ///< one value for each 0 in the corresponding level bv
    std::vector<BlockLAType> la_blocks;                ///< LA-blocks sorted according to the block tree in-order visit
    sdsl::sd_vector<> leaves;
    sdsl::sd_vector<>::rank_1_type leaves_rank;
    sdsl::sd_vector<>::select_1_type leaves_select;
//    sdsl::int_vector<> leaves;

    std::map<std::string, std::string> metadata;

    // Structures needed for rank
    std::vector<sdsl::int_vector<>> samples; ///< per-level samples. one sample for each child of an internal block

    class Pointer;

public:

    BlockEpsilonTree(std::vector<uint32_t> &data,
                     uint8_t branching_factor,
                     size_t leaf_size = 8,
                     size_t hint_starting_block_size = 1 << 22)
        : branching_factor(branching_factor) {
        assert(std::is_sorted(data.begin(), data.end()));

        // (1) Construct the pointer-based block tree
        std::basic_string<uint64_t> gap_string;
        gap_string.reserve(data.size() + hint_starting_block_size);
        gap_string.push_back(uint64_t(data.front()));
        for (auto it = data.begin() + 1; it != data.end(); ++it)
            gap_string.push_back(uint64_t(*it - *std::prev(it)));

        BlockTree bt(gap_string, branching_factor, leaf_size);
        bt.process_back_pointers(hint_starting_block_size);
        bt.clean_unnecessary_expansions();

        // (1.1) Find the first level, which is the one containing at least one back block
        std::vector<Block *> first_level = {bt.root_block_};
        while (true) {
            auto it = std::find_if(first_level.begin(), first_level.end(), [](auto &b) { return b->is_leaf(); });
            if (it != first_level.end())
                break;
            first_level = bt.next_level(first_level, false);
        }

        // (1.2) Assign a cost to the nodes of the block tree
        const auto n = data.size();
        const auto u = data.back();
        const auto log_n = BIT_WIDTH(n);
        const auto log_u = BIT_WIDTH(u);

        size_t leaf_elements = 0;
        for (auto &b : first_level)
            leaf_elements += b->count_leaf_elements();
        const auto leaf_approx_cost = std::ceil(ef_cost(leaf_elements, u - data.front()) / leaf_elements);

        auto cost = 0;
        for (auto &b : first_level)
            cost += b->la_visit(data.data(), log_n, log_u, 0, leaf_approx_cost);

        // (1.3) Assign IDs to LA blocks
        size_t la_blocks_count = 0;
        for (auto &b : first_level)
            la_blocks_count = b->number_la_blocks(la_blocks_count);

        la_blocks.resize(la_blocks_count);
        ptr_shift = la_blocks_count > 0 ? BIT_WIDTH(la_blocks_count - 1) : 0;
        if (log_n + ptr_shift + Pointer::selector_width > 64)
            throw std::overflow_error("Pointers too large");

        // (1.4) Compute stats on the pruned tree
        size_t depth = 0;
        long double depth_sum = 0;
        for (auto &b : first_level) {
            auto[subtree_max_depth, subtree_depth_sum] = b->compute_depth(0);
            depth = std::max<size_t>(depth, subtree_max_depth);
            depth_sum += subtree_depth_sum;
        }
        const auto average_depth = depth_sum / data.size();

        // (2) Level-wise compress the tree topology
        bitvectors_rank.reserve(depth);
        pointers.reserve(depth);
        bitvectors.reserve(depth);
        difference_values.reserve(depth);

        starting_block_size = first_level[0]->length();
        number_of_levels = 0;

        auto current_level = std::move(first_level);
        auto next_level = bt.next_level(current_level, true);

        sdsl::int_vector<> top_samples;
        top_samples.reserve(current_level.size());
        for (auto &b : current_level)
            top_samples.push_back(b->end_ < n ? data[b->end_] : u);
        sdsl::util::bit_compress(top_samples);
        samples.push_back(std::move(top_samples));

        while (!next_level.empty()) {
            // (2.1) Compress the level
            sdsl::bit_vector level_bv(current_level.size());
            sdsl::int_vector<> level_difference_values;
            sdsl::int_vector<> level_pointers;
            sdsl::int_vector<> level_samples;
            level_difference_values.reserve(current_level.size());
            level_pointers.reserve(current_level.size());
            level_samples.reserve(current_level.size());
            DEBUG_OUT("Level " << int(number_of_levels) << " samples: ");

            auto block_size = current_level.front()->length();
            for (size_t i = 0; i < current_level.size(); ++i) {
                auto *&block = current_level[i];
                block->level_start_pos_ = i;

                if (block->is_la_leaf()) {
                    assert(!block->is_leaf());

                    BlockLAType la;
                    if (!block->get_la(data.data(), la))
                        throw std::runtime_error("");
                    la_blocks[block->la_block_id_] = la;
                    DEBUG_OUTLN("Pruned " << block->start_ << " " << block->end_ << " pos=" << block->la_block_id_);

                    level_bv[i] = false;
                    level_pointers.push_back(Pointer::make(0b11, 0, block->la_block_id_, ptr_shift));
                    level_difference_values.push_back(0);
                } else if (block->is_leaf()) {
                    auto *&left = block->first_block_;
                    auto *&right = block->second_block_;
                    assert(right || block->offset_ == 0);

                    level_bv[i] = false;
                    level_difference_values.push_back(data[block->start_] - data[left->start_ + block->offset_]);

                    uint64_t pointer;
                    uint64_t position;
                    auto left_is_pruned = left->parent_->is_la_leaf();
                    auto right_is_pruned = (right && right->parent_->is_la_leaf()) || (left_is_pruned && !right);
                    auto selector = uint8_t((left_is_pruned << 1) | right_is_pruned);
                    switch (selector) {
                        case 0b00: {
                            position = left->level_start_pos_ * block_size + block->offset_;
                            pointer = Pointer::make(selector, position, Pointer::nil, ptr_shift);
                            break;
                        }

                        case 0b01: {
                            position = left->level_start_pos_ * block_size + block->offset_;
                            pointer = Pointer::make(selector, position, right->parent_->la_block_id_, ptr_shift);
                            break;
                        }

                        case 0b10: {
                            auto offset = int64_t(right->level_start_pos_ * block_size + block->offset_) - block_size;
                            auto shift = left->start_ - left->parent_->start_;
                            position = offset + shift;
                            pointer = Pointer::make(selector, position, left->parent_->la_block_id_, ptr_shift);
                            break;
                        }

                        case 0b11: {
                            position = block->offset_ + left->start_ - left->parent_->start_;
                            pointer = Pointer::make(selector, position, left->parent_->la_block_id_, ptr_shift);
                            break;
                        }
                    }

                    level_pointers.push_back(pointer);
                } else {
                    level_bv[i] = true;
                    for (auto &c : block->children_)
                        level_samples.push_back(c->end_ >= n ? u : data[c->end_]);
                }
            }

            // (2.2) Store the level
            DEBUG_OUTLN("BITV:\t" << level_bv << std::endl
                                  << "DIFF:\t" << level_difference_values << std::endl
                                  << "PTRS:\t" << level_pointers << std::endl
                                  << std::string(80, '-'));

            samples.emplace_back(std::move(level_samples));
            pointers.emplace_back(std::move(level_pointers));
            bitvectors.emplace_back(std::move(level_bv));
            difference_values.emplace_back(std::move(level_difference_values));
            sdsl::util::bit_compress(pointers.back());
            sdsl::util::bit_compress(difference_values.back());
            sdsl::util::bit_compress(samples.back());

            current_level = std::move(next_level);
            next_level = bt.next_level(current_level, true);
            ++number_of_levels;
        }

        // (3) Init auxiliary structures and prepare leaf string
        if (!samples.empty())
            samples.pop_back();

        for (auto &bv : bitvectors)
            bitvectors_rank.emplace_back(&bv);

        ++number_of_levels;

//        std::vector<uint32_t> tmp_leaves;
//        tmp_leaves.reserve(current_level.size() * bt.leaf_length_);
//        for (const auto &b: current_level) {
//            std::copy(data.begin() + b->start_,
//                      std::min(data.end(), data.begin() + b->end_ + 1),
//                      std::back_inserter(tmp_leaves));
//        }
//        auto leaves_count = tmp_leaves.size();
//        leaves = decltype(leaves)(tmp_leaves.begin(), tmp_leaves.end());

        size_t leaves_count = 0;
        for (const auto &b: current_level)
            leaves_count += b->actual_length();

        if (leaves_count) {
            leaves_shift = data[current_level.front()->start_];
            auto leaves_u = data[std::min<size_t>(data.size() - 1, current_level.back()->end_)] - leaves_shift;
            sdsl::sd_vector_builder builder(leaves_u + 1, leaves_count);
            for (const auto &b: current_level)
                for (size_t j = b->start_; j < std::min<size_t>(data.size(), b->end_ + 1); ++j)
                    builder.set(data[j] - leaves_shift);
            leaves = decltype(leaves)(builder);
        }
        sdsl::util::init_support(leaves_rank, &leaves);
        sdsl::util::init_support(leaves_select, &leaves);

        // (4) Compute stats / metadata
        size_t la_bytes = 0;
        for (auto &l: la_blocks)
            la_bytes += l.size_in_bytes();

        size_t internal_nodes_count = 0;
        size_t back_pointers_count = 0;
        for (size_t i = 0; i < bitvectors.size(); ++i) {
            auto num_ones = bitvectors_rank[i].rank(bitvectors[i].size());
            internal_nodes_count += num_ones;
            back_pointers_count += bitvectors[i].size() - num_ones;
        }
        back_pointers_count -= la_blocks_count;

        auto to_bpi = [&](auto bytes) { return std::to_string(bytes * 8. / n); };
        metadata["bpi"] = to_bpi(size_in_bytes());
        metadata["bitvectors_bpi"] = to_bpi(sdsl::size_in_bytes(bitvectors));
        metadata["bitvectors_rank_bpi"] = to_bpi(sdsl::size_in_bytes(bitvectors_rank));
        metadata["pointers_bpi"] = to_bpi(sdsl::size_in_bytes(pointers));
        metadata["difference_values_bpi"] = to_bpi(sdsl::size_in_bytes(difference_values));
        metadata["samples_bpi"] = to_bpi(sdsl::size_in_bytes(samples));
        metadata["leaves_bpi"] = to_bpi(sdsl::size_in_bytes(leaves));
        metadata["la_blocks_bpi"] = to_bpi(la_bytes);
        metadata["la_blocks_count"] = std::to_string(la_blocks_count);
        metadata["internal_nodes_count"] = std::to_string(internal_nodes_count);
        metadata["back_pointers_count"] = std::to_string(back_pointers_count);
        metadata["leaves_count"] = std::to_string(leaves_count);
        metadata["depth"] = std::to_string(depth);
        metadata["average_depth"] = std::to_string(average_depth);
        metadata["starting_block_size"] = std::to_string(starting_block_size);
    }

    std::map<std::string, std::string> get_metadata() { return metadata; }

    size_t size_in_bytes() const {
        size_t sum = 0;
        sum += sdsl::size_in_bytes(bitvectors);
        sum += sdsl::size_in_bytes(bitvectors_rank);
        sum += sdsl::size_in_bytes(pointers);
        sum += sdsl::size_in_bytes(difference_values);
        sum += sdsl::size_in_bytes(samples);
        for (auto &l: la_blocks)
            sum += l.size_in_bytes();
        sum += sdsl::size_in_bytes(leaves);
        return sum;
    }

    size_t rank(uint64_t x) const {
        if (samples.empty())
            return leaves_rank(std::max(x, leaves_shift) - leaves_shift);

        auto block_size = starting_block_size;
        auto &top_samples = samples.front();
        auto block = size_t(std::lower_bound(top_samples.begin(), top_samples.end(), x) - top_samples.begin());
        auto x_remapped = x;
        auto input_shift = block * block_size;

        for (auto level = 0; level < number_of_levels - 1; ++level) {
            if (!bitvectors[level][block]) {
                auto[diff_val, ptr] = get_leftward_data(level, block);
                x_remapped -= diff_val;

                switch (ptr.get_selector()) {
                    case 0b00: {
                        block = ptr.get_position() / block_size;
                        input_shift -= ptr.get_position() % block_size;
                        auto sample = samples[level][block];
                        if (x_remapped > sample) {
                            ++block;
                            input_shift += block_size;
                        }
                        break;
                    }

                    case 0b01: {
                        block = ptr.get_position() / block_size;
                        input_shift -= ptr.get_position() % block_size;
                        auto sample = samples[level][block];
                        if (x_remapped > sample) {
                            auto &la_block = la_blocks[ptr.get_right_la_block()];
                            return input_shift + block_size + la_block.rank(x_remapped, 0, block_size);
                        }
                        break;
                    }

                    case 0b10: {
                        auto &la_block = la_blocks[ptr.get_left_la_block()];
                        input_shift -= ptr.get_position() % block_size;
                        block = ptr.get_position() / block_size + 1 - (la_block.size() - block_size) / block_size;
                        auto sample = samples[level][block];
                        if (x_remapped <= sample) {
                            auto lo = la_block.size() - block_size;
                            auto hi = la_block.size();
                            return input_shift + la_block.rank(x_remapped, lo, hi) - (la_block.size() - block_size);
                        }
                        input_shift += block_size;
                        break;
                    }

                    case 0b11: {
                        auto k = ptr.get_position();
                        input_shift -= k;
                        auto &left_la_block = la_blocks[ptr.get_left_la_block()];
                        if (k < left_la_block.size()) {
                            auto lo = ptr.get_position();
                            auto hi = left_la_block.size();
                            return input_shift + left_la_block.rank(x_remapped, lo, hi);
                        }
                        auto &right_la_block = la_blocks[ptr.get_right_la_block()];
                        auto lo = ptr.get_position() - left_la_block.size();
                        auto hi = ptr.get_position() + block_size;
                        return input_shift + right_la_block.rank(x_remapped, lo, hi);
                    }
                }
            }

            block = bitvectors_rank[level].rank(block) * branching_factor;
            uint8_t child = 0;
            if (level != number_of_levels - 2) { // There are no samples at the last level because we use EF with rank
                for (; child < branching_factor - 1; ++child) {
                    if (block + child >= samples[level + 1].size())
                        break;
                    auto child_ub = samples[level + 1][block + child];
                    if (x_remapped <= child_ub)
                        break;
                }
            }
            block = block + child;
            block_size /= branching_factor;
            input_shift += child * block_size;
        }

        // Uncomment if there is no efficient rank support on the leaves
        //     auto start = block * block_size;
        //     auto j = start;
        //     uint64_t val = 0;
        //     for (; j < start + block_size; ++j) {
        //         val = leaves_shift + leaves_select(j + 1);
        //         if (val >= x_remapped)
        //             break;
        //     }
        //     return input_shift + j - start;

        auto leaves_rank_val = x_remapped <= leaves_shift ? 0 : leaves_rank(x_remapped - leaves_shift);
        return input_shift + leaves_rank_val - block * block_size;
    }

    uint64_t select(size_t i) const {
        assert(i > 0);
        return operator[](i - 1);
    }

    uint64_t operator[](size_t i) const {
        auto block_size = starting_block_size;
        auto block = i / block_size;
        auto offset = i % block_size;
        auto result = uint64_t(0);

        for (auto level = 0; level < number_of_levels - 1; ++level) {
            if (!bitvectors[level][block]) {
                auto[diff_val, ptr] = get_leftward_data(level, block);
                result += diff_val;

                switch (ptr.get_selector()) {
                    case 0b00:
                        offset += ptr.get_position() % block_size;
                        block = ptr.get_position() / block_size;
                        if (offset >= block_size) {
                            ++block;
                            offset -= block_size;
                        }
                        break;

                    case 0b01:
                        offset += ptr.get_position() % block_size;
                        block = ptr.get_position() / block_size;
                        if (offset >= block_size)
                            return result + la_blocks[ptr.get_right_la_block()][offset - block_size];
                        break;

                    case 0b10: {
                        offset += ptr.get_position() % block_size;
                        auto &la_block = la_blocks[ptr.get_left_la_block()];
                        if (offset < block_size)
                            return result + la_block[la_block.size() - block_size + offset];
                        block = ptr.get_position() / block_size + 1 - (la_block.size() - block_size) / block_size;
                        offset -= block_size;
                        break;
                    }

                    case 0b11:
                        offset += ptr.get_position();
                        auto &left_la_block = la_blocks[ptr.get_left_la_block()];
                        if (offset < left_la_block.size())
                            return result + left_la_block[offset];
                        auto &right_la_block = la_blocks[ptr.get_right_la_block()];
                        return result + right_la_block[offset - left_la_block.size()];
                }
            }
            block_size /= branching_factor;
            auto child = offset / block_size;
            block = bitvectors_rank[level].rank(block) * branching_factor + child;
            offset -= child * block_size;
        }

        return result + leaves_shift + leaves_select(block * block_size + offset + 1);
//         return result + leaves_shift + leaves[block * block_size + offset];
    }

private:

    class Pointer {
        uint8_t selector;
        size_t position;
        size_t la_block;

    public:

        static constexpr uint8_t selector_width = 2;
        static constexpr auto nil = std::numeric_limits<size_t>::max();

        Pointer(uint64_t value, uint8_t shift) {
            selector = value & sdsl::bits::lo_set[selector_width];
            if (selector) {
                la_block = (value >> selector_width) & sdsl::bits::lo_set[shift];
                la_block -= selector == 0b01;
                position = value >> (selector_width + shift);
            } else
                position = value >> selector_width;
        }

        static uint64_t make(uint8_t selector, size_t offset, size_t la_block, uint8_t shift) {
            if (selector == 0) {
                assert(la_block == nil);
                return offset << selector_width;
            }

            assert(selector <= 3);
            assert(la_block != nil);
            uint64_t ptr = 0;
            ptr |= selector;
            ptr |= la_block << selector_width;
            ptr |= offset << (selector_width + shift);
            return ptr;
        }

        uint8_t get_selector() const { return selector; }
        size_t get_position() const { return position; }
        size_t get_left_la_block() const { return la_block; }
        size_t get_right_la_block() const { return la_block + 1; }
    };

    std::pair<uint64_t, Pointer> get_leftward_data(size_t level, size_t block) const {
        auto rank0 = block - bitvectors_rank[level].rank(block);
        return {difference_values[level][rank0], Pointer(pointers[level][rank0], ptr_shift)};
    }

};