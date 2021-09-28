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

#include "la_vector/piecewise_linear_model.hpp"
#include <sdsl/int_vector.hpp>
#include <cstdint>
#include <memory>
#include <type_traits>

/** Computes (bits_per_correction > 0 ? 2^(bits_per_correction-1) - 1 : 0) without the conditional operator. */
#define BPC_TO_EPSILON(bits_per_correction) (((1ul << (bits_per_correction)) + 1) / 2 - 1)

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

/** Computes the smallest integral value not less than x / y, where x and y must be positive integers. */
#define CEIL_UINT_DIV(x, y) ((x) / (y) + ((x) % (y) != 0))

#pragma pack(push, 1)

template<typename K>
class LABlock {
    static constexpr uint8_t field_bits_for_bpc = 5;
    static constexpr uint8_t field_bits_for_n = 32 - field_bits_for_bpc;

    K slope_numerator;
    K slope_denominator;
    std::make_signed_t<K> intercept;
    uint32_t n: field_bits_for_n;
    uint8_t bpc: field_bits_for_bpc;
    std::unique_ptr<uint64_t[]> corrections;
    using larger_key_type = typename std::conditional_t<sizeof(K) <= 4, int64_t, __int128>;

public:

    LABlock() = default;

    LABlock(const LABlock &la)
        : slope_numerator(la.slope_numerator),
          slope_denominator(la.slope_denominator),
          intercept(la.intercept),
          n(la.n),
          bpc(la.bpc),
          corrections() {
        if (bpc) {
            auto s = CEIL_UINT_DIV(bpc * n, 64) + 1;
            corrections = std::make_unique<uint64_t[]>(s);
            std::copy(la.corrections.get(), la.corrections.get() + s, corrections.get());
        }
    }

    LABlock &operator=(const LABlock &la) {
        slope_numerator = la.slope_numerator;
        slope_denominator = la.slope_denominator;
        intercept = la.intercept;
        n = la.n;
        bpc = la.bpc;
        if (bpc) {
            auto s = CEIL_UINT_DIV(bpc * n, 64) + 1;
            corrections = std::make_unique<uint64_t[]>(s);
            std::copy(la.corrections.get(), la.corrections.get() + s, corrections.get());
        }
        return *this;
    }

    larger_key_type approximate(size_t i) const {
        return larger_key_type(slope_numerator * i) / slope_denominator + intercept;
    }

    K operator[](size_t i) const {
        if (!bpc)
            return approximate(i);
        auto epsilon = BPC_TO_EPSILON(bpc);
        auto j = i * bpc;
        auto correction = sdsl::bits::read_int(corrections.get() + (j >> 6u), j & 0x3F, bpc);
        return approximate(i) + correction - epsilon;
    }

    size_t rank(const K &value) const {
        auto[pos, bound] = approximate_position(value);
        auto lo = pos <= bound ? 0 : pos - bound;
        auto hi = std::min<size_t>(pos + bound + 1, n);

        while (lo < hi) {
            auto mid = lo + (hi - lo) / 2;
            if (operator[](mid) < value)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo;
    }

    size_t rank(const K &value, size_t lo_bound, size_t hi_bound) const {
        auto[pos, bound] = approximate_position(value);
        auto lo = std::max<size_t>(pos <= bound ? 0 : pos - bound, lo_bound);
        auto hi = std::min<size_t>(pos + bound + 1, hi_bound);

        while (lo < hi) {
            auto mid = lo + (hi - lo) / 2;
            if (operator[](mid) < value)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo;
    }

    std::pair<size_t, size_t> approximate_position(const K &value) const {
        auto numerator = std::max<larger_key_type>(1, slope_numerator);
        auto position = ((larger_key_type(value) - intercept) * slope_denominator) / numerator;
        auto epsilon = larger_key_type(BPC_TO_EPSILON(this->bpc));
        auto bound = 1 + (epsilon * slope_denominator) / numerator;
        return {std::clamp<larger_key_type>(position, 0, n), bound};
    }

    template<typename RandomIt>
    static bool make(RandomIt begin, RandomIt end, uint8_t bpc, LABlock &out) {
        const auto n = (size_t) std::distance(begin, end);
        if (BIT_WIDTH(n) > field_bits_for_n)
            throw std::overflow_error("increase bits assigned to n");
        if (BIT_WIDTH(bpc) > field_bits_for_bpc)
            throw std::overflow_error("increase bits assigned to bpc");

        const auto epsilon = BPC_TO_EPSILON(bpc);

        OptimalPiecewiseLinearModel<K, K> opt(epsilon);
        opt.add_point(0, begin[0]);

        for (size_t i = 1; i < n; ++i)
            if (!opt.add_point(i, begin[i]))
                return false;

        auto cs = opt.get_segment();
        auto max_slope = cs.rectangle[3] - cs.rectangle[1];
        auto intercept_numerator = cs.rectangle[3].x * cs.rectangle[1].y - cs.rectangle[1].x * cs.rectangle[3].y;
        out.slope_numerator = max_slope.dy;
        out.slope_denominator = max_slope.dx;
        out.intercept = max_slope.dx == 0 ? begin[0] : intercept_numerator / max_slope.dx;
        out.n = n;
        out.bpc = bpc;
        if (bpc == 0) {
            out.corrections = nullptr;
            return true;
        }

        out.corrections = std::make_unique<uint64_t[]>(CEIL_UINT_DIV(bpc * n, 64) + 1);
        for (size_t i = 0; i < n; ++i) {
            auto error = begin[i] - out.approximate(i);
            auto correction = uint64_t(error + epsilon);
            if (BIT_WIDTH(correction) > bpc)
                throw std::overflow_error("Segment correction too large");
            auto j = i * bpc;
            sdsl::bits::write_int(out.corrections.get() + (j >> 6), correction, j & 0x3F, bpc);
        }

        return true;
    }

    size_t size() const { return n; }

    uint8_t get_bpc() const { return bpc; }

    size_t size_in_bytes() const { return sizeof(*this) + n * bpc / 8; }
};

#pragma pack(pop)