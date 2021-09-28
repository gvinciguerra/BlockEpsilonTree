// This source file is adapted from Manuel Cáceres's https://github.com/elarielcl/BlockTrees, distributed under GPL-3.0

#pragma once

#include <string>

class HashString {
public:
    uint64_t hash_;
    const std::basic_string<uint64_t> &ws_;
    size_t init_;
    size_t end_;

    HashString(uint64_t, const std::basic_string<uint64_t> &, size_t, size_t);

    bool operator==(const HashString &) const;
};

namespace std {
template<>
struct hash<HashString> {
    std::size_t operator()(const HashString &hS) const { return hS.hash_; }
};
}

HashString::HashString(uint64_t hash, const std::basic_string<uint64_t> &s, size_t init, size_t end)
    : hash_(hash), ws_(s), init_(init), end_(end) {

}

bool HashString::operator==(const HashString &other) const {
    auto length = end_ - init_ + 1;
    if (length != other.end_ - other.init_ + 1)
        return false;

    for (int i = 0; i < length; ++i)
        if (ws_[init_ + i] != other.ws_[other.init_ + i])
            return false;
    return true;
}
