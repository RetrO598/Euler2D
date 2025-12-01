#pragma once

#include <cstddef>

namespace preprocess {
#ifndef RANS_USE_STRONG_INLINE
#define RANS_USE_STRONG_INLINE 1
#endif

#ifndef RANS_USE_ALWAYS_INLINE
#define RANS_USE_ALWAYS_INLINE 1
#endif

#ifdef USE_FLOAT
using real = float;
#else
using real = double;
#endif

// Index type and invalid sentinel
using Index = std::size_t;
constexpr Index INVALID_INDEX = static_cast<Index>(-1);

// Helper: check if index is valid
inline constexpr bool is_valid_index(Index i) noexcept {
  return i != INVALID_INDEX;
}

// Conversion: Index to int (e.g., for printing/file output expecting -1)
inline constexpr int index_to_int(Index i) noexcept {
  return is_valid_index(i) ? static_cast<int>(i) : -1;
}

// Conversion: int (possibly negative) to Index
inline constexpr Index int_to_index(int v) noexcept {
  return (v < 0) ? INVALID_INDEX : static_cast<Index>(v);
}

}  // namespace preprocess