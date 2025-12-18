#pragma once
#include <common/config.hpp>
namespace preprocess {
#if RANS_USE_STRONG_INLINE && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
#define RANS_STRONG_INLINE __forceinline
#else
#define RANS_STRONG_INLINE inline
#endif

#if RANS_USE_ALWAYS_INLINE && defined(__GNUC__)
#define RANS_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define RANS_ALWAYS_INLINE RANS_STRONG_INLINE
#endif
}  // namespace preprocess