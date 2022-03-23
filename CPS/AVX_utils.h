#pragma once

#ifdef USE_AVX
#include <immintrin.h>

inline double* to_double_ptr(__m256d& m)
{
    return (double*)&m;
}

inline const double* to_double_ptr(const __m256d& m)
{
    return (const double*)&m;
}

/// <summary>
/// https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx/49943540#49943540
/// </summary>
/// <param name="v"></param>
/// <returns>horizontal sum of v</returns>
inline double hsum_double_avx(const __m256d& v) {
    __m128d vlow = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);        // high 128
    vlow = _mm_add_pd(vlow, vhigh);                     // reduce down to 128

    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));    // reduce to scalar
}

/// <summary>
/// https://stackoverflow.com/questions/9795529/how-to-find-the-horizontal-maximum-in-a-256-bit-avx-vector
/// </summary>
/// <param name="x"></param>
/// <returns>horizontal maximum of x</returns>
inline double hmax_double_avx(const __m256d& x)
{
    __m256d y = _mm256_permute2f128_pd(x, x, 1); // permute 128-bit values
    __m256d m1 = _mm256_max_pd(x, y); // m1[0] = max(x[0], x[2]), m1[1] = max(x[1], x[3]), etc.
    __m256d m2 = _mm256_permute_pd(m1, 5); // set m2[0] = m1[1], m2[1] = m1[0], etc.
    __m256d m = _mm256_max_pd(m1, m2); // all m[0] ... m[3] contain the horizontal max(x[0], x[1], x[2], x[3])
    return *to_double_ptr(m);
}

#endif

