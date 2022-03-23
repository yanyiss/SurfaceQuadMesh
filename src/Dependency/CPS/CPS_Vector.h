#pragma once

#include <iostream>
#include <vector>
#include <array>
#include "AVX_utils.h"

namespace ClosestPointSearch
{
#ifdef USE_AVX

	class Vec3d
	{
	public:
		__m256d vec;
	public:
		Vec3d() {}

		Vec3d(const double& v0, const double& v1, const double& v2)
		{
			vec = _mm256_set_pd(0.0, v2, v1, v0);
		}

		Vec3d(const double* v)
		{
			vec = _mm256_set_pd(0.0, v[2], v[1], v[0]);
		}

		Vec3d(const Vec3d& v)
			:vec(v.vec)
		{}

		Vec3d(const __m256d& v)
			:vec(v)
		{}

		~Vec3d()
		{}

		double x()const
		{
			return to_double_ptr(vec)[0];
		}

		double y()const
		{
			return to_double_ptr(vec)[1];
		}

		double z()const
		{
			return to_double_ptr(vec)[2];
		}

		double& operator[](int dim)
		{
			return to_double_ptr(vec)[dim];
		}

		const double operator[](int dim)const
		{
			return to_double_ptr(vec)[dim];
		}

		Vec3d& operator=(const Vec3d& v)
		{
			vec = v.vec;
			return *this;
		}

		inline Vec3d operator-() const
		{
			__m256d zero = _mm256_setzero_pd();
			return Vec3d(_mm256_sub_pd(zero, vec));
		}

		inline Vec3d operator-(const Vec3d& rhs) const
		{
			return Vec3d(_mm256_sub_pd(vec, rhs.vec));
		}

		inline Vec3d operator+(const Vec3d& rhs) const
		{
			return Vec3d(_mm256_add_pd(vec, rhs.vec));
		}

		inline Vec3d operator*(const double rhs) const
		{
			__m256d m = _mm256_set_pd(rhs, rhs, rhs, rhs);
			return Vec3d(_mm256_mul_pd(vec, m));
		}

		inline Vec3d operator/(const double rhs) const
		{
			__m256d d = _mm256_set_pd(rhs, rhs, rhs, rhs);
			return Vec3d(_mm256_div_pd(vec, d));
		}

		inline bool operator==(const Vec3d& rhs) const
		{
			return _mm256_movemask_pd(_mm256_cmp_pd(vec, rhs.vec, _CMP_EQ_OS)) == 0xf;
		}

		inline bool operator!=(const Vec3d& rhs) const
		{
			return !(*this == rhs);
		}

		inline Vec3d cross(const Vec3d& rhs) const
		{
			// https://www.nersc.gov/assets/Uploads/Language-Impact-on-Vectorization-Vector-Programming-in-C++.pdf
			__m256d a201 = _mm256_permute4x64_pd(vec, _MM_SHUFFLE(3, 1, 0, 2));
			__m256d b201 = _mm256_permute4x64_pd(rhs.vec, _MM_SHUFFLE(3, 1, 0, 2));
			__m256d tmp = _mm256_fmsub_pd(rhs.vec, a201, _mm256_mul_pd(vec, b201));
			tmp = _mm256_permute4x64_pd(tmp, _MM_SHUFFLE(3, 1, 0, 2));
			tmp = _mm256_blend_pd(_mm256_setzero_pd(), tmp, 0x7); // put zero on 4th position
			return Vec3d(tmp);

			// https://gist.github.com/garrettsickles/85a9ab8385172bd0e762f38e4cfb045f
			/*
			return Vec3d(
				_mm256_permute4x64_pd(
					_mm256_sub_pd(
						_mm256_mul_pd(vec, _mm256_permute4x64_pd(rhs.vec, _MM_SHUFFLE(3, 0, 2, 1))),
						_mm256_mul_pd(rhs.vec, _mm256_permute4x64_pd(vec, _MM_SHUFFLE(3, 0, 2, 1)))
					), _MM_SHUFFLE(3, 0, 2, 1)
				));
			*/
		}

		inline double dot(const Vec3d& rhs) const
		{
			return hsum_double_avx(_mm256_mul_pd(vec, rhs.vec));
		}

		inline Vec3d normalized() const
		{
			double z = squaredNorm();
			return z > 0 ? (*this) / sqrt(z) : (*this);
		}

		inline double squaredNorm() const
		{
			return hsum_double_avx(_mm256_mul_pd(vec, vec));
		}

		inline double norm() const
		{
			return std::sqrt(squaredNorm());
		}

		inline double less_on(const int dim, const Vec3d& rhs)const
		{
			return to_double_ptr(vec)[dim] < to_double_ptr(rhs.vec)[dim];
		}

		inline double less_on(const int dim, const double rhs)const
		{
			return to_double_ptr(vec)[dim] < rhs;
		}
	};

	inline Vec3d operator*(const double lhs, const Vec3d& rhs)
	{
		return rhs * lhs;
	}
#endif	// USE_AVX

	typedef std::vector<Vec3d>		Vec3ds;
	typedef Vec3ds::iterator		Vec3dIter;
	typedef Vec3d*					Vec3dPtr;
	typedef std::vector<Vec3dPtr>	Vec3dPtrs;
	typedef Vec3dPtrs::iterator		Vec3dPtrIter;
}
