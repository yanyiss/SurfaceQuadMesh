#pragma once

#include "CPS_Vector.h"

namespace ClosestPointSearch
{
#ifdef USE_AVX
	/// <summary>
	/// Reference CGAL/Bbox_3.h
	/// </summary>
	class BoundingBox
	{
	protected:
		__m256d minBound, maxBound;
	public:
		BoundingBox()
		{
			minBound = _mm256_setzero_pd();
			maxBound = _mm256_setzero_pd();
		}

		BoundingBox(const Vec3d& minB, const Vec3d& maxB)
		{
			minBound = minB.vec;
			maxBound = maxB.vec;
		}

		BoundingBox operator+(const BoundingBox& b) const
		{
			BoundingBox result;

			result.minBound = _mm256_min_pd(minBound, b.minBound);
			result.maxBound = _mm256_max_pd(maxBound, b.maxBound);

			return result;
		}

		BoundingBox& operator+=(const BoundingBox& b)
		{
			minBound = _mm256_min_pd(minBound, b.minBound);
			maxBound = _mm256_max_pd(maxBound, b.maxBound);

			return *this;
		}

		int longest_axis() const
		{
			__m256d diff = _mm256_sub_pd(maxBound, minBound);
			double* d = to_double_ptr(diff);
			double dx = d[0], dy = d[1], dz = d[2];

			if (dx >= dy)
			{
				if (dx >= dz)
				{
					return 0;
				}
				else // dz>dx and dx>=dy
				{
					return 2;
				}
			}
			else // dy>dx
			{
				if (dy >= dz)
				{
					return 1;
				}
				else  // dz>dy and dy>dx
				{
					return 2;
				}
			}
		}

		bool do_intersect_sphere(const Vec3d& center, const double square_radius) const
		{

			__m256d FourZero = _mm256_setzero_pd();
			__m256d C = center.vec;

			__m256d T1, T2, T3, T4;
			T1 = _mm256_sub_pd(minBound, C);	// T1 = minBound - C
			T2 = _mm256_max_pd(T1, FourZero);	// T2 = max{T1, FourZero} 
			T3 = _mm256_mul_pd(T2, T2);			// T3 = T2 * T2

			T1 = _mm256_sub_pd(C, maxBound);	// T1 = C - maxBound
			T2 = _mm256_max_pd(T1, FourZero);	// T2 = max{T1, FourZero}
			T4 = _mm256_fmadd_pd(T2, T2, T3);	// T4 = T2 * T2

			return hsum_double_avx(T4) <= square_radius;
		}
	};
#endif // USE_AVX
}
