#pragma once

#include "CPS_Vector.h"
#include "CPS_BoundingBox.h"
#include "../src/MeshViewer/MeshDefinition.h"

namespace ClosestPointSearch
{
#ifdef USE_AVX
	class Triangle
	{
	public:
		Vec3d ver0, ver1, ver2;
		// used for accelarating.
		Vec3d face_normal;
		Vec3d edge_vec01, edge_vec02, edge_vec12;
		Vec3d edge_nor01, edge_nor02, edge_nor12;
		double edge_sqrlen01, edge_sqrlen02, edge_sqrlen12;
		bool is_plane_degenerate;
		bool has_obtuse_angle;
		// used for recording openmesh face handle
		OpenMesh::FaceHandle face_handle;
	public:
		Triangle() {}

		Triangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3)
		{
			Vec3d p12 = p2 - p1;
			Vec3d p23 = p3 - p2;
			Vec3d p31 = p1 - p3;
			// find obtuse angle, move obtuse angle point to ver0
			if (p12.dot(p31) > 0)
			{
				has_obtuse_angle = true;
				ver0 = p1;	ver1 = p2;	ver2 = p3;
			}
			else if (p12.dot(p23) > 0)
			{
				has_obtuse_angle = true;
				ver0 = p2;	ver1 = p3;	ver2 = p1;
			}
			else if (p23.dot(p31) > 0)
			{
				has_obtuse_angle = true;
				ver0 = p3;	ver1 = p1;	ver2 = p2;
			}
			else
			{
				has_obtuse_angle = false;
				ver0 = p1;	ver1 = p2;	ver2 = p3;
			}
			// calculate edge vector and squared length
			edge_vec01 = ver1 - ver0;
			edge_sqrlen01 = edge_vec01.squaredNorm();
			edge_vec02 = ver2 - ver0;
			edge_sqrlen02 = edge_vec02.squaredNorm();
			edge_vec12 = ver2 - ver1;
			edge_sqrlen12 = edge_vec12.squaredNorm();
			// calculate face normal
			face_normal = edge_vec01.cross(edge_vec02).normalized();
			is_plane_degenerate = face_normal == Vec3d(0, 0, 0);
			// calculate edge normal
			if (!is_plane_degenerate)
			{
				edge_nor01 = (edge_vec01).cross(face_normal).normalized();
				edge_nor02 = (-edge_vec02).cross(face_normal).normalized();
				edge_nor12 = (edge_vec12).cross(face_normal).normalized();
			}
		}

		Triangle(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3, const OpenMesh::FaceHandle& fh)
			:Triangle(p1, p2,p3)
		{
			face_handle = fh;
		}

		Vec3d normal()
		{
			return face_normal;
		}

		Vec3d normal() const
		{
			return face_normal;
		}

		Vec3d centroid()
		{
			return (ver0 + ver1 + ver2) / 3.0;
		}

		Vec3d centroid() const
		{
			return (ver0 + ver1 + ver2) / 3.0;
		}

		BoundingBox compute_bbox() const
		{
			Vec3d minBound(_mm256_min_pd(_mm256_min_pd(ver0.vec, ver1.vec), ver2.vec)),
				maxBound(_mm256_max_pd(_mm256_max_pd(ver0.vec, ver1.vec), ver2.vec));
			return BoundingBox(minBound, maxBound);
		}

		double closest_point(const Vec3d& query, Vec3d& result) const
		{
			result = closest_point(query);
			return (result - query).squaredNorm();
		}

		Vec3d closest_point(const Vec3d& query) const
		{
			if (is_plane_degenerate)
			{
				// If the plane is degenerate, then the triangle is degenerate, and
				// one tries to find to which segment it is equivalent.
				auto segment = find_longest_segment();
				if (is_segment_degerate(segment))
				{
					return *(segment.first);
				}

				return project_to_segment(segment, query);
			}

			if (!has_obtuse_angle)
			{
				Vec3d ver0_to_query = query - ver0;
				if (ver0_to_query.dot(edge_nor01) > 0.0)	// is outside edge 01 ?
				{
					double numerator = ver0_to_query.dot(edge_vec01);
					if (numerator < 0.0)
					{
						return ver0;
					}
					else if (numerator > edge_sqrlen01)
					{
						return ver1;
					}
					else
					{
						return ver0 + (numerator / edge_sqrlen01) * edge_vec01;
					}
				}
				else if (ver0_to_query.dot(edge_nor02) > 0.0)	// is outside edge 02 ?
				{
					double numerator = ver0_to_query.dot(edge_vec02);
					if (numerator < 0.0)
					{
						return ver0;
					}
					else if (numerator > edge_sqrlen02)
					{
						return ver1;
					}
					else
					{
						return ver0 + (numerator / edge_sqrlen02) * edge_vec02;
					}
				}
			}
			else
			{
				Vec3d ver0_to_query = query - ver0;
				if (ver0_to_query.dot(edge_nor01) > 0.0)	// is outside edge 01?
				{
					double numerator01 = ver0_to_query.dot(edge_vec01);
					if (numerator01 > edge_sqrlen01)
					{
						return ver1;
					}
					else if (numerator01 < 0.0)
					{
						// must also outside edge 02
						double numerator02 = ver0_to_query.dot(edge_vec02);
						if (numerator02 > edge_sqrlen02)
						{
							return ver2;
						}
						else if (numerator02 < 0.0)
						{
							return ver0;
						}
						else
						{
							return ver0 + (numerator02 / edge_sqrlen02) * edge_vec02;
						}
					}
					else
					{
						return ver0 + (numerator01 / edge_sqrlen01) * edge_vec01;
					}
				}
				else if (ver0_to_query.dot(edge_nor02) > 0.0)	// is outside edge 02?
				{
					double numerator = ver0_to_query.dot(edge_vec02);
					if (numerator > edge_sqrlen02)
					{
						return ver2;
					}
					else
					{
						return ver0 + (numerator / edge_sqrlen02) * edge_vec02;
					}
				}
			}
			Vec3d ver1_to_query = query - ver1;
			if (ver1_to_query.dot(edge_nor12) > 0.0)	// is outside edge 12 ?
			{
				double numerator = ver1_to_query.dot(edge_vec12);
				if (numerator < 0.0)
				{
					return ver1;
				}
				else if (numerator > edge_sqrlen12)
				{
					return ver2;
				}
				else
				{
					return ver1 + (numerator / edge_sqrlen12) * edge_vec12;
				}
			}
			double propotion = (ver1_to_query).dot(face_normal);	// is inside triangle!
			return query - propotion * face_normal;
		}
	private:
		bool is_segment_degerate(std::pair<const Vec3d*, const Vec3d*> segment) const
		{
			return *(segment.first) == *(segment.second);
		}

		std::pair<const Vec3d*, const Vec3d*> find_longest_segment()const
		{
			if (edge_sqrlen01 > edge_sqrlen12)
			{
				if (edge_sqrlen01 > edge_sqrlen02)
				{
					return std::pair<const Vec3d*, const Vec3d*>(&ver0, &ver1);
				}
				else
				{
					return std::pair<const Vec3d*, const Vec3d*>(&ver2, &ver0);
				}
			}
			else
			{
				if (edge_sqrlen12 > edge_sqrlen02)
				{
					return std::pair<const Vec3d*, const Vec3d*>(&ver1, &ver2);
				}
				else
				{
					return std::pair<const Vec3d*, const Vec3d*>(&ver2, &ver0);
				}
			}
		}

		Vec3d project_to_segment(const std::pair<const Vec3d*, const Vec3d*> segment, const Vec3d& query)const
		{
			// degeneration is a rare case, so I didn't optimize this function.
			Vec3d segment_to_vector = *segment.second - *segment.first;
			Vec3d source_to_query_x = query - *segment.first;

			double dot1 = segment_to_vector.dot(source_to_query_x);

			if (dot1 <= 0)
			{
				return *(segment.first);
			}
			else
			{
				Vec3d target_to_query = query - *segment.second;

				double dot2 = segment_to_vector.dot(target_to_query);

				if (dot2 >= 0)
				{
					return *(segment.second);
				}
				else
				{
					double num = dot1;
					double dom = segment_to_vector.squaredNorm();
					return Vec3d(*segment.first + (num / dom) * segment_to_vector);
				}
			}
		}
	};
#endif // USE_AVX
	typedef Triangle* TrianglePtr;
	typedef std::vector<Triangle> Triangles;
	typedef Triangles::iterator TriIter;
	typedef Triangles::const_iterator ConstTriIter;
}
