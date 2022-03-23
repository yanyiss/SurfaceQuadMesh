#ifndef CGAL_DEFINITION_H
#define CGAL_DEFINITION_H
#if 0
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_Inexact_Kernal;
typedef CGAL_Inexact_Kernal::Point_2 CGAL_2_Point;
typedef CGAL::Polygon_2<CGAL_Inexact_Kernal> CGAL_2_Polygon;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> CGAL_Simple_Cartesian_Kernal;
typedef CGAL_Simple_Cartesian_Kernal::Point_3 CGAL_double_3_Point;

typedef CGAL_Simple_Cartesian_Kernal::Triangle_3 CGAL_3_Triangle;
typedef std::vector<CGAL_3_Triangle>::iterator CGAL_Triangle_Iterator;
typedef CGAL::AABB_triangle_primitive<CGAL_Simple_Cartesian_Kernal, CGAL_Triangle_Iterator> CGAL_AABB_Triangle_Primitive;
typedef CGAL::AABB_traits<CGAL_Simple_Cartesian_Kernal, CGAL_AABB_Triangle_Primitive> CGAL_AABB_triangle_traits;
typedef CGAL::AABB_tree<CGAL_AABB_triangle_traits> CGAL_AABB_Tree;

typedef CGAL_Simple_Cartesian_Kernal::Segment_3 CGAL_3_Segment;
typedef std::vector<CGAL_3_Segment>::iterator CGAL_Segment_Iterator;
typedef CGAL::AABB_segment_primitive<CGAL_Simple_Cartesian_Kernal, CGAL_Segment_Iterator> CGAL_AABB_Segment_Primitive;
typedef CGAL::AABB_traits<CGAL_Simple_Cartesian_Kernal, CGAL_AABB_Segment_Primitive> CGAL_AABB_Segment_traits;
typedef CGAL::AABB_tree<CGAL_AABB_Segment_traits> CGAL_AABB_Segment_Tree;
#endif
#endif