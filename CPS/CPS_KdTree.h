#pragma once
#include <deque>
#include "CPS_Triangle.h"

namespace ClosestPointSearch
{
#ifdef USE_AVX
    class KdTree_Vec3d : public Vec3d
    {
    public:
        TriIter m_from;
    public:
        KdTree_Vec3d()
            :Vec3d()
        {}

        KdTree_Vec3d(const double& v0, const double& v1, const double& v2)
            :Vec3d(v0, v1, v2)
        {}

        KdTree_Vec3d(const double* v)
            :Vec3d(v)
        {}

        KdTree_Vec3d(const Vec3d& v)
            :Vec3d(v)
        {}

        KdTree_Vec3d(const Vec3d& v, TriIter tri)
            :Vec3d(v), m_from(tri)
        {}

        KdTree_Vec3d(const __m256d& v)
            :Vec3d(v)
        {}

        KdTree_Vec3d& operator=(const KdTree_Vec3d& v) = default;
    };

    typedef KdTree_Vec3d                    KdVec3d;
    typedef std::vector<KdVec3d>	    	KdVec3ds;
    typedef KdVec3ds::iterator		        KdVec3dIter;
    typedef KdVec3d* KdVec3dPtr;
    typedef std::vector<KdVec3dPtr>     	KdVec3dPtrs;
    typedef KdVec3dPtrs::iterator   		KdVec3dPtrIter;

    class KdTree_BBox : public BoundingBox
    {
    public:
        KdTree_BBox()
            :BoundingBox()
        {}

        KdTree_BBox(const KdVec3d& minB, const KdVec3d& maxB)
        {
            minBound = minB.vec;
            maxBound = maxB.vec;
        }

        inline double min_coord(int dim)
        {
            return to_double_ptr(minBound)[dim];
        }

        inline double max_coord(int dim)
        {
            return to_double_ptr(maxBound)[dim];
        }

        inline void set_min_coord(int dim, double val)
        {
            to_double_ptr(minBound)[dim] = val;
        }

        inline void set_max_coord(int dim, double val)
        {
            to_double_ptr(maxBound)[dim] = val;
        }

        void update_from_point_pointers(KdVec3dPtrIter begin, KdVec3dPtrIter end)
        {
            if (begin == end)
                return;

            minBound = (**begin).vec;
            maxBound = (**begin).vec;
            for (begin++; begin != end; begin++)
            {
                minBound = _mm256_min_pd(minBound, (**begin).vec);
                maxBound = _mm256_max_pd(maxBound, (**begin).vec);
            }
        }

        double distance_to_point(const KdVec3d& query, KdVec3d& dists)
        {
            __m256d FourZero = _mm256_setzero_pd();
            __m256d Q = query.vec;

            __m256d T1, T2, T3, T4, T5;
            T1 = _mm256_sub_pd(minBound, Q);	// T1 = minBound - Q
            T2 = _mm256_max_pd(T1, FourZero);	// T2 = max{T1, FourZero} 
            T3 = _mm256_mul_pd(T2, T2);			// T3 = T2 * T2

            T1 = _mm256_sub_pd(Q, maxBound);	// T1 = Q - maxBound
            T4 = _mm256_max_pd(T1, FourZero);	// T2 = max{T1, FourZero}
            T5 = _mm256_fmadd_pd(T4, T4, T3);	// T5 = T4 * T4 + T3

            dists.vec = _mm256_max_pd(T2, T4);
            return hsum_double_avx(T5);
        }
    };

    typedef KdTree_BBox     KdBBox;

    class Separator
    {
    public:
        int cut_dim;
        double cut_val;
    public:
        Separator()
        {}
        Separator(int cutdim, double cutval)
            :cut_dim(cutdim), cut_val(cutval)
        {}
    };

    // forward declaration
    class KdTree;

    class PointContainer
    {
    public:
        friend class KdTree;
    private:
        KdVec3dPtrIter m_begin, m_end;
        int m_build_dim;
        KdBBox m_bbox, m_tbox;
    public:
        PointContainer() { }

        PointContainer(KdVec3dPtrIter begin, KdVec3dPtrIter end)
            :m_begin(begin), m_end(end)
        {
            m_bbox = KdBBox(**begin, **begin);
            for (begin++; begin != end; begin++)
            {
                m_bbox += KdBBox(**begin, **begin);
            }
            m_tbox = m_bbox;
            m_build_dim = m_bbox.longest_axis();
        }

        inline std::size_t size()const
        {
            return std::distance(m_begin, m_end);
        }

        inline KdVec3dPtrIter begin()const
        {
            return m_begin;
        }

        inline KdVec3dPtrIter end()const
        {
            return m_end;
        }

        inline const KdBBox& bbox()const
        {
            return m_bbox;
        }

        inline const KdBBox& tbox()const
        {
            return m_tbox;
        }

        void split(PointContainer& c_low, Separator& sep, bool sliding = false)
        {
            c_low.m_bbox = m_bbox;
            m_build_dim = sep.cut_dim;
            c_low.m_build_dim = sep.cut_dim;

            KdVec3dPtrIter it = std::partition(m_begin, m_end, [&](KdVec3dPtr v) {return v->less_on(sep.cut_dim, sep.cut_val); });
            // now [begin,it) are lower and [it,end) are upper
            if (sliding)    // avoid empty lists
            {
                if (it == m_begin)
                {
                    KdVec3dPtrIter min_elt = std::min_element(m_begin, m_end, [&](KdVec3dPtr lhs, KdVec3dPtr rhs) {return lhs->less_on(sep.cut_dim, *rhs); });
                    if (min_elt != it)
                    {
                        std::iter_swap(min_elt, it);
                    }
                    sep.cut_val = (**it)[sep.cut_dim];
                    it++;
                }
                if (it == m_end)
                {
                    KdVec3dPtrIter max_elt = std::max_element(m_begin, m_end, [&](KdVec3dPtr lhs, KdVec3dPtr rhs) {return lhs->less_on(sep.cut_dim, *rhs); });
                    it--;
                    if (max_elt != it)
                    {
                        std::iter_swap(max_elt, it);
                    }
                    sep.cut_val = (**it)[sep.cut_dim];
                    it++;
                }
            }

            c_low.m_begin = m_begin;
            c_low.m_end = it;
            m_begin = it;
            // adjusting boxes
            m_bbox.set_min_coord(sep.cut_dim, sep.cut_val);
            m_tbox.update_from_point_pointers(m_begin, m_end);
            c_low.m_bbox.set_max_coord(sep.cut_dim, sep.cut_val);
            c_low.m_tbox.update_from_point_pointers(c_low.m_begin, c_low.m_end);
        }
    };

    class KdTree_Node
    {
    public:
        friend class KdTree;
        friend class OrthogonalNearestSeach;
    private:
        bool is_leaf;
    public:
        KdTree_Node(bool leaf)
            :is_leaf(leaf)
        {}
    };

    class KdTree_InternalNode : public KdTree_Node
    {
    public:
        typedef KdTree_Node     Node;
        typedef Node* NodePtr;

        friend class KdTree;
        friend class OrthogonalNearestSeach;
    private:
        int m_cut_dim;
        double m_cut_val;
        double m_lower_low_val, m_lower_high_val, m_upper_low_val, m_upper_high_val;
        NodePtr m_lower_ch, m_upper_ch;

    public:
        KdTree_InternalNode()
            :KdTree_Node(false)
        {}

        void set_separator(const Separator& sep)
        {
            m_cut_dim = sep.cut_dim;
            m_cut_val = sep.cut_val;
        }
    };

    class KdTree_LeafNode : public KdTree_Node
    {
    public:
        friend class KdTree;
        friend class OrthogonalNearestSeach;
    public:
        int n;
        KdVec3dIter data;
    public:
        KdTree_LeafNode()
            :KdTree_Node(true)
        {}

        KdTree_LeafNode(int n_)
            :KdTree_Node(true), n(n_)
        {}
    };

    class KdTree
    {
    public:
        typedef KdTree_Node             Node;
        typedef Node* NodePtr;
        typedef KdTree_InternalNode     InternalNode;
        typedef InternalNode* InternalPtr;
        typedef KdTree_LeafNode         LeafNode;
        typedef LeafNode* LeafPtr;
    public:
        std::deque<InternalNode>    internal_nodes;
        std::deque<LeafNode>        leaf_nodes;

        NodePtr                 tree_root;
        KdBBox                  bbox;
        KdVec3ds                pts;
        KdVec3dPtrs             data;

        int bucket_size = 10;
    public:
        KdTree(const Vec3ds& points, const std::vector<TriIter>& tri_iters)
        {
            insert(points, tri_iters);
            build();
        }


        void insert(const Vec3ds& points, const std::vector<TriIter>& tri_iters)
        {
            pts.reserve(points.size());
            for (int i = 0; i < points.size(); i++)
            {
                pts.emplace_back(points[i], tri_iters[i]);
            }
        }

        void build()
        {
            data.reserve(pts.size());
            for (int i = 0; i < pts.size(); i++)
            {
                data.emplace_back(&pts[i]);
            }

            PointContainer c(data.begin(), data.end());
            bbox = c.bbox();
            if (c.size() <= bucket_size)
            {
                tree_root = create_leaf_node(c);
            }
            else
            {
                tree_root = new_internal_node();
                create_internal_node(tree_root, c);
            }

            KdVec3ds pts_tmp;
            pts_tmp.reserve(pts.size());
            for (int i = 0; i < pts.size(); i++)
            {
                pts_tmp.emplace_back(*data[i]);
            }
            for (int i = 0; i < leaf_nodes.size(); i++)
            {
                std::ptrdiff_t tmp = leaf_nodes[i].data - pts.begin();
                leaf_nodes[i].data = pts_tmp.begin() + tmp;
            }
            pts.swap(pts_tmp);

            data.clear();
        }

        inline bool empty()
        {
            return pts.empty();
        }

        std::pair<Vec3d, TriIter> search_nearest_point(const Vec3d& query)
        {
            OrthogonalNearestSeach ons(this, KdVec3d(query));
            KdVec3dIter cp = ons.closest_point_iter();
            return std::pair<Vec3d, TriIter>(Vec3d(cp->vec), cp->m_from);
        }

    private:
        NodePtr create_leaf_node(PointContainer& c)
        {
            LeafNode node(c.size());
            std::ptrdiff_t tmp = c.begin() - data.begin();
            node.data = pts.begin() + tmp;

            leaf_nodes.emplace_back(node);
            return &(leaf_nodes.back());
        }

        NodePtr new_internal_node()
        {
            internal_nodes.emplace_back();
            return &(internal_nodes.back());
        }

        void create_internal_node(NodePtr n, PointContainer& c)
        {
            InternalPtr nh = static_cast<InternalPtr>(n);

            Separator sep;
            PointContainer c_low;

            split(sep, c, c_low);
            nh->set_separator(sep);

            handle_extended_node(nh, c, c_low);

            int low_size = c_low.size();
            if (c_low.size() > bucket_size)
            {
                nh->m_lower_ch = new_internal_node();
                create_internal_node(nh->m_lower_ch, c_low);
            }
            else
            {
                nh->m_lower_ch = create_leaf_node(c_low);
            }

            int high_size = c.size();
            if (c.size() > bucket_size)
            {
                nh->m_upper_ch = new_internal_node();
                create_internal_node(nh->m_upper_ch, c);
            }
            else
            {
                nh->m_upper_ch = create_leaf_node(c);
            }
        }

        void split(Separator& sep, PointContainer& c_origin, PointContainer& c_low)
        {
            int cutdim = c_origin.m_bbox.longest_axis();

            // Fix degenerated cases
            if (c_origin.m_tbox.min_coord(cutdim) != c_origin.m_tbox.max_coord(cutdim))
            {
                sep = Separator(cutdim,
                    (c_origin.m_bbox.max_coord(cutdim) + c_origin.m_bbox.min_coord(cutdim)) / 2.0);
            }
            else
            {
                cutdim = c_origin.m_tbox.longest_axis();
                sep = Separator(cutdim,
                    (c_origin.m_tbox.max_coord(cutdim) + c_origin.m_tbox.min_coord(cutdim)) / 2.0);
            }

            double max_span_upper = c_origin.m_tbox.max_coord(cutdim);
            double max_span_lower = c_origin.m_tbox.min_coord(cutdim);
            if (max_span_upper <= sep.cut_val)
            {
                sep.cut_val = max_span_upper;
            }
            if (max_span_lower >= sep.cut_val)
            {
                sep.cut_val = max_span_lower;
            }

            c_origin.split(c_low, sep, true);
        }

        void handle_extended_node(InternalPtr nh, PointContainer& c, PointContainer& c_low)
        {
            int cut_dim = nh->m_cut_dim;
            if (c_low.size() > 0)
            {
                nh->m_lower_low_val = c_low.m_tbox.min_coord(cut_dim);
                nh->m_lower_high_val = c_low.m_tbox.max_coord(cut_dim);
            }
            else
            {
                nh->m_lower_low_val = nh->m_cut_val;
                nh->m_lower_high_val = nh->m_cut_val;
            }
            if (c.size() > 0)
            {
                nh->m_upper_low_val = c.m_tbox.min_coord(cut_dim);
                nh->m_upper_high_val = c.m_tbox.max_coord(cut_dim);
            }
            else
            {
                nh->m_upper_low_val = nh->m_cut_val;
                nh->m_upper_high_val = nh->m_cut_val;
            }
        }
    private:
        class OrthogonalNearestSeach
        {
        private:
            KdVec3d m_dists;
            KdVec3d m_query;
            KdTree* m_tree;

            KdVec3dIter m_result;
            double m_square_distance;
        public:
            OrthogonalNearestSeach(KdTree* tree, const KdVec3d& query)
                :m_tree(tree), m_query(query), m_dists(0, 0, 0), m_square_distance(DBL_MAX)
            {
                if (m_tree->empty()) return;

                double distance_to_root = m_tree->bbox.distance_to_point(m_query, m_dists);
                compute_nearest_neightbors_orthogonally(m_tree->tree_root, distance_to_root);
            }

            KdVec3dIter closest_point_iter()
            {
                return m_result;
            }
        private:
            void compute_nearest_neightbors_orthogonally(NodePtr N, double rd)
            {
                if (N->is_leaf)
                {
                    LeafPtr node = static_cast<LeafPtr>(N);
                    if (node->n > 0)
                    {
                        search_nearest_in_leaf(node);
                    }
                }
                else
                {
                    InternalPtr node = static_cast<InternalPtr>(N);
                    int cut_dim = node->m_cut_dim;

                    NodePtr best_ch, other_ch;

                    double offset;
                    double val = m_query[cut_dim];
                    double diff1 = val - node->m_upper_low_val;
                    double diff2 = val - node->m_lower_high_val;

                    if (diff1 + diff2 < 0)
                    {
                        offset = diff1;
                        best_ch = node->m_lower_ch;
                        other_ch = node->m_upper_ch;
                    }
                    else
                    {
                        offset = diff2;
                        best_ch = node->m_upper_ch;
                        other_ch = node->m_lower_ch;
                    }

                    compute_nearest_neightbors_orthogonally(best_ch, rd);
                    double dst = m_dists[cut_dim];
                    double new_rd = new_square_distance(rd, dst, offset);
                    m_dists[cut_dim] = offset;
                    if (new_rd < m_square_distance)
                    {
                        compute_nearest_neightbors_orthogonally(other_ch, new_rd);
                    }
                    m_dists[cut_dim] = dst;
                }
            }

            void search_nearest_in_leaf(LeafPtr node)
            {
                for (KdVec3dIter begin = node->data, end = node->data + node->n; begin != end; begin++)
                {
                    double distance = ((*begin) - m_query).squaredNorm();
                    if (distance < m_square_distance)
                    {
                        m_result = begin;
                        m_square_distance = distance;
                    }
                }
            }

            double new_square_distance(double dist, double old_off, double new_off)
            {
                return dist + (new_off * new_off - old_off * old_off);
            }
        };
    };
#endif  // USE_AVX
}
