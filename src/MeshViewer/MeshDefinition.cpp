#include "MeshDefinition.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

double calc_mesh_ave_edge_length(Mesh* mesh_)
{
	double ave_len = 0.0; int ne = mesh_->n_edges();
	for(Mesh::EdgeIter e_it = mesh_->edges_begin(); e_it != mesh_->edges_end(); ++e_it)
	{
		Mesh::HalfedgeHandle heh = mesh_->halfedge_handle(e_it, 0);
		ave_len += (mesh_->point(mesh_->from_vertex_handle(heh)) - mesh_->point(mesh_->to_vertex_handle(heh)) ).norm();
	}

	ave_len /= ne; double s = 1.0 / ave_len;

	/*for(Mesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d p = mesh_->point(v_it);
		mesh_->set_point(v_it, p*s);
	}*/

	return (ave_len);
}

bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_)
{
	// boundary edges cannot be flipped
	if ( mesh_.is_boundary(eh) ) return false;

	Mesh::HalfedgeHandle hh = mesh_.halfedge_handle(eh, 0);
	Mesh::HalfedgeHandle oh = mesh_.halfedge_handle(eh, 1);

	// check if the flipped edge is already present
	// in the mesh

	Mesh::VertexHandle ah = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(hh));
	Mesh::VertexHandle bh = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(oh));

	if(ah == bh)   // this is generally a bad sign !!!
		return false;

	for (Mesh::ConstVertexVertexIter vvi = mesh_.vv_iter(ah); vvi; ++vvi)
		if( vvi.handle() == bh )
			return false;

	return true;
}

bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_)
{
	// CAUTION : Flipping a halfedge may result in
	// a non-manifold mesh, hence check for yourself
	// whether this operation is allowed or not!
	if( !is_flip_ok_openmesh(eh, mesh_) )
		return false;//let's make it sure it is actually checked
	//assert( is_flip_ok_openmesh(eh, mesh_ ) );

	Mesh::HalfedgeHandle a0 = mesh_.halfedge_handle(eh, 0);
	Mesh::HalfedgeHandle b0 = mesh_.halfedge_handle(eh, 1);

	Mesh::HalfedgeHandle a1 = mesh_.next_halfedge_handle(a0);
	Mesh::HalfedgeHandle a2 = mesh_.next_halfedge_handle(a1);

	Mesh::HalfedgeHandle b1 = mesh_.next_halfedge_handle(b0);
	Mesh::HalfedgeHandle b2 = mesh_.next_halfedge_handle(b1);

	Mesh::VertexHandle   va0 = mesh_.to_vertex_handle(a0);
	Mesh::VertexHandle   va1 = mesh_.to_vertex_handle(a1);

	Mesh::VertexHandle   vb0 = mesh_.to_vertex_handle(b0);
	Mesh::VertexHandle   vb1 = mesh_.to_vertex_handle(b1);

	Mesh::FaceHandle     fa  = mesh_.face_handle(a0);
	Mesh::FaceHandle     fb  = mesh_.face_handle(b0);

	mesh_.set_vertex_handle(a0, va1);
	mesh_.set_vertex_handle(b0, vb1);

	mesh_.set_next_halfedge_handle(a0, a2);
	mesh_.set_next_halfedge_handle(a2, b1);
	mesh_.set_next_halfedge_handle(b1, a0);

	mesh_.set_next_halfedge_handle(b0, b2);
	mesh_.set_next_halfedge_handle(b2, a1);
	mesh_.set_next_halfedge_handle(a1, b0);

	mesh_.set_face_handle(a1, fb);
	mesh_.set_face_handle(b1, fa);

	mesh_.set_halfedge_handle(fa, a0);
	mesh_.set_halfedge_handle(fb, b0);

	if(mesh_.halfedge_handle(va0) == b0)
		mesh_.set_halfedge_handle(va0, a1);
	if(mesh_.halfedge_handle(vb0) == a0)
		mesh_.set_halfedge_handle(vb0, b1);

	return true;
}

bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p)
{
	OpenMesh::Vec3d v1 = tri[1] - tri[0]; OpenMesh::Vec3d v2 = tri[2] - tri[0];
	OpenMesh::Vec3d n = OpenMesh::cross(v1, v2);
	double face_area = n.norm(); n.normalize(); double all_area = 0;
	for(unsigned i=0; i < tri.size(); ++i)
	{
		unsigned next_i = ( i+1 )%tri.size(); unsigned prev_i = ( i + tri.size() - 1 )%tri.size();
		v1 = tri[next_i] - p; v2 = tri[prev_i] - p;
		double area = OpenMesh::dot(OpenMesh::cross(v1, v2), n); all_area += area;
		if(area < 0)
		{
			return false;
		}
	}
	if(std::abs(all_area - face_area) < 1e-8) {return true;}
	else {return false;}
}

bool baryCoord( const OpenMesh::Vec3d& _p, const OpenMesh::Vec3d& _u, const OpenMesh::Vec3d& _v, const OpenMesh::Vec3d& _w, OpenMesh::Vec3d&_result )
{
	Mesh::Point  vu = _v - _u, wu = _w - _u, pu = _p - _u;

	// find largest absolute coodinate of normal
	Mesh::Scalar nx = vu[1]*wu[2] - vu[2]*wu[1],
		ny = vu[2]*wu[0] - vu[0]*wu[2],
		nz = vu[0]*wu[1] - vu[1]*wu[0],
		ax = fabs(nx),
		ay = fabs(ny),
		az = fabs(nz);


	unsigned char max_coord;

	if ( ax > ay )
	{
		if ( ax > az ) 
		{
			max_coord = 0;
		}
		else
		{
			max_coord = 2;
		}
	}
	else
	{
		if ( ay > az )
		{
			max_coord = 1;
		}
		else 
		{
			max_coord = 2;
		}
	}


	// solve 2D problem
	switch (max_coord)
	{
	case 0:
		{      
			if (1.0+ax == 1.0) return false;
			_result[1] = 1.0 + (pu[1]*wu[2]-pu[2]*wu[1])/nx - 1.0;
			_result[2] = 1.0 + (vu[1]*pu[2]-vu[2]*pu[1])/nx - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;

	case 1:
		{
			if (1.0+ay == 1.0) return false;
			_result[1] = 1.0 + (pu[2]*wu[0]-pu[0]*wu[2])/ny - 1.0;
			_result[2] = 1.0 + (vu[2]*pu[0]-vu[0]*pu[2])/ny - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;

	case 2:
		{
			if (1.0+az == 1.0) return false;
			_result[1] = 1.0 + (pu[0]*wu[1]-pu[1]*wu[0])/nz - 1.0;
			_result[2] = 1.0 + (vu[0]*pu[1]-vu[1]*pu[0])/nz - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;
	}

	return true;
}

void compute_point_area(Mesh* mesh_, std::vector<std::map<int,double>>& cornerArea, std::vector<double>& pointArea, bool use_np)
{
	pointArea.resize(mesh_->n_vertices(),0.0);
	cornerArea.resize(mesh_->n_faces());

	Mesh::FaceIter f_it = mesh_->faces_begin();
	int temp_face_index = 0;
	Mesh::FaceHalfedgeIter fhe_it;
	Mesh::Point e[3];
	int v[3];

	//#pragma omp parallel for
	for( f_it; f_it != mesh_->faces_end(); ++f_it)
	{
		temp_face_index = f_it.handle().idx();

		fhe_it = mesh_->fh_iter(f_it);
		if( !use_np )
		{
			e[0] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		}
		else
		{
			e[0] = mesh_->data( mesh_->to_vertex_handle(fhe_it) ).get_New_Pos() - mesh_->data( mesh_->from_vertex_handle(fhe_it) ).get_New_Pos();
		}
		v[2] = mesh_->to_vertex_handle(fhe_it).idx();
		++fhe_it;
		if( !use_np )
		{
			e[1] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		}
		else
		{
			e[1] = mesh_->data( mesh_->to_vertex_handle(fhe_it) ).get_New_Pos() - mesh_->data( mesh_->from_vertex_handle(fhe_it) ).get_New_Pos();
		}
		v[0] = mesh_->to_vertex_handle(fhe_it).idx();
		++fhe_it;
		if( !use_np )
		{
			e[2] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		}
		else
		{
			e[2] = mesh_->data( mesh_->to_vertex_handle(fhe_it) ).get_New_Pos() - mesh_->data( mesh_->from_vertex_handle(fhe_it) ).get_New_Pos();
		}
		v[1] = mesh_->to_vertex_handle(fhe_it).idx();

		double area = 0.5f * ( e[0]%e[1] ).norm();
		double l2[3] = { e[0].sqrnorm(), e[1].sqrnorm(), e[2].sqrnorm() };
		double ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
			l2[1] * (l2[2] + l2[0] - l2[1]),
			l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (ew[0] <= 0.0) 
		{
			cornerArea[temp_face_index][v[1]] = -0.25 * l2[2] * area / OpenMesh::dot(e[0],e[2]);
			cornerArea[temp_face_index][v[2]] = -0.25 * l2[1] * area / OpenMesh::dot(e[0],e[1]);
			cornerArea[temp_face_index][v[0]] = area - cornerArea[temp_face_index][v[1]] -cornerArea[temp_face_index][v[2]];
		}
		else if (ew[1] <= 0.0)
		{
			cornerArea[temp_face_index][v[2]] = -0.25 * l2[0] * area / OpenMesh::dot(e[1],e[0]);
			cornerArea[temp_face_index][v[0]] = -0.25 * l2[2] * area / OpenMesh::dot(e[1],e[2]);
			cornerArea[temp_face_index][v[1]] = area - cornerArea[temp_face_index][v[2]] - cornerArea[temp_face_index][v[0]];
		}
		else if (ew[2] <= 0.0f)
		{
			cornerArea[temp_face_index][v[0]] = -0.25 * l2[1] * area / OpenMesh::dot(e[2],e[1]);
			cornerArea[temp_face_index][v[1]] = -0.25 * l2[0] * area / OpenMesh::dot(e[2],e[0]);
			cornerArea[temp_face_index][v[2]] = area - cornerArea[temp_face_index][v[0]] - cornerArea[temp_face_index][v[1]];
		}
		else 
		{
			double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
			{
				cornerArea[temp_face_index][v[j]] = ewscale * (ew[(j+1)%3] + ew[(j+2)%3]);
			}
		}

		//#pragma omp atomic
		pointArea[v[0]] += cornerArea[temp_face_index][v[0]];
		//#pragma omp atomic
		pointArea[v[1]] += cornerArea[temp_face_index][v[1]];
		//#pragma omp atomic
		pointArea[v[2]] += cornerArea[temp_face_index][v[2]];
	}
}

void rot_coord_sys(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, const OpenMesh::Vec3d &new_norm, OpenMesh::Vec3d &new_u, OpenMesh::Vec3d &new_v)
{
	new_u = old_u;
	new_v = old_v;
	OpenMesh::Vec3d old_norm = OpenMesh::cross( old_u , old_v );
	double ndot = OpenMesh::dot( old_norm , new_norm);
	if ( ndot <= -1.0 )
	{
		new_u = -new_u;
		new_v = -new_v;
		return;
	}
	OpenMesh::Vec3d perp_old = new_norm - ndot * old_norm;
	OpenMesh::Vec3d dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
	new_u -= dperp * OpenMesh::dot(new_u , perp_old);
	new_v -= dperp * OpenMesh::dot(new_v , perp_old);
}

void proj_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, double old_ku, double old_kuv, double old_kv, const OpenMesh::Vec3d &new_u, const OpenMesh::Vec3d &new_v, double &new_ku, double &new_kuv, double &new_kv)
{
	OpenMesh::Vec3d r_new_u; OpenMesh::Vec3d r_new_v;
	rot_coord_sys(new_u, new_v, OpenMesh::cross(old_u, old_v), r_new_u, r_new_v);

	double u1 = OpenMesh::dot( r_new_u , old_u);
	double v1 = OpenMesh::dot( r_new_u , old_v);
	double u2 = OpenMesh::dot( r_new_v , old_u);
	double v2 = OpenMesh::dot( r_new_v , old_v);
	new_ku  = old_ku * u1*u1 + old_kuv * (2.0f  * u1*v1) + old_kv * v1*v1;
	new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
	new_kv  = old_ku * u2*u2 + old_kuv * (2.0f  * u2*v2) + old_kv * v2*v2;
}

void diagonalize_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, double ku, double kuv, double kv, const OpenMesh::Vec3d &new_norm, OpenMesh::Vec3d &pdir1, OpenMesh::Vec3d &pdir2, double &vk1, double &vk2)
{
	OpenMesh::Vec3d r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	double c = 1.0, s = 0.0, tt = 0.0;
	if(kuv != 0.0) 
	{
		// Jacobi rotation to diagonalize
		double h = 0.5 * (kv - ku) / kuv;
		tt = (h < 0.0) ?
			1.0 / (h - std::sqrt(1.0 + h*h)) :
		1.0 / (h + std::sqrt(1.0 + h*h));
		c = 1.0 / std::sqrt(1.0 + tt*tt);
		s = tt * c;
	}

	vk1 = ku - tt * kuv;
	vk2 = kv + tt * kuv;

	if (std::abs(vk1) >= std::abs(vk2))
	{
		pdir1 = c*r_old_u - s*r_old_v;
	} 
	else
	{
		std::swap(vk1, vk2);
		pdir1 = s*r_old_u + c*r_old_v;
	}
	pdir2 = OpenMesh::cross( new_norm , pdir1);
}

void compute_principal_curvature(Mesh* mesh_, 
								 std::vector<double>& K1, std::vector<double>& K2, 
								 std::vector<OpenMesh::Vec3d>& dir1,std::vector<OpenMesh::Vec3d>& dir2)
{
	if(!mesh_->has_vertex_normals())
	{
		mesh_->request_vertex_normals();
		//mesh_->update_vertex_normals();
		mesh_->update_normals();
	}
	int nv = mesh_->n_vertices();
	int nf = mesh_->n_faces();
	Mesh::FaceIter f_it;

	//compute vertex normal
	std::vector<OpenMesh::Vec3d> vertex_normal(nv,OpenMesh::Vec3d(0.0,0.0,0.0));
	Mesh::FaceVertexIter fv_it;
	OpenMesh::Vec3d v[3];
	int v_index[3];
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		fv_it = mesh_->fv_iter(f_it);
		int count = 0;
		for(fv_it;fv_it;++fv_it)
		{
			v[count] = mesh_->point(fv_it);
			v_index[count] = fv_it.handle().idx();
			++count;
		}
		OpenMesh::Vec3d a = v[0] - v[1];
		OpenMesh::Vec3d b = v[1] - v[2];
		OpenMesh::Vec3d c = v[2] - v[0];
		double l2a = a.sqrnorm(); double l2b = b.sqrnorm(); double l2c = c.sqrnorm();
		if (!l2a || !l2b || !l2c)
			continue;

		OpenMesh::Vec3d temp_fn = OpenMesh::cross(a,b);
		vertex_normal[v_index[0]] += temp_fn * (1.0 / (l2a * l2c));
		vertex_normal[v_index[1]] += temp_fn * (1.0 / (l2b * l2a));
		vertex_normal[v_index[2]] += temp_fn * (1.0 / (l2c * l2b));
	}

	K1.clear(); K1.resize(nv,0.0); K2.clear(); K2.resize(nv,0.0);
	dir1.clear(); dir1.resize(nv,OpenMesh::Vec3d(0,0,0));
	dir2.clear(); dir2.resize(nv,OpenMesh::Vec3d(0,0,0));
	std::vector<double> k12(nv,0.0);

	Mesh::FaceHalfedgeIter fhe_it;
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		fhe_it = mesh_->fh_iter(f_it);
		for(fhe_it;fhe_it;++fhe_it)
		{
			dir1[mesh_->from_vertex_handle(fhe_it).idx()] =
				mesh_->point( mesh_->to_vertex_handle(fhe_it) ) 
				- mesh_->point( mesh_->from_vertex_handle(fhe_it));
		}
	}

	Mesh::VertexIter v_it;
	int vertex_id; OpenMesh::Vec3d vn;
	for(v_it = mesh_->vertices_begin();v_it != mesh_->vertices_end(); ++v_it)
	{
		//vn = mesh_->normal(v_it); vn.normalize();
		vertex_normal[v_it.handle().idx()].normalize();
		vn = vertex_normal[v_it.handle().idx()];
		vertex_id = v_it.handle().idx();
		dir1[ vertex_id ] = OpenMesh::cross(dir1[ vertex_id ], vn);
		dir1[ vertex_id ].normalize();
		dir2[ vertex_id ] = OpenMesh::cross(vn, dir1[ vertex_id ]);
		dir2[ vertex_id ].normalize();
	}

	std::vector<std::map<int,double>> cornerArea;
	std::vector<double> pointArea;
	compute_point_area(mesh_, cornerArea, pointArea);

	OpenMesh::Vec3d ev[3]; 
	Mesh::HalfedgeHandle heh; Mesh::HalfedgeHandle temp_heh; 
	OpenMesh::Vec3d t; OpenMesh::Vec3d n;OpenMesh::Vec3d b;
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		// Edges
		fhe_it = mesh_->fh_iter(f_it); heh = fhe_it.handle();
		temp_heh = mesh_->next_halfedge_handle(heh);
		ev[0] = mesh_->point(mesh_->to_vertex_handle(temp_heh)) - mesh_->point(mesh_->from_vertex_handle(temp_heh));
		temp_heh = mesh_->prev_halfedge_handle(heh);
		ev[1] = mesh_->point(mesh_->to_vertex_handle(temp_heh)) - mesh_->point(mesh_->from_vertex_handle(temp_heh));
		ev[2] = mesh_->point(mesh_->to_vertex_handle(     heh)) - mesh_->point(mesh_->from_vertex_handle(     heh));
		
		// N-T-B coordinate system per face
		t = ev[0]; t.normalize();
		n = OpenMesh::cross(ev[0],ev[1]);
		b = OpenMesh::cross(n,t); b.normalize();

		// Estimate curvature based on variation of normals
		// along edges
		temp_heh = mesh_->next_halfedge_handle(heh);
		Eigen::Vector3d m(0.0,0.0,0.0);
		Eigen::Vector3d x;
		Eigen::Matrix3d w; w.setZero();
		for(int j=0;j<3;++j)
		{
			double u = OpenMesh::dot(ev[j],t);
			double v = OpenMesh::dot(ev[j],b);
			w(0,0) += u*u;
			w(0,1) += u*v;
			//w[1][1] += v*v + u*u; 
			//w[1][2] += u*v; 
			w(2,2) += v*v;
			/*OpenMesh::Vec3d dn = mesh_->normal(mesh_->to_vertex_handle(temp_heh)) 
			- mesh_->normal(mesh_->from_vertex_handle(temp_heh));*/
			OpenMesh::Vec3d dn = vertex_normal[mesh_->to_vertex_handle(temp_heh).idx()] 
							   - vertex_normal[mesh_->from_vertex_handle(temp_heh).idx()];
			double dnu = OpenMesh::dot(dn, t);
			double dnv = OpenMesh::dot(dn, b);
			m(0) += dnu*u;
			m(1) += dnu*v + dnv*u;
			m(2) += dnv*v;
			temp_heh = mesh_->next_halfedge_handle(temp_heh);
		}
		w(1,1) = w(0,0) + w(2,2);
		w(1,2) = w(0,1);
		w(2,1) = w(1,2);
		w(1,0) = w(0,1);
		//std::cout << w;
		x = w.fullPivHouseholderQr().solve(m);

		temp_heh = heh;
		for(int j=0;j<3;++j)
		{
			vertex_id = mesh_->from_vertex_handle(temp_heh).idx();
			double c1, c12, c2;
			proj_curv(t, b, x(0), x(1), x(2),
				dir1[vertex_id], dir2[vertex_id], c1, c12, c2);
			double wt = cornerArea[f_it.handle().idx()][vertex_id] / pointArea[vertex_id];
			K1[vertex_id]  += wt * c1;
			k12[vertex_id] += wt * c12;
			K2[vertex_id] += wt * c2;
			temp_heh = mesh_->next_halfedge_handle(temp_heh);
		}
	}

	for(v_it = mesh_->vertices_begin();v_it != mesh_->vertices_end(); ++v_it)
	{
		vertex_id = v_it.handle().idx();
		/*diagonalize_curv(dir1[vertex_id], dir2[vertex_id],
		k1[vertex_id], k12[vertex_id], k2[vertex_id],
		mesh_->normal(v_it), dir1[vertex_id], dir2[vertex_id],
		k1[vertex_id], k2[vertex_id]);*/
		diagonalize_curv(dir1[vertex_id], dir2[vertex_id],
			K1[vertex_id], k12[vertex_id], K2[vertex_id],
			vertex_normal[vertex_id], dir1[vertex_id], dir2[vertex_id],
			K1[vertex_id], K2[vertex_id]);
	}
}


#pragma region functions by yanyisheshou at GCL
#include "..\src\Toolbox\dprint.h"
double meshMinAngle(TriMesh &mesh)
{
	double angle = 4;
	for (auto th : mesh.halfedges())
	{
		angle = std::min(angle, mesh.calc_sector_angle(th));
	}
	return angle;
}

//the definition of single triangle quality in trimesh is from paper: Automatic and High-quality Surface Mesh Generation for CAD Models(Section 5.1)
void printMeshQuality(TriMesh &mesh)
{
	double c = 4 * std::sqrt(3);
	double minAngle = 4, maxAngle = 0, avgAngle = 0;
	double minQuality = 1, avgQuality = 0;
	for (auto tf : mesh.faces())
	{
		double h = 0, p = 0;
		for (auto &tfh : mesh.fh_range(tf))
		{
			double angle = mesh.calc_sector_angle(tfh);
			minAngle = std::min(angle, minAngle);
			maxAngle = std::max(angle, maxAngle);
			avgAngle += angle;

			double l = mesh.calc_edge_length(tfh);
			h = std::max(l, h);
			p += l;
		}
		double q = c * mesh.calc_face_area(tf) / (p*h);
		minQuality = std::min(q, minQuality);
		avgQuality += q;
	}
	dprint("mesh quality info:\nmin, max and avg angle:", minAngle, maxAngle, avgAngle / mesh.n_halfedges(),
		"\nmin and avg quality:", minQuality, avgQuality / mesh.n_faces());
}
#pragma endregion
