#include "crossField.h"
void crossField::calculateMeshFaceBase()
{
	mesh->update_face_normals();
	//faceBase.reserve(mesh->n_faces() * 2);
	faceBase.resize(3, mesh->n_faces() * 2);
	for (auto& tf : mesh->faces())
	{
		/*std::vector<OpenMesh::Vec3d> v;
		v.reserve(3);
		for (auto& tfv : mesh->fv_range(tf))
		{
			v.push_back(mesh->point(tfv));
		}*/
		/*faceBase.push_back((v[1] - v[0]).normalized());
		faceBase.push_back(cross(faceBase.back(), cross(v[2] - v[0], v[1] - v[0])).normalized());*/
		Eigen::Matrix3d v;
		int i = 0;
		for (auto &tfv : mesh->fv_range(tf))
		{
			v.col(i) = position.col(tfv.idx());
			++i;
		}
		i = tf.idx() * 2;
		faceBase.col(i) = (v.col(1) - v.col(0)).normalized();
		faceBase.col(i + 1) = faceBase.col(i).cross((v.col(2) - v.col(0)).cross(v.col(1) - v.col(0))).normalized();
	}
}

void crossField::runPolynomial()
{
	setCurvatureConstraint();

	using namespace std;

	auto fnum = mesh->n_faces();
	vector<int> status(fnum, 0);

	std::vector<COMPLEX> f_dir(fnum);

	/*for (int i = 0; i < constraintId.size(); i++)
	{
		int fid = constraintId[i];
		status[fid] = 1;
		OpenMesh::Vec3d cf = constraintVector[i].normalized();
		f_dir[fid] = std::pow(COMPLEX(dot(cf, faceBase[fid * 2]), dot(cf, faceBase[fid * 2 + 1])), 4);
	}*/
	for (int i = 0; i < constraintId.size(); ++i)
	{
		int fid = constraintId[i];
		status[fid] = 1;
		constraintVector.col(i).normalize();
		f_dir[fid] = std::pow(COMPLEX(constraintVector.col(i).dot(faceBase.col(fid * 2)), constraintVector.col(i).dot(faceBase.col(fid * 2 + 1))), 4);
	}
	vector<int> id2sln(fnum, -1);
	vector<int> sln2id(0);
	int count = 0;
	for (int i = 0; i < fnum; i++)
	{
		if (status[i] == 0)
		{
			sln2id.push_back(i);
			id2sln[i] = count;
			count++;
		}
	}

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
	Eigen::SparseMatrix<COMPLEX> A;
	Eigen::VectorXcd b_pre(mesh->n_edges());
	Eigen::VectorXcd b;
	b_pre.setZero();
	vector<Eigen::Triplet<COMPLEX>> tris;

	count = 0;
	for (auto &tf : mesh->faces())
	{
		int id_f = tf.idx();
		for (auto &tfh : mesh->fh_range(tf))
		{
			if (!mesh->is_boundary(tfh.edge()))
			{
				/*auto p1 = mesh->point(tfh.to());
				auto p2 = mesh->point(tfh.from());*/
				auto &p1 = position.col(tfh.to().idx());
				auto &p2 = position.col(tfh.from().idx());
				auto id_g = mesh->face_handle(tfh.opp()).idx();
				if (id_f < id_g)
				{
					if (status[id_f] == 1 && status[id_g] == 1)continue;
					auto e = (p2 - p1).normalized();

					/*COMPLEX e_f = COMPLEX(dot(e, faceBase[id_f * 2]), dot(e, faceBase[id_f * 2 + 1]));
					COMPLEX e_g = COMPLEX(dot(e, faceBase[id_g * 2]), dot(e, faceBase[id_g * 2 + 1]));*/
					COMPLEX e_f = COMPLEX(e.dot(faceBase.col(id_f * 2)), e.dot(faceBase.col(id_f * 2 + 1)));
					COMPLEX e_g = COMPLEX(e.dot(faceBase.col(id_g * 2)), e.dot(faceBase.col(id_g * 2 + 1)));

					COMPLEX e_f_c_4 = pow(conj(e_f), 4);
					COMPLEX e_g_c_4 = pow(conj(e_g), 4);

					if (status[id_f] == 0)
					{
						tris.emplace_back(count, id2sln[id_f], e_f_c_4);
					}
					else
					{
						b_pre[count] += -e_f_c_4 * f_dir[id_f];
					}
					if (status[id_g] == 0)
					{
						tris.emplace_back(count, id2sln[id_g], -e_g_c_4);
					}
					else
					{
						b_pre[count] += e_g_c_4 * f_dir[id_g];
					}

					count++;
				}
			}
		}

	}
	A.resize(count, sln2id.size());
	b.resize(count);
	b = b_pre.head(count);
	A.setFromTriplets(tris.begin(), tris.end());
	Eigen::SparseMatrix<COMPLEX> AT = A.adjoint();
	slu.compute(AT * A);
	Eigen::VectorXcd x = slu.solve(AT * b);
#if 0
	for (int i = 0; i < x.size(); i++)
	{
		dprint(std::sqrt(x(i).real()*x(i).real() + x(i).imag()*x(i).imag()));
	}
#endif
	for (int i = 0; i < fnum; i++)
	{
		if (status[i] == 0)
		{
			f_dir[i] = x(id2sln[i]);
		}
	}
#if 0
	ec4.resize(mesh->n_halfedges());
	for (auto &te : mesh->edges())
	{
		auto ev = position.col(te.h0().to().idx()) - position.col(te.h0().from().idx());
		auto fid = mesh->face_handle(te.h0()).idx() * 2;
		auto gid = mesh->face_handle(te.h1()).idx() * 2;
		auto ec0 = COMPLEX(faceBase.col(fid).dot(ev), faceBase.col(fid + 1).dot(ev)); ec0 *= ec0; ec0 *= ec0;
		auto ec1 = COMPLEX(faceBase.col(gid).dot(ev), faceBase.col(gid + 1).dot(ev)); ec1 *= ec1; ec1 *= ec1;
		ec4[te.h0().idx()] = ec0 / ec1;
		ec4[te.h1().idx()] = ec1 / ec0;
	}
	for (int i = 0; i < 100; i++)
	{
		std::vector<std::complex<double>> field4(mesh->n_faces(), std::complex<double>(0.0, 0.0));

		for (auto tf : mesh->faces())
		{
			for (auto tfh : mesh->fh_range(tf))
			{
				auto tg = mesh->opposite_face_handle(tfh);
				if (!tg.is_valid()) continue;

				//field4[f0.idx()] += field_complex4[f1.idx()] * e_w4[fh_h.idx()];
				field4[tf.idx()] += f_dir[tg.idx()] * ec4[tfh.idx()];
			}

			//field4[f0.idx()] += field_complex4[f0.idx()];
			//field4[f0.idx()] /= std::sqrt(std::norm(field4[f0.idx()]));
			field4[tf.idx()] += f_dir[tf.idx()];
			field4[tf.idx()] /= std::sqrt(std::norm(field4[tf.idx()]));
		}

		//field_complex4 = std::move(field4);
		f_dir = std::move(field4);
	}
#endif
	crossfield.resize(3, 4 * fnum);
	for (int i = 0; i < fnum; i++)
	{
		double arg = std::arg(f_dir[i]) * 0.25;
		for (int j = 0; j < 4; ++j)
		{
			crossfield.col(i * 4 + j) = faceBase.col(i * 2)*cos(arg + j * PI*0.5) + faceBase.col(i * 2 + 1)*sin(arg + j * PI*0.5);
		}
	}
	setNormal();
	setMatching();
	setSingularity();
}

//void crossField::runIteration(std::vector<OpenMesh::Vec3d>& crossfield)
//{
//	using namespace std;
//	using namespace Eigen;
//	typedef complex<double> COMPLEX;
//
//	size_t nf = mesh->n_faces();
//	size_t ne = mesh->n_edges();
//
//	SimplicialLDLT<SparseMatrix<double>> solver;
//	SparseMatrix<double> B_CCB, ATA, Hessian;
//	B_CCB.resize(ne * 2, nf * 2);
//	ATA.resize(nf * 2, nf * 2);
//	Hessian.resize(nf * 2, nf * 2);
//	VectorXd x(nf * 2);
//#if 0
//	x.setConstant(1.0);
//#else
//	runPolynomial(crossfield);
//	int i = 0;
//	for (auto dir : f_dir)
//	{
//		double sum = dir.real()*dir.real() + dir.imag()*dir.imag();
//		x(i) = dir.real() / sum;
//		x(i + nf) = dir.imag() / sum;
//		++i;
//	}
//#endif
//	VectorXd b(nf * 2);
//	vector<Triplet<double>> tri;
//
//	for (int itertimes = 0; itertimes < 10; ++itertimes)
//	{
//		dprint("itertimes:", itertimes);
//		int count = 0;
//		tri.reserve(8 * ne);
//		for (auto& th : mesh->halfedges())
//		{
//			if (mesh->is_boundary(th.edge()))
//				continue;
//			size_t id_f = th.face().idx();
//			size_t id_g = th.opp().face().idx();
//			if (id_f < id_g)
//			{
//				auto& p1 = mesh->point(th.to());
//				auto& p2 = mesh->point(th.from());
//				OpenMesh::Vec3d e = (p2 - p1).normalized();
//				COMPLEX ef = COMPLEX(dot(e, faceBase[id_f * 2]), -dot(e, faceBase[id_f * 2 + 1]));
//				ef *= ef; ef *= ef;
//				COMPLEX eg = COMPLEX(dot(e, faceBase[id_g * 2]), -dot(e, faceBase[id_g * 2 + 1]));
//				eg *= eg; eg *= -eg;
//
//				tri.emplace_back(count, id_f, ef.real());
//				tri.emplace_back(count + ne, id_f + nf, ef.real());
//				tri.emplace_back(count, id_f + nf, -ef.imag());
//				tri.emplace_back(count + ne, id_f, ef.imag());
//
//				tri.emplace_back(count, id_g, eg.real());
//				tri.emplace_back(count + ne, id_g + nf, eg.real());
//				tri.emplace_back(count, id_g + nf, -eg.imag());
//				tri.emplace_back(count + ne, id_g, eg.imag());
//
//				++count;
//			}
//		}
//
//		B_CCB.setZero();
//		B_CCB.setFromTriplets(tri.begin(), tri.end());
//		ATA = B_CCB.transpose().eval() * B_CCB;
//		tri.clear();
//		tri.reserve(4 * nf);
//		for (count = 0; count < nf; ++count)
//		{
//			double temp = 1.0 / (x(count) * x(count) + x(count + nf) * x(count + nf));
//			double t1 = temp * temp; double t2 = temp * t1;
//			tri.emplace_back(count, count, 2 * (1 - t1) + 8 * x(count) * x(count) * t2);
//			tri.emplace_back(count, count + nf, 8 * x(count) * x(count + nf) * t2);
//			tri.emplace_back(count + nf, count, 8 * x(count) * x(count + nf) * t2);
//			tri.emplace_back(count + nf, count + nf, 2 * (1 - t1) + 8 * x(count + nf) * x(count + nf) * t2);
//		}
//		Hessian.setZero();
//		Hessian.setFromTriplets(tri.begin(), tri.end());
//		Hessian += ATA;
//
//		for (count = 0; count < nf; ++count)
//		{
//			double temp = x(count) * x(count) + x(count + nf) * x(count + nf);
//			temp = 2 * (1 - 1.0 / temp * temp);
//			b(count) = x(count) * temp;
//			b(count + nf) = x(count + nf) * temp;
//		}
//		b += ATA * x;
//
//		if (!itertimes)
//		{
//			//solver.analyzePattern(Hessian);
//		}
//		solver.compute(Hessian);
//		//solver.factorize(Hessian);
//		x -= solver.solve(b);
//		for (count = 0; count < nf; ++count)
//		{
//			double temp = x(count) * x(count) + x(count + nf) * x(count + nf);
//			if (temp < 1.01)
//			{
//				temp = std::sqrt(temp);
//				x(count) /= temp*0.99;
//				x(count + nf) /= temp*0.99;
//			}
//		}
//	}
//	crossfield.resize(4 * nf);
//	for (int i = 0; i < nf; ++i)
//	{
//		//dprint(i);
//		double arg = std::arg(COMPLEX(x(i), x(i + nf))) * 0.25;
//		for (int j = 0; j < 4; j++)
//		{
//			crossfield[i * 4 + j] = faceBase[i * 2] * cos(arg + j * PI * 0.5) + faceBase[i * 2 + 1] * sin(arg + j * PI * 0.5);
//		}
//	}
//}

void crossField::setCurvatureConstraint()
{
	if (!mesh->has_vertex_normals())
	{
		mesh->request_vertex_normals();
		mesh->update_vertex_normals();
	}

	int n_constraints = std::max(20, (int)mesh->n_faces() / 2);
	StatisticsMostValues<double, true> max_cur(n_constraints);

	std::vector<Eigen::Matrix2d> e_dm(mesh->n_edges());

	OpenMesh::HPropHandleT<COMPLEX> h_local;
	mesh->add_property(h_local);
	for (auto &th : mesh->halfedges())
	{
		if (mesh->is_boundary(th))
			continue;
		auto fid = mesh->face_handle(th).idx() * 2;
		/*auto vec = mesh->calc_edge_vector(th);
		mesh->property(h_local, th) = COMPLEX(faceBase[fid].dot(vec), faceBase[fid + 1].dot(vec));*/
		auto vec = position.col(th.to().idx()) - position.col(th.from().idx());
		mesh->property(h_local, th) = COMPLEX(faceBase.col(fid).dot(vec), faceBase.col(fid + 1).dot(vec));
	}

	for (auto e_h : mesh->edges())
	{
		if (mesh->is_boundary(e_h)) continue;

		auto h0 = mesh->halfedge_handle(e_h, 0);
		auto h1 = mesh->halfedge_handle(e_h, 1);

		auto vec0 = mesh->calc_edge_vector(h0);
		auto vec1 = mesh->calc_edge_vector(h1);

		/*auto el0 = +h_local[h0];
		auto el1 = -h_local[h1];*/
		auto el0 = mesh->property(h_local, h0);
		auto el1 = -mesh->property(h_local, h1);

		el0 /= std::sqrt(std::norm(el0));
		el1 /= std::sqrt(std::norm(el1));

		Eigen::Matrix2d Re0, Re1;
		Re0 << el0.real(), -el0.imag(), el0.imag(), el0.real();
		Re1 << el1.real(), -el1.imag(), el1.imag(), el1.real();

		e_dm[e_h.idx()] = Re1 * Re0.transpose();
	}

	std::vector<Eigen::Matrix2d> f_II(mesh->n_faces());

	for (auto f_h : mesh->faces())
	{
		int f_id = f_h.idx();
		Eigen::Matrix3d mat_LMN;
		Eigen::Vector3d vec_LMN;

		int i_line = 0;
		for (auto fh_h : mesh->fh_range(f_h))
		{
			//auto vec2 = h_local[fh_h];
			auto vec2 = mesh->property(h_local, fh_h);
			auto dn = mesh->normal(mesh->to_vertex_handle(fh_h)) - mesh->normal(mesh->from_vertex_handle(fh_h));
			dn -= mesh->normal(f_h) * OpenMesh::dot(dn, mesh->normal(f_h));

			if (dn.norm() < 1.0e-9) dn = OpenMesh::Vec3d(0.0);

			mat_LMN(i_line, 0) = 1.0 * vec2.real() * vec2.real();
			mat_LMN(i_line, 1) = 2.0 * vec2.real() * vec2.imag();
			mat_LMN(i_line, 2) = 1.0 * vec2.imag() * vec2.imag();

			vec_LMN(i_line) = -OpenMesh::dot(dn, mesh->calc_edge_vector(fh_h));

			i_line++;
		}

		vec_LMN = mat_LMN.fullPivHouseholderQr().solve(vec_LMN);

		f_II[f_h.idx()] << vec_LMN(0), vec_LMN(1), vec_LMN(1), vec_LMN(2);
	}
	mesh->remove_property(h_local);

	for (int i = 0; i < 100; i++)
	{
		std::vector<Eigen::Matrix2d> f_II_new(mesh->n_faces());

		for (auto f0 : mesh->faces())
		{
			f_II_new[f0.idx()] = f_II[f0.idx()];
			double sum_weight = 1.0;

			for (auto fh_h : mesh->fh_range(f0))
			{
				auto f1 = mesh->opposite_face_handle(fh_h);
				if (!f1.is_valid()) continue;

				Eigen::Matrix2d dm = e_dm[fh_h.idx() / 2];
				Eigen::Matrix2d f1_II = f_II[f1.idx()];

				f1_II = (fh_h.idx() & 1) ? (Eigen::Matrix2d)(dm * f1_II * dm.transpose()) : (Eigen::Matrix2d)(dm.transpose() * f1_II * dm);

				f_II_new[f0.idx()] += f1_II;
				sum_weight += 1.0;
			}

			f_II_new[f0.idx()] /= sum_weight;
		}

		f_II = std::move(f_II_new);
	}

	//OpenMesh::FPropHandleT<OpenMesh::Vec3d> f0;
	OpenMesh::FPropHandleT<Eigen::Vector3d> f0;
	mesh->add_property(f0);

	for (auto f_h : mesh->faces())
	{
		int f_id = f_h.idx();
		Eigen::EigenSolver<Eigen::Matrix2d> cur_solver(f_II[f_h.idx()]);

		double k0 = std::abs(cur_solver.eigenvalues()[0].real());
		double k1 = std::abs(cur_solver.eigenvalues()[1].real());

		//std::cout << k0 << " " << k1 << std::endl;
		if (k0 <= 10.0 * k1 && k1 <= 10.0 * k0) continue;
		//if (k0 <= 1.0 * k1 && k1 <= 1.0 * k0) continue;

		Eigen::Matrix2d eigen_vectors = cur_solver.eigenvectors().real();

		int max_k = (k0 > k1) ? 0 : 1;
		//f0[f_h] = eigen_vectors(0, max_k) * e0[f_h] + eigen_vectors(1, max_k) * e1[f_h];
		//mesh->property(f0, f_h) = eigen_vectors(0, max_k)*faceBase[f_id * 2] + eigen_vectors(1, max_k)*faceBase[f_id * 2 + 1];
		mesh->property(f0, f_h) = eigen_vectors(0, max_k)*faceBase.col(f_id * 2) + eigen_vectors(1, max_k)*faceBase.col(f_id * 2 + 1);
		if(k0 > 1.0e-2 || k1 > 1.0e-2)
			max_cur.update(f_h.idx(), std::abs(k0 * k1));
	}
	std::cout << "#Constraints " << max_cur.size() << std::endl;

	auto id = max_cur.get_ids();
	constraintId.resize(id.size());
	constraintVector.resize(3, id.size());
	for (int i = 0; i < id.size(); ++i)
	{
		constraintId[i] = id[i];
		constraintVector.col(i) = mesh->property(f0, mesh->face_handle(id[i]));
	}

	mesh->remove_property(f0);
}

//void crossField::setCurvatureConstraint()
//{
//	constraintId.clear();
//	constraintVector.clear();
//
//	std::vector<double> K1;
//	std::vector<double> K2;
//	std::vector<OpenMesh::Vec3d> dir1;
//	std::vector<OpenMesh::Vec3d> dir2;
//	compute_principal_curvature(mesh, K1, K2, dir1, dir2);
//
//	//由顶点上的曲率张量估计面上曲率张量 https://zhuanlan.zhihu.com/p/386668016
//	Eigen::Matrix3d *C = new Eigen::Matrix3d[mesh->n_vertices()];
//	Eigen::Matrix3d P;
//	Eigen::Matrix3d D; D.setZero();
//	Eigen::Vector3d singular;
//	for (auto &tv : mesh->vertices())
//	{
//		auto id = tv.idx();
//		auto n = dir1[id].cross(dir2[id]);
//		P << dir1[id][0], dir2[id][0], n[0],
//			dir1[id][1], dir2[id][1], n[1],
//			dir1[id][2], dir2[id][2], n[2];
//		D(0, 0) = K1[id]; D(1, 1) = K2[id];
//		C[id] = P * D * P.transpose();
//	}
//
//	for (auto &tf : mesh->faces())
//	{
//		P.setZero();
//		for (auto &tfv : mesh->fv_range(tf))
//		{
//			P += C[tfv.idx()];
//		}
//		Eigen::JacobiSVD<Eigen::Matrix3d> svd(P * 0.3333333333333333333, Eigen::ComputeFullU | Eigen::ComputeFullV);
//		singular = svd.singularValues();
//		double k1 = std::fabs(singular(0, 0)); double k2 = std::fabs(singular(1, 1));
//
//		if (k1 < 1.0e-4 && k2 < 1.0e-4 || k1 < 30 * k2 && k2 < 30 * k1)
//			continue;
//		constraintId.push_back(tf.idx());
//		singular = svd.matrixU().col(0);
//		constraintVector.emplace_back(singular(0), singular(1), singular(2));
//	}
//	delete C;
//}

void crossField::setNormal()
{
	normal.resize(3, mesh->n_faces());
	for (auto& tf : mesh->faces())
	{
		std::vector<int> id; id.reserve(3);
		for (auto& tfv : mesh->fv_range(tf)) id.push_back(tfv.idx());
		normal.col(tf.idx()) = (position.col(id[1]) - position.col(id[0])).cross(position.col(id[2]) - position.col(id[1])).normalized();
	}
}

void crossField::setMatching()
{
	matching.resize(mesh->n_halfedges());
	//matching.setZero();
	double invHalfPI = 2.0 / PI;
	for (auto& te : mesh->edges())
	{
		if (te.is_boundary())
		{
			matching[te.h0().idx()] = 0;
			matching[te.h1().idx()] = 0;
			continue;
		}
		auto fid = mesh->face_handle(te.h0()).idx();
		auto gid = mesh->face_handle(te.h1()).idx();
		Eigen::Vector3d ev = position.col(te.h0().to().idx()) - position.col(te.h0().from().idx());
		auto fec = COMPLEX(faceBase.col(fid * 2).dot(ev), faceBase.col(fid * 2 + 1).dot(ev));
		auto gec = COMPLEX(faceBase.col(gid * 2).dot(ev), faceBase.col(gid * 2 + 1).dot(ev));
		auto fc = COMPLEX(crossfield.col(fid * 4).dot(faceBase.col(fid * 2)), crossfield.col(fid * 4).dot(faceBase.col(fid * 2 + 1)));
		auto gc = COMPLEX(crossfield.col(gid * 4).dot(faceBase.col(gid * 2)), crossfield.col(gid * 4).dot(faceBase.col(gid * 2 + 1)));
		int m = std::floor((std::arg(fc * gec / (gc * fec)) + PI * 0.25) * invHalfPI);
		matching[te.h0().idx()] = (m + 4) % 4;
		matching[te.h1().idx()] = (4 - m) % 4;
	}
}

 void crossField::setSingularity()
{
	 int count = 0;
	 for (auto& tv : mesh->vertices())
	 {
		 int sum = 0;
		 for (auto& tvoh : mesh->voh_range(tv))
		 {
			 sum += matching[tvoh.idx()];
		 }
		 if (sum % 4)
		 {
			 //dprint(tv.idx());
			 //singularity[count] = tv.idx();
			 singularity.push_back(tv.idx());
			 ++count;
		 }
	 }
}

 crossField::crossField(std::string& filename)
 {
	 std::ifstream file_reader;
	 file_reader.open(filename, std::ios::in);
	 char line[1024] = { 0 };
	 file_reader.getline(line, sizeof(line));
	 std::stringstream n(line);
	 int mark[5];
	 n >> mark[0] >> mark[1] >> mark[2] >> mark[3] >> mark[4];
	 crossfield.resize(3, mark[0]); mark[0] += 1;
	 normal.resize(3, mark[1]);     mark[1] += mark[0];
	 matching.resize(mark[2]);      mark[2] += mark[1];
	 position.resize(3, mark[3]);   mark[3] += mark[2];
	 singularity.resize(mark[4]);   mark[4] += mark[3];
	 int row = 1;
	 while (file_reader.getline(line, sizeof(line)))
	 {
		 std::stringstream num(line);
		 if (row < mark[0])
		 {
			 num >> crossfield(0, row - 1) >> crossfield(1, row - 1) >> crossfield(2, row - 1);
		 }
		 else if (row < mark[1])
		 {
			 num >> normal(0, row - mark[0]) >> normal(1, row - mark[0]) >> normal(2, row - mark[0]);
		 }
		 else if (row < mark[2])
		 {
			 num >> matching[row - mark[1]];
		 }
		 else if (row < mark[3])
		 {
			 num >> position(0, row - mark[2]) >> position(1, row - mark[2]) >> position(2, row - mark[2]);
		 }
		 else
		 {
			 num >> singularity[row - mark[3]];
		 }
		 ++row;
	 }
	 file_reader.close();
 }

 void crossField::write_field(std::string& field_file)
 {
	 std::ofstream file_writer;
	 file_writer.open(field_file);
	 if (file_writer.fail()) {
		 std::cout << "fail to open\n";
	 }
	 //crossfield, normal, matching, position, singularity
	 file_writer << mesh->n_faces() * 4 << " " << mesh->n_faces() << " " << mesh->n_halfedges() << " " << mesh->n_vertices() << " " << singularity.size() << "\n";
	 for (auto i = 0; i < crossfield.cols(); ++i)
		 file_writer << crossfield(0, i) << " " << crossfield(1, i) << " " << crossfield(2, i) << "\n";
	 for (auto i = 0; i < normal.cols(); ++i)
		 file_writer << normal(0, i) << " " << normal(1, i) << " " << normal(2, i) << "\n";
	 for (auto match : matching)
		 file_writer << match << "\n";
	 for (auto i = 0; i < position.cols(); ++i)
		 file_writer << position(0, i) << " " << position(1, i) << " " << position(2, i) << "\n";
	 for (auto sing : singularity)
		 file_writer << sing << "\n";

	 file_writer.close();
 }