#include "crossField.h"
namespace LoopGen
{
	void crossField::setPosition()
	{
		position.resize(3, mesh->n_vertices());
		for (auto& tv : mesh->vertices())
		{
			auto& p = mesh->point(tv);
			position.col(tv.idx()) << p[0], p[1], p[2];
		}
	}

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

	void crossField::setFaceBase()
	{
		mesh->update_face_normals();
		faceBase.resize(3, mesh->n_faces() * 2);
		for (auto& tf : mesh->faces())
		{
			Eigen::Matrix3d v;
			int i = 0;
			for (auto& tfv : mesh->fv_range(tf))
			{
				v.col(i) = position.col(tfv.idx());
				++i;
			}
			i = tf.idx() * 2;
			faceBase.col(i) = (v.col(1) - v.col(0)).normalized();
			faceBase.col(i + 1) = faceBase.col(i).cross((v.col(2) - v.col(0)).cross(v.col(1) - v.col(0))).normalized();
		}
	}

	void crossField::initMeshInfo()
	{
		setPosition();
		setNormal();
		setFaceBase();
	}

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
		for (auto& th : mesh->halfedges())
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
			mesh->property(f0, f_h) = eigen_vectors(0, max_k) * faceBase.col(f_id * 2) + eigen_vectors(1, max_k) * faceBase.col(f_id * 2 + 1);
			if (k0 > 1.0e-2 || k1 > 1.0e-2)
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

	void crossField::setOuterConstraint(std::deque<bool>& cons_flag, Eigen::Matrix3Xd& cons_direction)
	{
		int count = 0;
		for (int i = 0; i < mesh->n_faces(); ++i)
			if (cons_flag[i])
				++count;
		constraintId.clear();
		constraintId.reserve(count);
		constraintVector.resize(3, count);
		int cols = 0;
		for (int i = 0; i < mesh->n_faces(); ++i)
		{
			if (cons_flag[i])
			{
				constraintId.push_back(i);
				constraintVector.col(cols) = cons_direction.col(i);
				++cols;
			}
		}
	}

	void crossField::setField()
	{
		using namespace std;

		auto fnum = mesh->n_faces();
		vector<int> status(fnum, 0);

		std::vector<COMPLEX> f_dir(fnum);

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
		for (auto& tf : mesh->faces())
		{
			int id_f = tf.idx();
			for (auto& tfh : mesh->fh_range(tf))
			{
				if (!mesh->is_boundary(tfh.edge()))
				{
					auto& p1 = position.col(tfh.to().idx());
					auto& p2 = position.col(tfh.from().idx());
					auto id_g = mesh->face_handle(tfh.opp()).idx();
					if (id_f < id_g)
					{
						if (status[id_f] == 1 && status[id_g] == 1)
							continue;
						auto e = (p2 - p1).normalized();
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
		for (int i = 0; i < fnum; i++)
		{
			if (status[i] == 0)
			{
				f_dir[i] = x(id2sln[i]);
			}
		}
		crossfield.resize(3, 4 * fnum);
		for (int i = 0; i < fnum; i++)
		{
			double arg = std::arg(f_dir[i]) * 0.25;
			for (int j = 0; j < 4; ++j)
			{
				crossfield.col(i * 4 + j) = faceBase.col(i * 2) * cos(arg + j * PI * 0.5) + faceBase.col(i * 2 + 1) * sin(arg + j * PI * 0.5);
			}
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
		singularity.clear();
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

	//void crossField::setWeight(double alpha)
	//{
	//	double doublePI = 2.0 * PI;
	//	double halfPI = 0.5 * PI;
	//	double _PI = -1.0 * PI;
	//	double _halfPI = -0.5 * PI;
	//	weight.resize(4, mesh->n_halfedges());
	//	for (auto eitr = mesh->edges_begin(); eitr != mesh->edges_end(); ++eitr)
	//	{
	//		auto h0 = mesh->halfedge_handle(eitr.handle(), 0);
	//		auto h1 = mesh->halfedge_handle(eitr.handle(), 1);
	//		auto fid = mesh->face_handle(h0).idx();
	//		auto gid = mesh->face_handle(h1).idx();
	//		auto& fv = crossfield.col(fid * 4);
	//		auto& gv = crossfield.col(gid * 4 + matching[h0.idx()]);
	//		auto ev = position.col(mesh->to_vertex_handle(h0).idx()) - position.col(mesh->from_vertex_handle(h0).idx());
	//		double arc0 = atan2(ev.cross(fv).dot(normal.col(fid)), ev.dot(fv)); //arc0 += arc0 > 0 ? 0 : doublePI;
	//		double arc1 = atan2(ev.cross(gv).dot(normal.col(gid)), ev.dot(gv)); //arc1 += arc1 > 0 ? 0 : doublePI;
	//		double arc = atan2(sin(arc0) + sin(arc1), cos(arc0) + cos(arc1));
	//		auto& w0 = weight.col(h0.idx());
	//		auto& w1 = weight.col(h1.idx());
	//		double s = fabs(sin(arc));
	//		double c = fabs(cos(arc));
	//		if (s < c)
	//		{
	//			s = sqrt(alpha * s * s + 1);
	//			c = YYSS_INFINITE;
	//		}
	//		else
	//		{
	//			s = YYSS_INFINITE;
	//			c = sqrt(alpha * c * c + 1);
	//		}
	//		if (arc >= 0 && arc < halfPI)
	//			w0 << s, YYSS_INFINITE, YYSS_INFINITE, c;
	//		else if (arc >= halfPI && arc < PI)
	//			w0 << YYSS_INFINITE, YYSS_INFINITE, s, c;
	//		else if (arc >= _PI && arc < _halfPI)
	//			w0 << YYSS_INFINITE, c, s, YYSS_INFINITE;
	//		else
	//			w0 << s, c, YYSS_INFINITE, YYSS_INFINITE;
	//		switch (matching[h0.idx()])
	//		{
	//		case 0:
	//			w1 << w0(2), w0(3), w0(0), w0(1);
	//			break;
	//		case 1:
	//			w1 << w0(1), w0(2), w0(3), w0(0);
	//			break;
	//		case 2:
	//			w1 << w0(0), w0(1), w0(2), w0(3);
	//			break;
	//		case 3:
	//			w1 << w0(3), w0(0), w0(1), w0(2);
	//			break;
	//		}
	//	}
	//	dprint("Initialize Graph Weight Done!");
	//}

	void crossField::initFieldInfo()
	{
		setMatching();
		setSingularity();
		//setWeight();
	}

	void crossField::read_field()
	{
		std::ifstream file_reader;
		file_reader.open(file_name, std::ios::in);
		char line[1024] = { 0 };
		file_reader.getline(line, sizeof(line));
		std::stringstream n(line);
		int mark[6];
		n >> mark[0] >> mark[1] >> mark[2] >> mark[3] >> mark[4] >> mark[5];// >> mark[6];
		crossfield.resize(3, mark[0]); mark[0] += 1;
		normal.resize(3, mark[1]);     mark[1] += mark[0];
		matching.resize(mark[2]);      mark[2] += mark[1];
		position.resize(3, mark[3]);   mark[3] += mark[2];
		singularity.resize(mark[4]);   mark[4] += mark[3];
		faceBase.resize(3, mark[5]);   mark[5] += mark[4];
		//weight.resize(4, mark[6]);     mark[6] += mark[5];
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
			else if (row < mark[4])
			{
				num >> singularity[row - mark[3]];
			}
			else if (row < mark[5])
			{
				num >> faceBase(0, row - mark[4]) >> faceBase(1, row - mark[4]) >> faceBase(2, row - mark[4]);
			}
			/*else
			{
				num >> weight(0, row - mark[5]) >> weight(1, row - mark[5]) >> weight(2, row - mark[5]) >> weight(3, row - mark[5]);
			}*/
			++row;
		}
		file_reader.close();
	}

	void crossField::write_field()
	{
		std::ofstream file_writer;
		file_writer.open(file_name);
		if (file_writer.fail()) {
			std::cout << "fail to open\n";
		}
		//crossfield, normal, matching, position, singularity, faceBase, weight
		file_writer << mesh->n_faces() * 4 << " " << mesh->n_faces() << " " << mesh->n_halfedges() << " " << mesh->n_vertices()
			<< " " << singularity.size() << " " << mesh->n_faces() * 2 /*<< " " << mesh->n_halfedges() */<<"\n";
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
		for (auto i = 0; i < faceBase.cols(); ++i)
			file_writer << faceBase(0, i) << " " << faceBase(1, i) << " " << faceBase(2, i) << "\n";
		//for (auto i = 0; i < weight.cols(); ++i)
		//	file_writer << weight(0, i) << " " << weight(1, i) << " " << weight(2, i) << " " << weight(3, i) << "\n";

		file_writer.close();
	}
}