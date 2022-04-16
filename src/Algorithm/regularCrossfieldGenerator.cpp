#include "regularCrossfieldGenerator.h"
void regularCrossfieldGenerator::calculateMeshFaceBase()
{
	mesh->update_face_normals();
	faceBase.reserve(mesh->n_faces() * 2);
	for (auto& tf : mesh->faces())
	{
		std::vector<OpenMesh::Vec3d> v;
		v.reserve(3);
		for (auto& tfv : mesh->fv_range(tf))
		{
			v.push_back(mesh->point(tfv));
		}
		faceBase.push_back((v[1] - v[0]).normalized());
		faceBase.push_back(cross(faceBase.back(), cross(v[2] - v[0], v[1] - v[0])).normalized());
	}
}

void regularCrossfieldGenerator::run(std::vector<OpenMesh::Vec3d>& crossfield)
{
	using namespace std;
	using namespace Eigen;
	typedef complex<double> COMPLEX;

	auto nf = mesh->n_faces();
	auto ne = mesh->n_edges();
	int count = 0;

	/*Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
	Eigen::SparseMatrix<COMPLEX> A;
	Eigen::VectorXcd b(fnum);
	vector<Eigen::Triplet<COMPLEX>> tris;*/
	SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseMatrix<double> B_CCB, ATA, Hessian;
	B_CCB.resize(ne << 1, nf << 1);
	ATA.resize(nf << 1, nf << 1);
	Hessian.resize(nf << 1, nf << 1);
	VectorXd b(nf << 1);
	vector<Triplet<double>> tri;

	count = 0;
	for (auto& tf : mesh->faces())
	{
		int id_f = tf.idx();

		COMPLEX sum = 0;
		for (auto& tfh : mesh->fh_range(tf))
		{
			if (!mesh->is_boundary(tfh.edge()))
			{
				auto id_g = mesh->face_handle(tfh.opp()).idx();
				if (id_f < id_g)
				{
					auto p1 = mesh->point(tfh.to());
					auto p2 = mesh->point(tfh.from());
					OpenMesh::Vec3d e = (p2 - p1).normalized();
					COMPLEX ef = COMPLEX(dot(e, faceBase[id_f << 1]), -dot(e, faceBase[(id_f << 1) + 1]));
					ef *= ef; ef *= ef;
					COMPLEX eg = COMPLEX(dot(e, faceBase[id_g << 1]), -dot(e, faceBase[(id_g << 1) + 1]));
					eg *= eg; eg *= eg;

					/*OpenMesh::Vec3d e = (p2 - p1).normalized();

					COMPLEX e_f = COMPLEX(dot(e, faceBase[id_f * 2]), dot(e, faceBase[id_f * 2 + 1]));
					COMPLEX e_g = COMPLEX(dot(e, faceBase[id_g * 2]), dot(e, faceBase[id_g * 2 + 1]));

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
					}*/


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
		dprint(std::sqrt(x(i).real() * x(i).real() + x(i).imag() * x(i).imag()));
	}
#endif
	crossfield.resize(4 * fnum);
	for (int i = 0; i < fnum; i++)
	{
		double length = std::sqrt(f_dir[i].real() * f_dir[i].real() + f_dir[i].imag() * f_dir[i].imag());
		double arg = std::arg(f_dir[i]) / 4;
		for (int j = 0; j < 4; j++)
		{
			crossfield[i * 4 + j] = faceBase[i * 2] * length * cos(arg + j * PI * 0.5) + faceBase[i * 2 + 1] * length * sin(arg + j * PI * 0.5);
		}
	}
}