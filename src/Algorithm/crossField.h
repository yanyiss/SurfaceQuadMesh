#include "..\MeshViewer\MeshDefinition.h"
#include <complex>
#include <Eigen\Sparse>

#define PI 3.14159265358979323
#include "..\Toolbox\dprint.h"
#include "..\Toolbox\filesOperator.h"
class crossField
{
public:
	crossField(TriMesh *m) : mesh(m) { calculateMeshFaceBase(); };
	~crossField() {};

private:
	void calculateMeshFaceBase()
	{
		mesh->update_face_normals();
		faceBase.reserve(mesh->n_faces() * 2);
		for (auto &tf : mesh->faces())
		{
			std::vector<OpenMesh::Vec3d> v;
			v.reserve(3);
			for (auto &tfv : mesh->fv_range(tf))
			{
				v.push_back(mesh->point(tfv));
			}
			faceBase.push_back((v[1] - v[0]).normalized());
			faceBase.push_back(cross(faceBase.back(), cross(v[2] - v[0], v[1] - v[0])).normalized());
		}
	}
public:
	void crossfieldCreator(std::vector<OpenMesh::Vec3d>& crossfield, std::vector<size_t> &constraintId, std::vector<OpenMesh::Vec3d> &constraintVector)
	{
		using namespace std;
		typedef complex<double> COMPLEX;

		auto fnum = mesh->n_faces();
		vector<int> status(fnum, 0);
		vector<COMPLEX> f_dir(fnum);

		for (int i = 0; i < constraintId.size(); i++)
		{
			int fid = constraintId[i];
			status[fid] = 1;
			OpenMesh::Vec3d cf = constraintVector[i].normalized();
			f_dir[fid] = std::pow(COMPLEX(dot(cf, faceBase[fid * 2]), dot(cf, faceBase[fid * 2 + 1])), 4);
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
		for(auto &tf : mesh->faces())
		{
			int id_f = tf.idx();

			//COMPLEX sum = 0;
			for(auto &tfh : mesh->fh_range(tf))
			{
				if(!mesh->is_boundary(tfh.edge()))
				{
					auto p1 = mesh->point(tfh.to());
					auto p2 = mesh->point(tfh.from());
					auto id_g = mesh->face_handle(tfh.opp()).idx();
					if (id_f < id_g)
					{
						if (status[id_f] == 1 && status[id_g] == 1)continue;
						OpenMesh::Vec3d e = (p2 - p1).normalized();

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
			dprint(std::sqrt(x(i).real()*x(i).real()+x(i).imag()*x(i).imag()));
		}
#endif
		crossfield.resize(4 * fnum);
		for (int i = 0; i < fnum; i++)
		{
			if (status[i] == 0)
			{
				f_dir[i] = x(id2sln[i]);
			}
			/*double length = 1;
			double arg = std::arg(f_dir[i]) / 4;
			for (int j = 0; j < 4; j++)
			{
				crossfield[i * 4 + j] = faceBase[i * 2] * length * cos(arg + j * PI * 0.5) + faceBase[i * 2 + 1] * length * sin(arg + j * PI * 0.5);
			}*/
			double length = std::sqrt(f_dir[i].real()*f_dir[i].real() + f_dir[i].imag()*f_dir[i].imag());
			double arg = std::arg(f_dir[i]) / 4;
			for (int j = 0; j < 4; j++)
			{
				crossfield[i * 4 + j] = faceBase[i * 2] * cos(arg + j * PI * 0.5) + faceBase[i * 2 + 1] * sin(arg + j * PI * 0.5);
			}
		}
		/*auto cnorm = [&](COMPLEX c)
		{
			return std::sqrt(c.real()*c.real() + c.imag()*c.imag());
		};
		std::vector<double> data;
		data.reserve(mesh->n_edges());
		for (auto &te : mesh->edges())
		{
			if (mesh->is_boundary(te)) continue;
			auto v0 = f_dir[mesh->face_handle(te.h0()).idx()];
			auto v1 = f_dir[mesh->face_handle(te.h1()).idx()];
			data.push_back(cnorm(v0 / cnorm(v0) - v1 / cnorm(v1)));
		}
		WriteVector("data.txt", data);*/
	}

private:
	TriMesh* mesh;
	std::vector<OpenMesh::Vec3d> faceBase;

};
