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

	size_t nf = mesh->n_faces();
	size_t ne = mesh->n_edges();

	SimplicialLDLT<SparseMatrix<double>> solver;
	SparseMatrix<double> B_CCB, ATA, Hessian;
	B_CCB.resize(ne * 2, nf * 2);
	ATA.resize(nf * 2, nf * 2);
	Hessian.resize(nf * 2, nf * 2);
	VectorXd x(nf * 2); x.setConstant(1.0);
	VectorXd b(nf * 2);
	vector<Triplet<double>> tri;

	for(int itertimes = 0; itertimes < 10; ++itertimes)
	{
		dprint("itertimes:", itertimes);
		int count = 0;
		tri.reserve(8 * ne);
		for (auto& th : mesh->halfedges())
		{
			if (mesh->is_boundary(th.edge()))
				continue;
			size_t id_f = th.face().idx();
			size_t id_g = th.opp().face().idx();
			if (id_f < id_g)
			{
				auto& p1 = mesh->point(th.to());
				auto& p2 = mesh->point(th.from());
				OpenMesh::Vec3d e = (p2 - p1).normalized();
				COMPLEX ef = COMPLEX(dot(e, faceBase[id_f * 2]), -dot(e, faceBase[id_f * 2 + 1]));
				ef *= ef; ef *= ef;
				COMPLEX eg = COMPLEX(dot(e, faceBase[id_g * 2]), -dot(e, faceBase[id_g * 2 + 1]));
				eg *= eg; eg *= -eg;

				tri.emplace_back(count, id_f, ef.real());
				tri.emplace_back(count + ne, id_f + nf, ef.real());
				tri.emplace_back(count, id_f + nf, -ef.imag());
				tri.emplace_back(count + ne, id_f, ef.imag());

				tri.emplace_back(count, id_g, eg.real());
				tri.emplace_back(count + ne, id_g + nf, eg.real());
				tri.emplace_back(count, id_g + nf, -eg.imag());
				tri.emplace_back(count + ne, id_g, eg.imag());

				++count;
			}
		}

		B_CCB.setZero();
 		B_CCB.setFromTriplets(tri.begin(), tri.end());
		ATA = B_CCB.transpose().eval() * B_CCB;
		tri.clear();
		tri.reserve(4 * nf);
		for (count = 0; count < nf; ++count)
		{
			double temp = 1.0 / (x(count) * x(count) + x(count + nf) * x(count + nf));
			double t1 = temp * temp; double t2 = temp * t1;
			tri.emplace_back(count, count, 2 * (1 - t1) + 8 * x(count) * x(count) * t2);
			tri.emplace_back(count, count + nf, 8 * x(count) * x(count + nf) * t2);
			tri.emplace_back(count + nf, count, 8 * x(count) * x(count + nf) * t2);
			tri.emplace_back(count + nf, count + nf, 2 * (1 - t1) + 8 * x(count + nf) * x(count + nf) * t2);
		}
		Hessian.setZero();
		Hessian.setFromTriplets(tri.begin(), tri.end());
		Hessian += ATA;
		
		for (count = 0; count < nf; ++count)
		{
			double temp = x(count) * x(count) + x(count + nf) * x(count + nf);
			temp = 2 * (1 - 1.0 / temp * temp);
			b(count) = x(count) * temp;
			b(count + nf) = x(count + nf) * temp;
		}
		b += ATA * x;
		
		if (!itertimes)
		{
			solver.analyzePattern(Hessian);
		}
		//solver.compute(Hessian);
		solver.factorize(Hessian);
		x -= solver.solve(b);
		for (count = 0; count < nf; ++count)
		{
			double temp = x(count) * x(count) + x(count + nf) * x(count + nf);
			if (temp < 1.0)
			{
				temp = std::sqrt(temp);
				x(count) /= temp;
				x(count + nf) /= temp;
			}
		}
	}
	crossfield.resize(4 * nf);
	for (int i = 0; i < nf; ++i)
	{
		//dprint(i);
		double arg = std::arg(COMPLEX(x(i), x(i + nf))) * 0.25;
		for (int j = 0; j < 4; j++)
		{
			crossfield[i * 4 + j] = faceBase[i * 2] * cos(arg + j * PI * 0.5) + faceBase[i * 2 + 1] * sin(arg + j * PI * 0.5);
		}
	}
}