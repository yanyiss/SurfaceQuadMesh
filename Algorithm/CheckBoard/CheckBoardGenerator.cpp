#include "CheckBoardGenerator.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
using namespace Eigen;
using namespace OpenMesh;
#define updateConstraint

void CheckBoardGenerator::init()
{
	//assert mesh is closed 
	assert(isClosedMesh(*originalMesh));
	dprint("Closed Mesh!");

	//construct aabbtree
	aabbtree = new ClosestPointSearch::AABBTree(*originalMesh);
	dprint("construct aabbtree finished!");

	//init plane function of each triangle 
	auto nf = originalMesh->n_faces();
	cof.resize(3, nf);
	ct.resize(nf);
	for (OF tf : originalMesh->faces())
	{
		Vec3d normal = originalMesh->calc_face_normal(tf);
		Vec3d pos = originalMesh->point(originalMesh->fv_begin(tf));
		//cof.col(tf.idx()) << normal[0] / 2 , normal[1] / 2, normal[2] / 2;//here we divided by 2 because we replace (v_i+f_j)/2 with (v_i+f_j)
		cof.col(tf.idx()) << normal[0], normal[1], normal[2];
		ct(tf.idx()) = -normal.dot(pos);
	}
	dprint("init plane function finised!");

	//init primal mesh and dual mesh
	Mesh polymesh;
	OpenMesh::IO::read_mesh(polymesh, "../model/human_quad.obj");
	initMeshStatusAndNormal(polymesh);

	num[0] = polymesh.n_vertices();
	num[1] = polymesh.n_faces();
	num[2] = polymesh.n_edges();

	//set primal mesh info
	m3xd& primalV = v[0];                                 primalV.resize(3, num[0]);
	std::vector<std::vector<triple>>& primalAdj = adj[0]; primalAdj.resize(num[0]);
	for (OV tv : polymesh.vertices())
	{
		auto vid = tv.idx();
		//set adj
		std::vector<triple>& pAv = primalAdj[vid];
		pAv.reserve(polymesh.valence(tv));
		for (auto tvoh : polymesh.voh_range(tv))
		{
			pAv.emplace_back(polymesh.to_vertex_handle(tvoh).idx(), polymesh.opposite_face_handle(tvoh).idx(), polymesh.face_handle(tvoh).idx());
		}
		//set pos
		Vec3d pos = polymesh.point(tv);
		primalV.col(vid) << pos[0], pos[1], pos[2];
	}
	//set dual mesh info
	m3xd& dualV = v[1];                                   dualV.resize(3, num[1]);
	std::vector<std::vector<triple>>& dualAdj = adj[1];   dualAdj.resize(num[1]);
	for (OF tf : polymesh.faces())
	{
		auto fid = tf.idx();
		v3d pos(0, 0, 0);
		//set adj
		std::vector<triple>& dAv = dualAdj[fid];
		dAv.reserve(polymesh.valence(tf));
		for (auto tfh : polymesh.fh_range(tf))
		{
			dAv.emplace_back(polymesh.opposite_face_handle(tfh).idx(), polymesh.from_vertex_handle(tfh).idx(), polymesh.to_vertex_handle(tfh).idx());
			pos += primalV.col(dAv.back().y);
		}
		//set pos
		dualV.col(fid) = pos / polymesh.valence(tf);
	}
	dprint("init mesh info finised!");
}

#include <Eigen\Sparse>
#include <Eigen\SparseCholesky>
void CheckBoardGenerator::run()
{
	std::vector<Triplet<double>> triplet[2];
	SparseMatrix<double> M[2];
	vxd righthand[2];
	for (int i = 0; i < 2; ++i)
	{
		triplet[i].reserve(9 * (2 * num[2] + num[i]));
		M[i].resize(3 * num[i], 3 * num[i]);
		righthand[i].resize(3 * num[i]);
	}

	double lambda = 1.0;
	printCurrentInfo();
	int itertimes = 20;
	for (int i = 1; i < itertimes; ++i)
	{
		dprint("iteration times:", i);
		//model info
		m3xd& activeV = v[i % 2];
		const m3xd& fixedV = v[(i + 1) % 2];
		const std::vector<std::vector<triple>>& activeAdj = adj[i % 2];

		//matrix info
		std::vector<Triplet<double>>& activeTriplet = triplet[i % 2];
		activeTriplet.clear();
		SparseMatrix<double> &activeM = M[i % 2];
		activeM.setZero();
		vxd &activeRighthand = righthand[i % 2];

		//construct triplet with all non-zero elements in sparse matrix
		int n = num[i % 2];
		for (ui j = 0; j < n; ++j)
		{
			m3d matSum(3, 3);
			matSum.setZero();
			auto row = 3 * j;
			v3d vpos = activeV.col(j);
			v3d rh(0, 0, 0);
			for (const triple& aj : activeAdj[j])
			{
				auto col = aj.x * 3;
				const auto &f0 = aj.y;
				const auto &f1 = aj.z;

				//orthogonality energy
				m3d emat = -(fixedV.col(f0) - fixedV.col(f1))*(fixedV.col(f0) - fixedV.col(f1)).transpose() * lambda;
				for (int p = 0; p < 3; ++p)
				{
					for (int q = 0; q < 3; ++q)
					{
						activeTriplet.emplace_back(row + p, col + q, emat(p, q));
					}
				}
				matSum -= emat;

				//error energy
				v3d anchor = (vpos + fixedV.col(f0)) / 2;
				auto fid = aabbtree->closest_point_and_face_handle(anchor.data()).second.idx();
				matSum += cof.col(fid)*cof.col(fid).transpose()*0.25;

				//righthand side
				rh -= 0.5*cof.col(fid)*(0.5*cof.col(fid).dot(fixedV.col(f0)) + ct(fid));
			}
#ifdef updateConstraint
			matSum += m3d::Identity(3, 3);
			rh += vpos;
#endif
			//fill diagonal elements
			for (int p = 0; p < 3; ++p)
			{
				for (int q = 0; q < 3; ++q)
				{
					activeTriplet.emplace_back(row + p, row + q, matSum(p, q));
				}
				activeRighthand(row + p) = rh(p);
			}

			//lambda += 1.0e-8;
		}
		//solve sparse linear problem
		activeM.setFromTriplets(activeTriplet.begin(), activeTriplet.end());
		SimplicialCholesky<SparseMatrix<double>> solver(activeM);
		activeRighthand = solver.solve(activeRighthand);

		//update active mesh vertices position
		for (ui j = 0; j < n; ++j)
		{
			activeV.col(j) << activeRighthand(3 * j), activeRighthand(3 * j + 1), activeRighthand(3 * j + 2);
		}
		printCurrentInfo();
	}
}

void CheckBoardGenerator::getCheckBoard(m3xd &V, std::vector<vxu> &F)
{
	V.resize(3, 2 * num[2]);
	F.resize(num[2]);

	std::vector<std::map<ui, ui>> indexTransfer;
	indexTransfer.resize(num[0]);

	ui count = 0;
	for (int i = 0; i < num[0]; ++i)
	{
		for (const triple &aj : adj[0][i])
		{
			indexTransfer[i].emplace(aj.y, count);
			V.col(count++) = (v[0].col(i) + v[1].col(aj.y)) * 0.5;
		}
	}
	count = 0;
	for (int i = 0; i < num[0]; ++i)
	{
		for (const triple &aj : adj[0][i])
		{
			if (aj.x > i)
				continue;
			F[count].resize(4);
			F[count++] << indexTransfer[i][aj.y], indexTransfer[aj.x][aj.y], indexTransfer[aj.x][aj.z], indexTransfer[i][aj.z];
		}
	}
}

void CheckBoardGenerator::getPrimal(m3xd &V, std::vector<vxu> &F)
{
	V = v[0];
	F.resize(num[1]);
	for (int i = 0; i < num[1]; ++i)
	{
		F[i].resize(adj[1][i].size());
		for (int j = 0; j < F[i].size(); ++j)
		{
			F[i](j) = adj[1][i][j].y;
		}
	}
}

void CheckBoardGenerator::getDual(m3xd &V, std::vector<vxu> &F)
{
	V = v[1];
	F.resize(num[0]);
	for (int i = 0; i < num[0]; ++i)
	{
		F[i].resize(adj[0][i].size());
		for (int j = 0; j < F[i].size(); ++j)
		{
			F[i](j) = adj[0][i][j].y;
		}
	}
}

void CheckBoardGenerator::getMesh(Mesh &m, m3xd &V, std::vector<vxu> &F)
{
	m.clear();
	for (int i = 0; i < V.cols(); ++i)
	{
		v3d pos = V.col(i);
		m.add_vertex(Vec3d(pos(0), pos(1), pos(2)));
	}
	std::vector<OV> vh;
	vh.reserve(8);
	for (const vxu &f : F)
	{
		vh.clear();
		for (int i = 0; i < f.size(); ++i)
		{
			vh.push_back(m.vertex_handle(f(i)));
		}
		m.add_face(vh);
	}
	initMeshStatusAndNormal(m);
}

void CheckBoardGenerator::printCurrentInfo()
{
	double orthogonalityEnergy = 0, errorEnergy = 0;
	double orthogonalityMax = 0, errorMax = 0;

	for (int i = 0; i < num[0]; ++i)
	{
		for (const triple &ai : adj[0][i])
		{
			double temp = (v[0].col(i) - v[0].col(ai.x)).dot(v[1].col(ai.y) - v[1].col(ai.z));
			orthogonalityEnergy += temp * temp;
			/*v3d t0 = v[0].col(i) - v[0].col(ai.x);
			v3d t1 = v[1].col(ai.y) - v[1].col(ai.z);*/
			orthogonalityMax = std::max(orthogonalityMax, temp);
		}
	}
	orthogonalityEnergy *= 0.5;

	for (int i = 0; i < num[0]; ++i)
	{
		for (const triple &ai : adj[0][i])
		{
			v3d pos = 0.5*(v[0].col(i) + v[1].col(ai.y));
			double d = aabbtree->closest_distance(Vec3d(pos(0), pos(1), pos(2)));
			errorEnergy += d * d;
			errorMax = std::max(errorMax, d);
		}
	}

	dprint("current energy:", orthogonalityEnergy, errorEnergy, orthogonalityEnergy + errorEnergy);
	dprint("max single term:", orthogonalityMax, errorMax);
}

#include <deque>
#include <queue>
//void CheckBoardGenerator::setDiagonalMeshIndex()
//{
//	mesh->add_property(diagonalMeshIndex);
//
//	size_t nv = mesh->n_vertices();
//
//	std::queue<OV> vertexQueue;
//	vertexQueue.push(mesh->vertex_handle(0));
//	mesh->property(diagonalMeshIndex, mesh->vertex_handle(0)) = true;
//	std::deque<bool> visited(nv, false);
//	visited[0] = true;
//
//	while (!vertexQueue.empty())
//	{
//		OV tv = vertexQueue.front();
//		vertexQueue.pop();
//		bool flag = !mesh->property(diagonalMeshIndex, tv);
//		for (auto tvv : mesh->vv_range(tv))
//		{
//			if (visited[tvv.idx()])
//				continue;
//			mesh->property(diagonalMeshIndex, tvv) = flag;
//			visited[tvv.idx()] = true;
//			vertexQueue.push(tvv);
//		}
//	}
//}