#pragma once

#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <numeric>
#include "..\src\Dependency\CPS\CPS_AABBTree.h"
#include "..\src\Algorithm\SurfaceMesher\Generator\basic_def.h"

namespace CADMesher
{
	class TriangleMeshRemeshing
	{
	public:
		explicit TriangleMeshRemeshing(TriMesh *mesh_, double target_length = -1)
			:mesh(mesh_), expected_length(target_length)
		{
			if (expected_length <= 0)
			{
				expected_length = meshAverageLength(*mesh);
			}
			high = 1.33*expected_length;
			low = 0.8*expected_length;
			aabbtree = new ClosestPointSearch::AABBTree(*mesh);
		};
		TriangleMeshRemeshing(const TriangleMeshRemeshing &tmr) = delete;
		~TriangleMeshRemeshing() { 
			if (aabbtree) { delete aabbtree; aabbtree = nullptr; } 
		}

	public:
		void run();

	private:
		//main step
		void split();
		void collapse();
		void equalize_valence();
		void tangential_relaxation();
		//auxiliary step
		void adjustTargetLength();
		void processAngle();
		//geometry support
		O3d GravityPos(const OV &v);


	private:
		double high;
		double low;
		double expected_length;

		double lowerAngleBound = 0.05;
		TriMesh *mesh = nullptr;
		ClosestPointSearch::AABBTree *aabbtree = nullptr;

#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
	public:
		explicit TriangleMeshRemeshing(PolyMesh *mesh_, double target_length = -1);//, polymeshInput(true)
	private:
		bool polymeshInput = false;
		int boundaryNum;
		PolyMesh *polymesh;
	private:
		void assembleMesh();
#endif
	};
}

