#pragma once
#include "crossField.h"
namespace QuadLayout
{
	using namespace Eigen;

	class dualLoop
	{
	public:
		dualLoop(Mesh *m) :mesh(m) {
			cf = new crossField(mesh);
		};
		dualLoop(const dualLoop& dl) = delete;
		~dualLoop()
		{
			if (cf) { delete cf; cf = nullptr; };
		}

	private:
		void initGraph();

	private:
		MatrixXd graph[2];
		crossField* cf = nullptr;
		Mesh *mesh;
	};
}
