#pragma once
#include "../src/MeshViewer//MeshDefinition.h"
namespace LoopGen
{
	class AngleDefectMinimizer
	{
	public:
		AngleDefectMinimizer(Mesh* m) :mesh(m) {};
		~AngleDefectMinimizer() {};
	public:
		void run();
	private:
		Mesh* mesh;
	};
}
