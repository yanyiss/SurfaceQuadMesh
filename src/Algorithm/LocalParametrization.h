#include "crossField.h"
#include "LoopDef.h"
#define PRINT_DEBUG_INFO 0
namespace LoopGen
{
	class LocalParametrization
	{
	public:
		LocalParametrization() {};
		LocalParametrization(VertexLayer* vl_, cylinder &cy_, spread_info &sp_);
		~LocalParametrization() {};
	public:
		inline double GetRegularU(int vlid) { return cy->GetRegularU(vlid); }
		inline double GetU(int vlid) { return cy->GetU(vlid); }
		inline double GetV(int vlid) { return cy->GetV(vlid); }
		inline Eigen::VectorXd& GetU() { return cy->uv[0]; }
		inline Eigen::VectorXd& GetV() { return cy->uv[1]; }

		cylinder* cy;
		spread_info *sp;
		BoolVector cutv_flag;
		BoolVector cutf_flag;

		void run(const Eigen::Matrix3Xd &normal);
		void modify_cut();
	};
}