#include "crossField.h"
#include "LoopDef.h"
#define PRINT_DEBUG_INFO 0
namespace LoopGen
{
	class LocalParametrization
	{
	public:
		LocalParametrization() {};
		LocalParametrization(M4& m4_, VertexLayer* vl_);
		~LocalParametrization() {};
	public:
		inline double GetRegularU(int vlid) { return uv[0](vidmap[vlid]) - std::floor(uv[0](vidmap[vlid])); }
		inline double GetU(int vlid) { return uv[0](vidmap[vlid]); }
		inline double GetV(int vlid) { return uv[1](vidmap[vlid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }


		//private:
		M4* m4;

		std::vector<VertexLayer*> new_vertex;
		std::vector<FaceLayer*> new_face;
		std::vector<VertexLayer*> region_vertex;
		std::vector<FaceLayer*> region_face;
		std::vector<VertexLayer*> cut;

		BoolVector new_f_flag;
		BoolVector new_v_flag;
		BoolVector region_f_flag;
		BoolVector region_v_flag;
		BoolVector cutv_flag;
		BoolVector cutf_flag;
		std::vector<int> grow_dir;

		std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];
#if USE_NEW_SIMILARITY_ENERGY
		Eigen::Matrix3Xd normal_sampling;
		bool has_ns = false;
#else
		Eigen::VectorXd normal_similarity_angle;
		bool has_nsa = false;
#endif


		void run(const Eigen::Matrix3Xd &normal);
		void modify_cut();
	};
}