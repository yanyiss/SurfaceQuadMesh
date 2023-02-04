#include "crossField.h"
#include "LoopDef.h"
namespace LoopGen
{
	class LocalParametrization
	{
	public:
		LocalParametrization(M4& m4_, VertexLayer* vl_);
		~LocalParametrization() {};
	public:
		inline double GetRegularU(int vid) { return uv[0](vidmap[vid]) - std::floor(uv[0](vidmap[vid])); }
		inline double GetU(int vid) { return uv[0](vidmap[vid]); }
		inline double GetV(int vid) { return uv[1](vidmap[vid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }


		//private:
		M4* m4;

		std::vector<VertexLayer*> new_vertex;
		std::vector<FaceLayer*> new_face;
		std::vector<VertexLayer*> region_vertex;
		std::vector<FaceLayer*> region_face;
		std::vector<VertexLayer*> cut;

		std::deque<bool> new_f_flag;
		std::deque<bool> new_v_flag;
		std::deque<bool> region_f_flag;
		std::deque<bool> region_v_flag;
		std::deque<bool> cutv_flag;
		std::deque<bool> cutf_flag;
		std::vector<int> grow_dir;

		std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];
		Eigen::VectorXd normal_similarity_angle;
		bool has_nsa = false;

		void run(const Eigen::Matrix3Xd &normal);
	};
}