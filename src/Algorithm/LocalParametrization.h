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
		inline double GetV(int vlid) {
			if (vidmap[vlid] >= int(uv[1].size()))
			{
				dprint(vlid, vidmap[vlid], vidmap.size(), uv[1].size());
				dprint(new_v_flag[vlid], region_v_flag[vlid]);
				int p = 0;
			}
			return uv[1](vidmap[vlid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }


		//private:
		M4* m4;

		std::vector<VertexLayer*> new_vertex;
		std::vector<FaceLayer*> new_face;
		std::vector<VertexLayer*> region_vertex;
		std::vector<FaceLayer*> region_face;
		std::vector<VertexLayer*> cut;

		/*std::deque<bool> new_f_flag;
		std::deque<bool> new_v_flag;
		std::deque<bool> region_f_flag;
		std::deque<bool> region_v_flag;
		std::deque<bool> cutv_flag;
		std::deque<bool> cutf_flag;*/
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
		Eigen::VectorXd normal_similarity_angle;
		bool has_nsa = false;

		void run(const Eigen::Matrix3Xd &normal);
		void modify_cut();
	};
}