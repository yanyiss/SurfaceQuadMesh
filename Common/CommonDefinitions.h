#ifndef COMMON_DEFINITION_H
#define COMMON_DEFINITION_H

#include <limits>

namespace CommonDefinition
{
	enum COLOR_Definition{RED=0,BLUE,GREEN,PURPLE,ORANGE,N_COLOR};
	extern const int colorDisplay[5][3];
	extern const char Color_String[5][512];
	extern const int PointColor[2][3];
	extern const char PointColor_String[2][512];
	extern const int SingularityColor[12][3];

	enum TRI_DISTORTION_TYPE
	{
		CONFORMAL_D = 0,
		ISOMETRIC_D,
		LSCM_D,
	};
	extern const char distortion_type_name[3][128];

	enum TRI_MAPPING_ENERGY
	{
		MIPS_AREA_N =0,
		MIPS_AREA_G,
		LSCM,
		EXP_MIPS_AREA_G,
		EXP_MIPS_AREA_G_B,
		EXP_MIPS_AREA_N,
		EXP_MIPS_AREA_G_Untangling,
		MIPS_AREA_G_S,
		LSCM_S,
		EXP_MIPS_AREA_G_S,
		LCOT_ISOMRTEIC,
		DILATATION_G,
		DILATATION_N,
	};
	extern const char tri_mapping_energy[13][128];

	enum NODE_TYPE{ SMOOTH_NODE = 0, CREASE_NODE, CORNER_NODE};

	enum BOUND_SMOOTH_ERROR_METHOD
	{
		BSE_LSCM = 0,
		BSE_LSCM_BETA,
		BSE_LSCM_BOUND,
	};
	extern const char bound_smooth_error_method[3][128];
}

#include <OpenMesh/Core/Geometry/VectorT.hh>
struct local_frame 
{
	local_frame()
		:e_x(OpenMesh::Vec3d(1,0,0)), e_y(OpenMesh::Vec3d(0,1,0)), n(OpenMesh::Vec3d(0,0,1))
	{}
	~local_frame(){}

	void find_e_x_y()
	{
		if(std::abs(n[2]) >= std::abs(n[1]) && std::abs(n[2]) >= std::abs(n[0]))
		{
			e_x[0] = 1.0; e_x[1] = 1.0; e_x[2] = ( -n[0] - n[1] ) / n[2];
		}
		else if (std::abs(n[1]) >= std::abs(n[2]) && std::abs(n[1]) >= std::abs(n[0]))
		{
			e_x[0] = 1.0; e_x[2] = 1.0; e_x[1] = ( -n[0] - n[2] ) / n[1];
		}
		else 
		{
			e_x[1] = 1.0; e_x[2] = 1.0; e_x[0] = ( -n[2] - n[1] ) / n[0];
		}
		e_x.normalize();
		e_y = OpenMesh::cross(n, e_x);
	}

	OpenMesh::Vec3d e_x;
	OpenMesh::Vec3d e_y;
	OpenMesh::Vec3d n;
};

#endif