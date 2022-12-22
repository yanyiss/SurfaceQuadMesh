#include "LoopToolbox.h"
namespace LoopGen
{
	bool LoopGen::FieldAligned_PlanarLoop(VertexHandle v, std::vector<int>& loop, int shift)
	{
		// shift only 0 & 1
		// һ��ʼ����root��ջ������ΧһȦ��ջ��Ȼ���������Ϊdblmax���Ϳ���һֱ˳�ŷ�����
		// ���ɻ�
		// �ؼ��������ڣ�����ҵ��˷���ĳͷ�ֵ�������ϣ����ҵĻ�����

		int nv = mesh->n_vertices();
		int vid = v.idx();
		std::deque<bool> visited(nv, false);
		std::vector<double> distance(nv, DBL_MAX);
		std::vector<HalfedgeHandle> prev(nv);
		shift %= 2;

		struct VertexPQ
		{
			int id;
			int shift;
			HalfedgeHandle prev;
			double dist; 
			// ��˭�ľ���
			int count;
			VertexPQ() {}
			VertexPQ(int id_, int shift_, double dist_, int count_) :id(id_), shift(shift_), dist(dist_), count(count_) {}
			bool operator>(const VertexPQ& x) const { return dist > x.dist; }
		};

		std::priority_queue<VertexPQ, std::vector<VertexPQ>, std::greater<VertexPQ>> pq;

		std::vector<int> count(nv, 0);
		// ���count�Ǹ�ʲô�ģ���¼ÿ������ʵĴ�����

		Eigen::Vector3d plane_normal(0, 0, 0);
		auto& crossfield = cf->getCrossField();
		auto& matching = cf->getMatching();
		auto& position = cf->getPosition();

		double halfPI = PI * 0.5;
		double doublePI = PI * 2.0;
		double triple_halfPI = halfPI * 3.0;

		// ��ʼ������ v ��outgoing����㣬Ȼ�󣬸�һ����ʼֵ������shiftû����
		for (auto voh = mesh->voh_begin(v); voh != mesh->voh_end(v); ++voh)
		{
			// outgoing ��ߵĵ�
			int fid = mesh->face_handle(voh.handle()).idx();
			plane_normal += crossfield.col(4 * fid + shift + 1);
			// ��ֱ�Ӽӵ�ÿ����ĵ�һ���ߵڶ�����
			if (weight(shift, voh->idx()) < 2)
			{
				int toid = mesh->to_vertex_handle(voh.handle()).idx();
				distance[toid] = weight(shift, voh->idx());
				pq.emplace(toid, shift, distance[toid], ++count[toid]);
				prev[toid] = voh.handle();
			}
			shift += matching[voh->idx()];
			shift %= 4;
		}
		plane_normal.normalize();
		// ������������1�����ƽ���normal

		while (true)
		{
			VertexPQ vert;
			do
			{
				if (pq.empty())
				{ 
					return false;
				}
				vert = pq.top();
				pq.pop();
			} while (vert.count != count[vert.id]);
			// Ӧ���������Ҷ�ջ����ĵ�һ������

			//loop.push_back(vert.id); loop.push_back(mesh->from_vertex_handle(prev[vert.id]).idx());
			int fromid = vert.id;
			visited[fromid] = true;
			if (fromid == vid)
			{
				break;
			}
			auto voh = mesh->next_halfedge_handle(prev[fromid]);
			int valence = mesh->valence(mesh->vertex_handle(fromid)) - 1;
			shift = vert.shift;
			// ΪʲôshiftҪ�䣬��Ϊshift��������Ҫ���ĸ�����ʼ�ң�����Ҫȥ��һ����ĸ��������ֵ
			// �Ҿ�������������ǶԵģ���Ϊ��֮ǰ����ȫ�ֵķ���ȥ�ң��������ͻ�䣬�������൱���ھֲ��Һ��ʵķ��򣬲�ѡ����ʵı�
			for (int i = 0; i < valence; ++i)
			{
				double w = weight(shift, voh.idx()); w *= w;
				if (w < 2)
				{
					int toid = mesh->to_vertex_handle(voh).idx();
					//double dot_ = (position.col(toid) - position.col(fromid)).normalized().dot(plane_normal); dot_ *= dot_;
					//w = sqrt(w + dot_);
					w = sqrt(w * 900 + 1 - w);
					if (distance[fromid] + w < distance[toid])
					{
						distance[toid] = distance[fromid] + w;
						pq.emplace(toid, shift, distance[toid], ++count[toid]);
						prev[toid] = voh;
					}
				}
				shift += matching[voh.idx()]; shift %= 4;
				voh = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(voh));
			}
		}

		loop.clear();
		loop.push_back(vid);
		int previd = mesh->from_vertex_handle(prev[vid]).idx();
		while (previd != vid)
		{
			loop.push_back(previd);
			previd = mesh->from_vertex_handle(prev[previd]).idx();
		}
		loop.push_back(vid);
		// ���˸����ϱ����ߵĵ�Ľṹ�����·��
		return true;
	}

	void LoopGen::InitializeField()
	{
#if 0
		cf = new crossField(mesh);
#else
		cf = new crossField("..//resource//field//vase.field");
#endif
		dprint("Initialize Field Done!");
	}

	void LoopGen::InitializeGraphWeight()
	{
		auto& crossfield = cf->getCrossField();
		auto& matching = cf->getMatching();
		auto& normal = cf->getNormal();
		auto& position = cf->getPosition();

		double doublePI = 2.0 * PI;
		double halfPI = 0.5 * PI;
		double _PI = -1.0 * PI;
		double _halfPI = -0.5 * PI;

		weight.resize(4, mesh->n_halfedges());
		//�ĸ�ֵ��ʲô���ĸ������weight
		for (auto eitr = mesh->edges_begin(); eitr != mesh->edges_end(); ++eitr)
		{
			auto h0 = mesh->halfedge_handle(eitr.handle(), 0);
			auto h1 = mesh->halfedge_handle(eitr.handle(), 1);
			auto fid = mesh->face_handle(h0).idx();
			auto gid = mesh->face_handle(h1).idx();
			auto& fv = crossfield.col(fid * 4);
			auto& gv = crossfield.col(gid * 4 + matching[h0.idx()]);
			auto ev = position.col(mesh->to_vertex_handle(h0).idx()) - position.col(mesh->from_vertex_handle(h0).idx());

			// �����ÿ������corssfield��һ�������빫���ߵļнǣ�atan2 ���������ޣ������Ƕ�
			double arc0 = atan2(ev.cross(fv).dot(normal.col(fid)), ev.dot(fv)); //arc0 += arc0 > 0 ? 0 : doublePI;
			double arc1 = atan2(ev.cross(gv).dot(normal.col(gid)), ev.dot(gv)); //arc1 += arc1 > 0 ? 0 : doublePI;
			double arc = atan2(sin(arc0) + sin(arc1), cos(arc0) + cos(arc1));
			// ���arc��ɶ����ô��⣿��LC�ĽǶ���?�������淽�򳡵Ĳ�ֵ��

			auto& w0 = weight.col(h0.idx());
			auto& w1 = weight.col(h1.idx());
			double s = fabs(sin(arc));
			double c = fabs(cos(arc));

			/*
			if (eitr->idx() == 6097)
			{
				dprint(arc0/PI *180);
				dprint(arc1 / PI * 180);
				dprint(arc / PI * 180);
				dprint("eitr0");
				dprint(fv.dot(ev.normalized()), gv.dot(ev.normalized()), crossfield.col(4 * gid).dot(ev.normalized()));
				dprint(mesh->from_vertex_handle(h0).idx(), mesh->to_vertex_handle(h0).idx());
				dprint(matching[h0.idx()]);
				int p = 0;
			}
			*/

			if (s < c)
				c = DBL_MAX;
			else
				s = DBL_MAX;
			if (arc >= 0 && arc < halfPI)
			{
				w0 << s, DBL_MAX, DBL_MAX, c;
				//w1 << DBL_MAX, c, s, DBL_MAX;
			}
			else if (arc >= halfPI && arc < PI)
			{
				w0 << DBL_MAX, DBL_MAX, s, c;
				//w1 << s, c, DBL_MAX, DBL_MAX;
			}
			else if (arc >= _PI && arc < _halfPI)
			{
				w0 << DBL_MAX, c, s, DBL_MAX;
				//w1 << s, DBL_MAX, DBL_MAX, c;
			}
			else
			{
				w0 << s, c, DBL_MAX, DBL_MAX;
				//w1 << DBL_MAX, DBL_MAX, s, c;
			}
			switch (matching[h0.idx()])
			{
			case 0:
				w1 << w0(2), w0(3), w0(0), w0(1);
				break;
			case 1:
				w1 << w0(1), w0(2), w0(3), w0(0);
				break;
			case 2:
				w1 << w0(0), w0(1), w0(2), w0(3);
				break;
			case 3:
				w1 << w0(3), w0(0), w0(1), w0(2);
				break;
			}
		}
		dprint("Initialize Graph Weight Done!");
		// ������˼�������ĸ�����ģ���ͬ��ߵ�ֵ
		// �൱�ڰ���ӽ��ߵ��Ǹ�����������ˣ�����һ��ֵ����������doublemax
	}

	void LoopGen::InitializePQ()
	{
		//compute_principal_curvature(m, cur[0], cur[1], dir[0], dir[1]);
		for (int i = 0; i < 2; ++i)
		{
			//for (auto &k : cur[i])
				//k = fabs(k);
		}
		for (auto v : mesh->vertices())
		{
			auto vid = v.idx();
			EnergyOnVertex ev;
			ev.v = v;
			//if (cur[0][vid] < 1.0e-6 || cur[1][vid] < 1.0e-6)
			{
				ev.energy = DBL_MAX;
				pq.push(ev);
				continue;
			}


			//FieldAligned_PlanarLoop(v, dir[cur[0][vid] >= cur[1][vid] ? 0 : 1][vid]);
			LocalParametrization lp;
			ev.energy = 0;
			pq.push(ev);
		}
	}

#pragma region 
	void LoopGen::cal_M4_edge_Fd()
	{
		// ÿһ������ÿһ��Ĳ�ֵ�ķ��򳡣�ϣ��������edge.idx() + current_direction_index
		edge_fd.resize(3, mesh->n_edges() * 4);

		auto& crossfield = cf->getCrossField();
		auto& matching = cf->getMatching();

		for (auto eitr = mesh->edges_begin(); eitr != mesh->edges_end(); ++eitr)
		{
			auto eid = eitr.handle().idx();
			auto h0 = mesh->halfedge_handle(eitr.handle(), 0);
			auto h1 = mesh->halfedge_handle(eitr.handle(), 1);
			auto fid = mesh->face_handle(h0).idx();
			auto gid = mesh->face_handle(h1).idx();

			for (size_t direction = 0; direction < 4; direction++)
			{
				// ���ĸ�������о�ֵ ָ��ʱ�뷽�� ����ϣ����һ�µķ���ӵ�һ��
				auto& fv = crossfield.col(fid * 4 + (direction + face_M4layer_index_begin[fid] ) % 4);
				auto& gv = crossfield.col(gid * 4 + (direction + face_M4layer_index_begin[gid] ) % 4);
				auto temp_fd = (fv + gv) / 2;
				edge_fd.col(eid * 4 + direction) = temp_fd;
			}
		}
		/*
		for (size_t i = 0; i < 10; i++)
		{
			dprint(edge_fd.col(i));
		}
		*/
	}


	void LoopGen::cal_local_node() {
		// ��������ȵĹ�ȱ��������ÿ�����k-�����
		// ��ÿ�����dijkstra�����k���������·��������ȫ�ֵ�Pͼ(��������Ȩ��)
		// ��ʼ��Local_node,ȷ��ÿ�����k�����
		Local_node.clear();
		Local_node.resize(mesh->n_vertices());
		for (auto vitr = mesh->vertices_begin(); vitr != mesh->vertices_end(); ++vitr) {
			// ��ÿһ��v��������k���򣬴洢��Local_node
			Local_node[(*vitr).idx()] = cal_K_ring_node((*vitr));
			// ��ÿһ��v��dijkstra�����k���������·��
			cal_K_ring_dijkstra((*vitr), Local_node[(*vitr).idx()]);
		}
	}
	void LoopGen::cal_Weighted_LocalGraph(int current_direction_index)
	{
		// ��ÿһ��k-ring�����������ľ���Ȩ��
		this->current_direction_index = current_direction_index;
		for (auto vitr = mesh->vertices_begin(); vitr != mesh->vertices_end(); ++vitr) {
			for (auto& node : Local_node[(*vitr).idx()])
			{
				if (node.v_id == (*vitr).idx()) {
					continue;
				}
				auto newdistance = cal_node_distance((*vitr), node.v, Local_node[(*vitr).idx()]);
				node.distance = newdistance;
			}

		}
		dprint("cal_Weighted_LocalGraph done");
	}

	std::vector<LoopGen::GNode> LoopGen::get_global_node() {
		return Global_node;
	}

	std::vector<LoopGen::LNode> LoopGen::cal_K_ring_node(VertexHandle v) {
		// ����v��k-ring node��������һ����ʼ���õ�LNode vector
		std::vector<LoopGen::LNode> node_stack;
		
		// ������ڵ� 
		LoopGen::LNode n = {v,v.idx(),0,OpenMesh::SmartHalfedgeHandle(),false,0};
		int k_vertex_num = 0;
		node_stack.push_back(n);
		while (node_stack[k_vertex_num].deepth != k_default) {
			// ��ȡ��ǰnode
			auto tmp_node = node_stack[k_vertex_num];
			// �����node���������stack����node
			for (auto v_handle : mesh->vv_range(tmp_node.v)) {
				if (find_node_in_vector(v_handle, node_stack) == -1) {
					// û�ҵ������node
					node_stack.push_back(LoopGen::LNode({ v_handle,v_handle.idx(),DBL_MAX,OpenMesh::SmartHalfedgeHandle(),false,tmp_node.deepth + 1 }));
				}
				else {

				}
			}
			//���º��
			k_vertex_num++;
		}
		while (node_stack.size() != k_vertex_num) {
			node_stack.pop_back();	
		}
		// std::cout << "begin vertex" << v.idx() << "k-ring num:" << k_vertex_num << std::endl;
		return node_stack;
	}

	void LoopGen::cal_K_ring_dijkstra(VertexHandle root, std::vector<LoopGen::LNode>& node_vector) {
		// ��ָ����v ��������k-ring�ľ��룬��һ���Ծ���·��������£�֮�����������ľ���
		// ÿ����һ��ûȷ������̾���ĵ㣬ȷ������������������������ľ���
		// n_vectorֻ���˵㣬id �� deepth��halfedge��distance��fix���ǿ�

		int k_num = node_vector.size();
		for (int fix_num = 0; fix_num < k_num;  fix_num++) {
			//  ȡ�� tmp_node
			double tmp_dis = DBL_MAX;
			int tmp_min_id = 0;
			for (size_t i = 0; i < k_num; i++)
			{
				if (node_vector[i].have_fixed == false && node_vector[i].distance < tmp_dis) {
					tmp_dis = node_vector[i].distance;
					tmp_min_id = i;
				}
			}
			auto& tmp_node = node_vector[tmp_min_id];
			tmp_node.have_fixed = true;

			// ����tmp_node 1-ring ���� �µĵ�
			for (auto voh : mesh->voh_range(tmp_node.v)) {
				// ���v_handle��node_vector����¾��롢halfedge
				auto v_handle = voh.to();
				int update_index = find_node_in_vector(v_handle, node_vector);
				if (update_index != -1 && node_vector[update_index].have_fixed==false) {
					auto& update_node = node_vector[update_index];
					// �ȽϾ��� ������ ��ͶӰ
				    // double newdistance = cal_node_distance(root, node_vector, v_handle, voh);

					// �ȽϾ��� ��������
					if (update_node.distance >= tmp_node.distance+mesh->calc_edge_length(voh)) {
						update_node.distance = tmp_node.distance + mesh->calc_edge_length(voh);
						update_node.from_halfedge = voh;
					}
				}
			}
		}
	}

	double LoopGen::cal_node_distance(VertexHandle root, VertexHandle target, std::vector<LNode>& n_vector)
	{
		// ���ߵ�ֱ�߾����ͶӰ�����ѷ��򳡼ӹ�ȥ��crossfield��edge_fd���ǹ�һ���Ľ��
		Vec3d e_origin = mesh->point(target) - mesh->point(root);
		double length_e = e_origin.norm();
		Eigen::Matrix<double, 3, 1> e_fd = { 0,0,0 };
		int index = find_node_in_vector(target, n_vector);
		if (index == -1) {
			dprint("û�ҵ�");
			return DBL_MAX;
		}
		auto heh = n_vector[index].from_halfedge;
		do
		{
			auto v0 = heh.from();
			auto v1 = heh.to();
			Vec3d e01 = mesh->point(v1) - mesh->point(v0);
			double w = e01.dot(e_origin) / length_e;
			e_fd = e_fd + w * edge_fd.col(mesh->edge_handle(heh).idx()*4 + current_direction_index);
			if (v0.idx() == root.idx()) {
				break;
			}
			else{
				int tmp_index = find_node_in_vector(v0, n_vector);
				heh = n_vector[tmp_index].from_halfedge;
			}

		} while (true);

		Eigen::Matrix<double, 3, 1> e = {e_origin[0],e_origin[1],e_origin[2] };
		e_fd.normalize();
		// �жϷ�����������ˣ�����DBLMAX
		if (e_fd.dot(e)/(e_fd.norm()*e.norm()) < 0.70710) {
		//if(e_fd.dot(e) / (e_fd.norm() * e.norm()) < 0){
			// ϣ��������45����
			return DBL_MAX;
		}
		double alpha = 30;
		double distance = pow((e.transpose()*e_fd),2) + alpha* alpha * ((e.transpose()*e)-pow((e.transpose()*e_fd),2));
		distance = sqrt(distance);
		return distance;
	}

	std::vector<OpenMesh::SmartHalfedgeHandle> LoopGen::cal_k_ring_path(VertexHandle root, VertexHandle target, std::vector<LoopGen::LNode>& node_vector) {
		std::vector<OpenMesh::SmartHalfedgeHandle> path;
		int target_id = find_node_in_vector(target, node_vector);
		auto heh = node_vector[target_id].from_halfedge;
		do
		{
			path.push_back(heh);
			auto v0 = heh.from();
			if (v0.idx() == root.idx()) {
				break;
			}
			else {
				int index = find_node_in_vector(v0, node_vector);
				heh = node_vector[index].from_halfedge;
			}

		} while (true);
		return path;
	}

	int LoopGen::find_node_in_vector(VertexHandle v, std::vector<LoopGen::LNode>& n_vector) {
		int id = v.idx();
		int index = 0;
		for (auto tn : n_vector) {
			if (tn.v_id == id) {
				return index;
			}
			index++;
		}
		return -1;
	}

	void LoopGen::calM4Layers(VertexHandle root)
	{
		// һ������ʵ�����⣺������ҵ����ʼ�ķ��򣬺�ÿ����ķ���Ķ�Ӧ��ϵ��
		// �㻹�ǻ�����������������ϣ�Ӧ��ֻ���ͽ���Ƭ�ķ����·���Ĺ�ϵ
		 
		face_M4layer_index_begin.clear();
		// ��faceid=0��ʼ��ͨ��matchings��������������+��
		// ָ�������ĵ�һ�������Ӧ�ķ������ĸ�����
		int tmp_face_num = mesh->n_faces();
		face_M4layer_index_begin.resize(tmp_face_num);
		for (auto i : face_M4layer_index_begin) i = 0;

		auto begin_face = mesh->face_handle(mesh->voh_begin(root));
		int begin_id = begin_face.idx();
		face_M4layer_index_begin[begin_id] = 0;

		std::vector<bool> face_has_checked(tmp_face_num);
		for (size_t i = 0; i < face_has_checked.size(); i++) { face_has_checked[i] = false; }
		face_has_checked[begin_id] = true;

		std::vector<int> face_pool;
		face_pool.push_back(begin_id);
		auto matchings = cf->getMatching();

		while (face_pool.empty() != true) {
			// �ó���һ��
			int tmp_id = face_pool[0];
			face_pool.erase(face_pool.begin());
			// ������һ��
			auto fh = mesh->face_handle(tmp_id);
			for (auto it = mesh->fh_begin(fh); it != mesh->fh_end(fh); ++it) {
				auto he = it.handle();
				int new_face_index = mesh->face_handle(he.opp()).idx();

				if (face_has_checked[new_face_index] == false) {
					face_M4layer_index_begin[new_face_index] = (face_M4layer_index_begin[tmp_id] + matchings[he.idx()]) % 4;
					// ���ڴ���ǣ�ÿ�������һ����һ�µ��Ǹ���������
					face_pool.push_back(new_face_index);
					face_has_checked[new_face_index] = true;
				}
			}
			//Ϊ����������
		}
		/*
		for (auto m : face_M4layer_index_begin) {
			std::cout << m << std::endl;
		}
		*/
		return ;
	}

	std::vector<int>  LoopGen::get_face_M4layer_index_begin() {
		return face_M4layer_index_begin;
	}

	void LoopGen::cal_global_dijkstra(VertexHandle root, int current_direction_index)
	{
		// ����ȫ��dijkstra������С������浽global_node��
		
		//��ʼ��
		this->current_direction_index = current_direction_index;
		this->Global_node.clear();
		this->Global_node.resize(mesh->n_vertices());

		std::vector<bool> global_node_is_fixed;
		global_node_is_fixed.resize(mesh->n_vertices(), false);

		std::vector<double> global_distance;
		global_distance.resize(mesh->n_vertices(), DBL_MAX);
		global_distance[root.idx()] = 0;
		
		std::priority_queue<GNode, std::vector<GNode>, std::greater<GNode>> heap;// �����ȶ��У���������н��
		heap.push(GNode({ root.idx(),root.idx(),0 }));

		while (heap.size())
		{
			// ȡ��������С��δ�̶���Ԫ��
			auto tmp_gnode = heap.top();
			heap.pop();
			if (global_node_is_fixed[tmp_gnode.v_id]==true) continue;
			global_node_is_fixed[tmp_gnode.v_id] = true;

			// �浱ǰ�Ľ��
			this->Global_node[tmp_gnode.v_id]=tmp_gnode;
			global_distance[tmp_gnode.v_id] = tmp_gnode.distance;

			// ���ݵ�ǰ�����¶���
			for (auto lnode:Local_node[tmp_gnode.v_id])
			{
				if (global_node_is_fixed[lnode.v_id] == true) {
					continue;
				}
				// ����v_id��kring , ÿһ��lnode����tmp_gnode�ܵ��ĵ�
				if (global_distance[lnode.v_id] > tmp_gnode.distance + lnode.distance) {
					global_distance[lnode.v_id] = tmp_gnode.distance + lnode.distance;
					heap.push(GNode({ lnode.v_id,tmp_gnode.v_id,global_distance[lnode.v_id] }));
				}
			}
		}
		dprint("cal_global_dijkstra end");
	}

	std::vector<OpenMesh::SmartHalfedgeHandle> LoopGen::cal_global_path(VertexHandle root, VertexHandle target) {
		// ��Ҫ���� ȫ�ֵ�dijkstra������rootҪһ��
		std::vector<OpenMesh::SmartHalfedgeHandle> path;
		int target_id = target.idx();
		int root_id = root.idx();
		auto tmp_gnode = Global_node[target_id];
		if (tmp_gnode.v_id == 0 && tmp_gnode.last_v_id == 0) {
			return path;
			// bug can not find last v_id
		}

		while (tmp_gnode.v_id!=root_id)
		{
			// ���㵱ǰ��·�������ӵ���·����
			int tmp_target_id = tmp_gnode.v_id;
			int tmp_root_id = tmp_gnode.last_v_id;

			auto tmp_path = cal_k_ring_path(mesh->vertex_handle(tmp_root_id), mesh->vertex_handle(tmp_target_id),Local_node[tmp_root_id]);
			path.insert(path.end(), tmp_path.begin(), tmp_path.end());

			// �����µ�tmp_gnode
			tmp_gnode = Global_node[tmp_root_id];
		}
		return path;


	}
#pragma endregion
}


	