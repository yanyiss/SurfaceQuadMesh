#include "MeshField.h"
#include "..\src\Toolbox\dprinter\dprint.h"
#define grainSize 1024
#define epsilon 1.1e-15

MeshField::hierarchy* MeshField::h = nullptr;

#if 0
#pragma region MeshField
MeshField::hierarchy::hierarchy(const TriMesh &mesh, ui depth)
{
	//以 mesh 初始化该层
	init(mesh);
	if (depth < 2)
	{
		hierarchy_depth = 1;
		return;
	}

	//循环初始化剩下各层
	finer_hierarchy = nullptr;
	coarser_hierarchy = new hierarchy(this);
	hierarchy* h = coarser_hierarchy;
	h->finer_hierarchy = this;
	hierarchy_depth = depth;
	for (; depth > 2 && h->v.size() > 1; depth--, h = h->coarser_hierarchy)
	{
		h->coarser_hierarchy = new hierarchy(h);
		h->coarser_hierarchy->finer_hierarchy = h;
	}
	h->coarser_hierarchy = nullptr;
	hierarchy_depth -= depth - 2;

	//初始化各层 hierarchy_depth
	depth = hierarchy_depth;
	h = this->coarser_hierarchy;
	while (h)
	{
		h->hierarchy_depth = --depth;
		h = h->coarser_hierarchy;
	}
}

double MeshField::hierarchy::calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv)
{
	double area = 0;
	Vec3d pos = mesh.point(tv);
	std::vector<Link> &ad = adj[tv.idx()];
	ad.reserve(mesh.valence(tv));
	for (OH tvoh : mesh.voh_range(tv))
	{
		ad.push_back(mesh.to_vertex_handle(tvoh).idx());
		if (!mesh.face_handle(tvoh).is_valid())
			continue;
		Vec3d s = mesh.point(mesh.to_vertex_handle(tvoh));
		Vec3d t = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(tvoh)));
		Vec3d c = (pos + s + t) / 3.0;
		area += 0.25*((s - pos).cross(c - pos).norm() + (t - pos).cross(c - pos).norm());
	}
	return area;
}

void MeshField::hierarchy::init(const TriMesh &mesh)
{
	ui nv = mesh.n_vertices();
	v.reserve(nv);
	n.reserve(nv);
	a.reserve(nv);
	adj.resize(nv);
	for (OV tv : mesh.vertices())
	{
		ui id = tv.idx();
		v.push_back(mesh.point(tv));
		n.push_back(mesh.calc_normal(tv));
		a.push_back(calcVoronoiAndAdjacency(mesh, tv));
	}
	graph_color();
}

void MeshField::hierarchy::graph_color()
{
	ui nsizes = adj.size();
	std::vector<int> color(nsizes, -1);
	std::vector<ui> possible_color;

	ui ncolors = 0;
	for (size_t i = 0; i < nsizes; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1);
		for each (auto link in adj[i])
		{
			int c = color[link.id];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}

	phase.resize(ncolors);
	for (size_t i = 0; i < nsizes; i++)
	{
		phase[color[i]].push_back(i);
	}
}

#include <deque>
MeshField::hierarchy::hierarchy(hierarchy *fh)
{
	const vv3d &v_ = fh->v;
	const vv3d &n_ = fh->n;
	const std::vector<double> &a_ = fh->a;
	const vvL &adj_ = fh->adj;

	std::vector<edgeWeight> edge;
	ui size = 0;
	for (const std::vector<Link> &l : adj_)
		size += l.size();
	edge.reserve(size);

	size = 0;
	for (const std::vector<Link> &v0 : adj_)
	{
		for (Link v1 : v0)
			edge.emplace_back(size, v1.id, n_[size].dot(n_[v1.id])*std::max(a_[size] / a_[v1.id], a_[v1.id] / a_[size]));
		size++;
	}

	std::sort(edge.begin(), edge.end(), std::greater<edgeWeight>());
	std::deque<bool> mergeFlag(size, false);
	size = 0;
	for (const edgeWeight &e : edge)
	{
		if (mergeFlag[e.i] || mergeFlag[e.j])
			continue;
		mergeFlag[e.i] = true;
		mergeFlag[e.j] = true;
		edge[size++] = e;
	}

	ui nVertices = fh->v.size();
	ui leftVertices = nVertices - size;
	v.reserve(leftVertices);
	n.reserve(leftVertices);
	a.reserve(leftVertices);
	adj.resize(leftVertices);

	tofiner.reserve(leftVertices);
	std::vector<ui> &tc = fh->tocoarser;
	tc.resize(nVertices);

	for (ui id = 0; id < size; id++)
	{
		edgeWeight &e = edge[id];
		double a0 = a_[e.i];
		double a1 = a_[e.j];
		double area = a0 + a1;
		if (area > epsilon)
		{
			v.push_back((v_[e.i] * a0 + v_[e.j] * a1) / area);
		}
		else
		{
			v.push_back((v_[e.i] + v_[e.j]) / 2.0);
		}

		n.push_back((n_[e.i] * a0 + n_[e.j] * a1).normalized());
		if (isnan(n.back()[0]))
		{
			dprint(e.i, e.j);
			dprint(n_[e.i], n_[e.j]);
			int p = 0;
		}
		a.push_back(area);
		tofiner.emplace_back(e.i, e.j);
		tc[e.i] = id; tc[e.j] = id;
	}

	for (ui id = 0; id < nVertices; id++)
	{
		if (mergeFlag[id])
			continue;
		v.push_back(v_[id]);
		n.push_back(n_[id]);
		a.push_back(a_[id]);
		tofiner.emplace_back(id, nVertices);
		tc[id] = size;
		size++;
	}
	assert(size == leftVertices);
	for (ui id = 0; id < leftVertices; id++)
	{
		std::vector<Link> &oneRing = adj[id];
		//std::vector<Link> oneRing;
		//导入精细层网格的one-ring
		ui i = tofiner[id].first;
		for (const Link &var : adj_[i])
		{
			if (tc[var.id] == id) continue;
			oneRing.emplace_back(tc[var.id], var.weight);
		}
		i = tofiner[id].second;
		if (i != nVertices)
		{
			for (const Link &var : adj_[i])
			{
				if (tc[var.id] == id) continue;
				oneRing.emplace_back(tc[var.id], var.weight);
			}
		}
		if (oneRing.size() <= 1)
			continue;

		//删除oneRing重复元素
		std::sort(oneRing.begin(), oneRing.end());
		std::vector<Link>::iterator itr = oneRing.begin();
		itr++;
		for (; itr != oneRing.end(); itr++)
		{
			std::vector<Link>::iterator newItr = itr;
			newItr--;
			if (itr->id == newItr->id)
			{
				itr->weight += newItr->weight;
				oneRing.erase(newItr);
				itr--;
			}
		}
	}
	graph_color();
}

MeshField::hierarchy::~hierarchy()
{
	finer_hierarchy = nullptr;
	if (hierarchy_depth >> 1)
	{
		delete coarser_hierarchy;
		coarser_hierarchy = nullptr;
	}
}
#pragma endregion



#pragma region OrientationField
void OrientationField::loopDivision(TriMesh &mesh, ui divideTimes)
{
	if (divideTimes < 1)
		return;
	OpenMesh::Subdivider::Uniform::LoopT<TriMesh> loopDivider;
	loopDivider.attach(mesh);
	loopDivider(divideTimes);
	loopDivider.detach();
	mesh.update_normals();
}

void OrientationField::randomInit(TriMesh &mesh)
{
	loopDivision(mesh, 1);

	h = new hierarchy(mesh, DefaultDepth);
	hierarchy *hw = h;

	while (hw->coarser_hierarchy)
	{
		hw->o.resize(hw->v.size());
		hw = hw->coarser_hierarchy;
	}
	ui size = hw->n.size();
	hw->o.resize(size);

	for (ui i = 0; i < size; i++)
	{
		hw->calcXAxis(hw->n[i], hw->o[i]);
	}
	dprint("Hierarchy Build Done!\n\n");
}

void OrientationField::GSOptimize()
{
	auto calcDirection = [&](const Vec3d &o0, const Vec3d &n0, const Vec3d &o1, const Vec3d &n1)
	{
		Vec3d _P1[2] = { o0,n0.cross(o0) };
		Vec3d _P2[2] = { o1,n1.cross(o1) };
		double flag_max = 0;
		ui id1 = 0, id2 = 0;
		for (ui i = 0; i < 2; i++)
		{
			for (ui j = 0; j < 2; j++)
			{
				double dottemp = _P1[i].dot(_P2[j]);
				if (std::fabs(dottemp) > std::fabs(flag_max))
				{
					flag_max = dottemp;
					id1 = i;
					id2 = j;
				}
			}
		}
		return std::make_pair(_P1[id1], _P2[id2] * (flag_max > 0 ? 1 : -1));
	};

	hierarchy* hn = h;
	while (hn->coarser_hierarchy)
		hn = hn->coarser_hierarchy;

	const ui iterationTimes = 5;
	while (true)
	{
		//在当前层 Gauss-Seidel 迭代，得到优化 Orientation Field
		vv3d &o_ = hn->o;
		ui nv = o_.size();
		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
		const vv3d &n_ = hn->n;
		const vvL &adj_ = hn->adj;

		for (int i = 0; i < iterationTimes; i++)
		{
#ifdef RUN_MESHFIELD_PARALLEL
			for (const std::vector<ui> &ph : hn->phase)
			{
				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
				tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
				{
					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
					{
						ui id = ph[phaseIdx];
#else
			for (ui id = 0; id < nv; id++)
			{
#endif
						double sum = 0;
						Vec3d &oid = o_[id];
						const Vec3d &nid = n_[id];
						for (const Link &oneRing : adj_[id])
						{
							double w = oneRing.weight;
							//dprint("je.l", oid, nid, o_[oneRing.id], n_[oneRing.id]);
							std::pair<Vec3d, Vec3d> pair_ = calcDirection(oid, nid, o_[oneRing.id], n_[oneRing.id]);
							//dprint(pair_.first, pair_.second);
							oid = pair_.first*sum + pair_.second*w;
							//dprint(oid);
							sum += w;
							oid -= nid.dot(oid)*nid;
							//dprint(oid);
							oid.normalize();
							//dprint(oid);
						}
#ifdef RUN_MESHFIELD_PARALLEL
					}
				});
#endif
			}
			//dprint("iteration times:", i);
		}

		//将当前层 Orientation Field 传递到精细层
		if (!(hn->finer_hierarchy))
			break;

		const std::vector<std::pair<ui, ui>> &tofiner = hn->tofiner;
		vv3d &o_finer = hn->finer_hierarchy->o;
		const vv3d &n_finer = hn->finer_hierarchy->n;
		const ui finer_size = n_finer.size();

#ifdef RUN_MESHFIELD_PARALLEL
		tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
#else
		for (ui id = 0; id < nv; id++)
#endif
		{
			ui finer_id = tofiner[id].first;
			o_finer[finer_id] = (o_[id] - n_finer[finer_id].dot(o_[id]) * n_finer[finer_id]).normalized();
			finer_id = tofiner[id].second;
			//if (finer_id != invalid)//大坑!!! unsigned int 不能用 finer_id > -1
			if (finer_id != finer_size)
				o_finer[finer_id] = (o_[id] - n_finer[finer_id].dot(o_[id]) * n_finer[finer_id]).normalized();
		}
#ifdef RUN_MESHFIELD_PARALLEL
		);
#endif

		hn = hn->finer_hierarchy;
	}
	dprint("\nOrientation Field Build Done!\n\n");
}
#pragma endregion


#pragma region PositionField
void PositionField::randomInit(TriMesh &mesh)
{
	if (!h)
	{
		dprint("Orientation Field need to be build first !!!");
		exit(0);
	}
	hierarchy* hw = h;
	while (hw->coarser_hierarchy)
	{
		hw->p.resize(hw->v.size());
		hw = hw->coarser_hierarchy;
	}
	hw->p = hw->v;

	length = 0;
	for (OE te : mesh.edges())
		length += mesh.calc_edge_length(te);
	length /= mesh.n_edges() / scale;
}

void PositionField::GSOptimize()
{
	auto middlePoint = [&](const Vec3d &v0, const Vec3d &n0, const Vec3d &v1, const Vec3d &n1)
	{
		double n0n1 = n0.dot(n1);
		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
	};
	double inv_length = 1.0 / length;
	auto calcPosition = [&](const Vec3d &p0, const Vec3d &v0, const Vec3d &n0, const Vec3d &o0,
		const Vec3d &p1, const Vec3d &v1, const Vec3d &n1, const Vec3d &o1)
	{
		const Vec3d mP = middlePoint(v0, n0, v1, n1);
		Vec3d y0 = n0.cross(o0);
		Vec3d origin0 = p0 + length * (o0*std::floor(o0.dot(mP - p0) * inv_length) + y0 * std::floor(y0.dot(mP - p0) * inv_length));
		Vec3d y1 = n1.cross(o1);
		Vec3d origin1 = p1 + length * (o1*std::floor(o1.dot(mP - p1) * inv_length) + y1 * std::floor(y1.dot(mP - p1) * inv_length));

		/*std::pair<Vec3d, Vec3d> position, updatePosition; 
		double dis = DBL_MAX;
		ui index[2];
		for (ui i = 0; i < 4; i++)
		{
			updatePosition.first = origin0 + ((i & 1)*o0 + ((i & 2) >> 1)*y0)*length;
			for (ui j = 0; j < 4; j++)
			{
				updatePosition.second = origin1 + ((j & 1)*o1 + ((j & 2) >> 1)*y1)*length;
				double d = (updatePosition.second - updatePosition.first).norm();
				if (d < dis)
				{
					dis = d;
					position = updatePosition;
					index[0] = i; index[1] = j;
				}
			}
		}
		return position;*/
		vv3d updateP0(4), updateP1(4);
		for (ui i = 0; i < 4; i++)
		{
			updateP0[i] = origin0 + ((i & 1)*o0 + ((i & 2) >> 1)*y0)*length;
			updateP1[i] = origin1 + ((i & 1)*o1 + ((i & 2) >> 1)*y1)*length;
		}
		double dis = DBL_MAX;
		ui index[2];
		for (ui i = 0; i < 4; i++)
		{
			for (ui j = 0; j < 4; j++)
			{
				double d = (updateP0[i] - updateP1[j]).norm();
				if (d < dis)
				{
					dis = d;
					index[0] = i; index[1] = j;
				}
			}
		}
		updateP0 = { updateP0[index[0]],updateP1[index[1]] };

		return updateP0;
	};
	auto positionRound = [&](Vec3d &p, const Vec3d &v, const Vec3d &n, const Vec3d &o,  const double length)
	{
		const Vec3d y = n.cross(o);
		p += length * (std::floor((v - p).dot(o) / length + 0.5)*o + std::floor((v - p).dot(y) / length + 0.5)*y);
	};

	hierarchy* hn = h;
	while (hn->coarser_hierarchy)
		hn = hn->coarser_hierarchy;

	const ui iterTimes = 10;
	while (true)
	{
		//在当前层 Gauss-Seidel 迭代，得到优化 Position Field
		vv3d &p_ = hn->p;
		ui nv = p_.size();
		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
		const vv3d &v_ = hn->v;
		const vv3d &n_ = hn->n;
		const vv3d &o_ = hn->o;
		const vvL &adj_ = hn->adj;

		//for (ui id = 0; id < nv; id++)
		//{
		//	double sum = 0;
		//	Vec3d &pid = p_[id];
		//	const Vec3d &vid = v_[id];
		//	const Vec3d &nid = n_[id];
		//	const Vec3d &oid = o_[id];
		//	//auto temp = pid;
		//	for (const Link &oneRing : adj_[id])
		//	{
		//		double w = oneRing.weight;
		//		ui index = oneRing.id;
		//		std::pair<Vec3d, Vec3d> pair_ = calcPosition(pid, vid, nid, oid, p_[index], v_[index], n_[index], o_[index], length);
		//		pid = pair_.first*sum + pair_.second*w;
		//		sum += w;
		//		pid /= sum;
		//		pid -= nid.dot(pid - vid)*nid;
		//	}
		//	positionRound(pid, vid, nid, oid, length);
		//}

		for (int i = 0; i < /*hn->hierarchy_depth*/iterTimes; i++)
		{
#ifdef RUN_MESHFIELD_PARALLEL
			for (const std::vector<ui> &ph : hn->phase)
			{
				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
				tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
				{
					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
					{
						ui id = ph[phaseIdx];
#else
			for (ui id = 0; id < nv; id++)
			{
#endif
						double sum = 0;
						Vec3d &pid = p_[id];
						const Vec3d &vid = v_[id];
						const Vec3d &nid = n_[id];
						const Vec3d &oid = o_[id];
						for (const Link &oneRing : adj_[id])
						{
							double w = oneRing.weight;
							ui index = oneRing.id;
							vv3d pair_ = calcPosition(pid, vid, nid, oid, p_[index], v_[index], n_[index], o_[index]);
							pid = pair_[0]*sum + pair_[1]*w;
							/*auto pair_= calcPosition(pid, vid, nid, oid, p_[index], v_[index], n_[index], o_[index]);
							pid = pair_.first*sum + pair_.second*w;*/
							sum += w;
							pid /= sum;
							pid -= nid.dot(pid - vid)*nid;
						}
						positionRound(pid, vid, nid, oid, length);
#ifdef RUN_MESHFIELD_PARALLEL
					}
				});
#endif
			}
			//dprint("iteration times:", i);
		}
        
  
//for (ui id = 0; id < nv; id++)
//{
//	double sum = 0;
//	Vec3d &pid = p_[id];
//	const Vec3d &vid = v_[id];
//	const Vec3d &nid = n_[id];
//	const Vec3d &oid = o_[id];
//	//auto temp = pid;
//	for (const Link &oneRing : adj_[id])
//	{
//		double w = oneRing.weight;
//		ui index = oneRing.id;
//		std::pair<Vec3d, Vec3d> pair_ = calcPosition(pid, vid, nid, oid, p_[index], v_[index], n_[index], o_[index], length);
//		pid = pair_.first*sum + pair_.second*w;
//		sum += w;
//		pid /= sum;
//		pid -= nid.dot(pid - vid)*nid;
//	}
//	positionRound(pid, vid, nid, oid, length);
//}

		//将当前层 Position Field 传递到精细层
		if (!(hn->finer_hierarchy))
			break;

//		hierarchy* ht = hn;
//		while (ht->finer_hierarchy)
//		{
//			const std::vector<std::pair<ui, ui>> &tofiner = ht->tofiner;
//			const vv3d &v_finer = ht->finer_hierarchy->v;
//			const vv3d &n_finer = ht->finer_hierarchy->n;
//			vv3d &p_finer = ht->finer_hierarchy->p;
//			const ui finer_size = n_finer.size();
//
//#ifdef RUN_MESHFIELD_PARALLEL
//			tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//#else
//			for (ui id = 0; id < nv; id++)
//#endif
//			{
//				ui finer_id = tofiner[id].first;
//				p_finer[finer_id] = p_[id] - n_finer[finer_id].dot(p_[id] - v_finer[finer_id])*n_finer[finer_id];
//				finer_id = tofiner[id].second;
//				if (finer_id != finer_size)
//					p_finer[finer_id] = p_[id] - n_finer[finer_id].dot(p_[id] - v_finer[finer_id])*n_finer[finer_id];
//			}
//#ifdef RUN_MESHFIELD_PARALLEL
//			);
//#endif
//			ht = ht->finer_hierarchy;
//		}

		const std::vector<std::pair<ui, ui>> &tofiner = hn->tofiner;
		const vv3d &v_finer = hn->finer_hierarchy->v;
		const vv3d &n_finer = hn->finer_hierarchy->n;
		vv3d &p_finer = hn->finer_hierarchy->p;
		const ui finer_size = n_finer.size();

#ifdef RUN_MESHFIELD_PARALLEL
		tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
#else
		for (ui id = 0; id < nv; id++)
#endif
		{
			ui finer_id = tofiner[id].first;
			p_finer[finer_id] = p_[id] - n_finer[finer_id].dot(p_[id] - v_finer[finer_id])*n_finer[finer_id];
			finer_id = tofiner[id].second;
			if (finer_id != finer_size)
				p_finer[finer_id] = p_[id] - n_finer[finer_id].dot(p_[id] - v_finer[finer_id])*n_finer[finer_id];
		}
#ifdef RUN_MESHFIELD_PARALLEL
		);
#endif
		hn = hn->finer_hierarchy;
	}
	dprint("\nPosition Field Build Done!\n\n");
}

#include <set>
void PositionField::extractMesh(PolyMesh &polymesh, std::vector<Vec3d> &pos, std::vector<VertexValue> &polyFrame)
{
	computePolyFrame(polymesh, pos, polyFrame);
	auto rrr = polyFrame;

#if 1
	for (VertexValue &pF : polyFrame)
	{
		std::list<ui> &orientedAdj = pF.orientedAdj;
		for (const ui &va : pF.adj)
			orientedAdj.push_back(va);
		Vec3d x;
		h->calcXAxis(pF.normal, x);
		Vec3d y = pF.normal.cross(x);
		orientedAdj.sort([&](const ui &id0, const ui &id1) {
			Vec3d dir0 = polymesh.point(polymesh.vertex_handle(id0)) - polymesh.point(pF.vert);
			Vec3d dir1 = polymesh.point(polymesh.vertex_handle(id1)) - polymesh.point(pF.vert);
			return std::atan2(y.dot(dir0), x.dot(dir0)) > std::atan2(y.dot(dir1), x.dot(dir1));
		});
	}

	ui size = 0;
	for (VertexValue &pF : polyFrame)
	{
		//dprint(size);
		if (size == 1)
		{
			int p = 0;
		}
		ui sizeTT = 0;
		for (const ui &oh : pF.adj)
		{
			//dprint("sizeTT:", sizeTT);
			//if (size == 280 && sizeTT == 0)
			//{
			//	int p = 0; 
			//	/*polyFrame = rrr;
			//	return;*/
			//}
			std::vector<OV> loop;
			loop.reserve(4);
			loop.push_back(polymesh.vertex_handle(oh));
			ui previousId = size;
			do
			{
				std::list<ui> &currentOrientedAdj = polyFrame[loop.back().idx()].orientedAdj;
				std::list<ui>::iterator pointer = currentOrientedAdj.begin();
				while (*pointer != previousId)
					++pointer;
				previousId = loop.back().idx();
				pointer = ++pointer == currentOrientedAdj.end() ? currentOrientedAdj.begin() : pointer;
				polyFrame[previousId].adj.erase(*pointer);
				loop.push_back(polymesh.vertex_handle(*pointer));
				if (*pointer == size) 
				{
					break;
				}
			} while (true);
			/*for (auto tf : polymesh.faces())
			{
				dprint("face:", tf.idx());
				for (auto tfv : polymesh.fv_range(tf))
					dprint(tfv.idx());
			}*/
			polymesh.add_face(loop);
			/*if (size == 12 && sizeTT == 1)
			{
				for (auto l : loop)
					dprint(polymesh.point(l));
			}*/
			sizeTT++;
		}
		pF.adj.clear();
		++size;
	}
#endif

//#pragma region MyRegion
//	//利用广度优先搜索找到回路，该回路即为polymesh的面,其中每一条边都是贪心构造
//	auto searchLoop = [&](/*PolyMesh &polymesh, std::vector<VertexValue> &polyFrame, */ui currentVertex, ui targetVertex)
//	{
//	    typedef std::pair<ui, ui> puu;
//		std::vector<puu> tree;
//		tree.reserve(16);
//		tree.emplace_back(0, targetVertex);
//		tree.emplace_back(0, currentVertex);
//		//std::vector<puu>::iterator pointer = tree.begin() + 1;
//		ui pointer = 1;
//		ui presentLocation = 1;
//		while (true)
//		{
//			//ui previousVertex = tree[pointer->first].second;
//			//ui presentVertex = pointer->second;
//			ui previousVertex = tree[tree[pointer].first].second;
//			ui presentVertex = tree[pointer].second;
//			Vec3d centralPosition = polymesh.point(polyFrame[presentVertex].vert);
//			Vec3d referenceDirection = (centralPosition - polymesh.point(polyFrame[previousVertex].vert)).normalize();
//			Vec3d forwardDirection = polyFrame[presentVertex].normal.cross(referenceDirection);
//			bool ifHold = true;
//			for (ui va : polyFrame[presentVertex].adj)
//			{
//				if (va == previousVertex)
//					continue;
//				if (va == targetVertex)
//				{
//					tree.emplace_back(presentLocation, va);
//					goto goto20210924;
//				}
//				//Vec3d testDirection = (polymesh.point(polyFrame[va].vert) - centralPosition).normalize();
//				if ((polymesh.point(polyFrame[va].vert) - centralPosition).dot(forwardDirection) < 0)
//					continue;
//				tree.emplace_back(presentLocation, va);
//				ifHold = false;
//			}
//			if (ifHold)
//			{
//				for (ui va : polyFrame[presentVertex].adj)
//				{
//					if (va == previousVertex)
//						continue;
//					tree.emplace_back(presentLocation, va);
//				}
//			}
//			presentLocation++;
//			pointer++;
//		}
//	goto20210924:;
//#if 1
//		std::vector<OV> loop;
//		loop.reserve(4);
//		ui handle = tree.back().first;
//		while (handle)
//		{
//			loop.push_back(polymesh.vertex_handle(tree[handle].second));
//			handle = tree[handle].first;
//		}
//		loop.push_back(polymesh.vertex_handle(targetVertex));
//		std::reverse(loop.begin(), loop.end());
//		if (loop.size() > 4)
//		{
//			int p = 0;
//		}
//#else
//		std::vector<ui> loop;
//		loop.reserve(5);
//		loop.push_back(targetVertex);
//		ui handle = tree.back().first;
//		while (handle)
//		{
//			loop.push_back(tree[handle].second);
//			handle = tree[handle].first;
//		}
//		loop.push_back(targetVertex);
//		if (loop.size() != 4 && loop.size() != 5)
//		{
//			dprint(loop.size());
//			//system("pause");
//			loop = std::vector<ui>(loop.rbegin(), loop.rend());
//		}
//		else
//			std::swap(loop[1], loop[loop.size() - 2]);
//#endif
//		return loop;
//	};
//	
//	ui size = 0;
//#if 1
//	for (VertexValue &pF : polyFrame)
//	{
//		dprint("\n", size, "\n");
//		ui ddd = 0;
//		for (ui oh : pF.adj)
//		{
//			if (size == 56 && ddd == 0)
//			{
//				int lp = 0;
//			}
//#if 1
//			std::vector<OV> loop = searchLoop(/*polymesh, polyFrame, */oh, size);
//			/*if (loop.size() > 4)
//				dprint(ddd++, "loop size:", loop.size());
//			else
//				dprint(ddd++);*/
//
//			polymesh.add_face(loop);
//			loop.push_back(loop.front());
//			std::vector<OV>::const_iterator vend = loop.end() - 1;
//			for (std::vector<OV>::const_iterator v = loop.begin(); v != vend; v++)
//				polyFrame[v->idx()].adj.erase((v + 1)->idx());
//#else
//			std::vector<ui> loop = searchLoop(/*polymesh, polyFrame, */oh, size);
//			if (loop.empty())
//				continue;
//			std::vector<ui>::const_iterator vend = loop.end() - 1;
//			for (std::vector<ui>::const_iterator v = loop.begin(); v != vend; v++)
//			{
//				polyFrame[*v].adj.erase(*(v + 1));
//			}
//			switch (loop.size())
//			{
//			case 4:
//				polymesh.add_face(polyFrame[loop[0]].vert, polyFrame[loop[1]].vert, polyFrame[loop[2]].vert);
//				break;
//			case 5: {
//				auto rrr = polymesh.add_face(polyFrame[loop[0]].vert, polyFrame[loop[1]].vert, polyFrame[loop[2]].vert, polyFrame[loop[3]].vert); 
//				for (auto hgp : polymesh.fv_range(rrr))
//				{
//					//dprint(hgp.idx());
//				}
//				//dprint(rrr.idx());
//			}
//				break;
//			case 6:
//			{
//				//system("pause");
//				loop.pop_back();
//				std::vector<OpenMesh::VertexHandle> vhandle;
//				for (ui l : loop)
//					vhandle.push_back(polyFrame[l].vert);
//				polymesh.add_face(vhandle);
//				break;
//			}
//			default:
//				break;
//			}
//#endif
//		}
//		size++;
//	}
//#endif
//#pragma endregion

	polyFrame = rrr;
}

void PositionField::computeLinkProperty(vvlt &linkProperty)
{
	auto calcDirection = [&](const Vec3d &o0, const Vec3d &n0, const Vec3d &o1, const Vec3d &n1)
	{
		Vec3d _P1[2] = { o0,n0.cross(o0) };
		Vec3d _P2[2] = { o1,n1.cross(o1) };
		double flag_max = 0;
		ui id1 = 0, id2 = 0;
		for (ui i = 0; i < 2; i++)
		{
			for (ui j = 0; j < 2; j++)
			{
				float dottemp = _P1[i].dot(_P2[j]);
				if (std::fabs(dottemp) > std::fabs(flag_max))
				{
					flag_max = dottemp;
					id1 = i;
					id2 = j;
				}
			}
		}
		return std::make_pair(_P1[id1], _P2[id2] * (flag_max > 0 ? 1 : -1));
	};
	auto middlePoint = [&](const Vec3d &v0, const Vec3d &n0, const Vec3d &v1, const Vec3d &n1)
	{
		double n0n1 = n0.dot(n1);
		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
	};
	auto calcProperty = [&](const Vec3d &p0, const Vec3d &v0, const Vec3d &n0, const Vec3d &o0,
		const Vec3d &p1, const Vec3d &v1, const Vec3d &n1, const Vec3d &o1, const double length)
	{
		const Vec3d mP = middlePoint(v0, n0, v1, n1); 
		signed short linkP[4];
		Vec3d y0 = n0.cross(o0);
		linkP[0] = static_cast<signed short>(std::floor(o0.dot(mP - p0) / length));
		linkP[1] = static_cast<signed short>(std::floor(y0.dot(mP - p0) / length));
		Vec3d y1 = n1.cross(o1);
		linkP[2] = static_cast<signed short>(std::floor(o1.dot(mP - p1) / length));
		linkP[3] = static_cast<signed short>(std::floor(y1.dot(mP - p1) / length));
		double dis = DBL_MAX;
		signed short id0 = 4, id1 = 4;
		for (signed short i = 0; i < 4; i++)
		{
			Vec3d position0 = p0 + ((linkP[0] + (i & 1)) * o0 + (linkP[1] + ((i & 2) >> 1)) * y0)*length;
			for (signed short j = 0; j < 4; j++)
			{
				Vec3d position1 = p1 + ((linkP[2] + (j & 1)) * o1 + (linkP[3] + ((j & 2) >> 1)) * y1)*length;
				double d = (position0 - position1).norm();
				if (d < dis)
				{
					dis = d; id0 = i; id1 = j;
				}
			}
		}
		switch (std::abs(linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1)) +
			std::abs(linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1)))
		{
		case 0:
			return SAME;
		case 1:
			return UNIT;
		default:
			return DISCONNECTION;
		}
	};
	const vv3d &p_ = h->p;
	const vv3d &v_ = h->v;
	const vv3d &n_ = h->n;
	const vv3d &o_ = h->o;
	const vvL &adj_ = h->adj;
	
	ui nv = p_.size();
	linkProperty.resize(nv);
	for (ui id = 0; id < nv; id++)
		linkProperty[id].reserve(adj_[id].size());

	for(ui id = 0; id < nv; id++)
	{
		const Vec3d &pid = p_[id];
		const Vec3d &vid = v_[id];
		const Vec3d &nid = n_[id];
		const Vec3d &oid = o_[id];
		for (ui i = 0; i < adj_[id].size(); i++)
		{
			ui index = adj_[id][i].id;
			std::pair<Vec3d, Vec3d> value = calcDirection(oid, nid, o_[index], n_[index]);
			linkProperty[id].push_back(calcProperty(pid, vid, nid, value.first, p_[index], v_[index], n_[index], value.second, length));
		}
	}
}

void PositionField::computePolyFrame(PolyMesh &polymesh, std::vector<Vec3d> &pos, std::vector<VertexValue> &polyFrame)
{
	vvlt linkProperty;
	computeLinkProperty(linkProperty);
	lP = linkProperty;

	const vv3d &v_ = h->v;
	const vv3d &n_ = h->n;
	const vv3d &p_ = h->p;
	const std::vector<std::vector<Link>> &adj_ = h->adj;
	const ui nv = v_.size();

	//构造图，图的顶点是原网格顶点，边是带有UNIT标记的原网格边
	struct graphKnot
	{
		ui id;
		std::set<ui> adj;
		Vec3d weightPosition;
		Vec3d weightNormal;
		double weight;
		ui principleKnot;
		std::set<ui> minorKnot;
		graphKnot(ui i) : id(i), principleKnot(i) { }
		~graphKnot() {}
		inline bool deleted()
		{
			return principleKnot != id;
		}
	};

	std::vector<graphKnot*> graph(nv);
	//计算需要合并的点集C,将对应需要合并的边加入edge
	std::vector<edgeWeight> edge;
	ui size = 0;
	for (const std::vector<LinkType> &lt : linkProperty)
	{
		for (const LinkType &l : lt)
		{
			if (l == SAME)
			{
				size++;
			}
		}
	}
	edge.reserve(size);

	const double cof = -9.0 / length / length;
	for (ui id = 0; id < nv; id++)
	{
		ui lpSize = linkProperty[id].size();
		const Vec3d &pid = p_[id];
		const Vec3d &nid = n_[id];
		graphKnot* &gid = graph[id];
		gid = new graphKnot(id);
		gid->weight = exp(cof * (pid - v_[id]).sqrnorm());
		gid->weightPosition = gid->weight*pid;
		gid->weightNormal = gid->weight*nid;
		for (ui i = 0; i < lpSize; i++)
		{
			if (linkProperty[id][i] == SAME)
			{
				ui index = adj_[id][i].id;
				edge.emplace_back(id, index, (pid - p_[index]).norm());
			}
			else if (linkProperty[id][i] == UNIT)
			{
				gid->adj.insert(adj_[id][i].id);
			}
		}
	}

	std::sort(edge.begin(), edge.end(), std::less<edgeWeight>());
	size = edge.size() / 2;
	for (ui id = 1; id < size; id++)
	{
		std::swap(edge[id], edge[id << 1]);
	}
	edge.erase(edge.begin() + size, edge.end());

	//访问C,并将相应的顶点（在同一个block中）合并
	for (const edgeWeight &e : edge)
	{
		graphKnot* &gi = graph[graph[e.i]->principleKnot];
		graphKnot* &gj = graph[graph[e.j]->principleKnot];
		if (gi == gj)
			continue;
		gj->principleKnot = gi->id;
		for (ui gmK : gj->minorKnot)
		{
			graph[gmK]->principleKnot = gi->id;
		}
		gi->minorKnot.insert(gj->id);
		gi->minorKnot.insert(gj->minorKnot.begin(), gj->minorKnot.end());
		gi->weight += gj->weight;
		gi->weightPosition += gj->weightPosition;
		gi->weightNormal += gj->weightNormal;
		gi->adj.insert(gj->adj.begin(), gj->adj.end());
	}
	//至此合并完毕，其中principleKnot指向自身的graphKnot为block的代表元，这些graphKnot将成为四边形网格的顶点，其与adj将构成边

	//将上述几何拓扑关系整合
	{
		ui leftGKSize = 0;
		for (graphKnot* &gK : graph)
			if (!(gK->deleted()))
				leftGKSize++;
		ui leastBlockSize = std::max((ui)1, nv / (leftGKSize * 10));

		std::vector<ui> INDEX2index(nv);
		std::vector<ui> index2INDEX;
		index2INDEX.reserve(nv);
		for (graphKnot* &gK : graph)
		{
			if (!(gK->deleted()))
			{
				std::set<ui> newadj;
				for (ui ga : gK->adj)
				{
					newadj.insert(graph[ga]->principleKnot);
				}
				gK->adj = newadj;
				INDEX2index[gK->id] = index2INDEX.size();
				index2INDEX.push_back(gK->id);
			}
		}

		//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
		//std::vector<VertexValue> poly;
		ui polysize = index2INDEX.size();
		polyFrame.reserve(polysize);
		for (ui i = 0; i < polysize; i++)
		{
			ui id = index2INDEX[i];
			std::set<ui> &gadj = graph[id]->adj;
			polyFrame.emplace_back(polymesh.add_vertex(graph[id]->weightPosition / graph[id]->weight), graph[id]->weightNormal.normalize(), std::set<ui>());
			//polyFrame.emplace_back(OpenMesh::VertexHandle(), graph[id]->weightNormal.normalize(), std::set<ui>());
			pos.push_back(graph[id]->weightPosition / graph[id]->weight);
			std::set<ui> &padj = polyFrame.back().adj;
			for (std::set<ui>::const_iterator gitr = gadj.begin(); gitr != gadj.end(); gitr++)
			{
				padj.insert(INDEX2index[*gitr]);
			}
			padj.erase(i);
		}
	}

	for (graphKnot* &gK : graph)
	{
		delete gK;
		gK = nullptr;
	}

#pragma region post-processing for position field to remedy artifacts
	//处理过于靠近某条边的顶点，此时认为该点在这条边上，删除这条边
	size = 0;
	double minThreshold = 0.3*length;

	for (VertexValue &pFi : polyFrame)
	{
		Vec3d &posi = polymesh.point(pFi.vert);
		for (ui pFj : pFi.adj)
		{
			Vec3d &posj = polymesh.point(polymesh.vertex_handle(pFj));
			double dij = (posi - posj).norm();
			for (ui pFk : polyFrame[pFj].adj)
			{
				if (size == pFk)
					continue;
				Vec3d &posk = polymesh.point(polymesh.vertex_handle(pFk));
				double djk = (posj - posk).norm();
				double dki = (posk - posi).norm();
				if (djk < std::max(dij, dki))
					continue;
				double dcir = 0.5*(dij + djk + dki);
				if (sqrt(dcir*(dcir - dij)*(dcir - djk)*(dcir - dki)) * 2.0 < minThreshold * djk)
				{
					if (!pFi.adj.count(pFk))
						continue;
					polymesh.set_point(pFi.vert, 0.5*(posj + posk));
					polyFrame[pFj].adj.erase(pFk);
					polyFrame[pFk].adj.erase(pFj);
					goto goto20211009;
				}
			}
		}
	goto20211009:;
		++size;
	}


	//删除valence<3的顶点
	size = polyFrame.size() - 1;
	for (std::vector<VertexValue>::reverse_iterator vertex = polyFrame.rbegin(); vertex != polyFrame.rend(); ++vertex, --size)
	{
		if (vertex->adj.size() > 2)
			continue;
		ui pFsize = polyFrame.size() - 1;
		for (ui var : vertex->adj)
		{
			polyFrame[var].adj.erase(size);
		}
		vertex->normal = polyFrame.back().normal;
		vertex->adj = polyFrame.back().adj;
		for (ui var : vertex->adj)
		{
			polyFrame[var].adj.erase(pFsize);
			polyFrame[var].adj.insert(size);
		}
		polyFrame.pop_back();

		polymesh.set_point(polymesh.vertex_handle(size), polymesh.point(polymesh.vertex_handle(pFsize)));
		polymesh.delete_vertex(polymesh.vertex_handle(polymesh.n_vertices() - 1));
	}
	polymesh.garbage_collection();
#pragma endregion



}

#pragma endregion

#else
//#pragma region MeshField
//MeshField::hierarchy::hierarchy(const TriMesh &mesh, ui depth)
//{
//	//以 mesh 初始化该层
//	init(mesh);
//	if (depth < 2)
//	{
//		hierarchy_depth = 1;
//		return;
//	}
//
//	//循环初始化剩下各层
//	finer_hierarchy = nullptr;
//	coarser_hierarchy = new hierarchy(this);
//	hierarchy* h = coarser_hierarchy;
//	h->finer_hierarchy = this;
//	hierarchy_depth = depth;
//	for (; depth > 2 && h->v.size() > 1; depth--, h = h->coarser_hierarchy)
//	{
//		h->coarser_hierarchy = new hierarchy(h);
//		h->coarser_hierarchy->finer_hierarchy = h;
//	}
//	h->coarser_hierarchy = nullptr;
//	hierarchy_depth -= depth - 2;
//
//	//初始化各层 hierarchy_depth
//	depth = hierarchy_depth;
//	h = this->coarser_hierarchy;
//	while (h)
//	{
//		h->hierarchy_depth = --depth;
//		h = h->coarser_hierarchy;
//	}
//}
//
//double MeshField::hierarchy::calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv)
//{
//	double area = 0;
//	v3d pos = OpenMesh2EigenVector(mesh.point(tv));
//	std::vector<Link> &ad = adj[tv.idx()];
//	ad.reserve(mesh.valence(tv));
//	for (auto tvoh : mesh.voh_range(tv))
//	{
//		ad.push_back(mesh.to_vertex_handle(tvoh).idx());
//		if (!mesh.face_handle(tvoh).is_valid())
//			continue;
//		v3d s = OpenMesh2EigenVector(mesh.point(mesh.to_vertex_handle(tvoh)));
//		v3d t = OpenMesh2EigenVector(mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(tvoh))));
//		v3d c = (pos + s + t) / 3.0;
//		area += 0.25*((s - pos).cross(c - pos).norm() + (t - pos).cross(c - pos).norm());
//	}
//	return area;
//}
//
//void MeshField::hierarchy::init(const TriMesh &mesh)
//{
//	ui nv = mesh.n_vertices();
//	v.resize(3, nv);
//	n.resize(3, nv);
//	a.resize(nv);
//	adj.resize(nv);
////	ui id = 0;
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////	{
////		OV tv = mesh.vertex_handle(id);
////#else
////	for (OV tv : mesh.vertices())
////	{
////#endif // RUN_MESHFIELD_PARALLEL
////		v.col(id) = OpenMesh2EigenVector(mesh.point(tv));
////		n.col(id) = OpenMesh2EigenVector(mesh.calc_normal(tv));
////		a(id) = calcVoronoiAndAdjacency(mesh, tv);
////#ifdef RUN_MESHFIELD_PARALLEL
////	});
////#else
////		++id;
////    }
////#endif
//	SerialParallelBlock(
//		OV tv = mesh.vertex_handle(id);
//		v.col(id) = OpenMesh2EigenVector(mesh.point(tv));
//		n.col(id) = OpenMesh2EigenVector(mesh.calc_normal(tv));
//		a(id) = calcVoronoiAndAdjacency(mesh, tv);
//		)
//	graph_color();
//}
//
//void MeshField::hierarchy::graph_color()
//{
//	ui nsizes = adj.size();
//	std::vector<int> color(nsizes, -1);
//	std::vector<ui> possible_color;
//
//	ui ncolors = 0;
//	for (size_t i = 0; i < nsizes; i++)
//	{
//		std::fill(possible_color.begin(), possible_color.end(), 1);
//		for each (auto link in adj[i])
//		{
//			int c = color[link.id];
//			if (c >= 0)
//			{
//				possible_color[c] = 0;
//			}
//		}
//		int color_id = -1;
//		for (auto j = 0; j < possible_color.size(); j++)
//		{
//			if (possible_color[j] != 0)
//			{
//				color_id = j;
//				break;
//			}
//		}
//		if (color_id < 0)
//		{
//			color_id = ncolors++;
//			possible_color.resize(ncolors);
//		}
//		color[i] = color_id;
//	}
//
//	phase.resize(ncolors);
//	for (size_t i = 0; i < nsizes; i++)
//	{
//		phase[color[i]].push_back(i);
//	}
//}
//
//#include <deque>
//MeshField::hierarchy::hierarchy(hierarchy *fh)
//{
//	const m3xd &v_ = fh->v;
//	const m3xd &n_ = fh->n;
//	const vxd &a_ = fh->a;
//	const vvL &adj_ = fh->adj;
//	ui nv = v_.cols();
//
//
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::concurrent_vector<edgeWeight> edge;
//#else
//	std::vector<edgeWeight> edge;
//#endif
//	ui size = 0;
//	for (const std::vector<Link> &l : adj_)
//		size += l.size();
//	edge.reserve(size);
//
//	/*size = 0;
//	for (const std::vector<Link> &v0 : adj_)
//	{
//		for (Link v1 : v0)
//			edge.emplace_back(size, v1.id, n_.col(size).dot(n_.col(v1.id))*std::max(a_(size) / a_(v1.id), a_(v1.id) / a_(size)));
//		size++;
//	}*/
//
//	//size = adj_.size();
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for(ui id = 0; id < nv; id ++)
////#endif
////	{
////		const vL &v0 = adj_[id];
////		for (const Link &v1 : v0)
////		{
////			const ui &vid = v1.id;
////			edge.emplace_back(id, vid, n_.col(id).dot(n_.col(vid))*std::max(a_(id) / a_(vid), a_(vid) / a_(id)));
////		}
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//
//	SerialParallelBlock(
//		const vL &v0 = adj_[id];
//	    for (const Link &v1 : v0)
//	    {
//		    const ui &vid = v1.id;
//		    edge.emplace_back(id, vid, n_.col(id).dot(n_.col(vid))*std::max(a_(id) / a_(vid), a_(vid) / a_(id)));
//        }
//	)
//
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::parallel_sort(edge.begin(), edge.end(), std::greater<edgeWeight>());
//#else
//	std::sort(edge.begin(), edge.end(), std::greater<edgeWeight>());
//#endif
//	std::deque<bool> mergeFlag(nv, false);
//	size = 0;
//	for (const edgeWeight &e : edge)
//	{
//		if (mergeFlag[e.i] || mergeFlag[e.j])
//			continue;
//		mergeFlag[e.i] = true;
//		mergeFlag[e.j] = true;
//		edge[size++] = e;
//	}
//
//	ui leftVertices = nv - size;
//
//	v.resize(3, leftVertices);
//	n.resize(3, leftVertices);
//	a.resize(leftVertices);
//	adj.resize(leftVertices);
//
//	tofiner.resize(2, leftVertices);
//	vxu &tc = fh->tocoarser;
//	tc.resize(nv);
//
//	nv = size;
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; id++)
////#endif
////	{
////		edgeWeight &e = edge[id];
////		double a0 = a_(e.i);
////		double a1 = a_(e.j);
////		double area = a0 + a1;
////		if (area > epsilon)
////		{
////			v.col(id) = (v_.col(e.i)*a0 + v_.col(e.j)*a1) / area;
////		}
////		else
////		{
////			v.col(id) = (v_.col(e.i) + v_.col(e.j))*0.5;
////		}
////
////		n.col(id) = (n_.col(e.i)*a0 + n_.col(e.j)*a1).normalized();
////		a(id) = area;
////		tofiner.col(id) << e.i, e.j;
////		tc(e.i) = id; tc(e.j) = id;
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//
//	SerialParallelBlock(
//		edgeWeight &e = edge[id];
//	    double a0 = a_(e.i);
//	    double a1 = a_(e.j);
//	    double area = a0 + a1;
//	    if (area > epsilon)
//	    {
//		    v.col(id) = (v_.col(e.i)*a0 + v_.col(e.j)*a1) / area;
//   	    }
//	    else
//	    {
//		    v.col(id) = (v_.col(e.i) + v_.col(e.j))*0.5;
//	    }
//
//	    n.col(id) = (n_.col(e.i)*a0 + n_.col(e.j)*a1).normalized();
//	    a(id) = area;
//		tofiner(0, id) = e.i; tofiner(1, id) = e.j;
//	    tc(e.i) = id; tc(e.j) = id;
//	)
//
//	nv = v_.cols();
//	for (ui id = 0; id < nv; id++)
//	{
//		if (mergeFlag[id])
//			continue;
//		
//		v.col(size) = v_.col(id);
//		n.col(size) = n_.col(id);
//		a(size) = a_(id);
//		tofiner.col(size) << id, nv;
//		tc(id) = size;
//		size++;
//	}
//	assert(size == leftVertices);
//
//	std::swap(nv, size);
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; ++id)
////#endif
////	{
////		std::vector<Link> &oneRing = adj[id];
////		//std::vector<Link> oneRing;
////		//导入精细层网格的one-ring
////		ui i = tofiner(0, id);
////		for (const Link &var : adj_[i])
////		{
////			if (tc(var.id) == id) continue;
////			oneRing.emplace_back(tc(var.id), var.weight);
////		}
////		i = tofiner(1, id);
////		if (i != size)
////		{
////			for (const Link &var : adj_[i])
////			{
////				if (tc(var.id) == id) continue;
////				oneRing.emplace_back(tc(var.id), var.weight);
////			}
////		}
////		if (oneRing.size() > 1)
////		{
////			//删除oneRing重复元素
////			std::sort(oneRing.begin(), oneRing.end());
////			std::vector<Link>::iterator itr = oneRing.begin();
////			++itr;
////			for (; itr != oneRing.end(); ++itr)
////			{
////				std::vector<Link>::iterator newItr = itr;
////				--newItr;
////				if (itr->id == newItr->id)
////				{
////					itr->weight += newItr->weight;
////					itr = oneRing.erase(newItr);
////				}
////			}
////		}
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//	SerialParallelBlock(
//		std::vector<Link> &oneRing = adj[id];
//	    //std::vector<Link> oneRing;
//	    //导入精细层网格的one-ring
//	    ui i = tofiner(0, id);
//	    for (const Link &var : adj_[i])
//	    {
//		    if (tc(var.id) == id) continue;
//		    oneRing.emplace_back(tc(var.id), var.weight);
//	    }
//	    i = tofiner(1, id);
//	    if (i != size)
//	    {
//		    for (const Link &var : adj_[i])
//		    {
//			    if (tc(var.id) == id) continue;
//			    oneRing.emplace_back(tc(var.id), var.weight);
//		    }
//	    }
//	    if (oneRing.size() > 1)
//	    {
//		    //删除oneRing重复元素
//		    std::sort(oneRing.begin(), oneRing.end());
//		    std::vector<Link>::iterator itr = oneRing.begin();
//		    ++itr;
//		    for (; itr != oneRing.end(); ++itr)
//		    {
//			    std::vector<Link>::iterator newItr = itr;
//			    --newItr;
//			    if (itr->id == newItr->id)
//			    {
//				    itr->weight += newItr->weight;
//				    itr = oneRing.erase(newItr);
//			    }
//		    }
//	    }
//	)
//	graph_color();
//}
//
//MeshField::hierarchy::~hierarchy()
//{
//	finer_hierarchy = nullptr;
//	if (hierarchy_depth >> 1)
//	{
//		delete coarser_hierarchy;
//		coarser_hierarchy = nullptr;
//	}
//}
//#pragma endregion
//
//
//#pragma region OrientationField
//void OrientationField::loopDivision(TriMesh &mesh, ui divideTimes)
//{
//	if (divideTimes < 1)
//		return;
//	OpenMesh::Subdivider::Uniform::LoopT<TriMesh> loopDivider;
//	loopDivider.attach(mesh);
//	loopDivider(divideTimes);
//	loopDivider.detach();
//	mesh.update_normals();
//	OpenMesh::IO::write_mesh(mesh, "human2.obj");
//}
//
//void OrientationField::randomInit(TriMesh &mesh)
//{
//	loopDivision(mesh, 0);
//
//	h = new hierarchy(mesh, DefaultDepth);
//	hierarchy *hw = h;
//
//	while (hw->coarser_hierarchy)
//	{
//		hw->o.resize(3, hw->v.cols());
//		hw = hw->coarser_hierarchy;
//	}
//	m3xd &n_ = hw->n;
//	m3xd &o_ = hw->o;
//	ui cols = n_.cols();
//	o_.resize(3, cols);
//
//	for (ui i = 0; i < cols; i++)
//	{
//		o_.col(i) = hw->calcXAxis(n_.col(i));
//	}
//	dprint("Hierarchy Build Done!\n\n");
//}
//
//void OrientationField::GSOptimize()
//{
//	auto calcDirection = [&](const v3d &o0, const v3d &n0, const v3d &o1, const v3d &n1, v3d *direction)
//	{
//		v3d _P1[2] = { o0,n0.cross(o0) };
//		v3d _P2[2] = { o1,n1.cross(o1) };
//		double flag_max = 0;
//		ui id1 = 0, id2 = 0;
//		for (ui i = 0; i < 2; i++)
//		{
//			for (ui j = 0; j < 2; j++)
//			{
//				double dottemp = _P1[i].dot(_P2[j]);
//				if (std::fabs(dottemp) > std::fabs(flag_max))
//				{
//					flag_max = dottemp;
//					id1 = i;
//					id2 = j;
//				}
//			}
//		}
//		direction[0] = _P1[id1];
//		direction[1] = _P2[id2] * (flag_max > 0 ? 1 : -1);
//	};
//
//	hierarchy* hn = h;
//	while (hn->coarser_hierarchy)
//		hn = hn->coarser_hierarchy;
//
//	const ui iterationTimes = 5;
//	while (true)
//	{
//		//在当前层 Gauss-Seidel 迭代，得到优化 Orientation Field
//		m3xd &o_ = hn->o;
//		ui nv = o_.cols();
//		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
//		const m3xd &n_ = hn->n;
//		const vvL &adj_ = hn->adj;
//
//		for (int i = 0; i < iterationTimes; i++)
//		{
//#ifdef RUN_MESHFIELD_PARALLEL
//			for (const std::vector<ui> &ph : hn->phase)
//			{
//				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
//				tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
//				{
//					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
//					{
//						ui id = ph[phaseIdx];
//#else
//			for (ui id = 0; id < nv; id++)
//			{
//#endif
//				double sum = 0;
//				v3d oid = o_.col(id);
//				const v3d &nid = n_.col(id);
//				for (const Link &oneRing : adj_[id])
//				{
//					double w = oneRing.weight;
//					v3d direction[2];
//					calcDirection(oid, nid, o_.col(oneRing.id), n_.col(oneRing.id), direction);
//					oid = direction[0]*sum + direction[1]*w;
//					sum += w;
//					oid -= nid.dot(oid)*nid;
//					oid.normalize();
//				}
//				o_.col(id) = oid;
//#ifdef RUN_MESHFIELD_PARALLEL
//			}
//				});
//#endif
//			}
//	//dprint("iteration times:", i);
//		}
//    //将当前层 Orientation Field 传递到精细层
//    if (!(hn->finer_hierarchy))
//       break;
//
//	const m2xu &tofiner = hn->tofiner;
//    m3xd &o_finer = hn->finer_hierarchy->o;
//    const m3xd &n_finer = hn->finer_hierarchy->n;
//    const ui finer_size = n_finer.cols();
//
//#ifdef RUN_MESHFIELD_PARALLEL
//    tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//#else
//    for (ui id = 0; id < nv; id++)
//#endif
//    {
//		ui finer_id = tofiner(0, id);
//		o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();
//		finer_id = tofiner(1, id);
//		if (finer_id != finer_size)
//			o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();
//    }
//#ifdef RUN_MESHFIELD_PARALLEL
//    );
//#endif
//
//     hn = hn->finer_hierarchy;
//	}
//	dprint("\nOrientation Field Build Done!\n\n");
//}
//#pragma endregion
//
//
//#pragma region PositionField
//void PositionField::randomInit(TriMesh &mesh)
//{
//	clock_t start = clock();
//	if (!h)
//	{
//		dprint("Orientation Field need to be build first !!!");
//		exit(0);
//	}
//	hierarchy* hw = h;
//	while (hw->coarser_hierarchy)
//	{
//		hw->p.resize(3, hw->v.cols());
//		hw = hw->coarser_hierarchy;
//	}
//	hw->p = hw->v;
//
//	length = 0;
//	for (OE te : mesh.edges())
//		length += mesh.calc_edge_length(te);
//	length /= mesh.n_edges() / scale;
//	dprint("预计算", clock() - start);
//}
//
//void PositionField::GSOptimize()
//{
//	auto middlePoint = [&](const v3d &v0, const v3d &n0, const v3d &v1, const v3d &n1)
//	{
//		double n0n1 = n0.dot(n1);
//		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
//			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
//	};
//	double inv_length = 1.0 / length;
//	auto calcPosition = [&](const v3d &p0, const v3d &v0, const v3d &n0, const v3d &o0,
//		const v3d &p1, const v3d &v1, const v3d &n1, const v3d &o1, v3d *position)
//	{
//		const v3d mP = middlePoint(v0, n0, v1, n1);
//		v3d y0 = n0.cross(o0);
//		v3d origin0 = p0 + length * (o0*std::floor(o0.dot(mP - p0) * inv_length) + y0 * std::floor(y0.dot(mP - p0) * inv_length));
//		v3d y1 = n1.cross(o1);
//		v3d origin1 = p1 + length * (o1*std::floor(o1.dot(mP - p1) * inv_length) + y1 * std::floor(y1.dot(mP - p1) * inv_length));
//
//		v3d updateP0[4], updateP1[4];
//		for (ui i = 0; i < 4; i++)
//		{
//			updateP0[i] = origin0 + ((i & 1)*o0 + ((i & 2) >> 1)*y0)*length;
//			updateP1[i] = origin1 + ((i & 1)*o1 + ((i & 2) >> 1)*y1)*length;
//		}
//		double dis = DBL_MAX;
//		ui index[2];
//		for (ui i = 0; i < 4; i++)
//		{
//			for (ui j = 0; j < 4; j++)
//			{
//				double d = (updateP0[i] - updateP1[j]).squaredNorm();
//				if (d < dis)
//				{
//					dis = d;
//					index[0] = i; index[1] = j;
//				}
//			}
//		}
//		position[0] = updateP0[index[0]]; position[1] = updateP1[index[1]];
//	};
//	auto positionRound = [&](v3d &p, const v3d &v, const v3d &n, const v3d &o)
//	{
//		const v3d y = n.cross(o);
//		p += length * (std::floor((v - p).dot(o) * inv_length + 0.5)*o + std::floor((v - p).dot(y) * inv_length + 0.5)*y);
//	};
//
//	hierarchy* hn = h;
//	while (hn->coarser_hierarchy)
//		hn = hn->coarser_hierarchy;
//
//	const ui iterTimes = 10;
//	while (true)
//	{
//		//在当前层 Gauss-Seidel 迭代，得到优化 Position Field
//		m3xd &p_ = hn->p;
//		ui nv = p_.cols();
//		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
//		clock_t start = clock();
//		const m3xd &v_ = hn->v;
//		const m3xd &n_ = hn->n;
//		const m3xd &o_ = hn->o;
//		const vvL &adj_ = hn->adj;
//		for (int i = 0; i < /*hn->hierarchy_depth*/iterTimes; i++)
//		{
//#ifdef RUN_MESHFIELD_PARALLEL
//			for (const std::vector<ui> &ph : hn->phase)
//			{
//				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
//				tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
//				{
//					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
//					{
//						ui id = ph[phaseIdx];
//#else
//			        for (ui id = 0; id < nv; id++)
//			        {
//#endif
//						double sum = 0;
//						v3d pid = p_.col(id);
//						const v3d &vid = v_.col(id);
//						const v3d &nid = n_.col(id);
//						const v3d &oid = o_.col(id);
//						for (const Link &oneRing : adj_[id])
//						{
//							double w = oneRing.weight;
//							ui index = oneRing.id;
//							v3d position[2];
//							calcPosition(pid, vid, nid, oid, p_.col(index), v_.col(index), n_.col(index), o_.col(index), position);
//							pid = position[0] * sum + position[1] * w;
//							sum += w;
//							pid /= sum;
//							pid -= nid.dot(pid - vid)*nid;
//						}
//						positionRound(pid, vid, nid, oid);
//						p_.col(id) = pid;
//#ifdef RUN_MESHFIELD_PARALLEL
//			        }
//				});
//#endif
//			}
//		}
//
//		dprint("time0:", clock() - start);
//		start = clock();
//
//			//将当前层 Position Field 传递到精细层
//	    if (!(hn->finer_hierarchy))
//		   break;
//
//	    const m2xu &tofiner = hn->tofiner;
//	    const m3xd &v_finer = hn->finer_hierarchy->v;
//	    const m3xd &n_finer = hn->finer_hierarchy->n;
//	    m3xd &p_finer = hn->finer_hierarchy->p;
//	    const ui finer_size = n_finer.cols();
//
//#ifdef RUN_MESHFIELD_PARALLEL
//	    tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//#else
//	    for (ui id = 0; id < nv; id++)
//#endif
//	    {
//		    ui finer_id = tofiner(0, id);
//		    p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
//		    finer_id = tofiner(1, id);
//		    if (finer_id != finer_size)
//	 		    p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
//	    }
//#ifdef RUN_MESHFIELD_PARALLEL
//	    );
//#endif
//    	hn = hn->finer_hierarchy; dprint("time1:", clock() - start);
//	}
//	dprint("\nPosition Field Build Done!\n\n");
//}
//
//#include <set>
//void PositionField::extractMesh(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame)
//{
//	computePolyFrame(polymesh, pos, polyFrame);
//	ui nv = polymesh.n_vertices();
//	VertexValue *rrr = new VertexValue[nv];
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; id++)
////#endif
////	{
////		rrr[id] = polyFrame[id];
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//	SerialParallelBlock(rrr[id] = polyFrame[id];)
//
//	clock_t start = clock();
//	//#ifdef RUN_MESHFIELD_PARALLEL
//	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//	//#else
//	//	for(ui id = 0; id < nv; id++)
//	//#endif
//	//	{
//	//		VertexValue &pF = polyFrame[id];
//	//		std::list<ui> &orientedAdj = pF.orientedAdj;
//	//		for (const ui &va : pF.adj)
//	//			orientedAdj.push_back(va);
//	//		v3d x = h->calcXAxis(pF.normal);
//	//		v3d y = pF.normal.cross(x);
//	//		orientedAdj.sort([&](const ui &id0, const ui &id1) {
//	//			v3d dir0 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id0)) - polymesh.point(pF.vert));
//	//			v3d dir1 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id1)) - polymesh.point(pF.vert));
//	//			return std::atan2(y.dot(dir0), x.dot(dir0)) > std::atan2(y.dot(dir1), x.dot(dir1));
//	//		});
//	//	}
//	//#ifdef RUN_MESHFIELD_PARALLEL
//	//	);
//	//#endif
//#if 1
//	SerialParallelBlock(
//		VertexValue &pF = polyFrame[id];
//	    std::list<ui> &orientedAdj = pF.orientedAdj;
//		for (const ui &va : pF.adj) {
//			orientedAdj.push_back(va);
//		}
//	    v3d x = h->calcXAxis(pF.normal);
//	    v3d y = pF.normal.cross(x);
//	    orientedAdj.sort([&](const ui &id0, const ui &id1) {
//		    v3d dir0 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id0)) - polymesh.point(pF.vert));
//		    v3d dir1 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id1)) - polymesh.point(pF.vert));
//		    return std::atan2(y.dot(dir0), x.dot(dir0)) > std::atan2(y.dot(dir1), x.dot(dir1));
//	    });
//		)
//	dprint("oriente edges:", clock() - start, "ms");
//	start = clock();
//	//auto ppp = polyFrame[974];
//
//	for (ui id = 0; id < nv; id++)
//	{
//		dprint(id);
//		//if (id == 669)//human.obj
//		if(id==426)//vase.obj
//		{
//#if 1
//			int p = 0;
//#else
//			break;
//#endif
//		}
//		VertexValue &pF = polyFrame[id];
//		for (const ui &oh : pF.adj)
//		{
//			std::vector<OV> loop;
//			loop.reserve(4);
//			loop.push_back(polymesh.vertex_handle(oh));
//			ui previousId = id;
//			do
//			{
//				std::list<ui> &currentOrientedAdj = polyFrame[loop.back().idx()].orientedAdj;
//				std::list<ui>::iterator pointer = currentOrientedAdj.begin();
//				while (*pointer != previousId) {
//					++pointer;
//				}
//				previousId = loop.back().idx();
//				pointer = ++pointer == currentOrientedAdj.end() ? currentOrientedAdj.begin() : pointer;
//				polyFrame[previousId].adj.erase(*pointer);
//				loop.push_back(polymesh.vertex_handle(*pointer));
//				if (*pointer == id)
//				{
//					break;
//				}
//			} while (true);
//			polymesh.add_face(loop);
//		}
//		pF.adj.clear();
//	}
//	for (auto tf : polymesh.faces())
//	{
//		std::cout << "\n" << tf.idx() << ":   ";
//		for (auto tfv : polymesh.fv_range(tf))
//		{
//			std::cout << tfv.idx() << ",";
//		}
//	}
//	dprint("extract mesh:", clock() - start, "ms");
//#endif
//
//	delete[] polyFrame;
//	polyFrame = rrr;
//}
//
//void PositionField::computeLinkProperty(vvlt &linkProperty)
//{
//	double inv_length = 1.0 / length;
//	auto middlePoint = [&](const v3d &v0, const v3d &n0, const v3d &v1, const v3d &n1)
//	{
//		double n0n1 = n0.dot(n1);
//		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
//			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
//	};
//	auto calcDirection = [&](const v3d &o0, const v3d &n0, const v3d &o1, const v3d &n1, v3d *direction)
//	{
//		v3d _P1[2] = { o0,n0.cross(o0) };
//		v3d _P2[2] = { o1,n1.cross(o1) };
//		double flag_max = 0;
//		ui id1 = 0, id2 = 0;
//		for (ui i = 0; i < 2; i++)
//		{
//			for (ui j = 0; j < 2; j++)
//			{
//				double dottemp = _P1[i].dot(_P2[j]);
//				if (std::fabs(dottemp) > std::fabs(flag_max))
//				{
//					flag_max = dottemp;
//					id1 = i;
//					id2 = j;
//				}
//			}
//		}
//		direction[0] = _P1[id1]; direction[1] = _P2[id2] * (flag_max > 0 ? 1 : -1);
//	};
//	auto calcProperty = [&](const v3d &p0, const v3d &v0, const v3d &n0, const v3d &o0,
//		const v3d &p1, const v3d &v1, const v3d &n1, const v3d &o1)
//	{
//		const v3d mP = middlePoint(v0, n0, v1, n1);
//		signed short linkP[4];
//		v3d y0 = n0.cross(o0);
//		linkP[0] = static_cast<signed short>(std::floor(o0.dot(mP - p0) * inv_length));
//		linkP[1] = static_cast<signed short>(std::floor(y0.dot(mP - p0) * inv_length));
//		v3d y1 = n1.cross(o1);
//		linkP[2] = static_cast<signed short>(std::floor(o1.dot(mP - p1) * inv_length));
//		linkP[3] = static_cast<signed short>(std::floor(y1.dot(mP - p1) * inv_length));
//		v3d position0[4], position1[4];
//		for (signed short i = 0; i < 4; i++)
//		{
//			position0[i] = p0 + ((linkP[0] + (i & 1)) * o0 + (linkP[1] + ((i & 2) >> 1)) * y0)*length;
//			position1[i] = p1 + ((linkP[2] + (i & 1)) * o1 + (linkP[3] + ((i & 2) >> 1)) * y1)*length;
//		}
//		double dis = DBL_MAX;
//		signed short id0 = 4, id1 = 4;
//		for (signed short i = 0; i < 4; i++)
//		{
//			for (signed short j = 0; j < 4; j++)
//			{
//				double d = (position0[i] - position1[j]).squaredNorm();
//				if (d < dis)
//				{
//					dis = d; id0 = i; id1 = j;
//				}
//			}
//		}
//		switch (std::abs(linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1)) +
//			std::abs(linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1)))
//		{
//		case 0:
//			return SAME;
//		case 1:
//			return UNIT;
//		default:
//			return DISCONNECTION;
//		}
//	};
//	const m3xd &p_ = h->p;
//	const m3xd &v_ = h->v;
//	const m3xd &n_ = h->n;
//	const m3xd &o_ = h->o;
//	const vvL &adj_ = h->adj;
//
//	ui nv = p_.cols();
//	linkProperty.resize(nv);
//
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; id++)
////#endif
////	{
////		//linkProperty[id].reserve(adj_[id].size());
////		linkProperty[id].resize(adj_[id].size(),DISCONNECTION);
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//	SerialParallelBlock(linkProperty[id].resize(adj_[id].size(), DISCONNECTION);)
//
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; ++id)
////#endif
////	{
////		const v3d &pid = p_.col(id);
////		const v3d &vid = v_.col(id);
////		const v3d &nid = n_.col(id);
////		const v3d &oid = o_.col(id);
////		ui adjSize = adj_[id].size();
////		for (ui i = 0; i < adjSize; i++)
////		{
////			ui index = adj_[id][i].id;
////			//if (id > index)
////			//{
////			//	//linkProperty[id].push_back(UNDEFINED);
////			//	/*for (const Link &var : adj_[id])
////			//	{
////			//		if (var.id == id)
////			//		{
////			//			linkProperty[id].push_back(linkProperty[index][])
////			//		}
////			//	}*/
////			//	for (ui j = 0; j < adj_[index].size(); ++j)
////			//	{
////			//		if (adj_[index][j].id == id)
////			//		{
////			//			linkProperty[id].push_back(linkProperty[index][j]);
////			//			break;
////			//		}
////			//	}
////			//}
////			//else
////			//{
////			//	v3d direction[2];
////			//	calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
////			//	linkProperty[id].push_back(calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]));
////			//	//linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
////			//}
////
////			if (id > index)
////				continue;
////			v3d direction[2];
////			calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
////			//linkProperty[id].push_back(calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]));
////			linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
////			/*for (ui j = 0; j < adj_[index].size(); ++j)
////			{
////				if (adj_[index][j].id == id)
////				{
////					linkProperty[index][j] = linkProperty[id][i];
////					break;
////				}
////			}*/
////		}
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//	SerialParallelBlock(
//		const v3d &pid = p_.col(id);
//	    const v3d &vid = v_.col(id);
//	    const v3d &nid = n_.col(id);
//	    const v3d &oid = o_.col(id);
//	    ui adjSize = adj_[id].size();
//	    for (ui i = 0; i < adjSize; i++)
//	    {
//		    const ui &index = adj_[id][i].id;
//
//		    if (id > index)
//			    continue;
//		    v3d direction[2];
//		    calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
//		    linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
//	    }
//	)
//}
//
//void PositionField::computePolyFrame(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame)
//{		
//	vvlt linkProperty;
//	clock_t start = clock();
//	computeLinkProperty(linkProperty);
//	dprint("compute link property:", clock() - start, "ms");
//	//lP = linkProperty;
//	start = clock();
//	const m3xd &v_ = h->v;
//	const m3xd &n_ = h->n;
//	const m3xd &p_ = h->p;
//	const std::vector<std::vector<Link>> &adj_ = h->adj;
//	ui nv = v_.cols();
//
//	//构造图，图的顶点是原网格顶点，边是带有UNIT标记的原网格边
//	struct graphKnot
//	{
//		ui id;
//		v3d weightPosition;
//		v3d weightNormal;
//		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//		std::set<ui> adj;
//		double weight;
//		ui principleKnot;
//		std::set<ui> minorKnot;
//		graphKnot(){}
//		graphKnot(ui i) : id(i), principleKnot(i) { }
//		~graphKnot() {}
//		inline bool deleted()
//		{
//			return principleKnot != id;
//		}
//	};
//
//	std::vector<graphKnot> graph(nv);
//	//计算需要合并的点集C,将对应需要合并的边加入edge
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::concurrent_vector<edgeWeight> edge;
//#else
//	std::vector<edgeWeight> edge;
//#endif
//	ui size = 0;
//	for (const std::vector<LinkType> &lt : linkProperty)
//		for (const LinkType &l : lt)
//			if (l == SAME)
//				size++;
//	edge.reserve(size);
//
//	const double cof = -9.0 / (length * length);
////#ifdef RUN_MESHFIELD_PARALLEL
////	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
////#else
////	for (ui id = 0; id < nv; id++)
////#endif
////	{
////		ui lpSize = linkProperty[id].size();
////		/*const Vec3d &pid = p_[id];
////		const Vec3d &nid = n_[id];*/
////		const v3d &pid = p_.col(id);
////		const v3d &nid = n_.col(id);
////		graphKnot &gid = graph[id];
////		gid.id = id;
////		gid.principleKnot = id;
////		//gid = new graphKnot(id);
////		//gid->weight = exp(cof * (pid - v_.col(id)).sqrnorm());
////		gid.weight = exp(cof*(pid - v_.col(id)).squaredNorm());
////		gid.weightPosition = gid.weight*pid;
////		gid.weightNormal = gid.weight*nid;
////		for (ui i = 0; i < lpSize; i++)
////		{
////			ui index = adj_[id][i].id;
////			switch (linkProperty[id][i])
////			{
////			case SAME:
////				edge.emplace_back(id, index, (pid - p_.col(index)).norm());
////				break;
////			case UNIT:
////				gid.adj.insert(index);
////				graph[index].adj.insert(id);
////				break;
////			default:
////				break;
////			}
////		}
////	}
////#ifdef RUN_MESHFIELD_PARALLEL
////	);
////#endif
//
//	SerialParallelBlock(
//		ui lpSize = linkProperty[id].size();
//	    const v3d &pid = p_.col(id);
//	    const v3d &nid = n_.col(id);
//	    graphKnot &gid = graph[id];
//	    gid.id = id;
//	    gid.principleKnot = id;
//	    gid.weight = exp(cof*(pid - v_.col(id)).squaredNorm());
//	    gid.weightPosition = gid.weight*pid;
//	    gid.weightNormal = gid.weight*nid;
//	    for (ui i = 0; i < lpSize; i++)
//	    {
//		    ui index = adj_[id][i].id;
//		    switch (linkProperty[id][i])
//		    {
//		    case SAME:
//			    edge.emplace_back(id, index, (pid - p_.col(index)).norm());
//			    break;
//		    case UNIT:
//			    gid.adj.insert(index);
//			    graph[index].adj.insert(id);
//			    break;
//		    default:
//			    break;
//		    }
//	    }
//	)
//	dprint("set graph:", clock() - start, "ms");
//	start = clock();
//
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::parallel_sort(edge.begin(), edge.end(), std::less<edgeWeight>());
//#else
//	std::sort(edge.begin(), edge.end(), std::less<edgeWeight>());
//#endif
//
//	/*size = edge.size() / 2;
//	for (ui id = 1; id < size; id++)
//	{
//		std::swap(edge[id], edge[id << 1]);
//	}*/
//	//edge.erase(edge.begin() + size, edge.end());
//
//	//访问C,并将相应的顶点（在同一个block中）合并
//	for (const edgeWeight &e : edge)
//	//for(ui id = 0; id < size; id++)
//	{
//		//const edgeWeight &e = edge[id];
//		if (graph[e.i].principleKnot == graph[e.j].principleKnot)
//			continue;
//		graphKnot &gi = graph[graph[e.i].principleKnot];
//		graphKnot &gj = graph[graph[e.j].principleKnot];
//		gj.principleKnot = gi.id;
//		for (const ui &gmK : gj.minorKnot)
//		{
//			graph[gmK].principleKnot = gi.id;
//		}
//		gi.minorKnot.insert(gj.id);
//		gi.minorKnot.insert(gj.minorKnot.begin(), gj.minorKnot.end());
//		gi.weight += gj.weight;
//		gi.weightPosition += gj.weightPosition;
//		gi.weightNormal += gj.weightNormal;
//		gi.adj.insert(gj.adj.begin(), gj.adj.end());
//	}
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::concurrent_vector<edgeWeight>().swap(edge);
//#else
//	std::vector<edgeWeight>().swap(edge);
//#endif
//	//至此合并完毕，其中principleKnot指向自身的graphKnot为block的代表元，这些graphKnot将成为四边形网格的顶点，其与adj将构成边
//	dprint("merge:", clock() - start, "ms");
//	start = clock();
//	//将上述几何拓扑关系整合
//	{
//		ui leftGKSize = 0;
//		for (graphKnot &gK : graph)
//			if (!(gK.deleted()))
//				leftGKSize++;
//		ui leastBlockSize = std::max((ui)1, nv / (leftGKSize * 10));
//
//		std::vector<ui> INDEX2index(nv);
//		std::vector<ui> index2INDEX;
//		index2INDEX.reserve(nv);
//		for (graphKnot &gK : graph)
//		{
//			if (!(gK.deleted()))
//			{
//				std::set<ui> newadj;
//				for (const ui &ga : gK.adj)
//				{
//					newadj.insert(graph[ga].principleKnot);
//				}
//				gK.adj = newadj;
//				INDEX2index[gK.id] = index2INDEX.size();
//				index2INDEX.push_back(gK.id);
//			}
//		}
//
//		//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
//		//std::vector<VertexValue> poly;
//		nv = index2INDEX.size();
//		//polyFrame.reserve(polysize);
//		polyFrame = new VertexValue[nv];
//		for (ui id = 0; id < nv; id++)
//		{
//			ui i = index2INDEX[id];
//			std::set<ui> &gadj = graph[i].adj;
//			polyFrame[id].setValue(polymesh.add_vertex(Eigen2OpenMeshVector(graph[i].weightPosition / graph[i].weight)), graph[i].weightNormal);
//			std::set<ui> &padj = polyFrame[id].adj;
//			for (std::set<ui>::const_iterator gitr = gadj.begin(); gitr != gadj.end(); gitr++)
//			{
//				padj.insert(INDEX2index[*gitr]);
//			}
//			padj.erase(id);
//		}
//	}
//
//	/*for (graphKnot* &gK : graph)
//	{
//		delete gK;
//		gK = nullptr;
//	}*/
//
//	dprint("set polyFrame:", clock() - start, "ms");
//	start = clock();
//
//#pragma region post-processing for position field to remedy artifacts
//	//处理过于靠近某条边的顶点，此时认为该点在这条边上，删除这条边
//	//size = 0;
//	double minThreshold = 0.3*length;
//	nv = polymesh.n_vertices();
//	//for (VertexValue &pFi : polyFrame)
//	for (ui id = 0; id < nv; id++)
//	{
//		VertexValue &pFi = polyFrame[id];
//		Vec3d &posi = polymesh.point(pFi.vert);
//		for (ui pFj : pFi.adj)
//		{
//			Vec3d &posj = polymesh.point(polymesh.vertex_handle(pFj));
//			double dij = (posi - posj).sqrnorm();
//			for (ui pFk : polyFrame[pFj].adj)
//			{
//				if (id == pFk)
//					continue;
//				Vec3d &posk = polymesh.point(polymesh.vertex_handle(pFk));
//				double djk = (posj - posk).sqrnorm();
//				double dki = (posk - posi).sqrnorm();
//				if (djk < std::max(dij, dki))
//					continue;
//				double dcir = 0.5*(dij + djk + dki);
//				if (sqrt(dcir*(dcir - dij)*(dcir - djk)*(dcir - dki)) * 2.0 < minThreshold * djk)
//				{
//					if (!pFi.adj.count(pFk))
//						continue;
//					polymesh.set_point(pFi.vert, 0.5*(posj + posk));
//					polyFrame[pFj].adj.erase(pFk);
//					polyFrame[pFk].adj.erase(pFj);
//					goto goto20211009;
//				}
//			}
//		}
//	goto20211009:;
//		++id;
//	}
//
//
//	//删除valence<3的顶点
//	for (size = --nv; size != -1; --size)
//	{
//		VertexValue &vertex = polyFrame[size];
//		if (vertex.adj.size() > 2)
//			continue;
//		for (ui var : vertex.adj)
//		{
//			polyFrame[var].adj.erase(size);
//		}
//		vertex.normal = polyFrame[nv].normal;
//		vertex.adj = polyFrame[nv].adj;
//		for (ui var : vertex.adj)
//		{
//			polyFrame[var].adj.erase(nv);
//			polyFrame[var].adj.insert(size);
//		}
//		polymesh.set_point(polymesh.vertex_handle(size), polymesh.point(polymesh.vertex_handle(nv)));
//		polymesh.delete_vertex(polymesh.vertex_handle(nv));
//		--nv;
//	}
//	polymesh.garbage_collection();
//#pragma endregion
//	dprint("post processing:", clock() - start, "ms");
//}
//#pragma endregion
#endif

#pragma region MeshField
MeshField::hierarchy::hierarchy(const TriMesh &mesh, ui depth)
{
	//以 mesh 初始化该层
	init(mesh);
	if (depth < 2)
	{
		hierarchy_depth = 1;
		return;
	}

	//循环初始化剩下各层
	finer_hierarchy = nullptr;
	coarser_hierarchy = new hierarchy(this);
	hierarchy* h = coarser_hierarchy;
	h->finer_hierarchy = this;
	hierarchy_depth = depth;
	for (; depth > 2 && h->v.size() > 1; depth--, h = h->coarser_hierarchy)
	{
		h->coarser_hierarchy = new hierarchy(h);
		h->coarser_hierarchy->finer_hierarchy = h;
	}
	h->coarser_hierarchy = nullptr;
	hierarchy_depth -= depth - 2;

	//初始化各层 hierarchy_depth
	depth = hierarchy_depth;
	h = this->coarser_hierarchy;
	while (h)
	{
		h->hierarchy_depth = --depth;
		h = h->coarser_hierarchy;
	}
}

double MeshField::hierarchy::calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv)
{
	double area = 0;
	v3d pos = OpenMesh2EigenVector(mesh.point(tv));
	std::vector<Link> &ad = adj[tv.idx()];
	ad.reserve(mesh.valence(tv));
	for (auto tvoh : mesh.voh_range(tv))
	{
		ad.push_back(mesh.to_vertex_handle(tvoh).idx());
		if (!mesh.face_handle(tvoh).is_valid())
			continue;
		v3d s = OpenMesh2EigenVector(mesh.point(mesh.to_vertex_handle(tvoh)));
		v3d t = OpenMesh2EigenVector(mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(tvoh))));
		v3d c = (pos + s + t) / 3.0;
		area += 0.25*((s - pos).cross(c - pos).norm() + (t - pos).cross(c - pos).norm());
	}
	return area;
}

void MeshField::hierarchy::init(const TriMesh &mesh)
{
	ui nv = mesh.n_vertices();
	v.resize(3, nv);
	n.resize(3, nv);
	a.resize(nv);
	adj.resize(nv);
	//	ui id = 0;
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
	//	{
	//		OV tv = mesh.vertex_handle(id);
	//#else
	//	for (OV tv : mesh.vertices())
	//	{
	//#endif // RUN_MESHFIELD_PARALLEL
	//		v.col(id) = OpenMesh2EigenVector(mesh.point(tv));
	//		n.col(id) = OpenMesh2EigenVector(mesh.calc_normal(tv));
	//		a(id) = calcVoronoiAndAdjacency(mesh, tv);
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	});
	//#else
	//		++id;
	//    }
	//#endif
	SerialParallelBlock(
		OV tv = mesh.vertex_handle(id);
	v.col(id) = OpenMesh2EigenVector(mesh.point(tv));
	n.col(id) = OpenMesh2EigenVector(mesh.calc_normal(tv));
	a(id) = calcVoronoiAndAdjacency(mesh, tv);
	)
	graph_color();
}

void MeshField::hierarchy::graph_color()
{
	ui nsizes = adj.size();
	std::vector<int> color(nsizes, -1);
	std::vector<ui> possible_color;

	ui ncolors = 0;
	for (size_t i = 0; i < nsizes; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1);
		for each (auto link in adj[i])
		{
			int c = color[link.id];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}

	phase.resize(ncolors);
	for (size_t i = 0; i < nsizes; i++)
	{
		phase[color[i]].push_back(i);
	}
}

#include <deque>
MeshField::hierarchy::hierarchy(hierarchy *fh)
{
	const m3xd &v_ = fh->v;
	const m3xd &n_ = fh->n;
	const vxd &a_ = fh->a;
	const vvL &adj_ = fh->adj;
	ui nv = v_.cols();


#ifdef RUN_MESHFIELD_PARALLEL
	tbb::concurrent_vector<edgeWeight> edge;
#else
	std::vector<edgeWeight> edge;
#endif
	ui size = 0;
	for (const std::vector<Link> &l : adj_)
		size += l.size();
	edge.reserve(size);

	/*size = 0;
	for (const std::vector<Link> &v0 : adj_)
	{
		for (Link v1 : v0)
			edge.emplace_back(size, v1.id, n_.col(size).dot(n_.col(v1.id))*std::max(a_(size) / a_(v1.id), a_(v1.id) / a_(size)));
		size++;
	}*/

	//size = adj_.size();
//#ifdef RUN_MESHFIELD_PARALLEL
//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//#else
//	for(ui id = 0; id < nv; id ++)
//#endif
//	{
//		const vL &v0 = adj_[id];
//		for (const Link &v1 : v0)
//		{
//			const ui &vid = v1.id;
//			edge.emplace_back(id, vid, n_.col(id).dot(n_.col(vid))*std::max(a_(id) / a_(vid), a_(vid) / a_(id)));
//		}
//	}
//#ifdef RUN_MESHFIELD_PARALLEL
//	);
//#endif

	SerialParallelBlock(
		const vL &v0 = adj_[id];
	for (const Link &v1 : v0)
	{
		const ui &vid = v1.id;
		edge.emplace_back(id, vid, n_.col(id).dot(n_.col(vid))*max(a_(id) / a_(vid), a_(vid) / a_(id)));
	}
	)

#ifdef RUN_MESHFIELD_PARALLEL
	tbb::parallel_sort(edge.begin(), edge.end(), std::greater<edgeWeight>());
#else
	std::sort(edge.begin(), edge.end(), std::greater<edgeWeight>());
#endif
	std::deque<bool> mergeFlag(nv, false);
	size = 0;
	for (const edgeWeight &e : edge)
	{
		if (mergeFlag[e.i] || mergeFlag[e.j])
			continue;
		mergeFlag[e.i] = true;
		mergeFlag[e.j] = true;
		edge[size++] = e;
	}

	ui leftVertices = nv - size;

	v.resize(3, leftVertices);
	n.resize(3, leftVertices);
	a.resize(leftVertices);
	adj.resize(leftVertices);

	tofiner.resize(2, leftVertices);
	vxu &tc = fh->tocoarser;
	tc.resize(nv);

	nv = size;
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
	//#else
	//	for (ui id = 0; id < nv; id++)
	//#endif
	//	{
	//		edgeWeight &e = edge[id];
	//		double a0 = a_(e.i);
	//		double a1 = a_(e.j);
	//		double area = a0 + a1;
	//		if (area > epsilon)
	//		{
	//			v.col(id) = (v_.col(e.i)*a0 + v_.col(e.j)*a1) / area;
	//		}
	//		else
	//		{
	//			v.col(id) = (v_.col(e.i) + v_.col(e.j))*0.5;
	//		}
	//
	//		n.col(id) = (n_.col(e.i)*a0 + n_.col(e.j)*a1).normalized();
	//		a(id) = area;
	//		tofiner.col(id) << e.i, e.j;
	//		tc(e.i) = id; tc(e.j) = id;
	//	}
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	);
	//#endif

	SerialParallelBlock(
		edgeWeight &e = edge[id];
	double a0 = a_(e.i);
	double a1 = a_(e.j);
	double area = a0 + a1;
	if (area > epsilon)
	{
		v.col(id) = (v_.col(e.i)*a0 + v_.col(e.j)*a1) / area;
	}
	else
	{
		v.col(id) = (v_.col(e.i) + v_.col(e.j))*0.5;
	}

	n.col(id) = (n_.col(e.i)*a0 + n_.col(e.j)*a1).normalized();
	a(id) = area;
	tofiner(0, id) = e.i; tofiner(1, id) = e.j;
	tc(e.i) = id; tc(e.j) = id;
	)

		nv = v_.cols();
	for (ui id = 0; id < nv; id++)
	{
		if (mergeFlag[id])
			continue;

		v.col(size) = v_.col(id);
		n.col(size) = n_.col(id);
		a(size) = a_(id);
		tofiner.col(size) << id, nv;
		tc(id) = size;
		size++;
	}
	assert(size == leftVertices);

	std::swap(nv, size);
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
	//#else
	//	for (ui id = 0; id < nv; ++id)
	//#endif
	//	{
	//		std::vector<Link> &oneRing = adj[id];
	//		//std::vector<Link> oneRing;
	//		//导入精细层网格的one-ring
	//		ui i = tofiner(0, id);
	//		for (const Link &var : adj_[i])
	//		{
	//			if (tc(var.id) == id) continue;
	//			oneRing.emplace_back(tc(var.id), var.weight);
	//		}
	//		i = tofiner(1, id);
	//		if (i != size)
	//		{
	//			for (const Link &var : adj_[i])
	//			{
	//				if (tc(var.id) == id) continue;
	//				oneRing.emplace_back(tc(var.id), var.weight);
	//			}
	//		}
	//		if (oneRing.size() > 1)
	//		{
	//			//删除oneRing重复元素
	//			std::sort(oneRing.begin(), oneRing.end());
	//			std::vector<Link>::iterator itr = oneRing.begin();
	//			++itr;
	//			for (; itr != oneRing.end(); ++itr)
	//			{
	//				std::vector<Link>::iterator newItr = itr;
	//				--newItr;
	//				if (itr->id == newItr->id)
	//				{
	//					itr->weight += newItr->weight;
	//					itr = oneRing.erase(newItr);
	//				}
	//			}
	//		}
	//	}
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	);
	//#endif
	SerialParallelBlock(
		std::vector<Link> &oneRing = adj[id];
	//std::vector<Link> oneRing;
	//导入精细层网格的one-ring
	ui i = tofiner(0, id);
	for (const Link &var : adj_[i])
	{
		if (tc(var.id) == id) continue;
		oneRing.emplace_back(tc(var.id), var.weight);
	}
	i = tofiner(1, id);
	if (i != size)
	{
		for (const Link &var : adj_[i])
		{
			if (tc(var.id) == id) continue;
			oneRing.emplace_back(tc(var.id), var.weight);
		}
	}
	if (oneRing.size() > 1)
	{
		//删除oneRing重复元素
		std::sort(oneRing.begin(), oneRing.end());
		std::vector<Link>::iterator itr = oneRing.begin();
		++itr;
		for (; itr != oneRing.end(); ++itr)
		{
			std::vector<Link>::iterator newItr = itr;
			--newItr;
			if (itr->id == newItr->id)
			{
				itr->weight += newItr->weight;
				itr = oneRing.erase(newItr);
			}
		}
	}
	)
	graph_color();
}

MeshField::hierarchy::~hierarchy()
{
	finer_hierarchy = nullptr;
	if (hierarchy_depth >> 1)
	{
		delete coarser_hierarchy;
		coarser_hierarchy = nullptr;
	}
}
#pragma endregion


#pragma region OrientationField
void OrientationField::loopDivision(TriMesh &mesh, ui divideTimes)
{
	if (divideTimes < 1)
		return;
	OpenMesh::Subdivider::Uniform::LoopT<TriMesh> loopDivider;
	loopDivider.attach(mesh);
	loopDivider(divideTimes);
	loopDivider.detach();
	mesh.update_normals();
	OpenMesh::IO::write_mesh(mesh, "human2.obj");
}

void OrientationField::randomInit(TriMesh &mesh)
{
	loopDivision(mesh, 1);

	h = new hierarchy(mesh, DefaultDepth);
	hierarchy *hw = h;

	while (hw->coarser_hierarchy)
	{
		hw->o.resize(3, hw->v.cols());
		hw = hw->coarser_hierarchy;
	}
	m3xd &n_ = hw->n;
	m3xd &o_ = hw->o;
	ui cols = n_.cols();
	o_.resize(3, cols);

	for (ui i = 0; i < cols; i++)
	{
		o_.col(i) = hw->calcXAxis(n_.col(i));
	}
	dprint("Hierarchy Build Done!\n\n");

}

void OrientationField::GSOptimize()
{
	auto calcDirection = [&](const v3d *_P1, const v3d *_P2, ui *index)
	{
		double flag_max = 0;
		ui id1 = 0, id2 = 0;
		for (ui i = 0; i < 2; i++)
		{
			for (ui j = 0; j < 2; j++)
			{
				double dottemp = _P1[i].dot(_P2[j]);
				if (std::fabs(dottemp) > std::fabs(flag_max))
				{
					flag_max = dottemp;
					id1 = i;
					id2 = j;
				}
			}
		}
		index[0] = id1;
		index[1] = flag_max > 0 ? id2 : id2 + 2;
	};

	hierarchy* hn = h;
	while (hn->coarser_hierarchy)
		hn = hn->coarser_hierarchy;


	const ui iterationTimes = 5;
	while (true)
	{
		//在当前层 Gauss-Seidel 迭代，得到优化 Orientation Field
		m3xd &o_ = hn->o;
		ui nv = o_.cols();
		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
		const m3xd &n_ = hn->n;
		const vvL &adj_ = hn->adj;
		const m3xd &v_ = hn->v;


		//for (ui i = 0; i < nv; ++i)
		//{
		//	//method 1
		//	//o_.col(i) = h->calcXAxis(n_.col(i));
		//	//method 2
		//	v3d vi = v_.col(i);
		//	double t = std::max(std::fabs(vi(0)), std::fabs(vi(1)));
		//	t = std::max(t, std::fabs(vi(2)));
		//	vi = v3d(t,t,t);
		//	o_.col(i) = (vi - vi.dot(n_.col(i))*n_.col(i)).normalized();
		//	//method 3
		//	/*v3d ni = n_.col(i);
		//	double e = adj_[i].size();
		//	v3d vi = v3d(0, 0, 0);
		//	for (auto l : adj_[i])
		//	{
		//		ui j = l.id;
		//		vi += v_.col(j) / ((n_.col(j) - ni).norm() * 10 + e);
		//	}
		//	o_.col(i) = (vi - vi.dot(n_.col(i))*n_.col(i)).normalized();*/
		//}
#if 1
		for (int i = 0; i < iterationTimes; i++)
		{
#ifdef RUN_MESHFIELD_PARALLEL
			for (const std::vector<ui> &ph : hn->phase)
			{
				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
				tbb:: parallel_for(range, [&](const tbb::blocked_range<ui> &range)
				{
					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
					{
						ui id = ph[phaseIdx];
#else
			for (ui id = 0; id < nv; id++)
			{
#endif
				double sum = 0;
				v3d oid = o_.col(id);
				const v3d &nid = n_.col(id);
				for (const Link &oneRing : adj_[id])
				{
					double w = oneRing.weight;
					v3d _P1[2] = { oid, nid.cross(oid) };
					v3d _P2[2] = { o_.col(oneRing.id), n_.col(oneRing.id).cross(o_.col(oneRing.id)) };
					ui index[2];
					calcDirection(_P1, _P2, index);
					v3d direction[2] = { _P1[index[0]], _P2[index[1] & 1] * (( index[1] & 2) ? -1 : 1) };
					oid = direction[0] * sum + direction[1] * w;
					sum += w;
					oid -= nid.dot(oid)*nid;
					oid.normalize();
				}
				o_.col(id) = oid;
#ifdef RUN_MESHFIELD_PARALLEL
			}
				});
#endif
			}
	    //dprint("iteration times:", i);
		}
#endif

        //将当前层 Orientation Field 传递到精细层
        if (!(hn->finer_hierarchy))
            break;

        const m2xu &tofiner = hn->tofiner;
        m3xd &o_finer = hn->finer_hierarchy->o;
        const m3xd &n_finer = hn->finer_hierarchy->n;
        const ui finer_size = n_finer.cols();

#ifdef RUN_MESHFIELD_PARALLEL
        tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
#else
        for (ui id = 0; id < nv; id++)
#endif
        {
	        ui finer_id = tofiner(0, id);
	        o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();
	        finer_id = tofiner(1, id);
	        if (finer_id != finer_size)
	         	o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();
        }
#ifdef RUN_MESHFIELD_PARALLEL
        );
#endif
		//SerialParallelBlock(
		//	ui finer_id = tofiner.col(id)(0);// (0, id);
		//o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();
		//finer_id = tofiner.col(id)(1);// (1, id);
		//if (finer_id != finer_size)
		//	o_finer.col(finer_id) = (o_.col(id) - n_finer.col(finer_id).dot(o_.col(id))*n_finer.col(finer_id)).normalized();)

        hn = hn->finer_hierarchy;
    }

	dprint("\nOrientation Field Build Done!\n\n");

	//计算奇异点
#if 1
	//compute singularities
	m3xd &o_ = h->o;
	const m3xd &n_ = hn->n;
	const vvL &adj_ = h->adj;
	ui nv = o_.cols();
	ui table[2][4] = { 0,3,2,1,1,0,3,2 };
#ifdef RUN_MESHFIELD_PARALLEL
	tbb::concurrent_vector<ui> sing;
	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
#else
	std::vector<ui> sing;
	for (ui id = 0; id < nv; id++)
#endif
	{
		std::vector<ui> oneRing;
		ui as = adj_[id].size();
		oneRing.reserve(as + 1);
		for (const Link & or : adj_[id])
		{
			oneRing.push_back(or .id);
		}
		oneRing.push_back(oneRing.front());
		ui sum = 0;
		ui index[2];
		v3d _P1[2]; v3d _P2[2];
		for (ui i = 0; i < as; ++i)
		{
			_P1[0] = o_.col(oneRing[i]); _P1[1] = n_.col(oneRing[i]).cross(o_.col(oneRing[i]));
			_P2[0] = o_.col(oneRing[i + 1]); _P2[1] = n_.col(oneRing[i + 1]).cross(o_.col(oneRing[i + 1]));
			calcDirection(_P1, _P2, index);
			sum += table[index[0]][index[1]];
		}
		if (sum & 3)
			sing.push_back(id);
	}
#ifdef RUN_MESHFIELD_PARALLEL
	);
	for (auto &s : sing)
	{
		singularities.push_back(s);
	}
#else
	singularities = sing;
#endif
#endif

#if 0
	OV firstVertex = mesh_ptr->vertex_handle(7547);
	ui secondId;
	ui thirdId;
	{
		for (const ui &s : singularities)
		{
			if (mesh_ptr->find_halfedge(firstVertex, mesh_ptr->vertex_handle(s)).is_valid())
			{
				secondId = s;
				break;
			}
		}
		thirdId = mesh_ptr->to_vertex_handle(mesh_ptr->next_halfedge_handle(mesh_ptr->find_halfedge(firstVertex, mesh_ptr->vertex_handle(secondId)))).idx();
		bool flag = true;
		for (const ui &s : singularities)
		{
			if (s == thirdId)
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			thirdId = mesh_ptr->to_vertex_handle(mesh_ptr->next_halfedge_handle(mesh_ptr->find_halfedge(mesh_ptr->vertex_handle(secondId), firstVertex))).idx();
		}
	}
	ui vId[3] = { firstVertex.idx(),secondId,thirdId };
	double dif[3];
	ui index[2];
	v3d _P1[2]; v3d _P2[2];
	for (ui i = 0; i < 3; ++i)
	{
		_P1[0] = o_.col(vId[i]); _P1[1] = n_.col(vId[i]).cross(o_.col(vId[i]));
		_P2[0] = o_.col(vId[(i + 1) % 3]); _P2[1] = n_.col(vId[(i + 1) % 3]).cross(o_.col(vId[(i + 1) % 3]));
		calcDirection(_P1, _P2, index);
		v3d direction[2] = { _P1[index[0]], _P2[index[1] & 1] * ((index[1] & 2) ? -1 : 1) };
		dif[i] = acos(direction[0].dot(direction[1]));
	}
	ui fixedId = 6402;
	for (ui i = 0; i < 3; ++i)
	{
		//if (dif[i] < dif[(i + 1) % 3] && dif[i] < dif[(i + 2) % 3])
		if(fixedId == vId[(i + 2) % 3])
		{
			fixedId = vId[(i + 2) % 3];
			v3d onew = (o_.col(i) + o_.col((i + 1) % 3)) / 2.0;
			onew -= n_.col(fixedId).dot(onew)*n_.col(fixedId);
			onew.normalize();
			o_.col(fixedId) = onew;
			//o_.col(fixedId) << onew(0), onew(1), onew(2);
			break;
		}
	}
	for (int i = 0; i < iterationTimes; i++)
	{
#ifdef RUN_MESHFIELD_PARALLEL
		for (const std::vector<ui> &ph : hn->phase)
		{
			tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
			tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
			{
				for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
				{
					ui id = ph[phaseIdx];
#else
		for (ui id = 0; id < nv; id++)
		{
#endif
			double sum = 0;
			v3d oid = o_.col(id);
			const v3d &nid = n_.col(id);
			for (const Link &oneRing : adj_[id])
			{
				double w = oneRing.weight;
				v3d _P1[2] = { oid, nid.cross(oid) };
				v3d _P2[2] = { o_.col(oneRing.id), n_.col(oneRing.id).cross(o_.col(oneRing.id)) };
				ui index[2];
				calcDirection(_P1, _P2, index);
				v3d direction[2] = { _P1[index[0]], _P2[index[1] & 1] * ((index[1] & 2) ? -1 : 1) };
				oid = direction[0] * sum + direction[1] * w;
				sum += w;
				oid -= nid.dot(oid)*nid;
				oid.normalize();
			}
			bool f = true;
			for (ui i = 0; i < 3; ++i)
			{
				if (id == vId[i])
				{
					f = false;
					break;
				}
			}
			if (f)
				o_.col(id) = oid;
#ifdef RUN_MESHFIELD_PARALLEL
		}
			});
#endif
		}
	}
#endif

}
#pragma endregion


#pragma region PositionField
void PositionField::randomInit(TriMesh &mesh)
{
	clock_t start = clock();
	if (!h)
	{
		dprint("Orientation Field need to be build first !!!");
		exit(0);
	}
	hierarchy* hw = h;
	while (hw->coarser_hierarchy)
	{
		hw->p.resize(3, hw->v.cols());
		hw = hw->coarser_hierarchy;
	}
	hw->p = hw->v;

	length = 0;
	for (OE te : mesh.edges())
		length += mesh.calc_edge_length(te);

	length /= mesh.n_edges() / scale;
	dprint("预计算", clock() - start);
}

void PositionField::GSOptimize()
{
	auto middlePoint = [&](const v3d &v0, const v3d &n0, const v3d &v1, const v3d &n1)
	{
		double n0n1 = n0.dot(n1);
		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
	};
	double inv_length = 1.0 / length;
	auto calcPosition = [&](const v3d &p0, const v3d &v0, const v3d &n0, const v3d &o0,
		const v3d &p1, const v3d &v1, const v3d &n1, const v3d &o1, v3d *position)
	{
		const v3d mP = middlePoint(v0, n0, v1, n1);
		v3d y0 = n0.cross(o0);
		v3d origin0 = p0 + length * (o0*std::floor(o0.dot(mP - p0) * inv_length) + y0 * std::floor(y0.dot(mP - p0) * inv_length));
		v3d y1 = n1.cross(o1);
		v3d origin1 = p1 + length * (o1*std::floor(o1.dot(mP - p1) * inv_length) + y1 * std::floor(y1.dot(mP - p1) * inv_length));

		v3d updateP0[4], updateP1[4];
		for (ui i = 0; i < 4; i++)
		{
			updateP0[i] = origin0 + ((i & 1)*o0 + ((i & 2) >> 1)*y0)*length;
			updateP1[i] = origin1 + ((i & 1)*o1 + ((i & 2) >> 1)*y1)*length;
		}
		double dis = DBL_MAX;
		ui index[2];
		for (ui i = 0; i < 4; i++)
		{
			for (ui j = 0; j < 4; j++)
			{
				double d = (updateP0[i] - updateP1[j]).squaredNorm();
				if (d < dis)
				{
					dis = d;
					index[0] = i; index[1] = j;
				}
			}
		}
		position[0] = updateP0[index[0]]; position[1] = updateP1[index[1]];
	};
	auto positionRound = [&](v3d &p, const v3d &v, const v3d &n, const v3d &o)
	{
		const v3d y = n.cross(o);
		p += length * (std::floor((v - p).dot(o) * inv_length + 0.5)*o + std::floor((v - p).dot(y) * inv_length + 0.5)*y);
	};

	hierarchy* hn = h;
	while (hn->coarser_hierarchy)
		hn = hn->coarser_hierarchy;

	const ui iterTimes = 10;
	while (true)
	{
		//在当前层 Gauss-Seidel 迭代，得到优化 Position Field
		m3xd &p_ = hn->p;
		ui nv = p_.cols();
		dprint("Hierarchy Number:", hn->hierarchy_depth, "\tVertex Number:", nv);
		clock_t start = clock();
		const m3xd &v_ = hn->v;
		const m3xd &n_ = hn->n;
		const m3xd &o_ = hn->o;
		const vvL &adj_ = hn->adj;
		for (int i = 0; i < /*hn->hierarchy_depth*/iterTimes; i++)
		{
#ifdef RUN_MESHFIELD_PARALLEL
			for (const std::vector<ui> &ph : hn->phase)
			{
				tbb::blocked_range<ui> range(0, static_cast<ui>(ph.size()), grainSize);
				tbb::parallel_for(range, [&](const tbb::blocked_range<ui> &range)
				{
					for (ui phaseIdx = range.begin(); phaseIdx < range.end(); ++phaseIdx)
					{
						ui id = ph[phaseIdx];
#else
			for (ui id = 0; id < nv; id++)
			{
#endif
				double sum = 0;
				v3d pid = p_.col(id);
				const v3d &vid = v_.col(id);
				const v3d &nid = n_.col(id);
				const v3d &oid = o_.col(id);
				for (const Link &oneRing : adj_[id])
				{
					double w = oneRing.weight;
					ui index = oneRing.id;
					v3d position[2];
					calcPosition(pid, vid, nid, oid, p_.col(index), v_.col(index), n_.col(index), o_.col(index), position);
					pid = position[0] * sum + position[1] * w;
					sum += w;
					pid /= sum;
					pid -= nid.dot(pid - vid)*nid;
				}
				positionRound(pid, vid, nid, oid);
				p_.col(id) = pid;
#ifdef RUN_MESHFIELD_PARALLEL
			        }
				});
#endif
			}
		}

		dprint("time0:", clock() - start);
		start = clock();

		//将当前层 Position Field 传递到精细层
		if (!(hn->finer_hierarchy))
			break;

		const m2xu &tofiner = hn->tofiner;
		const m3xd &v_finer = hn->finer_hierarchy->v;
		const m3xd &n_finer = hn->finer_hierarchy->n;
		m3xd &p_finer = hn->finer_hierarchy->p;
		const ui finer_size = n_finer.cols();

//#ifdef RUN_MESHFIELD_PARALLEL
//        tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
//#else
//        for (ui id = 0; id < nv; id++)
//#endif
//		{
//			ui finer_id = tofiner(0, id);
//			p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
//			finer_id = tofiner(1, id);
//			if (finer_id != finer_size)
//				p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
//		}
//#ifdef RUN_MESHFIELD_PARALLEL
//        );
//#endif
		SerialParallelBlock(
			ui finer_id = tofiner.col(id)(0);// (0, id);
		p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
		finer_id = tofiner.col(id)(1);// (1, id);
		if (finer_id != finer_size)
			p_finer.col(finer_id) = p_.col(id) - n_finer.col(finer_id).dot(p_.col(id) - v_finer.col(finer_id))*n_finer.col(finer_id);
		)
        hn = hn->finer_hierarchy;
		dprint("time1:", clock() - start);
	}
	dprint("\nPosition Field Build Done!\n\n");
}

#include <set>
#include <OpenMesh\Tools\Subdivider\Uniform\CatmullClarkT.hh>
void PositionField::extractMesh(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame)
{
	timeRecorder tr;
	computePolyFrame(polymesh, pos, polyFrame);
	tr.out("compute polyframe:");
	ui nv = polymesh.n_vertices();
	VertexValue *rrr = new VertexValue[nv];
	SerialParallelBlock(rrr[id] = polyFrame[id];)

#if 1
	SerialParallelBlock(
		VertexValue &pF = polyFrame[id];
	std::list<ui> &orientedAdj = pF.orientedAdj;
	for (const ui &va : pF.adj) {
		orientedAdj.push_back(va);
	}
	v3d x = h->calcXAxis(pF.normal);
	v3d y = pF.normal.cross(x);
	orientedAdj.sort([&](const ui &id0, const ui &id1) {
		v3d dir0 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id0)) - polymesh.point(pF.vert));
		v3d dir1 = OpenMesh2EigenVector(polymesh.point(polymesh.vertex_handle(id1)) - polymesh.point(pF.vert));
		return std::atan2(y.dot(dir0), x.dot(dir0)) > std::atan2(y.dot(dir1), x.dot(dir1));
	});
	)
	tr.out("orient:");
	//auto ppp = polyFrame[974];

	for (ui id = 0; id < nv; id++)
	{
		//dprint(id);
		//if (id == 669)//human.obj
		if (id == 841)//vase.obj
		{
#if 1
			int p = 0;
#else
			break;
#endif
		}
		VertexValue &pF = polyFrame[id];
		for (const ui &oh : pF.adj)
		{
			std::vector<OV> loop;
			loop.reserve(4);
			loop.push_back(polymesh.vertex_handle(oh));
			ui previousId = id;
			do
			{
				std::list<ui> &currentOrientedAdj = polyFrame[loop.back().idx()].orientedAdj;
				std::list<ui>::iterator pointer = currentOrientedAdj.begin();
				while (*pointer != previousId) {
					++pointer;
				}
				previousId = loop.back().idx();
				pointer = ++pointer == currentOrientedAdj.end() ? currentOrientedAdj.begin() : pointer;
				polyFrame[previousId].adj.erase(*pointer);
				loop.push_back(polymesh.vertex_handle(*pointer));
				if (*pointer == id)
				{
					break;
				}
			} while (true);
			polymesh.add_face(loop);
		}
		pF.adj.clear();
	}
	/*for (auto tf : polymesh.faces())
	{
		std::cout << "\n" << tf.idx() << ":   ";
		for (auto tfv : polymesh.fv_range(tf))
		{
			std::cout << tfv.idx() << ",";
		}
	}*/
	//dprint("extract mesh:", clock() - start, "ms");
	tr.out("get face:");
#endif
	tr.sum("extract mesh:");
	delete[] polyFrame;
	polyFrame = rrr;

	OpenMesh::Subdivider::Uniform::CatmullClarkT<PolyMesh> loopDivider;
	loopDivider.attach(polymesh);
	loopDivider(1);
	loopDivider.detach();
	polymesh.update_normals();
}

void PositionField::computeLinkProperty(vvlt &linkProperty)
{
	double inv_length = 1.0 / length;
	auto middlePoint = [&](const v3d &v0, const v3d &n0, const v3d &v1, const v3d &n1)
	{
		double n0n1 = n0.dot(n1);
		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
	};
	auto calcDirection = [&](const v3d &o0, const v3d &n0, const v3d &o1, const v3d &n1, v3d *direction)
	{
		v3d _P1[2] = { o0,n0.cross(o0) };
		v3d _P2[2] = { o1,n1.cross(o1) };
		double flag_max = 0;
		ui id1 = 0, id2 = 0;
		for (ui i = 0; i < 2; i++)
		{
			for (ui j = 0; j < 2; j++)
			{
				double dottemp = _P1[i].dot(_P2[j]);
				if (std::fabs(dottemp) > std::fabs(flag_max))
				{
					flag_max = dottemp;
					id1 = i;
					id2 = j;
				}
			}
		}
		direction[0] = _P1[id1]; direction[1] = _P2[id2] * (flag_max > 0 ? 1 : -1);
	};
	auto calcProperty = [&](const v3d &p0, const v3d &v0, const v3d &n0, const v3d &o0,
		const v3d &p1, const v3d &v1, const v3d &n1, const v3d &o1)
	{
		const v3d mP = middlePoint(v0, n0, v1, n1);
		signed short linkP[4];
		v3d y0 = n0.cross(o0);
		linkP[0] = static_cast<signed short>(std::floor(o0.dot(mP - p0) * inv_length));
		linkP[1] = static_cast<signed short>(std::floor(y0.dot(mP - p0) * inv_length));
		v3d y1 = n1.cross(o1);
		linkP[2] = static_cast<signed short>(std::floor(o1.dot(mP - p1) * inv_length));
		linkP[3] = static_cast<signed short>(std::floor(y1.dot(mP - p1) * inv_length));
		v3d position0[4], position1[4];
		for (signed short i = 0; i < 4; i++)
		{
			position0[i] = p0 + ((linkP[0] + (i & 1)) * o0 + (linkP[1] + ((i & 2) >> 1)) * y0)*length;
			position1[i] = p1 + ((linkP[2] + (i & 1)) * o1 + (linkP[3] + ((i & 2) >> 1)) * y1)*length;
		}
		double dis = DBL_MAX;
		signed short id0 = 4, id1 = 4;
		for (signed short i = 0; i < 4; i++)
		{
			for (signed short j = 0; j < 4; j++)
			{
				double d = (position0[i] - position1[j]).squaredNorm();
				if (d < dis)
				{
					dis = d; id0 = i; id1 = j;
				}
			}
		}
		int fdsa = std::abs(linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1));
		int flgl = std::abs(linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1));
		if (std::abs(linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1)) > 1 ||
			std::abs(linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1)) > 1)
		{
			int pfr = 0;
		}
		switch (std::abs(linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1)) +
			std::abs(linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1)))
		{
		case 0:
			return SAME;
		case 1:
			return UNIT;
		default:
			return DISCONNECTION;
		}
	};
	const m3xd &p_ = h->p;
	const m3xd &v_ = h->v;
	const m3xd &n_ = h->n;
	const m3xd &o_ = h->o;
	const vvL &adj_ = h->adj;

	ui nv = p_.cols();
	linkProperty.resize(nv);

	//#ifdef RUN_MESHFIELD_PARALLEL
	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
	//#else
	//	for (ui id = 0; id < nv; id++)
	//#endif
	//	{
	//		//linkProperty[id].reserve(adj_[id].size());
	//		linkProperty[id].resize(adj_[id].size(),DISCONNECTION);
	//	}
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	);
	//#endif
	SerialParallelBlock(linkProperty[id].resize(adj_[id].size(), DISCONNECTION);)

		//#ifdef RUN_MESHFIELD_PARALLEL
		//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
		//#else
		//	for (ui id = 0; id < nv; ++id)
		//#endif
		//	{
		//		const v3d &pid = p_.col(id);
		//		const v3d &vid = v_.col(id);
		//		const v3d &nid = n_.col(id);
		//		const v3d &oid = o_.col(id);
		//		ui adjSize = adj_[id].size();
		//		for (ui i = 0; i < adjSize; i++)
		//		{
		//			ui index = adj_[id][i].id;
		//			//if (id > index)
		//			//{
		//			//	//linkProperty[id].push_back(UNDEFINED);
		//			//	/*for (const Link &var : adj_[id])
		//			//	{
		//			//		if (var.id == id)
		//			//		{
		//			//			linkProperty[id].push_back(linkProperty[index][])
		//			//		}
		//			//	}*/
		//			//	for (ui j = 0; j < adj_[index].size(); ++j)
		//			//	{
		//			//		if (adj_[index][j].id == id)
		//			//		{
		//			//			linkProperty[id].push_back(linkProperty[index][j]);
		//			//			break;
		//			//		}
		//			//	}
		//			//}
		//			//else
		//			//{
		//			//	v3d direction[2];
		//			//	calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
		//			//	linkProperty[id].push_back(calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]));
		//			//	//linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
		//			//}
		//
		//			if (id > index)
		//				continue;
		//			v3d direction[2];
		//			calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
		//			//linkProperty[id].push_back(calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]));
		//			linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
		//			/*for (ui j = 0; j < adj_[index].size(); ++j)
		//			{
		//				if (adj_[index][j].id == id)
		//				{
		//					linkProperty[index][j] = linkProperty[id][i];
		//					break;
		//				}
		//			}*/
		//		}
		//	}
		//#ifdef RUN_MESHFIELD_PARALLEL
		//	);
		//#endif
	SerialParallelBlock(
		const v3d &pid = p_.col(id);
	const v3d &vid = v_.col(id);
	const v3d &nid = n_.col(id);
	const v3d &oid = o_.col(id);
	ui adjSize = adj_[id].size();
	for (ui i = 0; i < adjSize; i++)
	{
		const ui &index = adj_[id][i].id;
		if (id > index)
			continue;
		v3d direction[2];
		calcDirection(oid, nid, o_.col(index), n_.col(index), direction);
		linkProperty[id][i] = calcProperty(pid, vid, nid, direction[0], p_.col(index), v_.col(index), n_.col(index), direction[1]);
	}
	)
}

void PositionField::computePolyFrame(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame)
{
	timeRecorder tr;
	vvlt linkProperty;
	computeLinkProperty(linkProperty);
	//dprint("compute link property:", clock() - start, "ms");
	tr.out("compute link property:");
	lP = linkProperty;
	const m3xd &v_ = h->v;
	const m3xd &n_ = h->n;
	const m3xd &p_ = h->p;
	const std::vector<std::vector<Link>> &adj_ = h->adj;
	ui nv = v_.cols();

	//构造图，图的顶点是原网格顶点，边是带有UNIT标记的原网格边
	struct graphKnot
	{
		ui id;
		v3d weightPosition;
		v3d weightNormal;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			std::set<ui> adj;
		double weight;
		ui principleKnot;
		std::set<ui> minorKnot;
		graphKnot() {}
		graphKnot(ui i) : id(i), principleKnot(i) { }
		~graphKnot() {}
		inline bool deleted()
		{
			return principleKnot != id;
		}
	};

	std::vector<graphKnot> graph(nv);
	//计算需要合并的点集C,将对应需要合并的边加入edge
#ifdef RUN_MESHFIELD_PARALLEL
	tbb::concurrent_vector<edgeWeight> edge;
#else
	std::vector<edgeWeight> edge;
#endif
	ui size = 0;
	for (const std::vector<LinkType> &lt : linkProperty)
		for (const LinkType &l : lt)
			if (l == SAME)
				size++;
	edge.reserve(size);

	const double cof = -9.0 / (length * length);
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id)
	//#else
	//	for (ui id = 0; id < nv; id++)
	//#endif
	//	{
	//		ui lpSize = linkProperty[id].size();
	//		/*const Vec3d &pid = p_[id];
	//		const Vec3d &nid = n_[id];*/
	//		const v3d &pid = p_.col(id);
	//		const v3d &nid = n_.col(id);
	//		graphKnot &gid = graph[id];
	//		gid.id = id;
	//		gid.principleKnot = id;
	//		//gid = new graphKnot(id);
	//		//gid->weight = exp(cof * (pid - v_.col(id)).sqrnorm());
	//		gid.weight = exp(cof*(pid - v_.col(id)).squaredNorm());
	//		gid.weightPosition = gid.weight*pid;
	//		gid.weightNormal = gid.weight*nid;
	//		for (ui i = 0; i < lpSize; i++)
	//		{
	//			ui index = adj_[id][i].id;
	//			switch (linkProperty[id][i])
	//			{
	//			case SAME:
	//				edge.emplace_back(id, index, (pid - p_.col(index)).norm());
	//				break;
	//			case UNIT:
	//				gid.adj.insert(index);
	//				graph[index].adj.insert(id);
	//				break;
	//			default:
	//				break;
	//			}
	//		}
	//	}
	//#ifdef RUN_MESHFIELD_PARALLEL
	//	);
	//#endif

	SerialParallelBlock(
		ui lpSize = linkProperty[id].size();
	const v3d &pid = p_.col(id);
	const v3d &nid = n_.col(id);
	graphKnot &gid = graph[id];
	gid.id = id;
	gid.principleKnot = id;
	gid.weight = exp(cof*(pid - v_.col(id)).squaredNorm());
	gid.weightPosition = gid.weight*pid;
	gid.weightNormal = gid.weight*nid;
	for (ui i = 0; i < lpSize; i++)
	{
		ui index = adj_[id][i].id;
		switch (linkProperty[id][i])
		{
		case SAME:
			edge.emplace_back(id, index, (pid - p_.col(index)).norm());
			break;
		case UNIT:
			gid.adj.insert(index);
			graph[index].adj.insert(id);
			break;
		default:
			break;
		}
	}
	)
	//dprint("set graph:", clock() - start, "ms");
	tr.out("set graph:");

#ifdef RUN_MESHFIELD_PARALLEL
	tbb::parallel_sort(edge.begin(), edge.end(), std::less<edgeWeight>());
#else
	std::sort(edge.begin(), edge.end(), std::less<edgeWeight>());
#endif

	/*size = edge.size() / 2;
	for (ui id = 1; id < size; id++)
	{
		std::swap(edge[id], edge[id << 1]);
	}*/
	//edge.erase(edge.begin() + size, edge.end());

	//访问C,并将相应的顶点（在同一个block中）合并
	for (const edgeWeight &e : edge)
		//for(ui id = 0; id < size; id++)
	{
		//const edgeWeight &e = edge[id];
		if (graph[e.i].principleKnot == graph[e.j].principleKnot)
			continue;
		graphKnot &gi = graph[graph[e.i].principleKnot];
		graphKnot &gj = graph[graph[e.j].principleKnot];
		gj.principleKnot = gi.id;
		for (const ui &gmK : gj.minorKnot)
		{
			graph[gmK].principleKnot = gi.id;
		}
		gi.minorKnot.insert(gj.id);
		gi.minorKnot.insert(gj.minorKnot.begin(), gj.minorKnot.end());
		gi.weight += gj.weight;
		gi.weightPosition += gj.weightPosition;
		gi.weightNormal += gj.weightNormal;
		gi.adj.insert(gj.adj.begin(), gj.adj.end());
	}
#ifdef RUN_MESHFIELD_PARALLEL
	tbb::concurrent_vector<edgeWeight>().swap(edge);
#else
	std::vector<edgeWeight>().swap(edge);
#endif
	//至此合并完毕，其中principleKnot指向自身的graphKnot为block的代表元，这些graphKnot将成为四边形网格的顶点，其与adj将构成边
	tr.out("merge:");
	//将上述几何拓扑关系整合
	{
		ui leftGKSize = 0;
		for (graphKnot &gK : graph)
			if (!(gK.deleted()))
				leftGKSize++;
		ui leastBlockSize = max((ui)1, nv / (leftGKSize * 10));

		std::vector<ui> INDEX2index(nv);
		std::vector<ui> index2INDEX;
		index2INDEX.reserve(nv);
		for (graphKnot &gK : graph)
		{
			if (!(gK.deleted()))
			{
				std::set<ui> newadj;
				for (const ui &ga : gK.adj)
				{
					newadj.insert(graph[ga].principleKnot);
				}
				gK.adj = newadj;
				INDEX2index[gK.id] = index2INDEX.size();
				index2INDEX.push_back(gK.id);
			}
		}

		//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
		//std::vector<VertexValue> poly;
		nv = index2INDEX.size();
		//polyFrame.reserve(polysize);
		polyFrame = new VertexValue[nv];
		for (ui id = 0; id < nv; id++)
		{
			ui i = index2INDEX[id];
			std::set<ui> &gadj = graph[i].adj;
			polyFrame[id].setValue(polymesh.add_vertex(Eigen2OpenMeshVector(graph[i].weightPosition / graph[i].weight)), graph[i].weightNormal);
			std::set<ui> &padj = polyFrame[id].adj;
			for (std::set<ui>::const_iterator gitr = gadj.begin(); gitr != gadj.end(); gitr++)
			{
				padj.insert(INDEX2index[*gitr]);
			}
			padj.erase(id);
		}
	}

	
	tr.out("set polyFrame:");

#pragma region post-processing for position field to remedy artifacts
	//处理过于靠近某条边的顶点，此时认为该点在这条边上，删除这条边
	//size = 0;
	double minThreshold = 0.3*length;
	nv = polymesh.n_vertices();
	//for (VertexValue &pFi : polyFrame)
	for (ui id = 0; id < nv; id++)
	{
		VertexValue &pFi = polyFrame[id];
		Vec3d &posi = polymesh.point(pFi.vert);
		for (ui pFj : pFi.adj)
		{
			Vec3d &posj = polymesh.point(polymesh.vertex_handle(pFj));
			double dij = (posi - posj).sqrnorm();
			for (ui pFk : polyFrame[pFj].adj)
			{
				if (id == pFk)
					continue;
				Vec3d &posk = polymesh.point(polymesh.vertex_handle(pFk));
				double djk = (posj - posk).sqrnorm();
				double dki = (posk - posi).sqrnorm();
				if (djk < max(dij, dki))
					continue;
				double dcir = 0.5*(dij + djk + dki);
				if (sqrt(dcir*(dcir - dij)*(dcir - djk)*(dcir - dki)) * 2.0 < minThreshold * djk)
				{
					if (!pFi.adj.count(pFk))
						continue;
					polymesh.set_point(pFi.vert, 0.5*(posj + posk));
					polyFrame[pFj].adj.erase(pFk);
					polyFrame[pFk].adj.erase(pFj);
					goto goto20211009;
				}
			}
		}
	goto20211009:;
		++id;
	}


	//删除valence<3的顶点
	for (size = --nv; size != -1; --size)
	{
		VertexValue &vertex = polyFrame[size];
		if (vertex.adj.size() > 2)
			continue;
		for (ui var : vertex.adj)
		{
			polyFrame[var].adj.erase(size);
		}
		vertex.normal = polyFrame[nv].normal;
		vertex.adj = polyFrame[nv].adj;
		for (ui var : vertex.adj)
		{
			polyFrame[var].adj.erase(nv);
			polyFrame[var].adj.insert(size);
		}
		polymesh.set_point(polymesh.vertex_handle(size), polymesh.point(polymesh.vertex_handle(nv)));
		polymesh.delete_vertex(polymesh.vertex_handle(nv));
		--nv;
	}
	polymesh.garbage_collection();
#pragma endregion
	
	tr.out("post processing:");
	//tr.sum("compute polyFrame time:");
}

void PositionField::searchBFSTree(TriMesh &trimesh)
{
#if 0
	double inv_length = 1.0 / length;
	const m3xd &v_ = h->v;
	const m3xd &n_ = h->n;
	const m3xd &o_ = h->o;
	const m3xd &p_ = h->p;

	typedef signed short ss;
	typedef unsigned short us;

	auto middlePoint = [&](const v3d &v0, const v3d &n0, const v3d &v1, const v3d &n1)
	{
		double n0n1 = n0.dot(n1);
		return 0.5*(v0 + v1) - 0.5 / (1 - n0n1 * n0n1 + 0.0001)
			*((n0 + n0n1 * n1).dot(v1 - v0)*n0 + (n1 + n0n1 * n0).dot(v0 - v1)*n1);
	};
	auto calcRotation = [&](const v3d &o0, /*const v3d &n0, */const v3d &o1, const v3d &y1)
	{
		us rot;
		v3d testP[4] = { o1, y1, -o1, -y1 };
		double flag_max = -1.0;
		for (unsigned short i = 0; i < 4; ++i)
		{
			double dot = o0.dot(testP[i]);
			if (dot > flag_max)
			{
				flag_max = dot;
				rot = i;
			}
		}
		return rot;
	};
	auto rotateOnVector = [&](const v3d &oi, const v3d &yi, const unsigned short &rotateTimes)
	{
		v3d outputVector = oi;
		if (rotateTimes & 1)
			outputVector = yi;
		if ((rotateTimes >> 1) & 1)
			outputVector *= -1.0;
		return outputVector;
	};
	auto calcShift = [&](const v3d &p0, const v3d &o0, const v3d &y0, const v3d &p1, const v3d &o1, const v3d &y1, const v3d &mP, ss *shift)
	{
		ss linkP[4];
		linkP[0] = static_cast<ss>(std::floor(o0.dot(mP - p0) * inv_length));
		linkP[1] = static_cast<ss>(std::floor(y0.dot(mP - p0) * inv_length));
		linkP[2] = static_cast<ss>(std::floor(o1.dot(mP - p1) * inv_length));
		linkP[3] = static_cast<ss>(std::floor(y1.dot(mP - p1) * inv_length));
		v3d position0[4], position1[4];
		for (ss i = 0; i < 4; i++)
		{
			position0[i] = p0 + ((linkP[0] + (i & 1)) * o0 + (linkP[1] + ((i & 2) >> 1)) * y0)*length;
			position1[i] = p1 + ((linkP[2] + (i & 1)) * o1 + (linkP[3] + ((i & 2) >> 1)) * y1)*length;
		}
		double dis = DBL_MAX;
		ss id0 = 4, id1 = 4;
		for (ss i = 0; i < 4; i++)
		{
			for (ss j = 0; j < 4; j++)
			{
				double d = (position0[i] - position1[j]).squaredNorm();
				if (d < dis)
				{
					dis = d; id0 = i; id1 = j;
				}
			}
		}
		shift[0] = linkP[0] + (id0 & 1) - linkP[2] - (id1 & 1);
		shift[1] = linkP[1] + ((id0 & 2) >> 1) - linkP[3] - ((id1 & 2) >> 1);
	};
	
	struct edgeData {
		ss dir[2];
		us rot;
		//bool sign;
		us reference;
		edgeData()/* : sign(false)*/ {}
	};
	ui ne = trimesh.n_edges();
	std::vector<edgeData> EdgesData(ne);
	for (auto &te : trimesh.edges())
	{
		edgeData &anedge = EdgesData[te.idx()];
		auto &edir = anedge.dir;
		us &erot = anedge.rot;

		ui ev[2] = { te.v0().idx(),te.v1().idx() };
		if (ev[0] > ev[1]) 
			std::swap(ev[0], ev[1]);

		const v3d &o0 = o_.col(ev[0]);
		const v3d &o1 = o_.col(ev[1]);
		const v3d &n0 = n_.col(ev[0]);
		const v3d &n1 = n_.col(ev[1]);
		v3d y0 = n0.cross(o0);
		v3d y1 = n1.cross(o1);
		erot = calcRotation(o0, o1, y1);
		//calcShift(p_.col(ev[0]), o0, y0, p_.col(ev[1]), o1, y1, middlePoint(v_.col(ev[0]), n0, v_.col(ev[1]), n1), edir);
		calcShift(p_.col(ev[0]), o0, y0, p_.col(ev[1]), rotateOnVector(o1, y1, erot), rotateOnVector(o1, y1, erot + 1), 
			middlePoint(v_.col(ev[0]), n0, v_.col(ev[1]), n1), edir);
	}

	struct faceData {
		ui vert[3];
		ui edge[3];
		//bool sign;
		faceData() {}
	};
	ui nf = trimesh.n_faces();
	std::vector<faceData> FacesData(nf);

	OpenMesh::FPropHandleT<bool> visited;
	trimesh.add_property(visited);
	for (auto &tf : trimesh.faces())
	{
		trimesh.property(visited, tf) = false;
	}

	std::queue<OF> tree;
	tree.push(trimesh.face_handle(0));
	trimesh.property(visited, trimesh.face_handle(0)) = true;
	//FacesData[0].sign = true;
	us id = 0;
	for (OH &tfh : trimesh.fh_range(trimesh.face_handle(0)))
	{
		FacesData[0].vert[id] = trimesh.from_vertex_handle(tfh).idx();
		FacesData[0].edge[id] = trimesh.edge_handle(tfh).idx();
		++id;
	}
	//EdgesData[FacesData[0].edge[0]].sign = true;
	EdgesData[FacesData[0].edge[0]].reference = 0;

	typedef Eigen::Matrix<ui, 4, 1> v4u;
	std::vector<v4u> edge_to_constraints(ne << 1, v4u(invalid, 0, invalid, 0));
	std::vector<int> initial(nf << 1, 0);
	while (!tree.empty())
	{
		OF &tf = tree.front();
		ui doubleFid = tf.idx() << 1;
		faceData &facedata = FacesData[tf.idx()];
		auto &fedge = facedata.edge;
		auto &fvert = facedata.vert;

		//us referenceRotation = 0;//与上一个面的相邻边的指标
		//while (!EdgesData[fedge[referenceRotation]].sign)
		//{
		//	++referenceRotation;
		//}
		//us reference = EdgesData[fedge[referenceRotation]].reference;
		us reference = EdgesData[fedge[0]].reference;
		
		us rank1 = EdgesData[fedge[0]].rot;
		us rank2 = EdgesData[fedge[2]].rot;

		us faceOrientation[3];
		{
			if (fvert[0] < fvert[1])
				faceOrientation[0] = reference + 2;
			else
				faceOrientation[0] = rank1 + reference + 2;
			/////////////////////////
			if (fvert[1] < fvert[2])
				faceOrientation[1] = rank1 + reference;
			else
				faceOrientation[1] = rank2 + reference;
			/////////////////////////
			if (fvert[0] < fvert[2])
				faceOrientation[2] = reference;
			else
				faceOrientation[2] = rank2 + reference + 2;
		}
		
		ui eid[2];
		for (id = 0; id < 3; ++id)
		{
			ui esign[2] = { 1,1 };
			if (faceOrientation[id] & 1)
			{
				eid[0] = (fedge[id] << 1) + 1;
				eid[1] = fedge[id] << 1;
				esign[0] = 0;
			}
			else
			{
				eid[0] = fedge[id] << 1;
				eid[1] = (fedge[id] << 1) + 1;
			}
			if ((faceOrientation[id] >> 1) & 1)
			{
				esign[0] = esign[0] ? 0 : 1;
				esign[1] = esign[1] ? 0 : 1;
			}
			
			for (us k = 0; k < 2; ++k)
			{
				ui &index = eid[k];
				ui equationId = doubleFid + k;
				v4u &ec = edge_to_constraints[index];
				if (ec[0] == invalid)
				{
					ec[0] = equationId;
					ec[1] = esign[k];
				}
				else
				{
					ec[2] = equationId;
					ec[3] = esign[k];
				}
				initial[equationId] += EdgesData[fedge[id]].dir[eid[k] % 2] * (esign[k] ? 1 : -1);
			}
		}

		

		OH tempH = trimesh.find_halfedge(trimesh.vertex_handle(fvert[0]), trimesh.vertex_handle(fvert[1]));
		OH tfh = tempH;
		do
		{
			OH adjH = trimesh.opposite_halfedge_handle(tfh);
			OF adjF = trimesh.face_handle(adjH);
			if (!adjF.is_valid()) continue;
			if (trimesh.property(visited, adjF)) continue;
			trimesh.property(visited, adjF) = true;
			//EdgesData[trimesh.edge_handle(tfh).idx()].sign = true;
			auto &adjEdge = FacesData[adjF.idx()].edge;
			auto &adjVert = FacesData[adjF.idx()].vert;
			for (id = 0; id < 3; ++id)
			{
				adjEdge[id] = trimesh.edge_handle(adjH).idx();
				adjVert[id] = trimesh.from_vertex_handle(adjH).idx();
				adjH = trimesh.next_halfedge_handle(adjH);
			}
			/*id = 0;
			for (auto &tfhfh : trimesh.fh_range(adjF))
			{
				adjEdge[id] = trimesh.edge_handle(tfhfh).idx();
				adjVert[id] = trimesh.from_vertex_handle(tfhfh).idx();
				++id;
			}*/
			tree.push(adjF);
			tfh = trimesh.next_halfedge_handle(tfh);
		} while (tfh != tempH);
		tree.pop();
	}
	trimesh.remove_property(visited);

	std::vector<triple<ui, ui, ss>> arcs;
	std::vector<ui> arcId;
	for (id = 0; id < edge_to_constraints.size(); ++id)
	{
		v4u ecid = edge_to_constraints[id];
		if (ecid[0] == invalid || ecid[2] == invalid)
			continue;
		if (ecid[1] + ecid[3] != 1)
			continue;
		ui v0 = ecid[0];
		ui v1 = ecid[2];
		if (!ecid[1])
			std::swap(v0, v1);
		arcs.emplace_back(v0, v1, EdgesData[id >> 1].dir[id % 2]);
		arcId.push_back(id);
	}

	int supply = 0;
	int demand = 0;
	int initial_size = initial.size();
	for (int i = 0; i < initial_size; ++i) {
		int init_val = initial[i];
		if (init_val > 0) {
			arcs.emplace_back(invalid, i, initial[i]);
			supply += init_val;
		}
		else if (init_val < 0) {
			arcs.emplace_back(i, initial_size, -init_val);
			demand -= init_val;
		}
	}
	using namespace qflow;
	std::unique_ptr<MaxFlowHelper> solver = nullptr;
	solver = std::make_unique<NetworkSimplexFlowHelper>();
	solver->resize(initial.size() + 2, arcId.size());

	std::set<int> ids;
	for (int i = 0; i < arcs.size(); ++i) {
		int v1 = arcs[i].first + 1;
		int v2 = arcs[i].second + 1;
		int c = arcs[i].second;
		if (v1 == 0 || v2 == initial.size() + 1) {
			solver->addEdge(v1, v2, c, 0, -1);
		}
		else {
			if (arcId[i] > 0)
				solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
					std::max(0, -c + edge_capacity), arcId[i] - 1);
			/*else {
				if (c > 0)
					solver->addEdge(v1, v2, std::max(0, c - 1),
						std::max(0, -c + edge_capacity), -1 - arc_ids[i]);
				else
					solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
						std::max(0, -c - 1), -1 - arc_ids[i]);
			}*/
		}
	}
	int flow_count = solver->compute();

	//solver->applyTo(EdgeDiff);
#endif
}
#pragma endregion