#include <QMouseEvent>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtCore>
#include <QUrl>

#include "InteractiveViewerWidget.h"
#include "..\src\Toolbox\dprint.h"

InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
	:MeshViewerWidget(parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
}

InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
:MeshViewerWidget(_fmt, _parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
}

InteractiveViewerWidget::~InteractiveViewerWidget()
{
	if(kdTree) delete kdTree;
}

void InteractiveViewerWidget::setMouseMode(int mm)
{
	if(mouse_mode_ != T2_MODE)
	{
		mouse_mode_ = mm;
		if( TRANS != mouse_mode_ )
		{ buildIndex(); }
		emit setMouseMode_signal(mm);
	}
}

void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mousePressEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE)
		{
			pick_point( _event->x(), _event->y() );
			if(mouse_mode_ == VERTEXPICK)
			{
				pick_vertex( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == FACEPICK)
			{
				pick_face( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == EDGEPICK)
			{
				pick_edge( _event->x(), _event->y() );
			}
			else if(mouse_mode_ == POINTPICK)
			{
			}
			else if( mouse_mode_ == MOVE )
			{
				pick_vertex( _event->x(), _event->y() );//set the selected handle
			}
			else if(mouse_mode_ == EDGECOLLAPSE)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( mesh.edge_handle(desired_edge), 0 );
					OpenMesh::Vec3d from_p = mesh.point(mesh.from_vertex_handle(heh));
					OpenMesh::Vec3d to_p = mesh.point(mesh.to_vertex_handle(heh));
					OpenMesh::Vec3d sp(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
					bool collapse_ok = true;
					if( (sp-from_p).sqrnorm() > (to_p-sp).sqrnorm() )
					{
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					else
					{
						heh = mesh.opposite_halfedge_handle(heh);
						if( mesh.is_collapse_ok(heh) )
						{
							mesh.collapse(heh);
						}
						else
						{
							collapse_ok = false;
							printf("[%d] Collapse Not OK!\n", desired_edge);
						}
					}
					if(collapse_ok)
					{
						mesh.garbage_collection();
						buildIndex();
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGEFLIP)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					if( is_flip_ok_openmesh(eh, mesh))
					{
						flip_openmesh(eh, mesh);
						if( mesh_vector.size() - 1 > mesh_vector_index )
						{
							mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
						}
						mesh_vector.push_back( mesh ); mesh_vector_index += 1;
						emit set_edit_undo_enable_viewer_signal( true );
						emit set_edit_redo_enable_viewer_signal( false );
					}
					else
					{
						printf("[%d] Flip Not OK!\n", desired_edge);
					}
					clearSelectedData();
				}
			}
			else if (mouse_mode_ == EDGESPLIT)
			{
				int desired_edge = find_edge_using_selected_point();
				if(desired_edge >= 0) 
				{
					Mesh::EdgeHandle eh = mesh.edge_handle(desired_edge);
					Mesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );
					Mesh::HalfedgeHandle heh_ = mesh.halfedge_handle( eh, 1 );
					Mesh::VertexHandle vh0 = mesh.to_vertex_handle(heh);
					Mesh::VertexHandle vh1 = mesh.to_vertex_handle(heh_);
					OpenMesh::Vec3d s = mesh.point( vh1 );
					OpenMesh::Vec3d e = mesh.point( vh0 );
					Mesh::VertexHandle vh = mesh.add_vertex( (s + e)*0.5 );
					std::vector<Mesh::VertexHandle> one_face(3);
					if(mesh.is_boundary(eh))
					{
						if(Mesh::InvalidFaceHandle != mesh.face_handle(heh))
						{
							Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						}
						else
						{
							Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
							mesh.delete_edge(eh, false); mesh.garbage_collection();
							one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
							one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
						}
					}
					else
					{
						Mesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
						Mesh::VertexHandle vh3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh_));
						mesh.delete_edge(eh, false); mesh.garbage_collection();
						one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh.add_face(one_face);
						one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh.add_face(one_face);
					}

					mesh.update_normals();
					buildIndex();
					clearSelectedData();

					if( mesh_vector.size() - 1 > mesh_vector_index )
					{
						mesh_vector.erase( mesh_vector.begin() + mesh_vector_index + 1, mesh_vector.end() );
					}
					mesh_vector.push_back( mesh ); mesh_vector_index += 1;
					emit set_edit_undo_enable_viewer_signal( true );
					emit set_edit_redo_enable_viewer_signal( false );
				}
			}
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if( mouse_mode_ != T2_MODE)
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				updateGL();
			}
		}
		else
		{
			
		}
		
	}
}

void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE )
		{
			if( mouse_mode_ == MOVE )
			{
				move_point_based_lastVertex( _event->x(), _event->y() );
				Mesh::Point P(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
				mesh.set_point( mesh.vertex_handle(lastestVertex), P );
				selectedVertex.clear();
				updateGL();
			}
		}
		else
		{
		}
	}
	
}

void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
{
	if(mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
	{
		MeshViewerWidget::wheelEvent(_event);
	}
}

void InteractiveViewerWidget::dragEnterEvent(QDragEnterEvent* event)
{
	if( event->mimeData()->hasFormat("text/uri-list") )
	{
		event->acceptProposedAction();
	}
}

void InteractiveViewerWidget::dropEvent(QDropEvent* event)
{
	QList<QUrl> urls = event->mimeData()->urls();
	if( urls.isEmpty() )
		return;
	QString fileName = urls.first().toLocalFile();
	if (fileName.isEmpty())
		return;

	if( fileName.endsWith(".off") || fileName.endsWith(".obj") || fileName.endsWith(".stl") || fileName.endsWith(".ply"))
	{
		if( openMesh(fileName.toLocal8Bit()))
		{
			emit(loadMeshOK(true,fileName));
			setDrawMode(FLAT_POINTS);
			setMouseMode(TRANS);
		}
		else
		{
			emit(loadMeshOK(false,"No Mesh"));
		}
	}
}

void InteractiveViewerWidget::pick_vertex(int x,int y)
{
	int r = find_vertex_using_selected_point();
	lastestVertex = r;
	//printf("Select Vertex : %d\n", r);
	dprint("Select Vertex:", r, mesh.voh_begin(mesh.vertex_handle(r)).handle().idx() / 2);
	dprint("Vertex Position:", mesh.point(mesh.vertex_handle(r)));
	if (lg)
	{/*
		dprint("pl info:", lg->m4.sing_flag[r], lg->pls[lg->InfoOnMesh[lg->m4.verticemap[r]].plid].size(),
			lg->m4.sing_flag[r] ? -1 : lg->pls[lg->InfoOnMesh[lg->m4.verticemap[r] + 1].plid].size());*/
		//dprint("vertex layer:", lg->m4.verticemap[r]);
		dprint("vertex layer:", r * 4);
		if (lg->uv_para[0].size()!=0 && lg->vertexidmap.size() != 0)
			dprint("uv para:", lg->uv_para[0](lg->vertexidmap[r]));
		if (!lg->cset.vertex_bound_index.empty())
			dprint("vertex bound index:", lg->cset.vertex_bound_index[r].first, lg->cset.vertex_bound_index[r].second);
	}
	//dprint("uv para:", lg->uv_para[0](lg->idmap[r]), lg->uv_para[1](lg->idmap[r]));

	std::vector<int>::iterator it;
	if( (it = std::find(selectedVertex.begin(),selectedVertex.end(), r)) == selectedVertex.end() )
	{
		selectedVertex.push_back(r);
	}
	else
	{
		selectedVertex.erase(it);
	}

	updateGL();
}
void InteractiveViewerWidget::pick_face(int x,int y)
{
	int desiredFace = find_face_using_selected_point();
	if(desiredFace < 0) return;
	lastestFace = desiredFace;
	printf("Select Face : %d\n", desiredFace);
	std::vector<int>::iterator it;
	if( (it = std::find(selectedFace.begin(),selectedFace.end(),desiredFace)) == selectedFace.end() )
	{
		selectedFace.push_back(desiredFace);
	}
	else
	{
		selectedFace.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_edge(int x,int y)
{
	int desiredEdge = find_edge_using_selected_point();
	if(desiredEdge < 0) return;
	lastestEdge = desiredEdge;
	dprint("Select Edge:", lastestEdge);

	std::vector<int>::iterator it;
	if( (it = std::find(selectedEdge.begin(),selectedEdge.end(),desiredEdge)) == selectedEdge.end() )
	{
		selectedEdge.push_back(desiredEdge);
	}
	else
	{
		selectedEdge.erase(it);
	}
	updateGL();
}
void InteractiveViewerWidget::pick_point(int x,int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = double(x);
	GLdouble winY = double( height() - y );
	GLfloat winZ = 0.0;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
{
	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
	{
		return;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	OpenMesh::Vec3d p = mesh.point(mesh.vertex_handle(lastestVertex));
	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

int InteractiveViewerWidget::find_vertex_using_selected_point()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	return nnIdx[0];
}

int InteractiveViewerWidget::find_face_using_selected_point()
{
	int rv = find_vertex_using_selected_point();
	Mesh::VertexFaceIter vf_it = mesh.vf_iter( mesh.vertex_handle(rv) );
	int desiredFace = -1; //double minLen = 10*radius();
	std::vector<OpenMesh::Vec3d> tri_p(3); int tri_count = 0;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for( vf_it; vf_it; ++vf_it )
	{
		tri_count = 0;
		for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle()); fv_it; ++fv_it)
		{
			tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
		}
		if( check_in_triangle_face(tri_p, resultP) )
		{
			desiredFace = vf_it.handle().idx(); break;
		}
	}
	if(desiredFace < 0)
	{
		for(Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			tri_count = 0;
			for(Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it.handle()); fv_it; ++fv_it)
			{
				tri_p[tri_count] = mesh.point(fv_it); ++tri_count;
			}
			if( check_in_triangle_face(tri_p, resultP) )
			{
				desiredFace = f_it.handle().idx(); break;
			}
		}
	}

	return  desiredFace;
}

int InteractiveViewerWidget::find_edge_using_selected_point()
{
	int desiredFace = find_face_using_selected_point(); if(desiredFace < 0) return -1;
	Mesh::FaceHandle fh = mesh.face_handle(desiredFace);
	double min_len= 1e30; int desiredEdge = -1;
	Mesh::Point resultP(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	for(Mesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh); fhe_it; ++fhe_it)
	{
		OpenMesh::Vec3d s = mesh.point( mesh.from_vertex_handle(fhe_it) );
		OpenMesh::Vec3d e = mesh.point( mesh.to_vertex_handle(fhe_it) );
		double dis = OpenMesh::cross(resultP - s, resultP - e).norm() / (s - e).norm();
		if(dis < min_len){ min_len = dis; desiredEdge = mesh.edge_handle(fhe_it.handle()).idx(); }
	}
	
	return desiredEdge;
}

void InteractiveViewerWidget::buildIndex()
{
	if(mesh.n_vertices() == 0)
		return;

	Mesh::VertexIter v_it(mesh.vertices_begin());
	Mesh::VertexIter v_end(mesh.vertices_end());
	Mesh::Point p;
	unsigned nv = mesh.n_vertices();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	int count = 0;
	for(; v_it != v_end; ++v_it)
	{
		p = mesh.point(v_it);
		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
		++count;
	}

	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
}

//with the first mesh
void InteractiveViewerWidget::draw_interactive_portion(int drawmode)
{
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	
	emit draw_from_out_signal();

	{
		//draw select vertex, face, edge.
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glPointSize(1);

		switch(mouse_mode_)
		{
		case POINTPICK:
			draw_selected_point();
			break;
		case VERTEXPICK:
			draw_selected_vertex();
			break;
		case FACEPICK:
			draw_selected_face();
			break;
		case EDGEPICK:
			draw_selected_edge();
			break;
		default:
			draw_selected_vertex();
			draw_selected_face();
			draw_selected_edge();
			//showIsotropicMesh();
			break;
		}
	}
	draw_field();
	if (if_draw_energy) draw_energy();
	if (if_draw_submesh) draw_submesh();
	//draw_plane();
	draw_planeloop();
	if(draw_new_mesh)
	{
		draw_scene_mesh(drawmode);
	}
}

//with the second mesh
void InteractiveViewerWidget::draw_interactive_portion_mesh2()
{
	return;
}

void InteractiveViewerWidget::draw_selected_point()
{
	glColor3f(1.0, 0.5, 0.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3d(selectedPoint[0],selectedPoint[1],selectedPoint[2]);
	glEnd();
	glPointSize(1);
}

void InteractiveViewerWidget::draw_selected_vertex()
{
	if( selectedVertex.size() > 0 )
	{
		Mesh::Point p;
		glColor3f(1.0, 0.1, 0.1);
		glPointSize(12);
		glBegin(GL_POINTS);
		for(unsigned int i=0;i<selectedVertex.size();++i)
		{
			p = mesh.point( mesh.vertex_handle(selectedVertex[i]) );
			glVertex3dv(p.data());
		}
		glEnd();
		glPointSize(1);
	}
}

void InteractiveViewerWidget::draw_selected_face()
{
	if( selectedFace.size() > 0 )
	{
		glColor3f(1.0, 0.5, 1.0);
		Mesh::Point p;
		Mesh::ConstFaceVertexIter fv_it;
		Mesh::FaceHandle f_handle;
		for( unsigned int i=0; i<selectedFace.size(); ++i )
		{
			f_handle = mesh.face_handle(selectedFace[i]);
			fv_it = mesh.fv_iter(f_handle);
			glBegin(GL_POLYGON);
			for( fv_it; fv_it; ++fv_it )
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_selected_edge()
{
	if( selectedEdge.size() > 0)
	{
		glLineWidth(2);
		glColor3f(0.1, 0.1, 0.1);
		Mesh::Point p1; Mesh::Point p2;
		Mesh::EdgeHandle e_handle;
		Mesh::HalfedgeHandle he_handle;
		for(unsigned int i=0;i<selectedEdge.size();++i)
		{
			e_handle = mesh.edge_handle(selectedEdge[i]);
			he_handle = mesh.halfedge_handle( e_handle, 0 );
			p1 = mesh.point( mesh.from_vertex_handle( he_handle ) );
			p2 = mesh.point( mesh.to_vertex_handle( he_handle ) );
			glBegin(GL_LINES);
			glVertex3dv( p1.data() );
			glVertex3dv( p2.data() );
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_field()
{
	if (if_draw_field)
	{
		glLineWidth(1);
		glColor3d(0.9, 0.1, 0.1);
		glBegin(GL_LINES);
		for (int i = 0; i < crossfield.cols(); i += 4)
		{
			Eigen::Vector3d dd = (crossfield.col(i) + crossfield.col(i + 2)) * 0.5;
#if 0
			glVertex3dv(dd.data());
			glVertex3dv(crossfield.col(i).data());
			glVertex3dv(dd.data());
			glVertex3dv(crossfield.col(i + 1).data());
#else
			glVertex3dv(dd.data()); glVertex3dv(crossfield.col(i).data());
			glVertex3dv(dd.data()); glVertex3dv(crossfield.col(i + 1).data());
			glVertex3dv(dd.data()); glVertex3dv(crossfield.col(i + 2).data());
			glVertex3dv(dd.data()); glVertex3dv(crossfield.col(i + 3).data());
#endif
		}

		glLineWidth(5);
		glColor3d(0.1, 0.1, 0.1);
		glVertex3dv(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(2 * 75815))).data());
		glVertex3dv(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(2 * 75815))).data());
		glEnd();

		auto& sing = lg->cf->getSingularity();
		glPointSize(6);
		glColor3d(0.1, 0.1, 0.1);
		glBegin(GL_POINTS);
		for (auto& s : sing)
		{
			glVertex3dv(mesh.point(mesh.vertex_handle(s)).data());
		}
		glEnd();


		glPointSize(12);
		glBegin(GL_POINTS);
		glColor3d(0.7, 0.1, 0.0);
		//glVertex3dv(mesh.point(mesh.vertex_handle(2836)).data());
		/*glColor3d(0.0, 0.7, 0.60);
		glVertex3dv(mesh.point(mesh.vertex_handle(37644)).data());
		glVertex3dv(mesh.point(mesh.vertex_handle(34504)).data());
		glVertex3dv(mesh.point(mesh.vertex_handle(37606)).data());
		glVertex3dv(mesh.point(mesh.vertex_handle(27276)).data());
		glVertex3dv(mesh.point(mesh.vertex_handle(40641)).data());
		glVertex3dv(mesh.point(mesh.vertex_handle(27270)).data());*/
		glEnd();
	}
}

void InteractiveViewerWidget::draw_scene(int drawmode)
{
	if (!mesh.n_vertices()) { return; }
	draw_interactive_portion_mesh2();
	draw_interactive_portion(drawmode);

	if( !draw_new_mesh )
	{
		MeshViewerWidget::draw_scene(drawmode);
	}
}

void InteractiveViewerWidget::render_text_slot(OpenMesh::Vec3d pos, QString str)
{
	/*GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	render_text(x,y,str);*/
	render_text(pos[0],pos[1],pos[2],str);
}

#include "../Algorithm/LoopGen.h"
void InteractiveViewerWidget::showField()
{
	if (!if_has_field)
	{
		if_has_field = true;
		lg = new LoopGen::LoopGen(mesh);
		lg->SetModelName(file_name);
		lg->InitializeField();
		crossfield = lg->cf->getCrossField();
		avgLen = 0.2 * calc_mesh_ave_edge_length(&mesh);
		for (auto tf : mesh.faces())
		{
			OpenMesh::Vec3d c = mesh.calc_centroid(tf);
			int i = tf.idx() * 4;
			Eigen::Vector3d vc(c[0], c[1], c[2]);
			Eigen::Vector3d temp = crossfield.col(i + 1);
#if 0
			crossfield.col(i) = vc + crossfield.col(i) * avgLen;
			crossfield.col(i + 1) = vc + crossfield.col(i + 2) * avgLen;
			crossfield.col(i + 2) = vc + temp * avgLen;
			crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
#else
			crossfield.col(i) = vc + crossfield.col(i) * avgLen;
			crossfield.col(i + 1) = vc + crossfield.col(i + 1) * avgLen;
			crossfield.col(i + 2) = vc + crossfield.col(i + 2) * avgLen;
			crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
#endif
		}
	}
	if_draw_field = !if_draw_field;
	setDrawMode(InteractiveViewerWidget::SOLID_FLAT);
	setMouseMode(InteractiveViewerWidget::TRANS);
}

void InteractiveViewerWidget::showLoop()
{
	setDrawMode(InteractiveViewerWidget::SOLID_FLAT);
	setMouseMode(InteractiveViewerWidget::TRANS);
	timeRecorder tr;
	if (!if_has_field)
	{
		showField();
	}
	if (!loop_gen_init)
	{
		loop_gen_init = true;

#if 0
		lg->m4.set_base(&mesh, lg->cf); lg->m4.update(); lg->m4.set_weight();
#else
#if 0
		lg->InitializePQ();
		lg->OptimizeLoop();
#else
		lg->ReLoop();
#endif
#endif
		/*std::ifstream file_reader;
		file_reader.open("..//resource//energy//vase.energy", std::ios::in);
		char line[1024] = { 0 };
		int row = 0; lg->eov.resize(mesh.n_vertices());
		while (file_reader.getline(line, sizeof(line)))
		{
			std::stringstream num(line);
			num >> lg->eov[row];
			++row;
}
		file_reader.close();*/

	}
	tr.out("all time:");

	/*if (selectedVertex.empty())
		return;
	selectedVertex = { selectedVertex.back() };
	plane_loop[0] = lg->InfoOnMesh[selectedVertex.back() * 2].pl;
	plane_loop[1] = lg->InfoOnMesh[selectedVertex.back() * 2 + 1].pl;*/
#if 1
#if 1
	if (selectedVertex.empty())
		return;
#else
	selectedVertex.push_back(601);//2836
#endif
	selectedVertex = { selectedVertex.back() };
	selectedEdge.clear();
	Eigen::VectorXd xyz[3];
	VertexHandle vnow = mesh.vertex_handle(selectedVertex.back());
	//LoopGen::VertexLayer* vl = &lg->m4.verticelayers[lg->m4.verticemap[vnow.idx()]];
	LoopGen::VertexLayer* vl = &lg->m4.verticelayers[vnow.idx() * 4];
	if (!lg->m4.sing_flag[vnow.idx()])
	{
		for (int i = 0; i < 2; ++i, ++vl)
		{
			//if (i == 0)
				//continue;
			if (lg->FieldAligned_PlanarLoop(vl, loop))
			{
				for (int j = 0; j < loop.size() - 1; ++j)
				{
					selectedEdge.push_back(mesh.find_halfedge(loop[j]->v, loop[j + 1]->v).idx() / 2);
				}
				lg->GetPositionFromLoop(loop, xyz);
				if_draw_plane = true;
				plane_loop[i].clear();
#if 1
				//dprint(lg->RefineLoopByPlanarity(loop, plane_loop[0], 0));
				dprint(lg->RefineLoopByPlanarity(loop, plane_loop[i]));
#else
				lg->GetPositionFromLoop(loop, xyz);
				double plane[4];
				//LoopGen::LeastSquarePlane(xyz, plane0);
				lg->LeastSquarePlane(xyz, plane0);
#endif
			}
		}
	}
#endif
	//plane_loop[0] = lg->pls[79180];
	//plane_loop[1] = lg->pls[80522];

	updateGL();
}

void InteractiveViewerWidget::showAnisotropicMesh()
{
	if_draw_energy = !if_draw_energy;
	setDrawMode(InteractiveViewerWidget::SOLID_FLAT);
	setMouseMode(InteractiveViewerWidget::TRANS);
}

#include "..\src\Toolbox\filesOperator.h"
void InteractiveViewerWidget::showDebugTest()
{
	if_draw_submesh = !if_draw_submesh;
	setDrawMode(InteractiveViewerWidget::SOLID_FLAT);
	setMouseMode(InteractiveViewerWidget::TRANS);
}

void InteractiveViewerWidget::draw_energy()
{
	double max_e = 0.25;
	double step_e =1.0/max_e;
	if (loop_gen_init)
	{
#if 0
		glLineWidth(10);
		glBegin(GL_LINES);
		//for(int i=0;i<lg->similarity_energy.size()/8;++i)
		for(int i=0;i<mesh.n_edges();++i)
		{
			double t = std::min(lg->similarity_energy[8 * i], lg->similarity_energy[8 * i + 1]);
			if (i == 75717)
			{
				//dprint("ee", 8*i, lg->similarity_energy[8 * i], lg->similarity_energy[8 * i + 1]);
			}
			if (t > max_e)
				glColor3d(0, 1, 0);
			else
			{
				double c = t * step_e;
				glColor3d(c, 0, 1 - c);
			}
			auto h = mesh.halfedge_handle(2 * i);
			glVertex3dv(mesh.point(mesh.from_vertex_handle(h)).data());
			glVertex3dv(mesh.point(mesh.to_vertex_handle(h)).data());
		}
#else
		glBegin(GL_POINTS);
		/*for (int i = 0; i < lg->m4.verticelayers.size(); ++i)
		{
			if (lg->pls[lg->InfoOnMesh[i].plid].empty())
			{
				glColor3d(0, 0, 0);
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}*/
		for (int i = 0; i < mesh.n_vertices(); ++i)
		{
			//int ee = lg->m4.verticemap[i];
			int ee = i * 4;
			//dprint(mesh.n_vertices(), lg->m4.sing_flag[ee], lg->InfoOnMesh[ee].plid, lg->pls.size());
			if (lg->m4.sing_flag[i])
				continue;
			
			//if (lg->InfoOnMesh[lg->m4.verticemap[i]].energy < 0.25||
			//	lg->InfoOnMesh[lg->m4.verticemap[i] + 1].energy < 0.25)
			if(lg->InfoOnMesh[i*4].energy<0.25||lg->InfoOnMesh[i*4+1].energy<0.25)
				glColor3d(0, 1, 0);
			else
				glColor3d(1, 0, 0);
			glVertex3dv(mesh.point(mesh.vertex_handle(i)).data());
		}
#endif
		glEnd();
	}

}

void InteractiveViewerWidget::draw_submesh()
{
	//画某个点
#if 1
	glColor3d(1, 0, 0);
	glPointSize(15);
	glBegin(GL_POINTS);
	glVertex3dv(mesh.point(mesh.vertex_handle(143946 /4)).data());
	glEnd();
#endif

	//画某条线
#if 1
	glLineWidth(3);
	glBegin(GL_LINES);
	glColor3d(0.8, 0.2, 0.1);

	for (auto &path : lg->all_vertice_path)
	{
		if (path.size() < 2)
			continue;

		if (path.front().hl)
			glVertex3dv(path.front().point(lg->m4).data());
		else
			glVertex3dv(mesh.point(mesh.vertex_handle(path.front().c)).data());
		OpenMesh::Vec3d p;
		for (int i = 1; i < path.size() - 1; ++i)
		{
			if (path[i].hl)
				p = path[i].point(lg->m4);
			else
				p = mesh.point(mesh.vertex_handle(path[i].c));
			glVertex3dv(p.data());
			glVertex3dv(p.data());
		}
		if (path.back().hl)
			glVertex3dv(path.back().point(lg->m4).data());
		else
			glVertex3dv(mesh.point(mesh.vertex_handle(path.back().c)).data());
	}
	glEnd();
#endif

	//画u参数
#if 0
	double vmi = 1.0e12;
	double vma = -1.0e12;
	for (int i = 0; i < lg->old_vert_flag.size(); ++i)
	{
		if (lg->old_vert_flag[i])
		{
			double c = lg->uv_para[0](lg->vertexidmap[i]);
			vmi = std::min(c, vmi);
			vma = std::max(c, vma);
		}
	}
	vma = 1.0 / (vma - vmi);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (int i = 0; i < lg->old_vert_flag.size(); ++i)
	{
		if (lg->old_vert_flag[i])
		{
			double c = lg->uv_para[0](lg->vertexidmap[i]);
			c = (c - vmi)*vma;
			glColor3d(c, c, c);
			glVertex3dv(mesh.point(mesh.vertex_handle(i)).data());
		}
	}
	glEnd();
#endif

	//画v参数
#if 0
	double vmi = 1.0e12;
	double vma = -1.0e12;
	for (int i = 0; i < lg->old_vert_flag.size(); ++i)
	{
		if (lg->old_vert_flag[i])
		{
			double c = lg->uv_para[1](lg->vertexidmap[i]);
			vmi = std::min(c, vmi);
			vma = std::max(c, vma);
		}
	}
	vma = 1.0 / (vma - vmi);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (int i = 0; i < lg->old_vert_flag.size(); ++i)
	{
		if (lg->old_vert_flag[i])
		{
			double c = lg->uv_para[1](lg->vertexidmap[i]);
			c = (c - vmi)*vma;
			glColor3d(c, c, c);
			glVertex3dv(mesh.point(mesh.vertex_handle(i)).data());
		}
	}
	glEnd();
#endif

	//画方向
#if 1
	double avgl = 0.2*calc_mesh_ave_edge_length(&mesh);
	glColor3d(1, 0, 0);
	glBegin(GL_LINES);
	for (int i = 0; i < lg->new_face_flag.size(); ++i)
	{
		if (lg->new_face_flag[i])
		{
			auto c = mesh.calc_face_centroid(mesh.face_handle(i));
			Eigen::Vector3d vc(c[0], c[1], c[2]);
			Eigen::Vector3d topoint = vc + lg->xaxis.col(i) * avgl;
			glVertex3dv(vc.data());
			glVertex3dv(topoint.data());
		}
	}
	glEnd();
#endif

	//画区域顶点
#if 1
	glColor3d(1, 1, 1);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (int i = 0; i < lg->new_vert_flag.size(); ++i)
	{
		if (lg->new_vert_flag[i])
		{
			glVertex3dv(mesh.point(mesh.vertex_handle(i)).data());
		}
	}
	glEnd();
#endif

	//画区域的面
#if 1
	glColor3d(0, 1, 0);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < lg->new_face_flag.size(); ++i)
	{
		if (lg->new_face_flag[i])
		{
			auto itr = mesh.fv_begin(mesh.face_handle(i));
			glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
			glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
			glVertex3dv(mesh.point(itr.handle()).data());
		}
		}
	glEnd();
#endif

	//画边界搜索路径
#if 1
	glLineWidth(3);
	glBegin(GL_LINES);
	glColor3d(0.8, 0.2, 0.1);
	int nv = mesh.n_vertices();
	for (auto &path : lg->all_vertice_path)
	{
		if (path.size() < 2)
			continue;
		for (int i = 0; i < path.size() - 1; ++i)
		{
			glVertex3dv(path[i].point(lg->m4).data());
			glVertex3dv(path[i + 1].point(lg->m4).data());
		}
		}
	glEnd();
#endif

	if (loop_gen_init)
	{
		//画区域顶点
#if 0
		glColor3d(1, 1, 1);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->old_vert_flag.size(); ++i)
		{
			if (lg->old_vert_flag[i])
			{
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif

		//画区域的面
#if 0
		glColor3d(0, 1, 0);
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < lg->old_face_flag.size(); ++i)
		{
			if (lg->old_face_flag[i])
			{
				auto itr = mesh.fv_begin(mesh.face_handle(i / 4));
				glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
				glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
				glVertex3dv(mesh.point(itr.handle()).data());
			}
		}
		glEnd();
#endif

		//画新搜索的顶点
#if 0
		glColor3d(1, 0, 0);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->new_vert_flag.size(); ++i)
		{
			if (lg->new_vert_flag[i])
			{
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif

		//画新搜索的面
#if 0
		glColor3d(0, 0, 1);
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < lg->new_face_flag.size(); ++i)
		{
			if (lg->new_face_flag[i])
			{
				auto itr = mesh.fv_begin(mesh.face_handle(i / 4));
				glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
				glVertex3dv(mesh.point(itr.handle()).data()); ++itr;
				glVertex3dv(mesh.point(itr.handle()).data());
			}
		}
		glEnd();
#endif

		//画所有柱体区域边界
#if 1
		glLineWidth(6);
		glBegin(GL_LINES);
		for (auto &cy : lg->cset.cylinders)
		{
			for (int i = 0; i < 2; ++i)
			{
				glColor3d(i, i, 1 - i);
				if (cy.bounds[i].empty())
					continue;
				for (auto hl : cy.bounds[i])
				{
					glVertex3dv(mesh.point(mesh.to_vertex_handle(hl->h)).data());
					glVertex3dv(mesh.point(mesh.from_vertex_handle(hl->h)).data());
				}
			}
		}
		glEnd();
#endif

		//画所有柱体顶点
#if 0
		glColor3d(0, 1, 1);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (auto &cy : lg->cset.cylinders)
		{
			for (auto vl : cy.vertices)
			{
				glVertex3dv(mesh.point(vl->v).data());
			}
		}
		glEnd();
#endif
		
		//画所有柱体的面
#if 0
		glColor3d(0, 0.8, 0);
		glBegin(GL_TRIANGLES);
		for (auto &cy : lg->cset.cylinders)
		{
			for (auto fl : cy.faces)
			{
				auto hl = fl->hl;
				glVertex3dv(mesh.point(lg->m4.verticelayers[hl->to].v).data());
				glVertex3dv(mesh.point(lg->m4.verticelayers[hl->next->to].v).data());
				glVertex3dv(mesh.point(lg->m4.verticelayers[hl->from].v).data());
			}
		}
		glEnd();
#endif

		//画所有柱体的cut
#if 0
		glColor3d(1, 0, 0);
		glBegin(GL_LINES);
		for (auto &cy : lg->cset.cylinders)
		{
			for (int i = 0; i < cy.cut.size() - 1; ++i)
			{
				glVertex3dv(mesh.point(cy.cut[i]->v).data());
				glVertex3dv(mesh.point(cy.cut[i + 1]->v).data());
			}
		}
		glEnd();
#endif

		//画所有柱体的seed
#if 0
		glColor3d(0, 1, 0);
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto &cy : lg->cset.cylinders)
		{
			glVertex3dv(mesh.point(cy.vertices.front()->v).data());
		}
		glEnd();
#endif

		//画所有柱体的u参数化
#if 0
		glPointSize(5);
		glBegin(GL_POINTS);
		for (auto &cy : lg->cset.cylinders)
		{
			for (auto vl : cy.vertices)
			{
				double c = cy.uv[0](cy.vidmap[vl->id]);
				c -= std::floor(c);
				glColor3d(c, c, c);
				glVertex3dv(mesh.point(vl->v).data());
			}
		}
		glEnd();
#endif

		//画所有柱体的v参数化
#if 0
		double vmi = 1.0e12;
		double vma = -1.0e12;
		for (auto &cy : lg->cset.cylinders)
		{
			for (auto vl : cy.vertices)
			{
				double c = cy.uv[1](cy.vidmap[vl->id]);
				vmi = std::min(vmi, c);
				vma = std::max(vma, c);
			}
		}
		vma = 1.0 / (vma - vmi);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (auto &cy : lg->cset.cylinders)
		{
			for (auto vl : cy.vertices)
			{
				double c = cy.uv[1](cy.vidmap[vl->id]);
				c = (c - vmi)*vma;
				glColor3d(c, c, c);
				glVertex3dv(mesh.point(vl->v).data());
			}
		}
		glEnd();
#endif

		//画所有柱体的边界搜索路径
#if 0
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3d(0.8, 0.2, 0.1);
		int nv = mesh.n_vertices();
		for (auto &path : lg->cset.all_path)
		{
			if (path.size() < 2)
				continue;
			for (int i = 0; i < path.size() - 1; ++i)
			{
				glVertex3dv(path[i].point(lg->m4).data());
				glVertex3dv(path[i + 1].point(lg->m4).data());
			}
		}
		glEnd();
#endif

		//画u参数
#if 0
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->old_vert_flag.size(); ++i)
		{
			if (lg->old_vert_flag[i] || lg->new_vert_flag[i])
			{
				double c = lg->uv_para[0](lg->vertexidmap[i]);
				c -= std::floor(c);
				glColor3d(c, c, c);
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif
		//画v参数
#if 0
		double vma = -1.0 * YYSS_FAIRLY_SMALL;
		double vmi = YYSS_INFINITE;
		for (int i = 0; i < lg->old_vert_flag.size(); ++i)
		{
			if (lg->old_vert_flag[i])
			{
				double c = lg->uv_para[1](lg->vertexidmap[i]);
				vma = std::max(vma, c);
				vmi = std::min(vmi, c);
			}
		}
		double step = 1.0 / (vma - vmi);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->old_vert_flag.size(); ++i)
		{
			if (lg->old_vert_flag[i])
			{
				double c = (lg->uv_para[1](lg->idmap[i]) - vmi) * step;
				glColor3d(c, c, c);
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif

		//画生长方向
#if 0
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->growDIR.size(); ++i)
		{
			if (lg->growDIR[i] == 0)
				glColor3d(1, 0, 0);
			else if (lg->growDIR[i] == 1)
				glColor3d(0, 1, 0);
			if (lg->growDIR[i] != -1)
			{
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif

		//画待加入的点
#if 0
		glColor3d(1, 0, 0);
		glPointSize(5);
		glBegin(GL_POINTS);
		/*for (int i = 0; i < lg->v_cache_flag.size(); ++i)
		{
			if (lg->v_cache_flag[i])
			{
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}*/
		for (auto vl : lg->v_cache)
		{
			glVertex3dv(mesh.point(vl->v).data());
		}
		glEnd();
#endif

		//画cut顶点
#if 0
		glColor3d(1, 1, 0);
		glPointSize(5);
		glBegin(GL_POINTS);
		for (int i = 0; i < lg->cut_vertex_flag.size(); ++i)
		{
			if (lg->cut_vertex_flag[i])
			{
				glVertex3dv(mesh.point(lg->m4.verticelayers[i].v).data());
			}
		}
		glEnd();
#endif

		//画新设的场
#if 0
		double avgl = 0.2*calc_mesh_ave_edge_length(&mesh);
		glColor3d(1, 0, 0);
		glBegin(GL_LINES);
		for (int i = 0; i < lg->new_face_flag.size(); ++i)
		{
			if (lg->new_face_flag[i])
			{
				auto &fl = lg->m4.facelayers[i];
				auto c = mesh.calc_face_centroid(fl.f);
				Eigen::Vector3d vc(c[0], c[1], c[2]);
				Eigen::Vector3d topoint = vc + lg->xaxis.col(fl.f.idx()) * avgl;
				glVertex3dv(vc.data());
				glVertex3dv(topoint.data());
			}
		}
		glEnd();
#endif

		//画loop
#if 0
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3d(0.8, 0.2, 0.1);
		int nv = mesh.n_vertices();
		for (int j = 0; j < lg->all_plane_loop.size(); ++j)
		{
			if (j % 5 != 0)
				continue;
			auto& pl = lg->all_plane_loop[j];
			int nn = pl.size();
			if (nn < 1)
				continue;
			auto poin = pl.front().point(lg->m4);
			glVertex3dv(poin.data());
			for (int i = 1; i < nn - 1; ++i)
			{
				poin = pl[i].point(lg->m4);
				glVertex3dv(poin.data());
				glVertex3dv(poin.data());
			}
			poin = pl.back().point(lg->m4);
			glVertex3dv(poin.data());
		}
		glEnd();
#endif

#if 1
		static bool pe = false;
		if (!pe)
		{
			pe = true;
			crossfield = lg->cf->getCrossField();
			avgLen = 0.2 * calc_mesh_ave_edge_length(&mesh);
			for (auto tf : mesh.faces())
			{
				OpenMesh::Vec3d c = mesh.calc_centroid(tf);
				int i = tf.idx() * 4;
				Eigen::Vector3d vc(c[0], c[1], c[2]);

				crossfield.col(i) = vc + crossfield.col(i) * avgLen;
				crossfield.col(i + 1) = vc + crossfield.col(i + 1) * avgLen;
				crossfield.col(i + 2) = vc + crossfield.col(i + 2) * avgLen;
				crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
			}
		}
#endif

	}
}

void InteractiveViewerWidget::draw_plane()
{
	if (if_draw_plane)
	{
		glColor3d(0.9, 0.9, 0.9);
		glBegin(GL_POLYGON);
		glVertex3d(2, 2, -(plane0[0] * 2 + plane0[1] * 2 + plane0[3]) / plane0[2]);
		glVertex3d(-2, 2, -(plane0[0] * -2 + plane0[1] * 2 + plane0[3]) / plane0[2]);
		glVertex3d(-2, -2, -(plane0[0] * -2 + plane0[1] * -2 + plane0[3]) / plane0[2]);
		glVertex3d(2, -2, -(plane0[0] * 2 + plane0[1] * -2 + plane0[3]) / plane0[2]);
		glEnd();

		glBegin(GL_POLYGON);
		glVertex3d(2, 2, -(plane1[0] * 2 + plane1[1] * 2 + plane1[3]) / plane1[2]);
		glVertex3d(-2, 2, -(plane1[0] * -2 + plane1[1] * 2 + plane1[3]) / plane1[2]);
		glVertex3d(-2, -2, -(plane1[0] * -2 + plane1[1] * -2 + plane1[3]) / plane1[2]);
		glVertex3d(2, -2, -(plane1[0] * 2 + plane1[1] * -2 + plane1[3]) / plane1[2]);
		glEnd();
	}
}

void InteractiveViewerWidget::draw_planeloop()
{
#if 1
	glColor3d(0.1, 0.8, 0.3);
	glBegin(GL_LINES);
	for (int i = 0; i < 2; ++i)
	{
		if (plane_loop[i].empty())
			continue;
		for (int j = 0; j < plane_loop[i].size() - 1; ++j)
		{
			auto& pl = plane_loop[i][j];
			auto& ps = plane_loop[i][j + 1];
			//glVertex3dv(((1 - pl.c) * mesh.point(mesh.to_vertex_handle(pl.h)) + pl.c * mesh.point(mesh.from_vertex_handle(pl.h))).data());
			//glVertex3dv(((1 - ps.c) * mesh.point(mesh.to_vertex_handle(ps.h)) + ps.c * mesh.point(mesh.from_vertex_handle(ps.h))).data());
			glVertex3dv(pl.point(lg->m4).data());
			glVertex3dv(ps.point(lg->m4).data());
		}
	}
	glEnd();
#else
	if (lg)
	{
		auto h = mesh.halfedge_handle(mesh.edges_begin().handle(), 0);
		if (lg->InfoOnMesh.empty())
			return;
		
		plane_loop[0] = lg->InfoOnMesh[mesh.from_vertex_handle(h).idx() * 2].pl;
		plane_loop[1] = lg->InfoOnMesh[mesh.to_vertex_handle(h).idx() * 2].pl;
		for (int i = 0; i < 2; ++i)
		{
			if (plane_loop[i].empty())
				continue;
			for (int j = 0; j < plane_loop[i].size() - 1; ++j)
			{
				auto& pl = plane_loop[i][j];
				auto& ps = plane_loop[i][j + 1];
				glVertex3dv(((1 - pl.c) * mesh.point(mesh.to_vertex_handle(pl.h)) + pl.c * mesh.point(mesh.from_vertex_handle(pl.h))).data());
				glVertex3dv(((1 - ps.c) * mesh.point(mesh.to_vertex_handle(ps.h)) + ps.c * mesh.point(mesh.from_vertex_handle(ps.h))).data());
			}
		}
		/*glColor3d(0.1, 0.8, 0.3);
		glBegin(GL_LINES);
		for (int i = 0; i < lg->loop0.cols(); ++i)
		{
			glVertex3dv(lg->loop0.col(i).data());
			glVertex3dv(lg->loop0.col((i + 1) % lg->loop0.cols()).data());
		}
		for (int i = 0; i < lg->loop1.cols(); ++i)
		{
			glVertex3dv(lg->loop1.col(i).data());
			glVertex3dv(lg->loop1.col((i + 1) % lg->loop1.cols()).data());
		}
		glEnd();*/

		/*glPointSize(15);
		glBegin(GL_POINTS);
		glColor3d(0.9, 0.9, 0.9);
		glVertex3dv(lg->fragment0.col(0).data());
		glVertex3dv(lg->fragment1.col(0).data());
		glColor3d(0.1, 0.1, 0.9);
		glVertex3dv(lg->fragment0.col(1).data());
		glVertex3dv(lg->fragment1.col(1).data());
		glColor3d(0.9, 0.1, 0.1);
		for (int i = 2; i < lg->fragment0.cols()-3; ++i)
		{
			glVertex3dv(lg->fragment0.col(i).data());
		}
		for (int i = 2; i < lg->fragment1.cols()-3; ++i)
		{
			glVertex3dv(lg->fragment1.col(i).data());
		}
		glEnd();*/
	}
#endif

}


