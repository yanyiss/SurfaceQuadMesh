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
	if(lg->simi_e.empty())
		dprint("Select Edge :", desiredEdge);// , "similarity energy:", std::min(lg->simi_e[desiredEdge * 2], lg->simi_e[desiredEdge * 2 + 1]));
	else
		dprint("Select Edge :", desiredEdge, "similarity energy:", std::min(lg->simi_e[desiredEdge * 2], lg->simi_e[desiredEdge * 2 + 1]));
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
	draw_energy();
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
			Eigen::Vector3d dd = (crossfield.col(i) + crossfield.col(i + 1)) * 0.5;
			glVertex3dv(dd.data());
			glVertex3dv(crossfield.col(i ).data());
		}
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

#include "../Algorithm/LoopToolbox.h"
void InteractiveViewerWidget::showField()
{
	if (!if_has_field)
	{
		if_has_field = true;
		lg = new LoopGen::LoopGen(mesh);
		lg->InitializeField();
		crossfield = lg->cf->getCrossField();
		avgLen = 0.2 * calc_mesh_ave_edge_length(&mesh);
		for (auto tf : mesh.faces())
		{
			OpenMesh::Vec3d c = mesh.calc_centroid(tf);
			int i = tf.idx() * 4;
			Eigen::Vector3d vc(c[0], c[1], c[2]);
			Eigen::Vector3d temp = crossfield.col(i + 1);
			crossfield.col(i) = vc + crossfield.col(i) * avgLen;
			crossfield.col(i + 1) = vc + crossfield.col(i + 2) * avgLen;
			crossfield.col(i + 2) = vc + temp * avgLen;
			crossfield.col(i + 3) = vc + crossfield.col(i + 3) * avgLen;
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
	if (!if_has_field)
	{
		showField();
	}
	if (!loop_gen_init)
	{
		loop_gen_init = true;
		lg->InitializeGraphWeight();
		lg->InitializePQ();

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
	selectedVertex.push_back(18442);
#endif
	selectedVertex = { selectedVertex.back() };
	selectedEdge.clear();
	Eigen::VectorXd xyz[3];
	if (lg->FieldAligned_PlanarLoop(mesh.vertex_handle(selectedVertex.back()), loop, 0))
	{
		for (int i = 0; i < loop.size() - 5; ++i)
		{
			selectedEdge.push_back(mesh.find_halfedge(loop[i], loop[i + 1]).idx() / 2);
		}
		lg->GetPositionFromLoop(loop, xyz);
		if_draw_plane = true;
		plane_loop[0].clear();
#if 1
		dprint(lg->RefineLoop(loop, plane_loop[0], 0));
#else
		lg->GetPositionFromLoop(loop, xyz);
		double plane[4];
		LoopGen::LeastSquarePlane(xyz, plane0);
#endif
	}
	if (lg->FieldAligned_PlanarLoop(mesh.vertex_handle(selectedVertex.back()), loop, 1))
	{
		for (int i = 0; i < loop.size() - 5; ++i)
		{
			selectedEdge.push_back(mesh.find_halfedge(loop[i], loop[i + 1]).idx() / 2);
		}
		lg->GetPositionFromLoop(loop, xyz);
		if_draw_plane = true;
		plane_loop[1].clear();
#if 1
		dprint(lg->RefineLoop(loop, plane_loop[1], 1));
#else
		lg->GetPositionFromLoop(loop, xyz);
		double plane[4];
		LoopGen::LeastSquarePlane(xyz, plane1);
#endif
	}
	/*lg->FieldAligned_PlanarLoop(mesh.vertex_handle(29054), loop, 0);
	for (int i = 0; i < loop.size() - 1; i += 2)
	{
		selectedEdge.push_back(mesh.find_halfedge(mesh.vertex_handle(loop[i]), mesh.vertex_handle(loop[i + 1])).idx() / 2);
	}*/
#endif


	updateGL();
}

void InteractiveViewerWidget::showAnisotropicMesh()
{
}

#include "..\src\Toolbox\filesOperator.h"
void InteractiveViewerWidget::showDebugTest()
{

}

void InteractiveViewerWidget::draw_energy()
{
#if 0
	double ie = DBL_MAX;
	double ae = 0;
	static bool ff = true;
	static double ll = 0;
	if (loop_gen_init && ff)
	{
		ff = false;
		for (auto e : lg->eov)
		{
			if (e > 10e10)
				continue;
			ie = std::min(ie, e);
			ae = std::max(ae, e);
		}
		ll = 1.0 / (ae - ie);
	}
	if (loop_gen_init)
	{
		ie = 0.00030153500000000002;
		ae = 0.22753200000000001;
		ll = 1.0 / (ae - ie);
		//double c = 0;
		/*glPointSize(5);
		for (auto v : mesh.vertices())
		{
			if (lg->eov[v.idx()] > 10e10)
			{
				glColor3d(0, 1, 0);
				glBegin(GL_POINTS);
				glVertex3dv(mesh.point(v).data());
				glEnd();
			}
			else
			{
				c = pow((lg->eov[v.idx()] - ie) * ll, 0.4);
				glColor3d(c, 0, 1 - c);
				glBegin(GL_POINTS);
				glVertex3dv(mesh.point(v).data());
				glEnd();
			}
		}*/
		for (auto f : mesh.faces())
		{
			double c = 0; VertexHandle v[3];
			auto fv = mesh.fv_begin(f); /*c += pow((lg->eov[fv->idx()] - ie) * ll, 0.4);*/ v[0] = fv.handle();
			++fv; /*c += pow((lg->eov[fv->idx()] - ie) * ll, 0.4);*/ v[1] = fv.handle();
			++fv; /*c += pow((lg->eov[fv->idx()] - ie) * ll, 0.4);*/ v[2] = fv.handle();
			for (int i = 0; i < 3; ++i)
			{
				if (lg->eov[v[i].idx()] > 10e10)
				{
					c = 10e11;
					break;
				}
				c += pow((lg->eov[v[i].idx()] - ie) * ll, 0.4);
			}
			if (c > 10e10)
			{
				glColor3d(0, 1, 0);
			}
			else
			{
				c *= 0.33333333333333333333;
				glColor3d(c, 0.2, 1 - c);
			}
			glBegin(GL_TRIANGLES);
			glVertex3dv(mesh.point(v[0]).data());
			glVertex3dv(mesh.point(v[2]).data());
			glVertex3dv(mesh.point(v[1]).data());
			glEnd();
		}
	}
#else
	double max_e = 2;
	double step_e =1.0/max_e;
	if (loop_gen_init)
	{
		glLineWidth(10);
		glBegin(GL_LINES);
		for(int i=0;i<lg->simi_e.size()/2;++i)
		{
			double t = std::min(lg->simi_e[2 * i], lg->simi_e[2 * i + 1]);
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
		glEnd();
	}

#endif
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
			glVertex3dv(((1 - pl.c) * mesh.point(mesh.to_vertex_handle(pl.h)) + pl.c * mesh.point(mesh.from_vertex_handle(pl.h))).data());
			glVertex3dv(((1 - ps.c) * mesh.point(mesh.to_vertex_handle(ps.h)) + ps.c * mesh.point(mesh.from_vertex_handle(ps.h))).data());
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


