#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent=0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize( int( rect.width()*0.15), rect.height() );
	}

private:
	QTabWidget* tabWidget;

signals:
	void print_info_signal();
	void show_crossfield_signal();
	void change_M4_signal();
	void change_Direction_signal();
	void test_matchings_signal();
	void cal_vertex_loop_signal();

private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;
	QPushButton* show_crossfield_button;
	QPushButton* change_M4_button;
	QPushButton* change_Direction_button;
	QPushButton* test_matchings;
	QPushButton* cal_vertex_loop_button;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};

#endif