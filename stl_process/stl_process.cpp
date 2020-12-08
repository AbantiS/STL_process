// stl_process.cpp : Defines the entry point for the console application.
//

// cad_i17.cpp : Defines the entry point for the console application.
//updated on 17.11.19
// updated on 25.11.19

#include "stdafx.h"
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <chrono> 
//#include "stdafx.h"
#include <sstream>
//#include "CSVFile.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <list>
#include <numeric>

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "spline.h"

//header for opengl pllotting starts
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\Dependencies\glew\glew.h"
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\Dependencies\freeglut\freeglut.h"
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\Dependencies\glfw\GLFW\glfw3.h"
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\Dependencies\glm\glm.hpp"
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\Dependencies\glm\gtc\matrix_transform.hpp"
#include<cmath> 

#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\stl_process\blockcad.h"
#include "D:\code_works\cad_i17_versions\v_18_32bits\stl_process\stl_process\interface_ogl3d.h"

using namespace std;
using namespace cadblock;
using namespace ogl_Point3D;
using namespace Eigen;
using Eigen::MatrixXd;
//using namespace std::chrono; 

// ine ~862  to line -- for visualization
struct color4opengl {
	float x, y, z;
};
vector<color4opengl> color_plot_mat, color_plot_mat2, color_plot_mat3, color_plot_mat3_l;
int plot_win_max_x, plot_win_min_x, plot_win_max_y, plot_win_min_y, plot_win_max_z, plot_win_min_z;
int plot_win2_max_x, plot_win2_min_x, plot_win2_max_y, plot_win2_min_y;
int plot_win3_max_x, plot_win3_min_x, plot_win3_max_y, plot_win3_min_y;
vector<cad3d::Point3D> target_plot, target_plot2, targetl1, targetl2;
vector<cad3d::Point3D> targetl31, targetl32, target_plot3;
vector<cad3d::Point3D> target_poly_3;

vector<color4opengl> color_plot_mat_poly3;
vector<color4opengl> color_plot_mat4, color_plot_mat4_l;
vector<cad3d::Point3D> targetl4_1, targetl4_2, target_plot4;

vector<cad3d::Point3D> targetl31_mod, targetl32_mod, target_plot3_mod;
cad3d blockplot;
int WindowID1, WindowID2, WindowID3, WindowID4, WindowID5, WindowID6, WindowID101, WindowID102, WindowID103;
void window_range_update(int x1min, int x1max, int y1min, int y1max, int x2min, int x2max, int y2min, int y2max) {
	plot_win_min_x = x1min;
	plot_win_max_x = x1max;
	plot_win_min_y = y1min;
	plot_win_max_y = y1max;
	plot_win2_min_x = x2min;
	plot_win2_max_x = x2max;
	plot_win2_min_y = y2min;
	plot_win2_max_y = y2max;
}
int init_opengl(int &argc, char* argv[])
{
	//updated 25.11.19
	//for window 1 
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	
	// 3 arguments all are 0.0 
	glClearColor(0.9, 0.9, 0.9, 1.0);

	return 0;
}
void display_targetplot_point()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POINTS);
	//float x, y;

	for (int i = 0; i < target_plot.size(); i++)
	{
		glVertex2i(target_plot[i].x, target_plot[i].y);
		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
	}
	glEnd();
	glFlush();
}
void display_targetplot_line()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_LINES);
	//float x, y;

	for (int i = 0; i < (target_plot.size() - 1); i++)
	{
		//glVertex2i(target_plot[i].x, target_plot[i].y);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glVertex2f(target_plot[i].x, target_plot[i].y);// , target_plot[i].z);
		glVertex2f(target_plot[i + 1].x, target_plot[i + 1].y);// , target_plot[i + 1].z);

		glColor3f(1.0, 1.0, 1.0);
	}
	glEnd();
	glFlush();
}
void display_targetplot_pointv2()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_LINES);
	//float x, y;

	for (int i = 0; i < target_plot.size(); i++)
	{

		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].x, target_plot[i].y);
		glVertex2i((target_plot[i].x + 5), target_plot[i].y);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
	}
	glEnd();
	glFlush();
}
void display_mt201() {
	glutSetWindow(WindowID2);
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(3.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;
	j = blockplot.find_point3d_max(target_plot2, 1);
	plot_win2_max_x = target_plot2[j].x + 10;
	//cout << plot_win2_max_x << ", ";
	j = blockplot.find_point3d_min(target_plot2, 1);
	plot_win2_min_x = target_plot2[j].x - 10;
	//cout << plot_win2_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot2, 2);
	plot_win2_max_y = target_plot2[j].y + 10;
	//cout << plot_win2_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot2, 2);
	plot_win2_min_y = target_plot2[j].y - 10;
	//cout << plot_win2_min_y << endl;
	gluOrtho2D(plot_win2_min_x, plot_win2_max_x, plot_win2_min_y, plot_win2_max_y);
	//gluOrtho2D(-200, 200, -200, 200);
	glBegin(GL_POINTS);
	//glBegin(GL_LINES);
	for (int i = 0; i < target_plot2.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glVertex2i((target_plot2[i].x + 5), target_plot2[i].y);
	}

	glEnd();
	glFlush();
}
void display_mt301() {
	//glutSetWindow(WindowID2);
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.3, 0.3, 0.3, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glLineWidth(1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	int j;
	j = blockplot.find_point3d_max(target_plot2, 1);
	plot_win2_max_x = target_plot2[j].x + 20;
	//cout << plot_win2_max_x << ", ";
	j = blockplot.find_point3d_min(target_plot2, 1);
	plot_win2_min_x = target_plot2[j].x - 20;
	//cout << plot_win2_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot2, 2);
	plot_win2_max_y = target_plot2[j].y + 20;
	//cout << plot_win2_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot2, 2);
	plot_win2_min_y = target_plot2[j].y - 20;
	//cout << plot_win2_min_y << endl;
	gluOrtho2D(plot_win2_min_x, plot_win2_max_x, plot_win2_min_y, plot_win2_max_y);
	//gluOrtho2D(-200, 200, -200, 200);
	//glBegin(GL_POINTS);
	glBegin(GL_LINES);
	for (int i = 0; i < targetl1.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		glVertex2f(targetl1[i].x, targetl1[i].y);
		glVertex2f((targetl2[i].x), targetl2[i].y);
	}
	glEnd();
	glBegin(GL_LINES);
	for (int i = 0; i < targetl1.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		glVertex2f(targetl1[i].x, targetl1[i].y);
		glVertex2f((targetl2[i].x), targetl2[i].y);
	}

	glEnd();
	glFlush();
}
void display_mt401() {
	glutSetWindow(WindowID3);
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;
	j = blockplot.find_point3d_max(target_plot3, 2);
	/*
	//cout << plot_win2_min_y << endl;
	*/
	gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	//gluOrtho2D(40, 120, -10, 30);
	glBegin(GL_POINTS);
	//glBegin(GL_LINES);
	for (int i = 0; i < target_plot3.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
		glVertex2i(target_plot3[i].y, target_plot3[i].z);
		//glVertex2i((target_plot2[i].x + 5), target_plot2[i].y);
	}

	glEnd();
	glFlush();

}
void display_mt402() {
	glutSetWindow(WindowID3);
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.8, 0.8, 0.8, 1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	glBegin(GL_POINTS);
	glPointSize(20.0);
	glColor3f(0.1, 0.1, 0.1);
	if (target_plot3.size()>0) {

		for (int i = 0; i < target_plot3.size(); i++)
		{

			glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
			glVertex2i(target_plot3[i].y, target_plot3[i].z);
		}

	}
	glEnd();
	glBegin(GL_LINES);
	//glPointSize(20.0);
	glLineWidth(10.0);
	glColor3f(0.1, 0.1, 0.1);
	if (targetl31.size()>0) {

		for (int i = 0; i < targetl31.size() - 1; i++)
		{
			glColor3f(color_plot_mat3_l[i].x, color_plot_mat3_l[i].y, color_plot_mat3_l[i].z);
			glVertex2i(targetl31[i].y, targetl31[i].z);
			glVertex2i(targetl32[i].y, targetl32[i].z);
			//for (int i = 0; i < show1.size()-2;i ++) {
			//std::cout << targetl31[i].y << ", " << targetl31[i].z << "~";
			//std::cout << targetl32[i].y << ", " << targetl32[i].z << endl;
			//}
		}
		//glEnd();
	}
	glEnd();
	glFlush();

}
void display_mt403() {
	//https://www.ntu.edu.sg/home/ehchua/programming/opengl/CG_examples.html
	glutSetWindow(WindowID3);
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.3, 0.3, 0.3, 1.0);
	//glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	//glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	//glShadeModel(GL_SMOOTH);   // Enable smooth shading
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glTranslatef(1.5f, 0.0f, 7.0f);  // Move right and into the screen

	// Enable perspective projection with fovy, aspect, zNear and zFar
	//gluPerspective(45.0f, aspect, 0.1f, 100.0f);

	gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 0.0, 1.0, 0.0);

	//gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	glBegin(GL_POINTS);
	glPointSize(20.0);
	glColor3f(0.1, 0.1, 0.1);
	if (target_plot3.size() > 0) {

		for (int i = 0; i < target_plot3.size(); i++)
		{

			glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
			//glVertex2i(target_plot3[i].y, target_plot3[i].z);
			glVertex3f(target_plot3[i].x, target_plot3[i].y, target_plot3[i].z);
		}

	}

	glEnd();
	glFlush();

}

void display_mt501() {
	glutSetWindow(WindowID4);
	//cross section display
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;

	//gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	gluOrtho2D(-10, 160, -10, 30);
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot4.size(); i++)
	{
		glColor3f(0.0, 0, 1.0);
		glVertex2f(target_plot4[i].y, target_plot4[i].z);
		//glVertex2i((target_plot2[i].x + 5), target_plot2[i].y);
	}
	//glBegin(GL_LINES);
	

	glEnd();
	glFlush();

}
void display_mt502() {
	glutSetWindow(WindowID4);
	//cross section display
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);

	//glLoadIdentity();

	std::cout << "inside window 3 function: " << target_plot4.size() << endl;
	//int j;

	glColor3f(0.0, 1.0, 1.0);
	glPointSize(4.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;

	j = blockplot.find_point3d_max(target_plot4, 1);
	plot_win3_max_x = target_plot4[j].x + 20;
	//cout << plot_win2_max_x << ", ";
	j = blockplot.find_point3d_min(target_plot4, 1);
	plot_win3_min_x = target_plot4[j].x - 20;
	//cout << plot_win2_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot4, 2);
	plot_win3_max_y = target_plot4[j].y + 20;
	//cout << plot_win2_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot4, 2);
	plot_win3_min_y = target_plot4[j].y - 20;
	//cout << plot_win2_min_y << endl;
	std::cout << "window scale count target plot 3" << endl;
	gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	//gluOrtho2D(-10, 160, -10, 30);
	
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot4.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat4[i].x, color_plot_mat4[i].y, color_plot_mat4[i].z);
		glVertex2f(target_plot4[i].x, target_plot4[i].y);
	}
	glEnd();
	if (targetl4_1.size() > 0) {
		glLineWidth(2.0);
		glBegin(GL_LINES);
		for (int i = 0; i < targetl4_1.size(); i++) {
			glColor3f(color_plot_mat4_l[i].x, color_plot_mat4_l[i].y, color_plot_mat4_l[i].z);
			glVertex2f(targetl4_1[i].x, targetl4_1[i].y);
			glVertex2f(targetl4_2[i].x, targetl4_2[i].y);

		}
		glEnd();
	}

	glFlush();

}
void display_mt503() {
	glutSetWindow(WindowID5);
	//cross section display
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;

	j = blockplot.find_point3d_max(target_plot3, 1);
	plot_win3_max_x = target_plot3[j].x + 20;
	//cout << plot_win2_max_x << ", ";
	j = blockplot.find_point3d_min(target_plot3, 1);
	plot_win3_min_x = target_plot3[j].x - 20;
	//cout << plot_win2_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot3, 2);
	plot_win3_max_y = target_plot3[j].y + 20;
	//cout << plot_win2_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot3, 2);
	plot_win3_min_y = target_plot3[j].y - 20;
	//cout << plot_win2_min_y << endl;
	std::cout << "window scale count target plot 3" << endl;
	//gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	gluOrtho2D(-10, 160, -10, 30);
	
	glBegin(GL_LINES);
	for (int i = 0; i < target_plot3.size(); i++)
	{
		//glVertex2i(target_plot2[i].x, target_plot2[i].y);
		//glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
		glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
		glVertex2i(target_plot3[i].y, target_plot3[i].z);
		//glVertex2i((target_plot2[i].x + 5), target_plot2[i].y);
	}

	glEnd();
	glFlush();

}
void display_mt504() {
	//created on 17.02.2020
	glutSetWindow(WindowID5);
	//display top plane with entroid 
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int j;

	//gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	gluOrtho2D(-10, 160, -10, 30);
	glBegin(GL_LINES);
	for (int i = 0; i < target_plot3.size(); i++)
	{
		glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
		glVertex2i(target_plot3[i].x, target_plot3[i].y);
	}

	glEnd();
	glFlush();

	glutSetWindow(WindowID6);
	//display top plane with entroid 
	//http://www.lighthouse3d.com/tutorials/glut-tutorial/rendering-to-multiple-subwindows/
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//int j;

	//gluOrtho2D(plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
	gluOrtho2D(-10, 160, -10, 30);
	glBegin(GL_LINES);
	for (int i = 0; i < target_plot3.size(); i++)
	{
		glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
		glVertex2i(target_plot3[i].x, target_plot3[i].y);
	}

	glEnd();
	glFlush();

}
void display_window1_v303()
{
	//updated on 26.11.19
	//glutSetWindow(WindowID1);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);

	glLoadIdentity();
	glLineWidth(1.0);
	std::cout << "inside window 2 function: " << target_plot2.size() << endl;
	int j;

	j = blockplot.find_point3d_max(target_plot2, 1);
	plot_win2_max_x = target_plot2[j].x + 20;
	//cout << plot_win2_max_x << ", ";
	j = blockplot.find_point3d_min(target_plot2, 1);
	plot_win2_min_x = target_plot2[j].x - 20;
	//cout << plot_win2_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot2, 2);
	plot_win2_max_y = target_plot2[j].y + 20;
	//cout << plot_win2_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot2, 2);
	plot_win2_min_y = target_plot2[j].y - 20;
	//cout << plot_win2_min_y << endl;
	std::cout << "window scale count" << endl;

	int a;
	//std::cout << "press any digit to continue: ";
	//std::cin >> a;

	gluOrtho2D(plot_win2_min_x, plot_win2_max_x, plot_win2_min_y, plot_win2_max_y);
	//gluOrtho2D(-200, 200, -200, 200);
	//glBegin(GL_POINTS);

	glBegin(GL_LINES);
	for (int i = 0; i < targetl1.size(); i++)
	{
		glColor3f(color_plot_mat2[i].x, color_plot_mat2[i].y, color_plot_mat2[i].z);
		glVertex2f(targetl1[i].x, targetl1[i].y);
		glVertex2f((targetl2[i].x), targetl2[i].y);
	}
	glEnd();
	std::cout << "window 2 end" << endl;

	// draw coordinate lines
	glColor3f(1.0, 0.0, 0.0); // red x
	glLineWidth(2.0);
	glBegin(GL_LINES);
	// x axis
	glVertex3f(plot_win2_min_x + 30, plot_win2_min_y + 25.0f, 0.0f);//5
	glVertex3f(plot_win2_min_x + 70.0, plot_win2_min_y + 25.0f, 0.0f);

	glVertex3f(plot_win2_min_x + 70.0, plot_win2_min_y + 25.0f, 0.0f);//40.0, 0.0f, 0.0f
	glVertex3f(plot_win2_min_x + 60.0, plot_win2_min_y + 30.0f, 0.0f);//35.0, 5.0f, 0.0f

	glVertex3f(plot_win2_min_x + 70.0, plot_win2_min_y + 25.0f, 0.0f);// 40.0, 0.0f, 0.0f
	glVertex3f(plot_win2_min_x + 60.0, plot_win2_min_y + 20.0f, 0.0f);// 35.0, -5.0f, 0.0f
	glEnd();

	// y axis
	glColor3f(0.0, 1.0, 0.0); // green y
	glBegin(GL_LINES);
	glVertex3f(plot_win2_min_x + 30.0, plot_win2_min_y + 25.0f, 0.0f);//5.0, 0.0f, 0.0f
	glVertex3f(plot_win2_min_x + 30.0, plot_win2_min_y + 45.0f, 0.0f);//5.0, 20.0f, 0.0f

	glVertex3f(plot_win2_min_x + 30.0, plot_win2_min_y + 45.0f, 0.0f);//5.0, 20.0f, 0.0f
	glVertex3f(plot_win2_min_x + 25.0, plot_win2_min_y + 40.0f, 0.0f);//0.0, 15.0f, 0.0f

	glVertex3f(plot_win2_min_x + 30.0, plot_win2_min_y + 45.0f, 0.0f);//5.0, 20.0f, 0.0f
	glVertex3f(plot_win2_min_x + 35.0, plot_win2_min_y + 40.0f, 0.0f);//10.0, 15.0f, 0.0f
	glEnd();



	glFlush();
}
GLfloat angleCube = 0.0f;     // Rotational angle for cube [NEW]
int refreshMills = 15;        // refresh interval in milliseconds [NEW]
void display3D_v101() {
	//initGL 3D
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
														/**/
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix

									// Render a color-cube consisting of 6 quads with different colors
	glLoadIdentity();                 // Reset the model-view matrix
	glTranslatef(-250.0f, -100.0f, -60.0f);  // Move right and into the screen
											 //gluLookAt(550.0, 200.0, 30.0, 80.0, 80.0, 10.0, 550.0, 20.0, 60.0);
											 //glTranslatef(-200.0f, -100.0f, -50.0f);
											 //glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
											 // Top face (y = 1.0f)
											 // Define vertices in counter-clockwise (CCW) order with normal pointing out
	glLineWidth(4.0);
	glBegin(GL_LINES);
	for (int i = 0; i < targetl31.size(); i++) {
		glColor3f(color_plot_mat3_l[i].x, color_plot_mat3_l[i].y, color_plot_mat3_l[i].z);
		glVertex3f(targetl31[i].x, targetl31[i].y, targetl31[i].z);
		glVertex3f(targetl32[i].x, targetl32[i].y, targetl32[i].z);
	}
	glEnd();
	if (target_plot3.size()>0) {
		glBegin(GL_POINTS);
		for (int i = 0; i < target_plot3.size(); i++) {
			glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
			glVertex3f(target_plot3[i].x, target_plot3[i].y, target_plot3[i].z);

		}
		glEnd();
	}

	// Render a pyramid consists of 4 triangles
	//glLoadIdentity();                  // Reset the model-view matrix
	//glTranslatef(-1.5f, 0.0f, -6.0f);  // Move left and into the screen


	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)

	angleCube -= 0.15f;
}

/*  Globals from https://www.davidwparker.com/2011/08/20/opengl-screencast-6-drawing-in-3d-part-1-glut-objects/*/
/*Globals for 3d view*/
double dim = 100.0; /* dimension of orthogonal box */
					/*  Various global state */
					/*  Toggles */
int toggleAxes = 0;   /* toggle axes on and off */
int toggleValues = 1; /* toggle values on and off */
					  /*  Display view */
int th = 15;  /* azimuth of view angle */// original are 0 0 
int ph = 45;  /* elevation of view angle */
int objId = 0;      /* object to draw */
float xminbase1, yminbase1, zminbase1, xmaxbase1, ymaxbase1, zmaxbase1;
void map_target_plot_into_3d(vector<cad3d::Point3D> &vin, cad3d newblock) {

	xminbase1 = -50 + 1;
	yminbase1 = -50 + 1;
	zminbase1 = -20 + 1;
	xmaxbase1 = -50 + 120;
	ymaxbase1 = -50 + 100;
	zmaxbase1 = -20 + 80;

	float xmin_o1, ymin_o1, zmin_o1, xmax_o1, ymax_o1, zmax_o1;
	int ii_temp;

	ii_temp = newblock.find_point3d_max(vin, 1);
	xmax_o1 = vin[ii_temp].x;
	ii_temp = newblock.find_point3d_min(vin, 1);
	xmin_o1 = vin[ii_temp].x;

	ii_temp = newblock.find_point3d_max(vin, 2);
	ymax_o1 = vin[ii_temp].y;
	ii_temp = newblock.find_point3d_min(vin, 2);
	ymin_o1 = vin[ii_temp].y;

	ii_temp = newblock.find_point3d_max(vin, 3);
	zmax_o1 = vin[ii_temp].z;
	ii_temp = newblock.find_point3d_min(vin, 3);
	zmin_o1 = vin[ii_temp].z;

	cad3d::Point3D ptemp1, ptemp2;
	for (int i = 0; i < targetl31.size(); i++) {
		ptemp1.x = (xminbase1 - xmin_o1) + (targetl31[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp1.y = (zminbase1 - zmin_o1) + (targetl31[i].z - zmin_o1)*(zmaxbase1 - zminbase1) / (zmax_o1 - zmin_o1);
		ptemp1.z = (yminbase1 - ymin_o1) + (targetl31[i].y - ymin_o1)*(ymaxbase1 - yminbase1) / (ymax_o1 - ymin_o1);
		targetl31_mod.push_back(ptemp1);
		ptemp2.x = (xminbase1 - xmin_o1) + (targetl32[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp2.y = (zminbase1 - zmin_o1) + (targetl32[i].z - zmin_o1)*(zmaxbase1 - zminbase1) / (zmax_o1 - zmin_o1);
		ptemp2.z = (yminbase1 - ymin_o1) + (targetl32[i].y - ymin_o1)*(ymaxbase1 - yminbase1) / (ymax_o1 - ymin_o1);
		targetl32_mod.push_back(ptemp2);
	}

}
void map_target_plot_into_3d_v2(vector<cad3d::Point3D> &vin, cad3d newblock) {
	// for mapping rectangular objects
	// something is missing with y and z conversion
	
	xminbase1 = -50; yminbase1 = -50; zminbase1 = -20;
	xmaxbase1 = -50 + 120; ymaxbase1 = -50 + 100; zmaxbase1 = -20 + 80;

	float xmin_o1, ymin_o1, zmin_o1, xmax_o1, ymax_o1, zmax_o1;
	int ii_temp;

	ii_temp = newblock.find_point3d_max(vin, 1);
	xmax_o1 = vin[ii_temp].x;
	ii_temp = newblock.find_point3d_min(vin, 1);
	xmin_o1 = vin[ii_temp].x;

	ii_temp = newblock.find_point3d_max(vin, 2);
	ymax_o1 = vin[ii_temp].y;
	ii_temp = newblock.find_point3d_min(vin, 2);
	ymin_o1 = vin[ii_temp].y;

	ii_temp = newblock.find_point3d_max(vin, 3);
	zmax_o1 = vin[ii_temp].z;
	ii_temp = newblock.find_point3d_min(vin, 3);
	zmin_o1 = vin[ii_temp].z;
	std::cout << "mapping for 3d plot in Cartesian" << endl;
	//std::cout << xmin_o1 << "," << ymin_o1 << "," << zmin_o1 << "~";
	//std::cout << xmax_o1 << ", " << ymax_o1 << ", " << zmax_o1 << endl;
	cad3d::Point3D ptemp1, ptemp2, ptemp3;
	for (int i = 0; i < targetl31.size(); i++) {
		/*ptemp1.x = (xmin_o1 - xminbase1) + ((targetl31[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1));
		ptemp1.y = (ymin_o1 - yminbase1) + ((targetl31[i].y - ymin_o1)*(ymaxbase1 - yminbase1) / (ymax_o1 - ymin_o1));
		ptemp1.z = (zmin_o1 - zminbase1) + ((targetl31[i].z - zmin_o1)*(zmaxbase1 - zminbase1) / (zmax_o1 - zmin_o1))*/;
	ptemp1.x = (xminbase1 - xmin_o1 + 0.0) + (targetl31[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
	ptemp1.y = (yminbase1 - zmin_o1 - 0.0) + (targetl31[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
	ptemp1.z = (zminbase1 - ymin_o1 + 0.0) + (targetl31[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
	targetl31_mod.push_back(ptemp1);

	ptemp2.x = (xminbase1 - xmin_o1 + 0.0) + (targetl32[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
	ptemp2.y = (yminbase1 - zmin_o1 - 0.0) + (targetl32[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
	ptemp2.z = (zminbase1 - ymin_o1 + 0.0) + (targetl32[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
	targetl32_mod.push_back(ptemp2);
	}
	std::cout << xminbase1 << "," << yminbase1 << "," << zminbase1 << "~";
	std::cout << xmaxbase1 << ", " << ymaxbase1 << ", " << zmaxbase1 << endl;
	//std::cout <<  << ", " << ymaxbase1 << ", " << zmaxbase1 << endl;
	// this mapping model works fine for rectangular block next test is for polar block
	for (int i = 0; i < target_plot3.size(); i++) {
		ptemp3.x = (xminbase1 - xmin_o1 + 0.0) + (target_plot3[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp3.y = (yminbase1 - zmin_o1 - 0.0) + (target_plot3[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
		ptemp3.z = (zminbase1 - ymin_o1 + 0.0) + (target_plot3[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
		target_plot3_mod.push_back(ptemp3);
	}

	ii_temp = newblock.find_point3d_min(target_plot3_mod, 1);
	std::cout << target_plot3_mod[ii_temp].x << ",";
	ii_temp = newblock.find_point3d_min(target_plot3_mod, 2);
	std::cout << target_plot3_mod[ii_temp].y << ",";
	ii_temp = newblock.find_point3d_min(target_plot3_mod, 3);
	std::cout << target_plot3_mod[ii_temp].z << endl;
}
void map_target_plot_into_3d_v3(vector<cad3d::Point3D> &vin, cad3d newblock) {

	//for mapping circular objects;
	xminbase1 = -50; yminbase1 = -50; zminbase1 = -20;
	xmaxbase1 = -50 + 120; ymaxbase1 = -50 + 100; zmaxbase1 = -20 + 80;

	float xmin_o1, ymin_o1, zmin_o1, xmax_o1, ymax_o1, zmax_o1;
	int ii_temp;

	ii_temp = newblock.find_point3d_max(vin, 1);
	xmax_o1 = vin[ii_temp].x;
	ii_temp = newblock.find_point3d_min(vin, 1);
	xmin_o1 = vin[ii_temp].x;

	ii_temp = newblock.find_point3d_max(vin, 2);
	ymax_o1 = vin[ii_temp].y;
	ii_temp = newblock.find_point3d_min(vin, 2);
	ymin_o1 = vin[ii_temp].y;

	ii_temp = newblock.find_point3d_max(vin, 3);
	zmax_o1 = vin[ii_temp].z;
	ii_temp = newblock.find_point3d_min(vin, 3);
	zmin_o1 = vin[ii_temp].z;
	std::cout << "mapping for 3d plot" << endl;
	//std::cout << xmin_o1 << "," << ymin_o1 << "," << zmin_o1 << "~";
	//std::cout << xmax_o1 << ", " << ymax_o1 << ", " << zmax_o1 << endl;
	cad3d::Point3D ptemp1, ptemp2, ptemp3;
	for (int i = 0; i < targetl31.size(); i++) {

		ptemp1.x = -50.0 + (targetl31[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp1.y = -50.0 + (targetl31[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
		ptemp1.z = -20.0 + (targetl31[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
		targetl31_mod.push_back(ptemp1);

		ptemp2.x = -50.0 + (targetl32[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp2.y = -50.0 + (targetl32[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
		ptemp2.z = -20.0 + (targetl32[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
		targetl32_mod.push_back(ptemp2);
	}
	std::cout << xminbase1 << "," << yminbase1 << "," << zminbase1 << "~";
	std::cout << xmaxbase1 << ", " << ymaxbase1 << ", " << zmaxbase1 << endl;

	for (int i = 0; i < target_plot3.size(); i++) {
		ptemp3.x = -50 + (target_plot3[i].x - xmin_o1)*(xmaxbase1 - xminbase1) / (xmax_o1 - xmin_o1);
		ptemp3.y = -50.0 + (target_plot3[i].z - zmin_o1)*(ymaxbase1 - yminbase1) / (zmax_o1 - zmin_o1);
		ptemp3.z = -20 + (target_plot3[i].y - ymin_o1)*(zmaxbase1 - zminbase1) / (ymax_o1 - ymin_o1);
		target_plot3_mod.push_back(ptemp3);
	}
	ii_temp = newblock.find_point3d_min(target_plot3_mod, 1);
	std::cout << target_plot3_mod[ii_temp].x << ",";
	ii_temp = newblock.find_point3d_min(target_plot3_mod, 2);
	std::cout << target_plot3_mod[ii_temp].y << ",";
	ii_temp = newblock.find_point3d_min(target_plot3_mod, 3);
	std::cout << target_plot3_mod[ii_temp].z << endl;
}
void drawAxes_v3()
{
	double len = 20.0;
	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(5);
	glBegin(GL_LINES);
	glVertex3d(xminbase1 - 1, yminbase1 - 1, zminbase1 - 1);
	glVertex3d(xminbase1 + len, yminbase1 - 1, zminbase1 - 1);
	glColor3f(0.0, 0.0, 1.0);
	glVertex3d(xminbase1 - 1, yminbase1 - 1, zminbase1 - 1);
	glVertex3d(xminbase1 - 1, yminbase1 + len, zminbase1 - 1);
	glColor3f(0.0, 1.0, 0.0);
	glVertex3d(xminbase1 - 1, yminbase1 - 1, zminbase1 - 1);
	glVertex3d(xminbase1 - 1, yminbase1 - 1, zminbase1 + len);
	glEnd();
	std::cout << "drawing axes" << endl;
	//}
}
void drawOutline() {
	//glutWireCone(1, 2, objSlices, objStacks);
	//glColor3f(1.0, 1.0, 0.0);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque
	glLineWidth(2);
	glBegin(GL_LINES);
	//glColor3f(color_plot_mat3_l[i].x, color_plot_mat3_l[i].y, color_plot_mat3_l[i].z);
	//back face in yellow
	glColor3f(0.1, 0.1, 0.1);
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 1));

	//glColor3f(0.4, 1.0, 0.4);
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));

	//glColor3f(0.4, 0.4, 1.0);
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 1), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));

	//glColor3f(0.4, 1.0, 0.4);
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 1));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 1), (zminbase1 + 80));
	glVertex3f((xminbase1 + 120), (yminbase1 + 100), (zminbase1 + 80));

	glEnd();

	//glFlush();
	//glutSwapBuffers();

}
void drawOutline_v2() {
	//glutWireCone(1, 2, objSlices, objStacks);
	//glColor3f(1.0, 1.0, 0.0);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque
	glLineWidth(2);
	glBegin(GL_LINES);
	//glColor3f(color_plot_mat3_l[i].x, color_plot_mat3_l[i].y, color_plot_mat3_l[i].z);
	//back face in yellow
	glColor3f(0.8, 0.8, 0.8);
	glVertex3f(xminbase1, yminbase1, zminbase1);
	glVertex3f(xmaxbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, ymaxbase1, zminbase1);
	glVertex3f(xminbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zminbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zmaxbase1);

	//glColor3f(0.4, 1.0, 0.4);
	glVertex3f(xminbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xminbase1, ymaxbase1, zminbase1);
	glVertex3f(xminbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xminbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xminbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zmaxbase1);

	//glColor3f(0.4, 0.4, 1.0);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zmaxbase1);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xminbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xminbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zmaxbase1);
	glVertex3f(xminbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);

	//glColor3f(0.4, 1.0, 0.4);
	glVertex3f(xminbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zmaxbase1);
	glVertex3f(xmaxbase1, yminbase1, zminbase1);
	glVertex3f(xminbase1, ymaxbase1, zminbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zminbase1);
	glVertex3f(xmaxbase1, yminbase1, zminbase1);
	glVertex3f(xmaxbase1, ymaxbase1, zminbase1);
	
	glEnd();

	//actual drawing
	
	//glFlush();
	//glutSwapBuffers();

}
void display3D_v102() {
	//initGL 3D
	std::cout << "function 3d view" << endl;
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque

	glLoadIdentity();

	/*  Set View Angle */
	glRotated(ph, 1, 0, 0);
	glRotated(th, 0, 1, 0);

	/*  Draw  */
	//drawAxes();//original

	drawAxes_v3();
	drawOutline_v2();
	std::cout << "Axes and outline drawn" << endl;
	//actual drawing
	if (target_plot3.size()>0) {
		glPointSize(4.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < target_plot3.size(); i++) {
			glColor3f(color_plot_mat3[i].x, color_plot_mat3[i].y, color_plot_mat3[i].z);
			glVertex3f(target_plot3_mod[i].x, target_plot3_mod[i].y, target_plot3_mod[i].z);

		}
		glEnd();
	}

	if (targetl31.size()>0) {
		glLineWidth(1.0);
		glBegin(GL_LINES);
		//glBegin(GL_POINTS);
		for (int i = 0; i < targetl31_mod.size(); i++) {
			glColor3f(color_plot_mat3_l[i].x, color_plot_mat3_l[i].y, color_plot_mat3_l[i].z);
			glVertex3f(targetl31_mod[i].x, targetl31_mod[i].y, targetl31_mod[i].z);
			glVertex3f(targetl32_mod[i].x, targetl32_mod[i].y, targetl32_mod[i].z);
		}
		glEnd();
	}

	if (target_poly_3.size()>0) {
		glBegin(GL_POLYGON);
		//glColor3f(1, 0, 0); glVertex3f(-0.6, -0.75, 0.5);
		//glColor3f(0, 1, 0); glVertex3f(0.6, -0.75, 0);
		//glColor3f(0, 0, 1); glVertex3f(0, 0.75, 0);
		//https://cs.lmu.edu/~ray/notes/openglexamples/
		for (int i = 0; i < target_poly_3.size(); i++) {
			glColor3f(color_plot_mat_poly3[i].x, color_plot_mat_poly3[i].y, color_plot_mat_poly3[i].z);
			glVertex3f(target_poly_3[i].x, target_poly_3[i].y, target_poly_3[i].z);
			//glVertex3f(targetl32_mod[i].x, targetl32_mod[i].y, targetl32_mod[i].z);
		}

		glEnd();
	}

	glFlush();
	//glutSwapBuffers();
}
void reshape3D_v101(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
													  // Compute aspect ratio of the new window
	if (height == 0) height = 1;                // To prevent divide by 0
												//width = 500;
												//height = 100;
	GLfloat aspect = (GLfloat)width / (GLfloat)height;
	//GLfloat aspect = 1.0;
	// Set the viewport to cover the new window
	glViewport(-10, -10, width, height);

	// Set the aspect ratio of the clipping volume to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
								  // Enable perspective projection with fovy, aspect, zNear and zFar
	gluPerspective(160.0f, aspect, 0.1f, 80.0f);
	//glOrtho(-5, 500, -5, 150,0.1f,50);
	//gluOrtho2D(-10, -10, 500, 160);
}
void reshape3D_v102(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
													  // Compute aspect ratio of the new window
	double w2h = (height>0) ? (double)width / height : 1;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	/* orthogonal projection */
	glOrtho(-dim*w2h, +dim*w2h, -dim, +dim, -dim, +dim);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
void reshape3D_v104(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
													  // Compute aspect ratio of the new window
													  //double w2h = (height>0) ? (double)width / height : 1;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	/* orthogonal projection */
	glOrtho(-dim, +dim, -dim, +dim, -dim, +dim);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void windowSpecial(int key, int x, int y)
{
	/*  Right arrow key - increase azimuth by 5 degrees */
	if (key == GLUT_KEY_RIGHT) th += 5;
	/*  Left arrow key - decrease azimuth by 5 degrees */
	else if (key == GLUT_KEY_LEFT) th -= 5;
	/*  Up arrow key - increase elevation by 5 degrees */
	else if (key == GLUT_KEY_UP) ph += 5;
	/*  Down arrow key - decrease elevation by 5 degrees */
	else if (key == GLUT_KEY_DOWN) ph -= 5;

	/*  Keep angles to +/-360 degrees */
	th %= 360;
	ph %= 360;

	glutPostRedisplay();
}
void windowKey(unsigned char key, int x, int y)
{
	//  Exit on ESC /
	if (key == 27) exit(0);
	else if (key == 'a' || key == 'A') toggleAxes = 1 - toggleAxes;
	else if (key == 'v' || key == 'V') toggleValues = 1 - toggleValues;

	///  Spacebar - Toggle through shapes /
	else if (key == 32) {
		if (objId == 18) objId = 0;
		else objId++;
	}

	glutPostRedisplay();
}
void windowMenu(int value)
{
	windowKey((unsigned char)value, 0, 0);
}

void reshape3D_v103(GLsizei width, GLsizei height) {
	// GLsizei for non-negative integer
	
	glViewport(-10, -10, width, height);
	/*glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	//glViewport(0, 0, (GLsizei)700, (GLsizei)600);
	*/
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-0.75, 0.75, -1.0, 1.0, 1.0, 600.0);
	glMatrixMode(GL_MODELVIEW);

}
void Timer3d(int iUnused)
{
	glutPostRedisplay();
	glutTimerFunc(30, Timer3d, 0);
}
static bool paused = false;
void handleKeypress(unsigned char key, //The key that was pressed                                                                                                           
	int x, int y) {    //The current mouse coordinates                                                                                  
	switch (key) {
	case 27: //Escape key                                                                                                                                       
		exit(0); //Exit the program                                                                                                                               
	case 'p':
		paused = !paused;
		break;
	}
}
void target_plot_color_update(float a, float b, float c, int n) {
	//int n_plot = color_plot_mat.size();
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat.push_back(dummy_c);
	}
}
void target_plot_update(vector<cad3d::Point3D> &v) {
	for (int i = 0; i < v.size(); i++) {
		target_plot.push_back(v[i]);
	}
}
void free_target_plot() {
	target_plot.resize(0);
}
void target_plot2_color_update(float a, float b, float c, int n) {
	//int n_plot = color_plot_mat.size();
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat2.push_back(dummy_c);
	}
}
void target_plot2_update(vector<cad3d::Point3D> &v) {
	for (int i = 0; i < v.size(); i++) {
		target_plot2.push_back(v[i]);
	}
}
void target_plot2line_update(vector<cad3d::Point3D> &v1, vector<cad3d::Point3D> &v2) {
	for (int i = 0; i < v1.size(); i++) {
		targetl1.push_back(v1[i]);
		targetl2.push_back(v2[i]);
	}
}

void target_plot3_update(vector<cad3d::Point3D> &v) {
	for (int i = 0; i < v.size(); i++) {
		target_plot3.push_back(v[i]);
	}
}
void target_plot3line_update(vector<cad3d::Point3D> &v1, vector<cad3d::Point3D> &v2) {
	for (int i = 0; i < v1.size(); i++) {
		targetl31.push_back(v1[i]);
		targetl32.push_back(v2[i]);
	}
}
void target_plot3line_color_update(float a, float b, float c, int n) {
	//int n_plot = color_plot_mat.size();
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat3_l.push_back(dummy_c);
	}
}
void target_plot3_color_update(float a, float b, float c, int n) {
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat3.push_back(dummy_c);
	}
}
void target_plot3poly_update(vector<cad3d::Point3D> &v1) {
	for (int i = 0; i < v1.size(); i++) {
		target_poly_3.push_back(v1[i]);
	}
}
void target_plot3poly_color_update(float a, float b, float c, int n) {
	//int n_plot = color_plot_mat.size();
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat_poly3.push_back(dummy_c);
	}
}

void target_plot4_update(vector<cad3d::Point3D> &v) {
	for (int i = 0; i < v.size(); i++) {
		target_plot4.push_back(v[i]);
	}
}

void target_plot4line_update(vector<cad3d::Point3D> &v1, vector<cad3d::Point3D> &v2) {
	for (int i = 0; i < v1.size(); i++) {
		targetl4_1.push_back(v1[i]);
		targetl4_2.push_back(v2[i]);
	}
}
void target_plot4line_color_update(float a, float b, float c, int n) {
	//int n_plot = color_plot_mat.size();
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat4_l.push_back(dummy_c);
	}
}
void target_plot4_color_update(float a, float b, float c, int n) {
	color4opengl dummy_c;
	dummy_c.x = a;
	dummy_c.y = b;
	dummy_c.z = c;
	for (int i = 0; i < n; i++) {
		color_plot_mat4.push_back(dummy_c);
	}
}

ogl_Point3D::RenderState rs;
int button, state, x, y;

void display_coord_temp()
{
	rs.display_coord();
}
void idleCallback_temp()
{
	rs.idleCallback();
}
void mouseCallback_temp(int button, int state, int x, int y)
{
	rs.mouseX = x;
	rs.mouseY = y;

	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			rs.mouseLeftDown = true;
		}
		else if (state == GLUT_UP)
			rs.mouseLeftDown = false;
	}
}
void mouseMotionCallback_temp(int x, int y)
{
	if (rs.mouseLeftDown)
	{
		rs.cameraY += (x - rs.mouseX);
		rs.cameraX += (y - rs.mouseY);
		rs.mouseX = x;
		rs.mouseY = y;
	}
}
void display_targetplot_pointv301()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_LINES);
	//float x, y;

	for (int i = 0; i < target_plot.size(); i++)
	{

		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].x, target_plot[i].y);
		glVertex2i((target_plot[i].x + 5), target_plot[i].y);
		//glVertex3i(target_plot[i].x, target_plot[i].y, target_plot[i].z);
	}
	glEnd();
	glFlush();
}

void drawCoordinates() {
	glClear(GL_COLOR_BUFFER_BIT);
	// draw some lines
	glColor3f(1.0, 0.0, 0.0); // red x
	glBegin(GL_LINES);
	// x aix

	glVertex3f(-4.0, 0.0f, 0.0f);
	glVertex3f(4.0, 0.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, 1.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, -1.0f, 0.0f);
	glEnd();

	// y 
	glColor3f(0.0, 1.0, 0.0); // green y
	glBegin(GL_LINES);
	glVertex3f(0.0, -4.0f, 0.0f);
	glVertex3f(0.0, 4.0f, 0.0f);

	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(1.0, 3.0f, 0.0f);

	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(-1.0, 3.0f, 0.0f);
	glEnd();

	// z 
	glColor3f(0.0, 0.0, 1.0); // blue z
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0f, -4.0f);
	glVertex3f(0.0, 0.0f, 4.0f);


	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, 1.0f, 3.0f);

	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, -1.0f, 3.0f);
	glEnd();
	glFlush();
}

void window_range_v67(int &x_max_i, int &x_min_i, int &y_max_i, int &y_min_i, int &z_max_i, int &z_min_i, cad3d &newblock, cad3d::stl_data &block_a) {
	x_max_i = newblock.find_point3d_max(block_a.vertices, 1);
	y_max_i = newblock.find_point3d_max(block_a.vertices, 2);
	z_max_i = newblock.find_point3d_max(block_a.vertices, 3);
	x_min_i = newblock.find_point3d_min(block_a.vertices, 1);
	y_min_i = newblock.find_point3d_min(block_a.vertices, 2);
	z_min_i = newblock.find_point3d_min(block_a.vertices, 3);
	plot_win_min_z = block_a.vertices[z_min_i].z;
	plot_win_max_z = block_a.vertices[z_max_i].z;
	std::cout << block_a.vertices[z_min_i].z << ", " << block_a.vertices[z_max_i].z;
}

void principle_axis_select(int argc, char* argv[], cad3d &newblock, cad3d::stl_data &block_a, vector<cad3d::Point3D> &vin) {

	int x_max_i, x_min_i, y_max_i, y_min_i, z_max_i, z_min_i;
	window_range_v67(x_max_i, x_min_i, y_max_i, y_min_i, z_max_i, z_min_i, newblock, block_a);

	int j;
	//j = blockplot.find_point3d_max(target_plot, 1);
	//window_range_update(x_min_i - 10, x_max_i + 10, y_min_i - 10, y_max_i + 10, z_min_i, slot_plane[ixmax].x, slot_plane[iymin].y, slot_plane[iymax].y);
	vector<cad3d::Point3D> zslice_points;
	float z_now;
	int dummy1;
	init_opengl(argc, argv);

	//along x axis
	//newblock.axis_wise_feature_plane_approximation(block_a, zslice_points, z_now, 1);
	newblock.feature_axis_selection(block_a, zslice_points, z_now, 1);
	target_plot_update(zslice_points);
	target_plot_color_update(0.8, 0.0, 0.0, zslice_points.size());
	glutSetWindow(WindowID101);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(3.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(block_a.vertices[y_min_i].y - 10, block_a.vertices[y_max_i].y + 10, block_a.vertices[z_min_i].z - 10, block_a.vertices[z_max_i].z + 10);
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot.size(); i++)
	{
		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].y, target_plot[i].z);
	}
	glEnd();
	glFlush();
	std::cout << "This is yz view, press any key to proceed" << endl;
	std::cin >> dummy1;
	zslice_points.resize(0); target_plot.resize(0); color_plot_mat.resize(0);

	//along y axis
	//newblock.axis_wise_feature_plane_approximation(block_a, zslice_points, z_now, 2);
	newblock.feature_axis_selection(block_a, zslice_points, z_now, 2);
	target_plot_update(zslice_points);
	target_plot_color_update(0.0, 0.8, 0.0, zslice_points.size());
	glutSetWindow(WindowID102);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(3.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluOrtho2D(block_a.vertices[x_min_i].x - 10, block_a.vertices[x_max_i].x + 10, block_a.vertices[y_min_i].y - 10, block_a.vertices[y_max_i].y + 10);
	gluOrtho2D(block_a.vertices[x_min_i].x - 10, block_a.vertices[x_max_i].x + 10, block_a.vertices[z_min_i].z, block_a.vertices[z_max_i].z + 10);
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot.size(); i++)
	{
		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].x, target_plot[i].z);
	}
	glEnd();
	glFlush();
	std::cout << "This is xz view, press any key to proceed" << endl;
	std::cin >> dummy1;
	zslice_points.resize(0); target_plot.resize(0); color_plot_mat.resize(0);

	//newblock.axis_wise_feature_plane_approximation(block_a, zslice_points, z_now, 3);
	newblock.feature_axis_selection(block_a, zslice_points, z_now, 3);
	target_plot_update(zslice_points);
	target_plot_color_update(0.0, 0.0, 0.8, zslice_points.size());
	glutSetWindow(WindowID102);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(3.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(block_a.vertices[x_min_i].x - 10, block_a.vertices[x_max_i].x + 10, block_a.vertices[y_min_i].y - 10, block_a.vertices[y_max_i].y + 10);
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot.size(); i++)
	{
		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].x, target_plot[i].y);
	}
	glEnd();
	glFlush();
	std::cout << "This is xy view, escape from opengl and proceed" << endl;
	//cin >> dummy1;
	zslice_points.resize(0); target_plot.resize(0); color_plot_mat.resize(0);



	glutKeyboardFunc(handleKeypress);
	glutMainLoop();
	//return 5;
}
void plot_dim_1(cad3d &newblock, cad3d::stl_data &block_a) {
	int x_max_i, x_min_i, y_max_i, y_min_i, z_max_i, z_min_i;
	window_range_v67(x_max_i, x_min_i, y_max_i, y_min_i, z_max_i, z_min_i, newblock, block_a);

	glutSetWindow(WindowID1);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.4, 0.4, 0.4, 1.0);
	glColor3f(0.0, 1.0, 1.0);
	glPointSize(3.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(block_a.vertices[y_min_i].y - 10, block_a.vertices[y_max_i].y + 10, block_a.vertices[z_min_i].z - 10, block_a.vertices[z_max_i].z + 10);
	glBegin(GL_POINTS);
	for (int i = 0; i < target_plot.size(); i++)
	{
		glColor3f(color_plot_mat[i].x, color_plot_mat[i].y, color_plot_mat[i].z);
		glVertex2i(target_plot[i].y, target_plot[i].z);
	}
	glEnd();
	glFlush();
	char dummy1;
	std::cout << "This is yz view, press any key to proceed" << endl;
	//cin >> dummy1;

}

int main(int argc, char* argv[])
{
	std::cout << "Start of processing article 2 works\n";
	
	string stl_base_file_name = "D:/code_works/cad_i17_versions/";//;//"F:/CAD/art2/"
	string model;
	model = "block2";//block1, block2, block4
	string stl_file_name = stl_base_file_name + model + ".STL";
	std::cout << stl_file_name << endl;
	cad3d newblock, newblock2;
	cad3d::stl_data block_a("new");
	cad3d::stl_data block_b("transform_r");

	//string result_quantify_file = "D:/phd_articles/art2/mario_edits/result_block1_part1.txt";

	int tx_required = 0;
	std::cout << "enter 1 if transform required or 0 if not required: ";
	std::cin >> tx_required;

	if (tx_required == 0) {
		ifstream stl_file(stl_file_name.c_str(), std::ios::in | std::ios::binary);
		char header_info[80] = "";
		stl_file.read(header_info, 80);
		char *output = NULL;
		char b[] = "solid CATIA";
		output = strstr(header_info, b);

		if (output) {
			std::cout << "String Found\n";
			newblock.stl_catia_reader(stl_file_name, block_a);
			/*target_plot_update(block_a.vertices);
			target_plot_color_update(0.0, 0.8, 0.0, block_a.vertices.size());*/
			//return 0;
		}
		else {
			std::cout << "catia string not Found\n";
			newblock.stl_reader(stl_file_name, block_a);
			std::cout << "reading performed" << endl;
		}
	}

	if (tx_required == 1) {
		ifstream stl_file(stl_file_name.c_str(), std::ios::in | std::ios::binary);
		char header_info[80] = "";
		stl_file.read(header_info, 80);
		char *output = NULL;
		char b[] = "solid CATIA";
		output = strstr(header_info, b);
		//std::cout << "file name: " << output << endl;
		if (output) {
			std::cout << "String Found\n";
			newblock.stl_catia_reader(stl_file_name, block_b);

		}
		else {
			std::cout << "catia string not Found\n";
			newblock.stl_reader(stl_file_name, block_b);
			std::cout << "reading performed" << endl;
		}
		newblock.stl_axis_tx_v2(block_a, block_b);
	}
	//finding the points onthe required plane
	vector<cad3d::Point3D> centered_block;
	vector<cad3d::Point3D> all_points;
	vector<cad3d::Point3D> zslice_points;
	vector<float> z_range;

	int x_max_i, x_min_i, y_max_i, y_min_i, z_max_i, z_min_i;
	float z_now;
	//newblock.feature_plane_approximation(block_a, zslice_points, z_now);

	int ixmax, iymax, ixmin, iymin;
	/**/
	char dd1;
	std::cout << endl << "press any character to start" << endl;
	std::cin >> dd1;
	//principle_axis_select(argc, argv, newblock, block_a, zslice_points);
	int dim = 3;
	//return 0;
	std::cout << "dimension (1, 2 or 3):" << endl;
	std::cin >> dim;

	//std::cout << "zslice_points size: " << zslice_points.size() << endl;

	//return 0;
	int cart_vs_pol;
	std::cout << "Enter 1 for Cartesean calculation and 2 for Polar calculation " << endl;
	std::cin >> cart_vs_pol;


	vector<cad3d::triangle_index> tslot_i, tbi_1;
	vector<cad3d::Point3D> slot_plane;
	vector<cad3d::Point3D>  v_ext_l1, v_ext_l2, feat_e1, feat_e2, v_ext1, v_ext2;

	vector<cad3d::Point3D> red_vcorner, red_vend_corner;
	vector<cad3d::Point3D> vcorner, vend_corner;

	vector<int> vi_all;
	//tracking use of all points
	vi_all.resize(block_a.vertices.size());
	newblock.initialize_0i(vi_all);

	vector<int> vi_temp551, vi_temp552;
	vector<cad3d::Point3D>  slot_e01, slot_e02, slot_e03, slot_e04;
	vector<cad3d::Point3D>  slot_e1, slot_e2, slot_e3, slot_e4;
	vector<cad3d::Point3D>  slot_ip1, slot_ip2;

	vector<cad3d::Point3D> v2, v1, v2_rev;
	vector<cad3d::Point3D> cross_l1, cross_l2, c_line, cent1, cent2;
	vector<cad3d::Point3D> vline1new, vline2new;
	vector<cad3d::Point3D> v_temp, v_line1, v_line2;
	vector<cad3d::Point3D> point_show_3d;

	int c_line_count = 0;
	std::cout << "press 2 for new centreline count, 1 for old centreline count and 0 for no processing: ";
	std::cin >> c_line_count;
	int cl_n;
	int characterize_2d = 0;
	// 0 no characterization, 1 path based characterization, 2 boundary based characterization
	std::cout << "press 0 for none, 1 for path based and 2 for 2D boundary based slot characterization: ";
	std::cin >> characterize_2d;

	int fine_segment_3d = 0;
	// 0 no no segment, 1 for fine segment processing 
	std::cout << "press 1 for fine 3D segment calucation along a 2d intersecting line, 0 for none: ";
	std::cin >> fine_segment_3d;

	//vector<cad3d::Point3D>  v2_ext_l1, v2_ext_l2, feat2_e1, feat2_e2, v2_ext1, v2_ext2;
	vector<cad3d::Point3D>  vc_line1, vc_line2, v1_ext_l1, v1_ext_l2, feat1_e1, feat1_e2, v1_ext1, v1_ext2;
	vector<cad3d::Point3D>  v_int_l1, v_int_l2;
	vector<cad3d::Point3D>  interpolated_edge_v1, interpolated_edge_v2, interpolated_border_pts;
	vector<cad3d::Point3D>  line_ext_top1, line_ext_top2;
	vector<cad3d::Point3D>   line_ext_low1, line_ext_low2;

	//variables for intersectionig line cross-section generation
	vector<cad3d::Point3D>   cross_line1, cross_line2;
	vector<cad3d::Point3D>   principal_edge1, principal_edge2;
	vector<cad3d::Point3D>   intersect_line1, intersect_line2, single_cs;

	//variables for sub-volume genration
	vector<cad3d::Point3D>   sub_vol_points, sub_vol_slice, t_edge1, t_edge2, c_line_temp;
	vector<cad3d::Point3D>   sub_vol_ext;
	vector<cad3d::Point3D>   f_boundary_line1, f_boundary_line2;
	vector<cad3d::Point3D>   boundary_3d_line1, boundary_3d_line2;
	vector<cad3d::Point3D>   sub_slice1, sub_slice2, v1_sub;
	vector<cad3d::Point3D>   slice_3d, slice3d_l1, slice3d_l2;
	vector<vector<cad3d::Point3D >> sub_vol_main_slices;
	vector<vector<cad3d::Point3D >> sub_vol_sub_slices;

	//vatiable for 2D and 3D characterization 
	int slot_type = 0;
	vector<cad3d::Point3D> normal_intersect, cline_temp, t_edge11, t_edge12, t_interp1, t_interp2, cross_out1, cross_out2;
	vector<cad3d::Point3D> ext_intersect_l1, ext_intersect_l2;

	if (cart_vs_pol == 1) {
		//z_now

		///...1. align and orient if necessary ...///
		newblock.feature_axis_selection(block_a, zslice_points, z_now, dim);
		newblock.cartesean_processing(block_a, slot_plane, tslot_i, v_ext_l1, v_ext_l2, feat_e1, feat_e2);

		//start of z border poiints
		float z_low;
		int i_temp;
		std::cout << "lowest Z: " << endl;
		i_temp = newblock.find_point3d_min(block_a.vertices, 3);
		newblock.find_point_3d_1feature(slot_plane, block_a.vertices, (block_a.vertices[i_temp].z + 0.5), 1.0, 3);

		slot_plane.resize(0); v_ext_l1.resize(0);  v_ext_l2.resize(0); ;  tslot_i.resize(0);
		feat_e1.resize(0); feat_e2.resize(0);// centered_block.resize(0);

		std::cout << "target Z: ";
		std::cin >> z_low;
		auto time_start = chrono::steady_clock::now();
		newblock.find_point_3d_1feature(slot_plane, block_a.vertices, z_low, 1.0, 3);
		/*if (slot_plane.size()>0) {
		target_plot3_update(slot_plane);
		target_plot3_color_update(0.2, 0.2, 0.8, slot_plane.size());
		}*/
		auto time_end = chrono::steady_clock::now();
		auto time_diff = time_end - time_start;
		std::cout << "top feature plane extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
		slot_plane.resize(0); v_ext_l1.resize(0);  v_ext_l2.resize(0); ;  tslot_i.resize(0);
		feat_e1.resize(0); feat_e2.resize(0); //centered_block.resize(0);
		std::cout << "highest Z: ";
		i_temp = newblock.find_point3d_max(block_a.vertices, 3);
		newblock.find_point_3d_1feature(slot_plane, block_a.vertices, (block_a.vertices[i_temp].z - 0.5), 1.0, 3);

		/*if (slot_plane.size()>0) {
		target_plot3_update(slot_plane);
		target_plot3_color_update(0.4, 0.4, 0.8, slot_plane.size());
		}*/
		slot_plane.resize(0); v_ext_l1.resize(0);  v_ext_l2.resize(0); ;  tslot_i.resize(0);
		feat_e1.resize(0); feat_e2.resize(0); centered_block.resize(0);
		//end of z border points

		//getting internal lines in new style =>// new triangle dependent style

		vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
		vector<int> i_line;
		vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
		vector<cad3d::Point3D>  feat01_e1, feat01_e2;
		time_start = chrono::steady_clock::now();
		newblock.cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
		//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
		newblock.triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);
		newblock.order_cartesian_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2);
		time_end = chrono::steady_clock::now();
		time_diff = time_end - time_start;
		std::cout << "cross model and feature edges extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
		int ixmax, iymax, ixmin, iymin;
		std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
		newblock.cartesian_separate_internal_regions_v2(feat1_e1, feat1_e2, v_int_l1, v_int_l2, 0, slot_plane1);

		std::cout << "Enter 1 for  processing semi-open slots and 2 for closed slots: ";
		std::cin >> slot_type;

		///...30i.a groove boundary extract ...///
		if (slot_type == 2) {

			//for closed region
			//if region type 1
			//newblock.cartesian_separate_internal_regions(feat1_e1, feat1_e2, v_int_l1, v_int_l2);
			newblock.single_closed_region_internal_boundary_extract_v2(block_a, feat1_e1, feat1_e2, v_int_l1, v1, slot_e01, slot_e02);

			//int closed_region_no = 0;
			//std::cout << "closed region number: ";
			//std::cin >> closed_region_no;
			int vn = 0;
			int closed_i_temp = 0;
			float it;
			/*if (feat1_e1.size()>0) {
			//target_plot3line_update(v_int_l1, v_int_l2);
			//target_plot3line_color_update(0.6, 0.6, 0.8, v_int_l1.size());
			target_plot3_update(v_int_l1);
			target_plot3_color_update(0.8, 0.2, 0.8, v_int_l1.size());
			target_plot3_update(v1);
			target_plot3_color_update(0.8,0.4,0.4,v1.size());
			target_plot2line_update(feat1_e1, feat1_e2);
			target_plot2_color_update(0.0, 0.8, 0.0, feat1_e1.size());
			target_plot2_update(v_int_l1);
			target_plot2_color_update(0.8, 0.4, 0.4, v_int_l1.size());
			};*/
			/*newblock.multiple_layer_cross_section_extract(block_a, sub_vol_points, c_line, normal_intersect, t_edge11, t_edge12);
			if (t_edge11.size()>0) {
			target_plot3line_update(t_edge11, t_edge12);
			target_plot3line_color_update(0.6, 0.6, 0.8, t_edge11.size());
			}
			newblock.cross_section_outline_generate_v1(t_edge11, t_edge12, cross_out1, cross_out2);
			if (cross_line1.size()>0) {
			target_plot3line_update(cross_out1, cross_out2);
			target_plot3line_color_update(0.8, 0.5, 0.5, cross_out1.size());
			}*/

			/*if (characterize_2d == 2) {
			std::cout << "processing boundary pick points based 2D characterization" << endl;
			int c_line_segment = 10;
			newblock.characterization_2d_4_closed_slot_v1(v1, c_line, cent1, cent2, c_line_segment);
			}*/
		}

		///...30i.b groove boundary extract ...///
		//if region type 2
		if (c_line_count == 1) {
			// old GE style
			std::wcout << "Old GE style: does no processing here" << endl;
		}
		if (c_line_count == 2) {

			//actually now for semi-open slots
			//for semi-open regions
			vector<cad3d::Point3D> pref_start, pref2;
			//finding corner points of semi open regions
			//needs a generalized function

			vector<cad3d::Point3D> open_slot_pts;
			//to be updated
			int model_group = 0;
			std::cout << "press 0 for GE models and 1 for Marco's models:  ";
			std::cin >> model_group;
			//for GE blocks
			if (model_group == 0) {
				ixmax = newblock.find_point3d_max(v_int_l1, 1);
				newblock.find_point_3d_1feature(open_slot_pts, v_int_l1, v_int_l1[ixmax].x, 0.05, 1);
				ixmax = newblock.find_point3d_max(v_int_l2, 1);
				newblock.find_point_3d_1feature(open_slot_pts, v_int_l2, v_int_l2[ixmax].x, 0.05, 1);
			}
			/**/
			//for Marco blocks
			if (model_group == 1) {
				ixmax = newblock.find_point3d_min(v_int_l1, 1);
				newblock.find_point_3d_1feature(open_slot_pts, v_int_l1, v_int_l1[ixmax].x, 0.05, 1);
				ixmax = newblock.find_point3d_min(v_int_l2, 1);
				newblock.find_point_3d_1feature(open_slot_pts, v_int_l2, v_int_l2[ixmax].x, 0.05, 1);
				std::cout << "corner point indices - " << endl;
			}

			for (int i = 0; i < open_slot_pts.size(); i++) {
				std::cout << i << ": " << open_slot_pts[i].x << ", ";
				std::cout << open_slot_pts[i].y << ", ";
				std::cout << open_slot_pts[i].z << endl;
			}
			newblock.single_slot_seperate_edges_v101(block_a, open_slot_pts, v1, v2, v_int_l1, v_int_l2, slot_e1, slot_e2, slot_e3, slot_e4);


			if (characterize_2d == 1) {
				std::cout << "processing path based 2D characterization" << endl;
				//newblock.single_slot_seperate_edges_v102(block_a, open_slot_pts, v1, v2, v_int_l1, v_int_l2, slot_e1, slot_e2, slot_e3, slot_e4);

				//interpolate edges if necessary
				// previous versions are with Lagrange, v5 is with spline header file
				//needs more than 6 points
				newblock.interpolate_edge_if_necessary_v5(slot_e1, slot_e2, slot_ip1, interpolated_edge_v1);
				newblock.interpolate_edge_if_necessary_v5(slot_e3, slot_e4, slot_ip2, interpolated_edge_v2);
				std::cout << "Edge interpolation completed" << endl;

				//newblock.calculate_fplane_slot_centreline(slot_e1, slot_e3, slot_ip1, c_line, cl_n);
				//newblock.calculate_fplane_slot_centreline_v3(interpolated_edge_v1, interpolated_edge_v2, c_line);
				//with prof's idea
				time_start = chrono::steady_clock::now();
				newblock.calculate_fplane_slot_centreline_v4(interpolated_edge_v1, interpolated_edge_v2, c_line);
				//with minimum distance 
				time_end = chrono::steady_clock::now();
				time_diff = time_end - time_start;
				std::cout << "centerline calculation time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
				if (c_line.size() > 0) {
					//newblock.view_scale_for_3d(block_a.vertices, c_line, point_show_3d, plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
					for (int i = 0; i < c_line.size() - 1; i++) {
						cent1.push_back(c_line[i]);
						cent2.push_back(c_line[i + 1]);
					}
				}


			}
			//std::cout << ti_all2.size() << ", " <<tslot_i1.size() << endl;
			/*if (intersect_line1.size()>0) {
			if (intersect_line1.size() == intersect_line2.size()) {
			target_plot3line_update(intersect_line1, intersect_line2);
			target_plot3line_color_update(0.9, 0.0, 0.0, intersect_line1.size());
			target_plot2line_update(intersect_line1, intersect_line2);
			target_plot2_color_update(0.9, 0.0, 0.0, intersect_line1.size());
			}
			}*/
			/*if (cross_line1.size()>0) {
			if (cross_line1.size() == cross_line2.size()) {
			//target_plot3line_update(cross_line1, cross_line2);
			//target_plot3line_color_update(0.9, 0.3, 0.3, cross_line1.size());
			target_plot2line_update(cross_line1, cross_line2);
			target_plot2_color_update(0.9, 0.3, 0.3, cross_line1.size());
			}
			}*/
		}

		///...401a. 2D boundary extract ...///
		int boundary_extract_2d = 0;
		std::wcout << "enter 1 to transform vertex boundary to line boundary: ";
		std::cin >> boundary_extract_2d;
		vector<cad3d::Point3D> dum1, dum2;
		if (boundary_extract_2d == 1) {
			newblock.vertex_2_boundary_generate_v1(v1, v2, f_boundary_line1, f_boundary_line2, cart_vs_pol, slot_type, dum1);
		}

		///...401b. sub_volume extract ...///
		int sub_vol_process = 0;
		std::wcout << "enter 1 to extract generalized sub volume regions: ";
		std::cin >> sub_vol_process;
		if (sub_vol_process == 1) {
			vector<cad3d::Point3D> slot_plane1;
			vector<cad3d::Point3D>   v1_next_slice, slice1_01, slice1_02;
			vector<cad3d::Point3D>   slice2_01, slice2_02;
			vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
			vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
			newblock.sub_volume_generate_v6(block_a, v1, slot_e01, slot_e02, sub_vol_ext, l1_next_slice, l2_next_slice, sub_vol_main_slices);
			std::cout << "sub_volume boundary layers: " << sub_vol_main_slices.size() << endl;
		}
		/*
		if (slot_type == 1) {
		int slot1_intersect = 0;
		std::wcout << "Enter 5* (not rec0mmended) if you want to extract centerline inside slot 1";
		if (slot1_intersect ==5) {
		//newblock.sub_volume_generate_v1(block_a,interpolated_edge_v1,interpolated_edge_v2,sub_vol_points);
		newblock.sub_volume_generate_v2(block_a, v1, sub_vol_points);
		//newblock.adaptive_2p5d_decomposition(sub_vol_points, sub_vol_slice, z_now, 3, 20, 0.5);

		//newblock.single_slot_seperate_edges_v103(block_a, sub_vol_points, t_edge1, t_edge2);
		//newblock.single_slot_seperate_edges_v103(block_a, sub_vol_slice, t_edge1, t_edge2);


		// find normal to centerline at a given point and intersection, start
		}
		}
		*/
		//std::cout << "v1, v2, slot_e1 size" << v1.size() << ", " << v2.size() << ", " << slot_e1.size() << endl;

		///...501.3D fine slot extract ...///
		if (fine_segment_3d == 1) {
			//if ((c_line.size()>0) & (sub_vol_main_slices.size()>0)) {

			///...501a Estimates intersect line only on feature plane ...///
			if (characterize_2d >= 1) {
				if (slot_type == 1) {
					//for semi-open slots
					if (c_line.size()>0) {
						std::cout << "step 501b" << endl;
						int single_cross = 0;
						std::cout << "2D characterization path vector size: " << c_line.size() << endl;
						std::cout << "Enter any of its positive integer index for intersection line estimation or zero for none" << endl;
						std::cin >> single_cross;
						//int non_planar_point_stat = 0;
						cad3d::Point3D ptemp3;
						if (single_cross > 0) {
							//interpolated_edge_v1, interpolated_edge_v2
							if (slot_type == 1) {
								time_start = chrono::steady_clock::now();
								single_cs.push_back(interpolated_edge_v1[single_cross - 1]);
								intersect_line1.push_back(interpolated_edge_v1[single_cross - 1]);
								//newblock.generate_single_cross_along_centreline_v1(block_a, slot_plane1, slot_ip1[single_cross], non_planar_point_stat, single_cs, ti_all2, slot_e1, single_cross);
								ptemp3 = newblock.find_correlated_point(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, single_cross - 1);
								intersect_line2.push_back(ptemp3);
								//find the cross correlated point
								//std::cout << "staep 501a" << endl;
								cross_line1.push_back(interpolated_edge_v1[single_cross]);
								i_temp = newblock.calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, single_cross);
								cross_line2.push_back(interpolated_edge_v2[i_temp]);
								time_end = chrono::steady_clock::now();
								time_diff = time_end - time_start;
								std::cout << "cross section gnereration time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
								//generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross)
								std::cout << "interpolated sizes: " << interpolated_edge_v1.size() << ", " << interpolated_edge_v2.size() << endl;
							}
							//newblock.generate_single_cross_along_centreline_v7(interpolated_edge_v1, interpolated_edge_v2, c_line, normal_intersect, intersect_line1, intersect_line2, 20, 2);
							//generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross)
							time_start = chrono::steady_clock::now();
							newblock.generate_extended_projection_line_v1(intersect_line1, intersect_line2, ext_intersect_l1, ext_intersect_l2, 1.5);
							time_end = chrono::steady_clock::now();
							time_diff = time_end - time_start;
							std::cout << "extended projection line gnereration time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
						}
					}
				}
				if (slot_type == 2) {
					//for closed slots
					std::cout << "processing boundary pick points based 2D characterization" << endl;
					int c_line_segment = 10;
					newblock.characterization_2d_4_closed_slot_v1(v1, c_line, cent1, cent2, c_line_segment);
					//newblock.characterization_2d_4_closed_slot_v1(v1, c_line, intersect_line1, intersect_line2, c_line_segment);
					intersect_line1.push_back(c_line[0]);
					intersect_line2.push_back(c_line[c_line.size() - 1]);
					time_start = chrono::steady_clock::now();
					newblock.generate_extended_projection_line_v1(intersect_line1, intersect_line2, ext_intersect_l1, ext_intersect_l2, 1.5);
					time_end = chrono::steady_clock::now();
					time_diff = time_end - time_start;
					std::cout << "extended projection line gnereration time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
				}
			}

			///...501b Estimates intersect line only on feature plane ...///
			//std::cout << "step 501b" << endl;
			if (sub_vol_main_slices.size()>0) {
				//newblock.generate_3D_sub_layer_v1(sub_vol_main_slices, boundary_3d_line1, boundary_3d_line2, 2);
				time_start = chrono::steady_clock::now();
				newblock.generate_3D_sub_layer_v2(sub_vol_main_slices, sub_vol_sub_slices, boundary_3d_line1, boundary_3d_line2, slot_type);
				time_end = chrono::steady_clock::now();
				time_diff = time_end - time_start;
				std::cout << "sub volume layer genereration time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			}

			///...501c Generate 3D fine slice ...///
			std::cout << "step 501c" << endl;
			time_start = chrono::steady_clock::now();
			newblock.generate_3D_fine_slice_v2(sub_vol_main_slices, ext_intersect_l1, ext_intersect_l2, slice3d_l1, slice3d_l2, slot_type);
			//newblock.generate_3D_fine_slice_v3(sub_vol_main_slices, ext_intersect_l1, ext_intersect_l2, slice3d_l1, slice3d_l2, slot_type);// 1 => slot_type
			time_end = chrono::steady_clock::now();
			time_diff = time_end - time_start;
			std::cout << "2D template genereration time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			std::cout << slice3d_l1.size() << ", " << slice3d_l2.size() << endl;
			/*else {
			std::cout << "Not processing for 3D fine segment, either 2D charcterization and/or sub-volume is missing " << endl;
			}*/
		}
		else {
			std::cout << "Not processing for 3D fine segment" << endl;
		}
	}

	vector<cad3d::Point3D> open_slot_pts, open_slot_pts_polar, open_slot_pts_temp;

	if (cart_vs_pol == 2) {
		float z_low;
		int i_temp;

		///...1. align and orient if necessary ...///
		//newblock.polar_processing(block_a, slot_plane, tslot_i, v_ext_l1, v_ext_l2, feat_e1, feat_e2, z_low, centered_block);
		newblock.centre_at_x0y0(block_a.vertices, centered_block);

		///...2. slice and predict feature plane...///
		newblock.feature_plane_approximation(block_a, zslice_points, z_low);

		i_temp = newblock.find_point3d_min(centered_block, 3);
		newblock.find_point_3d_1feature(slot_plane, centered_block, (centered_block[i_temp].z + 0.5), 1.0, 3);
		slot_plane.resize(0); v_ext_l1.resize(0);  v_ext_l2.resize(0); ;  tslot_i.resize(0);
		feat_e1.resize(0); feat_e2.resize(0);// centered_block.resize(0);

		i_temp = newblock.find_point3d_max(centered_block, 3);
		newblock.find_point_3d_1feature(slot_plane, centered_block, (centered_block[i_temp].z - 0.5), 1.0, 3);
		slot_plane.resize(0); v_ext_l1.resize(0);  v_ext_l2.resize(0); ;  tslot_i.resize(0);
		feat_e1.resize(0); feat_e2.resize(0); //centered_block.resize(0);

		auto time_start = chrono::steady_clock::now();
		newblock.feature_plane_approximation(block_a, zslice_points, z_now);
		auto time_end = chrono::steady_clock::now();
		auto time_diff = time_end - time_start;
		std::cout << "polar feature plane extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;

		time_start = chrono::steady_clock::now();
		newblock.polar_processing(block_a, slot_plane, tslot_i, v_ext_l1, v_ext_l2, feat_e1, feat_e2, z_now, centered_block);
		time_end = chrono::steady_clock::now();
		time_diff = time_end - time_start;
		std::cout << "polar processing time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;

		///...3. feature plane external boundary detection (embedded into)...//
		std::cout << "one_start_point_polar_block_v1" << endl;
		//3a. Getting start points (and ending) of all open or semi open slots in poiints mode
		cad3d::Point3D pref_start;
		pref_start.x = 0;
		newblock.one_start_point_polar_block_v1(block_a, centered_block, pref_start, vcorner, red_vcorner, vend_corner, z_now);


		if (c_line_count == 1) {
			int j = 1;
			std::cout << "slot index: ";
			std::cin >> j;

			v1.push_back(red_vcorner[j]);
			slot_e01.push_back(red_vcorner[j]);
			vi_temp551.resize(0);
			vi_temp552.resize(0);
			std::cout << " v1 first point " << endl;
			newblock.one_start_point_polar_block_v1_phase4(centered_block, red_vcorner[j], feat_e1, feat_e2, v1);
			slot_e02.push_back(v1[v1.size() - 1]);
			//newblock.get_internal_slot_v107_phase3(centered_block, slot_e01, slot_e02, feat_e1, feat_e2, v1);

			time_start = chrono::steady_clock::now();
			newblock.get_internal_slot_v107_phase3(centered_block, slot_e01, slot_e02, feat_e1, feat_e2, v1);
			//newblock.get_internal_slot_v107_phase4(centered_block, slot_e1, slot_e2, feat_e1, feat_e2, v1);
			newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
			auto time_end = chrono::steady_clock::now();
			auto time_diff = time_end - time_start;
			std::cout << "polar feature plane extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			std::cout << " v1 calculation completed: " << slot_e1.size() << endl;

			v2.push_back(vend_corner[j]);
			slot_e03.push_back(vend_corner[j]);
			vi_temp551.resize(0);
			vi_temp552.resize(0);
			//newblock.one_start_point_polar_block_v1_phase2(centered_block, vend_corner[j], feat_e1, feat_e2, v2);
			newblock.one_start_point_polar_block_v1_phase4(centered_block, vend_corner[j], feat_e1, feat_e2, v2);
			slot_e04.push_back(v2[v2.size() - 1]);
			time_start = chrono::steady_clock::now();
			newblock.get_internal_slot_v107_phase5(centered_block, slot_e03, slot_e04, feat_e1, feat_e2, v2);
			newblock.crop_slot_single_edge(slot_e3, slot_e4, slot_e03, slot_e04);
			time_end = chrono::steady_clock::now();
			time_diff = time_end - time_start;
			std::cout << "polar feature plane extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			std::cout << " v2 calculation completed: " << slot_e3.size() << endl;

			//for path symmetric structures
			for (int i = 2; i < slot_e1.size(); i++) {
				slot_ip1.push_back(slot_e1[i]);
				//i = i + 2;
			}
			//finding connected lines for slots
			cad3d::Point3D n0, n1, n2, n_temp, o_temp;
			//vector<cad3d::Point3D>  n0_ray1,n0_ray2, n1_ray1, n1_ray2, n2_ray1,n2_ray2;
			//cad3d::Point3D ncross_temp1, ncross_temp2, ncross1_dir, ncross2_dir;
			vector<cad3d::Point3D> intersect_pt, pt_of_int;
			//newblock.calculate_cross_rays_v3(slot_e1, slot_e3, slot_ip1, n1_ray1, n1_ray2, n2_ray1, n2_ray2);
			//for multiple point
			//for single point use v4
			//int i_vt = 1;
			//finding cross section feature plane lines
			//o_temp = newblock.find_correlated_point(slot_e1, slot_e3, slot_ip1, i_vt);
			//finding centrelines
			int cl_n;
			std::cout << "Input edge size: " << slot_ip1.size() << endl;
			std::cout << "Enter maximum index of input for centerline:";
			std::cin >> cl_n;

			time_start = chrono::steady_clock::now();
			newblock.calculate_fplane_slot_centreline(slot_e1, slot_e3, slot_ip1, c_line, cl_n);
			time_end = chrono::steady_clock::now();
			time_diff = time_end - time_start;
			std::cout << "ceneterline time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			if (c_line.size() > 0) {
				for (int i = 0; i < c_line.size() - 1; i++) {
					cent1.push_back(c_line[i]);
					cent2.push_back(c_line[i + 1]);
					//i = i + 2;
				}
			}
		}
		if (c_line_count == 2) {
			// new triangle dependent style

			///...4. Feature plane triangle data extract ...
			//(i.e., internal features and ...
			//internal connected features processing)...///
			//4a. triangle data extraction//
			vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
			vector<int> i_line;
			vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
			vector<cad3d::Point3D>  feat01_e1, feat01_e2;
			//4b. Relating triangle data with polar modified data//
			time_start = chrono::steady_clock::now();
			newblock.cross_model_generate_polar(block_a, centered_block, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
			auto time_end = chrono::steady_clock::now();
			auto time_diff = time_end - time_start;
			std::cout << "polar cross model generate time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;

			std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
			newblock.triangle2line(centered_block, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);
			// ***reference of points must be centered here

			///...5. Identifying feature regions ...
			//i.e., feature plane intrusion and extrusion finding ...///
			//5a. Useful for polar symmetric outward to inward paths//
			vector<cad3d::Point3D> feat1_e1_polar, feat1_e2_polar, feat01_e1_polar, feat01_e2_polar, slot_plane1_polar;
			newblock.convert_cart2pol(feat01_e1, feat01_e1_polar);
			newblock.convert_cart2pol(feat01_e2, feat01_e2_polar);
			newblock.order_polar_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2, feat01_e1_polar, feat01_e2_polar);
			std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
			//5b. Seperate internal feature lines where lines ...
			//include atleast one internal point than the external boundary lines  //
			newblock.polar_separate_internal_regions_v1(feat1_e1, feat1_e2, v_int_l1, v_int_l2, 0, slot_plane1);

			//5c. Feature regions quantification
			//needs an updated version for removing duplication of points//
			newblock.convert_cart2pol(slot_plane1, slot_plane1_polar);
			int ixmax = newblock.find_point3d_max(slot_plane1_polar, 1);
			vector<cad3d::Point3D> v_int_l1_polar, v_int_l2_polar;
			newblock.convert_cart2pol(v_int_l1, v_int_l1_polar);
			newblock.convert_cart2pol(v_int_l2, v_int_l2_polar);

			//5d. Sorting feature regions, needs update //
			newblock.find_point_3d_1feature(open_slot_pts_polar, v_int_l1_polar, slot_plane1_polar[ixmax].x, 0.1, 1);
			//possible change conditions
			//ixmax = newblock.find_point3d_max(v_int_l2_polar, 1);
			newblock.find_point_3d_1feature(open_slot_pts_polar, v_int_l2_polar, slot_plane1_polar[ixmax].x, 0.1, 1);
			std::cout << open_slot_pts_polar.size() << endl;
			float delz;
			int open_count = 0;
			newblock.convert_pol2cart(open_slot_pts_polar, open_slot_pts_temp);
			//5e Filtering feature regions only within feature plane//
			for (int i = 0; i < open_slot_pts_polar.size(); i++) {
				delz = abs(open_slot_pts_polar[i].z - z_now);
				if (delz < 0.6) {
					open_slot_pts.push_back(open_slot_pts_temp[i]);
					std::cout << open_count << ": " << open_slot_pts[open_count].x << ", ";
					std::cout << open_slot_pts[open_count].y << ", ";
					std::cout << open_slot_pts[open_count].z << endl;
					open_count++;
				}
			}

			//std::cout << open_slot_pts.size() << endl;
			//target_plot3_update(open_slot_pts);
			//target_plot3_color_update(0.8, 0.0, 0.0, open_slot_pts.size());
			//target_plot3_update(feat1_e2);
			//target_plot3_color_update(0.8, 0.0, 0.0, feat1_e2.size());
			//target_plot3line_update(feat1_e1, feat1_e2);
			//target_plot3line_color_update(0.2,0.2,0.2,feat1_e1.size());

			///...6. Sub region extraction in 2D ...///
			int model_group = 0;
			std::cout << "press 0 for GE models and 1 for Marco's models:  ";
			std::cin >> model_group;
			//for GE blocks
			//needs updates for GE blocks
			if (model_group == 0) {
				int i1, i2;

				std::cout << "corner point index 1:";
				std::cin >> i1;
				v1.push_back(open_slot_pts[i1]);//i_st1
				slot_e01.push_back(open_slot_pts[i1]);
				vi_temp551.resize(0);
				//vi_temp552.resize(0);
				newblock.one_start_point_polar_block_v1_phase4(centered_block, open_slot_pts[i1], v_int_l1, v_int_l2, v1);
				slot_e02.push_back(v1[v1.size() - 1]);
				std::cout << "v1 size: " << v1.size() << endl;
				newblock.get_internal_slot_v107_phase7(centered_block, slot_e01, slot_e02, feat_e1, feat_e2, v1);
				//newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);

				/*std::cout << "corner point index 2:";
				std::cin >> i2;
				v2.push_back(open_slot_pts[i2]);//i_st1
				slot_e03.push_back(open_slot_pts[i2]);
				vi_temp551.resize(0);
				//vi_temp552.resize(0);
				newblock.one_start_point_polar_block_v1_phase4(block_a.vertices, open_slot_pts[i2], v_int_l1, v_int_l2, v2);
				slot_e04.push_back(v2[v2.size() - 1]);
				newblock.get_internal_slot_v107_phase3(centered_block, slot_e03, slot_e04, feat_e1, feat_e2, v2);
				newblock.crop_slot_single_edge(slot_e3, slot_e4, slot_e03, slot_e04);
				*/
			}
			//for Marco blocks
			if (model_group == 1) {
				//6a getting the first line of one side of an open/ semi-open slot //
				int i1, i2;
				vector<cad3d::Point3D> v1pol;
				std::cout << "corner point index 1:";
				std::cin >> i1;
				v1.push_back(open_slot_pts[i1]);//i_st1
				slot_e01.push_back(open_slot_pts[i1]);
				vi_temp551.resize(0);
				newblock.one_start_point_polar_block_v1_phase4(centered_block, open_slot_pts[i1], v_int_l1, v_int_l2, v1);
				slot_e02.push_back(v1[v1.size() - 1]);

				//6b Getting inward points after the first line segment //
				time_start = chrono::steady_clock::now();
				newblock.get_internal_slot_v107_phase7(centered_block, slot_e01, slot_e02, v_int_l1, v_int_l2, v1);
				//std::cout << "v1 size: " << v1.size() << endl;
				newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
				auto time_end = chrono::steady_clock::now();
				auto time_diff = time_end - time_start;
				std::cout << "open edge1 detection time time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;

				//6c getting the first line of the other side of an open/ semi-open slot //
				std::cout << "corner point index 2:";
				std::cin >> i2;
				v2.push_back(open_slot_pts[i2]);//i_st1
				slot_e03.push_back(open_slot_pts[i2]);
				vi_temp551.resize(0);
				//6d Getting inward points after the first line segment //
				newblock.one_start_point_polar_block_v1_phase4(centered_block, open_slot_pts[i2], v_int_l1, v_int_l2, v2);
				slot_e04.push_back(v2[v2.size() - 1]);
				time_start = chrono::steady_clock::now();
				newblock.get_internal_slot_v107_phase7(centered_block, slot_e03, slot_e04, v_int_l1, v_int_l2, v2);
				newblock.crop_slot_single_edge(slot_e3, slot_e4, slot_e03, slot_e04);
				time_end = chrono::steady_clock::now();
				time_diff = time_end - time_start;
				std::cout << "open edge2 detection time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			}

			if (characterize_2d == 1) {
				///...7. Interpolate points if necessary for the edges ...///
				std::cout << "slot e01size:" << slot_e01.size() << endl;
				newblock.interpolate_edge_if_necessary_v5(slot_e1, slot_e2, slot_ip1, interpolated_edge_v1);
				newblock.interpolate_edge_if_necessary_v5(slot_e3, slot_e4, slot_ip2, interpolated_edge_v2);

				//...8. 2D centreline calculation on feature plane...//
				//newblock.calculate_fplane_slot_centreline_v3(interpolated_edge_v1, interpolated_edge_v2, c_line);
				//with prof's idea
				time_start = chrono::steady_clock::now();
				newblock.calculate_fplane_slot_centreline_v4(interpolated_edge_v1, interpolated_edge_v2, c_line);
				auto time_end = chrono::steady_clock::now();
				auto time_diff = time_end - time_start;
				std::cout << "centerline calculation time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
				//with minimum distance 
				if (c_line.size() > 0) {
					for (int i = 0; i < c_line.size() - 1; i++) {
						cent1.push_back(c_line[i]);
						cent2.push_back(c_line[i + 1]);
						//i = i + 2;
					}
				}

				//...9. 2D centreline cross-section calculation...//
				//std::cout << ti_all2.size() << ", " <<tslot_i1.size() << endl;
				int non_planar_point_stat = 0;
				int i_temp = 0;
				cad3d::Point3D ptemp3;
				std::cout << "centerline vector length: " << cent1.size() + 1 << endl;
				//int single_cross = 0;
				//std::cout << "Enter positive integer index for cross section extraction or zero for none" << endl;
				/*std::cin >> single_cross;
				if (single_cross > 0) {
				//single_cs.push_back(slot_ip1[single_cross]);
				//newblock.generate_single_cross_along_centreline_v1(block_a, slot_plane1, slot_ip1[single_cross], non_planar_point_stat, single_cs, ti_all2, slot_ip1, single_cross);
				single_cs.push_back(interpolated_edge_v1[single_cross]);
				intersect_line1.push_back(interpolated_edge_v1[single_cross]);
				//newblock.generate_single_cross_along_centreline_v1(block_a, slot_plane1, slot_ip1[single_cross], non_planar_point_stat, single_cs, ti_all2, slot_e1, single_cross);
				ptemp3 = newblock.find_correlated_point(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, single_cross);
				intersect_line2.push_back(ptemp3);
				//generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross)
				std::cout << "interpolated sizes: " << interpolated_edge_v1.size() << ", " << interpolated_edge_v2.size() << endl;
				//find the cross correlated point

				cross_line1.push_back(interpolated_edge_v1[single_cross]);
				i_temp = newblock.calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, single_cross);
				cross_line2.push_back(interpolated_edge_v2[i_temp]);

				//centreline calculation in phase 3 ends
				}
				*/
				if (intersect_line1.size() > 0) {
					/*if (intersect_line1.size() == intersect_line2.size()) {
					target_plot2line_update(intersect_line1, intersect_line2);
					target_plot2_color_update(0.9, 0.0, 0.0, intersect_line1.size());
					}*/

				}

				///...10. sub region and sub volume generate...///
				/*
				///...10. sub region and sub volume generate...///
				//identify differnce between cross-line and intersect_line variables ***
				//...2D centreline cross-section calculation ends...//


				std::cout << "Sub region and sub volume calculation starts" << endl;
				vector<cad3d::Point3D> sub_region_points;
				// becuause these are now open endges
				for (int i = 0; i < v1.size(); i++) {
				sub_region_points.push_back(v1[i]);
				}
				for (int i = 0; i < v2.size(); i++) {
				sub_region_points.push_back(v2[i]);
				}
				// we need the whole sub region in one vector for sub volume generation
				newblock.sub_volume_generate_v3(block_a, centered_block, sub_region_points, sub_vol_points);
				//target_plot3_update(sub_vol_points);
				//target_plot3_color_update(0.2, 0.2, 0.2, sub_vol_points.size());
				//target_plot3_update(sub_region_points);
				//target_plot3_color_update(0.8, 0.2, 0.2, sub_region_points.size());

				//...11. 3D centreline cross-section calculation...//

				//7.8.19 proceed from here
				int i_cent = 0;
				//11a. find normal to centerline at a given point and intersection, start
				std::cout << "centre line index: ";
				std::cin >> i_cent;
				vector<cad3d::Point3D> normal_intersect, near_edge_points;
				cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
				float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;
				float mt1, ct1, mt2, ct2;
				tan1 = newblock.copy_point3d(c_line[i_cent - 1]);
				tan2 = newblock.copy_point3d(c_line[i_cent + 1]);
				norm0 = newblock.copy_point3d(c_line[i_cent]);
				}
				*/
				/*int process_far = 0;
				std::cout << "press 1 to process further" << endl;
				std::cin >> process_far;
				if (process_far == 1) {

				}*/

			}
		}

		///...401a. 2D boundary extract ...///
		std::cout << "Enter 1 for  processing semi-open slots and 2 for closed slots, 3 for open slots: ";
		std::cin >> slot_type;
		int boundary_extract_2d = 0;
		std::wcout << "enter 1 to transform vertex boundary to line boundary: ";
		std::cin >> boundary_extract_2d;

		vector<cad3d::Point3D> dum1, dum2;
		if (boundary_extract_2d == 1) {
			newblock.vertex_2_boundary_generate_v1(v1, v2, f_boundary_line1, f_boundary_line2, cart_vs_pol, slot_type, dum1);
			newblock.update_point3d_vector(v1, principal_edge1);
			newblock.update_point3d_vector(v2, principal_edge2);
		}

		///...401b. sub_volume extract ...///
		int sub_vol_process = 0;
		std::wcout << "enter 1 to extract generalized sub volume regions: ";
		std::cin >> sub_vol_process;

		if (sub_vol_process == 1) {
			vector<cad3d::Point3D> slot_plane1;
			vector<cad3d::Point3D>   v1_next_slice, slice1_01, slice1_02;
			vector<cad3d::Point3D>   slice2_01, slice2_02;
			vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
			vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;

			if (v2.size()>0) {
				//newblock.update_point3d_vector(v2, v1);
				for (int i = v2.size() - 1; i >= 0; i--) {
					v1.push_back(v2[i]);
				}
				v2.resize(0);
			}
			time_start = chrono::steady_clock::now();
			newblock.sub_volume_generate_v7(block_a, v1, slot_e01, slot_e02, sub_vol_ext, l1_next_slice, l2_next_slice, sub_vol_main_slices);
			auto time_end = chrono::steady_clock::now();
			auto time_diff = time_end - time_start;
			std::cout << "sub volume extraction time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
			std::cout << "sub_volume boundary layers: " << sub_vol_main_slices.size() << endl;

		}

		///...501.3D fine slot extract ...///
		vector<cad3d::Point3D> near_edge_points;
		//newblock.generate_single_cross_along_centreline_v5(interpolated_edge_v1, interpolated_edge_v2, c_line, normal_intersect, dum1, dum2,15);
		/*std::cout << "3D slice processings " << endl;
		int i_cent = 0;
		std::cout << "centre line index: ";
		std::cin >> i_cent;
		//normal_intersect.push_back(c_line[i_cent - 1]);
		normal_intersect.push_back(c_line[i_cent]);
		//normal_intersect.push_back(c_line[i_cent + 1]);


		vector<float> cline_d_mat, cline_m_mat;
		int i_dmat, i_mmat;
		cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
		float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;
		float mt1, ct1, mt2, ct2, d11, d12;
		tan1 = newblock.copy_point3d(c_line[i_cent - 1]);
		tan2 = newblock.copy_point3d(c_line[i_cent + 1]);
		norm0 = newblock.copy_point3d(c_line[i_cent]);
		newblock.get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
		m_n = -1 / m_t;
		c_n = -m_n*norm0.x + norm0.y;

		std::cout << "tangent values (m_t, m_n): " << m_t << ", " << m_n << endl;
		newblock.find_close_points_from_edge(interpolated_edge_v1, near_edge_points, c_line[10], 15);
		//newblock.update_point3d_vector(near_edge_points, normal_intersect);
		for (int i = 0; i < near_edge_points.size(); i++) {
		mt1 = (near_edge_points[i].y - c_line[i_cent].y) / (near_edge_points[i].x - c_line[i_cent].x);
		//std::cout << i << ": "<<mt1 << ", " << newblock.dist201(near_edge_points[i], c_line[i_cent])<<endl;
		d11 = newblock.dist201(near_edge_points[i], c_line[i_cent]);
		cline_d_mat.push_back(mt1);
		cline_m_mat.push_back(d11);
		}
		i_dmat = newblock.find_1d_vector_min(cline_d_mat,1);
		normal_intersect.push_back(near_edge_points[i_dmat]);
		near_edge_points.resize(0);
		cline_d_mat.resize(0);
		cline_m_mat.resize(0);
		newblock.find_close_points_from_edge(interpolated_edge_v2, near_edge_points, c_line[10], 15);
		newblock.update_point3d_vector(near_edge_points, normal_intersect);
		near_edge_points.resize(0);*/
		if (fine_segment_3d == 1) {

			if ((c_line.size() > 0) & (sub_vol_main_slices.size() > 0)) {
				//equation of normal line
				//...11. 3D centreline cross-section calculation...//
				//int single_cross = 0;

				std::cout << "2D polar centerline vector size: " << c_line.size() << ", " << cent1.size() << endl;
				std::cout << "Enter any of its positive integer index for intersection line estimation or zero for none: " << endl;
				//std::cin >> single_cross;

				//7.8.19 proceed from here
				//11a. find normal to centerline at a given point and intersection, start
				//3D fine segment
				//newblock.generate_single_cross_along_centreline_v6(interpolated_edge_v1, interpolated_edge_v2, c_line, normal_intersect, intersect_line1, intersect_line2, 10, 2);

				newblock.generate_single_cross_along_centreline_v7(interpolated_edge_v1, interpolated_edge_v2, c_line, normal_intersect, intersect_line1, intersect_line2, 10, 2);

				time_start = chrono::steady_clock::now();
				newblock.generate_extended_projection_line_v1(intersect_line1, intersect_line2, ext_intersect_l1, ext_intersect_l2, 3);
				time_end = chrono::steady_clock::now();
				time_diff = time_end - time_start;
				std::cout << "extended projection line calculation time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;

				time_start = chrono::steady_clock::now();
				newblock.generate_3D_fine_slice_v2(sub_vol_main_slices, ext_intersect_l1, ext_intersect_l2, slice3d_l1, slice3d_l2, 3);// 3=> slot_type
				time_end = chrono::steady_clock::now();
				time_diff = time_end - time_start;
				std::cout << "2D template generation time: " << chrono::duration <double, milli>(time_diff).count() << " ms" << endl;
				std::cout << slice3d_l1.size() << ", " << slice3d_l2.size() << endl;
				//newblock.generate_3D_sub_layer_v1(sub_vol_main_slices, boundary_3d_line1, boundary_3d_line2, 2);
				newblock.generate_3D_sub_layer_v2(sub_vol_main_slices, sub_vol_sub_slices, boundary_3d_line1, boundary_3d_line2, 2);
				std::cout << "sub volume sub layer aize: " << sub_vol_sub_slices.size() << endl;
			}

		}
		else {
			std::cout << "Not processing for 3D fine segment, either 2D charcterization and/or sub-volume is missing " << endl;
		}

	}

	if (dim == 3) {
		//std::cout << "display started" << endl;
		ixmax = newblock.find_point3d_max(zslice_points, 1);
		iymax = newblock.find_point3d_max(zslice_points, 2);
		ixmin = newblock.find_point3d_min(zslice_points, 1);
		iymin = newblock.find_point3d_min(zslice_points, 2);
		//window_range_update(zslice_points[ixmin].x - 10, zslice_points[ixmax].x + 10, zslice_points[iymin].y - 10, zslice_points[iymax].y + 10, zslice_points[ixmin].x, zslice_points[ixmax].x, zslice_points[iymin].y, zslice_points[iymax].y);
		init_opengl(argc, argv);


		//25.11.19
		//visualize feature plane
		int output_stat = 0;
		std::cout << "Enter 1 for fplane, 2 for singulated 2d region, " << endl;;
		std::cout << "3 for feature lines, 4 for singulated 2D and feature lines, " << endl;
		std::cout << "5 for 2D characterization of sub-regions, 6 for sub-volume generation..." << endl;
		std::cout << "7 for fine 2D segment and 8* for detailed methodology steps" << endl;
		std::cout << "9* for sub-volume intermediate slices" << endl;
		std::cin >> output_stat;
		slot_plane.resize(0);
		vector<cad3d::triangle_index> tslot_i_plot;
		vector<cad3d::Point3D> v_temp, v_line1, v_line2;
		vector<int> i_line;
		vector<cad3d::Point3D> v_dummy, nline1, nline2;
		vector<cad3d::Point3D> main_slice1, main_slice2;

		if (cart_vs_pol == 1) {

			//major lines in 3D start
			slot_plane.resize(0);	tslot_i_plot.resize(0);
			v_line1.resize(0);	v_line2.resize(0);
			i_line.resize(0);	v_temp.resize(0);
			newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 20.0);//15
			newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
			newblock.triangle2line(block_a.vertices, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
			target_plot3line_update(v_line1, v_line2);
			target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
			slot_plane.resize(0);	tslot_i_plot.resize(0);
			v_line1.resize(0);	v_line2.resize(0);
			i_line.resize(0);	v_temp.resize(0);
			newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 7.0);
			newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
			newblock.triangle2line(block_a.vertices, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
			target_plot3line_update(v_line1, v_line2);
			target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
			slot_plane.resize(0);	tslot_i_plot.resize(0);
			v_line1.resize(0);	v_line2.resize(0);
			i_line.resize(0);	v_temp.resize(0);
			newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 0.4);
			newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
			newblock.triangle2line(block_a.vertices, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
			target_plot3line_update(v_line1, v_line2);
			target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
			slot_plane.resize(0);	tslot_i_plot.resize(0);
			v_line1.resize(0);	v_line2.resize(0);


			vector<cad3d::Point3D> dum1, dum2;

			newblock.find_point_3d_1feature(slot_plane, block_a.vertices, 14.8, 1.0, 3);
			target_plot3_update(block_a.vertices);
			target_plot3_color_update(0.2, 0.8, 0.2, block_a.vertices.size());
			target_plot2_update(feat1_e1);
			if (v1.size() > 0) {
				target_plot4_update(v1);
				target_plot4_color_update(0.8, 0.0, 0.8, v1.size());
				newblock.vertex_2_boundary_generate_v1(v1, v2, v_line1, v_line2, cart_vs_pol, slot_type, dum1);
				//if (slot_type == 1) { slot_type = 5; }
				target_plot4line_update(v_line1, v_line2);
				target_plot4line_color_update(0.4, 0.4, 0.4, v_line1.size());
				v_line1.resize(0);	v_line2.resize(0);
			}

			/*if (v_line1.size() > 0) {

			}*/

			//target_plot2_color_update(0.9, 0.9, 0.9, feat1_e1.size());


			if (output_stat >= 1) {
				target_plot3_update(slot_plane);
				target_plot3_color_update(0.2, 0.2, 0.8, slot_plane.size());
			}
			if (output_stat >= 3) {
				newblock.getplane_simplfied_block(block_a, slot_plane, tslot_i_plot);
				newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
				newblock.triangle2line(block_a.vertices, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
				target_plot2line_update(v_line1, v_line2);
				target_plot2_color_update(0.8, 0.8, 0.8, v_line1.size());
				if (feat1_e1.size()<1) {
					newblock.order_cartesian_points_v1(v_line1, v_line2, feat1_e1, feat1_e2);
				}
			}

			if (output_stat >= 2) {

				if (feat1_e1.size()>0) {
					target_plot2line_update(feat1_e1, feat1_e2);
					target_plot2_color_update(0.0, 1.01, 0.0, feat1_e1.size());
				}
				std::cout << "v1 and feature 1 sizes: " << v1.size() << ", " << feat_e1.size() << ", " << endl;
				/*if (feat1_e1.size()>0) {
				target_plot3line_update(feat1_e1, feat1_e2);
				target_plot3line_color_update(0.0, 1.01, 0.0, feat1_e1.size());
				}*/
			}

			if (output_stat >= 4) {
				if (v1.size() > 0) {
					target_plot3_update(v1);
					target_plot3_color_update(0.8, 0.0, 0.8, v1.size());
				}
				/*if (v2.size() > 0) {
				target_plot3_update(v2);
				target_plot3_color_update(0.8, 0.0, 0.8, v2.size());
				}*/
				if (slot_e1.size() > 0) {
					target_plot2line_update(slot_e1, slot_e2);
					target_plot2_color_update(0.2, 0.6, 0.2, slot_e1.size());
					//target_plot4_color_update(0.2, 0.6, 0.2, slot_e1.size());
				}
				if (f_boundary_line1.size() > 0) {

					newblock.vertex_2_boundary_generate_v1(v1, v2, f_boundary_line1, f_boundary_line2, cart_vs_pol, slot_type, dum1);
					if (slot_type == 1) { slot_type = 5; }
					target_plot2line_update(f_boundary_line1, f_boundary_line2);
					target_plot2_color_update(0.6, 0.2, 0.2, f_boundary_line1.size());
					//target_plot3line_update(f_boundary_line1, f_boundary_line2);
					//target_plot3line_color_update(1.0,0.8,0.8, f_boundary_line1.size());
				}
			}
			if ((output_stat == 5) || (output_stat >= 7)) {
				//target_plot3.push_back(point_show_3d[i]);
				//target_plot3.push_back(slot_e3[i]);
				//target_plot3.push_back(point_show_3d[size(point_show_3d)-1]);
				//target_plot3.push_back(slot_e3[size(slot_e3) - 1]);
				/*if (c_line.size()>0) {
				target_plot3_update(c_line);
				target_plot3_color_update(0.2, 0.2, 0.2, c_line.size());
				}
				if (cent1.size()>0) {
				target_plot3line_update(cent1, cent2);
				target_plot3line_color_update(1.2, 0.2, 0.6, cent1.size());
				target_plot2line_update(cent1, cent2);
				target_plot2_color_update(1.2, 0.2, 0.6, cent1.size());
				}*/
			}
			if (output_stat >= 6) {

				if (sub_vol_ext.size()>0) {
					std::cout << "plotting single layer" << endl;
					target_plot3_update(sub_vol_ext);
					target_plot3_color_update(0.8, 0.2, 0.2, sub_vol_ext.size());
				}
				for (int i = 0; i < sub_vol_main_slices.size(); i++) {
					newblock.vertex_2_boundary_generate_v1(sub_vol_main_slices[i], v_dummy, main_slice1, main_slice2, cart_vs_pol, slot_type, v1);//cart_vs_pol, slot_type
					target_plot3line_update(main_slice1, main_slice2);
					target_plot3line_color_update(0.6, 0.2, 0.2, main_slice1.size());
					main_slice1.resize(0);
					main_slice2.resize(0);
				}
				/*if (c_line.size()>0) {
				target_plot4_update(c_line);
				target_plot4_color_update(0.2, 0.2, 0.2, c_line.size());
				}*/
				if (cent1.size()>0) {
					target_plot4line_update(cent1, cent2);
					target_plot4line_color_update(1.2, 0.2, 0.6, cent1.size());
				}
				/*for (int i = 0; i < sub_vol_sub_slices.size(); i++) {
				//newblock.update_point3d_vector();
				newblock.vertex_2_boundary_generate_v1(sub_vol_sub_slices[i], v_dummy, sub_slice1, sub_slice2, cart_vs_pol, slot_type, v1_sub);//cart_vs_pol, slot_type
				target_plot3line_update(sub_slice1, sub_slice2);
				target_plot3line_color_update(0.5, 0.6, 1.0, sub_slice1.size());
				sub_slice1.resize(0);
				sub_slice2.resize(0);
				//target_plot3_update(sub_vol_sub_slices[i]);
				//target_plot3_color_update(0.7, 0.7, 1.0, sub_vol_sub_slices[i].size());
				}*/
			}
			if (output_stat == 7) {
				/*if (intersect_line1.size() > 0) {
				target_plot3line_update(intersect_line1, intersect_line2);
				target_plot3line_color_update(0.2, 0.2, 0.2, intersect_line1.size());
				//target_plot2line_update(intersect_line1, intersect_line2);
				//target_plot2_color_update(0.2, 0.2, 0.2, intersect_line1.size());
				}*/
				if (ext_intersect_l1.size() > 0) {
					target_plot2line_update(ext_intersect_l1, ext_intersect_l2);
					target_plot2_color_update(1.2, 0.0, 0.0, ext_intersect_l1.size());
					//target_plot3line_update(ext_intersect_l1, ext_intersect_l2);
					//target_plot3line_color_update(0.2, 0.2, 1.2, ext_intersect_l1.size());
				}
				if (slice3d_l1.size() > 0) {
					target_plot3_update(slice3d_l1);
					target_plot3_color_update(0, 1, 0.5, slice3d_l1.size());
					target_plot3line_update(slice3d_l1, slice3d_l2);
					target_plot3line_color_update(0.2, 0.2, 0.2, slice3d_l1.size());
					//target_plot3poly_update(slice3d_l1);
					//target_plot3poly_color_update(1.0, 0.7, 0.7, slice3d_l1.size());
				}
				if (slice3d_l2.size() > 0) {
					target_plot3_update(slice3d_l2);
					target_plot3_color_update(0.2, 0.2, 0.2, slice3d_l2.size());
					//target_plot3poly_update(slice3d_l1);
					//target_poly_3.push_back(slice3d_l2[slice3d_l2.size() - 1]);
					//target_plot3poly_color_update(0.8, 0.2, 0.2, slice3d_l1.size()+1);
				}
				/*if (boundary_3d_line1.size() > 0) {
				target_plot3line_update(boundary_3d_line1, boundary_3d_line2);
				target_plot3line_color_update(0.4, 0.4, 0.8, boundary_3d_line1.size());
				}*/

			}
			if (output_stat == 9) {

				target_plot3_update(sub_vol_ext);
				target_plot3_color_update(0.8, 0.2, 0.2, sub_vol_ext.size());
				/*for (int i = 0; i < sub_vol_main_slices.size(); i++) {
				newblock.vertex_2_boundary_generate_v1(sub_vol_main_slices[i], v_dummy, main_slice1, main_slice2, cart_vs_pol, slot_type, v1);//cart_vs_pol, slot_type
				target_plot3line_update(main_slice1, main_slice2);
				target_plot3line_color_update(1.0, 0.7, 0.7, main_slice1.size());
				main_slice1.resize(0);
				main_slice2.resize(0);
				}*/
			}
		}
		if (cart_vs_pol == 2) {
			target_plot3_update(centered_block);
			target_plot3_color_update(0.2, 0.8, 0.2, centered_block.size());
			target_plot2_update(centered_block);
			newblock.find_point_3d_1feature(slot_plane, centered_block, 14.8, 0.5, 3);

			vector<cad3d::Point3D> dum1, dum2;

			if (v1.size() > 0) {
				target_plot4_update(v1);
				target_plot4_color_update(0.8, 0.0, 0.8, v1.size());
				newblock.vertex_2_boundary_generate_v1(v1, v2, v_line1, v_line2, cart_vs_pol, slot_type, dum1);
				//if (slot_type == 1) { slot_type = 5; }
				target_plot4line_update(v_line1, v_line2);
				target_plot4line_color_update(0.4, 0.4, 0.4, v_line1.size());
				v_line1.resize(0);	v_line2.resize(0);
			}

			if (output_stat >= 1) {
				//major lines in 3D start
				slot_plane.resize(0);	tslot_i_plot.resize(0);
				v_line1.resize(0);	v_line2.resize(0);
				i_line.resize(0);	v_temp.resize(0);
				newblock.polar_processing(block_a, slot_plane, tslot_i_plot, v_ext_l1, v_ext_l2, feat_e1, feat_e2, 15.0, centered_block);
				//newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 15.0);
				newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
				newblock.triangle2line(centered_block, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
				target_plot3line_update(v_line1, v_line2);
				target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
				slot_plane.resize(0);	tslot_i_plot.resize(0);
				v_line1.resize(0);	v_line2.resize(0);
				i_line.resize(0);	v_temp.resize(0);
				newblock.polar_processing(block_a, slot_plane, tslot_i_plot, v_ext_l1, v_ext_l2, feat_e1, feat_e2, 7.0, centered_block);
				//newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 7.0);
				newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
				newblock.triangle2line(centered_block, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
				target_plot3line_update(v_line1, v_line2);
				target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
				slot_plane.resize(0);	tslot_i_plot.resize(0);
				v_line1.resize(0);	v_line2.resize(0);
				i_line.resize(0);	v_temp.resize(0);
				newblock.polar_processing(block_a, slot_plane, tslot_i_plot, v_ext_l1, v_ext_l2, feat_e1, feat_e2, 0.0, centered_block);
				//newblock.getplane_simplfied_block_v2(block_a, slot_plane, tslot_i_plot, 0.4);
				newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
				newblock.triangle2line(centered_block, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
				target_plot3line_update(v_line1, v_line2);
				target_plot3line_color_update(0.6, 0.6, 0.6, v_line1.size());
				slot_plane.resize(0);	tslot_i_plot.resize(0);
				v_line1.resize(0);	v_line2.resize(0);
				i_line.resize(0);	v_temp.resize(0);
				newblock.find_point_3d_1feature(v_temp, centered_block, 14.8, 0.5, 3);
				target_plot3_update(v_temp);
				target_plot3_color_update(0.0, 0.0, 1.0, v_temp.size());
				v_temp.resize(0);
				//major lines in 3d end

				//newblock.find_point_3d_1feature(slot_plane, centered_block, 7.1, 1.0, 3);//14.5
				target_plot3_update(slot_plane);
				target_plot3_color_update(0.2, 0.2, 0.8, slot_plane.size());
			}
			if (output_stat >= 3) {
				newblock.polar_processing(block_a, slot_plane, tslot_i_plot, v_ext_l1, v_ext_l2, feat_e1, feat_e2, z_now, centered_block);
				//newblock.getplane_simplfied_block(block_a, slot_plane, tslot_i_plot);
				newblock.triangle_normal_update(block_a.triangles, tslot_i_plot);
				newblock.triangle2line(centered_block, tslot_i_plot, v_line1, v_line2, i_line, v_temp);
				target_plot2line_update(v_line1, v_line2);
				target_plot2_color_update(0.8, 0.8, 0.8, v_line1.size());

			}
			if (output_stat >= 2) {
				if (v1.size() > 0) {
					target_plot3_update(v1);
					target_plot3_color_update(0.8, 0.0, 0.8, v1.size());
				}
				if (v2.size() > 0) {
					target_plot3_update(v2);
					target_plot3_color_update(0.8, 0.0, 0.8, v2.size());
				}
				//target_plot2line_update(feat1_e1, feat1_e2);
				//target_plot2_color_update(0.8, 0.0, 0.0, feat1_e1.size());
				if (feat1_e1.size()>0) {
					target_plot2line_update(feat1_e1, feat1_e2);
					target_plot2_color_update(0.0, 1.01, 0.0, feat1_e1.size());
				}
				if (slot_e1.size() > 0) {
					target_plot2line_update(slot_e1, slot_e2);
					target_plot2_color_update(0.8, 0.0, 0.8, slot_e1.size());
					//target_plot3line_update(slot_e1, slot_e2);
					//target_plot3line_color_update(0.0, 1.01, 0.0, slot_e1.size());
					//target_plot3line_update(slot_e3, slot_e4);
					//target_plot3line_color_update(0.0, 1.01, 0.0, slot_e3.size());
				}
				/*if (feat1_e1.size()>0) {
				target_plot3line_update(feat1_e1, feat1_e2);
				target_plot3line_color_update(0.0, 1.01, 0.0, feat1_e1.size());
				}*/

			}
			if (output_stat >= 4) {
				if (v1.size() > 0) {
					target_plot3_update(v1);
					target_plot3_color_update(0.8, 0.0, 0.8, v1.size());
				}
				if (slot_e1.size() > 0) {
					target_plot2line_update(slot_e1, slot_e2);
					target_plot2_color_update(0.2, 0.6, 0.2, slot_e1.size());
				}
				if (f_boundary_line1.size() > 0) {
					target_plot2line_update(f_boundary_line1, f_boundary_line2);
					target_plot2_color_update(0.6, 0.2, 0.2, f_boundary_line1.size());
					//target_plot3line_update(f_boundary_line1, f_boundary_line2);
					//target_plot3line_color_update(1.0, 0.8, 0.8, f_boundary_line1.size());
				}
			}
			if ((output_stat == 5) || (output_stat >= 7)) {
				/*if (c_line.size()>0) {
				//target_plot3_update(c_line);
				//target_plot3_color_update(0.8, 0.8, 0.2, c_line.size());
				target_plot3_update(cent1);
				target_plot3_color_update(0.8, 0.8, 0.2, cent1.size());
				}*/
				if (cent1.size()>0) {
					//target_plot3line_update(cent1, cent2);
					//target_plot3line_color_update(1.2, 0.2, 0.6, cent1.size());
					target_plot2line_update(cent1, cent2);
					target_plot2_color_update(1.2, 0.2, 0.6, cent1.size());
				}
			}
			if (output_stat >= 6) {
				slot_type = 4;
				target_plot3_update(sub_vol_ext);
				target_plot3_color_update(0.8, 0.2, 0.2, sub_vol_ext.size());
				for (int i = 0; i < sub_vol_main_slices.size(); i++) {
					newblock.vertex_2_boundary_generate_v1(sub_vol_main_slices[i], v_dummy, main_slice1, main_slice2, cart_vs_pol, 4, f_boundary_line1);//cart_vs_pol, slot_type//3,4
					target_plot3line_update(main_slice1, main_slice2);
					target_plot3line_color_update(0.6, 0.2, 0.2, main_slice1.size());
					main_slice1.resize(0);
					main_slice2.resize(0);
				}
				if (cent1.size()>0) {
					target_plot4line_update(cent1, cent2);
					target_plot4line_color_update(1.2, 0.2, 0.6, cent1.size());
				}
				/*for (int i = 0; i < sub_vol_sub_slices.size(); i++) {
				newblock.vertex_2_boundary_generate_v1(sub_vol_sub_slices[i], v_dummy, sub_slice1, sub_slice2, cart_vs_pol, 4, v1_sub);//cart_vs_pol, slot_type
				target_plot3line_update(sub_slice1, sub_slice2);
				target_plot3line_color_update(0.5, 0.6, 1.0, sub_slice1.size());
				target_plot3_update(sub_slice1);
				target_plot3_color_update(0.9, 0.9, 0.9, sub_slice1.size());
				sub_slice1.resize(0);
				sub_slice2.resize(0);
				}*/

			}
			if (output_stat == 7) {

				/*if (intersect_line1.size() > 0) {
				target_plot2line_update(intersect_line1, intersect_line2);
				target_plot2_color_update(0.2, 0.2, 0.2, intersect_line1.size());
				}*/
				/*if (normal_intersect.size()>0) {
				target_plot3_update(normal_intersect);
				target_plot3_color_update(0.2, 0.2, 0.2, normal_intersect.size());
				newblock.vertex_2_boundary_generate_v1(normal_intersect, v_dummy, nline1, nline2, 1, 1, f_boundary_line1);
				}*/
				/*if (intersect_line1.size() > 0) {
				target_plot3line_update(intersect_line1, intersect_line2);
				target_plot3line_color_update(0.2, 0.2, 0.2, intersect_line1.size());
				target_plot2line_update(intersect_line1, intersect_line2);
				target_plot2_color_update(0.2, 0.2, 0.2, intersect_line1.size());
				}*/
				/*if (boundary_3d_line1.size() > 0) {
				target_plot3line_update(boundary_3d_line1, boundary_3d_line2);
				target_plot3line_color_update(0.4, 0.4, 0.8, boundary_3d_line1.size());
				}*/
				if (ext_intersect_l1.size() > 0) {
					target_plot2line_update(ext_intersect_l1, ext_intersect_l2);
					target_plot2_color_update(1.2, 0.0, 0.0, ext_intersect_l1.size());
					//target_plot3line_update(ext_intersect_l1, ext_intersect_l2);
					//target_plot3line_color_update(0.2, 0.2, 1.2, ext_intersect_l1.size());
				}
				if (slice3d_l1.size() > 0) {
					target_plot3_update(slice3d_l1);
					target_plot3_color_update(0, 1, 0.5, slice3d_l1.size());
					target_plot3line_update(slice3d_l1, slice3d_l2);
					target_plot3line_color_update(0.2, 0.2, 0.2, slice3d_l1.size());
				}
				if (slice3d_l2.size() > 0) {
					target_plot3_update(slice3d_l2);
					target_plot3_color_update(0.2, 0.2, 0.2, slice3d_l2.size());
				}
			}
		}



		if ((targetl31.size() > 0) || (target_plot3.size()>0)) {
			if (cart_vs_pol == 1) {
				map_target_plot_into_3d_v2(block_a.vertices, newblock);
			}
			if (cart_vs_pol == 2) {
				map_target_plot_into_3d_v3(centered_block, newblock);
			}

			std::cout << "target_plot_3 size: " << target_plot3.size() << ", " << targetl31.size() << endl;
			//target_plot3_color_update(0.8, 0.0, 0.0, target_plot3.size());
			glutInitWindowPosition(50, 50);
			glutInitWindowSize(800, 600);

			WindowID3 = glutCreateWindow("Window 1: Complete3D view");
			//glutDisplayFunc(display_mt401); 
			//glutDisplayFunc(display_mt402); // for points
			std::cout << "3d view function call" << endl;
			glutDisplayFunc(display3D_v102);
			glutReshapeFunc(reshape3D_v104);       // Register callback handler for window re-size event
			glutKeyboardFunc(windowKey);
			glutSpecialFunc(windowSpecial);
			glutCreateMenu(windowMenu);
			glutAddMenuEntry("Toggle Axes [a]", 'a');
			glutAttachMenu(GLUT_RIGHT_BUTTON);
			std::cout << targetl1.size() << endl;


		}
		if (targetl1.size() > 0) {
			//WindowID1 = glutCreateWindow("Window 2: Single slot in 2D");
			glutInitWindowPosition(850, 50);//400,150
			glutInitWindowSize(400, 295);
			WindowID3 = glutCreateWindow("Window 2: Top view in 2D");
			//glViewport(0, 0, (GLsizei)400, (GLsizei)295);
			glutDisplayFunc(display_window1_v303);
			glutMouseFunc(mouseCallback_temp);
			glutMotionFunc(mouseMotionCallback_temp);
			//glutInitWindowPosition(850, 50);
			//glutInitWindowSize(400, 275);
			glMatrixMode(GL_PROJECTION);
			//glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			//glMatrixMode(GL_MODELVIEW);


		}

		if (target_plot3.size() > 0) {
			//std::cout << "cross section size" << target_plot4.size() << endl;
			//target_plot3_color_update(0.8, 0.0, 0.0, target_plot3.size());
			//target_plot4
			glutInitWindowPosition(850, 375);
			//glutInitWindowSize(400, 275);
			glutInitWindowSize(400, 275);//400,150
			WindowID4 = glutCreateWindow("Window 3: Localized top view");
			//glutDisplayFunc(display_mt201); // for only points
			glutDisplayFunc(display_mt502); // for points


		}


		glutKeyboardFunc(handleKeypress);
		glutMainLoop();
	}



	//generalized plot

	//opengl plot ends
	return 0;
}





