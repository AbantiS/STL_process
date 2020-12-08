//http://makble.com/draw-coordinate-lines-in-opengl
//#include "pch.h"
#include "stdafx.h"
#include <iostream>

#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\glew\glew.h"
#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\freeglut\freeglut.h"
//#include "F:\code_works\opengl6\Dependencies\"

#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\glfw\GLFW\glfw3.h"
#include <string>

//GLFWwindow* window;

#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\glm\glm.hpp"
#include <D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\glm\gtc\matrix_transform.hpp>
using namespace std;
//using namespace glm;

#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\Dependencies\common\shader_s.h"
//#include "F:\code_works\opengl6\Dependencies\oglp\loadshaders.h"
//#include "F:\code_works\opengl6\Dependencies\freeglut\glxext.h""
//#include "F:\code_works\opengl_versions\opengl65\opengl6\interface_ogl3d.h"
#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\cad_i17\interface_ogl3d.h"
#include<math.h> 




using namespace std;
using namespace ogl_Point3D;

/*void RenderState::init_coord(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_FLAT);

	RenderState* rsp = new RenderState();
	rs = *rsp;
}*/
RenderState::RenderState() {
	this->mouseX = 0;
	this->mouseY = 0;
	this->mouseLeftDown = false;
	this->mouseRightDown = false;
	this->cameraX = 0.0f;
	this->cameraY = 0.0f;
}
void RenderState::init_coord_v2(RenderState &rs)
{
	glClearColor(0.1, 0.1, 0.3, 0.0);
	glShadeModel(GL_FLAT);

	RenderState* rsp = new RenderState();
	rs = *rsp;
}
void RenderState::exit_v2(RenderState &rs) {
	delete &rs;
}
void RenderState::set_render_range(float xmin,float xmax, float ymin, float ymax) {
	arrrow_x_min = xmin;
	arrow_x_max = xmax;
	arrrow_y_min = ymin;
	arrow_y_max = ymax;
}
void RenderState::drawCoordinates() {

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

}
void RenderState::drawCoordinates2() {
	//updates the arrows according to input 3d
	glColor3f(1.0, 0.0, 0.0); // red x
	glBegin(GL_LINES);
	// x aix

	glVertex3f(arrrow_x_min, arrrow_y_min, 0.0f);
	glVertex3f(arrow_x_max, arrrow_y_min, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, 1.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, -1.0f, 0.0f);
	glEnd();

	// y 
	glColor3f(0.0, 1.0, 0.0); // green y
	glBegin(GL_LINES);
	glVertex3f(arrrow_x_min, arrrow_y_min, 0.0f);
	glVertex3f(arrrow_x_min, arrow_y_max, 0.0f);

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

}
//void RenderState::display_coord(RenderState &rs)

void RenderState::display_coord()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 1.0);
	glLoadIdentity();
	glTranslatef(0, 0, -15.0f);
	glRotatef(cameraX, 1, 0, 0);
	glRotatef(cameraY, 0, 1, 0);

	glScalef(1.0, 2.0, 1.0);
	glutWireCube(1.0);
	drawCoordinates();

	glFlush();

}
void RenderState::display_coord2()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 1.0);
	glLoadIdentity();
	glTranslatef(0, 0, -15.0f);
	glRotatef(cameraX, 1, 0, 0);
	glRotatef(cameraY, 0, 1, 0);

	glScalef(2.0, 2.0, 2.0);
	glutWireCube(1.0);
	drawCoordinates();

	glFlush();

}
void RenderState::display_coord3()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 1.0);
	glLoadIdentity();
	glTranslatef(0, 0, -15.0f);
	glRotatef(cameraX, 1, 0, 0);
	glRotatef(cameraY, 0, 1, 0);

	glScalef(5.0, 5.0, 5.0);
	glutWireCube(1.0);
	drawCoordinates2();

	glFlush();

}
//pmf = &RenderState::display_coord; // Assign the function to the pointer
/*void RenderState::pmf_displaycoord_generate(&rs) {
	 *pmf = rs.display_coord; // Assign the function to the pointer

}*/
void RenderState::mouseCallback(RenderState &rs, int button, int state, int x, int y)
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
void RenderState::mouseMotionCallback(RenderState &rs, int x, int y)
{
	if (rs.mouseLeftDown)
	{
		rs.cameraY += (x - rs.mouseX);
		rs.cameraX += (y - rs.mouseY);
		rs.mouseX = x;
		rs.mouseY = y;
	}
}
void RenderState::idleCallback()
{
	glutPostRedisplay();
}
