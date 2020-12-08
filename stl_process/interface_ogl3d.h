#pragma once

namespace ogl_Point3D {
	class RenderState {
	public:
		int x, y;
		float arrrow_x_min, arrow_x_max;
		float arrrow_y_min, arrow_y_max;
		float arrrow_z_min, arrow_z_max;
		float mouseX, mouseY, cameraX, cameraY;
		bool mouseLeftDown, mouseRightDown;
		//pmf = &RenderState::display_coord;
		RenderState(); /*{
			this->mouseX = 0;
			this->mouseY = 0;
			this->mouseLeftDown = false;
			this->mouseRightDown = false;
			this->cameraX = 0.0f;
			this->cameraY = 0.0f;
		}*/
		//struct RS_data {};

		//void init_coord();
		void init_coord_v2(RenderState &rs);
		void exit_v2(RenderState &rs);
		void drawCoordinates();
		void drawCoordinates2();
		//void drawCoordinates3();
		void set_render_range(float xmin, float xmax, float ymin, float ymax);
		
		void (RenderState::* pmf)(void); // pointer to member function
		//pmf = &RenderState::display_coord; // Assign the function to the pointer
		//void pmf_displaycoord_generate(&rs); // pointer to member function
		void display_coord();
		void display_coord2();
		void display_coord3();
		void mouseCallback(RenderState &rs, int button, int state, int x, int y);
		void mouseMotionCallback(RenderState &rs, int x, int y);
		void idleCallback();
	};
}