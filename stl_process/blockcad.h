#pragma once
#pragma once
#pragma once
#include <string>
#include <vector>
namespace cadblock {
	class cad3d {
	public:
		struct Point3D
		{
			float x, y, z; //in cartesion if polar its theta,r,z
						   //int i;//index in the original point cloud
			int r1, r2; // repetition status: 0,1,2
						//int orient;// 3 point orientatin measurement
		};
		struct triangle {
			Point3D normal;
			Point3D v1;
			Point3D v2;
			Point3D v3;
			triangle(Point3D normalp, Point3D v1p, Point3D v2p, Point3D v3p) :
				normal(normalp), v1(v1p), v2(v2p), v3(v3p) {}
		};
		struct triangle_index {
			int ti;//triangle index
			int v1i;
			int v2i;
			int v3i;
		};
		struct stl_data {
			std::string name;
			std::vector<triangle> triangles;
			std::vector<Point3D> vertices; //error corrected in header file with std::

			stl_data(std::string namep) : name(namep) {}
		};
		struct circle2d
		{
			float x, y; // center coordinates
			float r; //radius
		};
		/*struct color4opengl {
		float x, y, z;
		};*/
		Point3D p0;

		float parse_float(std::ifstream& s);
		Point3D parse_point(std::ifstream& s);
		void stl_reader(const std::string& f, stl_data& info);
		void point3d_write(std::string &f, std::vector<Point3D> &v1);

		void stl_catia_reader(const std::string& f, stl_data& info);
		Point3D split_char_array(const std::string& s, char delimiter, std::vector<std::string> &tokens);

		int find_point3d_max(std::vector<cad3d::Point3D> &v1, int dim);
		int find_point3d_min(std::vector<cad3d::Point3D> &v1, int dim);
		void find_point_3d_1feature(std::vector<Point3D> &vo, std::vector<Point3D> &vi, float test_value, float tol, int dim);
		void find_point3d_i_1feature(std::vector<int> &vo, std::vector<Point3D> &vi, float test_value, float tol, int dim);
		int find_1d_vector_max(std::vector<float> &v1, int dim);
		int find_1d_vector_min(std::vector<float> &v1, int dim);

		float dist201(Point3D p1, Point3D p2);
		float dist301(Point3D p1, Point3D p2);
		float del_theta(Point3D pol1, Point3D pol2);
		void find_Point3D_index(Point3D ref, std::vector<Point3D> &vall, std::vector<int> &vi_used, std::vector<int> &v_temp);
		void find_Point3D_indexv2(Point3D ref, std::vector<Point3D> &vall, std::vector<int> &v_temp);
		void find_Point3D_indexv3(Point3D ref, std::vector<Point3D> &vall, std::vector<int> &vi_used, std::vector<int> &v_temp);

		void exclude_2i_from_1i(std::vector<int> &v1, std::vector<int> &v2, std::vector<int> &vout);
		void initialize_0i(std::vector<int> &v);
		void show_index_frequency(std::vector<int> &v);
		void i_2_point3d(std::vector<int> &v_i, std::vector<Point3D> &v_all, std::vector<Point3D> &vo);

		void centre_at_x0y0(std::vector<Point3D> &vin, std::vector<Point3D> &vout);
		void triangle_centre_at_x0y0(std::vector<Point3D> &vin, std::vector<triangle> &tin, std::vector<triangle> &tout);
		void find_zrange(std::vector<Point3D> &vin, std::vector<float> &z_range);

		void triangle_vertexindex_init(std::vector<triangle> &t_all, std::vector<triangle_index> &ti_all);
		void triangle_allpoint_1feature(std::vector<triangle> &t_all, std::vector<triangle> &tout, float ref, float tol, int dim);
		void triangle_allpoint_i1feature(std::vector<Point3D> &p_all, std::vector<triangle_index> &ti_all, std::vector<triangle_index> &ti_out, float ref, float tol, int dim);
		void triangle2line(std::vector<Point3D> &centered_block, std::vector<triangle_index> &tslot_i2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<int> &i_line, std::vector<Point3D> &v_temp);
		void triangle_anypoint_1feature(std::vector<triangle> &t_all, std::vector<triangle> &tout, float ref, float tol, int dim);
		void triangle_anypoint_i1feature(std::vector<Point3D> &p_all, std::vector<triangle_index> &ti_all, std::vector<triangle_index> &ti_out, float ref, float tol, int dim);

		void convert_cart2pol(std::vector<Point3D> &vin, std::vector<Point3D> &vout);
		void convert_cart2polPoint3D(Point3D &vin, Point3D &vout);
		void convert_pol2cart(std::vector<Point3D> &vin, std::vector<Point3D> &vout);
		void convert_pol2cartPoint3D(Point3D &vin, Point3D &vout);

		void find_external_segmentsv4(std::vector<Point3D> &v_temp, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<int> &used_i_v_external, std::vector<triangle_index> &t_slot, std::vector<int> &viu);
		int get_next_segmentv2(std::vector<Point3D> &vtl1, std::vector<Point3D> &vtl2, Point3D ref, std::vector<int> &viu);
		void get_external_partv2(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &viu2, Point3D &ref, std::vector<int> &vline_i, std::vector<Point3D> &vselected, std::vector<int> &vi_used, std::vector<Point3D> &vall);
		void get_external_partv3(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &viu2, Point3D &ref, std::vector<int> &vline_i, std::vector<Point3D> &vselected, std::vector<int> &vi_used, std::vector<Point3D> &vall);
		void get_segment_numberv5(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected);
		void get_segment_numberv501(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected);
		
		void get_segment_numberv502(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected);
		void get_segment_numberv503(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected, std::vector<Point3D> &vcorner);
		void get_segment_numberv504(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected, std::vector<Point3D> &vcorner);
		
		void get_segment_numberv505(std::vector<Point3D> &vline1, std::vector<Point3D> &vline2, std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2, std::vector<int> &vlineusedindex, std::vector<int> &vi_all, std::vector<Point3D> &vall, std::vector<Point3D> &vselected, std::vector<Point3D> &vcorner);
		void get_end_corner_v101( std::vector<Point3D> &vlineordered1, std::vector<Point3D> &vlineordered2,  std::vector<Point3D> &vendcorner, std::vector<Point3D> &vcorner);
		
		int get_trianglesegment_wo_rv3(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_trianglesegment_wo_rv301(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_trianglesegment_wo_rv302(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		/*void mark_inner_radius(vector<int> &vi_used, vector<Point3D> &vall);*/
		int get_trianglesegment_wo_rv303(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_trianglesegment_wo_rv3031(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_1st_inward_point_v101(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_1st_inward_point_v102(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);

		circle2d compute_circle_4_point3d(Point3D p1, Point3D p2, Point3D p3);
		int get_slot_part(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_slot_part_v101(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);
		int get_trianglesegment_wo_rv304(std::vector<int> &vi_used, std::vector<Point3D> &vall, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int iref);

		void get_internal_slot_v101(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1);
		void get_internal_slot_v102(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1);
		void get_internal_slot_v103(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1);
		void get_internal_slot_v1031(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1);
		
		void relate_corner_w_line(std::vector<int> &vi_all, std::vector<Point3D> &centered_block,  std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &vcorner,std::vector<Point3D> &vcorner_new);
		
		void get_internal_slot_v104(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1s1, int vi_temp1s2);
		void get_internal_slot_v105(std::vector<int> &vi_all, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v1, int vi_temp1s1, int vi_temp1s2);
		void get_internal_slot_v106(stl_data &block_a, std::vector<int> &vi_all, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<int> &vi_temp552, std::vector<Point3D> &v1);

		void get_internal_slot_v107(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<int> &vi_temp552, std::vector<Point3D> &v1);
		void get_internal_slot_v107_phase2(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_internal_slot_v107_phase3(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_internal_slot_v107_phase4(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_internal_slot_v107_phase5(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_internal_slot_v107_phase6(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_internal_slot_v107_phase7(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		
		void get_internal_slot_v108(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<int> &vi_temp552, std::vector<Point3D> &v1);
		int point_repeat_check_v1(Point3D p1, std::vector<Point3D> &vin);
		int point_repeat_check_v2(Point3D p1, std::vector<Point3D> &vin);
		int point_superposition_check_v1(Point3D p1, std::vector<Point3D> &vin, float tol);

		void getplane_simplfied_block(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i);
		void getplane_simplfied_block_v2(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i,float z_temp);

		void one_start_point_simplified_block(stl_data &block_a, std::vector<int> &vi_all, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, Point3D pref_start, std::vector<int> &vi_temp551, std::vector<int> &vi_temp552);

		void cross_multiply_point3d(Point3D p1, Point3D p2, Point3D p3, Point3D &p_normal);
		void interpolate_lagrange_v5(std::vector<Point3D> &v1, int ii, int n1, Point3D &v2);
		float segment_length_measure(std::vector<Point3D> &vin);
		void interpolate_lagrange_v6(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present);
		void interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present);

		void triangle_normal_update(std::vector<triangle> &t_all, std::vector<triangle_index> &t_i);
		triangle generate_triangle(Point3D p1, Point3D p2, Point3D p3);
		Point3D copy_point3d(Point3D pin);
		void update_point3d_vector(std::vector<Point3D> v_in, std::vector<Point3D> &v_out);
		Point3D feature_edge_predict_phase1(Point3D p1, Point3D p2);
		
		int find_any_int(std::vector<int> &v_i, int test_val);
		int find_line_repete_index(Point3D p1, Point3D p2,std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int present_i, std::vector<int> &repete_i);
		void feature_edge_predict_phase2(std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2);
		void feature_edge_predict_phase0(std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2);

		void cartesean_processing(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2);
		int feature_plane_approximation(stl_data &block_a, std::vector<Point3D> &zslice_points, float &z_now);
		int axis_wise_feature_plane_approximation(stl_data &block_a, std::vector<Point3D> &zslice_points, float &z_now, int dim);
		void apply_interpolation(std::vector<Point3D> &v_in, std::vector<Point3D> &v_out, int point_no);
		int feature_axis_selection(stl_data &block_a, std::vector<Point3D> &zslice_points, float &z_now, int dim);
		int adaptive_2p5d_decomposition(std::vector<Point3D> &v_in, std::vector<Point3D> &zslice_points, float &z_now, int dim, int n_total, float z_tol);
		void fine_slice_calculation(std::vector<Point3D> &v_in, float &z_now, int dim, float z_tol);

		//x y dimensiono ones
		void cartesean_processing_extended(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2,int dim, float z_now);
		void getplane_simplfied_block_extended(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i,int dim,float z_now);

		Point3D find_cross_ray(Point3D &p1, Point3D &p2, Point3D &p3);
		Point3D return_unitv_pt3d_mode(Point3D p1);
		void calculate_cross_rays_v1(Point3D p1, Point3D p2, Point3D p3, Point3D &n0, Point3D &n1, Point3D &n2);
		void calculate_cross_rays_v2(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &n1_ray1, std::vector<Point3D> &n1_ray2, std::vector<Point3D> &n2_ray1, std::vector<Point3D> &n2_ray2);
		void calculate_cross_rays_v3(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &n1_ray1, std::vector<Point3D> &n1_ray2, std::vector<Point3D> &n2_ray1, std::vector<Point3D> &n2_ray2);
		void calculate_cross_rays_v4(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &n1_ray1, std::vector<Point3D> &n1_ray2, std::vector<Point3D> &n2_ray1, std::vector<Point3D> &n2_ray2,int i);

		Point3D test_n_correct_slope_direction(Point3D p0, Point3D n1, std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3);
		float DistancePointLine_modv1(Point3D p0, Point3D &LineStart, Point3D &LineEnd, float &dist);
		
		void polar_processing(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2, float z_now, std::vector<Point3D> &centered_block);
		void one_start_point_polar_block_v1(stl_data &block_a, std::vector<Point3D> &centered_block, Point3D pref_start, std::vector<Point3D> &vcorner, std::vector<Point3D> &red_vcorner, std::vector<Point3D> &vend_corner, float z_now);
		void one_start_point_polar_block_v1_phase2(std::vector<Point3D> &centered_block, Point3D pref_start, std::vector<int> &vi_all, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<int> &vi_temp551, std::vector<int> &vi_temp552, std::vector<Point3D> &v_out);
		void one_start_point_polar_block_v1_phase3(std::vector<Point3D> &centered_block, Point3D pref_start, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v_out);
		void one_start_point_polar_block_v1_phase4(std::vector<Point3D> &centered_block, Point3D pref_start, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &v_out);
		
		void find_feature_edge_break_points_v1(stl_data &block_a, std::vector<int> &vi_all, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &pref_start, std::vector<int> &vi_temp551, std::vector<int> &vi_temp552);
		void find_feature_edge_break_points_v2(stl_data &block_a, std::vector<int> &vi_all, std::vector<Point3D> &slot_line1, std::vector<Point3D> &slot_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &pref_start, std::vector<int> &vi_temp551, std::vector<int> &vi_temp552);
		void find_cartesian_external_points_v1(stl_data &block_a, std::vector<int> &vi_all, std::vector<Point3D> &b_line1, std::vector<Point3D> &b_line2, std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &pref_start, std::vector<int> &vi_temp551, std::vector<int> &vi_temp552);
		void crop_slot_single_edge(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_e01, std::vector<Point3D> &slot_e02);
		void crop_slot_single_edge_v2(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_e01, std::vector<Point3D> &slot_e02);

		Point3D find_correlated_point(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, int i_vt);
		Point3D find_correlated_point_v2(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, int i_vt);
		Point3D find_correlated_point_v3(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, int i_vt);
		void calculate_fplane_slot_centreline(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &c_line, int cl_n);
		void calculate_fplane_slot_centreline_v2(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &c_line, int cl_n);
		void calculate_fplane_slot_centreline_v3(std::vector<Point3D> &interpolated_edge_v1, std::vector<Point3D> &interpolated_edge_v2, std::vector<Point3D> &c_line);
		void calculate_fplane_slot_centreline_v4(std::vector<Point3D> &interpolated_edge_v1, std::vector<Point3D> &interpolated_edge_v2, std::vector<Point3D> &c_line);


		void stl_axis_tx_v1(stl_data &block_a, stl_data &block_b);
		void stl_axis_tx_v2(stl_data &block_a, stl_data &block_b);
		void cartesean_processing_v3(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &tslot_i, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2);
		
		void cross_model_generate(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &ti_all2, std::vector<triangle_index> &tslot_i1, std::vector<triangle_index> &tslot_i2, std::vector<triangle_index> &tsloti_feature);
		void cross_model_generate_v2(stl_data &block_a, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &ti_all2, std::vector<triangle_index> &tslot_i1, std::vector<triangle_index> &tslot_i2, std::vector<triangle_index> &tsloti_feature, float tol);
		void cross_model_generate_polar(stl_data &block_a, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_plane, std::vector<triangle_index> &ti_all2, std::vector<triangle_index> &tslot_i1, std::vector<triangle_index> &tslot_i2, std::vector<triangle_index> &tsloti_feature);
		void order_cartesian_points_v1(std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &vord1_line1, std::vector<Point3D> &vord1_line2);
		void order_polar_points_v1(std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, std::vector<Point3D> &vord1_line1, std::vector<Point3D> &vord1_line2, std::vector<Point3D> &vpolar_line1, std::vector<Point3D> &vpolar_line2);
		
		void cartesian_separate_internal_regions(std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2);
		void cartesian_separate_internal_regions_v2(std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2, int cp_mode, std::vector<Point3D> &slot_plane1);
		void single_closed_region_internal_boundary_extract(stl_data &block_a, std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v1, std::vector<Point3D> &slot_e01, std::vector<Point3D> &slot_e02);
		void single_closed_region_internal_boundary_extract_v2(stl_data &block_a, std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v1, std::vector<Point3D> &slot_e01, std::vector<Point3D> &slot_e02);
		void find_closed_slot_number_v1(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2, std::vector<Point3D> &v_used, int internal_v_length);
		float find_closed_slot_number_v2(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feat_e1, std::vector<Point3D> &feat_e2, std::vector<Point3D> &v_used);
		
		void polar_separate_internal_regions_v1(std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2, int cp_mode, std::vector<Point3D> &slot_plane1);

		void get_semi_open_slot_v101(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);
		void get_semi_open_slot_v102(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4, std::vector<Point3D> &feature_e1, std::vector<Point3D> &feature_e2, std::vector<Point3D> &v2);

		Point3D find_non_planar_connected_points(stl_data &block_a, std::vector<Point3D> &slot_plane1, std::vector<triangle_index> &tslot_i1, Point3D pin);
		Point3D find_non_planar_connected_points_v2(stl_data &block_a, std::vector<Point3D> &slot_plane1, std::vector<triangle_index> &tslot_i1, Point3D pin, int &status);
		Point3D find_non_planar_connected_points_v3(stl_data &block_a, std::vector<Point3D> &slot_plane1, std::vector<triangle_index> &tslot_i1, Point3D pin, int &status);
		Point3D find_un_connected_non_planar_points(stl_data &block_a, Point3D pin, int & status);
		Point3D find_un_connected_non_planar_points_v2(stl_data &block_a, std::vector<Point3D> &centered_block, Point3D pin, int & status);
		
		void generate_single_cross_along_centreline_v1(stl_data &block_a, std::vector<Point3D> &slot_plane1, Point3D pin, int &status_slot_plane1, std::vector<Point3D> &single_cs, std::vector<triangle_index> &ti_all2, std::vector<Point3D> &slot_ip1, int single_cross);
		void generate_single_cross_along_centreline_v101(stl_data &block_a, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_plane1, Point3D pin, int &status_slot_plane1, std::vector<Point3D> &single_cs, std::vector<triangle_index> &ti_all2, std::vector<Point3D> &slot_ip1, int single_cross);
		int find_lower_point(stl_data &block_a, std::vector<Point3D> &slot_plane1, std::vector<triangle_index> &ti_all2, Point3D &pout, Point3D  pin);
		int find_lower_point_v2(stl_data &block_a, std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_plane1, std::vector<triangle_index> &ti_all2, Point3D &pout, Point3D  pin);
		
		int calculate_intersectline_from_interpolated_curves_v1(std::vector<Point3D> &interpolated_edge_v1, std::vector<Point3D> &interpolated_edge_v2, int i_present);
		
		void sub_volume_generate_v1(stl_data &block_a, std::vector<Point3D> &interpolated_edge_v1, std::vector<Point3D> &interpolated_edge_v2,std::vector<Point3D> &sub_vol);
		void sub_volume_generate_v2(stl_data &block_a, std::vector<Point3D> &v_in,  std::vector<Point3D> &sub_vol);
		void sub_volume_generate_v3(stl_data &block_a, std::vector<Point3D> &v_all, std::vector<Point3D> &v_in, std::vector<Point3D> &sub_vol);
		void sub_volume_generate_v4(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<Point3D> &sub_vol);
		int sub_volume_generate_v5(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out);
		void sub_volume_generate_v6(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out, std::vector<std::vector<cad3d::Point3D >> &sub_vol_main_slices);
		int sub_volume_generate_v5_polar(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out);
		void sub_volume_generate_v7(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out, std::vector<std::vector<cad3d::Point3D >> &sub_vol_main_slices);
		
		void vertex_2_boundary_generate_v1(std::vector<Point3D> &v1, std::vector<Point3D> &v2, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, int cart_vs_pol, int slot_type, std::vector<Point3D> &v_feat_plane);
		void characterization_2d_4_closed_slot_v1(std::vector<Point3D> &v_border, std::vector<Point3D> &c_line, std::vector<Point3D> &cent1, std::vector<Point3D> &cent2, int c_line_count);

		void sub_volume_z_level_calculate_v1(std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<float> z_levels);
		int sub_volume_z_level_calculate_v2(std::vector<Point3D> &v_in, std::vector<Point3D> &l1_in, std::vector<Point3D> &l2_in, std::vector<float> z_levels);
		void get_next_sub_slice_v1(std::vector<Point3D> &v_in, std::vector<Point3D> &feat01_e1, std::vector<Point3D> &feat01_e2, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out);
		void get_next_sub_slice_v2(std::vector<Point3D> &v_in, std::vector<Point3D> &feat01_e1, std::vector<Point3D> &feat01_e2, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out, float tol);
		void get_next_sub_slice_v3(std::vector<Point3D> &v_in, std::vector<Point3D> &feat01_e1, std::vector<Point3D> &feat01_e2, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out, float tol);
		void get_next_sub_slice_v4(std::vector<Point3D> &centered_block, std::vector<Point3D> &slot_plane, std::vector<Point3D> &v_in, std::vector<Point3D> &feat01_e1, std::vector<Point3D> &feat01_e2, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out, float tol);
		
		void get_same_z_sub_slice_v1(std::vector<Point3D> &v_in, std::vector<Point3D> &feat01_e1, std::vector<Point3D> &feat01_e2, std::vector<Point3D> &sub_vol, std::vector<Point3D> &l1_out, std::vector<Point3D> &l2_out);
		void get_same_z_sub_slice_v2(stl_data &block_a, std::vector<Point3D> v_last_slice, std::vector<Point3D> l1_last_slice, std::vector<Point3D> l2_last_slice, std::vector<Point3D> &v_next_slice);
		void get_same_z_sub_slice_v3(stl_data &block_a, std::vector<Point3D> v_last_slice, std::vector<Point3D> l1_last_slice, std::vector<Point3D> l2_last_slice, std::vector<Point3D> &v_next_slice, std::vector<Point3D> l1_next_slice, std::vector<Point3D> l2_next_slice);
		void get_same_z_sub_slice_v4(stl_data &block_a, std::vector<Point3D> v_last_slice, std::vector<Point3D> l1_last_slice, std::vector<Point3D> l2_last_slice, std::vector<Point3D> &v_next_slice, std::vector<Point3D> l1_next_slice, std::vector<Point3D> l2_next_slice);
		void get_same_z_sub_slice_v5(stl_data &block_a, std::vector<Point3D> v_last_slice, std::vector<Point3D> l1_last_slice, std::vector<Point3D> l2_last_slice, std::vector<Point3D> &v_next_slice, std::vector<Point3D> l1_next_slice, std::vector<Point3D> l2_next_slice);
		int test_if_in_shadow_v1(Point3D v_in, std::vector<Point3D> v_ref, float tol);
		int test_if_in_shadow_v2(Point3D v_in, std::vector<Point3D> v_ref, float tol);


		void view_scale_for_3d(std::vector<Point3D> &v_all, std::vector<Point3D> &v_in, std::vector<Point3D> &v_out, int &xmin, int &x_max, int &y_min, int &y_max);
		void view_scale_for_3d_v2(std::vector<Point3D> &v_all, int &xmin, int &x_max, int &y_min, int &y_max);
		void view_scale_for_3d_v3(std::vector<Point3D> &v_all, std::vector<Point3D> &v_in, std::vector<Point3D> &v_out);
		void calculate_3d_coordinate(std::vector<Point3D> &v_all, std::vector<Point3D> &show_coord_3d);
		void line_view_update_3d();

		void single_slot_seperate_edges_v101(stl_data &block_a, std::vector<Point3D> &open_slot_pts, std::vector<Point3D> &v1, std::vector<Point3D> &v2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2, std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4);
		void single_slot_seperate_edges_v102(stl_data &block_a, std::vector<Point3D> &open_slot_pts, std::vector<Point3D> &v1, std::vector<Point3D> &v2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2, std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_e3, std::vector<Point3D> &slot_e4);
		void single_slot_seperate_edges_v103(stl_data &block_a, std::vector<Point3D> &v_in, std::vector<Point3D> &v1, std::vector<Point3D> &v2);
		
		void interpolate_edge_if_necessary(std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1);
		void interpolate_edge_if_necessary_v2(std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1);
		void interpolate_edge_if_necessary_v3(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1);
		void interpolate_edge_if_necessary_v4(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1);
		void interpolate_edge_if_necessary_v5(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1);
		void interpolate_edge_if_necessary_v6(std::vector<Point3D> &slot_e1, std::vector<Point3D> &slot_e2, std::vector<Point3D> &slot_ip1, std::vector<Point3D> &interpolated_edge_v1, float tol);
		void interpolate_edge_if_necessary_v7(std::vector<Point3D> &slot_e1, std::vector<Point3D> &interpolated_edge_v1, float tol);
		void interpolate_edge_if_necessary_v8(std::vector<Point3D> &slot_e1, std::vector<Point3D> &interpolated_edge_v1, float tol);

		float find_rotation();
		void rotate_translate_point3d(std::vector<cad3d::Point3D> &vi, std::vector<cad3d::Point3D> &vo, float theta);
		void align_with_robot_CAD_local_coords(std::vector<Point3D> &v_in, std::vector<Point3D> &v_out);

		void generate_border_lines(stl_data &block_a, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, float z_level);
		void cartesian_separate_external_regions_v1(std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_int_l1, std::vector<Point3D> &v_int_l2, int cp_mode, std::vector<Point3D> &slot_plane1);
		void generate_border_line_polar(stl_data &block_a, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2);
		void polar_separate_external_regions_v1(std::vector<Point3D> &feat1_e1, std::vector<Point3D> &feat1_e2, std::vector<Point3D> &v_ext_l1, std::vector<Point3D> &v_ext_l2, int cp_mode, std::vector<Point3D> &slot_plane1);

		void find_close_points_from_edge(std::vector<Point3D> &v_in, std::vector<Point3D> &v_ouy, Point3D p_base, float d_tol);
		void get_straight_line_2d(float &m, float &c, Point3D tan1, Point3D tan2, Point3D norm0);
		Point3D get_intersect_point_phase3_v1(Point3D norm0, float m1, float c1, float m2, float c2);
		void generate_single_cross_along_centreline_v2(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &c_line, std::vector<cad3d::Point3D> &normal_intersect);
		void generate_single_cross_along_centreline_v3(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &c_line, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2);
		void generate_single_cross_along_centreline_v4(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &c_line, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2);
		
		void generate_single_cross_along_centreline_v5(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &normal_intersect, std::vector<Point3D> &c_line, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, float tol);
		void generate_single_cross_along_centreline_v6(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &normal_intersect, std::vector<Point3D> &c_line, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, float tol, int mode);
		void generate_single_cross_along_centreline_v7(std::vector<Point3D> &edge1, std::vector<Point3D> &edge2, std::vector<Point3D> &normal_intersect, std::vector<Point3D> &c_line, std::vector<Point3D> &intersect_line1, std::vector<Point3D> &intersect_line2, float tol, int mode);
		
		void multiple_layer_cross_section_extract(stl_data &block_a, std::vector<Point3D> &sub_vol_points, std::vector<Point3D> &c_line, std::vector<Point3D> &normal_intersect, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2);
		void cross_section_outline_generate_v1(std::vector<Point3D> &vin1, std::vector<Point3D> &vin2, std::vector<Point3D> &vout1, std::vector<Point3D> &vout2);

		void fine_3D_segment_extraction_polar(std::vector<Point3D> &c_line, int i_cent, std::vector<Point3D> &slot_e1);
		void generate_3D_sub_layer_v1(std::vector<std::vector<Point3D>> &sub_vol_main_slices, std::vector<Point3D> &boundary_3d_line1, std::vector<Point3D> &boundary_3d_line2, int mode);
		void generate_3D_sub_layer_v2(std::vector<std::vector<Point3D>> &sub_vol_main_slices, std::vector<std::vector<Point3D>> &sub_vol_sub_slices, std::vector<Point3D> &boundary_3d_line1, std::vector<Point3D> &boundary_3d_line2, int mode);

		void generate_3D_fine_slice_v1(std::vector<std::vector<Point3D>> &sub_vol_main_slices, std::vector<std::vector<Point3D>> &sub_vol_sub_slices, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, std::vector<Point3D> &slice3d_l1, std::vector<Point3D> &slice3d_l2, int mode);
		void generate_3D_fine_slice_v2(std::vector<std::vector<Point3D>> &sub_vol_slices, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, std::vector<Point3D> &slice3d_l1, std::vector<Point3D> &slice3d_l2, int mode);
		void generate_3D_fine_slice_v3(std::vector<std::vector<Point3D>> &sub_vol_slices, std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, std::vector<Point3D> &slice3d_l1, std::vector<Point3D> &slice3d_l2, int mode);
		void generate_extended_projection_line_v1(std::vector<Point3D> &intersect_l1, std::vector<Point3D> &intersect_l2, std::vector<Point3D> &ext_intersect_l1, std::vector<Point3D> &ext_intersect_l2, float stretch_val);
	};
}