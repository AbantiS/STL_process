//#include "pch.h"
#include "stdafx.h"

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include <sstream>
//#include "CSVFile.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <list>
#include <numeric>

#include "D:\code_works\cad_i17_versions\v9_4_dell\cad_i17\cad_i17\blockcad.h"
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "spline.h"

using namespace std;
using namespace cadblock;
using namespace Eigen;
using Eigen::MatrixXd;
//using Eigen::Rotation2Df;

struct Point3D
{
	float x, y, z; //in cartesion if polar its theta,r,z
				   //int i;//index in the original point cloud
	int r1, r2; // repetition status: 0,1,2
				//int orient;// 3 point orientatin measurement
};
// A global point needed for  sorting points with reference 
// to the first point. Used in compare function of qsort() 
Point3D p0;
// A global point3d vector for plotting with void display function
vector<Point3D> target_plot, target_plot2, targetl1, targetl2;
int mat_x, mat_y;
float parse_float(std::ifstream& s) {
	char f_buf[sizeof(float)];
	s.read(f_buf, 4);
	float* fptr = (float*)f_buf;
	return *fptr;
}
Point3D parse_point(std::ifstream& s) {
	Point3D a;
	a.x = parse_float(s);
	a.y = parse_float(s);
	a.z = parse_float(s);
	return a;
}
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
	vector<Point3D> vertices;

	stl_data(std::string namep) : name(namep) {}
};
struct circle2d
{
	float x, y; // center coordinates
	float r; //radius
};

float cad3d::parse_float(std::ifstream& s) {
	char f_buf[sizeof(float)];
	s.read(f_buf, 4);
	float* fptr = (float*)f_buf;
	return *fptr;
}
cad3d::Point3D cad3d::parse_point(std::ifstream& s) {
	Point3D a;
	a.x = cad3d::parse_float(s);
	a.y = cad3d::parse_float(s);
	a.z = cad3d::parse_float(s);
	return a;
}
void cad3d::stl_reader(const std::string& f, stl_data& info) {
	ifstream stl_file(f.c_str(), std::ios::in | std::ios::binary);
	//ifstream stl_file(f,ios::binary);
	std::cout << "opened" << endl;
	if (!stl_file) {
		std::cout << "ERROR: COULD NOT READ FILE" << std::endl;
		assert(false);
	}
	char header_info[80] = "";
	char n_triangles[4];
	stl_file.read(header_info, 80);
	//cout << header_info << endl;
	stl_file.read(n_triangles, 4);
	std::string h(header_info);
	//stl_data info(h);
	//cout << info.name;

	unsigned int* r = (unsigned int*)n_triangles;
	unsigned int num_triangles = *r;
	char dummy[2];
	std::cout << num_triangles << endl;
	for (unsigned int i = 0; i < num_triangles; i++) {
		auto normal = cad3d::parse_point(stl_file);
		auto v1 = cad3d::parse_point(stl_file);
		auto v2 = cad3d::parse_point(stl_file);
		auto v3 = cad3d::parse_point(stl_file);
		info.triangles.push_back(triangle(normal, v1, v2, v3));
		//char dummy[2];
		info.vertices.push_back(v1);
		info.vertices.push_back(v2);
		info.vertices.push_back(v3);
		stl_file.read(dummy, 2);
	}

	std::cout << info.vertices.size() << endl;
	stl_file.close();
	//return info;
	//https://github.com/dillonhuff/stl_parser/blob/master/main.cpp
}
void cad3d::point3d_write(string &f, vector<Point3D> &v1) {
	ofstream vout(f);
	for (int i = 0; i < v1.size(); i++) {
		vout << v1[i].x << " " << v1[i].y << " " << v1[i].z << endl;
	}
	vout.close();
}

void cad3d::stl_catia_reader(const std::string& f, stl_data& info) {
	ifstream stl_file(f.c_str(), std::ios::in | std::ios::binary);
	//ifstream stl_file(f,ios::binary);
	std::cout << "opened" << endl;
	if (!stl_file) {
		std::cout << "ERROR: COULD NOT READ FILE" << std::endl;
		assert(false);
	}
	std::string line;
	char stl_start[] = "solid CATIA";
	char normal_start[] = "facet normal";
	char *output_facetnormal = NULL;
	char *output_outerloop = NULL;
	char *output_endloop = NULL;
	char *output_endfacet = NULL;
	int dummy;
	char dummy_read[2];
	char line_char[80];
	char * pch;
	double a, b, c;
	string dd1, dd2;
	string token;
	//std::getline(stl_file, line,' ');
	vector<string> tokens;
	Point3D n1, p1, p2, p3;
	line = "";
	int triangle_no = 0;
	//for (int i = 0; i < 5; i++) {
	for (int i = 0; i < 5000; i++) {
		//while (std::getline(stl_file, line)) {
		std::getline(stl_file, line);
		// line contains the current line
		std::copy(line.begin(), line.end(), line_char);
		output_facetnormal = strstr(line_char, normal_start);
		if (!stl_file.eof()) {
			triangle_no++;
			//std::cout << triangle_no << ": " << endl;
			if ((i != 0) & (output_facetnormal > 0)) {
				n1 = split_char_array(line, ' ', tokens);//j = 0;
				line = "";
				std::getline(stl_file, line); //j = 1;
				//split_char_array(line, ' ', tokens); 
				line = "";
				std::getline(stl_file, line); //j = 2;v1
				p1 = split_char_array(line, ' ', tokens);
				line = "";
				std::getline(stl_file, line); //j = 3;v2
				p2 = split_char_array(line, ' ', tokens);
				line = "";
				std::getline(stl_file, line); //j = 4;v3
				p3 = split_char_array(line, ' ', tokens);
				line = "";
				std::getline(stl_file, line); //j = 5;
				line = "";
				std::getline(stl_file, line); //j = 6;
				line = "";
				info.triangles.push_back(triangle(n1, p1, p2, p3));
				//char dummy[2];
				info.vertices.push_back(p1);
				info.vertices.push_back(p2);
				info.vertices.push_back(p3);
			}
		}
		
		dummy = 0;
		//output_start = NULL;
		line = "";
	
	}
	std::cout << "triangle_no:" << triangle_no <<endl;
	//return info;
	//https://github.com/dillonhuff/stl_parser/blob/master/main.cpp
}
cad3d::Point3D cad3d::split_char_array(const std::string& s, char delimiter, std::vector<std::string> &tokens) {
	//std::vector<std::string> tokens;
	tokens.resize(0);
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter))
	{
		tokens.push_back(token);
		std::cout << token<<", ";
	}
	std::cout << " splitted ";
	Point3D p;
	p.z = stof(tokens[tokens.size() - 1]);
	//p.y = stof(tokens[tokens.size() - 3]);
	//p.x = 2.0;
	
	//
	
	int jx, jy, jz;
	
	if (p.z >= 0) {
		jy = 3;//<< endl;}
	}
	else if (p.z < 0)
	{
		jy = 2;
	}
	p.y = stof(tokens[tokens.size() - jy]);
	
	if (p.y >= 0) {
		p.x =  stof(tokens[tokens.size() - (jy+2)]);//<< endl;}
	}
	else if(p.y < 0)
	{
		p.x =  stof(tokens[tokens.size() - (jy + 1)]);
	}
	/*std::cout << p.x << "* ";//px
	std::cout << stof(tokens[tokens.size() - jy]) << "~ ";//py
	std::cout << stof(tokens[tokens.size() - 1]) << ", ";//pz
	std::cout << endl;*/
	return p;
}

int cad3d::find_point3d_max(vector<cad3d::Point3D> &v1, int dim) { // version 2 with both direction option
																   //return the first index of maximum value
	int max_index = 0;
	float iValue = 0.0;
	float newValue = 0.0;
	if (dim == 1) //x values 
	{

		iValue = v1[0].x;// an arbtrary value
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].x;
			//cout << newValue << endl;
			if (newValue > iValue) {
				max_index = i;
				iValue = newValue;

			}
		}


	}
	if (dim == 2) {//if search in y dimension
				   //for (int i = 0; i < 100; i++) {
		iValue = v1[0].y;
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].y;
			//cout << newValue << endl;
			if (newValue > iValue) {
				max_index = i;
				iValue = newValue;

			}
		}
	}
	if (dim == 3) {//if search in y dimension
				   //for (int i = 0; i < 100; i++) {
		iValue = v1[0].z;
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].z;
			//cout << newValue << endl;
			if (newValue > iValue) {
				max_index = i;
				iValue = newValue;

			}
		}
	}
	//cout << max_index << endl;
	return max_index;
}
int cad3d::find_point3d_min(vector<cad3d::Point3D> &v1, int dim) { // version 2 with both direction option
																   //return the first index of maximum value
	int min_index = 0;
	float iValue = 0.0;
	float newValue = 0.0;
	if (dim == 1) //x values 
	{

		iValue = v1[0].x;// an arbtrary value
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].x;
			//cout << newValue << endl;
			if (newValue < iValue) {
				min_index = i;
				iValue = newValue;

			}
		}


	}
	if (dim == 2) {//if search in y dimension
				   //for (int i = 0; i < 100; i++) {
		iValue = v1[0].y;
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].y;
			//cout << newValue << endl;
			if (newValue < iValue) {
				min_index = i;
				iValue = newValue;

			}
		}
	}
	if (dim == 3) {//if search in y dimension
				   //for (int i = 0; i < 100; i++) {
		iValue = v1[0].z;
		for (int i = 0; i < v1.size(); i++) {
			newValue = v1[i].z;
			//cout << newValue << endl;
			if (newValue < iValue) {
				min_index = i;
				iValue = newValue;

			}
		}
	}
	///cout << max_index << endl;
	return min_index;
}
void cad3d::find_point_3d_1feature(vector<Point3D> &vo, vector<Point3D> &vi, float test_value, float tol, int dim) {
	//performs frequency measurement within a certain tolerance 
	// n being thetotal number of rows

	float newValue;
	int frequency = 0;
	if (dim == 1) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].x;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vi[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}
	if (dim == 2) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].y;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vi[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}
	if (dim == 3) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].z;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vi[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}

	//std::cout << frequency << endl;
}
void cad3d::find_point3d_i_1feature(vector<int> &vo, vector<Point3D> &vi, float test_value, float tol, int dim) {
	//this one returns an index file of integers
	//used furthur for 2 feature comparison 

	float newValue;
	int frequency = 0;

	if (dim == 1) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].x;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo[i] = 1;
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			else
				vo[i] = 0;
			// At the end of this loop the value you want is indexToReturn
		}
	}

	if (dim == 2) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].y;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo[i] = 1;
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			else
				vo[i] = 0;
			// At the end of this loop the value you want is indexToReturn
		}
	}

	if (dim == 3) { // for z dimension
		for (int i = 0; i < vi.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vi[i].z;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo[i] = 1;
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			else
				vo[i] = 0;
			// At the end of this loop the value you want is indexToReturn
		}
	}

	//std::cout << frequency << endl;
}
int cad3d::find_1d_vector_max(vector<float> &v1, int dim) { // version 2 with both direction option
													 //return the first index of maximum value
	int max_index = 0;
	float iValue;// = 0.0;
	float newValue = 0.0;
	if (dim == 1) //x values 
	{
		iValue = v1[0];// an arbtrary value
		for (int i = 1; i < v1.size(); i++) {
			newValue = v1[i];
			//cout << newValue << endl;
			if (newValue > iValue) {
				max_index = i;
				iValue = newValue;

			}
		}
	}
	return max_index;
}
int cad3d::find_1d_vector_min(vector<float> &v1, int dim) { // version 2 with both direction option
													 //return the first index of maximum value
	int min_index = 0;
	float iValue;// = 0.0;
	float newValue = 0.0;
	if (dim == 1) //x values 
	{
		iValue = v1[0];// an arbtrary value
		for (int i = 1; i < v1.size(); i++) {
			newValue = v1[i];
			//cout << newValue << endl;
			if (newValue < iValue) {
				min_index = i;
				iValue = newValue;

			}
		}
	}
	return min_index;
}

float cad3d::dist201(Point3D p1, Point3D p2) // added by abanti
{
	return (p1.x - p2.x)*(p1.x - p2.x) +
		(p1.y - p2.y)*(p1.y - p2.y);
}
float cad3d::dist301(Point3D p1, Point3D p2) // added by abanti
{
	float res = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	return res;
}
float cad3d::del_theta(Point3D pol1, Point3D pol2) // works on polar angle difference
{
	//return (sqrt((pol2.y - pol1.y)^2));
	return ((pol2.y - pol1.y)*(pol2.y - pol1.y));
}
void cad3d::find_Point3D_index(Point3D ref, vector<Point3D> &vall, vector<int> &vi_used, vector<int> &v_temp) {
	v_temp.resize(0);
	float ds;
	for (int i = 1; i < vall.size(); i++) {
		ds = dist201(ref, vall[i]);
		if (ds == 00) {
			v_temp.push_back(i);
			//std::cout << i << ", ";
			vi_used[i] = vi_used[i] + 1;
		}

	}
	//std::cout << endl;
}
void cad3d::find_Point3D_indexv2(Point3D ref, vector<Point3D> &vall, vector<int> &v_temp) {
	v_temp.resize(0);
	float ds;
	for (int i = 1; i < vall.size(); i++) {
		ds = dist201(ref, vall[i]);
		if (ds == 00) {
			v_temp.push_back(i);
			//index of all the points who have 0 distance with the concerned point
			//std::cout << i << ", ";
			//vi_used[i] = vi_used[i] + 1;
		}

	}
	//std::cout << endl;
}
void cad3d::find_Point3D_indexv3(Point3D ref, vector<Point3D> &vall, vector<int> &vi_used, vector<int> &v_temp) {
	v_temp.resize(0);
	float ds;
	for (int i = 0; i < vall.size(); i++) {
		ds = dist201(ref, vall[i]);
		if (ds < 0.1) {
			v_temp.push_back(i);
			//index of all the points who have 0 distance with the concerned point
			//std::cout << i << ", ";
			//vi_used[i] = vi_used[i] + 1;
		}

	}
	//std::cout << endl;
}

void cad3d::exclude_2i_from_1i(vector<int> &v1, vector<int> &v2, vector<int> &vout) {
	//takes 2 input int vector and returns one without 2i
	// first copy and then remove;
	// can be used cumulatively
	for (int i = 0; i < v1.size(); i++) {
		if ((v1[i] == 1) & v2[i] == 1)
			vout[i] = 0;

		//use of else makes it noncumulative
	}
}
void cad3d::initialize_0i(vector<int> &v) {
	for (int i = 0; i < v.size(); i++) {
		v[i] = 0;
	}
}
void cad3d::show_index_frequency(vector<int> &v) {
	// takes an index vector with 1 and 0, 
	//if 1 incrmenets frequency and shows in the end
	int frequency = 0.0;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] == 1)
			frequency++;
	}
	std::cout << "frequency: " << frequency << endl;
}
void cad3d::i_2_point3d(vector<int> &v_i, vector<Point3D> &v_all, vector<Point3D> &vo) {
	for (int i = 0; i < v_all.size(); i++) {
		if (v_i[i] == 1)
			vo.push_back(v_all[i]);
	}
	//
}

void cad3d::centre_at_x0y0(vector<Point3D> &vin, vector<Point3D> &vout) {

	int xmaxi = find_point3d_max(vin, 1);
	int xmini = find_point3d_min(vin, 1);
	int ymaxi = find_point3d_max(vin, 2);
	int ymini = find_point3d_min(vin, 2);
	//int zmaxi = find_point3d_max(vin, 1);
	//int zmini = find_point3d_min(vin, 1);
	float tol = 0.5;
	float xmax = vin[xmaxi].x;
	float xmin = vin[xmini].x;
	float ymax = vin[ymaxi].y;
	float ymin = vin[ymini].y;

	float center_x = (xmin + xmax) / 2;
	float center_y = (ymin + ymax) / 2;
	vout.resize(vin.size());
	for (int i = 0; i < vin.size(); i++) {
		vout[i].x = vin[i].x - center_x;
		vout[i].y = vin[i].y - center_y;
		vout[i].z = vin[i].z;
	}
	cout << "2D shifted to Cartesean (0,0)" << endl;

	//window set up 
	/*plot_win_max_x = vout[xmaxi].x + 10;
	plot_win_min_x = vout[xmini].x - 10;
	plot_win_max_y = vout[ymaxi].y + 10;
	plot_win_min_y = vout[ymini].y - 10;*/
}
void cad3d::triangle_centre_at_x0y0(vector<Point3D> &vin, vector<triangle> &tin, vector<triangle> &tout) {

	int xmaxi = find_point3d_max(vin, 1);
	int xmini = find_point3d_min(vin, 1);
	int ymaxi = find_point3d_max(vin, 2);
	int ymini = find_point3d_min(vin, 2);
	//int zmaxi = find_point3d_max(vin, 1);
	//int zmini = find_point3d_min(vin, 1);
	float tol = 0.5;
	float xmax = vin[xmaxi].x;
	float xmin = vin[xmini].x;
	float ymax = vin[ymaxi].y;
	float ymin = vin[ymini].y;

	float center_x = (xmin + xmax) / 2;
	float center_y = (ymin + ymax) / 2;
	//tout.resize(tin.size());
	Point3D ptemp;
	ptemp.x = 0;
	ptemp.y = 0;
	ptemp.z = 0;
	triangle t_temp(ptemp,ptemp,ptemp,ptemp);
	Point3D p1, p2, p3;
	for (int i = 0; i < tin.size(); i++) {
	//for (int i = 0; i < 5; i++) {
		p1 = copy_point3d(tin[i].v1);
		p2 = copy_point3d(tin[i].v2);
		p3 = copy_point3d(tin[i].v3);
		//generate_triangle(Point3D p1, Point3D p2, Point3D p3);
		
		//std::cout << "p1: "<<p1.x <<", "<<p1.y<< "p2 " <<p2.x << endl;
		
		p1.x = tin[i].v1.x - center_x;
		p1.y = tin[i].v1.y - center_y;
		//tout[i].v1.z = tin[i].v1.z;

		p2.x = tin[i].v2.x - center_x;
		p2.y = tin[i].v2.y - center_y;
		//tout[i].v2.z = tin[i].v2.z;

		p3.x = tin[i].v3.x - center_x;
		p3.y = tin[i].v3.y - center_y;
		//tout[i].v3.z = tin[i].v3.z;

		t_temp = generate_triangle(tin[i].v1, tin[i].v2, tin[i].v3);
		tout.push_back(t_temp);
		//std::cout << p1.x << ", "<< t_temp.v2.x <<endl;
	}
	//cout << "2D shifted to Cartesean (0,0)" << endl;

	//window set up 
	/*plot_win_max_x = vout[xmaxi].x + 10;
	plot_win_min_x = vout[xmini].x - 10;
	plot_win_max_y = vout[ymaxi].y + 10;
	plot_win_min_y = vout[ymini].y - 10;*/
}
void cad3d::find_zrange(vector<Point3D> &vin, vector<float> &z_range) {
	int xmaxi = find_point3d_max(vin, 1);
	int ymaxi = find_point3d_max(vin, 2);
	//int zmaxi = find_point3d_max(vin, 1);
	//int zmini = find_point3d_min(vin, 1);
	float tol = 0.5;
	float xmax = vin[xmaxi].x;
	float ymax = vin[ymaxi].y;

	int zcount = 0;
	for (int i = 0; i < vin.size(); i++) {
		if ((vin[i].x > (xmax - tol)) || (vin[i].y > (ymax - tol))) {
			zcount++;
			z_range.push_back(vin[i].z);
			//vout.push_back(vin[i]);
		}
	}
	int zmaxvi = find_1d_vector_max(z_range, 1);
	int zminvi = find_1d_vector_min(z_range, 1);
	cout << "z range: " << z_range.size() <<
		", " << z_range[zmaxvi] << ", " << z_range[zminvi] << endl;
}

void convert_cart2pol(vector<Point3D> &vin, vector<Point3D> &vout) {
	const float pi = 3.14159265;
	vout.resize(vin.size());
	for (int i = 0; i < vin.size(); i++) {
		vout[i].x = sqrt(pow(vin[i].x, 2) + pow(vin[i].y, 2));
		vout[i].y = atan2(vin[i].y, vin[i].x) * 180 / pi;;
		vout[i].z = vin[i].z;
	}
}
void convert_cart2polPoint3D(Point3D &vin, Point3D &vout) {
	const float pi = 3.14159265;

	vout.x = sqrt(pow(vin.x, 2) + pow(vin.y, 2));
	vout.y = atan2(vin.y, vin.x) * 180 / pi;;
	vout.z = vin.z;

}
void convert_pol2cart(vector<Point3D> &vin, vector<Point3D> &vout) {
	vout.resize(vin.size());
	const float pi = 3.14159265;
	for (int i = 0; i < vin.size(); i++) {
		vout[i].x = vin[i].x*cos(vin[i].y*pi / 180);
		vout[i].y = vin[i].x*sin(vin[i].y*pi / 180);
		vout[i].z = vin[i].z;
	}
}
void convert_pol2cartPoint3D(Point3D &vin, Point3D &vout) {
	//vout.resize(vin.size());
	const float pi = 3.14159265;
	//for (int i = 0; i < vin.size(); i++) {
	vout.x = vin.x*cos(vin.y*pi / 180);
	vout.y = vin.x*sin(vin.y*pi / 180);
	vout.z = vin.z;
	//}
}
void find_point_3d_1feature_polcart(vector<Point3D> &vo, vector<Point3D> &vipol, vector<Point3D> &vicart, float test_value, float tol, int dim) {
	//performs frequency measurement within a certain tolerance 
	// n being thetotal number of rows

	float newValue;
	int frequency = 0;

	if (dim == 3) { // for z dimension
		for (int i = 0; i < vipol.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vipol[i].z;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vicart[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}
	if (dim == 1) { // for x dimension
		for (int i = 0; i < vipol.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vipol[i].x;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vicart[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}
	if (dim == 2) { // for y dimension
		for (int i = 0; i < vipol.size(); i++) {
			//for (int i = 0; i < 10; i++) {
			newValue = vipol[i].y;
			//cout << newValue << endl;
			if (newValue >= (test_value - tol) && newValue <= (test_value + tol))
			{
				vo.push_back(vicart[i]);
				//mat1[i] = 1;// stores the index number of vertex having required features
				frequency++;
			}
			// At the end of this loop the value you want is indexToReturn
		}
	}

	//std::cout << frequency << endl;
}
int find_2feature_occurance_polcart_v1(vector<Point3D> &vpol, vector<Point3D> &vcart, float test_valuec, float test_valuep, float tol, int dim1, int dim2) {
	// this one will be sequence sensitive
	// which satisfies the first feature as well as the 2nd feature
	//both vectors should be of same length
	bool a, b;
	int i_count = 0;
	int j_send = 0;
	if (dim1 == 1 & dim2 == 1) //x value and 
	{
		for (int i = 0; i < vpol.size(); i++) {
			a = (vpol[i].x < test_valuep + tol) & (vpol[i].x > test_valuep - tol);
			b = (vcart[i].x < test_valuec + tol) & (vcart[i].x > test_valuec - tol);
			if (a == 1 & b == 1) {
				cout << i << endl;
				i_count++;
				j_send = i;
			}
		}

	}
	if (i_count == 1)
		return j_send;
	else
		return 7000;//change to a number thant does not exit

}

void cad3d::triangle_vertexindex_init(vector<triangle> &t_all, vector<triangle_index> &ti_all) {
	ti_all.resize(t_all.size());
	for (int i = 0; i < t_all.size(); i++) {
		ti_all[i].ti = i;
		ti_all[i].v1i = 3 * i;
		ti_all[i].v2i = 3 * i + 1;
		ti_all[i].v3i = 3 * i + 2;
	}
}
void cad3d::triangle_allpoint_1feature(vector<triangle> &t_all, vector<triangle> &tout, float ref, float tol, int dim) {
	int a, b, c;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < t_all.size(); i++) {
		if (dim == 1) {
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))a = 1;
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))b = 1;
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))c = 1;
			if (a&b&c) tout.push_back(t_all[i]);
		}
		if (dim == 2) {
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))a = 1;
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))b = 1;
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))c = 1;
			if (a&b&c) tout.push_back(t_all[i]);
		}
		if (dim == 3) {
			if ((t_all[i].v1.z <= ref + tol) || (t_all[i].v1.z > ref - tol))a = 1;
			if ((t_all[i].v1.z <= ref + tol) || (t_all[i].v1.z > ref - tol))b = 1;
			if ((t_all[i].v1.z <= ref + tol) || (t_all[i].v1.z > ref - tol))c = 1;
			if (a&b&c) tout.push_back(t_all[i]);
		}

	}

}
void cad3d::triangle_allpoint_i1feature(vector<Point3D> &p_all, vector<triangle_index> &ti_all, vector<triangle_index> &ti_out, float ref, float tol, int dim) {
	//reas the vertex index from triangle index file
	//checks if all*** the points fullfill the criteria
	//if comples, pushes the index to tout index file
	int a, b, c, i1, i2, i3;
	Point3D p1, p2, p3;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < ti_all.size(); i++) {
		//for (int i = 0; i < 20; i++) {

		i1 = ti_all[i].v1i; i2 = ti_all[i].v2i; i3 = ti_all[i].v3i;
		p1 = p_all[i1]; p2 = p_all[i2]; p3 = p_all[i3];

		if (dim == 1) {
			if ((p1.x <= ref + tol) & (p1.x > ref - tol))a = 1;
			if ((p2.x <= ref + tol) & (p2.x > ref - tol))b = 1;
			if ((p3.x <= ref + tol) & (p3.x > ref - tol))c = 1;
			if (a&b&c) { ti_out.push_back(ti_all[i]); }
		}
		if (dim == 2) {
			if ((p1.y <= ref + tol) & (p1.y > ref - tol))a = 1;
			if ((p2.y <= ref + tol) & (p2.y > ref - tol))b = 1;
			if ((p3.y <= ref + tol) & (p3.y > ref - tol))c = 1;
			if (a&b&c) { ti_out.push_back(ti_all[i]); }
		}
		if (dim == 3) {
			if ((p1.z <= ref + tol) & (p1.z > ref - tol))a = 1;
			if ((p2.z <= ref + tol) & (p2.z > ref - tol))b = 1;
			if ((p3.z <= ref + tol) & (p3.z > ref - tol))c = 1;
			if (a&b&c) {
				ti_out.push_back(ti_all[i]);
				//std::cout << i << ": " << p1.z << endl;
			}
		}
		a = 0; b = 0; c = 0;
	}

}
void cad3d::triangle2line(vector<Point3D> &centered_block, vector<triangle_index> &tslot_i2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<int> &i_line, vector<Point3D> &v_temp) {
	vector<int>it_temp;
	//vector<int>i_line, i_line_nr;
	vector<int>i_line_frequency;// saves frequency if a line is repeated more than once
								//vector<Point3D> v_temp;
	Point3D  ref1, ref2;
	//vector<Point3D> v_line1, v_line2;
	int nt1, nt2, nt3, ilinet, tempfrequency;
	ilinet = 0;
	float dd1, dd2;
	for (int i = 0; i < tslot_i2.size(); i++) {
		//for (int i = 10; i < 50; i++) {
		nt1 = tslot_i2[i].v1i;
		v_temp.push_back(centered_block[nt1]);
		nt2 = tslot_i2[i].v2i;
		v_temp.push_back(centered_block[nt2]);
		nt3 = tslot_i2[i].v3i;
		v_temp.push_back(centered_block[nt3]);
		if (centered_block[nt1].x <= centered_block[nt2].x) {
			v_line1.push_back(centered_block[nt1]); v_line2.push_back(centered_block[nt2]);
		}
		else {
			v_line1.push_back(centered_block[nt2]);
			v_line2.push_back(centered_block[nt1]);
		}
		i_line.push_back(ilinet);
		ilinet++;
		if (centered_block[nt2].x <= centered_block[nt3].x) {
			v_line1.push_back(centered_block[nt2]);
			v_line2.push_back(centered_block[nt3]);
		}
		else {
			v_line1.push_back(centered_block[nt3]);
			v_line2.push_back(centered_block[nt2]);
		}
		i_line.push_back(ilinet);
		ilinet++;
		if (centered_block[nt3].x <= centered_block[nt1].x) {
			v_line1.push_back(centered_block[nt3]);
			v_line2.push_back(centered_block[nt1]);
		}
		else {
			v_line1.push_back(centered_block[nt1]);
			v_line2.push_back(centered_block[nt3]);
		}
		i_line.push_back(ilinet);
		ilinet++;

		/*std::cout << "(" << centered_block[nt1].x << "," << centered_block[nt1].y << ":";
		std::cout << centered_block[nt2].x << "," << centered_block[nt2].y << "),";
		std::cout << "(" << centered_block[nt2].x << "," << centered_block[nt2].y << ":";
		std::cout << centered_block[nt3].x << "," << centered_block[nt3].y << "),";
		std::cout << "(" << centered_block[nt3].x << "," << centered_block[nt3].y << ":";
		std::cout << centered_block[nt1].x << "," << centered_block[nt1].y << ")"<<endl;*/
		//std::cout << centered_block[nt1].x << "," << centered_block[nt2].x << "," << centered_block[nt3].x << endl;
	}
}
void cad3d::triangle_anypoint_1feature(vector<triangle> &t_all, vector<triangle> &tout, float ref, float tol, int dim) {
	int a, b, c;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < t_all.size(); i++) {
		if (dim == 1) {
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))a = 1;
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))b = 1;
			if ((t_all[i].v1.x <= ref + tol) || (t_all[i].v1.x > ref - tol))c = 1;
			if (a&b&c) tout.push_back(t_all[i]);
		}
		if (dim == 2) {
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))a = 1;
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))b = 1;
			if ((t_all[i].v1.y <= ref + tol) || (t_all[i].v1.y > ref - tol))c = 1;
			if (a&b&c) tout.push_back(t_all[i]);
		}
		if (dim == 3) {
			if ((t_all[i].v1.z <= ref + tol) & (t_all[i].v1.z > ref - tol))a = 1;
			if ((t_all[i].v2.z <= ref + tol) & (t_all[i].v2.z > ref - tol))b = 1;
			if ((t_all[i].v3.z <= ref + tol) & (t_all[i].v3.z > ref - tol))c = 1;
			if (a||b||c) tout.push_back(t_all[i]);
		}

	}

}
void cad3d::triangle_anypoint_i1feature(vector<Point3D> &p_all, vector<triangle_index> &ti_all, vector<triangle_index> &ti_out, float ref, float tol, int dim) {
	//reas the vertex index from triangle index file
	//checks if all*** the points fullfill the criteria
	//if comples, pushes the index to tout index file
	int a, b, c, i1, i2, i3;
	Point3D p1, p2, p3;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < ti_all.size(); i++) {
		//for (int i = 0; i < 20; i++) {

		i1 = ti_all[i].v1i; i2 = ti_all[i].v2i; i3 = ti_all[i].v3i;
		p1 = p_all[i1]; p2 = p_all[i2]; p3 = p_all[i3];

		if (dim == 1) {
			if ((p1.x <= ref + tol) & (p1.x > ref - tol))a = 1;
			if ((p2.x <= ref + tol) & (p2.x > ref - tol))b = 1;
			if ((p3.x <= ref + tol) & (p3.x > ref - tol))c = 1;
			if (a&b&c) { ti_out.push_back(ti_all[i]); }
		}
		if (dim == 2) {
			if ((p1.y <= ref + tol) & (p1.y > ref - tol))a = 1;
			if ((p2.y <= ref + tol) & (p2.y > ref - tol))b = 1;
			if ((p3.y <= ref + tol) & (p3.y > ref - tol))c = 1;
			if (a&b&c) { ti_out.push_back(ti_all[i]); }
		}
		if (dim == 3) {
			if ((p1.z <= ref + tol) & (p1.z > ref - tol))a = 1;
			if ((p2.z <= ref + tol) & (p2.z > ref - tol))b = 1;
			if ((p3.z <= ref + tol) & (p3.z > ref - tol))c = 1;
			if (a||b||c) {
				ti_out.push_back(ti_all[i]);
				//std::cout << i << ": " << p1.z << endl;
			}
		}
		a = 0; b = 0; c = 0;
	}

}

void cad3d::convert_cart2pol(vector<Point3D> &vin, vector<Point3D> &vout) {
	const float pi = 3.14159265;
	vout.resize(vin.size());
	for (int i = 0; i < vin.size(); i++) {
		vout[i].x = sqrt(pow(vin[i].x, 2) + pow(vin[i].y, 2));
		vout[i].y = atan2(vin[i].y, vin[i].x) * 180 / pi;;
		vout[i].z = vin[i].z;
	}
}
void cad3d::convert_cart2polPoint3D(Point3D &vin, Point3D &vout) {
	const float pi = 3.14159265;

	vout.x = sqrt(pow(vin.x, 2) + pow(vin.y, 2));
	vout.y = atan2(vin.y, vin.x) * 180 / pi;;
	vout.z = vin.z;

}
void cad3d::convert_pol2cart(vector<Point3D> &vin, vector<Point3D> &vout) {
	vout.resize(vin.size());
	const float pi = 3.14159265;
	for (int i = 0; i < vin.size(); i++) {
		vout[i].x = vin[i].x*cos(vin[i].y*pi / 180);
		vout[i].y = vin[i].x*sin(vin[i].y*pi / 180);
		vout[i].z = vin[i].z;
	}
}
void cad3d::convert_pol2cartPoint3D(Point3D &vin, Point3D &vout) {
	//vout.resize(vin.size());
	const float pi = 3.14159265;
	//for (int i = 0; i < vin.size(); i++) {
	vout.x = vin.x*cos(vin.y*pi / 180);
	vout.y = vin.x*sin(vin.y*pi / 180);
	vout.z = vin.z;
	//}
}

void cad3d::find_external_segmentsv4(vector<Point3D> &v_temp, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<int> &used_i_v_external, vector<triangle_index> &t_slot, vector<int> &viu) {
	//finds the points that have higher radius than the threshold value
	// i.e., external radius  tolerance
	vector<Point3D> vpol_temp, vpol_line1, vpol_line2;
	convert_cart2pol(v_temp, vpol_temp);
	convert_cart2pol(v_line1, vpol_line1);
	convert_cart2pol(v_line2, vpol_line2);
	int i_rho_max = find_point3d_max(vpol_temp, 1);
	bool a, b;
	vector<int> i_temp1;
	/*viu.resize(t_slot.size()*3);
	initialize_0i(viu);*/
	int vt1, vt2;
	//index of vertices to keep track of how many times they have been used
	//for (int i = 0; i < v_line1.size(); i++) {
	for (int i = 0; i < t_slot.size(); i++) {
		//for (int i = 0; i < 500; i++) {
		for (int j = 0; j < 3; j++) {
			a = vpol_line1[i * 3 + j].x > (vpol_temp[i_rho_max].x - 0.5);
			b = vpol_line2[i * 3 + j].x > (vpol_temp[i_rho_max].x - 0.5);
			if (a & b) {
				//if both points are on external radius
				//i_temp1.push_back(i);
				//vl1_ui1.push_back(i);
				v_ext_l1.push_back(v_line1[i * 3 + j]);
				v_ext_l2.push_back(v_line2[i * 3 + j]);
				//viu[i * 3 + j]=
				/*t1 = tslot_i2[i].v1i;vi_used[t1] = 1;
				t1 = tslot_i2[i].v2i;vi_used[t1] = 1;
				t1 = tslot_i2[i].v3i;vi_used[t1] = 1;*/
				if (j == 0) { vt1 = t_slot[i].v1i; vt2 = t_slot[i].v2i; }
				if (j == 1) { vt1 = t_slot[i].v2i; vt2 = t_slot[i].v3i; }
				if (j == 2) { vt1 = t_slot[i].v3i; vt2 = t_slot[i].v1i; }
				viu[vt1] = viu[vt1]++;
				viu[vt2] = viu[vt2]++;
				used_i_v_external.push_back(vt1);
				used_i_v_external.push_back(vt2);
				/*std::cout << "vt1, vt2 //";
				std::cout << vt1 << ", " << vt2 << "// ";*/
			}
		}

	}
	std::cout << "slot points size" << (t_slot.size() * 3) << endl;
}
int cad3d::get_next_segmentv2(vector<Point3D> &vtl1, vector<Point3D> &vtl2, Point3D ref, vector<int> &viu) {
	//getting the external segment
	float dd;//distance
	int j = -1;
	for (int i = 0; i < vtl1.size(); i++) {
		dd = dist201(ref, vtl1[i]);
		//std::cout << i << ": " << dd << endl;

		//if (dd == 0.0) {
		if (dd == 0.0) {
			j = i;
			/*v2tl1.push_back(vtl1[j]);
			v2tl2.push_back(vtl2[j]);*/
			//ref = vtl2[j];
		}

	}
	return j;
}
void cad3d::get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall) {
	int j;
	vector<int> v_temp;
	for (int i = 1; i < vline1.size(); i++) {
		//j = get_next_segmentv1(vtl1, vtl2, ref);
		j = get_next_segmentv2(vline1, vline2, ref, viu2);
		if (vline_i[i] != 1) {
			if (j < 0)break;
			if (j >= 0) {
				//std::cout << j << endl;
				vlineordered1.push_back(vline1[j]);
				find_Point3D_index(vline1[j], vall, vi_used, v_temp);
				vlineordered2.push_back(vline2[j]);
				ref = vline2[j];
				viu2.push_back(vlineusedindex[j * 2]);
				viu2.push_back(vlineusedindex[j * 2 + 1]);
				//seg_ui2.push_back(vl1_ui1[j]);
				//jend = j;
				vline_i[j] = vline_i[j] + 1;
				vselected.push_back(vline1[j]);
			}
		}
	}
}
void cad3d::get_external_partv3(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall) {
	//this one will not save in the used index
	int j;
	vector<int> v_temp;
	for (int i = 1; i < vline1.size(); i++) {
		//j = get_next_segmentv1(vtl1, vtl2, ref);
		j = get_next_segmentv2(vline1, vline2, ref, viu2);
		if (vline_i[i] != 1) {
			if (j < 0)break;
			if (j >= 0) {
				//std::cout << j << endl;
				vlineordered1.push_back(vline1[j]);
				find_Point3D_index(vline1[j], vall, vi_used, v_temp);
				vlineordered2.push_back(vline2[j]);
				ref = vline2[j];
				viu2.push_back(vlineusedindex[j * 2]);
				viu2.push_back(vlineusedindex[j * 2 + 1]);
				//seg_ui2.push_back(vl1_ui1[j]);
				//jend = j;
				vline_i[j] = vline_i[j] + 1;
				vselected.push_back(vline1[j]);
			}
		}
	}
}
/*void get_trianglesegment_wo_r(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	std::cout << "points and points index size: ";
	std::cout << vall.size() << "," << vi_used.size() << endl;
	show_index_frequency(vi_used);
	float d1, d2;
	vector<Point3D>  vt_line1, vt_line2;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist2(vall[iref], v_line1[i]);
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);
			if (ivt1.size()>0) {
				for (int j = 0; j < ivt1.size(); j++) {
					std::cout << "clock:" << ivt1[j] << endl;
					std::cout << ": " << vi_used[ivt1[j]] << endl;
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
					}
				}
			}

		}
		d2 = dist2(vall[iref], v_line2[i]);
		if (d2 == 0) {
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			if (ivt2.size()>0) {
				for (int j = 0; j < ivt2.size(); j++) {
					std::cout << "counter clock:" << ivt2[j];
					std::cout << ": " << vi_used[ivt2[j]] << endl;
					if (vi_used[ivt2[j]] < 1) {
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

		}
	}
	target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.8, 0.0, 0.0, vt_line1.size());
}*/
void cad3d::get_segment_numberv5(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points
	/*vector<Point3D> vpol_line1, vpol_line2;
	convert_cart2pol(v_line1, vpol_line1);
	convert_cart2pol(v_line2, vpol_line2);*/
	int segment_no = 0;
	vector<int> viu2;
	//vector<int> vi_all;

	vector<int> vline_i;//index vector for line count?
	vline_i.resize(vline1.size());
	initialize_0i(vline_i);

	Point3D ref;
	vector<int> vi_temp;
	vector<Point3D> v_temp;//used for?

	int iymax2 = find_point3d_max(vline1, 2);
	//finding highest , left most point, for initializing count
	vlineordered1.push_back(vline1[iymax2]);
	vlineordered2.push_back(vline2[iymax2]);
	find_Point3D_index(vline1[iymax2], vall, vi_all, vi_temp);// for marking use of points
	std::cout << endl;
	vselected.push_back(vline1[iymax2]);
	viu2.push_back(vlineusedindex[iymax2 * 2]);
	viu2.push_back(vlineusedindex[iymax2 * 2 + 1]);
	vline_i[iymax2] = vline_i[iymax2] + 1;
	if (vline1.size() == (vlineusedindex.size() / 2)) {
		std::cout << "index length matched" << endl;
	}
	//seg_ui2.push_back(vl1_ui1[iymax2]);

	float dd;
	ref = vline2[iymax2];
	int j, jstart, jend;
	jstart = iymax2;

	int i1;
	for (int k = 0; k < 15; k++) {
		//for (int k = 0; ; k++) {
		//get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, ...
		//vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, ...
		//vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall)
		std::cout << "k-" << k << ": ";
		get_external_partv2(vline1, vline2, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
		//find next point with  minimum distance and line repitation within external segment
		find_Point3D_index(ref, vall, vi_all, vi_temp);
		// this is the last point/ slot corner start 
		std::cout << "corner indices: ";
		for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
			v_temp.push_back(vall[vi_temp[ii2]]);
			std::cout << vi_temp[ii2] << ",";
		}
		std::cout << endl;
		//std::cout << "slot corner vertex repitation:" << vi_temp.size() << endl;
		/*target_plot_update(v_temp);
		target_plot_color_update(0.8, 0.0, 0.0, v_temp.size());*/
		v_temp.resize(0);

		//get the slot border inside


		//look for next external next external
		int next_index = -1;
		float iValue;// = 0.0;
		float newValue = 20000.0;
		for (int i = 0; i < vline1.size(); i++) {
			iValue = dist201(ref, vline1[i]);
			if ((iValue < newValue) & (vline_i[i] != 1)) {
				next_index = i;
				newValue = iValue;
			}
		}
		if (next_index > 0) {
			vselected.push_back(vline1[next_index]);
			//vselected.push_back(vline2[next_index]);
			vline_i[next_index] = vline_i[next_index] + 1;
			vlineordered1.push_back(vline1[next_index]);
			vlineordered2.push_back(vline2[next_index]);
			ref = vline2[next_index];

		}
		else if (next_index < 0) break;
		segment_no++;
	}
	std::cout << "segment number: " << segment_no << endl;
	int vt;
	std::cout << "total points: " << vi_all.size() << endl;
	std::cout << "viu2 size: " << viu2.size() << endl;
	std::cout << "v2tl1,2 size: " << vlineordered1.size() << endl;
	/*for (int i = 0; i < vline_i.size(); i++) {
	std::cout << i << ": " << vline_i[i] << endl;
	}*/
	for (int i = 0; i < vlineusedindex.size(); i++) {
		vt = vlineusedindex[i];
		vi_all[vt] = vi_all[vt] + 1;
		//std::cout << i << ": " << vt << endl;
	}

	//find next point
	//get the corner of the slot, which is ref, find all the segments with ref
	std::cout << "vi_used: " << viu2.size() << ", " << vlineordered1.size() << endl;
}
void cad3d::get_segment_numberv501(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new;
	convert_cart2pol(vline1, vpol_line1);
	convert_cart2pol(vline2, vpol_line2);
	for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}
	int segment_no = 0;
	vector<int> viu2;
	//vector<int> vi_all;

	vector<int> vline_i;//index vector for line count?
	vline_i.resize(vline1.size());
	initialize_0i(vline_i);

	Point3D ref;
	vector<int> vi_temp;
	vector<Point3D> v_temp;//used for?

	int iymax2 = find_point3d_max(vline1, 2);
	//finding highest , left most point, for initializing count
	vlineordered1.push_back(vline1[iymax2]);
	vlineordered2.push_back(vline2[iymax2]);
	std::cout << "starting angles: " << vpol_line1[iymax2].y;
	std::cout << " " << vpol_line2[iymax2].y << endl;
	find_Point3D_index(vline1[iymax2], vall, vi_all, vi_temp);// for marking use of points
	std::cout << endl;
	vselected.push_back(vline1[iymax2]);
	viu2.push_back(vlineusedindex[iymax2 * 2]);// why?
	viu2.push_back(vlineusedindex[iymax2 * 2 + 1]);
	vline_i[iymax2] = vline_i[iymax2] + 1;
	if (vline1.size() == (vlineusedindex.size() / 2)) {
		std::cout << "index length matched" << endl;
	}
	//seg_ui2.push_back(vl1_ui1[iymax2]);

	float dd;
	ref = vline2[iymax2];
	int j, jstart, jend;
	jstart = iymax2;
	/*string f = "D:/code_works/cad_i17_versions/big_1_corners_all.txt";
	ofstream vout(f);*/

	int i1;
	// k for segment number
	for (int k = 0; k < 26; k++) {
		//for (int k = 0; ; k++) {
		//get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, ...
		//vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, ...
		//vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall)
		std::cout << "k-" << k << ": ";
		get_external_partv2(vline1new, vline2new, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
		//find next point with  minimum distance and line repitation within external segment
		find_Point3D_index(ref, vall, vi_all, vi_temp);
		// this is the last point/ slot corner start 
		std::cout << "corner indices: ";
		for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
			v_temp.push_back(vall[vi_temp[ii2]]);
			//shows corner indices
			//std::cout << vi_temp[ii2] << ",";
			//vout << vi_temp[ii2] << ",";
			/*vout << vall[vi_temp[ii2]].x << " ";
			vout << vall[vi_temp[ii2]].y << " ";
			vout << vall[vi_temp[ii2]].z << endl;*/
		}
		std::cout << endl;
		//vout << endl;
		//std::cout << "slot corner vertex repitation:" << vi_temp.size() << endl;
		/*target_plot_update(v_temp);
		target_plot_color_update(0.8, 0.0, 0.0, v_temp.size());*/
		v_temp.resize(0);

		//get the slot border inside


		//look for next external next external
		int next_index = -1;
		float iValue;// = 0.0;
		float newValue = 20000.0;
		for (int i = 0; i < vline1.size(); i++) {
			iValue = dist201(ref, vline1[i]);
			if ((iValue < newValue) & (vline_i[i] != 1)) {
				next_index = i;
				newValue = iValue;
			}
		}
		if (next_index > 0) {
			vselected.push_back(vline1[next_index]);
			//vselected.push_back(vline2[next_index]);
			vline_i[next_index] = vline_i[next_index] + 1;
			vlineordered1.push_back(vline1[next_index]);
			vlineordered2.push_back(vline2[next_index]);
			ref = vline2[next_index];

		}
		else if (next_index < 0) break;
		segment_no++;
	}
	//vout.close();
	std::cout << "segment number: " << segment_no << endl;
	int vt;
	std::cout << "total points: " << vi_all.size() << endl;
	std::cout << "viu2 size: " << viu2.size() << endl;
	std::cout << "v2tl1,2 size: " << vlineordered1.size() << endl;
	/*for (int i = 0; i < vline_i.size(); i++) {
	std::cout << i << ": " << vline_i[i] << endl;
	}*/
	for (int i = 0; i < vlineusedindex.size(); i++) {
		vt = vlineusedindex[i];
		vi_all[vt] = vi_all[vt] + 1;
		//std::cout << i << ": " << vt << endl;
	}

	//find next point
	//get the corner of the slot, which is ref, find all the segments with ref
	std::cout << "vi_used: " << viu2.size() << ", " << vlineordered1.size() << endl;
}
void cad3d::get_segment_numberv502(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points and also corresponding ending
	std::cout << "function: get_segment_numberv502" << endl;
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new;
	convert_cart2pol(vline1, vpol_line1);
	convert_cart2pol(vline2, vpol_line2);
	//vline1new and 2 are all and only outside external segments wqhich is calculated in the next for loop
	for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}
	int segment_no = 0;
	vector<int> viu2;
	//vector<int> vi_all;

	vector<int> vline_i;//index vector for line count?
	vline_i.resize(vline1.size());
	initialize_0i(vline_i);

	Point3D ref;
	vector<int> vi_temp, viendcorner;
	vector<Point3D> v_temp;//used for?

	int iymax2 = find_point3d_max(vline1, 2);
	//finding highest , left most point, for initializing count
	vlineordered1.push_back(vline1[iymax2]);
	vlineordered2.push_back(vline2[iymax2]);
	std::cout << "starting angles: " << vpol_line1[iymax2].y;
	std::cout << " " << vpol_line2[iymax2].y << endl;
	find_Point3D_index(vline1[iymax2], vall, vi_all, vi_temp);// for marking use of points
	std::cout << endl;
	vselected.push_back(vline1[iymax2]);
	viu2.push_back(vlineusedindex[iymax2 * 2]);// why? most probably because line to point conversion ans all
	viu2.push_back(vlineusedindex[iymax2 * 2 + 1]);
	vline_i[iymax2] = vline_i[iymax2] + 1;
	if (vline1.size() == (vlineusedindex.size() / 2)) {
		std::cout << "index length matched" << endl;
	}
	//seg_ui2.push_back(vl1_ui1[iymax2]);

	float dd;
	ref = vline2[iymax2];
	int j, jstart, jend;
	jstart = iymax2;
	string f = "D:/code_works/cad_i17_versions/big_11_doublecorners_all.txt";
	ofstream vout(f);

	int i1;
	// k for segment number
	for (int k = 0; k < 26; k++) {
		//for (int k = 0; ; k++) {
		//get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, ...
		//vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, ...
		//vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall)
		//std::cout << "k-" << k << ": ";
		get_external_partv2(vline1new, vline2new, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
		//find next point with  minimum distance and line repitation within external segment
		find_Point3D_index(ref, vall, vi_all, vi_temp);
		// this is the last point/ slot corner start 
		//std::cout << "corner indices: ";
		for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
			v_temp.push_back(vall[vi_temp[ii2]]);
			//shows corner indices
			//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
			vout << vi_temp[ii2] << ", ";

		}
		//std::cout << endl;
		vout << endl;
		//std::cout << "slot corner vertex repitation:" << vi_temp.size() << endl;
		/*target_plot_update(v_temp);
		target_plot_color_update(0.8, 0.0, 0.0, v_temp.size());*/
		v_temp.resize(0);

		//get the slot border inside


		//look for next external next external
		int next_index = -1;
		float iValue;// = 0.0;
		float newValue = 20000.0;
		for (int i = 0; i < vline1.size(); i++) {
			iValue = dist201(ref, vline1[i]);
			if ((iValue < newValue) & (vline_i[i] != 1)) {
				next_index = i;
				newValue = iValue;
			}
		}
		if (next_index > 0) {
			vselected.push_back(vline1[next_index]);
			//find_Point3D_index(ref, vall, vi_all, vi_temp);
			find_Point3D_indexv2(vline1[next_index], vall, viendcorner);
			vout << "end corner: " << endl;
			for (int ii3 = 0; ii3 < viendcorner.size(); ii3++) {

				//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
				vout << viendcorner[ii3] << ", ";

			}
			vout << endl;
			//vselected.push_back(vline2[next_index]);
			vline_i[next_index] = vline_i[next_index] + 1;
			vlineordered1.push_back(vline1[next_index]);
			vlineordered2.push_back(vline2[next_index]);
			ref = vline2[next_index];

		}
		else if (next_index < 0) break;
		segment_no++;
	}
	vout.close();
	std::cout << "segment number: " << segment_no << endl;
	int vt;
	std::cout << "total points: " << vi_all.size() << endl;
	std::cout << "viu2 size: " << viu2.size() << endl;
	std::cout << "v2tl1,2 size: " << vlineordered1.size() << endl;
	/*for (int i = 0; i < vline_i.size(); i++) {
	std::cout << i << ": " << vline_i[i] << endl;
	}*/
	for (int i = 0; i < vlineusedindex.size(); i++) {
		vt = vlineusedindex[i];
		vi_all[vt] = vi_all[vt] + 1;
		//std::cout << i << ": " << vt << endl;
	}

	//find next point
	//get the corner of the slot, which is ref, find all the segments with ref
	std::cout << "vi_used: " << viu2.size() << ", " << vlineordered1.size() << endl;
}
void cad3d::get_segment_numberv503(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected, vector<Point3D> &vcorner) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points and also corresponding ending
	std::cout << "function: get_segment_numberv503" << endl;
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new;
	convert_cart2pol(vline1, vpol_line1);
	convert_cart2pol(vline2, vpol_line2);
	//vline1new and 2 are all and only outside external segments wqhich is calculated in the next for loop
	for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}
	int segment_no = 0;
	vector<int> viu2;
	//vector<int> vi_all;

	vector<int> vline_i;//index vector for line count?
	vline_i.resize(vline1.size());
	initialize_0i(vline_i);

	Point3D ref;
	vector<int> vi_temp, viendcorner;
	vector<Point3D> v_temp;//used for?

	int iymax2 = find_point3d_max(vline1, 2);
	//finding highest , left most point, for initializing count
	vlineordered1.push_back(vline1[iymax2]);
	vlineordered2.push_back(vline2[iymax2]);
	std::cout << "starting angles: " << vpol_line1[iymax2].y;
	std::cout << " " << vpol_line2[iymax2].y << endl;
	find_Point3D_index(vline1[iymax2], vall, vi_all, vi_temp);// for marking use of points
	std::cout << endl;
	vselected.push_back(vline1[iymax2]);
	viu2.push_back(vlineusedindex[iymax2 * 2]);// why? most probably because line to point conversion ans all
	viu2.push_back(vlineusedindex[iymax2 * 2 + 1]);
	vline_i[iymax2] = vline_i[iymax2] + 1;
	if (vline1.size() == (vlineusedindex.size() / 2)) {
		std::cout << "index length matched" << endl;
	}
	//seg_ui2.push_back(vl1_ui1[iymax2]);

	float dd;
	ref = vline2[iymax2];
	int j, jstart, jend;
	jstart = iymax2;
	//string f = "D:/code_works/cad_i17_versions/big_11_doublecorners_all.txt";
	//ofstream vout(f);

	int i1;
	// k for segment number
	for (int k = 0; k < 26; k++) {
		//for (int k = 0; ; k++) {
		//get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, ...
		//vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, ...
		//vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall)
		std::cout << "k-" << k << ": ";
		get_external_partv2(vline1new, vline2new, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
		//find next point with  minimum distance and line repitation within external segment
		find_Point3D_index(ref, vall, vi_all, vi_temp);
		// this is the last point/ slot corner start 
		std::cout << "corner indices: ";
		for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
			v_temp.push_back(vall[vi_temp[ii2]]);
			vcorner.push_back(vall[vi_temp[ii2]]);
			//shows corner indices
			//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
			//vout << vi_temp[ii2] << ", ";

		}
		std::cout << endl;
		//vout << endl;
		//std::cout << "slot corner vertex repitation:" << vi_temp.size() << endl;
		/*target_plot_update(v_temp);
		target_plot_color_update(0.8, 0.0, 0.0, v_temp.size());*/
		v_temp.resize(0);

		//get the slot border inside


		//look for next external next external
		int next_index = -1;
		float iValue;// = 0.0;
		float newValue = 20000.0;
		for (int i = 0; i < vline1.size(); i++) {
			iValue = dist201(ref, vline1[i]);
			if ((iValue < newValue) & (vline_i[i] != 1)) {
				next_index = i;
				newValue = iValue;
			}
		}
		if (next_index > 0) {
			vselected.push_back(vline1[next_index]);
			//find_Point3D_index(ref, vall, vi_all, vi_temp);
			find_Point3D_indexv2(vline1[next_index], vall, viendcorner);
			//vout << "end corner: " << endl;
			for (int ii3 = 0; ii3 < viendcorner.size(); ii3++) {
				vcorner.push_back(vall[viendcorner[ii3]]);
				//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
				//vout << viendcorner[ii3] << ", ";

			}
			//vout << endl;
			//vselected.push_back(vline2[next_index]);
			vline_i[next_index] = vline_i[next_index] + 1;
			vlineordered1.push_back(vline1[next_index]);
			vlineordered2.push_back(vline2[next_index]);
			ref = vline2[next_index];

		}
		else if (next_index < 0) break;
		segment_no++;
	}
	//vout.close();
	std::cout << "segment number: " << segment_no << endl;
	int vt;
	std::cout << "total points: " << vi_all.size() << endl;
	std::cout << "viu2 size: " << viu2.size() << endl;
	std::cout << "v2tl1,2 size: " << vlineordered1.size() << endl;
	/*for (int i = 0; i < vline_i.size(); i++) {
	std::cout << i << ": " << vline_i[i] << endl;
	}*/
	for (int i = 0; i < vlineusedindex.size(); i++) {
		vt = vlineusedindex[i];
		vi_all[vt] = vi_all[vt] + 1;
		//std::cout << i << ": " << vt << endl;
	}

	//find next point
	//get the corner of the slot, which is ref, find all the segments with ref
	std::cout << "vi_used: " << viu2.size() << ", " << vlineordered1.size() << endl;
}
void cad3d::get_segment_numberv504(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected, vector<Point3D> &vcorner) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points and also corresponding ending
	std::cout << "function: get_segment_numberv504" << endl;
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new;
	convert_cart2pol(vline1, vpol_line1);
	convert_cart2pol(vline2, vpol_line2);
	//vline1new and 2 are all and only outside external segments wqhich is calculated in the next for loop
	//ordering the points angle wise (very important)
	for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}// getting only exernal only segments
	int segment_no = 0;
	vector<int> viu2;
	//vector<int> vi_all;

	vector<int> vline_i;//index vector for line count?
	vline_i.resize(vline1.size());
	initialize_0i(vline_i);

	Point3D ref;
	vector<int> vi_temp, viendcorner;
	vector<Point3D> v_temp;//used for?

	int iymax2 = find_point3d_min(vpol_line1, 2);// minding point with minimum angle
	//finding highest , left most point, for initializing count
	vlineordered1.push_back(vline1[iymax2]);
	vlineordered2.push_back(vline2[iymax2]);
	std::cout << "starting angles: " << vpol_line1[iymax2].y;
	std::cout << ", " << vpol_line2[iymax2].y << endl;
	find_Point3D_index(vline1[iymax2], vall, vi_all, vi_temp);// for marking use of points
	std::cout << endl;
	vselected.push_back(vline1[iymax2]);
	viu2.push_back(vlineusedindex[iymax2 * 2]);// why? most probably because line to point conversion ans all
	viu2.push_back(vlineusedindex[iymax2 * 2 + 1]);
	vline_i[iymax2] = vline_i[iymax2] + 1;
	if (vline1.size() == (vlineusedindex.size() / 2)) {
		std::cout << "index length matched" << endl;
	}
	//seg_ui2.push_back(vl1_ui1[iymax2]);

	float dd;
	ref = vline2[iymax2];
	int j, jstart, jend;
	jstart = iymax2;
	get_external_partv2(vline1new, vline2new, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
	std::cout << "corner indices: ";
	for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
		v_temp.push_back(vall[vi_temp[ii2]]);
		vcorner.push_back(vall[vi_temp[ii2]]);
	}
	std::cout << endl;
	int next_index = -1;
	float iValue;// = 0.0;
	float newValue = 20000.0;
	for (int i = 0; i < vline1.size(); i++) {
		iValue = dist201(ref, vline1[i]);
		if ((iValue < newValue) & (vline_i[i] != 1)) {
			next_index = i;
			newValue = iValue;
		}
	}
	if (next_index > 0) {
		vselected.push_back(vline1[next_index]);
		//find_Point3D_index(ref, vall, vi_all, vi_temp);
		find_Point3D_indexv2(vline1[next_index], vall, viendcorner);
		//vout << "end corner: " << endl;
		for (int ii3 = 0; ii3 < viendcorner.size(); ii3++) {
			vcorner.push_back(vall[viendcorner[ii3]]);
			//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
			//vout << viendcorner[ii3] << ", ";

		}
		//vout << endl;
		//vselected.push_back(vline2[next_index]);
		vline_i[next_index] = vline_i[next_index] + 1;
		vlineordered1.push_back(vline1[next_index]);
		vlineordered2.push_back(vline2[next_index]);
		ref = vline2[next_index];

	}

	v_temp.resize(0);
	//string f = "D:/code_works/cad_i17_versions/big_11_doublecorners_all.txt";
	//ofstream vout(f);

	int i1;
	// k for segment number
	/*for (int k = 0; k < 26; k++) {
		//for (int k = 0; ; k++) {
		//get_external_partv2(vector<Point3D> &vline1, vector<Point3D> &vline2, ...
		//vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, ...
		//vector<int> &vlineusedindex, vector<int> &viu2, Point3D &ref, vector<int> &vline_i, vector<Point3D> &vselected, vector<int> &vi_used, vector<Point3D> &vall)
		std::cout << "k-" << k << ": ";
		get_external_partv2(vline1new, vline2new, vlineordered1, vlineordered2, vlineusedindex, viu2, ref, vline_i, vselected, vi_all, vall);
		//find next point with  minimum distance and line repitation within external segment
		find_Point3D_index(ref, vall, vi_all, vi_temp);
		// this is the last point/ slot corner start
		std::cout << "corner indices: ";
		for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
			v_temp.push_back(vall[vi_temp[ii2]]);
			vcorner.push_back(vall[vi_temp[ii2]]);
			//shows corner indices
			//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
			//vout << vi_temp[ii2] << ", ";

		}
		std::cout << endl;
		//vout << endl;
		//std::cout << "slot corner vertex repitation:" << vi_temp.size() << endl;

		v_temp.resize(0);

		//get the slot border inside


		//look for next external next external
		int next_index = -1;
		float iValue;// = 0.0;
		float newValue = 20000.0;
		for (int i = 0; i < vline1.size(); i++) {
			iValue = dist201(ref, vline1[i]);
			if ((iValue<newValue) & (vline_i[i] != 1)) {
				next_index = i;
				newValue = iValue;
			}
		}
		if (next_index > 0) {
			vselected.push_back(vline1[next_index]);
			//find_Point3D_index(ref, vall, vi_all, vi_temp);
			find_Point3D_indexv2(vline1[next_index], vall, viendcorner);
			//vout << "end corner: " << endl;
			for (int ii3 = 0; ii3 < viendcorner.size(); ii3++) {
				vcorner.push_back(vall[viendcorner[ii3]]);
				//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
				//vout << viendcorner[ii3] << ", ";

			}
			//vout << endl;
			//vselected.push_back(vline2[next_index]);
			vline_i[next_index] = vline_i[next_index] + 1;
			vlineordered1.push_back(vline1[next_index]);
			vlineordered2.push_back(vline2[next_index]);
			ref = vline2[next_index];

		}
		else if (next_index < 0) break;
		segment_no++;
	}*/
	//vout.close();
	std::cout << "segment number: " << segment_no << endl;
	int vt;
	std::cout << "total points: " << vi_all.size() << endl;
	std::cout << "viu2 size: " << viu2.size() << endl;
	std::cout << "v2tl1,2 size: " << vlineordered1.size() << endl;

	for (int i = 0; i < vlineusedindex.size(); i++) {
		vt = vlineusedindex[i];
		vi_all[vt] = vi_all[vt] + 1;
		//std::cout << i << ": " << vt << endl;
	}

	//find next point
	//get the corner of the slot, which is ref, find all the segments with ref
	std::cout << "vi_used: " << viu2.size() << ", " << vlineordered1.size() << endl;
}
void cad3d::get_segment_numberv505(vector<Point3D> &vline1, vector<Point3D> &vline2, vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2, vector<int> &vlineusedindex, vector<int> &vi_all, vector<Point3D> &vall, vector<Point3D> &vselected, vector<Point3D> &vcorner) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points and also corresponding ending
	//std::cout << "function: get_segment_numberv505 with starting corners" << endl;
	std::cout << "Open slot segment number calculation for circular boundary" << endl;
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new;
	convert_cart2pol(vline1, vpol_line1);
	convert_cart2pol(vline2, vpol_line2);
	//vline1new and 2 are oriented in counter clock wise direction
	for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			vlineordered1.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			vlineordered1.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}
	//vlineordered1 = vline1new;
	//vlineordered2 = vline2new;
	// for finding v start corner
	int zero_count = 0;
	float dd1, dd2;
	int k = 0;// k  = segment number
	Point3D ref;
	vector<int> vi_temp, viendcorner;
	for (int i = 0; i < vline2new.size(); i++) {
		//dd = dist201(vline2new[i], vline1new[5]);
		for (int j = 0; j < vline1new.size(); j++) {
			dd1 = dist201(vline2new[i], vline1new[j]);
			if (dd1 == 0)
				zero_count = zero_count + 1;
			//below one for removing the starting point
			if (i != j) {
				dd2 = dist201(vline2new[i], vline2new[j]);
				if (dd2 == 0)
					zero_count = zero_count + 1;
			}
		}
		if (!zero_count) {
			vcorner.push_back(vline2new[i]);
			//std::cout << "k-" << k << ": ";
			ref = vline2new[i];
			find_Point3D_index(ref, vall, vi_all, vi_temp);
			// this is the last point/ slot corner start 
			//std::cout << "start corner indices: ";
			for (int ii2 = 0; ii2 < vi_temp.size(); ii2++) {
				//v_temp.push_back(vall[vi_temp[ii2]]);
				vcorner.push_back(vall[vi_temp[ii2]]);
				//shows corner indices
				//vout << vall[vi_temp[ii2]].x << " " << vall[vi_temp[ii2]].y << " "<< vall[vi_temp[ii2]].z << endl;
				//std::cout << vi_temp[ii2] << ", ";

			}
			//std::cout << endl;
			k++;
		}

		zero_count = 0;

	}
	std::cout << "Open slot number: " <<k << endl;
}
void cad3d::get_end_corner_v101(vector<Point3D> &vlineordered1, vector<Point3D> &vlineordered2,  vector<Point3D> &vendcorner, vector<Point3D> &vcorner) {
	//vline1 and 2 are input as external segments
	//finds out the starting corner points and also corresponding ending
	std::cout << "get_end_corner_v101" << endl;
	vector<Point3D> vpol_line1, vpol_line2, vline1new, vline2new,vcornerpol;
	convert_cart2pol(vlineordered1, vpol_line1);
	convert_cart2pol(vlineordered1, vpol_line2);
	convert_cart2pol(vcorner, vcornerpol);
	//vline1new and 2 are oriented in counter clock wise direction
	/*for (int i = 0; i < vline1.size(); i++) {
		if (vpol_line1[i].y >= vpol_line2[i].y) {
			vline1new.push_back(vline1[i]);
			//vlineordered1.push_back(vline1[i]);
			vline2new.push_back(vline2[i]);
		}
		else {
			vline1new.push_back(vline2[i]);
			//vlineordered1.push_back(vline2[i]);
			vline2new.push_back(vline1[i]);
		}

	}*/
	
	// for finding v end corner
	int zero_count = 0;
	float dd1, dd2;
	int k = 0;// k  = segment number
	Point3D ref;
	int i_temp;
	vector<int> vi_temp, viendcorner;
	for (int i = 0; i < vcorner.size(); i++) {
		//dd = dist201(vline2new[i], vline1new[5]);
		dd1 = dist201(vcorner[i], vlineordered1[0]);
		i_temp = 0;//index of vlineordered1
		for (int j = 1; j < vlineordered1.size(); j++) {
			dd2 = dist201(vcorner[i], vlineordered1[j]);
			if ((dd2 < dd1) & (dd2 != 0) & (vcornerpol[i].y > vpol_line1[j].y)) {
				dd1 = dd2;
				i_temp = j;
			}
			
		}
		vendcorner.push_back(vlineordered1[i_temp]);

	}
	
}

int cad3d::get_trianglesegment_wo_rv3(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	std::cout << "points and points index size: ";
	//std::cout << vall.size() << "," << vi_used.size() << endl;
	//show_index_frequency(vi_used);
	float d1, d2, d3, dmin;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2;
	convert_cart2pol(vall, vallpol);
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);
			if (ivt1.size() > 0) {
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " << vi_used[ivt1[j]] << endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					d3 = dist201(vall[iref], vinit2[k]);
					if (d3 < dmin) {
						dmin = d3;
						iref_out = inext[k];
					}
				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {
					//std::cout << "counter clock:" << ivt2[j];
					//std::cout << ": " << vi_used[ivt2[j]] << endl;
					//std::cout << "d2: " << ivt2[j] << ", ";
					//std::cout << vi_used[ivt2[j]] << ", ";
					//std::cout << vallpol[ivt2[j]].x << ", ";
					//std::cout << vallpol[ivt2[j]].y << endl;
					if (vi_used[ivt2[j]] < 1) {
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					d3 = dist201(vall[iref], vinit2[k]);
					if (d3 < dmin) {
						dmin = d3;
						iref_out = inext[k];
					}
				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
		}
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);
	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	return iref_out;
}
int cad3d::get_trianglesegment_wo_rv301(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	std::cout << "points and points index size: ";
	//std::cout << vall.size() << "," << vi_used.size() << endl;
	//show_index_frequency(vi_used);
	float d1, d2, d3, dmin;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2;
	convert_cart2pol(vall, vallpol);
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					d3 = dist201(vall[iref], vinit2[k]);
					if (d3 < dmin) {
						dmin = d3;
						iref_out = inext[k];
					}
				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {
					//dtheta = del_theta(v_line2[i], v_line1[i]);
					//std::cout << "counter clock:" << ivt2[j];
					//std::cout << ": " << vi_used[ivt2[j]] << endl;
					//std::cout << "d2: " << ivt2[j] << ", ";
					//std::cout << vi_used[ivt2[j]] << ", ";
					//std::cout << vallpol[ivt2[j]].x << ", ";
					//std::cout << vallpol[ivt2[j]].y << endl;
					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
							//error
							//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					d3 = dist201(vall[iref], vinit2[k]);
					if (d3 < dmin) {
						dmin = d3;
						iref_out = inext[k];
					}
				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
		}
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);
	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	return iref_out;
}
int cad3d::get_trianglesegment_wo_rv302(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	std::cout << "points and points index size: ";

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
	//std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta;
	int idump;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					//std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// ";
					//if (vinit2pol[k].x >=275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)) {
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];
						/*if (vinit2pol[k].x < 275.5)
						std::cout <<" error " ;*/
						//continue;
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
						//if (r_temp > rmin ){
							//- 0.005) || r_temp >(rmin + 0.005)) {
							dmin = d3;
							iref_out = inext[k];
						//}
							//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}

		}
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);
	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
cad3d::circle2d cad3d::compute_circle_4_point3d(Point3D p1, Point3D p2, Point3D p3) {
	circle2d Crcl;
	float ax = p1.x;
	float ay = p1.y;
	float bx = p2.x;
	float by = p2.y;
	float cx = p3.x;
	float cy = p3.y;

	float x1 = (bx + ax) / 2;
	float y1 = (by + ay) / 2;
	float dy1 = bx - ax;
	float dx1 = -(by - ay);

	float x2 = (bx + cx) / 2;
	float y2 = (by + cy) / 2;
	float dy2 = cx - bx;
	float dx2 = -(cy - by);

	float ox = (y1 * dx1 * dx2 + x2 * dx1 * dy2 - x1 * dy1 * dx2 - y2 * dx1 * dx2) / (dx1 * dy2 - dy1 * dx2);
	float oy = (ox - x1) * dy1 / dx1 + y1;

	float dx = ox - ax;
	float dy = oy - ay;
	float radius = (int)sqrt(dx * dx + dy * dy);

	Crcl.x = ox;
	Crcl.y = oy;
	Crcl.r = radius;

	/*return (p1.x - p2.x)*(p1.x - p2.x) +
	(p1.y - p2.y)*(p1.y - p2.y);*/
	return Crcl;
	//http://forums.codeguru.com/showthread.php?442320-draw-circle-from-3-points
}
/*void cad3d::mark_inner_radius(vector<int> &vi_used, vector<Point3D> &vall) {

}*/
int cad3d::get_trianglesegment_wo_rv303(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge
	//closed edge
	// outward edge

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					//std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// ";
					//if (vinit2pol[k].x >=275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)) {
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];
						/*if (vinit2pol[k].x < 275.5)
						std::cout <<" error " ;*/
						//continue;
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
					//if (r_temp > rmin ){
					//- 0.005) || r_temp >(rmin + 0.005)) {
					dmin = d3;
					iref_out = inext[k];
					//}
					//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}

		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
int cad3d::get_trianglesegment_wo_rv3031(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	//convert_cart2polPoint3D();
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);
			//if the point has been used and what are the connected points
			
			//if there are connected points, pass them to vinit1
			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}
			
			//if only one point, trouble reduces
			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}

			////if more than one point, follow the sub routine
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vallpol[iref], vinit2pol[k]);
					std::cout << dtheta << "//";
					//if ((d3 < dmin) & (vinit2pol[k].x < vallpol[iref].x) & (abs(dtheta) < 5)) {
					if ((d3 < dmin) & (vinit2pol[k].x < vallpol[iref].x) ) {
						dmin = d3;
						iref_out = inext[k];
						/*if (vinit2pol[k].x < 275.5)
						std::cout <<" error " ;*/
						//continue;
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << vinit2pol[iref_out].y << endl;
			}
			//std::cout << iref_out << ": " << vinit2pol[inext[iref_out]].y << endl;
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vallpol[iref], vinit2pol[k]);
					r_temp = vinit2pol[k].x;
				
					std::cout << dtheta << "//";
					//if ((d3 < dmin) & (vinit2pol[k].x < vallpol[iref].x)& (abs(dtheta) < 5)) {
					if ((d3 < dmin) & (vinit2pol[k].x < vallpol[iref].x)) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						//if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}
			//std::cout << iref_out << ": " << vinit2pol[inext[iref_out]].y << endl;
		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);
		//dtheta = del_theta(vt1, vt2);
		//std::cout << dtheta << endl;

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
int cad3d::get_1st_inward_point_v101(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge
	

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	//convert_cart2polPoint3D();
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);
			//if the point has been used and what are the connected points

			//if there are connected points, pass them to vinit1
			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			//if only one point, trouble reduces
			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
				std::cout <<"single"<<endl;
			}

			////if more than one point, follow the sub routine
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);

					if ((d3 < dmin) &  (d3 != 0)) {
						if (vinit2pol[k].x < vallpol[iref].x)  {
							if (vinit2pol[k].y <= vallpol[iref].y) {
								//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
								if (vinit2pol[k].x < 275.9) continue;
								dmin = d3;
								iref_out = inext[k];
								/*if (vinit2pol[k].x < 275.5)
								std::cout <<" error " ;*/
								//continue;
								idump = 0;
							}
						}
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
				std::cout << "single" << endl;
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
					//if (r_temp > rmin ){
					//- 0.005) || r_temp >(rmin + 0.005)) {
					dmin = d3;
					iref_out = inext[k];
					//}
					//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					//if (d3 < dmin) {
					/*if ((d3 < dmin) & (vinit2pol[k].x < vallpol[iref].x) & (d3 != 0)) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}*/
					//}
					if ((d3 < dmin) &  (d3 != 0)) {
						if (vinit2pol[k].x < vallpol[iref].x) {
							if (vinit2pol[k].y <= vallpol[iref].y) {
								//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
								if (vinit2pol[k].x < 275.9) continue;
								dmin = d3;
								iref_out = inext[k];
								/*if (vinit2pol[k].x < 275.5)
								std::cout <<" error " ;*/
								//continue;
								idump = 0;
							}
						}
					}
				}

			}

		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
int cad3d::get_1st_inward_point_v102(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge
	std::cout << "get_1st_inward_point_v102" << endl;

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	//convert_cart2polPoint3D();
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump,imin;
	//std::cout << "i		used st		r		theta" << endl;
	dmin = dist201(vall[iref], vall[0]);
	imin = 0;
	//d1 = dist201(vall[iref], vall[i]);
	for (int i = 1; i < vall.size(); i++) {
		d1 = dist201(vall[iref], vall[i]);
		//for anti clockwise  calculation
		if ((d1 < dmin) & (d1 != 0)) {
			//if ((vallpol[i].x < vallpol[iref].x) & (vi_used[i] != 0)) {
				//if (vallpol[i].y < vallpol[iref].y) {
			if (vallpol[i].x < (vallpol[iref].x-1)) {
					//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {

					dmin = d1;
					imin = i;
					/*if (vinit2pol[k].x < 275.5)
					std::cout <<" error " ;*/
					//continue;
					//idump = 0;
				}

			//}
			/*find_Point3D_indexv2(v_line2[i], vall, ivt1);
			//if the point has been used and what are the connected points

			//if there are connected points, pass them to vinit1
			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			//if only one point, trouble reduces
			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
				std::cout << "single" << endl;
			}

			////if more than one point, follow the sub routine
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);

					if ((d3 < dmin) &  (d3 != 0)) {
						if (vinit2pol[k].x < vallpol[iref].x) {
							if (vinit2pol[k].y <= vallpol[iref].y) {
								//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
								if (vinit2pol[k].x < 275.9) continue;
								dmin = d3;
								iref_out = inext[k];
								/*if (vinit2pol[k].x < 275.5)
								std::cout <<" error " ;*/
		}

	}
	iref_out = imin;
	
	return iref_out;
}

int cad3d::get_slot_part(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge
	//closed edge
	// outward edge
	//std::cout << "get_slot_part" << endl;
	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;
	//circle2d c_temp;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					//std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// ";
					//if (vinit2pol[k].x >=275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)) {
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];
						/*if (vinit2pol[k].x < 275.5)
						std::cout <<" error " ;*/
						//continue;
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
					//if (r_temp > rmin ){
					//- 0.005) || r_temp >(rmin + 0.005)) {
					dmin = d3;
					iref_out = inext[k];
					//}
					//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}

		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		/*if (i >= 2) {
			c_temp = compute_circle_4_point3d(vt_line1[i-2], vt_line1[i - 1], vt_line1[i]);
			std::cout << "r: "<<c_temp.r << endl;
		}*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
int cad3d::get_slot_part_v101(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge

	//std::cout << "function: get_slot_part v101" << endl;
	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;

	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {

					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						std::cout << "d1: " << " (cw angle difference: " << dtheta << ")";
						std::cout << " selected" << endl;
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// " << endl;
					//if (vinit2pol[k].x >=275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)) {
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;*/
			std::cout << iref_out << endl;
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						std::cout << "d2: " << vi_used[ivt2[j]] << " selected" << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// " << endl;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
					//if (r_temp > rmin ){
					//- 0.005) || r_temp >(rmin + 0.005)) {
					dmin = d3;
					iref_out = inext[k];
					//}
					//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}
			std::cout << iref_out << ": " << endl;
		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		/*if (i >= 2) {
		c_temp = compute_circle_4_point3d(vt_line1[i-2], vt_line1[i - 1], vt_line1[i]);
		std::cout << "r: "<<c_temp.r << endl;
		}*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}
int cad3d::get_trianglesegment_wo_rv304(vector<int> &vi_used, vector<Point3D> &vall, vector<Point3D> &v_line1, vector<Point3D> &v_line2, int iref) {
	// inward edge
	//closed edge
	// outward edge

	float d1, d2, d3, dmin, dnew;
	vector<Point3D>  vt_line1, vt_line2, vallpol;
	vector<Point3D> vinit1, vinit2, vinit1pol, vinit2pol;
	convert_cart2pol(vall, vallpol);
	float rmin = find_point3d_min(vallpol, 1); //since this is the index
											   //std::cout << "minimum radius: " <<rmin<< endl;
	float r_temp;
	Point3D vt1, vt2;
	vector<int> ivt1, ivt2, ivdump, inext;
	int iref_out = -1;
	float dtheta, dr, rad1;
	int idump;
	//std::cout << "i		used st		r		theta" << endl;
	for (int i = 0; i < v_line1.size(); i++) {
		d1 = dist201(vall[iref], v_line1[i]);
		//for anti clockwise  calculation
		if (d1 == 0) {
			find_Point3D_indexv2(v_line2[i], vall, ivt1);

			if (ivt1.size() > 0) {
				dtheta = del_theta(vall[iref], v_line2[i]);
				for (int j = 0; j < ivt1.size(); j++) {
					//std::cout << "clock:" << ivt1[j] << endl;
					//std::cout << "d1: " << ivt1[j] << ", ";

					//std::cout << " (cw angle difference: " << dtheta << ")";
					if (vi_used[ivt1[j]] < 1) {
						vt_line1.push_back(v_line1[i]);
						vt_line2.push_back(v_line2[i]);
						vinit1.push_back(v_line1[i]);
						vinit2.push_back(v_line2[i]);
						find_Point3D_index(v_line2[i], vall, vi_used, ivdump);
						//iref_out = ivt1[j];
						inext.push_back(ivt1[j]);
						//std::cout << "d1: " <<  " (cw angle difference: " << dtheta << ")"<<endl;
						//std::cout << " selected";
					}
				}
			}

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				dmin = dist201(vall[iref], vinit2[0]);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {
					r_temp = vinit2pol[k].x;
					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					//std::cout << inext[k] << ": distance: " << d3 << "// " << "angle: " << dtheta << "// ";
					//if (vinit2pol[k].x >=275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)) {
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];
						/*if (vinit2pol[k].x < 275.5)
						std::cout <<" error " ;*/
						//continue;
						idump = 0;
					}
					//}

				}
				//std::cout << iref_out << ": " << d3 << "// ";
			}
			/*;
			std::cout << iref_out << ": " << dnew << endl;*/
			/*dnew = dist201(vall[iref], vall[iref_out]);
			dtheta = del_theta(vall[iref], vall[iref_out]);
			std::cout << iref_out << ": " << dnew << "// " << "angle: " << dtheta << endl;*/
		}
		d2 = dist201(vall[iref], v_line2[i]);
		if (d2 == 0) {
			//for clock wise calculation
			find_Point3D_indexv2(v_line1[i], vall, ivt2);
			//dtheta = del_theta(vall[iref], v_line1[i]);
			//std::cout << " (ccw angle difference: " << dtheta << ")";
			if (ivt2.size() > 0) {
				for (int j = 0; j < ivt2.size(); j++) {

					dtheta = del_theta(vall[iref], vall[ivt2[j]]);
					if (vi_used[ivt2[j]] < 1) {
						//if ((vi_used[ivt2[j]] < 1)&(dtheta<70)) {
						//error
						//std::cout << " (ccw angle difference: " << dtheta << ")";
						vt_line1.push_back(v_line2[i]);
						vt_line2.push_back(v_line1[i]);
						vinit1.push_back(v_line2[i]);
						vinit2.push_back(v_line1[i]);
						find_Point3D_index(v_line1[i], vall, vi_used, ivdump);
						iref_out = ivt2[j];
						inext.push_back(ivt2[j]);
						//std::cout << "d2: " << " (ccw angle difference: " << dtheta << ")" << endl;
						//std::cout << "d2: " << vi_used[ivt2[j]] << endl;
						//std::cout << " selected";
					}
				}
			}
			//std::cout << v_line1[i].x << ", " << v_line2[i].x << endl;

			if (vinit1.size() == 1) {
				vt1 = vinit1[0];
				vt2 = vinit2[0];
				iref_out = inext[0];
			}
			if (vinit1.size() > 1) {
				dmin = dist201(vall[iref], vinit2[0]);
				convert_cart2pol(vinit1, vinit1pol);
				convert_cart2pol(vinit2, vinit2pol);
				iref_out = inext[0];
				for (int k = 1; k < vinit1.size(); k++) {

					d3 = dist201(vall[iref], vinit2[k]);
					dtheta = del_theta(vall[iref], vinit2[k]);
					r_temp = vinit2pol[k].x;
					/*if (d3 < dmin) {
					//if ((d3 < dmin) & (r_temp > (rmin+0.5))){
					//if (r_temp > rmin ){
					//- 0.005) || r_temp >(rmin + 0.005)) {
					dmin = d3;
					iref_out = inext[k];
					//}
					//continue;

					}*/
					//if (vinit2pol[k].x >= 275.8) {
					//if ((vinit2pol[k].x >= 275.8) & (dtheta < 100)){
					if (d3 < dmin) {
						//if ((d3 < dmin) & (r_temp > (rmin + 0.5))) {
						if (vinit2pol[k].x < 275.9) continue;
						dmin = d3;
						iref_out = inext[k];

						idump = 0;
					}
					//}
				}

			}

		}
		/*dr = vallpol[iref].x - vallpol[iref_out].x;
		std::cout << " del r: " << dr;*/
		vt_line1.push_back(vt1);
		vt_line2.push_back(vt2);

	}
	/*target_plot2line_update(vt_line1, vt_line2);
	target_plot2_color_update(0.4, 0.4, 0.8, vt_line1.size());*/
	//dnew = dist201(vall[iref], vall[iref_out]);
	//std::cout << iref_out << ": " << dnew << endl;
	return iref_out;
}

void cad3d::get_internal_slot_v101(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1) {
	slot_line1.push_back(centered_block[vi_temp1]);

	int vi_temp2 = get_trianglesegment_wo_rv301(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	slot_line2.push_back(centered_block[vi_temp2]);
	for (int i = 0; i < 200; i++) {
		slot_line1.push_back(centered_block[vi_temp2]);
		vi_temp2 = get_trianglesegment_wo_rv301(vi_all, centered_block, v_line1, v_line2, vi_temp2);
		if (vi_temp2 < 0) {
			slot_line1.pop_back();
			break;
		}
		slot_line2.push_back(centered_block[vi_temp2]);
		cout << i << ": " << vi_temp2 << endl << endl;
		v1.push_back(centered_block[vi_temp2]);
	}
}
void cad3d::get_internal_slot_v102(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1) {
	slot_line1.push_back(centered_block[vi_temp1]);
	vector<Point3D> vallpol;
	convert_cart2pol(centered_block, vallpol);
	int vi_temp2 = get_trianglesegment_wo_rv301(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	slot_line2.push_back(centered_block[vi_temp2]);
	float dtheta;
	//for (int i = 0; i < 200; i++) {
	for (int i = 0; i < 60; i++) {
		slot_line1.push_back(centered_block[vi_temp2]);
		vi_temp2 = get_trianglesegment_wo_rv302(vi_all, centered_block, v_line1, v_line2, vi_temp2);
		if (vi_temp2 < 0) {
			slot_line1.pop_back();
			break;
		}
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp2]);
		//cout << i << ": " << vi_temp2 << ", del theta: " << dtheta << " radius:  "<< vallpol[vi_temp2].x <<endl;
		std::cout << i << ": del theta : " << dtheta << " radius : " << vallpol[vi_temp2].x << endl;
		/*if (dtheta > 130) {
			slot_line1.pop_back();
			continue;
		}*/
		slot_line2.push_back(centered_block[vi_temp2]);

		v1.push_back(centered_block[vi_temp2]);
	}
}
void cad3d::get_internal_slot_v103(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1) {
	slot_line1.push_back(centered_block[vi_temp1]);
	vector<Point3D> vallpol;
	convert_cart2pol(centered_block, vallpol);
	int vi_temp2 = get_trianglesegment_wo_rv303(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	slot_line2.push_back(centered_block[vi_temp2]);
	float dtheta;
	//for (int i = 0; i < 200; i++) {
	for (int i = 0; i < 60; i++) {
		slot_line1.push_back(centered_block[vi_temp2]);
		vi_temp2 = get_trianglesegment_wo_rv303(vi_all, centered_block, v_line1, v_line2, vi_temp2);
		if (vi_temp2 < 0) {
			slot_line1.pop_back();
			break;
		}
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp2]);
		//cout << i << ": " << vi_temp2 << ", del theta: " << dtheta << " radius:  "<< vallpol[vi_temp2].x <<endl;
		std::cout << i << ": del theta : " << dtheta << " radius : " << vallpol[vi_temp2].x << endl;
		/*if (dtheta > 130) {
		slot_line1.pop_back();
		continue;
		}*/
		slot_line2.push_back(centered_block[vi_temp2]);

		v1.push_back(centered_block[vi_temp2]);
	}
}
void cad3d::get_internal_slot_v1031(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1) {
	slot_line1.push_back(centered_block[vi_temp1]);
	v1.push_back(centered_block[vi_temp1]);
	vector<Point3D> vallpol;
	convert_cart2pol(centered_block, vallpol);
	//int vi_temp2 = get_trianglesegment_wo_rv3031(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	int vi_temp2 = get_1st_inward_point_v102(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	slot_line2.push_back(centered_block[vi_temp2]);
	float dtheta;
	v1.push_back(centered_block[vi_temp2]);
	
	slot_line1.push_back(centered_block[vi_temp2]);
	vi_temp2 = get_trianglesegment_wo_rv301(vi_all, centered_block, v_line1, v_line2, vi_temp2); (vi_all, centered_block, v_line1, v_line2, vi_temp2);
	v1.push_back(centered_block[vi_temp2]);
	slot_line2.push_back(centered_block[vi_temp2]);

	/*slot_line1.push_back(centered_block[vi_temp2]);
	vi_temp2 = get_trianglesegment_wo_rv3031(vi_all, centered_block, v_line1, v_line2, vi_temp2); (vi_all, centered_block, v_line1, v_line2, vi_temp2);
	v1.push_back(centered_block[vi_temp2]);
	slot_line2.push_back(centered_block[vi_temp2]);*/
	//vi_temp2 = get_1st_inward_point_v102(vi_all, centered_block, v_line1, v_line2, vi_temp1);
	//for (int i = 0; i < 200; i++) {
	/*for (int i = 1; i < 20; i++) {
		slot_line1.push_back(centered_block[vi_temp2]);
		vi_temp2 = get_trianglesegment_wo_rv3031(vi_all, centered_block, v_line1, v_line2, vi_temp2);
		if (vi_temp2 < 0) {
			slot_line1.pop_back();
			break;
		}
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp2]);
		//cout << i << ": " << vi_temp2 << ", del theta: " << dtheta << " radius:  "<< vallpol[vi_temp2].x <<endl;
		//std::cout << i << ": del theta : " << dtheta << " radius : " << vallpol[vi_temp2].x << endl;
		
		slot_line2.push_back(centered_block[vi_temp2]);

		v1.push_back(centered_block[vi_temp2]);
		dtheta = 0;
	}*/
}

void cad3d::relate_corner_w_line(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &vcorner, vector<Point3D> &vcorner_new) {
	std::cout << vcorner.size() <<endl;
	int itemp1 = -1;
	int itemp2 = -1;
	float d1, d2,d;
	int n = 0;
	//vector<Point3D> vcorner_new;
	//reduce vconrer to segment nubmber
	vcorner_new.push_back(vcorner[0]);
	for (int i = 1; i < vcorner.size(); i++) {
		d = dist201(vcorner_new[n], vcorner[i]);
		if (d != 0) {
			vcorner_new.push_back(vcorner[i]);
			n++;
		}
	}
	
	std::cout << vcorner.size() << ", "<< vcorner_new.size()<<endl;
	
	/*for (int i = 0; i < vcorner_new.size(); i++) {
		std::cout << i << ": ";
		for (int j = 0; j < v_line1.size(); j++) {
			d1 = dist201(vcorner_new[i], v_line1[j]);
			d2 = dist201(vcorner_new[i], v_line2[j]);
			if (d1 == 0) {
				std::cout << "d1: " << j << endl;
			}
			if (d2 == 0) {
				std::cout << "d2: " << j << endl;
			}
		}
	}*/

}

void cad3d::get_internal_slot_v104(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1s1, int vi_temp1s2) {
	slot_line1.push_back(centered_block[vi_temp1s1]);
	vector<Point3D> vallpol;
	Point3D temp_slotline1pol;
	std::cout << "function: get_internal_slot_v104" << endl;
	convert_cart2pol(centered_block, vallpol);
	int rmin = find_point3d_min(vallpol, 1); //

	//first side or edge 1
	int vi_temp21 = get_trianglesegment_wo_rv304(vi_all, centered_block, v_line1, v_line2, vi_temp1s1);
	slot_line2.push_back(centered_block[vi_temp21]);
	float dtheta, drr;
	circle2d c_temp;
	std::cout << "e1:: r,   del r,    del theta,   del_r*del_theta  " << endl;
	//for (int i = 0; i < 200; i++) {
	//one side
	for (int i = 0; i < 40; i++) {
		slot_line1.push_back(centered_block[vi_temp21]);
		convert_cart2polPoint3D(centered_block[vi_temp21], temp_slotline1pol);
		vi_temp21 = get_slot_part(vi_all, centered_block, v_line1, v_line2, vi_temp21);
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp21]);
		drr = temp_slotline1pol.x - vallpol[vi_temp21].x;
		std::cout << i << ": " << vallpol[vi_temp21].x << ", " << drr << ", ";
		std::cout << dtheta << ", " << dtheta * drr << endl;
		if (vi_temp21 < 0) {
			slot_line1.pop_back();
			break;
		}
		/*else if (vallpol[vi_temp21].x < (vallpol[rmin].x + 0.5)) {
			slot_line1.pop_back();
			break;
		}*/
		/*else if (drr < 0) {
			slot_line1.pop_back();
			//slot_line1.pop_back();
			//slot_line2.pop_back();
			//v1.pop_back();
			//break;
		}*/
		slot_line2.push_back(centered_block[vi_temp21]);
		v1.push_back(centered_block[vi_temp21]);
	}

	//other side or edge 2

	int vi_temp22 = get_trianglesegment_wo_rv304(vi_all, centered_block, v_line1, v_line2, vi_temp1s2);
	slot_line2.push_back(centered_block[vi_temp22]);
	//float dtheta, drr;
	//circle2d c_temp;
	std::cout << "e2:: r,   del r,    del theta,   del_r*del_theta  " << endl;
	for (int i = 0; i < 40; i++) {
		slot_line1.push_back(centered_block[vi_temp22]);
		convert_cart2polPoint3D(centered_block[vi_temp22], temp_slotline1pol);
		vi_temp22 = get_slot_part(vi_all, centered_block, v_line1, v_line2, vi_temp22);
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp22]);
		drr = temp_slotline1pol.x - vallpol[vi_temp22].x;
		std::cout << i << ": " << vallpol[vi_temp22].x << ", " << drr << ", ";
		std::cout << dtheta << ", " << dtheta * drr << endl;
		if (vi_temp22 < 0) {
			slot_line1.pop_back();
			break;
		}
		else if (vallpol[vi_temp22].x < (vallpol[rmin].x + 0.5)) {
			slot_line1.pop_back();
			break;
		}
		else if (drr < 0) {
			slot_line1.pop_back();
			slot_line1.pop_back();
			slot_line2.pop_back();
			v1.pop_back();
			break;
		}
		slot_line2.push_back(centered_block[vi_temp22]);
		v1.push_back(centered_block[vi_temp22]);
	}
}
void cad3d::get_internal_slot_v105(vector<int> &vi_all, vector<Point3D> &centered_block, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v1, int vi_temp1s1, int vi_temp1s2) {
	slot_line1.push_back(centered_block[vi_temp1s1]);
	vector<Point3D> vallpol;
	Point3D temp_slotline1pol;
	std::cout << "function: get_internal_slot_v105" << endl;
	convert_cart2pol(centered_block, vallpol);
	int rmin = find_point3d_min(vallpol, 1); //

	//find_Point3D_index(Point3D ref, vector<Point3D> &vall, vector<int> &vi_used, vector<int> &v_temp);

	//first side or edge 1
	/*int vi_temp21 = get_trianglesegment_wo_rv304(vi_all, centered_block, v_line1, v_line2, vi_temp1s1);
	std::cout << "phase 2::" << vi_temp21 << endl;
	if (vi_temp21 > 0) {
		slot_line2.push_back(centered_block[vi_temp21]);
	}

	float dtheta, drr;
	circle2d c_temp;
	//atd::cout << before fir
	std::cout << "e1:: r,   del r,    del theta,   del_r*del_theta  " << endl;
	//for (int i = 0; i < 200; i++) {
	//one side
	for (int i = 0; i < 60; i++) {
		slot_line1.push_back(centered_block[vi_temp21]);
		convert_cart2polPoint3D(centered_block[vi_temp21], temp_slotline1pol);
		vi_temp21 = get_slot_part_v101(vi_all, centered_block, v_line1, v_line2, vi_temp21);
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp21]);
		drr = temp_slotline1pol.x - vallpol[vi_temp21].x;
		std::cout << i << ": " << vallpol[vi_temp21].x << ", " << drr << ", ";
		std::cout << dtheta << ", " << dtheta * drr << endl;
		if (vi_temp21 < 0) {
			slot_line1.pop_back();
			break;
		}

		slot_line2.push_back(centered_block[vi_temp21]);
		v1.push_back(centered_block[vi_temp21]);
	}*/

	//other side or edge 2

	int vi_temp22 = get_trianglesegment_wo_rv304(vi_all, centered_block, v_line1, v_line2, vi_temp1s2);
	slot_line2.push_back(centered_block[vi_temp22]);
	//float dtheta, drr;
	//circle2d c_temp;
	/*std::cout << "e2:: r,   del r,    del theta,   del_r*del_theta  " << endl;
	for (int i = 0; i < 40; i++) {
		slot_line1.push_back(centered_block[vi_temp22]);
		convert_cart2polPoint3D(centered_block[vi_temp22], temp_slotline1pol);
		vi_temp22 = get_slot_part_v101(vi_all, centered_block, v_line1, v_line2, vi_temp22);
		dtheta = del_theta(slot_line1[slot_line1.size() - 1], centered_block[vi_temp22]);
		drr = temp_slotline1pol.x - vallpol[vi_temp22].x;
		std::cout << i << ": " << vallpol[vi_temp22].x << ", " << drr << ", ";
		std::cout << dtheta << ", " << dtheta*drr << endl;
		if (vi_temp22 < 0) {
			slot_line1.pop_back();
			break;
		}
		else if (vallpol[vi_temp22].x < (vallpol[rmin].x + 0.5)) {
			slot_line1.pop_back();
			break;
		}
		else if (drr < 0) {
			slot_line1.pop_back();
			slot_line1.pop_back();
			slot_line2.pop_back();
			v1.pop_back();
			break;
		}
		slot_line2.push_back(centered_block[vi_temp22]);
		v1.push_back(centered_block[vi_temp22]);
	}*/
}

void cad3d::get_internal_slot_v106(stl_data &block_a, vector<int> &vi_all, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<int> &vi_temp552, vector<Point3D> &v1) {
	cad3d::circle2d c1;
	int vn;
	int vi_temp1s1, vi_temp1s2, vi_temp2, vi_temp11;
	vi_temp2 = get_trianglesegment_wo_rv3031(vi_all, block_a.vertices, v_line1, v_line2, vi_temp552[0]);
	slot_line1.push_back(block_a.vertices[vi_temp552[0]]);
	slot_line2.push_back(block_a.vertices[vi_temp2]);
	v1.push_back(block_a.vertices[vi_temp2]);
	vn = v1.size();
	for (int j = 0; j < 120; j++) {
		slot_line1.push_back(block_a.vertices[vi_temp2]);
		//newblock.convert_cart2polPoint3D(centered_block[vi_temp2], vtpol1);
		vi_temp2 = get_trianglesegment_wo_rv3031(vi_all, block_a.vertices, v_line1, v_line2, vi_temp2);
		if (vi_temp2 < 0) {
			slot_line1.pop_back();
			break;
		}
		slot_line2.push_back(block_a.vertices[vi_temp2]);

		v1.push_back(block_a.vertices[vi_temp2]);
		vn++;
		c1 = compute_circle_4_point3d(v1[vn - 1], v1[vn - 2], v1[vn - 3]);
		//std::cout << c1.r << endl;
		if (c1.r < 50) {
			slot_line1.pop_back();
			slot_line2.pop_back();
			v1.pop_back();
			break;
		}
	}
}

void cad3d::get_internal_slot_v107(vector<Point3D> &centered_block, vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<int> &vi_temp552, vector<Point3D> &v1) {
	//initial ordering is done in outer radius in edge 2
	//inner radius in edge 3
	int vi_temp1s1, vi_temp1s2, vi_temp2, vi_temp11;
	float d1, d2;
	//second point for starting
	for (int i = 0; i < feat_e1.size(); i++) {
		d2 = dist201(centered_block[vi_temp552[0]], feat_e2[i]);
		//for anti clockwise  calculation
		if (d2 < sqrt(0.01)) {
			std::cout << "d2: " << i << endl;
			vi_temp2 = i;
			v1.push_back(feat_e1[i]);
			slot_e2.push_back(feat_e2[i]);
			slot_e1.push_back(feat_e1[i]);
		}
		//algorithm, going from outside to inside

	}
	float del_r;
	Point3D p_ext, p_int;
	int point_found = 0;
	/*newblock.convert_cart2polPoint3D(feat_e2[vi_temp2], aa);
	std::cout << "r, theta:" << aa.x << ", " << aa.y << endl;
	newblock.convert_cart2polPoint3D(feat_e1[vi_temp2], aa);
	std::cout << "r, theta:" << aa.x << ", " << aa.y << endl;*/
	//next points till internal 
	for (int j = 0; j < 100; j++) {
		point_found = 0;
		for (int i = 0; i < feat_e1.size(); i++) {

			d2 = dist201(v1[v1.size() - 1], feat_e2[i]);

			//for anti clockwise  calculation and outer to inner
			if (d2 < sqrt(0.01)) {
				convert_cart2polPoint3D(v1[v1.size() - 1], p_ext);
				convert_cart2polPoint3D(feat_e1[i], p_int);
				del_r = p_ext.x - p_int.x;
				/*if (del_r < 0)
					break;*/
					//std::cout << "d2: " << i << endl;
				vi_temp2 = i;
				v1.push_back(feat_e1[i]);
				slot_e2.push_back(feat_e2[i]);
				slot_e1.push_back(feat_e1[i]);
				//std::cout << "del r:" << j <<", "<< del_r << endl;
				point_found = 1;
				break;
			}
			/*else (d2 < sqrt(0.01)) {
				convert_cart2polPoint3D(v1[v1.size() - 1], p_ext);
				convert_cart2polPoint3D(feat_e1[i], p_int);
				del_r = p_ext.x - p_int.x;

				std::cout << "d2: " << i << endl;
				vi_temp2 = i;
				v1.push_back(feat_e1[i]);
				slot_e2.push_back(feat_e2[i]);
				slot_e1.push_back(feat_e1[i]);
				std::cout << "del r:" << j << ", " << del_r << endl;
			}*/
			//d2 = -1;
			/*if (point_found == 0)
				break;*/
		}
		//std::cout << point_found << endl;
		if (del_r == 0) {
			std::cout << "(point_found == 0: " << j << endl;
			break;
		}
		/*
		else if (d2 < 0) {
			std::cout << "end of line @: " << j << endl;
			break;
		}*/
	}
}
void cad3d::get_internal_slot_v108(vector<Point3D> &centered_block, vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<int> &vi_temp552, vector<Point3D> &v1) {
	//initial ordering is done in outer radius in edge 2
	//inner radius in edge 3
	int vi_temp1s1, vi_temp1s2, vi_temp2, vi_temp11;
	float d1, d2;
	//second point for starting
	for (int i = 0; i < feat_e1.size(); i++) {
		d2 = dist201(centered_block[vi_temp552[0]], feat_e2[i]);
		//for anti clockwise  calculation
		if (d2 < sqrt(0.01)) {
			std::cout << "d2: " << i << endl;
			vi_temp2 = i;
			v1.push_back(feat_e1[i]);
			slot_e2.push_back(feat_e2[i]);
			slot_e1.push_back(feat_e1[i]);
			std::cout << "found second point" << endl;
			//break;
		}
		//algorithm, going from outside to inside

	}
	
	float del_r;
	Point3D p_ext, p_int;
	int point_found = 0;
	/*newblock.convert_cart2polPoint3D(feat_e2[vi_temp2], aa);
	std::cout << "r, theta:" << aa.x << ", " << aa.y << endl;
	newblock.convert_cart2polPoint3D(feat_e1[vi_temp2], aa);
	std::cout << "r, theta:" << aa.x << ", " << aa.y << endl;*/
	//next points till internal 
	for (int j = 0; j < 60; j++) {
		point_found = 0;
		for (int i = 0; i < feat_e1.size(); i++) {

			d2 = dist201(v1[v1.size() - 1], feat_e2[i]);
			d1 = dist201(v1[v1.size() - 1], feat_e1[i]);
			
			//for anti clockwise  calculation and outer to inner
			if (d2 < sqrt(0.01)) {
				/*convert_cart2polPoint3D(v1[v1.size() - 1], p_ext);
				convert_cart2polPoint3D(feat_e1[i], p_int);
				del_r = p_ext.x - p_int.x;*/
				/*if (del_r < 0)
					break;*/
					//std::cout << "d2: " << i << endl;
				vi_temp2 = i;
				v1.push_back(feat_e1[i]);
				slot_e2.push_back(feat_e2[i]);
				slot_e1.push_back(feat_e1[i]);
				//std::cout << "del r:" << j <<", "<< del_r << endl;
				/*point_found = 1;
				break;*/
			}
			

			//for anti clockwise  calculation and outer to inner
			if (d1 < sqrt(0.01)) {
				
					//std::cout << "d2: " << i << endl;
				vi_temp2 = i;
				v1.push_back(feat_e2[i]);
				slot_e2.push_back(feat_e1[i]);
				slot_e1.push_back(feat_e2[i]);
				//std::cout << "del r:" << j <<", "<< del_r << endl;
				
			}
			//d2 = -1;
			
		}
		//std::cout << point_found << endl;
		
		
	}
}
void cad3d::get_internal_slot_v107_phase2(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);

			if (d1 < 0.01) {
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
				jr1.push_back(j);
			}
			if (d2 < 0.01) {
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
				jrmat2.push_back(p2_con.x);
				jr2.push_back(j);
				slot_e3.push_back(feat_e2[j]);
				v2.push_back(feat_e1[j]);
				slot_e4.push_back(feat_e1[j]);
			}
		}
		if (jr2.size() < 0) {
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	a = jr1.size() == 1;
	b = jr2.size() == 1;
	/*if ((jr1.size() == 1) & jr2.size() == 1) {
		if (p1_con.x <= p2_con.x) {
			v2.push_back(feat_e2[jr1[0]]);
		}
		else {
			v2.push_back(feat_e1[jr2[0]]);
		}
	}*/
	if (jr2.size() == 1) {
		v2.push_back(feat_e1[jr2[0]]);
	}
	/*for (int i = 0; i < jr2.size(); i++) {
	}*/
	/*int jd2min = find_1d_vector_min(jrmat2, 1);
	v_out.push_back(v_line1[jr2[jd2min]]);*/
	/*int vi_temp1s1, vi_temp1s2, vi_temp2, vi_temp11;
	float d1, d2;
	d2 = 1;
	//second point for starting
	for (int i = 0; i < feat_e1.size(); i++) {
		d1 = dist201(v1[v1.size() - 1], feat_e1[i]);
		//for anti clockwise  calculation
		if (d1 < sqrt(0.01)) {
			std::cout << "d2: " << i << endl;
			vi_temp2 = i;
			v2.push_back(feat_e1[i]);
			slot_e4.push_back(feat_e1[i]);
			slot_e3.push_back(feat_e2[i]);
		}

	}
	float del_r;
	Point3D p_ext, p_int;
	int point_found = 0;*/


}
void cad3d::get_internal_slot_v107_phase3(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	for (int i = 0; i < 40; i++) {
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);

			if (d1 < 0.01) {
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
				jrmat1.push_back(p1_con.x);
				slot_e3.push_back(feat_e1[j]);
				v2.push_back(feat_e1[j]);
				slot_e4.push_back(feat_e2[j]);
			}
			if (d2 < 0.01) {
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				slot_e3.push_back(feat_e2[j]);
				v2.push_back(feat_e1[j]);
				slot_e4.push_back(feat_e1[j]);
			}
		}
		if (jr1.size() < 0) {
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	a = jr1.size() == 1;
	b = jr2.size() == 1;

	if (jr1.size() == 1) {
		v2.push_back(feat_e2[jr1[0]]);
	}

	/*if (jr2.size() == 1) {
		v2.push_back(feat_e1[jr2[0]]);
	}
	if ((jr1.size() == 1) & jr2.size() == 1) {
		if (p1_con.x <= p2_con.x) {
			v2.push_back(feat_e2[jr1[0]]);
		}
		else {
			v2.push_back(feat_e1[jr2[0]]);
		}
	}*/
}
void cad3d::get_internal_slot_v107_phase4(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	for (int i = 0; i < 100; i++) {
		std::cout << i << ": ";
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);
			
			if (d1 < 0.1) {
				std::cout <<" d1- " << d1 ;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
				jrmat1.push_back(p1_con.x);
				slot_e3.push_back(feat_e1[j]);
				v2.push_back(feat_e2[j]);
				slot_e4.push_back(feat_e2[j]);
			}
			if (d2 < 0.1) {
				std::cout << " d2- " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				slot_e3.push_back(feat_e2[j]);
				v2.push_back(feat_e1[j]);
				slot_e4.push_back(feat_e1[j]);
			}
		}
		std::cout << "~ "<<jr1.size() << ", " << jr2.size() << endl;
		if ((jr1.size() < 0)|| (jr2.size() < 0)) {
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	a = jr1.size() == 1;
	b = jr2.size() == 1;

	if (jr1.size() == 1) {
		v2.push_back(feat_e2[jr1[0]]);
	}

}
void cad3d::get_internal_slot_v107_phase5(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	//function has issues of when to stop
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	int vn = v2.size();
	
	for (int i = 0; i < 40; i++) {
		std::cout << i << ": " <<v2.size() ;
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);
			
			if (d1 < 0.1) {
				std::cout << " j d1- " << j<<", "<<d1;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
				//point_repeat_check_v1(feat_e2[j], v2);
				jrmat1.push_back(feat_e2[j].x);
				slot_e3.push_back(feat_e1[j]);
				v2.push_back(feat_e2[j]);
				slot_e4.push_back(feat_e2[j]);
				//break;
				//vn++;
			}
			else if (d2 < 0.1) {
				std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				//point_repeat_check_v1(feat_e1[j], v2);
				slot_e3.push_back(feat_e2[j]);
				v2.push_back(feat_e1[j]);
				slot_e4.push_back(feat_e1[j]);
				//break;
				//vn++;
			}
			vn++;
		}
		std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		/*if ((jr1.size() < 1) || (jr2.size() < 1)) {
			break;
		}*/
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	/*a = jr1.size() == 1;
	b = jr2.size() == 1;

	if (jr1.size() == 1) {
		v2.push_back(feat_e2[jr1[0]]);
	}
	else if (jr2.size() == 1) {
		v2.push_back(feat_e1[jr2[0]]);
	}*/

}
void cad3d::get_internal_slot_v107_phase6(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	//function has issues of when to stop
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	//Point3D p_start_dif;
	int vn = v2.size();
	int point_reap_val = 0;
	vector<float> int_index_check;
	//for (int i = 0; i < 100; i++) {
	for (int i = 0; i < 40; i++) {
		//std::cout << i << ": " << v2.size();
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);

			if (d1 < 0.01) {
				//std::cout << " j d1- " << j << ", " << d1;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				std::cout << "d1: " << j <<".."<<p1_con.x<< endl;
				jrmat1.push_back(p1_con.x);
				point_reap_val = point_repeat_check_v1(feat_e2[j], v2);
				if (point_reap_val==0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					v2.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e2[j]);
					break;
				}
				//
				//vn++;
			}
			else if (d2 < 0.01) {
				//std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				std::cout << "d2: " << j << ".." << p2_con.x << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				point_reap_val = point_repeat_check_v1(feat_e1[j], v2);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					v2.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e1[j]);
					break;
				}
				//vn++;
			}
			vn++;
			point_reap_val = 0;
		}
		//std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		//std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		if ((jr1.size() < 1) & (jr2.size() < 1)) {
			std::cout << "discontinuity" << endl;
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	int i_temp = find_1d_vector_max(int_index_check,1);
	std::cout << "discontinuity max: " << int_index_check [i_temp]<<endl;
	//end closed region
	slot_e3.push_back(v2[v2.size()-1]);
	//v2.push_back(feat_e2[j]);
	slot_e4.push_back(v2[0]);

}
void cad3d::get_internal_slot_v107_phase7(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	std::cout << "function internal slot\n ";
	//update done from phase 3
	//the direction will always be polar inward
	//still based on tolerance level of d1 and d2, need to figure it out
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con, pbase_con;
	float delz;
	for (int i = 0; i < 40; i++) {
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v2[v2.size() - 1], feat_e2[j]);
			convert_cart2polPoint3D(v2[v2.size() - 1], pbase_con);

			if (d1 < 0.1) {
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				delz = abs(v2[v2.size() - 1].z - feat_e1[j].z);
				//std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
				//if (p1_con.x < pbase_con.x) {
				if ((p1_con.x < pbase_con.x) &(delz < 0.6)){
					
					//std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
					v2.push_back(feat_e2[j]);
					jr1.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e2[j]);
				}
				/*
				

				jrmat1.push_back(p1_con.x);
				*/
			}
			if (d2 < 0.1) {
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				delz = abs(v2[v2.size() - 1].z - feat_e2[j].z);
				//std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
				//if (p2_con.x < pbase_con.x) {
				if ((p2_con.x < pbase_con.x) &(delz < 0.6)) {
					
					//std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
					v2.push_back(feat_e1[j]);
					jr2.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e1[j]);
				}

				/*
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				*/
			}

			/*if ((jr1.size() < 1) & (jr2.size() < 1)){
				break;
				//v2.push_back(feat_e2[jr1[0]]);
			}*/

		}
		/*if ((jr1.size() < 1) & (jr2.size() < 1)) {
			break;
			//v2.push_back(feat_e2[jr1[0]]);
		}*/
		jr1.resize(0);
		jr2.resize(0);
	}
	/*if (jr1.size() < 0) {
				break;
		}
		
		jrmat1.resize(0);
		jrmat2.resize(0);
		}
		bool a, b;
		a = jr1.size() == 1;
		b = jr2.size() == 1;

		if (jr1.size() == 1) {
			v2.push_back(feat_e2[jr1[0]]);
		}


	}*/
}

int cad3d::point_repeat_check_v1(Point3D p1, vector<Point3D> &vin) {
	//at 2d level repeat
	int repeat_value = 0;
	Point3D p_dist;
	for (int i = 0; i < vin.size(); i++) {
		p_dist = feature_edge_predict_phase1(p1,vin[i]);
		if (p_dist.x <0.1 & p_dist.y <0.1) {
			//std::cout << i << "repeat found"  << endl;
			repeat_value = 1;
			break;
		}
	}
	return repeat_value;
}

int cad3d::point_repeat_check_v2(Point3D p1, vector<Point3D> &vin) {
	//at 3d level repeat
	int repeat_value = 0;
	Point3D p_dist;
	
	float dist;
	bool a = 0;
	bool b = 0;
	bool c = 0;
	for (int i = 0; i < vin.size(); i++) {
		p_dist = feature_edge_predict_phase1(p1, vin[i]);
		//dist = dist201(p1, vin[i]);
		a = p_dist.x < 0.1;
		b = p_dist.y < 0.1;
		c = p_dist.z < 0.1;
		if (a & b & c) {
			//std::cout << i << "repeat found"  << endl;
			repeat_value = 1;
			break;
		}
	}
	/*for (int i = 0; i < vin.size(); i++) {
		p_dist = feature_edge_predict_phase1(p1, vin[i]);
		if (p_dist.x <0.1 & p_dist.y <0.1) {
			//std::cout << i << "repeat found"  << endl;
			repeat_value = 1;
			break;
		}
	}*/
	return repeat_value;
}
int cad3d::point_superposition_check_v1(Point3D p1, vector<Point3D> &vin, float tol) {
	// at 3d level repeat
	//at 2d level repeat
	int repeat_value = 0;
	Point3D p_dist;
	float dist;
	for (int i = 0; i < vin.size(); i++) {
		p_dist = feature_edge_predict_phase1(p1, vin[i]);
		dist = dist201(p1, vin[i]);
		if (p_dist.x <tol & p_dist.y <tol) {
		//if (dist < tol) {
			//std::cout << i << "repeat found"  << endl;
			repeat_value = 1;
			break;
		}
	}
	return repeat_value;
	/*int repeat_value = 0;
	Point3D p_dist;
	float dist;
	
	for (int i = 0; i < vin.size(); i++) {
		//p_dist = feature_edge_predict_phase1(p1, vin[i]);
		dist = dist201(p1, vin[i]);
		
		if (dist < tol) {
			//std::cout << i << "repeat found"  << endl;
			repeat_value = 1;
			break;
		}
	}
	return repeat_value;*/
}

void cad3d::getplane_simplfied_block(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i) {
	int x_max_i = find_point3d_max(block_a.vertices, 1);
	int y_max_i = find_point3d_max(block_a.vertices, 2);
	int z_max_i = find_point3d_max(block_a.vertices, 3);
	int x_min_i = find_point3d_min(block_a.vertices, 1);
	int y_min_i = find_point3d_min(block_a.vertices, 2);
	int z_min_i = find_point3d_min(block_a.vertices, 3);
	std::cout << "z_max_value: " << block_a.vertices[z_max_i].z << endl;

	
	vector<int> slot_plane_i;
	//centered_block = block_a.vertices;
	find_point_3d_1feature(slot_plane, block_a.vertices, block_a.vertices[z_max_i].z, 0.05, 3);
	slot_plane_i.resize(block_a.vertices.size());
	find_point3d_i_1feature(slot_plane_i, block_a.vertices, block_a.vertices[z_max_i].z, 0.05, 3);
	show_index_frequency(slot_plane_i);

	vector<cad3d::triangle_index> t_slot, ti_all,  tslot_i1, tslot_i2, tslot_i3;
	triangle_vertexindex_init(block_a.triangles, ti_all);
	vector<cad3d::triangle> tslot;
	std::cout << "total triangle index " << ti_all.size() << endl;
	//triangle and line coordination
	float zmid = (block_a.vertices[z_max_i].z + block_a.vertices[z_min_i].z) / 2;
	triangle_allpoint_1feature(block_a.triangles, tslot, block_a.vertices[z_max_i].z, 0.5, 3);
	std::cout << "tslot size" << tslot.size() << endl;
	triangle_allpoint_i1feature(block_a.vertices, ti_all, tslot_i, block_a.vertices[z_max_i].z, 0.05, 3);
}
void cad3d::getplane_simplfied_block_v2(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i, float z_temp) {
	/*int x_max_i = find_point3d_max(block_a.vertices, 1);
	int y_max_i = find_point3d_max(block_a.vertices, 2);
	int z_max_i = find_point3d_max(block_a.vertices, 3);
	int x_min_i = find_point3d_min(block_a.vertices, 1);
	int y_min_i = find_point3d_min(block_a.vertices, 2);
	int z_min_i = find_point3d_min(block_a.vertices, 3);
	std::cout << "z_max_value: " << block_a.vertices[z_max_i].z << endl;*/


	vector<int> slot_plane_i;
	//centered_block = block_a.vertices;
	find_point_3d_1feature(slot_plane, block_a.vertices, z_temp, 0.05, 3);
	slot_plane_i.resize(block_a.vertices.size());
	find_point3d_i_1feature(slot_plane_i, block_a.vertices, z_temp, 0.05, 3);
	show_index_frequency(slot_plane_i);

	vector<cad3d::triangle_index> t_slot, ti_all, tslot_i1, tslot_i2, tslot_i3;
	triangle_vertexindex_init(block_a.triangles, ti_all);
	vector<cad3d::triangle> tslot;
	std::cout << "total triangle index " << ti_all.size() << endl;
	//triangle and line coordination
	//float zmid = (block_a.vertices[z_max_i].z + block_a.vertices[z_min_i].z) / 2;
	triangle_allpoint_1feature(block_a.triangles, tslot, z_temp, 1.0, 3);
	std::cout << "tslot size" << tslot.size() << endl;
	triangle_allpoint_i1feature(block_a.vertices, ti_all, tslot_i, z_temp, 1.0, 3);
}
void cad3d::one_start_point_simplified_block(stl_data &block_a, vector<int> &vi_all, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, Point3D pref_start, vector<int> &vi_temp551, vector<int> &vi_temp552) {
	
	//vector<cad3d::Point3D> slot_line1, slot_line2, v1, vinward1;
	vector<int> jd1, jd2;
	find_Point3D_index(pref_start, block_a.vertices, vi_all, vi_temp551);
	std::cout << vi_temp551.size() << "- " << vi_temp551[0] << endl;
	float d1, d2;
	vector<float> jdmat1, jdmat2;
	for (int j = 0; j < v_line1.size(); j++) {
		d1 = dist201(pref_start, v_line1[j]);//v_ext_l2[0]
		d2 = dist201(pref_start, v_line2[j]);
		if (d1 == 0) {
			std::cout << "d1: " << j << "- " << dist201(pref_start, v_line2[j]) << endl;
			jd1.push_back(j);
			jdmat1.push_back(dist201(pref_start, v_line2[j]));
		}
		if (d2 == 0) {
			//find the points in which they are connected
			std::cout << "d2: " << j << "- " << dist201(pref_start, v_line1[j]) << endl;
			jd2.push_back(j);
			jdmat2.push_back(dist201(pref_start, v_line1[j]));
		}
	}
	std::cout << jd2.size() << ", " << jdmat2.size() << endl;
	//corner detection, not optimum for single basic block
	if (jd2.size() > 1) {
		if (jdmat2[1] < jdmat2[0]) {
			find_Point3D_index(v_line1[jd2[1]], block_a.vertices, vi_all, vi_temp552);
		}
		else
		{
			find_Point3D_index(v_line1[jd2[0]], block_a.vertices, vi_all, vi_temp552);
		}



	}
	if (jd1.size() == 1 & jd2.size() == 1) {
		if (jdmat1[0] < jdmat2[0]) {
			find_Point3D_index(v_line2[jd1[0]], block_a.vertices, vi_all, vi_temp552);
		}
		else
		{
			find_Point3D_index(v_line1[jd2[0]], block_a.vertices, vi_all, vi_temp552);
		}
	}
	if (jd2.size() == 1) {

		find_Point3D_index(v_line1[jd2[0]], block_a.vertices, vi_all, vi_temp552);

	}
	if (jd1.size() == 1) {

		find_Point3D_index(v_line2[jd1[0]], block_a.vertices, vi_all, vi_temp552);

	}
}

void cad3d::cross_multiply_point3d(Point3D p1, Point3D p2, Point3D p3, Point3D &p_normal) {
	//https://www.sanfoundry.com/cpp-program-compute-cross-product-two-vectors/
	//u1i +u2j + u3k = (p1-p2);
	//v1i +v2j + v3k = (p2-p3);
	float u1, u2, u3, v1, v2, v3;
	u1 = p1.x - p2.x; u2 = p1.y - p2.y; u3 = p1.z - p2.z;
	v1 = p2.x - p3.x; v2 = p2.y - p3.y; v3 = p2.z - p3.z;
	p_normal.x = u2 * v3 - v2 * u3;
	p_normal.y = v1 * u3 - u1 * v3;
	p_normal.y = u1 * v2 - v1 * u2;

}
void cad3d::interpolate_lagrange_v5(std::vector<Point3D> &v1, int ii, int n1, Point3D &v2) {
	/*modify according to the 3d point structures
	takes an array of 4 points and x
	returns the point 3d where y in interpolated
	n1 is the number of point size in here
	this one performs index wise interpolation using header file spline.h*/

	//here ii works as xi, so index based

	float xo = 0; // Initialize result 
	float yo = 0;
	float zo = 0;

	for (int i = 0; i < n1; i++) {
		// Compute individual terms of above formula 
		float termx = v1[i].x;
		float termy = v1[i].y;
		float termz = v1[i].z;
		//std::cout << termx << ", " << termy<< endl;
		for (int j = 0; j < n1; j++) {
			if (j != i) {
				termx = termx * float(ii - j) / double(i - j);
				termy = termy * float(ii - j) / double(i - j);
				termz = termz * float(ii - j) / double(i - j);
			}
		}

		// Add current term to result 
		xo += termx;
		yo += termy;
		zo += termz;
	}

	v2.x = xo;
	v2.y = yo;
	v2.z = zo;
	//std::cout << xo << ", " << yo << ", ";
	//return result;
	//https://www.geeksforgeeks.org/lagranges-interpolation/
}
float cad3d::segment_length_measure(vector<Point3D> &vin) {
	//assumed that vin is already ordered
	float total_dist = 0;
	float temp_d;
	Point3D v1;
	for (int i = 1; i < vin.size(); i++) {
		v1 = feature_edge_predict_phase1(vin[i - 1], vin[i]);
		temp_d = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
		total_dist = total_dist + temp_d;
	}
	return total_dist;
}
void cad3d::interpolate_lagrange_v6(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present) {
	/* float ii is the comparative position
	// test with n1 for 3 points
	v1 is input array
	*/

	//here ii works as xi, so index based

	//initializing output result
	float xo = 0; // Initialize result 
	float yo = 0;
	float zo = 0;

	// debug, i should be present index ans should not start from 0 only
	for (int i = i_present; i <( i_present+n1); i++) {
		// Compute individual terms of above formula 
		float termx = v1[i].x;
		float termy = v1[i].y;
		float termz = v1[i].z;
		//std::cout << termx << ", " << termy<< endl;
		for (int j = 0; j < n1; j++) {
			if (j != i) {
				termx = termx * float(ii - j) / double(i - j);
				termy = termy * float(ii - j) / double(i - j);
				termz = termz * float(ii - j) / double(i - j);
			}
		}

		// Add current term to result 
		xo += termx;
		yo += termy;
		zo += termz;
	}

	v2.x = xo;
	v2.y = yo;
	v2.z = zo;
	//std::cout << xo << ", " << yo << ", ";
	//return result;
	//https://www.geeksforgeeks.org/lagranges-interpolation/
}
void cad3d::interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present) {
	/* float ii is the comparative position
	// test with n1 for 3 points
	v1 is input array
	i_present at which it will start taking input
	ii index at which interpolation needs to be found
	*/
	//std::cout << "function spline:" << i_present<< ", " << ii <<", " << n1<<endl;
	//here ii works as xi, so index based

	//initializing output result
	float xo = 0; // Initialize result 
	float yo = 0;
	float zo = 0;
	//vector<double> xx(5), yy(5), zz(5), i_index(5);
	vector<double> xx, yy, zz, i_index;
	/*td::vector<double> X(5), Y(5);
   X[0]=0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=2.0;
   Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;

   tk::spline s;
   s.set_points(X,Y);    // currently it is required that X is already sorted

   double x=1.5;

   printf("spline at %f is %f\n", x, s(x));*/
	// debug, i should be present index ans should not start from 0 only
	//for (int i = i_present; i < (i_present + n1); i++) {
	for (int i = 0; i < n1; i++) {
		/*i_index[i] = i_present + n1;
		xx[i] = v1[i + i_present].x;
		yy[i] = v1[i + i_present].y;
		zz[i] = v1[i + i_present].z;*/
		i_index.push_back(i_present + i);
		xx.push_back(v1[i + i_present].x);
		yy.push_back(v1[i + i_present].y);
		zz.push_back(v1[i + i_present].z);
	}

	tk::spline s1, s2, s3;
	s1.set_points(i_index,xx);
	s2.set_points(i_index, yy);
	s3.set_points(i_index, zz);
	//std::cout << "setting points done" << endl;
	v2.x = s1(ii);
	v2.y = s2(ii);
	v2.z = s3(ii);
	//std::cout << ii << endl;
	//std::cout << v2.x << ", " << v2.y << endl;
	//return result;
	//https://www.geeksforgeeks.org/lagranges-interpolation/
}

void cad3d::triangle_normal_update(vector<triangle> &t_all, std::vector<triangle_index> &t_i) {
	Point3D p1, p2, p3, pn;
	int ii;
	for (int i = 0; i < t_i.size();i++) {
		ii = t_i[i].ti;
		p1 = t_all[ii].v1;
		p2 = t_all[ii].v2;
		p3 = t_all[ii].v3;
		cross_multiply_point3d(p1, p2, p3, pn);
		t_all[ii].normal = pn;
	}
}
cad3d::triangle cad3d::generate_triangle(Point3D p1, Point3D p2, Point3D p3) {
	//triangle t_out;
	Point3D pn;
	cross_multiply_point3d(p1, p2, p3, pn);
	triangle t_out(pn, p1, p2, p3);
	t_out.v1 = copy_point3d(p1);
	t_out.v2 = copy_point3d(p2);
	t_out.v3 = copy_point3d(p3);
	t_out.normal = pn;
	return t_out;
}
cad3d::Point3D cad3d::copy_point3d(Point3D pin) {
	Point3D p_out;
	p_out.x = pin.x;
	p_out.y = pin.y;
	p_out.z = pin.z;

	return p_out;

}
void cad3d::update_point3d_vector(vector<Point3D> v_in, vector<Point3D> &v_out) {
	// for updating a point3d matrix into another one 
	//09.10.19
	for (int i = 0; i < v_in.size(); i++) {
		v_out.push_back(v_in[i]);
	}
}
cad3d::Point3D cad3d::feature_edge_predict_phase1(Point3D p1, Point3D p2) {
	Point3D line_values;
	line_values.x = abs(p1.x - p2.x);
	line_values.y = abs(p1.y - p2.y);
	line_values.z = abs(p1.z - p2.z);
	return line_values;
}
int cad3d::find_line_repete_index(Point3D p1, Point3D p2,std::vector<Point3D> &v_line1, std::vector<Point3D> &v_line2, int present_i, std::vector<int> &repete_i) {
	//if returns 0 no repeat if >0 repeat occured cannot be pushed in 
	Point3D line_value1, line_value2;
	int total_repeat=0;
	bool a, b, c, d;
	//for (int i = present_i+1; i < v_line1.size(); i++) {
	for (int i = 0; i < v_line1.size(); i++) {
		if (i != present_i) {
			line_value1 = feature_edge_predict_phase1(p1, v_line1[i]);
			line_value2 = feature_edge_predict_phase1(p2, v_line2[i]);
			a = line_value1.x == 0.0;
			b = line_value2.x == 0.0;
			c = line_value1.y == 0.0;
			d = line_value2.y == 0.0;
			if ((a & b)&(c & d)) {
				total_repeat++;
				repete_i.push_back(i);
			}
		}
	}
	return total_repeat;
}

int cad3d::find_any_int(std::vector<int> &v_i, int test_val) {
	int val_exist = 0;
	int del_i=0;
	for (int i = 0; i < v_i.size(); i++) {
		del_i == test_val - v_i[i];
		if (del_i == 0) {
			val_exist++;
		}
	}
	return val_exist;
}
void cad3d::feature_edge_predict_phase2(vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2) {
	//takes input of ordered lines vline1,2
	//gives output to feature edges without any repetition
	//performs xy planer 2d comparison
	int repeat_status =0;
	int find_status = 0;
	vector<int> repeat_temp, repeat_all;
	repeat_status = find_line_repete_index(v_line1[0], v_line2[0], v_line1, v_line2, 0, repeat_temp);
	if (repeat_status == 0) {
		slot_line1.push_back(v_line1[0]);
		slot_line2.push_back(v_line2[0]);
		//std::cout << 0<<": NR\n";
	}
	else {
		repeat_all.push_back(0);
		for (int j = 0; j < repeat_temp.size(); j++) {
			repeat_all.push_back(repeat_temp[j]);
		}
		repeat_temp.resize(0);
	}
	for (int i = 1; i < v_line1.size(); i++) {
		repeat_status = find_line_repete_index(v_line1[i], v_line2[i], v_line1, v_line2, i, repeat_temp);
		if (repeat_status == 0) {
			//std::cout << i << ": NR\n";
			slot_line1.push_back(v_line1[i]);
			slot_line2.push_back(v_line2[i]);
			find_status = find_any_int(repeat_all,i);
			/*if (find_status == 0) {
				slot_line1.push_back(v_line1[i]);
				slot_line2.push_back(v_line2[i]);
			}*/
		}
		else {
			repeat_all.push_back(i);
			for (int j = 0; j < repeat_temp.size(); j++) {
				repeat_all.push_back(repeat_temp[j]);
			}
			repeat_temp.resize(0);
		}
		repeat_status = 0;
		find_status = 0;
	}
}

void cad3d::feature_edge_predict_phase0(vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2) {
	cad3d::Point3D line_value;
	vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	bool a, b;
	//ordering lines
	//for rectangular block
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = v_line2[i].x - v_line1[i].x;
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = v_line2[i].y - v_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}

	for (int i = 0; i < vord1_line1.size(); i++) {
		//for (int i = 0; i < tslot_i.size(); i++) {
		//for (int i = 0; i < 100; i++) {
		line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
		/*a = line_value.x < 20.0;
		b = line_value.y < 20.0;
		//b = 0;
		if (a || b) {
			if (sqrt(line_value.x*line_value.x + line_value.y*line_value.y) < 100) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}

			//std::cout << sqrt(line_value.x*line_value.x + line_value.y*line_value.y) << endl;
		}*/
		if (line_value.x == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
		if (line_value.y == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}

	}
	feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);
	
}
void cad3d::order_cartesian_points_v1(vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &vord1_line1, vector<Point3D> &vord1_line2) {
	cad3d::Point3D line_value;
	//vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	bool a, b;
	//ordering lines
	//for rectangular block
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = v_line2[i].x - v_line1[i].x;
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = v_line2[i].y - v_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}

	/*for (int i = 0; i < vord1_line1.size(); i++) {
		//for (int i = 0; i < tslot_i.size(); i++) {
		//for (int i = 0; i < 100; i++) {
		line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
		
		if (line_value.x == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
		if (line_value.y == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}

	}
	feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);*/

}
void cad3d::order_polar_points_v1(vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &vord1_line1, vector<Point3D> &vord1_line2, vector<Point3D> &vpolar_line1, vector<Point3D> &vpolar_line2) {
	cad3d::Point3D line_value;
	//vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	bool a, b;
	//ordering lines
	//for rectangular block
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = vpolar_line2[i].x - vpolar_line1[i].x;
		//radius comparison
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = vpolar_line2[i].y - vpolar_line1[i].y;
			//angle comparison
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}

}
void cad3d::cartesean_processing(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2) {
	getplane_simplfied_block(block_a, slot_plane, tslot_i);
	triangle_normal_update(block_a.triangles, tslot_i);
	int ii;
	bool a, b;
	
	vector<cad3d::Point3D> v_temp, v_line1, v_line2;
	vector<int> i_line;
	triangle2line(block_a.vertices, tslot_i, v_line1, v_line2, i_line, v_temp);
	//newblock.triangle2line(block_a.vertices, tbi_1, v_line1, v_line2, i_line, v_temp);
	/*std::cout << "triangle to line completed " << endl;
	std::cout << "tslot_i size: " << tslot_i.size() << endl;
	std::cout << "v_line size:/3 " << v_line1.size() / 3 << endl;*/
	cad3d::Point3D line_value;

	vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	//ordering lines
	//for rectangular block
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = v_line2[i].x - v_line1[i].x;
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = v_line2[i].y - v_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}

	for (int i = 0; i < vord1_line1.size(); i++) {
		//for (int i = 0; i < tslot_i.size(); i++) {
		//for (int i = 0; i < 100; i++) {
		line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
		a = line_value.x < 20.0;
		b = line_value.y < 20.0;
		//b = 0;
		if (a || b) {
			if (sqrt(line_value.x*line_value.x + line_value.y*line_value.y) < 100) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}

			//std::cout << sqrt(line_value.x*line_value.x + line_value.y*line_value.y) << endl;
		}
		if (line_value.x == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
		if (line_value.y == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}

	}
	feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);
}
int cad3d::feature_plane_approximation(stl_data &block_a, vector<Point3D> &zslice_points, float &z_now) {
	//returns 1 if feature plane found or returns 0 or -1;

	int z_min_i = find_point3d_min(block_a.vertices, 3);
	int z_max_i = find_point3d_max(block_a.vertices, 3);

	float z_init = block_a.vertices[z_min_i].z;
	float del_z = 1;
	//should be minimum to maximum points
	for (int i = 0; i < 100; i++) {
		if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].z)
			break;
		find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 3);
		std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
		zslice_points.resize(0);
	}

	std::cout << "z value:";
	cin >> z_now;

	find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 3);
	std::cout << "points in the slice plane: " << zslice_points.size() << endl;
	if (zslice_points.size() < 1) {
		std::cout << "No points centering this value\n";
		return -1;
	}
	else
		return 1;
	//std::cout << "y_min_value: " << block_a.vertices[y_min_i].y << endl;
	/*target_plot_update(block_a.vertices);
	target_plot_color_update(0.0, 0.8, 0.0, block_a.vertices.size());*/
}
int cad3d::axis_wise_feature_plane_approximation(stl_data &block_a, vector<Point3D> &zslice_points, float &z_now, int dim) {
	//returns 1 if feature plane found or returns 0 or -1;

	if (dim == 1) {int z_min_i = find_point3d_min(block_a.vertices, 1);
	int z_max_i = find_point3d_max(block_a.vertices, 1);

	float z_init = block_a.vertices[z_min_i].x;
	float del_z = 1;
	//should be minimum to maximum points
	for (int i = 0; i < 100; i++) {
		if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].x)
			break;
		find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 1);
		std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
		zslice_points.resize(0);
	}

	std::cout << "x value:";
	cin >> z_now;

	find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 1);
	std::cout << "points in the slice plane: " << zslice_points.size() << endl;
	if (zslice_points.size() < 1) {
		std::cout << "No points centering this value\n";
		return -1;
	}
	else
		return 1;
	}
	
	else if (dim ==2) {
		int z_min_i = find_point3d_min(block_a.vertices, 2);
		int z_max_i = find_point3d_max(block_a.vertices, 2);

		float z_init = block_a.vertices[z_min_i].y;
		float del_z = 1;
		//should be minimum to maximum points
		for (int i = 0; i < 100; i++) {
			if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].y)
				break;
			find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 2);
			std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
			zslice_points.resize(0);
		}

		std::cout << "y value:";
		cin >> z_now;

		find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 1);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}

	else if (dim==3) {
		int z_min_i = find_point3d_min(block_a.vertices, 3);
		int z_max_i = find_point3d_max(block_a.vertices, 3);

		float z_init = block_a.vertices[z_min_i].z;
		float del_z = 1;
		//should be minimum to maximum points
		for (int i = 0; i < 100; i++) {
			if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].z)
				break;
			find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 3);
			std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
			zslice_points.resize(0);
		}

		std::cout << "z value:";
		cin >> z_now;

		find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 3);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}
	else
		return -1;
	//std::cout << "y_min_value: " << block_a.vertices[y_min_i].y << endl;
	/*target_plot_update(block_a.vertices);
	target_plot_color_update(0.0, 0.8, 0.0, block_a.vertices.size());*/
}
void cad3d::apply_interpolation(vector<Point3D> &v_in, vector<Point3D> &v_out, int point_no) {
	//using lagrange v3 and x,y = f(i)
	point_no = 40; // number of points in centerline
	//v_out.resize(point_no + 1);
	//1st point
	v_out.push_back(v_in[0]);
	//v_out[0] = v_in[0];
	int kk = 1;
	int i_tol = 4;
	vector<Point3D> temp1, temp2;
	vector<Point3D> templ2, tempu2;// , low_int1, up_int1;
	Point3D t1, t2, p1, p2, pdist;
	float del_s = segment_length_measure(v_in) / point_no;
	std::cout << "function: apply interpoaltion" << endl;
	//std::cout << "number of points, total length: " << v_in.size();
	//std::cout <<", "<< segment_length_measure(v_in) <<endl;
	float d_temp;
	p1 = v_in[0];
	int jinit = 0;
	int jtemp, j_present;
	//second point
	jtemp = jinit;
	for (int j_i = 0; j_i < v_in.size(); j_i++) {
		//for (int j_i = 6; j_i < 13; j_i++) {
		p1 = v_in[jinit];
		std::cout << jinit << ": ";
		for (int jj = 1; jj <= 4; jj++) {
			jtemp = jinit + jj;
			j_present = jj;
			if (jinit + jj == v_in.size()) {
				//std::cout << jinit + jj << endl;
				break;
			}
			jtemp = jinit + jj;
			templ2.push_back(v_in[jinit + jj]);
			p2 = v_in[jinit + jj];
			pdist = feature_edge_predict_phase1(p1, p2);
			//d_temp = sqrt(pdist.x*pdist.x + pdist.y*pdist.y + pdist.z*pdist.z);

		}
		if (jtemp == v_in.size()) {
			break;
		}
		d_temp = sqrt(pdist.x*pdist.x + pdist.y*pdist.y + pdist.z*pdist.z);
		//std::cout << d_temp << "~ 8" << endl;
		jinit = jtemp;
		//std::cout << templ2.size();
		interpolate_lagrange_v5(templ2, round(templ2.size() / 2), templ2.size(), t1);
		templ2.resize(0);
		v_out.push_back(t1);
	}
	interpolate_lagrange_v5(templ2, round(templ2.size() / 2), templ2.size(), t1);
	v_out.push_back(t1);
	templ2.resize(0);
	v_out.push_back(v_in[v_in.size() - 1]);

}

int cad3d::feature_axis_selection(stl_data &block_a, vector<Point3D> &zslice_points, float &z_now, int dim) {
	//returns 1 if feature plane found or returns 0 or -1;
	
	vector<float> zslice_point_count;
	vector<float> zslice_value;
	int zslice_max;
	if (dim == 1) {
		int z_min_i = find_point3d_min(block_a.vertices, 1);
		int z_max_i = find_point3d_max(block_a.vertices, 1);

		float z_init = block_a.vertices[z_min_i].x;
		float del_z = 1;
		//should be minimum to maximum points
		for (int i = 0; i < 100; i++) {
			if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].x)
				break;
			find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 1);
			std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
			zslice_point_count.push_back(zslice_points.size());
			zslice_value.push_back(z_init + 0.5 + i * del_z);
			zslice_points.resize(0);
		}

		
		zslice_max = find_1d_vector_max(zslice_point_count,1);
		z_now = zslice_value[zslice_max];
		//cin >> z_now;
		std::cout << "x value:" << z_now << endl;

		find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 1);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}

	else if (dim == 2) {
		int z_min_i = find_point3d_min(block_a.vertices, 2);
		int z_max_i = find_point3d_max(block_a.vertices, 2);

		float z_init = block_a.vertices[z_min_i].y;
		float del_z = 1;
		//should be minimum to maximum points
		for (int i = 0; i < 100; i++) {
			if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].y)
				break;
			find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 2);
			std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
			zslice_point_count.push_back(zslice_points.size());
			zslice_value.push_back(z_init + 0.5 + i * del_z);
			zslice_points.resize(0);
		}

		//cin >> z_now;
		zslice_max = find_1d_vector_max(zslice_point_count, 1);
		z_now = zslice_value[zslice_max];
		std::cout << "y value:" << z_now << endl;

		find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 1);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}

	else if (dim == 3) {
		int z_min_i = find_point3d_min(block_a.vertices, 3);
		int z_max_i = find_point3d_max(block_a.vertices, 3);

		float z_init = block_a.vertices[z_min_i].z;
		float del_z = 1;
		//should be minimum to maximum points
		for (int i = 0; i < 100; i++) {
			if (z_init + 0.5 + i * del_z > block_a.vertices[z_max_i].z)
				break;
			find_point_3d_1feature(zslice_points, block_a.vertices, z_init + 0.5 + i * del_z, 0.5, 3);
			std::cout << "points in the slice plane: " << z_init + 0.5 + i * del_z << ", " << zslice_points.size() << endl;
			zslice_point_count.push_back(zslice_points.size());
			zslice_value.push_back(z_init + 0.5 + i * del_z);
			zslice_points.resize(0);
		}
		std::cout << "z_now: ";
		cin >> z_now;
		/*zslice_max = find_1d_vector_max(zslice_point_count, 1);
		z_now = zslice_value[zslice_max];*/
		std::cout << "z value:" << z_now << endl;

		find_point_3d_1feature(zslice_points, block_a.vertices, z_now, 0.5, 3);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}
	else
		return -1;
	//std::cout << "y_min_value: " << block_a.vertices[y_min_i].y << endl;
	/*target_plot_update(block_a.vertices);
	target_plot_color_update(0.0, 0.8, 0.0, block_a.vertices.size());*/
}
int cad3d::adaptive_2p5d_decomposition(vector<Point3D> &v_in, vector<Point3D> &zslice_points, float &z_now, int dim, int n_total, float z_tol) {
	//returns 1 if feature plane found or returns 0 or -1;

	vector<float> zslice_point_count;
	vector<float> zslice_value;
	int zslice_max;
	float z_temp;
	
	if (dim == 3) {
		int z_min_i = find_point3d_min(v_in, 3);
		int z_max_i = find_point3d_max(v_in, 3);

		float z_init = v_in[z_min_i].z;
		//float del_z = 1;
		float del_z = (v_in[z_max_i].z- v_in[z_min_i].z)/n_total;
		//should be minimum to maximum points
		for (int i = 0; i < n_total; i++) {
			if (z_init  + i * del_z > v_in[z_max_i].z)
				break;
			find_point_3d_1feature(zslice_points, v_in, z_init+ i * del_z, 0.5, 3);
			if (zslice_points.size()>0) {
				std::cout << "points in the slice plane: " << z_init + (i * del_z) << ", " << zslice_points.size() << endl;
				z_temp = z_init + (i * del_z);
				fine_slice_calculation(v_in,z_temp,3, del_z);
			}
			zslice_point_count.push_back(zslice_points.size());
			zslice_value.push_back(z_init + i * del_z);
			zslice_points.resize(0);
		}
		std::cout << "z_now: ";
		cin >> z_now;
		/*zslice_max = find_1d_vector_max(zslice_point_count, 1);
		z_now = zslice_value[zslice_max];*/
		std::cout << "z value:" << z_now << endl;

		find_point_3d_1feature(zslice_points, v_in, z_now, 0.1, 3);
		std::cout << "points in the slice plane: " << zslice_points.size() << endl;
		if (zslice_points.size() < 1) {
			std::cout << "No points centering this value\n";
			return -1;
		}
		else
			return 1;
	}
	else
		return -1;
	//std::cout << "y_min_value: " << block_a.vertices[y_min_i].y << endl;
	/*target_plot_update(block_a.vertices);
	target_plot_color_update(0.0, 0.8, 0.0, block_a.vertices.size());*/
}
void cad3d::fine_slice_calculation(vector<Point3D> &v_in, float &z_now, int dim, float z_tol) {
	float del_z = 2*z_tol/10;
	float z_init = z_now - z_tol;
	vector<Point3D> zslice_points;
	//should be minimum to maximum points
	for (int i = 0; i < 10; i++) {
		//if (z_init + 0.1 + i * del_z > v_in[z_max_i].z)
			//break;
		find_point_3d_1feature(zslice_points, v_in, z_init + i * del_z, del_z, 3);
		std::cout << "points in the fine slice plane: " << z_init + (i * del_z) << ", " << zslice_points.size() << endl;
		
		zslice_points.resize(0);
	}
}

void cad3d::cartesean_processing_extended(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2,int dim, float z_now) {

	if (dim == 1)  {
		getplane_simplfied_block_extended(block_a, slot_plane, tslot_i,1,z_now);
		triangle_normal_update(block_a.triangles, tslot_i);
		int ii;
		bool a, b;

		vector<cad3d::Point3D> v_temp, v_line1, v_line2;
		vector<int> i_line;
		triangle2line(block_a.vertices, tslot_i, v_line1, v_line2, i_line, v_temp);
		//newblock.triangle2line(block_a.vertices, tbi_1, v_line1, v_line2, i_line, v_temp);
		std::cout << "triangle to line completed " << endl;
		std::cout << "tslot_i size: " << tslot_i.size() << endl;
		std::cout << "v_line size:/3 " << v_line1.size() / 3 << endl;
		cad3d::Point3D line_value;

		vector<cad3d::Point3D>  vord1_line1, vord1_line2;
		cad3d::Point3D ordt1, ordt2;
		float delz1, dely1;
		//ordering lines
		//for rectangular block
		for (int i = 0; i < v_line1.size(); i++) {
			dely1 = v_line2[i].y - v_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
			else {
				delz1 = v_line2[i].z - v_line1[i].z;
				if (delz1 > 0) {
					ordt2 = v_line2[i];
					ordt1 = v_line1[i];
				}
				else if (delz1 < 0) {
					ordt2 = v_line1[i];
					ordt1 = v_line2[i];
				}
			}
			vord1_line1.push_back(ordt1);
			vord1_line2.push_back(ordt2);

		}

		for (int i = 0; i < vord1_line1.size(); i++) {
			//for (int i = 0; i < tslot_i.size(); i++) {
			//for (int i = 0; i < 100; i++) {
			line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
			a = line_value.x < 20.0;
			b = line_value.y < 20.0;
			//b = 0;
			if (a || b) {
				if (sqrt(line_value.x*line_value.x + line_value.y*line_value.y + line_value.z*line_value.z) < 100) {
					v_ext_l1.push_back(vord1_line1[i]);
					v_ext_l2.push_back(vord1_line2[i]);
				}

				//std::cout << sqrt(line_value.x*line_value.x + line_value.y*line_value.y) << endl;
			}
			if (line_value.x == 0.0) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}
			if (line_value.y == 0.0) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}

		}
	}
	//feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);
	if (dim == 3) {
		//getplane_simplfied_block_extended(block_a, slot_plane, tslot_i, 3, z_now);
		//triangle_normal_update(block_a.triangles, tslot_i);
		getplane_simplfied_block(block_a, slot_plane, tslot_i);
		triangle_normal_update(block_a.triangles, tslot_i);
		int ii;
		bool a, b;

		vector<cad3d::Point3D> v_temp, v_line1, v_line2;
		vector<int> i_line;
		triangle2line(block_a.vertices, tslot_i, v_line1, v_line2, i_line, v_temp);
		//newblock.triangle2line(block_a.vertices, tbi_1, v_line1, v_line2, i_line, v_temp);
		std::cout << "triangle to line completed " << endl;
		std::cout << "tslot_i size: " << tslot_i.size() << endl;
		std::cout << "v_line size:/3 " << v_line1.size() / 3 << endl;
	}
}
void cad3d::getplane_simplfied_block_extended(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i,int dim, float z_now) {
	if (dim ==1) {
		/*int x_max_i = find_point3d_max(block_a.vertices, 1);
		int y_max_i = find_point3d_max(block_a.vertices, 2);
		int z_max_i = find_point3d_max(block_a.vertices, 3);
		int x_min_i = find_point3d_min(block_a.vertices, 1);
		int y_min_i = find_point3d_min(block_a.vertices, 2);
		int z_min_i = find_point3d_min(block_a.vertices, 3);
		std::cout << "x_max_value: " << block_a.vertices[z_now].x << endl;*/


		vector<int> slot_plane_i;
		//centered_block = block_a.vertices;
		find_point_3d_1feature(slot_plane, block_a.vertices, z_now, 0.5, dim);
		slot_plane_i.resize(block_a.vertices.size());
		find_point3d_i_1feature(slot_plane_i, block_a.vertices, z_now, 0.5, 3);
		show_index_frequency(slot_plane_i);

		vector<cad3d::triangle_index> t_slot, ti_all, tslot_i1, tslot_i2, tslot_i3;
		triangle_vertexindex_init(block_a.triangles, ti_all);
		vector<cad3d::triangle> tslot;
		std::cout << "total triangle index " << ti_all.size() << endl;
		//triangle and line coordination
		float zmid = (z_now-0.5 + z_now +0.5) / 2;
		triangle_allpoint_1feature(block_a.triangles, tslot, z_now, 0.5, dim);
		std::cout << "tslot size" << tslot.size() << endl;
		triangle_allpoint_i1feature(block_a.vertices, ti_all, tslot_i, z_now, 0.5, dim);
	}
}

/*void cad3d::return_unit_vector(Vector3d &v) {
	float v_mod = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v = v / v_mod;
}*/
cad3d::Point3D cad3d::find_cross_ray(Point3D &p1, Point3D &p2, Point3D &p3) {
	//p2 is th center point along which the cross is generated
	Point3D cross_pt;
	//https://eigen.tuxfamily.org/dox/group__TutorialMatrixArithmetic.html

	float dx, dy, dz;
	dx = p2.x - p1.x; dy = p2.y - p1.y; dz = p2.z - p1.z;
	Vector3d v_r1(dx, dy, dz);

	dx = p3.x - p2.x;	dy = p3.y - p2.y;	dz = p3.z - p2.z;
	Vector3d v_r2(dx, dy, dz);

	Vector3d n0 = v_r1.cross(v_r2);
	//return_unit_vector(n0);
	/*float n_mod = sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
	n0 = n0 / n_mod;
	std::cout << "n0 (unit vector): (" << n0[0] << ", ";
	std::cout << n0[1]  << ", "<< n0[2]  << ")"<<endl;*/

	cross_pt.x = n0[0];
	cross_pt.y = n0[1];
	cross_pt.z = n0[2];
	return cross_pt;
	/*Vector3d v_cross1 = v_r1.cross(n0);
	return_unit_vector(v_cross1);
	Vector3d v_cross2 = v_r2.cross(n0);
	return_unit_vector(v_cross2);
	std::cout << "v_cross1 (unit vector): (" << v_cross1[0] << ", ";
	std::cout << v_cross1[1] << ", " << v_cross1[2] << ")" << endl;
	std::cout << "v_cross2 (unit vector): (" << v_cross2[0] << ", ";
	std::cout << v_cross2[1] << ", " << v_cross2[2] << ")" << endl;*/
}
cad3d::Point3D cad3d::return_unitv_pt3d_mode(Point3D p1) {
	Point3D p_out;
	float v_mod = sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
	p_out.x = p1.x / v_mod;
	p_out.y = p1.y / v_mod;
	p_out.z = p1.z / v_mod;

	return p_out;
}
void cad3d::calculate_cross_rays_v1(Point3D p1, Point3D p2, Point3D p3, Point3D &n0, Point3D &n1, Point3D &n2) {
	// for a single segment
	//normal vector perpendicular t the slot plane
	n0 = find_cross_ray(p1, p2, p3);

	Point3D n_temp;
	n_temp.x = p2.x + n0.x; n_temp.y = p2.y + n0.y; n_temp.z = p2.z + n0.z;

	n1 = find_cross_ray(p1, p2, n_temp);
	//ncross_temp1.x = p2.x + n1.x; ncross_temp1.y = p2.y + n1.y; ncross_temp1.z = p2.z + n1.z;

	n2 = find_cross_ray(p3, p2, n_temp);
}
void cad3d::calculate_cross_rays_v2(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, vector<Point3D> &n1_ray1, vector<Point3D> &n1_ray2, vector<Point3D> &n2_ray1, vector<Point3D> &n2_ray2) {
	// for a multiple segments
	//need to add one more start point and one more end pointca

	//Point3D p1, Point3D p2, Point3D p3, Point3D &n0, Point3D &n1, Point3D &n2
	cad3d::Point3D n0, n1, n2, n_temp;
	//vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;
	cad3d::Point3D ncross_temp1, ncross_temp2, ncross1_dir, ncross2_dir;

	for (int i = 1; i < slot_ip1.size() - 2; i++) {
		//for (int i = 1; i < 6 - 1; i++) {

		calculate_cross_rays_v1(slot_ip1[i - 1], slot_ip1[i], slot_ip1[i + 1], n0, n1, n2);

		ncross1_dir = test_n_correct_slope_direction(slot_ip1[4], n1, slot_e1, slot_e3);
		ncross2_dir = test_n_correct_slope_direction(slot_ip1[4], n2, slot_e1, slot_e3);

		n1_ray1.push_back(slot_ip1[i]);
		ncross_temp1.x = slot_ip1[i].x + n1.x;
		ncross_temp1.y = slot_ip1[i].y + n1.y*ncross1_dir.y;
		ncross_temp1.z = slot_ip1[i].z + n1.z;
		n1_ray2.push_back(ncross_temp1);

		n2_ray1.push_back(slot_ip1[i]);
		ncross_temp2.x = slot_ip1[i].x + n2.x;
		ncross_temp2.y = slot_ip1[i].y + n2.y*ncross2_dir.y;
		ncross_temp2.z = slot_ip1[i].z + n2.z;
		n2_ray2.push_back(ncross_temp2);
	}
}
void cad3d::calculate_cross_rays_v3(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, vector<Point3D> &n1_ray1, vector<Point3D> &n1_ray2, vector<Point3D> &n2_ray1, vector<Point3D> &n2_ray2) {
	// for a multiple segments
	//need to add one more start point and one more end pointca

	//Point3D p1, Point3D p2, Point3D p3, Point3D &n0, Point3D &n1, Point3D &n2
	cad3d::Point3D n0, n1, n2, n1_0, n2_0, n_temp;
	//vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;
	cad3d::Point3D ncross_temp1, ncross_temp2, ncross1_dir, ncross2_dir;

	for (int i = 1; i < slot_ip1.size() - 2; i++) {
		//for (int i = 1; i < 6 - 1; i++) {

		calculate_cross_rays_v1(slot_ip1[i - 1], slot_ip1[i], slot_ip1[i + 1], n0, n1_0, n2_0);

		n1 = return_unitv_pt3d_mode(n1_0);
		n2 = return_unitv_pt3d_mode(n2_0);

		/*if (abs(n1.z) > 0) {
			n1.z = 0.0;
			n1 = return_unitv_pt3d_mode(n1);
		}

		if (abs(n2.z) > 0) {
			n2.z = 0.0;
			n2 = return_unitv_pt3d_mode(n2);
		}*/

		std::cout << i << "- n1: " << n1.x << ", " << n1.y << ", " << n1.z << " ~ ";
		std::cout << i << "- n2: " << n2.x << ", " << n2.y << ", " << n2.z << endl;

		ncross1_dir = test_n_correct_slope_direction(slot_ip1[4], n1, slot_e1, slot_e3);
		ncross2_dir = test_n_correct_slope_direction(slot_ip1[4], n2, slot_e1, slot_e3);
		/*https://stackoverflow.com/questions/1800138/given-a-start-and-end-point-and-a-distance-calculate-a-point-along-a-line*/
		float orig2dist = 40.0;

		n1_ray1.push_back(slot_ip1[i]);
		ncross_temp1.x = slot_ip1[i].x + orig2dist * n1.x;
		ncross_temp1.y = slot_ip1[i].y + orig2dist * n1.y*ncross1_dir.y;
		ncross_temp1.z = slot_ip1[i].z + orig2dist * n1.z;
		n1_ray2.push_back(ncross_temp1);

		n2_ray1.push_back(slot_ip1[i]);
		ncross_temp2.x = slot_ip1[i].x + orig2dist * n2.x;
		ncross_temp2.y = slot_ip1[i].y + orig2dist * n2.y*ncross2_dir.y;
		ncross_temp2.z = slot_ip1[i].z + orig2dist * n2.z;
		n2_ray2.push_back(ncross_temp2);
	}
}
void cad3d::calculate_cross_rays_v4(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, vector<Point3D> &n1_ray1, vector<Point3D> &n1_ray2, vector<Point3D> &n2_ray1, vector<Point3D> &n2_ray2, int i) {
	// for a testing single segment
	//need to add one more start point and one more end pointca

	//Point3D p1, Point3D p2, Point3D p3, Point3D &n0, Point3D &n1, Point3D &n2
	cad3d::Point3D n0, n1, n2, n1_0, n2_0, n_temp;
	//vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;
	cad3d::Point3D ncross_temp1, ncross_temp2, ncross1_dir, ncross2_dir;

	//int i = 3;

	//for (int i = 1; i < slot_ip1.size() - 2; i++) {
		//for (int i = 1; i < 6 - 1; i++) {

	calculate_cross_rays_v1(slot_ip1[i - 1], slot_ip1[i], slot_ip1[i + 1], n0, n1_0, n2_0);

	n1 = return_unitv_pt3d_mode(n1_0);
	n2 = return_unitv_pt3d_mode(n2_0);

	if ((n1.z == 0) && (n2.z == 0)) {
		//std::cout << i << "- n1: " << n1.x << ", " << n1.y << ", " << n1.z << " ~ ";
		//std::cout << i << "- n2: " << n2.x << ", " << n2.y << ", " << n2.z << endl;

		ncross1_dir = test_n_correct_slope_direction(slot_ip1[4], n1, slot_e1, slot_e3);
		ncross2_dir = test_n_correct_slope_direction(slot_ip1[4], n2, slot_e1, slot_e3);
		/*https://stackoverflow.com/questions/1800138/given-a-start-and-end-point-and-a-distance-calculate-a-point-along-a-line*/
		float orig2dist = 40.0;

		n1_ray1.push_back(slot_ip1[i]);
		ncross_temp1.x = slot_ip1[i].x + orig2dist * n1.x;
		ncross_temp1.y = slot_ip1[i].y + orig2dist * n1.y*ncross1_dir.y;
		ncross_temp1.z = slot_ip1[i].z + orig2dist * n1.z;
		n1_ray2.push_back(ncross_temp1);

		n2_ray1.push_back(slot_ip1[i]);
		ncross_temp2.x = slot_ip1[i].x + orig2dist * n2.x;
		ncross_temp2.y = slot_ip1[i].y + orig2dist * n2.y*ncross2_dir.y;
		ncross_temp2.z = slot_ip1[i].z + orig2dist * n2.z;
		n2_ray2.push_back(ncross_temp2);
	}
}
cad3d::Point3D cad3d::test_n_correct_slope_direction(Point3D p0, Point3D n1, vector<Point3D> &slot_e1, vector<Point3D> &slot_e3) {
	Point3D p_dir;//return direction values in signed 0,1 
	//and to be multiplied with the input normal in point3D format
	Point3D n_start; //open edge direction
	n_start.x = slot_e3[0].x - slot_e1[0].x;
	n_start.y = slot_e3[0].y - slot_e1[0].y;
	n_start.z = slot_e3[0].z - slot_e1[0].z;

	if (n1.x*n_start.x > 0)
		p_dir.x = 1;
	else if (n1.x*n_start.x < 0)
		p_dir.x = -1;
	else
		p_dir.x = 0;

	if (n1.y*n_start.y > 0)
		p_dir.y = 1;
	else if (n1.y*n_start.y < 0)
		p_dir.y = -1;
	else
		p_dir.y = 0;

	if (n1.z*n_start.z > 0)
		p_dir.z = 1;
	else if (n1.z*n_start.z < 0)
		p_dir.z = -1;
	else
		p_dir.z = 0;

	return p_dir;
}

float cad3d::DistancePointLine_modv1(Point3D p0, Point3D &LineStart, Point3D &LineEnd, float &dist)
{
	//compares the unit line vectors
	//original from http://paulbourke.net/geometry/pointlineplane/source.c
	float LineMag;
	float U;

	Point3D new_line, new_line_norm;
	new_line.x = p0.x - LineStart.x;
	new_line.y = p0.y - LineStart.y;
	new_line.z = p0.z - LineStart.z;
	new_line_norm = return_unitv_pt3d_mode(new_line);

	Point3D line_org, line_org_norm;
	line_org.x = LineEnd.x - LineStart.x;
	line_org.y = LineEnd.y - LineStart.y;
	line_org.z = LineEnd.z - LineStart.z;
	line_org_norm = return_unitv_pt3d_mode(line_org);
	//std::cout << new_line_norm.x <<", "<<line_org_norm.x << endl;
	//std::cout << dist201(new_line_norm, line_org_norm) << endl;
	float dd = dist201(new_line_norm, line_org_norm);
	return dd;
}

void cad3d::polar_processing(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, float z_now, vector<Point3D> &centered_block) {
	//output is slot_plane, tslot_i, v_ext_l1, v_ext_l2, feat_e1, feat_e2
	//vector<cad3d::Point3D> centered_block;
	vector<int> slot_plane_i;
	centre_at_x0y0(block_a.vertices, centered_block);
	find_point_3d_1feature(slot_plane, centered_block, z_now, 0.5, 3);
	slot_plane_i.resize(centered_block.size());
	find_point3d_i_1feature(slot_plane_i, centered_block, z_now, 0.5, 3);
	show_index_frequency(slot_plane_i);

	vector<cad3d::triangle_index> ti_all;//, tslot_i, tslot_i1, tslot_i2, tslot_i3;
	vector<cad3d::triangle> tslot;
	triangle_normal_update(block_a.triangles, tslot_i);
	triangle_vertexindex_init(block_a.triangles, ti_all);

	int ii;
	bool a, b;
	triangle_vertexindex_init(block_a.triangles, ti_all);

	triangle_allpoint_1feature(block_a.triangles, tslot, z_now, 0.5, 3);
	std::cout << "tslot size" << tslot.size() << endl;
	triangle_allpoint_i1feature(centered_block, ti_all, tslot_i, z_now, 0.5, 3);
	std::cout << "tslot i size" << tslot_i.size() << endl;

	vector<cad3d::Point3D> v_temp, v_line1, v_line2;
	vector<int> i_line;
	triangle2line(block_a.vertices, tslot_i, v_line1, v_line2, i_line, v_temp);
	std::cout << "triangle to line completed " << endl;
	std::cout << "v_line size/3: " << v_line1.size() / 3 << endl;

	//ordering lines
	//for polar block
	vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	//for polar block
	vector<cad3d::Point3D> vpol_line1, vpol_line2;
	convert_cart2pol(v_line1, vpol_line1);
	convert_cart2pol(v_line2, vpol_line2);
	//iner radis to puter radius
	//lpwer angle to higher angle
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = vpol_line2[i].x - vpol_line1[i].x;
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = vpol_line2[i].y - vpol_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}
	cad3d::Point3D line_value;
	for (int i = 0; i < vord1_line1.size(); i++) {
		//for (int i = 0; i < tslot_i.size(); i++) {
		//for (int i = 0; i < 100; i++) {
		line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
		a = line_value.x < 25.0;
		b = line_value.y < 25.0;
		//b = 0;
		if (a || b) {
			if (sqrt(line_value.x*line_value.x + line_value.y*line_value.y) < 40) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}

			//std::cout << sqrt(line_value.x*line_value.x + line_value.y*line_value.y) << endl;
		}
		if (line_value.x == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
		if (line_value.y == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
	}
	feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);
}
void cad3d::one_start_point_polar_block_v1(stl_data &block_a, vector<Point3D> &centered_block, Point3D pref_start, vector<Point3D> &vcorner, vector<Point3D> &red_vcorner, vector<Point3D> &vend_corner, float z_now) {
	//vector<cad3d::Point3D> slot_line1, slot_line2, v_line1, v_line2;
	vector<int> vi_all, vi_temp;//tracking use of all points
	vector<int> vi_temp551, vi_temp552;
	vector<int> used_i_v_external, seg_ui2, seg_ui3, vi_ext;//, vl1_ui2, vl2_ui2; // index reference to vline1 and 2
	vector<cad3d::Point3D> vv1, vv2, v2tl1, v2tl2, v2poltl1, v2poltl2, v3ext, v3tl1, v3tl2;
	vi_ext.resize(centered_block.size());
	initialize_0i(vi_ext);
	vector<cad3d::Point3D> v_temp, v_line1, v_line2;
	vector<int>i_line, i_line_nr, i_line_r, i_line_frequency;

	vector<triangle_index> tslot_i, ti_all;
	vector<cad3d::triangle> tslot;
	triangle_vertexindex_init(block_a.triangles, ti_all);
	triangle_allpoint_1feature(block_a.triangles, tslot, z_now, 0.5, 3);
	//std::cout << "tslot size" << tslot.size() << endl;
	triangle_allpoint_i1feature(centered_block, ti_all, tslot_i, z_now, 0.5, 3);
	//std::cout << "tslot i size" << tslot_i.size() << endl;

	triangle2line(centered_block, tslot_i, v_line1, v_line2, i_line, v_temp);
	//std::cout << "triangle to line completed" << endl;
	//std::cout << "line array size: " << v_line2.size() << endl;

	vi_ext.resize(centered_block.size());
	initialize_0i(vi_ext);
	find_external_segmentsv4(v_temp, v_line1, v_line2, vv1, vv2, used_i_v_external, tslot_i, vi_ext);
	//newblock.find_external_segmentsv4(v_temp, v_ext_l1, v_ext_l2, vv1, vv2, used_i_v_external, tslot_i, vi_ext);
	//std::cout << "external points size" << vv1.size() << endl;
	vi_all.resize(centered_block.size());
	initialize_0i(vi_all);
	vi_temp.resize(centered_block.size());
	initialize_0i(vi_temp);
	get_segment_numberv505(vv1, vv2, v2tl1, v2tl2, used_i_v_external, vi_all, centered_block, v3ext, vcorner);

	vector<cad3d::Point3D> slot_line1, slot_line2, v1, vinward1;
	int vi_temp1s1, vi_temp1s2, vi_temp2, vi_temp11;
	//vector<int> vi_temp551, vi_temp552;

	relate_corner_w_line(vi_all, centered_block, v_line1, v_line2, vcorner, red_vcorner);
	get_end_corner_v101(v2tl1, v2tl2, vend_corner, red_vcorner);
	//std::cout << "end corner size: " << vend_corner.size() << endl;
	//newblock.get_internal_slot_v1031(vi_all, centered_block, slot_line1, slot_line2, v_line1, v_line2, v1, vi_temp1s1);

	//std::cout << red_vcorner.size();
}
void cad3d::one_start_point_polar_block_v1_phase2(vector<Point3D> &centered_block, Point3D pref_start, vector<int> &vi_all, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<int> &vi_temp551, vector<int> &vi_temp552, vector<Point3D> &v_out) {
	/* vector<Point3D> &v_line1, vector<Point3D> &v_line2, */
	vector<int> jd1, jd2;
	std::cout << "pref_start: " << pref_start.x << ", ";
	std::cout << pref_start.y << endl;
	find_Point3D_index(pref_start, centered_block, vi_all, vi_temp551);
	float d1, d2;
	vector<float> jdmat1, jdmat2;
	for (int j = 0; j < v_line1.size(); j++) {
		d1 = dist201(pref_start, v_line1[j]);//v_ext_l2[0]
		d2 = dist201(pref_start, v_line2[j]);
		//if (d1 == 0) {
		if (d1 < 0.1) {
			std::cout << "d1: " << j << "- " << dist201(pref_start, v_line1[j]) << endl;
			jd1.push_back(j);
			jdmat1.push_back(dist201(pref_start, v_line2[j]));
		}
		//if (d2 == 0) {
		if (d2 < 0.1) {
			//find the points in which they are connected
			std::cout << "d2: " << j << "- " << dist201(pref_start, v_line2[j]) << endl;
			jd2.push_back(j);
			jdmat2.push_back(dist201(pref_start, v_line1[j]));
		}
	}
	std::cout << jd1.size() << ", " << jdmat1.size() << endl;
	std::cout << jd2.size() << ", " << jdmat2.size() << endl;
	//corner detection, not optimum for single basic block
	if (jd2.size() > 1) {
		if (jdmat2[1] < jdmat2[0]) {
			find_Point3D_indexv3(v_line1[jd2[1]], centered_block, vi_all, vi_temp552);
		}
		else
		{
			find_Point3D_indexv3(v_line1[jd2[0]], centered_block, vi_all, vi_temp552);
		}



	}
	if (jd1.size() == 1 & jd2.size() == 1) {
		if (jdmat1[0] < jdmat2[0]) {
			find_Point3D_indexv3(v_line2[jd1[0]], centered_block, vi_all, vi_temp552);
		}
		else
		{
			find_Point3D_indexv3(v_line1[jd2[0]], centered_block, vi_all, vi_temp552);
		}
	}
	else {
		if (jd2.size() == 1) {

			find_Point3D_indexv3(v_line1[jd2[0]], centered_block, vi_all, vi_temp552);

		}
		if (jd1.size() == 1) {

			find_Point3D_indexv3(v_line2[jd1[0]], centered_block, vi_all, vi_temp552);

		}
	}
	std::cout << vi_temp552.size();
}
void cad3d::one_start_point_polar_block_v1_phase3(vector<Point3D> &centered_block, Point3D pref_start, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v_out) {
	/* vector<Point3D> &v_line1, vector<Point3D> &v_line2, */
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	for (int j = 0; j < v_line1.size(); j++) {
		d1 = dist201(pref_start, v_line1[j]);//v_ext_l2[0]
		d2 = dist201(pref_start, v_line2[j]);
		Point3D p1, p2, p2_con;
		if (d1 < 0.01) {
			convert_cart2polPoint3D(v_line1[j], p1);
			//std::cout << "d1: "<< j <<", r: "<<p1.x << endl;
		}
		if (d2 < 0.01) {
			jr2.push_back(j);
			convert_cart2polPoint3D(v_line1[j], p2_con);
			std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
			jrmat2.push_back(p2_con.x);
		}
	}
	/*for (int i = 0; i < jr2.size(); i++) {
	}*/
	int jd2min = find_1d_vector_min(jrmat2, 1);
	v_out.push_back(v_line1[jr2[jd2min]]);
}
void cad3d::one_start_point_polar_block_v1_phase4(vector<Point3D> &centered_block, Point3D pref_start, vector<Point3D> &v_line1, vector<Point3D> &v_line2, vector<Point3D> &v_out) {
	/* vector<Point3D> &v_line1, vector<Point3D> &v_line2, */
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	for (int j = 0; j < v_line1.size(); j++) {
		d1 = dist201(pref_start, v_line1[j]);//v_ext_l2[0]
		d2 = dist201(pref_start, v_line2[j]);
		Point3D p1_con, p2, p2_con;
		if (d1 < 0.01) {
			jr1.push_back(j);
			convert_cart2polPoint3D(v_line2[j], p1_con);
			//std::cout << "d1: " << j << ", r: " << p1_con.x << endl;
			jrmat1.push_back(p1_con.x);
		}
		if (d2 < 0.01) {
			jr2.push_back(j);
			convert_cart2polPoint3D(v_line1[j], p2_con);
			//std::cout << "d2: " << j << ", r: " << p2_con.x << endl;
			jrmat2.push_back(p2_con.x);
		}
	}
	/*for (int i = 0; i < jr2.size(); i++) {
	}*/
	/*int */
	int jdmin;
	if (jrmat1.size() > 0) {
		std::cout << "inward edge \n";
		jdmin = find_1d_vector_min(jrmat1, 1);
		v_out.push_back(v_line2[jr1[jdmin]]);
	}
	else if (jrmat2.size() > 0) {
		std::cout << "outward edge \n";
		jdmin = find_1d_vector_min(jrmat2, 1);
		v_out.push_back(v_line1[jr2[jdmin]]);
	}
}

void cad3d::find_feature_edge_break_points_v1(stl_data &block_a, vector<int> &vi_all, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, std::vector<Point3D>  &pref_start, vector<int> &vi_temp551, vector<int> &vi_temp552) {
	//vector<cad3d::Point3D> slot_line1, slot_line2, v1, vinward1;
	bool a1, b1, c1, d1, e, f;
	bool a2, b2, c2, d2;

	//corner point estmation 
	int ixmax, iymax, ixmin, iymin;
	ixmax = find_point3d_max(block_a.vertices, 1);
	ixmin = find_point3d_min(block_a.vertices, 1);
	iymax = find_point3d_max(block_a.vertices, 2);
	iymin = find_point3d_min(block_a.vertices, 2);

	a1 = 0; b1 = 0; c1 = 0; d1 = 0;
	a2 = 0; b2 = 0; c2 = 0; d2 = 0;
	e = 0; f = 0;

	for (int i = 0; i < v_line1.size(); i++) {
		a1 = abs(v_line1[i].x - block_a.vertices[ixmax].x) < 0.1;
		b1 = abs(v_line1[i].x - block_a.vertices[ixmin].x) < 0.1;
		c1 = abs(v_line1[i].y - block_a.vertices[iymax].y) < 0.1;
		d1 = abs(v_line1[i].y - block_a.vertices[iymin].y) < 0.1;

		e = (a1 || b1) || (c1 || d1);

		a2 = abs(v_line2[i].x - block_a.vertices[ixmax].x) < 0.1;
		b2 = abs(v_line2[i].x - block_a.vertices[ixmin].x) < 0.1;
		c2 = abs(v_line2[i].y - block_a.vertices[iymax].y) < 0.1;
		d2 = abs(v_line2[i].y - block_a.vertices[iymin].y) < 0.1;

		f = (a2 || b2) || (c2 || d2);
		/*if (e == 1)
			std::cout << i << "e" << endl;
		if (f == 1)
			std::cout << i << "f" << endl;*/
		if (e*(!f) == 1) {
			//std::cout << i << "e.!f" << endl;
			pref_start.push_back(v_line1[i]);
		}
			
		if (f*(!e) == 1) {
			//std::cout << i << "!e.f" << endl;
			pref_start.push_back(v_line1[i]);
		}
			
		a1 = 0; b1 = 0; c1 = 0; d1 = 0;
		a2 = 0; b2 = 0; c2 = 0; d2 = 0;
		e = 0; f = 0;
	}

	
	//open slot corner line inspect
}
void cad3d::find_feature_edge_break_points_v2(stl_data &block_a, vector<int> &vi_all, vector<Point3D> &slot_line1, vector<Point3D> &slot_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, std::vector<Point3D>  &pref_start, vector<int> &vi_temp551, vector<int> &vi_temp552) {
	//vector<cad3d::Point3D> slot_line1, slot_line2, v1, vinward1;
	bool a1, b1, c1, d1, e, f;
	bool a2, b2, c2, d2;

	//corner point estmation 
	int ixmax, iymax, ixmin, iymin;
	ixmax = find_point3d_max(block_a.vertices, 1);
	ixmin = find_point3d_min(block_a.vertices, 1);
	iymax = find_point3d_max(block_a.vertices, 2);
	iymin = find_point3d_min(block_a.vertices, 2);

	a1 = 0; b1 = 0; c1 = 0; d1 = 0;
	a2 = 0; b2 = 0; c2 = 0; d2 = 0;
	e = 0; f = 0;

	for (int i = 0; i < v_line1.size(); i++) {
		a1 = abs(v_line1[i].x - block_a.vertices[ixmax].x) < 0.1;
		b1 = abs(v_line1[i].x - block_a.vertices[ixmin].x) < 0.1;
		c1 = abs(v_line1[i].y - block_a.vertices[iymax].y) < 0.1;
		d1 = abs(v_line1[i].y - block_a.vertices[iymin].y) < 0.1;

		e = (a1 || b1) || (c1 || d1);

		a2 = abs(v_line2[i].x - block_a.vertices[ixmax].x) < 0.1;
		b2 = abs(v_line2[i].x - block_a.vertices[ixmin].x) < 0.1;
		c2 = abs(v_line2[i].y - block_a.vertices[iymax].y) < 0.1;
		d2 = abs(v_line2[i].y - block_a.vertices[iymin].y) < 0.1;

		f = (a2 || b2) || (c2 || d2);
		/*if (e == 1)
			std::cout << i << "e" << endl;
		if (f == 1)
			std::cout << i << "f" << endl;*/
		if (e*(!f) == 1) {
			std::cout << i << "e.!f" << endl;
			pref_start.push_back(v_line1[i]);
		}

		if (f*(!e) == 1) {
			std::cout << i << "!e.f" << endl;
			pref_start.push_back(v_line1[i]);
		}

		a1 = 0; b1 = 0; c1 = 0; d1 = 0;
		a2 = 0; b2 = 0; c2 = 0; d2 = 0;
		e = 0; f = 0;
	}


	//open slot corner line inspect
}
void cad3d::find_cartesian_external_points_v1(stl_data &block_a, vector<int> &vi_all, vector<Point3D> &b_line1, vector<Point3D> &b_line2, vector<Point3D> &v_line1, vector<Point3D> &v_line2, std::vector<Point3D>  &pref_start, vector<int> &vi_temp551, vector<int> &vi_temp552) {
	//vector<cad3d::Point3D> slot_line1, slot_line2, v1, vinward1;
	bool a1, b1, c1, d1, e, f;
	bool a2, b2, c2, d2;

	//corner point estmation 
	int ixmax, iymax, ixmin, iymin;
	ixmax = find_point3d_max(block_a.vertices, 1);
	ixmin = find_point3d_min(block_a.vertices, 1);
	iymax = find_point3d_max(block_a.vertices, 2);
	iymin = find_point3d_min(block_a.vertices, 2);

	a1 = 0; b1 = 0; c1 = 0; d1 = 0;
	a2 = 0; b2 = 0; c2 = 0; d2 = 0;
	e = 0; f = 0;

	for (int i = 0; i < v_line1.size(); i++) {
		a1 = abs(v_line1[i].x - block_a.vertices[ixmax].x) < 0.1;
		b1 = abs(v_line1[i].x - block_a.vertices[ixmin].x) < 0.1;
		c1 = abs(v_line1[i].y - block_a.vertices[iymax].y) < 0.1;
		d1 = abs(v_line1[i].y - block_a.vertices[iymin].y) < 0.1;

		e = (a1 || b1) || (c1 || d1);

		a2 = abs(v_line2[i].x - block_a.vertices[ixmax].x) < 0.1;
		b2 = abs(v_line2[i].x - block_a.vertices[ixmin].x) < 0.1;
		c2 = abs(v_line2[i].y - block_a.vertices[iymax].y) < 0.1;
		d2 = abs(v_line2[i].y - block_a.vertices[iymin].y) < 0.1;

		f = (a2 || b2) || (c2 || d2);
		/*if (e == 1)
			std::cout << i << "e" << endl;
		if (f == 1)
			std::cout << i << "f" << endl;*/
		if (e*f == 1) {
			//std::cout << i << "e.f" << endl;
			b_line1.push_back(v_line1[i]);
			b_line2.push_back(v_line2[i]);
		}

		/*if (f*(!e) == 1) {
			std::cout << i << "!e.f" << endl;
			pref_start.push_back(v_line1[i]);
		}*/

		a1 = 0; b1 = 0; c1 = 0; d1 = 0;
		a2 = 0; b2 = 0; c2 = 0; d2 = 0;
		e = 0; f = 0;
	}


	//open slot corner line inspect
}
void cad3d::crop_slot_single_edge(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_e01, vector<Point3D> &slot_e02) {
	cad3d::circle2d c_temp;
	slot_e1.push_back(slot_e01[0]); slot_e2.push_back(slot_e02[0]);
	slot_e1.push_back(slot_e01[1]); slot_e2.push_back(slot_e02[1]);
	for (int i = 2; i < slot_e01.size(); i++) {
		c_temp = compute_circle_4_point3d(slot_e01[i], slot_e01[i - 1], slot_e01[i - 2]);
		//std::cout << i << ": " << c_temp.r << endl;
		if (abs(c_temp.r) > 50) {
			slot_e1.push_back(slot_e01[i]);
			slot_e2.push_back(slot_e02[i]);
		}
		else if (abs(c_temp.r) <= 100) {
			break;
		}


	}
}
void cad3d::crop_slot_single_edge_v2(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_e01, vector<Point3D> &slot_e02) {
	cad3d::circle2d c_temp;
	slot_e1.push_back(slot_e01[0]); slot_e2.push_back(slot_e02[0]);
	slot_e1.push_back(slot_e01[1]); slot_e2.push_back(slot_e02[1]);
	for (int i = 2; i < slot_e01.size(); i++) {
		c_temp = compute_circle_4_point3d(slot_e01[i], slot_e01[i - 1], slot_e01[i - 2]);
		//std::cout << i << ": " << c_temp.r << endl;
		if (abs(c_temp.r) > 50) {
			slot_e1.push_back(slot_e01[i]);
			slot_e2.push_back(slot_e02[i]);
		}
		else if (abs(c_temp.r) <= 100) {
			break;
		}


	}
}

cad3d::Point3D cad3d::find_correlated_point(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, int i_vt) {

	vector<cad3d::Point3D>  intersect_pt, pt_of_int;
	vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;

	cad3d::Point3D oppo_p;
	oppo_p.x = 0.0; oppo_p.y = 0.0; oppo_p.z = 0.0;

	//int i_vt = 18;
	//base edge
	//intersect_pt.push_back(slot_ip1[i_vt]);
	calculate_cross_rays_v4(slot_e1, slot_e3, slot_ip1, n1_ray1, n1_ray2, n2_ray1, n2_ray2, i_vt);
	//for single point
	/*std::cout << "n1 and n2 rays" << endl;
	std::cout << n1_ray1[0].x << ", " << n1_ray2[0].x <<", ";
	std::cout << n2_ray1[0].x << ", "<< n2_ray2[0].x << endl;*/
	int point_on_line = 0;
	float dtemp, dd;
	vector<float> dmat1;
	//unit vector distance minimize for first ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n1_ray1[0], n1_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	int i_d_min = find_1d_vector_min(dmat1, 1);
	dmat1.resize(0);
	//intersect_pt.push_back(slot_ip1[3]);
	intersect_pt.push_back(slot_e3[i_d_min]);

	//unit vector distance minimize for second ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n2_ray1[0], n2_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	i_d_min = find_1d_vector_min(dmat1, 1);
	intersect_pt.push_back(slot_e3[i_d_min]);

	for (int i = 0; i < intersect_pt.size(); i++) {
		oppo_p.x = oppo_p.x + intersect_pt[i].x;
		oppo_p.y = oppo_p.y + intersect_pt[i].y;
		oppo_p.z = oppo_p.z + intersect_pt[i].z;
	}
	oppo_p.x = oppo_p.x / intersect_pt.size();
	oppo_p.y = oppo_p.y / intersect_pt.size();
	oppo_p.z = oppo_p.z / intersect_pt.size();
	//std::cout << "" << intersect
	return oppo_p;
}
cad3d::Point3D cad3d::find_correlated_point_v3(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, int i_vt) {

	vector<cad3d::Point3D>  intersect_pt, pt_of_int;
	vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;

	cad3d::Point3D oppo_p;
	oppo_p.x = 0.0; oppo_p.y = 0.0; oppo_p.z = 0.0;

	
	calculate_cross_rays_v4(slot_e1, slot_e3, slot_ip1, n1_ray1, n1_ray2, n2_ray1, n2_ray2, i_vt);
	//for single point
	//std::cout << "n1 and n2 rays~" << n1_ray1.size() << ", " << n1_ray2.size() << endl;
	bool a, b;
	a = n1_ray1.size() == 0;
	b = n1_ray2.size() == 0;
	if (a||b) {
		return oppo_p;
	}
	//std::cout << n1_ray1[0].x << ", " << n1_ray2[0].x <<", ";
	//std::cout << n2_ray1[0].x << ", "<< n2_ray2[0].x << endl;
	int point_on_line = 0;
	float dtemp, dd;
	vector<float> dmat1;
	//unit vector distance minimize for first ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n1_ray1[0], n1_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	int i_d_min = find_1d_vector_min(dmat1, 1);
	dmat1.resize(0);
	//intersect_pt.push_back(slot_ip1[3]);
	intersect_pt.push_back(slot_e3[i_d_min]);

	//unit vector distance minimize for second ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n2_ray1[0], n2_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	i_d_min = find_1d_vector_min(dmat1, 1);
	intersect_pt.push_back(slot_e3[i_d_min]);
	//std::cout << "intersect size: " << intersect_pt.size() << endl;
	for (int i = 0; i < intersect_pt.size(); i++) {
		oppo_p.x = oppo_p.x + intersect_pt[i].x;
		oppo_p.y = oppo_p.y + intersect_pt[i].y;
		oppo_p.z = oppo_p.z + intersect_pt[i].z;
	}
	oppo_p.x = oppo_p.x / intersect_pt.size();
	oppo_p.y = oppo_p.y / intersect_pt.size();
	oppo_p.z = oppo_p.z / intersect_pt.size();
	
	return oppo_p;
}
cad3d::Point3D cad3d::find_correlated_point_v2(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, int i_vt) {

	vector<cad3d::Point3D>  intersect_pt, pt_of_int;
	vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;

	cad3d::Point3D oppo_p;
	oppo_p.x = 0.0; oppo_p.y = 0.0; oppo_p.z = 0.0;

	//int i_vt = 18;
	//base edge
	//intersect_pt.push_back(slot_ip1[i_vt]);
	calculate_cross_rays_v4(slot_e1, slot_e3, slot_ip1, n1_ray1, n1_ray2, n2_ray1, n2_ray2, i_vt);
	//for single point
	std::cout << "n1 and n2 ray size" << endl;
	std::cout << n1_ray1.size() << ", " << n1_ray2.size()<< ", ";
	std::cout << n2_ray1.size() << ", " << n2_ray2.size() << endl;
	bool a, b, c, d;
	a = n1_ray1.size() == 0;
	b = n1_ray2.size() == 0;
	c = n2_ray1.size() == 0;
	d = n2_ray2.size() == 0;
	int point_on_line = 0;
	if ((a||b)||(c||d)) {
		return oppo_p;
	}
	float dtemp, dd;
	vector<float> dmat1;
	//unit vector distance minimize for first ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n1_ray1[0], n1_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	int i_d_min = find_1d_vector_min(dmat1, 1);
	dmat1.resize(0);
	//intersect_pt.push_back(slot_ip1[3]);
	intersect_pt.push_back(slot_e3[i_d_min]);

	//unit vector distance minimize for second ray
	for (int i = 0; i < slot_e3.size(); i++) {
		dtemp = DistancePointLine_modv1(slot_e3[i], n2_ray1[0], n2_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	i_d_min = find_1d_vector_min(dmat1, 1);
	intersect_pt.push_back(slot_e3[i_d_min]);

	for (int i = 0; i < intersect_pt.size(); i++) {
		oppo_p.x = oppo_p.x + intersect_pt[i].x;
		oppo_p.y = oppo_p.y + intersect_pt[i].y;
		oppo_p.z = oppo_p.z + intersect_pt[i].z;
	}
	oppo_p.x = oppo_p.x / intersect_pt.size();
	oppo_p.y = oppo_p.y / intersect_pt.size();
	oppo_p.z = oppo_p.z / intersect_pt.size();

	return oppo_p;
}
void cad3d::calculate_fplane_slot_centreline(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, vector<Point3D> &c_line, int cl_n) {
	//need to make more robust and flexible
	vector<cad3d::Point3D> intersect_pt, pt_of_int;
	vector<cad3d::Point3D> cross_l1, cross_l2;
	cad3d::Point3D o_temp, c_temp, temp0;
	c_temp.x = 0.0; c_temp.y = 0.0; c_temp.z = 0.0;
	o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
	temp0.x = 0.0; temp0.y = 0.0; temp0.z = 0.0;
	//starting point
	//starting from the open slot
	c_temp.x = (slot_e1[0].x + slot_e3[0].x) / 2;
	c_temp.y = (slot_e1[0].y + slot_e3[0].y) / 2;
	c_temp.z = (slot_e1[0].z + slot_e3[0].z) / 2;
	c_line.push_back(c_temp);
	float d1;
	for (int i_vt = 1; i_vt < cl_n; i_vt++) {
		//intersect_pt.push_back(slot_ip1[i_vt]);
		//cross_l1.push_back(slot_ip1[i_vt]);
		o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
		o_temp = find_correlated_point(slot_e1, slot_e3, slot_ip1, i_vt);
		//intersect_pt.push_back(o_temp);
		//cross_l2.push_back(o_temp);
		d1 = dist201(o_temp, temp0);
		//std::cout << i_vt << ": " << d1 << endl;
		if (d1 != 0) {
			c_temp.x = (slot_ip1[i_vt].x + o_temp.x) / 2;
			c_temp.y = (slot_ip1[i_vt].y + o_temp.y) / 2;
			c_temp.z = (slot_ip1[i_vt].z + o_temp.z) / 2;
			c_line.push_back(c_temp);
		}
	}

	int l1, l3;
	l1 = slot_e1.size();
	l3 = slot_e3.size();
	c_temp.x = (slot_e1[l1 - 1].x + slot_e3[l3 - 1].x) / 2;
	c_temp.y = (slot_e1[l1 - 1].y + slot_e3[l3 - 1].y) / 2;
	c_temp.z = (slot_e1[l1 - 1].z + slot_e3[l3 - 1].z) / 2;
	c_line.push_back(c_temp);
}
void cad3d::calculate_fplane_slot_centreline_v2(vector<Point3D> &slot_e1, vector<Point3D> &slot_e3, vector<Point3D> &slot_ip1, vector<Point3D> &c_line, int cl_n) {
	//for path symmetric structures
	/*for (int i = 2; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}*/
	//finding connected lines for slots
	Point3D n0, n1, n2, n_temp, o_temp;
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
	//int cl_n;
	std::cout << "Input edge size: " << slot_ip1.size() << endl;
	std::cout << "Enter maximum index of input for centerline:";
	std::cin >> cl_n;
	calculate_fplane_slot_centreline(slot_e1, slot_e3, slot_ip1, c_line, cl_n);
	
}
void cad3d::calculate_fplane_slot_centreline_v3(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line) {
	//fplane slot v13 starts here
	vector<cad3d::Point3D> intersect_pt, pt_of_int;
	vector<cad3d::Point3D> cross_l1, cross_l2;
	cad3d::Point3D o_temp, c_temp, temp0;
	vector<cad3d::Point3D> slot_ip11, slot_ip12;
	int l1, l3;
	l1 = interpolated_edge_v1.size();
	l3 = interpolated_edge_v2.size();
	for (int i = 0; i < interpolated_edge_v1.size(); i++) {
		slot_ip11.push_back(interpolated_edge_v1[i]);
		//i = i + 2;
	}
	std::cout << "slot_ip11 calculation completed" << endl;

	c_temp.x = 0.0; c_temp.y = 0.0; c_temp.z = 0.0;
	o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
	temp0.x = 0.0; temp0.y = 0.0; temp0.z = 0.0;
	//starting point
	c_temp.x = (interpolated_edge_v1[0].x + interpolated_edge_v2[0].x) / 2;
	c_temp.y = (interpolated_edge_v1[0].y + interpolated_edge_v2[0].y) / 2;
	c_temp.z = (interpolated_edge_v1[0].z + interpolated_edge_v2[0].z) / 2;
	c_line.push_back(c_temp);
	std::cout << "1st centeline calculated" << endl;
	float d1, d2;;
	for (int i_vt = 1; i_vt < (l1 - 2); i_vt++) {
		//intersect_pt.push_back(slot_ip1[i_vt]);
		//cross_l1.push_back(slot_ip1[i_vt]);
		o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
		o_temp = find_correlated_point_v3(interpolated_edge_v1, interpolated_edge_v2, slot_ip11, i_vt);
		//intersect_pt.push_back(o_temp);
		//cross_l2.push_back(o_temp);
		d1 = dist201(o_temp, temp0);
		//std::cout << i_vt << ": " << sqrt(d1) << "~" << o_temp.z << endl;
		if (d1 != 0) {
			c_temp.x = (interpolated_edge_v1[i_vt].x + o_temp.x) / 2;
			c_temp.y = (interpolated_edge_v1[i_vt].y + o_temp.y) / 2;
			c_temp.z = (interpolated_edge_v1[i_vt].z + o_temp.z) / 2;
			c_line.push_back(c_temp);
		}
		d2 = 0;
	}
	std::cout << "main centreline calculation completed" << endl;
	c_temp.x = (interpolated_edge_v1[l1 - 1].x + interpolated_edge_v2[l3 - 1].x) / 2;
	c_temp.y = (interpolated_edge_v1[l1 - 1].y + interpolated_edge_v2[l3 - 1].y) / 2;
	c_temp.z = (interpolated_edge_v1[l1 - 1].z + interpolated_edge_v2[l3 - 1].z) / 2;
	c_line.push_back(c_temp);
	std::cout << "end centeline calculation completed" << endl;
	
	//fplane slot v3 ends here 
}
void cad3d::calculate_fplane_slot_centreline_v4(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line) {
	//fplane slot v13 starts here
	vector<cad3d::Point3D> intersect_pt, pt_of_int;
	vector<cad3d::Point3D> cross_l1, cross_l2;
	cad3d::Point3D o_temp, c_temp, temp0;
	vector<cad3d::Point3D> slot_ip11, slot_ip12;
	int l1, l3;
	l1 = interpolated_edge_v1.size();
	l3 = interpolated_edge_v2.size();
	/*for (int i = 0; i < interpolated_edge_v1.size(); i++) {
		slot_ip11.push_back(interpolated_edge_v1[i]);
		//i = i + 2;
	}*/
	std::cout << "centerline calculation module" << endl;

	c_temp.x = 0.0; c_temp.y = 0.0; c_temp.z = 0.0;
	o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
	temp0.x = 0.0; temp0.y = 0.0; temp0.z = 0.0;
	//starting point
	c_temp.x = (interpolated_edge_v1[0].x + interpolated_edge_v2[0].x) / 2;
	c_temp.y = (interpolated_edge_v1[0].y + interpolated_edge_v2[0].y) / 2;
	c_temp.z = (interpolated_edge_v1[0].z + interpolated_edge_v2[0].z) / 2;
	c_line.push_back(c_temp);
	std::cout << "1st centeline calculated" << endl;
	float d1, d2;
	int i_temp = 0;
	for (int i_vt = 1; i_vt < (l1 - 2); i_vt++) {
		//intersect_pt.push_back(slot_ip1[i_vt]);
		//cross_l1.push_back(slot_ip1[i_vt]);
		o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
		i_temp = calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, i_vt);
		//std::cout << i_vt << ": " << sqrt(d1) << "~" << o_temp.z << endl;
		
		c_temp.x = (interpolated_edge_v1[i_vt].x + interpolated_edge_v2[i_temp].x) / 2;
		c_temp.y = (interpolated_edge_v1[i_vt].y + interpolated_edge_v2[i_temp].y) / 2;
		c_temp.z = (interpolated_edge_v1[i_vt].z + interpolated_edge_v2[i_temp].z) / 2;
		c_line.push_back(c_temp);
		//}
		//d2 = 0;
	}
	std::cout << "main centreline calculation completed" << endl;
	c_temp.x = (interpolated_edge_v1[l1 - 1].x + interpolated_edge_v2[l3 - 1].x) / 2;
	c_temp.y = (interpolated_edge_v1[l1 - 1].y + interpolated_edge_v2[l3 - 1].y) / 2;
	c_temp.z = (interpolated_edge_v1[l1 - 1].z + interpolated_edge_v2[l3 - 1].z) / 2;
	c_line.push_back(c_temp);
	std::cout << "end centeline calculation completed" << endl;

	//fplane slot v3 ends here 
}
/*void cad3d::calculate_fplane_slot_centreline_v4(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line) {
	//fplane slot v13 starts here
	vector<cad3d::Point3D> intersect_pt, pt_of_int;
	vector<cad3d::Point3D> cross_l1, cross_l2;
	vector<cad3d::Point3D> intersect_line1, intersect_line2;
	cad3d::Point3D o_temp, c_temp, temp0;
	vector<cad3d::Point3D> slot_ip11, slot_ip12;
	int non_planar_point_stat = 0;
	cad3d::Point3D ptemp3;
	int single_cross = 0;
	//vector<cad3d::Point3D>  intersect_pt, pt_of_int;
	vector<cad3d::Point3D>  n0_ray1, n0_ray2, n1_ray1, n1_ray2, n2_ray1, n2_ray2;
	std::cout << "Enter positive integer index for cross section extraction or zero for none" << endl;
	std::cin >> single_cross;
	//ptemp3 = newblock.find_correlated_point(slot_e1, slot_e3, slot_ip1, single_cross);
	ptemp3 = find_correlated_point(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, single_cross);
	//single_cs.push_back(ptemp3);
	intersect_line2.push_back(ptemp3);
	//generate_single_cross_along_centreline_v1(block_a, slot_plane1, ptemp3, non_planar_point_stat, single_cs, ti_all2, slot_ip2, single_cross);
	//centreline calculation in phase 3
	std::cout << "interpolated sizes: " << interpolated_edge_v1.size() << ", " << interpolated_edge_v2.size() << endl;
	//cross_line1.push_back(interpolated_edge_v1[single_cross]);
	//ptemp3 = newblock.find_correlated_point_v3(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, single_cross);
	//cross_line2.push_back(interpolated_edge_v2[single_cross]);
	calculate_cross_rays_v4(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, n1_ray1, n1_ray2, n2_ray1, n2_ray2, single_cross);
	int point_on_line = 0;
	float dtemp, dd;
	vector<float> dmat1;
	//unit vector distance minimize for first ray
	for (int i = 0; i < interpolated_edge_v2.size(); i++) {
		dtemp = DistancePointLine_modv1(interpolated_edge_v2[i], n1_ray1[0], n1_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	cross_l1.push_back(interpolated_edge_v1[single_cross]);
	int i_d_min = find_1d_vector_min(dmat1, 1);
	cross_l2.push_back(interpolated_edge_v2[i_d_min]);
	dmat1.resize(0);
	//intersect_pt.push_back(slot_ip1[3]);
	intersect_pt.push_back(interpolated_edge_v2[i_d_min]);
	//unit vector distance minimize for second ray
	for (int i = 0; i < interpolated_edge_v2.size(); i++) {
		dtemp = DistancePointLine_modv1(interpolated_edge_v2[i], n2_ray1[0], n2_ray2[0], dd);
		//std::cout << i << ": ";
		//point_on_line = 0;
		dmat1.push_back(dtemp);
	}
	cross_l1.push_back(interpolated_edge_v1[single_cross]);
	i_d_min = find_1d_vector_min(dmat1, 1);
	cross_l2.push_back(interpolated_edge_v2[i_d_min]);
	intersect_pt.push_back(interpolated_edge_v2[i_d_min]);
	std::cout << "cross size: " << n1_ray1.size() << ", " << n1_ray2.size() << ", " << n2_ray1.size() << ", " << n2_ray2.size() << endl;
	
}*/
void cad3d::stl_axis_tx_v1(stl_data &block_a, stl_data &block_b) {
	float temp1 = 0.0;
	int dir;
	std::cout << "enter 1*, 2 or 3 for convenient direction ttansformation" << endl;
	cin >> dir;
	if (dir ==1) {
		for (int i = 0; i < block_b.triangles.size(); i++) {
			block_a.triangles.push_back(block_b.triangles[i]);
		}
		block_a.vertices.resize(block_b.vertices.size());
		for (int i = 0; i < block_b.vertices.size(); i++) {
			block_a.vertices[i].z = block_b.vertices[i].x;
			block_a.vertices[i].x = block_b.vertices[i].y;
			block_a.vertices[i].y = block_b.vertices[i].z;
		}
	}
}
void cad3d::stl_axis_tx_v2(stl_data &block_a, stl_data &block_b) {
	float temp1 = 0.0;
	int dir;
	std::cout << "enter 1*, 2 or 3 for convenient direction ttansformation" << endl;
	cin >> dir;
	int ii;
	vector<cad3d::triangle_index> ti_all;
	
	Point3D p1, p2, p3,pn;
	
	if (dir == 1) {
		//block_a.triangles.resize(block_b.triangles.size());
		//triangle_vertexindex_init(block_b.triangles, ti_all);
		for (int i = 0; i < block_b.triangles.size(); i++) {
			p1 = block_b.triangles[i].v1;
			p2 = block_b.triangles[i].v2;
			p3 = block_b.triangles[i].v3;
			block_a.triangles.push_back(block_b.triangles[i]);
			//ii = ti_all[i].ti;
			block_a.triangles[i].v1.z = p1.x;
			block_a.triangles[i].v1.x = p1.y;
			block_a.triangles[i].v1.y = p1.z;
			
			block_a.triangles[i].v2.z = p2.x;
			block_a.triangles[i].v2.x = p2.y;
			block_a.triangles[i].v2.y = p2.z;

			block_a.triangles[i].v3.z = p3.x;
			block_a.triangles[i].v3.x = p3.y;
			block_a.triangles[i].v3.y = p3.z;
		
		}
		block_a.vertices.resize(block_b.vertices.size());
		for (int i = 0; i < block_b.vertices.size(); i++) {
			block_a.vertices[i].z = block_b.vertices[i].x;
			block_a.vertices[i].x = block_b.vertices[i].y;
			block_a.vertices[i].y = block_b.vertices[i].z;
		}
	}
}
void cad3d::cartesean_processing_v3(stl_data &block_a, vector<Point3D> &slot_plane, vector<triangle_index> &tslot_i, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2) {
	//for crossection and closed region identification
	//getplane_simplfied_block(block_a, slot_plane, tslot_i);
	float z_now = 0.0;
	std::cin >> z_now;
	vector<int> slot_plane_i;
	find_point_3d_1feature(slot_plane, block_a.vertices, z_now, 0.5, 3);
	slot_plane_i.resize(block_a.vertices.size());
	find_point3d_i_1feature(slot_plane_i, block_a.vertices, z_now, 0.5, 3);
	show_index_frequency(slot_plane_i);
	/*triangle_normal_update(block_a.triangles, tslot_i);
	int ii;
	bool a, b;

	vector<cad3d::Point3D> v_temp, v_line1, v_line2;
	vector<int> i_line;
	triangle2line(block_a.vertices, tslot_i, v_line1, v_line2, i_line, v_temp);
	//newblock.triangle2line(block_a.vertices, tbi_1, v_line1, v_line2, i_line, v_temp);
	
	cad3d::Point3D line_value;

	vector<cad3d::Point3D>  vord1_line1, vord1_line2;
	cad3d::Point3D ordt1, ordt2;
	float delx1, dely1;
	//ordering lines
	//for rectangular block
	for (int i = 0; i < v_line1.size(); i++) {
		delx1 = v_line2[i].x - v_line1[i].x;
		if (delx1 > 0) {
			ordt2 = v_line2[i];
			ordt1 = v_line1[i];
		}
		else if (delx1 < 0) {
			ordt2 = v_line1[i];
			ordt1 = v_line2[i];
		}
		else {
			dely1 = v_line2[i].y - v_line1[i].y;
			if (dely1 > 0) {
				ordt2 = v_line2[i];
				ordt1 = v_line1[i];
			}
			else if (dely1 < 0) {
				ordt2 = v_line1[i];
				ordt1 = v_line2[i];
			}
		}
		vord1_line1.push_back(ordt1);
		vord1_line2.push_back(ordt2);

	}

	for (int i = 0; i < vord1_line1.size(); i++) {
		//for (int i = 0; i < tslot_i.size(); i++) {
		//for (int i = 0; i < 100; i++) {
		line_value = feature_edge_predict_phase1(vord1_line1[i], vord1_line2[i]);
		a = line_value.x < 20.0;
		b = line_value.y < 20.0;
		//b = 0;
		if (a || b) {
			if (sqrt(line_value.x*line_value.x + line_value.y*line_value.y) < 100) {
				v_ext_l1.push_back(vord1_line1[i]);
				v_ext_l2.push_back(vord1_line2[i]);
			}

			//std::cout << sqrt(line_value.x*line_value.x + line_value.y*line_value.y) << endl;
		}
		if (line_value.x == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}
		if (line_value.y == 0.0) {
			v_ext_l1.push_back(vord1_line1[i]);
			v_ext_l2.push_back(vord1_line2[i]);
		}

	}
	feature_edge_predict_phase2(v_ext_l1, v_ext_l2, feat_e1, feat_e2);*/
}

void cad3d::cross_model_generate(stl_data &block_a, vector<Point3D> &slot_plane1, vector<triangle_index> &ti_all2, vector<triangle_index> &tslot_i1, vector<triangle_index> &tslot_i2,vector<triangle_index> &tsloti_feature) {
	float z_int;
	std::cout << "z of target layer: " << endl;
	std::cin >> z_int;
	find_point_3d_1feature(slot_plane1, block_a.vertices, z_int, 0.5, 3);
	//std::cout << "points in the intermediate slice plane: " << slot_plane1.size() << endl;

	triangle_vertexindex_init(block_a.triangles, ti_all2);
	triangle_anypoint_i1feature(block_a.vertices, ti_all2, tslot_i1, z_int, 0.5, 3);
	triangle_allpoint_i1feature(block_a.vertices, ti_all2, tslot_i2, z_int, 0.5, 3);
	//std::cout << "triangle size diff: " << tslot_i1.size() << ", " << tslot_i2.size() << endl;
	triangle_normal_update(block_a.triangles, tslot_i1);

	int ii;
	cad3d::Point3D p1, p2, p3;
	int a, b, c;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = block_a.triangles[ii].v1;
		p2 = block_a.triangles[ii].v2;
		p3 = block_a.triangles[ii].v3;
		//std::cout << i << ": " << p1.z << ", ";
		//std::cout << p2.z << ", " << p3.z << endl;
		a = abs(p1.z - z_int) < 0.5;
		b = abs(p2.z - z_int) < 0.5;
		c = abs(p3.z - z_int) < 0.5;
		if (b*c*(!a) == 1) {
			//std::cout << ii << ": !a.b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (c*a*(!b) == 1) {
			//std::cout << ii << ": a.!b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (a*b*(!c) == 1) {
			//std::cout << ii << ": a.b.!c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		a = 0; b = 0; c = 0;
	}
}
void cad3d::cross_model_generate_v2(stl_data &block_a, vector<Point3D> &slot_plane1, vector<triangle_index> &ti_all2, vector<triangle_index> &tslot_i1, vector<triangle_index> &tslot_i2, vector<triangle_index> &tsloti_feature, float tol) {
	float z_int;
	std::cout << "z of target layer: " << endl;
	std::cin >> z_int;
	find_point_3d_1feature(slot_plane1, block_a.vertices, z_int, tol, 3);
	//std::cout << "points in the intermediate slice plane: " << slot_plane1.size() << endl;

	triangle_vertexindex_init(block_a.triangles, ti_all2);
	triangle_anypoint_i1feature(block_a.vertices, ti_all2, tslot_i1, z_int, tol, 3);
	//triangle_allpoint_i1feature(block_a.vertices, ti_all2, tslot_i2, z_int, 0.5, 3);
	//std::cout << "triangle size diff: " << tslot_i1.size() << ", " << tslot_i2.size() << endl;
	triangle_normal_update(block_a.triangles, tslot_i1);

	int ii;
	cad3d::Point3D p1, p2, p3;
	int a, b, c;
	a = 0; b = 0; c = 0;
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = block_a.triangles[ii].v1;
		p2 = block_a.triangles[ii].v2;
		p3 = block_a.triangles[ii].v3;
		//std::cout << i << ": " << p1.z << ", ";
		//std::cout << p2.z << ", " << p3.z << endl;
		a = abs(p1.z - z_int) < 0.5;
		b = abs(p2.z - z_int) < 0.5;
		c = abs(p3.z - z_int) < 0.5;
		if (b*c*(!a) == 1) {
			//std::cout << ii << ": !a.b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (c*a*(!b) == 1) {
			//std::cout << ii << ": a.!b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (a*b*(!c) == 1) {
			//std::cout << ii << ": a.b.!c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		a = 0; b = 0; c = 0;
	}
}
void cad3d::cross_model_generate_polar(stl_data &block_a, vector<Point3D> &centered_block, vector<Point3D> &slot_plane1, vector<triangle_index> &ti_all2, vector<triangle_index> &tslot_i1, vector<triangle_index> &tslot_i2, vector<triangle_index> &tsloti_feature) {
	//updated 29.11.19
	float z_int;
	std::cout << "z_int: " << endl;
	std::cin >> z_int;
	find_point_3d_1feature(slot_plane1, centered_block, z_int, 0.5, 3);
	std::cout << "points in the intermediate slice plane: " << slot_plane1.size() << endl;

	triangle_vertexindex_init(block_a.triangles, ti_all2);
	triangle_anypoint_i1feature(block_a.vertices, ti_all2, tslot_i1, z_int, 15, 3);
	triangle_allpoint_i1feature(block_a.vertices, ti_all2, tslot_i2, z_int, 15, 3);
	std::cout << "triangle size diff: " << tslot_i1.size() << ", " << tslot_i2.size() << endl;
	triangle_normal_update(block_a.triangles, tslot_i1);

	int ii;
	cad3d::Point3D p1, p2, p3;
	int a, b, c;
	a = 0; b = 0; c = 0;
	//centre_at_x0y0(block_a.vertices, centered_block);
	//triangle points transfer to centre
	vector<triangle> t_centred;
	triangle_centre_at_x0y0(block_a.vertices, block_a.triangles, t_centred);
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = copy_point3d(t_centred[ii].v1);
		p2 = copy_point3d(t_centred[ii].v2);
		p3 = copy_point3d(t_centred[ii].v3);
		//std::cout << i << ": " << p1.z << ", ";
		//std::cout << p2.z << ", " << p3.z << endl;
		a = abs(p1.z - z_int) < 0.5;
		b = abs(p2.z - z_int) < 0.5;
		c = abs(p3.z - z_int) < 0.5;
		if (b*c*(!a) == 1) {
			//std::cout << ii << ": !a.b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (c*a*(!b) == 1) {
			//std::cout << ii << ": a.!b.c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		if (a*b*(!c) == 1) {
			//std::cout << ii << ": a.b.!c" << endl;
			tsloti_feature.push_back(tslot_i1[i]);
		}
		a = 0; b = 0; c = 0;
	}
}

void cad3d::cartesian_separate_internal_regions(vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_int_l1, vector<Point3D> &v_int_l2) {
	int int_boundary = 0;
	std::cout << "internal region number: ";
	std::cin >> int_boundary;

	bool a1, b1, c1, d1, e, f;
	bool a2, b2, c2, d2;
	//corner point estmation 
	int ixmax, iymax, ixmin, iymin;
	ixmax = find_point3d_max(feat1_e1, 1);
	ixmin = find_point3d_min(feat1_e1, 1);
	iymax = find_point3d_max(feat1_e1, 2);
	iymin = find_point3d_min(feat1_e1, 2);
	for (int i = 0; i < feat1_e1.size(); i++) {
		a1 = abs(feat1_e1[i].x - feat1_e1[ixmax].x) < 0.1;
		b1 = abs(feat1_e1[i].x - feat1_e1[ixmin].x) < 0.1;
		c1 = abs(feat1_e1[i].y - feat1_e1[iymax].y) < 0.1;
		d1 = abs(feat1_e1[i].y - feat1_e1[iymin].y) < 0.1;

		e = (a1 || b1) || (c1 || d1);

		a2 = abs(feat1_e2[i].x - feat1_e1[ixmax].x) < 0.1;
		b2 = abs(feat1_e2[i].x - feat1_e1[ixmin].x) < 0.1;
		c2 = abs(feat1_e2[i].y - feat1_e1[iymax].y) < 0.1;
		d2 = abs(feat1_e2[i].y - feat1_e1[iymin].y) < 0.1;

		f = (a2 || b2) || (c2 || d2);

		if (!e * !f == 1) {
			//std::cout << i << "e.f" << endl;
			v_int_l1.push_back(feat1_e1[i]);
			v_int_l2.push_back(feat1_e2[i]);
		}


		a1 = 0; b1 = 0; c1 = 0; d1 = 0;
		a2 = 0; b2 = 0; c2 = 0; d2 = 0;
		e = 0; f = 0;
	}
	std::cout << "v_int_l1 size: " << v_int_l1.size() << endl;
}
void cad3d::cartesian_separate_internal_regions_v2(vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_int_l1, vector<Point3D> &v_int_l2, int cp_mode, vector<Point3D> &slot_plane1) {
	//cp mode is Cartesian 0 or polar 1 mode
	int int_boundary = 0;
	std::cout << "internal region number: ";
	std::cin >> int_boundary;

	if (cp_mode == 0) {
		bool a1, b1, c1, d1, e, f;
		bool a2, b2, c2, d2;
		//corner point estmation 
		int ixmax, iymax, ixmin, iymin;
		
		ixmax = find_point3d_max(slot_plane1, 1);
		ixmin = find_point3d_min(slot_plane1, 1);
		iymax = find_point3d_max(slot_plane1, 2);
		iymin = find_point3d_min(slot_plane1, 2);
		//std::cout << slot_plane1[ixmax].x << ", " << slot_plane1[ixmin].x << ", " << slot_plane1[iymax].y << ", " << slot_plane1[iymin].y << endl;
		for (int i = 0; i < feat1_e1.size(); i++) {
			a1 = abs(feat1_e1[i].x - slot_plane1[ixmax].x) < 0.1;
			b1 = abs(feat1_e1[i].x - slot_plane1[ixmin].x) < 0.1;
			c1 = abs(feat1_e1[i].y - slot_plane1[iymax].y) < 0.1;
			d1 = abs(feat1_e1[i].y - slot_plane1[iymin].y) < 0.1;

			e = (a1 || b1) || (c1 || d1); // at any corner point e = 1

			a2 = abs(feat1_e2[i].x - slot_plane1[ixmax].x) < 0.1;
			b2 = abs(feat1_e2[i].x - slot_plane1[ixmin].x) < 0.1;
			c2 = abs(feat1_e2[i].y - slot_plane1[iymax].y) < 0.1;
			d2 = abs(feat1_e2[i].y - slot_plane1[iymin].y) < 0.1;
			
			f = (a2 || b2) || (c2 || d2); // at any corner point f = 1
			//std::cout << i << ": " << e << ", " << f;
			//if any of the point is internal, include in internal feature edges
			if ((e == 0) || (f == 0)) {
				//std::cout << " internal " << feat1_e1[i].x << ", "<<feat1_e1[i].y;
				//std::cout << " ~ " << feat1_e2[i].x << ", " << feat1_e2[i].y;
				v_int_l1.push_back(feat1_e1[i]);
				v_int_l2.push_back(feat1_e2[i]);
			}
			/*else if (!f ==1) {
				v_int_l1.push_back(feat1_e1[i]);
				v_int_l2.push_back(feat1_e2[i]);
			}*/
			//std::cout << endl;
			a1 = 0; b1 = 0; c1 = 0; d1 = 0;
			a2 = 0; b2 = 0; c2 = 0; d2 = 0;
			e = 0; f = 0;
		}
	}
	std::cout << "v_int_l1 size: " << v_int_l1.size() << endl;
}
void cad3d:: single_closed_region_internal_boundary_extract(stl_data &block_a, vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_int_l1, vector<Point3D> &v1, vector<Point3D> &slot_e01, vector<Point3D> &slot_e02) {
	std::cout << "internal vertex array size: " << v_int_l1.size()<< endl;
	int i_st1 = find_point3d_max(v_int_l1, 2);
	std::cout << i_st1 << endl;
	int i_start = 0;
	vector<int> vi_temp551, vi_temp552;
	std::cout << "internal edges starting index: " << endl;
	std::cin >> i_start;
	std::cout << i_start;
	v1.push_back(v_int_l1[i_start]);//i_st1
	slot_e01.push_back(v_int_l1[i_start]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, v_int_l1[i_start], feat1_e1, feat1_e2, v1);
	slot_e02.push_back(v1[v1.size() - 1]);
	//newblock.get_internal_slot_v106(block_a, vi_all, slot_e1, slot_e2, feat_e1, feat_e2, vi_temp552, v1);
	//For GE blocks
	//newblock.get_internal_slot_v107_phase4(block_a.vertices, slot_e01, slot_e02, feat1_e1, feat1_e2, v1);

	get_internal_slot_v107_phase6(block_a.vertices, slot_e01, slot_e02, feat1_e1, feat1_e2, v1);
	//newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	std::cout << "slot size: " << slot_e01.size() << ", " << slot_e02.size() << endl;
}
void cad3d::single_closed_region_internal_boundary_extract_v2(stl_data &block_a, vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_int_l1, vector<Point3D> &v1, vector<Point3D> &slot_e01, vector<Point3D> &slot_e02) {
	std::cout << "internal vertex array size: " << v_int_l1.size() << endl;
	//int i_st1 = find_point3d_max(v_int_l1, 2);
	//std::cout << i_st1 << endl;
	int i_start = 0;
	vector<int> vi_temp551, vi_temp552;
	std::cout << "internal edges starting index: " << endl;
	std::cin >> i_start;
	//std::cout << i_start;
	v1.push_back(v_int_l1[i_start]);//i_st1
	slot_e01.push_back(v_int_l1[i_start]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, v_int_l1[i_start], feat1_e1, feat1_e2, v1);
	slot_e02.push_back(v1[v1.size() - 1]);
	//newblock.get_internal_slot_v106(block_a, vi_all, slot_e1, slot_e2, feat_e1, feat_e2, vi_temp552, v1);
	//For GE blocks
	//newblock.get_internal_slot_v107_phase4(block_a.vertices, slot_e01, slot_e02, feat1_e1, feat1_e2, v1);

	get_internal_slot_v107_phase6(block_a.vertices, slot_e01, slot_e02, feat1_e1, feat1_e2, v1);
	//newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	std::cout << "slot size: " << slot_e01.size() << ", " << slot_e02.size() << endl;
}
void cad3d::find_closed_slot_number_v1(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v_used, int internal_v_length) {
	//function has issues of when to stop
	std::cout << "function internal slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	//Point3D p_start_dif;
	int vn = 0;
	int point_reap_val = 0;
	vector<float> int_index_check;
	
	for (int i = 0; i < 100; i++) {
		//std::cout << i << ": " << v_used.size();
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v_used[v_used.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v_used[v_used.size() - 1], feat_e2[j]);

			if (d1 < 0.01) {
				//std::cout << " j d1- " << j << ", " << d1;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << endl;
				jrmat1.push_back(p1_con.x);
				point_reap_val = point_repeat_check_v1(feat_e2[j], v_used);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					v_used.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e2[j]);
					break;
				}
				//
				//vn++;
			}
			else if (d2 < 0.01) {
				//std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				point_reap_val = point_repeat_check_v1(feat_e1[j], v_used);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					v_used.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e1[j]);
					break;
				}
				//vn++;
			}
			vn++;
			point_reap_val = 0;
		}
		//std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		//std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		if ((jr1.size() < 1) & (jr2.size() < 1)) {
			std::cout << "discontinuity" << endl;
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
		bool a, b;
		int i_temp = find_1d_vector_max(int_index_check, 1);
		std::cout << "discontinuity max: " << int_index_check[i_temp] << endl;
		//end closed region
		slot_e3.push_back(v_used[v_used.size() - 1]);
		//v2.push_back(feat_e2[j]);
		slot_e4.push_back(v_used[0]);
		std::cout << dist201(slot_e3[0], v_used[v_used.size() - 1]) << ", " << dist201(slot_e4[0], v_used[v_used.size() - 1]) << endl;
	}
	//return int_index_check[i_temp];
}
float cad3d::find_closed_slot_number_v2(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v_used) {
	//function has issues of when to stop
	std::cout << "function closed slot count\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	//Point3D p_start_dif;
	int vn = 0;
	int point_reap_val = 0;
	vector<float> int_index_check;
	int i_temp = 0;
	for (int i = 0; i < 100; i++) {
		//std::cout << i << ": " << v_used.size();
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist201(v_used[v_used.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist201(v_used[v_used.size() - 1], feat_e2[j]);

			if (d1 < 0.01) {
				//std::cout << " j d1- " << j << ", " << d1;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << endl;
				jrmat1.push_back(p1_con.x);
				point_reap_val = point_repeat_check_v1(feat_e2[j], v_used);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					v_used.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e2[j]);
					break;
				}
				//
				//vn++;
			}
			else if (d2 < 0.01) {
				//std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				point_reap_val = point_repeat_check_v1(feat_e1[j], v_used);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					v_used.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e1[j]);
					break;
				}
				//vn++;
			}
			vn++;
			point_reap_val = 0;
		}
		//std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		//std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		if ((jr1.size() < 1) & (jr2.size() < 1)) {
			std::cout << "discontinuity" << endl;
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
		bool a, b;
		
		
	}
	i_temp = find_1d_vector_max(int_index_check, 1);
	std::cout << "discontinuity max: " << int_index_check[i_temp] << endl;
	//end closed region
	slot_e3.push_back(v_used[v_used.size() - 1]);
	//v2.push_back(feat_e2[j]);
	slot_e4.push_back(v_used[0]);
	std::cout << dist201(slot_e3[0], v_used[v_used.size() - 1]) << ", " << dist201(slot_e4[0], v_used[v_used.size() - 1]) << endl;
	return int_index_check[i_temp];
}
void cad3d::polar_separate_internal_regions_v1(vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_int_l1, vector<Point3D> &v_int_l2, int cp_mode, vector<Point3D> &slot_plane1) {
	//cp mode is Cartesian 0 or polar 1 mode
	
	vector<cad3d::Point3D> feat1_e1_polar, feat1_e2_polar,slot_plane1_polar;
	convert_cart2pol(feat1_e1, feat1_e1_polar);
	convert_cart2pol(feat1_e2, feat1_e2_polar);
	convert_cart2pol(slot_plane1, slot_plane1_polar);

	int int_boundary = 0;
	std::cout << "internal region number: ";
	std::cin >> int_boundary;

	if (cp_mode == 0) {
		bool a1, b1, c1, d1, e, f;
		bool a2, b2, c2, d2;
		//corner point estmation 
		int ixmax,ixmax1, ixmax2;

		ixmax = find_point3d_max(feat1_e1_polar, 1);
		ixmax1 = find_point3d_max(slot_plane1_polar, 1);
		ixmax2 = find_point3d_max(feat1_e2_polar, 1);
		std::cout << "rcomp: "<< slot_plane1_polar[ixmax1].x << ","<< feat1_e2_polar[ixmax2].x << endl;
		/*ixmin = find_point3d_min(slot_plane1, 1);
		iymax = find_point3d_max(slot_plane1, 2);
		iymin = find_point3d_min(slot_plane1, 2);*/
		//std::cout << slot_plane1[ixmax].x << ", " << slot_plane1[ixmin].x << ", " << slot_plane1[iymax].y << ", " << slot_plane1[iymin].y << endl;
		for (int i = 0; i < feat1_e1.size(); i++) {
			a1 = abs(feat1_e1_polar[i].x - slot_plane1_polar[ixmax1].x) < 0.1;
			// if e1 lies on extrnal point
			/*b1 = abs(feat1_e1[i].x - slot_plane1[ixmin].x) < 0.1;
			c1 = abs(feat1_e1[i].y - slot_plane1[iymax].y) < 0.1;
			d1 = abs(feat1_e1[i].y - slot_plane1[iymin].y) < 0.1;*/

			//e = (a1 || b1) || (c1 || d1); // at any corner point e = 1

			a2 = abs(feat1_e2_polar[i].x - slot_plane1_polar[ixmax1].x) < 0.1;
			//if e2 lies on external point
			/*b2 = abs(feat1_e2[i].x - slot_plane1[ixmin].x) < 0.1;
			c2 = abs(feat1_e2[i].y - slot_plane1[iymax].y) < 0.1;
			d2 = abs(feat1_e2[i].y - slot_plane1[iymin].y) < 0.1;*/

			//f = (a2 || b2) || (c2 || d2); // at any corner point f = 1
			//std::cout << i << ": " << e << ", " << f;
			//if any of the point is internal, include in internal feature edges
			if ((a1 == 0) || (a2 == 0)) {
				//std::cout << " internal " << feat1_e1[i].x << ", "<<feat1_e1[i].y;
				//std::cout << " ~ " << feat1_e2[i].x << ", " << feat1_e2[i].y;
				v_int_l1.push_back(feat1_e1[i]);
				v_int_l2.push_back(feat1_e2[i]);
			}
			/*else if (!f ==1) {
				v_int_l1.push_back(feat1_e1[i]);
				v_int_l2.push_back(feat1_e2[i]);
			}*/
			//std::cout << endl;
			a1 = 0; //b1 = 0; c1 = 0; d1 = 0;
			a2 = 0; //b2 = 0; c2 = 0; d2 = 0;
			//e = 0; f = 0;
		}
	}
	std::cout << "v_int_l1 size: " << v_int_l1.size() << endl;
}

void cad3d::get_semi_open_slot_v101(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	//function has issues of when to stop
	std::cout << "function semi open slot\n ";
	//this one is for inward direction from end edge
	float d1, d2;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	//Point3D p_start_dif;
	int vn = v2.size();
	int point_reap_val = 0;
	vector<float> int_index_check;
	//for (int i = 0; i < feat_e1.size(); i++) {
	for (int i = 0; i < 100; i++) {
		//std::cout << i << ": " << v2.size();
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist301(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist301(v2[v2.size() - 1], feat_e2[j]);

			if (d1 < 0.1) {
				//std::cout << " j d1- " << j << ", " << d1;
				jr1.push_back(j);
				convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << endl;
				jrmat1.push_back(p1_con.x);
				point_reap_val = point_repeat_check_v1(feat_e2[j], v2);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					v2.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e2[j]);
					break;
				}
				//
				//vn++;
			}
			else if (d2 < 0.1) {
				//std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << endl;
				jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				point_reap_val = point_repeat_check_v1(feat_e1[j], v2);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					v2.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e1[j]);
					break;
				}
				//vn++;
			}
			vn++;
			point_reap_val = 0;
		}
		//std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		//std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		if ((jr1.size() < 1) & (jr2.size() < 1)) {
			std::cout << "discontinuity" << endl;
			break;
		}
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
	}
	bool a, b;
	int i_temp = find_1d_vector_max(int_index_check, 1);
	std::cout << "discontinuity max: " << int_index_check[i_temp] << endl;
	

}
void cad3d::get_semi_open_slot_v102(vector<Point3D> &centered_block, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4, vector<Point3D> &feat_e1, vector<Point3D> &feat_e2, vector<Point3D> &v2) {
	//function has issues of when to stop
	std::cout << "function semi open slot\n ";
	//this one is for inward direction from end edge
	float d1, d2, d1z, d2z;
	vector<int> jr1, jr2;
	vector<float> jrmat1, jrmat2;
	Point3D p1_con, p2_con;
	//Point3D p_start_dif;
	int vn = v2.size();
	int point_reap_val = 0;
	vector<float> int_index_check;
	for (int i = 0; i < feat_e1.size(); i++) {
	//for (int i = 0; i < 100; i++) {
		//std::cout << i << ": " << v2.size();
		for (int j = 0; j < feat_e1.size(); j++) {
			d1 = dist301(v2[v2.size() - 1], feat_e1[j]);//v_ext_l2[0]
			d2 = dist301(v2[v2.size() - 1], feat_e2[j]);
			//to check if the other point of the line is on the same plane
			d1z = abs(v2[v2.size() - 1].z - feat_e2[j].z);
			d2z = abs(v2[v2.size() - 1].z - feat_e1[j].z);
			if ((d1 < 0.1) & (d1z <0.1)) {
				//std::cout << " j d1- " << j << ", " << d1;
				jr1.push_back(j);
				//convert_cart2polPoint3D(feat_e2[j], p1_con);
				//std::cout << "d1: " << j << endl;
				//jrmat1.push_back(p1_con.x);
				point_reap_val = point_repeat_check_v1(feat_e2[j], v2);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e1[j]);
					v2.push_back(feat_e2[j]);
					slot_e4.push_back(feat_e2[j]);
					break;
				}
				//
				//vn++;
			}
			else if ((d2 < 0.1) & (d2z < 0.1)) {
				//std::cout << " j d2- " << j << ", " << d2;
				jr2.push_back(j);
				//convert_cart2polPoint3D(feat_e1[j], p2_con);
				//std::cout << "d2: " << j << endl;
				//jrmat2.push_back(p2_con.x);
				//jr2.push_back(j);
				point_reap_val = point_repeat_check_v1(feat_e1[j], v2);
				if (point_reap_val == 0) {
					int_index_check.push_back(j);
					slot_e3.push_back(feat_e2[j]);
					v2.push_back(feat_e1[j]);
					slot_e4.push_back(feat_e1[j]);
					break;
				}
				//vn++;
			}
			vn++;
			point_reap_val = 0;
		}
		//std::cout << endl;
		//std::cout << "~ " << d1 << ", " << d2 << endl;
		//std::cout << "~ " << jr1.size() << ", " << jr2.size() << endl;
		/*if ((jr1.size() < 1) & (jr2.size() < 1)) {
			std::cout << "discontinuity" << endl;
			break;
		}*/
		jr1.resize(0);
		jr2.resize(0);
		jrmat1.resize(0);
		jrmat2.resize(0);
		//std::cout << i << ": " << v2[v2.size() - 1].x << ", " << v2[v2.size() - 1].y << " "<< v2[v2.size() - 1].z<< endl;
	}
	bool a, b;
	//int i_temp = find_1d_vector_max(int_index_check, 1);
	//std::cout << "discontinuity max: " << int_index_check[i_temp] << endl;


}

cad3d::Point3D cad3d::find_non_planar_connected_points(stl_data &block_a, vector<Point3D> &slot_plane1, vector<triangle_index> &tslot_i1, Point3D pin) {
	int ii;
	Point3D p1, p2, p3;
	float dd1, dd2, dd3;
	bool a, b, c;
	Point3D pout;
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = block_a.triangles[ii].v1;
		p2 = block_a.triangles[ii].v2;
		p3 = block_a.triangles[ii].v3;
		
		dd1 = dist301(p1, pin);
		dd2 = dist301(p2, pin);
		dd3 = dist301(p3, pin);
		a = abs(dd1) < 0.5;
		b = abs(dd2) < 0.5;
		c = abs(dd3) < 0.5;
		if (a || b || c) {
			if ((pin.z - p1.z) > 1) {
				pout = copy_point3d(p1);
				//std::cout << i << ": " << p1.z << endl;
				break;
			}
			else if ((pin.z - p2.z) > 1) {
				pout = copy_point3d(p2);
				//std::cout << i << ": " << p2.z << endl;
				break;
			}
			else if ((pin.z - p3.z) > 1) {
				pout = copy_point3d(p3);
				//std::cout << i << ": " << p3.z << endl;
				break;
			}
		}
	}
	return pout;
}
cad3d::Point3D cad3d::find_non_planar_connected_points_v2(stl_data &block_a, vector<Point3D> &slot_plane1, vector<triangle_index> &tslot_i1, Point3D pin, int & status) {
	int ii;
	Point3D p1, p2, p3;
	float dd1, dd2, dd3;
	bool a, b, c;
	Point3D pout;
	//std::cout << pin.x << ", "<<pin.z << endl;
	pout = copy_point3d(pin);
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = block_a.triangles[ii].v1;
		p2 = block_a.triangles[ii].v2;
		p3 = block_a.triangles[ii].v3;

		dd1 = dist301(p1, pin);
		dd2 = dist301(p2, pin);
		dd3 = dist301(p3, pin);
		a = abs(dd1) < 0.5;
		b = abs(dd2) < 0.5;
		c = abs(dd3) < 0.5;
		if (a || b || c) {
			//std::cout << i;
			if ((pin.z - p1.z) > 1) {
				pout = copy_point3d(p1);
				//std::cout << i << ": " << p1.z << endl;
				status = 1;
				break;
			}
			else if ((pin.z - p2.z) > 1) {
				pout = copy_point3d(p2);
				//std::cout << i << ": " << p2.z << endl;
				status = 1;
				break;
			}
			else if ((pin.z - p3.z) > 1) {
				pout = copy_point3d(p3);
				//std::cout << i << ": " << p3.z << endl;
				status = 1;
				break;
			}
		}
	}
	return pout;
}
cad3d::Point3D cad3d::find_non_planar_connected_points_v3(stl_data &block_a,  vector<Point3D> &slot_plane1, vector<triangle_index> &tslot_i1, Point3D pin, int & status) {
	int ii;
	Point3D p1, p2, p3;
	float dd1, dd2, dd3;
	bool a, b, c;
	Point3D pout;
	//std::cout << pin.x << ", "<<pin.z << endl;
	pout = copy_point3d(pin);
	vector<triangle> t_centred;
	triangle_centre_at_x0y0(block_a.vertices, block_a.triangles, t_centred);
	for (int i = 0; i < tslot_i1.size(); i++) {
		//for (int i = 0; i < 20; i++) {
		ii = tslot_i1[i].ti;
		p1 = t_centred[ii].v1;
		p2 = t_centred[ii].v2;
		p3 = t_centred[ii].v3;

		dd1 = dist301(p1, pin);
		dd2 = dist301(p2, pin);
		dd3 = dist301(p3, pin);
		a = abs(dd1) < 0.5;
		b = abs(dd2) < 0.5;
		c = abs(dd3) < 0.5;
		if (a || b || c) {
			//std::cout << i;
			if ((pin.z - p1.z) > 1) {
				pout = copy_point3d(p1);
				//std::cout << i << ": " << p1.z << endl;
				status = 1;
				break;
			}
			else if ((pin.z - p2.z) > 1) {
				pout = copy_point3d(p2);
				//std::cout << i << ": " << p2.z << endl;
				status = 1;
				break;
			}
			else if ((pin.z - p3.z) > 1) {
				pout = copy_point3d(p3);
				//std::cout << i << ": " << p3.z << endl;
				status = 1;
				break;
			}
		}
	}
	return pout;
}
cad3d::Point3D cad3d::find_un_connected_non_planar_points(stl_data &block_a, Point3D pin, int & status) {
	int ii;
	Point3D p1, p2, p3;
	float dd1, dd2, dd3;
	bool a, b, c;
	Point3D pout;
	pout = copy_point3d(pin);
	for (int i = 0; i < block_a.vertices.size(); i++) {
		
		p1 = block_a.vertices[i];
		dd1 = sqrt(dist201(p1, pin));
		if (dd1 < 1) {
			if ((pin.z - p1.z) > 0) {
				//std::cout << i << ": " << p1.z << endl;
				pout = copy_point3d(p1);
				status = 1;
				break;
				
			}
			
		}
	}
	//std::cout << pout.x << pout.y << pout.z << endl;
	return pout;
}
cad3d::Point3D cad3d::find_un_connected_non_planar_points_v2(stl_data &block_a, vector<Point3D> &centered_block, Point3D pin, int & status) {
	int ii;
	Point3D p1, p2, p3;
	float dd1, dd2, dd3;
	bool a, b, c;
	Point3D pout;
	pout = copy_point3d(pin);
	for (int i = 0; i < centered_block.size(); i++) {

		p1 = centered_block[i];
		dd1 = sqrt(dist201(p1, pin));
		if (dd1 < 1) {
			if ((pin.z - p1.z) > 0) {
				//std::cout << i << ": " << p1.z << endl;
				pout = copy_point3d(p1);
				status = 1;
				break;

			}

		}
	}
	//std::cout << pout.x << pout.y << pout.z << endl;
	return pout;
}

void cad3d::generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) {
	//one side
	Point3D ptemp1,ptemp3;//for one side and the other
	int non_planar_stat = 1;
	while (non_planar_stat==1) {
		non_planar_stat = find_lower_point(block_a, slot_plane1, ti_all2, ptemp1, pin);
		if (non_planar_stat == 1) {
			pin = copy_point3d(ptemp1);
			single_cs.push_back(ptemp1);
		}
		else
			break;
	}
	//ptemp3 = find_correlated_point(slot_e1, slot_e3, slot_ip1, single_cross);
	
	
}
void cad3d::generate_single_cross_along_centreline_v101(stl_data &block_a, vector<Point3D> &centered_block, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) {
	//one side
	Point3D ptemp1, ptemp3;//for one side and the other
	int non_planar_stat = 1;
	while (non_planar_stat == 1) {
		non_planar_stat = find_lower_point_v2(block_a, centered_block, slot_plane1, ti_all2, ptemp1, pin);
		if (non_planar_stat == 1) {
			pin = copy_point3d(ptemp1);
			single_cs.push_back(ptemp1);
		}
		else
			break;
	}
	//ptemp3 = find_correlated_point(slot_e1, slot_e3, slot_ip1, single_cross);


}
int cad3d::find_lower_point(stl_data &block_a, vector<Point3D> &slot_plane1, vector<triangle_index> &ti_all2,  Point3D  &pout, Point3D  pin) {
	Point3D ptemp1;//for one side and the other
	int non_planar_point_stat = 0;
	ptemp1 = find_non_planar_connected_points_v2(block_a, slot_plane1, ti_all2, pin, non_planar_point_stat);
	if (non_planar_point_stat == 1) {
		pout = copy_point3d(ptemp1);
		//std::cout << ptemp1.x << ", " << ptemp1.y;
		//std::cout << ", " << ptemp1.z << endl;
		non_planar_point_stat = 1;
	}
	else {
		ptemp1 = find_un_connected_non_planar_points(block_a, pin, non_planar_point_stat);
		if (non_planar_point_stat == 1) {
			pout = copy_point3d(ptemp1);
			//std::cout << ptemp1.x << ", " << ptemp1.y;
			//std::cout << ", " << ptemp1.z << endl;
			//target_plot3.push_back(ptemp3);
			non_planar_point_stat = 1;
		}
	}
	return non_planar_point_stat;
}
int cad3d::find_lower_point_v2(stl_data &block_a, vector<Point3D> &centered_block, vector<Point3D> &slot_plane1, vector<triangle_index> &ti_all2, Point3D  &pout, Point3D  pin) {
	Point3D ptemp1;//for one side and the other
	int non_planar_point_stat = 0;
	ptemp1 = find_non_planar_connected_points_v3(block_a, slot_plane1, ti_all2, pin, non_planar_point_stat);
	if (non_planar_point_stat == 1) {
		pout = copy_point3d(ptemp1);
		//std::cout << ptemp1.x << ", " << ptemp1.y;
		//std::cout << ", " << ptemp1.z << endl;
		non_planar_point_stat = 1;
	}
	else {
		ptemp1 = find_un_connected_non_planar_points_v2(block_a, centered_block, pin, non_planar_point_stat);
		if (non_planar_point_stat == 1) {
			pout = copy_point3d(ptemp1);
			//std::cout << ptemp1.x << ", " << ptemp1.y;
			//std::cout << ", " << ptemp1.z << endl;
			//target_plot3.push_back(ptemp3);
			non_planar_point_stat = 1;
		}
	}
	return non_planar_point_stat;
}
int cad3d::calculate_intersectline_from_interpolated_curves_v1(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, int i_present) {
	vector<float> dmat;
	float dtemp;
	for (int i = 0; i < interpolated_edge_v2.size(); i++) {
		dtemp = dist301(interpolated_edge_v1[i_present],interpolated_edge_v2[i]);
		dmat.push_back(dtemp);
	}
	int i_min_dtemp = 0;
	i_min_dtemp = find_1d_vector_min(dmat, 1);
	std::cout << dmat[i_min_dtemp]<<endl;

	return i_min_dtemp;
}


void cad3d::sub_volume_generate_v2(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &sub_vol) {
	// sub volume generate from single border
	cad3d::Point3D ptemp3;
	int i_temp;
	vector<cad3d::Point3D>   cross_line1, cross_line2;
	int l1, l3;
	l1 = v_in.size();
	//l3 = interpolated_edge_v2.size();


	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
	cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);
	/*newblock.order_cartesian_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2);
	std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
	int ixmax, iymax, ixmin, iymin;
	newblock.cartesian_separate_internal_regions_v2(feat1_e1, feat1_e2, v_int_l1, v_int_l2, 0, slot_plane1);*/

	//get depth data
	int non_planar_point_stat = 0;
	for (int i_vt = 0; i_vt < l1; i_vt++) {
		//for (int i_vt = 1; i_vt < 7; i_vt++) {

			//get depth data
		non_planar_point_stat = 0;
		sub_vol.push_back(v_in[i_vt]);
		generate_single_cross_along_centreline_v1(block_a, slot_plane1, v_in[i_vt], non_planar_point_stat, sub_vol, ti_all2, v_in, i_vt);
		//void cad3d::generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) 

	}
}
void cad3d::sub_volume_generate_v1(stl_data &block_a, vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &sub_vol) {

		cad3d::Point3D ptemp3;
		int i_temp;
		vector<cad3d::Point3D>   cross_line1, cross_line2;
		int l1, l3;
		l1 = interpolated_edge_v1.size();
		l3 = interpolated_edge_v2.size();


		vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
		vector<int> i_line;
		vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
		vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
		cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
		//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
		triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);
		/*newblock.order_cartesian_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2);
		std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
		int ixmax, iymax, ixmin, iymin;
		newblock.cartesian_separate_internal_regions_v2(feat1_e1, feat1_e2, v_int_l1, v_int_l2, 0, slot_plane1);*/

		//get intersect line
		//get depth data
		int non_planar_point_stat = 0;
		for (int i_vt = 0; i_vt < l1; i_vt++) {
			//for (int i_vt = 1; i_vt < 7; i_vt++) {

			//get intersect line
			cross_line1.push_back(interpolated_edge_v1[i_vt]);
			i_temp = calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, i_vt);
			cross_line2.push_back(interpolated_edge_v2[i_temp]);

			//get depth data
			//1st side
			non_planar_point_stat = 0;
			generate_single_cross_along_centreline_v1(block_a, slot_plane1, interpolated_edge_v1[i_vt], non_planar_point_stat, sub_vol, ti_all2, interpolated_edge_v1, i_vt);
			//void cad3d::generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) 
			//2nd side
			non_planar_point_stat = 0;
			generate_single_cross_along_centreline_v1(block_a, slot_plane1, interpolated_edge_v2[i_temp], non_planar_point_stat, sub_vol, ti_all2, interpolated_edge_v2, i_temp);
		}
	//i_temp = calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, single_cross);
	//cross_line2.push_back(interpolated_edge_v2[i_temp]);
	
}
void cad3d::sub_volume_generate_v3(stl_data &block_a, vector<Point3D> &centred_block, vector<Point3D> &v_in, vector<Point3D> &sub_vol) {
	// sub volume generate from single border
	cad3d::Point3D ptemp3;
	int i_temp;
	vector<cad3d::Point3D>   cross_line1, cross_line2;
	int l1, l3;
	l1 = v_in.size();
	//l3 = interpolated_edge_v2.size();


	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
	/*cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);*/

	cross_model_generate_polar(block_a, centred_block, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(centred_block, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);


	//get intersect line
	//get depth data
	int non_planar_point_stat = 0;
	for (int i_vt = 0; i_vt < l1; i_vt++) {
		//for (int i_vt = 1; i_vt < 7; i_vt++) {

		//get depth data
		non_planar_point_stat = 0;
		sub_vol.push_back(v_in[i_vt]);
		generate_single_cross_along_centreline_v101(block_a, centred_block, slot_plane1, v_in[i_vt], non_planar_point_stat, sub_vol, ti_all2, v_in, i_vt);
		//void cad3d::generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) 

	}
}
void cad3d::sub_volume_generate_v4(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in,vector<Point3D> &sub_vol) {
	//last version before October, 2019
	// sub volume generate from single border
	cad3d::Point3D ptemp3;
	int i_temp;
	
	int l1, l3;
	l1 = v_in.size();
	//l3 = interpolated_edge_v2.size();


	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
	cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);

	//get intersect line
	//get depth data
	int non_planar_point_stat = 0;

	std::cout << "sub volume generation v4: downwards" << endl;
	vector<float> sub_vol_z_levels;
	// finding how many levels the top border is connected with
	//initialize sub vol z levels
	sub_vol_z_levels.push_back(l1_in[0].z);
	if (l1_in[0].z != l2_in[0].z) {
		sub_vol_z_levels.push_back(l2_in[0].z);
		std::cout << "l1 z 0" << l1_in[0].z << "~";
		std::cout << "l2 z 0" << l2_in[0].z << endl;
	}
	
	int repeat_stat = 0;
	bool a = 0;
	bool b = 0;
	vector<cad3d::Point3D>  next_slice;
	for (int i_vt = 1; i_vt < l1_in.size(); i_vt++) {
		for (int j = 0; j < sub_vol_z_levels.size(); j++) {
			if (l1_in[i_vt].z == sub_vol_z_levels[j]) a = 1;
			if (l2_in[i_vt].z == sub_vol_z_levels[j]) b = 1;
			if (a || b == 1) repeat_stat = repeat_stat + 1;
			if (repeat_stat == 0) {
				if (a == 0) sub_vol_z_levels.push_back(l1_in[i_vt].z);
				if (b == 0) sub_vol_z_levels.push_back(l2_in[i_vt].z);
			}
		}
		/*if (repeat_stat == 0) {
			if(a==0) sub_vol_z_levels.push_back(l1_in[0].z);
			if (b == 0) sub_vol_z_levels.push_back(l2_in[0].z);
		}*/
		/*if (repeat_stat == 0) {
			if (a == 0) sub_vol_z_levels.push_back(l1_in[i_vt].z);
			if (b == 0) sub_vol_z_levels.push_back(l2_in[i_vt].z);
		}*/
		repeat_stat = 0;
		a = 0;
		b = 0;
		std::cout << "l1 z i" << l1_in[i_vt].z << "~";
		std::cout << "l2 z i" << l2_in[i_vt].z << endl;
	}
	std::cout << sub_vol_z_levels.size() << endl;
	std::cout << "connected z levels: " << endl;
	for (int j = 0; j < sub_vol_z_levels.size(); j++) {
		std::cout << "" << sub_vol_z_levels[j] << endl;
	}
	
	//find the points connected to this these points or these level
	// transfer to these points
	float z_current;
	std::cout << "check for connectivity level:" << endl;
	std::cin >> z_current;
	// check for connected points at the same level
	repeat_stat = 0;
	int c = 0;
	int d = 0;
	int e = 0;
	int f = 0;
	for (int j = 0; j < feat01_e1.size(); j++) {
		//if (l1_in[j].z == z_current) {
		//if ((l1_in[j].z <= z_current+0.5) & (l1_in[j].z <= z_current + 0.5)) {
		if ((feat01_e1[j].z <= z_current + 0.5) & (feat01_e1[j].z <= z_current + 0.5)) {
			//l1_in[j]
			//c = point_repeat_check_v1(l1_in[j], feat01_e1);
			//d = point_repeat_check_v1(l1_in[j], feat01_e2);
			sub_vol.push_back(feat01_e1[j]);
		}
		if ((feat01_e2[j].z <= z_current + 0.5) & (feat01_e2[j].z <= z_current + 0.5)) {
			//l1_in[j]
			//c = point_repeat_check_v1(l1_in[j], feat01_e1);
			//d = point_repeat_check_v1(l1_in[j], feat01_e2);
			sub_vol.push_back(feat01_e2[j]);
		}
		
	}
	/*//get depth data
	int non_planar_point_stat = 0;
	for (int i_vt = 0; i_vt < l1; i_vt++) {
		//for (int i_vt = 1; i_vt < 7; i_vt++) {

			//get depth data
		non_planar_point_stat = 0;
		sub_vol.push_back(v_in[i_vt]);
		generate_single_cross_along_centreline_v1(block_a, slot_plane1, v_in[i_vt], non_planar_point_stat, sub_vol, ti_all2, v_in, i_vt);
		//void cad3d::generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross) 

	}*/
}
int cad3d::sub_volume_generate_v5(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out) {
	//first version before October, 2019
	//updated on 11.10.2019
	//updated on 13.10.2019
	//updated on 22.10.2019
	//updated on 18.11.2019

	// sub volume generate from single border
	int slice_stat = 0;
	cad3d::Point3D ptemp3;
	int i_temp;

	int l1, l3;
	l1 = v_in.size();
	//l3 = interpolated_edge_v2.size();


	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
	cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);

	//get intersect line
	//get depth data
	int non_planar_point_stat = 0;

	std::cout << "sub volume generation v5: downwards" << endl;
	vector<float> sub_vol_z_levels;
	int z_level_status;
	//get depth data
	// finding how many levels the top border is connected with
	//initialize sub vol z levels
	
	//sub_volume_z_level_calculate_v1(v_in, feat01_e1, feat01_e2, sub_vol_z_levels);
	z_level_status = sub_volume_z_level_calculate_v2(v_in, feat01_e1, feat01_e2, sub_vol_z_levels);
	// need a function for knowing z levels

	int level_dif_show = 0;
	std::cout << "Enter 1 to show l1 and l2 z levels";
	std::cin >> level_dif_show;
	if (level_dif_show > 0) {
		for (int i = 0; i < l1_in.size(); i++) {
			std::cout << i << ": " << l1_in[i].z << ", " << l2_in[i].z << endl;
		}
	}
	//will run in branch-loop from here
	if (z_level_status>0) {
		std::cout << "a connected lower plane found" << endl;
		get_next_sub_slice_v2(v_in, feat01_e1, feat01_e2, sub_vol, l1_out, l2_out, 8);
	}
	if (z_level_status == 0) {
		std::cout << "no connected points on a lower plane" << endl;
		//get_same_z_sub_slice_v1(v_in, feat01_e1, feat01_e2, sub_vol, l1_out, l2_out);
		//get_same_z_sub_slice_v4(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
		std::cout << "input same level points size:" << v_in.size() << endl;
		//get_same_z_sub_slice_v3(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
		get_same_z_sub_slice_v5(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
	}
	bool a = 0;
	bool b = 0;
	bool c = 0;

	if (sub_vol.size() > 0) a = 1;
	if (l1_out.size() > 0) b = 1;
	if (l2_out.size() > 0) c = 1;
	
	if (a == 1) {
		if (b == 1 & c == 1) {
			slice_stat = 2;
		}
		else
			slice_stat = 1;
	}
	
	/**/
	return slice_stat;
}
void cad3d::sub_volume_generate_v6(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out, vector<vector<cad3d::Point3D >> &sub_vol_main_slices) {
	//last version 9 October, 2019
	// updated on 13.10.2019
	//updated on 22.10.2019
	//updated on 11.09.2019
	//updated on 3.12.2019

	// sub volume generate from single border in multiple layers
	//this will run in a loop unless a status 0 is received 
	// next work: include status variables for better control of the loop 
	//include line points in get_same_z_sub_slice_v2
	//include withon dialated boundary pointsboundary 
	vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	int i = 0;
	update_point3d_vector(v_in, v_last_slice);//update_point3d_vector(vector<Point3D> v_in, vector<Point3D> &v_out);
	update_point3d_vector(l1_in, l1_last_slice);
	update_point3d_vector(l2_in, l2_last_slice);
	sub_vol_main_slices.push_back(v_last_slice);
	//updating outputs
	update_point3d_vector(v_in, sub_vol);
	update_point3d_vector(l1_in, l1_out);
	update_point3d_vector(l2_in, l2_out);
	int slice_stat = 0;
	int user_next_level_stat = 1;
	while (i < 10) {
		
		std::cout << "enter positive integers if you want to extract a deeper layer: ";
		std::cin >> user_next_level_stat;
		if (user_next_level_stat < 1) break;
		std::cout << "sub_volume iteration step: " << i << endl;
		slice_stat = sub_volume_generate_v5(block_a, v_last_slice, l1_last_slice, l2_last_slice, v_next_slice, l1_next_slice, l2_next_slice);
		std::cout << "slice status: " << slice_stat<< endl;
		if (slice_stat < 1) break;
		//updating outputs
		update_point3d_vector(v_next_slice, sub_vol);
		update_point3d_vector(l1_next_slice, l1_out);
		update_point3d_vector(l2_next_slice, l2_out);
		sub_vol_main_slices.push_back(v_next_slice);

		//updating next iternation inputs
		v_last_slice.resize(0);
		update_point3d_vector(v_next_slice, v_last_slice);
		l1_last_slice.resize(0);
		update_point3d_vector(l1_next_slice, l1_last_slice);
		l2_last_slice.resize(0);
		update_point3d_vector(l2_next_slice, l2_last_slice);
		i++;
		v_next_slice.resize(0);
		l1_next_slice.resize(0);
		l2_next_slice.resize(0);
	}
	std::cout << "slice levels found: " <<i+1 << endl;
}
int cad3d::sub_volume_generate_v5_polar(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out) {
	//first version 29.11.19
	//updated on 4.12.19
	
	vector<Point3D> centered_block;
	centre_at_x0y0(block_a.vertices, centered_block);
	// sub volume generate from single border
	int slice_stat = 0;
	cad3d::Point3D ptemp3;
	int i_temp;

	int l1, l3;
	l1 = v_in.size();
	//l3 = interpolated_edge_v2.size();


	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  feat01_e1, feat01_e2, v_temp;
	triangle_vertexindex_init(block_a.triangles, ti_all2);
	triangle_normal_update(block_a.triangles, tsloti_feature);
	triangle_vertexindex_init(block_a.triangles, ti_all2);

	int ii;
	bool a = 0;
	bool b = 0;
	
	//cross_model_generate(block_a, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	cross_model_generate_polar(block_a, centered_block, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	//std::cout << "possible feature triangle no: " << tsloti_feature.size() << endl;
	triangle2line(centered_block, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);

	//get intersect line
	//get depth data
	int non_planar_point_stat = 0;

	std::cout << "sub volume generation v5: downwards" << endl;
	vector<float> sub_vol_z_levels;
	int z_level_status;
	//get depth data
	// finding how many levels the top border is connected with
	//initialize sub vol z levels

	//sub_volume_z_level_calculate_v1(v_in, feat01_e1, feat01_e2, sub_vol_z_levels);
	z_level_status = sub_volume_z_level_calculate_v2(v_in, feat01_e1, feat01_e2, sub_vol_z_levels);
	// need a function for knowing z levels

	int level_dif_show = 0;
	std::cout << "Enter 1 to show l1 and l2 z levels";
	std::cin >> level_dif_show;
	if (level_dif_show > 0) {
		for (int i = 0; i < l1_in.size(); i++) {
			std::cout << i << ": " << l1_in[i].z << ", " << l2_in[i].z << endl;
		}
	}
	//will run in branch-loop from here
	if (z_level_status>0) {
		std::cout << "a connected lower plane found" << endl;
		//get_next_sub_slice_v2(v_in, feat01_e1, feat01_e2, sub_vol, l1_out, l2_out, 10);//original 8
		get_next_sub_slice_v4(centered_block, slot_plane1, v_in, feat01_e1, feat01_e2, sub_vol, l1_out, l2_out, 10);//original 8
	}
	if (z_level_status == 0) {
		std::cout << "no connected points on a lower plane" << endl;
		//get_same_z_sub_slice_v1(v_in, feat01_e1, feat01_e2, sub_vol, l1_out, l2_out);
		//get_same_z_sub_slice_v4(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
		std::cout << "input same level points size:" << v_in.size() << endl;
		//get_same_z_sub_slice_v3(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
		get_same_z_sub_slice_v5(block_a, v_in, l1_in, l2_in, sub_vol, l1_out, l2_out);
	}
	bool c = 0;

	if (sub_vol.size() > 0) a = 1;
	if (l1_out.size() > 0) b = 1;
	if (l2_out.size() > 0) c = 1;

	if (a == 1) {
		if (b == 1 & c == 1) {
			slice_stat = 2;
		}
		else
			slice_stat = 1;
	}

	/**/
	return slice_stat;
}
void cad3d::sub_volume_generate_v7(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out, vector<vector<cad3d::Point3D >> &sub_vol_main_slices) {
	//updated on 11.09.2019
	//updated on 3.12.2019

	// sub volume generate from single border in multiple layers for polar models
	//this will run in a loop unless a status 0 is received 
	// next work: include status variables for better control of the loop 
	//include line points in get_same_z_sub_slice_v2
	//include withon dialated boundary pointsboundary 
	
	
	//vector<Point3D> centered_block;
	
	vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	int i = 0;
	update_point3d_vector(v_in, v_last_slice);//update_point3d_vector(vector<Point3D> v_in, vector<Point3D> &v_out);
	update_point3d_vector(l1_in, l1_last_slice);
	update_point3d_vector(l2_in, l2_last_slice);
	sub_vol_main_slices.push_back(v_last_slice);
	//updating outputs
	update_point3d_vector(v_in, sub_vol);
	update_point3d_vector(l1_in, l1_out);
	update_point3d_vector(l2_in, l2_out);
	int slice_stat = 0;
	int user_next_level_stat = 1;
	while (i < 10) {

		std::cout << "enter positive integers if you want to extract a deeper layer: ";
		std::cin >> user_next_level_stat;
		if (user_next_level_stat < 1) break;
		std::cout << "sub_volume iteration step: " << i << endl;
		slice_stat = sub_volume_generate_v5_polar(block_a, v_last_slice, l1_last_slice, l2_last_slice, v_next_slice, l1_next_slice, l2_next_slice);
		std::cout << "slice status: " << slice_stat << endl;
		if (slice_stat < 1) break;
		//updating outputs
		update_point3d_vector(v_next_slice, sub_vol);
		update_point3d_vector(l1_next_slice, l1_out);
		update_point3d_vector(l2_next_slice, l2_out);
		sub_vol_main_slices.push_back(v_next_slice);

		//updating next iternation inputs
		v_last_slice.resize(0);
		update_point3d_vector(v_next_slice, v_last_slice);
		l1_last_slice.resize(0);
		update_point3d_vector(l1_next_slice, l1_last_slice);
		l2_last_slice.resize(0);
		update_point3d_vector(l2_next_slice, l2_last_slice);
		i++;
		v_next_slice.resize(0);
		l1_next_slice.resize(0);
		l2_next_slice.resize(0);
	}
	//std::cout << "slice levels found: " << i + 1 << endl;
}

void cad3d::vertex_2_boundary_generate_v1(vector<Point3D> &v1, vector<Point3D> &v2, vector<Point3D> &f_boundary_line1, vector<Point3D> &f_boundary_line2, int cart_vs_pol, int slot_type, vector<Point3D> &v_feat_plane){
	//created on 03.12.2019
	//updated on 6.12.12
	//updated on 6.12.12
	//v1 => v_in
	vector<cad3d::Point3D> v1_red,v1_sel;
	int repeat_stat = 0;
	if (cart_vs_pol == 1) {
		if (slot_type == 1) {
			// semi open slot
			// created on 1.12.19
			//feature boundry generate fr open slot in Cartesian
			std::cout << "slot type semi open" << endl;
			for (int i = 0; i < v1.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				f_boundary_line1.push_back(v1[i]);
				f_boundary_line2.push_back(v1[i + 1]);
				//std::cout << i << ": " << v1[i].x << ", ";
				//std::cout << v1[i].y << ", " << v1[i].z << endl;
			}
		}
		if (slot_type == 2) {
			// created on 1.12.19
			// closed slot
			//// updated on 6.12.19
			//feature boundry generate fr closed slot in Cartesian
			//f_boundary_line1.push_back(v1[0]);
			//f_boundary_line2.push_back(v1[1]);
			for (int i = 0; i < v1.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				f_boundary_line1.push_back(v1[i]);
				f_boundary_line2.push_back(v1[i + 1]);
				//std::cout << i << ": " << v1[i].x << ", ";
				//std::cout << v1[i].y << ", " << v1[i].z << endl;
			}
			f_boundary_line1.push_back(v1[v1.size() - 1]);
			f_boundary_line2.push_back(v1[0]);
		}
		/*if (slot_type == 4) {
			// created on 6.12.19
			vector<cad3d::Point3D> f_line1, f_line2;
			std::wcout << "ordering boundaries" << endl;
			for (int i = 0; i < v1.size(); i++) {
				repeat_stat = point_repeat_check_v2(v1[i], v1_red);
				if (repeat_stat == 0) {
					v1_red.push_back(v1[i]);
				}
				repeat_stat = 0;
			}
			//from lines to points
			for (int i = 0; i < v1_red.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				//f_boundary_line1.push_back(v1_red[i]);
				//f_boundary_line2.push_back(v1_red[i + 1]);
				f_line1.push_back(v1_red[i]);
				f_line2.push_back(v1_red[i + 1]);
				//std::cout << i << ": " << v1_red[i].x << ", ";
				//std::cout << v1_red[i].y << ", " << v1_red[i].z << endl;
			}
			vector<int> rep_i_mat;
			int a = 0;
			int b = 0;
			float d1, d2, d3,d4;
			// check if they correspond with the feature plane boundary
			for (int i = 0; i < f_line1.size(); i++) {
				//repeat_stat = find_line_repete_index(f_line1[i],f_line2[i],ref_line1, ref_line2,3,rep_i_mat);
				for (int j = 0; j < ref_line1.size(); j++) {
					d1 = dist201(f_line1[i], ref_line1[j]);
					d2 = dist201(f_line2[i], ref_line2[j]);
					if ((d1 == 0) & (d2 == 0)) {
						f_boundary_line1.push_back(f_line1[i]);
						f_boundary_line2.push_back(f_line2[i]);
					}
					d3 = dist201(f_line1[i], ref_line2[j]);
					d4 = dist201(f_line2[i], ref_line1[j]);
					if ((d3 == 0) & (d4 == 0)) {
						f_boundary_line1.push_back(f_line2[i]);
						f_boundary_line2.push_back(f_line1[i]);
					}

				}

			}
		}*/
		if (slot_type == 5) {
			// created on 6.12.19
			vector<cad3d::Point3D> f_line1, f_line2;
			std::wcout << "ordering boundaries" << endl;
			vector<int> rep_i_mat;
			int a = 0;
			int b = 0;
			float d1, d2, d3, d4;
			for (int i = 0; i < v1.size(); i++) {
				repeat_stat = point_repeat_check_v2(v1[i], v1_red);
				if (repeat_stat == 0) {
					v1_red.push_back(v1[i]);
				}
				repeat_stat = 0;
			}
			for (int i = 0; i < v_feat_plane.size(); i++) {
				for (int j = 0; j < v1_red.size(); j++) {
					d1 = sqrt(dist201(v_feat_plane[i], v1_red[j]));
					if (d1 < 0.05) {
						v1_sel.push_back(v1_red[j]);
					}
				}
			}
			//from points to lines
			for (int i = 0; i < v1_sel.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				f_boundary_line1.push_back(v1_sel[i]);
				f_boundary_line2.push_back(v1_sel[i + 1]);
				//f_line1.push_back(v1_red[i]);
				//f_line2.push_back(v1_red[i + 1]);
				//std::cout << i << ": " << v1_red[i].x << ", ";
				//std::cout << v1_red[i].y << ", " << v1_red[i].z << endl;
			}


		}
	}
	if (cart_vs_pol == 2) {
		if (slot_type == 3) {
			std::wcout << "open slot re order" << endl;
			// created on 1.12.19
			//feature boundry generate fr open slot in Polar
			for (int i = 0; i < v1.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				f_boundary_line1.push_back(v1[i]);
				f_boundary_line2.push_back(v1[i + 1]);
				//std::cout << i << ": " << v1[i].x << ", ";
				//std::cout << v1[i].y << ", " << v1[i].z << endl;
			}

			if (v2.size() > 0) {
				std::wcout << "v2 starting" << endl;
				f_boundary_line1.push_back(v1[v1.size() - 1]);
				f_boundary_line2.push_back(v2[v2.size() - 1]);
				for (int i = v2.size() - 1; i > 0; i--) {
					//for (int i = 0; i < v2.size() - 1; i++) {
						//for (int i = 1; i < 10; i++) {
						f_boundary_line1.push_back(v2[i]);
						f_boundary_line2.push_back(v2[i - 1]);
						//std::cout << i << ": " << v2[i].x << ", ";
						//std::cout << v2[i].y << ", " << v2[i].z << endl;
					}
				//f_boundary_line1.push_back(v1[v1.size() - 1]);
				//f_boundary_line2.push_back(v2[v2.size() - 1]);
			}
		}

		if (slot_type == 4) {
			// created on 4.12.19
			//post process for Cartesian semi open slot
			//feature boundry generate fr open slot in Polar
			std::wcout << "ordering boundaries" << endl;
			for (int i = 0; i < v1.size(); i++) {
				repeat_stat = point_repeat_check_v2(v1[i], v1_red);
				if (repeat_stat == 0) {
					v1_red.push_back(v1[i]);
				}
				repeat_stat = 0;
			}
			for (int i = 0; i < v1_red.size() - 1; i++) {
				//for (int i = 1; i < 10; i++) {
				f_boundary_line1.push_back(v1_red[i]);
				f_boundary_line2.push_back(v1_red[i + 1]);
				//std::cout << i << ": " << v1_red[i].x << ", ";
				//std::cout << v1_red[i].y << ", " << v1_red[i].z << endl;
			}
		}
	}

}
void cad3d::characterization_2d_4_closed_slot_v1(vector<Point3D> &v_border, vector<Point3D> &c_line, vector<Point3D> &cent1, vector<Point3D> &cent2,int c_line_count) {
	std::cout << "pick external boundary point indices - " << endl;
	int i_cl1, i_cl2;
	for (int i = 0; i < v_border.size(); i++) {
		std::cout << i << ": " << v_border[i].x << ", ";
		std::cout << v_border[i].y << ", " << v_border[i].z << endl;
	}
	std::cout << "i1: ";
	std::cin >> i_cl1;
	std::cout << "i2: ";
	std::cin >> i_cl2;
	//c_line_count = 10;
	c_line.push_back(v_border[i_cl1]);
	cad3d::Point3D ptemp;
	float del_x = (v_border[i_cl2].x - v_border[i_cl1].x) / c_line_count;
	float del_y = (v_border[i_cl2].y - v_border[i_cl1].y) / c_line_count;
	for (int i = 1; i < c_line_count - 2; i++) {
		ptemp.x = c_line[0].x + del_x*i;
		ptemp.y = c_line[0].y + del_y*i;
		ptemp.z = c_line[0].z;
		c_line.push_back(ptemp);
	}
	c_line.push_back(v_border[i_cl2]);

	if (c_line.size() > 0) {
		//newblock.view_scale_for_3d(block_a.vertices, c_line, point_show_3d, plot_win3_min_x, plot_win3_max_x, plot_win3_min_y, plot_win3_max_y);
		for (int i = 0; i < c_line.size() - 1; i++) {
			cent1.push_back(c_line[i]);
			cent2.push_back(c_line[i + 1]);
		}
	}
}

void cad3d::sub_volume_z_level_calculate_v1(vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<float> z_levels) {
	//3.10.19
	//only checking z levels
	int repeat_stat = 0;
	int a = 0;
	int b = 0;
	z_levels.push_back(v_in[0].z);
	//vector<cad3d::Point3D>  next_slice;
	for (int i = 0; i < l1_in.size(); i++) {
		for (int j = 0; j < z_levels.size(); j++) {
			if (l1_in[i].z == z_levels[j]) {
				a++;
			}
			if (l2_in[i].z == z_levels[j]) {
				b++;
			}
			/*if (a || b == 1) repeat_stat = repeat_stat + 1;
			if (repeat_stat == 0) {
				
				if (b == 0) sub_vol_z_levels.push_back(l2_in[i_vt].z);
			}*/
		}
		if (a == 0) z_levels.push_back(l1_in[i].z);
		if (b == 0) z_levels.push_back(l2_in[i].z);
		
		repeat_stat = 0;
		a = 0;
		b = 0;

	}
	std::cout << "z levels of input slice: " <<z_levels.size() << endl;
	std::cout << "connected z levels: " << endl;
	for (int j = 0; j < z_levels.size(); j++) {
		std::cout << j<<": " << z_levels[j] << endl;
	}
}
int cad3d::sub_volume_z_level_calculate_v2(vector<Point3D> &v_in, vector<Point3D> &l1_in, vector<Point3D> &l2_in, vector<float> z_levels) {
	int z_level_status = 0;
	//3.10.19
	//checking point repitations
	int repeat_stat = 0;
	int a = 0;
	int b = 0;
	int point_repeat1 = 0;
	int point_repeat2 = 0;
	z_levels.push_back(v_in[0].z);
	//vector<cad3d::Point3D>  temp_slice;
	//for (int i = 0; i < l1_in.size(); i++) {
	for (int i = 0; i < l1_in.size(); i++) {
		point_repeat1 = point_repeat_check_v2(l1_in[i], v_in);
		point_repeat2 = point_repeat_check_v2(l2_in[i], v_in);
		if ((point_repeat1>0) || (point_repeat2>0)) {
			for (int j = 0; j < z_levels.size(); j++) {
				if (l1_in[i].z == z_levels[j]) {
					a++;
				}
				if (l2_in[i].z == z_levels[j]) {
					b++;
				}
				/*if (a || b == 1) repeat_stat = repeat_stat + 1;
				if (repeat_stat == 0) {

				if (b == 0) sub_vol_z_levels.push_back(l2_in[i_vt].z);
				}*/
			}
			if (a == 0) z_levels.push_back(l1_in[i].z);
			if (b == 0) z_levels.push_back(l2_in[i].z);
		}

		repeat_stat = 0;
		a = 0;
		b = 0;
		point_repeat1 = 0;
		point_repeat2 = 0;

	}
	std::cout << "z levels of input slice: " << z_levels.size() << endl;
	std::cout << "connected z levels: " << endl;
	for (int j = 0; j < z_levels.size(); j++) {
		std::cout << j << ": " << z_levels[j] << endl;
		if (z_levels[j] < v_in[0].z) z_level_status++;
	}
	return z_level_status;
}
void cad3d::get_next_sub_slice_v1(vector<Point3D> &v_in,vector<Point3D> &feat01_e1, vector<Point3D> &feat01_e2,vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out) {
	// or finding lower z level slice value
	float z_current;
	std::cout << "next level z:" << endl;
	std::cin >> z_current;
	int a = 0;
	int b = 0;
	//to check if any point connected with v_in
	for (int i = 0; i < feat01_e1.size(); i++) {
		a = point_repeat_check_v2(feat01_e2[i], v_in);
		b = point_repeat_check_v2(feat01_e1[i], v_in);
		if ((a > 0) & (b == 0)) {
			if ((feat01_e1[i].z <= z_current + 0.5) & (feat01_e1[i].z >= z_current - 0.5)) {
				sub_vol.push_back(feat01_e1[i]);
				l1_out.push_back(feat01_e1[i]);
				l2_out.push_back(feat01_e2[i]);
			}
		}

		//if (b>0) {
		if ((a == 0) & (b > 0)) {
			if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
				sub_vol.push_back(feat01_e2[i]);
				l1_out.push_back(feat01_e2[i]);
				l2_out.push_back(feat01_e1[i]);
			}
		}
	}
}
void cad3d::get_next_sub_slice_v2(vector<Point3D> &v_in, vector<Point3D> &feat01_e1, vector<Point3D> &feat01_e2, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out, float tol) {
	//first version 22.10.2019
	//updated on 4.12.19

	// for finding lower z level slice value withiin a tolerance shadow level
	float z_current;
	std::cout << "next level z:" << endl;
	std::cin >> z_current;
	int a = 0;
	int b = 0;
	int shadow_a = 0;
	//int shadow_b = 0;
	//to check if any point connected with v_in
	vector<Point3D> feat_pol1, feat_pol2;
	convert_cart2pol(feat01_e1, feat_pol1);
	convert_cart2pol(feat01_e2, feat_pol2);
	for (int i = 0; i < feat01_e1.size(); i++) {
		//std::cout << i << ": " << feat01_e1[i].x << ", " << feat01_e1[i].y << "~";
		//std::cout << feat01_e2[i].x << ", " << feat01_e2[i].y << "~";
		a = point_repeat_check_v2(feat01_e2[i], v_in);
		shadow_a = test_if_in_shadow_v1(feat01_e2[i], v_in,tol);
		b = point_repeat_check_v2(feat01_e1[i], v_in);
		if (shadow_a>0) {
			if ((a > 0) & (b == 0)) {
				if ((feat01_e1[i].z <= z_current + 1.0) & (feat01_e1[i].z >= z_current - 1.0)) {
					sub_vol.push_back(feat01_e1[i]);
					l1_out.push_back(feat01_e1[i]);
					l2_out.push_back(feat01_e2[i]);
				}
				//std::cout << i << ": " << feat_pol1[i].x << ", " << feat_pol1[i].y << "~";
				//std::cout << feat_pol2[i].x << ", " << feat_pol2[i].y << "~";
				//std::cout << i << ": " <<feat_pol1[i].x - feat_pol2[i].x << "~";
				//std::cout << a << ", " << b << endl;
			}

			//if (b>0) {
			if ((a == 0) & (b > 0)) {
				if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
					sub_vol.push_back(feat01_e2[i]);
					l1_out.push_back(feat01_e2[i]);
					l2_out.push_back(feat01_e1[i]);
				}
				//std::cout << i << ": " << feat_pol1[i].x << ", " << feat_pol1[i].y << "~";
				//std::cout << feat_pol2[i].x << ", " << feat_pol2[i].y << "~";
				//std::cout << i << ": " << feat_pol1[i].x - feat_pol2[i].x << "~";
				//std::cout << a << ", " << b << endl;
			}
		}
		//std::cout << endl;
		a = 0;
		b = 0;
		shadow_a = 0;
	}
}
void cad3d::get_next_sub_slice_v3(vector<Point3D> &v_in, vector<Point3D> &feat01_e1, vector<Point3D> &feat01_e2, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out, float tol) {
	//first version 4.12.2019
	//updated on ...
	//for polar mode lower points

	// for finding lower z level slice value withiin a tolerance shadow level
	float z_current;
	std::cout << "next level z:" << endl;
	std::cin >> z_current;
	int a = 0;
	int b = 0;
	int c = 0;
	int d = 0;
	int shadow_1 = 0;
	int shadow_2 = 0;
	vector<cad3d::Point3D>  shadow_used;
	//int shadow_b = 0;
	//to check if any point connected with v_in
	vector<Point3D> feat_pol1, feat_pol2;
	convert_cart2pol(feat01_e1, feat_pol1);
	convert_cart2pol(feat01_e2, feat_pol2);
	for (int j = 0; j < v_in.size(); j++) {
		for (int i = 0; i < feat01_e1.size(); i++) {
			//std::cout << i << ": " << feat01_e1[i].x << ", " << feat01_e1[i].y << "~";
			//std::cout << feat01_e2[i].x << ", " << feat01_e2[i].y << "~";
			a = point_repeat_check_v2(v_in[j], shadow_used);
			//b = point_repeat_check_v2(feat01_e1[i], v_in);
			shadow_1 = test_if_in_shadow_v1(feat01_e1[i], v_in, 10);
			shadow_2 = test_if_in_shadow_v1(feat01_e2[i], v_in, 10);
			if (shadow_1 > 0 && shadow_2 > 0) {
				if ((feat01_e1[i].z <= z_current + 1.0) & (feat01_e1[i].z >= z_current - 1.0)) {
					sub_vol.push_back(feat01_e1[i]);
					if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
						l1_out.push_back(feat01_e1[i]);
						l2_out.push_back(feat01_e2[i]);
					}
				}
				if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
					sub_vol.push_back(feat01_e2[i]);
				}
			}
			/*shadow_a = test_if_in_shadow_v1(feat01_e2[i], v_in, tol);
			
			if (shadow_a>0) {
			if ((a > 0) & (b == 0)) {
			if ((feat01_e1[i].z <= z_current + 1.0) & (feat01_e1[i].z >= z_current - 1.0)) {
			sub_vol.push_back(feat01_e1[i]);
			l1_out.push_back(feat01_e1[i]);
			l2_out.push_back(feat01_e2[i]);
			}
			//std::cout << i << ": " << feat_pol1[i].x << ", " << feat_pol1[i].y << "~";
			//std::cout << feat_pol2[i].x << ", " << feat_pol2[i].y << "~";
			std::cout << i << ": " << feat_pol1[i].x - feat_pol2[i].x << "~";
			std::cout << a << ", " << b << endl;
			}

			//if (b>0) {
			if ((a == 0) & (b > 0)) {
			if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
			sub_vol.push_back(feat01_e2[i]);
			l1_out.push_back(feat01_e2[i]);
			l2_out.push_back(feat01_e1[i]);
			}
			//std::cout << i << ": " << feat_pol1[i].x << ", " << feat_pol1[i].y << "~";
			//std::cout << feat_pol2[i].x << ", " << feat_pol2[i].y << "~";
			std::cout << i << ": " << feat_pol1[i].x - feat_pol2[i].x << "~";
			std::cout << a << ", " << b << endl;
			}
			}
			//std::cout << endl;*/
			a = 0;
			b = 0;
			shadow_1 = 0;
			shadow_2 = 0;
		}
	}
	
}
void cad3d::get_next_sub_slice_v4(vector<Point3D> &centered_block, vector<Point3D> &slot_plane, vector<Point3D> &v_in, vector<Point3D> &feat01_e1, vector<Point3D> &feat01_e2, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out, float tol) {
	//first version 5.12.2019
	//updated on 6.12.2019
	//for polar mode lower points

	// for finding lower z level slice value withiin a tolerance shadow level
	float z_current;
	std::cout << "next level z:" << endl;
	std::cin >> z_current;
	int a = 0;
	int b = 0;
	int c = 0;
	int d = 0;
	int shadow_1 = 0;
	int shadow_2 = 0;
	float d11, d12;
	float d21, d22;
	vector<cad3d::Point3D>  shadow_used;
	//std::cout << "input plane points: " << feat01_e1.size() << ", " << feat01_e2.size() << endl;
	vector<cad3d::Point3D>  v1_int;
	vector<int> vi_temp551, vi_temp552;

	vector<cad3d::Point3D>  v1, v1_polar;
	vector<cad3d::Point3D>  feat1_e1_polar, feat1_e2_polar, feat01_e1_polar, feat01_e2_polar, slot_plane1_polar;
	vector<cad3d::Point3D>  feat1_e1, feat1_e2, v_int_l1, v_int_l2;
	convert_cart2pol(feat01_e1, feat01_e1_polar);
	convert_cart2pol(feat01_e2, feat01_e2_polar);
	order_polar_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2, feat01_e1_polar, feat01_e2_polar);
	//5b
	polar_separate_internal_regions_v1(feat1_e1, feat1_e2, v_int_l1, v_int_l2, 0, slot_plane);
	std::cout << "v_int_l1 points: " << v_int_l1.size() << ", " << v_int_l1.size() << endl;
	//update_point3d_vector(feat01_e1, sub_vol);
	//int shadow_b = 0;
	//to check if any point connected with v_in
	
	//update_point3d_vector(open_slot_pts_temp, sub_vol);
	/*for (int j = 0; j < v_in.size;j++) {
	//for (int j = 0; j < 13; j++) {
		for (int i = 0; i < feat01_e1.size(); i++) {
			d11 = sqrt(dist201(v_in[j], feat01_e1[i]));
			d12 = sqrt(dist201(v_in[j], feat01_e2[i]));
			if (d11 < 1) {
				v1.push_back(feat01_e1[i]);
				v1.push_back(feat01_e2[i]);
				if ((feat01_e1[i].z <= z_current + 1.0) & (feat01_e1[i].z >= z_current - 1.0)) {
					
					if (sub_vol.size()> 0) {
						d21 = sqrt(dist201(sub_vol[sub_vol.size()-1], feat01_e1[i]));
						std::cout << j << ": " << d21 << endl;
					}
					
					v1_int.push_back(feat01_e1[i]);
				}
				if ((feat01_e2[i].z <= z_current + 1.0) & (feat01_e2[i].z >= z_current - 1.0)) {
					
					if (sub_vol.size()> 0) {
						d21 = sqrt(dist201(sub_vol[sub_vol.size() - 1], feat01_e2[i]));
						std::cout << j << ": " << d21 << endl;
					}
				}
						
			}
			if (d12 < 1) {
				v1.push_back(feat01_e1[i]);
				v1.push_back(feat01_e2[i]);
				if ((feat01_e1[i].z <= z_current + 1.0) & (feat01_e1[i].z >= z_current - 1.0)) {
					
					if (sub_vol.size()> 0) {
						d21 = sqrt(dist201(sub_vol[sub_vol.size() - 1], feat01_e1[i]));
						std::cout << j << ": " << d21 << endl;
					}
					v1_int.push_back(feat01_e1[i]);
				}
				if ((feat01_e2[i].z <= z_current + 1.0) & (feat01_e2[i].z >= z_current - 1.0)) {
					
					if (sub_vol.size()> 0) {
						d21 = sqrt(dist201(sub_vol[sub_vol.size() - 1], feat01_e2[i]));
						std::cout << j << ": " << d21 << endl;
					}
				
			}
		}
		
	}*/
	std::cout << "teasting shadow..." <<  endl;
	//convert_cart2pol(v1, v1_polar);
	//sub_vol.push_back(v1[0]);
	//int k = 1;// suze of new sub_vol
	//int dummy = 0;
	/*for (int i = 1; i < v1.size() - 1; i++) {
		
		if (dist201(v1[i], v1[i-1])!=0) {
			//if (abs(v1_polar[i].z - v1_polar[i - 1].z)<9) {
			if((v1[i].z <= z_current + 1.0) & (v1[i].z >= z_current - 1.0)){
				//if (abs(v1_polar[i].y - v1_polar[i - 1].y) < 25) {
					sub_vol.push_back(v1[i]);
					std::cout << i << ", " << v1_polar[i].x - v1_polar[i - 1].x << ", ";
					std::cout << v1_polar[i].y - v1_polar[i - 1].y << ", ";
					std::cout << v1_polar[i].z - v1_polar[i - 1].z << endl;
					k++;
				//}
			
			}
			
		}
	}*/
	//update_point3d_vector(v1, sub_vol);
	vector<cad3d::Point3D>  v1_pseudo, v_plane;
	update_point3d_vector(v_in, v1_pseudo);
	for (int i = 0; i < v_in.size(); i++) {
		v1_pseudo[i].z = z_current;
	}
	find_point_3d_1feature(v_plane,centered_block, z_current, 1,3);
	int rep_stat = 0;;
	for (int i = 0; i < v1_pseudo.size(); i++) {
		//d21 = dist201()
		rep_stat = point_repeat_check_v2(v1_pseudo[i], v_plane);
		if (rep_stat == 1) {
			sub_vol.push_back(v1_pseudo[i]);
		}
		rep_stat = 0;
	}
	std::cout << "next level found: " << sub_vol.size() << endl;
}

void cad3d::get_same_z_sub_slice_v1(vector<Point3D> &v_in, vector<Point3D> &feat01_e1, vector<Point3D> &feat01_e2, vector<Point3D> &sub_vol, vector<Point3D> &l1_out, vector<Point3D> &l2_out) {
	//same level sub slice
	//3.10.19
	float z_current;
	std::cout << "same level z:" << endl;
	std::cin >> z_current;
	int a = 0;
	int b = 0;
	float d;
	//to check if any point connected with v_in
	for (int i = 0; i < feat01_e1.size(); i++) {
		
		a = point_repeat_check_v2(feat01_e2[i], v_in);
		b = point_repeat_check_v2(feat01_e1[i], v_in);
		if ((a > 0) & (b == 0)) {
			//std::cout << feat01_e2[i].x << ", " << feat01_e2[i].y << "~";
			//std::cout << feat01_e1[i].x << ", " << feat01_e1[i].y << endl;
			d = dist301(feat01_e2[i], feat01_e1[i]);
			//if ((feat01_e1[i].z <= z_current + 0.5) & (feat01_e1[i].z >= z_current - 0.5)) {
			if (feat01_e1[i].z <= z_current + 0.5)  {
				std::cout << d << endl;
				sub_vol.push_back(feat01_e1[i]);
				l1_out.push_back(feat01_e2[i]);
				l2_out.push_back(feat01_e1[i]);
			}
		}

		//if (b>0) {
		if ((a == 0) & (b > 0)) {
			//std::cout << feat01_e1[i].x << ", " << feat01_e1[i].y << "~";
			//std::cout << feat01_e2[i].x << ", " << feat01_e2[i].y << endl;
			d = dist301(feat01_e2[i], feat01_e1[i]);
			//if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
			if (feat01_e2[i].z <= z_current + 0.5){
				std::cout << d << endl;
				sub_vol.push_back(feat01_e2[i]);
				l1_out.push_back(feat01_e2[i]);
				l2_out.push_back(feat01_e1[i]);
			}
		}
	}
}
void cad3d::get_same_z_sub_slice_v2(stl_data &block_a, vector<cad3d::Point3D> v_last_slice, vector<Point3D> l1_last_slice, vector<Point3D> l2_last_slic, vector<cad3d::Point3D> &v_next_slice) {
	//same level sub slice
	//9.10.19

	vector<cad3d::Point3D> slot_plane1;
	vector<cad3d::Point3D>   sub_vol_ext, v1_next_slice, slice1_01, slice1_02;
	vector<cad3d::Point3D>   slice2_01, slice2_02;
	//vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	//vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	vector<int> i_line;
	
	//newblock.sub_volume_generate_v4(block_a, v1, slot_e01, slot_e02,sub_vol_ext); 
	vector<cad3d::Point3D>  line_e1, line_e2, vv_temp, slot_plane_t;
	vector<cad3d::triangle_index> tsi_feature, ts_i1, ts_i2, ti_temp2;
	cross_model_generate_v2(block_a, slot_plane_t, ti_temp2, ts_i1, ts_i2, tsi_feature, 15);
	// gets point depending on atolerance level, for getting multiple z slice levels
	i_line.resize(0);
	
	triangle2line(block_a.vertices, tsi_feature, line_e1, line_e2, i_line, vv_temp);
	
	float z_current;
	std::cout << "current level z:" << endl;
	std::cin >> z_current;

	// gets all the lines on  and connected to the same plane
	//step 1: get all points
	int a = 0;	int b = 0;
	int c = 0;	int d = 0;
	int e = 0;	int f = 0;
	int p = 0;	int q = 0;
	for(int i = 0; i < line_e1.size(); i++) {
		a = point_repeat_check_v2(line_e1[i], v_last_slice);
		b = point_repeat_check_v2(line_e2[i], v_last_slice);
		p = point_repeat_check_v2(line_e1[i], v_next_slice);
		q = point_repeat_check_v2(line_e2[i], v_next_slice);
		if (!a || !b) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				if (!p & !q) {
					v_next_slice.push_back(line_e1[i]);
					v_next_slice.push_back(line_e2[i]);
				}
			}
		}
		a = 0;	b = 0;
		c = 0;	d = 0;
		p = 0;	q = 0;
	}

	// need to keep points only on the same plane and not repeated
	
}
void cad3d::get_same_z_sub_slice_v3(stl_data &block_a, vector<cad3d::Point3D> v_last_slice, vector<Point3D> l1_last_slice, vector<Point3D> l2_last_slic, vector<cad3d::Point3D> &v_out_slice, vector<Point3D> l1_next_slice, vector<Point3D> l2_next_slice) {
	//same level sub slice
	// has troubles getting the lines on the same plane 
	//11.10.19
	// updated: 13.10.19
	// updated: 13.11.19

	vector<cad3d::Point3D> slot_plane1;
	vector<cad3d::Point3D>   sub_vol_ext, v1_next_slice, slice1_01, slice1_02;
	vector<cad3d::Point3D>   slice2_01, slice2_02;
	//vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	//vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	vector<int> i_line;

	vector<cad3d::Point3D> v_next_slice;

	//newblock.sub_volume_generate_v4(block_a, v1, slot_e01, slot_e02,sub_vol_ext); 
	vector<cad3d::Point3D>  line_e1, line_e2, vv_temp, slot_plane_t;
	vector<cad3d::triangle_index> tsi_feature, ts_i1, ts_i2, ti_temp2;
	cross_model_generate_v2(block_a, slot_plane_t, ti_temp2, ts_i1, ts_i2, tsi_feature, 15);
	// gets point depending on atolerance level, for getting multiple z slice levels
	i_line.resize(0);

	triangle2line(block_a.vertices, tsi_feature, line_e1, line_e2, i_line, vv_temp);

	float z_current;
	std::cout << "current level z:" << endl;
	std::cin >> z_current;

	// gets all the lines on  and connected to the same plane
	//step 1: get all points
	int a = 0;	int b = 0;
	int c = 0;	int d = 0;
	int e = 0;	int f = 0;
	int p = 0;	int q = 0;
	//shadow_a = test_if_in_shadow_v1(feat01_e2[i], v_in,tol);
	for (int i = 0; i < line_e1.size(); i++) {
		a = point_repeat_check_v2(line_e1[i], v_last_slice);
		b = point_repeat_check_v2(line_e2[i], v_last_slice);
		p = point_repeat_check_v2(line_e1[i], v_next_slice);
		q = point_repeat_check_v2(line_e2[i], v_next_slice);
		/*if (!a || !b) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				if (!p & !q) {
					v_next_slice.push_back(line_e1[i]);
					v_next_slice.push_back(line_e2[i]);
				}
			}
		}*/
		if ((a > 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
			//if (c || d) {
				//if (!p & !q) {
				//if (!p) {
					v_next_slice.push_back(line_e1[i]);
					//v_next_slice.push_back(line_e2[i]);
					//l1_next_slice.push_back(line_e1[i]);
					//l2_next_slice.push_back(line_e2[i]);
				}
			//}
		}
		if ((a == 0) & (b > 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
			//if (c || d) {
				//if (!p & !q) {
				//if (!q) {
					v_next_slice.push_back(line_e2[i]);
					//v_next_slice.push_back(line_e2[i]);
					//l1_next_slice.push_back(line_e1[i]);
					//l2_next_slice.push_back(line_e2[i]);
				//}
			}
		}
		if ((a == 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
			//if (c || d) {
				//if (!p & !q) {
				//if (!q) {
					v_next_slice.push_back(line_e2[i]);
					v_next_slice.push_back(line_e2[i]);
					//l1_next_slice.push_back(line_e1[i]);
					//l2_next_slice.push_back(line_e2[i]);
				//}
			}
		}
		a = 0;	b = 0;
		c = 0;	d = 0;
		p = 0;	q = 0;
	}
	// for points to line
	int m = 0;	int n = 0;
	/*for (int i = 0; i < v_next_slice.size(); i++) {
		m = point_repeat_check_v2(line_e1[i], v_last_slice);
		//n = point_repeat_check_v2(line_e2[i], v_next_slice);
		if ((m ==0){
			ppo
		}
		m = 0;
		//n = 0;
	}*/
	int superposition_check = 0;
	Point3D p_pop;
	//std:cout << "last slice size: " << v_next_slice.size() << endl;
	std::cout << "input boundary border size [xmin, xmax, ymin, ymax]: ";
	vector<cad3d::Point3D>  bound_box, v_test;
	int i_temp = 0;
	i_temp = find_point3d_min(v_last_slice, 1);
	std::cout << v_last_slice[i_temp].x << ", " ;
	bound_box.push_back(v_last_slice[i_temp]);
	i_temp = find_point3d_max(v_last_slice, 1);
	std::cout << v_last_slice[i_temp].x << ", "  ;
	i_temp = find_point3d_min(v_last_slice, 2);
	std::cout << v_last_slice[i_temp].y << ", "  ;
	i_temp = find_point3d_max(v_last_slice, 2);
	std::cout << v_last_slice[i_temp].y << endl ;

	for (int i = 0; i < v_next_slice.size(); i++) {
		m = point_superposition_check_v1(v_next_slice[i], v_last_slice, 50.0);
		//m = point_repeat_check_v1(v_next_slice[i], v_last_slice);
		//n = point_repeat_check_v2(line_e2[i], v_next_slice);
		//if ((m > 0) || (n > 0)) {
		if(m > 0){
			v_out_slice.push_back(v_next_slice[i]);
			//l1_next_slice.push_back(line_e1[i]);
			//l2_next_slice.push_back(line_e2[i]);
		}
		m = 0;
	}
	std::cout << "super position test (in & out): " << v_next_slice.size() <<", " << v_out_slice.size()<< endl;
	/*
	

	a = point_repeat_check_v2(feat01_e2[i], v_in);
	b = point_repeat_check_v2(feat01_e1[i], v_in);
	

	//if (b>0) {
	if ((a == 0) & (b > 0)) {
	//std::cout << feat01_e1[i].x << ", " << feat01_e1[i].y << "~";
	//std::cout << feat01_e2[i].x << ", " << feat01_e2[i].y << endl;
	d = dist301(feat01_e2[i], feat01_e1[i]);
	//if ((feat01_e2[i].z <= z_current + 0.5) & (feat01_e2[i].z >= z_current - 0.5)) {
	if (feat01_e2[i].z <= z_current + 0.5){
	std::cout << d << endl;
	sub_vol.push_back(feat01_e2[i]);
	l1_out.push_back(feat01_e2[i]);
	l2_out.push_back(feat01_e1[i]);
	}
	}
	}
	*/
	// need to keep points only on the same plane and not repeated

}
void cad3d::get_same_z_sub_slice_v4(stl_data &block_a, vector<cad3d::Point3D> v_last_slice, vector<Point3D> l1_last_slice, vector<Point3D> l2_last_slic, vector<cad3d::Point3D> &v_out_slice, vector<Point3D> l1_next_slice, vector<Point3D> l2_next_slice) {
	//same level sub slice
	//has troubles getting the lines on the same plane 
	//written on 22.10.2019
	// updated: ...
	vector<cad3d::Point3D> slot_plane1;
	vector<cad3d::Point3D>   sub_vol_ext, v1_next_slice, slice1_01, slice1_02;
	vector<cad3d::Point3D>   slice2_01, slice2_02;
	//vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	//vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	vector<int> i_line;

	vector<cad3d::Point3D> v_next_slice;

	//newblock.sub_volume_generate_v4(block_a, v1, slot_e01, slot_e02,sub_vol_ext); 
	vector<cad3d::Point3D>  line_e1, line_e2, vv_temp, slot_plane_t;
	vector<cad3d::triangle_index> tsi_feature, ts_i1, ts_i2, ti_temp2;
	cross_model_generate_v2(block_a, slot_plane_t, ti_temp2, ts_i1, ts_i2, tsi_feature, 15);
	// gets point depending on atolerance level, for getting multiple z slice levels
	i_line.resize(0);

	triangle2line(block_a.vertices, tsi_feature, line_e1, line_e2, i_line, vv_temp);

	float z_current;
	std::cout << "current level z:" << endl;
	std::cin >> z_current;

	// gets all the lines on  and connected to the same plane
	//step 1: get all points
	int a = 0;	int b = 0;
	int c = 0;	int d = 0;
	int e = 0;	int f = 0;
	int p = 0;	int q = 0;
	int shadow_a = 0;
	int shadow_b = 0;
	//shadow_a = test_if_in_shadow_v1(feat01_e2[i], v_in,tol);
	std::cout << line_e1.size() << endl;
	for (int i = 0; i < line_e1.size(); i++) {
		a = point_repeat_check_v2(line_e1[i], v_last_slice);
		b = point_repeat_check_v2(line_e2[i], v_last_slice);
		shadow_a = test_if_in_shadow_v1(line_e1[i], v_last_slice, 10.0);
		shadow_b = test_if_in_shadow_v1(line_e2[i], v_last_slice, 10.0);
		
		p = point_repeat_check_v2(line_e1[i], v_next_slice);
		q = point_repeat_check_v2(line_e2[i], v_next_slice);
		if ((shadow_a > 0) || (shadow_b > 0)) {
			if (!a || !b) {
				std::cout << i << ": " << shadow_a << ", " << shadow_b << endl;
				v_next_slice.push_back(line_e1[i]);
			}
			
			/*if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) v_next_slice.push_back(line_e1[i]);*/
		}
		/*if ((a > 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				
				v_next_slice.push_back(line_e1[i]);
				
			}
			//}
		}
		if ((a == 0) & (b > 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				
				v_next_slice.push_back(line_e2[i]);
				
			}
		}
		if ((a == 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				
				v_next_slice.push_back(line_e2[i]);
				v_next_slice.push_back(line_e2[i]);
				
			}
		}*/
		a = 0;	b = 0;
		c = 0;	d = 0;
		p = 0;	q = 0;
		shadow_a = 0; shadow_b = 0;
	}
	// for points to line
	/*int m = 0;	int n = 0;
	
	int superposition_check = 0;
	Point3D p_pop;
	std:cout << "last slice size: " << v_next_slice.size() << endl;
	for (int i = 0; i < v_next_slice.size(); i++) {
		m = point_superposition_check_v1(v_next_slice[i], v_last_slice, 50.0);
		
		//if ((m > 0) || (n > 0)) {
		if (m > 0) {
			v_out_slice.push_back(v_next_slice[i]);
			//l1_next_slice.push_back(line_e1[i]);
			//l2_next_slice.push_back(line_e2[i]);
		}
		m = 0;
	}
	std::cout << "super position test (in & out): " << v_next_slice.size() << ", " << v_out_slice.size() << endl;
	*/
	// need to keep points only on the same plane and not repeated

}
void cad3d::get_same_z_sub_slice_v5(stl_data &block_a, vector<cad3d::Point3D> v_last_slice, vector<Point3D> l1_last_slice, vector<Point3D> l2_last_slic, vector<cad3d::Point3D> &v_out_slice, vector<Point3D> l1_next_slice, vector<Point3D> l2_next_slice) {
	//same level sub slice
	// has troubles getting the lines on the same plane 
	//18.11.19

	vector<cad3d::Point3D> slot_plane1;
	vector<cad3d::Point3D>   sub_vol_ext, v1_next_slice, slice1_01, slice1_02;
	vector<cad3d::Point3D>   slice2_01, slice2_02;
	//vector<cad3d::Point3D>   v_last_slice, l1_last_slice, l2_last_slice;
	//vector<cad3d::Point3D>   v_next_slice, l1_next_slice, l2_next_slice;
	vector<int> i_line;

	vector<cad3d::Point3D> v_next_slice;

	//newblock.sub_volume_generate_v4(block_a, v1, slot_e01, slot_e02,sub_vol_ext); 
	vector<cad3d::Point3D>  line_e1, line_e2, vv_temp, slot_plane_t;
	vector<cad3d::triangle_index> tsi_feature, ts_i1, ts_i2, ti_temp2;
	cross_model_generate_v2(block_a, slot_plane_t, ti_temp2, ts_i1, ts_i2, tsi_feature, 15);
	// gets point depending on atolerance level, for getting multiple z slice levels
	i_line.resize(0);

	triangle2line(block_a.vertices, tsi_feature, line_e1, line_e2, i_line, vv_temp);

	float z_current;
	std::cout << "current level z:" << endl;
	std::cin >> z_current;

	// gets all the lines on  and connected to the same plane
	//step 1: get all points
	int a = 0;	int b = 0;
	int c = 0;	int d = 0;
	int e = 0;	int f = 0;
	int p = 0;	int q = 0;
	//shadow_a = test_if_in_shadow_v1(feat01_e2[i], v_in,tol);
	for (int i = 0; i < line_e1.size(); i++) {
		a = point_repeat_check_v2(line_e1[i], v_last_slice);
		b = point_repeat_check_v2(line_e2[i], v_last_slice);
		p = point_repeat_check_v2(line_e1[i], v_next_slice);
		q = point_repeat_check_v2(line_e2[i], v_next_slice);
		
		if ((a > 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				//if (c || d) {
				//if (!p & !q) {
				//if (!p) {
				v_next_slice.push_back(line_e1[i]);
				//v_next_slice.push_back(line_e2[i]);
				//l1_next_slice.push_back(line_e1[i]);
				//l2_next_slice.push_back(line_e2[i]);
			}
			//}
		}
		if ((a == 0) & (b > 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				//if (c || d) {
				//if (!p & !q) {
				//if (!q) {
				v_next_slice.push_back(line_e2[i]);
				//v_next_slice.push_back(line_e2[i]);
				//l1_next_slice.push_back(line_e1[i]);
				//l2_next_slice.push_back(line_e2[i]);
				//}
			}
		}
		if ((a == 0) & (b == 0)) {
			if (line_e1[i].z == z_current) c = 1;
			if (line_e2[i].z == z_current) d = 1;
			if (c & d) {
				//if (c || d) {
				//if (!p & !q) {
				//if (!q) {
				v_next_slice.push_back(line_e2[i]);
				v_next_slice.push_back(line_e2[i]);
				//l1_next_slice.push_back(line_e1[i]);
				//l2_next_slice.push_back(line_e2[i]);
				//}
			}
		}
		a = 0;	b = 0;
		c = 0;	d = 0;
		p = 0;	q = 0;
	}
	// for points to line
	int m = 0;	int n = 0;
	
	int superposition_check = 0;
	Point3D p_pop;
	//std:cout << "last slice size: " << v_next_slice.size() << endl;
	std::cout << "input boundary border size [xmin, xmax, ymin, ymax]: ";
	vector<cad3d::Point3D>  bound_box, v_test;
	int i_temp = 0;
	i_temp = find_point3d_min(v_last_slice, 1);
	std::cout << v_last_slice[i_temp].x << ", ";
	bound_box.push_back(v_last_slice[i_temp]);//xmin
	i_temp = find_point3d_max(v_last_slice, 1);
	std::cout << v_last_slice[i_temp].x << ", ";
	bound_box.push_back(v_last_slice[i_temp]);//xmax
	i_temp = find_point3d_min(v_last_slice, 2);
	std::cout << v_last_slice[i_temp].y << ", ";
	bound_box.push_back(v_last_slice[i_temp]);//ymin
	i_temp = find_point3d_max(v_last_slice, 2);
	std::cout << v_last_slice[i_temp].y << endl;
	bound_box.push_back(v_last_slice[i_temp]);//ymax

	for (int i = 0; i < v_next_slice.size(); i++) {
		
		if ((v_next_slice[i].x >= bound_box[0].x) & (v_next_slice[i].x <= bound_box[1].x)) {
			if ((v_next_slice[i].y >= bound_box[2].y) & (v_next_slice[i].y <= bound_box[3].y)) {
				v_out_slice.push_back(v_next_slice[i]);
			}
		}
	}
	std::cout << "super position test (in & out): " << v_next_slice.size() << ", " << v_out_slice.size() << endl;
	

}
int cad3d::test_if_in_shadow_v1(Point3D v_in, vector<Point3D> v_ref, float tol)
{
	int shadow_stat = 0;
	float test_dist = 0;
	for (int i = 0; i < v_ref.size(); i++) {
		test_dist = dist201(v_in, v_ref[i]);
		if (test_dist < tol) {
			shadow_stat = 1;
			break;
		}
	}
	return shadow_stat;
}
int cad3d::test_if_in_shadow_v2(Point3D v_in, vector<Point3D> v_ref, float tol)
{
	int shadow_stat = 0;
	float test_dist = 0;
	for (int i = 0; i < v_ref.size(); i++) {
		test_dist = dist201(v_in, v_ref[i]);
		if (test_dist < tol) {
			shadow_stat = 1;
			break;
		}
	}
	return shadow_stat;
}

void cad3d::view_scale_for_3d(vector<Point3D> &v_all, vector<Point3D> &v_in, vector<Point3D> &v_out, int &x_min,int &x_max, int &y_min, int &y_max ) {
	// for window size and all in block
	float x_range, y_range, z_range, z_min, z_max;
	
	int j = find_point3d_min(v_all, 1);
	x_min = v_all[j].x - 10;
	std::cout <<x_min << ", ";
	j = find_point3d_max(v_all, 1);
	x_max = v_all[j].x + 10;
	std::cout << x_max << ", ";
	x_range = x_max - x_min;

	j = find_point3d_min(v_all, 2);
	y_min = v_all[j].y - 10;
	std::cout << y_min << ", ";
	j = find_point3d_max(v_all, 2);
	y_max = v_all[j].y + 10;
	std::cout << y_max << ", ";
	y_range = y_max - y_min;

	j = find_point3d_min(v_all, 3);
	z_min = v_all[j].z - 10;
	std::cout << z_min << ", ";
	j = find_point3d_max(v_all, 3);
	z_max = v_all[j].z + 10;
	std::cout << z_max << ", ";
	z_range = z_max - z_min;
	//v_out.resize(size(v_in));
	for (int i = 0; i < v_in.size(); i++) {
		v_out.push_back(v_in[i]);
		v_out[size(v_out) - 1].x = v_in[i].x;
		v_out[size(v_out) - 1].y = v_in[i].y + y_range*(v_in[i].x - x_min) / x_range;
		v_out[size(v_out) - 1].z = v_in[i].z + z_range * (v_in[i].x - x_min) / x_range;
		
		/*std::cout << i <<": "<< v_out[i].x << ", ";
		std::cout << v_out[i].y<< ", " << v_out[i].z << endl;*/
		//some how left and right flipped, why? need to correct
	}
	x_min = y_min;
	y_min = z_min;
	x_max = y_min + 2 * y_range;
	y_max = z_min + 2 * z_range;
	/*plot_win3_min_x = target_plot3[j].y - 10;
	std::cout << plot_win3_min_x << ", ";
	j = blockplot.find_point3d_max(target_plot3, 3);
	plot_win3_max_y = target_plot2[j].z + 10;
	std::cout << plot_win3_max_y << ", ";
	j = blockplot.find_point3d_min(target_plot3, 3);
	plot_win3_min_y = target_plot2[j].z - 10;
	std::cout << plot_win3_min_y << endl;*/
}
void cad3d::view_scale_for_3d_v2(vector<Point3D> &v_all, int &x_min, int &x_max, int &y_min, int &y_max) {
	// for window size and coordinate lines at x,y, z min
	float x_range, y_range, z_range, z_min, z_max;

	int j = find_point3d_min(v_all, 1);
	x_min = v_all[j].x - 0.0;
	std::cout << x_min << ", ";
	j = find_point3d_max(v_all, 1);
	x_max = v_all[j].x + 0.0;
	std::cout << x_max << ", ";
	x_range = x_max - x_min;

	j = find_point3d_min(v_all, 2);
	y_min = v_all[j].y - 0.0;
	std::cout << y_min << ", ";
	j = find_point3d_max(v_all, 2);
	y_max = v_all[j].y + 0.0;
	std::cout << y_max << ", ";
	y_range = y_max - y_min;

	j = find_point3d_min(v_all, 3);
	z_min = v_all[j].z - 0.0;
	std::cout << z_min << ", ";
	j = find_point3d_max(v_all, 3);
	z_max = v_all[j].z + 0.0;
	std::cout << z_max << ", ";
	z_range = z_max - z_min;
	
	x_min = y_min;
	y_min = z_min;
	x_max = y_min + 2 * y_range;
	y_max = z_min + 2 * z_range;
	
}
void cad3d::view_scale_for_3d_v3(vector<Point3D> &v_all, vector<Point3D> &v_in, vector<Point3D> &v_out) {
	// for transforming lines and points
	float x_min, x_max, y_min, y_max, x_range, y_range, z_range, z_min, z_max;

	int j = find_point3d_min(v_all, 1);
	x_min = v_all[j].x - 0.0;
	//std::cout << x_min << ", ";
	j = find_point3d_max(v_all, 1);
	x_max = v_all[j].x + 0.0;
	//std::cout << x_max << ", ";
	x_range = x_max - x_min;

	j = find_point3d_min(v_all, 2);
	y_min = v_all[j].y - 0.0;
	//std::cout << y_min << ", ";
	j = find_point3d_max(v_all, 2);
	y_max = v_all[j].y + 0.0;
	//std::cout << y_max << ", ";
	y_range = y_max - y_min;

	j = find_point3d_min(v_all, 3);
	z_min = v_all[j].z - 0.0;
	//std::cout << z_min << ", ";
	j = find_point3d_max(v_all, 3);
	z_max = v_all[j].z + 0.0;
	//std::cout << z_max << ", ";
	z_range = z_max - z_min;
	//v_out.resize(size(v_in));
	for (int i = 0; i < v_in.size(); i++) {
		v_out.push_back(v_in[i]);
		v_out[size(v_out) - 1].x = v_in[i].x;
		v_out[size(v_out) - 1].y = v_in[i].y + y_range * (v_in[i].x - x_min) / x_range;
		v_out[size(v_out) - 1].z = v_in[i].z + z_range * (v_in[i].x - x_min) / x_range;

		//some how left and right flipped, why? need to correct
	}
	

}
void cad3d::calculate_3d_coordinate(vector<Point3D> &v_all, vector<Point3D> &show_coord_3d) {
	//output sequence is centre, xlinemax, ylinemax, zlinemax
	cad3d::Point3D ptemp;
	float x_range, y_range, z_range, z_min, z_max, y_min, y_max, x_min, x_max;

	int j = find_point3d_min(v_all, 1); x_min = v_all[j].x;
	j = find_point3d_max(v_all, 1); x_max = v_all[j].x;
	x_range = x_max - x_min;

	j = find_point3d_min(v_all, 2); y_min = v_all[j].y;
	j = find_point3d_max(v_all, 2); y_max = v_all[j].y ;
	y_range = y_max - y_min;

	j = find_point3d_min(v_all, 3); z_min = v_all[j].z;
	j = find_point3d_max(v_all, 3); z_max = v_all[j].z;
	z_range = z_max - z_min;
	
	//centre
	ptemp.x = x_min+5.0; ptemp.y = y_min +5.0; ptemp.z = z_min +5.0;
	show_coord_3d.push_back(ptemp);
	//x coordinate
	ptemp.x = x_min +0.0 +(x_range/8); ptemp.y = y_min + 5.0; ptemp.z = z_min + 5.0;
	show_coord_3d.push_back(ptemp);
	//y coordinate
	ptemp.x = x_min + 5.0; ptemp.y = y_min -0.0 + (y_range / 8); ptemp.z = z_min +5.0;
	show_coord_3d.push_back(ptemp);
	//z coordinate
	ptemp.x = x_min +5.0; ptemp.y = y_min +5.0 ; ptemp.z = z_min -0.0 + (z_range / 4);
	show_coord_3d.push_back(ptemp);/**/
}
void cad3d::line_view_update_3d() {
}


void cad3d::single_slot_seperate_edges_v101(stl_data &block_a, vector<Point3D> &open_slot_pts, vector<Point3D> &v1, vector<Point3D> &v2, vector<Point3D> &v_int_l1, vector<Point3D> &v_int_l2, vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4) {
	//Cartesian Processing
	vector<int> vi_temp551, vi_temp552;
	vector<cad3d::Point3D>  slot_e01, slot_e02, slot_e03, slot_e04;
	//vector<cad3d::Point3D>  slot_e1, slot_e2, slot_e3, slot_e4;
	vector<cad3d::Point3D>  slot_ip1, slot_ip2;
	
	int i1, i2;
	std::cout << "corner point index 1:";
	std::cin >> i1;
	v1.push_back(open_slot_pts[i1]);//i_st1
	slot_e01.push_back(open_slot_pts[i1]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, open_slot_pts[i1], v_int_l1, v_int_l2, v1);
	slot_e02.push_back(v1[v1.size() - 1]);
	get_semi_open_slot_v102(block_a.vertices, slot_e01, slot_e02, v_int_l1, v_int_l2, v1);
	crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	//newblock.get_internal_slot_v107_phase6(block_a.vertices, slot_e01, slot_e02, v_int_l1, v_int_l2, v1);
	//newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	std::cout << "slot e1 size: " << slot_e01.size() << "~ " << slot_e1.size() << endl;

	std::cout << "corner point index 2:";
	std::cin >> i2;
	v2.push_back(open_slot_pts[i2]);//i_st1
	slot_e03.push_back(open_slot_pts[i2]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, open_slot_pts[i2], v_int_l1, v_int_l2, v2);
	slot_e04.push_back(v2[v2.size() - 1]);
	get_semi_open_slot_v102(block_a.vertices, slot_e03, slot_e04, v_int_l1, v_int_l2, v2);
	crop_slot_single_edge(slot_e3, slot_e4, slot_e03, slot_e04);
	std::cout << "slot e1 size: " << slot_e03.size() << "~ " << slot_e3.size() << endl;

	//for path symmetric structures
	/*for (int i = 2; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}*/
}
void cad3d::single_slot_seperate_edges_v102(stl_data &block_a, vector<Point3D> &open_slot_pts, vector<Point3D> &v1, vector<Point3D> &v2, vector<Point3D> &v_int_l1, vector<Point3D> &v_int_l2, vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_e3, vector<Point3D> &slot_e4) {
	//Cartesian Processing
	vector<int> vi_temp551, vi_temp552;
	vector<cad3d::Point3D>  slot_e01, slot_e02, slot_e03, slot_e04;
	//vector<cad3d::Point3D>  slot_e1, slot_e2, slot_e3, slot_e4;
	vector<cad3d::Point3D>  slot_ip1, slot_ip2;

	int i1, i2;
	std::cout << "corner point index 1:";
	std::cin >> i1;
	v1.push_back(open_slot_pts[i1]);//i_st1
	slot_e01.push_back(open_slot_pts[i1]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, open_slot_pts[i1], v_int_l1, v_int_l2, v1);
	slot_e02.push_back(v1[v1.size() - 1]);
	get_semi_open_slot_v102(block_a.vertices, slot_e01, slot_e02, v_int_l1, v_int_l2, v1);
	crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	//newblock.get_internal_slot_v107_phase6(block_a.vertices, slot_e01, slot_e02, v_int_l1, v_int_l2, v1);
	//newblock.crop_slot_single_edge(slot_e1, slot_e2, slot_e01, slot_e02);
	std::cout << "slot e1 size: " << slot_e01.size() << "~ " << slot_e1.size() << endl;

	std::cout << "corner point index 2:";
	std::cin >> i2;
	v2.push_back(open_slot_pts[i2]);//i_st1
	slot_e03.push_back(open_slot_pts[i2]);
	vi_temp551.resize(0);
	//vi_temp552.resize(0);
	one_start_point_polar_block_v1_phase4(block_a.vertices, open_slot_pts[i2], v_int_l1, v_int_l2, v2);
	slot_e04.push_back(v2[v2.size() - 1]);
	get_semi_open_slot_v102(block_a.vertices, slot_e03, slot_e04, v_int_l1, v_int_l2, v2);
	crop_slot_single_edge(slot_e3, slot_e4, slot_e03, slot_e04);
	std::cout << "slot e1 size: " << slot_e03.size() << "~ " << slot_e3.size() << endl;

	//for path symmetric structures
	/*for (int i = 2; i < slot_e1.size(); i++) {
	slot_ip1.push_back(slot_e1[i]);
	//i = i + 2;
	}*/
}
void cad3d::single_slot_seperate_edges_v103(stl_data &block_a, vector<Point3D> &v_in, vector<Point3D> &v1, vector<Point3D> &v2) {
	//used for split edges at sub-z level of sub-volume
	// get all the slice points
	cad3d::circle2d c_temp;
	vector<Point3D> zslice_points;
	float sub_z;
	std::cout << "sub z slice level:" ;
	std::cin >> sub_z;
	find_point_3d_1feature(zslice_points, v_in, sub_z, 0.1, 3);
	std::cout << "points in the slice plane: " << zslice_points.size() << endl;
	
	v1.push_back(zslice_points[0]);
	v1.push_back(zslice_points[1]);
	for (int i = 2; i < zslice_points.size(); i++) {
		c_temp = compute_circle_4_point3d(zslice_points[i], zslice_points[i - 1], zslice_points[i - 2]);
		std::cout << i << ": " << c_temp.r << endl;
		v1.push_back(zslice_points[i]);
		if (abs(c_temp.r) > 50) {
			//v1.push_back(v1[i]);
			v1.push_back(zslice_points[i]);
		}
		else if (abs(c_temp.r) <= 100) {
			break;
		}
	}
	v2.push_back(zslice_points[zslice_points.size()-1]);
	v2.push_back(zslice_points[zslice_points.size() - 1]);
	for (int i = zslice_points.size() - 3; i >-1; i--) {
		c_temp = compute_circle_4_point3d(zslice_points[i], zslice_points[i + 1], zslice_points[i + 2]);
		//std::cout << i << ": " << c_temp.r << endl;
		if (abs(c_temp.r) > 50) {
			//v1.push_back(v1[i]);
			v2.push_back(zslice_points[i]);
		}
		else if (abs(c_temp.r) <= 100) {
			break;
		}
	}
	
	
}

void cad3d::interpolate_edge_if_necessary(vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1) {
	//interpolate edges if necessary
	float edge_point_dist;
	Point3D intp_pt;
	for (int i = 0; i < slot_ip1.size() - 4; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		std::cout << i << ": " << edge_point_dist << endl;
		if (abs(edge_point_dist) > 5) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.25), 3, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 3, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
		}
	}
}
void cad3d::interpolate_edge_if_necessary_v2(vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1) {
	//interpolate edges if necessary
	float edge_point_dist;
	Point3D intp_pt;
	//interpolated_edge_v1.push_back(slot_ip1[0]);
	for (int i = 0; i < slot_ip1.size() - 5; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		//std::cout << i << ": " << edge_point_dist <<"~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);
			
		}*/
		if ((abs(edge_point_dist) > 5) & (abs(edge_point_dist) < 10)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i+1], intp_pt) << ", ";
			/*interpolate_lagrange_v6(slot_ip1, (i + 0.66), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);*/
		}
		else if ((abs(edge_point_dist) > 10) & (abs(edge_point_dist) < 15)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.33), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.66), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			/*interpolate_lagrange_v6(slot_ip1, (i + 0.75), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);*/
		}
		else if (abs(edge_point_dist) > 15) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.25), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 4, intp_pt,i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
		}
		//std::cout << endl;
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size()-5]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v3(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1) {
	//interpolate edges if necessary incartesian mode
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}
	slot_ip1.push_back(slot_e2[slot_e2.size() - 1]);
	std::cout << "Input edge size: " << slot_ip1.size() << endl;
	
	float edge_point_dist;
	Point3D intp_pt;
	//interpolated_edge_v1.push_back(slot_ip1[0]);
	for (int i = 0; i < slot_ip1.size() - 4; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		std::cout << i << ": " << edge_point_dist <<"~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);

		}*/
		if ((abs(edge_point_dist) > 4) & (abs(edge_point_dist) < 8)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i+1], intp_pt) << ", ";
			/*interpolate_lagrange_v6(slot_ip1, (i + 0.66), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);*/
		}
		else if ((abs(edge_point_dist) > 8) & (abs(edge_point_dist) < 12)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.33), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.66), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			/*interpolate_lagrange_v6(slot_ip1, (i + 0.75), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);*/
		}
		else if (abs(edge_point_dist) > 12) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.25), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 3, intp_pt, i);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
		}
		std::cout << endl;
	}
	//interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 5]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v4(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1) {
	//interpolate edges if necessary incartesian mode
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}
	slot_ip1.push_back(slot_e2[slot_e2.size() - 1]);
	std::cout << "Input edge size: " << slot_ip1.size() << endl;

	float edge_point_dist;
	Point3D intp_pt;
	interpolated_edge_v1.push_back(slot_ip1[0]);
	edge_point_dist = dist301(slot_ip1[0], slot_ip1[1]);
	
	if (abs(edge_point_dist) > 5) {
		interpolate_lagrange_v6(slot_ip1, (0 + 0.33), 3, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
		//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
		interpolate_lagrange_v6(slot_ip1, (0 + 0.66), 3, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
	}
	// test with interpolate all points
	for (int i = 1; i < slot_ip1.size()-2; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		std::cout << i << ": " << edge_point_dist << "~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);

		}*/
		if ((abs(edge_point_dist) > 4) & (abs(edge_point_dist) < 8)) {
			//interpolate_lagrange_v6(slot_ip1, (i + 0.5), slot_ip1.size(), intp_pt, 0);
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 4, intp_pt, i-1);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i+1], intp_pt) << ", ";
			/*interpolate_lagrange_v6(slot_ip1, (i + 0.66), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), 4, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);*/
		}
		//else if ((abs(edge_point_dist) > 8) & (abs(edge_point_dist) < 12)) {
		else if (abs(edge_point_dist) > 8) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.33), 4, intp_pt, i-1);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.66), 4, intp_pt, i-1);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			
		}
		/*else if (abs(edge_point_dist) > 12) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.25), slot_ip1.size(), intp_pt, 0);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), slot_ip1.size(), intp_pt, 0);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
			interpolate_lagrange_v6(slot_ip1, (i + 0.75), slot_ip1.size(), intp_pt, 0);
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i + 1], intp_pt) << ", ";
		}*/
		std::cout << endl;
	}
	//interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 5]);
	/*interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);*/
	
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 3], slot_ip1[slot_ip1.size() - 2]);
	if (abs(edge_point_dist) > 5) {
		interpolate_lagrange_v6(slot_ip1, (slot_ip1.size() - 2.5), 5, intp_pt, (slot_ip1.size() - 6));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v5(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1) {
	//interpolate edges if necessary with spline
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}
	slot_ip1.push_back(slot_e2[slot_e2.size() - 1]);
	//std::cout << "Input edge size: " << slot_ip1.size() << endl;

	float edge_point_dist;
	Point3D intp_pt;
	interpolated_edge_v1.push_back(slot_ip1[0]);
	edge_point_dist = dist301(slot_ip1[0], slot_ip1[1]);

	if (abs(edge_point_dist) > 4) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (0.5), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
		
	}
	// test with interpolate all points
	for (int i = 1; i < slot_ip1.size() - 5; i++) {
	//for (int i = 1; i < 5; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		//std::cout << i << ": " << edge_point_dist << "~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
			interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
			interpolated_edge_v1.push_back(intp_pt);

		}*/
		if (abs(edge_point_dist) > 4)  {
			//interpolate_lagrange_v6(slot_ip1, (i + 0.5), slot_ip1.size(), intp_pt, 0);
			//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
			interpolate_spline_v1(slot_ip1, (i + 0.5), 4, intp_pt, (i-1));
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i+1], intp_pt) << ", ";
			
		}
		
		//std::cout << endl;
	}

	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 5]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 5], slot_ip1[slot_ip1.size() - 4]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 4.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 4], slot_ip1[slot_ip1.size() - 3]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 3.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 3], slot_ip1[slot_ip1.size() - 2]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 2], slot_ip1[slot_ip1.size() - 1]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v6(vector<Point3D> &slot_e1, vector<Point3D> &slot_e2, vector<Point3D> &slot_ip1, vector<Point3D> &interpolated_edge_v1,float tol) {
	//interpolate edges if necessary with spline and with tol distance value
	// slot e1 and e2 are main input
	//slot ip1 is temporary point array
	//interpolated edge v1 is the output
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}
	slot_ip1.push_back(slot_e2[slot_e2.size() - 1]);
	//std::cout << "Input edge size: " << slot_ip1.size() << endl;

	float edge_point_dist;
	Point3D intp_pt;
	interpolated_edge_v1.push_back(slot_ip1[0]);
	edge_point_dist = dist301(slot_ip1[0], slot_ip1[1]);

	if (abs(edge_point_dist) > 4) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (0.5), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);

	}
	// test with interpolate all points
	for (int i = 1; i < slot_ip1.size() - 5; i++) {
		//for (int i = 1; i < 5; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		//std::cout << i << ": " << edge_point_dist << "~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
		interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
		interpolated_edge_v1.push_back(intp_pt);

		}*/
		if (abs(edge_point_dist) > 4) {
			//interpolate_lagrange_v6(slot_ip1, (i + 0.5), slot_ip1.size(), intp_pt, 0);
			//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
			interpolate_spline_v1(slot_ip1, (i + 0.5), 4, intp_pt, (i - 1));
			interpolated_edge_v1.push_back(intp_pt);
			//std::cout << dist301(slot_ip1[i+1], intp_pt) << ", ";

		}

		//std::cout << endl;
	}

	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 5]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 5], slot_ip1[slot_ip1.size() - 4]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 4.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 4], slot_ip1[slot_ip1.size() - 3]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 3.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 3], slot_ip1[slot_ip1.size() - 2]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 2], slot_ip1[slot_ip1.size() - 1]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v7(vector<Point3D> &slot_e1, vector<Point3D> &interpolated_edge_v1, float tol) {
	vector<Point3D> slot_ip1;
	//interpolate edges if necessary with spline
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}
	

	float edge_point_dist;
	Point3D intp_pt;
	interpolated_edge_v1.push_back(slot_ip1[0]);
	edge_point_dist = dist301(slot_ip1[0], slot_ip1[1]);

	if (abs(edge_point_dist) > 4) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (0.33), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (0.66), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
	}
	// test with interpolate all points
	for (int i = 1; i < slot_ip1.size() - 3; i++) {
		//for (int i = 1; i < 5; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		//std::cout << i << ": " << edge_point_dist << "~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
		interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
		interpolated_edge_v1.push_back(intp_pt);

		}*/
		if (abs(edge_point_dist) > 4) {
			
			//interpolate_spline_v1(slot_ip1, (i + 0.5), 4, intp_pt, (i - 1));
			//interpolated_edge_v1.push_back(intp_pt);
			interpolate_spline_v1(slot_ip1, (i + 0.33), 2, intp_pt, (i - 1));
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_spline_v1(slot_ip1, (i + 0.66), 2, intp_pt, (i - 1));
			interpolated_edge_v1.push_back(intp_pt);
			

		}

		//std::cout << endl;
	}

	/*interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 5]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 5], slot_ip1[slot_ip1.size() - 4]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 4.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 4]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 4], slot_ip1[slot_ip1.size() - 3]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 3.5), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}*/
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 3], slot_ip1[slot_ip1.size() - 2]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.66), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.33), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 2], slot_ip1[slot_ip1.size() - 1]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.66), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.33), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}
void cad3d::interpolate_edge_if_necessary_v8(vector<Point3D> &slot_e1, vector<Point3D> &interpolated_edge_v1, float tol) {
	vector<Point3D> slot_ip1;
	//interpolate edges if necessary with spline
	for (int i = 0; i < slot_e1.size(); i++) {
		slot_ip1.push_back(slot_e1[i]);
		//i = i + 2;
	}


	float edge_point_dist;
	Point3D intp_pt;
	interpolated_edge_v1.push_back(slot_ip1[0]);
	edge_point_dist = dist301(slot_ip1[0], slot_ip1[1]);

	if (abs(edge_point_dist) > 4) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (0.33), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (0.66), 4, intp_pt, 0);
		interpolated_edge_v1.push_back(intp_pt);
	}
	// test with interpolate all points
	for (int i = 1; i < slot_ip1.size() - 3; i++) {
		//for (int i = 1; i < 5; i++) {
		edge_point_dist = dist301(slot_ip1[i], slot_ip1[i + 1]);
		interpolated_edge_v1.push_back(slot_ip1[i]);
		//std::cout << i << ": " << edge_point_dist << "~ ";
		/*if ((abs(edge_point_dist) > 2) & (abs(edge_point_dist) < 5)) {
		interpolate_lagrange_v6(slot_ip1, (i + 0.5), 3, intp_pt);
		interpolated_edge_v1.push_back(intp_pt);

		}*/
		if (abs(edge_point_dist) > tol) {

			//interpolate_spline_v1(slot_ip1, (i + 0.5), 4, intp_pt, (i - 1));
			//interpolated_edge_v1.push_back(intp_pt);
			interpolate_spline_v1(slot_ip1, (i + 0.33), 2, intp_pt, (i - 1));
			interpolated_edge_v1.push_back(intp_pt);
			interpolate_spline_v1(slot_ip1, (i + 0.66), 2, intp_pt, (i - 1));
			interpolated_edge_v1.push_back(intp_pt);


		}

		//std::cout << endl;
	}

	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 3]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 3], slot_ip1[slot_ip1.size() - 2]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.66), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 2.33), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 2]);
	edge_point_dist = dist301(slot_ip1[slot_ip1.size() - 2], slot_ip1[slot_ip1.size() - 1]);
	if (abs(edge_point_dist) > 5) {
		//interpolate_spline_v1(std::vector<Point3D> &v1, float ii, int n1, Point3D &v2, int i_present)
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.66), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
		interpolate_spline_v1(slot_ip1, (slot_ip1.size() - 1.33), 4, intp_pt, (slot_ip1.size() - 5));
		interpolated_edge_v1.push_back(intp_pt);
	}
	interpolated_edge_v1.push_back(slot_ip1[slot_ip1.size() - 1]);
}

float cad3d::find_rotation() {
	float theta=0.0;
	return theta;
}
void cad3d::rotate_translate_point3d(vector<cad3d::Point3D> &vi, vector<cad3d::Point3D> &vo, float theta) {
	vo.resize(vi.size());

	MatrixXd rot(2, 2);
	rot(0, 0) = cos(theta);
	rot(0, 1) = -sin(theta);
	rot(1, 0) = sin(theta);
	rot(1, 1) = cos(theta);
	//https://en.wikipedia.org/wiki/Rotation_matrix
	//https://stackoverflow.com/questions/37443689/eigen-rotation2d-applied-to-a-2-vector

	MatrixXd v_temp(2, 1);
	MatrixXd v_ot(2, 1);
	float delx, dely;
	//measure the initial displacement
	v_temp(0, 0) = vi[0].x;
	v_temp(1, 0) = vi[0].y;
	v_ot = rot * v_temp;
	delx = vi[0].x - v_ot(0, 0);
	dely = vi[0].y - v_ot(1, 0);
	vo[0].x = v_ot(0, 0) + delx;
	vo[0].y = v_ot(1, 0) + dely;
	vo[0].z = vi[0].z;

	float dxr, dyr, dzr;

	for (int i = 1; i < vi.size(); i++) {
		v_temp(0, 0) = vi[i].x;
		v_temp(1, 0) = vi[i].y;
		v_ot = rot * v_temp;
		vo[i].x = v_ot(0, 0) + delx;
		vo[i].y = v_ot(1, 0) + dely;
		vo[i].z = vi[i].z;
	}
}
void cad3d::align_with_robot_CAD_local_coords(vector<Point3D> &v_in, vector<Point3D> &v_out) {
	Point3D cad_local_top_corner1, cad_local_top_corner2, cad_local_top_corner3, cad_local_top_corner4;
	vector<cad3d::Point3D> cl_rot;
	/*float theta = 0.0;
	float theta_r = theta * 3.14159 / 180.0;
	rotate_translate_point3d(cl4, cl_rot, theta_r);

	string robot_f = "F:/rt_works/Impeller_int/RC1/Program/centerline.txt";
	ofstream vout(robot_f);
	vout << cl_rot.size() + 1 << endl;
	//needs use of rotation matrix, but check for now
	for (int i = 0; i < cl_rot.size(); i++) {
		vout << "( " << robot_tx[0] - cl_rot[i].y << ", ";
		vout << robot_tx[1] - cl_rot[i].x << ", ";
		//vout << robot_tx[2] + cl3[i].z << ", ";
		vout << robot_tx[2] << ", ";
		vout << robot_tx[3] << ", ";
		vout << robot_tx[4] << ", ";
		vout << robot_tx[5] << ") " << endl;
	}
	vout.close();
	string c3_file_name = "F:/matlab_works/cad_bologna/cl3_54_cad15_v10.txt";
	newblock.point3d_write(c3_file_name, cl3);*/
}

void cad3d::generate_border_lines(stl_data &block_a, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, float z_level) {
	vector<cad3d::triangle_index> tsloti_feature3, tslot_i31, tslot_i32, ti_all3;
	vector<cad3d::Point3D>  v_temp, feat01_e1, feat01_e2,slot_plane1, feat_e1_border, feat_e2_border;
	vector<int> i_line;

	cross_model_generate(block_a, slot_plane1, ti_all3, tslot_i31, tslot_i32, tsloti_feature3);
	std::cout << "possible feature triangle no: " << tsloti_feature3.size() << endl;
	triangle2line(block_a.vertices, tsloti_feature3, feat01_e1, feat01_e2, i_line, v_temp);
	order_cartesian_points_v1(feat01_e1, feat01_e2, feat_e1_border, feat_e2_border);
	//std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
	int ixmax, iymax, ixmin, iymin;
	cartesian_separate_external_regions_v1(feat_e1_border, feat_e2_border, v_ext_l1, v_ext_l2, 0, slot_plane1);
}
void cad3d::cartesian_separate_external_regions_v1(vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, int cp_mode, vector<Point3D> &slot_plane1) {
	//cp mode is Cartesian 0 or polar 1 mode
	int int_boundary = 0;
	std::cout << "internal region number: ";
	std::cin >> int_boundary;

	if (cp_mode == 0) {
		bool a1, b1, c1, d1, e, f;
		bool a2, b2, c2, d2;
		bool m, n;
		//corner point estmation 
		int ixmax, iymax, ixmin, iymin;

		ixmax = find_point3d_max(slot_plane1, 1);
		ixmin = find_point3d_min(slot_plane1, 1);
		iymax = find_point3d_max(slot_plane1, 2);
		iymin = find_point3d_min(slot_plane1, 2);
		//std::cout << slot_plane1[ixmax].x << ", " << slot_plane1[ixmin].x << ", " << slot_plane1[iymax].y << ", " << slot_plane1[iymin].y << endl;
		for (int i = 0; i < feat1_e1.size(); i++) {
			a1 = abs(feat1_e1[i].x - slot_plane1[ixmax].x) < 0.1;
			b1 = abs(feat1_e1[i].x - slot_plane1[ixmin].x) < 0.1;
			c1 = abs(feat1_e1[i].y - slot_plane1[iymax].y) < 0.1;
			d1 = abs(feat1_e1[i].y - slot_plane1[iymin].y) < 0.1;

			e = (a1 || b1) || (c1 || d1); // at any corner point and at same z plane e = 1

			a2 = abs(feat1_e2[i].x - slot_plane1[ixmax].x) < 0.1;
			b2 = abs(feat1_e2[i].x - slot_plane1[ixmin].x) < 0.1;
			c2 = abs(feat1_e2[i].y - slot_plane1[iymax].y) < 0.1;
			d2 = abs(feat1_e2[i].y - slot_plane1[iymin].y) < 0.1;


			f = (a2 || b2) || (c2 || d2); // at any corner point f = 1
			//std::cout << i << ": " << e << ", " << f;
			//if any of the point is internal, include in internal feature edges
			
			m = abs(feat1_e1[i].z - feat1_e2[i].z) ;

			if ((e == 1) & (f == 1) & (m < 1)) {
				//std::cout << " internal " << feat1_e1[i].x << ", "<<feat1_e1[i].y;
				//std::cout << " ~ " << feat1_e1[i].z << ", " << feat1_e2[i].z << "=>"<<m << endl;
				v_ext_l1.push_back(feat1_e1[i]);
				v_ext_l2.push_back(feat1_e2[i]);
			}
			/*else if (!f ==1) {
				v_int_l1.push_back(feat1_e1[i]);
				v_int_l2.push_back(feat1_e2[i]);
			}*/
			//std::cout << endl;
			a1 = 0; b1 = 0; c1 = 0; d1 = 0;
			a2 = 0; b2 = 0; c2 = 0; d2 = 0;
			e = 0; f = 0;
			m = 0;
		}
	}
	//std::cout << "v_int_l1 size: " << v_ext_l1.size() << endl;
}
void cad3d::generate_border_line_polar(stl_data &block_a, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2) {
	int a = 1;
	//feature_plane_approximation(block_a, zslice_points, z_now);
	//polar_processing(block_a, slot_plane, tslot_i, v_ext_l1, v_ext_l2, feat_e1, feat_e2, z_now, centered_block);
	vector<cad3d::Point3D> slot_plane1, slot_plane2, v_used, st1, st2, v_temp;
	vector<int> i_line;
	vector<cad3d::triangle_index> tsloti_feature, tslot_i1, tslot_i2, ti_all2;
	vector<cad3d::Point3D>  centered_block,feat01_e1, feat01_e2, feat1_e1, feat1_e2;
	centre_at_x0y0(block_a.vertices, centered_block);
	cross_model_generate_polar(block_a, centered_block, slot_plane1, ti_all2, tslot_i1, tslot_i2, tsloti_feature);
	std::cout << "cross model generated: " << slot_plane1.size() <<endl;
	triangle2line(block_a.vertices, tsloti_feature, feat01_e1, feat01_e2, i_line, v_temp);
	vector<cad3d::Point3D> feat01_e1_polar, feat01_e2_polar, slot_plane1_polar;
	std::cout << "triangle to line transferred: "  << feat01_e1.size() << endl;
	convert_cart2pol(feat01_e1, feat01_e1_polar);
	convert_cart2pol(feat01_e2, feat01_e2_polar);
	order_polar_points_v1(feat01_e1, feat01_e2, feat1_e1, feat1_e2, feat01_e1_polar, feat01_e2_polar);
	//std::cout << "feat1_e1 size: " << feat1_e1.size() << endl;
	polar_separate_external_regions_v1(feat1_e1, feat1_e2, v_ext_l1, v_ext_l2, 0, slot_plane1);
}
void cad3d::polar_separate_external_regions_v1(vector<Point3D> &feat1_e1, vector<Point3D> &feat1_e2, vector<Point3D> &v_ext_l1, vector<Point3D> &v_ext_l2, int cp_mode, vector<Point3D> &slot_plane1) {
	vector<cad3d::Point3D> feat1_e1_polar, feat1_e2_polar,slot_plane1_polar;
	
	convert_cart2pol(feat1_e1, feat1_e1_polar);
	convert_cart2pol(feat1_e2, feat1_e2_polar);
	convert_cart2pol(slot_plane1, slot_plane1_polar);

	bool a1, b1, c1, d1, e, f;
	bool a2, b2, c2, d2;
	//corner point estmation 
	int ixmax,ixmax1, ixmax2;
	
	ixmax = find_point3d_max(feat1_e1_polar, 1);
	ixmax1 = find_point3d_max(slot_plane1_polar, 1);
	ixmax2 = find_point3d_max(feat1_e2_polar, 1);
	
	for (int i = 0; i < feat1_e1.size(); i++) {
		a1 = abs(feat1_e1_polar[i].x - slot_plane1_polar[ixmax1].x) < 0.1;
		// if e1 lies on extrnal point
		
		a2 = abs(feat1_e2_polar[i].x - slot_plane1_polar[ixmax1].x) < 0.1;
		//if e2 lies on external point
		
		//if aboth the points are external, include in external feature edges
		if ((a1 == 1) & (a2 == 1)) {
			//std::cout << " internal " << feat1_e1[i].x << ", "<<feat1_e1[i].y;
			//std::cout << " ~ " << feat1_e2[i].x << ", " << feat1_e2[i].y;
			v_ext_l1.push_back(feat1_e1[i]);
			v_ext_l2.push_back(feat1_e2[i]);
		}
		
		//std::cout << endl;
		a1 = 0; //b1 = 0; c1 = 0; d1 = 0;
		a2 = 0; //b2 = 0; c2 = 0; d2 = 0;
		//e = 0; f = 0;
	}
}


void cad3d::find_close_points_from_edge(vector<Point3D> &v_in, vector<Point3D> &v_out,Point3D p_base,float d_tol) {
	//find distance from edges close to 70% of width
	float dtemp;
	for (int i = 0; i < v_in.size(); i++) {
		dtemp = dist301(p_base, v_in[i]);
		if (dtemp < 0.7 * d_tol) {
			v_out.push_back(v_in[i]);
			//normal_intersect.push_back(interpolated_edge_v1[i]);
		}
	}

}
void cad3d::get_straight_line_2d(float &m, float &c, Point3D tan1, Point3D tan2, Point3D norm0) {
	// y = mx + c 
	//returns m and c
	//https://www.ugrad.math.ubc.ca/coursedoc/math100/notes/zoo/eqline.html
	vector<cad3d::Point3D> normal_intersect, near_edge_points;
	//cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
	float dy, dx, m_n, c_n, slot_width, dtemp;

	/*tan1 = copy_point3d(c_line[i_cent - 1]);
	tan2 = newblock.copy_point3d(c_line[i_cent + 1]);
	norm0 = newblock.copy_point3d(c_line[i_cent]);*/

	dy = tan2.y - tan1.y;
	dx = tan2.x - tan1.x;
	m = dy/dx;
	c = -m*norm0.x + norm0.y;
	//equation of normal line
}
cad3d::Point3D cad3d::get_intersect_point_phase3_v1(Point3D norm0, float m1, float c1, float m2, float c2) {
	Point3D intersect_pt; 
	//float  m1, c1, m2, c2;
	intersect_pt.z = norm0.z;
	if ((m1-m2)==0) {
		intersect_pt.x = 0.0;
		intersect_pt.y = 0.0;
		return intersect_pt;
	}
	intersect_pt.x = (c2-c1) / (m1-m2);
	intersect_pt.y = m1*intersect_pt.x + c1;
	return intersect_pt;
}
void cad3d::generate_single_cross_along_centreline_v2(vector<Point3D> &edge1, vector<Point3D> &edge2, vector<Point3D> &c_line, vector<cad3d::Point3D> &normal_intersect) {
	//with normal_intersect, it should send only two pi=points at every edge 
	// find normal to centerline at a given point and intersection, start
	int i_cent = 0;
	std::cout << "centre line index: ";
	std::cin >> i_cent;
	vector<cad3d::Point3D> near_edge_points;
	cad3d::Point3D tan1, tan2, norm1, norm2, norm0, intersect_p1, intersect_p2;
	float dy, dx, m_n, c_n, slot_width, dtemp;
	float m_t, c_t, mt1, ct1, mt2, ct2;
	std::cout << "slot width: ";
	std::cin >> slot_width;

	tan1	= copy_point3d(c_line[i_cent - 1]);
	tan2	= copy_point3d(c_line[i_cent + 1]);
	norm0	= copy_point3d(c_line[i_cent]);


	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
	m_n = -1 / m_t;
	c_n = -m_n*norm0.x + norm0.y;

	find_close_points_from_edge(edge1, near_edge_points, c_line[i_cent], slot_width);
	//from edge 1
	if (near_edge_points.size()>0) {
		get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p1 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt1, ct1);
		std::cout << near_edge_points.size() << ", " << intersect_p1.x << ", " << intersect_p1.y << endl;
		normal_intersect.push_back(intersect_p1);
	}
	near_edge_points.resize(0);

	//from edge 2
	find_close_points_from_edge(edge2, near_edge_points, c_line[i_cent], slot_width);
	if (near_edge_points.size() > 0) {
		get_straight_line_2d(mt2, ct2, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p2 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt2, ct2);
		normal_intersect.push_back(intersect_p2);
		std::cout << near_edge_points.size() << ", " << intersect_p2.x << ", " << intersect_p2.y << endl;
		near_edge_points.resize(0);
		normal_intersect.push_back(c_line[i_cent]);
		std::cout << normal_intersect.size() << endl;
	}
	

}
void cad3d::generate_single_cross_along_centreline_v3(vector<Point3D> &edge1, vector<Point3D> &edge2, vector<Point3D> &c_line, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2) {
	//with normal_intersect, it should send only two pi=points at every edge 
	// find normal to centerline at a given point and intersection, start
	int i_cent = 0;
	std::cout << "centre line index: ";
	std::cin >> i_cent;
	vector<cad3d::Point3D> near_edge_points;
	cad3d::Point3D tan1, tan2, norm1, norm2, norm0, intersect_p1, intersect_p2;
	float dy, dx, m_n, c_n, slot_width, dtemp;
	float m_t, c_t, mt1, ct1, mt2, ct2;
	std::cout << "slot width: ";
	std::cin >> slot_width;

	tan1 = copy_point3d(c_line[i_cent - 1]);
	tan2 = copy_point3d(c_line[i_cent + 1]);
	norm0 = copy_point3d(c_line[i_cent]);


	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
	m_n = -1 / m_t;
	c_n = -m_n*norm0.x + norm0.y;

	
	//from edge 1
	find_close_points_from_edge(edge1, near_edge_points, c_line[i_cent], slot_width);
	if (near_edge_points.size()>0) {
		intersect_l1.push_back(c_line[i_cent]);
		get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p1 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt1, ct1);
		std::cout << near_edge_points.size() << ", " << intersect_p1.x << ", " << intersect_p1.y << ", "<<intersect_p1.z << endl;
		//normal_intersect.push_back(intersect_p1);
		intersect_l2.push_back(intersect_p1);
	}
	near_edge_points.resize(0);

	
	//from edge 2
	find_close_points_from_edge(edge2, near_edge_points, c_line[i_cent], slot_width);
	if (near_edge_points.size()>0) {
		intersect_l1.push_back(c_line[i_cent]);
		get_straight_line_2d(mt2, ct2, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p2 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt2, ct2);
		//normal_intersect.push_back(intersect_p2);
		intersect_l2.push_back(intersect_p2);
		std::cout << near_edge_points.size() << ", " << intersect_p2.x << ", " << intersect_p2.y << endl;
	}
	near_edge_points.resize(0);
	//normal_intersect.push_back(c_line[i_cent]);
	//std::cout << normal_intersect.size() << endl;

}
void cad3d::generate_single_cross_along_centreline_v4(vector<Point3D> &edge1, vector<Point3D> &edge2, vector<Point3D> &c_line, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2) {
	//with normal_intersect, it should send only two pi=points at every edge 
	// find normal to centerline at a given point and intersection, start
	int i_cent = 0;
	std::cout << "centre line index: ";
	std::cin >> i_cent;
	vector<cad3d::Point3D> near_edge_points;
	cad3d::Point3D tan1, tan2, norm1, norm2, norm0, intersect_p1, intersect_p2;
	float dy, dx, m_n, c_n, slot_width, dtemp;
	float m_t, c_t, mt1, ct1, mt2, ct2;
	std::cout << "slot width: ";
	std::cin >> slot_width;

	tan1 = copy_point3d(c_line[i_cent - 1]);
	tan2 = copy_point3d(c_line[i_cent + 1]);
	norm0 = copy_point3d(c_line[i_cent]);


	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
	m_n = -1 / m_t;
	c_n = -m_n*norm0.x + norm0.y;


	//from edge 1
	find_close_points_from_edge(edge1, near_edge_points, c_line[i_cent], slot_width);
	if (near_edge_points.size()>0) {
		//intersect_l1.push_back(c_line[i_cent]);
		get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p1 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt1, ct1);
		std::cout << near_edge_points.size() << ", " << intersect_p1.x << ", " << intersect_p1.y << ", " << intersect_p1.z << endl;
		//normal_intersect.push_back(intersect_p1);
		//intersect_l2.push_back(intersect_p1);
		intersect_l1.push_back(intersect_p1);
	}
	near_edge_points.resize(0);


	//from edge 2
	find_close_points_from_edge(edge2, near_edge_points, c_line[i_cent], slot_width);
	if (near_edge_points.size()>0) {
		//intersect_l1.push_back(c_line[i_cent]);
		get_straight_line_2d(mt2, ct2, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		intersect_p2 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt2, ct2);
		//normal_intersect.push_back(intersect_p2);
		intersect_l2.push_back(intersect_p2);
		std::cout << near_edge_points.size() << ", " << intersect_p2.x << ", " << intersect_p2.y << endl;
	}
	near_edge_points.resize(0);
	//normal_intersect.push_back(c_line[i_cent]);
	//std::cout << normal_intersect.size() << endl;

}

void cad3d::generate_single_cross_along_centreline_v5(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line, vector<Point3D> &normal_intersect, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, float tol) {
	//created on 11.12.19

	std::cout << "3D slice processings " << endl;
	int i_cent = 0;
	std::cout << "centre line index: ";
	std::cin >> i_cent;
	//normal_intersect.push_back(c_line[i_cent - 1]);
	
	//normal_intersect.push_back(c_line[i_cent + 1]);

	vector<cad3d::Point3D> near_edge_points;
	vector<float> cline_d_mat, cline_m_mat;
	int i_dmat, i_mmat;
	cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
	float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;
	float mt1, ct1, mt2, ct2, d11, d12;
	tan1 = copy_point3d(c_line[i_cent - 1]);
	tan2 = copy_point3d(c_line[i_cent + 1]);
	norm0 = copy_point3d(c_line[i_cent]);
	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
	m_n = -1 / m_t;
	c_n = -m_n*norm0.x + norm0.y;

	std::cout << "tangent values (m_t, m_n): " << m_t << ", " << m_n << endl;
	find_close_points_from_edge(interpolated_edge_v1, near_edge_points, c_line[i_cent], tol);
	//newblock.update_point3d_vector(near_edge_points, normal_intersect);
	for (int i = 0; i < near_edge_points.size(); i++) {
		mt1 = (near_edge_points[i].y - c_line[i_cent].y) / (near_edge_points[i].x - c_line[i_cent].x);
		//std::cout << i << ": "<<mt1 << ", " << newblock.dist201(near_edge_points[i], c_line[i_cent])<<endl;
		d11 = dist201(near_edge_points[i], c_line[i_cent]);
		cline_d_mat.push_back(d11);
		cline_m_mat.push_back(mt1);
	}
	i_dmat = find_1d_vector_min(cline_d_mat, 1);
	normal_intersect.push_back(near_edge_points[i_dmat]);
	near_edge_points.resize(0);
	cline_d_mat.resize(0);
	cline_m_mat.resize(0);
	
	normal_intersect.push_back(c_line[i_cent]);

	find_close_points_from_edge(interpolated_edge_v2, near_edge_points, c_line[i_cent], tol);
	//newblock.update_point3d_vector(near_edge_points, normal_intersect);
	for (int i = 0; i < near_edge_points.size(); i++) {
		mt1 = (near_edge_points[i].y - c_line[i_cent].y) / (near_edge_points[i].x - c_line[i_cent].x);
		//std::cout << i << ": "<<mt1 << ", " << newblock.dist201(near_edge_points[i], c_line[i_cent])<<endl;
		d11 = dist201(near_edge_points[i], c_line[i_cent]);
		cline_d_mat.push_back(d11);
		cline_m_mat.push_back(mt1);
	}
	i_dmat = find_1d_vector_min(cline_d_mat, 1);
	normal_intersect.push_back(near_edge_points[i_dmat]);
	near_edge_points.resize(0);
	cline_d_mat.resize(0);
	cline_m_mat.resize(0);
	
}
void cad3d::generate_single_cross_along_centreline_v6(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line, vector<Point3D> &normal_intersect, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, float tol, int mode) {
	//created on 22.01.2020
	//updated on 23.01.2020
	//updated on 03.02.2020
	
	vector<cad3d::Point3D> near_edge_points;
	int i_cent = 0;
	std::cin >> i_cent;
	
	//polar mode
	if (mode == 2) {
		vector<float> cline_d_mat, cline_m_mat;
		float mt1, ct1, mt2, ct2, d11, d12;
		int i_dmat, i_mmat;
		cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
		cad3d::Point3D intersect_p1, intersect_p2;
		float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;

		tan1 = copy_point3d(c_line[i_cent - 1]);
		tan2 = copy_point3d(c_line[i_cent + 1]);
		norm0 = copy_point3d(c_line[i_cent]);
		get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
		m_n = -1 / m_t;
		c_n = -m_n*norm0.x + norm0.y;

		std::cout << "tangent values (m_t, m_n): " << m_t << ", " << m_n << endl;

		float edge_d = dist201(interpolated_edge_v1[0], interpolated_edge_v2[0]);
		//std::cout << "distance" << edge_d << endl;


		find_close_points_from_edge(interpolated_edge_v1, near_edge_points, c_line[i_cent], 10);
		std::cout << "near edge points from edge 1: " << near_edge_points.size() << endl;
		for (int i = 0; i < near_edge_points.size(); i++) {
			mt1 = (near_edge_points[i].y - c_line[i_cent].y) / (near_edge_points[i].x - c_line[i_cent].x);
			std::cout << i << ": "<<mt1 << ", " << dist201(near_edge_points[i], c_line[i_cent])<<endl;
			d11 = dist201(near_edge_points[i], c_line[i_cent]);
			cline_d_mat.push_back(d11);
			cline_m_mat.push_back(mt1);
		}
		i_dmat = find_1d_vector_min(cline_d_mat, 1);
		std::cout << "i_dmat" << i_dmat << endl;
		//mt1 = (near_edge_points[i_dmat].y - c_line[i_cent].y) / (near_edge_points[i_dmat].x - c_line[i_cent].x);
		if (near_edge_points.size()>0) {
			get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[i_dmat], near_edge_points[i_dmat]);
			//intersect_p1 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt1, ct1);
			intersect_p1.x = norm0.x + 50;
			intersect_p1.y = m_n*intersect_p1.x + c_n;
			intersect_p1.z = norm0.z;
			//std::cout << near_edge_points.size() << ", " << intersect_p1.x << ", " << intersect_p1.y << endl;
			
		}
		//normal_intersect.push_back(near_edge_points[i_dmat]);
		//intersect_l1.push_back(near_edge_points[i_dmat]);
		normal_intersect.push_back(intersect_p1);
		intersect_l1.push_back(intersect_p1);
		intersect_l2.push_back(c_line[i_cent]);
		//newblock.update_point3d_vector(near_edge_points, normal_intersect);
		near_edge_points.resize(0);
		cline_d_mat.resize(0);
		cline_m_mat.resize(0);

		find_close_points_from_edge(interpolated_edge_v2, near_edge_points, c_line[i_cent], 10);
		std::cout << "near edge points from edge 2: " << near_edge_points.size() << endl;
		for (int i = 0; i < near_edge_points.size(); i++) {
			mt1 = (near_edge_points[i].y - c_line[i_cent].y) / (near_edge_points[i].x - c_line[i_cent].x);
			std::cout << i << ": "<<mt1 << ", " << dist201(near_edge_points[i], c_line[i_cent])<<endl;
			d11 = dist201(near_edge_points[i], c_line[i_cent]);
			cline_d_mat.push_back(d11);
			cline_m_mat.push_back(mt1);
		}
		i_dmat = find_1d_vector_min(cline_d_mat, 1);
		if (near_edge_points.size()>0) {
			get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
			intersect_p1 = get_intersect_point_phase3_v1(norm0, m_n, c_n, mt1, ct1);
			//std::cout << near_edge_points.size() << ", " << intersect_p1.x << ", " << intersect_p1.y << endl;
		}
		std::cout << "i_dmat" << i_dmat << endl;
		//normal_intersect.push_back(near_edge_points[i_dmat]);
		
		//newblock.update_point3d_vector(near_edge_points, normal_intersect);
		normal_intersect.push_back(c_line[i_cent]);
		normal_intersect.push_back(intersect_p1);
		intersect_l1.push_back(c_line[i_cent]);
		intersect_l2.push_back(intersect_p1);
		//intersect_l2.push_back(near_edge_points[i_dmat]);
	}
}
void cad3d::generate_single_cross_along_centreline_v7(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line, vector<Point3D> &normal_intersect, vector<Point3D> &intersect_line1, vector<Point3D> &intersect_line2, float tol, int mode) {
	//created on 04.02.2020
	//updated on 23.01.2020
	vector<cad3d::Point3D> near_edge_points;
	int i_cent = 0;
	std::cout << "enter path index: " ;
	std::cin >> i_cent;

	//updated on 03.02.2020
	//polar mode
	if (mode == 2) {
		intersect_line1.push_back(interpolated_edge_v1[0]);
		intersect_line2.push_back(c_line[0]);
		intersect_line1.push_back(c_line[0]);
		intersect_line2.push_back(interpolated_edge_v2[0]);
		cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
		float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;
		vector<float> cline_mn_mat, cline_mt_mat, del_mt_mat;

		tan1 = copy_point3d(c_line[0]);
		tan2 = copy_point3d(c_line[1]);
		norm0 = copy_point3d(c_line[0]);

		get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
		m_n = -1 / m_t;
		c_n = -m_n*norm0.x + norm0.y;
		//std::cout << "0:  (m_t, m_n) at start => " << m_t << ", " << m_n << endl;
		cline_mt_mat.push_back(m_t);

		tan1 = copy_point3d(interpolated_edge_v1[0]);
		tan2 = copy_point3d(interpolated_edge_v2[0]);
		norm0 = copy_point3d(c_line[0]);
		get_straight_line_2d(m_n, c_n, tan1, tan2, norm0);
		//std::cout << "intersection values at start (m_n) at start: " << m_n << endl;
		cline_mn_mat.push_back(m_n);

		float mt1, ct1, mt2, ct2, d11, d12, del_mt;
		vector<float> cline_d_mat, cline_m_mat;
		int i_dmat, i_mmat;
		cad3d::Point3D intersect_p1, intersect_p2;

		for (int i = 1; i < c_line.size()-1; i++) {
		//for (int i = 1; i < 6; i++) {
			
			//measure tangent and normal along centerline
			tan1 = copy_point3d(c_line[i - 1]);
			tan2 = copy_point3d(c_line[i]);
			norm0 = copy_point3d(c_line[i]);
			get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
			//m_n = -1 / m_t;
			//c_n = -m_n*norm0.x + norm0.y;
			cline_mt_mat.push_back(m_t);
			del_mt = cline_mt_mat[i] - cline_mt_mat[i - 1];
			del_mt_mat.push_back(del_mt);
			//std::cout << i << ": " << m_t << ", " << del_mt << endl;
			cline_mn_mat.push_back(cline_mn_mat[i - 1] + del_mt);
			//cline_mn_mat.push_back(cline_mn_mat[i]- del_mt_mat[i - 1]);
			c_n = -cline_mn_mat[i] * c_line[i].x + c_line[i].y;

			//find line for edge 1
			find_close_points_from_edge(interpolated_edge_v1, near_edge_points, c_line[i], 20);
			//std::cout << i <<": near edge points from edge 1 => " << near_edge_points.size() << endl;
			for (int j = 0; j < near_edge_points.size(); j++) {
				mt1 = (near_edge_points[j].y - c_line[i].y) / (near_edge_points[j].x - c_line[i].x);
				//std::cout << i << ": " << mt1 << ", " << dist201(near_edge_points[i], c_line[i_cent]) << endl;
				d11 = dist201(near_edge_points[j], c_line[i]);
				cline_d_mat.push_back(d11);
				cline_m_mat.push_back(mt1);
			}
			i_dmat = find_1d_vector_min(cline_d_mat, 1);
			if (near_edge_points.size()>0) {
				get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[i_dmat], near_edge_points[i_dmat]);
				intersect_p1 = get_intersect_point_phase3_v1(c_line[i], cline_mn_mat[i], c_n, mt1, ct1);
				//intersect_p1.x = norm0.x + 50;
				//intersect_p1.y = m_n*intersect_p1.x + c_n;
				intersect_p1.z = norm0.z;
				//std::cout << "intersect_p1 " << intersect_p1.x << ", " << intersect_p1.y << "~ ";

			}
			
			near_edge_points.resize(0);
			cline_d_mat.resize(0);
			cline_m_mat.resize(0);

			//find line for edge 2
			find_close_points_from_edge(interpolated_edge_v2, near_edge_points, c_line[i], 30);
			//std::cout << "near edge points from edge 1: " << near_edge_points.size() << endl;
			for (int j = 0; j < near_edge_points.size(); j++) {
				mt1 = (near_edge_points[j].y - c_line[i].y) / (near_edge_points[j].x - c_line[i].x);
				//std::cout << i << ": " << mt1 << ", " << dist201(near_edge_points[i], c_line[i_cent]) << endl;
				d11 = dist201(near_edge_points[j], c_line[i]);
				cline_d_mat.push_back(d11);
				cline_m_mat.push_back(mt1);
			}
			i_dmat = find_1d_vector_min(cline_d_mat, 1);
			if (near_edge_points.size()>0) {
				get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[i_dmat], near_edge_points[i_dmat]);
				intersect_p2 = get_intersect_point_phase3_v1(c_line[i], cline_mn_mat[i], c_n, mt1, ct1);
				//intersect_p1.x = norm0.x + 50;
				//intersect_p1.y = m_n*intersect_p1.x + c_n;
				intersect_p2.z = norm0.z;
				//std::cout << "intersect_p2 " << intersect_p2.x << ", " << intersect_p2.y << endl;

			}
			
			if (i == i_cent) {
				if (!isnan(intersect_p1.x) && !isnan(intersect_p1.y)) {
					if (!isnan(intersect_p2.x) && !isnan(intersect_p2.y)) {
						normal_intersect.push_back(intersect_p1);
						intersect_line1.push_back(intersect_p1);
						intersect_line2.push_back(c_line[i]);
						normal_intersect.push_back(intersect_p2);
						intersect_line1.push_back(intersect_p2);
						intersect_line2.push_back(c_line[i]);
					}
				}
			}
			
			near_edge_points.resize(0);
			cline_d_mat.resize(0);
			cline_m_mat.resize(0);
		}
	}
}
void cad3d::multiple_layer_cross_section_extract(stl_data &block_a, vector<Point3D> &sub_vol_points, vector<Point3D> &c_line, vector<Point3D> &normal_intersect, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2) {
	//updated on 13.10.2020
	float z_now, z_shift;
	int extract_status = 1;
	vector<cad3d::Point3D> sub_vol_slice, cline_temp,t_edge1, t_edge2, t_interp1,t_interp2;
	
	while (extract_status ==1) {
		adaptive_2p5d_decomposition(sub_vol_points, sub_vol_slice, z_now, 3, 20, 0.5);
		
		t_edge1.resize(0);	t_edge2.resize(0);
		t_interp1.resize(0); t_interp2.resize(0);
		cline_temp.resize(0);
		single_slot_seperate_edges_v103(block_a, sub_vol_slice, t_edge1, t_edge2);
		if (t_edge1.size()>0) {
			interpolate_edge_if_necessary_v7(t_edge1, t_interp1, 4.0);
			interpolate_edge_if_necessary_v7(t_edge2, t_interp2, 4.0);
			//interpolate_edge_if_necessary_v8(t_edge1, t_interp1, 2.0);
			//interpolate_edge_if_necessary_v8(t_edge2, t_interp2, 2.0);
			z_shift = t_edge1[0].z;
			for (int i = 0; i < c_line.size(); i++) {
				cline_temp.push_back(c_line[i]);
				cline_temp[i].z = z_shift;
			}

			generate_single_cross_along_centreline_v2(t_interp1, t_interp2, cline_temp, normal_intersect);
			//generate_single_cross_along_centreline_v3(t_interp1, t_interp2, cline_temp, intersect_l1, intersect_l2);
			generate_single_cross_along_centreline_v4(t_interp1, t_interp2, cline_temp, intersect_l1, intersect_l2);
			float z_shift = t_edge1[0].z;
		}
	
		std::cout <<"Press 1 to extract another intersection layer and 0 to none: ";
		std::cin >> extract_status;
		sub_vol_slice.resize(0);
	}


}
void cad3d::cross_section_outline_generate_v1(vector<Point3D> &vin1, vector<Point3D> &vin2, vector<Point3D> &vout1, vector<Point3D> &vout2) {
	Point3D ptemp;
	ptemp.x = 0.0;
	ptemp.y = 0;
	ptemp.z = 0.0;
	if (vin1.size() == vin2.size()) {
		for (int i = 0; i < vin1.size()-1; i++) {
			vout1.push_back(vin1[i]);
			vout2.push_back(vin1[i+1]);
			vout1.push_back(vin2[i]);
			vout2.push_back(vin2[i + 1]);
		}
		vout1.push_back(vin1[vin1.size() - 1]);
		vout2.push_back(vin2[vin2.size() - 1]);
	}
	else {
		vout1.push_back(ptemp);
		vout2.push_back(ptemp);
	}
	

}

void cad3d::fine_3D_segment_extraction_polar(vector<Point3D> &c_line, int i_cent, vector<Point3D> &slot_e1) {
	vector<cad3d::Point3D> normal_intersect, near_edge_points;
	cad3d::Point3D tan1, tan2, norm1, norm2, norm0;
	float dy, dx, m_t, c_t, m_n, c_n, slot_width, dtemp;
	float mt1, ct1, mt2, ct2;
	tan1 = copy_point3d(c_line[i_cent - 1]);
	tan2 = copy_point3d(c_line[i_cent + 1]);
	norm0 = copy_point3d(c_line[i_cent]);
	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);
	m_n = -1 / m_t;
	c_n = -m_n*norm0.x + norm0.y;
	float d_temp;
	find_close_points_from_edge(slot_e1, near_edge_points, c_line[i_cent], 60);
}
void cad3d::generate_3D_sub_layer_v1(vector<vector<Point3D>> &sub_vol_main_slices, vector<Point3D> &boundary_3d_line1, vector<Point3D> &boundary_3d_line2, int mode) {
	
	// calculate sub volume border lines, Omega_lines
	//created on 21.02.2020
	//updated on ...

	cad3d::Point3D v_temp1;
	cad3d::Point3D v_temp2;
	vector<cad3d::Point3D> present_slice, next_slice;
	float p2p, del_plane_z;
	int close_point_stat = 0;
	for (int i = 0; i < sub_vol_main_slices.size() - 1; i++) {
		close_point_stat == 0;
		std::cout << "calculating  3D sub layers: " << i << "=>" << sub_vol_main_slices[i].size() << endl;
		del_plane_z = sub_vol_main_slices[i][1].z - sub_vol_main_slices[i + 1][1].z;
		std::cout << del_plane_z << "enter 0 to continue";
		std::cin >> close_point_stat;
		if (close_point_stat == 1) {
			break;
		}
		for (int jj = 0; jj < sub_vol_main_slices[i].size(); jj++) {
			present_slice.push_back(sub_vol_main_slices[i][jj]);
			v_temp1 = copy_point3d(present_slice[jj]);
			
			
			//continue if same plane
			
			/*if (abs(del_plane_z) < 0.1) {
				std::cout << "not counting for points on the same plane" << endl;
				continue;
			}*/
			// find distance
			std::cout << jj << "=>" << endl;
			del_plane_z = v_temp1.z - sub_vol_main_slices[i + 1][1].z;
			for (int k = 0; k < sub_vol_main_slices[i + 1].size(); k++) {
				v_temp2.x = sub_vol_main_slices[i + 1][k].x;
				v_temp2.y = sub_vol_main_slices[i + 1][k].y;
				v_temp2.z = sub_vol_main_slices[i + 1][k].z;
				p2p = dist201(v_temp1, v_temp2);
				
				//std::cout << i << "---" << abs(del_plane_z) << "---";
				if (p2p < 1) {
					std::cout << k << ": " << p2p << "~" << endl;
					close_point_stat = 1;
					boundary_3d_line1.push_back(v_temp1);
					boundary_3d_line2.push_back(v_temp2);
					break;
				}
				/*if (abs(del_plane_z) > 1) {
					//std::cout << "not counting for points on the same plane" << endl;
					//continue;
					
				}
				else {
					std::cout << "not counting for points on the same plane" << endl;
					close_point_stat = 1;
				}*/
				
			}
			
		}
		//std::cout << present_slice.size() << endl;
	}
}
void cad3d::generate_3D_sub_layer_v2(vector<vector<Point3D>> &sub_vol_main_slices, vector<vector<Point3D>> &sub_vol_sub_slices, vector<Point3D> &boundary_3d_line1, vector<Point3D> &boundary_3d_line2, int mode) {

	// calculate sub volume border lines, Omega_lines
	//created on 21.02.2020
	//updated on 25.02.2020

	cad3d::Point3D v_temp1;
	cad3d::Point3D v_temp2;
	cad3d::Point3D v_sub;
	vector<cad3d::Point3D> present_slice, next_slice, temp_sub_slice;
	float p2p, del_plane_z, z_present;
	int close_point_stat = 0;
	int z_step;
	
	if (mode == 2) {
		z_step = 5;// total iteration will be step size -1
		for (int i = 0; i < sub_vol_main_slices.size() - 1; i++) {
			close_point_stat == 0;
			std::cout << "calculating  3D sub layers: " << i << "=>" << sub_vol_main_slices[i].size() << endl;
			del_plane_z = sub_vol_main_slices[i][0].z - sub_vol_main_slices[i + 1][0].z;
			std::cout << "diff in z: " << del_plane_z << " ~ enter 0 to continue" << endl;
			std::cin >> close_point_stat;

			if (close_point_stat == 0) {
				// get vertical/ along z axes boundary lines
				for (int jj = 0; jj < sub_vol_main_slices[i].size(); jj++) {
					present_slice.push_back(sub_vol_main_slices[i][jj]);
					v_temp1 = copy_point3d(present_slice[jj]);
					// find distance
					std::cout << jj << "=>" << endl;
					del_plane_z = v_temp1.z - sub_vol_main_slices[i + 1][1].z;
					for (int k = 0; k < sub_vol_main_slices[i + 1].size(); k++) {
						v_temp2.x = sub_vol_main_slices[i + 1][k].x;
						v_temp2.y = sub_vol_main_slices[i + 1][k].y;
						v_temp2.z = sub_vol_main_slices[i + 1][k].z;
						p2p = dist201(v_temp1, v_temp2);
						//std::cout << i << "---" << abs(del_plane_z) << "---";
						if (p2p < 1) {
							std::cout << k << ": " << p2p << "~" << endl;
							close_point_stat = 1;
							boundary_3d_line1.push_back(v_temp1);
							boundary_3d_line2.push_back(v_temp2);
							break;
						}
					}
				}
				// split in sub volume sub borderpoints
				std::cout << "sub step size along Z: 4 ";
				//std::cin >> z_step;
				for (int i1 = 0; i1<4; i1++) {
					z_present = sub_vol_main_slices[i][1].z - (i1 + 1)*del_plane_z / z_step;
					for (int i2 = 0; i2 < boundary_3d_line1.size(); i2++) {
						v_sub.x = (boundary_3d_line1[i2].x + boundary_3d_line2[i2].x) / 2;
						v_sub.y = (boundary_3d_line1[i2].y + boundary_3d_line2[i2].y) / 2;
						v_sub.z = z_present;
						temp_sub_slice.push_back(v_sub);
					}
					sub_vol_sub_slices.push_back(temp_sub_slice);
					std::cout << "temp_sub_slice size:" << temp_sub_slice.size() << endl;
					temp_sub_slice.resize(0);
				}
				//std::cout << present_slice.size() << endl;
				present_slice.resize(0);
			}
		}
	}
	else {
		for (int i = 0; i < sub_vol_main_slices.size() - 1; i++) {
			close_point_stat == 0;
			std::cout << "calculating  3D sub layers: " << i << "=>" << sub_vol_main_slices[i].size() << endl;
			del_plane_z = sub_vol_main_slices[i][1].z - sub_vol_main_slices[i + 1][1].z;
			std::cout << "diff in z: " << del_plane_z << " ~ enter 0 to continue";
			std::cin >> close_point_stat;

			if (close_point_stat == 1) {
				break;
			}
			// get vertical/ along z axes boundary lines
			for (int jj = 0; jj < sub_vol_main_slices[i].size(); jj++) {
				present_slice.push_back(sub_vol_main_slices[i][jj]);
				v_temp1 = copy_point3d(present_slice[jj]);
				// find distance
				std::cout << jj << "=>" << endl;
				del_plane_z = v_temp1.z - sub_vol_main_slices[i + 1][1].z;
				for (int k = 0; k < sub_vol_main_slices[i + 1].size(); k++) {
					v_temp2.x = sub_vol_main_slices[i + 1][k].x;
					v_temp2.y = sub_vol_main_slices[i + 1][k].y;
					v_temp2.z = sub_vol_main_slices[i + 1][k].z;
					p2p = dist201(v_temp1, v_temp2);
					//std::cout << i << "---" << abs(del_plane_z) << "---";
					if (p2p < 1) {
						std::cout << k << ": " << p2p << "~" << endl;
						close_point_stat = 1;
						boundary_3d_line1.push_back(v_temp1);
						boundary_3d_line2.push_back(v_temp2);
						break;
					}
				}
			}

			// split in sub volume sub borderpoints
			std::cout << "sub step size along Z: 4 ";
			//std::cin >> z_step;

			//std::cout << present_slice.size() << endl;
		}
		z_step = 5;// total iteration will be step size -1
		for (int i1 = 0; i1<4; i1++) {
			z_present = sub_vol_main_slices[0][1].z - (i1 + 1)*del_plane_z / z_step;
			for (int i2 = 0; i2 < boundary_3d_line1.size(); i2++) {
				v_sub.x = (boundary_3d_line1[i2].x + boundary_3d_line2[i2].x) / 2;
				v_sub.y = (boundary_3d_line1[i2].y + boundary_3d_line2[i2].y) / 2;
				v_sub.z = z_present;
				temp_sub_slice.push_back(v_sub);
			}
			sub_vol_sub_slices.push_back(temp_sub_slice);
			std::cout << "temp_sub_slice size:" << temp_sub_slice.size() << endl;
			temp_sub_slice.resize(0);
		}
	}
}

void cad3d::generate_3D_fine_slice_v1(vector<vector<Point3D>> &sub_vol_main_slices, vector<vector<Point3D>> &sub_vol_sub_slices, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, vector<Point3D> &slice3d_l1, vector<Point3D> &slice3d_l2, int mode) {
	//created on 26.02.2020
	std::cout << "Generating fine volume segment in 3D" << endl;

	//step 1: reorder slices
	// for the same level, compare with enclosing rectangular boundary box

	//step 2: get projection line on all slices from intersect line 
}
void cad3d::generate_3D_fine_slice_v2(vector<vector<Point3D>> &sub_vol_slices, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, vector<Point3D> &slice3d_l1, vector<Point3D> &slice3d_l2, int mode) {
	//created on 26.02.2020
	std::cout << "Generating fine volume segment in 3D" << endl;

	//step 1: reorder slices
	// for the same level, compare with enclosing rectangular boundary box

	//step 2: get projection line on all slices from intersect line 
	int itemp1, itemp2;
	vector<float> d_mat1, d_mat2;
	vector<Point3D> slice3d_v1, slice3d_v2;
	vector<Point3D> int_slice3d_l1, int_slice3d_l2;
	Point3D ext_v1, ext_v2;
	Point3D intersect_p1, intersect_p2;
	float m_int, m1, m2, c_int, c1, c2, crit_dist;
	
	for (int i = 0; i < sub_vol_slices.size(); i++) {
		std::cout << i << endl << endl;
		ext_v1 = copy_point3d(intersect_l1[0]);
		ext_v1.z = sub_vol_slices[i][0].z;
		ext_v2 = copy_point3d(intersect_l2[0]);
		ext_v2.z = sub_vol_slices[i][0].z;
		get_straight_line_2d(m_int, c_int, ext_v1, ext_v2, ext_v1);

		for (int j = 0; j < sub_vol_slices[i].size(); j++) {
			d_mat1.push_back(dist301(ext_v1, sub_vol_slices[i][j]));
			d_mat2.push_back(dist301(ext_v2, sub_vol_slices[i][j]));
			std::cout << j << ": " << d_mat1[j] << "~ " << d_mat2[j] << endl;

		}
		itemp1 = find_1d_vector_min(d_mat1, 1);
		itemp2 = find_1d_vector_min(d_mat2, 1);
		if (mode == 1) { crit_dist = 2.0; }
		if (mode == 2) { crit_dist = 1.1; }//1.1
		if (mode == 3) { crit_dist = 1.3; }
		for (int j = 0; j < sub_vol_slices[i].size(); j++) {
			if ((d_mat1[j] < crit_dist * d_mat1[itemp1])) {
				slice3d_v1.push_back(sub_vol_slices[i][j]);
			}
			if (d_mat2[j] < crit_dist * d_mat2[itemp2]) {
				slice3d_v2.push_back(sub_vol_slices[i][j]);
			}
		}
		//get_intersect_point_phase3_v1(Point3D norm0, float m1, float c1, float m2, float c2)
		get_straight_line_2d(m1, c1, slice3d_v1[0], slice3d_v1[slice3d_v1.size() - 1], slice3d_v1[0]);
		intersect_p1 = get_intersect_point_phase3_v1(ext_v1, m_int, c_int, m1, c1);
		intersect_p1.z = ext_v1.z;
		int_slice3d_l1.push_back(intersect_p1);
		std::cout << "intersect 1:" << intersect_p1.x << "," << intersect_p1.y << endl;

		get_straight_line_2d(m2, c2, slice3d_v2[0], slice3d_v2[slice3d_v2.size() - 1], slice3d_v2[0]);
		intersect_p2 = get_intersect_point_phase3_v1(ext_v2, m_int, c_int, m2, c2);
		intersect_p2.z = ext_v1.z;
		int_slice3d_l2.push_back(intersect_p2);
		std::cout << "intersect 2:" << intersect_p2.x << "," << intersect_p2.y << endl;

		//update_point3d_vector(slice3d_v1, slice3d_l1);
		//update_point3d_vector(slice3d_v2, slice3d_l2);
		slice3d_v1.resize(0);
		slice3d_v2.resize(0);
		d_mat1.resize(0);
		d_mat2.resize(0);
	}
	
	//reorder slice 3d as intersection 
	for (int i = 0; i < int_slice3d_l1.size()-1; i++) {
		
		slice3d_l1.push_back(int_slice3d_l1[i]);
		slice3d_l2.push_back(int_slice3d_l1[i+1]);

		slice3d_l1.push_back(int_slice3d_l2[i]);
		slice3d_l2.push_back(int_slice3d_l2[i + 1]);
	}
	slice3d_l1.push_back(int_slice3d_l1[int_slice3d_l1.size()-1]);
	slice3d_l2.push_back(int_slice3d_l2[int_slice3d_l2.size() - 1]);
}
void cad3d::generate_3D_fine_slice_v3(vector<vector<Point3D>> &sub_vol_slices, vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, vector<Point3D> &slice3d_l1, vector<Point3D> &slice3d_l2, int mode) {
	//created on 13.10.2020
	std::cout << "Generating fine volume segment in 3D" << endl;

	//step 1: reorder slices
	// for the same level, compare with enclosing rectangular boundary box

	//step 2: get projection line on all slices from intersect line 
	int itemp1, itemp2;
	vector<float> d_mat1, d_mat2;
	vector<Point3D> slice3d_v1, slice3d_v2;
	vector<Point3D> int_slice3d_l1, int_slice3d_l2;
	Point3D ext_v1, ext_v2;
	Point3D intersect_p1, intersect_p2;
	float m_int, m1, m2, c_int, c1, c2, crit_dist;

	//vector<Point3D> interpolated_slice;
	std::cout << "sub-volume size: " << sub_vol_slices.size() << endl;
	for (int i = 0; i < sub_vol_slices.size(); i++) {
		std::cout << i <<  endl;
		//interpolate_edge_if_necessary_v8(t_edge1, t_interp1, 2.0);
		//interpolate_edge_if_necessary_v8(sub_vol_slices[i], interpolated_slice, 2.0);
		ext_v1 = copy_point3d(intersect_l1[0]);
		ext_v1.z = sub_vol_slices[i][0].z;
		ext_v2 = copy_point3d(intersect_l2[0]);
		ext_v2.z = sub_vol_slices[i][0].z;
		get_straight_line_2d(m_int, c_int, ext_v1, ext_v2, ext_v1);
		//std::cout << "done with starting points" << endl;
		for (int j = 0; j < sub_vol_slices[i].size(); j++) {
			d_mat1.push_back(dist301(ext_v1, sub_vol_slices[i][j]));
			d_mat2.push_back(dist301(ext_v2, sub_vol_slices[i][j]));
			std::cout << j << ": " << d_mat1[j] << "~ " << d_mat2[j] << endl;

		}
		//std::cout << "done with phase 2" << endl;
		itemp1 = find_1d_vector_min(d_mat1, 1);
		itemp2 = find_1d_vector_min(d_mat2, 1);
		if (mode == 1) { crit_dist = 20.00; }//2.0
		if (mode == 2) { crit_dist = 10.00; }//1.1
		if (mode == 3) { crit_dist = 1.5; }
		for (int j = 0; j < sub_vol_slices[i].size(); j++) {
			//if ((d_mat1[j] < crit_dist * d_mat1[itemp1])) {
			if ((d_mat1[j] < crit_dist )) {
				slice3d_v1.push_back(sub_vol_slices[i][j]);
			}
			//if (d_mat2[j] < crit_dist * d_mat2[itemp2]) {
			if ((d_mat2[j] < crit_dist )) {
				slice3d_v2.push_back(sub_vol_slices[i][j]);
			}
		}
		std::cout << "mode and crit_dist" <<mode << ", "<<crit_dist<< endl;
		//get_intersect_point_phase3_v1(Point3D norm0, float m1, float c1, float m2, float c2)
		std::cout << "slice_3d sizes: " << slice3d_v1.size() << ", " << slice3d_v2.size()<< endl;
		get_straight_line_2d(m1, c1, slice3d_v1[0], slice3d_v1[slice3d_v1.size() - 1], slice3d_v1[0]);
		intersect_p1 = get_intersect_point_phase3_v1(ext_v1, m_int, c_int, m1, c1);
		intersect_p1.z = ext_v1.z;
		int_slice3d_l1.push_back(intersect_p1);
		std::cout << "intersect 1:" << intersect_p1.x << "," << intersect_p1.y << endl;

		get_straight_line_2d(m2, c2, slice3d_v2[0], slice3d_v2[slice3d_v2.size() - 1], slice3d_v2[0]);
		intersect_p2 = get_intersect_point_phase3_v1(ext_v2, m_int, c_int, m2, c2);
		intersect_p2.z = ext_v1.z;
		int_slice3d_l2.push_back(intersect_p2);
		std::cout << "intersect 2:" << intersect_p2.x << "," << intersect_p2.y << endl;

		//update_point3d_vector(slice3d_v1, slice3d_l1);
		//update_point3d_vector(slice3d_v2, slice3d_l2);
		//interpolated_slice.resize(0);
		slice3d_v1.resize(0);
		slice3d_v2.resize(0);
		d_mat1.resize(0);
		d_mat2.resize(0);
	}

	//reorder slice 3d as intersection 
	for (int i = 0; i < int_slice3d_l1.size() - 1; i++) {

		slice3d_l1.push_back(int_slice3d_l1[i]);
		slice3d_l2.push_back(int_slice3d_l1[i + 1]);

		slice3d_l1.push_back(int_slice3d_l2[i]);
		slice3d_l2.push_back(int_slice3d_l2[i + 1]);
	}
	slice3d_l1.push_back(int_slice3d_l1[int_slice3d_l1.size() - 1]);
	slice3d_l2.push_back(int_slice3d_l2[int_slice3d_l2.size() - 1]);
}

/*newblock.interpolate_edge_if_necessary_v5(slot_e1, slot_e2, slot_ip1, interpolated_edge_v1);
				newblock.interpolate_edge_if_necessary_v5(slot_e3, slot_e4, slot_ip2, interpolated_edge_v2);*/
void cad3d::generate_extended_projection_line_v1(vector<Point3D> &intersect_l1, vector<Point3D> &intersect_l2, vector<Point3D> &ext_intersect_l1, vector<Point3D> &ext_intersect_l2, float stretch_val) {
	//created on 26.02.2020
	//updated on 27.02.2020
	float delx, dely;
	delx = intersect_l2[0].x - intersect_l1[0].x;
	dely = intersect_l2[0].y - intersect_l1[0].y;
	cad3d::Point3D v_temp1;
	cad3d::Point3D v_temp2;
	std::cout << "Generating extended projection line: " <<delx<<", "<<dely<< endl;
	//measure m and c
	cad3d::Point3D tan1, tan2, norm0;
	float m_t, c_t;
	tan1 = copy_point3d(intersect_l1[intersect_l1.size() - 1]);
	tan2 = copy_point3d(intersect_l2[intersect_l1.size() - 1]);
	norm0 = copy_point3d(intersect_l1[intersect_l1.size() - 1]);
	get_straight_line_2d(m_t, c_t, tan1, tan2, norm0);

	if (abs(dely) >= abs(delx)) {
		std::cout << "abs(dely) >= abs(delx)" << endl;
		if (dely >=0) {
			v_temp1.y = intersect_l1[intersect_l1.size()-1].y - stretch_val*dely;
			v_temp1.x = (v_temp1.y - c_t) / m_t;
			v_temp1.z = intersect_l1[intersect_l1.size() - 1].z;

			v_temp2.y = intersect_l2[intersect_l1.size() - 1].y + stretch_val*dely;
			v_temp2.x = (v_temp2.y - c_t) / m_t;
			v_temp2.z = intersect_l2[intersect_l1.size() - 1].z;
		}
		else
		{
			v_temp1.y = intersect_l1[intersect_l1.size() - 1].y + stretch_val*dely;
			v_temp1.x = (v_temp1.y - c_t) / m_t;
			v_temp1.z = intersect_l1[intersect_l1.size() - 1].z;

			v_temp2.y = intersect_l2[intersect_l1.size() - 1].y - stretch_val*dely;
			v_temp2.x = (v_temp2.y - c_t) / m_t;
			v_temp2.z = intersect_l2[intersect_l1.size() - 1].z;
		}
	}
	else {
		//do nothing for now
	}
	ext_intersect_l1.push_back(v_temp1);
	ext_intersect_l2.push_back(v_temp2);
	ext_intersect_l1.push_back(v_temp1);
	ext_intersect_l2.push_back(v_temp2);
}

/*void cad3d::generate_single_cross_along_centreline_v7(vector<Point3D> &interpolated_edge_v1, vector<Point3D> &interpolated_edge_v2, vector<Point3D> &c_line, vector<Point3D> &normal_intersect, vector<Point3D> &intersect_line1, vector<Point3D> &intersect_line2, float tol, int mode) {
	//created on 04.02.2020
	//updated on 23.01.2020
	vector<cad3d::Point3D> near_edge_points;
	}*/
/*void commented_lines_fine_segment_polar() {
	if (slot_type == 3) {
		
		//ptemp3 = newblock.find_correlated_point(c_line, v1, c_line, i_cent);
		//normal_intersect.push_back(ptemp3);
		//cad3d::Point3D o_temp, c_temp, temp0;
		//o_temp.x = 0.0; o_temp.y = 0.0; o_temp.z = 0.0;
		//int i_intersect;
		//i_intersect = newblock.calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, c_line, i_cent);
		//normal_intersect.push_back(interpolated_edge_v1[i_intersect]);

		
		//newblock.get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		//std::cout << "slot width: " << newblock.dist201(interpolated_edge_v1[0],interpolated_edge_v2[0])<<endl;
		//newblock.update_point3d_vector(near_edge_points,normal_intersect);
		//find_close_points_from_edge(edge1, near_edge_points, c_line[i_cent], slot_width);
		//from edge 1
		cad3d::Point3D intersect_p1;*/
		/*for (int i = 0; i < v1.size(); i++) {
		d_temp = newblock.dist201(c_line[i_cent], v1[i]);
 
		if (d_temp<80) {
		std::cout << i << ": " << d_temp << endl;
		normal_intersect.push_back(v1[i]);
		}
		}
		*/
		/*cad3d::Point3D ptemp3;
		newblock.find_close_points_from_edge(interpolated_edge_v1, near_edge_points, c_line[i_cent], 14.0);
		newblock.get_straight_line_2d(mt1, ct1, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);

		intersect_line1.push_back(c_line[i_cent]);
		intersect_line2.push_back(near_edge_points[0]);
		near_edge_points.resize(0);
		//target_plot3_update(near_edge_points);
		//target_plot3_color_update(0.2, 0.2, 0.2, near_edge_points.size());


		newblock.find_close_points_from_edge(interpolated_edge_v2, near_edge_points, c_line[i_cent], 14.0);
		newblock.get_straight_line_2d(mt2, ct2, near_edge_points[0], near_edge_points[near_edge_points.size() - 1], near_edge_points[0]);
		//target_plot3_update(near_edge_points);
		//target_plot3_color_update(0.2, 0.2, 0.2, near_edge_points.size());
		intersect_line1.push_back(c_line[i_cent]);
		intersect_line2.push_back(near_edge_points[0]);
		near_edge_points.resize(0);
		near_edge_points.resize(0);

		single_cs.push_back(interpolated_edge_v1[i_cent - 1]);
		//intersect_line1.push_back(interpolated_edge_v1[i_cent - 1]);
		//newblock.generate_single_cross_along_centreline_v1(block_a, slot_plane1, slot_ip1[single_cross], non_planar_point_stat, single_cs, ti_all2, slot_e1, single_cross);
		ptemp3 = newblock.find_correlated_point(interpolated_edge_v1, interpolated_edge_v2, interpolated_edge_v1, i_cent);
		//intersect_line2.push_back(ptemp3);

		//generate_single_cross_along_centreline_v1(stl_data &block_a, vector<Point3D> &slot_plane1, Point3D pin, int & status_slot_plane1, vector<Point3D> &single_cs, vector<triangle_index> &ti_all2, vector<Point3D> &slot_ip1, int single_cross)
		std::cout << "interpolated sizes: " << interpolated_edge_v1.size() << ", " << interpolated_edge_v2.size() << endl;
		//find the cross correlated point


		//cross_line1.push_back(interpolated_edge_v1[single_cross]);
		//i_temp = newblock.calculate_intersectline_from_interpolated_curves_v1(interpolated_edge_v1, interpolated_edge_v2, single_cross);
		//cross_line2.push_back(interpolated_edge_v2[i_temp]);
}*/

