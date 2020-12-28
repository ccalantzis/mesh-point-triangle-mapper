#ifndef MAPPER_H
#define MAPPER_H

#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <chrono>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writePLY.h>
#include <igl/voxel_grid.h>
#include <igl/barycentric_coordinates.h>
#include <fstream>
#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <igl/point_mesh_squared_distance.h>
#include <chrono>
#include <random>
#include <cstdlib>
#include <ctime>

class mapper
{
public:
    mapper(std::string filename);
    bool test_point_to_triangle(int &test_triangle, Eigen::RowVector3d &test_point);
    void display_point_and_triangle(int test_triangle, Eigen::RowVector3d test_point);
    int point_to_triangle(const Eigen::RowVector3d &point);

private:
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::RowVector3i sides;           // number of cubes on each side of the voxel grid
    double csl;                         // cube side length
    Eigen::MatrixXd voxel_cubes;        // a matrix of pairs containing the min and max vertices of each cube in the voxel grid
    Eigen::RowVector3d vg_min;          // minimum vertex of the voxel grid
    Eigen::RowVector3d vg_max;          // maximum vertex of the voxel grid
    std::vector<std::vector<int>> map;  // mapping of voxels to triangles

    bool is_point_within_range(const Eigen::RowVector3d &point);
    int point_to_voxel(const Eigen::RowVector3d &point);
    Eigen::MatrixXd get_voxel_cubes(const Eigen::MatrixXd &cube_centres);
    std::vector<int> triangle_to_voxels(const Eigen::Matrix<double,3,3> &triangle);
    std::vector<std::vector<int>> voxel_to_triangle_mapping();    
    double randMToN(double M, double N);
    std::string point_to_string(Eigen::Vector3d point);
    void set_vertex_colours(Eigen::MatrixXd &colours, const int &test_triangle);
    void print_voxel_info(std::vector<std::vector<int>> map);
};

#endif // MAPPER_H
