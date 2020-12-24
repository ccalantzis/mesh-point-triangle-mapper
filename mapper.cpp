#include "mapper.h"

mapper::mapper(std::string filename)
{
    if(filename.substr(filename.length()-3) == "obj")
    {
        igl::readOBJ(filename, V, F);
    }
    if(filename.substr(filename.length()-3) == "ply")
    {
        igl::readPLY(filename, V, F);
    }

    Eigen::MatrixXd GV;
    int cubes_longest_side = 80;
    Eigen::AlignedBox<double, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());
    igl::voxel_grid(box, cubes_longest_side, 0, GV, sides);
    csl = (GV.row(0)-GV.row(1)).cwiseAbs().maxCoeff();
    voxel_cubes = get_voxel_cubes(GV);
    vg_min = voxel_cubes.colwise().minCoeff();
    vg_max = voxel_cubes.colwise().maxCoeff();
    map = voxel_to_triangle_mapping();
}

Eigen::MatrixXd mapper::get_voxel_cubes(const Eigen::MatrixXd &cube_centres)
{
    Eigen::MatrixXd cube_V(cube_centres.rows()*2,3);
    double chsl = csl/2;
    for (int i = 0; i < cube_centres.rows(); ++i)
    {
        auto tcv = cube_centres.row(i); //centre of cube
        cube_V.row(2*i) << tcv(0) - chsl, tcv(1) - chsl, tcv(2) - chsl;
        cube_V.row(2*i + 1) << tcv(0) + chsl, tcv(1) + chsl, tcv(2) + chsl;
    }
    return cube_V;
}

// checking that the specified point is inside the bounding box of the object
bool mapper::is_point_within_range(const Eigen::RowVector3d &point)
{
    for(int i = 0; i <3; ++i)
    {
        if(point(i) < vg_min(i) || point(i) > vg_max(i))
        {
            std::cout << point_to_string(point) << " not contained in bounding box of object\n";
            return false;
        }
    }
    return true;
}

int mapper::point_to_voxel(const Eigen::RowVector3d &point)
{
    if(!is_point_within_range(point)) return -1;
    int sum1 = std::floor((std::abs((point(0))-vg_min(0)))/csl)+1;
    int sum2 = sides(0)*(std::floor((std::abs(point(1)-vg_min(1)))/csl));
    int sum3 = sides(0)*sides(1)*(std::floor((std::abs(point(2)-vg_min(2)))/csl));
    int voxel = sum1 + sum2 + sum3;
    return voxel;
}

// this method takes a triangle and returns the voxels that the triangle is present in
std::vector<int> mapper::triangle_to_voxels(const Eigen::Matrix<double,3,3> &triangle)
{
    std::vector<int> cubes;
    // min and max coordinates of the triangle
    Eigen::RowVector3d min = triangle.colwise().minCoeff();
    Eigen::RowVector3d max = triangle.colwise().maxCoeff();
    Eigen::RowVector3d max_xy_min_z = {max(0), max(1), min(2)};
    int min_voxel = point_to_voxel(min);
    int max_voxel = point_to_voxel(max);
    int max_xy_min_z_voxel = point_to_voxel(max_xy_min_z);
    if(min_voxel == -1 || max_voxel == -1 || max_xy_min_z_voxel == -1)
    {
        std::cout << "triangle out of range\n";
        std::vector<int> empty;
        empty.push_back(-1);
        return empty;
    }
    if(min_voxel != max_voxel){
        int mult = std::floor((max_xy_min_z_voxel - min_voxel)/(sides(0)));
        int add = (max_xy_min_z_voxel - min_voxel)%sides(0);
        int z_mult = std::floor((max_voxel - min_voxel)/(sides(0)*sides(1)));
        for(int i = 0; i <= z_mult; ++i){
            for(int j = 0; j <= mult; ++j)
            {
                for(int k = 0; k <= add; ++k)
                {
                    cubes.push_back(min_voxel + j*sides(0) + k + i*sides(0)*sides(1));
                }
            }
        }
    }
    else
    {
        cubes.push_back(min_voxel);
    }
    return cubes;
}

// each vector represents a voxel and the triangles in that voxel
std::vector<std::vector<int>> mapper::voxel_to_triangle_mapping()
{
    std::vector<std::vector<int>> mapping;
    for(int i = 0; i < sides.prod(); ++i)
    {
        std::vector<int> v;
        mapping.push_back(v);
    }

    for(int i = 0; i < F.rows(); ++i)
    {
        Eigen::Matrix<double,3,3> triangle;
        triangle << V.row(F.row(i)(0)), V.row(F.row(i)(1)), V.row(F.row(i)(2));
        std::vector<int> voxels = triangle_to_voxels(triangle);
        for(int voxel : voxels)
        {
            mapping[voxel].push_back(i);
        }
    }
    return mapping;
}

// find the triangle that 'point' lies on, given the triangles in a voxel
int mapper::point_to_triangle(const Eigen::RowVector3d &point)
{
    Eigen::MatrixXd point_d;
    point_d = point;
    int voxel = point_to_voxel(point);
    if(voxel == -1) return -1;
    std::vector<int> triangles = map[voxel];
    for(int i = 0; i < triangles.size(); ++i)
    {
        int F_index = triangles[i];
        Eigen::MatrixXd bc;
        igl::barycentric_coordinates(point, V.row(F.row(F_index)(0)), V.row(F.row(F_index)(1)), V.row(F.row(F_index)(2)), bc);
        if(bc(0) >=0 && bc(1) >= 0 && bc(2) >= 0)
        {
            Eigen::VectorXd sqrD;
            Eigen::VectorXi I;
            Eigen::MatrixXd C;
            Eigen::Matrix<int,1,3> Fd;
            Fd = F.row(F_index);
            igl::point_mesh_squared_distance(point_d, V, Fd, sqrD, I, C);
            if(sqrD.sum() < 1e-7)
            {
                return F_index;
            }
        }
    }
    std::cout << point_to_string(point) << " does not lie on a face of the mesh\n";
    return -1;
}

double mapper::randMToN(double M, double N)
{
    srand (static_cast <unsigned> (time(0)));
    return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}

std::string mapper::point_to_string(Eigen::Vector3d point)
{
    return std::to_string(point(0)) + ", " + std::to_string(point(1)) + ", " + std::to_string(point(2));
}

// choose a random triangle, and then a random point on that triangle
// write into test_triangle the value returned by point_to_triangle
// return true if test_triangle matches the randomly chosen triangle
// test_triangle and test_point can also be used to visually check the result
bool mapper::test_point_to_triangle(int &test_triangle, Eigen::RowVector3d &test_point)
{
    double u = randMToN(0.1, 0.8);
    double v = randMToN(0.1, 0.9-u);
    double w = 1-u-v;
    int F_index = rand()%(F.rows()+1);
    test_point << u*V.row(F.row(F_index)(0)) + v*V.row(F.row(F_index)(1)) + w*V.row(F.row(F_index)(2));
    test_triangle = F_index;
    int triangle = point_to_triangle(test_point);
    return triangle == F_index;
}

// set every face colour to white except the test_triangle, which is green
void mapper::set_vertex_colours(Eigen::MatrixXd &colours, const int &test_triangle)
{
    colours = Eigen::MatrixXd(colours.rows(),3);
    for(int i = 0; i < colours.rows(); ++i)
    {
        colours(i,0) = 1;
        colours(i,1) = 1;
        colours(i,2) = 1;
    }
    for(int i = 0; i < 3; ++i)
    {
        colours(test_triangle, 0) = 0;
        colours(test_triangle, 1) = 1;
        colours(test_triangle, 2) = 0;
    }
}

// this displays the mesh with every triangle coloured white, except test_triangle
// test_point is displayed as a red circle
void mapper::display_point_and_triangle(int test_triangle, Eigen::RowVector3d test_point)
{
    Eigen::MatrixXd colours = Eigen::MatrixXd(F.rows(),3);
    set_vertex_colours(colours, test_triangle);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(colours);

//    Eigen::RowVector3d v1 = V.row(F(test_triangle,1)) - V.row(F(test_triangle,0));
//    Eigen::RowVector3d v2 = V.row(F(test_triangle,2)) - V.row(F(test_triangle,0));
//    double triangle_area = v1.cross(v2).norm();
//    viewer.data().point_size = triangle_area*3;

    viewer.data().point_size = 10;
    viewer.data().add_points(test_point, Eigen::RowVector3d(1,0,0));
    viewer.data().set_face_based(true);
    viewer.launch();
}

void mapper::print_voxel_info(std::vector<std::vector<int>> map)
{
    int sum = 0;
    int sum_no0 = 0;
    int count_no0 = 0;
    for(std::vector<int> v : map){
        sum += v.size();
        if(v.size() != 0)
        {
            sum_no0 += v.size();
            count_no0 += 1;
        }
    }
    std::cout << "average no. of triangles per voxel: " << sum/map.size() <<'\n';
    std::cout << "voxels with 0 triangles: " << map.size() - count_no0 << '\n';
    std::cout << "average no. of triangles per voxel excluding zeros: " << sum_no0/count_no0 <<'\n';
}

