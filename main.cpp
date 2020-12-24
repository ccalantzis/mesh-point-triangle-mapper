#include "mapper.h"

void point_provided(mapper &m)
{
    Eigen::RowVector3d point = {0,0,0};
    int triangle = m.point_to_triangle(point);
    if(triangle == -1) return;
    m.display_point_and_triangle(triangle, point);
}

void random_point_generated(mapper &m)
{
    int triangle;
    Eigen::RowVector3d point;
    bool passed = m.test_point_to_triangle(triangle, point);
    if(triangle == -1) return;
    m.display_point_and_triangle(triangle, point);
}

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        std::cout << "required command line argument: path to mesh file\n";
        return -1;
    }

    std::string mesh_filename = argv[1];    
    mapper m = mapper(mesh_filename);

    random_point_generated(m);
}
