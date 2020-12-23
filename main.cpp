#include "mapper.h"

int main(int argc, char *argv[])
{
    srand (static_cast <unsigned> (time(0)));
    std::cout << std::boolalpha;;
//    mapper m = mapper("inspired_mesh.obj");
    mapper m = mapper("head.ply");

    int test_triangle = 0;
    Eigen::RowVector3d test_point;
    bool passed = m.test_point_to_triangle(test_triangle, test_point);
    std::cout << "test passed: " << passed << "\n\n";

    m.display_point_and_triangle(test_triangle, test_point);
}
