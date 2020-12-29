# mesh-point-triangle-mapper
Given a 3D mesh of triangular faces and a 3D point, this program will return the face of the mesh that the point lies on, if any.

## Dependencies

- [libigl](https://github.com/libigl/libigl.git)

## Overview

The method `point_to_triangle` returns the row ID of the face that the given point lies on, or -1 if the point does not lie on any face.

`display_point_and_triangle` takes a point and face as input and opens a window displaying the mesh, with the point as a red circle, and the face coloured green. All other faces are coloured white. This method can be used to confirm the output of `point_to_triangle`: the red circle should lie on the green face.

main.cpp contains examples on how to use this program. `point_provided` allows the user to provide a point, which may or may not lie on one of the faces of the mesh. If the point does lie on a face of the mesh, `display_point_and_triangle` is called to allow the user to confirm the result visually.

`random_point_generated` will choose a random triangle on the mesh, generate a point on this triangle, and then pass the point to `test_point_to_triangle`. `test_point_to_triangle` checks that the output of `point_to_triangle` matches the randomly generated triangle, which the point lies on. `display_point_to_triangle` is then called to confirm the result visually.

## Screenshots

This is an example using the head.ply mesh, which can be found in the build folder. The first image is zoomed out, showing the point in red. After zooming in we can see the green triangle, confirming that the output of `point_to_triangle` is correct.

![head.ply zoomed out](https://github.com/ccalantzis/mesh-point-triangle-mapper/blob/main/demo_images/head-ply-zoomed-out.png)

![head.ply zoomed in](https://github.com/ccalantzis/mesh-point-triangle-mapper/blob/main/demo_images/head-ply-zoomed-in.png)

## Technical details

Testing if a point lies on a triangle involves calculating the barycentric coordinates of the point with respect to the triangle. If all 3 of these coordinates are positive, the point lies within the triangle, if any are negative, the point is outside the triangle. However, positive barycentric coordinates do not guarantee the point will lie on the triangle, we also need to check the distance between the point and the triangle. If it is less than some threshold (I use 1e-7), we say the point lies on the triangle. Rather than repeating this process for each triangle in the mesh, a simpler solution is to split the mesh into a number of cubic sections, and map each section to the triangles it contains. First we check which of these sections contains the point, and then we only have to test the triangles contained in that section. 

The axis aligned bounding box of a mesh is the minimum rectangular prism which contains the entire mesh, with the constraint that the sides of the box are aligned with the x, y, and z axes. The bounding box can then be split into cubes, and is called a voxel grid. Each instance of the mapper class stores a matrix containing the minimum and maximum corners of each cube in the voxel grid. This matrix can be used to find which voxel cube the point is present in, if any. Each instance of the mapper class also stores a mapping of voxel cubes to triangles. This map is used to find a subset of triangles to test with the point, using the method described above.



