//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_GEO_UTIL_H
#define TLC_GEO_UTIL_H

#include <Eigen/Core>
#include <vector>

#include "corecrt_math_defines.h"

typedef Eigen::Vector2d Point;

// return the vector after rotating 90 degree counter clockwise
Eigen::Vector2d rotate_90deg(const Eigen::Vector2d &vec);

// compute the angle mod 2pi
// return: angle in [0, 2pi)
double angle_mod_2PI(double a);

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
double compute_vector_angle(const Point &p);

// compute the counter clockwise rotation angle from a1 to a2
// return: angle in [0, 2pi)
double compute_rotation_angle(double a1, double a2);

// for a ray sweeping counter-clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_ccw(double a, double a1, double a2);

// for a ray sweeping clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_cw(double a, double a1, double a2);

// test if angle a is between a1(included) and a2(not included)
bool is_angle_between(double a, double a1, double a2);

// compute total signed area of a polygon
double compute_total_signed_area(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges);

// compute total signed area of a polygon and its gradient wrt. polygon vertices
double compute_total_signed_area_with_gradient(const std::vector<Point> &vertices,
                                               const std::vector<std::pair<size_t,size_t>> &edges,
                                               Eigen::Matrix2Xd &dArea_dv);

// compute triangle area using (robust) Heron's formula
// input: squared edge lengths of the triangle
double compute_Heron_tri_area(double d1, double d2, double d3);

// computed total signed area of a triangle mesh
double compute_total_signed_mesh_area(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F);

// compute minimum signed area of a triangle mesh
double compute_min_signed_mesh_area(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F);

// compute signed area of all triangles in a mesh
void compute_signed_tri_areas(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
                                Eigen::VectorXd &areaList);

// compute signed area of a triangle
// input: 2D coordinates of 3 points of triangle
// return: signed area of the triangle
double compute_tri_signed_area(const Point &p1, const Point &p2, const Point &p3);

// compute total unsigned area of a triangle mesh
double compute_total_unsigned_area(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

//
void compute_squared_edge_Length(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
                                 Eigen::Matrix3Xd &D);

//
void extract_mesh_boundary_edges(const Eigen::Matrix3Xi &faces,
                                 std::vector<std::pair<size_t,size_t>> &boundary_edges);

// compute angle between from vector1 to vector2
// return: angle in (-pi, pi]
double compute_vec_vec_angle(const Eigen::Vector2d &vec1, const Eigen::Vector2d &vec2);



// compute indices of over-winded interior vertices
// input: 2D triangle mesh (V, F); is_boundary_vertex
// output: indices of winded interior vertices
void compute_winded_interior_vertices(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
                                      const std::vector<bool> &is_boundary_vertex,
                                      std::vector<size_t> &winded_vertices);

#endif //TLC_GEO_UTIL_H

