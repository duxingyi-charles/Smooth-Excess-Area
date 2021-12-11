//
// Created by Charles Du on 4/28/21.
//

#include "Arc_Occupancy.h"
#include "Arrangement.h"
#include "timing.h"

double Arc_Occupancy::compute_arc_occupancy(const std::vector<Point> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges) {

    return compute_arc_occupancy(vertices, edges, param_theta);
}

double Arc_Occupancy::compute_arc_loop_area(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges)
{
    std::vector<std::pair<size_t, size_t>> raw_edges;
    raw_edges.reserve(edges.size());
    for (const auto & e : edges) {
        raw_edges.emplace_back(e.id1, e.id2);
    }

    // signed area of the polygon
    double polygon_area = compute_total_signed_area(pts, raw_edges);

    // area of arc segments
    double arc_seg_area = 0;
    for (const auto & e: edges) {
        arc_seg_area += e.arc.get_segment_area();
    }

    //
    return polygon_area + arc_seg_area;
}

double Arc_Occupancy::compute_arc_occupancy(const std::vector<Point> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges, double theta) {
    // create arc edges
    std::vector<Arc_Edge> arc_edges;
    arc_edges.reserve(edges.size());

    for (const auto & e : edges) {
        arc_edges.emplace_back(Arc_Edge{e.first, e.second,
                                        Circular_Arc(vertices[e.first], vertices[e.second], theta)});
    }

    // compute arrangement and winding numbers
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    std::vector<int>  arc1_of_intersection;
    std::vector<int>  arc2_of_intersection;
    std::vector<std::vector<SubArc_Edge>> edges_of_cell;
    std::vector<int> windings;
//    Arrangement::compute_arrangement(vertices,arc_edges,
//                                     pts,pEdges,is_intersection_point,
//                                     arc1_of_intersection,arc2_of_intersection,
//                                     edges_of_cell,windings);
    Arrangement::compute_multi_arrangement(vertices,arc_edges,
                                     pts,pEdges,is_intersection_point,
                                     arc1_of_intersection,arc2_of_intersection,
                                     edges_of_cell,windings);
    has_found_intersection = (pts.size() > vertices.size());

    // compute occupancy
    double occupancy = 0;
    for (int i = 0; i < windings.size(); ++i) {
        if (windings[i] > 0) {
            occupancy += compute_arc_loop_area(pts, edges_of_cell[i]);
        }
    }

    //
    return occupancy;
}

double Arc_Occupancy::compute_arc_occupancy_with_gradient(const std::vector<Point> &vertices,
                                                          const std::vector<std::pair<size_t, size_t>> &edges,
                                                          double theta, Eigen::Matrix2Xd &grad) {

    // create arc edges
    std::vector<Arc_Edge> arc_edges;
    arc_edges.reserve(edges.size());

    for (const auto & e : edges) {
        arc_edges.emplace_back(Arc_Edge{e.first, e.second,
                                        Circular_Arc(vertices[e.first], vertices[e.second], theta)});
    }

    // compute arc's derivatives wrt. input vertices
    for (auto & e : arc_edges) {
        e.arc.update_derivatives();
    }

    // compute arrangement and winding numbers
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    std::vector<int>  arc1_of_intersection;
    std::vector<int>  arc2_of_intersection;
    std::vector<std::vector<SubArc_Edge>> edges_of_cell;
    std::vector<int> windings;
//    Arrangement::compute_arrangement(vertices,arc_edges,
//                                     pts,pEdges,is_intersection_point,
//                                     arc1_of_intersection, arc2_of_intersection,
//                                     edges_of_cell,windings);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Arrangement::compute_multi_arrangement(vertices,arc_edges,
                                     pts,pEdges,is_intersection_point,
                                     arc1_of_intersection, arc2_of_intersection,
                                     edges_of_cell,windings);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    global_arc_arrangement_time += std::chrono::duration_cast<Time_duration>(t2 - t1);
    has_found_intersection = (pts.size() > vertices.size());

    // compute occupancy and its derivatives wrt. pts and edge radii
    Eigen::Matrix2Xd dOccu_dp = Eigen::Matrix2Xd::Zero(2, pts.size());
    Eigen::VectorXd dOccu_dr2 = Eigen::VectorXd::Zero(edges.size());
    double occupancy = 0;
    for (int i = 0; i < windings.size(); ++i) {
        if (windings[i] > 0) {
            occupancy += compute_arc_loop_area_with_gradient(pts, edges_of_cell[i], dOccu_dp, dOccu_dr2);
        }
    }

    // apply chain rule to compute occupancy's derivatives wrt. input vertices
    grad = Eigen::Matrix2Xd::Zero(2, vertices.size());

    //// gradient from dOccu_dr2
    for (int ei = 0; ei < edges.size(); ++ei) {
        auto v1 = edges[ei].first;
        auto v2 = edges[ei].second;
        grad.col(v1) += dOccu_dr2(ei) * arc_edges[ei].arc.get_dr2_dP1();
        grad.col(v2) += dOccu_dr2(ei) * arc_edges[ei].arc.get_dr2_dP2();
    }

    //// gradient from dOccu_dp
    std::vector<bool> is_active_point(pts.size(), false);
    for (int i = 0; i < windings.size(); ++i) {
        if (windings[i] > 0) {
            for (const auto & e : edges_of_cell[i]) {
                is_active_point[e.id1] = true;
                is_active_point[e.id2] = true;
            }
        }
    }

    for (int i = 0; i < pts.size(); ++i) {
        if (!is_active_point[i]) continue;
        if (is_intersection_point[i]) {
            const auto &p = pts[i];
            auto arc1 = arc1_of_intersection[i];
            auto arc2 = arc2_of_intersection[i];
            auto o1 = arc_edges[arc1].arc.get_center();
            auto o2 = arc_edges[arc2].arc.get_center();
            Eigen::Matrix2d dp_do1, dp_do2;
            Eigen::Vector2d dp_dr1s, dp_dr2s;
            Circular_Arc::compute_intersection_gradient(p,o1,o2,dp_do1,dp_do2,dp_dr1s,dp_dr2s);
            // gradient from arc1 of intersection
            auto arc1_v1 = edges[arc1].first;
            auto arc1_v2 = edges[arc1].second;
            grad.col(arc1_v1) += (dOccu_dp.col(i).transpose() * dp_do1 * arc_edges[arc1].arc.get_dO_dP1()).transpose()
                            + (dOccu_dp.col(i).dot(dp_dr1s)) * arc_edges[arc1].arc.get_dr2_dP1();
            grad.col(arc1_v2) += (dOccu_dp.col(i).transpose() * dp_do1 * arc_edges[arc1].arc.get_dO_dP2()).transpose()
                            + (dOccu_dp.col(i).dot(dp_dr1s)) * arc_edges[arc1].arc.get_dr2_dP2();
            // gradient from arc2 of intersection
            auto arc2_v1 = edges[arc2].first;
            auto arc2_v2 = edges[arc2].second;
            grad.col(arc2_v1) += (dOccu_dp.col(i).transpose() * dp_do2 * arc_edges[arc2].arc.get_dO_dP1()).transpose()
                              + (dOccu_dp.col(i).dot(dp_dr2s)) * arc_edges[arc2].arc.get_dr2_dP1();
            grad.col(arc2_v2) += (dOccu_dp.col(i).transpose() * dp_do2 * arc_edges[arc2].arc.get_dO_dP2()).transpose()
                              + (dOccu_dp.col(i).dot(dp_dr2s)) * arc_edges[arc2].arc.get_dr2_dP2();
        } else {
            // pts[i] is not an intersection point, i.e. it equals vertices[i]
            grad.col(i) += dOccu_dp.col(i);
        }
    }


    //
    return occupancy;
}

double
Arc_Occupancy::compute_arc_loop_area_with_gradient(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges,
                                                   Eigen::Matrix2Xd &dArea_dp, Eigen::VectorXd &dArea_dr2) {
// note:
// we don't check the size of dArea_dp and dArea_dr2;
// we don't reset dArea_dp and dArea_dr2 to 0 at the beginning, we only accumulate on them.
// dArea_dr2[i] is the derivative wrt. squared radius of the i-th parent edge.

    std::vector<std::pair<size_t, size_t>> raw_edges;
    raw_edges.reserve(edges.size());
    for (const auto & e : edges) {
        raw_edges.emplace_back(e.id1, e.id2);
    }

    // signed area of the polygon and its gradient
    double polygon_area = compute_total_signed_area_with_gradient(pts, raw_edges, dArea_dp);

    // area of arc segments
    double arc_seg_area = 0;
    for (const auto & e: edges) {
        arc_seg_area += e.arc.get_segment_area();
    }

    // gradient associated with area of arc segments
    for (const auto & e : edges) {
        auto i = e.id1;
        auto j = e.id2;
        auto e_vec = pts[i] - pts[j];
        //auto e_len_squared = e_vec.squaredNorm();
        auto theta = e.arc.get_arc_angle();
        auto radius = e.arc.get_radius();
        Eigen::Vector2d da_dp = (tan(theta/2)/2) * e_vec;
        //
        dArea_dp.col(i) += da_dp;
        dArea_dp.col(j) -= da_dp;
        //dArea_dr2(e.parent_id) += (theta-sin(theta))/2 - tan(theta/2) * e_len_squared / (radius*radius*4);
        dArea_dr2(e.parent_id) += (theta - sin(theta)) / 2 - tan(theta / 2) * sin(theta/2) * sin(theta/2);
    }

    //
    return polygon_area + arc_seg_area;
}

double Arc_Occupancy::compute_arc_curve_segment_area(const std::vector<Point> &vertices,
                                                     const std::vector<std::pair<size_t, size_t>> &edges) const {
    double arc_seg_area = 0;

    double theta_factor = 0.25 * (param_theta - sin(param_theta)) / (1 - cos(param_theta));

    for (const auto & e: edges) {
        const auto & p1 = vertices[e.first];
        const auto & p2 = vertices[e.second];
        arc_seg_area += theta_factor * (p1-p2).squaredNorm();
    }

    return arc_seg_area;
}

double Arc_Occupancy::compute_arc_curve_segment_area_with_gradient(const std::vector<Point> &vertices,
                                                                   const std::vector<std::pair<size_t, size_t>> &edges,
                                                                   Eigen::Matrix2Xd &grad) const {
    double arc_seg_area = 0;
    grad = Eigen::Matrix2Xd::Zero(2, vertices.size());

    // compute arc curve segment area and gradient
    double theta_factor = 0.25 * (param_theta - sin(param_theta)) / (1 - cos(param_theta));
    double theta_factor2 = 2 * theta_factor;

    for (const auto & e: edges) {
        const auto & p1 = vertices[e.first];
        const auto & p2 = vertices[e.second];
        const auto & e_vec = p2 - p1;
        arc_seg_area += theta_factor * e_vec.squaredNorm();
        grad.col(e.second) += theta_factor2 * e_vec;
        grad.col(e.first)  -= theta_factor2 * e_vec;
    }

    return arc_seg_area;
}

double Arc_Occupancy::compute_arc_occupancy_with_gradient(const std::vector<Point> &vertices,
                                                          const std::vector<std::pair<size_t, size_t>> &edges,
                                                          Eigen::Matrix2Xd &grad) {
    return compute_arc_occupancy_with_gradient(vertices, edges, param_theta, grad);
}

double Arc_Occupancy::compute_arc_curve_segment_area(const std::vector<Point> &vertices,
                                                     const std::vector<std::pair<size_t, size_t>> &edges,
                                                     Eigen::VectorXd &segment_area_list) const {
    segment_area_list.resize(edges.size());

    double theta_factor = 0.25 * (param_theta - sin(param_theta)) / (1 - cos(param_theta));

    for (int i = 0; i < edges.size(); ++i) {
        const auto & p1 = vertices[edges[i].first];
        const auto & p2 = vertices[edges[i].second];
        segment_area_list(i) = theta_factor * (p1-p2).squaredNorm();
    }

    return segment_area_list.sum();
}

double Arc_Occupancy::compute_arc_curve_segment_area_with_gradient(const std::vector<Point> &vertices,
                                                                   const std::vector<std::pair<size_t, size_t>> &edges,
                                                                   Eigen::VectorXd &segment_area_list,
                                                                   Eigen::Matrix2Xd &grad) const {
//    double arc_seg_area = 0;
    segment_area_list.resize(edges.size());
    grad = Eigen::Matrix2Xd::Zero(2, vertices.size());

    // compute arc curve segment area and gradient
    double theta_factor = 0.25 * (param_theta - sin(param_theta)) / (1 - cos(param_theta));
    double theta_factor2 = 2 * theta_factor;

    for (int i = 0; i < edges.size(); ++i) {
        const auto & p1 = vertices[edges[i].first];
        const auto & p2 = vertices[edges[i].second];
        const auto & e_vec = p2 - p1;
        segment_area_list(i) = theta_factor * e_vec.squaredNorm();
        grad.col(edges[i].second) += theta_factor2 * e_vec;
        grad.col(edges[i].first)  -= theta_factor2 * e_vec;
    }

    return segment_area_list.sum();
}
