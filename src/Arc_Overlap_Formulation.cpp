//
// Created by Charles Du on 4/29/21.
//

#include "Arc_Overlap_Formulation.h"
#include <iostream>
#include <utility>

Arc_Overlap_Formulation::Arc_Overlap_Formulation(const MatrixXd &rest_vertices, Matrix2Xd init_vertices,
                                                 Matrix3Xi faces, const VectorXi &handles,
                                                 const std::string &form, double alphaRatio,
                                                 double alpha, double theta) :
                                                 F(std::move(faces)), V(std::move(init_vertices)),
                                                 arcOccupancy(theta),
                                                 // we assume there are intersection at the beginning
                                                 has_found_intersection(true)
                                                 {
    // compute freeI: indices of free vertices
    int nV = rest_vertices.cols();
    int vDim = 2;

    std::vector<bool> freeQ(nV, true);
    for (auto i = 0; i < handles.size(); ++i) {
        freeQ[handles(i)] = false;
    }
    freeI.resize(nV - handles.size());
    int ii = 0;
    for (int i = 0; i < nV; ++i) {
        if (freeQ[i]) {
            freeI[ii] = i;
            ++ii;
        }
    }
    std::sort(freeI.data(), freeI.data() + freeI.size());


    // compute indexDict and F_free
    indexDict = VectorXi::Constant(nV, -1);
    for (auto i = 0; i < freeI.size(); ++i) {
        indexDict(freeI(i)) = i;
    }

    F_free.resize(F.rows(), F.cols());
    for (auto i = 0; i < F.cols(); ++i) {
        for (auto j = 0; j < F.rows(); ++j) {
            F_free(j, i) = indexDict(F(j, i));
        }
    }

    // compute alpha if it's not explicitly given
    double param_alpha;
    if (alpha >= 0) param_alpha = alpha;
    else {  //alpha < 0. Deduce alpha from alphaRatio
        param_alpha = computeAlpha(rest_vertices, V, F, form, alphaRatio);
    }
    std::cout << "form: " << form << std::endl;
    std::cout << "alphaRatio: " << alphaRatio << std::endl;
    std::cout << "alpha: " << param_alpha << std::endl;
    std::cout << "theta: " << theta << std::endl;

    // compute x0 from initV
    x0.resize(vDim * freeI.size());
    for (auto i = 0; i < freeI.size(); ++i) {
        int vi = freeI(i);
        for (int j = 0; j < vDim; ++j) {
            x0(i * vDim + j) = V(j, vi);
        }
    }

    // initialize TLC
    tlc.initialize(rest_vertices, F, form, param_alpha);

    // extract boundary edges
    extract_mesh_boundary_edges(F, boundary_edges);

}

double Arc_Overlap_Formulation::computeAlpha(const MatrixXd &restV, const Matrix2Xd &initV, const Matrix3Xi &Faces,
                                             const std::string &form, double alphaRatio) {
    unsigned nF = Faces.cols();
    double rest_measure;
    double init_measure;
    // tri
    init_measure = compute_total_signed_mesh_area(initV, Faces);
    if (form == "harmonic") {
        rest_measure = compute_total_unsigned_area(restV, Faces);
    } else { // Tutte form
        rest_measure = nF * sqrt(3) / 4;
    }

    // alpha
    return alphaRatio * fabs(init_measure) / rest_measure;
}

double Arc_Overlap_Formulation::compute_energy(const VectorXd &x) {
    update_V(x);

    // TLC
    double tlc_energy = tlc.compute_total_lifted_content(V);

    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    double arc_segment_area = arcOccupancy.compute_arc_curve_segment_area(vertices, boundary_edges);
    double arc_occupancy = arcOccupancy.compute_arc_occupancy(vertices, boundary_edges);
    has_found_intersection = arcOccupancy.has_found_intersection;

    return tlc_energy + arc_segment_area - arc_occupancy;
}

void Arc_Overlap_Formulation::update_V(const VectorXd &x) {
    int vDim = 2;
    for (auto i = 0; i < freeI.size(); ++i) {
        for (auto j = 0; j < vDim; ++j) {
            V(j, freeI(i)) = x[i * vDim + j];
        }
    }
}

double Arc_Overlap_Formulation::compute_energy_with_gradient(const VectorXd& x, VectorXd& grad) {
    update_V(x);

    // TLC
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Matrix2Xd tlc_grad;
    double tlc_energy = tlc.compute_total_lifted_content_with_gradient(V, tlc_grad);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    global_TLC_time += std::chrono::duration_cast<Time_duration>(t2 - t1);
    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }

    // arc segment area
    t1 = std::chrono::steady_clock::now();
    Matrix2Xd arc_seg_grad;
    double arc_segment_area = arcOccupancy.compute_arc_curve_segment_area_with_gradient(vertices, boundary_edges,
        arc_seg_grad);
    t2 = std::chrono::steady_clock::now();
    global_arc_seg_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // arc occupancy
    t1 = std::chrono::steady_clock::now();
    Matrix2Xd arc_occupancy_grad;
    double arc_occupancy = arcOccupancy.compute_arc_occupancy_with_gradient(vertices, boundary_edges,
        arc_occupancy_grad);
    t2 = std::chrono::steady_clock::now();
    global_arc_occupancy_time += std::chrono::duration_cast<Time_duration>(t2 - t1);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // full gradient
    Matrix2Xd m_grad = tlc_grad + arc_seg_grad - arc_occupancy_grad;

    // gradient of free vertices
    grad.resize(freeI.size() * m_grad.rows());
    for (int i = 0; i < freeI.size(); ++i) {
        grad(2 * i) = m_grad(0, freeI(i));
        grad(2 * i + 1) = m_grad(1, freeI(i));
    }

    // energy
    return tlc_energy + arc_segment_area - arc_occupancy;
}


double Arc_Overlap_Formulation::compute_energy_with_gradient_approxProjectedHessian(const VectorXd &x, VectorXd &grad,
                                                                                    SpMat &Hess) {
    update_V(x);

    // TLC energy, TLC full gradient,
    // and PSD projected Hessian of (TLC - total signed area) on free vertices
    Matrix2Xd tlc_grad;
    double tlc_energy = tlc.compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(V, freeI, F_free,
                                                                                                 tlc_grad, Hess);

    // arc segment area - arc occupancy
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    Matrix2Xd arc_seg_grad;
    double arc_segment_area =arcOccupancy.compute_arc_curve_segment_area_with_gradient(vertices, boundary_edges,
                                                                                       arc_seg_grad);
    Matrix2Xd arc_occupancy_grad;
    double arc_occupancy = arcOccupancy.compute_arc_occupancy_with_gradient(vertices, boundary_edges,
                                                                            arc_occupancy_grad);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // full gradient
    Matrix2Xd m_grad = tlc_grad + arc_seg_grad - arc_occupancy_grad;

    // gradient of free vertices
    grad.resize(freeI.size() * m_grad.rows());
    for (int i = 0; i < freeI.size(); ++i) {
        grad(2*i) = m_grad(0,freeI(i));
        grad(2*i+1) = m_grad(1,freeI(i));
    }

    // energy
    return tlc_energy + arc_segment_area - arc_occupancy;
}

double Arc_Overlap_Formulation::compute_energy(const VectorXd &x, VectorXd &energy_list) {
    update_V(x);

    // TLC
    VectorXd lifted_content_list;
    double tlc_energy = tlc.compute_total_lifted_content(V,lifted_content_list);

    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    VectorXd arc_curve_segment_area_list;
    double arc_segment_area = arcOccupancy.compute_arc_curve_segment_area(vertices, boundary_edges,
                                                                          arc_curve_segment_area_list);
    double arc_occupancy = arcOccupancy.compute_arc_occupancy(vertices, boundary_edges);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // fill the energy decomposition into energy_list
    energy_list.resize(lifted_content_list.size()+arc_curve_segment_area_list.size()+1);
    int ii = 0;
    for (int i = 0; i < lifted_content_list.size(); ++i) {
        energy_list(ii) = lifted_content_list(i);
        ++ii;
    }
    for (int i = 0; i < arc_curve_segment_area_list.size(); ++i) {
        energy_list(ii) = arc_curve_segment_area_list(i);
        ++ii;
    }
    energy_list(ii) = -arc_occupancy;


    return tlc_energy + arc_segment_area - arc_occupancy;
}

double
Arc_Overlap_Formulation::compute_energy_with_gradient_approxProjectedHessian(const VectorXd &x, VectorXd &energy_list,
                                                                             VectorXd &grad, SpMat &Hess) {
    update_V(x);

    // TLC energy, TLC full gradient,
    // and PSD projected Hessian of (TLC - total signed area) on free vertices
    VectorXd lifted_content_list;
    Matrix2Xd tlc_grad;
    double tlc_energy = tlc.compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(V, freeI, F_free,
                                                                                                 lifted_content_list,
                                                                                                 tlc_grad, Hess);

    // arc segment area - arc occupancy
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    VectorXd  arc_curve_segment_area_list;
    Matrix2Xd arc_seg_grad;
    double arc_segment_area =arcOccupancy.compute_arc_curve_segment_area_with_gradient(vertices, boundary_edges,
                                                                                       arc_curve_segment_area_list,
                                                                                       arc_seg_grad);
    Matrix2Xd arc_occupancy_grad;
    double arc_occupancy = arcOccupancy.compute_arc_occupancy_with_gradient(vertices, boundary_edges,
                                                                            arc_occupancy_grad);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // full gradient
    Matrix2Xd m_grad = tlc_grad + arc_seg_grad - arc_occupancy_grad;

    // gradient of free vertices
    grad.resize(freeI.size() * m_grad.rows());
    for (int i = 0; i < freeI.size(); ++i) {
        grad(2*i) = m_grad(0,freeI(i));
        grad(2*i+1) = m_grad(1,freeI(i));
    }

    // fill the energy decomposition into energy_list
    energy_list.resize(lifted_content_list.size()+arc_curve_segment_area_list.size()+1);
    int ii = 0;
    for (int i = 0; i < lifted_content_list.size(); ++i) {
        energy_list(ii) = lifted_content_list(i);
        ++ii;
    }
    for (int i = 0; i < arc_curve_segment_area_list.size(); ++i) {
        energy_list(ii) = arc_curve_segment_area_list(i);
        ++ii;
    }
    energy_list(ii) = -arc_occupancy;

    // energy
    return tlc_energy + arc_segment_area - arc_occupancy;
}


