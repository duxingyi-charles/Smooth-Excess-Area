//
// Created by Charles Du on 4/29/21.
//

#ifndef TLC_ARC_OVERLAP_FORMULATION_H
#define TLC_ARC_OVERLAP_FORMULATION_H

#include "Total_Lifted_Content.h"
#include "Arc_Occupancy.h"

#include "timing.h"


using namespace Eigen;

class Arc_Overlap_Formulation {
public:
    Arc_Overlap_Formulation(const MatrixXd &rest_vertices, Matrix2Xd init_vertices,
                            Matrix3Xi faces,
                            const VectorXi &handles, const std::string &form,
                            double alphaRatio, double alpha, double theta);

    ~Arc_Overlap_Formulation() = default;


    double compute_energy(const Eigen::VectorXd &x);

    // compute arc overlap energy,
    // record the decomposition of the energy in energy_list
    // decomposition: lifted content per triangle, arc segment area per edge, arc occupancy area
    double compute_energy(const Eigen::VectorXd &x, Eigen::VectorXd &energy_list);

    double compute_energy_with_gradient(const Eigen::VectorXd &x,
                                        Eigen::VectorXd &grad);

    // approximate Hessian using the Hessian of (TLC - total signed area)
    double compute_energy_with_gradient_approxProjectedHessian(const Eigen::VectorXd &x,
                                                               Eigen::VectorXd &grad, SpMat &Hess);

    // approximate Hessian using the Hessian of (TLC - total signed area)
    double compute_energy_with_gradient_approxProjectedHessian(const Eigen::VectorXd &x,
                                                               Eigen::VectorXd &energy_list,
                                                               Eigen::VectorXd &grad, SpMat &Hess);


    static double computeAlpha(const MatrixXd &restV,
                               const Matrix2Xd &initV,
                               const Matrix3Xi &Faces,
                               const std::string &form, double alphaRatio);

    VectorXd get_x0() const { return x0; }

    std::vector<std::pair<size_t,size_t>> get_boundary_edges() const { return boundary_edges; }

    VectorXi get_freeI() const { return freeI; }

    const Matrix2Xd& get_V() const { return V; }

    // whether arc-arc intersection is found in the most recent call of
    // compute_energy() and compute_energy_with_gradient()
    bool has_found_intersection;

private:
    // x = Flatten(freeV)
    void update_V(const Eigen::VectorXd &x);

private:
    // V indices of triangles
    Matrix3Xi F;
    // indices of free vertices
    VectorXi freeI;
    // map: V index --> freeV index. If V(i) is not free, indexDict(i) = -1.
    VectorXi indexDict;
    // freeV indices of triangles.
    Matrix3Xi F_free;

    // boundary edges
    std::vector<std::pair<size_t,size_t>> boundary_edges;

    // TLC energy
    Total_Lifted_Content tlc;
    // arc occupancy
    Arc_Occupancy arcOccupancy;


    // current V of target mesh
    Matrix2Xd V;

    // initial variable vector
    VectorXd x0;


};


#endif //TLC_ARC_OVERLAP_FORMULATION_H
