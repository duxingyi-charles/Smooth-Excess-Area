//
// Created by Charles Du on 4/26/21.
//

#ifndef TLC_ARRANGEMENT_H
#define TLC_ARRANGEMENT_H

#include "Circular_Arc.h"

struct Arc_Edge {
    // first vertex index
    size_t id1;
    // second vertex index
    size_t id2;
    // circular arc geometry
    Circular_Arc arc;
};

struct SubArc_Edge {
    // first vertex index
    size_t id1;
    // second vertex index
    size_t id2;
    // parent edge index
    size_t parent_id;
    // circular arc geometry
    Circular_Arc arc;
};

class Arrangement {
public:
    Arrangement() = default;
    ~Arrangement() = default;

    // compute self-arrangement of a circular arc loop
    // return: arc edges of cell arrangements, winding numbers of cell arrangements
    static void compute_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                    // output
                                    std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
                                    std::vector<bool> &is_intersection_point,
                                    std::vector<int>  &arc1_of_intersection,
                                    std::vector<int>  &arc2_of_intersection,
                                    std::vector<std::vector<SubArc_Edge>> &edges_of_cell,
                                    std::vector<int> &windings);

    // compute self-arrangement of a set of arc loops
    // return: arc edges of cell arrangements, winding number of cell arrangements
    static void compute_multi_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                    // output
                                    std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
                                    std::vector<bool> &is_intersection_point,
                                    std::vector<int>  &arc1_of_intersection,
                                    std::vector<int>  &arc2_of_intersection,
                                    std::vector<std::vector<SubArc_Edge>> &edges_of_cell,
                                    std::vector<int> &windings);

    // subdivide input arcs by their intersections
    static void subdivide_polyArc_by_intersection(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                             // output
                                             std::vector<Point> &pts,
                                             std::vector<SubArc_Edge> &pEdges,
                                             std::vector<bool> &is_intersection_point,
                                             std::vector<int>  &arc1_of_intersection,
                                             std::vector<int>  &arc2_of_intersection);

    // decompose 2D plane into arrangement cells
    // input: subdivided polyArc {pts, pEdges}
    // output: cells, each cell represented by a list of boundary half-edges
    // notes:
    // * index convention for half-edge:
    // ** i : index of original edge e
    // ** 2i+1: same direction as e, on the left side of e
    // ** 2i  : opposite direction of e, on the right side of e
    static void decompose_into_cells(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &pEdges,
                                     // output
                                     std::vector<std::vector<size_t>> &eIn,
                                     std::vector<std::vector<size_t>> &eOut,
                                     std::vector<std::vector<size_t>> &cells);

    // compute winding numbers for arrangement cells
    static void compute_cell_windings(const std::vector<SubArc_Edge> &pEdges,
                                      const std::vector<std::vector<size_t>> &eIn,
                                      const std::vector<std::vector<size_t>> &eOut,
                                      const std::vector<std::vector<size_t>> &cells,
                                      // output
                                      std::vector<int> &windings);

    // compute arrangement cells and their winding numbers from chains
    static void compute_cells_and_windings(const std::vector<SubArc_Edge> &pEdges,
                                           const std::vector<std::vector<size_t>> &eIn,
                                           const std::vector<std::vector<size_t>> &eOut,
                                           const std::vector<std::vector<size_t>> &chains,
                                           // output
                                           std::vector<std::vector<size_t>> &cells,
                                           std::vector<int> &windings);

private:
    // trace half-edges into chains
    // input: next_hEdge[i] is the index of the half-edge after half-edge i.
    // output: chains
    // - each chain is a connected component of half-edges
    static void trace_chains(const std::vector<size_t> &next_hEdge,
                      // output
                      std::vector<std::vector<size_t>> &chains);

    // find the unbounded cell
    // return: index of the unbounded cell
    static size_t find_the_unbounded_cell(const std::vector<SubArc_Edge> &pEdges,
                                          const std::vector<std::vector<size_t>> &eIn,
                                          const std::vector<std::vector<size_t>> &eOut,
                                          const std::vector<std::vector<size_t>> &cells,
                                          const std::vector<size_t> &cell_of_hE);

    // find the unbounded cell given the left-most edge
    // return: index of the unbounded cell
    static size_t find_the_unbounded_cell(const std::vector<SubArc_Edge> &pEdges,
                                                const std::vector<std::vector<size_t>> &eIn,
                                                const std::vector<std::vector<size_t>> &eOut,
                                                const std::vector<std::vector<size_t>> &cells,
                                                const std::vector<size_t> &cell_of_hE,
                                                Point_Arc_Location most_left_location,
                                                size_t most_left_edge_id);
};


#endif //TLC_ARRANGEMENT_H
