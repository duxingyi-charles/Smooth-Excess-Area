//
// Created by Charles Du on 4/26/21.
//

#include <cassert>
#include <map>
#include <queue>
#include <algorithm>
#include <limits>
#include "Arrangement.h"
#include "timing.h"

void Arrangement::compute_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                      std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
                                      std::vector<bool> &is_intersection_point,
                                      std::vector<int>  &arc1_of_intersection,
                                      std::vector<int>  &arc2_of_intersection,
                                      std::vector<std::vector<SubArc_Edge>> &edges_of_cell,
                                      std::vector<int> &windings)
{
    // step 1: subdivide input arcs by their intersections
    subdivide_polyArc_by_intersection(vertices, edges,pts, pEdges,is_intersection_point,
                                      arc1_of_intersection,arc2_of_intersection);

    // step 2: get all arrangement cells
    std::vector<std::vector<size_t>> eIn;
    std::vector<std::vector<size_t>> eOut;
    std::vector<std::vector<size_t>> cells;
    decompose_into_cells(pts, pEdges, eIn, eOut, cells);

    // convert each cell from a list of half-edge indices to a list arc edges
    edges_of_cell.clear();
    edges_of_cell.reserve(cells.size());
    for (const auto & hEdges : cells) {
        edges_of_cell.emplace_back();
        auto &es = edges_of_cell.back();
        for (auto hE : hEdges) {
            // create SubArc_Edge for half-edge
            // here, SubArc_Edge.parent_id record the index of the input edge in [std::vector<Arc_Edge> edges]
            if (hE%2 == 0) {
                const auto &pEdge = pEdges[hE/2];
                es.emplace_back(SubArc_Edge{pEdge.id2, pEdge.id1, pEdge.parent_id,
                                            Circular_Arc::reverse(pEdge.arc)});
            } else {
                es.emplace_back(pEdges[(hE-1)/2]);
            }
        }
    }

    // step 3: compute winding numbers for each cell
    compute_cell_windings(pEdges, eIn, eOut, cells, windings);

}


void Arrangement::subdivide_polyArc_by_intersection(
        const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
        std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
        std::vector<bool> &is_intersection_point,
        std::vector<int>  &arc1_of_intersection,
        std::vector<int>  &arc2_of_intersection)
{
    // init
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    pts = vertices;
    pEdges.clear();

    is_intersection_point.clear();
    is_intersection_point.resize(pts.size(), false);

    arc1_of_intersection.clear();
    arc1_of_intersection.resize(pts.size(), -1); // arc_of_intersection[i] = -1 means pts[i] is not an intersection
    arc2_of_intersection.clear();
    arc2_of_intersection.resize(pts.size(), -1);

    typedef std::pair<size_t, double> PID_Angle_Pair;
    std::vector<std::vector<PID_Angle_Pair>> edge_intersection_list(edges.size()); // record intersections on each input arc edge

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    global_arc_intersection_init_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // filter potentially intersecting arcs by bounding box
    t1 = std::chrono::steady_clock::now();
    std::vector<Rectangle> bbox_list;
    bbox_list.resize(edges.size());
    for (const auto& e : edges) {
        bbox_list.emplace_back(e.arc.get_bounding_box());
    }

    std::vector<bool> is_degenerate_edge(edges.size(),false);
    for (size_t i = 0; i < edges.size(); i++)
    {
        if (edges[i].arc.get_start_point() == edges[i].arc.get_end_point()) {
            is_degenerate_edge[i] = true;
        }
    }
    std::vector<std::pair<int,int>> potential_pair_list;
    for (int i = 0; i < edges.size(); ++i) {
        if (is_degenerate_edge[i]) {
            // skip degenerate arcs
            continue;
        }
        for (int j = i+1; j < edges.size(); ++j) {
            if (!is_degenerate_edge[j] && Rectangle::is_intersect(bbox_list[i], bbox_list[j])) {
                potential_pair_list.emplace_back(i, j);
            }
        }
    }
    t2 = std::chrono::steady_clock::now();
    global_arc_bbox_intersection_test_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // compute intersections
    t1 = std::chrono::steady_clock::now();
    // parallel version 1: separate parallel part from non-parallel part.
    std::vector<std::vector<Intersection_Point>> intersection_results(potential_pair_list.size());
#pragma omp parallel
    {
#pragma omp for nowait
        for (int pi = 0; pi < potential_pair_list.size(); pi++)
        {
            const auto& pair = potential_pair_list[pi];
            int i = pair.first;
            int j = pair.second;
            Circular_Arc::compute_intersection(edges[i].arc, edges[j].arc, intersection_results[pi]);
        }
    }
    

    for (size_t pi = 0; pi < potential_pair_list.size(); pi++)
    {
        const auto& pair = potential_pair_list[pi];
        int i = pair.first;
        int j = pair.second;
        const auto& result = intersection_results[pi];
        if (!result.empty()) {
            for (const auto& intersection : result) {
                pts.emplace_back(intersection.location);
                is_intersection_point.push_back(true);
                edge_intersection_list[i].emplace_back(pts.size() - 1, intersection.angle1);
                edge_intersection_list[j].emplace_back(pts.size() - 1, intersection.angle2);
                // intersectInfo
                arc1_of_intersection.push_back(i);
                arc2_of_intersection.push_back(j);
            }
        }
    }

    // parallel version 2: using omp critical
    // note: in this version, the order of points in pts is determined only at run-time. 
    //       In consequence, the result may differ from time to time.
//#pragma omp parallel
//    {
//#pragma omp for
//        for (int pi = 0; pi < potential_pair_list.size(); ++pi)
//        {
//            const auto& pair = potential_pair_list[pi];
//            int i = pair.first;
//            int j = pair.second;
//            std::vector<Intersection_Point> result;
//            Circular_Arc::compute_intersection(edges[i].arc, edges[j].arc, result);
//            if (!result.empty()) {
//#pragma omp critical
//                {
//                    for (const auto& intersection : result) {
//                        pts.emplace_back(intersection.location);
//                        is_intersection_point.push_back(true);
//                        edge_intersection_list[i].emplace_back(pts.size() - 1, intersection.angle1);
//                        edge_intersection_list[j].emplace_back(pts.size() - 1, intersection.angle2);
//                        // intersectInfo
//                        arc1_of_intersection.push_back(i);
//                        arc2_of_intersection.push_back(j);
//                    }
//                }
//            }
//            
//        }
//    }

    // non-parallel version:
    //for (const auto & pair : potential_pair_list) {
    //    int i = pair.first;
    //    int j = pair.second;
    //    std::vector<Intersection_Point> result;
    //    Circular_Arc::compute_intersection(edges[i].arc, edges[j].arc, result);
    //    if (!result.empty()) {
    //        for (const auto & intersection : result) {
    //            pts.emplace_back(intersection.location);
    //            is_intersection_point.push_back(true);
    //            edge_intersection_list[i].emplace_back(pts.size()-1, intersection.angle1);
    //            edge_intersection_list[j].emplace_back(pts.size()-1, intersection.angle2);
    //            // intersectInfo
    //            arc1_of_intersection.push_back(i);
    //            arc2_of_intersection.push_back(j);
    //        }
    //    }
    //}

    t2 = std::chrono::steady_clock::now();
    global_arc_arc_intersection_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // subdivide input arcs
    t1 = std::chrono::steady_clock::now();
    auto PID_angle_less_than = [](const PID_Angle_Pair &left, const PID_Angle_Pair &right)
    { return left.second < right.second; };
    auto PID_angle_greater_than = [](const PID_Angle_Pair &left, const PID_Angle_Pair &right)
    { return left.second > right.second; };

    for (size_t ei = 0; ei < edges.size(); ++ei) {
        if (edge_intersection_list[ei].empty()) {
            // no intersection point in ei, copy the input edge
            pEdges.emplace_back(SubArc_Edge{edges[ei].id1, edges[ei].id2, ei, edges[ei].arc});
        } else {
            double theta1 = edges[ei].arc.get_start_angle();
            double theta2 = edges[ei].arc.get_end_angle();
            double theta  = edges[ei].arc.get_arc_angle();

            std::vector<PID_Angle_Pair> sorted_intersections;
            if (theta > 0) { // ccw rotation from theta1 to theta2
                for (const auto &pid_angle : edge_intersection_list[ei]) {
                    sorted_intersections.emplace_back(pid_angle.first, compute_rotation_angle(theta1, pid_angle.second));
                }
                std::sort(sorted_intersections.begin(), sorted_intersections.end(), PID_angle_less_than);
            } else { // cw rotation from theta1 to theta2
                for (const auto &pid_angle : edge_intersection_list[ei]) {
                    sorted_intersections.emplace_back(pid_angle.first,
                                                      theta + compute_rotation_angle(theta2, pid_angle.second));
                }
                std::sort(sorted_intersections.begin(), sorted_intersections.end(), PID_angle_greater_than);
            }

            sorted_intersections.emplace_back(edges[ei].id2, theta);
            size_t last_vId = edges[ei].id1;
            double last_angle = 0;
            for (const auto & pid_angle : sorted_intersections) {
                pEdges.emplace_back(SubArc_Edge{last_vId, pid_angle.first, ei,
                                                Circular_Arc(pts[last_vId],pts[pid_angle.first],
                                                             pid_angle.second - last_angle,
                                                             edges[ei].arc.get_center(), edges[ei].arc.get_radius(),
                                                             theta1 + last_angle,
                                                             theta1 + pid_angle.second)});

                last_vId = pid_angle.first;
                last_angle = pid_angle.second;
            }
        }
    }

    t2 = std::chrono::steady_clock::now();
    global_arc_subdivision_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    //done
}

void Arrangement::decompose_into_cells(const std::vector<Point> &vertices, const std::vector<SubArc_Edge> &edges,
                                       std::vector<std::vector<size_t>> &eIn,
                                       std::vector<std::vector<size_t>> &eOut,
                                       std::vector<std::vector<size_t>> &cells)
{
    // find in-coming and out-going edges for each vertex
    eIn.clear();
    eIn.resize(vertices.size());
    eOut.clear();
    eOut.resize(vertices.size());


    for (int ei = 0; ei < edges.size(); ++ei) {
        eIn[edges[ei].id2].push_back(ei);
        eOut[edges[ei].id1].push_back(ei);
    }

    // find successor for each half-edge
    std::vector<size_t> next_hEdge(edges.size()*2, 0);

    struct Incident_Edge_Data {
        // edge index
        size_t id;
        // is the edge an in-coming edge for the vertex?
        bool is_in;
        // incident angle: angle of the vector tangent to the arc at the vertex
        double angle;
    };
    auto Incident_Edge_Data_greater_than = [](const Incident_Edge_Data &left, const Incident_Edge_Data &right)
    { return left.angle > right.angle; };

    for (int i = 0; i < vertices.size(); ++i) {
        const auto &in_edges = eIn[i];
        const auto &out_edges = eOut[i];
        if (in_edges.empty()) continue;

        if (in_edges.size() == 1) {
            // degree-2 vertex
            next_hEdge[2*in_edges[0]+1] = 2*out_edges[0] + 1;
            next_hEdge[2*out_edges[0]] = 2*in_edges[0];
        } else {
            // intersection vertex
            std::vector<Incident_Edge_Data> incident_edge_list;
            for (const auto id : in_edges) {
                double arc_angle = edges[id].arc.get_arc_angle();
                double angle2 = edges[id].arc.get_end_angle();
                double incident_angle = angle2 + ((arc_angle > 0) ? (-M_PI_2) : M_PI_2);
                incident_angle = angle_mod_2PI(incident_angle);
                incident_edge_list.emplace_back(Incident_Edge_Data{id, true, incident_angle});
            }
            for (const auto id : out_edges) {
                double arc_angle = edges[id].arc.get_arc_angle();
                double angle1 = edges[id].arc.get_start_angle();
                double incident_angle = angle1 + ((arc_angle > 0) ? M_PI_2 : (-M_PI_2));
                incident_angle = angle_mod_2PI(incident_angle);
                incident_edge_list.emplace_back(Incident_Edge_Data{id, false, incident_angle});
            }

            // sort incident edges by their incident angle (clockwise order)
            std::sort(incident_edge_list.begin(), incident_edge_list.end(), Incident_Edge_Data_greater_than);

            // find successor for each incident half-edge
            size_t n_incident = incident_edge_list.size();
            for (int j = 0; j < n_incident; ++j) {
                const auto &e = incident_edge_list[j];
                const auto &f = incident_edge_list[(j+1)%n_incident];
                size_t he = e.is_in ? (2*e.id+1) : 2*e.id;
                size_t hf = f.is_in ? 2*f.id : (2*f.id+1);
                next_hEdge[he] = hf;
            }
        }
    }

    // trace half-edges into chains (each chain bounds an arrangement cell)
    trace_chains(next_hEdge, cells);

}

void Arrangement::trace_chains(const std::vector<size_t> &next_hEdge, std::vector<std::vector<size_t>> &chains)
{
    size_t n_hEdge = next_hEdge.size();
    std::vector<bool> is_visited(n_hEdge, false);

    // trace chains
    chains.clear();
    for (size_t i = 0; i < n_hEdge; ++i) {
        if (!is_visited[i]) {
            is_visited[i] = true;
            chains.emplace_back();
            auto &chain = chains.back();
            chain.push_back(i);
            auto next = next_hEdge[i];
            while (!is_visited[next]) {
                is_visited[next] = true;
                chain.push_back(next);
                next = next_hEdge[next];
            }
        }
    }
}

void Arrangement::compute_cell_windings(const std::vector<SubArc_Edge> &pEdges,
                                        const std::vector<std::vector<size_t>> &eIn,
                                        const std::vector<std::vector<size_t>> &eOut,
                                        const std::vector<std::vector<size_t>> &cells,
                                        std::vector<int> &windings)
{
    // map: half-edge index -> cell index
    std::vector<size_t> cell_of_hE(2*pEdges.size(),0);
    for (int i = 0; i < cells.size(); ++i) {
        for (const auto h : cells[i]) {
            cell_of_hE[h] = i;
        }
    }

    // find the unbounded cell
    size_t unbound_cell_id = find_the_unbounded_cell(pEdges, eIn, eOut, cells, cell_of_hE);
    windings.clear();
    windings.resize(cells.size(), 0);
    windings[unbound_cell_id] = 0;

    // build cell adjacency list
    // adjacent_cells[i][j] is the half-edge between cell i and cell j
    std::vector<std::map<size_t,size_t>> adjacent_cells(cells.size());
    for (int i = 0; i < pEdges.size(); ++i) {
        adjacent_cells[cell_of_hE[2*i+1]][cell_of_hE[2*i]] = 2*i + 1;
        adjacent_cells[cell_of_hE[2*i]][cell_of_hE[2*i+1]] = 2*i;
    }

    // propagate winding numbers on the cell adjacency graph
    std::queue<size_t> Q;
    std::vector<bool> is_visited(cells.size(), false);

    Q.push(unbound_cell_id);
    is_visited[unbound_cell_id] = true;
    while (!Q.empty()) {
        auto cell_id = Q.front();
        Q.pop();
        for (const auto &c : adjacent_cells[cell_id]) {
            auto c_id = c.first;
            auto hEdge_id = c.second;
            if (!is_visited[c_id]) {
                is_visited[c_id] = true;
                if (hEdge_id%2 == 0) {
                    windings[c_id] = windings[cell_id] + 1;
                } else {
                    windings[c_id] = windings[cell_id] - 1;
                }
                Q.push(c_id);
            }
        }
    }

}


size_t Arrangement::find_the_unbounded_cell(const std::vector<SubArc_Edge> &pEdges,
                                            const std::vector<std::vector<size_t>> &eIn,
                                            const std::vector<std::vector<size_t>> &eOut,
                                            const std::vector<std::vector<size_t>> &cells,
                                            const std::vector<size_t> &cell_of_hE)
{  //todo: implement using the other find_the_unbounded_cell function
    // find the arc farthest to the left
    Point most_left_point(std::numeric_limits<double>::infinity(), 0);
    size_t most_left_edge_id;
    Point_Arc_Location most_left_location;

    for (int i = 0; i < pEdges.size(); ++i) {
        std::pair<Point,Point_Arc_Location> left_point_info = pEdges[i].arc.get_most_left_point();
        if (left_point_info.first.x() < most_left_point.x()) {
            most_left_point = left_point_info.first;
            most_left_edge_id = i;
            most_left_location = left_point_info.second;
        }
    }

    // pick out the unbounded cell
    size_t unbounded_cell_id;
    double theta1 = pEdges[most_left_edge_id].arc.get_start_angle();
    double theta2 = pEdges[most_left_edge_id].arc.get_end_angle();

    if (most_left_location == Middle) {
        if (theta1 < theta2) {
            // counter clockwise arc
            unbounded_cell_id = cell_of_hE[2*most_left_edge_id];
        } else {
            unbounded_cell_id = cell_of_hE[2*most_left_edge_id+1];
        }
    } else {
        size_t most_left_vert_id = (most_left_location == Start) ?
                                   pEdges[most_left_edge_id].id1 :
                                   pEdges[most_left_edge_id].id2;
        assert(eIn[most_left_vert_id].size() == 1);

        auto in_edge_id = eIn[most_left_vert_id][0];
        auto out_edge_id = eOut[most_left_vert_id][0];

        auto in_vec = pEdges[in_edge_id].arc.get_out_tangent_vector();
        auto out_vec = pEdges[out_edge_id].arc.get_in_tangent_vector();

        double cross_product = in_vec.x()*out_vec.y() - in_vec.y()*out_vec.x();
        if (cross_product > 0) {
            // left turn
            unbounded_cell_id = cell_of_hE[2*in_edge_id];
        } else if (cross_product < 0) {
            // right turn
            unbounded_cell_id = cell_of_hE[2*in_edge_id+1];
        } else {
            // co-linear
            if ((theta2-theta1)*(in_vec.dot(out_vec)) > 0) {
                unbounded_cell_id = cell_of_hE[2*in_edge_id];
            } else {
                unbounded_cell_id = cell_of_hE[2*in_edge_id+1];
            }
        }
    }

    return unbounded_cell_id;
}

size_t Arrangement::find_the_unbounded_cell(const std::vector<SubArc_Edge> &pEdges,
                                            const std::vector<std::vector<size_t>> &eIn,
                                            const std::vector<std::vector<size_t>> &eOut,
                                            const std::vector<std::vector<size_t>> &cells,
                                            const std::vector<size_t> &cell_of_hE,
                                            Point_Arc_Location most_left_location,
                                            size_t most_left_edge_id)
{
    // pick out the unbounded cell
    size_t unbounded_cell_id;
    double theta1 = pEdges[most_left_edge_id].arc.get_start_angle();
    double theta2 = pEdges[most_left_edge_id].arc.get_end_angle();

    if (most_left_location == Middle) {
        if (theta1 < theta2) {
            // counter clockwise arc
            unbounded_cell_id = cell_of_hE[2*most_left_edge_id];
        } else {
            unbounded_cell_id = cell_of_hE[2*most_left_edge_id+1];
        }
    } else {
        size_t most_left_vert_id = (most_left_location == Start) ?
                                   pEdges[most_left_edge_id].id1 :
                                   pEdges[most_left_edge_id].id2;
        assert(eIn[most_left_vert_id].size() == 1);

        auto in_edge_id = eIn[most_left_vert_id][0];
        auto out_edge_id = eOut[most_left_vert_id][0];

        auto in_vec = pEdges[in_edge_id].arc.get_out_tangent_vector();
        auto out_vec = pEdges[out_edge_id].arc.get_in_tangent_vector();

        double cross_product = in_vec.x()*out_vec.y() - in_vec.y()*out_vec.x();
        if (cross_product > 0) {
            // left turn
            unbounded_cell_id = cell_of_hE[2*in_edge_id];
        } else if (cross_product < 0) {
            // right turn
            unbounded_cell_id = cell_of_hE[2*in_edge_id+1];
        } else {
            // co-linear
            if ((theta2-theta1)*(in_vec.dot(out_vec)) > 0) {
                unbounded_cell_id = cell_of_hE[2*in_edge_id];
            } else {
                unbounded_cell_id = cell_of_hE[2*in_edge_id+1];
            }
        }
    }

    return unbounded_cell_id;
}

void Arrangement::compute_multi_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                            std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
                                            std::vector<bool> &is_intersection_point,
                                            std::vector<int> &arc1_of_intersection,
                                            std::vector<int> &arc2_of_intersection,
                                            std::vector<std::vector<SubArc_Edge>> &edges_of_cell,
                                            std::vector<int> &windings) {
    // step 1: subdivide input arcs by their intersections
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    subdivide_polyArc_by_intersection(vertices, edges,pts, pEdges,is_intersection_point,
                                      arc1_of_intersection,arc2_of_intersection);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    global_subdivide_polyArc_by_intersection_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // step 2: get all arrangement chains (a cell can have many bounding chains)
    t1 = std::chrono::steady_clock::now();
    std::vector<std::vector<size_t>> eIn;
    std::vector<std::vector<size_t>> eOut;
    std::vector<std::vector<size_t>> chains;
    decompose_into_cells(pts, pEdges, eIn, eOut, chains);
    t2 = std::chrono::steady_clock::now();
    global_decompose_into_cells_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // step 3: compute arrangement cells and their winding numbers
    t1 = std::chrono::steady_clock::now();
    std::vector<std::vector<size_t>> cells;
    compute_cells_and_windings(pEdges,eIn,eOut,chains,cells,windings);
//    compute_cell_windings(pEdges, eIn, eOut, cells, windings);
    t2 = std::chrono::steady_clock::now();
    global_compute_cells_and_windings_time += std::chrono::duration_cast<Time_duration>(t2 - t1);

    // convert each cell from a list of half-edge indices to a list arc edges
    edges_of_cell.clear();
    edges_of_cell.reserve(cells.size());
    for (const auto & hEdges : cells) {
        edges_of_cell.emplace_back();
        auto &es = edges_of_cell.back();
        for (auto hE : hEdges) {
            // create SubArc_Edge for half-edge
            // here, SubArc_Edge.parent_id record the index of the input edge in [std::vector<Arc_Edge> edges]
            if (hE%2 == 0) {
                const auto &pEdge = pEdges[hE/2];
                es.emplace_back(SubArc_Edge{pEdge.id2, pEdge.id1, pEdge.parent_id,
                                            Circular_Arc::reverse(pEdge.arc)});
            } else {
                es.emplace_back(pEdges[(hE-1)/2]);
            }
        }
    }

}

void Arrangement::compute_cells_and_windings(const std::vector<SubArc_Edge> &pEdges,
                                             const std::vector<std::vector<size_t>> &eIn,
                                             const std::vector<std::vector<size_t>> &eOut,
                                             const std::vector<std::vector<size_t>> &chains,
                                             std::vector<std::vector<size_t>> &cells, std::vector<int> &windings) {
    // map: half-edge index -> chain index
    std::vector<size_t> chain_of_hE(2*pEdges.size(),0);
    for (int i = 0; i < chains.size(); ++i) {
        for (const auto h : chains[i]) {
            chain_of_hE[h] = i;
        }
    }

    // step 1: build chain adjacency list
    // adjacent_chains[i][j] is the half-edge between chain i and chain j
    std::vector<std::map<size_t,size_t>> adjacent_chains(chains.size());
    for (int i = 0; i < pEdges.size(); ++i) {
        adjacent_chains[chain_of_hE[2*i+1]][chain_of_hE[2*i]] = 2*i + 1;
        adjacent_chains[chain_of_hE[2*i]][chain_of_hE[2*i+1]] = 2*i;
    }

    // step 2: identify connected components of chain adjacency graph
    std::vector<int> component_of_chain(chains.size(),-1);  // -1 for unknown
    int cur_component = -1;
    for (size_t i = 0; i < chains.size(); ++i) {
        if (component_of_chain[i] == -1) { // find new component
            ++cur_component;
            std::queue<size_t> Q;
            Q.push(i);
            component_of_chain[i] = cur_component;
            while (!Q.empty()) {
                auto chain_id = Q.front();
                Q.pop();
                for (const auto &c : adjacent_chains[chain_id]) {
                    auto c_id = c.first;
                    if (component_of_chain[c_id] == -1) {
                        component_of_chain[c_id] = cur_component;
                        Q.push(c_id);
                    }
                }
            }
        }
    }

    std::vector<std::vector<size_t>> components(cur_component+1);
    for (int i = 0; i < chains.size(); ++i) {
        components[component_of_chain[i]].push_back(i);
    }


    // step 3: find the left-most edge of each component

    std::vector<std::pair<Point, Point_Arc_Location>> left_point_info_list;
    left_point_info_list.reserve(pEdges.size());
    for (const auto &pE: pEdges) {
        left_point_info_list.emplace_back(pE.arc.get_most_left_point());
    }

    std::vector<size_t> left_most_edge_of_component(components.size());

    typedef std::pair<size_t ,double> ID_xMin;
    std::vector<ID_xMin> component_id_xmin_list;
    component_id_xmin_list.reserve(components.size());

    for (size_t i = 0; i < components.size(); ++i) {
        // find the arc farthest to the left
        double x_min = std::numeric_limits<double>::infinity();
        size_t most_left_edge_id;
        //
        for (auto ci : components[i]) {
            for (auto he : chains[ci]) {
                size_t e_id = he/2;
                if (left_point_info_list[e_id].first.x() < x_min) {
                    x_min = left_point_info_list[e_id].first.x();
                    most_left_edge_id = e_id;
                }
            }
        }
        left_most_edge_of_component[i] = most_left_edge_id;
        component_id_xmin_list.emplace_back(i,x_min);
    }

    // step 4: sort components in left to right order
    auto ID_xMin_less_than = [](const ID_xMin &left, const ID_xMin &right)
    { return left.second < right.second; };
    std::sort(component_id_xmin_list.begin(), component_id_xmin_list.end(), ID_xMin_less_than);

    std::vector<std::vector<size_t>> sorted_components;
    sorted_components.reserve(components.size());
    std::vector<size_t> left_most_edge_of_sorted_component;
    for (const auto &id_xmin : component_id_xmin_list) {
        sorted_components.emplace_back(components[id_xmin.first]);
        left_most_edge_of_sorted_component.push_back(left_most_edge_of_component[id_xmin.first]);
    }

    // step 5: find unbounded chain of each component
    std::vector<size_t> unbounded_chain_of_sorted_component(sorted_components.size());
    for (int i = 0; i < sorted_components.size(); ++i) {
        unbounded_chain_of_sorted_component[i] =
                find_the_unbounded_cell(pEdges,eIn,eOut,chains,chain_of_hE,
                                        left_point_info_list[left_most_edge_of_sorted_component[i]].second,
                                        left_most_edge_of_sorted_component[i]);
    }

    // step 6: compute arrangement cells and their adjacency

    std::vector<std::vector<size_t>> new_chains = chains;

    std::vector<size_t> final_component = sorted_components[0];
    size_t unbounded_chain_id = unbounded_chain_of_sorted_component[0];

    for (int i = 1; i < sorted_components.size(); ++i) {
        const auto &component = sorted_components[i];

        // pick a point on the current component
        auto he = chains[component[0]][0];
        size_t e_id = he/2;
        Point p = pEdges[e_id].arc.get_start_point();

        // find the chain that encloses the current component
        bool enclosing_chain_found = false;
        size_t enclosing_chain_id;
        for (auto ci : final_component) {
            double total_signed_angle = 0;
            for (auto half_edge : new_chains[ci]) {
                size_t edge_id = half_edge / 2;
                double angle = pEdges[edge_id].arc.get_view_angle(p);
                if (half_edge%2 == 0) { // opposite direction half-edge
                    angle *= -1;
                }
                total_signed_angle += angle;
            }
            if (round(0.5*total_signed_angle/M_PI) == 1) {
                enclosing_chain_id = ci;
                enclosing_chain_found = true;
                break;
            }
        }
        if (!enclosing_chain_found) {
            enclosing_chain_id = unbounded_chain_id;
        }

        // merge the current component into the final component
        auto cur_unbounded_chain_id = unbounded_chain_of_sorted_component[i];
        new_chains[enclosing_chain_id].insert(new_chains[enclosing_chain_id].end(),
                                              chains[cur_unbounded_chain_id].begin(),
                                              chains[cur_unbounded_chain_id].end());

        for (const auto &ci : component) {
            if (ci != cur_unbounded_chain_id) {
                final_component.push_back(ci);
            }
        }

        for (const auto &cId_he :adjacent_chains[cur_unbounded_chain_id]) {
            size_t c_id = cId_he.first;
            size_t half_edge1 = cId_he.second;
            size_t half_edge2 = adjacent_chains[c_id][cur_unbounded_chain_id];
            // new edges from the enclosing chain
            adjacent_chains[enclosing_chain_id][c_id] = half_edge1;
            // new edges to the enclosing chain
            adjacent_chains[c_id].erase(cur_unbounded_chain_id);
            adjacent_chains[c_id][enclosing_chain_id] = half_edge2;
        }

    }

    // step 7: compute cells and their winding numbers

    std::vector<int> chain_windings(new_chains.size(), 0);
    chain_windings[unbounded_chain_id] = 0;

    // propagate winding numbers on the cell adjacency graph
    std::queue<size_t> Q;
    std::vector<bool> is_visited(new_chains.size(), false);

    Q.push(unbounded_chain_id);
    is_visited[unbounded_chain_id] = true;
    while (!Q.empty()) {
        auto chain_id = Q.front();
        Q.pop();
        for (const auto &c : adjacent_chains[chain_id]) {
            auto c_id = c.first;
            auto hEdge_id = c.second;
            if (!is_visited[c_id]) {
                is_visited[c_id] = true;
                if (hEdge_id%2 == 0) {
                    chain_windings[c_id] = chain_windings[chain_id] + 1;
                } else {
                    chain_windings[c_id] = chain_windings[chain_id] - 1;
                }
                Q.push(c_id);
            }
        }
    }

    cells.clear();
    for (const auto & chain_id: final_component) {
        cells.emplace_back(new_chains[chain_id]);
    }

    windings.clear();
    for (const auto & chain_id: final_component) {
        windings.push_back(chain_windings[chain_id]);
    }

}
