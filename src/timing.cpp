//
// Created by Charles Du on 8/15/21.
//

#include "timing.h"
#include <iostream>

Time_duration global_TLC_time;
Time_duration global_arc_seg_time;
Time_duration global_arc_occupancy_time;
Time_duration global_arc_arrangement_time;

Time_duration global_subdivide_polyArc_by_intersection_time;
Time_duration global_arc_intersection_init_time;
Time_duration global_arc_bbox_intersection_test_time;
Time_duration global_arc_arc_intersection_time;
Time_duration global_arc_subdivision_time;

Time_duration global_decompose_into_cells_time;
Time_duration global_compute_cells_and_windings_time;

void reset_timings() {
    global_TLC_time = Time_duration::zero();
    global_arc_seg_time = Time_duration::zero();
    global_arc_occupancy_time = Time_duration::zero();
    global_arc_arrangement_time = Time_duration::zero();

    global_subdivide_polyArc_by_intersection_time = Time_duration::zero();
    global_arc_intersection_init_time = Time_duration::zero();
    global_arc_bbox_intersection_test_time = Time_duration::zero();
    global_arc_arc_intersection_time = Time_duration::zero();
    global_arc_subdivision_time = Time_duration::zero();

    global_decompose_into_cells_time = Time_duration::zero();
    global_compute_cells_and_windings_time = Time_duration::zero();
}

void print_timings() {
    std::cout << "TLC time: " << global_TLC_time.count() << " seconds." << std::endl;
    std::cout << "arc segment time: " << global_arc_seg_time.count() << " seconds." << std::endl;
    std::cout << "arc occupancy time: " << global_arc_occupancy_time.count() << " seconds." << std::endl;
    std::cout << "- arc arrangement time: " << global_arc_arrangement_time.count() << " seconds." << std::endl;
    std::cout << "** subdivide polyArc time: " << global_subdivide_polyArc_by_intersection_time.count() << " seconds." << std::endl;
    std::cout << "### arc intersection init time: " << global_arc_intersection_init_time.count() << " seconds." << std::endl;
    std::cout << "### arc bbox intersection test time: " << global_arc_bbox_intersection_test_time.count() << " seconds." << std::endl;
    std::cout << "### arc-arc intersection time: " << global_arc_arc_intersection_time.count() << " seconds." << std::endl;
    std::cout << "### arc subdivision time: " << global_arc_subdivision_time.count() << " seconds." << std::endl;
    std::cout << "** decompose into cells time: " << global_decompose_into_cells_time.count() << " seconds." << std::endl;
    std::cout << "** compute cells windings time: " << global_compute_cells_and_windings_time.count() << " seconds." << std::endl;
}
