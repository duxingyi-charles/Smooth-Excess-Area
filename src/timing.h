//
// Created by Charles Du on 8/15/21.
//

#ifndef TLC_TIMING_H
#define TLC_TIMING_H

#include <chrono>
typedef std::chrono::duration<double> Time_duration;

// computing TLC func/grad
extern Time_duration global_TLC_time;
// computing arc segment area func/grad
extern Time_duration global_arc_seg_time;
// computing arc occupancy func/grad
extern Time_duration global_arc_occupancy_time;

// computing arrangement by boundary arcs (this is part of global_arc_occupancy_time)
extern Time_duration global_arc_arrangement_time;

extern Time_duration global_subdivide_polyArc_by_intersection_time;
extern Time_duration global_arc_intersection_init_time;
extern Time_duration global_arc_bbox_intersection_test_time;
extern Time_duration global_arc_arc_intersection_time;
extern Time_duration global_arc_subdivision_time;

extern Time_duration global_decompose_into_cells_time;
extern Time_duration global_compute_cells_and_windings_time;


void reset_timings();

void print_timings();



#endif //TLC_TIMING_H
