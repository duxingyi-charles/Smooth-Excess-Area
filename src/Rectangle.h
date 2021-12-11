//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_RECTANGLE_H
#define TLC_RECTANGLE_H

#include "geo_util.h"

class Rectangle {
public:
    // default: unit square
    Rectangle() : p_min(0,0), p_max(1,1) {};
    Rectangle(const Point &p1, const Point &p2) : p_min(p1), p_max(p2) {};

    ~Rectangle() = default;

public:
    static bool is_intersect(const Rectangle &rect1, const Rectangle &rect2);


private:
    Point p_min;
    Point p_max;

};


#endif //TLC_RECTANGLE_H
