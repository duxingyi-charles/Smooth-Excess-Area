//
// Created by Charles Du on 4/21/21.
//

#include "Rectangle.h"

bool Rectangle::is_intersect(const Rectangle &rect1, const Rectangle &rect2) {
    return (rect1.p_max.x() >= rect2.p_min.x()) && (rect1.p_min.x() <= rect2.p_max.x())
            && (rect1.p_max.y() >= rect2.p_min.y()) && (rect1.p_min.y() <= rect2.p_max.y());
}