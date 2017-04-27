/* Copyright (C) 2016-2018, Stanford University
 * This file is part of MESH
 * Written by Kaifeng Chen (kfchen@stanford.edu)
 *
 * MESH is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MESH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef _GEOMETRY_FUNCTIONS
#define _GEOMETRY_FUNCTIONS

#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#define INF 1e10
/*======================================================*/
// Functions related to grating
/*=======================================================*/
inline double getGratingArea(const double a){
  return a;
}
// check whether center2 is contained in center1
inline bool isContainedInGrating(const double center1, const double center2, const double width){
  return std::abs(center2 - center1) <= width / 2;
}
/*======================================================*/
// Functions related to rectangle
/*=======================================================*/
inline double getRectangleArea(const double width1, const double width2){
  return width1 * width2;
}
inline bool isContainedInRectangle(const double center1[2], const double center2[2], const double width[2]){
  return (std::abs(center2[0] - center1[0]) <= width[0] / 2) && (std::abs(center2[1] - center1[1]) <= width[1] / 2);
}
/*======================================================*/
// Functions related to circle
/*=======================================================*/
inline double getCircleArea(const double r){
  return M_PI * r * r;
}
inline bool isContainedInCircle(const double center1[2], const double center2[2], const double r){
  return std::pow(center2[0] - center1[0], 2) + std::pow(center2[1] - center1[1], 2) <= std::pow(r, 2);
}
/*======================================================*/
// Functions related to ellipse
/*=======================================================*/
inline double getEllipseArea(const double a, const double b){
  return M_PI * a * b;
}

inline bool isContainedInEllipse(const double center1[2], const double center2[2], const double a, const double b){
  return std::pow((center2[0] - center1[0])/a, 2) + std::pow((center2[1] - center1[1])/b, 2) <= 1;
}
/*======================================================*/
// Functions related to polygon
/*=======================================================*/
typedef std::pair<double, double> Point;
typedef std::vector< Point > EdgeList;

inline double crossProduct(const Point& point1, const Point& point2){
  return point1.first * point2.second - point1.second * point2.first;
}

inline double getPolygonArea(const EdgeList& edgeList){
  double area = 0;
  for(size_t i = 0; i < edgeList.size() - 1; i++){
    area += crossProduct(edgeList[i], edgeList[i + 1]);
  }
  area += crossProduct(edgeList[edgeList.size() - 1], edgeList[0]);
  return std::abs(area) / 2;
}

// function checking whether q lies on segment p-r, assuming pqr are colinear
inline bool onSegment(const Point& p, const Point& q, const Point& r){
  if(q.first <= std::max(p.first, r.first)
    && q.first >= std::min(p.first, r.first)
    && q.second <= std::max(p.second, r.second)
    && q.second >= std::min(p.second, r.second)){
      return true;
  }
  return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int orientation(const Point& p, const Point& q, const Point& r){
  double val = (q.second - p.second) * (r.first - q.first) - (q.first - p.first) * (r.second - q.second);
  if(val == 0.0) return 0;
  return (val > 0)? 1:2;
}

// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
inline bool doIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2){
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
    return false; // Doesn't fall in any of the above cases
}
/*
inline bool isContainedInPolygon(const double center1[2], const double center2[2], const EdgeList& edgeList){
  Point pointToCheck = std::make_pair(center2[0] - center1[0], center2[1] - center1[1]);
  Point pointExtreme = std::make_pair(INF, pointToCheck.second);
  int count = 0, i = 0, n = edgeList.size();
  do{
    int next = (i+1) % n;
    if(doIntersect(edgeList[i], edgeList[next], pointToCheck, pointExtreme)){
      std::cout << orientation(edgeList[i], pointToCheck, edgeList[next]) << std::endl;
        if(orientation(edgeList[i], pointToCheck, edgeList[next]) == 0){
          return onSegment(edgeList[i], pointToCheck, edgeList[next]);
        }
        count++;
    }
    i = next;
  } while(i != 0);

  return count & 1;
}
*/
inline bool isContainedInPolygon(const double center1[2], const double center2[2], const EdgeList& edgeList){
  bool isInside = false;
  Point pointToCheck = std::make_pair(center2[0] - center1[0], center2[1] - center1[1]);
  for(size_t i = 0, j = edgeList.size() - 1; i < edgeList.size(); j = i++){
    if( ((edgeList[i].second > pointToCheck.second) != (edgeList[j].second > pointToCheck.second)) &&
      (pointToCheck.first < (edgeList[j].first - edgeList[i].first) * (pointToCheck.second - edgeList[i].second)/(edgeList[j].second - edgeList[i].second) + edgeList[i].first) ){
        isInside = !isInside;
      }
  }
  return isInside;
}

#endif