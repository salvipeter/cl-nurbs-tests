#ifndef POINT_HH
#define POINT_HH

#include <cmath>
#include <vector>

struct Point {
  double x;
  double y;
  double z;
  Point(double a, double b) { x = a; y = b; z = 0.0; }
  Point(double a, double b, double c) { x = a; y = b; z = c; }
};
typedef std::vector<Point> PointVector;

inline double sqr(double x) { return x * x; }

inline double PointDistance(Point a, Point b)
{
  return std::sqrt(sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z));
}

#endif	// POINT_HH
