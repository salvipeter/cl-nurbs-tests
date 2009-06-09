#ifndef INTERPOLATE_HH
#define INTERPOLATE_HH

#include "point.hh"

typedef std::vector<double> DoubleVector;

struct BSpline {
  int degree;
  DoubleVector knots;
  PointVector points;
};

BSpline BSplineInterpolate(PointVector const &points, int degree);

#endif	// INTERPOLATE_HH
