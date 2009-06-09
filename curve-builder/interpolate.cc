#include <vector>

#include "lu.hh"
#include "interpolate.hh"

DoubleVector
GenerateParameters(PointVector const &v)
{
  size_t const length = v.size();
  DoubleVector result;

  double d = 0.0;
  for(size_t k = 1; k < length; ++k)
    d += PointDistance(v[k - 1], v[k]);
  d = 1.0 / d;

  result.resize(length);
  result[0] = 0.0;
  result[length - 1] = 1.0;

  for(size_t k = 1; k < length - 1; ++k)
    result[k] = result[k - 1] + PointDistance(v[k - 1], v[k]) * d;

  return result;
}

DoubleVector
GenerateKnots(DoubleVector const &parameters, int p)
{
  size_t const plength = parameters.size();
  size_t const klength = plength + p + 1;
  DoubleVector result;

  result.resize(klength);
  for(size_t k = 0; k <= p; ++ k) {
    result[k] = 0.0;
    result[klength - k - 1] = 1.0;
  }

  for(size_t j = 1; j < plength - p; ++j) {
    double sum = 0.0;
    for(size_t i = j; i <= j + p - 1; ++i)
      sum += parameters[i];
    result[j + p] = sum / p;
  }

  return result;
}

int
BSplineSpan(DoubleVector const &knots, int p, double u)
{
  size_t m = knots.size();
  size_t n = m - p - 1;
  if(u == knots[n])
    return n - 1;

  size_t low = p, high = m;
  size_t mid = (low + high) / 2;
  while(u < knots[mid] || u >= knots[mid + 1]) {
    if(u < knots[mid])
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2;
  }
  return mid;
}

DoubleVector
BSplineBasis(DoubleVector const &knots, int i, int p, double u)
{
  DoubleVector left, right, result;

  result.resize(p + 1); result[0] = 1.0;

  for(int j = 1; j <= p; ++j) {
    left.push_back(u - knots[i + 1 - j]);
    right.push_back(knots[i + j] - u);
    double saved = 0.0;
    for(int r = 0; r < j; ++r) {
      double temp = result[r] / (right[r] + left[j - r - 1]);
      result[r] = saved + right[r] * temp;
      saved = left[j - r - 1] * temp;
    }
    result[j] = saved;
  }

  return result;
}

PointVector
GenerateControlPoints(PointVector const &points, DoubleVector const &parameters,
		      DoubleVector const &knots, int p)
{
  size_t const n = parameters.size();
  DoubleVector a;

  a.resize(n * n);
  for(size_t i = 0; i < n; ++i) {
    for(size_t j = 0; j < n; ++j)
      a[i * n + j] = 0.0;
    int span = BSplineSpan(knots, p, parameters[i]);
    DoubleVector basis = BSplineBasis(knots, span, p, parameters[i]);
    for(size_t j = 0; j <= p; ++j)
      a[i * n + span - p + j] = basis[j];
  }

  std::vector<int> index;
  LUDecompose(a, n, index);

  DoubleVector temp;
  for(size_t i = 0; i < n; ++i) temp.push_back(points[i].x);
  DoubleVector x = LUBackSubstitution(a, n, index, temp);
  for(size_t i = 0; i < n; ++i) temp[i] = points[i].y;
  DoubleVector y = LUBackSubstitution(a, n, index, temp);
  for(size_t i = 0; i < n; ++i) temp[i] = points[i].z;
  DoubleVector z = LUBackSubstitution(a, n, index, temp);

  PointVector result;
  for(size_t i = 0; i < n; ++i) result.push_back(Point(x[i], y[i], z[i]));

  return result;
}

BSpline
BSplineInterpolate(PointVector const &points, int degree)
{
  BSpline result;

  DoubleVector parameters = GenerateParameters(points);
  result.degree = degree;
  result.knots = GenerateKnots(parameters, degree);
  result.points = GenerateControlPoints(points, parameters, result.knots, degree);

  return result;
}
