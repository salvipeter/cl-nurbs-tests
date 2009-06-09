#include <cmath>
#include <iostream>
#include <limits>

#include "lu.hh"

double
LUDecompose(std::vector<double> &a, int n, std::vector<int> &indx)
{
  double temp, d = 1.0;
  std::vector<double> vv;

  indx.resize(n);
  vv.resize(n);
  for(int i = 0; i < n; ++i) {
    double big = 0.0;
    for(int j = 0; j < n; ++j)
      if ((temp = std::abs(a[i * n + j])) > big) big = temp;
    if (big == 0.0) { std::cerr << "Singular matrix!" << std::endl; return 0; }
    vv[i] = 1.0 / big;
  }
  for(int j = 0; j < n; ++j) {
    for(int i = 0; i < j; ++i) {
      double sum = a[i * n + j];
      for(int k = 0; k < i; ++k) sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
    }
    double big = 0.0; int imax;
    for(int i = j; i < n; ++i) {
      double sum = a[i * n + j];
      for(int k = 0; k < j; ++k)
	sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
      if((temp = vv[i] * std::abs(sum)) >= big) {
	big = temp;
	imax = i;
      }
    }
    if(j != imax) {
      for(int k = 0; k < n; ++k)
	std::swap(a[imax * n + k], a[j * n + k]);
      std::swap(vv[imax], vv[j]);
      d = -d;
    }
    indx[j] = imax;
    if(a[j * n + j] == 0.0) a[j * n + j] = std::numeric_limits<double>::min();
    if(j != n - 1) {
      double dum = 1.0 / (a[j * n + j]);
      for(int i = j + 1; i < n; ++i) a[i * n + j] *= dum;
    }
  }
  return d;
}

std::vector<double>
LUBackSubstitution(std::vector<double> const &a, int n,
		   std::vector<int> const &indx, std::vector<double> b)
{
  for(int i = 0, ii = -1; i < n; ++i) {
    int ip = indx[i];
    double sum = b[ip];
    b[ip] = b[i];
    if(ii >= 0)
      for(int j = ii; j <= i-1; ++j) sum -= a[i * n + j] * b[j];
    else if(sum) ii = i;
    b[i] = sum;
  }
  for(int i = n - 1; i >= 0; --i) {
    double sum = b[i];
    for(int j = i+1; j < n; ++j) sum -= a[i * n + j] * b[j];
    b[i] = sum / a[i * n + i];
  }
  return b;
}

std::vector<double>
LUSolveLEQS(std::vector<double> const &a, std::vector<double> const &b)
{
  std::vector<double> m = a;
  std::vector<int> i;
  LUDecompose(m, b.size(), i);
  return LUBackSubstitution(m, b.size(), i, b);
}
