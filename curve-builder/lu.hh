#ifndef LU_HH
#define LU_HH

#include <vector>

double
LUDecompose(std::vector<double> &a, int n, std::vector<int> &indx);

std::vector<double>
LUBackSubstitution(std::vector<double> const &a, int n,
		   std::vector<int> const &indx, std::vector<double> b);

std::vector<double>
LUSolveLEQS(std::vector<double> const &a, std::vector<double> const &b);

#endif	// LU_HH
