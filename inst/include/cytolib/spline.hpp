/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * spline.hpp
 *
 *  Created on: Apr 26, 2012
 *      Author: wjiang2
 */


#ifndef SPLINE_HPP_

#define SPLINE_HPP_
#include <vector>
using namespace std;

namespace cytolib
{
#define NATURAL 2
/*TODO:change it to c++ version
 * n :number of input points
 * x,y:coordinates of input points
 * b,c,d:output with the same length of x
 */
void natural_spline_C(int n, double *x, double *y, double *b, double *c, double *d);

/*
 * *method=2 is natural method
 * nu is the length of input vector
 * u: is the input vector x
 * v:is the output vector y
 * n,x,y, b,c,d: from spline_coef
 */
void spline_eval_C(int *method, int *nu, double *u, double *v,
		 int *n, double *x, double *y, double *b, double *c, double *d);
/*
 * vector version
 */

void natural_spline(vector<double>x, vector<double> y, vector<double>& b,vector<double>& c,vector<double>& d);
void spline_eval(int method, double* u,int nSize,
		  const vector<double> & x, const vector<double> & y, const vector<double> & b, const vector<double> & c, const vector<double> & d);
};
#endif /* SPLINE_HPP_ */
