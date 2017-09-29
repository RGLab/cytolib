/*
 * spline.hpp
 *
 *  Created on: Apr 26, 2012
 *      Author: wjiang2
 */


#ifndef SPLINE_HPP_

#define SPLINE_HPP_

#include "config.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <cstring>
#define NATURAL 2

using namespace std;
void natural_spline_C(int n, double *x, double *y, double *b, double *c, double *d);

void spline_eval_C(int *method, int *nu, double *u, double *v,
		 int *n, double *x, double *y, double *b, double *c, double *d);


void natural_spline(vector<double>x, vector<double> y, vector<double>& b,vector<double>& c,vector<double>& d);
void spline_eval(int method, double* u,int nSize,const vector<double> & x, const vector<double> & y, const vector<double> & b, const vector<double> & c, const vector<double> & d);
#endif /* SPLINE_HPP_ */
