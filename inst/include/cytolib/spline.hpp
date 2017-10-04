/*
 * spline.hpp
 *
 *  Created on: Apr 26, 2012
 *      Author: wjiang2
 */


#ifndef SPLINE_HPP_

#define SPLINE_HPP_
#ifdef ROUT
#define COUT cout //TODO: originally it was Rout, to get around Rcpp, we will have to use Rprintf/printf to rewrite all the console print
#endif


#ifndef ROUT
#define COUT cout
#endif

#include "global.hpp"
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
