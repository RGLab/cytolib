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

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <valarray>
#define NATURAL 2

using namespace std;
void natural_spline_C(int n, double *x, double *y, double *b, double *c, double *d);

void spline_eval_C(int *method, int *nu, double *u, double *v,
		 int *n, double *x, double *y, double *b, double *c, double *d);


void natural_spline(valarray<double>x, valarray<double> y, valarray<double>& b,valarray<double>& c,valarray<double>& d);
void spline_eval(int method, valarray<double> u, valarray<double> & v,const valarray<double> & x, const valarray<double> & y, const valarray<double> & b, const valarray<double> & c, const valarray<double> & d);
#endif /* SPLINE_HPP_ */
