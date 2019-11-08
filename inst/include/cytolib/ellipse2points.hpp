/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * ellipse2points.hpp
 *
 *  Created on: Feb 16, 2017
 *      Author: wjiang2
 */

#ifndef INCLUDE_ELLIPSE2POINTS_HPP_
#define INCLUDE_ELLIPSE2POINTS_HPP_

#include <vector>
using namespace std;
namespace cytolib
{
struct ellipse_parsed{

  float mu_x, mu_y, a, b, alpha;
//  vector<float>x,y;

};
struct matrix{
	vector<float> x;
	vector<float> y;
};
ellipse_parsed parseEllipse(vector<float> x, vector<float> y);
/**
 * translated from flowClust:::.ellipsePoints R code
 * @param res
 * @param n
 * @return
 */
matrix toPoly(ellipse_parsed res, int n);
};


#endif /* INCLUDE_ELLIPSE2POINTS_HPP_ */
