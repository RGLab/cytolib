/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * calibrationTable.hpp
 *
 *  Created on: May 14, 2012
 *      Author: wjiang2
 */

#ifndef CALIBRATIONTABLE_HPP_
#define CALIBRATIONTABLE_HPP_

#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include "spline.hpp"
using namespace std;
#include <cytolib/GatingSet.pb.h>

namespace cytolib
{
struct Spline_Coefs{
	map<string,vector<double> > coefs;
	int method;
	string type;//to be deprecated
};

class calibrationTable{
private:
	vector<double> x,y,b,c,d;
	int spline_method;
	string caltype;//TODO:move this to transformation class
	bool flag;
public:
	vector<double> getX();
	vector<double> getY();
	void setY(vector<double> _y);
	void setX(vector<double> _x);
	vector<double> getB();
	vector<double> getC();
	vector<double> getD();
	void setCaltype(string _caltype);
	string getCaltype();
	void setMethod(int _spline_method);
	int getMethod();
	void setInterpolated(bool _flag);
	bool isInterpolated();
	calibrationTable();
	calibrationTable(string _caltype,int _spline_method);
	void interpolate();
	/*
	 * API provided for Rcpp to access calibration table
	 */
	Spline_Coefs getSplineCoefs();
	void transforming(double * input, int nSize);
	void convertToPb(pb::calibrationTable & cal_pb);

	calibrationTable(const pb::calibrationTable & cal_pb);
};
};


#endif /* CALIBRATIONTABLE_HPP_ */
