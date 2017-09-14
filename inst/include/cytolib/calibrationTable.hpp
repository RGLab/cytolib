/*
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
#include <boost/config.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/foreach.hpp>
#include "GatingSet.pb.h"

struct Spline_Coefs{
	map<string,vector<double> > coefs;
	int method;
	string type;//to be deprecated
};

class calibrationTable{
	friend std::ostream & operator<<(std::ostream &os, const calibrationTable &gh);
private:
	vector<double> x,y,b,c,d;
	int spline_method;
	string caltype;//TODO:move this to transformation class
	bool flag;
public:
	calibrationTable();
	calibrationTable(string _caltype,int _spline_method);
	void interpolate();
	void transforming(double * input, int nSize);
	Spline_Coefs getSplineCoefs();
	vector<double> getX(){return x;};
	vector<double> getY(){return y;};
	void setY(vector<double> _y){
			y=_y;
			};
	void setX(vector<double> _x){
				x=_x;
				};
	vector<double> getB(){return b;};
	vector<double> getC(){return c;};
	vector<double> getD(){return d;};
	void setCaltype(string _caltype){caltype=_caltype;};
	string getCaltype(){return caltype;};
	void setMethod(int _spline_method){spline_method=_spline_method;};
	int getMethod(){return spline_method;};
	void setInterpolated(bool _flag){flag=_flag;};
	bool isInterpolated(){return flag;}
	void convertToPb(pb::calibrationTable & cal_pb);
	calibrationTable(const pb::calibrationTable & cal_pb);
};


#endif /* CALIBRATIONTABLE_HPP_ */
