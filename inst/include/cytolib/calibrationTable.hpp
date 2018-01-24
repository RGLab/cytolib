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
#include <GatingSet.pb.h>

namespace cytolib
{
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
	calibrationTable(){
		flag=false;
	}

	calibrationTable(string _caltype,int _spline_method){
	//								type=CALTBL;
			caltype=_caltype;
			spline_method=_spline_method;
			flag=false;
	}



	void interpolate(){

	//	PRINT("entering interpolate\n");


		if(!flag)
		{
			b.resize(x.size());
			c.resize(x.size());
			d.resize(x.size());
			natural_spline(x, y, b, c, d);
			flag=true;
		}


	}
	/*
	 * API provided for Rcpp to access calibration table
	 */
	Spline_Coefs getSplineCoefs(){

		map<string,vector<double> > coefs;


		coefs["x"]=x;
		coefs["y"]=y;
		coefs["b"]=b;
		coefs["c"]=c;
		coefs["d"]=d;

		Spline_Coefs res;
		res.coefs=coefs;
		res.method=spline_method;
		res.type=caltype;

		return res;
	}
	void transforming(double * input, int nSize){


		int imeth=2;

		spline_eval(imeth,input, nSize, x, y, b, c, d);

	}


	void convertToPb(pb::calibrationTable & cal_pb){
		if(!isInterpolated())
			interpolate();
		for(unsigned i = 0; i < x.size(); i++){
			cal_pb.add_x(x[i]);
			cal_pb.add_y(y[i]);
			cal_pb.add_b(b[i]);
			cal_pb.add_c(c[i]);
			cal_pb.add_d(d[i]);
		}
		cal_pb.set_spline_method(spline_method);
		cal_pb.set_caltype(caltype);
		cal_pb.set_flag(flag);
	}

	calibrationTable(const pb::calibrationTable & cal_pb){
		int nSize = cal_pb.x_size();
		x.resize(nSize);
		y.resize(nSize);
		b.resize(nSize);
		c.resize(nSize);
		d.resize(nSize);
		for(int i = 0; i < nSize; i++){
			x[i] = cal_pb.x(i);
			y[i] = cal_pb.y(i);
			b[i] = cal_pb.b(i);
			c[i] = cal_pb.c(i);
			d[i] = cal_pb.d(i);
		}
		spline_method = cal_pb.spline_method();
		caltype = cal_pb.caltype();
		flag = cal_pb.flag();
	}

};
};


#endif /* CALIBRATIONTABLE_HPP_ */
