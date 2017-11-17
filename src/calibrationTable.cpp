/*
 * calibrationTable.cpp
 *
 *  Created on: May 14, 2012
 *      Author: wjiang2
 */

#include <cytolib/calibrationTable.hpp>
#include <fstream>
calibrationTable::calibrationTable(){
	flag=false;
}

calibrationTable::calibrationTable(string _caltype,int _spline_method){
//								type=CALTBL;
		caltype=_caltype;
		spline_method=_spline_method;
		flag=false;
}



void calibrationTable::interpolate(){

//	COUT<<"entering interpolate"<<endl;


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
Spline_Coefs calibrationTable::getSplineCoefs(){

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
void calibrationTable::transforming(double * input, int nSize){


	int imeth=2;

	spline_eval(imeth,input, nSize, x, y, b, c, d);

}


void calibrationTable::convertToPb(pb::calibrationTable & cal_pb){
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

calibrationTable::calibrationTable(const pb::calibrationTable & cal_pb){
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
