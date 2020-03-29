/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * transformation.hpp
 *
 *  Created on: Apr 24, 2012
 *      Author: wjiang2
 */

#ifndef TRANSFORMATION_HPP_
#define TRANSFORMATION_HPP_
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include "calibrationTable.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/compare.hpp>


#include "in_polygon.hpp"

using namespace std;

namespace cytolib
{
#define CALTBL 0
#define LOG 1
#define LIN 2
#define FLIN 3
#define FASINH 4
#define BIEXP 5
#define LOGICLE 6
#define LOGGML2 7
#define SCALE 8
/* case insensitive compare predicate*/
struct ciLessBoost : std::binary_function<std::string, std::string, bool>
{
    bool operator() (const std::string & s1, const std::string & s2) const {
        return lexicographical_compare(s1, s2, boost::is_iless());
    }
};

typedef map<std::string, std::string, ciLessBoost> CHANNEL_MAP;

struct coordinate : public cytolib::CYTO_POINT
{
	coordinate(EVENT_DATA_TYPE _x,EVENT_DATA_TYPE _y):CYTO_POINT(_x, _y){};//{x=_x;y=_y;};
	coordinate(){};
	void convertToPb(pb::coordinate & coor_pb){
		coor_pb.set_x(x);
		coor_pb.set_y(y);
	};
	coordinate(const pb::coordinate & coor_pb):CYTO_POINT(coor_pb.x(),coor_pb.y()){};
};

class transformation;
typedef shared_ptr<transformation> TransPtr;
class transformation{
protected:
	calibrationTable calTbl; //no longer have to store calTbl since we now can compute it from biexp even for mac workspace
	bool isGateOnly;
	unsigned short type;//could have been avoided if it is not required by R API getTransformation that needs to extract concrete transformation
	string name;
	string channel;
	bool isComputed;//this flag allow lazy computCalTbl/interpolation
public:
	/*
		 * if it is pure transformation object,then assume calibration is directly read from ws
		 * so there is no need to compute calibration
		 */

	transformation();
	transformation(bool _isGate, unsigned short _type);
	virtual ~transformation(){};
	virtual void transforming(EVENT_DATA_TYPE * input, int nSize);

	virtual void computCalTbl();//dummy routine that does nothing
	virtual Spline_Coefs getSplineCoefs();
	virtual void setCalTbl(calibrationTable _tbl);

	virtual calibrationTable getCalTbl();
	virtual void interpolate();
	virtual bool isInterpolated();
	virtual bool gateOnly();
	virtual void setGateOnlyFlag(bool _flag);
	virtual bool computed();
	virtual void setComputeFlag(bool _flag);
	virtual string getName();
	virtual void setName(string _name);
	virtual string getChannel();
	virtual void setChannel(string _channel);
	virtual unsigned short getType();
	virtual unsigned short getType(string &ctype);
	virtual void setType(unsigned short _type);
	virtual TransPtr clone() const;
	transformation(const pb::transformation & trans_pb);
	virtual void convertToPb(pb::transformation & trans_pb);
	virtual TransPtr  getInverseTransformation();
	virtual void setTransformedScale(int scale);
	virtual int getTransformedScale();
	virtual int getRawScale();
};
typedef shared_ptr<transformation> TransPtr;


class biexpTrans:public transformation{
public:
	int channelRange;
	EVENT_DATA_TYPE pos, neg, widthBasis, maxValue;
public:

	biexpTrans();

	biexpTrans(int channelRange_,EVENT_DATA_TYPE pos_, EVENT_DATA_TYPE neg_, EVENT_DATA_TYPE widthBasis_, EVENT_DATA_TYPE maxValue_);
	/*
	 * directly translated from java routine from tree star
	 */
	EVENT_DATA_TYPE logRoot(EVENT_DATA_TYPE b, EVENT_DATA_TYPE w);

	 /*
	  * directly translated from java routine from tree star
	  * potential segfault risk: the inappropriate biexp parameters can cause
	  * the indexing of vector out of the boundary.
	  */

	void computCalTbl();
	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	biexpTrans(const pb::transformation & trans_pb);
	void setTransformedScale(int scale);
	int getTransformedScale();
	int getRawScale();

};

class fasinhTrans:public transformation{
public:
	EVENT_DATA_TYPE maxRange;
	EVENT_DATA_TYPE length;//unused at this moment
	EVENT_DATA_TYPE T, A, M;
public:
	fasinhTrans();

	fasinhTrans(EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _T, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M);
	/*
	 * implementation copied from flowCore
	 */
	virtual void transforming(EVENT_DATA_TYPE * input, int nSize);

	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	fasinhTrans(const pb::transformation & trans_pb);
	TransPtr  getInverseTransformation();

	void setTransformedScale(int scale);
	int getTransformedScale();
	int getRawScale();
};
/*
 * inverse transformation of fasinhTrans
 */
class fsinhTrans:public fasinhTrans{

public:

	fsinhTrans();

	fsinhTrans(EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _T, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M);

	void  transforming(EVENT_DATA_TYPE * input, int nSize);
	TransPtr getInverseTransformation();
};

/*
 * TODO:right now set two flags to TRUE in the contructor to avoid doing cal table stuff,
 * we should consider redesign the classes so that logTrans does not share this extra feature from parent class
 */
class logTrans:public transformation{
public:
		EVENT_DATA_TYPE offset;//i.e. min in flowjo
		EVENT_DATA_TYPE decade; //i.e. log(max) - log(min) in flowjo
		unsigned scale;
		unsigned T; //(deprecated)top value; derived from keyword $PnR for each channel
public:
	logTrans();

	logTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T);

	void transforming(EVENT_DATA_TYPE * input, int nSize);
	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	logTrans(const pb::transformation & trans_pb);

	TransPtr  getInverseTransformation();

	void setTransformedScale(int _scale);
	int getTransformedScale();
	int getRawScale();
};

class logInverseTrans:public logTrans{
public:
	logInverseTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T);
	void transforming(EVENT_DATA_TYPE * input, int nSize);
};

/*
 * Separate class created to faithfully represent GML 2.0 specification
 * for "parametrized logarithmic transformation -- flog" (Section 6.2.1)
 *
 * Names of data members intentionally mirror the specification
 */
class logGML2Trans:public transformation{
public:
		EVENT_DATA_TYPE T; // top of scale value
		EVENT_DATA_TYPE M; // number of logarithmic decades
public:
	logGML2Trans();

	logGML2Trans(EVENT_DATA_TYPE _T,EVENT_DATA_TYPE _M);

	void transforming(EVENT_DATA_TYPE * input, int nSize);
	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	logGML2Trans(const pb::transformation & trans_pb);

	TransPtr  getInverseTransformation();

//	void setTransformedScale(int _scale);
	int getTransformedScale();
	int getRawScale();
};

class logGML2InverseTrans:public logGML2Trans{
public:
	logGML2InverseTrans(EVENT_DATA_TYPE _T,EVENT_DATA_TYPE _M);
	void transforming(EVENT_DATA_TYPE * input, int nSize);
};



class linTrans:public transformation{
public:
	linTrans();
	void transforming(EVENT_DATA_TYPE * input, int nSize);

		TransPtr clone() const;
        void convertToPb(pb::transformation & trans_pb);
        linTrans(const pb::transformation & trans_pb);
        TransPtr getInverseTransformation();
        void setTransformedScale(int scale);
};

/*
 * This class is dedicated to scale the EllipsoidGate
 */
class scaleTrans:public linTrans{
public:
	int t_scale; //transformed scale
	int r_scale; // raw scale
	EVENT_DATA_TYPE scale_factor;

	scaleTrans();
	scaleTrans(int _t_scale, int _r_scale);
	scaleTrans(EVENT_DATA_TYPE _scale_factor);
	void transforming(EVENT_DATA_TYPE * input, int nSize);

	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	scaleTrans(const pb::transformation & trans_pb);

	TransPtr getInverseTransformation();

	void setTransformedScale(int _scale);

};



class flinTrans:public transformation{
	EVENT_DATA_TYPE min;
	EVENT_DATA_TYPE max;
public:
	flinTrans();
	flinTrans(EVENT_DATA_TYPE _minRange, EVENT_DATA_TYPE _maxRange);

	EVENT_DATA_TYPE flin(EVENT_DATA_TYPE x);


	void transforming(EVENT_DATA_TYPE * input, int nSize);

	TransPtr clone() const;

	void convertToPb(pb::transformation & trans_pb);
	flinTrans(const pb::transformation & trans_pb);
	TransPtr getInverseTransformation();
	void setTransformedScale(int scale);

};

/**
 * abridged from flowCore version
 */
struct logicle_params
{
	double T, W, M, A;

	double a, b, c, d, f;
	double w, x0, x1, x2;

	double xTaylor;
	vector<double> taylor;

	int bins;
	bool isInverse;
};

struct sfun_info{
	double b,w;
};
class logicleTrans:public transformation
{
//	const double DEFAULT_DECADES = 4.5;

	const double LN_10 = log(10.);
	const double EPSILON = std::numeric_limits<double>::epsilon();
	const double NaN = std::numeric_limits<double>::quiet_NaN();

	const int TAYLOR_LENGTH = 16;

	logicle_params p;
	bool isGml2;

public:
	logicle_params get_params();
	logicleTrans (double T, double W, double M, double A, bool _isGml2, int bins= 0, bool isInverse = false);
	void init();


	// f(w,b) = 2 * (ln(d) - ln(b)) + w * (b + d)
	static double logicle_fn(double x,void*info);

	/*
	 * root finder routines are copied from stats/src/zeroin.c
	 */
	double R_zeroin(			/* An estimate of the root */
	    double ax,				/* Left border | of the range	*/
	    double bx,				/* Right border| the root is seeked*/
	    double (*f)(double x, void *info),	/* Function under investigation	*/
	    void *info,				/* Add'l info passed on to f	*/
	    double *Tol,			/* Acceptable tolerance		*/
	    int *Maxit)				/* Max # of iterations */
	;

	/* R_zeroin2() is faster for "expensive" f(), in those typical cases where
	 *             f(ax) and f(bx) are available anyway : */

	double R_zeroin2(			/* An estimate of the root */
	    double ax,				/* Left border | of the range	*/
	    double bx,				/* Right border| the root is seeked*/
	    double fa, double fb,		/* f(a), f(b) */
	    double (*f)(double x, void *info),	/* Function under investigation	*/
	    void *info,				/* Add'l info passed on to f	*/
	    double *Tol,			/* Acceptable tolerance		*/
	    int *Maxit)				/* Max # of iterations */
	;
	/*
	 * use R built-in root finder API :R_zeroin
	 */
	double solve (double b, double w);

	double slope (double scale) const;

	double seriesBiexponential (double scale) const;

	double scale (double value) const;
	double inverse (double scale) const;

	virtual void transforming(EVENT_DATA_TYPE * input, int nSize);
	TransPtr clone() const;
	void convertToPb(pb::transformation & trans_pb);
	logicleTrans(const pb::transformation & trans_pb);
	TransPtr  getInverseTransformation();

	void setTransformedScale(int scale);
	int getTransformedScale();
	int getRawScale();
};

};

#endif /* TRANSFORMATION_HPP_ */
