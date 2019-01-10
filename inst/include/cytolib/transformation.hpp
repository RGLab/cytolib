/*
 * transformation.hpp
 *
 *  Created on: Apr 24, 2012
 *      Author: wjiang2
 */

#ifndef TRANSFORMATION_HPP_
#define TRANSFORMATION_HPP_
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include "calibrationTable.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/compare.hpp>


#include "global.hpp"

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
/* case insensitive compare predicate*/
struct ciLessBoost : std::binary_function<std::string, std::string, bool>
{
    bool operator() (const std::string & s1, const std::string & s2) const {
        return lexicographical_compare(s1, s2, boost::is_iless());
    }
};

typedef map<std::string, std::string, ciLessBoost> CHANNEL_MAP;

struct coordinate
{
//	friend class boost::serialization::access;

	EVENT_DATA_TYPE x,y;
	coordinate(EVENT_DATA_TYPE _x,EVENT_DATA_TYPE _y){x=_x;y=_y;};
	coordinate(){};
	void convertToPb(pb::coordinate & coor_pb){
		coor_pb.set_x(x);
		coor_pb.set_y(y);
	};
	coordinate(const pb::coordinate & coor_pb):x(coor_pb.x()),y(coor_pb.y()){};
};
inline bool compare_x(coordinate i, coordinate j) { return i.x<j.x; }
inline bool compare_y(coordinate i, coordinate j) { return i.y<j.y; }

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

	transformation():isGateOnly(false),type(CALTBL),isComputed(true){}
	transformation(bool _isGate, unsigned short _type):isGateOnly(_isGate),type(_type),isComputed(true){}
	virtual ~transformation(){};
	virtual void transforming(EVENT_DATA_TYPE * input, int nSize){
		if(!calTbl.isInterpolated()){
			 /* calculate calibration table from the function
			 */
			if(!computed())
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("computing calibration table...\n");
				computCalTbl();
			}

			if(!isInterpolated())
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("spline interpolating...\n");
				interpolate();
			}
		}

		calTbl.transforming(input, nSize);

}

	virtual void computCalTbl(){};//dummy routine that does nothing
	virtual Spline_Coefs getSplineCoefs(){return calTbl.getSplineCoefs();};
	virtual void setCalTbl(calibrationTable _tbl){
		calTbl=_tbl;
	}

	virtual calibrationTable getCalTbl(){return calTbl;};
	virtual void interpolate(){calTbl.interpolate();};
	virtual bool isInterpolated(){return calTbl.isInterpolated();}
	virtual bool gateOnly(){return isGateOnly;};
	virtual void setGateOnlyFlag(bool _flag){isGateOnly=_flag;};
	virtual bool computed(){return isComputed;};
	virtual void setComputeFlag(bool _flag){isComputed=_flag;};
	virtual string getName(){return name;};
	virtual void setName(string _name){name=_name;};
	virtual string getChannel(){return channel;};
	virtual void setChannel(string _channel){channel=_channel;};
	virtual unsigned short getType(){return type;};
	virtual void setType(unsigned short _type){type=_type;};
	virtual TransPtr clone() const{return TransPtr(new transformation(*this));};
	transformation(const pb::transformation & trans_pb){
		isComputed = trans_pb.iscomputed();
		isGateOnly = trans_pb.isgateonly();
		type = trans_pb.type();
		name = trans_pb.name();
		channel = trans_pb.channel();
		/* For PB_BIEXP, caltbl was not saved during archiving and thus could be 0x0
		  and thus trans_pb.caltbl() call is supposed to fall back to return the one from default_instance, which should not be null.
		  However strange enough, in some circumstances (e.g. after save_gs() call), this default_instance_->caltbl_ does become null
		  which leads this pb auto generated accessor function to be unsafe to be invoked.
		  So we have to skip it in case of biexp. (we don't need to do it anyway)
		*/
		if(trans_pb.trans_type() != pb::PB_BIEXP)
			calTbl = calibrationTable(trans_pb.caltbl());
	}
	virtual void convertToPb(pb::transformation & trans_pb){

		trans_pb.set_isgateonly(isGateOnly);
		trans_pb.set_type(type);
		trans_pb.set_name(name);
		trans_pb.set_channel(channel);
		/*skip saving calibration table to save disk space, which means it needs to be always recalculated when load it back
		 	 Setting the flag to FALSE can only make sure it is recomputed properly for the APIs where the flag is checked first
			but it is not sufficient to prevent the segfault because the pointer to the caltbl in pb object is unset and somehow the
			default_instance_->caltbl_ can also be null due to some previous operations (i.e.`save_gs`). So ::pb::calibrationTable& transformation::caltbl
			becomes unsafe to call.
		*/
		if(type == BIEXP)
			trans_pb.set_iscomputed(false);
		else
		{
			trans_pb.set_iscomputed(isComputed);
			pb::calibrationTable * cal_pb = trans_pb.mutable_caltbl();
			calTbl.convertToPb(*cal_pb);
		}



	}

	virtual TransPtr  getInverseTransformation(){
		if(!calTbl.isInterpolated()){
			 /* calculate calibration table from the function
			 */
			if(!computed())
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("computing calibration table...\n");
				computCalTbl();
			}

			if(!isInterpolated())
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("spline interpolating...\n");
				interpolate();
			}
		}

		//clone the existing trans
		TransPtr  inverse = TransPtr(new transformation(*this));
		//make sure to reset type to avoid type-discrepancy because
		//it returns the base transformation type instead of the original one (e.g. biexp)
		inverse->type = CALTBL;

		//swap the x, y vectors in calTbl

		inverse->calTbl.setX(this->calTbl.getY());
		inverse->calTbl.setY(this->calTbl.getX());

		//re-interpolate the inverse calibration tbl
		inverse->calTbl.setInterpolated(false);
		if(g_loglevel>=POPULATION_LEVEL)
				PRINT("spline interpolating...\n");
		inverse->interpolate();
		return inverse;
	}

	virtual void setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};
	virtual int getTransformedScale(){throw(domain_error("getTransformedScale function not defined!"));};
	virtual int getRawScale(){throw(domain_error("getRawScale function not defined!"));};
};
typedef shared_ptr<transformation> TransPtr;
/*
 * directly translated from java routine from tree star
 */
inline EVENT_DATA_TYPE logRoot(EVENT_DATA_TYPE b, EVENT_DATA_TYPE w)
{
	EVENT_DATA_TYPE xLo = 0;
	EVENT_DATA_TYPE xHi = b;
	EVENT_DATA_TYPE d = (xLo + xHi) / 2;
	EVENT_DATA_TYPE dX = abs((long) (xLo - xHi));
	EVENT_DATA_TYPE dXLast = dX;
	EVENT_DATA_TYPE fB = -2 * log(b) + w * b;
	EVENT_DATA_TYPE f = 2. * log(d) + w * b + fB;
	EVENT_DATA_TYPE dF = 2 / d + w;
	if (w == 0) return b;
	for (long i = 0; i < 100; i++)
	{
		if (((d - xHi) * dF - f) * ((d - xLo) * dF - f) >= 0 ||
				abs((long) (2 * f)) > abs((long) (dXLast * dF)))
		{
			dX = (xHi - xLo) / 2;
			d = xLo + dX;
			if (d == xLo)
				return d;
		}
		else
		{
			dX = f / dF;
			EVENT_DATA_TYPE t = d;
			d -= dX;
			if (d == t)
				return d;
		}
		if (abs((long) dX) < 1.0e-12)
			return d;
		dXLast = dX;
		f = 2 * log(d) + w * d + fB;
		dF = 2 / d + w;
		if (f < 0) xLo = d;
		else xHi = d;
	}
	return d;
}


class biexpTrans:public transformation{
public:
	int channelRange;
	EVENT_DATA_TYPE pos, neg, widthBasis, maxValue;
public:

	biexpTrans():transformation(false, BIEXP),channelRange(4096), pos(4.5), neg(0), widthBasis(-10),maxValue(262144){
		setComputeFlag(false);
		calTbl.setInterpolated(false);
	}

	biexpTrans(int channelRange_,EVENT_DATA_TYPE pos_, EVENT_DATA_TYPE neg_, EVENT_DATA_TYPE widthBasis_, EVENT_DATA_TYPE maxValue_):transformation(false, BIEXP),channelRange(channelRange_), pos(pos_), neg(neg_), widthBasis(widthBasis_),maxValue(maxValue_){
			setComputeFlag(false);
			calTbl.setInterpolated(false);
		}

	 /*
	  * directly translated from java routine from tree star
	  * potential segfault risk: the inappropriate biexp parameters can cause
	  * the indexing of vector out of the boundary.
	  */

	void computCalTbl(){
		/*
		 * directly translated from java routine from tree star
		 */

		EVENT_DATA_TYPE ln10 = log(10.0);
		EVENT_DATA_TYPE decades = pos;
		EVENT_DATA_TYPE lowScale = widthBasis;
		EVENT_DATA_TYPE width = log10(-lowScale);

		if (width < 0.5 || width > 3) width = 0.5;
		decades -= width / 2;
		EVENT_DATA_TYPE extra = neg;
		if (extra < 0) extra = 0;
		extra += width / 2;

		int zeroChan = (int)(extra * channelRange / (extra + decades));
		zeroChan = min(zeroChan, channelRange / 2);

		if (zeroChan > 0) decades = extra * channelRange / zeroChan;
		width /= 2 * decades;        // 1.1

		EVENT_DATA_TYPE maximum = maxValue;
		EVENT_DATA_TYPE positiveRange = ln10 * decades;
		EVENT_DATA_TYPE minimum = maximum / exp(positiveRange);
		EVENT_DATA_TYPE negativeRange = logRoot(positiveRange, width);

		EVENT_DATA_TYPE maxChannlVal = channelRange + 1;
		int nPoints = maxChannlVal;//4097;//fix the number of points so that it won't lost the precision when scale is set to 256 (i.e. channelRange = 256)

		vector<EVENT_DATA_TYPE> positive(nPoints), negative(nPoints), vals(nPoints);
		EVENT_DATA_TYPE step = (maxChannlVal-1)/(EVENT_DATA_TYPE)(nPoints -1);
		for (int j = 0; j < nPoints; j++)
		{
			vals[j] = j * step;
			positive[j] = exp((float)(j) / (float)(nPoints) * positiveRange);
			negative[j] = exp((float)(j) / (float)(nPoints) * (-negativeRange));
		}



		EVENT_DATA_TYPE s = exp((positiveRange + negativeRange) * (width + extra / decades));
		for(int j = 0; j < nPoints; j++)
			negative[j] *= s;

		//ensure it is not out-of-bound
		if(zeroChan<0||zeroChan>=nPoints)
			throw(logic_error("invalid zeroChan: " + std::to_string(zeroChan)));

		s = positive[zeroChan] - negative[zeroChan];
		for (int j = zeroChan; j < nPoints; j++)
			positive[j] = minimum * (positive[j] - negative[j] - s);
		for (int j = 0; j < zeroChan; j++){
			int m = 2 * zeroChan - j;
			if(m<0||m>=nPoints)
					throw(logic_error("invalid value from '2 * zeroChan - j': " + std::to_string(m)));
			positive[j] = -positive[m];
		}


		/*
		 * save the calibration table
		 */
		calTbl.setCaltype("flowJo");
		calTbl.setMethod(2);
		calTbl.setX(positive);
		calTbl.setY(vals);

		isComputed=true;


	}

	TransPtr clone() const{return TransPtr(new biexpTrans(*this));};
	void convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_BIEXP);
		pb::biexpTrans * bt_pb = trans_pb.mutable_bt();
		bt_pb->set_channelrange(channelRange);
		bt_pb->set_maxvalue(maxValue);
		bt_pb->set_neg(neg);
		bt_pb->set_pos(pos);
		bt_pb->set_widthbasis(widthBasis);
	}
	biexpTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		if(!trans_pb.has_bt())
			throw(domain_error("biexpTrans field not found in pb::transformation!"));
		const pb::biexpTrans & bt_pb = trans_pb.bt();

		channelRange = bt_pb.channelrange();
		maxValue = bt_pb.maxvalue();
		neg = bt_pb.neg();
		pos = bt_pb.pos();
		widthBasis = bt_pb.widthbasis();

		//make sure to always recompute caltbl (regardless of compute flag) since it was not saved for the sake of space
	//	computCalTbl(); //now we do lazy-compute since some global trans (in the legacy ws and these trans are not actually used by any samples) may not have the valid parameters which will fail the pb unarchiving process
	}
	void setTransformedScale(int scale){
		channelRange = scale;
		//recompute cal table
		computCalTbl();
		interpolate();
	};
	int getTransformedScale(){return channelRange;};
	int getRawScale(){return maxValue;};

};

class fasinhTrans:public transformation{
public:
	EVENT_DATA_TYPE maxRange;
	EVENT_DATA_TYPE length;//unused at this moment
	EVENT_DATA_TYPE T, A, M;
public:
	fasinhTrans():transformation(false,FASINH),maxRange(262144),length(256), T(262144),A(0),M(4.5){
		calTbl.setInterpolated(true);
	}

	fasinhTrans(EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _T, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M):transformation(false, FASINH),maxRange(_maxRange),length(_length), T(_T),A(_A),M(_M){
		calTbl.setInterpolated(true);
	}
	/*
	 * implementation copied from flowCore
	 */
	virtual void transforming(EVENT_DATA_TYPE * input, int nSize){


		for(int i=0;i<nSize;i++){
			input[i] = length * (asinh(input[i] * sinh(M * log(10)) / T) + A * log(10)) / ((M + A) * log(10));
		}
	//		EVENT_DATA_TYPE myB = (M + A) * log(10);
	//		EVENT_DATA_TYPE myC = A * log(10);
	//		EVENT_DATA_TYPE myA = T / sinh(myB - myC);
	//		input = input / myA;
	//
	//		// This formula for the arcsinh loses significance when x is negative
	//		//Therefore we take advantage of the fact that sinh is an odd function
	//		input = abs(input);
	//
	//		input = log(input + sqrt(input * input + 1));
	//		result = rep(NA, times=length(asinhx))
	//		result[negative] = (myC - asinhx[negative]) / myB
	//		result[!negative] = (asinhx[!negative] + myC) / myB
	//		result








	}


	TransPtr clone() const{return TransPtr(new fasinhTrans(*this));};
	void convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_FASIGNH);
		pb::fasinhTrans * ft_pb = trans_pb.mutable_ft();
		ft_pb->set_a(A);
		ft_pb->set_length(length);
		ft_pb->set_m(M);
		ft_pb->set_maxrange(maxRange);
		ft_pb->set_t(T);
	}
	fasinhTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::fasinhTrans & ft_pb = trans_pb.ft();
		length = ft_pb.length();
		maxRange = ft_pb.maxrange();
		T = ft_pb.t();
		A = ft_pb.a();
		M = ft_pb.m();
	}
	TransPtr  getInverseTransformation();

	void setTransformedScale(int scale){maxRange = scale;};
	int getTransformedScale(){return maxRange;};
	int getRawScale(){return T;};
};
/*
 * inverse transformation of fasinhTrans
 */
class fsinhTrans:public fasinhTrans{

public:

	fsinhTrans():fasinhTrans(){}

	fsinhTrans(EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _T, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M):fasinhTrans(_length,_maxRange, _T, _A, _M){}

	void  transforming(EVENT_DATA_TYPE * input, int nSize){
		for(int i=0;i<nSize;i++)
			input[i] = sinh(((M + A) * log(10)) * input[i]/length - A * log(10)) * T / sinh(M * log(10));

	}
	TransPtr getInverseTransformation(){throw(domain_error("inverse function not defined!"));};
//	fsinhTrans * clone(){return new fasinhTrans(*this);};
//	void convertToPb(pb::transformation & trans_pb);
//	fasinhTrans(const pb::transformation & trans_pb);
//	TransPtr getInverseTransformation();
//	void setTransformedScale(int scale){length = scale;};
};

inline TransPtr  fasinhTrans::getInverseTransformation(){
		return TransPtr(new fsinhTrans(length, maxRange, T, A , M));
	}

/*
 * TODO:right now set two flags to TRUE in the contructor to avoid doing cal table stuff,
 * we should consider redesign the classes so that logTrans does not share this extra feature from parent class
 */
class logTrans:public transformation{
public:
		EVENT_DATA_TYPE offset;
		EVENT_DATA_TYPE decade;
		unsigned scale;
		unsigned T; //top value; derived from keyword $PnR for each channel
public:
	logTrans():transformation(false,LOG),offset(0),decade(1), scale(1),T(262144){
		calTbl.setInterpolated(true);
	}


	logTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T):transformation(false,LOG),offset(_offset),decade(_decade), scale(_scale),T(_T){
		calTbl.setInterpolated(true);
	}


	/*
	 *
	 *now we switch back to zero imputation instead of min value since
	 *when convert to R version of transformation function, the data is
	 *no available anymore, thus no way to specify this minvalue
	 *
	 */
	EVENT_DATA_TYPE flog(EVENT_DATA_TYPE x,EVENT_DATA_TYPE T,EVENT_DATA_TYPE _min) {

		EVENT_DATA_TYPE M=decade;
		return x>0?(log10(x/T)/M+offset):_min;
	//	return x>0?(log10((x+offset)/T)/M):_min;

	}

	void transforming(EVENT_DATA_TYPE * input, int nSize){


	//		EVENT_DATA_TYPE thisMax=input.max();//max val must be globally determined during xml parsing
			EVENT_DATA_TYPE thisMin=0;//input.min();

			for(int i=0;i<nSize;i++){
				input[i]=flog(input[i],T,thisMin) * scale;
			}

	}
	TransPtr clone() const{return TransPtr(new logTrans(*this));};
	void convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LOG);
		pb::logTrans * lt_pb = trans_pb.mutable_lt();
		lt_pb->set_decade(decade);
		lt_pb->set_offset(offset);
		lt_pb->set_t(T);
	}
	logTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::logTrans & lt_pb = trans_pb.lt();
		decade = lt_pb.decade();
		offset = lt_pb.offset();
		T = lt_pb.t();
	}

	TransPtr  getInverseTransformation();

	void setTransformedScale(int _scale){scale = _scale;};
	int getTransformedScale(){return scale;};
	int getRawScale(){return T;};
};

class logInverseTrans:public logTrans{
public:
	logInverseTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T):logTrans(_offset, _decade, _scale, _T){};
	void transforming(EVENT_DATA_TYPE * input, int nSize){


	//		EVENT_DATA_TYPE thisMax=input.max();
	//		EVENT_DATA_TYPE thisMin=0;//input.min();

			for(int i=0;i<nSize;i++){
				input[i]= pow(10, (input[i]/scale - 1) * decade) * T;
			}


	//		input=log10(input);

	}

};

inline TransPtr  logTrans::getInverseTransformation(){
		return TransPtr(new logInverseTrans(offset, decade,scale, T));
	}

class linTrans:public transformation{
public:
	linTrans():transformation(true,LIN){
	        calTbl.setInterpolated(true);
	}
	void transforming(EVENT_DATA_TYPE * input, int nSize){
		for(int i=0;i<nSize;i++)
	                input[i]*=64;
	}

		TransPtr clone() const{return TransPtr(new linTrans(*this));};
        void convertToPb(pb::transformation & trans_pb){
        	transformation::convertToPb(trans_pb);
        	trans_pb.set_trans_type(pb::PB_LIN);

        }
        linTrans(const pb::transformation & trans_pb):transformation(trans_pb){}
        TransPtr getInverseTransformation(){throw(domain_error("inverse function not defined!"));};
        void setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};

};

/*
 * This class is dedicated to scale the EllipsoidGate
 */
class scaleTrans:public linTrans{
	int t_scale; //transformed scale
	int r_scale; // raw scale

public:

	scaleTrans():linTrans(),t_scale(256), r_scale(262144){isGateOnly = true;}
	scaleTrans(int _t_scale, int _r_scale):linTrans(),t_scale(_t_scale), r_scale(_r_scale){isGateOnly = true;}

	void transforming(EVENT_DATA_TYPE * input, int nSize){
		for(int i=0;i<nSize;i++)
			input[i]*=(t_scale/(EVENT_DATA_TYPE)r_scale);
	}

	TransPtr clone() const{return TransPtr(new scaleTrans(*this));};

	TransPtr getInverseTransformation(){
		return TransPtr(new scaleTrans(r_scale, t_scale));//swap the raw and trans scale
	}

	void setTransformedScale(int _scale){t_scale = _scale;};

};



class flinTrans:public transformation{
	EVENT_DATA_TYPE min;
	EVENT_DATA_TYPE max;
public:
	flinTrans():transformation(false,FLIN),min(0),max(0){
		calTbl.setInterpolated(true);
	}
	flinTrans(EVENT_DATA_TYPE _minRange, EVENT_DATA_TYPE _maxRange):transformation(false,FLIN),min(_minRange),max(_maxRange){
		calTbl.setInterpolated(true);
	}

	EVENT_DATA_TYPE flin(EVENT_DATA_TYPE x){
		EVENT_DATA_TYPE T=max;
		EVENT_DATA_TYPE A=min;
		return (x+A)/(T+A);
	}


	void transforming(EVENT_DATA_TYPE * input, int nSize){

		for(int i=0;i<nSize;i++){
			input[i]=flin(input[i]);
		}

	}

	TransPtr clone() const{return TransPtr(new flinTrans(*this));};

	void convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_FLIN);
		pb::flinTrans * ft_pb = trans_pb.mutable_flt();
		ft_pb->set_max(max);
		ft_pb->set_min(min);
	}
	flinTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::flinTrans & ft_pb = trans_pb.flt();
		max = ft_pb.max();
		min = ft_pb.min();
	}
	TransPtr getInverseTransformation(){throw(domain_error("inverse function not defined!"));};
	void setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};

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
	const double DEFAULT_DECADES = 4.5;

	const double LN_10 = log(10.);
	const double EPSILON = std::numeric_limits<double>::epsilon();
	const double NaN = std::numeric_limits<double>::quiet_NaN();

	const int TAYLOR_LENGTH = 16;

	logicle_params p;
	bool isGml2;

public:
	logicle_params get_params(){return p;}
	logicleTrans (double T, double W, double M, double A, bool _isGml2, int bins = 0, bool isInverse = false):transformation(false,LOGICLE)
	{
		calTbl.setInterpolated(true);

	  	if (T <= 0)
			throw domain_error("IllegalParameter: T is not positive");
			//throw IllegalParameter("T is not positive");
		if (W < 0)
	        throw domain_error("IllegalParameter: W is not positive");
	        //throw IllegalParameter("W is not positive");
		if (M <= 0)
			throw domain_error("IllegalParameter: M is not positive");
	        //throw IllegalParameter("M is not positive");
		if (2 * W > M)
			throw domain_error("IllegalParameter: W is too large");
	        //throw IllegalParameter("W is too large");
		if (-A > W || A + W > M - W)
	        throw domain_error("IllegalParameter: A is too large");
	        //throw IllegalParameter("A is too large");

		// if we're going to bin the data make sure that
		// zero is on a bin boundary by adjusting A
		if (bins > 0)
		{
			double zero = (W + A) / (M + A);
			zero = floor(zero * bins + .5) / bins;
			A = (M * zero - W) / (1 - zero);
		}
		isGml2 = _isGml2;
		// standard parameters
		p.T = T;
		p.M = M;
		p.W = W;
		p.A = A;
		p.isInverse = isInverse;

		// actual parameters
		// formulas from biexponential paper
		p.w = W / (M + A);
		p.x2 = A / (M + A);
		p.x1 = p.x2 + p.w;
		p.x0 = p.x2 + 2 * p.w;
		p.b = (M + A) * LN_10;
		p.d = solve(p.b, p.w);
		double c_a = exp(p.x0 * (p.b + p.d));
		double mf_a = exp(p.b * p.x1) - c_a / exp(p.d * p.x1);
		p.a = T / ((exp(p.b) - mf_a) - c_a / exp(p.d));
		p.c = c_a * p.a;
		p.f = -mf_a * p.a;

		// use Taylor series near x1, i.e., data zero to
		// avoid round off problems of formal definition
		p.xTaylor = p.x1 + p.w / 4;
		// compute coefficients of the Taylor series
		double posCoef = p.a * exp(p.b * p.x1);
		double negCoef = -p.c / exp(p.d * p.x1);
		// 16 is enough for full precision of typical scales
		p.taylor.resize(TAYLOR_LENGTH);
		for (int i = 0; i < TAYLOR_LENGTH; ++i)
		{
			posCoef *= p.b / (i + 1);
			negCoef *= -p.d / (i + 1);
			(p.taylor)[i] = posCoef + negCoef;
		}
		p.taylor[1] = 0; // exact result of Logicle condition
	}


	// f(w,b) = 2 * (ln(d) - ln(b)) + w * (b + d)
	static double logicle_fn(double x,void*info) {
		struct sfun_info *p = (struct sfun_info *)info;
		double B = 2 * (log(x) - log(p->b)) + p->w * (p->b + x);
	    return (B);
	}

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
	{
	    double fa = (*f)(ax, info);
	    double fb = (*f)(bx, info);
	    return R_zeroin2(ax, bx, fa, fb, f, info, Tol, Maxit);
	}

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
	{
	    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
	    double tol;
	    int maxit;

	    a = ax;  b = bx;
	    c = a;   fc = fa;
	    maxit = *Maxit + 1; tol = * Tol;

	    /* First test if we have found a root at an endpoint */
	    if(fa == 0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return a;
	    }
	    if(fb ==  0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return b;
	    }

	    while(maxit--)		/* Main iteration loop	*/
	    {
		double prev_step = b-a;		/* Distance from the last but one
						   to the last approximation	*/
		double tol_act;			/* Actual tolerance		*/
		double p;			/* Interpolation step is calcu- */
		double q;			/* lated in the form p/q; divi-
						 * sion operations is delayed
						 * until the last moment	*/
		double new_step;		/* Step at this iteration	*/

		if( fabs(fc) < fabs(fb) )
		{				/* Swap data for b to be the	*/
		    a = b;  b = c;  c = a;	/* best approximation		*/
		    fa=fb;  fb=fc;  fc=fa;
		}
		tol_act = 2*EPSILON*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if( fabs(new_step) <= tol_act || fb == (double)0 )
		{
		    *Maxit -= maxit;
		    *Tol = fabs(c-b);
		    return b;			/* Acceptable approx. is found	*/
		}

		/* Decide if the interpolation can be tried	*/
		if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
		    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
						 * Interpolation may be tried	*/
		    double t1,cb,t2;
		    cb = c-b;
		    if( a==c ) {		/* If we have only two distinct	*/
						/* points linear interpolation	*/
			t1 = fb/fa;		/* can only be applied		*/
			p = cb*t1;
			q = 1.0 - t1;
		    }
		    else {			/* Quadric inverse interpolation*/

			q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
			p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
			q = (q-1.0) * (t1-1.0) * (t2-1.0);
		    }
		    if( p>(double)0 )		/* p was calculated with the */
			q = -q;			/* opposite sign; make p positive */
		    else			/* and assign possible minus to	*/
			p = -p;			/* q				*/

		    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
			&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
			new_step = p/q;			/* it is accepted
							 * If p/q is too large then the
							 * bisection procedure can
							 * reduce [b,c] range to more
							 * extent */
		}

		if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
		    if( new_step > (double)0 )	/* than tolerance		*/
			new_step = tol_act;
		    else
			new_step = -tol_act;
		}
		a = b;	fa = fb;			/* Save the previous approx. */
		b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
		    /* Adjust c for it to have a sign opposite to that of b */
		    c = a;  fc = fa;
		}

	    }
	    /* failed! */
	    *Tol = fabs(c-b);
	    *Maxit = -1;
	    return b;
	}

	/*
	 * use R built-in root finder API :R_zeroin
	 */
	double solve (double b, double w)
	{
		// w == 0 means its really arcsinh
		if (w == 0)
			return b;

		// precision is the same as that of b
		double tolerance = 2 * b * EPSILON;
		struct sfun_info params;
		params.b=b;
		params.w=w;

		// bracket the root
		double d_lo = 0;
		double d_hi = b;


		int MaxIt = 20;
		double d ;
		d= R_zeroin(d_lo,d_hi,logicle_fn,(void*)&params,&tolerance,&MaxIt);
		return d;
	}

	double slope (double scale) const
	{
		// reflect negative scale regions
		if (scale < p.x1)
			scale = 2 * p.x1 - scale;

		// compute the slope of the biexponential
		return p.a * p.b * exp(p.b * scale) + p.c * p.d / exp(p.d * scale);
	}

	double seriesBiexponential (double scale) const
	{
		// Taylor series is around x1
		double x = scale - p.x1;
		// note that taylor[1] should be identically zero according
		// to the Logicle condition so skip it here
		double sum = p.taylor[TAYLOR_LENGTH - 1] * x;
		for (int i = TAYLOR_LENGTH - 2; i >= 2; --i)
			sum = (sum + p.taylor[i]) * x;
		return (sum * x + p.taylor[0]) * x;
	}

	double scale (double value) const
	{
		// handle true zero separately
		if (value == 0)
			return p.x1;

		// reflect negative values
		bool negative = value < 0;
		if (negative)
			value = -value;

		// initial guess at solution
		double x;
		if (value < p.f)
			// use linear approximation in the quasi linear region
			x = p.x1 + value / p.taylor[0];
		else
			// otherwise use ordinary logarithm
			x = log(value / p.a) / p.b;

		// try for double precision unless in extended range
		double tolerance = 3 * EPSILON;
		if (x > 1)
			tolerance = 3 * x * EPSILON;

		for (int i = 0; i < 20; ++i)
		{
			// compute the function and its first two derivatives
			double ae2bx = p.a * exp(p.b * x);
			double ce2mdx = p.c / exp(p.d * x);
			double y;
			if (x < p.xTaylor)
				// near zero use the Taylor series
				y = seriesBiexponential(x) - value;
			else
				// this formulation has better roundoff behavior
				y = (ae2bx + p.f) - (ce2mdx + value);
			double abe2bx = p.b * ae2bx;
			double cde2mdx = p.d * ce2mdx;
			double dy = abe2bx + cde2mdx;
			double ddy = p.b * abe2bx - p.d * cde2mdx;

			// this is Halley's method with cubic convergence
			double delta = y / (dy * (1 - y * ddy / (2 * dy * dy)));
			x -= delta;

			// if we've reached the desired precision we're done
			if (std::abs(delta) < tolerance) {
				// handle negative arguments
				if (negative)
					return 2 * p.x1 - x;
				else
					return x;
	        }
		}

	     throw "DidNotConverge: scale() didn't converge";
		//throw DidNotConverge("scale() didn't converge");
	}

	double inverse (double scale) const
	{
		// reflect negative scale regions
		bool negative = scale < p.x1;
		if (negative)
			scale = 2 * p.x1 - scale;

		// compute the biexponential
		double inverse;
		if (scale < p.xTaylor)
			// near x1, i.e., data zero use the series expansion
			inverse = seriesBiexponential(scale);
		else
			// this formulation has better roundoff behavior
			inverse = (p.a * exp(p.b * scale) + p.f) - p.c / exp(p.d * scale);

		// handle scale for negative values
		if (negative)
			return -inverse;
		else
			return inverse;
	}

	virtual void transforming(EVENT_DATA_TYPE * input, int nSize){
		float m = isGml2?1:p.M;//set scale to (0,1) for Gml2 version
		if(p.isInverse)
			for(int i=0;i<nSize;i++)
				input[i] = inverse(input[i]/m);
		else
			for(int i=0;i<nSize;i++)
				input[i] = scale(input[i]) * m;

		}


	TransPtr clone() const{return TransPtr(new logicleTrans(*this));};
	void convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LOGICLE);
		pb::logicleTrans * lt_pb = trans_pb.mutable_lgt();
		lt_pb->set_a(p.A);
		lt_pb->set_w(p.W);
		lt_pb->set_m(p.M);
		lt_pb->set_bins(p.bins);
		lt_pb->set_t(p.T);
		lt_pb->set_isgml2(isGml2);
		//no need to store isInverse flag assuming the inverse won't be used/stored directly by gs
	}
	logicleTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::logicleTrans & lt_pb = trans_pb.lgt();
		p.W = lt_pb.w();
		p.bins = lt_pb.bins();
		p.T = lt_pb.t();
		p.A = lt_pb.a();
		p.M = lt_pb.m();
		isGml2 = lt_pb.isgml2();
	}
	TransPtr  getInverseTransformation(){
		logicleTrans tt = *this;
		tt.p.isInverse = !tt.p.isInverse;
		return TransPtr(new logicleTrans(tt));
	}

	void setTransformedScale(int scale){
		if(isGml2)
			throw(logic_error("can't set scale for logicleGml2!"));
		p.M = scale;};
	int getTransformedScale(){return isGml2?1:p.M;};
	int getRawScale(){return p.T;};
};

};

#endif /* TRANSFORMATION_HPP_ */
