// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/transformation.hpp>
#include <cytolib/global.hpp>

namespace cytolib
{

	transformation::transformation():isGateOnly(false),type(CALTBL),isComputed(true){}
	transformation::transformation(bool _isGate, unsigned short _type):isGateOnly(_isGate),type(_type),isComputed(true){}
	void transformation::transforming(EVENT_DATA_TYPE * input, int nSize){
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

	void transformation::computCalTbl(){};//dummy routine that does nothing
	Spline_Coefs transformation::getSplineCoefs(){return calTbl.getSplineCoefs();};
	void transformation::setCalTbl(calibrationTable _tbl){
		calTbl=_tbl;
	}

	calibrationTable transformation::getCalTbl(){return calTbl;};
	void transformation::interpolate(){calTbl.interpolate();};
	bool transformation::isInterpolated(){return calTbl.isInterpolated();}
	bool transformation::gateOnly(){return isGateOnly;};
	void transformation::setGateOnlyFlag(bool _flag){isGateOnly=_flag;};
	bool transformation::computed(){return isComputed;};
	void transformation::setComputeFlag(bool _flag){isComputed=_flag;};
	string transformation::getName(){return name;};
	void transformation::setName(string _name){name=_name;};
	string transformation::getChannel(){return channel;};
	void transformation::setChannel(string _channel){channel=_channel;};
	unsigned short transformation::getType(){return type;};
	unsigned short transformation::getType(string &ctype){
		switch(type)
		{
			case  CALTBL:
				ctype = "CALTBL";
				break;
			case  LOG:
				ctype = "LOG";
				break;
			case  LIN:
				ctype = "LIN";
				break;
			case  FLIN :
				ctype = "FLIN";
				break;
			case  FASINH :
				ctype = "FASINH";
				break;
			case  BIEXP :
				ctype = "BIEXP";
				break;
			case  LOGICLE :
				ctype = "LOGICLE";
				break;
			case  LOGGML2 :
				ctype = "LOGGML2";
				break;
			case  SCALE :
				ctype = "SCALE";
				break;
			default:
				throw(domain_error("unknown trans type id: " + to_string(type)));
		}

		return type;
	};
	void transformation::setType(unsigned short _type){type=_type;};
	TransPtr transformation::clone() const{return TransPtr(new transformation(*this));};
	transformation::transformation(const pb::transformation & trans_pb){
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
	void transformation::convertToPb(pb::transformation & trans_pb){

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

	TransPtr  transformation::getInverseTransformation(){
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

	void transformation::setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};
	int transformation::getTransformedScale(){throw(domain_error("getTransformedScale function not defined!"));};
	int transformation::getRawScale(){throw(domain_error("getRawScale function not defined!"));};
	EVENT_DATA_TYPE biexpTrans::logRoot(EVENT_DATA_TYPE b, EVENT_DATA_TYPE w)
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


	biexpTrans::biexpTrans():transformation(false, BIEXP),channelRange(4096), pos(4.5), neg(0), widthBasis(-10),maxValue(262144){
		setComputeFlag(false);
		calTbl.setInterpolated(false);
	}

	biexpTrans::biexpTrans(int channelRange_,EVENT_DATA_TYPE pos_, EVENT_DATA_TYPE neg_, EVENT_DATA_TYPE widthBasis_, EVENT_DATA_TYPE maxValue_):transformation(false, BIEXP),channelRange(channelRange_), pos(pos_), neg(neg_), widthBasis(widthBasis_),maxValue(maxValue_){
			setComputeFlag(false);
			calTbl.setInterpolated(false);
		}

	 /*
	  * directly translated from java routine from tree star
	  * potential segfault risk: the inappropriate biexp parameters can cause
	  * the indexing of vector out of the boundary.
	  */

	void biexpTrans::computCalTbl(){
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

	TransPtr biexpTrans::clone() const{return TransPtr(new biexpTrans(*this));};
	void biexpTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_BIEXP);
		pb::biexpTrans * bt_pb = trans_pb.mutable_bt();
		bt_pb->set_channelrange(channelRange);
		bt_pb->set_maxvalue(maxValue);
		bt_pb->set_neg(neg);
		bt_pb->set_pos(pos);
		bt_pb->set_widthbasis(widthBasis);
	}
	biexpTrans::biexpTrans(const pb::transformation & trans_pb):transformation(trans_pb){
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
	void biexpTrans::setTransformedScale(int scale){
		channelRange = scale;
		//recompute cal table
		computCalTbl();
		interpolate();
	};
	int biexpTrans::getTransformedScale(){return channelRange;};
	int biexpTrans::getRawScale(){return maxValue;};

	fasinhTrans::fasinhTrans():transformation(false,FASINH),maxRange(262144),length(256), T(262144),A(0),M(4.5){
		calTbl.setInterpolated(true);
	}

	fasinhTrans::fasinhTrans(EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _T, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M):transformation(false, FASINH),maxRange(_maxRange),length(_length), T(_T),A(_A),M(_M){
		calTbl.setInterpolated(true);
	}
	/*
	 * implementation copied from flowCore
	 */
	void fasinhTrans::transforming(EVENT_DATA_TYPE * input, int nSize){


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


	TransPtr fasinhTrans::clone() const{return TransPtr(new fasinhTrans(*this));};
	void fasinhTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_FASIGNH);
		pb::fasinhTrans * ft_pb = trans_pb.mutable_ft();
		ft_pb->set_a(A);
		ft_pb->set_length(length);
		ft_pb->set_m(M);
		ft_pb->set_maxrange(maxRange);
		ft_pb->set_t(T);
	}
	fasinhTrans::fasinhTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::fasinhTrans & ft_pb = trans_pb.ft();
		length = ft_pb.length();
		maxRange = ft_pb.maxrange();
		T = ft_pb.t();
		A = ft_pb.a();
		M = ft_pb.m();
	}
	void fasinhTrans::setTransformedScale(int scale){maxRange = scale;};
	int fasinhTrans::getTransformedScale(){return maxRange;};
	int fasinhTrans::getRawScale(){return T;};
/*
 * inverse transformation of fasinhTrans
 */
	fsinhTrans::fsinhTrans():fasinhTrans(){}

	fsinhTrans::fsinhTrans(EVENT_DATA_TYPE _maxRange, EVENT_DATA_TYPE _length, EVENT_DATA_TYPE _T
			, EVENT_DATA_TYPE _A, EVENT_DATA_TYPE _M):fasinhTrans(_maxRange, _length, _T, _A, _M){}

	void  fsinhTrans::transforming(EVENT_DATA_TYPE * input, int nSize){
		for(int i=0;i<nSize;i++)
			input[i] = sinh(((M + A) * log(10)) * input[i]/length - A * log(10)) * T / sinh(M * log(10));

	}
	TransPtr fsinhTrans::getInverseTransformation(){throw(domain_error("inverse function not defined!"));};

	TransPtr  fasinhTrans::getInverseTransformation(){
		return TransPtr(new fsinhTrans(maxRange, length, T, A , M));
	}

	logTrans::logTrans():transformation(false,LOG),offset(0),decade(1), scale(1),T(262144){
		calTbl.setInterpolated(true);
	}


	logTrans::logTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T):transformation(false,LOG),offset(_offset),decade(_decade), scale(_scale),T(_T){
		calTbl.setInterpolated(true);
	}


	/*
	 *
	 *now we switch back to zero imputation instead of min value since
	 *when convert to R version of transformation function, the data is
	 *no available anymore, thus no way to specify this minvalue
	 *EDIT:flowjo log formula
	 * scale* (log(x) - log(min))/(log(max) - log(min))
	 * min is recorded as offset
	 * log(max) - log(min) is recorded as decade
	 * scale is 256 when store ellipsoid gate
	 */

	void logTrans::transforming(EVENT_DATA_TYPE * input, int nSize){

			for(int i=0;i<nSize;i++){
				auto & x = input[i];

				x = x>0?((log10(x)-log10(offset))/decade)* scale:0;
			}

	}
	TransPtr logTrans::clone() const{return TransPtr(new logTrans(*this));};
	void logTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LOG);
		pb::logTrans * lt_pb = trans_pb.mutable_lt();
		lt_pb->set_decade(decade);
		lt_pb->set_offset(offset);
		lt_pb->set_t(T);
		lt_pb->set_scale(scale);
	}
	logTrans::logTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::logTrans & lt_pb = trans_pb.lt();
		decade = lt_pb.decade();
		offset = lt_pb.offset();
		T = lt_pb.t();
		scale = lt_pb.scale();
	}

	void logTrans::setTransformedScale(int _scale){scale = _scale;};
	int logTrans::getTransformedScale(){return scale;};
	int logTrans::getRawScale(){return T;};
	logInverseTrans::logInverseTrans(EVENT_DATA_TYPE _offset,EVENT_DATA_TYPE _decade, unsigned _scale, unsigned _T):logTrans(_offset, _decade, _scale, _T){};
	void logInverseTrans::transforming(EVENT_DATA_TYPE * input, int nSize){

			for(int i=0;i<nSize;i++){
				input[i]= pow(10, (input[i]* decade/scale + log10(offset)));
			}

	}


 TransPtr  logTrans::getInverseTransformation(){
		return TransPtr(new logInverseTrans(offset, decade,scale, T));
	}

	logGML2Trans::logGML2Trans():transformation(false,LOGGML2),T(262144),M(1){
		calTbl.setInterpolated(true);
	}


	logGML2Trans::logGML2Trans(EVENT_DATA_TYPE _T,EVENT_DATA_TYPE _M):transformation(false,LOGGML2),T(_T),M(_M){
		calTbl.setInterpolated(true);
	}

	void logGML2Trans::transforming(EVENT_DATA_TYPE * input, int nSize){
		EVENT_DATA_TYPE min = 0.0;
		// Find smallest positive value
		for(int i=0; i< nSize; i++)
		{
			EVENT_DATA_TYPE x = input[i];
			min = (x > 0.0 && ( x < min || !min ))?x:min;
		}
		if(!min){
		  // For nSize == 1, just move up to the lower limit of the transform (0.0)
		  // this allows transformation of negative lower bound of transform range
		  if(nSize > 1)
		    throw(domain_error("All data values are negative. Cannot impute minimum value for GML2 log transform."));
		}
			
		for(int i=0;i<nSize;i++){
			auto & x = input[i];
			// Non GML2-standard imputation logic
			// Bring any negative values up to the smallest
			// positive value
			x = x>0.0?((log10(x)-log10(T))/M)+1:min;
		}
	}

	TransPtr logGML2Trans::clone() const{return TransPtr(new logGML2Trans(*this));};
	void logGML2Trans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LOGGML2);
		pb::logGML2Trans * lgml2t_pb = trans_pb.mutable_lgml2t();
		lgml2t_pb->set_t(T);
		lgml2t_pb->set_m(M);
	}
	logGML2Trans::logGML2Trans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::logGML2Trans & lgml2t_pb = trans_pb.lgml2t();
		T = lgml2t_pb.t();
		M = lgml2t_pb.m();
	}

//	void logTrans::setTransformedScale(int _scale){scale = _scale;};
	int logGML2Trans::getTransformedScale(){return 1;};
	int logGML2Trans::getRawScale(){return T;};
	logGML2InverseTrans::logGML2InverseTrans(EVENT_DATA_TYPE _T,EVENT_DATA_TYPE _M):logGML2Trans(_T, _M){};
	void logGML2InverseTrans::transforming(EVENT_DATA_TYPE * input, int nSize){

			for(int i=0;i<nSize;i++){
				// T*10^(M(x-1))
				input[i]= pow(10, (input[i]*M - M + log10(T)));
			}

	}

	TransPtr logGML2Trans::getInverseTransformation(){
		return TransPtr(new logGML2InverseTrans(T,M));
	}

 linTrans::linTrans():transformation(true,LIN){
	        calTbl.setInterpolated(true);
	}
	void linTrans::transforming(EVENT_DATA_TYPE * input, int nSize){
	for(int i=0;i<nSize;i++)
				input[i]*=64;
	}

	TransPtr linTrans::clone() const{return TransPtr(new linTrans(*this));};
	void linTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LIN);

	}
	linTrans::linTrans(const pb::transformation & trans_pb):transformation(trans_pb){}
	TransPtr linTrans::getInverseTransformation(){throw(domain_error("inverse function not defined!"));};
	void linTrans::setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};



	scaleTrans::scaleTrans():linTrans(),t_scale(256), r_scale(262144){
		isGateOnly = true;
		scale_factor = t_scale/(EVENT_DATA_TYPE)r_scale;
	}
	scaleTrans::scaleTrans(int _t_scale, int _r_scale):linTrans(),t_scale(_t_scale), r_scale(_r_scale){
		isGateOnly = true;
		if((_r_scale == 0) || (_t_scale == 0))
			throw(domain_error("Illegal arguments provided to scaleTrans constructor: t_scale and r_scale must be nonzero"));
		scale_factor = t_scale/(EVENT_DATA_TYPE)r_scale;
	}

	// Defines a scaleTrans solely on float scale factor (not as a ratio of ints)
	scaleTrans::scaleTrans(EVENT_DATA_TYPE _scale_factor):linTrans(),t_scale(0), r_scale(0), scale_factor(_scale_factor){isGateOnly = true;}

	void scaleTrans::transforming(EVENT_DATA_TYPE * input, int nSize){
		for(int i=0;i<nSize;i++)
			input[i]*=scale_factor;
	}

	TransPtr scaleTrans::clone() const{return TransPtr(new scaleTrans(*this));};

	void scaleTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_SCALE);
		// PICKUP and fix this for scaleTrans
		pb::scaleTrans * st_pb = trans_pb.mutable_st();
		st_pb->set_scale_factor(scale_factor);
		st_pb->set_t_scale(t_scale);
		st_pb->set_r_scale(r_scale);
	}

	scaleTrans::scaleTrans(const pb::transformation & trans_pb):linTrans(trans_pb){
		const pb::scaleTrans & st_pb = trans_pb.st();
		scale_factor = st_pb.scale_factor();
		t_scale = st_pb.t_scale();
		r_scale = st_pb.r_scale();
	}

	TransPtr scaleTrans::getInverseTransformation(){
		if(t_scale == 0)
			return TransPtr(new scaleTrans(1.0/scale_factor)); //just invert the scale factor
		else
			return TransPtr(new scaleTrans(r_scale, t_scale));//swap the raw and trans scale
	}

	void scaleTrans::setTransformedScale(int _scale){
		t_scale = _scale;
		scale_factor = t_scale/(EVENT_DATA_TYPE)r_scale;
	};

	flinTrans::flinTrans():transformation(false,FLIN),min(0),max(0){
		calTbl.setInterpolated(true);
	}
	flinTrans::flinTrans(EVENT_DATA_TYPE _minRange, EVENT_DATA_TYPE _maxRange):transformation(false,FLIN),min(_minRange),max(_maxRange){
		calTbl.setInterpolated(true);
	}

	EVENT_DATA_TYPE flinTrans::flin(EVENT_DATA_TYPE x){
		EVENT_DATA_TYPE T=max;
		EVENT_DATA_TYPE A=min;
		return (x+A)/(T+A);
	}


	void flinTrans::transforming(EVENT_DATA_TYPE * input, int nSize){

		for(int i=0;i<nSize;i++){
			input[i]=flin(input[i]);
		}

	}

	TransPtr flinTrans::clone() const{return TransPtr(new flinTrans(*this));};

	void flinTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_FLIN);
		pb::flinTrans * ft_pb = trans_pb.mutable_flt();
		ft_pb->set_max(max);
		ft_pb->set_min(min);
	}
	flinTrans::flinTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::flinTrans & ft_pb = trans_pb.flt();
		max = ft_pb.max();
		min = ft_pb.min();
	}
	TransPtr flinTrans::getInverseTransformation(){throw(domain_error("inverse function not defined!"));};
	void flinTrans::setTransformedScale(int scale){throw(domain_error("setTransformedScale function not defined!"));};

	logicle_params logicleTrans::get_params(){return p;}
	logicleTrans::logicleTrans (double T, double W, double M, double A, bool _isGml2
			, int bins, bool isInverse):transformation(false,LOGICLE)
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
		init();
	}
	void logicleTrans::init(){

		// actual parameters
		// formulas from biexponential paper
		p.w = p.W / (p.M + p.A);
		p.x2 = p.A / (p.M + p.A);
		p.x1 = p.x2 + p.w;
		p.x0 = p.x2 + 2 * p.w;
		p.b = (p.M + p.A) * LN_10;
		p.d = solve(p.b, p.w);
		double c_a = exp(p.x0 * (p.b + p.d));
		double mf_a = exp(p.b * p.x1) - c_a / exp(p.d * p.x1);
		p.a = p.T / ((exp(p.b) - mf_a) - c_a / exp(p.d));
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
	double logicleTrans::logicle_fn(double x,void*info) {
		struct sfun_info *p = (struct sfun_info *)info;
		double B = 2 * (log(x) - log(p->b)) + p->w * (p->b + x);
	    return (B);
	}

	/*
	 * root finder routines are copied from stats/src/zeroin.c
	 */
	double logicleTrans::R_zeroin(			/* An estimate of the root */
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

	double logicleTrans::R_zeroin2(			/* An estimate of the root */
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
	double logicleTrans::solve (double b, double w)
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

	double logicleTrans::slope (double scale) const
	{
		// reflect negative scale regions
		if (scale < p.x1)
			scale = 2 * p.x1 - scale;

		// compute the slope of the biexponential
		return p.a * p.b * exp(p.b * scale) + p.c * p.d / exp(p.d * scale);
	}

	double logicleTrans::seriesBiexponential (double scale) const
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

	double logicleTrans::scale (double value) const
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

	double logicleTrans::inverse (double scale) const
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

	void logicleTrans::transforming(EVENT_DATA_TYPE * input, int nSize){
		float m = isGml2?1:p.M;//set scale to (0,1) for Gml2 version
		if(p.isInverse)
			for(int i=0;i<nSize;i++)
				input[i] = inverse(input[i]/m);
		else
			for(int i=0;i<nSize;i++)
				input[i] = scale(input[i]) * m;

		}


	TransPtr logicleTrans::clone() const{return TransPtr(new logicleTrans(*this));};
	void logicleTrans::convertToPb(pb::transformation & trans_pb){
		transformation::convertToPb(trans_pb);
		trans_pb.set_trans_type(pb::PB_LOGICLE);
		pb::logicleTrans * lt_pb = trans_pb.mutable_lgt();
		lt_pb->set_a(p.A);
		lt_pb->set_w(p.W);
		lt_pb->set_m(p.M);
		lt_pb->set_bins(p.bins);
		lt_pb->set_t(p.T);
		lt_pb->set_isgml2(isGml2);
		lt_pb->set_isinverse(p.isInverse);

	}
	logicleTrans::logicleTrans(const pb::transformation & trans_pb):transformation(trans_pb){
		const pb::logicleTrans & lt_pb = trans_pb.lgt();
		p.W = lt_pb.w();
		p.bins = lt_pb.bins();
		p.T = lt_pb.t();
		p.A = lt_pb.a();
		p.M = lt_pb.m();
		isGml2 = lt_pb.isgml2();
		p.isInverse = lt_pb.isinverse();

		init();
	}
	TransPtr  logicleTrans::getInverseTransformation(){
		logicleTrans tt = *this;
		tt.p.isInverse = !tt.p.isInverse;
		return TransPtr(new logicleTrans(tt));
	}

	void logicleTrans::setTransformedScale(int scale){
		if(isGml2)
			throw(logic_error("can't set scale for logicleGml2!"));
		p.M = scale;};
	int logicleTrans::getTransformedScale(){return isGml2?1:p.M;};
	int logicleTrans::getRawScale(){return p.T;};

};

