/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * gate.hpp
 *
 *  Created on: Mar 16, 2012
 *      Author: wjiang2
 */

#ifndef GATE_HPP_
#define GATE_HPP_
#include "MemCytoFrame.hpp"
#include "trans_group.hpp"


using namespace std;


namespace cytolib
{
struct BOOL_GATE_OP{
	deque<string> path;
	char op;
	bool isNot;
	void convertToPb(pb::BOOL_GATE_OP & BOOL_GATE_OP_pb);
	BOOL_GATE_OP(){};
	BOOL_GATE_OP(const pb::BOOL_GATE_OP & BOOL_GATE_OP_pb);

} ;

const EVENT_DATA_TYPE pi = 3.1415926535897;



#define POLYGONGATE 1
#define RANGEGATE 2
#define BOOLGATE 3
#define ELLIPSEGATE 4
#define RECTGATE 5
#define LOGICALGATE 6
#define CURLYQUADGATE 7
#define CLUSTERGATE 8
#define QUADGATE 9

#define AND 1
#define OR 2
#define ANDNOT 3
#define ORNOT 4




class vertices_vector{
public:
	vector<EVENT_DATA_TYPE> x;
	vector<EVENT_DATA_TYPE> y;
public:
	void resize(unsigned nSize);
	vertices_vector(){};
	vertices_vector(vector<coordinate> vertices);
	//dummy api for backward compatibility
	vertices_vector toVector() const{return *this;};
	void print();
};


class paramRange
{

private:

	string name;
	EVENT_DATA_TYPE min, max;
public:
	paramRange(EVENT_DATA_TYPE _min,EVENT_DATA_TYPE _max,string _name){min=_min;max=_max;name=_name;};
	paramRange(){};
	vertices_vector toVector() const;
	void setName(string _n){name=_n;};
	void update_channels(const CHANNEL_MAP & chnl_map);
	string getName(){return name;}
	vector<string> getNameArray() const;
	EVENT_DATA_TYPE getMin(){return min;};
	void setMin(EVENT_DATA_TYPE _v){min=_v;};
	EVENT_DATA_TYPE getMax(){return max;};
	void setMax(EVENT_DATA_TYPE _v){max=_v;};
	void convertToPb(pb::paramRange & paramRange_pb){paramRange_pb.set_name(name);paramRange_pb.set_max(max);paramRange_pb.set_min(min);};
	paramRange(const pb::paramRange & paramRange_pb):name(paramRange_pb.name()),min(paramRange_pb.min()),max(paramRange_pb.max()){};
};
class paramPoly
{
private:


	vector<string> params;//params[0] is x, params[1] is y axis
	vector<coordinate> vertices;
public:
	vector<coordinate> getVertices() const{return vertices;};
	void setVertices(vector<coordinate> _v){vertices=_v;};
	vector<string>  getNameArray() const{return params;};
	void setName(vector<string> _params){params=_params;};
	void update_channels(const CHANNEL_MAP & chnl_map);
	vertices_vector toVector() const;
	string xName(){return params[0];};
	string yName(){return params[1];};
	paramPoly(){};
	void convertToPb(pb::paramPoly & paramPoly_pb);
	paramPoly(const pb::paramPoly & paramPoly_pb);
};
class gate;
typedef shared_ptr<gate> gatePtr;

/*
 * TODO:possibly implement getCentroid,getMajorAxis,getMinorAxis for all gate types
 */
/*
 * Important:
 *
 * now that nodePorperties class has customized copy constructor that uses clone member function
 * form gate class. Thus it is necessary to define clone function for each derived gate class
 * in order to avoid the dispatching to parent method and thus degraded to the parent gate object
 */
/**
 * \class gate
 * \brief the base gate class
 *
 * It is an abstract class that is inherited by other concrete gate types.
 */
class gate {
protected:
	bool neg;
	bool isTransformed;
	bool isGained;

public:
	/*
	 * exact string returned by std::type_info::name() is compiler-dependent
	 * so we can't rely on RTTI. instead we return the gate type by API
	 * However it is against the motivation for nodeProperty to use base gate pointer
	 * the very reason of this gate abstraction was to make gatingheirarhcy being agnostic
	 * about the gate type. The reason we are doing it is a compromise to the needs of R API getGate
	 */
	gate():neg(false),isTransformed(false),isGained(false){};
	gate(const pb::gate & gate_pb):neg(gate_pb.neg()),isTransformed(gate_pb.istransformed()),isGained(gate_pb.isgained()){}
	virtual void convertToPb(pb::gate & gate_pb);

	virtual ~gate(){};
	virtual unsigned short getType() const=0;
	virtual vector<BOOL_GATE_OP> getBoolSpec() const{throw(domain_error("undefined getBoolSpec function!"));};
	virtual INDICE_TYPE gating(MemCytoFrame &, INDICE_TYPE &){throw(domain_error("undefined gating function!"));};
	virtual void extend(MemCytoFrame &,float){throw(domain_error("undefined extend function!"));};
	virtual void extend(float,float){throw(domain_error("undefined extend function!"));};
	virtual void gain(map<string,float> &){throw(domain_error("undefined gain function!"));};
	virtual vector<string> getParamNames() const{throw(domain_error("undefined getParam function!"));};
	virtual vertices_vector getVertices() const{throw(domain_error("undefined getVertices function!"));};
	virtual void transforming(trans_local &){throw(domain_error("undefined transforming function!"));};
	virtual void update_channels(const CHANNEL_MAP & chnl_map){throw(domain_error("undefined update_channels function!"));};
	virtual gatePtr clone() const=0;
	virtual bool isNegate() const{return neg;};
	virtual bool gained() const{return isGained;};
	virtual void setNegate(bool _neg){neg=_neg;};
	virtual bool Transformed() const{return isTransformed;};
	virtual void setTransformed(bool _isTransformed){isTransformed=_isTransformed;};
	virtual void setShift(vector<EVENT_DATA_TYPE>) {throw(domain_error("undefined setShift function!"));};
	virtual vector<EVENT_DATA_TYPE> getShift() const{throw(domain_error("undefined getShift function!"));};
	virtual void shiftGate() {throw(domain_error("undefined shiftGate function!"));}
};


class rangeGate:public gate {
private:
	paramRange param;
	vector<EVENT_DATA_TYPE> shift;
public:
	rangeGate():gate(), shift(vector<EVENT_DATA_TYPE>{0.0}){}
	rangeGate(const pb::gate & gate_pb):gate(gate_pb),param(paramRange(gate_pb.rg().param())), shift(vector<EVENT_DATA_TYPE>{0.0}){}
	void convertToPb(pb::gate & gate_pb);
	unsigned short getType() const{return RANGEGATE;}
	void transforming(trans_local & trans);
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd);

	void extend(MemCytoFrame & fdata,float extend_val);
	void extend(float extend_val, float extend_to);
	void gain(map<string,float> & gains);
	paramRange getParam() const{return param;};
	vector<string> getParamNames() const{return param.getNameArray();};
	void setParam(paramRange _param){param=_param;};
	void update_channels(const CHANNEL_MAP & chnl_map){param.update_channels(chnl_map);};
	vertices_vector getVertices() const{return param.toVector();};
	gatePtr clone() const{return gatePtr(new rangeGate(*this));};
	void setShift(vector<EVENT_DATA_TYPE> _shift) {shift=_shift;};
	vector<EVENT_DATA_TYPE> getShift() const{return shift;};
	void shiftGate();
};


/*
 * TODO:using #include <boost/multi_array.hpp> instead to make it easier to convert to R data structure hopefully.
 *
 */
/**
 * \class polygonGate
 * \brief polygon shaped gate
 *
 * It is the most common gate type used in gating.
 */
class polygonGate:public gate {
protected:
	paramPoly param;
	vector<EVENT_DATA_TYPE> shift;
public:
	polygonGate():gate(), shift(vector<EVENT_DATA_TYPE>{0.0, 0.0}){};
	virtual unsigned short getType() const{return POLYGONGATE;}
	/*
	 * when the original gate vertices are at the threshold
	 * it is likely that the gates were truncated in flowJo xml
	 * currently what we can do is to extend it to the real data range to avoid losing
	 * the data points that are below this theshold range
	 * to cut data range)
	 */
	virtual void extend(MemCytoFrame & fdata,float extend_val);

	virtual void extend(float extend_val, float extend_to);
	void gain(map<string,float> & gains);
	 /*
	 *
	 *  reimplement c++ version of inPolygon_c
	 *  indices are allocated within gating function, so it is up to caller to free it
	 *  and now it is freed in destructor of its owner "nodeProperties" object
	 */
	virtual INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd);

	/*
	 * a wrapper that calls transforming(TransPtr , TransPtr )
	 */
	virtual void transforming(trans_local & trans);
	/*
	 * the actual transforming logic for polygonGate, that is shared by polyonGate and ellipsoidGate(due to the special scale)
	 */
	virtual void transforming(TransPtr trans_x, TransPtr trans_y);
	virtual vertices_vector getVertices() const{return param.toVector();};
	void setParam(paramPoly _param){param=_param;};
	void update_channels(const CHANNEL_MAP & chnl_map){param.update_channels(chnl_map);};
	virtual paramPoly getParam() const{return param;};
	virtual vector<string> getParamNames() const{return param.getNameArray();};
	virtual gatePtr clone() const{return gatePtr(new polygonGate(*this));};
	void convertToPb(pb::gate & gate_pb);
	polygonGate(const pb::gate & gate_pb):gate(gate_pb),param(paramPoly(gate_pb.pg().param())),shift(vector<EVENT_DATA_TYPE>{0.0,0.0}){}
	void setShift(vector<EVENT_DATA_TYPE> _shift) {shift=_shift;};
	vector<EVENT_DATA_TYPE> getShift() const{return shift;};
	void shiftGate();
};

enum QUAD{
	Q1 = 1,//-+
	Q2 = 2,//++
	Q3 = 3,//+-
	Q4 = 4//--

};

/*
 * rectgate is a special polygon requires simpler gating routine
 * it doesn't overload getType member function, which means it is exposed to R
 * as a regular polygonGate
 */
/**
 * \class rectGate
 * \brief rectangle gate
 *
 * It is a special polygonGate and has the simpler(faster) gating calculation.
 */
class rectGate:public polygonGate {
	bool is_quad;
	QUAD quadrant;//it is only valid when is_quad is true
public:
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd)
	{
			vector<coordinate> vertices=param.getVertices();
			unsigned nVertex=vertices.size();
			if(nVertex!=2)
				throw(domain_error("invalid number of vertices for rectgate!"));
			string x=param.xName();
			string y=param.yName();
			EVENT_DATA_TYPE * xdata = fdata.get_data_memptr(x, ColType::channel);
			EVENT_DATA_TYPE * ydata =fdata.get_data_memptr(y, ColType::channel);

			int nEvents=parentInd.size();
			INDICE_TYPE res;
			res.reserve(nEvents);

			/*
			 * actual gating
			 */
			for(auto i : parentInd)
			{
				bool inX=false,inY=false;
				EVENT_DATA_TYPE xMin=vertices[0].x;
				EVENT_DATA_TYPE yMin=vertices[0].y;

				EVENT_DATA_TYPE xMax=vertices[1].x;
				EVENT_DATA_TYPE yMax=vertices[1].y;

				if(xMin>xMax||yMin>yMax)
					throw(domain_error("invalid vertices for rectgate!"));

				if(is_quad)
				{
					//avoid the edge cells counted multiple times
					switch(quadrant)
					{
						case Q1:
						{
							inX=xdata[i]<=xMax&&xdata[i]>=xMin;
							inY=ydata[i]<=yMax&&ydata[i]>yMin;
							break;
						}
						case Q2:
						{
							inX=xdata[i]<=xMax&&xdata[i]>xMin;
							inY=ydata[i]<=yMax&&ydata[i]>=yMin;
							break;
						}
						case Q3:
						{
							inX=xdata[i]<=xMax&&xdata[i]>=xMin;
							inY=ydata[i]<yMax&&ydata[i]>=yMin;
							break;
						}
						case Q4:
						{
							inX=xdata[i]<xMax&&xdata[i]>=xMin;
							inY=ydata[i]<=yMax&&ydata[i]>=yMin;
							break;
						}
					}
				}
				else
				{
					inX=xdata[i]<=xMax&&xdata[i]>=xMin;
					inY=ydata[i]<=yMax&&ydata[i]>=yMin;
				}


				bool isIn = inX&&inY;
				if(isIn != neg)
				res.push_back(i);
			}

			return res;

		}
	unsigned short getType() const{return RECTGATE;}
	gatePtr clone() const{return gatePtr(new rectGate(*this));};
	void convertToPb(pb::gate & gate_pb);
	rectGate(const pb::gate & gate_pb):polygonGate(gate_pb), is_quad(false), quadrant(Q1){};;
	rectGate():polygonGate(), is_quad(false), quadrant(Q1){};
	void set_quadrant(QUAD _quadrant)
	{
		is_quad = true;
		quadrant = _quadrant;
	}
};

/**
 * \class ellipseGate
 * \brief ellipse gate
 *
 * It actually no longer needs to inherit polygonGate since we are now doing the gating
 * without interpolating it into polygon. But for backward compatibility (the legacy archive), we preserve this class definition.
 */
class ellipseGate:public polygonGate {
protected:
	vector<coordinate> antipodal_vertices; //four antipodal points of ellipse (to be deprecated)
	coordinate mu;// center point
	vector<coordinate> cov;//covariance matrix
	EVENT_DATA_TYPE dist; //size of ellipse
public:
	ellipseGate(){dist = 1; setShift(vector<EVENT_DATA_TYPE>{0.0,0.0});};
	vector<coordinate> getCovarianceMat() const;
	coordinate getMu() const;
	EVENT_DATA_TYPE getDist() const;
	vector<coordinate> getAntipodalVerts() const;
	virtual unsigned short getType() const{return ELLIPSEGATE;}
	ellipseGate(coordinate _mu, vector<coordinate> _cov, EVENT_DATA_TYPE _dist);

	ellipseGate(vector<coordinate> _antipodal, vector<string> _params);

	void extend(MemCytoFrame & fdata,float extend_val);
	void extend(float extend_val, float extend_to);
	void gain(map<string,float> & gains);
	/*
	 * covert antipodal points to covariance matrix and mean
	 * antipodal points must be transformed first.
	 */
	void computeCov();
	/*
	 * convert covariance matrix and mean to antipodal vertices
	 */
	void computeAntipodalVerts();
	/*
	 * translated from flowCore::%in% method for ellipsoidGate
	 */
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd);
	gatePtr clone() const{return gatePtr(new ellipseGate(*this));};
	void convertToPb(pb::gate & gate_pb);
	ellipseGate(const pb::gate & gate_pb);

	/*
	 * interpolation has to be done on the transformed original 4 coordinates
	 * otherwise, the interpolation results will be wrong
	 */
	void toPolygon(unsigned nVertices);

	/*
	 * Alternative to toPolygon using same logic as flowCore
	 * ellipsoidGate->polygonGate coercion method using Cholesky
	 * decomposition
	 */
	void interpolatePolygon(unsigned nVertices);

	// Need to shift mu and antipodal verts
	// to accomplish shift
	void shiftGate();

	void transforming(trans_local & trans);
};

/*
 * the purpose of having this class is to do the special scaling to the gate coordinates
 * due to the historical FlowJo's implementation (win/vX) of the ellipsoid gate that the foci, distance, and edge points are expressed in 256 x 256 display coordinates
 * to scale back to data space , for linear channel, the scaling factor is max_val/256
 * for non-linear channel, we need to
 * 1. Interpolate it to polygon
 * 2. inverse transform polygon back to raw scale
 * 3. then transform it to data scale
 * Thus we still need to preserve the inheritance to the polygonGate
 */
class ellipsoidGate:public ellipseGate {
public:
	ellipsoidGate():ellipseGate(){};
	ellipsoidGate(vector<coordinate> _antipodal, vector<string> _params):ellipseGate(_antipodal,_params)
	{
		/*
		 * interpolate to polygon gate
		 */
		toPolygon(100);
	}

	gatePtr clone() const{return gatePtr(new ellipsoidGate(*this));};
	void convertToPb(pb::gate & gate_pb);
	ellipsoidGate(const pb::gate & gate_pb);
	/*
	 *
	 * we moved the interpolation to polygonGate form gating method to here because
	 * gating may not be called when only gates to be extracted
	 *
	 *
	 * ellipsoidGate does not follow the regular transforming process
	 * for historical reason, it is defined in 256 * 256 scale.
	 * For linear channel, we simply linear scale it back to raw scale
	 * For non-linear channel, We need to first inverse transform it back to raw scale
	 * before transforming to the ultimate appropriate data scale.
	 */
	void transforming(trans_local & trans);
	/*
	 * ellipsoidGate can't use ellipseGate gating function due to its special treatment of the scale
	 */
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){
		return polygonGate::gating(fdata, parentInd);
	}
	unsigned short getType() const{return POLYGONGATE;}//expose it to R as polygonGate since the original antipodal points can't be used directly anyway

	//ellipsoidGate needs to shift its interpolated points
	void shiftGate();
	/*
	 * overload ellipseGate::extend function
	 * to skip extend logic (to avoid exception) since the gates no longer
	 * seems to be truncated from the latest flowJo wsp
	 * we may phase out the extend logic for other gates eventually
	 *
	 */
	void extend(MemCytoFrame & fdata,float extend_val){
		//do nothing
		};

	void extend(float extend_val, float extend_to){
		//do nothing
		};

};

/*
 *.
 * And gate classes are sits in more abstract level than GatingHierarchy in the C++ class tree,
 * thus GatingHierarchy data structure should be invisible to gate.
 */
/**
 * \class boolGate
 * \brief boolean gate
 *
 * It is not the geometric gate but the boolean combination of the other reference gates.
 * So instead of defining the gating function in this class, the actual gating logic for boolGate is defined
 * in GatingHierarchy::gating function because it needs the indices from the reference nodes which are only accessible at GatingHierarchy object.
 */
class boolGate:public gate {
public:
	boolGate():gate(){};
	vector<BOOL_GATE_OP> boolOpSpec;//the gatePaths with the their logical operators
public:
	vector<BOOL_GATE_OP> getBoolSpec() const{return boolOpSpec;};
	unsigned short getType() const{return BOOLGATE;}
	gatePtr clone() const{return gatePtr(new boolGate(*this));};
	void convertToPb(pb::gate & gate_pb);
	boolGate(const pb::gate & gate_pb);

};
/**
 * \class logicalGate
 * \brief a special boolGate
 * (Now deprecated by the dedicated clusterGate
 * This is mainly used to deal with the situation where the gating algorithm (typically clustering based gating) doesn't generate any type of gate object.
 * In order still be able to record the gating results (i.e. the logical indices), this logicalGate can be used as the dummy gate to be added to the node.
 * Because nodeProperties requires a population node to have a gate to be associated with.
 *
 */
class logicalGate:public boolGate {
private:
	unsigned short getType() const{return LOGICALGATE;}
	gatePtr clone() const{return gatePtr(new logicalGate(*this));};

public:
	void convertToPb(pb::gate & gate_pb);
	logicalGate(const pb::gate & gate_pb):boolGate(gate_pb){};

	logicalGate():boolGate(){};
};

class clusterGate:public boolGate {
private:
	string cluster_method_name_;
	unsigned short getType() const{return CLUSTERGATE; }
	gatePtr clone() const{return gatePtr(new clusterGate(*this));};

public:
	string get_cluster_method_name(){return cluster_method_name_;}
	void convertToPb(pb::gate & gate_pb);
	clusterGate(const pb::gate & gate_pb):boolGate(gate_pb), cluster_method_name_(gate_pb.cg().cluster_method()){};

	clusterGate(string cluster_method_name):boolGate(),cluster_method_name_(cluster_method_name){};
};

/*
 * Before interpolation, the intersection points are stored as the first element of param in polygonGate
 */
class CurlyQuadGate:public polygonGate{
	bool interpolated;
protected:
	QUAD quadrant;
public:
	CurlyQuadGate(paramPoly _inter, QUAD _quad):polygonGate(),interpolated(false),quadrant(_quad){
		param = _inter;
	};
	void transforming(trans_local & trans);



	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd);


	void interpolate(trans_local & trans);
	virtual unsigned short getType() const{return CURLYQUADGATE;}
	gatePtr clone() const{return gatePtr(new CurlyQuadGate(*this));};

};

class quadGate : public polygonGate{
	string uid_;//used to uniquely tag the quadrant gate so that they can be reassembled/grouped
	QUAD quadrant;
public:
	quadGate(paramPoly _intersect, string uid, QUAD _quadrant):polygonGate(), uid_(uid), quadrant(_quadrant)
	{
		param = _intersect;

	};
	coordinate get_intersection() const{return param.getVertices()[0];};
	QUAD get_quadrant() const{return quadrant;};
	void transforming(trans_local & trans){
		polygonGate::transforming(trans);
	}
	rectGate to_rectgate()
	{
		rectGate g;
		paramPoly params;
		vector<coordinate> verts(2);
		auto p = get_intersection();
		switch(quadrant)
		{
			case Q1:
			{
				verts[0] = {-numeric_limits<double>::infinity(), p.y};//bottom left
				verts[1] = {p.x, numeric_limits<double>::infinity()};//top right
				break;
			}
			case Q2:
			{
				verts[0] = {p.x, p.y};//bottom left
				verts[1] = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};//top right
				break;
			}
			case Q3:
			{
				verts[0] = {p.x, -numeric_limits<double>::infinity()};//bottom left
				verts[1] = {numeric_limits<double>::infinity(), p.y};//top right
				break;
			}
			case Q4:
			{
				verts[0] = {-numeric_limits<double>::infinity(), -numeric_limits<double>::infinity()};//bottom left
				verts[1] = {p.x, p.y};//top right
				break;
			}
		}
		params.setName(param.getNameArray());
		params.setVertices(verts);
		g.setParam(params);
		g.set_quadrant(quadrant);
		return g;
	}
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd)
	{
		//construct rect gate on the fly and do the quadrant-specific gating
		return to_rectgate().gating(fdata, parentInd);
	}
	virtual unsigned short getType() const{return QUADGATE;}
	gatePtr clone() const{return gatePtr(new quadGate(*this));};
	string get_uid(){return uid_;}
	void convertToPb(pb::gate & gate_pb){
		polygonGate::convertToPb(gate_pb);
		//cp nested gate
		pb::polygonGate * p_pb = gate_pb.mutable_pg();//already created by previous call, just get its pointer
		pb::quadGate * q_pb = p_pb->mutable_qg();
		gate_pb.set_type(pb::QUAD_GATE);
		q_pb->set_uid(uid_);
		q_pb->set_quadrant(static_cast<pb::QUADRANT>(quadrant));
	}
	quadGate(const pb::gate & gate_pb):polygonGate(gate_pb)
	{

		auto qg = gate_pb.pg().qg();
		uid_ = qg.uid();
		quadrant = static_cast<cytolib::QUAD>(qg.quadrant());

	};

};
};
#endif /* GATE_HPP_ */
