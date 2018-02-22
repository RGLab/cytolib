/*
 * gate.hpp
 *
 *  Created on: Mar 16, 2012
 *      Author: wjiang2
 */

#ifndef GATE_HPP_
#define GATE_HPP_
#include "MemCytoFrame.hpp"
#include "transformation.hpp"
#include "compensation.hpp"
#include "ellipse2points.hpp"


using namespace std;


namespace cytolib
{
struct BOOL_GATE_OP{
	vector<string> path;
	char op;
	bool isNot;
	void convertToPb(pb::BOOL_GATE_OP & BOOL_GATE_OP_pb){
		BOOL_GATE_OP_pb.set_isnot(isNot);
		BOOL_GATE_OP_pb.set_op(op);
		for(unsigned i = 0; i < path.size(); i++){
			 BOOL_GATE_OP_pb.add_path(path[i]);
		}
	};
	BOOL_GATE_OP(){};
	BOOL_GATE_OP(const pb::BOOL_GATE_OP & BOOL_GATE_OP_pb){
		op = BOOL_GATE_OP_pb.op();
		isNot = BOOL_GATE_OP_pb.isnot();
		for(int i = 0; i < BOOL_GATE_OP_pb.path_size(); i++)
			path.push_back(BOOL_GATE_OP_pb.path(i));
	};
	template<class Archive>
				    void serialize(Archive &ar, const unsigned int version)
				    {

						ar & BOOST_SERIALIZATION_NVP(path);
						ar & BOOST_SERIALIZATION_NVP(op);
						ar & BOOST_SERIALIZATION_NVP(isNot);
				    }

} ;

const EVENT_DATA_TYPE pi = 3.1415926535897;



#define POLYGONGATE 1
#define RANGEGATE 2
#define BOOLGATE 3
#define ELLIPSEGATE 4
#define RECTGATE 5
#define LOGICALGATE 6
#define CURLYQUADGATE 7

#define AND 1
#define OR 2
#define ANDNOT 3
#define ORNOT 4

typedef vector<unsigned> INDICE_TYPE;


class vertices_vector{
public:
	vector<EVENT_DATA_TYPE> x;
	vector<EVENT_DATA_TYPE> y;
public:
	void resize(unsigned nSize){
		x.resize(nSize);
		y.resize(nSize);
	}
	vertices_vector(){};
	vertices_vector(vector<coordinate> vertices){

			unsigned nSize=vertices.size();
			resize(nSize);
			for(unsigned i=0;i<nSize;i++)
			{
				x[i]=vertices[i].x;
				y[i]=vertices[i].y;
			}

	};
	//dummy api for backward compatibility
	vertices_vector toVector(){return *this;};
	void print(){
		PRINT("x:");
		for(unsigned i=0;i<x.size();i++)
				PRINT(to_string(x[i])+",");
//		PRINT("x:");
//		for(unsigned i=0;i<x.size();i++)
//				PRINT(x[i]+",");

	}
};


class paramRange
{

private:

	string name;
	EVENT_DATA_TYPE min, max;
public:
	paramRange(EVENT_DATA_TYPE _min,EVENT_DATA_TYPE _max,string _name){min=_min;max=_max;name=_name;};
	paramRange(){};
	vertices_vector toVector(){

		vertices_vector res;
		res.resize(2);
		res.x[0]=min;
		res.x[1]=max;

		return res;
	}
	void setName(string _n){name=_n;};
	void update_channels(const CHANNEL_MAP & chnl_map){

			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(name);
			if(itChnl!=chnl_map.end())
				name = itChnl->second;
	};
	string getName(){return name;}
	vector<string> getNameArray(){
			vector<string> res;
			res.push_back(name);
			return res;
		};
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
	vector<coordinate> getVertices(){return vertices;};
	void setVertices(vector<coordinate> _v){vertices=_v;};
	vector<string>  getNameArray(){return params;};
	void setName(vector<string> _params){params=_params;};
	void update_channels(const CHANNEL_MAP & chnl_map){

			for(vector<string>::iterator it = params.begin(); it != params.end(); it++)
			{
				string curName = *it;

				CHANNEL_MAP::const_iterator itChnl = chnl_map.find(curName);
				if(itChnl!=chnl_map.end())
					*it = itChnl->second;
			}
		};
	vertices_vector toVector(){

		vertices_vector res;
		unsigned nSize=vertices.size();
		res.resize(nSize);
		for(unsigned i=0;i<nSize;i++)
		{
			res.x[i]=vertices[i].x;
			res.y[i]=vertices[i].y;
		}
		return res;
	}

	string xName(){return params[0];};
	string yName(){return params[1];};
	paramPoly(){};
	void convertToPb(pb::paramPoly & paramPoly_pb){
		BOOST_FOREACH(vector<string>::value_type & it, params){
			paramPoly_pb.add_params(it);
		}
		BOOST_FOREACH(vector<coordinate>::value_type & it, vertices){
			pb::coordinate * coor_pb = paramPoly_pb.add_vertices();
			it.convertToPb(*coor_pb);
		}
	};
	paramPoly(const pb::paramPoly & paramPoly_pb){
		for(int i = 0; i < paramPoly_pb.params_size(); i++){
			params.push_back(paramPoly_pb.params(i));
		}
		for(int i = 0; i < paramPoly_pb.vertices_size(); i++){
			vertices.push_back(coordinate(paramPoly_pb.vertices(i)));
		}
	};
};


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
	virtual void convertToPb(pb::gate & gate_pb){
		//cp basic members
			gate_pb.set_istransformed(isTransformed);
			gate_pb.set_neg(neg);
			gate_pb.set_isgained(isGained);
	}

	virtual ~gate(){};
	virtual unsigned short getType()=0;
	virtual vector<BOOL_GATE_OP> getBoolSpec(){throw(domain_error("undefined getBoolSpec function!"));};
	virtual INDICE_TYPE gating(MemCytoFrame &, INDICE_TYPE &){throw(domain_error("undefined gating function!"));};
	virtual void extend(MemCytoFrame &,float){throw(domain_error("undefined extend function!"));};
	virtual void extend(float,float){throw(domain_error("undefined extend function!"));};
	virtual void gain(map<string,float> &){throw(domain_error("undefined gain function!"));};
	virtual vector<string> getParamNames(){throw(domain_error("undefined getParam function!"));};
	virtual vertices_vector getVertices(){throw(domain_error("undefined getVertices function!"));};
	virtual void transforming(trans_local &){throw(domain_error("undefined transforming function!"));};
	virtual void update_channels(const CHANNEL_MAP & chnl_map){throw(domain_error("undefined update_channels function!"));};
	virtual gate * clone()=0;
	virtual bool isNegate(){return neg;};
	virtual bool gained(){return isGained;};
	virtual void setNegate(bool _neg){neg=_neg;};
	virtual bool Transformed(){return isTransformed;};
	virtual void setTransformed(bool _isTransformed){isTransformed=_isTransformed;};
};


class rangeGate:public gate {
private:
	paramRange param;
public:
	rangeGate():gate(){}
	rangeGate(const pb::gate & gate_pb):gate(gate_pb),param(paramRange(gate_pb.rg().param())){}
	void convertToPb(pb::gate & gate_pb){
		gate::convertToPb(gate_pb);
		gate_pb.set_type(pb::RANGE_GATE);
		//cp nested gate
		pb::rangeGate * g_pb = gate_pb.mutable_rg();
		//cp its unique member
		pb::paramRange * pr_pb = g_pb->mutable_param();
		param.convertToPb(*pr_pb);
	}
	unsigned short getType(){return RANGEGATE;}
	void transforming(trans_local & trans){
		if(!Transformed())
		{
			EVENT_DATA_TYPE vert[2] = {param.getMin(),param.getMax()};

			transformation * curTrans=trans.getTran(param.getName());
			if(curTrans!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming "+param.getName()+"\n");

				curTrans->transforming(vert, 2);
				param.setMin(vert[0]);
				param.setMax(vert[1]);
			}
			isTransformed=true;
		}

	}

	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){

		EVENT_DATA_TYPE * data_1d = fdata.get_data_memptr(param.getName(), ColType::channel);

		int nEvents=parentInd.size();
		INDICE_TYPE res;
		res.reserve(nEvents);
		for(auto i : parentInd){
			bool isIn = data_1d[i]<=param.getMax()&&data_1d[i]>=param.getMin();
			if(isIn != neg)
			res.push_back(i);
		}

		return res;
	}

	void extend(MemCytoFrame & fdata,float extend_val){
		string pName=param.getName();
		EVENT_DATA_TYPE * data_1d = fdata.get_data_memptr(pName, ColType::channel);
		int nSize = fdata.n_rows();
		/*
		 * get R_min
		 */

		EVENT_DATA_TYPE xMin= *min_element(data_1d, data_1d + nSize);
		if(param.getMin()<=extend_val)
		{
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("extending "+pName+"from "+to_string(param.getMin())+" to :"+to_string(xMin)+"\n");
			param.setMin(min(xMin, param.getMin()));
		}


	}
	void extend(float extend_val, float extend_to){
		string pName=param.getName();


		EVENT_DATA_TYPE xMin= extend_to;
		if(param.getMin()<=extend_val)
		{
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("extending "+pName+"from "+to_string(param.getMin())+" to :"+to_string(xMin)+"\n");
			param.setMin(min(xMin, param.getMin()));
		}


	}
	void gain(map<string,float> & gains){
		if(!isGained)
		{
			vertices_vector vert(getVertices());

			map<string,float>::iterator it=gains.find(param.getName().c_str());
			if(it!=gains.end())
			{
				float this_gain = it->second;

				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("adjusting "+param.getName()+"\n");

				param.setMin(param.getMin()/this_gain);
				param.setMax(param.getMax()/this_gain);
			}
			isGained=true;
		}
	}
	paramRange getParam(){return param;};
	vector<string> getParamNames(){return param.getNameArray();};
	void setParam(paramRange _param){param=_param;};
	void update_channels(const CHANNEL_MAP & chnl_map){param.update_channels(chnl_map);};
	vertices_vector getVertices(){return param.toVector();};
	rangeGate * clone(){return new rangeGate(*this);};
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
public:
	polygonGate():gate(){};
	virtual unsigned short getType(){return POLYGONGATE;}
	/*
	 * when the original gate vertices are at the threshold
	 * it is likely that the gates were truncated in flowJo xml
	 * currently what we can do is to extend it to the real data range to avoid losing
	 * the data points that are below this theshold range
	 * to cut data range)
	 */
	virtual void extend(MemCytoFrame & fdata,float extend_val){
		string x=param.xName();
		string y=param.yName();
		EVENT_DATA_TYPE* xdata(fdata.get_data_memptr(x, ColType::channel));
		EVENT_DATA_TYPE* ydata(fdata.get_data_memptr(y, ColType::channel));
		int nSize = fdata.n_rows();
		vector<coordinate> v=param.getVertices();
		/*
		 * get R_min
		 */
		EVENT_DATA_TYPE xMin=*min_element(xdata, xdata + nSize);
		EVENT_DATA_TYPE yMin=*min_element(ydata, ydata + nSize);
		for(unsigned i=0;i<v.size();i++)
		{
			if(v[i].x<=extend_val)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("extending " + x + "from " + to_string(v[i].x)+" to :"+to_string(xMin)+"\n");
				v[i].x=min(xMin, v[i].x);
			}
			if(v[i].y<=extend_val)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("extending " + y + "from " + to_string(v[i].y)+" to :"+to_string(yMin)+"\n");
				v[i].y=min(yMin, v[i].y);

			}
		}
		param.setVertices(v);
	}

	virtual void extend(float extend_val, float extend_to){
		string x=param.xName();
		string y=param.yName();

		vector<coordinate> v=param.getVertices();
		/*
		 * get R_min
		 */
		EVENT_DATA_TYPE xMin=extend_to;
		EVENT_DATA_TYPE yMin=extend_to;
		for(unsigned i=0;i<v.size();i++)
		{
			if(v[i].x<=extend_val)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("extending " + x + "from " + to_string(v[i].x)+" to :"+to_string(xMin)+"\n");
				v[i].x=min(xMin,v[i].x);
			}
			if(v[i].y<=extend_val)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("extending " + y + "from " + to_string(v[i].y)+" to :"+to_string(yMin)+"\n");
				v[i].y=min(yMin, v[i].y);

			}
		}
		param.setVertices(v);
	}
	void gain(map<string,float> & gains){

		if(!isGained)
			{
				vector<coordinate> vertices=param.getVertices();
				/*
				 * get channel names to select respective transformation functions
				 */
				string channel_x=param.xName();
				string channel_y=param.yName();



				map<string,float>::iterator it=gains.find(channel_x);
				if(it!=gains.end())
				{
					float this_gain = it->second;
					if(g_loglevel>=POPULATION_LEVEL)
						PRINT("adjusting: "+channel_x+"\n");

					for(unsigned i=0;i<vertices.size();i++)
						vertices[i].x=vertices[i].x/this_gain;
				}

				it=gains.find(channel_y);
				if(it!=gains.end())
				{
					float this_gain = it->second;
					if(g_loglevel>=POPULATION_LEVEL)
						PRINT("adjusting: "+channel_y+"\n");

					for(unsigned i=0;i<vertices.size();i++)
						vertices[i].y=vertices[i].y/this_gain;
				}


				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("\n");
				param.setVertices(vertices);
				isGained=true;
			}



	}
	 /*
	 *
	 *  reimplement c++ version of inPolygon_c
	 *  indices are allocated within gating function, so it is up to caller to free it
	 *  and now it is freed in destructor of its owner "nodeProperties" object
	 */
	virtual INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){




		vector<coordinate> vertices=param.getVertices();
		unsigned nVertex=vertices.size();

		string x=param.xName();
		string y=param.yName();
		EVENT_DATA_TYPE * xdata = fdata.get_data_memptr(x, ColType::channel);
		EVENT_DATA_TYPE * ydata = fdata.get_data_memptr(y, ColType::channel);


		unsigned counter;
		EVENT_DATA_TYPE xinters;
		EVENT_DATA_TYPE p1x, p2x, p1y, p2y, p_y_min, p_y_max;

		int nEvents=parentInd.size();
		INDICE_TYPE res;
		res.reserve(nEvents);
		for(auto i : parentInd)
		{//iterate over points
		p1x=vertices[0].x;
		p1y=vertices[0].y;
		counter=0;
		for(unsigned j=1; j <= nVertex; j++)
		{// iterate over vertices
		  /*p1x,p1y and p2x,p2y are the endpoints of the current vertex*/
		  if (j == nVertex)
		  {//the last vertice must "loop around"
			p2x = vertices[0].x;
			p2y = vertices[0].y;
		  }
		  else
		  {
			p2x = vertices[j].x;
			p2y = vertices[j].y;
		  }
		  /*if horizontal ray is in range of vertex find the x coordinate where
			ray and vertex intersect*/

		  p_y_min = p1y;
		  p_y_max = p2y;
		  if(p1y > p2y)
			  swap(p_y_min, p_y_max);
		  if(ydata[i] >= p_y_min && ydata[i] <= p_y_max &&xdata[i] <= max(p1x, p2x))
		  {
			  xinters = (ydata[i]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x;
			/*if intersection x coordinate == point x coordinate it lies on the
			  boundary of the polygon, which means "in"*/
			if(xinters==xdata[i])
			{
			  counter=1;
			  break;
			}
			/*count how many vertices are passed by the ray*/
			if (xinters > xdata[i])counter++;
		  }
		  p1x=p2x;
		  p1y=p2y;
		}
		/*uneven number of vertices passed means "in"*/

		 bool isIn =((counter % 2) != 0);
		 if(isIn != neg)
			res.push_back(i);
		}

		return res;
	}


	/*
	 * a wrapper that calls transforming(transformation * , transformation * )
	 */
	virtual void transforming(trans_local & trans){

			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();


			/*
			 * do the actual transformations
			 */
			transformation * trans_x=trans.getTran(channel_x);
			transformation * trans_y=trans.getTran(channel_y);

			transforming(trans_x, trans_y);
	}

	/*
	 * the actual transforming logic for polygonGate, that is shared by polyonGate and ellipsoidGate(due to the special scale)
	 */
	virtual void transforming(transformation * trans_x, transformation * trans_y){
		if(!Transformed())
		{
			vector<coordinate> vertices=param.getVertices();
			int nSize = vertices.size();
			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();

			//get vertices in array format
			vertices_vector vert(vertices);


			if(trans_x!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_x+"\n");;
				trans_x->transforming(&vert.x[0], nSize);
				for(int i=0;i<nSize;i++)
					vertices[i].x=vert.x[i];
			}
			if(trans_y!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_y+"\n");;
				trans_y->transforming(&vert.y[0], nSize);
				for(int i=0;i<nSize;i++)
					vertices[i].y=vert.y[i];
			}
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("\n");
			param.setVertices(vertices);
			isTransformed=true;
		}
	}
	virtual vertices_vector getVertices(){return param.toVector();};
	void setParam(paramPoly _param){param=_param;};
	void update_channels(const CHANNEL_MAP & chnl_map){param.update_channels(chnl_map);};
	virtual paramPoly getParam(){return param;};
	virtual vector<string> getParamNames(){return param.getNameArray();};
	virtual polygonGate * clone(){return new polygonGate(*this);};
	void convertToPb(pb::gate & gate_pb){
		gate::convertToPb(gate_pb);

		gate_pb.set_type(pb::POLYGON_GATE);
		//cp nested gate
		pb::polygonGate * g_pb = gate_pb.mutable_pg();
		//cp its unique member
		pb::paramPoly * pr_pb = g_pb->mutable_param();
		param.convertToPb(*pr_pb);
	}
	polygonGate(const pb::gate & gate_pb):gate(gate_pb),param(paramPoly(gate_pb.pg().param())){}

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
public:



	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){

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
			bool inX,inY;
			EVENT_DATA_TYPE xMin=vertices[0].x;
			EVENT_DATA_TYPE yMin=vertices[0].y;

			EVENT_DATA_TYPE xMax=vertices[1].x;
			EVENT_DATA_TYPE yMax=vertices[1].y;

			if(xMin>xMax||yMin>yMax)
				throw(domain_error("invalid vertices for rectgate!"));

			inX=xdata[i]<=xMax&&xdata[i]>=xMin;
			inY=ydata[i]<=yMax&&ydata[i]>=yMin;
			bool isIn = inX&&inY;
			if(isIn != neg)
			res.push_back(i);
		}

		return res;

	}

	unsigned short getType(){return RECTGATE;}
	rectGate * clone(){return new rectGate(*this);};
	void convertToPb(pb::gate & gate_pb)
	{
		polygonGate::convertToPb(gate_pb);
			gate_pb.set_type(pb::RECT_GATE);
	}
	rectGate(const pb::gate & gate_pb):polygonGate(gate_pb){};;
	rectGate():polygonGate(){};
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
	ellipseGate(){dist = 1;};
	vector<coordinate> getCovarianceMat(){
		if(!Transformed())
			throw(domain_error("EllipseGate has not been transformed so covariance matrix is unavailable!"));
		return cov;};
	coordinate getMu(){
		if(!Transformed())
				throw(domain_error("EllipseGate has not been transformed so mu is unavailable!"));
		return mu;};
	EVENT_DATA_TYPE getDist(){
		if(!Transformed())
			throw(domain_error("EllipseGate has not been transformed so dist is unavailable!"));
		return dist;};
	virtual unsigned short getType(){return ELLIPSEGATE;}
	ellipseGate(coordinate _mu, vector<coordinate> _cov, EVENT_DATA_TYPE _dist):mu(_mu),cov(_cov), dist(_dist){
		isTransformed = true;
		isGained = true;
		neg = false;
	}

	ellipseGate(vector<coordinate> _antipodal, vector<string> _params):antipodal_vertices(_antipodal),dist(1){
		isTransformed = false;
		isGained = false;
		neg = false;

		/*
		 * init the dummy vertices for base class
		 * (this deprecated inheritance exists for the sake of legacy archive)
		 */
		param.setName(_params);

	}

	void extend(MemCytoFrame & fdata,float extend_val){

		/*
		 * get R_min
		 */
		vector<coordinate> v=param.getVertices();
		for(unsigned i=0;i<v.size();i++)
		{
			if((v[i].x<=extend_val)|(v[i].y<=extend_val))
			{
				throw(domain_error("try to extend the coordinates for ellipse gate!"));
			}

		}

	}
	void extend(float extend_val, float extend_to){

		/*
		 * get R_min
		 */
		vector<coordinate> v=param.getVertices();
		for(unsigned i=0;i<v.size();i++)
		{
			if((v[i].x<=extend_val)|(v[i].y<=extend_val))
			{
				throw(domain_error("try to extend the coordinates for ellipse gate!"));
			}

		}

	}
	void gain(map<string,float> & gains){
		if(!isGained)
		{
			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();


			map<string,float>::iterator it=gains.find(channel_x);
			if(it!=gains.end())
			{
				float this_gain = it->second;
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("adjusting: "+channel_x+"\n");;
				for(unsigned i=0;i<antipodal_vertices.size();i++)
					antipodal_vertices[i].x=antipodal_vertices[i].x/this_gain;
			}
			it=gains.find(channel_y);
			if(it!=gains.end())
			{
				float this_gain = it->second;
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("adjusting: "+channel_y+"\n");;
				for(unsigned i=0;i<antipodal_vertices.size();i++)
					antipodal_vertices[i].y=antipodal_vertices[i].y/this_gain;
			}
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("\n");

			isGained=true;
		}
	}

	/*
	 * covert antipodal points to covariance matrix and mean
	 * antipodal points must be transformed first.
	 */
	void computeCov(){
		if(!Transformed())
			throw(domain_error("antipodal points of ellipseGate must be transformed before computing covariance matrix!"));

		vector<coordinate> v=antipodal_vertices;
		unsigned short nSize = v.size();
		if (nSize != 4)
			throw(domain_error("invalid number of antipodal points"));

		/*
		 * get center and set mu
		 */
		mu.x=0;
		mu.y=0;
		for(vector<coordinate>::iterator it=v.begin();it!=v.end();it++)
		{
			mu.x+=it->x;
			mu.y+=it->y;
		}
		mu.x=mu.x/nSize;
		mu.y=mu.y/nSize;

		//center the antipods
		for(vector<coordinate>::iterator it=v.begin();it!=v.end();it++)
		{
			it->x = it->x - mu.x;
			it->y = it->y - mu.y;
		}

		/*
		 * find the four positions of four antipodals
		 */

		//far right point
		vector<coordinate>::iterator R_it=max_element(v.begin(),v.end(),compare_x);
		coordinate R = *R_it;

		//far left point
		vector<coordinate>::iterator L_it=min_element(v.begin(),v.end(),compare_x);
		coordinate L = *L_it;

		// calculate the a length
		EVENT_DATA_TYPE a = hypot(L.x-R.x,L.y-R.y)/2;

		//use the rest of two points for computing b
		vector<coordinate> Q;
		for(vector<coordinate>::iterator it = v.begin();it!= v.end();it++){
			if(it != R_it && it != L_it)
				Q.push_back(*it);
		}
		coordinate V1 = Q[0];
		coordinate V2 = Q[1];
		EVENT_DATA_TYPE b = hypot(V1.x-V2.x,V1.y-V2.y)/2;

		EVENT_DATA_TYPE a2 = a * a ;
		EVENT_DATA_TYPE b2 = b * b ;


		//normailize R and V1 first
		EVENT_DATA_TYPE L_norm = hypot(L.x, L.y);
		EVENT_DATA_TYPE x1 = L.x/L_norm;
		EVENT_DATA_TYPE y1 = L.y/L_norm;

		EVENT_DATA_TYPE V1_norm = hypot(V1.x, V1.y);
		EVENT_DATA_TYPE x2 = V1.x/V1_norm;
		EVENT_DATA_TYPE y2 = V1.y/V1_norm;

		coordinate p1;
		p1.x = x1 * x1 * a2 + x2 * x2 * b2;
		p1.y = x1 * y1 * a2 + x2 * y2 * b2;

		coordinate p2;
		p2.x = p1.y;
		p2.y = y1 * y1 * a2 + y2 * y2 * b2;


		//set cov
		cov.push_back(p1);
		cov.push_back(p2);

		//set distance (in this calculation should always be 1)
		dist = 1;
	}

	/*
	 * translated from flowCore::%in% method for ellipsoidGate
	 */
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){


		// get data

		EVENT_DATA_TYPE * xdata = fdata.get_data_memptr(param.xName(), ColType::channel);
		EVENT_DATA_TYPE * ydata = fdata.get_data_memptr(param.yName(), ColType::channel);


		//inverse the cov matrix
		/*
		 * 	| a,b |
			| c,d | --> | aa, bb |
						| cc, dd |
		 */
		EVENT_DATA_TYPE a , b, c, d;
		if(cov.size()!=2)
			throw(domain_error("invalid cov matrix!"));
		a = cov[0].x;
		b = cov[0].y;
		c = cov[1].x;
		d = cov[1].y;

		EVENT_DATA_TYPE det = a* d - b* c;
		EVENT_DATA_TYPE aa, bb, cc, dd;
		aa = d/det;
		bb = -b/det;
		cc = -c/det;
		dd = a/det;

		// if inside of the ellipse
		int nEvents=parentInd.size();
		INDICE_TYPE res;
		res.reserve(nEvents);
		for(auto i : parentInd){
			//center the data

			EVENT_DATA_TYPE x = xdata[i] - mu.x;
			EVENT_DATA_TYPE y = ydata[i] - mu.y;
			bool isIn = (x * x * aa + x* y * cc + x* y * bb + y * y * dd) <= pow(dist, 2);
			if(isIn != neg)
			res.push_back(i);
		}

		return res;
	}
	ellipseGate * clone(){return new ellipseGate(*this);};
	void convertToPb(pb::gate & gate_pb)
	{
		polygonGate::convertToPb(gate_pb);

			gate_pb.set_type(pb::ELLIPSE_GATE);
			//cp nested gate
			pb::ellipseGate * g_pb = gate_pb.mutable_eg();
			//cp its unique member
			g_pb->set_dist(dist);
			pb::coordinate * coor_pb = g_pb->mutable_mu();
			mu.convertToPb(*coor_pb);
			for(unsigned i = 0; i < cov.size(); i++){
				pb::coordinate * coor_pb = g_pb->add_cov();
				cov[i].convertToPb(*coor_pb);
			}
			for(unsigned i = 0; i < antipodal_vertices.size(); i++){
				pb::coordinate * coor_pb = g_pb->add_antipodal_vertices();
				antipodal_vertices[i].convertToPb(*coor_pb);
			}
	}
	ellipseGate(const pb::gate & gate_pb):polygonGate(gate_pb),mu(coordinate(gate_pb.eg().mu())),dist(gate_pb.eg().dist()){
		const pb::ellipseGate & eg_pb = gate_pb.eg();
		for(int i = 0; i < eg_pb.antipodal_vertices_size(); i++){
			antipodal_vertices.push_back(coordinate(eg_pb.antipodal_vertices(i)));
		}
		for(int i = 0; i < eg_pb.cov_size(); i++){
			cov.push_back(coordinate(eg_pb.cov(i)));
		}
	}

	/*
	 * interpolation has to be done on the transformed original 4 coordinates
	 * otherwise, the interpolation results will be wrong
	 */
	void toPolygon(unsigned nVertices){




		/*
		 * using 4 vertices to fit polygon points
		 */
		vector<coordinate> v=antipodal_vertices;
		vector<coordinate> vertices=param.getVertices();
		vertices.clear();//reset the vertices
		vertices.resize(nVertices);
		/*
		 * fit the polygon points
		 */

		vector<float> x, y;
		for(auto & i : antipodal_vertices)
		{
			x.push_back(i.x);
			y.push_back(i.y);
		}

		ellipse_parsed res = parseEllipse(x, y);
		matrix mat = toPoly(res,nVertices);
		for(unsigned short i=0;i<nVertices;i++)
		{
			vertices[i].x = mat.x[i];
			vertices[i].y = mat.y[i];
		}

		param.setVertices(vertices);

	}

	void transforming(trans_local & trans){
		if(!Transformed())
		{
			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();

			//get vertices in valarray format
			vertices_vector vert(antipodal_vertices);
			int nSize = antipodal_vertices.size();
			/*
			 * do the actual transformations
			 */
			transformation * trans_x=trans.getTran(channel_x);
			transformation * trans_y=trans.getTran(channel_y);


			if(trans_x!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_x+"\n");;

				trans_x->transforming(&vert.x[0],nSize);
				for(int i=0;i<nSize;i++)
					antipodal_vertices[i].x=vert.x[i];
			}
			if(trans_y!=NULL)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_y+"\n");;

				trans_y->transforming(&vert.y[0],nSize);
				for(int i=0;i<nSize;i++)
					antipodal_vertices[i].y=vert.y[i];
			}
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("\n");
			isTransformed=true;

			//compute the covariance matrix after transformed
			computeCov();

		}
	}

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

	ellipsoidGate * clone(){return new ellipsoidGate(*this);};
	void convertToPb(pb::gate & gate_pb)
	{
		ellipseGate::convertToPb(gate_pb);
			gate_pb.set_type(pb::ELLIPSOID_GATE);
	}
	ellipsoidGate(const pb::gate & gate_pb):ellipseGate(gate_pb){
		//deal with legacy archive that did not interpolate ellipsoidGate
		if(param.getVertices().size() == 0)
			toPolygon(100);
	}
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
	void transforming(trans_local & trans){

		if(!Transformed())
		{
			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();



			transformation * trans_x=trans.getTran(channel_x);
			transformation * trans_y=trans.getTran(channel_y);


			/*
			 * re-construct the trans object that was used by flowJo to transform ellipsoid gate to 256 scale
			 */
			unique_ptr<transformation> trans_gate_x,trans_gate_y;
			if(trans_x == NULL)
				throw(domain_error("ellipsoidGate::transforming can't find transformation for " + channel_x));
	//			trans_gate_x.reset(new scaleTrans()); //create default scale trans for linear, assuming the max value for linear scale is always 262144
	//		else
				trans_gate_x.reset(trans_x->clone()); //copy existing trans_x for non-linear
	//
			if(trans_y == NULL)
				throw(domain_error("ellipsoidGate::transforming can't find transformation for " + channel_y));
	//			trans_gate_y.reset(new scaleTrans()); //create default scale trans for linear
	//		else
				trans_gate_y.reset(trans_y->clone()); //copy existing trans_y for non-linear

			//set to scale 256
			trans_gate_x->setTransformedScale(256);
			trans_gate_y->setTransformedScale(256);

			//get its inverse
			boost::shared_ptr<transformation> inverseTrans_x = trans_gate_x->getInverseTransformation();
			boost::shared_ptr<transformation> inverseTrans_y = trans_gate_y->getInverseTransformation();


			/*
			 * transform the polygon from 256 to raw
			 */
			polygonGate::transforming(inverseTrans_x.get(), inverseTrans_y.get());



			/*
			 * transform the raw to the actual data scale (for non-linear channel)
			 */
			isTransformed = false;//reset transform flag otherwise the transforming won't get executed
			polygonGate::transforming(trans_x, trans_y);

			isTransformed=true;
		}

	}
	/*
	 * ellipsoidGate can't use ellipseGate gating function due to its special treatment of the scale
	 */
	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){
		return polygonGate::gating(fdata, parentInd);
	}
	unsigned short getType(){return POLYGONGATE;}//expose it to R as polygonGate since the original antipodal points can't be used directly anyway
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
	vector<BOOL_GATE_OP> getBoolSpec(){return boolOpSpec;};
	unsigned short getType(){return BOOLGATE;}
	boolGate * clone(){return new boolGate(*this);};
	void convertToPb(pb::gate & gate_pb){
		gate::convertToPb(gate_pb);

		gate_pb.set_type(pb::BOOL_GATE);
		//cp nested gate
		pb::boolGate * g_pb = gate_pb.mutable_bg();
		//cp its unique member
		for(unsigned i = 0; i < boolOpSpec.size(); i++){
			pb::BOOL_GATE_OP * gop_pb = g_pb->add_boolopspec();
			boolOpSpec[i].convertToPb(*gop_pb);
		}

	}
	boolGate(const pb::gate & gate_pb):gate(gate_pb){
		const pb::boolGate & bg_pb = gate_pb.bg();
		for(int i = 0; i < bg_pb.boolopspec_size(); i++){
			const pb::BOOL_GATE_OP & thisOP_pb = bg_pb.boolopspec(i);
			BOOL_GATE_OP thisOP = BOOL_GATE_OP(thisOP_pb);
			boolOpSpec.push_back(thisOP);


		}
	}

};
/**
 * \class logicalGate
 * \brief a special boolGate
 *
 * This is mainly used to deal with the situation where the gating algorithm (typically clustering based gating) doesn't generate any type of gate object.
 * In order still be able to record the gating results (i.e. the logical indices), this logicalGate can be used as the dummy gate to be added to the node.
 * Because nodeProperties requires a population node to have a gate to be associated with.
 *
 */
class logicalGate:public boolGate {
private:
	unsigned short getType(){return LOGICALGATE;}
	logicalGate * clone(){return new logicalGate(*this);};

public:
	void convertToPb(pb::gate & gate_pb){
		boolGate::convertToPb(gate_pb);
		gate_pb.set_type(pb::LOGICAL_GATE);
	}
	logicalGate(const pb::gate & gate_pb):boolGate(gate_pb){};

	logicalGate():boolGate(){};
};

enum QUAD{
	Q1,//-+
	Q2,//++
	Q3,//+-
	Q4//--

};
/*
 * Before interpolation, the intersection points are stored as the first element of param in polygonGate
 */
class CurlyGuadGate:public polygonGate{
	bool interpolated;
	QUAD quadrant;
public:
	CurlyGuadGate(paramPoly _inter, QUAD _quad):polygonGate(),interpolated(false),quadrant(_quad){
		param = _inter;
	};
	void transforming(trans_local & trans){
		if(interpolated)
			polygonGate::transforming(trans);
		else
			throw(logic_error("CurlyGuadGate can't not be transformed before interpolation!"));
	};



	INDICE_TYPE gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){
		if(interpolated)
		{
			return polygonGate::gating(fdata, parentInd);
		}
		else
		{
			throw(logic_error("CurlyQuad gate has not been converted to polygonGate yet!"));
		}


	}


	void interpolate(trans_local & trans){

		string x_chnl = param.xName();
		string y_chnl = param.yName();
		/*
		 * transform intersect back to raw
		 */

		transformation * trans_x = trans.getTran(x_chnl);
		transformation * trans_y = trans.getTran(y_chnl);


		/*
		 * and rescale raw to 256 space
		 */
		unique_ptr<transformation> trans_gate_x,trans_gate_y;
		if(trans_x == NULL)
			trans_gate_x.reset(new scaleTrans()); //create default scale trans for linear, assuming the max value for linear scale is always 262144
		else
			trans_gate_x.reset(trans_x->clone()); //copy existing trans_x for non-linear

		if(trans_y == NULL)
			trans_gate_y.reset(new scaleTrans()); //create default scale trans for linear
		else
			trans_gate_y.reset(trans_y->clone()); //copy existing trans_y for non-linear

		//set to scale 256
		int displayScale = 255;
		trans_gate_x->setTransformedScale(displayScale);
		trans_gate_y->setTransformedScale(displayScale);
		polygonGate::transforming(trans_gate_x.get(), trans_gate_y.get());

	//	/*
	//	 * directly map from log scale to 225 space to make the curve smoother
	//	 */
	//	int displayScale = 255;
	//	scaleTrans tx(displayScale, trans_x->getRawScale());
	//	scaleTrans ty(displayScale, trans_y->getRawScale());
	//	scaleTrans *trans_gate_x = &tx;
	//	scaleTrans *trans_gate_y = &ty;
	//	polygonGate::transforming(trans_gate_x, trans_gate_y);



		setTransformed(false);//reset flag so that it won't interfere the next transforming

		coordinate center = param.getVertices()[0];
		EVENT_DATA_TYPE x_mu = center.x;
		EVENT_DATA_TYPE y_mu = center.y;
		//locate the a value
		EVENT_DATA_TYPE multiplier = 0.001;


		/*
		 * interpolate two curves
		 */
		int nLen = 40;
		vector<coordinate> curve1(nLen), curve2(nLen);
		//curve1: round(multiplier * (x - x.mu) ^ 2) + y.mu (horizontal)
		EVENT_DATA_TYPE x_max = displayScale;//xdata.max();
		EVENT_DATA_TYPE y_max = displayScale;//ydata.max();
		EVENT_DATA_TYPE nStep = (x_max - x_mu) / nLen;
		EVENT_DATA_TYPE delta;
		for(auto i = 0; i < nLen; i++){
			delta = nStep * i;
			curve1[i].x = x_mu + delta;
			curve1[i].y = multiplier * pow(delta, 2) + y_mu;
		}
		//curve2:  (vertical)
		nStep = (y_max - y_mu) / nLen;
		for(auto i = 0; i < nLen; i++){
			delta = nStep * i;
			curve2[i].y = y_mu + delta;
			curve2[i].x = multiplier * pow(delta, 2) + x_mu;
		}

		vector<coordinate> polyVert; //the interpolated vertices for polygon
		EVENT_DATA_TYPE x_min = -4e3;//-numeric_limits<EVENT_DATA_TYPE>::max();//xdata.min();
		EVENT_DATA_TYPE y_min = -4e3;//-numeric_limits<EVENT_DATA_TYPE>::max();//ydata.min();


		/*
		 * add the other edges
		 */
		switch(quadrant)
		{
		case Q1:
		{
			//start with curv2
			polyVert = curve2;
			//top left
			polyVert.push_back(coordinate(x_min, y_max));
			//bottom left
			polyVert.push_back(coordinate(x_min, y_mu));
			//bottom right
			polyVert.push_back(curve2.front());
		}
			break;
		case Q2:
		{
			//start with curv1
			polyVert = curve1;
			//top right
			polyVert.push_back(coordinate(x_max, y_max));
			//top left
			polyVert.push_back(curve2.back());
			//add curve2 reversely
			unsigned len = polyVert.size();
			polyVert.resize(len+curve2.size());
			reverse_copy(curve2.begin(), curve2.end(), polyVert.begin()+len);
		}
			break;
		case Q3:
		{
			polyVert = curve1;
			//bottom right

			polyVert.push_back(coordinate(x_max,y_min));
			//bottom left
			polyVert.push_back(coordinate(x_mu,y_min));
			//top left
			polyVert.push_back(center);
		}
			break;
		case Q4://quadrant 4 is actually a rectangle
		{
			polyVert.push_back(center);

			polyVert.push_back(coordinate(x_mu, y_min));
			polyVert.push_back(coordinate(x_min, y_min));
			polyVert.push_back(coordinate(x_min, y_mu));
			polyVert.push_back(center);
		}
			break;
		default:
			throw(logic_error("invalid quadrant"));
		}

		param.setVertices(polyVert);

		/*
		 * scale back to the raw scale
		 */
		boost::shared_ptr<transformation> inverseGate_x,inverseGate_y;
		if(trans_gate_x){
			inverseGate_x = trans_gate_x->getInverseTransformation();
		}
		if(trans_gate_y!=NULL){
			inverseGate_y = trans_gate_y->getInverseTransformation();
		}
		polygonGate::transforming(inverseGate_x.get(), inverseGate_y.get());
		setTransformed(false);
		interpolated = true;
	}
	virtual unsigned short getType(){return CURLYQUADGATE;}
	CurlyGuadGate * clone(){return new CurlyGuadGate(*this);};

};
};
#endif /* GATE_HPP_ */
