// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/gate.hpp>
#include <cytolib/global.hpp>

#include <cytolib/ellipse2points.hpp>
#include <boost/foreach.hpp>


namespace cytolib
{
	void BOOL_GATE_OP::convertToPb(pb::BOOL_GATE_OP & BOOL_GATE_OP_pb){
		BOOL_GATE_OP_pb.set_isnot(isNot);
		BOOL_GATE_OP_pb.set_op(op);
		for(unsigned i = 0; i < path.size(); i++){
			 BOOL_GATE_OP_pb.add_path(path[i]);
		}
	};
	BOOL_GATE_OP::BOOL_GATE_OP(const pb::BOOL_GATE_OP & BOOL_GATE_OP_pb){
		op = BOOL_GATE_OP_pb.op();
		isNot = BOOL_GATE_OP_pb.isnot();
		for(int i = 0; i < BOOL_GATE_OP_pb.path_size(); i++)
			path.push_back(BOOL_GATE_OP_pb.path(i));
	};



	void vertices_vector::resize(unsigned nSize){
		x.resize(nSize);
		y.resize(nSize);
	}
	vertices_vector::vertices_vector(vector<coordinate> vertices){

			unsigned nSize=vertices.size();
			resize(nSize);
			for(unsigned i=0;i<nSize;i++)
			{
				x[i]=vertices[i].x;
				y[i]=vertices[i].y;
			}

	};
	void vertices_vector::print(){
		PRINT("x:");
		for(unsigned i=0;i<x.size();i++)
				PRINT(to_string(x[i])+",");
//		PRINT("x:");
//		for(unsigned i=0;i<x.size();i++)
//				PRINT(x[i]+",");

	}

	vertices_vector paramRange::toVector() const{

		vertices_vector res;
		res.resize(2);
		res.x[0]=min;
		res.x[1]=max;

		return res;
	}
	void paramRange::update_channels(const CHANNEL_MAP & chnl_map){

			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(name);
			if(itChnl!=chnl_map.end())
				name = itChnl->second;
	};
	vector<string> paramRange::getNameArray() const{
			vector<string> res;
			res.push_back(name);
			return res;
		};
	void paramPoly::update_channels(const CHANNEL_MAP & chnl_map){

			for(vector<string>::iterator it = params.begin(); it != params.end(); it++)
			{
				string curName = *it;

				CHANNEL_MAP::const_iterator itChnl = chnl_map.find(curName);
				if(itChnl!=chnl_map.end())
					*it = itChnl->second;
			}
		};
	vertices_vector paramPoly::toVector() const{

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

	void paramPoly::convertToPb(pb::paramPoly & paramPoly_pb){
		BOOST_FOREACH(vector<string>::value_type & it, params){
			paramPoly_pb.add_params(it);
		}
		BOOST_FOREACH(vector<coordinate>::value_type & it, vertices){
			pb::coordinate * coor_pb = paramPoly_pb.add_vertices();
			it.convertToPb(*coor_pb);
		}
	};
	paramPoly::paramPoly(const pb::paramPoly & paramPoly_pb){
		for(int i = 0; i < paramPoly_pb.params_size(); i++){
			params.push_back(paramPoly_pb.params(i));
		}
		for(int i = 0; i < paramPoly_pb.vertices_size(); i++){
			vertices.push_back(coordinate(paramPoly_pb.vertices(i)));
		}
	};
	void gate::convertToPb(pb::gate & gate_pb){
		//cp basic members
			gate_pb.set_istransformed(isTransformed);
			gate_pb.set_neg(neg);
			gate_pb.set_isgained(isGained);
	}

	void rangeGate::convertToPb(pb::gate & gate_pb){
		gate::convertToPb(gate_pb);
		gate_pb.set_type(pb::RANGE_GATE);
		//cp nested gate
		pb::rangeGate * g_pb = gate_pb.mutable_rg();
		//cp its unique member
		pb::paramRange * pr_pb = g_pb->mutable_param();
		param.convertToPb(*pr_pb);
	}
	void rangeGate::transforming(trans_local & trans){
		if(!Transformed())
		{
			EVENT_DATA_TYPE vert[2] = {param.getMin(),param.getMax()};

			TransPtr curTrans=trans.getTran(param.getName());
			if(curTrans)
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

	void rangeGate::shiftGate(){
		param.setMin(param.getMin() + shift[0]);
		param.setMax(param.getMax() + shift[0]);
	}

	INDICE_TYPE rangeGate::gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){

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

	void rangeGate::extend(MemCytoFrame & fdata,float extend_val){
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
	void rangeGate::extend(float extend_val, float extend_to){
		string pName=param.getName();


		EVENT_DATA_TYPE xMin= extend_to;
		if(param.getMin()<=extend_val)
		{
			if(g_loglevel>=POPULATION_LEVEL)
				PRINT("extending "+pName+"from "+to_string(param.getMin())+" to :"+to_string(xMin)+"\n");
			param.setMin(min(xMin, param.getMin()));
		}


	}
	void rangeGate::gain(map<string,float> & gains){
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
/*
	 * when the original gate vertices are at the threshold
	 * it is likely that the gates were truncated in flowJo xml
	 * currently what we can do is to extend it to the real data range to avoid losing
	 * the data points that are below this theshold range
	 * to cut data range)
	 */
	void polygonGate::extend(MemCytoFrame & fdata,float extend_val){
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

	void polygonGate::extend(float extend_val, float extend_to){
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
	void polygonGate::gain(map<string,float> & gains){

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

	void polygonGate::shiftGate(){
		vector<coordinate> vertices=param.getVertices();
		for(unsigned i=0;i<vertices.size();i++){
			vertices[i].x += shift[0];
			vertices[i].y += shift[1];
		}
		param.setVertices(vertices);
	}

	 /*
	 *
	 *  reimplement c++ version of inPolygon_c
	 *  indices are allocated within gating function, so it is up to caller to free it
	 *  and now it is freed in destructor of its owner "nodeProperties" object
	 */
	INDICE_TYPE polygonGate::gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){




		vector<coordinate> vertices=param.getVertices();


		string x=param.xName();
		string y=param.yName();
		EVENT_DATA_TYPE * xdata = fdata.get_data_memptr(x, ColType::channel);
		EVENT_DATA_TYPE * ydata = fdata.get_data_memptr(y, ColType::channel);

		int nEvents=parentInd.size();
		INDICE_TYPE res;
		res.reserve(nEvents);
		unsigned nVert = vertices.size();
		vector<cytolib::CYTO_POINT> points(nVert);
		for(unsigned i = 0; i < nVert; i++)
			points[i] = vertices[i];
		cytolib::in_polygon(xdata, ydata, points, parentInd, neg, res);
		return res;
	}

	/*
	 * a wrapper that calls transforming(TransPtr , TransPtr )
	 */
	void polygonGate::transforming(trans_local & trans){

			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();


			/*
			 * do the actual transformations
			 */
			TransPtr trans_x=trans.getTran(channel_x);
			TransPtr trans_y=trans.getTran(channel_y);

			transforming(trans_x, trans_y);
	}

	/*
	 * the actual transforming logic for polygonGate, that is shared by polyonGate and ellipsoidGate(due to the special scale)
	 */
	void polygonGate::transforming(TransPtr trans_x, TransPtr trans_y){
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


			if(trans_x)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_x+"\n");;
				trans_x->transforming(&vert.x[0], nSize);
				for(int i=0;i<nSize;i++)
					vertices[i].x=vert.x[i];
			}
			if(trans_y)
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
	void polygonGate::convertToPb(pb::gate & gate_pb){
		gate::convertToPb(gate_pb);

		gate_pb.set_type(pb::POLYGON_GATE);
		//cp nested gate
		pb::polygonGate * g_pb = gate_pb.mutable_pg();
		//cp its unique member
		pb::paramPoly * pr_pb = g_pb->mutable_param();
		param.convertToPb(*pr_pb);
	}


	void rectGate::convertToPb(pb::gate & gate_pb)
	{
		polygonGate::convertToPb(gate_pb);
			gate_pb.set_type(pb::RECT_GATE);
	}
	vector<coordinate> ellipseGate::getCovarianceMat() const{
		if(!Transformed())
			throw(domain_error("EllipseGate has not been transformed so covariance matrix is unavailable!"));
		return cov;};
	coordinate ellipseGate::getMu() const{
		if(!Transformed())
				throw(domain_error("EllipseGate has not been transformed so mu is unavailable!"));
		return mu;};
	EVENT_DATA_TYPE ellipseGate::getDist() const{
		if(!Transformed())
			throw(domain_error("EllipseGate has not been transformed so dist is unavailable!"));
		return dist;};
	vector<coordinate> ellipseGate::getAntipodalVerts() const{
		return antipodal_vertices;
	}
	ellipseGate::ellipseGate(coordinate _mu, vector<coordinate> _cov, EVENT_DATA_TYPE _dist):mu(_mu),cov(_cov), dist(_dist){
		isTransformed = true;
		isGained = true;
		neg = false;
		setShift(vector<EVENT_DATA_TYPE>{0.0,0.0});
	}

	ellipseGate::ellipseGate(vector<coordinate> _antipodal, vector<string> _params):antipodal_vertices(_antipodal),dist(1){
		isTransformed = false;
		isGained = false;
		neg = false;
		setShift(vector<EVENT_DATA_TYPE>{0.0,0.0});
		/*
		 * init the dummy vertices for base class
		 * (this deprecated inheritance exists for the sake of legacy archive)
		 */
		param.setName(_params);

	}

	void ellipseGate::extend(MemCytoFrame & fdata,float extend_val){

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
	void ellipseGate::extend(float extend_val, float extend_to){

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
	void ellipseGate::gain(map<string,float> & gains){
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
	void ellipseGate::computeCov(){
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
		vector<coordinate>::iterator R_it=max_element(v.begin(),v.end(),[](coordinate i, coordinate j) { return i.x<j.x; });
		coordinate R = *R_it;

		//far left point
		vector<coordinate>::iterator L_it=min_element(v.begin(),v.end(),[](coordinate i, coordinate j) { return i.x<j.x; });
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
		cov.clear();
		cov.push_back(p1);
		cov.push_back(p2);

		//set distance (in this calculation should always be 1)
		dist = 1;
	}

	void ellipseGate::computeAntipodalVerts(){
		antipodal_vertices.clear();
		/*
		 * Calculate the antipodal vertices by eigendecomposition
		 * of the covariance matrix
		 */
		Mat<EVENT_DATA_TYPE> covmat = {{cov[0].x, cov[0].y} , {cov[1].x,cov[1].y} };
		Col<EVENT_DATA_TYPE> evals;
		Mat<EVENT_DATA_TYPE> evecs;
		// Covariance matrix should be enforced symmetric by cytolib/flowCore
		// but otherwise could use eig_gen
		eig_sym(evals, evecs, covmat);

		// It sorts the evals, evecs in ascending mag order by default
		EVENT_DATA_TYPE a2 = dist*dist*evals(1);
		EVENT_DATA_TYPE a = sqrt(a2); // major semi-axis length
		EVENT_DATA_TYPE b2 = dist*dist*evals(0);
		EVENT_DATA_TYPE b = sqrt(b2); // minor semi-axis length
		EVENT_DATA_TYPE c = sqrt(a2-b2); // semi-focal length

		Col<EVENT_DATA_TYPE> a_vec(a*evecs.col(1));
		Col<EVENT_DATA_TYPE> b_vec(b*evecs.col(0));
		Col<EVENT_DATA_TYPE> c_vec(c*evecs.col(1));
		// Add the 4 antipodal vertices in the order FlowJo expects:
		// vertices on same axis are adjacent in vector
		antipodal_vertices.push_back(coordinate(mu.x+a_vec(0), mu.y+a_vec(1)));
		antipodal_vertices.push_back(coordinate(mu.x-a_vec(0), mu.y-a_vec(1)));
		antipodal_vertices.push_back(coordinate(mu.x+b_vec(0), mu.y+b_vec(1)));
		antipodal_vertices.push_back(coordinate(mu.x-b_vec(0), mu.y-b_vec(1)));
	}

	void ellipseGate::shiftGate(){
		mu.x += shift[0];
		mu.y += shift[1];
		computeAntipodalVerts();
	}

	/*
	 * translated from flowCore::%in% method for ellipsoidGate
	 */
	INDICE_TYPE ellipseGate::gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){


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
	void ellipseGate::convertToPb(pb::gate & gate_pb)
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
	ellipseGate::ellipseGate(const pb::gate & gate_pb):polygonGate(gate_pb),mu(coordinate(gate_pb.eg().mu())),dist(gate_pb.eg().dist()){
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
	void ellipseGate::toPolygon(unsigned nVertices){




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

	void ellipseGate::interpolatePolygon(unsigned nVertices = 50){
		Mat<EVENT_DATA_TYPE> covmat = {{cov[0].x, cov[0].y} , {cov[1].x,cov[1].y} };
		Mat<EVENT_DATA_TYPE> chol_lower = chol(covmat, "lower");
		Row<EVENT_DATA_TYPE> theta = linspace<Row<EVENT_DATA_TYPE>>(0.0, 2.0*datum::pi, nVertices);
		Row<EVENT_DATA_TYPE> xvals = dist*cos(theta);
		Row<EVENT_DATA_TYPE> yvals = dist*sin(theta);
		Mat<EVENT_DATA_TYPE> result = join_cols(xvals, yvals);
		result = chol_lower * result;
		vector<coordinate> vertices=param.getVertices();
		vertices.clear();
		vertices.resize(nVertices);
		for(unsigned short i=0;i<nVertices;i++)
		{
			vertices[i].x = result(0, i) + mu.x;
			vertices[i].y = result(1, i) + mu.y;
		}
		param.setVertices(vertices);
	}

	void ellipseGate::transforming(trans_local & trans){
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
			TransPtr trans_x=trans.getTran(channel_x);
			TransPtr trans_y=trans.getTran(channel_y);


			if(trans_x)
			{
				if(g_loglevel>=POPULATION_LEVEL)
					PRINT("transforming: "+channel_x+"\n");;

				trans_x->transforming(&vert.x[0],nSize);
				for(int i=0;i<nSize;i++)
					antipodal_vertices[i].x=vert.x[i];
			}
			if(trans_y)
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

	void ellipsoidGate::convertToPb(pb::gate & gate_pb)
	{
		ellipseGate::convertToPb(gate_pb);
			gate_pb.set_type(pb::ELLIPSOID_GATE);
	}
	ellipsoidGate::ellipsoidGate(const pb::gate & gate_pb):ellipseGate(gate_pb){
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
	void ellipsoidGate::transforming(trans_local & trans){

		if(!Transformed())
		{
			/*
			 * get channel names to select respective transformation functions
			 */
			string channel_x=param.xName();
			string channel_y=param.yName();



			TransPtr trans_x=trans.getTran(channel_x);
			TransPtr trans_y=trans.getTran(channel_y);


			/*
			 * re-construct the trans object that was used by flowJo to transform ellipsoid gate to 256 scale
			 */
			TransPtr trans_gate_x,trans_gate_y;
			if(!trans_x)
				throw(domain_error("ellipsoidGate::transforming can't find transformation for " + channel_x));
	//			trans_gate_x.reset(new scaleTrans()); //create default scale trans for linear, assuming the max value for linear scale is always 262144
	//		else
				trans_gate_x = trans_x->clone(); //copy existing trans_x for non-linear
	//
			if(!trans_y)
				throw(domain_error("ellipsoidGate::transforming can't find transformation for " + channel_y));
	//			trans_gate_y.reset(new scaleTrans()); //create default scale trans for linear
	//		else
				trans_gate_y = trans_y->clone(); //copy existing trans_y for non-linear

			//set to scale 256
			trans_gate_x->setTransformedScale(256);
			trans_gate_y->setTransformedScale(256);

			//get its inverse
			TransPtr inverseTrans_x = trans_gate_x->getInverseTransformation();
			TransPtr inverseTrans_y = trans_gate_y->getInverseTransformation();


			/*
			 * transform the polygon from 256 to raw
			 */
			polygonGate::transforming(inverseTrans_x, inverseTrans_y);



			/*
			 * transform the raw to the actual data scale (for non-linear channel)
			 */
			isTransformed = false;//reset transform flag otherwise the transforming won't get executed
			polygonGate::transforming(trans_x, trans_y);

			isTransformed=true;
		}

	}

	void ellipsoidGate::shiftGate(){
		//Let polygonGate::shiftGate() shift the interpolated points
		polygonGate::shiftGate();

		//Still need to shift the antipodal vertices
		for(auto & i : antipodal_vertices)
		{
			i.x += shift[0];
			i.y += shift[1];
		}
	}

	void boolGate::convertToPb(pb::gate & gate_pb){
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
	boolGate::boolGate(const pb::gate & gate_pb):gate(gate_pb){
		const pb::boolGate & bg_pb = gate_pb.bg();
		for(int i = 0; i < bg_pb.boolopspec_size(); i++){
			const pb::BOOL_GATE_OP & thisOP_pb = bg_pb.boolopspec(i);
			BOOL_GATE_OP thisOP = BOOL_GATE_OP(thisOP_pb);
			boolOpSpec.push_back(thisOP);


		}
	}
	void logicalGate::convertToPb(pb::gate & gate_pb){
		boolGate::convertToPb(gate_pb);
		gate_pb.set_type(pb::LOGICAL_GATE);
	}


	void clusterGate::convertToPb(pb::gate & gate_pb){
		boolGate::convertToPb(gate_pb);
		gate_pb.set_type(pb::CLUSTER_GATE);
		//cp nested gate
		pb::clusterGate * g_pb = gate_pb.mutable_cg();

		g_pb->set_cluster_method(cluster_method_name_);

	}


	void CurlyQuadGate::transforming(trans_local & trans){
		if(interpolated)
			polygonGate::transforming(trans);
		else
			throw(logic_error("CurlyQuadGate can't not be transformed before interpolation!"));
	};



	INDICE_TYPE CurlyQuadGate::gating(MemCytoFrame & fdata, INDICE_TYPE & parentInd){
		if(interpolated)
		{
			return polygonGate::gating(fdata, parentInd);
		}
		else
		{
			throw(logic_error("CurlyQuad gate has not been converted to polygonGate yet!"));
		}


	}


	void CurlyQuadGate::interpolate(trans_local & trans){

		string x_chnl = param.xName();
		string y_chnl = param.yName();
		/*
		 * transform intersect back to raw
		 */

		TransPtr trans_x = trans.getTran(x_chnl);
		TransPtr trans_y = trans.getTran(y_chnl);


		/*
		 * and rescale raw to 256 space
		 */
		TransPtr trans_gate_x,trans_gate_y;
		if(!trans_x)
			trans_gate_x.reset(new scaleTrans()); //create default scale trans for linear, assuming the max value for linear scale is always 262144
		else
			trans_gate_x = trans_x->clone(); //copy existing trans_x for non-linear

		if(!trans_y)
			trans_gate_y.reset(new scaleTrans()); //create default scale trans for linear
		else
			trans_gate_y = trans_y->clone(); //copy existing trans_y for non-linear

		//set to scale 256
		int displayScale = 255;
		trans_gate_x->setTransformedScale(displayScale);
		trans_gate_y->setTransformedScale(displayScale);
		polygonGate::transforming(trans_gate_x, trans_gate_y);

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
		TransPtr inverseGate_x,inverseGate_y;
		if(trans_gate_x){
			inverseGate_x = trans_gate_x->getInverseTransformation();
		}
		if(trans_gate_y){
			inverseGate_y = trans_gate_y->getInverseTransformation();
		}
		polygonGate::transforming(inverseGate_x, inverseGate_y);
		setTransformed(false);
		interpolated = true;
	}

};
