/*
 * MemCytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#include "CytoFrame.hpp"
#include "readFCSdata.hpp"

/**
 * the container that stores the different FCS parse arguments
 */
struct FCS_READ_PARAM{
	FCS_READ_HEADER_PARAM header;
	FCS_READ_DATA_PARAM data;
};

/**
 * The class represents the in-memory version of CytoFrame, which stores and owns the events data
 */
class MemCytoFrame: public CytoFrame{
	EVENT_DATA_VEC data;

public:

	/**
	 * Constructor from a generic CytoFrame object
	 * @param frm a reference to CytoFrame
	 */
	MemCytoFrame(const CytoFrame & frm)
	{
		pd = frm.getPData();
		keys.setPairs(frm.getKeywords());
		params = frm.getParams();
		buildHash();
		data = frm.getData();
	}
	/**
	 * Constructor from the FCS file
	 *
	 * @param filename FCS file path
	 * @param config the parse arguments.
	 * @param onlyTxt flag indicates whether to only parse text segment (which contains the keywords)
	 */
	MemCytoFrame(const string &filename, FCS_READ_PARAM & config,  bool onlyTxt = false){
		ifstream in(filename, ios::in|ios::binary);

		if(!in.is_open())
			throw(domain_error("can't open the file: " + filename + "\nPlease check if the path is normalized to be recognized by c++!"));
		FCS_Header header;
		readHeaderAndText(in, header, keys, params, config.header);

		buildHash();

		if(!onlyTxt)
		{
	//		double start = clock();
			//parse the data section
			readFCSdata(in, data,header, keys, params, config.data);

	//		cout << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << endl;
			//update min and max
			/*
			 * ## set transformed flag and fix the PnE and the Datatype keywords
			 */
			keys["FILENAME"] = filename;
			if(config.data.isTransformed)
			{
				keys["transformation"] ="applied";
				keys["$DATATYPE"] = "F";
			}
			for(int i = 0; i < params.size(); i++)
			{

				string pid = to_string(i+1);

				//insert our own PnR fields
				if(config.data.isTransformed)
				{
					keys["$P" + pid + "E"] = "0,0";
					params[i].PnE[0] = 0;
					params[i].PnE[1] = 0;
					keys["flowCore_$P" + pid + "Rmax"] = to_string(static_cast<int>(params[i].max + 1));
					keys["flowCore_$P" + pid + "Rmin"] = to_string(static_cast<int>(params[i].min));
				}
				else
					params[i].max--;
			}



			//GUID
			string oldguid;
			if(keys.find("GUID")!=keys.end()){
				oldguid = keys["GUID"];
				keys["ORIGINALGUID"] = oldguid;
			}
			else
				oldguid = filename;
			//strip dir
			vector<string> paths;
			boost::split(paths, oldguid, boost::is_any_of("/\\"));
			keys["GUID"] = paths.back();

		}



	}
	int nRow() const{

		return data.size()/nCol();
	}


	EVENT_DATA_VEC getData() const{
		return data;
	}
	EVENT_DATA_VEC getData(const string & colname, ColType type) const{
		int idx = getColId(colname, type);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		int nEvents = nRow();
		EVENT_DATA_VEC res(nEvents);
		memcpy(&res[0], &data[0] + idx * nEvents, nEvents*sizeof(EVENT_DATA_TYPE));
		return res ;
	}



	void  compensate(const compensation &){

	}



	/**
	 * return the pointer of a particular data column
	 *
	 * The reason to return the pointer is to keep it backward compatible
	 * TODO:modify gating apis to use iterator instead of raw pointer
	 * @param colname
	 * @param type
	 * @return
	 */
	EVENT_DATA_TYPE * subset(const string & colname, ColType type = ColType::unknown){
		int idx = getColId(colname, type);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		return data.data() + idx * nRow();
	}


};

#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
