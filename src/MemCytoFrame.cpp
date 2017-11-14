/*
 * MemCytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/MemCytoFrame.hpp"

MemCytoFrame::MemCytoFrame(const string &filename, FCS_READ_PARAM & config,  bool onlyTxt = false){
	ifstream in(filename, ios::in|ios::binary);

	if(!in.is_open())
		throw(domain_error("can't open the file: " + filename + "\nPlease check if the path is normalized to be recognized by c++!"));
	FCS_Header header;
	readHeaderAndText(in, header, keys, params, config.header);



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
int MemCytoFrame::nRow(){

	return data.size()/nCol();
}


EVENT_DATA_VEC MemCytoFrame::getData(){
	return data;
}
EVENT_DATA_VEC MemCytoFrame::getData(const string & colname, ColType type){
	int idx = getColId(colname, type);
	if(idx<0)
		throw(domain_error("colname not found: " + colname));
	int nEvents = nRow();
	EVENT_DATA_VEC res(nEvents);
	memcpy(&res[0], &data[0] + idx * nEvents, nEvents*sizeof(EVENT_DATA_TYPE));
	return res ;
}



void MemCytoFrame:: compensate(const compensation &){

}


