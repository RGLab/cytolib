/*
 * MemCytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/MemCytoFrame.hpp"

MemCytoFrame::MemCytoFrame(const string &filename, FCS_READ_PARAM & config,  bool onlyTxt = false){
	ifstream in(filename, ios::in|ios::binary);


	FCS_Header header;
	readHeaderAndText(in, header, keys, params, config.header);



	if(!onlyTxt)
	{
		//parse the data section
		data = move(readFCSdata(in, header, keys, params, nEvents, config.data));

		//update min and max
		for(int i = 0; i < params.size(); i++)
		{

			string pid = to_string(i);


			if(keys.find("transformation")!=keys.end() &&  keys["transformation"] == "custom")
				params[i].min = boost::lexical_cast<EVENT_DATA_TYPE>(keys["flowCore_$P" + pid + "Rmin"]);
			else
			{

				EVENT_DATA_TYPE * vec = getData(params[i].channel, ColType::channel);
				auto absMin  = *min_element(vec, vec + nRow());
				auto zeroVals = params[i].PnE.second;
				params[i].min = min(static_cast<EVENT_DATA_TYPE>(zeroVals), max(config.data.min_limit, absMin));

			}


			params[i].max--;
		}




	}
}


EVENT_DATA_TYPE * MemCytoFrame::getData(){
	return data.get();
}
EVENT_DATA_TYPE * MemCytoFrame::getData(const string & colname, ColType type){
	int idx = getColId(colname, type);
	return data.get() + idx * nRow();
}



void MemCytoFrame:: compensate(const compensation &){

}

void MemCytoFrame:: save(const string & filename, FrameType type){

}
