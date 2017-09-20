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
	readHeaderAndText(in, header, CytoFrame::keys, config.header);



	if(!onlyTxt)
	{
		//parse the data section
		data = move(readFCSdata(in, header, CytoFrame::keys, config.data));
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
