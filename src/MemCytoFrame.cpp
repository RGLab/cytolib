/*
 * MemCytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/MemCytoFrame.hpp"

MemCytoFrame::~MemCytoFrame(){
	delete [] data;
}

MemCytoFrame::MemCytoFrame(const string & filename){

}


EVENT_DATA_TYPE * MemCytoFrame::getData(){
	return data;
}
EVENT_DATA_TYPE * MemCytoFrame::getData(const string & colname, ColType type){
	int idx = getColId(colname, type);
	return data + idx * nRow();
}

void MemCytoFrame::setChannel(const string &, const string &){

}
void MemCytoFrame::setMarker(const string &, const string &){

}

