/*
 * MemCytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/MemCytoFrame.hpp"

MemCytoFrame::~MemCytoFrame(){
	if(data)
		delete [] data;
}
void readFCSHeader(ifstream in, FCS_Header & header){
	/*
		 * parse the header
		 */

		//parse version
		char version[7];
		in.get(version, 6);
		version[6]='\n';

	    if(strcmp(version, "FCS2.0")&&strcmp(version, "FCS3.0")&&strcmp(version, "FCS3.1"))
		     throw(domain_error("This does not seem to be a valid FCS2.0, FCS3.0 or FCS3.1 file"));

	    header.FCSversion = stof(version + 3);

	    char tmp[5];
	    tmp[4]='\n';
		in.get(tmp, 4);
		if(strcmp(tmp, "    "))
			 throw(domain_error("This does not seem to be a valid FCS header"));

		//parse offset
		char tmp1[9];
		tmp1[8]='\n';
		in.get(tmp1, 8);
	    header.textstart = stof(t)
		        coffs[i] <- readChar(con=con, nchars=8)
	//
	//	    ioffs <- c(as.double(version), as.integer(coffs), as.integer(start))
	//	    names(ioffs) <- c("FCSversion", "textstart", "textend", "datastart",
	//	                      "dataend", "anastart", "anaend", "additional")
	//	    ioffs[2:7] <- ioffs[2:7]+ioffs[8]
	//
	//	    if(all(is.na(ioffs[2:5]) || ioffs[2:5]==""))
	//	        stop("Missing header information to start parsing the binary ",
	//	             "section of the file")

}
MemCytoFrame::MemCytoFrame(const string &filename, bool emptyValue = false, int nDataset = 1, bool scale =false, double decades=0, double min_limit=-111, bool truncate_max_range = true, bool onlyTxt = false){
	ifstream in(filename, ios::in|ios::binary);
	FCS_Header header;
	readFCSHeader(in, header);
		if(onlyTxt)
		data = NULL;
	else
	{
		//parse the data section
	}
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

