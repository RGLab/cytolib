/*
 * readfcs.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */
#include <cytolib/MemCytoFrame.hpp>
int main(void)
{
	float a = 0;
	string filename = "/shared/silo_researcher/Gottardo_R/mike_working/flowCore_misc/sample_1071.001";
	FCS_READ_PARAM config;
	MemCytoFrame cytofrm(filename.c_str(), config,false);
	for(auto p : cytofrm.getKeywords())
		cout << p.first << ": " << p.second << endl;
	cout << cytofrm.nCol() << " " << cytofrm.nRow()<<endl;
}



