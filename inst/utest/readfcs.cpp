/*
 * readfcs.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */
#include <cytolib/MemCytoFrame.hpp>
int main(void)
{
	string filename = "/shared/silo_researcher/Gottardo_R/mike_working/flowCore_misc/sample_1071.001";
	MemCytoFrame cytofrm(filename.c_str());
	cout << cytofrm.nCol() << " " << cytofrm.nRow()<<endl;
}



