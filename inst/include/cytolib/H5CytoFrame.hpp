/*
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/CytoFrame.hpp>


/**
 * The class represents the H5 version of cytoFrame
 * It doesn't store and own the event data in memory.
 * Instead, data is read from H5 file on demand, which is more memory efficient.
 */
class H5CytoFrame:public CytoFrame{
protected:
	string filename;
	/*
	 * these H5 handlers remain open during the life cycle of H5CytoFrame
	 * for faster accessing the data
	 */
	H5File file;
	DataSet dataset;
	DataSpace dataspace;
	hsize_t dims[2];              // dataset dimensions
public:
	~H5CytoFrame(){};
	void compensate(const compensation &);
	int nRow();
	/**
	 * constructor from the H5 format of FCS
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & _filename);
//	void save(const string & filename, FrameType type);
	EVENT_DATA_VEC getData();
	EVENT_DATA_VEC getData(const string &, ColType);
//	void setKeyword(const string &, const string &);
//	void setChannel(const string &, const string &);
//	void setMarker(const string &, const string &);

};




#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
