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
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
namespace cytolib
{
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
	EVENT_DATA_VEC data_;//col-major

	// below are cached for fcs parsing, should be of no usage once the data section is parsed
	string filename_;
	FCS_READ_PARAM config_;
	FCS_Header header_;
	ifstream in_;//because of this member, the class needs to explicitly define copy/assignment constructor

	void parse_fcs_header(ifstream &in, int nOffset = 0);
	void string_to_keywords(string txt, bool emptyValue);
	void parse_fcs_text_section(ifstream &in, bool emptyValue);
	void open_fcs_file();

public:
//	void close_h5(){};
	MemCytoFrame(){}
	MemCytoFrame(const MemCytoFrame & frm);
	MemCytoFrame(MemCytoFrame && frm);
	MemCytoFrame & operator=(const MemCytoFrame & frm);
	MemCytoFrame & operator=(MemCytoFrame && frm);
	/**
	 * Constructor from a generic CytoFrame object
	 * @param frm a reference to CytoFrame
	 */
	MemCytoFrame(const CytoFrame & frm):CytoFrame(frm)
	{
		data_ = frm.get_data();
	}
	/**
	 * Constructor from the FCS file
	 *
	 * @param filename FCS file path
	 * @param config the parse arguments.
	 * @param onlyTxt flag indicates whether to only parse text segment (which contains the keywords)
	 */
	MemCytoFrame(const string &filename, const FCS_READ_PARAM & config):filename_(filename),config_(config){
		set_pheno_data("name", path_base_name(filename));
	}

	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const;
	unsigned n_rows() const;

	void read_fcs();

	void read_fcs_data();

	/**
	 * parse the data segment of FCS
	 *
	 * @param in (input) file stream object opened from FCS file
	 * @param config (input) the parsing arguments for data
	 */
	void read_fcs_data(ifstream &in, const FCS_READ_DATA_PARAM & config);
	void read_fcs_header();
	/**
	 * parse the FCS header and Text segment
	 *
	 * @param in (input) the file stream object opened from FCS file
	 * @param config (input) FCS_READ_HEADER_PARAM object gives the parsing arguments for header
	 */
	void read_fcs_header(ifstream &in, const FCS_READ_HEADER_PARAM & config);

	CytoFramePtr copy(const string & h5_filename = "") const
	{
		CytoFramePtr res(new MemCytoFrame(*this));
		return res;
	}

	/**
	 * realize in place
	 * @param row_idx
	 * @param col_idx
	 * @param h5_filename
	 * @return
	 */
	void realize_(uvec row_idx, uvec col_idx);

	CytoFramePtr copy_realized(uvec row_idx, uvec col_idx, const string & h5_filename = "") const;

	/**
 * Caller will receive a copy of data
 * @return
 */

	EVENT_DATA_VEC get_data() const
	{
		return data_;
	}


	EVENT_DATA_VEC get_data(uvec col_idx) const;
	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{
		check_write_permission();
		data_ = _data;
	}
	/**
	 * move setter
	 * @param _data
	 */
	void set_data(EVENT_DATA_VEC && _data)
	{
		check_write_permission();
		swap(data_, _data);
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
	EVENT_DATA_TYPE * get_data_memptr(const string & colname, ColType type);
	string get_h5_file_path() const{
		return "";
	}

};

#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif
};




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
