/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * MemCytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#include "CytoFrame.hpp"
#include "trans_group.hpp"
#include "readFCSdata.hpp"

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
	MemCytoFrame(const string &filename, const FCS_READ_PARAM & config);

	void convertToPb(pb::CytoFrame & fr_pb
			, const string & h5_filename
			, CytoFileOption h5_opt
			, const CytoCtx & ctx = CytoCtx()) const;
	unsigned n_rows() const;

	void read_fcs();

	void read_fcs_data();

	FileFormat get_backend_type() const
	{
		return FileFormat::MEM;
	};
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

	CytoFramePtr copy(const string & h5_filename = "", bool overwrite = false) const
	{
		CytoFramePtr res(new MemCytoFrame(*this));
		res->set_readonly(false);
		return res;
	}
	CytoFramePtr copy(uvec idx, bool is_row_indexed, const string & h5_filename = "", bool overwrite = false) const
	{
		unique_ptr<MemCytoFrame> ptr(new MemCytoFrame(*this));
		ptr->set_readonly(false);
		ptr->realize_(idx, is_row_indexed);

		return CytoFramePtr(ptr.release());
	}

	CytoFramePtr copy(uvec row_idx, uvec col_idx, const string & h5_filename = "", bool overwrite = false) const
	{
		unique_ptr<MemCytoFrame> ptr(new MemCytoFrame(*this));
		ptr->set_readonly(false);
		ptr->realize_(row_idx, col_idx);

		return CytoFramePtr(ptr.release());
	}

	/**
	 * realize in place
	 * @param row_idx
	 * @param col_idx
	 * @param h5_filename
	 * @return
	 */
	void realize_(uvec row_idx, uvec col_idx);
	void realize_(uvec idx, bool is_row_indexed);

	/**
 * Caller will receive a copy of data
 * @return
 */

	EVENT_DATA_VEC get_data() const
	{
		return data_;
	}
	EVENT_DATA_VEC & get_data_ref()
	{
		return data_;
	}


	EVENT_DATA_VEC get_data(uvec idx, bool is_col) const
	{
		if(is_col)
			return data_.cols(idx);
		else
			return data_.rows(idx);
	}
	EVENT_DATA_VEC get_data(uvec row_idx, uvec col_idx) const
	{
		return data_.submat(row_idx, col_idx);
	}

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data)
	{
		data_ = _data;
	}
	/**
	 * move setter
	 * @param _data
	 */
	void set_data(EVENT_DATA_VEC && _data)
	{
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
	string get_uri() const{
		return "";
	}

	void append_data_columns(const EVENT_DATA_VEC & new_cols);

	void transform_data(const trans_local & trans);
};


};




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
