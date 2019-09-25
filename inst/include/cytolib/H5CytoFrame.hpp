/*
 * H5CytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_
#include <cytolib/MemCytoFrame.hpp>


namespace cytolib
{
/**
 * The class represents the H5 version of cytoFrame
 * It doesn't store and own the event data in memory.
 * Instead, data is read from H5 file on demand, which is more memory efficient.
 */
class H5CytoFrame:public CytoFrame{
protected:
	string filename_;
	/*TODO: We may not want to maintain these handlers, instead treat each IO as atomic operations
	 * Because it is not easy to accomplish the resource sharing among multiple H5CytoFrame objects solely depending on H5's mechanisms.
	 * e.g. a second openFile call with H5F_ACC_RDONLY will NOT overwrite the previous H5F_ACC_RDWR , thus cause the unexpected data tampering
	 * these H5 handlers remain open during the life cycle of H5CytoFrame
	 * for faster accessing the data
	 */
	H5File file;
	DataSet dataset;
	DataSpace dataspace;
	hsize_t dims[2];              // dataset dimensions
	//flags indicating if cached meta data needs to be flushed to h5
	bool is_dirty_params;
	bool is_dirty_keys;
	bool is_dirty_pdata;
	EVENT_DATA_VEC read_data(uvec col_idx) const;
public:
	~H5CytoFrame(){
		/*
		 * catch the exception to prevent the destructor from throwing, which could crash the application
		 */
//		string msg = "Warning: failed to flush the meta data to h5!Changes to meta are unsaved.";
//
//		try{
//			flush_meta();
//		}catch(const H5::DataSetIException &e){
//			PRINT(e.getDetailMsg() + "\n");
//			PRINT(msg);
//		}catch(...){
//			PRINT(msg);
//		}

	};
	/*
	 * for simplicity, we don't want to handle the object that has all the h5 handler closed
	 * because it will require lots of validity checks before each disk IO operations
	 */
//	void close_h5(){
//		dataspace.close();
//		dataset.close();
//		file.close();
//	}
	void flush_meta();
	void flush_params();

	void flush_keys();
	void flush_pheno_data();

	H5CytoFrame(const H5CytoFrame & frm);
	H5CytoFrame(H5CytoFrame && frm);
	H5CytoFrame & operator=(const H5CytoFrame & frm);
	H5CytoFrame & operator=(H5CytoFrame && frm);
	unsigned n_rows() const{
				return dims[1];
		}

	void set_params(const vector<cytoParam> & _params)
	{
		CytoFrame::set_params(_params);
		is_dirty_params = true;

	}
	void set_keywords(const KEY_WORDS & keys){
		CytoFrame::set_keywords(keys);
		is_dirty_keys = true;

	}
	void set_keyword(const string & key, const string & value){
		CytoFrame::set_keyword(key, value);
		is_dirty_keys = true;

	}
	void set_channel(const string & oldname, const string &newname, bool is_update_keywords = true)
	{
		CytoFrame::set_channel(oldname, newname, is_update_keywords);
		is_dirty_params = true;
		if(is_update_keywords)
			is_dirty_keys = true;
	}
	void set_marker(const string & channelname, const string & markername)
	{
		CytoFrame::set_marker(channelname, markername);
		is_dirty_params = true;
	}
	void set_range(const string & colname, ColType ctype, pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> new_range, bool is_update_keywords = true){
		CytoFrame::set_range(colname, ctype, new_range, is_update_keywords);
		is_dirty_params = true;
		if(is_update_keywords)
			is_dirty_keys = true;

	}
	void set_pheno_data(const string & name, const string & value){
		CytoFrame::set_pheno_data(name, value);
		is_dirty_pdata = true;
	}
	void set_pheno_data(const PDATA & _pd)
	{
		CytoFrame::set_pheno_data(_pd);
		is_dirty_pdata = true;
	}
	void del_pheno_data(const string & name){
		CytoFrame::del_pheno_data(name);
		is_dirty_pdata = true;
	}
	/**
	 * constructor from FCS
	 * @param fcs_filename
	 * @param h5_filename
	 */
	H5CytoFrame(const string & fcs_filename, FCS_READ_PARAM & config, const string & h5_filename
			, bool readonly = false);
	/**
	 * constructor from the H5
	 * @param _filename H5 file path
	 */
	H5CytoFrame(const string & h5_filename, bool readonly = true);
	/**
	 * abandon the changes to the meta data in cache by reloading them from disk
	 */
	void load_meta();;;

	string get_h5_file_path() const{
		return filename_;
	}

	void convertToPb(pb::CytoFrame & fr_pb, const string & h5_filename, H5Option h5_opt) const;

	/**
	 * Read data from disk.
	 * The caller will directly receive the data vector without the copy overhead thanks to the move semantics supported for vector container in c++11
	 * @return
	 */


	EVENT_DATA_VEC get_data() const;
	/**
	 * Partial IO
	 * @param col_idx
	 * @return
	 */
	EVENT_DATA_VEC get_data(uvec col_idx) const
	{
		return read_data(col_idx);
	}

	CytoFramePtr copy(const string & h5_filename = "") const;
	CytoFramePtr copy_realized(uvec row_idx, uvec col_idx, const string & h5_filename = "") const;

	/**
	 * copy setter
	 * @param _data
	 */
	void set_data(const EVENT_DATA_VEC & _data);

	void set_data(EVENT_DATA_VEC && _data)
	{
		check_write_permission();
		set_data(_data);
	}
};

};



#endif /* INST_INCLUDE_CYTOLIB_H5CYTOFRAME_HPP_ */
