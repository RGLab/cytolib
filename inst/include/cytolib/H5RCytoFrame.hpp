/*
 * H5RCytoFrame.hpp
 *
 *  Created on: Mar 4, 2020
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_H5RCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_H5RCYTOFRAME_HPP_

#include <cytolib/H5CytoFrame.hpp>
#include <H5FDros3.h>
namespace cytolib
{
/**
 * The class represents the H5 version of cytoFrame
 * It doesn't store and own the event data in memory.
 * Instead, data is read from H5 file on demand, which is more memory efficient.
 */
class H5RCytoFrame:public H5CytoFrame{
	string access_key_id_;
	string access_key_;
	string region_;
public:
	void flush_meta(){throw(domain_error("Writing to the s3 based H5CytoFrame object is not supported yet!"));};
	void flush_params(){throw(domain_error("Writing to the s3 based H5CytoFrame object is not supported yet!"));};
	void flush_keys(){throw(domain_error("Writing to the s3 based H5CytoFrame object is not supported yet!"));};
	void set_data(const EVENT_DATA_VEC & _data){throw(domain_error("Writing to the s3 based H5CytoFrame object is not supported yet!"));};
	void set_data(EVENT_DATA_VEC && _data){throw(domain_error("Writing to the s3 based H5CytoFrame object is not supported yet!"));};
	H5RCytoFrame(const string & h5_filename, bool readonly = true
				, const string & secret_id = ""
				, const string & secret_key = ""
				, const string & aws_region = "us-west-1"):H5CytoFrame(h5_filename, readonly, false)
				, access_key_id_(secret_id), access_key_(secret_key), region_(aws_region){
		//set fapl for s3
//		access_plist_ = FileAccPropList::DEFAULT;
//		access_plist_.get
		hid_t fapl_id = access_plist_.getId();
		H5FD_ros3_fapl_t fa;
		fa.version = 1;
		fa.authenticate = secret_id != "";
		strcpy(fa.aws_region, aws_region.c_str());
		if(fa.authenticate)
		{
			strcpy(fa.secret_id, secret_id.c_str());
			strcpy(fa.secret_key, secret_key.c_str());

		}
		H5Pset_fapl_ros3(fapl_id, &fa);
		//reset flag
		default_flags = H5F_ACC_RDONLY;
		init_load();
	}
};

}

#endif /* INST_INCLUDE_CYTOLIB_H5RCYTOFRAME_HPP_ */
