/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * GatingSet.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: wjiang2
 */



#ifndef GATINGSET_HPP_
#define GATINGSET_HPP_
#include "GatingHierarchy.hpp"
#include <cytolib/CytoFrameView.hpp>
#include <string>
#include <cytolib/delimitedMessage.hpp>
#include <cytolib/global.hpp>

using namespace std;

namespace cytolib
{
#define ARCHIVE_TYPE_BINARY 0
#define ARCHIVE_TYPE_TEXT 1
#define ARCHIVE_TYPE_XML 2
#define PB true
#define BS false

/**
 * \class GatingSet
 * \brief A container class that stores multiple GatingHierarchy objects.
 *
 *
 */
class GatingSet{
	typedef unordered_map<string, GatingHierarchyPtr> ghMap;
	ghMap ghs_;
	vector<string> sample_names_;
	GatingHierarchyPtr get_first_gh() const;
	string uid_;
public:
	typedef typename ghMap::iterator iterator;
	typedef typename ghMap::const_iterator const_iterator;

	/*
	 * forwarding APIs
	 */
	 size_t size() const{return ghs_.size();}
	 iterator end(){return ghs_.end();}
	 iterator begin(){return ghs_.begin();}
	 const_iterator end() const{return ghs_.end();}
	 const_iterator begin() const{return ghs_.begin();}
	 iterator find(const string &sample_uid){
			 return ghs_.find(sample_uid);
	 }
	 const_iterator find(const string &sample_uid) const{
					 return ghs_.find(sample_uid);
			 }
	 size_t erase ( const string& k ){return ghs_.erase(k);}

	string get_uid(){return uid_;}
	void set_uid(const string & uid){uid_ = uid;}

	GatingSet(){
		uid_ = generate_uid();
	};
	bool is_cytoFrame_only() const{return size() == 0||get_first_gh()->is_cytoFrame_only();};
	/**
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param is_skip_data whether to skip writing cytoframe data to pb. It is typically remain as default unless for debug purpose (e.g. re-writing gs that is loaded from legacy pb archive without actual data associated)
	 */
	void serialize_pb(string path, H5Option h5_opt, bool is_skip_data = false);
	/**
	 * constructor from the archives (de-serialization)
	 * @param path
	 * @param is_skip_data whether to skip loading cytoframe data from h5. It should typically remain as default unless for debug purpose (e.g. legacy pb archive)
	 */
	GatingSet(string path, bool is_skip_data = false, bool readonly = true);
	/**
	 * constructor from the legacy archives (de-serialization)
	 * @param filename
	 * @param format
	 * @param isPB
	 */
	GatingSet(string pb_file, const GatingSet & gs_data);

	/*
	 * up to caller to free the memory
	 */
	GatingSet copy(bool is_copy_data = true, bool is_realize_data = true
			, const string & new_h5_dir = fs_tmp_path()) const;

	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet(const GatingHierarchy & gh_template,const GatingSet & cs, bool execute = true);
	/**
	 * assign the flow data from the source gs
	 * @param gs typically it is a root-only GatingSet that only carries cytoFrames
	 */
	void set_cytoset(const GatingSet & gs);

	/**
	 * Extract the ungated data
	 * @param node_path
	 * @return
	 */
	GatingSet get_cytoset();

	/**
	 * extract gated data
	 * @param node_path
	 * @return a root-only GatingSet that carries the subsetted cytoframes
	 */
	GatingSet get_cytoset(string node_path);
	string generate_h5_folder(string h5_dir) const;

	/**
	 * Retrieve the GatingHierarchy object from GatingSet by sample name.
	 *
	 * @param sampleName a string providing the sample name as the key
	 * @return a pointer to the GatingHierarchy object
	 */
	GatingHierarchyPtr getGatingHierarchy(string sample_uid) const;
	CytoFrameView & get_cytoframe_view_ref(string sample_uid){return getGatingHierarchy(sample_uid)->get_cytoframe_view_ref();}
	CytoFrameView get_cytoframe_view(string sample_uid) const{return getGatingHierarchy(sample_uid)->get_cytoframe_view();}
	/**
	 * insert an empty GatingHierarchy
	 * @param sn
	 */
	GatingHierarchyPtr add_GatingHierarchy(string sample_uid){
		GatingHierarchyPtr gh(new GatingHierarchy());
		add_GatingHierarchy(gh, sample_uid);
		return gh;
	}
	GatingHierarchyPtr add_GatingHierarchy(GatingHierarchyPtr gh, string sample_uid);

	 /**
	  * forward to the first element's getChannels
	  */
	vector<string> get_channels() const{return get_first_gh()->get_channels();};
	/**
	 * modify the channels for each individual fh
	 * @param _old
	 * @param _new
	 */
	void set_channel(const string & _old, const string & _new){
		for(auto & p : ghs_)
			p.second->set_channel(_old, _new);
	};

	/**
	 *
	 * update channel information stored in GatingSet
	 * @param chnl_map the mapping between the old and new channel names
	 */
	void set_channels(const CHANNEL_MAP & chnl_map);
//* forward to the first element's getChannels
	vector<string> get_markers(){return get_first_gh()->get_markers();};

	void set_marker(const string & _channel, const string & _marker){
		for(auto & p : ghs_)
			p.second->set_marker(_channel, _marker);
	};

	int n_cols(){return get_first_gh()->n_cols();}

	/**
	 * Subset by samples
	 * @param sample_uids
	 * @return
	 */
	GatingSet sub_samples(const vector<string> & sample_uids) const;

	/**
	 * Subset by samples (in place)
	 * @param sample_uids
	 * @return
	 */
	void sub_samples_(const vector<string> & sample_uids);
	/**
	 * Subet set by columns (in place)
	 * @param colnames
	 * @param col_type
	 * @return
	 */
	void cols_(vector<string> colnames, ColType col_type);


	/**
	 * Add sample
	 * @param sample_uid
	 * @param frame_ptr
	 */
	CytoFrameView & add_cytoframe_view(string sample_uid, const CytoFrameView & frame_view);

	/**
	 * update sample (move)
	 * @param sample_uid
	 * @param frame_ptr
	 */
	void update_cytoframe_view(string sample_uid, const CytoFrameView & frame_view);

	/**
	 * Constructor from FCS files
	 * @param file_paths
	 * @param config
	 * @param is_h5
	 * @param h5_dir
	 * @param is_add_root whether add root node. When false, gs is used as cytoset without gating tree
	 */
	GatingSet(const vector<string> & file_paths, const FCS_READ_PARAM & config= FCS_READ_PARAM()
			, bool is_h5 = true, string h5_dir = fs_tmp_path());

	GatingSet(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config = FCS_READ_PARAM()
			, bool is_h5 = true, string h5_dir = fs_tmp_path()):GatingSet()
	{
		add_fcs(sample_uid_vs_file_path, config, is_h5, h5_dir);
	}

	void add_fcs(const vector<pair<string,string>> & sample_uid_vs_file_path
			, const FCS_READ_PARAM & config, bool is_h5, string h5_dir, bool readonly = false);
	/**
	 * Update sample id
	 * @param _old
	 * @param _new
	 */
	void set_sample_uid(const string & _old, const string & _new);

	vector<string> get_sample_uids() const{
		return sample_names_;

	};

};


/**
 * validity checks on the frame to see if its data structure is consistent with cytoset
 *
 * @param frm
 * @return
 */
//	int isNotValidFrame(const CytoFrame & frm){
//		//validity check the channels against the existing frms
//		if(nCol() != frm.nCol())
//			return -1;
//
//		//check channel in linear time(taking advantage of the hash map)
//		auto frm1 = begin()->second;
//		for(const auto & c : frm.getChannels())
//		{
//			if(frm1.getColId(c, ColType::channel) <0 )
//				return -2;
//		}
//
//		//check the pdata
//		const auto & pd1 = frm1.getPData();
//		const auto & pd2 = frm.getPData();
//		if(pd1.size()!=pd2.size())
//			return -3;
//
//		for(const auto & p : pd2)
//		{
//			if(pd1.find(p.first)==pd1.end())
//				return -4;
//		}
//		return 0;
//	}


};


#endif /* GATINGSET_HPP_ */
