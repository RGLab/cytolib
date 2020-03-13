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
	 * @param select_sample_idx samples to load
	 */

	GatingSet(string path, bool is_skip_data = false, bool readonly = true, vector<string> select_samples = {}, bool print_lib_ver = false)
	{

		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		fs::path gs_pb_file;
		unordered_set<string> h5_samples;
		unordered_set<string> pb_samples;
		//search for h5
		for(auto & e : fs::directory_iterator(path))
		{
			fs::path p = e;
			string ext = p.extension().string();
			string fn = p.stem().string();
			if(ext == ".h5")
			{
				h5_samples.insert(fn);
			}
			else if(ext != ".pb")
				throw(domain_error(errmsg + "File not recognized: " + p.string()));

		}
		//search for pb file
		for(auto & e : fs::directory_iterator(path))
		{
			fs::path p = e;
			string ext = p.extension().string();
			string fn = p.stem().string();
			if(ext == ".pb")
			{
				if(h5_samples.find(fn)==h5_samples.end())
				{
					if(gs_pb_file.empty())
						gs_pb_file = p;
					else
					{
						errmsg += " Can't determine the pb file for gs since both .pb files do not match to any h5 sample files!";
						errmsg += gs_pb_file.string() + ", " +  p.string();
						throw(domain_error(errmsg));

					}

				}
				else
					pb_samples.insert(fn);


			}

		}

		if(gs_pb_file.empty())
		  throw(domain_error(errmsg + "No .pb file found for gs!"));

		if(pb_samples.size()==0)
		{
			cout << path + " seems to be the legacy archive and it is recommended to convert to the new format by saving it to the new folder!" << endl;
			deserialize_legacy(path, is_skip_data, readonly, select_samples, print_lib_ver);
		}
		else
		{
			if(pb_samples.size()<h5_samples.size())
			{
				for(auto sn : h5_samples)
				{
					if(pb_samples.find(sn)==pb_samples.end())
						  throw(domain_error(errmsg + "No .pb file found for sample " + sn + ".h5"));
				}
			}


			GOOGLE_PROTOBUF_VERIFY_VERSION;
			ifstream input(gs_pb_file.string(), ios::in | ios::binary);
			if (!input) {
				throw(invalid_argument("File not found.." + gs_pb_file.string()));
			} else{
				 pb::GatingSet pbGS;
				 //read entire file into buffer since message-lite doesn't support iostream
				 input.seekg (0, input.end);
				 int length = input.tellg();
				 input.seekg (0, input.beg);
				 vector<char> buf(length);
				 input.read(buf.data(), length);

				 google::protobuf::io::ArrayInputStream raw_input(buf.data(), length);
				 //read gs message
				 bool success = readDelimitedFrom(raw_input, pbGS);

				if (!success) {
					throw(domain_error("Failed to parse GatingSet."));
				}
				if(print_lib_ver)
				{
					PRINT("The GatingSet was archived by:\n");
					PRINT("cytolib: ");
					string cv = pbGS.cytolib_verion();
					if(cv=="")
						cv = "unknown";
					PRINT(cv);
					PRINT("\n");
					PRINT("protobuf: ");
					string pv = pbGS.pb_verion();
					if(pv=="")
						pv = "unknown";
					PRINT(pv);
					PRINT("\n");
					PRINT("HDF5: ");
					string hv = pbGS.h5_verion();
					if(hv=="")
						hv = "unknown";
					PRINT(hv);
					PRINT("\n");

				}
				uid_ = pbGS.guid();

				auto nSelect = select_samples.size();
				auto nTotal = pbGS.samplename_size();
				unordered_map<string,bool> sn_hash;

				//prescan select and update hash
				if(nSelect>0)
				{
					for(int i = 0; i < nTotal; i++){
						string sn = pbGS.samplename(i);
						sn_hash[sn] = false;
					}
					for(unsigned i = 0; i < nSelect; i++)
					{
						auto sel = select_samples[i];
						auto it = sn_hash.find(sel);
						if(it == sn_hash.end())
							throw(domain_error("sample selection is out of boundary: " + sel));
						it->second = true;

					}
				}
				//read gating hierarchy messages
				for(int i = 0; i < nTotal; i++){
					string sn = pbGS.samplename(i);

					//conditional add gh
					if(nSelect==0||sn_hash.find(sn)->second)
					{
						string gh_pb_file = (fs::path(path) / sn).string() + ".pb";
						ifstream input(gh_pb_file, ios::in | ios::binary);
						if (!input) {
							throw(invalid_argument("File not found.." + gh_pb_file));
						}
						else
						{
							pb::GatingHierarchy gh_pb;
							 //read entire file into buffer since message-lite doesn't support iostream
							 input.seekg (0, input.end);
							 int length = input.tellg();
							 input.seekg (0, input.beg);
							 vector<char> buf(length);
							 input.read(buf.data(), length);

							 google::protobuf::io::ArrayInputStream raw_input(buf.data(), length);
							 //read gs message
							 bool success = readDelimitedFrom(raw_input, gh_pb);

							if (!success) {
								throw(domain_error("Failed to parse GatingHierarchy " + sn));
							}

						pb::CytoFrame fr = *gh_pb.mutable_frame();
						string h5_filename = (fs::path(path) / (sn + ".h5")).string();

						add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(gh_pb, h5_filename, is_skip_data, readonly)), sn);
						}
					}
				}


				//reorder view based on select
				if(nSelect>0)
				{
					uid_ = generate_uid();
					sample_names_ = select_samples;
				}
			}
		}
	}
	/**
	 * legacy de-serialization for single pb file
	 * @param path
	 * @param is_skip_data whether to skip loading cytoframe data from h5. It should typically remain as default unless for debug purpose (e.g. legacy pb archive)
	 */
	void deserialize_legacy(string path, bool is_skip_data = false, bool readonly = true, vector<string> select_samples = {}, bool print_lib_ver = false)
	{
		GatingSet gs;
		fs::path pb_file;
		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		for(auto & e : fs::directory_iterator(path))
		{
			fs::path p = e;
			string ext = p.extension().string();
			if(ext == ".pb")
			{
				if(pb_file.empty())
					pb_file = p;
				else
				  throw(domain_error(errmsg + "Multiple .pb files found!"));
			}
			else if(ext != ".h5")
				throw(domain_error(errmsg + "File not recognized: " + p.string()));

		}

		if(pb_file.empty())
		  throw(domain_error(errmsg + "No .pb file found!"));


		GOOGLE_PROTOBUF_VERIFY_VERSION;
		ifstream input(pb_file.string(), ios::in | ios::binary);
		if (!input) {
			throw(invalid_argument("File not found.." ));
		} else{
			 pb::GatingSet pbGS;
			 //read entire file into buffer since message-lite doesn't support iostream
			 input.seekg (0, input.end);
			 int length = input.tellg();
			 input.seekg (0, input.beg);
			 vector<char> buf(length);
			 input.read(buf.data(), length);

			 google::protobuf::io::ArrayInputStream raw_input(buf.data(), length);
			 //read gs message
			 bool success = readDelimitedFrom(raw_input, pbGS);

			if (!success) {
				throw(domain_error("Failed to parse GatingSet."));
			}
			if(print_lib_ver)
			{
				PRINT("The GatingSet was archived by:\n");
				PRINT("cytolib: ");
				string cv = pbGS.cytolib_verion();
				if(cv=="")
					cv = "unknown";
				PRINT(cv);
				PRINT("\n");
				PRINT("protobuf: ");
				string pv = pbGS.pb_verion();
				if(pv=="")
					pv = "unknown";
				PRINT(pv);
				PRINT("\n");
				PRINT("HDF5: ");
				string hv = pbGS.h5_verion();
				if(hv=="")
					hv = "unknown";
				PRINT(hv);
				PRINT("\n");

			}
			uid_ = pbGS.guid();

			auto nSelect = select_samples.size();
			auto nTotal = pbGS.samplename_size();
			unordered_map<string,bool> sn_hash;

			//prescan select and update hash
			if(nSelect>0)
			{
				for(int i = 0; i < nTotal; i++){
					string sn = pbGS.samplename(i);
					sn_hash[sn] = false;
				}
				for(unsigned i = 0; i < nSelect; i++)
				{
					auto sel = select_samples[i];
					auto it = sn_hash.find(sel);
					if(it == sn_hash.end())
						throw(domain_error("sample selection is out of boundary: " + sel));
					it->second = true;

				}
			}
			//read gating hierarchy messages
			for(int i = 0; i < nTotal; i++){
				string sn = pbGS.samplename(i);
//				all_samples[i] =sn;
				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}
				//conditional add gh (thus avoid to load h5)
				if(nSelect==0||sn_hash.find(sn)->second)
				{
					pb::CytoFrame fr = *gh_pb.mutable_frame();
					string h5_filename = (fs::path(path) / (sn + ".h5")).string();

					add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(gh_pb, h5_filename, is_skip_data, readonly)), sn);
				}
			}


			//reorder view based on select
			if(nSelect>0)
			{
				uid_ = generate_uid();
				sample_names_ = select_samples;
			}
		}

	}
	/**
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
	/**
	 *
	 * @param gh
	 * @param sample_uid
	 * @param validity_check only set it to false when adding the legacy gh that has no data associated
	 * @return
	 */
	GatingHierarchyPtr add_GatingHierarchy(GatingHierarchyPtr gh, string sample_uid, bool validity_check = true)
	{
			if(ghs_.find(sample_uid)!=ghs_.end())
				throw(domain_error("Can't add new sample since it already exists for: " + sample_uid));
			if(validity_check)
			{
				auto view = channel_consistency_check<GatingSet, CytoFrameView>(*this, gh->get_cytoframe_view(), sample_uid);
				gh->set_cytoframe_view(view);//update potentially reordered view
			}
			ghs_[sample_uid] = gh;
			sample_names_.push_back(sample_uid);
			return ghs_[sample_uid];
	}

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
