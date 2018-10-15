/*
 * GatingSet.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: wjiang2
 */



#ifndef GATINGSET_HPP_
#define GATINGSET_HPP_
#include "GatingHierarchy.hpp"
#include <cytolib/H5CytoFrame.hpp>
#include <cytolib/MemCytoFrame.hpp>
#include <cytolib/CytoFrameView.hpp>
#include <string>
#include <cytolib/delimitedMessage.hpp>

using namespace std;

namespace cytolib
{
#define ARCHIVE_TYPE_BINARY 0
#define ARCHIVE_TYPE_TEXT 1
#define ARCHIVE_TYPE_XML 2
#define PB true
#define BS false

#define GS_ID_LEN 20
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
	GatingHierarchyPtr get_first_gh() const
	{
		if(size() == 0)
			throw(range_error("Empty GatingSet!"));
		return begin()->second;
	}

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
		uid_ = generate_uid(GS_ID_LEN);
	};
	bool is_cytoFrame_only() const{return size() == 0||get_first_gh()->is_cytoFrame_only();};
	/**
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param filename
	 */
	void serialize_pb(string path, H5Option h5_opt)
	{
		/*
		 * validity check for path
		 */
		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		if(fs::exists(path))
		{
			if(!fs::is_empty(fs::path(path)))
			{
//				if(is_overwrite)
//				{
//					fs::remove_all(path);
//					fs::create_directory(path);
//				}
//				else
//				{
					fs::path pb_file;
					unordered_set<string> h5_samples;
					for(auto & e : fs::directory_iterator(path))
					{
						fs::path p = e.path();
						string ext = p.extension();
						if(ext == ".pb")
						{
							if(pb_file.empty())
								pb_file = p;
							else
							  throw(domain_error(errmsg + "Multiple .pb files found!"));
						}
						else if(ext == ".h5")
						{
							string sample_uid = p.stem();
							if(find(sample_uid) == end())
							  throw(domain_error(errmsg + "h5 file not matched to any sample in GatingSet: " + p.string()));
							else
								h5_samples.insert(p.stem());
						}
						else
						  throw(domain_error(errmsg + "File not recognized: " + p.string()));

					}

					if(!pb_file.empty())
						if(pb_file.stem() != uid_)
							throw(domain_error(errmsg + "The pb file doesn't match to the uid of GatingSet!"));
					for(const auto & it : ghs_)
					{
						if(h5_samples.find(it.first) == h5_samples.end())
							throw(domain_error(errmsg + "h5 file missing for sample: " + it.first));
					}
//				}
			}
		}
		else
			fs::create_directory(path);

		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;
		//init the output stream
		string filename = (fs::path(path) / uid_).string() + ".pb";
		ofstream output(filename.c_str(), ios::out | ios::trunc | ios::binary);
		google::protobuf::io::OstreamOutputStream raw_output(&output);

		//empty message for gs
		pb::GatingSet gs_pb;

		//uid
		gs_pb.set_guid(uid_);

		/*
		 *save sn and gh
		 *gh message is stored in the same order as sample_names_ vector
		 */
		for(auto & sn : sample_names_)
		{
			gs_pb.add_samplename(sn);
		}
		//write gs message to stream
		bool success = writeDelimitedTo(gs_pb, raw_output);

		if (!success){
			google::protobuf::ShutdownProtobufLibrary();
			throw(domain_error("Failed to write GatingSet."));
		}

		//write each gh as a separate message to stream due to the pb message size limit
		for(auto & sn : sample_names_)
		{
			string h5_filename = (fs::path(path) / sn).string() + ".h5";

			pb::GatingHierarchy pb_gh;
			getGatingHierarchy(sn)->convertToPb(pb_gh, h5_filename, h5_opt);


			bool success = writeDelimitedTo(pb_gh, raw_output);
			if (!success)
				throw(domain_error("Failed to write GatingHierarchy."));
		}

				// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}

	/**
	 * constructor from the archives (de-serialization)
	 * @param path
	 * @param is_skip_data whether to skip loading cytoframe data from h5. It should typically remain as default unless for debug purpose (e.g. legacy pb archive)
	 */
	GatingSet(string path, bool is_skip_data = false)
	{
		fs::path pb_file;
		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		for(auto & e : fs::directory_iterator(path))
		{
			fs::path p = e.path();
			string ext = p.extension();
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
		ifstream input(pb_file.c_str(), ios::in | ios::binary);
		if (!input) {
			throw(invalid_argument("File not found.." ));
		} else{
			 pb::GatingSet pbGS;

			 google::protobuf::io::IstreamInputStream raw_input(&input);
			 //read gs message
			 bool success = readDelimitedFrom(raw_input, pbGS);

			if (!success) {
				throw(domain_error("Failed to parse GatingSet."));
			}

			uid_ = pbGS.guid();


			//read gating hierarchy messages
			for(int i = 0; i < pbGS.samplename_size(); i++){
				string sn = pbGS.samplename(i);
				sample_names_.push_back(sn);
				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}

				pb::CytoFrame fr = *gh_pb.mutable_frame();
				string h5_filename = fs::path(path) / (sn + ".h5");

				ghs_[sn].reset(new GatingHierarchy(gh_pb, h5_filename, is_skip_data));
			}

		}

	}
	/**
	 * constructor from the legacy archives (de-serialization)
	 * @param filename
	 * @param format
	 * @param isPB
	 */
	GatingSet(string pb_file, const GatingSet & gs_data):GatingSet()
	{
		GOOGLE_PROTOBUF_VERIFY_VERSION;
		ifstream input(pb_file.c_str(), ios::in | ios::binary);
		if (!input) {
			throw(invalid_argument("File not found.." ));
		} else{
			 pb::GatingSet pbGS;

			 google::protobuf::io::IstreamInputStream raw_input(&input);
			 //read gs message
			 bool success = readDelimitedFrom(raw_input, pbGS);

			if (!success) {
				throw(domain_error("Failed to parse GatingSet."));
			}

			//parse global trans tbl from message
			map<intptr_t, TransPtr> trans_tbl;

			for(int i = 0; i < pbGS.trans_tbl_size(); i++){
				const pb::TRANS_TBL & trans_tbl_pb = pbGS.trans_tbl(i);
				const pb::transformation & trans_pb = trans_tbl_pb.trans();
				intptr_t old_address = (intptr_t)trans_tbl_pb.trans_address();

				/*
				 * first two global trans do not need to be restored from archive
				 * since they use the default parameters
				 * simply add the new address
				 */

				switch(i)
				{
				case 0:
					trans_tbl[old_address] = TransPtr(new biexpTrans());
					break;
				case 1:
					trans_tbl[old_address] = TransPtr(new linTrans());
					break;
				default:
					{
						switch(trans_pb.trans_type())
						{
						case pb::PB_CALTBL:
							trans_tbl[old_address] = TransPtr(new transformation(trans_pb));
							break;
						case pb::PB_BIEXP:
							trans_tbl[old_address] = TransPtr(new biexpTrans(trans_pb));
							break;
						case pb::PB_FASIGNH:
							trans_tbl[old_address] = TransPtr(new fasinhTrans(trans_pb));
							break;
						case pb::PB_FLIN:
							trans_tbl[old_address] = TransPtr(new flinTrans(trans_pb));
							break;
						case pb::PB_LIN:
							trans_tbl[old_address] = TransPtr(new linTrans(trans_pb));
							break;
						case pb::PB_LOG:
							trans_tbl[old_address] = TransPtr(new logTrans(trans_pb));
							break;
	//					case pb::PB_SCALE:
	//						trans_tbl[old_address] = new scaleTrans(trans_pb);
	//						break;
						default:
							throw(domain_error("unknown type of transformation archive!"));
						}
					}
				}

			}
			/*
			 * recover the trans_global
			 */

//			for(int i = 0; i < pbGS.gtrans_size(); i++){
//				const pb::trans_local & trans_local_pb = pbGS.gtrans(i);
//				gTrans.push_back(trans_global(trans_local_pb, trans_tbl));
//			}
			//read gating hierarchy messages
			for(int i = 0; i < pbGS.samplename_size(); i++){
				string sn = pbGS.samplename(i);
				sample_names_.push_back(sn);

				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}


				ghs_[sn].reset(new GatingHierarchy(gh_pb, trans_tbl));
			}

			if(gs_data.size()>0)
			{
				set_cytoset(gs_data);
				sample_names_ = gs_data.get_sample_uids();//update sample view from the data
			}

		}

	}


	/*
	 * up to caller to free the memory
	 */
	GatingSet copy(bool is_copy_data = true, bool is_realize_data = true, const string & new_h5_dir = fs::temp_directory_path().string()){
		GatingSet gs;
		gs.sample_names_ = sample_names_;
		fs::path h5_dir = gs.generate_h5_folder(fs::path(new_h5_dir));
		for(auto & it : ghs_)
		{
			string sn = it.first;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+sn+"... \n");
			gs.ghs_[sn] = it.second->copy(is_copy_data, is_realize_data, h5_dir/sn);

		}

		return gs;
	}

	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet(const GatingHierarchy & gh_template,vector<string> sample_uids):GatingSet(){

		vector<string>::iterator it;
		for(it=sample_uids.begin();it!=sample_uids.end();it++)
		{
			string sn=*it;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start cloning GatingHierarchy for: "+sn+"... \n");


			ghs_[sn] = gh_template.copy(false, false, "");

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("Gating hierarchy cloned: "+sn+"\n");
		}
		sample_names_ = sample_uids;
	}

	/**
	 * assign the flow data from the source gs
	 * @param gs typically it is a root-only GatingSet that only carries cytoFrames
	 */
	void set_cytoset(const GatingSet & gs){
		if(!gs.is_cytoFrame_only())
			throw(domain_error("The input gs is not data-only object! "));

		for(const string & sn : sample_names_)
		{
			const auto & it = gs.find(sn);
			if(it==gs.end())
				throw(domain_error("Sample '"  + sn + "' is missing from the data to be assigned!"));
			GatingHierarchy & gh = *ghs_[sn];
			gh.set_cytoFrame_view(it->second->get_cytoframe_view_ref());
		}
	};

	/**
	 * Extract the ungated data
	 * @param node_path
	 * @return
	 */
	GatingSet get_cytoset(){
		GatingSet gs;
		gs.sample_names_ = sample_names_;
		for(auto & it : ghs_)
		{
			string sn = it.first;
			GatingHierarchyPtr gh = it.second;

			gs.add_cytoframe_view(sn, gh->get_cytoframe_view_ref());
		}
		return gs;
	}

	/**
	 * extract gated data
	 * @param node_path
	 * @return a root-only GatingSet that carries the subsetted cytoframes
	 */
	GatingSet get_cytoset(string node_path){
		GatingSet gs;
		gs.sample_names_ = sample_names_;
		for(auto & it : ghs_)
		{
			//copy data
			string sn = it.first;
			GatingHierarchyPtr gh = it.second;
			CytoFrameView & fr = gs.add_cytoframe_view(sn, gh->get_cytoframe_view());
			//subset by node
			nodeProperties & node = gh->getNodeProperty(gh->getNodeID(node_path));
			fr.rows_(node.getIndices_u());
		}
		return gs;
	}
	fs::path generate_h5_folder(fs::path h5_dir)
	{
		h5_dir /= uid_;
		const string sh5_dir = h5_dir.string();
		if(fs::exists(h5_dir))
			throw(domain_error(sh5_dir + " already exists!"));
		if(!create_directories(h5_dir))
			throw(domain_error("Failed to create directory: " + sh5_dir));
		return h5_dir;
	}
	/**
	 * Retrieve the GatingHierarchy object from GatingSet by sample name.
	 *
	 * @param sampleName a string providing the sample name as the key
	 * @return a pointer to the GatingHierarchy object
	 */
	GatingHierarchyPtr getGatingHierarchy(string sample_uid) const
	{

		const_iterator it=ghs_.find(sample_uid);
		if(it==ghs_.end())
			throw(domain_error(sample_uid + " not found!"));
		else
			return it->second;
	}
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
	GatingHierarchyPtr add_GatingHierarchy(GatingHierarchyPtr gh, string sample_uid){
			if(ghs_.find(sample_uid)!=ghs_.end())
				throw(domain_error("Can't add new GatingHierarchy since it already exists for: " + sample_uid));
			ghs_[sample_uid] = gh;
			sample_names_.push_back(sample_uid);
			return ghs_[sample_uid];
	}

	 /**
	  * forward to the first element's getChannels
	  */
	vector<string> get_channels(){return get_first_gh()->get_channels();};
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
	void set_channels(const CHANNEL_MAP & chnl_map)
	{

		//update gh
		for(auto & it : ghs_){

				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\nupdate channels for GatingHierarchy:"+it.first+"\n");
				it.second->set_channels(chnl_map);
				//comp
			}

	}
//* forward to the first element's getChannels
	vector<string> get_markers(){return get_first_gh()->get_markers();};

	void set_marker(const string & _old, const string & _new){
		for(auto & p : ghs_)
			p.second->set_marker(_old, _new);
	};

	int n_cols(){return get_first_gh()->n_cols();}

	/**
	 * Subset by samples
	 * @param sample_uids
	 * @return
	 */
	GatingSet sub_samples(const vector<string> & sample_uids) const
	{
		GatingSet res(*this);
		res.sub_samples_(sample_uids);
		res.uid_ = generate_uid(GS_ID_LEN);
		return res;
	}

	/**
	 * Subset by samples (in place)
	 * @param sample_uids
	 * @return
	 */
	void sub_samples_(const vector<string> & sample_uids)
	{
		ghMap ghs_new;
		//validity check
		for(const auto & uid : sample_uids)
		{
			const auto & it = find(uid);
			if(it==end())
				throw(domain_error("The data to be assigned is missing sample: " + uid));
			else
				ghs_new[uid] = it->second;
		}
		sample_names_ = sample_uids;
		ghs_ = ghs_new;
	}

	/**
	 * Subet set by columns (in place)
	 * @param colnames
	 * @param col_type
	 * @return
	 */
	void cols_(vector<string> colnames, ColType col_type)
	{

		for(auto & it : ghs_)
		{
			if(!it.second->is_cytoFrame_only())
				throw(domain_error("Can't subset by cols when gh is not data-only object! "));
			it.second->get_cytoframe_view_ref().cols_(colnames, col_type);

		}
	}


	/**
	 * Add sample
	 * @param sample_uid
	 * @param frame_ptr
	 */
	CytoFrameView & add_cytoframe_view(string sample_uid, const CytoFrameView & frame_view){
		if(!is_cytoFrame_only())
			throw(domain_error("Can't add cytoframes to gs when it is not data-only object! "));
		GatingHierarchyPtr gh = add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(frame_view)), sample_uid);
		return gh->get_cytoframe_view_ref();

	}

	/**
	 * update sample (move)
	 * @param sample_uid
	 * @param frame_ptr
	 */
	void update_cytoframe_view(string sample_uid, const CytoFrameView & frame_view){
		if(find(sample_uid) == end())
			throw(domain_error("Can't update the cytoframe since it doesn't exists: " + sample_uid));
		ghs_[sample_uid].reset(new GatingHierarchy(frame_view));
	}

	/**
	 * Constructor from FCS files
	 * @param file_paths
	 * @param config
	 * @param is_h5
	 * @param h5_dir
	 */
	GatingSet(const vector<string> & file_paths, const FCS_READ_PARAM & config, bool is_h5, string h5_dir):GatingSet()
	{
		vector<pair<string,string>> map(file_paths.size());
		transform(file_paths.begin(), file_paths.end(), map.begin(), [](string i){return make_pair(path_base_name(i), i);});
		add_fcs(map, config, is_h5, h5_dir);


	}

	GatingSet(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config, bool is_h5, string h5_dir):GatingSet()
	{
		add_fcs(sample_uid_vs_file_path, config, is_h5, h5_dir);
	}

	void add_fcs(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config, bool is_h5, string h5_dir)
	{

		fs::path h5_path(h5_dir);
		for(const auto & it : sample_uid_vs_file_path)
		{

			string h5_filename = (h5_path/it.first).string() + ".h5";
			CytoFramePtr fr_ptr(new MemCytoFrame(it.second,config));
			//set pdata
			fr_ptr->set_pheno_data("name", path_base_name(it.second));

			dynamic_cast<MemCytoFrame&>(*fr_ptr).read_fcs();
			if(is_h5)
			{
				fr_ptr->write_h5(h5_filename);
				fr_ptr.reset(new H5CytoFrame(h5_filename));
			}

			add_cytoframe_view(it.first, CytoFrameView(fr_ptr));

		}
	}

	/**
	 * Update sample id
	 * @param _old
	 * @param _new
	 */
	void set_sample_uid(const string & _old, const string & _new){
		if(_old.compare(_new) != 0)
		{
			auto it = find(_new);
			if(it!=end())
				throw(range_error(_new + " already exists!"));
			it = find(_old);
			if(it==end())
				throw(range_error(_old + " not found!"));

			ghs_[_new] = it->second;
			erase(_old);
		}

	};

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
