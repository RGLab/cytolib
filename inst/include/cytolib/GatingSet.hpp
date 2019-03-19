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
		uid_ = generate_uid();
	};
	bool is_cytoFrame_only() const{return size() == 0||get_first_gh()->is_cytoFrame_only();};
	/**
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param is_skip_data whether to skip writing cytoframe data to pb. It is typically remain as default unless for debug purpose (e.g. re-writing gs that is loaded from legacy pb archive without actual data associated)
	 */
	void serialize_pb(string path, H5Option h5_opt, bool is_skip_data = false)
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
		string buf;
		google::protobuf::io::StringOutputStream raw_output(&buf);

		//empty message for gs
		pb::GatingSet gs_pb;

		//uid
		gs_pb.set_guid(uid_);

		/*
		 *save sn and gh
		 *gh message is stored in the same order as sample_names_ vector
		 */
		const vector<string> sample_names = get_sample_uids();
		for(auto & sn : sample_names)
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
		for(auto & sn : sample_names)
		{
			string h5_filename = (fs::path(path) / sn).string() + ".h5";

			pb::GatingHierarchy pb_gh;
			getGatingHierarchy(sn)->convertToPb(pb_gh, h5_filename, h5_opt, is_skip_data);


			bool success = writeDelimitedTo(pb_gh, raw_output);
			if (!success)
				throw(domain_error("Failed to write GatingHierarchy."));
		}

				// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
		output.write(&buf[0], buf.size());
	}

	/**
	 * constructor from the archives (de-serialization)
	 * @param path
	 * @param is_skip_data whether to skip loading cytoframe data from h5. It should typically remain as default unless for debug purpose (e.g. legacy pb archive)
	 */
	GatingSet(string path, bool is_skip_data = false, unsigned int h5_flags = H5F_ACC_RDONLY)
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

			uid_ = pbGS.guid();


			//read gating hierarchy messages
			for(int i = 0; i < pbGS.samplename_size(); i++){
				string sn = pbGS.samplename(i);

				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}

				pb::CytoFrame fr = *gh_pb.mutable_frame();
				string h5_filename = fs::path(path) / (sn + ".h5");

				add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(gh_pb, h5_filename, is_skip_data, h5_flags)), sn);
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
						case pb::PB_LOGICLE:
							trans_tbl[old_address] = TransPtr(new logicleTrans(trans_pb));
							break;
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

				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}
				//only add the sample that is present in gs_data(in case fs was subsetted when gs was archived in legacy pb)
				if(gs_data.find(sn)!=gs_data.end())
					add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(gh_pb, trans_tbl)), sn);
			}

			if(gs_data.size()>0)
			{
				set_cytoset(gs_data);
				/*
				 * update sample view from the data, this direct operation on sample_names_is typically not recommended
				 * and herer only done for legacy gs where fs was separately stored and its sample order can be different from pb
				 */
				sample_names_ = gs_data.get_sample_uids();
			}

		}

	}


	/*
	 * up to caller to free the memory
	 */
	GatingSet copy(bool is_copy_data = true, bool is_realize_data = true, const string & new_h5_dir = fs::temp_directory_path().string()){
		GatingSet gs;
		fs::path h5_dir;
		if(is_copy_data)
			h5_dir = gs.generate_h5_folder(fs::path(new_h5_dir));
		for(const string & sn : get_sample_uids())
		{
			GatingHierarchyPtr gh = getGatingHierarchy(sn);

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+sn+"... \n");
			gs.add_GatingHierarchy(gh->copy(is_copy_data, is_realize_data, h5_dir/sn), sn);

		}

		return gs;
	}

	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet(const GatingHierarchy & gh_template,const GatingSet & cs, unsigned int flags = H5F_ACC_RDWR):GatingSet(){

		fs::path h5_dir = generate_h5_folder(fs::temp_directory_path());

		for(const string & sn : cs.get_sample_uids())
		{

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start cloning GatingHierarchy for: "+sn+"... \n");
			GatingHierarchyPtr gh = add_GatingHierarchy(gh_template.copy(false, false, ""), sn);
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... load flow data: "+sn+"... \n");
			MemCytoFrame fr = MemCytoFrame(*(cs.get_cytoframe_view(sn).get_cytoframe_ptr()));
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... compensate: "+sn+"... \n");
			gh->compensate(fr);
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... transform_data: "+sn+"... \n");
			gh->transform_data(fr);
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... gating: "+sn+"... \n");
			gh->gating(fr, 0, true, true);
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... save flow data: "+sn+"... \n");
			string h5_filename = (h5_dir/sn).string() + ".h5";
			fr.write_h5(h5_filename);
			//attach to gh
			gh->set_cytoframe_view(CytoFrameView(CytoFramePtr(new H5CytoFrame(h5_filename, flags))));
		}

	}

	/**
	 * assign the flow data from the source gs
	 * @param gs typically it is a root-only GatingSet that only carries cytoFrames
	 */
	void set_cytoset(const GatingSet & gs){
		if(!gs.is_cytoFrame_only())
			throw(domain_error("The input gs is not data-only object! "));

		for(const string & sn : get_sample_uids())
		{
			const auto & it = gs.find(sn);
			if(it==gs.end())
				throw(domain_error("Sample '"  + sn + "' is missing from the data to be assigned!"));
			GatingHierarchy & gh = *ghs_[sn];
			gh.set_cytoframe_view(it->second->get_cytoframe_view_ref());
		}
	};

	/**
	 * Extract the ungated data
	 * @param node_path
	 * @return
	 */
	GatingSet get_cytoset(){
		GatingSet gs;

		for(const string & sn : get_sample_uids())
		{
			GatingHierarchyPtr gh = getGatingHierarchy(sn);

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
		for(const string & sn : get_sample_uids())
		{
			GatingHierarchyPtr gh = getGatingHierarchy(sn);

			CytoFrameView & fr = gs.add_cytoframe_view(sn, gh->get_cytoframe_view());
			//subset by node
			auto u = gh->getNodeID(node_path);
			gh->check_ungated_bool_node(u);
			nodeProperties & node = gh->getNodeProperty(u);
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
		res.uid_ = generate_uid();
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
	 * @param is_add_root whether add root node. When false, gs is used as cytoset without gating tree
	 */
	GatingSet(const vector<string> & file_paths, const FCS_READ_PARAM & config= FCS_READ_PARAM(), bool is_h5 = true, string h5_dir = fs::temp_directory_path()):GatingSet()
	{
		vector<pair<string,string>> map(file_paths.size());
		transform(file_paths.begin(), file_paths.end(), map.begin(), [](string i){return make_pair(path_base_name(i), i);});
		add_fcs(map, config, is_h5, h5_dir);
	}

	GatingSet(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config = FCS_READ_PARAM(), bool is_h5 = true, string h5_dir = fs::temp_directory_path()):GatingSet()
	{
		add_fcs(sample_uid_vs_file_path, config, is_h5, h5_dir);
	}

	void add_fcs(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config, bool is_h5, string h5_dir, unsigned int flags = H5F_ACC_RDWR)
	{

		fs::path h5_path;
		if(is_h5)
			h5_path = generate_h5_folder(h5_dir);
		for(const auto & it : sample_uid_vs_file_path)
		{


			CytoFramePtr fr_ptr(new MemCytoFrame(it.second,config));
			//set pdata
			fr_ptr->set_pheno_data("name", path_base_name(it.second));

			dynamic_cast<MemCytoFrame&>(*fr_ptr).read_fcs();
			if(is_h5)
			{
				string h5_filename = (h5_path/it.first).string() + ".h5";
				fr_ptr->write_h5(h5_filename);
				fr_ptr.reset(new H5CytoFrame(h5_filename, flags));
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

			//update sample view
			auto it1 = std::find(sample_names_.begin(), sample_names_.end(),_new);
			if(it1 != sample_names_.end())
				throw(range_error(_new + " already exists!"));
			it1 = std::find(sample_names_.begin(), sample_names_.end(),_old);
			if(it1==sample_names_.end())
				throw(range_error(_old + " not found!"));
			*it1 = _new;

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
