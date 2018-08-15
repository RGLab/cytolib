/*
 * GatingSet.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: wjiang2
 */



#ifndef GATINGSET_HPP_
#define GATINGSET_HPP_
#include "GatingHierarchy.hpp"
#include "CytoSet.hpp"
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

	typedef unordered_map<string, GatingHierarchyPtr> gh_map;
	gh_map ghs;
	CytoSet cytoset_;
	string guid_;
public:
	typedef typename gh_map::iterator iterator;
	typedef typename gh_map::const_iterator const_iterator;

	/*
	 * forwarding APIs
	 */
	 size_t size(){return ghs.size();}
	 iterator end(){return ghs.end();}
	 iterator begin(){return ghs.begin();}
	 iterator find(const string &sample_uid){
			 return ghs.find(sample_uid);
	 }
	 size_t erase ( const string& k ){return ghs.erase(k);}

	 /**
	  * forward to the first element's getChannels
	  */
	vector<string> get_channels(){return cytoset_.get_channels();};
	/**
	 * modify the channels for each individual frame
	 * @param _old
	 * @param _new
	 */
	void set_channel(const string & _old, const string & _new){
		cytoset_.set_channel(_old, _new);
	};

	//* forward to the first element's getChannels
	vector<string> get_markers(){return cytoset_.get_markers();};

	void set_marker(const string & _old, const string & _new){
		cytoset_.set_marker(_old, _new);
	};

	int n_cols(){return cytoset_.n_cols();}

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
//todo:
//	CytoSet subset(const vector<string> & sampleNames);
//	vector<PDATA> getPData();
//	void setPData(const pData & pd);

	 GatingSet sub_samples(const vector<string> & sample_uids){
		 	 GatingSet gs;
		 	 if(cytoset_.size() > 0)//avoid subsetting cs when gs is not gated
		 		 gs.cytoset_ = cytoset_.sub_samples(sample_uids);
		 	 for(const auto & s :sample_uids)
		 	 {
		 		 auto it = find(s);
		 		 if(it == end())
		 			 throw(domain_error("Sample not found in GatingSet: " + s));

		 		gs.ghs[s] = it->second;
		 	 }

			 return gs;
	   }

	string get_gatingset_id(){return guid_;}
	void set_gatingset_id(const string & guid){guid_ = guid;}
	/**
	 * iterate through hash map to extract sample names
	 * @return
	 */
	vector<string> get_sample_uids() const{
		vector<string> res;
		for(const auto & f : ghs)
			res.push_back(f.first);
		return res;

	};
	/**
	 * modify the name of one sample , which involves delete/insert the existing frame
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

			auto gh = it->second;
			erase(_old);
			ghs[_new] = gh;
		}
		cytoset_.set_sample_uid(_old, _new);
	};

	GatingSet(){
		guid_ = generate_guid(20);
	};

	/**
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param filename
	 */
	void serialize_pb(string path, bool is_overwrite, H5Option h5_opt)
	{
		/*
		 * validity check for path
		 */
		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		if(fs::exists(path))
		{
			if(is_overwrite)
			{
				fs::remove_all(path);
				fs::create_directory(path);
			}
			else
			{
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

				if(pb_file.empty())
				  throw(domain_error(errmsg + "No .pb file found!"));
				else
					if(pb_file.stem() != guid_)
						throw(domain_error(errmsg + "The pb file doesn't match to the guid of GatingSet!"));
		        for(const auto & it : ghs)
		        {
		        	if(h5_samples.find(it.first) == h5_samples.end())
		        		throw(domain_error(errmsg + "h5 file missing for sample: " + it.first));
		        }
			}

		}
		else
			fs::create_directory(path);

		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;
		//init the output stream
		string filename = (fs::path(path) / guid_).string() + ".pb";
		ofstream output(filename.c_str(), ios::out | ios::trunc | ios::binary);
		google::protobuf::io::OstreamOutputStream raw_output(&output);

		//empty message for gs
		pb::GatingSet gs_pb;

		//guid
		gs_pb.set_guid(guid_);

		//cs
		pb::CytoSet * cs_pb = gs_pb.mutable_cs();
		cytoset_.convertToPb(*cs_pb, path, h5_opt);

		//add sample name
		for(auto & it : ghs){
				string sn = it.first;
				gs_pb.add_samplename(sn);
		}

		//write gs message to stream
		bool success = writeDelimitedTo(gs_pb, raw_output);

		if (!success){
			google::protobuf::ShutdownProtobufLibrary();
			throw(domain_error("Failed to write GatingSet."));
		}else
		{
			/*
			 * write pb message for each sample
			 */

			for(auto & it :ghs){
					string sn = it.first;
					GatingHierarchy & gh =  *(it.second);

					pb::GatingHierarchy pb_gh;
					gh.convertToPb(pb_gh);
					//write the message
					bool success = writeDelimitedTo(pb_gh, raw_output);
					if (!success)
						throw(domain_error("Failed to write GatingHierarchy."));
			}

		}

		// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}
	/**
	 * constructor from the archives (de-serialization)
	 * @param filename
	 * @param format
	 * @param isPB
	 */
	GatingSet(string path)
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

			guid_ = pbGS.guid();


			//read gating hierarchy messages
			for(int i = 0; i < pbGS.samplename_size(); i++){
				string sn = pbGS.samplename(i);
				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}


				ghs[sn].reset(new GatingHierarchy(gh_pb));
			}

			cytoset_ = CytoSet(pbGS.cs(), path);


		}

	}

	/**
	 * Retrieve the GatingHierarchy object from GatingSet by sample name.
	 *
	 * @param sampleName a string providing the sample name as the key
	 * @return a pointer to the GatingHierarchy object
	 */
	GatingHierarchyPtr getGatingHierarchy(string sample_uid)
	{

		iterator it=ghs.find(sample_uid);
		if(it==ghs.end())
			throw(domain_error(sample_uid + " not found in gating set!"));
		else
			return it->second;
	}


	/*
	 * up to caller to free the memory
	 */
	GatingSet copy(const string & new_h5_dir = ""){
		GatingSet gs;

		gs.cytoset_ = cytoset_.copy(new_h5_dir);

		for(auto & it : ghs)
		{
			string curSampleName = it.first;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+curSampleName+"... \n");
			gs.ghs[curSampleName] = it.second->copy();

		}

		return gs;
	}

	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet(const GatingHierarchy & gh_template,vector<string> sample_uids){

		vector<string>::iterator it;
		for(it=sample_uids.begin();it!=sample_uids.end();it++)
		{
			string curSampleName=*it;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start cloning GatingHierarchy for: "+curSampleName+"... \n");


			ghs[curSampleName] = gh_template.copy();

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("Gating hierarchy cloned: "+curSampleName+"\n");
		}
	}

	GatingSet(const CytoSet & cytoset):GatingSet()
	{
		for(const auto & it : cytoset.get_sample_uids())
		{

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start adding GatingHierarchy for: "+ it +"... \n");


			GatingHierarchyPtr gh=add_GatingHierarchy(it);
			gh->addRoot();//add default root


		}
		set_cytoset(cytoset);
	}

	void set_cytoset(const CytoSet & cytoset){cytoset_ = cytoset;};
	void set_cytoset(CytoSet && cytoset){swap(cytoset_ , cytoset);};
	const CytoSet & get_cytoset(){
		return cytoset_;
	}

	CytoSet get_cytoset(string node_path){
		CytoSet cs = cytoset_;
		for(auto & it : ghs)
		{
			GatingHierarchyPtr gh = it.second;
			nodeProperties & node = gh->getNodeProperty(gh->getNodeID(node_path));
			cs.get_cytoframe_view_ref(it.first).rows_(node.getIndices_u());
		}
		return cs;
	}
	fs::path generate_h5_folder(fs::path h5_dir)
	{
		h5_dir /= guid_;
		const string sh5_dir = h5_dir.string();
		if(fs::exists(h5_dir))
			throw(domain_error(sh5_dir + " already exists!"));
		if(!create_directories(h5_dir))
			throw(domain_error("Failed to create directory: " + sh5_dir));
		return h5_dir;
	}


	/**
	 * insert an empty GatingHierarchy
	 * @param sn
	 */
	GatingHierarchyPtr add_GatingHierarchy(string sample_uid){
		if(ghs.find(sample_uid)!=ghs.end())
			throw(domain_error("Can't add new GatingHierarchy since it already exists for: " + sample_uid));
		return ghs[sample_uid];
	}
	void add_GatingHierarchy(GatingHierarchyPtr gh, string sample_uid){
			if(ghs.find(sample_uid)!=ghs.end())
				throw(domain_error("Can't add new GatingHierarchy since it already exists for: " + sample_uid));
			ghs[sample_uid] = gh;
	}

	/**
	 *
	 * update channel information stored in GatingSet
	 * @param chnl_map the mapping between the old and new channel names
	 */
	void update_channels(const CHANNEL_MAP & chnl_map)
	{

		//update gh
		for(auto & it : ghs){

				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\nupdate channels for GatingHierarchy:"+it.first+"\n");
				it.second->update_channels(chnl_map);
				//comp
			}
		//update flow data
		cytoset_.update_channels(chnl_map);

	}


};

};


#endif /* GATINGSET_HPP_ */
