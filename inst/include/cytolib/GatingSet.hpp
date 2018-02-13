/*
 * GatingSet.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: wjiang2
 */



#ifndef GATINGSET_HPP_
#define GATINGSET_HPP_
#include "GatingHierarchy.hpp"
#include <string>
#include <cytolib/delimitedMessage.hpp>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

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
 *  It also stores the global transformations.
 *
 */
class GatingSet{

	biexpTrans globalBiExpTrans; //default bi-exponential transformation functions
	linTrans globalLinTrans;
	trans_global_vec gTrans;//parsed from xml workspace

	typedef unordered_map<string,GatingHierarchy> gh_map;
	gh_map ghs;
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
	vector<string> getChannels(){return begin()->second.getChannels();};
	/**
	 * modify the channels for each individual frame
	 * @param _old
	 * @param _new
	 */
	void setChannel(const string & _old, const string & _new){
		for(auto & p : ghs)
			p.second.setChannel(_old, _new);
	};

	//* forward to the first element's getChannels
	vector<string> getMarkers(){return begin()->second.getMarkers();};

	void setMarker(const string & _old, const string & _new){
		for(auto & p : ghs)
			p.second.setMarker(_old, _new);
	};

	int nCol(){return begin()->second.nCol();}

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
	 /**
	  * insert
	  * @param sampleName
	  * @return
	  */
	 GatingHierarchy & operator [](const string & sample_uid){
			 return ghs[sample_uid];
	   }

	string get_gatingset_id(){return guid_;}
	void set_gatingset_id(const string & guid){guid_ = guid;}
	/**
	 * iterate through hash map to extract sample names
	 * @return
	 */
	vector<string> get_sample_uids(){
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
	};

	~GatingSet()
	{
		if(g_loglevel>=GATING_SET_LEVEL)
			PRINT("\nstart to free GatingSet...\n");


		for(trans_global_vec::iterator it=gTrans.begin();it!=gTrans.end();it++)
		{
			trans_map curTrans=it->getTransMap();
			if(g_loglevel>=GATING_SET_LEVEL)
				PRINT("\nstart to free transformatioin group:"+it->getGroupName()+"\n");
			for(trans_map::iterator it1=curTrans.begin();it1!=curTrans.end();it1++)
			{
				transformation * curTran=it1->second;
				if(curTran!=NULL)
				{
					if(g_loglevel>=GATING_SET_LEVEL)
							PRINT("free transformatioin:"+curTran->getChannel()+"\n");

					delete curTran;
					curTran = NULL;
				}

			}

		}

	}

	GatingSet(){
		guid_ = generate_guid(10);
	};

	/**
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param filename
	 */
	void serialize_pb(string filename)
	{
		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;
		//init the output stream
		ofstream output(filename.c_str(), ios::out | ios::trunc | ios::binary);
		google::protobuf::io::OstreamOutputStream raw_output(&output);

		//empty message for gs
		pb::GatingSet gs_pb;

		/*
		 * add messages for trans
		 */

		//save the address of global biexp (as 1st entry)
		pb::TRANS_TBL * trans_tbl_pb = gs_pb.add_trans_tbl();
		intptr_t address = (intptr_t)&globalBiExpTrans;
		trans_tbl_pb->set_trans_address(address);


		//save the address of global lintrans(as 2nd entry)
		trans_tbl_pb = gs_pb.add_trans_tbl();
		address = (intptr_t)&globalLinTrans;
		trans_tbl_pb->set_trans_address(address);


		// cp trans group
		BOOST_FOREACH(trans_global_vec::value_type & it, gTrans){
			pb::trans_local * tg = gs_pb.add_gtrans();
			it.convertToPb(*tg, gs_pb);
		}

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
					GatingHierarchy & gh =  it.second;

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
	GatingSet(string filename)
	{
		GOOGLE_PROTOBUF_VERIFY_VERSION;
		ifstream input(filename.c_str(), ios::in | ios::binary);
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
			map<intptr_t, transformation *> trans_tbl;

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
					trans_tbl[old_address] = &globalBiExpTrans;
					break;
				case 1:
					trans_tbl[old_address] = &globalLinTrans;
					break;
				default:
					{
						switch(trans_pb.trans_type())
						{
						case pb::PB_CALTBL:
							trans_tbl[old_address] = new transformation(trans_pb);
							break;
						case pb::PB_BIEXP:
							trans_tbl[old_address] = new biexpTrans(trans_pb);
							break;
						case pb::PB_FASIGNH:
							trans_tbl[old_address] = new fasinhTrans(trans_pb);
							break;
						case pb::PB_FLIN:
							trans_tbl[old_address] = new flinTrans(trans_pb);
							break;
						case pb::PB_LIN:
							trans_tbl[old_address] = new linTrans(trans_pb);
							break;
						case pb::PB_LOG:
							trans_tbl[old_address] = new logTrans(trans_pb);
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

			for(int i = 0; i < pbGS.gtrans_size(); i++){
				const pb::trans_local & trans_local_pb = pbGS.gtrans(i);
				gTrans.push_back(trans_global(trans_local_pb, trans_tbl));
			}

			//read gating hierarchy messages
			for(int i = 0; i < pbGS.samplename_size(); i++){
				string sn = pbGS.samplename(i);
				//gh message is stored as the same order as sample name vector in gs
				pb::GatingHierarchy gh_pb;
				bool success = readDelimitedFrom(raw_input, gh_pb);

				if (!success) {
					throw(domain_error("Failed to parse GatingHierarchy."));
				}


				ghs[sn] = GatingHierarchy(gh_pb, trans_tbl);
			}
		}

	}

	/**
	 * Retrieve the GatingHierarchy object from GatingSet by sample name.
	 *
	 * @param sampleName a string providing the sample name as the key
	 * @return a pointer to the GatingHierarchy object
	 */
	GatingHierarchy & getGatingHierarchy(string sample_uid)
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
	GatingSet* clone(vector<string> sample_uids){
		GatingSet * gs=new GatingSet();

		//deep copying trans_global_vec
		trans_global_vec new_gtrans;
		map<transformation *, transformation * > trans_tbl;

		for(auto & it : gTrans){
			//copy trans_global
			trans_global newTransGroup;

			newTransGroup.setGroupName(it.getGroupName());
			newTransGroup.setSampleIDs(it.getSampleIDs());

			//clone trans_map
			trans_map newTmap;
			for(auto & v : it.getTransMap()){
				transformation * orig = v.second;
				transformation * trans = orig->clone();
				newTmap[v.first] = trans;
				//maintain a map between orig and new trans
				trans_tbl[orig] =  trans;
			}
			//save map into trans_global
			newTransGroup.setTransMap(newTmap);
			//save new trans_global into trans_global_vec
			new_gtrans.push_back(newTransGroup);
		}
		gs->set_gTrans(new_gtrans);

		//copy gh

		for(auto & it : sample_uids)
		{
			string curSampleName = it;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+curSampleName+"... \n");
			//except trans, the entire gh object is copiable
			GatingHierarchy curGh = getGatingHierarchy(curSampleName);

			//update trans_local with new trans pointers
			trans_map newTmap = curGh.getLocalTrans().getTransMap();
			for(auto & v : newTmap){
				transformation * orig = v.second;
				auto trans_it = trans_tbl.find(orig);
				if(trans_it == trans_tbl.end())
					throw(logic_error("the current transformation is not found in global trans table!"));

				transformation * trans = trans_it->second;
				v.second = trans;
			}
			curGh.addTransMap(newTmap);

			gs->ghs[curSampleName]=curGh;//add to the map

		}

		return gs;
	}
	/*
	 *
	 * copy gating trees from self to the dest gs
	 */
	void add(GatingSet & gs,vector<string> sample_uids){


		/*
		 * copy trans from gh_template into gtrans
		 * involve deep copying of transformation pointers
		 */
	//	if(g_loglevel>=GATING_SET_LEVEL)
	//		PRINT("copying transformation from gh_template...\n");
	//	trans_global newTransGroup;

	//	trans_map newTmap=gh_template->getLocalTrans().cloneTransMap();
	//	newTransGroup.setTransMap(newTmap);
	//	gTrans.push_back(newTransGroup);

		/*
		 * use newTmap for all other new ghs
		 */
		vector<string>::iterator it;
		for(it=sample_uids.begin();it!=sample_uids.end();it++)
		{
			string curSampleName=*it;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+curSampleName+"...\n ");


			GatingHierarchy &toCopy=gs.getGatingHierarchy(curSampleName);

			ghs[curSampleName]=toCopy.clone();


		}
	}
	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet(const GatingHierarchy & gh_template,vector<string> sample_uids){



		/*
		 * copy trans from gh_template into gtrans
		 * involve deep copying of transformation pointers
		 */
		if(g_loglevel>=GATING_SET_LEVEL)
			PRINT("copying transformation from gh_template...\n");
		trans_global newTransGroup;
	//	gh_template->printLocalTrans();
		trans_map newTmap=gh_template.getLocalTrans().cloneTransMap();
		newTransGroup.setTransMap(newTmap);
		gTrans.push_back(newTransGroup);

		/*
		 * use newTmap for all other new ghs
		 */
		vector<string>::iterator it;
		for(it=sample_uids.begin();it!=sample_uids.end();it++)
		{
			string curSampleName=*it;
			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start cloning GatingHierarchy for: "+curSampleName+"... \n");


			ghs[curSampleName] = gh_template.clone(newTmap,&gTrans);
	//		curGh->setNcPtr(&nc);//make sure to assign the global nc pointer to  GatingHierarchy

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("Gating hierarchy cloned: "+curSampleName+"\n");
		}
	}

	GatingSet(vector<pair<string,string>> sample_uid_vs_file_path, const FCS_READ_PARAM & config, string h5_dir):GatingSet()
	{
		fs::path h5_path = generate_h5_folder(fs::path(h5_dir));
		for(const auto & it : sample_uid_vs_file_path)
		{

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start adding GatingHierarchy for: "+ it.first +"... \n");


			GatingHierarchy & gh=addGatingHierarchy(it.first);
			gh.addRoot();//add default root

			string h5_filename = (h5_path/it.first).string() + ".h5";
			MemCytoFrame fr(it.second,config);
			fr.read_fcs();
			fr.writeH5(h5_filename);
			gh.set_cytoframe(h5_filename);

		}
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


	/*
	 * It is used by R API to add global transformation object during the auto gating.
	 * @param tName transformation name (usually channel name)
	 * @param tm trans_map
	 */
	void addTransMap(string gName,trans_map tm){
		/*
		 * assumption is there is either none transformation group exists before adding
		 */
		if(gTrans.size() == 0){

			trans_global thisTransGroup = trans_global();
			thisTransGroup.setGroupName(gName);
			thisTransGroup.setTransMap(tm);
			gTrans.push_back(thisTransGroup);
		}
		else
			throw(domain_error("transformation group already exists!Can't add the second one."));

	}
	/**
	 * insert an empty GatingHierarchy
	 * @param sn
	 */
	GatingHierarchy & addGatingHierarchy(string sample_uid){
		if(ghs.find(sample_uid)!=ghs.end())
			throw(domain_error("Can't add new GatingHierarchy since it already exists for: " + sample_uid));
		return ghs[sample_uid];
	}
	void addGatingHierarchy(const GatingHierarchy & gh, string sample_uid){
			if(ghs.find(sample_uid)!=ghs.end())
				throw(domain_error("Can't add new GatingHierarchy since it already exists for: " + sample_uid));
			ghs[sample_uid] = gh;
	}

	/**
	 *
	 * update channel information stored in GatingSet
	 * @param chnl_map the mapping between the old and new channel names
	 */
	void updateChannels(const CHANNEL_MAP & chnl_map)
	{
		//update trans
		for(trans_global_vec::iterator it=gTrans.begin();it!=gTrans.end();it++)
		{

			if(g_loglevel>=GATING_SET_LEVEL)
				PRINT("\nupdate channels for transformation group:" + it->getGroupName()+ "\n");
			it->updateChannels(chnl_map);

		}

		//update gh
		for(auto & it : ghs){

				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\nupdate channels for GatingHierarchy:"+it.first+"\n");
				it.second.updateChannels(chnl_map);
				//comp
			}
	}

	void set_gTrans(const trans_global_vec & _gTrans){gTrans = _gTrans;};
	biexpTrans * get_globalBiExpTrans(){return &globalBiExpTrans;};
	linTrans * get_globalLinTrans(){return &globalLinTrans;};

};

};


#endif /* GATINGSET_HPP_ */
