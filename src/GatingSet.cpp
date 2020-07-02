// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/GatingSet.hpp>
#include <cytolib/H5CytoFrame.hpp>
#include <cytolib/MemCytoFrame.hpp>
#include <cytolib/cytolibConfig.h>


namespace cytolib
{
	GatingHierarchyPtr GatingSet::get_first_gh() const
	{
		if(size() == 0)
			throw(range_error("Empty GatingSet!"));
		return begin()->second;
	}
	GatingSet::GatingSet(string path, bool is_skip_data
			, bool readonly
			, vector<string> select_samples
			, bool print_lib_ver
			, CytoCtx ctx):ctx_(ctx)
	{

		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		fs::path gs_pb_file;
		unordered_set<string> cf_samples;
		unordered_set<string> pb_samples;
		FileFormat fmt;
		bool is_first_sample = true;
		//search for h5
		CytoVFS vfs(ctx);
		for(auto & e : vfs.ls(path))
		{
			fs::path p(e);
			string ext = p.extension().string();
			string fn = p.stem().string();
			if(ext == ".h5"||ext == ".tile")
			{
				//init fmt for the first sample file
				if(is_first_sample)
				{
					if(ext == ".h5")
						fmt = FileFormat::H5;
					else
						fmt = FileFormat::TILE;
					is_first_sample = false;
				}
				else
				{
					//consistency check
					if((fmt == FileFormat::H5&&ext==".tile")||(fmt == FileFormat::TILE&&ext==".h5"))
						throw(domain_error(errmsg + "Multiple file formats found!"));
				}


				cf_samples.insert(fn);
			}
			else if(ext == ".pb")
			{
				pb_samples.insert(fn);

			}
			else if(ext == ".gs")
			{
				if(gs_pb_file.empty())
					gs_pb_file = p;
				else
					throw(domain_error(errmsg + "Multiple .gs files found for the same gs object!"));
			}
			else
				throw(domain_error(errmsg + "File not recognized: " + p.string()));
		}

		bool is_legacy = false;
		if(gs_pb_file.empty())
		{
			if(pb_samples.size()==1)
			{
				auto id = *(pb_samples.begin());//check if pb file seems like a guid of a legacy gs
				if(cf_samples.find(id)==cf_samples.end())
				{
					is_legacy = true;
				}
				else
				  throw(domain_error(errmsg + "No .gs file found!"));
			}
			else
			  throw(domain_error(errmsg + "No .gs file found!"));
		}

		if(is_legacy)
		{
			cout << path + " seems to be the legacy archive and it is recommended to convert to the new format by saving it to the new folder!" << endl;
			deserialize_legacy(path, is_skip_data, readonly, select_samples, print_lib_ver);
		}
		else
		{

			for(auto sn : cf_samples)
			{
				if(pb_samples.find(sn)==pb_samples.end())
					  throw(domain_error(errmsg + "No .pb file matched for sample " + sn));
			}

			for(auto sn : pb_samples)
			{
				if(cf_samples.find(sn)==cf_samples.end())
					  throw(domain_error(errmsg + "No cytoframe file matched for sample " + sn + ".pb"));
			}

			GOOGLE_PROTOBUF_VERIFY_VERSION;


			auto buf = vfs.read_buf(gs_pb_file.string());
			 pb::GatingSet pbGS;
			 //read entire file into buffer since message-lite doesn't support iostream



			 google::protobuf::io::ArrayInputStream raw_input(buf.data(), buf.size());
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
					if(is_remote_path(path))
						PRINT("loading GatingHierarchy: " + sn + " \n");
					string gh_pb_file = (fs::path(path) / sn).string() + ".pb";

					pb::GatingHierarchy gh_pb;
					auto buf = vfs.read_buf(gh_pb_file);

					 google::protobuf::io::ArrayInputStream raw_input(buf.data(), buf.size());
					 //read gs message
					 bool success = readDelimitedFrom(raw_input, gh_pb);

					if (!success) {
						throw(domain_error("Failed to parse GatingHierarchy " + sn));
					}

					pb::CytoFrame fr = *gh_pb.mutable_frame();

					string uri;
					auto cf_ext = "." + fmt_to_str(fmt);
					uri = (fs::path(path) / (sn + cf_ext)).string();
					add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(ctx_, gh_pb, uri, is_skip_data, readonly)), sn);

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
	 * separate filename from dir to avoid to deal with path parsing in c++
	 * @param path the dir of filename
	 * @param is_skip_data whether to skip writing cytoframe data to pb. It is typically remain as default unless for debug purpose (e.g. re-writing gs that is loaded from legacy pb archive without actual data associated)
	 */
	void GatingSet::serialize_pb(string path
			, CytoFileOption cf_opt
			, bool is_skip_data
			, const CytoCtx & ctx)
	{
		/*
		 * validity check for path
		 */
		if(is_remote_path(path))
		{
			if(begin()->second->get_cytoframe_view_ref().get_backend_type()!=FileFormat::TILE)
				throw(logic_error("Only tiledb backend supports saving to remote path!"));
		}

		string errmsg = "Not a valid GatingSet archiving folder! " + path + "\n";
		CytoVFS vfs(ctx);
		if(vfs.is_dir(path))
		{
			auto files = vfs.ls(path);
			if(files.size()>0)
			{
				fs::path gs_pb_file;
				unordered_set<string> cf_samples;
				unordered_set<string> pb_samples;
				for(auto & e : files)
				{
					fs::path p(e);
					string ext = p.extension().string();
					if(ext == ".gs")
					{
						string fn = p.stem().string();
						if(fn == uid_)
						{
							if(gs_pb_file.empty())
								gs_pb_file = p;
							else
								throw(domain_error(errmsg + "Multiple .gs files found for the same gs object!"));
						}
						else
						{
							  throw(domain_error(errmsg + "gs file not matched to GatingSet uid: " + p.string()));

						}

					}
					else if(ext == ".tile"||ext == ".h5"||ext == ".pb")
					{
						string sample_uid = p.stem().string();
						if(find(sample_uid) == end())
						  throw(domain_error(errmsg + "file not matched to any sample in GatingSet: " + p.string()));
						else
						{
							if(ext == ".pb")
								pb_samples.insert(p.stem().string());
							else
								cf_samples.insert(p.stem().string());
						}

					}
					else
					  throw(domain_error(errmsg + "File not recognized: " + p.string()));

				}

				if(gs_pb_file.empty())
				{
					if(pb_samples.size()>0)
					{
						throw(domain_error(errmsg + "gs file missing"));
					}
				}

				for(const auto & it : ghs_)
				{
					auto sn = it.first;
					if(cf_samples.find(sn) == cf_samples.end())
						throw(domain_error(errmsg + "cytoframe file missing for sample: " + sn));
					/*
					 * when no pb present, treat it as valid dest folder where h5 were pre-written
					 * e.g. when save_gs with cdf = skip
					 */
					if(pb_samples.size()>0)
					{
						if(pb_samples.find(sn) == pb_samples.end())
							throw(domain_error(errmsg + "pb file missing for sample: " + sn));
					}
				}
			}
		}
		else
			vfs.create_dir(path);

		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;

//		ofstream output(filename.c_str(), ios::out | ios::trunc | ios::binary);
		string buf;
		google::protobuf::io::StringOutputStream raw_output(&buf);

		//empty message for gs
		pb::GatingSet gs_pb;

		//ver
		gs_pb.set_cytolib_verion(CYTOLIB_VERSION);
		gs_pb.set_pb_verion(std::to_string(GOOGLE_PROTOBUF_VERSION));
		unsigned h5major,h5minor,h5rel;
		H5::H5Library::getLibVersion(h5major,h5minor,h5rel);
		string h5ver = std::to_string(h5major) + "." + std::to_string(h5minor) + "." +std::to_string(h5rel);
		gs_pb.set_h5_verion(h5ver);
		//uid
		gs_pb.set_guid(uid_);

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
		//init the output stream for gs
		string gs_pb_file = (fs::path(path) / uid_).string() + ".gs";

		vfs.write_buf(gs_pb_file, buf);
		//write each gh as a separate message to stream due to the pb message size limit
		//we now go one step further to save each message to individual file
		//due to the single string buffer used by lite-message won't be enough to hold the all samples for large dataset
		for(auto & sn : sample_names)
		{
			auto gh = getGatingHierarchy(sn);
			auto src_uri = gh->get_cytoframe_view_ref().get_uri();
			if(is_remote_path(path)||is_remote_path(src_uri))
				PRINT("saving GatingHierarchy: " + sn + " \n");
			string cf_filename = (fs::path(path) / sn).string();
			string buf;
			google::protobuf::io::StringOutputStream raw_output(&buf);


			pb::GatingHierarchy pb_gh;
			gh->convertToPb(pb_gh, cf_filename, cf_opt, is_skip_data, ctx);


			bool success = writeDelimitedTo(pb_gh, raw_output);
			if (!success)
			{
				google::protobuf::ShutdownProtobufLibrary();
				throw(domain_error("Failed to write GatingHierarchy."));
			}
			//init the output stream for gs
			string gh_pb_file = (fs::path(path) / sn).string() + ".pb";


			vfs.write_buf(gh_pb_file, buf);


		}

				// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}

	 /* constructor from the legacy archives (de-serialization)
	 * @param filename
	 * @param format
	 * @param isPB
	 */
	GatingSet::GatingSet(string pb_file, const GatingSet & gs_data):GatingSet()
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
						case pb::PB_SCALE:
							trans_tbl[old_address] = TransPtr(new scaleTrans(trans_pb));
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
					add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(gh_pb, trans_tbl)), sn, false);
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
	GatingSet GatingSet::copy(bool is_copy_data, bool is_realize_data, const string & new_cf_dir) const{
		GatingSet gs;
		fs::path cf_dir;
		if(is_copy_data)
			cf_dir = gs.generate_cytoframe_folder(fs::path(new_cf_dir).string());
		for(const string & sn : get_sample_uids())
		{
			GatingHierarchyPtr gh = getGatingHierarchy(sn);

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... copying GatingHierarchy: "+sn+"... \n");
			gs.add_GatingHierarchy(gh->copy(is_copy_data, is_realize_data, (cf_dir/sn).string()), sn, is_copy_data);

		}

		return gs;
	}

	/*Defunct
	 * TODO:current version of this contructor is based on gating template ,simply copying
	 * compensation and transformation,more options can be allowed in future like providing different
	 * comp and trans
	 */
	GatingSet::GatingSet(const GatingHierarchy & gh_template,const GatingSet & cs, bool execute):GatingSet(){
		auto samples = cs.get_sample_uids();
		for(const string & sn : samples)
		{

			if(g_loglevel>=GATING_HIERARCHY_LEVEL)
				PRINT("\n... start cloning GatingHierarchy for: "+sn+"... \n");
			auto gh = gh_template.copy(false, false, "");
			auto cfv = cs.get_cytoframe_view(sn);
			string cf_filename = cfv.get_uri();
			if(cf_filename=="")
				throw(logic_error("in-memory version of cs is not supported!"));
			if(execute)
			{
				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\n... load flow data: "+sn+"... \n");
				MemCytoFrame fr = MemCytoFrame(*(cfv.get_cytoframe_ptr()));
				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\n... compensate: "+sn+"... \n");
				gh->compensate(fr);
				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\n... transform_data: "+sn+"... \n");
				// fr.scale_time_channel();
				gh->transform_data(fr);
				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\n... gating: "+sn+"... \n");
				gh->gating(fr, 0, true, true);
				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\n... save flow data: "+sn+"... \n");
				cfv.set_params(fr.get_params());
				cfv.set_keywords(fr.get_keywords());
				cfv.set_data(fr.get_data());
			}
			//attach to gh
			gh->set_cytoframe_view(cfv);
			add_GatingHierarchy(gh, sn, false);

		}

	}

	/**
	 * assign the flow data from the source gs
	 * @param gs typically it is a root-only GatingSet that only carries cytoFrames
	 */
	void GatingSet::set_cytoset(const GatingSet & gs){
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
	GatingSet GatingSet::get_cytoset(){
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
	GatingSet GatingSet::get_cytoset(string node_path){
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

	/**
	 * Retrieve the GatingHierarchy object from GatingSet by sample name.
	 *
	 * @param sampleName a string providing the sample name as the key
	 * @return a pointer to the GatingHierarchy object
	 */
	GatingHierarchyPtr GatingSet::getGatingHierarchy(string sample_uid) const
	{

		const_iterator it=ghs_.find(sample_uid);
		if(it==ghs_.end())
			throw(domain_error(sample_uid + " not found!"));
		else
			return it->second;
	}

	/**
	 *
	 * update channel information stored in GatingSet
	 * @param chnl_map the mapping between the old and new channel names
	 */
	void GatingSet::set_channels(const CHANNEL_MAP & chnl_map)
	{

		//update gh
		for(auto & it : ghs_){

				if(g_loglevel>=GATING_HIERARCHY_LEVEL)
					PRINT("\nupdate channels for GatingHierarchy:"+it.first+"\n");
				it.second->set_channels(chnl_map);
				//comp
			}

	}

	/**
	 * Subset by samples
	 * @param sample_uids
	 * @return
	 */
	GatingSet GatingSet::sub_samples(const vector<string> & sample_uids) const
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
	void GatingSet::sub_samples_(const vector<string> & sample_uids)
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
	void GatingSet::cols_(vector<string> colnames, ColType col_type)
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
	CytoFrameView & GatingSet::add_cytoframe_view(string sample_uid, const CytoFrameView & frame_view){
		if(!is_cytoFrame_only())
			throw(domain_error("Can't add cytoframes to gs when it is not data-only object! "));
		//validity check
		auto res = channel_consistency_check<GatingSet, CytoFrameView>(*this, frame_view, sample_uid);
		GatingHierarchyPtr gh = add_GatingHierarchy(GatingHierarchyPtr(new GatingHierarchy(res)), sample_uid);
		return gh->get_cytoframe_view_ref();

	}

	/**
	 * update sample (move)
	 * @param sample_uid
	 * @param frame_ptr
	 */
	void GatingSet::update_cytoframe_view(string sample_uid, const CytoFrameView & frame_view){
		if(find(sample_uid) == end())
			throw(domain_error("Can't update the cytoframe since it doesn't exists: " + sample_uid));
		auto res = channel_consistency_check<GatingSet, CytoFrameView>(*this, frame_view, sample_uid);
		ghs_[sample_uid]->set_cytoframe_view(res);
	}


	/**
	 * Update sample id
	 * @param _old
	 * @param _new
	 */
	void GatingSet::set_sample_uid(const string & _old, const string & _new){
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


};


