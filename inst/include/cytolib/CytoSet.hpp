/*
 * CytoSet.hpp
 *
 *  Created on: Feb 16, 2018
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOSET_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOSET_HPP_
#include <cytolib/H5CytoFrame.hpp>
#include <cytolib/MemCytoFrame.hpp>
namespace cytolib
{

	class CytoSet
	{
		typedef unordered_map<string,CytoFramePtr> FrameMap;
		FrameMap frames_;

		CytoFramePtr & get_first_frame_ptr()
		{
			if(size() == 0)
				throw(range_error("Empty CytoSet!"));
			return begin()->second;
		}
	public:
		typedef typename FrameMap::iterator iterator;
		typedef typename FrameMap::const_iterator const_iterator;

		/*
		 * forwarding APIs
		 */
		 size_t size() const{return frames_.size();}
		 iterator end(){return frames_.end();}
		 iterator begin(){return frames_.begin();}
		 const_iterator end() const{return frames_.end();}
		 const_iterator begin() const{return frames_.begin();}
		 iterator find(const string &sample_uid){
				 return frames_.find(sample_uid);
		 }
		 size_t erase ( const string& k ){return frames_.erase(k);}

		 CytoSet(const pb::CytoSet & cs_pb, const string & path)
		 {
			 fs::path dir(path);
			 for(int i = 0; i < cs_pb.samplename_size(); i++)
			 {
				 string uid = cs_pb.samplename(i);
				 const pb::CytoFrame & fr = cs_pb.frame(i);
				 string h5_filename = dir / (uid + ".h5");
				 if(fs::exists(h5_filename))
				 {
					 frames_[uid].reset(new H5CytoFrame(h5_filename));
					 const CytoFrame & ref = *(frames_[uid]);
					 if(!fr.is_h5())
						 frames_[uid].reset(new MemCytoFrame(ref));
				 }
				 else
					 throw(domain_error("H5 file missing for sample: " + uid));

			 }
		 }

		 void convertToPb(pb::CytoSet & cs_pb, const string & path, H5Option h5_opt)
		 {
			 for(const auto & it : frames_)
			 {
				 cs_pb.add_samplename(it.first);
				 pb::CytoFrame * fr_pb = cs_pb.add_frame();
				 string h5_filename = (fs::path(path) / it.first).string() + ".h5";
				 it.second->convertToPb(*fr_pb, h5_filename, h5_opt);

			 }
		 }
		 /**
		  * forward to the first element's getChannels
		  */
		vector<string> get_channels(){return get_first_frame_ptr()->get_channels();};
		/**
		 * modify the channels for each individual frame
		 * @param _old
		 * @param _new
		 */
		void set_channel(const string & _old, const string & _new){
			for(auto & p : frames_)
				p.second->set_channel(_old, _new);
		};

		//* forward to the first element's getChannels
		vector<string> get_markers(){return get_first_frame_ptr()->get_markers();};

		void set_marker(const string & _old, const string & _new){
			for(auto & p : frames_)
				p.second->set_marker(_old, _new);
		};

		int n_cols(){return get_first_frame_ptr()->n_cols();}

		/**
		 * Subset by samples
		 * @param sample_uids
		 * @return
		 */
		CytoSet sub_samples(vector<string> sample_uids)
		{
			CytoSet res;
			for(const auto & uid : sample_uids)
				res.frames_[uid] = this->get_cytoframe(uid);
			return res;
		}

		/**
		 * Subset by samples (in place)
		 * @param sample_uids
		 * @return
		 */
		void sub_samples_(vector<string> sample_uids)
		{
			FrameMap tmp;

			for(const auto & uid : sample_uids)
				tmp[uid] = this->get_cytoframe(uid);
			swap(tmp, frames_);
		}

		/**
		 * Subet set by columns (in place)
		 * @param colnames
		 * @param col_type
		 * @return
		 */
		void cols_(vector<string> colnames, ColType col_type)
		{

			for(auto & it : frames_)
				it.second->cols_(colnames, col_type);
		}

		/**
		 * Extract the single CytoFrame
		 * @param sample_uid
		 * @return
		 */
		CytoFramePtr get_cytoframe(string sample_uid)
		{

			iterator it=find(sample_uid);
			if(it==end())
				throw(domain_error(sample_uid + " not found in gating set!"));
			else
				return it->second;
		}

		/**
		 * Add sample
		 * @param sample_uid
		 * @param frame_ptr
		 */
		void add_cytoframe(string sample_uid, const CytoFramePtr & frame_ptr){
			if(find(sample_uid) != end())
				throw(domain_error("Can't add new cytoframe since it already exists for: " + sample_uid));
			frames_[sample_uid] = frame_ptr;
		}

		/**
		 * update sample (move)
		 * @param sample_uid
		 * @param frame_ptr
		 */
		void update_cytoframe(string sample_uid, const CytoFramePtr & frame_ptr){
			if(find(sample_uid) == end())
				throw(domain_error("Can't update the cytoframe since it doesn't exists: " + sample_uid));
			frames_[sample_uid] = frame_ptr;
		}

		CytoSet(){}

		/**
		 *
		 * @param cs
		 */
		CytoSet copy(const string & new_h5_dir = "")
		{
			CytoSet cs;
			for(const auto & it : frames_)
			{
				string new_h5 = "";
				if(new_h5_dir != "")
				{
					new_h5 = fs::path(new_h5_dir) / (it.first + ".h5");
				}
				cs.frames_[it.first] = it.second->copy(new_h5);
			}

			return cs;
		}
		/**
		 * Constructor from FCS files
		 * @param file_paths
		 * @param config
		 * @param is_h5
		 * @param h5_dir
		 */
		CytoSet(const vector<string> & file_paths, const FCS_READ_PARAM & config, bool is_h5, string h5_dir)
		{
			vector<pair<string,string>> map(file_paths.size());
			transform(file_paths.begin(), file_paths.end(), map.begin(), [](string i){return make_pair(path_base_name(i), i);});
			add_fcs(map, config, is_h5, h5_dir);


		}

		CytoSet(const vector<pair<string,string>> & sample_uid_vs_file_path, const FCS_READ_PARAM & config, bool is_h5, string h5_dir)
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

				add_cytoframe(it.first, fr_ptr);

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

				frames_[_new] = it->second;
				erase(_old);
			}

		};

		vector<string> get_sample_uids() const{
			vector<string> res;
			for(const auto & f : frames_)
				res.push_back(f.first);
			return res;

		};

		void update_channels(const CHANNEL_MAP & chnl_map){
			//update gh
			for(auto & it : frames_){
					it.second->update_channels(chnl_map);
					//comp
				}

		}
	};
}




#endif /* INST_INCLUDE_CYTOLIB_CYTOSET_HPP_ */
