/*
 * CytoSet.hpp
 *
 *  Created on: Feb 16, 2018
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_CYTOSET_HPP_
#define INST_INCLUDE_CYTOLIB_CYTOSET_HPP_
#include <cytolib/H5CytoFrame.hpp>

namespace cytolib
{
	class CytoSet
	{
		typedef unordered_map<string,H5CytoFrame> FrameMap;
		FrameMap frames_;
	public:
		typedef typename FrameMap::iterator iterator;
		typedef typename FrameMap::const_iterator const_iterator;

		/*
		 * forwarding APIs
		 */
		 size_t size(){return frames_.size();}
		 iterator end(){return frames_.end();}
		 iterator begin(){return frames_.begin();}
		 iterator find(const string &sample_uid){
				 return frames_.find(sample_uid);
		 }
		 size_t erase ( const string& k ){return frames_.erase(k);}


	 /**
		  * forward to the first element's getChannels
		  */
		vector<string> get_channels(){return begin()->second.get_channels();};
		/**
		 * modify the channels for each individual frame
		 * @param _old
		 * @param _new
		 */
		void set_channel(const string & _old, const string & _new){
			for(auto & p : frames_)
				p.second.set_channel(_old, _new);
		};

		//* forward to the first element's getChannels
		vector<string> get_markers(){return begin()->second.get_markers();};

		void set_marker(const string & _old, const string & _new){
			for(auto & p : frames_)
				p.second.set_marker(_old, _new);
		};

		int n_cols(){return begin()->second.n_cols();}

		string get_h5_file_path(){return path_dir_name(begin()->second.get_h5_file_path());}

		H5CytoFrame & get_cytoframe(string sample_uid)
		{

			iterator it=find(sample_uid);
			if(it==end())
				throw(domain_error(sample_uid + " not found in gating set!"));
			else
				return it->second;
		}
		void add_cytoframe(string sample_uid, const H5CytoFrame & frame){
			if(find(sample_uid) != end())
				throw(domain_error("Can't add new cytoframe since it already exists for: " + sample_uid));
			frames_[sample_uid] = frame;
		}

		CytoSet(){}

		CytoSet(vector<pair<string,string>> sample_uid_vs_file_path, const FCS_READ_PARAM & config, string h5_dir)
		{
			fs::path h5_path(h5_dir);
			for(const auto & it : sample_uid_vs_file_path)
			{

				string h5_filename = (h5_path/it.first).string() + ".h5";
				MemCytoFrame fr(it.second,config);
				//set pdata
				fr.set_pheno_data("name", path_base_name(it.second));
				fr.read_fcs();
				fr.write_h5(h5_filename);
				add_cytoframe(it.first, H5CytoFrame(h5_filename));

			}
		}

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
				frames_[_new] = gh;
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
					it.second.update_channels(chnl_map);
					//comp
				}

		}
	};
}




#endif /* INST_INCLUDE_CYTOLIB_CYTOSET_HPP_ */
