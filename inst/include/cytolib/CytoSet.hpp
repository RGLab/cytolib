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
		typedef unordered_map<string,unique_ptr<CytoFrame>> FrameMap;
		FrameMap frames_;

		unique_ptr<CytoFrame> & get_first_frame_ptr()
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
		CytoSet operator[](vector<string> sample_uids)
		{
			CytoSet res;
			for(const auto & uid : sample_uids)
				res.add_cytoframe(uid, this->get_cytoframe(uid));
			return res;
		}

		/**
		 * Subset by columns
		 * @param colnames
		 * @param col_type
		 * @return
		 */
		CytoSet cols(vector<string> colnames, ColType col_type)
		{
			CytoSet res = *this;
			res.cols_(colnames, col_type);
			return res;
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
		CytoFrame & get_cytoframe(string sample_uid)
		{

			iterator it=find(sample_uid);
			if(it==end())
				throw(domain_error(sample_uid + " not found in gating set!"));
			else
				return *(it->second);
		}

		/**
		 * Add sample (move)
		 * @param sample_uid
		 * @param frame_ptr
		 */
		void add_cytoframe(string sample_uid, unique_ptr<CytoFrame> && frame_ptr){
			if(find(sample_uid) != end())
				throw(domain_error("Can't add new cytoframe since it already exists for: " + sample_uid));
			swap(frames_[sample_uid], frame_ptr);
		}

		/**
		 * Add sample (copy)
		 * @param sample_uid
		 * @param frame
		 */
		void add_cytoframe(string sample_uid, const CytoFrame & frame){
			add_cytoframe(sample_uid, unique_ptr<CytoFrame>(frame.copy()));
		}
		/**
		 * update sample (move)
		 * @param sample_uid
		 * @param frame_ptr
		 */
		void update_cytoframe(string sample_uid, unique_ptr<CytoFrame> && frame_ptr){
			if(find(sample_uid) == end())
				throw(domain_error("Can't update the cytoframe since it doesn't exists: " + sample_uid));
			swap(frames_[sample_uid], frame_ptr);
		}

		/**
		 * update sample (copy)
		 * @param sample_uid
		 * @param frame
		 */
		void update_cytoframe(string sample_uid, const CytoFrame & frame){
			update_cytoframe(sample_uid, unique_ptr<CytoFrame>(frame.copy()));

		}

		CytoSet(){}

		/**
		 * Copy constructor (shallow copy)
		 * @param cs
		 */
		CytoSet(const CytoSet & cs)
		{
			for(const auto & it : cs)
				this->frames_[it.first] = unique_ptr<CytoFrame>(it.second->copy());
		}

		/**
		 * Assignment operator (shallow copy)
		 * @param cs
		 * @return
		 */
		CytoSet & operator=(const CytoSet & cs)
		{
			for(const auto & it : cs)
				this->frames_[it.first] = unique_ptr<CytoFrame>(it.second->copy());
			return *this;
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
				unique_ptr<CytoFrame> fr_ptr;
				unique_ptr<MemCytoFrame> fr(new MemCytoFrame(it.second,config));
				//set pdata
				fr->set_pheno_data("name", path_base_name(it.second));
				fr->read_fcs();
				if(is_h5)
				{
					fr->write_h5(h5_filename);
					fr_ptr.reset(new H5CytoFrame(h5_filename));
				}
				else
					fr_ptr = move(fr);

				add_cytoframe(it.first, move(fr_ptr));

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

				frames_[_new] =  move(it->second);
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
