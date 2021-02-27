// Copyright 2019 Fred Hutchinson Cancer Research Center
// See the included LICENSE file for details on the licence that is granted to the user of this software.
#include <cytolib/compensation.hpp>
#include <boost/foreach.hpp>

namespace cytolib
{
	compensation::compensation(mat spillMat, vector<string> _markers){
		if(spillMat.n_cols!=spillMat.n_rows)
			throw(domain_error("Non-square shaped spillover matrix only possible if both marker names and detector names are provided!"));
		if(spillMat.n_cols!=_markers.size())
			throw(domain_error("The number of columns of the spillover matrix is not the same as the markers!"));

		marker = _markers;
		detector = _markers;
		mat B = spillMat.t();

		spillOver.resize(spillMat.n_elem);
		memcpy(spillOver.data(), B.mem, sizeof(double)*spillMat.n_elem);
	}

	compensation::compensation(mat spillMat, vector<string> _markers, vector<string> _detectors){
		if(spillMat.n_cols > spillMat.n_rows)
			throw(domain_error("Invalid spillover matrix (less channels than labels)!"));
		if(spillMat.n_cols!=_markers.size())
			throw(domain_error("The number of columns of the spillover matrix is not the same as the markers!"));

		marker = _markers;
		detector = _detectors;
		mat B = spillMat.t();

		spillOver.resize(spillMat.n_elem);
		memcpy(spillOver.data(), B.mem, sizeof(double)*spillMat.n_elem);
	}
	/**
	 * convert spillover matrix from row-majored std::vector
	 * to col-majored arma::ma
	 * @return
	 */
	mat compensation::get_spillover_mat ()const
	{

		unsigned nMarker = marker.size();
		unsigned nDetector = detector.size();
		if(spillOver.size()!=nDetector * nMarker)
			throw(domain_error("invalid spillover matrix!"));
		mat B(spillOver.data(), nDetector, nMarker);
		return B.t();
	}
	void compensation::update_channels(const CHANNEL_MAP & chnl_map){

		for(vector<string>::iterator it = marker.begin(); it != marker.end(); it++)
		{
			string curName = *it;

			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(curName);
			if(itChnl!=chnl_map.end())
				*it = itChnl->second;
		}
		for(vector<string>::iterator it = detector.begin(); it != detector.end(); it++)
		{
			string curName = *it;

			CHANNEL_MAP::const_iterator itChnl = chnl_map.find(curName);
			if(itChnl!=chnl_map.end())
				*it = itChnl->second;
		}
	}
	void compensation::convertToPb(pb::COMP & comp_pb){
		comp_pb.set_cid(cid);
		comp_pb.set_name(name);
		comp_pb.set_prefix(prefix);
		comp_pb.set_suffix(suffix);
		comp_pb.set_comment(comment);
		BOOST_FOREACH(vector<double>::value_type & it, spillOver){
			comp_pb.add_spillover(it);
		}
		BOOST_FOREACH(vector<string>::value_type & it, marker){
				comp_pb.add_marker(it);
			}
		BOOST_FOREACH(vector<string>::value_type & it, detector){
				comp_pb.add_detector(it);
			}
	}
	compensation::compensation(const pb::COMP & comp_pb):cid(comp_pb.cid()),prefix(comp_pb.prefix()),suffix(comp_pb.suffix()),name(comp_pb.name()),comment(comp_pb.comment()){
		// Newer pb files will have the detectors saved out
		if(comp_pb.detector_size()){
			for(int i = 0; i < comp_pb.marker_size(); i++)
				marker.push_back(comp_pb.marker(i));
			for(int i = 0; i < comp_pb.detector_size(); i++)
				detector.push_back(comp_pb.detector(i));
		}else{
			// Older pb files assume detectors == markers
			for(int i = 0; i < comp_pb.marker_size(); i++){
				marker.push_back(comp_pb.marker(i));
				detector.push_back(comp_pb.marker(i));
			}
		}

		for(int i = 0; i < comp_pb.spillover_size(); i++)
				spillOver.push_back(comp_pb.spillover(i));
	}


};

