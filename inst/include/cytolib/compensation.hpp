/* Copyright 2019 Fred Hutchinson Cancer Research Center
 * See the included LICENSE file for details on the license that is granted to the
 * user of this software.
 * compensation.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: wjiang2
 */

#ifndef INCLUDE_COMPENSATION_HPP_
#define INCLUDE_COMPENSATION_HPP_
#include <armadillo>
using namespace arma;

#include "transformation.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <queue>
#include "global.hpp"
namespace cytolib
{
class compensation{
public:
	string cid;
	string prefix;
	string suffix;
	string name;
	string comment;// store "Acquisition-defined" when the spillOver matrix is not supplied and cid==-1
	vector<string> marker;
	vector<double> spillOver;
	compensation(){};
	/**
	 * constructor from the existing spillover matrix
	 * @spillover arma:mat representing a col-major matrix
	 * @_markers string vector
	 */
	compensation(mat spillMat, vector<string> _markers);
	/**
	 * construct spillover matrix from string value of FCS keyword
	 * @param val
	 * @return
	 */
	compensation (const string & val)
	{

		vector<string> valVec;
		boost::split(valVec, val, boost::is_any_of(","));
		int n = boost::lexical_cast<unsigned int>(valVec[0]);
		// Pre-check size of valVec -- this is strictly imposed by the FCS spec
		if(valVec.size() != 1+n+n*n){
			PRINT("$SPILLOVER keyword value is of improper size for number of spillover channels!\n");
			// This is just a simple check. If something other than . was used as decimal
			// separator but used consistently, the size should be 1+n+2*(n*n)
			if(valVec.size() == 1+n+2*n*n){
				vector<string> valVec_fixed(1+n+n*n);
				// Copy over number of markers and marker names
				for(int i = 0; i < 1+n; i++)
					valVec_fixed[i] = valVec[i];
				// For the actual values, collapse them down by adjacent values, assuming prior improper split
				for(int i = 1+n, j = 1+n; i < valVec_fixed.size(); i++, j+=2)
					valVec_fixed[i] = valVec[j] + "." + valVec[j+1];

				valVec = std::move(valVec_fixed);
				PRINT("$SPILLOVER size discrepancy was resolved by assuming non-standard decimal separator.\n");
			}else
				throw(domain_error("$SPILLOVER size discrepancy could not be resolved!\n"));

		}
		unordered_map<string, queue<int>> chnls;
		if(n > 0)
		{
			spillOver.resize(n*n);
			marker.resize(n);
			bool isDuplicate = false;
			for(int i = 0; i < n; i++)//param name
			{
				marker[i] = valVec[i+1];

				// Keep track of where this marker has appeared
				auto found = chnls.find(marker[i]);
				if( found == chnls.end()){
					chnls[marker[i]] = queue<int>();
					chnls[marker[i]].push(i);
				}else{
					isDuplicate = true;
					found->second.push(i);
				}
			}

			// Disambiguate duplicates by appending -<N>
			if(isDuplicate){
				PRINT("channel_alias: Duplicate channel names in spillover matrix!\n"
				"Integer suffixes added to disambiguate channels.\n"
				"It is recommended to verify correct mapping of spillover matrix columns.\n");
				for ( auto chnl : chnls ){
					if( chnl.second.size() > 1 ){
						int dup_idx;
						int suffix = 1;
						while( !chnl.second.empty()){
							dup_idx = chnl.second.front();
							chnl.second.pop();
							marker[dup_idx] = marker[dup_idx] + "-" + std::to_string(suffix++);
						}
					}
				}
			}
			// This is not safe without the size pre-checks. Additionally, vector[] will not even
			// yield an out_of_range exception so it just results in an ugly segfault if size is off
			for(unsigned i = n + 1; i < valVec.size(); i++)//param name
			{
				//Removing any formatting padding from spillover matrix entries
				boost::algorithm::trim(valVec[i]);
				spillOver[i-n-1] = boost::lexical_cast<double>(valVec[i]);
			}

		}

	}
	/**
	 * convert spillover matrix to FCS keyword string
	 * @return
	 */
	string to_string() const
	{
		auto n = marker.size();
		string res = std::to_string(n);
		for(auto c : marker)
		{
			res += "," + c;
		}
		for(auto i : spillOver)
			res += "," + std::to_string(i);
		return res;
	}
	/**
	 * convert spillover matrix from row-majored std::vector
	 * to col-majored arma::ma
	 * @return
	 */
	mat get_spillover_mat ()const;
	void update_channels(const CHANNEL_MAP & chnl_map);
	void convertToPb(pb::COMP & comp_pb);
	compensation(const pb::COMP & comp_pb);
	bool empty() const{return marker.size() == 0;}

};

};

#endif /* INCLUDE_COMPENSATION_HPP_ */
