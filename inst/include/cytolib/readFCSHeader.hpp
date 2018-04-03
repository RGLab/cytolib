/*
 * readFCSHeader.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_
#define INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_

#include <fstream>
#include <cstring>
#include "global.hpp"
#include <vector>
#include <numeric>
#include <unordered_map>
#include <iostream>
#include <algorithm>
using namespace std;


namespace cytolib
{
enum class TransformType {none, linearize, scale,linearize_with_PnG_scaling};
enum class endianType {big, small, mixed};
struct FCS_Header{
	float FCSversion;
	int textstart, textend, datastart, dataend, anastart, anaend, additional;
};

//typedef unordered_map <string, string> KEY_WORDS;
typedef vector<pair <string, string>> KW_PAIR;

/**
 * this class mimic the map behavior so that the same code
 * can be used for both map and vector based container
 */
class vec_kw_constainer{
 KW_PAIR kw;
public:
 typedef KW_PAIR::iterator iterator;
 typedef KW_PAIR::const_iterator const_iterator;
 void resize(size_t n){kw.resize(n);}
 size_t size(){return kw.size();}
 const KW_PAIR & getPairs() const{return kw;}
 void setPairs(const KW_PAIR & _kw){kw = _kw;}
 iterator end() {return kw.end();}
 const_iterator end() const{return kw.end();}
 iterator begin(){return kw.begin();}
 const_iterator begin() const{return kw.begin();}
 iterator find(const string &key){
         return std::find_if(kw.begin(), kw.end(), [key](const pair<string, string> & p){return p.first == key;});
 }
 const_iterator find(const string &key) const{
          return std::find_if(kw.begin(), kw.end(), [key](const pair<string, string> & p){return p.first == key;});
  }
 string & operator [](const string & key){
         iterator it = find(key);
         if(it==end())
         {
                 kw.push_back(pair<string, string>(key, ""));
                 return kw.back().second;
         }
         else
                 return it->second;
   }
 pair <string, string> & operator [](const int & n){
	 return kw[n];
 }
};


typedef vec_kw_constainer KEY_WORDS;


/**
 * the struct that stores the FCS  header parse arguments
 */
struct FCS_READ_HEADER_PARAM{
	bool is_fix_slash_in_channel_name;
	bool isEmptyKeyValue; //whether allow the keyword value to be empty. When true, then double delimiter will be considered as empty keyword value.
	bool ignoreTextOffset; //whether to ignore the offset recorded in the text segment of FCS
	int nDataset; // which data set to be parsed when multi-data segments are detected
	FCS_READ_HEADER_PARAM(){
		is_fix_slash_in_channel_name = false;
		isEmptyKeyValue = false;
		ignoreTextOffset = false;
		nDataset = 1;
	};

};
/**
 * parameters parsed from keywords, but may be accessed frequently later thus saved as the copy
 * in this struct for fast and easy query
 */
struct cytoParam{
	string channel, marker;
	EVENT_DATA_TYPE min, max, PnG;
//	pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> PnE;//replace pair with simple array since it is not clear how to create compound type for pair
	EVENT_DATA_TYPE PnE[2];
	int PnB;
};


};

#endif /* INST_INCLUDE_CYTOLIB_READFCSHEADER_HPP_ */
