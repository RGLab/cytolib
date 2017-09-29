/*
 * CytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/CytoFrame.hpp"



KEY_WORDS CytoFrame::getKeywords(){
	return keys;
}
string CytoFrame::getKeyword(const string & key){
	return keys[key];
}
void CytoFrame::setKeyword(const string & key, const string & value){
	keys[key] = value;
}
int CytoFrame::nCol(){
	return params.size();
}
int CytoFrame::nRow(){
	return nEvents;
}

bool CytoFrame::isHashed(){
	return channel_vs_idx.size()==nCol();
}
void CytoFrame::buildHash(){
	if(!isHashed())
	{
		for(auto i = 0; i < nCol(); i++)
		{
			channel_vs_idx[params[i].channel] = i;
			marker_vs_idx[params[i].marker] = i;
		}
	}
}
vector<string> CytoFrame::getChannels(){
	vector<string> res(nCol());
	for(auto i = 0; i < nCol(); i++)
		res[i] = params[i].channel;
	return res;
}
vector<string> CytoFrame::getMarkers(){
	vector<string> res(nCol());
		for(auto i = 0; i < nCol(); i++)
			res[i] = params[i].marker;
	return res;
}
void CytoFrame::setChannel(const string & oldname, const string & newname){
	int id = getColId(oldname, ColType::channel);
	params[id].channel=newname;
	channel_vs_idx.erase(oldname);
	channel_vs_idx[newname] = id;
}
void CytoFrame::setMarker(const string & oldname, const string & newname){
	int id = getColId(oldname, ColType::marker);
	params[id].marker=newname;
	marker_vs_idx.erase(oldname);
	marker_vs_idx[newname] = id;
}
int CytoFrame::getColId(const string & colname, ColType type = ColType::unknown){
	if(!isHashed())
		buildHash();

	switch(type)
	{
	case ColType::channel:
		{
			return channel_vs_idx[colname];
		}
	case ColType::marker:
		{
			return marker_vs_idx[colname];
		}
	case ColType::unknown:
	{
		unordered_map<string, int>::iterator it1 = channel_vs_idx.find(colname);
		unordered_map<string, int>::iterator it2 = marker_vs_idx.find(colname);
		if(it1==channel_vs_idx.end()&&it2==marker_vs_idx.end())
			throw(domain_error("colname not found: " + colname));
		else if(it1!=channel_vs_idx.end()&&it2!=marker_vs_idx.end())
			throw(domain_error("ambiguous colname without colType: " + colname ));
		else if(it1!=channel_vs_idx.end())
			return it1->second;
		else
			return it2->second;
	}
	default:
		throw(domain_error("invalid col type"));
	}

}

pair<EVENT_DATA_TYPE, EVENT_DATA_TYPE> CytoFrame::getRange(const string & colname, ColType ctype, RangeType type){

	switch(type)
	{
	case RangeType::data:
		{

			EVENT_DATA_TYPE * vec = getData(colname, ctype);
			auto res = minmax_element(vec, vec + nRow());
			return make_pair(*res.first, *res.second);
		}
	case RangeType::instrument:
	{
		int idx = getColId(colname, ctype);
		return make_pair(params[idx].min, params[idx].max);
	}
	default:
		throw(domain_error("invalid range type"));
	}
}
