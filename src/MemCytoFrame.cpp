/*
 * MemCytoFrame.cpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#include "cytolib/MemCytoFrame.hpp"

MemCytoFrame::~MemCytoFrame(){
	if(data)
		delete [] data;
}
void readFCSHeader(ifstream &in, FCS_Header & header, int start = 0){
	/*
		 * parse the header
		 */

		//parse version
		char version[7];
		in.get(version, 7);

	    if(strcmp(version, "FCS2.0")!=0&&strcmp(version, "FCS3.0")!=0&&strcmp(version, "FCS3.1")!=0)
		     throw(domain_error("This does not seem to be a valid FCS2.0, FCS3.0 or FCS3.1 file"));

	    header.FCSversion = boost::lexical_cast<float>(version+3);

	    char tmp[5];
	    in.get(tmp, 5);
		if(strcmp(tmp, "    "))
			 throw(domain_error("This does not seem to be a valid FCS header"));

		//parse offset
		char tmp1[9];
		in.get(tmp1, 9);
	    header.textstart = stoi(tmp1) + start;
	    in.get(tmp1, 9);
		header.textend = stof(tmp1) + start;
		in.get(tmp1, 9);
		header.datastart = stof(tmp1) + start;
		in.get(tmp1, 9);
		header.dataend = stof(tmp1) + start;
		in.get(tmp1, 9);
		header.anastart = stof(tmp1) + start;
		in.get(tmp1, 9);
		header.anaend = stof(tmp1) + start;

		header.additional = start;
}

void fcsTextParse(string txt, KEY_WORDS & pairs, bool emptyValue){

		/*
		 * get the first character as delimiter
		 */
		std::string delimiter = txt.substr(0,1);

		/*
		 * check if string ends with delimiter
		 */
		bool isDelimiterEnd = txt.substr(txt.size()-1, 1) == delimiter;

//		regexes require double-escaping (*sigh*)
//		if(delimiter == "\\" || delimiter == "|")
			delimiter = "\\" + delimiter;



		std::string doubleDelimiter,magicString;
		doubleDelimiter = delimiter + delimiter;
		magicString = "\\0QuickAndDirty\\0";
//		std::cout << doubleDelimiter << ":" << magicString <<std::endl;
		unsigned i = 0; //counter
		string key;
		/*
		 *	when empty value is allowed, we have to take the assumption that there is no double delimiters in any keys or values,
		 */
		if(!emptyValue)//replace the double delimiter with a magic strings
			txt = boost::regex_replace(txt, boost::regex(doubleDelimiter), magicString);//somehow boost::replace_all won't do the job for \\\\
		std::cout << txt << std::endl;

		/*
		 * then split by single delimiter
		 */
		boost::sregex_token_iterator token_begin(txt.begin() + 1, txt.end(), boost::regex(delimiter), -1), token_end;
		while(token_begin != token_end){
			i++;
			std::string token = *token_begin++;
//			std::cout << token << " ";
			if(!emptyValue){
				/*
				 * restore double delimiter when needed
				 * (this slows down things quite a bit, but still a lot faster than R version,
				 *  and this double delimiter logic is not normally invoked anyway)
				 */
				token = boost::regex_replace(token, boost::regex(magicString), doubleDelimiter);
//				std::cout << token;
			}
//			std::cout << std::endl;

			if((i)%2 == 1)
			{
				if(token.empty())
					// Rcpp::stop (temporarily switch from stop to range_error due to a bug in Rcpp 0.12.8)
					throw std::range_error("Empty keyword name detected!If it is due to the double delimiters in keyword value, please set emptyValue to FALSE and try again!");
				boost::trim(token);
				key = token;//set key
			}
			else{
				pairs[key] = token;//set value
			}


		}

		/*
		 * check if kw and value are paired
		 */
		 if(i%2 == 1){
			 if(isDelimiterEnd){
			   // Rcpp::stop
			   std::string serror = "uneven number of tokens: ";
			   serror.append(boost::lexical_cast<std::string>(i-1));
			   throw std::range_error(serror.c_str());
			 }
			 else
				 cout << "the text section does not end with delimiter: " << delimiter << ". The last keyword is dropped." << std::endl;;
		 }
}

void readFCStext(ifstream &in, const FCS_Header & header, KEY_WORDS & pairs, bool emptyValue){
	 in.seekg(header.textstart);
	    /**
	     *  Certain software (e.g. FlowJo 8 on OS X) likes to put characters into
	    files that readChar can't read, yet readBin, rawToChar and iconv can
	     handle just fine.
	     */
//	    txt <- readBin(con,"raw", offsets["textend"]-offsets["textstart"]+1)
//	    txt <- iconv(rawToChar(txt), "", "latin1", sub="byte")
	 int nTxt = header.textend - header.textstart + 1;
	 char * tmp = new char[nTxt + 1];
	 in.get(tmp, nTxt + 1);
	 string txt(tmp);
	 delete [] tmp;
     fcsTextParse(txt, pairs, emptyValue);

	if(pairs.find("FCSversion")==pairs.end())
	  pairs["FCSversion"] = boost::lexical_cast<string>(header.FCSversion);

}
MemCytoFrame::MemCytoFrame(const string &filename, bool emptyValue = false, int nDataset = 1, bool scale =false, double decades=0, double min_limit=-111, bool truncate_max_range = true, bool onlyTxt = false){
	ifstream in(filename, ios::in|ios::binary);

	FCS_Header header;
	readFCSHeader(in, header);
	readFCStext(in, header, CytoFrame::keys, emptyValue);
	if(onlyTxt)
		data = NULL;
	else
	{
		//parse the data section
	}
}


EVENT_DATA_TYPE * MemCytoFrame::getData(){
	return data;
}
EVENT_DATA_TYPE * MemCytoFrame::getData(const string & colname, ColType type){
	int idx = getColId(colname, type);
	return data + idx * nRow();
}



void MemCytoFrame:: compensate(const compensation &){

}

void MemCytoFrame:: save(const string & filename, FrameType type){

}
