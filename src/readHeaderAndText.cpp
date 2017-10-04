/*
 * readHeaderAndText.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#include "cytolib/readFCSHeader.hpp"
void readFCSHeader(ifstream &in, FCS_Header & header, int nOffset = 0){
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
	    header.textstart = stoi(tmp1) + nOffset;
	    in.get(tmp1, 9);
		header.textend = stof(tmp1) + nOffset;
		in.get(tmp1, 9);
		header.datastart = stof(tmp1) + nOffset;
		in.get(tmp1, 9);
		header.dataend = stof(tmp1) + nOffset;
		in.get(tmp1, 9);
		header.anastart = stof(tmp1) + nOffset;
		in.get(tmp1, 9);
		header.anaend = stof(tmp1) + nOffset;

		header.additional = nOffset;

}




void fcsTextParse(string txt, KEY_WORDS & pairs, bool emptyValue){
	/*
	 * get the first character as delimiter
	 */
	char delimiter = txt[0];

	/*
	 * check if string ends with delimiter
	 */
	bool isDelimiterEnd = txt[txt.size()-1] == delimiter;



	std::string doubleDelimiter,magicString;
	doubleDelimiter.push_back(delimiter);
	doubleDelimiter.push_back(delimiter);
	//search for the first unused odd char as replacememnt for double delimiter
	//FCS 3.1 states only 0-126 ASCII are legal delimiter, but we can't assume the file always follows the standard
	//also the TEXT main contain some special characters , thus we want to make sure the replacement char is not used anywhere in FCS TEXT
	char oddChar = 127;
	for(; oddChar < 256; oddChar++)
	{

		if(oddChar==delimiter||txt.find(oddChar)!=std::string::npos)
			continue;
		else
			break;
	}
	if(oddChar==256)
		throw(domain_error("Can't find the unused odd character from ASCII(127-255) in FSC TEXT section!"));

	std::string soddChar;
	soddChar.push_back(oddChar);
	/*
	 *	when empty value is allowed, we have to take the assumption that there is no double delimiters in any keys or values,
	 */
	if(!emptyValue)//replace the double delimiter with the odd char
		boost::replace_all(txt, doubleDelimiter, soddChar);
	std::vector<std::string> tokens;
	boost::split(tokens, txt, [delimiter](char c){return c == delimiter;});

	unsigned j = isDelimiterEnd?tokens.size()-2:tokens.size()-1;//last token, skip the last empty one when end with delimiter
	string key;
	for(unsigned i = 1; i <= j; i++){//counter, start with 1 to skip the first empty tokens
		std::string token = tokens[i];
//			std::PRINT( token + " ");
		if(!emptyValue){
			/*
			 * restore double delimiter when needed
			 * (this slows down things quite a bit, but still a lot faster than R version,
			 *  and this double delimiter logic is not normally invoked anyway)
			 */
			boost::replace_all(token, soddChar, doubleDelimiter);
//				std::PRINT( token;
		}
//			std::PRINT( std::endl;

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
	 if(j%2 == 1){
		 std::string serror = "uneven number of tokens: \n";
	     serror.append(boost::lexical_cast<std::string>(j));
	     PRINT(serror);
	     PRINT("The last keyword is dropped.\n");
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


void readHeaderAndText(ifstream &in,FCS_Header & header, KEY_WORDS & keys, vector<cytoParam> & params, const FCS_READ_HEADER_PARAM & config){


	 //search the stream for the header and txt of the nth DataSet
	int nOffset = 0, nNextdata = 0;
	int n = config.nDataset <=0?1:config.nDataset;
	 //non C-style index: starting from 1
	for(int i = 1; i <= n; i++)
	{
		nOffset += nNextdata;
		readFCSHeader(in,header, nOffset);//read the header
		readFCStext(in, header, keys, config.isEmptyKeyValue);//read the txt section

		if(keys.find("$NEXTDATA")!=keys.end())
			nNextdata = boost::lexical_cast<int>(keys["$NEXTDATA"]);
		else
		{
			if(i<n)
				throw(domain_error("Can't find " + boost::lexical_cast<string>(n) + "th dataset in FCS!"));
			break;
		}


	}


	  if(config.nDataset <=0 && nNextdata >0)
	  {
		  PRINT("The file contains additional data segment%s.\n The default is to read the first segment only.\nPlease consider setting the 'dataset' argument.");

	  }

	/*
	 * checkOffset:Fix the offset when its values recorded in header and TEXT don't agree
	 */
	//##for DATA segment exceeding 99,999,999 byte.
	 if(header.FCSversion >= 3)
	 {
	   int datastart_h = header.datastart - header.additional;
	   int dataend_h = header.dataend - header.additional;

	//
	//   # Let's not be too strick here as unfortunatelly, some files exported from FlowJo
	//   # are missing the $BEGINDATA and $ENDDATA keywords and we still need to read those
	   int datastart, dataend;
	   if(keys.find("$BEGINDATA")==keys.end())
	   {
	     if (datastart_h != 0)
	     {
	       datastart = datastart_h;
	       PRINT("warning:Missing the required $BEGINDATA keyword! Reading data based on information in the FCS HEADER only.\n");
	     } else
	       throw(domain_error("Don't know where the data segment begins, there was no $BEGINDATA keyword and the FCS HEADER does not say it either."));
	   }
	   else
		   datastart = boost::lexical_cast<int>(keys["$BEGINDATA"]);

	   if(keys.find("$ENDDATA")==keys.end())
	   {
		 if (dataend_h != 0) {
			 dataend = dataend_h;
			 PRINT("warning:Missing the required $ENDDATA keyword! Reading data based on information in the FCS HEADER only.\n");
		 } else
		   throw(domain_error("Don't know where the data segment ends, there was no $ENDDATA keyword and the FCS HEADER does not say it either."));
	   }
	   else
		   dataend = boost::lexical_cast<int>(keys["$ENDDATA"]);

	//   # when both are present and they don't agree with each other
	   if(datastart_h != datastart)
	   {
		  if(datastart_h== 0) //#use the TEXT when header is 0
			  header.datastart =  datastart + header.additional;
		  else
		  {//#trust the header when it is non-zero
			  string msg = "The HEADER and the TEXT segment define different starting point (";
			  msg.append(boost::lexical_cast<string>(header.datastart) + ":" + boost::lexical_cast<string>(datastart) + ") to read the data.");
			 if(config.ignoreTextOffset)
			 {
				 msg.append(" The values in TEXT are ignored!\n");
				 PRINT(msg);
			 }
			 else
			   throw(domain_error(msg));
		  }
	   }
	//   #both are present and they don't agree
	   if(dataend_h != dataend)
	   {
		   if(dataend_h== 0 || dataend_h== 99999999)//#use TEXT when either header is 0 or TEXT is 99999999
			header.dataend = dataend + header.additional;
		   else
			{//#otherwise trust the header
				string msg = "The HEADER and the TEXT segment define different ending point (";
				msg.append(boost::lexical_cast<string>(header.dataend) + ":" + boost::lexical_cast<string>(dataend) + ") to read the data.");
				if(config.ignoreTextOffset)
				 {
					 msg.append(" The values in TEXT are ignored!\n");
					 PRINT(msg);
				 }
				else
				   throw(domain_error(msg));
			}
	   }

	}

	 //parse important params from keys
	 int nrpar = boost::lexical_cast<int>(keys["$PAR"]);
	params.resize(nrpar);
	KEY_WORDS::iterator it;
	for(int i = 1; i <= nrpar; i++)
	{
		string pid = to_string(i);
		string range_str;
		if( keys.find("transformation")!=keys.end() &&  keys["transformation"] == "custom")
				range_str = keys["flowCore_$P" + pid + "Rmax"];
		else
				range_str = keys["$P" + pid + "R"];

		params[i-1].max = boost::lexical_cast<EVENT_DATA_TYPE>(range_str);


		params[i-1].PnB = boost::lexical_cast<int>(keys["$P" + pid + "B"]);

		it = keys.find("$P" + pid + "E");
		if(it==keys.end())
			params[i-1].PnE = make_pair<int,int>(0,0);
		else
		{
			vector<string> tokens;
			boost::split(tokens, it->second, boost::is_any_of(","));
			params[i-1].PnE = make_pair<int,int>(boost::lexical_cast<int>(tokens[0]),boost::lexical_cast<int>(tokens[1]));
		}

		it = keys.find("$P" + pid + "G");
		if(it==keys.end())
			params[i-1].PnG = 1;
		else
		{

			params[i-1].PnG = boost::lexical_cast<EVENT_DATA_TYPE>(it->second);
		}

		params[i-1].channel = keys["$P" + pid + "N"];
		params[i-1].marker = keys["$P" + pid + "S"];
		boost::trim(params[i-1].marker);
	}




//	    origRange <- range
//	    range <- rbind(realMin,range-1)


}
