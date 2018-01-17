/*
 * MemCytoFrame.hpp
 *
 *  Created on: Sep 18, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#define INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_
#include "CytoFrame.hpp"
#include "readFCSdata.hpp"

/**
 * the container that stores the different FCS parse arguments
 */
struct FCS_READ_PARAM{
	FCS_READ_HEADER_PARAM header;
	FCS_READ_DATA_PARAM data;
};

/**
 * The class represents the in-memory version of CytoFrame, which stores and owns the events data
 */
class MemCytoFrame: public CytoFrame{
	EVENT_DATA_VEC data;

public:
	MemCytoFrame(){};
	/**
	 * Constructor from a generic CytoFrame object
	 * @param frm a reference to CytoFrame
	 */
	MemCytoFrame(const CytoFrame & frm)
	{

		keys.setPairs(frm.getKeywords());
		params = frm.getParams();
		buildHash();
		data = frm.getData();
	}
	/**
	 * Constructor from the FCS file
	 *
	 * @param filename FCS file path
	 * @param config the parse arguments.
	 * @param onlyTxt flag indicates whether to only parse text segment (which contains the keywords)
	 */
	MemCytoFrame(const string &filename, FCS_READ_PARAM & config,  bool onlyTxt = false){
		ifstream in(filename, ios::in|ios::binary);

		if(!in.is_open())
			throw(domain_error("can't open the file: " + filename + "\nPlease check if the path is normalized to be recognized by c++!"));
		FCS_Header header;
		readHeaderAndText(in, header, keys, params, config.header);

		buildHash();

		if(!onlyTxt)
		{
	//		double start = clock();
			//parse the data section
			readFCSdata(in, data,header, keys, params, config.data);

	//		cout << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << endl;
			//update min and max
			/*
			 * ## set transformed flag and fix the PnE and the Datatype keywords
			 */
			keys["FILENAME"] = filename;
			if(config.data.isTransformed)
			{
				keys["transformation"] ="applied";
				keys["$DATATYPE"] = "F";
			}
			for(int i = 0; i < params.size(); i++)
			{

				string pid = to_string(i+1);

				//insert our own PnR fields
				if(config.data.isTransformed)
				{
					keys["$P" + pid + "E"] = "0,0";
					params[i].PnE[0] = 0;
					params[i].PnE[1] = 0;
					keys["flowCore_$P" + pid + "Rmax"] = to_string(static_cast<int>(params[i].max + 1));
					keys["flowCore_$P" + pid + "Rmin"] = to_string(static_cast<int>(params[i].min));
				}
				else
					params[i].max--;
			}



			//GUID
			string oldguid;
			if(keys.find("GUID")!=keys.end()){
				oldguid = keys["GUID"];
				keys["ORIGINALGUID"] = oldguid;
			}
			else
				oldguid = filename;
			//strip dir
			vector<string> paths;
			boost::split(paths, oldguid, boost::is_any_of("/\\"));
			keys["GUID"] = paths.back();

		}



	}
	/**
	 * parse the data segment of FCS
	 *
	 * @param in (input) file stream object opened from FCS file
	 * @param output (output) the container to store the parsed data
	 * @param header (input) the header parsed by readHeaderAndText
	 * @param keys (input&output) the keywords parsed by readHeaderAndText
	 * @param params (input&output) the params parsed by readHeaderAndText
	 * @param config (input) the parsing arguments for data
	 */
	void readFCSdata(ifstream &in, EVENT_DATA_VEC & output, const FCS_Header & header,KEY_WORDS & keys, vector<cytoParam> & params,  FCS_READ_DATA_PARAM & config)
	{
		//## transform or scale data?
		  bool fcsPnGtransform = false, isTransformation, scale;
		  if(config.transform == TransformType::linearize)
		  {
			isTransformation =  true;
			scale = false;
		  } else if (config.transform == TransformType::scale) {
			  isTransformation = true;
			scale = true;
		  } else if (config.transform == TransformType::linearize_with_PnG_scaling) {
			  isTransformation =  true;
				scale = false;
			fcsPnGtransform = true;
		  } else if (config.transform == TransformType::none) {
			  isTransformation =  false;
			scale = false;
		  }

		  if (fcsPnGtransform)
			  keys["flowCore_fcsPnGtransform"] = "linearize-with-PnG-scaling";
		  bool transDefinedinKeys = keys.find("transformation")!=keys.end();
		 if(transDefinedinKeys)
			 if(keys["transformation"] == "applied"||keys["transformation"] ==  "custom")
				 isTransformation =  false;


		string byte_order = keys["$BYTEORD"];
		endianType endian;
		if(byte_order == "4,3,2,1" || byte_order == "2,1")
			endian = endianType::big;
		else if(byte_order == "1,2,3,4" || byte_order == "1,2")
			endian = endianType::small;
		else
			endian = endianType::mixed;



		 string dattype = keys["$DATATYPE"];
		 if(dattype!="I"&&dattype!="F"&&dattype!="D")
			 throw(domain_error("Don't know how to deal with $DATATYPE"));


	    if (keys["$MODE"] != "L")
	    	throw(domain_error("Don't know how to deal with $MODE " + keys["$MODE"]));


		fcsPnGtransform = keys.find("flowCore_fcsPnGtransform")!= keys.end() && keys["flowCore_fcsPnGtransform"] == "linearize-with-PnG-scaling";

	//	int nrowTotal= boost::lexical_cast<int>(keys["$TOT"]);
		int nCol = params.size();

		bool multiSize = false;
		for(int i = 1; i < nCol; i++)
		{
			if(params[i].PnB%8)
				throw(domain_error("Sorry, C parser doesn't support odd bitwidth!"));
			if(params[i].PnB!=params[0].PnB)
			{
				multiSize = true;
				break;
			}
		}
		if(dattype!="I"&&multiSize)
			throw(domain_error("Sorry, Numeric data type expects the same bitwidth for all parameters!"));

		bool splitInt;
		if(dattype=="I"){

	//	  vector<unsigned> range(nCol);
	//
	//	  for(int i = 1; i <= nCol; i++)
	//		  range[i-1] = boost::lexical_cast<unsigned>(range_str[i-1]);
		  if(multiSize)
			splitInt = false;
		  else
			  splitInt = params[0].PnB == 32;


		}
		else
		{
		  splitInt = false;

		}



		if(!multiSize){
		  if(params[0].PnB ==10){
			  string sys = keys["$SYS"];
			  transform(sys.begin(), sys.end(), sys.begin(),::tolower);
			if(sys !=  "cxp")
			  PRINT("Invalid bitwidth specification.\nThis is a known bug in Beckman Coulter's CPX software.\nThe data might be corrupted if produced by another software.\n");
			else
				PRINT("Beckma Coulter CPX data.\nCorrected for invalid bitwidth 10.\n");

			for(auto &p : params)
				p.PnB = 16;
		  }
		}

		bool isSigned;
		if(multiSize){

		  isSigned = false; // #dummy. not used in mutliSize logic.
		}else{


		  //# since signed = FALSE is not supported by readBin when size > 2
		  //# we set it to TRUE automatically then to avoid warning flooded by readBin
		  //# It shouldn't cause data clipping since we haven't found any use case where datatype is unsigned integer with size > 16bits
		 isSigned = !(params[0].PnB == 8 ||params[0].PnB == 16);
		}



	  in.seekg(header.datastart);

	  int nBytes = header.dataend - header.datastart + 1;


	//	vector<BYTE>bytes(nBytes);
	//	in.read((char *)&bytes[0], nBytes);
	//
	//	if(dattype != "i")
	//	  throw(domain_error("we don't support different bitwdiths for numeric data type!"));
		//total bits for each row
	  	size_t nRowSize = accumulate(params.begin(), params.end(), 0, [](size_t i, cytoParam p){return i + p.PnB;});

	  	unsigned nrow = nBytes * 8/nRowSize;

	  	vector<int>which_lines = config.which_lines;
	  	int nSelected = which_lines.size();
	  	if(nSelected>0){
	  		if(nSelected >= nrow)
	  			throw(domain_error("total number of which.lines exceeds the total number of events: " + to_string(nrow)));

	  		sort(which_lines.begin(), which_lines.end());
	  		nrow = nSelected;
	  		nBytes = nrow * nRowSize/8;
	  	}
	  	unique_ptr<char []> buf(new char[nBytes]);//we need to rearrange dat from row-major to col-major thus need a separate buf anyway (even for float)
	  	char * bufPtr = buf.get();
	  	if(nSelected>0)
	  	{
	  		char * thisBufPtr = bufPtr;
	  		auto nRowSizeBytes = nRowSize/8;
	  		for(auto i : which_lines)
	  		{
	  			auto pos =  header.datastart + i * nRowSizeBytes;
	  			if(pos > header.dataend || pos < header.datastart)
	  				throw(domain_error("the index of which.lines exceeds the data boundary: " + to_string(i)));
	  			in.seekg(pos);
	  			in.read(thisBufPtr, nRowSizeBytes);
	  			thisBufPtr += nRowSizeBytes;
	  		}
	  	}
	  	else
	  	{
	  		//load entire data section with one disk IO

			in.read(bufPtr, nBytes); //load the bytes from file
	  	}
	//	nEvents = nrow;
		//how many element to return
		auto nElement = nrow * nCol;
	//	EVENT_DATA_PTR output(new EVENT_DATA_TYPE[nElement]);
		output.resize(nElement);

	//	char *p = buf.get();//pointer to the current beginning byte location of the processing data element in the byte stream
		float decade = pow(10, config.decades);
		/*
		 * mixed endian parsing could be more efficiently
		 * integrated into the subsequent main loop of data parsing
		 * but since it is a rare case (legacy data), we do the separate simple
		 * preprocessing here to avoid adding extra overhead into the main loop
		 */
		if(endian == endianType::mixed)
		{
			  if(multiSize)
				throw(domain_error("Cant't handle diverse bitwidths while endian is mixed: " + byte_order));

			  vector<string> byteOrd;
			  boost::split(byteOrd, byte_order, boost::is_any_of(","));
			  int elementSize = byteOrd.size();

			  if(params[0].PnB/8 != elementSize)
				throw(domain_error("Byte order is not consistent with bidwidths!"));

			  vector<int> iByteOrd(elementSize);
			  for(auto i = 0; i < elementSize; i++)
			  {
				  iByteOrd[i] = boost::lexical_cast<int>(byteOrd[i])-1;
			  }
			  char * tmp = new char[elementSize];
			  for(auto ind = 0; ind < nElement; ind++){

				  memcpy(tmp, bufPtr + ind * elementSize, elementSize);

			     for(auto i = 0; i < elementSize; i++){
			       auto j = iByteOrd.at(i);

			       auto pos_old = ind * elementSize + i;
			       auto pos_new = ind * elementSize + j;
			 //       if(ind<=10)
			 //         Rcpp::Rcout << pos_old <<":" << pos_new << std::endl;


			       bufPtr[pos_new] = tmp[i];

			     }

			   }
			  	delete [] tmp;

			    endian = endianType::small;
		}

		bool isbyteswap = false;

		if((is_host_big_endian()&&endian==endianType::small)||(!is_host_big_endian()&&endian==endianType::big))
			isbyteswap = true;

		/**
		 * cp raw bytes(row-major) to a 2d mat (col-major) represented as 1d array(with different byte width for each elements)
		 */
	//	double start = omp_get_wtime();//clock();
	#ifdef _OPENMP
		omp_set_num_threads(config.num_threads);
	#endif

	    #pragma omp parallel for
		for(auto c = 0; c < nCol; c++)
		 {
			string pid = to_string(c+1);
			EVENT_DATA_TYPE realMin = numeric_limits<EVENT_DATA_TYPE>::max();
			size_t element_offset = nrow * c;
			size_t bits_offset = accumulate(params.begin(), params.begin() + c, 0, [](size_t i, cytoParam p){return i + p.PnB;});
			cytoParam & param = params[c];
			int usedBits = ceil(log2(param.max));
		    uint64_t base = static_cast<uint64_t>(1)<<usedBits;
		    auto thisSize = params[c-1].PnB;
			for(auto r = 0; r < nrow; r++)
		    {
		      //convert each element
				  auto thisSize = param.PnB;
				  size_t idx = element_offset + r;
				  EVENT_DATA_TYPE & outElement = output[idx];
				  size_t idx_bits = r * nRowSize + bits_offset;
				  char *p = bufPtr + idx_bits/8;
				  thisSize/=8;
				  if(isbyteswap)
					  std::reverse(p, p + thisSize);

				  if(dattype == "I")
				  {
					  switch(thisSize)
					  {
					  case sizeof(BYTE)://1 byte
						{
						  outElement = static_cast<EVENT_DATA_TYPE>(*p);
						}

						  break;
					  case sizeof(unsigned short): //2 bytes
						{

						  outElement = static_cast<EVENT_DATA_TYPE>(*reinterpret_cast<unsigned short *>(p));
						}

						break;
					  case sizeof(unsigned)://4 bytes
						{
						  outElement = static_cast<EVENT_DATA_TYPE>(*reinterpret_cast<unsigned *>(p));
						}

						break;
					  case sizeof(uint64_t)://8 bytes
						{
						  outElement = static_cast<EVENT_DATA_TYPE>(*reinterpret_cast<uint64_t *>(p));
						}

					  break;
					  default:
						  {
							  std::string serror = "unsupported byte width :";
							  serror.append(std::to_string(thisSize));
							  throw std::range_error(serror.c_str());
						  }
					  }
					  // apply bitmask for integer data

					 if(param.max > 0)
					 {

					   if(usedBits < param.PnB)
						   outElement = static_cast<uint64_t>(outElement) % base;
					}


				}
				else
				{
				  switch(thisSize)
				  {
				  case sizeof(float):
					{
					  outElement = *reinterpret_cast<float *>(p);
					}

					break;
				  case sizeof(double):
					{
					  outElement = static_cast<EVENT_DATA_TYPE>(*reinterpret_cast<double *>(p));
					}

					break;
				  default:
					std::string serror ="Unsupported bitwidths for numerical data type:";
					serror.append(std::to_string(thisSize));
					throw std::range_error(serror.c_str());
				  }
				}

			  // truncate data at range
				if(!transDefinedinKeys)
				{
					if(config.truncate_max_range&&outElement > param.max)
						outElement = param.max;

					if(config.truncate_min_val&&outElement < config.min_limit)
						outElement = config.min_limit;
				}



	//				## Transform or scale if necessary
	//				# J.Spidlen, Nov 13, 2013: added the flowCore_fcsPnGtransform keyword, which is
	//				# set to "linearize-with-PnG-scaling" when transformation="linearize-with-PnG-scaling"
	//				# in read.FCS(). This does linearization for log-stored parameters and also division by
	//				# gain ($PnG value) for linearly stored parameters. This is how the channel-to-scale
	//				# transformation should be done according to the FCS specification (and according to
	//				# Gating-ML 2.0), but lots of software tools are ignoring the $PnG division. I added it
	//				# so that it is only done when specifically asked for so that read.FCS remains backwards
	//				# compatible with previous versions.


				if(isTransformation)
				{


				  if(param.PnE[0] > 0)
				  {
		//				 # J.Spidlen, Nov 5, 2013: This was a very minor bug. The linearization transformation
		//				 # for $PnE != "0,0" is defined as:
		//				 # For $PnR/r/, r>0, $PnE/f,0/, f>0: n is a logarithmic parameter with channel values
		//				 # from 0 to r-1. A channel value xc is converted to a scale value xs as xs=10^(f*xc/r).
		//				 # Note the "r" instead of the "r-1" in the formula (which would admitedly make more sense)
		//				 # However, this is the standard that apparently has been followed by BD and other companies
		//				 # "forever" and it is therefore addoped as such by the ISAC DSTF (see FCS 3.1 specification)
		//				 # To bring this to compliance, I am just changing
		//				 # dat[,i] <- 10^((dat[,i]/(range[i]-1))*ampli[i,1])
		//				 # to
		//				 # dat[,i] <- 10^((dat[,i]/range[i])*ampli[i,1])
					 outElement = pow(10,((outElement/param.max)*param.PnE[0]));
	//				 param.max = pow(10,param.PnE.first);
				  }
				  else if (fcsPnGtransform && param.PnG != 1) {
					  outElement = outElement / param.PnG;
	//				  param.max = (param.max-1) / param.PnG;
				  }
	//			  else
	//				  param.max--;
				}
				if(scale)
				{
					if(param.PnE[0] > 0)
					{
						outElement = decade*((outElement-1)/(param.max-1));
	//					param.max = decade*(param.max/param.max-1);
					}
					else
					{
						outElement = decade*((outElement)/(param.max));
	//					param.max = decade;
					}
				}
				realMin = realMin > outElement?outElement:realMin;

		    }

			if(keys.find("transformation")!=keys.end() &&  keys["transformation"] == "custom")
				param.min = boost::lexical_cast<EVENT_DATA_TYPE>(keys["flowCore_$P" + pid + "Rmin"]);
			else
			{

				auto zeroVals = param.PnE[1];
				param.min = min(zeroVals, max(config.min_limit, realMin));

			}

		 }
	//	 cout << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << endl;
	//	 cout << (omp_get_wtime() - start) / (double)(CLOCKS_PER_SEC / 1000) << endl;

		//update params
		for(auto &p : params)
		{

			if(isTransformation)
			{


			  if(p.PnE[0] > 0)
			  {
				 p.min = pow(10,p.min/(p.max-1) * p.PnE[0]);
				 p.max = pow(10,p.PnE[0]);
			  }
			  else if (fcsPnGtransform && p.PnG != 1)
				  p.max = (p.max-1) / p.PnG;
			  else
				  p.max--;
			}
			if(scale)
			{
				if(p.PnE[0] > 0)
					p.max = decade*(p.max/p.max-1);
				else
					p.max = decade;

			}


		}
		config.isTransformed = isTransformation;


	}
	void readFCSHeader(ifstream &in, FCS_Header & header, int nOffset = 0){
		/*
			 * parse the header
			 */

			in.seekg(nOffset);
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
			string tmp2(" ",8);
			in.get(tmp1, 9);
			//skip whitespaces
			copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
		    header.textstart = boost::lexical_cast<int>(tmp2) + nOffset;

		    in.get(tmp1, 9);
		    tmp2.resize(8);
		    copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
			header.textend = boost::lexical_cast<int>(tmp2) + nOffset;

			in.get(tmp1, 9);
			tmp2.resize(8);
			copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
			header.datastart = boost::lexical_cast<int>(tmp2) + nOffset;

			in.get(tmp1, 9);
			tmp2.resize(8);
			copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
			header.dataend = boost::lexical_cast<int>(tmp2) + nOffset;

			in.get(tmp1, 9);
			tmp2.resize(8);
			copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
			if(tmp2.size()>0)
				header.anastart = boost::lexical_cast<int>(tmp2) + nOffset;

			in.get(tmp1, 9);
			tmp2.resize(8);
			copy(tmp1, tmp1+8, tmp2.begin());
			boost::trim(tmp2);
			if(tmp2.size()>0)
				header.anaend = boost::lexical_cast<int>(tmp2) + nOffset;

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

	/**
	 * parse the FCS header and Text segment
	 *
	 * @param in (input) the file stream object opened from FCS file
	 * @param header (output) the FCS_Header object to store the result of FCS header
	 * @param keys (output) the KEY_WORDS object to store parsed keywords
	 * @param params (output) cytoParam container to store the parsed parameters
	 * @param config (input) FCS_READ_HEADER_PARAM object gives the parsing arguments for header
	 */
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

			if(keys.find("$NEXTDATA")!=keys.end()){
				string nd = keys["$NEXTDATA"];
				boost::trim(nd);
				if(nd.size()==0)
					throw(domain_error("empty value in $NEXTDATA"));
				else
				 nNextdata = boost::lexical_cast<int>(nd);
			}
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
		   unsigned long datastart, dataend;
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
			   datastart = stoi(keys["$BEGINDATA"]);

		   if(keys.find("$ENDDATA")==keys.end())
		   {
			 if (dataend_h != 0) {
				 dataend = dataend_h;
				 PRINT("warning:Missing the required $ENDDATA keyword! Reading data based on information in the FCS HEADER only.\n");
			 } else
			   throw(domain_error("Don't know where the data segment ends, there was no $ENDDATA keyword and the FCS HEADER does not say it either."));
		   }
		   else
			   dataend = stoul(keys["$ENDDATA"]);

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
		 int nrpar = stoi(keys["$PAR"]);
		params.resize(nrpar);
		KEY_WORDS::iterator it;
		for(int i = 1; i <= nrpar; i++)
		{
			string pid = to_string(i);
			string range_str;
			if( keys.find("transformation")!=keys.end() &&  keys["transformation"] == "custom")
				range_str = "flowCore_$P" + pid + "Rmax";
			else
				range_str = "$P" + pid + "R";
			it = keys.find(range_str);
			if(it==keys.end())
				throw(domain_error(range_str + " not contained in Text section!"));
			else
				params[i-1].max = boost::lexical_cast<EVENT_DATA_TYPE>(it->second);


			params[i-1].PnB = stoi(keys["$P" + pid + "B"]);

			it = keys.find("$P" + pid + "E");
			if(it==keys.end())
			{
				params[i-1].PnE[0] = 0;
				params[i-1].PnE[1] = 0;
			}
			else
			{
				vector<string> tokens;
				boost::split(tokens, it->second, boost::is_any_of(","));
				params[i-1].PnE[0] = stof(tokens[0]);
				params[i-1].PnE[1] = stof(tokens[1]);
			}

			it = keys.find("$P" + pid + "G");
			if(it==keys.end())
				params[i-1].PnG = 1;
			else
			{

				params[i-1].PnG = boost::lexical_cast<EVENT_DATA_TYPE>(it->second);
			}

			params[i-1].channel = keys["$P" + pid + "N"];
			it = keys.find("$P" + pid + "S");
			if(it!=keys.end())
				params[i-1].marker = keys["$P" + pid + "S"];
			boost::trim(params[i-1].marker);
		}




	//	    origRange <- range
	//	    range <- rbind(realMin,range-1)


	}


	int nRow() const{
		if(nCol()==0)
			return 0;
		else
			return data.size()/nCol();
	}

/**
 * Caller will receive a copy of data
 * @return
 */
	EVENT_DATA_VEC getData() const{
		return data;
	}
	EVENT_DATA_VEC getData(const string & colname, ColType type) const{
		int idx = getColId(colname, type);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		int nEvents = nRow();
		EVENT_DATA_VEC res(nEvents);
		memcpy(&res[0], &data[0] + idx * nEvents, nEvents*sizeof(EVENT_DATA_TYPE));
		return res ;
	}



	void  compensate(const compensation &){

	}



	/**
	 * return the pointer of a particular data column
	 *
	 * The reason to return the pointer is to keep it backward compatible
	 * TODO:modify gating apis to use iterator instead of raw pointer
	 * @param colname
	 * @param type
	 * @return
	 */
	EVENT_DATA_TYPE * subset(const string & colname, ColType type = ColType::unknown){
		int idx = getColId(colname, type);
		if(idx<0)
			throw(domain_error("colname not found: " + colname));
		return data.data() + idx * nRow();
	}


};

#ifdef _OPENMP
#define gettime() omp_get_wtime()
#else
#define gettime() clock()/(double)(CLOCKS_PER_SEC / 1000)
#endif




#endif /* INST_INCLUDE_CYTOLIB_MEMCYTOFRAME_HPP_ */
