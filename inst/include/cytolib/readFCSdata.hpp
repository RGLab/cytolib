/*
 * readFCSdata.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: wjiang2
 */

#ifndef INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#define INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_
#include "readFCSHeader.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif
typedef unsigned char BYTE;


/**
 * The struct stores all the parsing arguments for events data
 */
struct FCS_READ_DATA_PARAM{
	 bool scale, truncate_max_range, truncate_min_val;
	 EVENT_DATA_TYPE decades, min_limit;
	 TransformType transform;
	 int num_threads; //number of cores to be used for parallel-read of data (channel / core)
	 vector<int> which_lines; //select rows to be read in
	 bool isTransformed;//record the outcome after parsing
	 FCS_READ_DATA_PARAM(){
		 scale = false;
		 truncate_max_range = true;
		 truncate_min_val = false;
		 decades = 0;
		 min_limit=-111;
		 transform =  TransformType::linearize;
		 num_threads = 1;
		 isTransformed = false;
	 }


};
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





#endif /* INST_INCLUDE_CYTOLIB_READFCSDATA_HPP_ */
