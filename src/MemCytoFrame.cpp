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
EVENT_DATA_TYPE * readFCSdata(ifstream &in, const FCS_Header & header,const KEY_WORDS & pairs,  bool scale, double decades, double min_limit, bool truncate_max_range){
//	byte_order <- readFCSgetPar(x, "$BYTEORD")
//	    endian <- switch(byte_order,
//	                     "4,3,2,1" = "big",
//	                     "2,1" = "big",
//	                     "1,2" = "little",
//	                     "1,2,3,4" = "little",
//	                     "mixed")
//
//	    dattype <- switch(readFCSgetPar(x, "$DATATYPE"),
//	                      "I" = "integer",
//	                      "F" = "numeric",
//	                      "D" = "numeric",
//	                      stop(paste("Don't know how to deal with $DATATYPE",
//	                                 readFCSgetPar(x, "$DATATYPE"))))
//
//	    if (readFCSgetPar(x, "$MODE") != "L")
//	        stop(paste("Don't know how to deal with $MODE",
//	                   readFCSgetPar(x, "$MODE")))

//	    nrpar    <- as.integer(readFCSgetPar(x, "$PAR"))
//	    nrowTotal <- as.integer(readFCSgetPar(x, "$TOT"))
//
//	    if( "transformation" %in% names(x) &&  x[["transformation"]] == "custom"){
//	      range_str <- sapply(seq_len(nrpar),function(k){
//	                x[[sprintf("flowCore_$P%sRmax", k)]]
//	             })
//	    } else {
//	        range_str <- readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep=""))
//	    }
//
//	    bitwidth_vec <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
//	    bitwidth <- unique(bitwidth_vec)
//	    multiSize <- length(bitwidth) > 1
//
//	    if(dattype=="numeric"&&multiSize)
//	      stop("Sorry, Numeric data type expects the same bitwidth for all parameters!")
//
//
//	    if(dattype=="integer"){
//	      suppressWarnings(range <- as.integer(range_str))
//
//	      #when any channel has range > 2147483647 (i.e. 2^31-1)
//	      if(any(is.na(range)))
//	        range <- as.numeric(range_str)
//
//	      if(any(is.na(range)))#throws if still fails
//	        stop("$PnR", range_str[is.na(range)][1], "is larger than R's integer limit:", .Machine$integer.max)
//	      else if(any(range>2^32)){
//	        #check if larger than C's uint32 limit ,which should be 2^32-1
//	        #but strangely(and inaccurately) these flow data uses 2^32 to specifiy the upper bound of 32 uint
//	        #we try to tolerate this and hopefully there is no such extreme value exsiting in the actual data section
//	        stop("$PnR", range_str[range>2^32][1], "is larger than C's uint32 limit:", 2^32-1)
//	      }
//
//	      if(multiSize){
//	        splitInt <- FALSE
//	      }else
//	      {
//	        splitInt <- bitwidth == 32
//	      }
//
//	    }
//	    else{
//	      splitInt <- FALSE
//	      range <- as.numeric(range_str)
//	      if(any(is.na(range)))
//	        stop("$PnR", range_str[is.na(range)][1], "is larger than R's numeric limit:", .Machine$double.xmax)
//	    }
//
//
//
//	    if(!multiSize){
//	      if(bitwidth==10){
//	        if(!gsub(" " ,"", tolower(readFCSgetPar(x, "$SYS"))) ==  "cxp")
//	          warning("Invalid bitwidth specification.\nThis is a known bug in Beckman ",
//	                  "Coulter's CPX software.\nThe data might be corrupted if produced ",
//	                  "by another software.", call.=FALSE)
//	        else
//	          warning("Beckma Coulter CPX data.\nCorrected for invalid bitwidth 10.",
//	                  call.=FALSE)
//	        bitwidth <- 16
//	      }
//	    }
//	    # multiSize <- T
//
//	    if(multiSize){
//	      size <- bitwidth_vec/8
//	      signed <- FALSE #dummy. not used in mutliSize logic.
//	    }else{
//	      size <- bitwidth/8
//
//	      # since signed = FALSE is not supported by readBin when size > 2
//	      # we set it to TRUE automatically then to avoid warning flooded by readBin
//	      # It shouldn't cause data clipping since we haven't found any use case where datatype is unsigned integer with size > 16bits
//	      signed <- !(size%in%c(1,2))
//	    }
//
//
//	    nwhichLines <- length(which.lines)
//	    ##Read all reports
//	    if(is.null(which.lines) || (nwhichLines >  nrowTotal)){
//	        if (nwhichLines >  nrowTotal){
//	            cat("Warning: the number of lines specified (",nwhichLines,
//	                ") is greater
//	                 than the number of collected events (",nrowTotal,
//	                "). All the events have been read. \n")
//	        }
//	      seek(con, offsets["datastart"])
//
//	      nBytes <- offsets["dataend"]-offsets["datastart"]+1
//	      if(nBytes > .Machine$integer.max)
//	        stop("Total number of bytes (", nBytes, ") in data segment exceeds the R integer limits!Please read the subset of FCS by specifying 'which.lines'")
//	      nBytes <- as.integer(nBytes)
//		    if(multiSize){
//	#	      if(splitInt&&dattype=="integer")
//	#	        stop("Mutliple bitwidths with big integer are not supported!")
//		      if(endian == "mixed")
//		        stop("Cant't handle diverse bitwidths while endian is mixed: ", byte_order)
//
//		      bytes <- readBin(con=con, what="raw",n = nBytes, size = 1)
//		      # browser()
//		      if(dattype == "numeric" && length(unique(size)) > 1)
//		        stop("we don't support different bitwdiths for numeric data type!")
//		      dat <- convertRawBytes(bytes, isInt = dattype == "integer", colSize = size, ncol = nrpar, isBigEndian = endian == "big")
//
//
//		    }else{
//		      dat <- .readFCSdataRaw(con, dattype
//		                             , count= nBytes/size
//		                             , size= size
//		                             , signed=signed, endian=endian, splitInt = splitInt, byte_order = byte_order)
//		    }
//
//
//	    }else {
//	      if(multiSize)
//	        stop("'which.lines' can not be used when bitwidths are different across parameters!")
//	      ##Read n lines with or without sampling
//	        if(length(which.lines)==1)
//	            which.lines <- sample(seq_len(nrowTotal), which.lines)
//	        which.lines <- sort(which.lines)
//	        outrange <- length(which(which.lines > nrowTotal))
//	        if(outrange!=0)
//	            stop("Some or all the line indices specified are greater that the",
//	                 "number of collected events.\n")
//	        dat <- c()
//	        for (i in 1:length(which.lines)){
//	            startP <- offsets["datastart"] + (which.lines[i]-1) * nrpar * size
//	            endP   <-  startP + nrpar * size
//	            seek(con, startP)
//				temp <- .readFCSdataRaw(con, dattype, count= as.integer(endP - startP+1)/size,
//						size=size, signed=signed, endian=endian, splitInt = splitInt, byte_order = byte_order)
//
//	            dat <- c(dat, temp)
//	        }
//	    }
//	    ## stopifnot(length(dat)%%nrpar==0)
//	    ## Do we want the function to bail out when the above condition is TRUE?
//	    ## Might be better to assume the data was ok up to this point and
//	    ## exit gracefully with a warning as done in the following lines...
//	    ld <- length(dat)
//	    if(ld %% nrpar != 0){
//	        dat <- dat[1:(ld - (ld %% nrpar))]
//	        warning("Error in reading data stream for file '",
//	                summary(con)$description, "'/nData may be truncated!")
//	    }
//
//
//	    ## apply bitmask for integer data
//	    if(dattype=="integer"){
//	        if(length(unique(range))==1)
//	        {
//
//	            if(range[1]>0){
//	              usedBits <- ceiling(log2(range[1]))
//	              if(usedBits<bitwidth)
//	                dat <- dat %% (2^usedBits)
//	            }
//
//	            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
//	        }
//	        else
//	        {
//	            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
//	            for(i in 1:ncol(dat))
//	            {
//
//	                if(range[i] > 0){
//	                  usedBits <- ceiling(log2(range[i]))
//	                  if(usedBits<bitwidth_vec[i])
//	                      dat[,i] <- dat[,i] %% (2^usedBits)
//	                }
//	            }
//	        }
//	    }
//	    else
//	    {
//	        dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
//	    }
//
//	    cn  <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
//	    cn <- if(alter.names)  structure(make.names(cn),names=names(cn)) else cn
//	    dimnames(dat) <- list(NULL, cn)
//	    ## truncate data at max range
//	    if(is.na(x["transformation"]))
//	    {
//	        if(truncate_max_range){
//	          for(i in seq_len(ncol(dat)))
//	            dat[dat[,i]>range[i],i] <- range[i]
//	        }
//
//	        if(!is.null(min.limit))
//	            dat[dat<min.limit] <- min.limit
//	    }
//
//	    ## Transform or scale if necessary
//	    # J.Spidlen, Nov 13, 2013: added the flowCore_fcsPnGtransform keyword, which is
//	    # set to "linearize-with-PnG-scaling" when transformation="linearize-with-PnG-scaling"
//	    # in read.FCS(). This does linearization for log-stored parameters and also division by
//	    # gain ($PnG value) for linearly stored parameters. This is how the channel-to-scale
//	    # transformation should be done according to the FCS specification (and according to
//	    # Gating-ML 2.0), but lots of software tools are ignoring the $PnG division. I added it
//	    # so that it is only done when specifically asked for so that read.FCS remains backwards
//	    # compatible with previous versions.
//	    fcsPnGtransform <- FALSE
//	    flowCore_fcsPnGtransform <- readFCSgetPar(x, "flowCore_fcsPnGtransform", strict=FALSE)
//	    if(!is.na(flowCore_fcsPnGtransform) && flowCore_fcsPnGtransform == "linearize-with-PnG-scaling") fcsPnGtransform <- TRUE
//	    if(transformation)
//	    {
//	       ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""),
//	             strict=FALSE)
//	       noPnE <- is.na(ampliPar)
//	       if(any(noPnE))
//	       {
//	          warning("No '$PnE' keyword available for the following channels: ",
//	                paste(which(noPnE), collapse=", "), "\nUsing '0,0' as default.",
//	                call.=FALSE)
//	          ampliPar[noPnE] <- "0,0"
//	       }
//	       ampli <- do.call(rbind,lapply(ampliPar, function(x)
//	                   as.integer(unlist(strsplit(x,",")))))
//	       PnGPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "G", sep=""), strict=FALSE)
//	       noPnG <- is.na(PnGPar)
//	       if(any(noPnG)) PnGPar[noPnG] <- "1"
//	       PnGPar = as.numeric(PnGPar)
//
//	       for (i in 1:nrpar){
//	          if(ampli[i,1] > 0){
//	             # J.Spidlen, Nov 5, 2013: This was a very minor bug. The linearization transformation
//	             # for $PnE != "0,0" is defined as:
//	             # For $PnR/r/, r>0, $PnE/f,0/, f>0: n is a logarithmic parameter with channel values
//	             # from 0 to r-1. A channel value xc is converted to a scale value xs as xs=10^(f*xc/r).
//	             # Note the "r" instead of the "r-1" in the formula (which would admitedly make more sense)
//				 # However, this is the standard that apparently has been followed by BD and other companies
//	             # "forever" and it is therefore addoped as such by the ISAC DSTF (see FCS 3.1 specification)
//	             # To bring this to compliance, I am just changing
//	             # dat[,i] <- 10^((dat[,i]/(range[i]-1))*ampli[i,1])
//	             # to
//	             # dat[,i] <- 10^((dat[,i]/range[i])*ampli[i,1])
//	             dat[,i] <- 10^((dat[,i]/range[i])*ampli[i,1])
//	             range[i] <- 10^ampli[i,1]
//	          }
//	          else if (fcsPnGtransform && PnGPar[i] != 1) {
//	             dat[,i] <- dat[,i] / PnGPar[i]
//	             range[i] <- (range[i]-1) / PnGPar[i]
//	          }
//	          else
//	             range[i] <- range[i]-1
//
//	       }
//	    }
//	    if(scale){
//	        d = 10^decades
//	        for(i in 1:nrpar)
//	            if(ampli[i,1] > 0){
//	               dat[,i] <- d*((dat[,i]-1)/(range[i]-1))
//	               range[i] <- d*(range[i]/range[i]-1)
//	            } else{
//	                dat[,i] <- d*((dat[,i])/(range[i]))
//	                range[i] <- d
//	            }
//	    }
//	    attr(dat, "ranges") <- range
//	    return(dat)

}

void readHeaderAndText(ifstream &in,FCS_Header & header, KEY_WORDS & keys, bool isEmptyKeyValue, bool ignoreTextOffset,  int nDataset){


	 //search the stream for the header and txt of the nth DataSet
	int nOffset = 0, nNextdata = 0;
	int n = nDataset <=0?1:nDataset;
	 //non C-style index: starting from 1
	for(int i = 1; i <= n; i++)
	{
		nOffset += nNextdata;
		readFCSHeader(in,header, nOffset);//read the header
		readFCStext(in, header, keys, isEmptyKeyValue);//read the txt section

		if(keys.find("$NEXTDATA")!=keys.end())
			nNextdata = boost::lexical_cast<int>(keys["$NEXTDATA"]);
		else
		{
			if(i<n)
				throw(domain_error("Can't find " + boost::lexical_cast<string>(n) + "th dataset in FCS!"));
			break;
		}


	}


	  if(nDataset <=0 && nNextdata >0)
	  {
			cout << "The file contains additional data segment%s.\n The default is to read the first segment only.\nPlease consider setting the 'dataset' argument." << endl;

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
	       cout << "warning:Missing the required $BEGINDATA keyword! Reading data based on information in the FCS HEADER only." << endl;
	     } else
	       throw(domain_error("Don't know where the data segment begins, there was no $BEGINDATA keyword and the FCS HEADER does not say it either."));
	   }
	   else
		   datastart = boost::lexical_cast<int>(keys["$BEGINDATA"]);

	   if(keys.find("$ENDDATA")==keys.end())
	   {
		 if (dataend_h != 0) {
			 dataend = dataend_h;
		   cout << "warning:Missing the required $ENDDATA keyword! Reading data based on information in the FCS HEADER only." << endl;
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
			 if(ignoreTextOffset)
			   cout <<msg << " The values in TEXT are ignored!";
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
				if(ignoreTextOffset)
				   cout <<msg << " The values in TEXT are ignored!";
				else
				   throw(domain_error(msg));
			}
	   }

	}
}
MemCytoFrame::MemCytoFrame(const string &filename
		, bool isEmptyKeyValue = false, int nDataset = 1
		, bool scale =false, double decades=0, double min_limit=-111
		, bool truncate_max_range = true, bool ignoreTextOffset = false, bool onlyTxt = false){
	ifstream in(filename, ios::in|ios::binary);


	FCS_Header header;
	readHeaderAndText(in, header, CytoFrame::keys, isEmptyKeyValue, ignoreTextOffset, nDataset);



	if(onlyTxt)
		data = NULL;
	else
	{
		//parse the data section
		data = readFCSdata(in, header, CytoFrame::keys, scale, decades, min_limit, truncate_max_range);
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
