/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSProcessor.h"
#include "unistd.h"



const AliHLTComponentDataType AliHLTPHOSProcessor::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array


AliHLTPHOSProcessor::AliHLTPHOSProcessor():AliHLTProcessor(), 
					   AliHLTPHOSBase(), 
					   fPhosEventCount(0), 
					   fModuleID(2), 
					   fPrintInfoModule(0), 
					   fPrintInfoFrequncyModule(1000), 
					   fRunNumber(0)
{
  //  ScanRunNumberFromFile();
}



AliHLTPHOSProcessor::~AliHLTPHOSProcessor()
{

}

 
bool 
AliHLTPHOSProcessor::CheckFileLog(const char *origin,   const char *filename,  const char *opt)
{
  sprintf(fFilepath, "%s/%s", getenv("PWD"), filename);
  FILE *fp = fopen(filename, opt);

  if(fp == 0)
    {
      //      if( (opt == "w")  || (opt == "a")) \\OD
      if( (!strcmp(opt,"w"))  || (!strcmp(opt,"a")))
	{
	  sprintf(fMessage, "for writing  please check that the directory exists and that you have write access to it"  );
	}
      else
	{
	  sprintf(fMessage, "for reading  please check that the directory exists and that you have read access to it"  );
	}
     Logging(kHLTLogFatal, origin , "cannot open file" , "Was not able to open file %s  %s", fFilepath, fMessage);
     return false;
    }
  else
    {
      //      if( (opt == "w")  || (opt == "a")) \\OD
      if( (!strcmp(opt,"w"))  || (!strcmp(opt,"a")))
	{
	  sprintf(fMessage, "for writing" );
	}
      else
	{
	  sprintf(fMessage, "for reading");
	}
      //    Logging(kHLTLogInfo, origin , "opening file" , "Sucessfully opening %s  %s", fFilepath, fMessage);
      fclose(fp); 
      return true;
    }
  
}

 
void 
AliHLTPHOSProcessor::DoneWritingLog(const char *origin, const char *filename)
{
  //  char filepath[1024];
  sprintf(fFilepath, "%s/%s", getenv("PWD"), filename);
  Logging(kHLTLogInfo, origin , "finnished writing file" , "wrote file %s", fFilepath);
}


void 
AliHLTPHOSProcessor::ScanRunNumberFromFile()
{
  char tmpDirectory[512];
  char tmpFileName[512];
  sprintf(tmpDirectory, "%s", getenv("HOME"));  

  //TODO, remove hardcoded file path
 
  
  sprintf(tmpFileName, "%s%s", tmpDirectory, "/hlt/rundir/runNumber.txt");
  int tmp = 0;
  if(CheckFileLog( __FILE__ , tmpFileName , "r")== true) 
    { 
      FILE *fp = fopen(tmpFileName, "r");
      tmp = fscanf(fp, "%d", &fRunNumber);
      fclose(fp);
    }

  


 }

const char*
AliHLTPHOSProcessor::IntToChar(int number)
{
  sprintf(lineNumber,"%d", number);
  return lineNumber;
}



int
AliHLTPHOSProcessor::ScanArgumentsModule(int argc, const char** argv)
{
  fPrintInfoModule = kFALSE;
  int iResult=0;
  TString argument="";

  for(int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      
      if (argument.IsNull()) 
	{
	  continue;
	}
                         
    if (argument.CompareTo("-printinfo") == 0) 
      {
	if(i+1 <= argc)
	  {
	    argument=argv[i+1];
	    fPrintInfoFrequncyModule = atoi(argv[i+1]);
	    fPrintInfoModule = kTRUE;
	  }
	else
	  {
	    Logging(kHLTLogWarning, __FILE__ , "inconsistency during init" , "asking for event info, but no update frequency is specified, option is ignored");
	  }
      }
 
    }
  return 0;
}

