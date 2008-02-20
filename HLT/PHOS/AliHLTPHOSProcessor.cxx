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


AliHLTPHOSProcessor::AliHLTPHOSProcessor():AliHLTProcessor(), AliHLTPHOSBase(), fPhosEventCount(0), fModuleID(0), fPrintInfo(0), fPrintInfoFrequncy(1000), fRunNumber(0)
{
  ScanRunNumberFromFile();
}



AliHLTPHOSProcessor::~AliHLTPHOSProcessor()
{

}


void 
AliHLTPHOSProcessor::ScanRunNumberFromFile()
{
  char tmpDirectory[512];
  char tmpFileName[512];
  sprintf(tmpDirectory, "%s", getenv("HOME"));  
  sprintf(tmpFileName, "%s%s", tmpDirectory, "/rundir/runNumber.txt");

  if(CheckFile(tmpFileName, "r") == true)
    {
      FILE *fp = fopen(tmpFileName, "r");
      fscanf(fp, "%d", &fRunNumber);
      fclose(fp);
    }

  else
    {
      cout << "ERROR, could not find file  " << tmpFileName <<endl;
    }
}
int
AliHLTPHOSProcessor::ScanArguments(int argc, const char** argv)
{
  fPrintInfo = kFALSE;
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
	    fPrintInfoFrequncy = atoi(argv[i+1]);
	    fPrintInfo = kTRUE;
	  }
	else
	  {
	    cout << "WARNING: asking for event info, but no update frequency is specified, option is ignored" << endl;
	  }
      }
 
    }
  return 0;
}
