/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSRcuProcessor.h"

AliHLTPHOSRcuProcessor::AliHLTPHOSRcuProcessor(): fkEquippmentID(0)
{

}

AliHLTPHOSRcuProcessor::~AliHLTPHOSRcuProcessor()
{

}


const AliHLTUInt16_t
AliHLTPHOSRcuProcessor::GetEquippmentID() const
{
  return fkEquippmentID;
}

void 
AliHLTPHOSRcuProcessor::SetEquippmentID(AliHLTUInt16_t id)
{
  AliHLTUInt16_t  &ref = const_cast<AliHLTUInt16_t&>(fkEquippmentID); 
  ref = id;
}


int
AliHLTPHOSRcuProcessor::ScanArguments(int argc, const char** argv)
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
                         
    if (argument.CompareTo("-equipmentID") == 0) 
	{
	  cout << "AliHLTPHOSProcessor:DoInit  argument = -equipmentID   "  <<endl;  
	  if(i+1 <= argc)
	    {
	      SetEquippmentID((AliHLTUInt16_t)atoi(argv[i+1]));
	      cout << "AliHLTPHOSRawAnalyzerComponent:DoInit  setting equippment ID to  " << fkEquippmentID <<endl;
	      SetCoordinates(fkEquippmentID);
	      fIsSetEquippmentID = kTRUE;
	      cout << " fIsSetEquippmentID = kTRUE"<< endl;
	    }
	  else
	    {
	       iResult= -1;
	       Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
			"The argument -equippmentID expects a number");
	       return  iResult;   
	    }
	}
    
    
    if (argument.CompareTo("-printinfo") == 0) 
      {
	if(i+1 <= argc)
	  {
	    argument=argv[i+1];
	    fPrintInfoFrequncy = atoi(argv[i+1]);
	    fPrintInfo = kTRUE;
	    cout << "AliHLTPHOSRawAnalyzerComponent::DoIni  setting printinfo = kTRUE, with update frequency every  "<< fPrintInfoFrequncy << "th event" <<endl; 
	  }
	else
	  {
	    cout << "WARNING: asking for event info, but no update frequency is specified, option is ignored" << endl;
	  }
      }
 
    }


  if(fIsSetEquippmentID == kFALSE)
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -3; 
    }
  return iResult;
}

void 
AliHLTPHOSRcuProcessor::SetCoordinates(AliHLTUInt16_t /*equippmentID*/)
{
  int rcuIndex =  (fkEquippmentID - 1792)%N_RCUS_PER_MODULE;
  fModuleID = (fkEquippmentID  -1792 -rcuIndex)/N_RCUS_PER_MODULE;
  
  if(rcuIndex == 0)
    {
      fRcuX = 0; 
      fRcuZ = 0;
    }

  if(rcuIndex == 1)
    {
      fRcuX = 0; 
      fRcuZ = 1;
    }
 
  if(rcuIndex == 2)
    {
      fRcuX = 1; 
      fRcuZ = 0;
    }

  if(rcuIndex == 3)
    {
      fRcuX = 1; 
      fRcuZ = 1;
    }

  fRcuZOffset =  N_ZROWS_RCU*fRcuZ;
  fRcuXOffset =  N_XCOLUMNS_RCU*fRcuX;

  cout <<"********InitInfo************"<< endl;
  cout <<"AliHLTPHOSRawAnalyzerComponent::SetCoordinate casted"<< endl;
  cout <<"Equpippment ID =\t"<< fkEquippmentID <<endl;
  cout <<"Module ID =\t"<<  (int)fModuleID<<endl;
  cout <<"RCUX =\t\t" << (int)fRcuX << endl;
  cout <<"RCUZ =\t\t" << (int)fRcuZ << endl;
  cout <<"RcuZOffset = \t" <<  (int)fRcuZOffset << endl;
  cout <<"RcuXOffset = \t" <<  (int)fRcuXOffset << endl << endl;

}
