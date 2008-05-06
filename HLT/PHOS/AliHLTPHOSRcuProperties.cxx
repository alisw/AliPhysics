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

#include "AliHLTPHOSRcuProperties.h"
//#include "AliHLTPHOSConstants.h"


//#include ""


AliHLTPHOSRcuProperties::AliHLTPHOSRcuProperties() :fTest(-2),
						    fkEquippmentID(0),
						    fModID(2), 
						    fRcuID(0),
						    fRcuX(0), 
						    fRcuZ(0),
						    fRcuZOffset(0),
						    fRcuXOffset(0),
						    fPrintInfo(false),
						    fPrintInfoFrequncy(10000),
						    fIsSetEquippmentID(kFALSE),
						    fLog(),
						    fIsInitialized(false)
{

}

AliHLTPHOSRcuProperties::~AliHLTPHOSRcuProperties()
{

}


const AliHLTUInt16_t
AliHLTPHOSRcuProperties::GetEquippmentID() const
{
  return fkEquippmentID;
}


const int 
AliHLTPHOSRcuProperties::GetRCUID() const
{
  cout << "AliHLTPHOSRcuProperties::GetRCUI returning " <<  fRcuID  << endl;
  if(fIsInitialized == true)
    {
      cout <<"AliHLTPHOSRcuProperties::GetRCUID fIsInitialized = TRUE"  << endl;
    }
  else
    {
      cout <<"AliHLTPHOSRcuProperties::GetRCUID fIsInitialized = FALSE"  << endl;
    }

 return fRcuID;
  // cout << "AliHLTPHOSRcuProperties::GetRCUI returning " <<  fRcuID  << endl;

  // int tmpRcuId = 0;


  /*
  if(fRcuZ ==0 && fRcuX == 0)
    {
      tmpRcuId = 0; 
    }

  if(fRcuZ ==1 && fRcuX == 0)
    {
      tmpRcuId = 1; 
    }
  
  if(fRcuZ ==1 && fRcuX == 1)
    {
      tmpRcuId = 2; 
    }

    if(fRcuZ ==0 && fRcuX == 1)
    {
      tmpRcuId = 3; 
    }
  */   

  //  return tmpRcuId;
}



void 
AliHLTPHOSRcuProperties::SetEquippmentID(AliHLTUInt16_t id)
{
  AliHLTUInt16_t  &ref = const_cast<AliHLTUInt16_t&>(fkEquippmentID); 
  ref = id;
}


int
AliHLTPHOSRcuProperties::ScanArguments(int argc, const char** argv)
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
	  cout << "AliHLTPHOSRcuProperties:ScanArguments  argument = -equipmentID   "  <<endl;  
	  if(i+1 <= argc)
	    {
	      SetEquippmentID((AliHLTUInt16_t)atoi(argv[i+1]));
	      cout <<"AliHLTPHOSRcuProperties:ScanArguments  equipmentID =  " << fkEquippmentID  <<endl;
	      fLog.Logging(kHLTLogInfo, __FILE__ , "Init info", "  setting equippment ID to  %lu ", fkEquippmentID); 
	      SetCoordinates(fkEquippmentID);
	      fIsSetEquippmentID = kTRUE;
	    }
	  else
	    {
	       iResult= -1;
	       fLog.Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuProperties::ScanArguments( int argc, const char** argv )", "Missing argument",
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
	    fLog.Logging(kHLTLogInfo, __FILE__ , "Info output", " setting printinfo = kTRUE, with update frequency every %lu th event ", fPrintInfoFrequncy);
	  }
	else
	  {
	    //	    cout << "WARNING: asking for event info, but no update frequency is specified, option is ignored" << endl;
	    fLog.Logging(kHLTLogWarning, __FILE__ , "Invalid request", " asking for event info, but no update frequency is specified, request  ignored");
	  }
      }
 
    }


  if(fIsSetEquippmentID == kFALSE)
    {
      fLog.Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuProperties::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -3; 
    }
  return iResult;
}


void 
AliHLTPHOSRcuProperties::SetCoordinates(AliHLTUInt16_t /*equippmentID*/)
{
  // int rcuIndex =  (fkEquippmentID - 1792)%N_RCUS_PER_MODULE;

  fRcuID =  (fkEquippmentID - 1792)%N_RCUS_PER_MODULE;
  //  fModID = (fkEquippmentID  -1792 -rcuIndex)/N_RCUS_PER_MODULE;
  fModID = (fkEquippmentID  -1792 - fRcuID)/N_RCUS_PER_MODULE;
 
  if( fRcuID  == 0)
    {
      fRcuX = 0; 
      fRcuZ = 0;
    }

  if( fRcuID  == 1)
    {
      fRcuX = 0; 
      fRcuZ = 1;
    }
 
  if( fRcuID == 2)
    {
      fRcuX = 1; 
      fRcuZ = 0;
    }

  if( fRcuID == 3)
    {
      fRcuX = 1; 
      fRcuZ = 1;
    }

  fRcuZOffset =  N_ZROWS_RCU*fRcuZ;
  fRcuXOffset =  N_XCOLUMNS_RCU*fRcuX;
  
  // cout <<"********InitInfo************"<< endl;
//   cout <<"AliHLTPHOSRcuProperties::SetCoordinate casted"<< endl;
//   cout <<"rcuIndex"<< fRcuID  <<endl;
//   cout <<"Equpippment ID =\t"<< fkEquippmentID <<endl;
//   cout <<"Mod ID =\t"<<  (int)fModID<<endl;
//   cout <<"RCUX =\t\t" << (int)fRcuX << endl;
//   cout <<"RCUZ =\t\t" << (int)fRcuZ << endl;
//   cout <<"RcuZOffset = \t" <<  (int)fRcuZOffset << endl;
//   cout <<"RcuXOffset = \t" <<  (int)fRcuXOffset << endl << endl;
  fIsInitialized = true;
}

