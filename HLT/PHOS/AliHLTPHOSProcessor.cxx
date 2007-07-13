#include "AliHLTPHOSProcessor.h"
//#include "AliHLTPHOSCommonDefs.h"
//#include "TString.h"

const AliHLTComponentDataType AliHLTPHOSProcessor::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array

AliHLTPHOSProcessor::AliHLTPHOSProcessor():AliHLTProcessor(), fkEquippmentID(0), 
					   fModuleID(0), fRcuX(0), fRcuZ(0),fRcuZOffset(0), fRcuXOffset(0),
					   fPrintInfo(kFALSE), fIsSetEquippmentID(kFALSE), fPrintInfoFrequncy(1000)
{

}


AliHLTPHOSProcessor::~AliHLTPHOSProcessor()
{

}


AliHLTPHOSProcessor::AliHLTPHOSProcessor(const AliHLTPHOSProcessor & ):AliHLTProcessor(), fkEquippmentID(0), 
								       fModuleID(0), fRcuX(0), fRcuZ(0),fRcuZOffset(0), fRcuXOffset(0),
								       fPrintInfo(kFALSE), fIsSetEquippmentID(kFALSE), fPrintInfoFrequncy(1000)
{

}

const AliHLTUInt16_t
AliHLTPHOSProcessor::GetEquippmentID() const
{
  return fkEquippmentID;
}

void 
AliHLTPHOSProcessor::SetEquippmentID(AliHLTUInt16_t id)
{
  AliHLTUInt16_t  &ref = const_cast<AliHLTUInt16_t&>(fkEquippmentID); 
  ref = id;
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
      iResult = -2; 
    }
  return iResult;
}


void 
AliHLTPHOSProcessor::SetCoordinates(AliHLTUInt16_t equippmentID)
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
