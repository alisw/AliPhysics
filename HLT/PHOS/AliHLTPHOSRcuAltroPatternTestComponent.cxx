// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/




#include "AliHLTPHOSRcuAltroPatternTestComponent.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSValidCellDataStruct.h" 
#include "AliHLTPHOSRcuAltroPatternTest.h"

AliHLTPHOSRcuAltroPatternTestComponent gAliHLTPHOSRcuAltroPatternTestComponent;

AliHLTPHOSRcuAltroPatternTestComponent:: AliHLTPHOSRcuAltroPatternTestComponent() : AliHLTPHOSRcuProcessor(), 
										    fPatternTestPtr(0),
										    fShmPtr(0),
										    fNTotalPatterns(0),
										    fNWrongPatterns(0),
										    fNTotalSamples(0), 
										    fNWrongSamples(0)
{
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
} 


AliHLTPHOSRcuAltroPatternTestComponent::~ AliHLTPHOSRcuAltroPatternTestComponent()
{
  //Destructor
}


int 
AliHLTPHOSRcuAltroPatternTestComponent::Deinit()
{
  //See html documentation of base class
  return 0;
}


const char* 
AliHLTPHOSRcuAltroPatternTestComponent::GetComponentID()
{
  //See html documentation of base class
  return "AltroPatternTester";
}


void
AliHLTPHOSRcuAltroPatternTestComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See html documentation of base class
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) 
    {
      list.push_back(*pType);
      pType++;
    }
}


AliHLTComponentDataType 
AliHLTPHOSRcuAltroPatternTestComponent::GetOutputDataType()
{
  //See html documentation of base class  
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}


void
AliHLTPHOSRcuAltroPatternTestComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //See html documentation of base class
  constBase = 30;
  inputMultiplier = 1;
}



int  AliHLTPHOSRcuAltroPatternTestComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
						      AliHLTComponentTriggerData& /*trigData*/,  AliHLTUInt8_t* outputPtr, 
						      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& /*outputBlocks */)
{
  //See html documentation of base class
  //cout << "AliHLTPHOSRcuAltroPatternTestComponent::DoEvent, processing event "<< fPhosEventCount  << endl;

  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  unsigned long ndx       = 0;
  //  UInt_t offset           = 0; 
  //  UInt_t mysize           = 0;
  UInt_t tSize            = 0;
  const AliHLTComponentBlockData* iter = NULL;   
  AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr;
  AliHLTUInt8_t* outBPtr;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkCellEnergyDataType)
	{
	  continue;
	}
      
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      fShmPtr->SetMemory(cellDataPtr);
      currentChannel = fShmPtr->NextChannel();
      
      Int_t* tmpRawPtr = 0;
     
      int tmp =0;      
      while(currentChannel != 0)
	{

	  int tmpZ =  currentChannel->fZ;
	  int tmpX =  currentChannel->fX;
	  int tmpGain =  currentChannel->fGain;
	  int tmpSamples = 0;
	  tmpRawPtr = fShmPtr->GetRawData(tmpSamples); 
	  
	  
	  if( (tmpZ > NZROWSRCU) || (tmpX > NXCOLUMNSRCU) || (tmpGain > NGAINS))
	    {
	      cout <<" ERROR parameters out of range z = "<< tmpZ <<" x = "<< tmpX<< " gain = " << tmpGain<<" nSamples = " <<  tmpSamples  <<endl;
	    }
	  else
	    {
	      //    cout  <<"  all parameters in range" << endl;
	    }
	  
	  //	  cout << "analyzing channelnr " << tmp <<" of event " << fPhosEventCount << endl;
	  //	  cout << "  z = " << currentChannel->fZ << "  x = "  << currentChannel->fX  << "  gain = " <<  currentChannel->fGain  << "  samples =" << currentChannel->fNSamples  <<endl;
	  tmp ++;

	  fPatternTestPtr->AddPattern(tmpRawPtr, currentChannel->fZ, currentChannel->fX,  currentChannel->fGain, tmpSamples, fNPresamples);
	  int ret =  fPatternTestPtr->ValidateAltroPattern(tmpRawPtr, tmpSamples);
	  
	  if(ret >= 0)
	    {
	      
	      fNTotalSamples += tmpSamples;
	      fNTotalPatterns ++;   

	      if(ret > 0)
		{
		  fNWrongSamples += ret;
		  fNWrongPatterns ++;
		  
		  if(fPhosEventCount%100 == 0)
		    { 
		      /*
		      cout << "warning: incorrect pattern found for event " << fPhosEventCount <<" X = " << currentChannel->fX <<"  Z = " << currentChannel->fZ ;
		      cout << "number of correct aptterns is " << fNTotalPatterns  << "number of wrong patterns is " << fNWrongPatterns<<endl;
		      float percent = 100*((float)fNWrongPatterns)/( float(fNTotalPatterns));  
		      cout <<"The corrpution rate is " << percent << " percent " <<endl;
		      cout <<" ERROR: incorrect parameters" << endl;
		      */
		    }
		

		  if(fPhosEventCount%100 == 0 && fPhosEventCount !=0)
		    {
		      fPatternTestPtr->countAllPatterns(fNSamples);  
		      fPatternTestPtr->PrintStatistics();
		      
		    }
		}
	    }
	  else
	    {
	      cout <<" ERROR: incorrect parameters" << endl;
	    }
	  currentChannel = fShmPtr->NextChannel();
	  
	}
    }

 
  outBPtr = outputPtr;
 
  //  mysize += sizeof(AliHLTPHOSRcuCellAccumulatedEnergyDataStruct);

  /*
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType;
  bd.fSpecification = 0xFFFFFFFF;
  outputBlocks.push_back( bd );
  tSize += mysize;
  outBPtr += mysize;
  */

  if( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuAltroPatternTestComponent::DoEvent", "Too much data",
	       "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
	       , tSize, size );
      return EMSGSIZE;
    }

  fPhosEventCount++; 

 
  return 0;
}//end DoEvent





int
AliHLTPHOSRcuAltroPatternTestComponent::DoInit(int argc, const char** argv )
{
  //See html documentation of base class
  cout << "AliHLTPHOSRcuAltroPatternTestComponent::DoInit, argc =" << argc << endl;
  char patternFilename[2000];
  bool isSetPartternArg = false;
  fPrintInfo = kFALSE;
  int iResult=0;
  TString argument="";
  iResult = ScanArguments(argc, argv);
  //  fPatternTestPtr = new  AliHLTPHOSRcuAltroPatternTest(fModuleID, fRcuX, fRcuZ);

  if(iResult < 0)
    {
      return iResult;
    }

  if(fIsSetEquippmentID == kFALSE)
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuAltroPatternTestComponent::DoInit( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      //  iResult = -2; 
      return -5;
    }

  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-patternfile", argv[i]))
	{
	  if(argc < (i+1))
	    {
	      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuAltroPatternTestComponent::DoInit( int argc, const char** argv )", "Missing argument",
		       "The argument  -patternfile requires that a filename follows: set it with a component argumet like this:  -patternfile   <filename>");
	      return -6;  
	    }
	  else
	    {
	      sprintf(patternFilename, "%s", argv[i+1]);
	      isSetPartternArg = true;
	    }
	}
    }

  if(isSetPartternArg == false)
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuAltroPatternTestComponent::DoInit( int argc, const char** argv )", "Missing argument",
	       "The argument patternfile is not set: set it with a component argumet like this: -patternfile  <filename>");
      return -7;
      
    }
  else  if(fUtilitiesPtr->CheckFile(patternFilename, "r") == false)
    {
      char tmpMessage[1024];
      sprintf(tmpMessage, "the file %s could not be found, please check that file exist and that you have read access to it", patternFilename);
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuAltroPatternTestComponent::DoInit( int argc, const char** argv )", "File not found",tmpMessage);

     return -8;
     
    }
  else
    {
      int tmpPattern[ALTROMAXSAMPLES]; 
      ScanPatternFromFile(patternFilename, tmpPattern, ALTROMAXSAMPLES) ;
      fPatternTestPtr = new  AliHLTPHOSRcuAltroPatternTest(fModuleID, fRcuX, fRcuZ, tmpPattern, ALTROMAXSAMPLES);
    }

  //  return iResult; 
  return 0;
}


AliHLTComponent*
AliHLTPHOSRcuAltroPatternTestComponent::Spawn()
{
  //See html documentation of base class
  return new AliHLTPHOSRcuAltroPatternTestComponent;
}


void 
AliHLTPHOSRcuAltroPatternTestComponent::ScanPatternFromFile(const char *filename, int *pattern, const int /*length*/) const
{
  FILE *fp = fopen(filename, "r");

  //  int tmpPattern[ALTRO_MAX_SAMPLES];
  int dummy = 0;
  int res = 0;
  for(int i=0; i<ALTROMAXSAMPLES; i++)
    {
       res = fscanf(fp,"w 0x%X 0x%X\n", &dummy, &pattern[i]);
      //      cout << tmpPattern[i] << endl;
    }

  
  //  fPatternTestPtr->SetAltroPattern(tmpPattern);
}
  
//w 0x6801 0x1
