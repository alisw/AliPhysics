/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//#include "AliHLTPHOSFileWriterComponent.h"


#include "AliHLTPHOSFileWriterComponent.h"
#include <iostream>
#include "stdio.h"
#include <TObjString.h>

#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDataHeaderStruct.h"
#include "AliHLTDataTypes.h"


const AliHLTComponentDataType AliHLTPHOSFileWriterComponent::fInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array


AliHLTPHOSFileWriterComponent gAliHLTPHOSFileWriterComponent;

AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent():AliHLTDataSink(), fCellEnergiesFileWriterPtr(0), fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEventCount(0)
  //AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent():AliHLTDataSink()
{
  for(int i=0; i<N_DATATYPES; i++)
    {
      fDataTypesToFile[i] = kAliHLTVoidDataType;
    }

  fCellEnergiesFileWriterPtr = new AliHLTPHOSCellEnergiesFileWriter();
  fDDLPackedFileWriterPtr    = new AliHLTPHOSDDLPackedFileWriter();

} 

AliHLTPHOSFileWriterComponent::~AliHLTPHOSFileWriterComponent()
{
  delete fCellEnergiesFileWriterPtr;
  delete fDDLPackedFileWriterPtr;
}


AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent(const AliHLTPHOSFileWriterComponent & ):AliHLTDataSink(), fCellEnergiesFileWriterPtr(0), fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEventCount(0)
  //AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent(const AliHLTPHOSFileWriterComponent & ):AliHLTDataSink()
{

}

int 
AliHLTPHOSFileWriterComponent::AddDataType(string dataType)
{
  int ret = -1;
  int tmpCnt = 0;
  for(int i=0; i< N_DATATYPES; i++)
    {
      if( fDataTypesToFile[i] != kAliHLTVoidDataType)
	{
	  tmpCnt ++;
	}
    }
  
  string cmpString("gkCellEnergyDataType");

  if(dataType.compare("gkCellEnergyDataType") == 0)
    {
      cout << "!!!!!!!!!!!!!!AliHLTPHOSFileWriterComponent::AddDataType"<< dataType << endl;
      //    fDataTypesToFilePtr[N_DATATYPES] = new AliHLTPHOSDefinitions::gkCellEnergyDataType;
      fDataTypesToFile[tmpCnt] = AliHLTPHOSDefinitions::gkCellEnergyDataType; 
      cout <<"regsitring dataType for filewriting: fDataTypesToFile[" << tmpCnt <<"]"<<endl; 
    } 
  else if(dataType.compare("gkDDLPackedRawDataType") == 0)
    {
 	  fDataTypesToFile[tmpCnt] = AliHLTPHOSDefinitions::gkDDLPackedRawDataType; 
    }

  cout << "dataType.compare(cmpString) = " <<dataType.compare(cmpString)<<endl;
  return ret;
}




int 
AliHLTPHOSFileWriterComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSFileWriterComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSFileWriterComponen DoDeinit");
  return 0;
}

const char* 
AliHLTPHOSFileWriterComponent::GetComponentID()
{
  return "PhosFileWriter";
}

AliHLTComponent*
AliHLTPHOSFileWriterComponent::Spawn()
{
  return new AliHLTPHOSFileWriterComponent;
}


void
AliHLTPHOSFileWriterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=fInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSFileWriterComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}

void
AliHLTPHOSFileWriterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 0.1;
}

int 
AliHLTPHOSFileWriterComponent::DumpEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData )
{
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;
  Int_t tmpChannelCnt     = 0;
  UInt_t offset           = 0;
  const AliHLTComponentDataType *tmpDataType; 
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  AliHLTPHOSDataHeaderStruct  dataHeader;
  
  dataHeader.fSize = sizeof(dataHeader); 
  dataHeader.fEventID =  evtData.fEventID;
  cout << "analyzing event: " << fEventCount << endl;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;
      tmpDataType = &(iter->fDataType);

      if(IsRegisteredDataType(*tmpDataType))
	{
	  if(*tmpDataType == AliHLTPHOSDefinitions::gkCellEnergyDataType)
	    {

	      fCellEnergiesFileWriterPtr->WriteFile(evtData, blocks, trigData, fEventCount); 
	      cout <<"AliHLTPHOSFileWriterComponen: data type = is  gkCellEnergyDataType. block index = "<< ndx\
		   <<" EventCount =" << fEventCount  << "Event ID"<<evtData.fEventID << endl;
	    }
	  else if(*tmpDataType == AliHLTPHOSDefinitions::gkDDLPackedRawDataType)
	    {

	       fDDLPackedFileWriterPtr->WriteFile(evtData, blocks, trigData, fEventCount); 
	    }
	}

      //  cout <<"AliHLTPHOSFileWriterComponen: data type = is  gkCellEnergyDataType. block index = "<< ndx\
	  //		   <<" EventCount =" << fEventCount  << "Event ID"<<evtData.fEventID << endl;

      fEventCount++;
    } 
  return 0;
}//end DumpEvent


int
AliHLTPHOSFileWriterComponent::DoInit( int argc, const char** argv )
{
  int iResult=0;
  TString argument="";
  Bool_t dirSet = kFALSE;
  Bool_t dataSet = kFALSE;
  string dataType="";
  int bMissingParam=0;

  //  fFilename.assign(256,0);
  //   fFilename.assign(256,0);

  for(int i=0; i<argc; i++)
    {
      argument=argv[i];
 
      if(argument.CompareTo("-directory")==0) 
	{
	  if ((bMissingParam=(++i>=argc))) 
	    {
	      break;
	    }
	  fDirectory.assign(argv[i]);
	  fCellEnergiesFileWriterPtr->SetDirectory(fDirectory);
	  fDDLPackedFileWriterPtr->SetDirectory(fDirectory) ;

	  fFilename.insert(0, fDirectory);
	  dirSet = kTRUE;
	  
	  cout << "fDirectory=" << fDirectory << endl;
	}

      if(argument.CompareTo("-datatype")==0) 
	{
	  if ((bMissingParam=(++i>=argc))) break;
	  cout << "datatype = " << argv[i] << endl;
	  dataType = argv[i];
	  
	  AddDataType(dataType);
	  dataSet = kTRUE;
	}
      
      cout << "argv[" << i <<"] = " << argv[i] << endl;

    }

  /*
   * We dont start the component if we don know what data to write
   * or where to store it
   */
  if((dataSet != kTRUE || dataSet != kTRUE))
    {
      iResult = -1;
      HLTFatal(" either direcory or datatype is not set, usage -datatype <datatype>  -driectory <directory>");
    }
  else
    {
      iResult = 0;
    }
  return iResult;
}


Bool_t 
AliHLTPHOSFileWriterComponent::IsRegisteredDataType(const AliHLTComponentDataType& dataType)
{
  Bool_t tmp = kFALSE;
  for(int i =0; i<N_DATATYPES; i++)
    {
      if((fDataTypesToFile[i] == dataType) && (dataType !=  kAliHLTVoidDataType))
	{
	  tmp = kTRUE;
	}
    }

  return tmp;
}

