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

#include "AliHLTPHOSFileWriterComponent.h"
#include <iostream>
#include <TObjString.h>
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDataHeaderStruct.h"
#include "AliHLTDataTypes.h"



class AliRawReaderMemory;
class AliCaloRawStream;
class AliHLTPHOSRcuCellEnergyDataStruct;


const AliHLTComponentDataType AliHLTPHOSFileWriterComponent::fInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array


AliHLTPHOSFileWriterComponent gAliHLTPHOSFileWriterComponent;

AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent():  AliHLTFileWriter(), fCellEnergiesFileWriterPtr(0) \
  , fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEventCount(0)
{
  // see header file for documentation
  for(int i=0; i<N_DATATYPES; i++)
    {
      fDataTypesToFile[i] = kAliHLTVoidDataType;
    }

  fCellEnergiesFileWriterPtr = new AliHLTPHOSCellEnergiesFileWriter();
  fDDLPackedFileWriterPtr    = new AliHLTPHOSDDLPackedFileWriter();

} 

AliHLTPHOSFileWriterComponent::~AliHLTPHOSFileWriterComponent()
{
  // see header file for documentation
  delete fCellEnergiesFileWriterPtr;
  delete fDDLPackedFileWriterPtr;
}



AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent(const AliHLTPHOSFileWriterComponent & ): AliHLTFileWriter(), fCellEnergiesFileWriterPtr(0), \
  fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEventCount(0)
{
  // see header file for documentation

}

int 
AliHLTPHOSFileWriterComponent::AddDataType(string dataType)
{
  // see header file for documentation 
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
  // see header file for documentation
  return 0;
}

int 
AliHLTPHOSFileWriterComponent::DoDeinit()
{
 // see header file for documentation
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSFileWriterComponen DoDeinit");
  return 0;
}

const char* 
AliHLTPHOSFileWriterComponent::GetComponentID()
{
 // see header file for documentation
  return "PhosFileWriter";
}

AliHLTComponent*
AliHLTPHOSFileWriterComponent::Spawn()
{
  // see header file for documentation
  return new AliHLTPHOSFileWriterComponent;
}


void
AliHLTPHOSFileWriterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for documentation
  const AliHLTComponentDataType* pType=fInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSFileWriterComponent::GetOutputDataType()
{
  // see header file for documentation
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}

void
AliHLTPHOSFileWriterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  // see header file for documentation
  constBase = 30;
  inputMultiplier = 0.1;
}

Bool_t 
AliHLTPHOSFileWriterComponent::IsRegisteredDataType(const AliHLTComponentDataType& dataType)
{
  // see header file for documentation
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
