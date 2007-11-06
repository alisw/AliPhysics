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

#include "AliHLTPHOSDefinitions.h"
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
#include "AliHLTPHOSDDLPackedFileWriter.h" 
#include "AliHLTPHOSCellEnergiesFileWriter.h"


const AliHLTComponentDataType AliHLTPHOSFileWriterComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array


AliHLTPHOSFileWriterComponent gAliHLTPHOSFileWriterComponent;

//____________________________________________________________________________________
AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent():  AliHLTFileWriter(), fCellEnergiesFileWriterPtr(0) \
  , fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEvtCnt(0)
{
  /* 
   *Default constructor
   */
  for(int i=0; i<N_DATATYPES; i++)
    {
      fDataTypesToFile[i] = kAliHLTVoidDataType;
    }

  fCellEnergiesFileWriterPtr = new AliHLTPHOSCellEnergiesFileWriter();
  fDDLPackedFileWriterPtr    = new AliHLTPHOSDDLPackedFileWriter();

} 

//____________________________________________________________________________________
AliHLTPHOSFileWriterComponent::~AliHLTPHOSFileWriterComponent()
{
  //comment
  delete fCellEnergiesFileWriterPtr;
  delete fDDLPackedFileWriterPtr;
}


//____________________________________________________________________________________
AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent(const AliHLTPHOSFileWriterComponent & ): AliHLTFileWriter(), fCellEnergiesFileWriterPtr(0), \
  fDDLPackedFileWriterPtr(0), fDirectory(""),fFilename(""), fEvtCnt(0)
{
  //comment
}


//____________________________________________________________________________________
int 
AliHLTPHOSFileWriterComponent::AddDataType(string dataType)
{
  //comment
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
      fDataTypesToFile[tmpCnt] = AliHLTPHOSDefinitions::fgkCellEnergyDataType; 
      cout <<"regsitring dataType for filewriting: fDataTypesToFile[" << tmpCnt <<"]"<<endl; 
    } 
  else if(dataType.compare("gkDDLPackedRawDataType") == 0)
    {
 	  fDataTypesToFile[tmpCnt] = AliHLTPHOSDefinitions::fgkDDLPackedRawDataType; 
    }

  cout << "dataType.compare(cmpString) = " <<dataType.compare(cmpString)<<endl;
  return ret;
}

//____________________________________________________________________________________
int 
AliHLTPHOSFileWriterComponent::Deinit()
{
  //comment
  return 0;
}

//____________________________________________________________________________________
int 
AliHLTPHOSFileWriterComponent::DoDeinit()
{
  //comment
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSFileWriterComponen DoDeinit");
  return 0;
}


//____________________________________________________________________________________
const char* 
AliHLTPHOSFileWriterComponent::GetComponentID()
{
  //comment
  return "PhosFileWriter";
}


//____________________________________________________________________________________
AliHLTComponent*
AliHLTPHOSFileWriterComponent::Spawn()
{  
  //comment
  return new AliHLTPHOSFileWriterComponent;
}


//____________________________________________________________________________________
void
AliHLTPHOSFileWriterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

//____________________________________________________________________________________
AliHLTComponentDataType 
AliHLTPHOSFileWriterComponent::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}

//____________________________________________________________________________________
void
AliHLTPHOSFileWriterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  // see header file for documentation
  constBase = 30;
  inputMultiplier = 0.1;
}

//____________________________________________________________________________________
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
