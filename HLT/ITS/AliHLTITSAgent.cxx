// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTITSAgent.cxx
    @author Matthias Richter
    @date   25.08.2008
    @brief  Agent of the libAliHLTITS library
*/

#include <cassert>
#include "AliHLTITSAgent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTOUT.h"
//#include "AliDAQ.h"

// header files of library components

// header file of the module preprocessor
#include "AliHLTITSCompressRawDataSDDComponent.h"
#include "AliHLTITSClusterFinderSPDComponent.h"
#include "AliHLTITSClusterFinderSSDComponent.h"

/** global instance for agent registration */
AliHLTITSAgent gAliHLTITSAgent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSAgent)

AliHLTITSAgent::AliHLTITSAgent()
  :
  AliHLTModuleAgent("ITS"),
  fRawDataHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTITSAgent::~AliHLTITSAgent()
{
  // see header file for class documentation
}

int AliHLTITSAgent::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,
					    AliRawReader* /*rawReader*/,
					    AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return 0;
}

const char* AliHLTITSAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						       AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation

  return "";
}

const char* AliHLTITSAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return "";
}

int AliHLTITSAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTITSCompressRawDataSDDComponent);
  pHandler->AddComponent(new AliHLTITSClusterFinderSPDComponent);
  pHandler->AddComponent(new AliHLTITSClusterFinderSSDComponent);

  return 0;
}

AliHLTModulePreprocessor* AliHLTITSAgent::GetPreprocessor()
{
  // see header file for class documentation
  return NULL;
}

int AliHLTITSAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					     AliHLTUInt32_t spec,
					     AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  // Handlers for ITS raw data. Even though there are 3 detectors
  // everything is handled in one module library and one HLTOUT handler.
  // This assumes that the data blocks are sent with data type
  // {DDL_RAW :ISDD} and the bit set in the specification corresponding.
  // to detector DDL id.
  // An HLTOUT handler is implemented to extract the equipment id from
  // the specification.
  // Note: Future versions of the framework will provide a default handler
  // class with that functionality.
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITSSDD)) {
      desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
      HLTInfo("module %s handles data block type %s specification %d (0x%x)", 
	      GetModuleId(), AliHLTComponent::DataType2Text(dt).c_str(), spec, spec);
      return 1;
  }
  return 0;
}

AliHLTOUTHandler* AliHLTITSAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITSSDD)) {
    // use the default handler
    if (!fRawDataHandler) {
      fRawDataHandler=new AliHLTOUTSDDRawDataHandler;
    }
    return fRawDataHandler;
  }
  return NULL;
}

int AliHLTITSAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  if (pInstance==fRawDataHandler) {
    delete fRawDataHandler;
    fRawDataHandler=NULL;
    return 0;
  }

  delete pInstance;
  return 0;
}

int AliHLTITSAgent::AliHLTOUTSDDRawDataHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  static int errorCount=0;
  const int maxErrorCount=10;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITSSDD)) {
      int ddlOffset=256;//AliDAQ::DdlIDOffset("ITSSDD");
      int numberOfDDLs=24;//AliDAQ::NumberOfDdls("ITSSDD");
      int ddlNo=0;
      for (;ddlNo<32 && ddlNo<numberOfDDLs; ddlNo++) {
	if (spec&(0x1<<ddlNo)) break;
      }
      if (ddlNo>=32 || ddlNo>=numberOfDDLs) {
	HLTError("invalid specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
	iResult=-ENODEV;
      } else if (spec^(0x1<<ddlNo)) {
	iResult=-EEXIST;
	HLTError("multiple links set in specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
      } else {
	iResult=ddlOffset+ddlNo;
      }
    } else {
      if (errorCount++<10) {
	HLTError("wrong data type: expecting %s, got %s; %s",
		 AliHLTComponent::DataType2Text(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITSSDD).c_str(),
		 AliHLTComponent::DataType2Text(dt).c_str(),
		   errorCount==maxErrorCount?"suppressing further error messages":"");
      }
      iResult=-EFAULT;
    }
  }
  return iResult;
}
