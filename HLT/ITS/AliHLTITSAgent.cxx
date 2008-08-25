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

// header files of library components

// header file of the module preprocessor

// raw data handler of HLTOUT data
#include "AliHLTOUTHandlerEquId.h"

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
  //pHandler->AddComponent(new ...);

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
  // {DDL_RAW :ITS } and the equipment id as specification
  // The default behavior of AliHLTOUTHandlerEquId is used.
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITS)) {
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
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginITS)) {
    // use the default handler
    if (!fRawDataHandler) {
      fRawDataHandler=new AliHLTOUTHandlerEquId;
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
