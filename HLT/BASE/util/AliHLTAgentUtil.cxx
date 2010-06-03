// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTAgentUtil.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTUtil library
*/

#include <cassert>
#include "AliHLTAgentUtil.h"
#include "AliHLTOUTHandlerChain.h"

// header files of library components
#include "AliHLTDataGenerator.h"
#include "AliHLTRawReaderPublisherComponent.h"
#include "AliHLTLoaderPublisherComponent.h"
#include "AliHLTRootFileStreamerComponent.h"
#include "AliHLTRootFileWriterComponent.h"
#include "AliHLTRootFilePublisherComponent.h"
#include "AliHLTRootSchemaEvolutionComponent.h"
//#include "AliHLTMCGeneratorComponent.h"
#include "AliHLTESDMCEventPublisherComponent.h"
#include "AliHLTFileWriter.h"
#include "AliHLTFilePublisher.h"
#include "AliHLTBlockFilterComponent.h"
#include "AliHLTMonitoringRelay.h"
#include "AliHLTEsdCollectorComponent.h"
#include "AliHLTReadoutListDumpComponent.h"
#include "AliHLTOUTPublisherComponent.h"
#include "AliHLTCompStatCollector.h"

/** global instance for agent registration */
AliHLTAgentUtil gAliHLTAgentUtil;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAgentUtil)

AliHLTAgentUtil::AliHLTAgentUtil()
  : AliHLTModuleAgent("Util")
  , fCompStatDataHandler(NULL)
  , fStreamerInfoDataHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAgentUtil::~AliHLTAgentUtil()
{
  // see header file for class documentation
}

int AliHLTAgentUtil::CreateConfigurations(AliHLTConfigurationHandler* handler,
					  AliRawReader* /*rawReader*/,
					  AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  if (!handler) return 0;

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // a kChain HLTOUT configuration for processing of {'COMPSTAT':'PRIV'} data blocks
  // produces a TTree object of the component statistics and writes it to disc

  // publisher component
  handler->CreateConfiguration("UTIL-hltout-compstat-publisher", "AliHLTOUTPublisher"   , NULL, "-disable-component-stat");

  // collector configuration
  handler->CreateConfiguration("UTIL-compstat-converter", "StatisticsCollector", "UTIL-hltout-compstat-publisher", "-disable-component-stat -arraysize 10000");

  // writer configuration
  handler->CreateConfiguration("UTIL-compstat-writer", "ROOTFileWriter", "UTIL-compstat-converter", "-datafile HLT.statistics.root -concatenate-events -overwrite");

  return 0;
}

const char* AliHLTAgentUtil::GetReconstructionChains(AliRawReader* /*rawReader*/,
						     AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  return NULL;
}

const char* AliHLTAgentUtil::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return NULL;
}

int AliHLTAgentUtil::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  assert(pHandler);
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTDataGenerator);
  pHandler->AddComponent(new AliHLTRawReaderPublisherComponent);
  pHandler->AddComponent(new AliHLTLoaderPublisherComponent);
  pHandler->AddComponent(new AliHLTRootFileStreamerComponent);
  pHandler->AddComponent(new AliHLTRootFileWriterComponent);
  pHandler->AddComponent(new AliHLTRootFilePublisherComponent);
  pHandler->AddComponent(new AliHLTRootSchemaEvolutionComponent);
  //  pHandler->AddComponent(new AliHLTMCGeneratorComponent);
  pHandler->AddComponent(new AliHLTESDMCEventPublisherComponent);
  pHandler->AddComponent(new AliHLTFileWriter);
  pHandler->AddComponent(new AliHLTFilePublisher);
  pHandler->AddComponent(new AliHLTBlockFilterComponent);
  pHandler->AddComponent(new AliHLTMonitoringRelay);
  pHandler->AddComponent(new AliHLTEsdCollectorComponent);
  pHandler->AddComponent(new AliHLTReadoutListDumpComponent);
  pHandler->AddComponent(new AliHLTOUTPublisherComponent);
  pHandler->AddComponent(new AliHLTCompStatCollector);
  return 0;
}

int AliHLTAgentUtil::GetHandlerDescription(AliHLTComponentDataType dt,
					   AliHLTUInt32_t /*spec*/,
					   AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  // handler for the component statistics data blocks {'COMPSTAT':'PRIV'}
  if (dt==kAliHLTDataTypeComponentStatistics ||
      dt==kAliHLTDataTypeComponentTable) {
      desc=AliHLTOUTHandlerDesc(kChain, dt, GetModuleId());
      return 1;
  }

  // handler for the component statistics data blocks {'ROOTSTRI':'HLT '}
  if (dt==kAliHLTDataTypeStreamerInfo) {
      desc=AliHLTOUTHandlerDesc(kProprietary, dt, GetModuleId());
      return 1;
  }
  return 0;
}

AliHLTOUTHandler* AliHLTAgentUtil::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  // handler for the component statistics data blocks {'COMPSTAT':'PRIV'}
  if (dt==kAliHLTDataTypeComponentStatistics ||
      dt==kAliHLTDataTypeComponentTable) {
    if (fCompStatDataHandler==NULL)
      fCompStatDataHandler=new AliHLTOUTHandlerChain("chains=UTIL-compstat-writer");
    return fCompStatDataHandler;
  }

  // handler for the component statistics data blocks {'ROOTSTRI':'HLT '}
  if (dt==kAliHLTDataTypeStreamerInfo) {
    if (fStreamerInfoDataHandler==NULL)
      fStreamerInfoDataHandler=new AliHLTStreamerInfoHandler;
    return fStreamerInfoDataHandler;
  }

  return NULL;
}

int AliHLTAgentUtil::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  if (pInstance==fCompStatDataHandler) {
    delete fCompStatDataHandler;
    fCompStatDataHandler=NULL;
  }

  if (pInstance==fStreamerInfoDataHandler) {
    delete fStreamerInfoDataHandler;
    fStreamerInfoDataHandler=NULL;
  }
  return 0;
}
