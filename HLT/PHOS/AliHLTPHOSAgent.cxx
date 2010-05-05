
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Oystein Djuvsland                                     *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTPHOSAgent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  Agent of the libAliHLTPHOS library
*/

#include "AliHLTPHOSAgent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTOUT.h"
#include "AliHLTOUTHandlerChain.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"


// #include "AliHLTPHOSConstant.h"
#include "AliHLTPHOSConstants.h"


#include "AliHLTPHOSMapper.h"

/** global instance for agent registration */
AliHLTPHOSAgent gAliHLTPHOSAgent;

// component headers
#include "AliHLTPHOSCalibrationComponent.h"
#include "AliHLTPHOSClusterAnalyserComponent.h"
#include "AliHLTPHOSClusterizerComponent.h"
#include "AliHLTPHOSDigitMakerComponent.h"
#include "AliHLTPHOSESDEntriesMakerComponent.h"
#include "AliHLTPHOSHistogramProducerComponent.h"
#include "AliHLTPHOSModuleCalibrationProcessorComponent.h"
#include "AliHLTPHOSMonitorTriggerComponent.h"
#include "AliHLTPHOSRawAnalyzerComponent.h"
#include "AliHLTPHOSRawAnalyzerCrudeComponentv2.h"
#include "AliHLTPHOSRawAnalyzerPeakFinderComponent.h"
#include "AliHLTPHOSRcuCalibrationProcessorComponent.h"
#include "AliHLTPHOSRcuDAComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTPHOSAgent)

AliHLTPHOSAgent::AliHLTPHOSAgent()
  :
  AliHLTModuleAgent("PHOS"),
  fRawDataHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTPHOSAgent::~AliHLTPHOSAgent()
{
  // see header file for class documentation
}

int AliHLTPHOSAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					  AliRawReader* /*rawReader*/,
					 AliRunLoader* /*runloader*/) const
{
  // see header file for class documentation
  if (handler) 
    {
      //      const char* cdbEntry="PHOS/Calib/Parameters";
      //      AliCDBManager* pMan=AliCDBManager::Instance();
//       AliPHOSParam* pPHOSParam=NULL;
//       if (pMan) 
// 	{
// 	  AliCDBEntry *pEntry = pMan->Get(cdbEntry);
// 	  if (pEntry && 
// 	      pEntry->GetObject() &&
// 	      (pPHOSParam=dynamic_cast<AliPHOSParam*>(pEntry->GetObject()))) 
// 	    {
// 	    } else 
// 	    {
// 	      HLTWarning("can not load AliPHOSParam from CDB entry %s", cdbEntry);
// 	    }
// 	}


      int moduleStart = 0;
      //   int moduleEnd = PhosHLTConst::NMODULES - 1;
      int moduleEnd = 5 - 1;

      int rcuStart = 0;
      //   int rcuEnd = PhosHLTConst::NRCUSPERMODULE - 1;
      int rcuEnd = PhosHLTConst::4 - 1;

      TString mergerInput;
      TString sinkClusterInput;
      TString emInput;
      for (int module = moduleStart; module < moduleEnd; module++) 
	{
	  TString clInput;

	  for(int rcu = rcuStart; rcu < rcuEnd; rcu++) 
	    {
	      TString arg, publisher, ra, dm;
	      // raw data publisher components
	      publisher.Form("PHOS-RP_%02d_%d", module, rcu);
	
	      //     arg.Form("-datatype 'DDL_RAW ' 'PHOS'  -dataspec 0x ", 0x1 << (module*PhosHLTConst::NRCUSPERMODULE + rcu));
	      arg.Form("-datatype 'DDL_RAW ' 'PHOS'  -dataspec 0x ", 0x1 << (module*4 + rcu));

	      handler->CreateConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

	      // Raw analyzer
	      arg = "";
	      ra.Form("PHOS-RA_%02d_%d", module, rcu);
	      handler->CreateConfiguration(ra.Data(), "PhosRawCrudev2", publisher.Data(), arg.Data());

	      // digit maker components
	      dm.Form("PHOS-DM_%02d_%d", module, rcu);
	      arg="";
	      arg.Form("-sethighgainfactor 0.005 -setlowgainfactor 0.08 -setdigitthresholds 0.005 0.002");
	      handler->CreateConfiguration(dm.Data(), "PhosDigitMaker", ra.Data(), arg.Data());

	      if(clInput.Length() > 0) clInput += " ";
	      clInput+=dm;
	    }
	  TString arg, cl, ca;

	  cl.Form("PHOS-CF_%02d", module);
	  arg = "";
	  arg.Form("-digitthreshold 0.005 -recpointthreshold 0.1 -modulemode");
	  handler->CreateConfiguration(cl.Data(), "PhosClusterizer", clInput.Data(), arg.Data());
	
	  ca.Form("PHOS-CA_%02d", module);
	  arg = " ";
	  handler->CreateConfiguration(ca.Data(), "PhosClusterAnalyser", cl.Data(), arg.Data());

	  if(emInput.Length() > 0) emInput += " ";
	  emInput += ca;
	}
      
      TString arg, em;
  
      // tracker finder components
      em.Form("PHOS-EM");
      arg = " ";
      handler->CreateConfiguration(em.Data(), "PhosEsdEntriesMaker", emInput.Data(), " ");
    }
  return 0;
}

const char* AliHLTPHOSAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						    AliRunLoader* runloader) const
{
  // see header file for class documentation
  if (runloader) {
    // reconstruction chains for AliRoot simulation
    // Note: run loader is only available while running embedded into
    // AliRoot simulation
    if (runloader->GetLoader("PHOSLoader") != NULL)
      return "PHOS-EM";
  }
  return NULL;
}

const char* AliHLTPHOSAgent::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return NULL;
}

int AliHLTPHOSAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
  // see header file for class documentation
  if (!pHandler) return -EINVAL;
  pHandler->AddComponent(new AliHLTPHOSRawAnalyzerCrudeComponentv2);
  pHandler->AddComponent(new AliHLTPHOSDigitMakerComponent);
  pHandler->AddComponent(new AliHLTPHOSClusterizerComponent);
  pHandler->AddComponent(new AliHLTPHOSClusterAnalyserComponent);			 
  pHandler->AddComponent(new AliHLTPHOSESDEntriesMakerComponent);

  return 0;
}

int AliHLTPHOSAgent::GetHandlerDescription(AliHLTComponentDataType dt,
					  AliHLTUInt32_t spec,
					  AliHLTOUTHandlerDesc& desc) const
{
  // see header file for class documentation

  AliHLTPHOSMapper mapper;
  mapper.InitDDLSpecificationMapping();
  
  // raw data blocks to be fed into offline reconstruction
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginPHOS)) 
    {
      if(mapper.GetDDLFromSpec(spec) >= 0)
	{
	  desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
	  return 1;
	} 
      else 
	{
	  HLTWarning("Handler can not process data inconsistent with a single PHOS DDL from specification % d", spec);
	  return 0;
	}
    }
  return 0;
}

AliHLTOUTHandler* AliHLTPHOSAgent::GetOutputHandler(AliHLTComponentDataType dt,
						   AliHLTUInt32_t /*spec*/)
{
  // see header file for class documentation

  // raw data blocks to be fed into offline reconstruction
  if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginPHOS)) 
    {
      if (!fRawDataHandler) 
	{
	  fRawDataHandler = new AliHLTPHOSAgent::AliHLTPHOSRawDataHandler;
	}
      return fRawDataHandler;
    }

  return NULL;
}

int AliHLTPHOSAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // see header file for class documentation
  if (pInstance==NULL) return -EINVAL;

  if (pInstance==fRawDataHandler) {
    delete fRawDataHandler;
    fRawDataHandler=NULL;
  }
  return 0;
}

AliHLTPHOSAgent::AliHLTPHOSRawDataHandler::AliHLTPHOSRawDataHandler()
{
  // see header file for class documentation
}

AliHLTPHOSAgent::AliHLTPHOSRawDataHandler::~AliHLTPHOSRawDataHandler()
{
  // see header file for class documentation
}

int AliHLTPHOSAgent::AliHLTPHOSRawDataHandler::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;

  AliHLTPHOSMapper mapper;
  mapper.InitDDLSpecificationMapping();

  AliHLTComponentDataType dt = kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult = pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) 
    {
      int ddl = -1;
      if((ddl = mapper.GetDDLFromSpec(spec)) >=0)
	{
	  iResult = ddl;
	}
    } 
  else 
    {
      HLTError("Handler can not process data inconsistent with a single PHOS DDL from specification % d", spec);
      iResult=-EBADMSG;
    }
  return iResult;
}
