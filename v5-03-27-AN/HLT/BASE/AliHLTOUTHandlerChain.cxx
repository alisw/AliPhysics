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

/** @file   AliHLTOUTHandlerChain.cxx
    @author Matthias Richter
    @date   24.06.2008
    @brief  HLTOUT handler of type kChain.
*/

#include "AliHLTOUTHandlerChain.h"
#include "AliHLTOUT.h"
#include "AliHLTSystem.h"
#include "AliHLTOUTTask.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandlerChain)

AliHLTOUTHandlerChain::AliHLTOUTHandlerChain(const char* arguments)
  :
  fChains(),
  fOptions(),
  fpSystem(NULL),
  fbHaveOutput(false)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (arguments) {
    TString args=arguments;
    TObjArray* pTokens=args.Tokenize(" ");
    if (pTokens) {
      int iEntries=pTokens->GetEntries();
      for (int i=0; i<iEntries; i++) {
	TString token=(((TObjString*)pTokens->At(i))->GetString());
	if (token.Contains("chains=")) {
	  TString param=token.ReplaceAll("chains=", "");
	  fChains=param.ReplaceAll(",", " ");
	} else {
	  if (!fOptions.IsNull()) fOptions+=" ";
	  fOptions+=token;
	}
      }
      delete pTokens;
    }
  } 
}

AliHLTOUTHandlerChain::~AliHLTOUTHandlerChain()
{
  // see header file for class documentation
  if (fpSystem) {
    // TODO: the EOR is currenttly not send because the reconstruction
    // chian is not stopped. Trying it here gives an error, there is
    // some state mismatch in AliHLTSystem. Probably due to the fact,
    // that the AliHLTSystem::Reconstruct method is not used
//     if (!fpSystem->CheckStatus(AliHLTSystem::kError)) {
//       // send specific 'event' to execute the stop sequence
//       fpSystem->Reconstruct(0, NULL, NULL);
//     }
    delete fpSystem;
  }
}

int AliHLTOUTHandlerChain::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  int iResult=0;

  if (CheckStatus(kHandlerError)) {
    HLTWarning("kChain handler '%s' in error state, skipping processing of associated HLTOUT blocks", fChains.Data());
    return -EPERM;
  }

  if (!fpSystem && (iResult=InitSystem())<0) {
    return iResult;
  }

  if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
    HLTWarning("kChain handler '%s': system in error state, skipping processing of associated HLTOUT blocks", fChains.Data());
    return -EACCES;
  }

  // run one event and do not stop the chain
  {
    AliHLTOUT::AliHLTOUTGlobalInstanceGuard g(pData);
    if ((iResult=fpSystem->Run(1,0))>=0) {
      // sub-collection is going to be reset from the
      // parent HLTOUT collection
      AliHLTOUTTask* pTask=fpSystem->GetHLTOUTTask();

      // either have the task or none of the chains controlled by the chain
      // handler has output
      assert(pTask || !fbHaveOutput);
      if (pTask) {
      AliHLTOUT* pSubCollection=dynamic_cast<AliHLTOUT*>(pTask);
      pSubCollection->Init();

      // filter out some data blocks which should not be processed
      // in the next stage:
      // 1. we are not interested in the component statistics
      //    produced in the HLTOUT handler chain
      for (iResult=pSubCollection->SelectFirstDataBlock();
	   iResult>=0;
	   iResult=pSubCollection->SelectNextDataBlock()) {
	AliHLTComponentDataType dt=kAliHLTVoidDataType;
	AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
	pSubCollection->GetDataBlockDescription(dt, spec);
	if (dt==kAliHLTDataTypeComponentStatistics) {
	  pSubCollection->MarkDataBlockProcessed();
	}
      }
      pData->AddSubCollection(pSubCollection);
      } else if (fbHaveOutput) {
	// this is an error condition since task has been created and should
	// be available
	HLTError("can not get instance of HLTOUT task from HLT system %p", fpSystem);
      }
    }
  }

  return iResult;
}

int AliHLTOUTHandlerChain::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/)
{
  //default implementation, nothing to do
  return 0;
}

int AliHLTOUTHandlerChain::InitSystem()
{
  // see header file for class documentation
  int iResult=0;
  if (!fpSystem) {
    // init AliHLTSystem
    TString systemName="kChain_"; systemName+=fChains;
    systemName.ReplaceAll(" ", "_");
    fpSystem = new AliHLTSystem(GetGlobalLoggingLevel(), systemName);
    if (fpSystem) {
      if ((iResult=fpSystem->ScanOptions(fOptions.Data()))>=0) {
	// load configurations if not specified by external macro
	if (!fOptions.Contains("config="))
	  iResult=CreateConfigurations(fpSystem->GetConfigurationHandler());

	if (iResult>=0) {
	  iResult=fpSystem->BuildTaskList(fChains.Data());
	}

	// add AliHLTOUTTask on top of the configuartions in order to
	// collect the data
	// remember if task has been created (result>0)
	fbHaveOutput=((iResult=fpSystem->AddHLTOUTTask(fChains.Data()))>0);
      }
    } else {
      iResult=-ENOMEM;
    }
    if (iResult<0) {
      SetStatusFlag(kHandlerError);
      if (fpSystem) delete fpSystem; fpSystem=NULL;
    }
  }
  return iResult;
}
