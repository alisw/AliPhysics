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

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandlerChain)

AliHLTOUTHandlerChain::AliHLTOUTHandlerChain(const char* arguments)
  :
  fChains(),
  fOptions(),
  fpSystem(NULL),
  fpTask(NULL)
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
      AliHLTOUT* pSubCollection=dynamic_cast<AliHLTOUT*>(fpTask);
      pSubCollection->Init();
      pData->AddSubCollection(pSubCollection);
    } else {
      fpTask->Reset();
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
    fpSystem = new AliHLTSystem(GetGlobalLoggingLevel());
    if (fpSystem) {
      if ((iResult=fpSystem->ScanOptions(fOptions.Data()))>=0) {
	// load configurations if not specified by external macro
	if (!fOptions.Contains("config="))
	  iResult=CreateConfigurations(fpSystem->fpConfigurationHandler);

	if (iResult>=0) {
	  iResult=fpSystem->BuildTaskList(fChains.Data());
	}

	// add AliHLTOUTTask on top of the configuartions in order to
	// collect the data
	fpTask=new AliHLTOUTTask(fChains.Data());
      }
    } else {
      iResult=-ENOMEM;
    }
    if (iResult>=0 && fpSystem && fpTask) {
      if (fpTask->GetConf() && fpTask->GetConf()->SourcesResolved()>=0) {
	iResult=fpSystem->InsertTask(fpTask);
      } else {
	HLTError("HLTOUT task (%s) sources not resolved", fpTask->GetName());
	iResult=-ENOENT;
      }
    }
    if (iResult<0 || !fpTask) {
      SetStatusFlag(kHandlerError);
      if (fpSystem) delete fpSystem; fpSystem=NULL;
      if (fpTask) delete fpTask; fpTask=NULL;
    }
  }
  return iResult;
}
