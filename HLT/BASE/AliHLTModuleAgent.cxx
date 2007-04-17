// @(#) $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTModuleAgent.cxx
    @author Matthias Richter
    @date   
    @brief  Agent helper class for component libraries.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTModuleAgent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTModuleAgent)

AliHLTModuleAgent::AliHLTModuleAgent()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  Register(this);
}

AliHLTModuleAgent::AliHLTModuleAgent(const AliHLTModuleAgent&)
  :
  TObject(),
  AliHLTLogging()
{
  // see header file for function documentation
}

AliHLTModuleAgent& AliHLTModuleAgent::operator=(const AliHLTModuleAgent&)
{
  // see header file for function documentation
  return *this;
}

TList AliHLTModuleAgent::fgAgentList;
TObjLink* AliHLTModuleAgent::fgCurrentLnk=NULL;

AliHLTModuleAgent::~AliHLTModuleAgent()
{
  // see header file for function documentation
  Unregister(this);
}

void AliHLTModuleAgent::PrintStatus(const char* agent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (agent) {
    TObject* pAgent=fgAgentList.FindObject(agent);
    if (pAgent) {
      log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		  "agent %s available", pAgent->GetName());
    } else {
      log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		  "agent %s not found", agent);
    }
  } else {
    TObjLink* lnk=fgAgentList.FirstLink();
    log.Logging(kHLTLogInfo, "AliHLT", "", "-----------------------");
    log.Logging(kHLTLogInfo, "AliHLT", "", "available module agents");
    if (lnk==NULL) 
      log.Logging(kHLTLogInfo, "AliHLT", "", "   none");
    while (lnk) {
      TString msg;
      msg.Form("   %s : %p", ((AliHLTModuleAgent*)lnk->GetObject())->GetName(), lnk->GetObject());
      log.Logging(kHLTLogInfo, "AliHLT", "", msg.Data());
      lnk=lnk->Next();
    }
    log.Logging(kHLTLogInfo, "AliHLT", "", "-----------------------");
  }
}

int AliHLTModuleAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
					    AliRunLoader* runloader) const
{
  // default method, nothing to be done, child classes can overload
  if (handler==NULL && runloader==NULL) {
    // get rid of 'unused parameter' warning
  }
  return 0;
}

const char* AliHLTModuleAgent::GetLocalRecConfigurations(AliRunLoader* runloader) const
{
  // default method, nothing to be done, child classes can overload
  if (runloader==NULL) {
    // get rid of 'unused parameter' warning
  }
  return NULL;
}

const char* AliHLTModuleAgent::GetEventRecConfigurations(AliRunLoader* runloader) const
{
  // default method, nothing to be done, child classes can overload
  if (runloader==NULL) {
    // get rid of 'unused parameter' warning
  }
  return NULL;
}

const char* AliHLTModuleAgent::GetRequiredComponentLibraries() const
{
  // default method, nothing to be done, child classes can overload
  return NULL;
}

int AliHLTModuleAgent::RegisterComponents(AliRunLoader* runloader) const
{
  if (runloader==NULL) {
    // get rid of 'unused parameter' warning
  }
  // default method, nothing to be done, child classes can overload
  return 0;
}

AliHLTModuleAgent* AliHLTModuleAgent::GetFirstAgent()
{
  // see header file for function documentation
  fgCurrentLnk=fgAgentList.FirstLink();
  if (fgCurrentLnk==NULL) return NULL;
  return (AliHLTModuleAgent*)fgCurrentLnk->GetObject();
}

AliHLTModuleAgent* AliHLTModuleAgent::GetNextAgent()
{
  // see header file for function documentation
  if (fgCurrentLnk==NULL) return NULL;
  fgCurrentLnk=fgCurrentLnk->Next();
  if (fgCurrentLnk==NULL) return NULL;
  return (AliHLTModuleAgent*)fgCurrentLnk->GetObject();
}

int AliHLTModuleAgent::Register(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  if (fgAgentList.FindObject(pAgent)==NULL) {
    log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Register", "", "module agent %p registered", pAgent);
    fgAgentList.Add(pAgent);
  }
  return 0;
}

int AliHLTModuleAgent::Unregister(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  if (fgAgentList.FindObject(pAgent)!=NULL) {
    log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Unregister", "", "module agent %p removed", pAgent);
    fgAgentList.Remove(pAgent);
  } else {
  }
  return 0;
}
