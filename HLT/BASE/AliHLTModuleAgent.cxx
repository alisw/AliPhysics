// @(#) $Id$

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
  :
  fpNext(NULL),
  fpComponentHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  Register(this);
}

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
   AliHLTModuleAgent* pCurrent=fAnchor;
   while (pCurrent!=NULL && strcmp(pCurrent->GetName(), agent)!=0) pCurrent=pCurrent->fpNext;
   if (pCurrent) {
     log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		 "agent %s available", pCurrent->GetName());
   } else {
     log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		 "agent %s not found", agent);
   }
  } else {
   AliHLTModuleAgent* pCurrent=fAnchor;
   log.Logging(kHLTLogInfo, "AliHLT", "", "-----------------------");
   log.Logging(kHLTLogInfo, "AliHLT", "", "available module agents");
   if (pCurrent==NULL)
     log.Logging(kHLTLogInfo, "AliHLT", "", "   none");
   while (pCurrent) {
     TString msg;
     msg.Form("   %s : %p", pCurrent->GetName(), pCurrent);
     log.Logging(kHLTLogInfo, "AliHLT", "", msg.Data());
     pCurrent=pCurrent->fpNext;
   }
   log.Logging(kHLTLogInfo, "AliHLT", "", "-----------------------");
  }
}

int AliHLTModuleAgent::CreateConfigurations(AliHLTConfigurationHandler* /*handler*/,
					    AliRawReader* /*rawReader*/,
					    AliRunLoader* /*runloader*/) const
{
  // default method, nothing to be done, child classes can overload
  return 0;
}

const char* AliHLTModuleAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
						       AliRunLoader* /*runloader*/) const
{
  // default method, nothing to be done, child classes can overload
  return NULL;
}

const char* AliHLTModuleAgent::GetRequiredComponentLibraries() const
{
  // default method, nothing to be done, child classes can overload
  return NULL;
}

int AliHLTModuleAgent::ActivateComponentHandler(AliHLTComponentHandler* pHandler)
{
  int iResult=0;
  if (pHandler==NULL) {
    if (fpComponentHandler!=NULL) {
      // reset and think about deregistration
      fpComponentHandler=NULL;
      HLTWarning("deregistration of components not yet implemented");
    }
    return 0;
  }
  if (fpComponentHandler!=NULL) {
    if (pHandler!=fpComponentHandler) {
      HLTError("only one component handler can be activated per agent");
      return -EINVAL;
    }
    return 0;
  }
  if ((iResult=RegisterComponents(pHandler))>=0) {
    fpComponentHandler=pHandler;
  }
  return iResult;
}

int AliHLTModuleAgent::RegisterComponents(AliHLTComponentHandler* /*pHandler*/) const
{
  // default method, nothing to be done, child classes can overload
  return 0;
}

AliHLTModuleAgent* AliHLTModuleAgent::fAnchor=NULL;
AliHLTModuleAgent* AliHLTModuleAgent::fCurrent=NULL;
int AliHLTModuleAgent::fCount=0;

AliHLTModuleAgent* AliHLTModuleAgent::GetFirstAgent()
{
  // see header file for function documentation
  fCurrent=fAnchor;
  return fAnchor;
}

AliHLTModuleAgent* AliHLTModuleAgent::GetNextAgent()
{
  // see header file for function documentation
  if (fCurrent!=NULL) fCurrent=fCurrent->fpNext;
  return fCurrent;
}

int AliHLTModuleAgent::Register(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  if (fAnchor==NULL) {
    fAnchor=pAgent;
  } else {
    pAgent->fpNext=fAnchor;
    fAnchor=pAgent;
  }
  //  log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Register", "", "module agent %p registered", pAgent);
  fCount++;
  return 0;	
}

int AliHLTModuleAgent::Unregister(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  fCurrent=NULL;
  AliHLTModuleAgent* prev=NULL;
  AliHLTModuleAgent* handler=fAnchor;
  while (handler!=NULL && handler!=pAgent) {
    prev=handler;
    handler=handler->fpNext;
  }
  if (handler) {
    if (prev==NULL) {
      fAnchor=handler->fpNext;
    } else {
      prev->fpNext=handler->fpNext;
    }
    //log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Unregister", "", "module agent %p removed", pAgent);
    fCount--;
  }
  return 0;
}
