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

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModuleAgent.h"
#include "AliHLTOUTHandler.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTModuleAgent)

AliHLTModuleAgent::AliHLTModuleAgent(const char* id)
  :
  fpNext(NULL),
  fpComponentHandler(NULL),
  fModuleId(id)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  Register(this);
}

const AliHLTModuleAgent::AliHLTOUTHandlerDesc AliHLTModuleAgent::fgkVoidHandlerDesc;

AliHLTModuleAgent::~AliHLTModuleAgent()
{
  // see header file for function documentation
  Unregister(this);
}

const char* AliHLTModuleAgent::GetModuleId() const
{
  // see header file for function documentation
  return fModuleId.Data();
}

void AliHLTModuleAgent::PrintStatus(const char* agent)
{
  // see header file for function documentation
  AliHLTLogging log;
 if (agent) {
   AliHLTModuleAgent* pCurrent=fgAnchor;
   while (pCurrent!=NULL && strcmp(pCurrent->GetName(), agent)!=0) pCurrent=pCurrent->fpNext;
   if (pCurrent) {
     log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		 "agent %s available", pCurrent->GetName());
   } else {
     log.Logging(kHLTLogInfo, "AliHLTModuleAgent::PrintStatus", "module agents", 
		 "agent %s not found", agent);
   }
  } else {
   AliHLTModuleAgent* pCurrent=fgAnchor;
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

int AliHLTModuleAgent::GetHandlerDescription(AliHLTComponentDataType /*dt*/,
					     AliHLTUInt32_t /*spec*/,
					     AliHLTOUTHandlerDesc& /*desc*/) const
{
  // default method, nothing to be done, child classes can overload
  return 0;
}

AliHLTOUTHandler* AliHLTModuleAgent::GetOutputHandler(AliHLTComponentDataType /*dt*/,
						      AliHLTUInt32_t /*spec*/)
{
  // default method, nothing to be done, child classes can overload
  return NULL;
}

int AliHLTModuleAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
  // default method, simply deletes object
  if (pInstance) return -EINVAL;
  delete pInstance;
  return 0;
}


// likely to be moved to AliHLTOUTHandler
// AliRawStream* AliHLTModuleAgent::GetRawStream(AliHLTComponentDataType /*dt*/,
// 					      AliHLTUInt32_t /*spec*/,
// 					      const AliHLTOUT* /*pData*/) const
// {
//   // default method, nothing to be done, child classes can overload
//   return NULL;
// }

int AliHLTModuleAgent::ActivateComponentHandler(AliHLTComponentHandler* pHandler)
{
  // see header file for function documentation
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

AliHLTModulePreprocessor* AliHLTModuleAgent::GetPreprocessor()
{
  // default method, nothing to be done, child classes can overload
  return NULL;
}

AliHLTModuleAgent* AliHLTModuleAgent::fgAnchor=NULL;
AliHLTModuleAgent* AliHLTModuleAgent::fgCurrent=NULL;
int AliHLTModuleAgent::fgCount=0;

AliHLTModuleAgent* AliHLTModuleAgent::GetFirstAgent()
{
  // see header file for function documentation
  fgCurrent=fgAnchor;
  return fgAnchor;
}

AliHLTModuleAgent* AliHLTModuleAgent::GetNextAgent()
{
  // see header file for function documentation
  if (fgCurrent!=NULL) fgCurrent=fgCurrent->fpNext;
  return fgCurrent;
}

string AliHLTModuleAgent::GetAgentIds()
{
  // see header file for function documentation
  string ids;
  ids.clear();
  for (AliHLTModuleAgent* pCurrent=fgAnchor;
       pCurrent;
       pCurrent=pCurrent->fpNext) {
    if (ids.size()>0) ids+=" ";
    ids+=pCurrent->GetModuleId();
  }

  return ids;
}

int AliHLTModuleAgent::Register(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  // The following check is for extra protection. In some strange cases the agent might
  // try to register itself more than once. So we need to check for that and prevent it.
  // Otherwise we create a cycle in our linked list and go into an infinite loop.
  AliHLTModuleAgent* current=fgAnchor;
  while (current!=NULL) {
    if (current == pAgent) return 0;
    current = current->fpNext;
  }
  if (fgAnchor==NULL) {
    fgAnchor=pAgent;
  } else {
    pAgent->fpNext=fgAnchor;
    fgAnchor=pAgent;
  }
  //  log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Register", "", "module agent %p registered", pAgent);
  fgCount++;
  return 0;	
}

int AliHLTModuleAgent::Unregister(AliHLTModuleAgent* pAgent)
{
  // see header file for function documentation
  AliHLTLogging log;
  if (!pAgent) return -EINVAL;
  fgCurrent=NULL;
  AliHLTModuleAgent* prev=NULL;
  AliHLTModuleAgent* handler=fgAnchor;
  while (handler!=NULL && handler!=pAgent) {
    prev=handler;
    handler=handler->fpNext;
  }
  if (handler) {
    if (prev==NULL) {
      fgAnchor=handler->fpNext;
    } else {
      prev->fpNext=handler->fpNext;
    }
    //log.Logging(kHLTLogDebug, "AliHLTModuleAgent::Unregister", "", "module agent %p removed", pAgent);
    fgCount--;
  }
  return 0;
}
