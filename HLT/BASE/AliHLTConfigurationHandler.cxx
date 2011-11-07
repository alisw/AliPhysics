// $Id$
// splitted from AliHLTConfiguration.cxx,v 1.25 2007/10/12 13:24:47
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

/** @file   AliHLTConfigurationHandler.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of HLT tasks.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include <iostream>
#include <string>
#include "AliHLTConfigurationHandler.h"
#include "AliHLTConfiguration.h"
#include "AliHLTErrorGuard.h"
#include "TMap.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTConfigurationHandler)

AliHLTConfigurationHandler::AliHLTConfigurationHandler()
  : AliHLTLogging()
  , fgListConfigurations()
  , fgListScheduledRegistrations()
  , fFlags(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  SetLocalLoggingLevel(kHLTLogInfo);
}

AliHLTConfigurationHandler::~AliHLTConfigurationHandler()
{
  // see header file for function documentation
  TObjLink* lnk=NULL;
  while ((lnk=fgListConfigurations.FirstLink())!=NULL) {
    AliHLTConfiguration* pConf=(AliHLTConfiguration*)lnk->GetObject();
    HLTDebug("delete configuration \"%s\"", pConf->GetName());
    fgListConfigurations.Remove(lnk);
    delete pConf;
  }
  fgListScheduledRegistrations.Delete();
}

AliHLTConfigurationHandler* AliHLTConfigurationHandler::fgpInstance=NULL;
int AliHLTConfigurationHandler::fgNofInstances=0;
TMap* AliHLTConfigurationHandler::fgpSubstitutions=NULL;

AliHLTConfigurationHandler* AliHLTConfigurationHandler::CreateHandler()
{
  // see header file for class documentation
  if (!fgpInstance) fgpInstance=new AliHLTConfigurationHandler;
  fgNofInstances++;
  return fgpInstance;
}

int AliHLTConfigurationHandler::Destroy()
{
  // see header file for class documentation
  int nofInstances=0;
  if (fgpInstance==this) {
    nofInstances=--fgNofInstances;
  }
  if (fgNofInstances==0) {
    fgpInstance = NULL;
    if (fgpSubstitutions) delete fgpSubstitutions;
    fgpSubstitutions=NULL;
  }
  if (nofInstances==0) delete this;
  return nofInstances;
}


int AliHLTConfigurationHandler::RegisterConfiguration(AliHLTConfiguration* pConf)
{
  // see header file for function documentation
  int iResult=0;
  if (pConf) {
    AliHLTConfiguration* pClone=new AliHLTConfiguration(*pConf);
    if (IsActive()) {      
      AliHLTConfiguration* pExisting=NULL;
      if ((pExisting=FindConfiguration(pConf->GetName())) == NULL) {
        fgListConfigurations.Add(pClone);
        HLTDebug("configuration \"%s\" (%p) registered from %p", pClone->GetName(), pClone, pConf);

        // mark all configurations with unresolved dependencies for re-evaluation
        TObjLink* lnk=fgListConfigurations.FirstLink();
        while (lnk) {
          AliHLTConfiguration* pSrc=(AliHLTConfiguration*)lnk->GetObject();
          if (pSrc && pSrc!=pClone && pSrc->SourcesResolved()!=1) {
            pSrc->InvalidateSources();
          }
          lnk=lnk->Next();
        }
      } else {
        if ((*pExisting)!=(*pConf)) {
          iResult=-EEXIST;
          HLTWarning("configuration \"%s\" already registered with different properties", pConf->GetName());
        }
      }
    } else if (IsScheduling()) {
      fgListScheduledRegistrations.Add(pClone);
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTConfigurationHandler::CreateConfiguration(const char* id, const char* component, const char* sources, const char* arguments)
{
  // see header file for function documentation
  int iResult=0;
  // if this handler is the global instance the configuration is added
  // automatically in the creation of the AliHLTConfiguration object
  // the global instance must be deactivated otherwise in order to just create
  // the object and then add it to THIS handler
  bool bIamGlobal=fgpInstance==this;
  if (!bIamGlobal && fgpInstance) {
    // deactivate the automatic registration in the global handler
    fgpInstance->Deactivate(false);
  }
  AliHLTConfiguration* pConf= new AliHLTConfiguration(id, component, sources, arguments);
  if (pConf) {
    if (bIamGlobal) {
    // the configuration will be registered automatically, if this failes the configuration
    // is missing -> delete it
    if (FindConfiguration(id)==NULL) {
      delete pConf;
      pConf=NULL;
      iResult=-EEXIST;
    }
    } else {
      RegisterConfiguration(pConf);
    }
  } else {
    HLTError("system error: object allocation failed");
    iResult=-ENOMEM;
  }
  if (!bIamGlobal && fgpInstance) {
    // deactivate the automatic registration in the global handler
    fgpInstance->Activate();
  }
  return iResult;
}

void AliHLTConfigurationHandler::PrintConfigurations()
{
  // see header file for function documentation
  HLTLogKeyword("configuration listing");
  HLTMessage("registered configurations:");
  TObjLink *lnk = fgListConfigurations.FirstLink();
  while (lnk) {
    TObject *obj = lnk->GetObject();
    HLTMessage("  %s", obj->GetName());
    lnk = lnk->Next();
  }
}

void AliHLTConfigurationHandler::Print(const char* option)
{
  // print info
  TString argument(option);
  if (argument.BeginsWith("treeroot=")) {
    argument.ReplaceAll("treeroot=", "");
    if (argument.IsNull()) {
      cout << "invalid argument to option 'treeroot=', please specify configuration" << endl;
      return;
    }
    // TODO: add functionality to print a dependency tree beginning from a root configuration
    // add also option to limit the depth
    cout << "need to implement option 'treeview', argument " << argument << endl;
    return;
  }

  // default: print all
  PrintConfigurations();
}

int AliHLTConfigurationHandler::RemoveConfiguration(const char* id)
{
  // see header file for function documentation
  int iResult=0;
  if (id) {
    AliHLTConfiguration* pConf=NULL;
    if ((pConf=FindConfiguration(id))!=NULL) {
      iResult=RemoveConfiguration(pConf);
      delete pConf;
      pConf=NULL;
    } else {
      HLTWarning("can not find configuration \"%s\"", id);
      iResult=-ENOENT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTConfigurationHandler::RemoveConfiguration(AliHLTConfiguration* pConf)
{
  // see header file for function documentation
  int iResult=0;
  if (pConf) {
    // remove the configuration from the list
    HLTDebug("remove configuration \"%s\"", pConf->GetName());
    fgListConfigurations.Remove(pConf);
    // remove cross links in the remaining configurations
    TObjLink* lnk=fgListConfigurations.FirstLink();
    while (lnk && iResult>=0) {
      AliHLTConfiguration* pRem=(AliHLTConfiguration*)lnk->GetObject();
      if (pRem) {
	pRem->InvalidateSource(pConf);
      } else {
	iResult=-EFAULT;
      }
      lnk=lnk->Next();
    }
  }
  return iResult;
}

AliHLTConfiguration* AliHLTConfigurationHandler::FindConfiguration(const char* id)
{
  // see header file for function documentation
  AliHLTConfiguration* pConf=NULL;
  if (id) {
    pConf=(AliHLTConfiguration*)fgListConfigurations.FindObject(id); 
  }
  return pConf;
}

int AliHLTConfigurationHandler::Deactivate(bool schedule) {
  // see header file for function documentation
  fFlags|=kInactive;
  if (schedule)
    fFlags|=kScheduling;
  return 0;
}

int AliHLTConfigurationHandler::Activate() {
  // see header file for function documentation
  fFlags&=~kInactive;
  if (IsScheduling()) {
    fFlags&=~kScheduling;
    TObjLink *lnk = fgListScheduledRegistrations.FirstLink();
    while (lnk) {
      RegisterConfiguration((AliHLTConfiguration*)lnk->GetObject());
      lnk = lnk->Next();
    }
    ClearScheduledRegistrations();
  }
  return 0;
}

int AliHLTConfigurationHandler::MissedRegistration(const char* name)
{
  /// indicate a failed attempt to register because of unavailable global instance

  /// everything fine if global instance is inactive
  if (fgpInstance) {
    if (fgpInstance->IsActive()) {
      static AliHLTErrorGuard g("AliHLTConfigurationHandler", "MissedRegistration",
				"internal error, global instance available but registration of configuration failed");
      (++g).Throw(1);
    }
    return 0;
  }
  TString message("Missing configuration handler, failed to register configuration");
  if (name) {message+=" '"; message+=name;}
  message+="'\n AliHLTSystem and configuration handler can be initialized by adding the line";
  message+="\n    AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();";
  message+="\n to the macro before the first AliHLTConfiguration definition. Suppressing further messages.\n";
  static AliHLTErrorGuard g("AliHLTConfigurationHandler", "MissedRegistration", message.Data());
  (++g).Throw(1);
  return 1;
}

int AliHLTConfigurationHandler::AddSubstitution(const char* componentId, const AliHLTConfiguration& subst)
{
  /// add component substitution for components of specified id
  if (!componentId) return -EINVAL;
  if (!fgpSubstitutions) fgpSubstitutions=new TMap;
  if (!fgpSubstitutions) return -ENOMEM;
  fgpSubstitutions->SetOwnerKeyValue(kTRUE);

  fgpSubstitutions->Add(new TObjString(componentId), new AliHLTConfiguration(subst));

  return 0;  
}

int AliHLTConfigurationHandler::AddSubstitution(const AliHLTConfiguration& conf , const AliHLTConfiguration& subst)
{
  /// add component substitution for components of specified id
  if (!fgpSubstitutions) fgpSubstitutions=new TMap;
  if (!fgpSubstitutions) return -ENOMEM;
  fgpSubstitutions->SetOwnerKeyValue(kTRUE);

  fgpSubstitutions->Add(new AliHLTConfiguration(conf), new AliHLTConfiguration(subst));

  return 0;  
}

const AliHLTConfiguration* AliHLTConfigurationHandler::FindSubstitution(const AliHLTConfiguration& conf)
{
  /// find component substitution for a configuration
  if (!fgpSubstitutions) return NULL;
  TObject* value=NULL;

  // check for specific configuration
  value=fgpSubstitutions->GetValue(conf.GetName());
  if (value) return dynamic_cast<AliHLTConfiguration*>(value);

  // check for component Id
  value=fgpSubstitutions->GetValue(conf.GetComponentID());
  if (value) return dynamic_cast<AliHLTConfiguration*>(value);

  return NULL;
}
