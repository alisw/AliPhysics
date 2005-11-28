// $Id$

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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLT configuration handling                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include <cerrno>
#include "AliHLTConfiguration.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include <iostream>
#include <string>

ClassImp(AliHLTConfiguration)

/* the global configuration handler which is used to automatically register the configuration
 */
AliHLTConfigurationHandler* AliHLTConfiguration::fConfigurationHandler=NULL;

AliHLTConfiguration::AliHLTConfiguration()
{ 
  fID=NULL;
  fComponent=NULL;
  fStringSources=NULL;
  fNofSources=-1;
  fArguments=NULL;
  fArgc=-1;
  fArgv=NULL;
  fListSrcElement=fListSources.begin();
}

AliHLTConfiguration::AliHLTConfiguration(const char* id, const char* component, const char* sources, const char* arguments)
{
  fArgc=-1;
  fArgv=NULL;
  
  if (id && component) {
    fID=id;
    fComponent=component;
    fStringSources=sources;
    fNofSources=-1;
    fArguments=arguments;
    fListSrcElement=fListSources.begin();
    if (fConfigurationHandler) {
      fConfigurationHandler->RegisterConfiguration(this);
    } else {
      HLTError("no configuration handler set, abort registration");
    }
  }
}

AliHLTConfiguration::~AliHLTConfiguration()
{
  if (fConfigurationHandler) {
    if (fConfigurationHandler->FindConfiguration(fID)!=NULL) {
      fConfigurationHandler->RemoveConfiguration(this);
    }
  }
  if (fArgv != NULL) {
    if (fArgc>0) {
      for (int i=0; i<fArgc; i++) {
	delete[] fArgv[i];
      }
    }
    delete[] fArgv;
    fArgv=NULL;
  }
}

int AliHLTConfiguration::GlobalInit(AliHLTConfigurationHandler* pHandler)
{
  int iResult=0;
  if (fConfigurationHandler!=NULL) {
    fConfigurationHandler->Logging(kHLTLogWarning, "AliHLTConfiguration::GlobalInit", HLT_DEFAULT_LOG_KEYWORD, "configuration handler already initialized, overriding object %p", fConfigurationHandler);
  }
  fConfigurationHandler=pHandler;
  return iResult;
}

int AliHLTConfiguration::GlobalDeinit()
{
  int iResult=0;
  fConfigurationHandler=NULL;
  return iResult;
}

const char* AliHLTConfiguration::GetName() const {
  if (fID)
    return fID;
  return TObject::GetName();
}

AliHLTConfiguration* AliHLTConfiguration::GetSource(const char* id)
{
  AliHLTConfiguration* pSrc=NULL;
  if (id) {
    // first check the current element
    if (fListSrcElement!=fListSources.end() && strcmp(id, (*fListSrcElement)->GetName())==0) {
      pSrc=*fListSrcElement;
      } else {
      // check the list

      pSrc=GetFirstSource();
      while (pSrc) {
	if (strcmp(id, pSrc->GetName())==0)
	  break;
	pSrc=GetNextSource();
      }
    }
  }
  return pSrc;
}

AliHLTConfiguration* AliHLTConfiguration::GetFirstSource()
{
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>=0 || ExtractSources()) {
    fListSrcElement=fListSources.begin();
    if (fListSrcElement!=fListSources.end()) pSrc=*fListSrcElement;
  } 
  return pSrc;
}

AliHLTConfiguration* AliHLTConfiguration::GetNextSource()
{
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>0) {
    if (fListSrcElement!=fListSources.end() && (++fListSrcElement)!=fListSources.end()) 
      pSrc=*fListSrcElement;
  } 
  return pSrc;
}

int AliHLTConfiguration::SourcesResolved(int bAuto) 
{
  int iResult=0;
  if (fNofSources>=0 || bAuto && (iResult=ExtractSources())>=0) {
    //HLTDebug("fNofSources=%d", fNofSources);
    //HLTDebug("list size = %d", fListSources.size());
    iResult=fNofSources==(int)fListSources.size();
  }
  return iResult;
}

int AliHLTConfiguration::InvalidateSource(AliHLTConfiguration* pConf)
{
  int iResult=0;
  if (pConf) {
    vector<AliHLTConfiguration*>::iterator element=fListSources.begin();
    while (element!=fListSources.end()) {
      if (*element==pConf) {
	fListSources.erase(element);
	fListSrcElement=fListSources.end();
	// there is no need to re-evaluate until there was a new configuration registered
	// -> postpone the invalidation, its done in AliHLTConfigurationHandler::RegisterConfiguration
	//InvalidateSources();
	break;
      }
      element++;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTConfiguration::PrintStatus()
{
  HLTLogKeyword("configuration status");
  HLTMessage("status of configuration \"%s\" (%p)", GetName(), this);
  if (fComponent) HLTMessage("  - component: \"%s\"", fComponent);
  else HLTMessage("  - component string invalid");
  if (fStringSources) HLTMessage("  - sources: \"%s\"", fStringSources);
  else HLTMessage("  - no sources");
  if (SourcesResolved(1)<=0)
    HLTMessage("    there are unresolved sources");
  AliHLTConfiguration* pSrc=GetFirstSource();
  while (pSrc) {
    HLTMessage("    source \"%s\" (%p) resolved", pSrc->GetName(), pSrc);
    pSrc=GetNextSource();
  }
}

int AliHLTConfiguration::GetArguments(int* pArgc, const char*** pArgv)
{
  int iResult=0;
  if (pArgc && pArgv) {
    *pArgc=fArgc;
    *pArgv=(const char**)fArgv;
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}


int AliHLTConfiguration::ExtractSources()
{
  int iResult=0;
  fNofSources=0;
  if (fStringSources!=NULL) {
    vector<char*> tgtList;
    fListSources.clear();
    if ((iResult=InterpreteString(fStringSources, tgtList))>=0) {
      fNofSources=tgtList.size();
      vector<char*>::iterator element=tgtList.begin();
      while ((element=tgtList.begin())!=tgtList.end()) {
	if (fConfigurationHandler) {
	  AliHLTConfiguration* pConf=fConfigurationHandler->FindConfiguration(*element);
	  if (pConf) {
	    HLTDebug("source \"%s\" inserted", pConf->GetName());
	    fListSources.push_back(pConf);
	  } else {
	    HLTError("can not find source \"%s\"", (*element));
	    iResult=-ENOENT;
	  }
	} else if (iResult>=0) {
	  iResult=-EFAULT;
	  HLTFatal("global configuration handler not initialized, can not resolve sources");
	}
	delete[] (*element);
	tgtList.erase(element);
      }
      fListSrcElement=fListSources.begin();
    }
  }
  return iResult;
}

int AliHLTConfiguration::ExtractArguments()
{
  int iResult=0;
  if (fArguments!=NULL) {
    vector<char*> tgtList;
    if ((iResult=InterpreteString(fArguments, tgtList))>=0) {
      fArgc=tgtList.size();
      //HLTDebug("found %d arguments", fArgc);
      if (fArgc>0) {
	fArgv = new char*[fArgc];
	if (fArgv) {
	  vector<char*>::iterator element=tgtList.begin();
	  int i=0;
	  while (element!=tgtList.end()) {
	    //HLTDebug("assign arguments %d (%s)", i, *element);
	    fArgv[i++]=(*element);
	    element++;
	  }
	} else {
	  iResult=-ENOMEM;
	}
      }
    }
  }
  return iResult;
}

int AliHLTConfiguration::InterpreteString(const char* arg, vector<char*>& argList)
{
  int iResult=0;
  if (arg) {
    //HLTDebug("interprete \"%s\"", arg);
    int i=0;
    int prec=-1;
    do {
      if (arg[i]==0 || arg[i]==' ') {
	if (prec>=0) {
	  char* pEntry= new char[i-prec+1];
	  if (pEntry) {
	    strncpy(pEntry, &arg[prec], i-prec);
	    pEntry[i-prec]=0; // terminate string
	    //HLTDebug("create string \"%s\", insert at %d", pEntry, argList.size());
	    argList.push_back(pEntry);
	  } else 
	    iResult=-ENOMEM;
	  prec=-1;
	}
      } else if (prec==-1) prec=i;
    } while (arg[i++]!=0 && iResult>=0); 
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTConfiguration::FollowDependency(const char* id, TList* pTgtList)
{
  int iResult=0;
  if (id) {
    AliHLTConfiguration* pDep=NULL;
    if ((pDep=GetSource(id))!=NULL) {
      if (pTgtList) pTgtList->Add(pDep);
      iResult++;
    } else {
      pDep=GetFirstSource();
      while (pDep && iResult==0) {
	if ((iResult=pDep->FollowDependency(id, pTgtList))>0) {
	  if (pTgtList) pTgtList->AddFirst(pDep);
	  iResult++;
	}
	pDep=GetNextSource();
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliHLTTask)

AliHLTTask::AliHLTTask()
{
  fpConfiguration=NULL;
  fpComponent=NULL;
  fpBlockDataArray=NULL;
}

AliHLTTask::AliHLTTask(AliHLTConfiguration* fConf, AliHLTComponentHandler* pCH)
{
  fpConfiguration=NULL;
  fpComponent=NULL;
  fpBlockDataArray=NULL;
  Init(fConf, pCH);
}

AliHLTTask::~AliHLTTask()
{
  if (fpComponent) delete fpComponent;
  fpComponent=NULL;
  if (fpBlockDataArray) delete[] fpBlockDataArray;
  fpBlockDataArray=NULL;
}

int AliHLTTask::Init(AliHLTConfiguration* fConf, AliHLTComponentHandler* pCH)
{
  int iResult=0;
  if (fConf) {
    fpConfiguration=fConf;
    if (pCH) {
      int argc=0;
      const char** argv=NULL;
      if ((iResult=fConf->GetArguments(&argc, &argv))>=0) {
	iResult=pCH->CreateComponent(fConf->GetComponentID(), NULL, argc, argv, fpComponent);
	if (fpComponent) {
	} else {
	  HLTError("can not find component \"%s\"", fConf->GetComponentID());
	}
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

const char *AliHLTTask::GetName() const
{
  if (fpConfiguration)
    return fpConfiguration->GetName();
  return TObject::GetName();
}

AliHLTConfiguration* AliHLTTask::GetConf()
{
  return fpConfiguration;
}

AliHLTComponent* AliHLTTask::GetComponent()
{
  return fpComponent;
}

AliHLTTask* AliHLTTask::FindDependency(const char* id)
{
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListDependencies.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::FollowDependency(const char* id, TList* pTgtList)
{
  int iResult=0;
  if (id) {
    AliHLTTask* pDep=NULL;
    if ((pDep=(AliHLTTask*)fListDependencies.FindObject(id))!=NULL) {
      if (pTgtList) pTgtList->Add(pDep);
      iResult++;
    } else {
      TObjLink* lnk=fListDependencies.FirstLink();
      while (lnk && iResult==0) {
	pDep=(AliHLTTask*)lnk->GetObject();
	if (pDep) {
	  if ((iResult=pDep->FollowDependency(id, pTgtList))>0) {
	    if (pTgtList) pTgtList->AddFirst(pDep);
	    iResult++;
	  }
	} else {
	  iResult=-EFAULT;
	}
	lnk=lnk->Next();
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTTask::PrintDependencyTree(const char* id, int bFromConfiguration)
{
  HLTLogKeyword("task dependencies");
  int iResult=0;
  TList tgtList;
  if (bFromConfiguration) {
    if (fpConfiguration)
      iResult=fpConfiguration->FollowDependency(id, &tgtList);
    else
      iResult=-EFAULT;
  } else
    iResult=FollowDependency(id, &tgtList);
  if (iResult>0) {
    HLTMessage("     task \"%s\": dependency level %d ", GetName(), iResult);
    TObjLink* lnk=tgtList.FirstLink();
    int i=iResult;
    char* pSpace = new char[iResult+1];
    if (pSpace) {
      memset(pSpace, 32, iResult);
      pSpace[i]=0;
      while (lnk) {
	TObject* obj=lnk->GetObject();
	HLTMessage("     %s^-- %s ", &pSpace[i--], obj->GetName());
	lnk=lnk->Next();
      }
      delete [] pSpace;
    } else {
      iResult=-ENOMEM;
    }
  }
}

int AliHLTTask::InsertBlockData(AliHLTComponent_BlockData* pBlock, AliHLTTask* pSource)
{
  int iResult=0;
  return iResult;
}

int AliHLTTask::SetDependency(AliHLTTask* pDep)
{
  int iResult=0;
  if (pDep) {
    if (FindDependency(pDep->GetName())==NULL) {
      fListDependencies.Add(pDep);
    } else {
      iResult=-EEXIST;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::CheckDependencies()
{
  int iResult=0;
  AliHLTConfiguration* pSrc=fpConfiguration->GetFirstSource();
  while (pSrc) {
    if (FindDependency(pSrc->GetName())==NULL) {
      //HLTDebug("dependency \"%s\" unresolved", pSrc->GetName());
      iResult++;
    }
    pSrc=fpConfiguration->GetNextSource();
  }
  return iResult;
}


int AliHLTTask::Depends(AliHLTTask* pTask)
{
  int iResult=0;
  if (pTask) {
    if (fpConfiguration) {
      iResult=fpConfiguration->GetSource(pTask->GetName())!=NULL;
      if (iResult>0) {
	//HLTDebug("task  \"%s\" depends on \"%s\"", GetName(), pTask->GetName());
      } else {
	//HLTDebug("task  \"%s\" independend of \"%s\"", GetName(), pTask->GetName());
      }
    } else {
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

AliHLTTask* AliHLTTask::FindTarget(const char* id)
{
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListTargets.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::SetTarget(AliHLTTask* pTgt)
{
  int iResult=0;
  if (pTgt) {
    if (FindTarget(pTgt->GetName())==NULL) {
      fListTargets.Add(pTgt);
    } else {
      iResult=-EEXIST;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::BuildBlockDataArray(AliHLTComponent_BlockData*& pTgt)
{
  int iResult=0;
  return iResult;
}


// int AliHLTTask::ProcessTask(...)
// {
//   int iResult=0;
//   return iResult;
// }


int AliHLTTask::ClearSourceBlocks()
{
  int iResult=0;
  return iResult;
}

void AliHLTTask::PrintStatus()
{
  HLTLogKeyword("task properties");
  if (fpComponent) {
    HLTMessage("     component: %s (%p)", fpComponent->GetComponentID(), fpComponent);
  } else {
    HLTMessage("     no component set!");
  }
  if (fpConfiguration) {
    AliHLTConfiguration* pSrc=fpConfiguration->GetFirstSource();
    while (pSrc) {
      const char* pQualifier="unresolved";
      if (FindDependency(pSrc->GetName()))
	pQualifier="resolved";
      HLTMessage("     source: %s (%s)", pSrc->GetName(), pQualifier);
      pSrc=fpConfiguration->GetNextSource();
    }
    TObjLink* lnk = fListTargets.FirstLink();
    while (lnk) {
      TObject *obj = lnk->GetObject();
      HLTMessage("     target: %s", obj->GetName());
      lnk = lnk->Next();
    }
  } else {
    HLTMessage("     task \"%s\" not initialized", GetName());
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TList AliHLTConfigurationHandler::fListConfigurations;
TList AliHLTConfigurationHandler::fListDynamicConfigurations;

ClassImp(AliHLTConfigurationHandler)

AliHLTConfigurationHandler::AliHLTConfigurationHandler()
{
}

AliHLTConfigurationHandler::~AliHLTConfigurationHandler()
{
  TObjLink* lnk=fListDynamicConfigurations.FirstLink();
  while (lnk) {
    TObject* obj=lnk->GetObject();
    if (fListConfigurations.FindObject(obj->GetName())==NULL) {
      HLTDebug("delete dynamic configuration \"%s\"", obj->GetName());
      delete obj;
    }
    lnk=lnk->Next();
  }
}

int AliHLTConfigurationHandler::RegisterConfiguration(AliHLTConfiguration* pConf)
{
  int iResult=0;
  if (pConf) {
    if (FindConfiguration(pConf->GetName()) == NULL) {
      fListConfigurations.Add(pConf);
      //HLTDebug("configuration \"%s\" registered", pConf->GetName());

      // mark all configurations with unresolved dependencies for re-evaluation
      TObjLink* lnk=fListConfigurations.FirstLink();
      while (lnk) {
	AliHLTConfiguration* pSrc=(AliHLTConfiguration*)lnk->GetObject();
	if (pSrc && pSrc!=pConf && pSrc->SourcesResolved()!=1) {
	  pSrc->InvalidateSources();
	}
	lnk=lnk->Next();
      }
    } else {
      iResult=-EEXIST;
      HLTWarning("configuration \"%s\" already registered", pConf->GetName());
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTConfigurationHandler::CreateConfiguration(const char* id, const char* component, const char* sources, const char* arguments)
{
  int iResult=0;
  AliHLTConfiguration* pConf= new AliHLTConfiguration(id, component, sources, arguments);
  if (pConf) {
    // the configuration will be registered automatically, if this failes the configuration
    // is missing -> delete it
    if (FindConfiguration(id)==NULL) {
      delete pConf;
      pConf=NULL;
      iResult=-EEXIST;
    } else {
      fListDynamicConfigurations.Add(pConf);
    }
  } else {
    HLTError("system error: object allocation failed");
    iResult=-ENOMEM;
  }
  return iResult;
}

void AliHLTConfigurationHandler::PrintConfigurations()
{
  HLTLogKeyword("configuration listing");
  HLTMessage("registered configurations:");
  TObjLink *lnk = fListConfigurations.FirstLink();
  while (lnk) {
    TObject *obj = lnk->GetObject();
    HLTMessage("  %s", obj->GetName());
    lnk = lnk->Next();
  }
}

int AliHLTConfigurationHandler::RemoveConfiguration(const char* id)
{
  int iResult=0;
  if (id) {
    AliHLTConfiguration* pConf=NULL;
    if ((pConf=FindConfiguration(id))!=NULL) {
      iResult=RemoveConfiguration(pConf);
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
  int iResult=0;
  if (pConf) {
    // remove the configuration from the list
    fListConfigurations.Remove(pConf);
    // remove cross links in the remaining configurations
    TObjLink* lnk=fListConfigurations.FirstLink();
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
  AliHLTConfiguration* pConf=NULL;
  if (id) {
    pConf=(AliHLTConfiguration*)fListConfigurations.FindObject(id); 
  }
  return pConf;
}

