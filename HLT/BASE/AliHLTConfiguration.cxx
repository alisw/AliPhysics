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

/** @file   AliHLTConfiguration.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of HLT configuration handler.
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
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTTask.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include <iostream>
#include <string>
#include "TList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTConfiguration)

AliHLTConfiguration::AliHLTConfiguration()
  :
  fID(""),
  fComponent(""),
  fStringSources(""),
  fNofSources(-1),
  fListSources(),
  fListSrcElement(),
  fArguments(""),
  fArgc(-1),
  fArgv(NULL)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  fListSrcElement=fListSources.begin();
}

AliHLTConfiguration::AliHLTConfiguration(const char* id, const char* component, const char* sources, const char* arguments)
  :
  fID(id),
  fComponent(component),
  fStringSources(sources),
  fNofSources(-1),
  fListSources(),
  fListSrcElement(),
  fArguments(arguments),
  fArgc(-1),
  fArgv(NULL)
{
  // see header file for function documentation
  fListSrcElement=fListSources.begin();
  if (id && component) {
    if (fgConfigurationHandler) {
      fgConfigurationHandler->RegisterConfiguration(this);
    } else {
      HLTError("no configuration handler set, abort registration");
    }
  }
}

AliHLTConfiguration::AliHLTConfiguration(const AliHLTConfiguration& src)
  :
  TObject(),
  AliHLTLogging(),
  fID(src.fID),
  fComponent(src.fComponent),
  fStringSources(src.fStringSources),
  fNofSources(-1),
  fListSources(),
  fListSrcElement(),
  fArguments(src.fArguments),
  fArgc(-1),
  fArgv(NULL)
{ 
  // see header file for function documentation
  fListSrcElement=fListSources.begin();
}

AliHLTConfiguration& AliHLTConfiguration::operator=(const AliHLTConfiguration& src)
{ 
  // see header file for function documentation
  fID=src.fID;
  fComponent=src.fComponent;
  fStringSources=src.fStringSources;
  fNofSources=-1;
  fArguments=src.fArguments;
  fArgc=-1;
  fArgv=NULL;
  return *this;
}

AliHLTConfiguration::~AliHLTConfiguration()
{
  // see header file for function documentation
  if (fgConfigurationHandler) {
    if (fgConfigurationHandler->FindConfiguration(fID.Data())!=NULL) {
      fgConfigurationHandler->RemoveConfiguration(this);
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

  vector<AliHLTConfiguration*>::iterator element=fListSources.begin();
  while (element!=fListSources.end()) {
    fListSources.erase(element);
    element=fListSources.begin();
  }
}

/* the global configuration handler which is used to automatically register the configuration
 */
AliHLTConfigurationHandler* AliHLTConfiguration::fgConfigurationHandler=NULL;

int AliHLTConfiguration::GlobalInit(AliHLTConfigurationHandler* pHandler)
{
  // see header file for function documentation
  int iResult=0;
  if (fgConfigurationHandler!=NULL) {
    fgConfigurationHandler->Logging(kHLTLogWarning, "AliHLTConfiguration::GlobalInit", HLT_DEFAULT_LOG_KEYWORD, "configuration handler already initialized, overriding object %p with %p", fgConfigurationHandler, pHandler);
  }
  fgConfigurationHandler=pHandler;
  return iResult;
}

int AliHLTConfiguration::GlobalDeinit(AliHLTConfigurationHandler* pHandler)
{
  // see header file for function documentation
  int iResult=0;
  if (pHandler!=NULL && fgConfigurationHandler!=pHandler) {
    fgConfigurationHandler->Logging(kHLTLogWarning, "AliHLTConfiguration::GlobalDeinit", HLT_DEFAULT_LOG_KEYWORD, "handler %p is not set, skip ...", pHandler);
    return -EBADF;
  }
  fgConfigurationHandler=NULL;
  return iResult;
}

const char* AliHLTConfiguration::GetName() const 
{
  // see header file for function documentation
  if (!fID.IsNull())
    return fID.Data();
  return TObject::GetName();
}

AliHLTConfiguration* AliHLTConfiguration::GetSource(const char* id)
{
  // see header file for function documentation
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
  // see header file for function documentation
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>=0 || ExtractSources()) {
    fListSrcElement=fListSources.begin();
    if (fListSrcElement!=fListSources.end()) pSrc=*fListSrcElement;
  } 
  return pSrc;
}

AliHLTConfiguration* AliHLTConfiguration::GetNextSource()
{
  // see header file for function documentation
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>0) {
    if (fListSrcElement!=fListSources.end() && (++fListSrcElement)!=fListSources.end()) 
      pSrc=*fListSrcElement;
  } 
  return pSrc;
}

int AliHLTConfiguration::SourcesResolved(int bAuto) 
{
  // see header file for function documentation
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
  // see header file for function documentation
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
  // see header file for function documentation
  HLTLogKeyword("configuration status");
  HLTMessage("status of configuration \"%s\" (%p)", GetName(), this);
  if (!fComponent.IsNull()) HLTMessage("  - component: \"%s\"", fComponent.Data());
  else HLTMessage("  - component string invalid");
  if (!fStringSources.IsNull()) HLTMessage("  - sources: \"%s\"", fStringSources.Data());
  else HLTMessage("  - no sources");
  if (SourcesResolved(1)<=0)
    HLTMessage("    there are unresolved sources");
  AliHLTConfiguration* pSrc=GetFirstSource();
  while (pSrc) {
    HLTMessage("    source \"%s\" (%p) resolved", pSrc->GetName(), pSrc);
    pSrc=GetNextSource();
  }
}

int AliHLTConfiguration::GetArguments(const char*** pArgv)
{
  // see header file for function documentation
  int iResult=0;
  if (pArgv) {
    if (fArgc==-1) {
      if ((iResult=ExtractArguments())<0) {
	HLTError("error extracting arguments for configuration %s", GetName());
	fArgc=-EINVAL;
      }
    } else if (fArgc<0) {
      HLTError("previous argument extraction failed");
    }
    //HLTDebug("%s fArgc %d", GetName(), fArgc);
    iResult=fArgc;
    *pArgv=(const char**)fArgv;
  } else {
    HLTError("invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}


int AliHLTConfiguration::ExtractSources()
{
  // see header file for function documentation
  int iResult=0;
  fNofSources=0;
  if (!fStringSources.IsNull()) {
    vector<char*> tgtList;
    fListSources.clear();
    if ((iResult=InterpreteString(fStringSources.Data(), tgtList))>=0) {
      fNofSources=tgtList.size();
      vector<char*>::iterator element=tgtList.begin();
      while ((element=tgtList.begin())!=tgtList.end()) {
	if (fgConfigurationHandler) {
	  AliHLTConfiguration* pConf=fgConfigurationHandler->FindConfiguration(*element);
	  if (pConf) {
	    //HLTDebug("configuration %s (%p): source \"%s\" (%p) inserted", GetName(), this, pConf->GetName(), pConf);
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
  // see header file for function documentation
  int iResult=0;
  if (!fArguments.IsNull()) {
    vector<char*> tgtList;
    if ((iResult=InterpreteString(fArguments, tgtList))>=0) {
      fArgc=tgtList.size();
      //HLTDebug("configuration %s: extracted %d arguments from \"%s\"", GetName(), fArgc, fArguments);
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
  } else {
    // there are zero arguments
    fArgc=0;
  }
  return iResult;
}

int AliHLTConfiguration::InterpreteString(const char* arg, vector<char*>& argList)
{
  // see header file for function documentation
  int iResult=0;
  if (arg) {
    //HLTDebug("interprete \"%s\"", arg);
    int i=0;
    int prec=-1;
    int bQuote=0;
    do {
      //HLTDebug("%d %x", i, arg[i]);
      if (arg[i]=='\'' && bQuote==0) {
	bQuote=1;
      } else if (arg[i]==0 || 
		 (arg[i]==' ' && bQuote==0) ||
		 (arg[i]=='\'' && bQuote==1)) {
	bQuote=0;
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
  // see header file for function documentation
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

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTask)

AliHLTTask::AliHLTTask()
  :
  fpConfiguration(NULL),
  fpComponent(NULL),
  fpDataBuffer(NULL),
  fListTargets(),
  fListDependencies(),
  fBlockDataArray()
{
  // see header file for function documentation
}

AliHLTTask::AliHLTTask(AliHLTConfiguration* pConf)
  :
  fpConfiguration(pConf),
  fpComponent(NULL),
  fpDataBuffer(NULL),
  fListTargets(),
  fListDependencies(),
  fBlockDataArray()
{
  // see header file for function documentation
}

AliHLTTask::~AliHLTTask()
{
  // see header file for function documentation
  TObjLink* lnk=fListDependencies.FirstLink();

  while (lnk!=NULL) {
    AliHLTTask* pTask=(AliHLTTask*)lnk->GetObject();
    pTask->UnsetTarget(this);
    lnk=lnk->Next();
  }
  lnk=fListTargets.FirstLink();

  while (lnk!=NULL) {
    AliHLTTask* pTask=(AliHLTTask*)lnk->GetObject();
    pTask->UnsetDependency(this);
    lnk=lnk->Next();
  }

  if (fpComponent) delete fpComponent;
  fpComponent=NULL;
}

int AliHLTTask::Init(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH)
{
  // see header file for function documentation
  int iResult=0;
  if (fpConfiguration!=NULL && pConf!=NULL && fpConfiguration!=pConf) {
    HLTWarning("overriding existing reference to configuration object %p (%s) by %p",
	       fpConfiguration, GetName(), pConf);
  }
  if (pConf!=NULL) fpConfiguration=pConf;
  if (fpConfiguration) {
    if (pCH) {
      int argc=0;
      const char** argv=NULL;
      if ((iResult=fpConfiguration->GetArguments(&argv))>=0) {
	argc=iResult; // just to make it clear
	// TODO: we have to think about the optional environment parameter,
	// currently just set to NULL. 
	iResult=pCH->CreateComponent(fpConfiguration->GetComponentID(), NULL, argc, argv, fpComponent);
	if (fpComponent || iResult<=0) {
	  //HLTDebug("component %s (%p) created", fpComponent->GetComponentID(), fpComponent); 
	} else {
	  HLTError("can not find component \"%s\" (%d)", fpConfiguration->GetComponentID(), iResult);
	}
      } else {
	HLTError("can not get argument list for configuration %s (%s)", fpConfiguration->GetName(), fpConfiguration->GetComponentID());
	iResult=-EINVAL;
      }
    } else {
      HLTError("component handler instance needed for task initialization");
      iResult=-EINVAL;
    }
  } else {
    HLTError("configuration object instance needed for task initialization");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Deinit()
{
  // see header file for function documentation
  int iResult=0;
  AliHLTComponent* pComponent=GetComponent();
  fpComponent=NULL;
  if (pComponent) {
    //HLTDebug("delete component %s (%p)", pComponent->GetComponentID(), pComponent); 
    pComponent->Deinit();
    delete pComponent;
  } else {
    HLTWarning("task %s (%p) doesn't seem to be in initialized", GetName(), this);
  }
  return iResult;
}

const char *AliHLTTask::GetName() const
{
  // see header file for function documentation
  if (fpConfiguration)
    return fpConfiguration->GetName();
  return TObject::GetName();
}

AliHLTConfiguration* AliHLTTask::GetConf() const
{
  // see header file for function documentation
  return fpConfiguration;
}

AliHLTComponent* AliHLTTask::GetComponent() const
{
  // see header file for function documentation
  return fpComponent;
}

AliHLTTask* AliHLTTask::FindDependency(const char* id)
{
  // see header file for function documentation
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListDependencies.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::FollowDependency(const char* id, TList* pTgtList)
{
  // see header file for function documentation
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
  // see header file for function documentation
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

int AliHLTTask::SetDependency(AliHLTTask* pDep)
{
  // see header file for function documentation
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

int AliHLTTask::UnsetDependency(AliHLTTask* pDep)
{
  // see header file for function documentation
  fListDependencies.Remove(pDep);
  if (fpConfiguration) {
    fpConfiguration->InvalidateSources();
  }
  return 0;
}

int AliHLTTask::CheckDependencies()
{
  // see header file for function documentation
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
  // see header file for function documentation
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
  // see header file for function documentation
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListTargets.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::SetTarget(AliHLTTask* pTgt)
{
  // see header file for function documentation
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

int AliHLTTask::UnsetTarget(AliHLTTask* pTarget)
{
  fListTargets.Remove(pTarget);
  return 0;
}

int AliHLTTask::StartRun()
{
  // see header file for function documentation
  int iResult=0;
  int iNofInputDataBlocks=0;
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent) {
    // determine the number of input data blocks provided from the source tasks
    TObjLink* lnk=fListDependencies.FirstLink();
    while (lnk && iResult>=0) {
      AliHLTTask* pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	if ((iResult=pSrcTask->GetNofMatchingDataTypes(this))>0) {
	  iNofInputDataBlocks+=iResult;
	} else if (iResult==0) {
	  HLTWarning("source task %s (%p) does not provide any matching data type for task %s (%p)", pSrcTask->GetName(), pSrcTask, GetName(), this);
	} else {
	  HLTError("task %s (%p): error getting matching data types for source task %s (%p)", GetName(), this, pSrcTask->GetName(), pSrcTask);
	  iResult=-EFAULT;
	}
      }
      lnk=lnk->Next();
    }
    if (iResult>=0) {
      if (fBlockDataArray.size()>0) {
	HLTWarning("block data array for task %s (%p) was not cleaned", GetName(), this);
	fBlockDataArray.resize(0);
      }

      // component init
      // the initialization of the component is done by the ComponentHandler after creation
      // of the component.
      //iResult=Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );

      // allocate internal task variables for bookkeeping aso.
      // we allocate the BlockData array with at least one member
      if (iNofInputDataBlocks==0) iNofInputDataBlocks=1;
      AliHLTComponentBlockData init;
      memset(&init, 0, sizeof(AliHLTComponentBlockData));
      fBlockDataArray.resize(iNofInputDataBlocks, init);

      // allocate the data buffer, which controls the output buffer and subscriptions
      if (iResult>=0) {
	fpDataBuffer=new AliHLTDataBuffer;
	if (fpDataBuffer!=NULL) {
	  HLTDebug("created data buffer %p for task %s (%p)", fpDataBuffer, GetName(), this);
	  TObjLink* lnk=fListTargets.FirstLink();
	  while (lnk && iResult>=0) {
	    AliHLTTask* pTgtTask=(AliHLTTask*)lnk->GetObject();
	    if (pTgtTask) {
	      if ((iResult=fpDataBuffer->SetConsumer(pTgtTask->GetComponent()))>=0) {
	      }
	    } else {
	      break;
	      iResult=-EFAULT;
	    }
	    lnk=lnk->Next();
	  }
	} else {
	  HLTFatal("can not create data buffer object, memory allocation failed");
	  iResult=-ENOMEM;
	}
      }
    }
  } else {
    HLTError("task %s (%p) does not have a component", GetName(), this);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTask::EndRun()
{
  // see header file for function documentation
  int iResult=0;
  if (fBlockDataArray.size()>0) {
    fBlockDataArray.resize(0);
  } else {
    HLTWarning("task %s (%p) doesn't seem to be in running mode", GetName(), this);
  }
  if (fpDataBuffer) {
    AliHLTDataBuffer* pBuffer=fpDataBuffer;
    fpDataBuffer=NULL;
    delete pBuffer;
  }
  return iResult;
}

int AliHLTTask::ProcessTask(Int_t eventNo)
{
  // see header file for function documentation
  int iResult=0;
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent && fpDataBuffer) {
    HLTDebug("Processing task %s (%p) fpDataBuffer %p", GetName(), this, fpDataBuffer);
    fpDataBuffer->Reset();
    int iSourceDataBlock=0;
    int iInputDataVolume=0;

    AliHLTTask* pSrcTask=NULL;
    TList subscribedTaskList;
    TObjLink* lnk=fListDependencies.FirstLink();

    // subscribe to all source tasks
    while (lnk && iResult>=0) {
      pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	int iMatchingDB=pSrcTask->GetNofMatchingDataBlocks(this);
	if (iMatchingDB>=0 && static_cast<unsigned int>(iMatchingDB)>fBlockDataArray.size()-iSourceDataBlock) {
	  AliHLTComponentBlockData init;
	  memset(&init, 0, sizeof(AliHLTComponentBlockData));
	  fBlockDataArray.resize(iSourceDataBlock+iMatchingDB, init);
	} else {
	  if (iMatchingDB<0) {
	    HLTError("task %s (%p): error getting no of matching data blocks from task %s (%p), error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iMatchingDB);
	    iResult=iMatchingDB;
	    break;
	  } else if (iMatchingDB==0) {
	    HLTDebug("source task %s (%p) does not provide any matching data type for task %s (%p)", pSrcTask->GetName(), pSrcTask, GetName(), this);
	  }
	}
	if ((iResult=pSrcTask->Subscribe(this, &fBlockDataArray[iSourceDataBlock],fBlockDataArray.size()-iSourceDataBlock))>=0) {
	  for (int i=0; i<iResult; i++) {
	    iInputDataVolume+=fBlockDataArray[i+iSourceDataBlock].fSize;
	    // put the source task as many times into the list as it provides data blocks
	    // makes the bookkeeping for the data release easier
	    subscribedTaskList.Add(pSrcTask);
	  }
	  HLTDebug("Task %s (%p) successfully subscribed to %d data block(s) of task %s (%p)", GetName(), this, iResult, pSrcTask->GetName(), pSrcTask);
	  iSourceDataBlock+=iResult;	  
	  iResult=0;
	} else {
	  HLTError("Task %s (%p): subscription to task %s (%p) failed with error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iResult);
	  iResult=-EFAULT;
	}
      } else {
	HLTFatal("fatal internal error in ROOT list handling");
	iResult=-EFAULT;
      }
      lnk=lnk->Next();
    }

    // process the event
    int iNofTrial=0; // repeat processing if component returns -ENOSPC
    AliHLTUInt32_t size=0;
    if (iResult>=0) {
    do {
      long unsigned int iConstBase=0;
      double fInputMultiplier=0;
      if (pComponent->GetComponentType()!=AliHLTComponent::kSink)
	pComponent->GetOutputDataSize(iConstBase, fInputMultiplier);
      if (fInputMultiplier<0) {
	HLTWarning("ignoring negative input multiplier");
	fInputMultiplier=0;
      }
      long unsigned int iOutputDataSize=int(fInputMultiplier*iInputDataVolume) + iConstBase;
      //HLTDebug("task %s: reqired output size %d", GetName(), iOutputDataSize);
      if (iNofTrial>0) {
	// dont process again if the buffer size is the same
	if (size==iOutputDataSize) break;
	HLTInfo("processing task %s again with buffer size %d", GetName(), iOutputDataSize);
      }
      AliHLTUInt8_t* pTgtBuffer=NULL;
      if (iOutputDataSize>0) pTgtBuffer=fpDataBuffer->GetTargetBuffer(iOutputDataSize);
      //HLTDebug("provided raw buffer %p", pTgtBuffer);
      AliHLTComponentEventData evtData;
      AliHLTComponent::FillEventData(evtData);
      evtData.fEventID=(AliHLTEventID_t)eventNo;
      evtData.fBlockCnt=iSourceDataBlock;
      AliHLTComponentTriggerData trigData;
      size=iOutputDataSize;
      AliHLTUInt32_t outputBlockCnt=0;
      AliHLTComponentBlockData* outputBlocks=NULL;
      AliHLTComponentEventDoneData* edd;
      if (pTgtBuffer!=NULL || iOutputDataSize==0) {
	iResult=pComponent->ProcessEvent(evtData, &fBlockDataArray[0], trigData, pTgtBuffer, size, outputBlockCnt, outputBlocks, edd);
	HLTDebug("task %s: component %s ProcessEvent finnished (%d): size=%d blocks=%d", GetName(), pComponent->GetComponentID(), iResult, size, outputBlockCnt);
	if (iResult>=0 && pTgtBuffer && outputBlocks) {
	  iResult=fpDataBuffer->SetSegments(pTgtBuffer, outputBlocks, outputBlockCnt);
	  delete [] outputBlocks; outputBlocks=NULL; outputBlockCnt=0;
	} else {
	  fpDataBuffer->Reset();
	}
      } else {
	HLTError("task %s: no target buffer available", GetName());
	iResult=-EFAULT;
      }
    } while (iResult==-ENOSPC && iNofTrial++<1);
    }

    // now release all buffers which we have subscribed to
    iSourceDataBlock=0;
    lnk=subscribedTaskList.FirstLink();
    while (lnk) {
      pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	int iTempRes=0;
	if ((iTempRes=pSrcTask->Release(&fBlockDataArray[iSourceDataBlock], this))>=0) {
	  HLTDebug("Task %s (%p) successfully released segment of task %s (%p)", GetName(), this, pSrcTask->GetName(), pSrcTask);
	} else {
	  HLTError("Task %s (%p): realease of task %s (%p) failed with error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iTempRes);
	}
      } else {
	HLTFatal("task %s (%p): internal error in ROOT list handling", GetName(), this);
	if (iResult>=0) iResult=-EFAULT;
      }
      subscribedTaskList.Remove(lnk);
      lnk=subscribedTaskList.FirstLink();
      iSourceDataBlock++;
    }
    if (subscribedTaskList.GetSize()>0) {
      HLTError("task %s (%p): could not release all data buffers", GetName(), this);
    }
  } else {
    HLTError("task %s (%p): internal failure (not initialized component %p, data buffer %p)", GetName(), this, fpComponent, fpDataBuffer);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTask::GetNofMatchingDataBlocks(const AliHLTTask* pConsumerTask) const
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->FindMatchingDataBlocks(pConsumerTask->GetComponent(), NULL);
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::GetNofMatchingDataTypes(const AliHLTTask* pConsumerTask) const
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    AliHLTComponent* pComponent=GetComponent();
    if (!pComponent) {
      // init ?
      HLTError("component not initialized");
      iResult=-EFAULT;
    }
    if (pComponent) {
      iResult=pComponent->FindMatchingDataTypes(pConsumerTask->GetComponent(), NULL);
    } else {
      HLTFatal("task initialization failed");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Subscribe(const AliHLTTask* pConsumerTask, AliHLTComponentBlockData* pBlockDesc, int iArraySize)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->Subscribe(pConsumerTask->GetComponent(), pBlockDesc, iArraySize);
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTTask* pConsumerTask)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask && pBlockDesc) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->Release(pBlockDesc, pConsumerTask->GetComponent());
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTTask::PrintStatus()
{
  // see header file for function documentation
  HLTLogKeyword("task properties");
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent) {
    HLTMessage("     component: %s (%p)", pComponent->GetComponentID(), pComponent);
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

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTConfigurationHandler)

AliHLTConfigurationHandler::AliHLTConfigurationHandler()
  :
  fgListConfigurations()
{
  // see header file for function documentation
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
}

int AliHLTConfigurationHandler::RegisterConfiguration(AliHLTConfiguration* pConf)
{
  // see header file for function documentation
  int iResult=0;
  if (pConf) {
    if (FindConfiguration(pConf->GetName()) == NULL) {
      AliHLTConfiguration* pClone=new AliHLTConfiguration(*pConf);
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
  // see header file for function documentation
  int iResult=0;
  AliHLTConfiguration* pConf= new AliHLTConfiguration(id, component, sources, arguments);
  if (pConf) {
    // the configuration will be registered automatically, if this failes the configuration
    // is missing -> delete it
    if (FindConfiguration(id)==NULL) {
      delete pConf;
      pConf=NULL;
      iResult=-EEXIST;
    }
  } else {
    HLTError("system error: object allocation failed");
    iResult=-ENOMEM;
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

