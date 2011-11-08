// $Id$

///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
///*                  for The ALICE HLT Project.                            *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          *
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************

/// @file   AliHLTConfiguration.cxx
/// @author Matthias Richter
/// @date   2007
/// @brief  HLT configuration description for a single component.
/// @note   The class is used in Offline (AliRoot) context

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
  fListSrcElementIdx(-1),
  fArguments(""),
  fArgc(-1),
  fArgv(NULL),
  fBufferSize(-1)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTConfiguration::AliHLTConfiguration(const char* id, const char* component, const char* sources,
					 const char* arguments, const char* bufsize)
  :
  fID(id),
  fComponent(component),
  fStringSources(sources),
  fNofSources(-1),
  fListSources(),
  fListSrcElementIdx(-1),
  fArguments(arguments),
  fArgc(-1),
  fArgv(NULL),
  fBufferSize(-1)
{
  // see header file for function documentation
  if (bufsize) fBufferSize=ConvertSizeString(bufsize);
  if (id && component) {
    if (AliHLTConfigurationHandler::Instance()) {
      AliHLTConfigurationHandler::Instance()->RegisterConfiguration(this);
    } else {
      AliHLTConfigurationHandler::MissedRegistration(id);
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
  fListSrcElementIdx(-1),
  fArguments(src.fArguments),
  fArgc(-1),
  fArgv(NULL),
  fBufferSize(src.fBufferSize)
{ 
  // see header file for function documentation
}

AliHLTConfiguration& AliHLTConfiguration::operator=(const AliHLTConfiguration& src)
{ 
  // see header file for function documentation
  if (this==&src) return *this;

  fID=src.fID;
  fComponent=src.fComponent;
  fStringSources=src.fStringSources;
  fNofSources=-1;
  fArguments=src.fArguments;
  fArgc=-1;
  fArgv=NULL;
  fBufferSize=src.fBufferSize;
  return *this;
}

AliHLTConfiguration::~AliHLTConfiguration()
{
  // see header file for function documentation
  if (AliHLTConfigurationHandler::Instance()) {
    if (AliHLTConfigurationHandler::Instance()->FindConfiguration(fID.Data())!=NULL) {
      // remove the configuration from the handler if it exists
      // but DO NOT remove the clone configuration
      AliHLTConfigurationHandler::Instance()->RemoveConfiguration(this);
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
    if (fListSrcElementIdx>=0 && fListSrcElementIdx<(int)fListSources.size() &&
	strcmp(id, (fListSources[fListSrcElementIdx])->GetName())==0) {
      pSrc=fListSources[fListSrcElementIdx];
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

AliHLTConfiguration* AliHLTConfiguration::GetFirstSource() const
{
  // see header file for function documentation
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>0) {
    const_cast<AliHLTConfiguration*>(this)->fListSrcElementIdx=-1;
    pSrc=GetNextSource();
  } 
  return pSrc;
}

AliHLTConfiguration* AliHLTConfiguration::GetNextSource() const
{
  // see header file for function documentation
  AliHLTConfiguration* pSrc=NULL;
  if (fNofSources>0) {
    if (fListSrcElementIdx+1<(int)fListSources.size()) {
      const_cast<AliHLTConfiguration*>(this)->fListSrcElementIdx++;
      pSrc=fListSources[fListSrcElementIdx];
    }
  } 
  return pSrc;
}

int AliHLTConfiguration::SourcesResolved() const
{
  // see header file for function documentation
  int iResult=0;
  if (fNofSources>=0) {
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
	fListSrcElementIdx=fListSources.size();
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

void AliHLTConfiguration::PrintStatus() const
{
  // see header file for function documentation
  HLTLogKeyword("configuration status");
  HLTMessage("status of configuration \"%s\" (%p)", GetName(), this);
  if (!fComponent.IsNull()) HLTMessage("  - component: \"%s\"", fComponent.Data());
  else HLTMessage("  - component string invalid");
  if (!fStringSources.IsNull()) HLTMessage("  - sources: \"%s\"", fStringSources.Data());
  else HLTMessage("  - no sources");
  if (SourcesResolved()!=1)
    HLTMessage("    there are unresolved sources");
  AliHLTConfiguration* pSrc=GetFirstSource();
  while (pSrc) {
    HLTMessage("    source \"%s\" (%p) resolved", pSrc->GetName(), pSrc);
    pSrc=GetNextSource();
  }
}

void AliHLTConfiguration::Print(const char* option) const
{
  // print information
  if (option && strcmp(option, "status")==0) {
    PrintStatus();
    return;
  }
  HLTLogKeyword("configuration");
  HLTMessage("configuration %s: component %s, sources %s, arguments %s",
	     GetName(),
	     GetComponentID(),
	     GetSourceSettings(),
	     GetArgumentSettings()
	     );
}

int AliHLTConfiguration::GetArguments(const char*** pArgv) const
{
  // see header file for function documentation
  int iResult=0;
  if (pArgv) {
    if (fArgc==-1) {
      if ((iResult=const_cast<AliHLTConfiguration*>(this)->ExtractArguments())<0) {
	HLTError("error extracting arguments for configuration %s", GetName());
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


int AliHLTConfiguration::ExtractSources(AliHLTConfigurationHandler* pHandler)
{
  // see header file for function documentation
  int iResult=0;
  fNofSources=0; // indicates that the function was called, there are either n or 0 sources
  fListSources.clear();
  if (!pHandler) {
    HLTError("configuration handler missing, can not resolve sources");
    return -EFAULT;
  }
  if (!fStringSources.IsNull()) {
    vector<char*> tgtList;
    if ((iResult=InterpreteString(fStringSources.Data(), tgtList))>=0) {
      fNofSources=tgtList.size();
      vector<char*>::iterator element=tgtList.begin();
      while ((element=tgtList.begin())!=tgtList.end()) {
	  AliHLTConfiguration* pConf=pHandler->FindConfiguration(*element);
	  if (pConf) {
	    //HLTDebug("configuration %s (%p): source \"%s\" (%p) inserted", GetName(), this, pConf->GetName(), pConf);
	    fListSources.push_back(pConf);
	  } else {
	    HLTError("can not find source \"%s\"", (*element));
	    iResult=-ENOENT;
	  }
	delete[] (*element);
	tgtList.erase(element);
      }
    }
  }
  fListSrcElementIdx=-1;
  return iResult<0?iResult:SourcesResolved();
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
  if (iResult<0) fArgc=iResult;
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

int AliHLTConfiguration::ConvertSizeString(const char* strSize) const
{
  // see header file for function documentation
  int size=0;
  if (!strSize) return -1;

  char* endptr=NULL;
  size=strtol(strSize, &endptr, 10);
  if (size>=0) {
    if (endptr) {
      if (endptr==strSize) {
	HLTWarning("ignoring unrecognized buffer size '%s'", strSize);
	size=-1;
      } else if (*endptr==0) {
	// no unit specifier
      } else if (*endptr=='k') {
	size*=1014;
      } else if (*endptr=='M') {
	size*=1024*1024;
      } else {
	HLTWarning("ignoring buffer size of unknown unit '%c'", endptr[0]);
      }
    } else {
      HLTWarning("ignoring negative buffer size specifier '%s'", strSize);
      size=-1;
    }
  }
  return size;
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


