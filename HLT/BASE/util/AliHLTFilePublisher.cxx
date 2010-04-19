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

/** @file   AliHLTFilePublisher.cxx
    @author Matthias Richter
    @date   
    @brief  HLT file publisher component implementation. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTFilePublisher.h"
#include "AliLog.h"
//#include <TObjString.h>
#include <TMath.h>
#include <TFile.h>


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFilePublisher)

AliHLTFilePublisher::AliHLTFilePublisher()
  :
  AliHLTDataSource(),
  fpCurrent(NULL),
  fEvents(),
  fMaxSize(0),
  fOpenFilesAtStart(false),
  fOutputDataTypes(),
  fIsRaw(kTRUE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // make the lists owners of their objects in order to automatically
  // de-allocate the objects
  fEvents.SetOwner();
}

AliHLTFilePublisher::~AliHLTFilePublisher()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

const char* AliHLTFilePublisher::GetComponentID()
{
  // see header file for class documentation
  return "FilePublisher";
}

AliHLTComponentDataType AliHLTFilePublisher::GetOutputDataType()
{
  // see header file for class documentation
  if (fOutputDataTypes.size()==0) return kAliHLTVoidDataType;
  else if (fOutputDataTypes.size()==1) return fOutputDataTypes[0];
  return kAliHLTMultipleDataType;
}

int AliHLTFilePublisher::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.assign(fOutputDataTypes.begin(), fOutputDataTypes.end());
  HLTInfo("%s %p provides %d output data types", GetComponentID(), this, fOutputDataTypes.size());
  return fOutputDataTypes.size();
}

void AliHLTFilePublisher::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=fMaxSize;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTFilePublisher::Spawn()
{
  // see header file for class documentation
  return new AliHLTFilePublisher;
}

int AliHLTFilePublisher::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  //HLTDebug("%d %s", argc, argv[0]);
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  int bHaveDatatype=0;
  int bHaveSpecification=0;
  fOpenFilesAtStart = false;
  AliHLTComponentDataType currDataType=kAliHLTVoidDataType;
  AliHLTUInt32_t          currSpecification=kAliHLTVoidDataSpec;
  EventFiles*             pCurrEvent=NULL;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -datafile
    if (argument.CompareTo("-datafile")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if (!bHaveDatatype) {
	HLTWarning("no data type available so far, please set data type and specification before the file name. The first available data type will be set for all files preceding it");
      }
      FileDesc* pDesc=new FileDesc(argv[i], currDataType, currSpecification, fIsRaw);
      if (pDesc) {
	iResult=InsertFile(pCurrEvent, pDesc);
      } else {
	iResult=-ENOMEM;
      }

      // -datafilelist
    } else if (argument.CompareTo("-datafilelist")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      HLTWarning("-datafilelist option not yet implemented");

      // -datatype
    } else if (argument.CompareTo("-datatype")==0) {
      currDataType=kAliHLTVoidDataType;
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&currDataType.fID, argv[i], TMath::Min(kAliHLTComponentDataTypefIDsize, (Int_t)strlen(argv[i])));
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&currDataType.fOrigin, argv[i], TMath::Min(kAliHLTComponentDataTypefOriginSize, (Int_t)strlen(argv[i])));

      // add all different data types to the list
      AliHLTComponentDataTypeList::iterator element=fOutputDataTypes.begin();
      while (element!=fOutputDataTypes.end() && *element!=currDataType) element++;
      if (element==fOutputDataTypes.end()) fOutputDataTypes.push_back(currDataType);

      if (bHaveDatatype==0 && pCurrEvent && iResult>=0) {
	// this is a workaround to make old tutorials working which contain
	// the arguments in the wrong sequence
	TList& files=*pCurrEvent; // type conversion operator defined
	TObjLink *flnk=files.FirstLink();
	while (flnk) {
	  FileDesc* pFileDesc=dynamic_cast<FileDesc*>(flnk->GetObject());
	  if (pFileDesc) {
	    pFileDesc->SetDataType(currDataType);
	  }
	  flnk=flnk->Next();
	}
      }
      bHaveDatatype=1;

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	currSpecification=(AliHLTUInt32_t)parameter.Atoi();
      } else if (parameter.BeginsWith("0x") &&
		 parameter.Replace(0,2,"",0).IsHex()) {
	sscanf(parameter.Data(),"%x", &currSpecification);
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      if (bHaveSpecification==0 && pCurrEvent && iResult>=0) {
	// this is a workaround to make old tutorials working which contain
	// the arguments in the wrong sequence
	TList& files=*pCurrEvent; // type conversion operator defined
	TObjLink *flnk=files.FirstLink();
	while (flnk) {
	  FileDesc* pFileDesc=dynamic_cast<FileDesc*>(flnk->GetObject());
	  if (pFileDesc) {
	    pFileDesc->SetSpecification(currSpecification);
	  }
	  flnk=flnk->Next();
	}
      }
      bHaveSpecification=1;
      // -nextevent
    } else if (argument.CompareTo("-nextevent")==0) {
      InsertEvent(pCurrEvent);
    } else if (argument.CompareTo("-open_files_at_start")==0) {
      fOpenFilesAtStart = true;
    } else {
      if ((iResult=ScanArgument(argc-i, &argv[i]))==-EINVAL) {
	HLTError("unknown argument %s", argument.Data());
	break;
      } else if (iResult==-EPROTO) {
	bMissingParam=1;
	break;
      } else if (iResult>=0) {
	i+=iResult;
	iResult=0;
      }
    }
  }
  InsertEvent(pCurrEvent);

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (fEvents.GetSize()==0) {
    HLTError("the publisher needs at least one file argument");
    iResult=-EINVAL;
  }
  if (iResult>=0) iResult=OpenFiles(fOpenFilesAtStart);
  if (iResult<0) {
    fEvents.Clear();
  }
  return iResult;
}

int AliHLTFilePublisher::InsertFile(EventFiles* &pCurrEvent, FileDesc* pDesc)
{
  // see header file for class documentation
  int iResult=0;
  if (pDesc) {
    if (pCurrEvent==NULL) {
      pCurrEvent=new EventFiles;
      if (pCurrEvent!=NULL) {
      } else {
	iResult=-ENOMEM;
      }
    }
    if (iResult>=0 && pCurrEvent!=NULL) {
      HLTDebug("Insert file %p to event %p", pDesc, pCurrEvent);
      pCurrEvent->Add(pDesc);
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTFilePublisher::InsertEvent(EventFiles* &pEvent)
{
  // see header file for class documentation
  int iResult=0;
  if (pEvent) {
    HLTDebug("Inserted event %p", pEvent);
    fEvents.Add(pEvent);
    pEvent=NULL;
  }
  return iResult;
}

int AliHLTFilePublisher::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation

  // there are no other arguments than the standard ones
  if (argc==0 && argv==NULL) {
    // this is just to get rid of the warning "unused parameter"
  }
  return -EINVAL;
}

int AliHLTFilePublisher::OpenFiles(bool keepOpen)
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fEvents.FirstLink();
  while (lnk && iResult>=0) {
    EventFiles* pEventDesc=dynamic_cast<EventFiles*>(lnk->GetObject());
    if (pEventDesc) {
      HLTDebug("open files for event %p", pEventDesc);
      TList& files=*pEventDesc; // type conversion operator defined
      TObjLink *flnk=files.FirstLink();
      int eventSize=0;
      while (flnk && iResult>=0) {
	FileDesc* pFileDesc=dynamic_cast<FileDesc*>(flnk->GetObject());
	if (pFileDesc) {
	  int size=pFileDesc->OpenFile();
	  if (not keepOpen) pFileDesc->CloseFile();
	  if (size<0) {
	    iResult=size;
	    HLTError("can not open file %s", pFileDesc->GetName());
	  } else {
	    eventSize+=size;
	  }
	}
	flnk=flnk->Next();
      }
      HLTDebug("event %p size %d", pEventDesc, eventSize);
      if (fMaxSize<eventSize) fMaxSize=eventSize;
    } else {
      HLTError("can not get event descriptor for TObjLink");
    }
    lnk = lnk->Next();
  }

  return iResult;
}

int AliHLTFilePublisher::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  fEvents.Clear();
  return iResult;
}

int AliHLTFilePublisher::GetEvent( const AliHLTComponentEventData& /*evtData*/,
				   AliHLTComponentTriggerData& /*trigData*/,
				   AliHLTUInt8_t* outputPtr, 
				   AliHLTUInt32_t& size,
				   AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation

  // process data events only
  if (!IsDataEvent()) return 0;
  AliHLTUInt32_t capacity=size;
  size=0;

  int iResult=0;
  TObjLink *lnk=fpCurrent;
  if (lnk==NULL) lnk=fEvents.FirstLink();
  fpCurrent=lnk;
  if (lnk) {
    EventFiles* pEventDesc=dynamic_cast<EventFiles*>(lnk->GetObject());
    if (pEventDesc) {
      HLTDebug("publishing files for event %p", pEventDesc);
      TList& files=*pEventDesc; // type conversion operator defined
      TObjLink *flnk=files.FirstLink();
      int iTotalSize=0;
      while (flnk && iResult>=0) {
	FileDesc* pFileDesc=dynamic_cast<FileDesc*>(flnk->GetObject());
	if (not fOpenFilesAtStart) pFileDesc->OpenFile();
	TFile* pFile=NULL;
	if (pFileDesc && (pFile=*pFileDesc)!=NULL) {
	  int iCopy=pFile->GetSize();
	  pFile->Seek(0);
	  if (iCopy+iTotalSize<=(int)capacity) {
	    if (pFile->ReadBuffer((char*)outputPtr+iTotalSize, iCopy)!=0) {
	      // ReadBuffer returns 1 in case of failure and 0 in case of success
	      iResult=-EIO;
	    } else {
	      AliHLTComponentBlockData bd;
	      FillBlockData(bd);
	      bd.fPtr=outputPtr;
	      bd.fOffset=iTotalSize;
	      bd.fSize=iCopy;
	      bd.fDataType=*pFileDesc;      // type conversion operator defined
	      bd.fSpecification=*pFileDesc; // type conversion operator defined
	      outputBlocks.push_back(bd);
	      iTotalSize+=iCopy;
// 	      TString msg;
// 	      msg.Form("get file %s ", pFile->GetName());
// 	      msg+="data type \'%s\'";
// 	      PrintDataTypeContent(bd.fDataType, msg.Data());
	    }
	  } else {
	    // output buffer too small, update GetOutputDataSize for the second trial
	    fMaxSize=iCopy;
	    iResult=-ENOSPC;
	  }
	  if (not fOpenFilesAtStart) pFileDesc->CloseFile();
	} else {
	  HLTError("no file available");
	  iResult=-EFAULT;
	}
	flnk=flnk->Next();
      }
      size=iTotalSize;
    } else {
      HLTError("can not get event descriptor from list link");
      iResult=-EFAULT;
    }
  } else {
    iResult=-ENOENT;
  }
  if (iResult>=0 && fpCurrent) fpCurrent=fpCurrent->Next();

  return iResult;
}

// AliHLTComponentDataType AliHLTFilePublisher::GetCurrentDataType() const
// {
//   return kAliHLTVoidDataType;
// }

// AliHLTUInt32_t          AliHLTFilePublisher::GetCurrentSpecification() const
// {
//   return 0;
// }

AliHLTFilePublisher::FileDesc::FileDesc(const char* name, AliHLTComponentDataType dt, AliHLTUInt32_t spec, Bool_t isRaw)
  :
  TObject(),
  fIsRaw(isRaw),
  fName(name),
  fpInstance(NULL),
  fDataType(dt),
  fSpecification(spec)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTFilePublisher::FileDesc::~FileDesc()
{
  // see header file for class documentation
  CloseFile();
}

void AliHLTFilePublisher::FileDesc::CloseFile()
{
  // see header file for class documentation
  if (fpInstance)
  {
    // Unfortunately had to use AliLog mechanisms rather that AliHLTLogging because
    // AliHLTFilePublisher::FileDesc does not derive from AliHLTLogging. It would
    // become a rather heavy class if it did.
#ifdef __DEBUG
    AliDebugGeneral("AliHLTFilePublisher::FileDesc",
      2, Form("File %s has been closed.", fName.Data())
    );
#endif
    delete fpInstance;
  }
  fpInstance=NULL;
}

int AliHLTFilePublisher::FileDesc::OpenFile()
{
  // see header file for class documentation
  int iResult=0;
  
  TString fullFN="";

  if ( fIsRaw ) fullFN = fName + "?filetype=raw";
  else fullFN = fName;

  fpInstance = TFile::Open(fullFN);
  if (fpInstance) {
    if (fpInstance->IsZombie()==0) {
      iResult=fpInstance->GetSize();
#ifdef __DEBUG
      AliDebugGeneral("AliHLTFilePublisher::FileDesc",
        2, Form("File %s has been opened.", fName.Data())
      );
#endif
    } else {
      iResult=-ENOENT;
    }
  }
  return iResult;
}
