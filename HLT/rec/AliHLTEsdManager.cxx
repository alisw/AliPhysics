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

/** @file   AliHLTEsdManager.cxx
    @author Matthias Richter
    @date   
    @brief  Manager for merging and writing of HLT ESDs
*/

#include "AliHLTEsdManager.h"
#include "AliHLTComponent.h"
#include "AliESDEvent.h"
#include "AliHLTMessage.h"
#include "AliESDEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "TClass.h"
#include "TObject.h"
#include "TObjectTable.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEsdManager)

AliHLTEsdManager::AliHLTEsdManager()
  :
  fESDs()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTEsdManager::~AliHLTEsdManager()
{
  // see header file for class documentation
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i]) {
      delete fESDs[i];
    }
    fESDs[i]=NULL;
  }
}

AliHLTEsdManager::AliHLTEsdListEntry* AliHLTEsdManager::Find(AliHLTComponentDataType dt) const
{
  // see header file for class documentation
  AliHLTEsdListEntry* pEntry=NULL;
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i] && *(fESDs[i])==dt) {
      pEntry=const_cast<AliHLTEsdListEntry*>(fESDs[i]);
    }
  }
  return pEntry;
}

int AliHLTEsdManager::WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size,
			       AliHLTComponentDataType dt, AliESDEvent* tgtesd, int eventno)
{
  // see header file for class documentation
  if (!pBuffer && size<=0) return -EINVAL;
  int iResult=0;
  AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)pBuffer);
  if (firstWord==size-sizeof(AliHLTUInt32_t)) {
    HLTDebug("create object from block %s size %d", AliHLTComponent::DataType2Text(dt).c_str(), size);
    AliHLTMessage msg(const_cast<AliHLTUInt8_t*>(pBuffer), size);
    TClass* objclass=msg.GetClass();
    TObject* pObj=msg.ReadObject(objclass);
    if (pObj && objclass) {
      HLTDebug("object %p type %s created", pObj, objclass->GetName());
      AliESDEvent* pESD=dynamic_cast<AliESDEvent*>(pObj);
      TTree* pTree=NULL;
      if (!pESD) {
	pTree=dynamic_cast<TTree*>(pObj);
	if (pTree) {
	  pESD=new AliESDEvent;
	  pESD->CreateStdContent();
	  if (pTree->GetEntries()>0) {
	    if (pTree->GetEntries()>1) {
	      HLTWarning("only one entry allowed for ESD embedded into tree, data block %s contains tree with %d entires, taking first entry",
			 AliHLTComponent::DataType2Text(dt).c_str(), pTree->GetEntries());
	    }
	    pESD->ReadFromTree(pTree);
	    pTree->GetEvent(0);
	  }
	} else {
	  HLTWarning("tree of data block %s has no events, skipping data block", AliHLTComponent::DataType2Text(dt).c_str());
	}
      }
      if (pESD) {
	AliHLTEsdListEntry* entry=Find(dt);
	if (!entry) {
	  AliHLTEsdListEntry* newEntry=new AliHLTEsdListEntry(dt);
	  fESDs.push_back(newEntry);
	}
	if (tgtesd) {
	}
	entry=Find(dt);
	if (entry) {
	  entry->WriteESD(pESD, eventno);
	} else {
	  HLTError("internal mismatch, can not create list entry");
	  iResult=-ENOMEM;
	}
      } else {
	HLTWarning("data block %s is not of class type AliESDEvent, ignoring ...", AliHLTComponent::DataType2Text(dt).c_str());
      }
      if (pTree) {
	// ESD has been created and must be cleaned up
	delete pESD;
	pESD=NULL;
      }
      delete pObj;
      pObj=NULL;
    } else {
    }
  }
  return iResult;
}

AliHLTEsdManager::AliHLTEsdListEntry::AliHLTEsdListEntry(AliHLTComponentDataType dt)
  :
  fName(),
  fpFile(NULL),
  fpTree(NULL),
  fpEsd(NULL),
  fDt(dt)
{
  // see header file for class documentation
}

AliHLTEsdManager::AliHLTEsdListEntry::~AliHLTEsdListEntry()
{
  // see header file for class documentation
  if (fpTree) {
    fpTree->GetUserInfo()->Clear();
    delete fpTree;
    fpTree=NULL;
  }

  // due to the Root garbage collection the ESD object might already be
  // deleted since the pTree->GetUserInfo()->Add(pESD) adds the ESD object to
  // an internal list which is cleaned when the tree is deleted
  if (fpEsd && gObjectTable->PtrIsValid(fpEsd)) {
    delete fpEsd;
  }
  fpEsd=NULL;

  if (fpFile) {
    fpFile->Close();
    delete fpFile;
    fpFile=NULL;
  }
}

bool AliHLTEsdManager::AliHLTEsdListEntry::operator==(AliHLTComponentDataType dt) const
{
  // see header file for class documentation
  return fDt==dt;
}

int AliHLTEsdManager::AliHLTEsdListEntry::WriteESD(AliESDEvent* pESD, int eventno)
{
  // see header file for class documentation
  if (!pESD) return -EINVAL;
  int iResult=0;
  if (!fpFile) {
    TString origin;
    origin.Insert(0, fDt.fOrigin, kAliHLTComponentDataTypefOriginSize);
    origin.Remove(TString::kTrailing, ' ');
    origin.ToUpper();
    fName="AliHLT"; fName+=origin;
    if (fDt!=kAliHLTDataTypeESDObject &&
	fDt!=kAliHLTDataTypeESDTree) {

      HLTWarning("non-standard ESD type %s", AliHLTComponent::DataType2Text(fDt).c_str());
      TString id;
      id.Insert(0, fDt.fID, kAliHLTComponentDataTypefIDsize);
      id.Remove(TString::kTrailing, ' ');
      id.ToUpper();
      fName+="_"; fName+=id; fName+=".root";
    } else {
      fName+="ESDs.root";
      fpFile=new TFile(fName, "RECREATE");
    }
  }
  if (fpFile && !fpFile->IsZombie() && iResult>=0) {
    if (!fpTree) {
      fpTree=new TTree("esdTree", "Tree with HLT ESD objects");
      if (!fpTree) {
	iResult=-ENOMEM;
      } else {
	fpTree->SetDirectory(0);
      }
    }
    if (fpTree && iResult>=0) {
      if (!fpEsd) {
	// create the ESD structure for filling into the tree
	fpEsd=new AliESDEvent;
	if (fpEsd) {
	  fpEsd->CreateStdContent();
	  fpTree->GetUserInfo()->Add(fpEsd);
	} else {
	  iResult=-ENOMEM;
	}
      } else {
	fpEsd->ResetStdContent();
      }
      if (eventno>=0) {
	// synchronize and add empty events
	for (int i=fpTree->GetEntries(); i<eventno; i++) {
	  fpTree->Fill();
	}
	if (fpTree->GetEntries()>eventno) {
	  HLTWarning("event %d ESD of type %s already written, skipping additional data block", eventno, AliHLTComponent::DataType2Text(fDt).c_str());
	}
      }
      if (fpEsd) {
	// we need to copy the ESD, I did not find an approptiate
	// method, the workaround is to save the ESD in a temporary
	// tree, read the content back into the ESD structure
	// used for filling
	TTree* dummy=new TTree("dummy","dummy");
	if (dummy) {
	  /*
	  dummy->SetDirectory(0);
	  pESD->WriteToTree(dummy);
	  dummy->Fill();
	  dummy->GetUserInfo()->Add(pESD);
	  fpEsd->ReadFromTree(dummy);
	  dummy->GetEvent(0);
	  */
	  fpEsd->WriteToTree(fpTree);
	  fpTree->Fill();
	  dummy->GetUserInfo()->Clear();
	  delete dummy;
	} else {
	  iResult=-ENOMEM;
	}
      } else {
	iResult=-ENOMEM;
      }
      if (iResult>=0) {
	fpTree->Write("",TObject::kOverwrite);
      }
    }
  }
  return iResult;
}
