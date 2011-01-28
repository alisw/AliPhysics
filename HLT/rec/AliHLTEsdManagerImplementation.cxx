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

/// @file   AliHLTEsdManagerImplementation.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Manager for merging and writing of HLT ESDs
///         This is an implementation of the abstract interface AliHLTEsdManager

#include "AliHLTEsdManagerImplementation.h"
#include "AliHLTComponent.h"
#include "AliESDEvent.h"
#include "AliHLTMessage.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "AliESDCaloCells.h"
#include "AliMultiplicity.h"
#include "AliESDACORDE.h"
#include "TFile.h"
#include "TTree.h"
#include "TClass.h"
#include "TObject.h"
#include "TList.h"

const float kAliESDVZERODefaultTime = -1024.;
const float kAliESDVZERODefaultTimeGlitch = 1e-6;
const float kAliESDZDCDefaultEMEnergy = 0.;
const float kAliESDZDCDefaultEMEnergyGlitch = 1e-6;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEsdManagerImplementation)

AliHLTEsdManagerImplementation::AliHLTEsdManagerImplementation()
  :
  fESDs()
  , fDirectory()
  , fTreeName()
  , fWriteLocal(false)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  CheckClassConditions();
}

AliHLTEsdManagerImplementation::~AliHLTEsdManagerImplementation()
{
  // see header file for class documentation
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i]) {
      delete fESDs[i];
    }
    fESDs[i]=NULL;
  }
}

int AliHLTEsdManagerImplementation::SetOption(const char* option)
{
  // see header file for class documentation
  int iResult=0;
  TString strOptions=option;
  TObjArray* pTokens=strOptions.Tokenize(" ");
  if (pTokens) {
    if (pTokens->GetEntriesFast()>0) {
      for (int n=0; n<pTokens->GetEntriesFast(); n++) {
	TString data=((TObjString*)pTokens->At(n))->GetString();
	if (data.IsNull()) continue;

	if (data.CompareTo("-writelocal")==0) {
	  fWriteLocal=true;
	} else if (data.Contains("-directory=")) {
	  data.ReplaceAll("-directory=", "");
	  SetDirectory(data.Data());
	} else if (data.Contains("-treename=")) {
	  data.ReplaceAll("-treename=", "");
	  fTreeName=data;
	} else {
	  HLTError("unknown argument %s", data.Data());
	  iResult=-EINVAL;
	  break;
	}
      }
    }
    delete pTokens;
  }
  return iResult;
}

AliHLTEsdManagerImplementation::AliHLTEsdListEntry* AliHLTEsdManagerImplementation::Find(AliHLTComponentDataType dt) const
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

int AliHLTEsdManagerImplementation::WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size,
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
	  if ((entry=new AliHLTEsdListEntry(dt))!=NULL) {
	    if (!fDirectory.IsNull()) {
	      entry->SetDirectory(fDirectory);
	    }
	    if (!fTreeName.IsNull()) {
	      entry->SetTreeName(fTreeName);
	    }
	    fESDs.push_back(entry);
	  }
	}
	if (entry) {
	  if (tgtesd) {
	    Merge(tgtesd, pESD);
	  }

	  // Matthias 2009-06-06: writing of individual ESD files for the different origins was a
	  // first attempt when functionality was missing in the AliRoot framework and remained as
	  // debugging feature. ESD merging is now implemented and data written to the hltEsd, so
	  // the feature is now disabled by default because it causes increasing memory consumption.
	  // Presumably not because of a memory leak but the way the internal TTree is used and kept
	  // in memory.
	  // Writing of local files can be optionally switched on as e.g. by the EsdCollector component.
	  if (fWriteLocal) entry->WriteESD(pESD, eventno);
	} else {
	  HLTError("internal mismatch, can not create list entry");
	  iResult=-ENOMEM;
	}
      } else {
	HLTWarning("data block %s is not of class type AliESDEvent, ignoring ...", AliHLTComponent::DataType2Text(dt).c_str());
      }
      if (pTree) {
	// ESD has been created and must be cleaned up
	pESD->Reset();
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

int AliHLTEsdManagerImplementation::PadESDs(int eventno)
{
  // see header file for class documentation
  int iResult=0;
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i]) {
      int res=fESDs[i]->WriteESD(NULL, eventno);
      if (res<0 && iResult>=0) iResult=res;
    }
  }
  return iResult;
}

void AliHLTEsdManagerImplementation::SetDirectory(const char* directory)
{
  // see header file for class documentation
  if (!directory) return;
  fDirectory=directory;
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i]) {
      fESDs[i]->SetDirectory(directory);
    }
  }
}

TString AliHLTEsdManagerImplementation::GetFileNames(AliHLTComponentDataType dt) const
{
  TString result;
  for (unsigned int i=0; i<fESDs.size(); i++) {
    if (fESDs[i] && *(fESDs[i])==dt) {
      if (!result.IsNull()) result+=" ";
      result+=fESDs[i]->GetFileName();
    }
  }
  return result;
}

TTree* AliHLTEsdManagerImplementation::EmbedIntoTree(AliESDEvent* pESD, const char* name, const char* title)
{
  // see header file for class documentation
  int iResult=0;
  TTree* pTree=new TTree(name, title);
  if (pTree) {
    pESD->WriteToTree(pTree);
    pTree->Fill();
    pTree->GetUserInfo()->Add(pESD);
  } else {
    iResult=-ENOMEM;
  }

  if (iResult<0) {
    pTree->GetUserInfo()->Clear();
    delete pTree;
  }

  return pTree;
}

AliHLTEsdManagerImplementation::AliHLTEsdListEntry::AliHLTEsdListEntry(AliHLTComponentDataType dt)
  :
  fName(),
  fDirectory(),
  fDt(dt),
  fpFile(NULL),
  fpTree(NULL),
  fpEsd(NULL),
  fPrefix(),
  fTreeName("esdTree")
{
  // see header file for class documentation
}

AliHLTEsdManagerImplementation::AliHLTEsdListEntry::~AliHLTEsdListEntry()
{
  // see header file for class documentation
  if (fpEsd) delete fpEsd;
  fpEsd=NULL;

  if (fpTree) delete fpTree;
  fpTree=NULL;

  if (fpFile) {
    fpFile->Close();
    delete fpFile;
  }
  fpFile=NULL;
}

bool AliHLTEsdManagerImplementation::AliHLTEsdListEntry::operator==(AliHLTComponentDataType dt) const
{
  // see header file for class documentation
  return fDt==dt;
}

int AliHLTEsdManagerImplementation::AliHLTEsdListEntry::WriteESD(AliESDEvent* pSrcESD, int eventno)
{
  // see header file for class documentation
  int iResult=0;

  if (fName.IsNull()) {
    // this is the first event, create the file name
    fName="";
    if (!fDirectory.IsNull()) {
      fName+=fDirectory; fName+="/";
    }
    fName+="Ali"; fName+=GetPrefix();
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
    }

    fpFile=new TFile(fName, "RECREATE");
    fpTree=new TTree(fTreeName, "Tree with HLT ESD objects");
    fpTree->SetDirectory(0);
    fpEsd=new AliESDEvent;
    if (fpEsd) {
      fpEsd->CreateStdContent();
      *fpEsd=*pSrcESD;
      if (fpTree) {
	fpEsd->WriteToTree(fpTree);
      }
    }
  }

  if (fpFile && fpTree && fpEsd) {
    // synchronize and add empty events
    fpEsd->Reset();
    int nofCurrentEvents=fpTree->GetEntries();
    if (nofCurrentEvents<eventno) {
      iResult=1; // indicate tree to be written
      HLTDebug("adding %d empty events to file %s", eventno-nofCurrentEvents, fName.Data());
      for (int i=nofCurrentEvents; i<eventno; i++) {
	fpTree->Fill();
      }
    }

    if (iResult>=0 && pSrcESD) {
      int nofObjects=fpEsd->GetList()->GetEntries();
      *fpEsd=*pSrcESD;
      if (nofObjects!=fpEsd->GetList()->GetEntries()) {
	// The source ESD contains object not present in the target ESD
	// before. Those objects will not be written to the tree since
	// the branch layout has been created earlier.
	// Create new tree with the additional branches, copy the entries
	// of the current tree into the new tree, and continue.
	TTree* pNewTree=new TTree(fTreeName, "Tree with HLT ESD objects");
	pNewTree->SetDirectory(0);
	AliESDEvent* readESD=new AliESDEvent;
	readESD->CreateStdContent();
	readESD->ReadFromTree(fpTree);
	fpEsd->Reset();
	fpEsd->WriteToTree(pNewTree);
	HLTDebug("cloning tree with %d entries", fpTree->GetEntries());
	for (int event=0; event<fpTree->GetEntries(); event++) {
	  fpTree->GetEntry(event);
	  *fpEsd=*readESD;
	  pNewTree->Fill();
	  fpEsd->Reset();
	}
	fpFile->Close();
	delete fpFile;
	delete readESD;
	delete fpTree;
	fpFile=new TFile(fName, "RECREATE");
	fpTree=pNewTree;
	*fpEsd=*pSrcESD;
	HLTDebug("new ESD with %d objects", fpEsd->GetList()->GetEntries());
      }
      fpTree->Fill();
      iResult=1; // indicate tree to be written
    }

    if (iResult>0) {
      fpFile->cd();
      fpTree->GetUserInfo()->Add(fpEsd);
      fpTree->Write(fpTree->GetName(),TObject::kOverwrite);
      fpTree->GetUserInfo()->Clear();
    }
  }

  return iResult;
}

void AliHLTEsdManagerImplementation::AliHLTEsdListEntry::SetDirectory(const char* directory)
{
  // see header file for class documentation
  if (!directory) return;
  if (!fName.IsNull()) {
    HLTWarning("ESD entry already in writing mode (%s), ignoring directory", fName.Data());
    return;
  }
  fDirectory=directory;
}

const char* AliHLTEsdManagerImplementation::AliHLTEsdListEntry::GetFileName() const
{
  // see header file for class documentation
  return fName.Data();
}

const char* AliHLTEsdManagerImplementation::AliHLTEsdListEntry::GetPrefix()
{
  // see header file for class documentation
  if (fPrefix.IsNull()) {
    fPrefix.Insert(0, fDt.fOrigin, kAliHLTComponentDataTypefOriginSize);
    fPrefix.Remove(TString::kTrailing, ' ');
    fPrefix.ToUpper();
    if (!fPrefix.Contains("HLT")) {
      fPrefix.Insert(0, "HLT");
    }
  }
  return fPrefix.Data();
}

int AliHLTEsdManagerImplementation::Merge(AliESDEvent* pTgt, AliESDEvent* pSrc) const
{
  // see header file for class documentation
  int iResult=0;
  if (!pTgt || !pSrc) return -EINVAL;

  TIter next(pSrc->GetList());
  TObject* pSrcObject=NULL;
  static int warningCount=0;
  while ((pSrcObject=next())) {
    if(!pSrcObject->InheritsFrom("TCollection")){
      // simple objects
      // for every type of object we have to find out whether it is empty or not
      // objects are only copied if non-empty, otherwise valid content would be
      // overridden by empty objects during the merging
      bool copy=false;
      TString name=pSrcObject->GetName();
      if(pSrcObject->InheritsFrom("AliHLTTriggerDecision")){
	copy=true;
      } else if (pSrcObject->IsA()==AliESDRun::Class()) {
	AliESDRun* pESDRun=dynamic_cast<AliESDRun*>(pSrcObject);
	// zero might be a valid run no in simulation, so we check also whether the CTP trigger classes are set
	copy=(pESDRun && (pESDRun->GetRunNumber()>0 || !pESDRun->GetActiveTriggerClasses().IsNull()));
      } else if (pSrcObject->IsA()==AliESDHeader::Class()) {
	AliESDHeader* pESDHeader=dynamic_cast<AliESDHeader*>(pSrcObject);
	copy=(pESDHeader && pESDHeader->GetTriggerMask()!=0);
      } else if (pSrcObject->IsA()==AliESDVertex::Class()) {
	AliESDVertex* pESDVertex=dynamic_cast<AliESDVertex*>(pSrcObject);
	copy=(pESDVertex && pESDVertex->GetNContributors()>0);
      } else if (pSrcObject->IsA()==AliESDTZERO::Class()) {
	AliESDTZERO* pESDTZero=dynamic_cast<AliESDTZERO*>(pSrcObject);
	copy=(pESDTZero && (pESDTZero->GetT0zVertex()!=0.0 || pESDTZero->GetT0()!=0.0));
      } else if (pSrcObject->IsA()==AliESDVZERO::Class()) {
	AliESDVZERO* pESDVZERO=dynamic_cast<AliESDVZERO*>(pSrcObject);
	// object contains data if one of the times exceeds the default value
	copy=pESDVZERO &&
	  (pESDVZERO->GetV0ATime() > kAliESDVZERODefaultTime+kAliESDVZERODefaultTimeGlitch ||
	   pESDVZERO->GetV0CTime() > kAliESDVZERODefaultTime+kAliESDVZERODefaultTimeGlitch);
      } else if (pSrcObject->IsA()==AliESDFMD::Class()) {
	AliESDFMD* pESDFMD=dynamic_cast<AliESDFMD*>(pSrcObject);
	copy=(pESDFMD && false); // have to find an easy valid condition
      } else if (pSrcObject->IsA()==AliESDZDC::Class()) {
	AliESDZDC* pESDZDC=dynamic_cast<AliESDZDC*>(pSrcObject);
	// object contains data if the EM energies are different from the default value
	copy=pESDZDC &&
	  (TMath::Abs(pESDZDC->GetZDCEMEnergy(0)-kAliESDZDCDefaultEMEnergy) > kAliESDZDCDefaultEMEnergyGlitch ||
	   TMath::Abs(pESDZDC->GetZDCEMEnergy(1)-kAliESDZDCDefaultEMEnergy) > kAliESDZDCDefaultEMEnergyGlitch);
      } else if (pSrcObject->IsA()==AliMultiplicity::Class()) {
	AliMultiplicity* pMultiplicity=dynamic_cast<AliMultiplicity*>(pSrcObject);
	copy=(pMultiplicity && pMultiplicity->GetNumberOfTracklets()>0);
      } else if (pSrcObject->IsA()==AliESDCaloTrigger::Class()) {
	AliESDCaloTrigger* pESDCaloTrigger=dynamic_cast<AliESDCaloTrigger*>(pSrcObject);
	copy=(pESDCaloTrigger && false); // have to find an easy valid condition
      } else if (pSrcObject->IsA()==AliESDCaloCells::Class()) {
	AliESDCaloCells* pESDCaloCells=dynamic_cast<AliESDCaloCells*>(pSrcObject);
	copy=(pESDCaloCells && false); // have to find an easy valid condition
      } else if (pSrcObject->IsA()==AliESDACORDE::Class()) {
	AliESDACORDE* pESDACORDE=dynamic_cast<AliESDACORDE*>(pSrcObject);
	copy=(pESDACORDE && false); // have to find an easy valid condition
      } else if (!AliHLTESDEventHelper::IsStdContent(name)) {
	// this is likely to be ok as long as it is not any object of the std content.
	copy=true;
      } else {
	HLTError("no merging implemented for object %s, omitting", name.Data());
      }
      if (copy) {
	//pSrcObject->Print();
	TObject* pTgtObject=pTgt->FindListObject(name);
	if (pTgtObject) {
	  pSrcObject->Copy(*pTgtObject);
	} else {
	  pTgt->AddObject(pSrcObject->Clone());
	}
      }
    } else if(pSrcObject->InheritsFrom("TClonesArray")){
      TClonesArray* pTClA=dynamic_cast<TClonesArray*>(pSrcObject);
      if (pTClA!=NULL && pTClA->GetEntriesFast()>0) {
	TString name=pTClA->GetName();
	TObject* pTgtObject=pTgt->GetList()->FindObject(name);
	TClonesArray* pTgtArray=NULL;
	if (pTgtObject!=NULL && pTgtObject->InheritsFrom("TClonesArray")){
	  pTgtArray=dynamic_cast<TClonesArray*>(pTgtObject);
	  if (pTgtArray) {
	    TString classType=pTClA->Class()->GetName();
	    if (classType.CompareTo(pTgtArray->Class()->GetName())==0) {
	      if (pTgtArray->GetEntries()==0) {
		pTgtArray->ExpandCreate(pTClA->GetEntries());
		for(int i=0; i<pTClA->GetEntriesFast(); ++i){
		  (*pTClA)[i]->Copy(*((*pTgtArray)[i]));
		}
	      } else {
		if (warningCount++<10) {
		  HLTWarning("TClonesArray \"%s\"  in target ESD %p is already filled with %d entries",
			     name.Data(), pTgt, pTgtArray->GetEntries());
		}
		iResult=-EBUSY;
	      }
	    } else {
	      if (warningCount++<10) {
		HLTWarning("TClonesArray \"%s\" exists in target ESD %p, but describes incompatible class type %s instead of %s",
			   name.Data(), pTgt, pTgtArray->GetClass()->GetName(), pTClA->GetClass()->GetName());
	      }
	      iResult=-EBUSY;
	    }
	  } else {
	    if (warningCount++<10) {
	      HLTError("internal error: dynamic cast failed for object %s %p", pTgtObject->GetName(), pTgtObject);
	    }
	    iResult=-EBUSY;
	  }
	} else if (pTgtObject) {
	  if (warningCount++<10) {
	    HLTWarning("object \"%s\" does already exist in target ESD %p, but is %s rather than TClonesArray"
		       " skipping data",
		       name.Data(), pTgt, pTgtObject->Class()->GetName());
	  }
	  iResult=-EBUSY;
	} else {
	  if (warningCount++<10) {
	    HLTWarning("object \"%s\" does not exist in target ESD %p, data can not be copied because it will be lost when filling the tree",
		       name.Data(), pTgt);
	  }
	  iResult=-ENOENT;
	}
      }
    }
  }
  return iResult;
}

bool AliHLTEsdManagerImplementation::AliHLTESDEventHelper::IsStdContent(const char* key)
{
  // check if the key denotes a std object
  TString needle=key;
  for (int i=0; i<kESDListN; i++) {
    if (needle.CompareTo(fgkESDListName[i])==0) return true;
  }
  return false;
}

int AliHLTEsdManagerImplementation::CheckClassConditions() const
{
  // this is a helper method which checks if some thing in the default
  // object initialization changes during the evolution of the source code
  // which is necessary for checking the validity of data in the object

  // check AliESDVZERO
  AliESDVZERO vzerodummy;
  if (TMath::Abs(vzerodummy.GetV0ATime()-kAliESDVZERODefaultTime) > kAliESDVZERODefaultTimeGlitch ||
      TMath::Abs(vzerodummy.GetV0CTime()-kAliESDVZERODefaultTime) > kAliESDVZERODefaultTimeGlitch) {
    HLTWarning("initialization of AliESDVZERO has changed, default time values now %f/%f, "
	       "revision of condition might be needed",
	       vzerodummy.GetV0ATime(), vzerodummy.GetV0CTime());
  }

  // check AliESDZDC
  AliESDZDC zdcdummy;
  if (TMath::Abs(zdcdummy.GetZDCEMEnergy(0)-kAliESDZDCDefaultEMEnergy) > kAliESDZDCDefaultEMEnergyGlitch ||
      TMath::Abs(zdcdummy.GetZDCEMEnergy(1)-kAliESDZDCDefaultEMEnergy) > kAliESDZDCDefaultEMEnergyGlitch) {
    HLTWarning("initialization of AliESDZDC has changed, default em energy values now %f/%f, "
	       "revision of condition might be needed",
	       zdcdummy.GetZDCEMEnergy(0), zdcdummy.GetZDCEMEnergy(1));
  }

  return 0;
}

TObject* AliHLTEsdManagerImplementation::CreateEsdEvent(bool bCreateStdContent) const
{
  // create ESDEvent object and optionally initialize standard content
  AliESDEvent* pESD=new AliESDEvent;
  if (pESD && bCreateStdContent) pESD->CreateStdContent();
  return pESD;
}

int AliHLTEsdManagerImplementation::DestroyEsdEvent(TObject* pESDInstance) const
{
  // destroy specified ESD object, pointer is invalid afterwords
  if (!pESDInstance) return -EINVAL;
  AliESDEvent* pESD=dynamic_cast<AliESDEvent*>(pESDInstance);
  if (!pESD) return -EINVAL;
  delete pESD;
  return 0;
}

int AliHLTEsdManagerImplementation::AddObject(TObject* pESDInstance, const TObject* pObject, const char* branchname) const
{
  // add object to the list of ESD contributors
  if (!pESDInstance || !pObject) return -EINVAL;
  AliESDEvent* pESD=dynamic_cast<AliESDEvent*>(pESDInstance);
  if (!pESD) return -EINVAL;
  TObject* pESDObject=pESD->FindListObject(branchname);
  if (pESDObject) {
    // copy the content to the already existing object
    pObject->Copy(*pESDObject);
  } else {
    // add a new object
    pESD->AddObject(pObject->Clone());
  }
  return 0;
}

int AliHLTEsdManagerImplementation::ResetEsdEvent(TObject* pESDInstance) const
{
  // reset the specified ESD object
  if (!pESDInstance) return -EINVAL;
  AliESDEvent* pESD=dynamic_cast<AliESDEvent*>(pESDInstance);
  if (!pESD) return -EINVAL;
  pESD->Reset();
  return 0;
}
