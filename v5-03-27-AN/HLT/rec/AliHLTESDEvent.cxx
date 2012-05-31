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

/// @file   AliHLTESDEvent.cxx
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  A streamlined container class for AliESDEvent.
/// @note   

#include "AliHLTESDEvent.h"
#include "AliHLTESDtrack.h"
#include "AliHLTOnlineESDtrack.h"
#include "AliHLTLogging.h"
#include "TList.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TClass.h"
#include <cerrno>
#include <memory>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTESDEvent)

AliHLTESDEvent::AliHLTESDEvent()
  : AliESDEvent()
  , fTemplateEsd(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTESDEvent::AliHLTESDEvent(const AliHLTESDEvent& src)
  : AliESDEvent(src)
  , fTemplateEsd(NULL)
{
  // copy contructor
  if (src.fTemplateEsd) {
    fTemplateEsd=new AliESDEvent;
    if (fTemplateEsd) *fTemplateEsd=*src.fTemplateEsd;
  }
}

AliHLTESDEvent& AliHLTESDEvent::operator=(const AliHLTESDEvent& src)
{
  // assignment operator
  if (this==&src) return *this;

  AliESDEvent::operator=(src);
  if (src.fTemplateEsd) {
    fTemplateEsd=new AliESDEvent;
    if (fTemplateEsd) *fTemplateEsd=*src.fTemplateEsd;
  }

  return *this;
}

AliHLTESDEvent::~AliHLTESDEvent()
{
  // destructor
  if (fTemplateEsd) delete fTemplateEsd;
  fTemplateEsd=NULL;
}

void AliHLTESDEvent::Print(const char* options) const
{
  /// overloaded from TObject, print info
  AliESDEvent::Print(options);
}

void AliHLTESDEvent::Dump() const
{
  /// overloaded from TObject, more crude data dump
  AliESDEvent::Dump();
}

void AliHLTESDEvent::Clear(Option_t * option)
{
  /// overloaded from TObject, clear object
  
  AliESDEvent::Clear(option);
}

TObject * AliHLTESDEvent::Clone(const char *newname) const
{
  /// overloaded from TObject, clone object

  return AliESDEvent::Clone(newname);
}

void AliHLTESDEvent::Copy(TObject &object) const
{
  /// overloaded from TObject, copy object

  AliESDEvent::Copy(object);
}

int AliHLTESDEvent::LoadTemplate(const char* /*cdbPath*/)
{
  /// load a template from OCDB or create the default template

  if (fTemplateEsd) delete fTemplateEsd;
  fTemplateEsd=NULL;

  AliHLTLogging log;

  // default list of skiped ESD objects
  TString skipObjects=
    // "AliESDRun,"
    // "AliESDHeader,"
    // "AliESDZDC,"
    "AliESDFMD,"
    // "AliESDVZERO,"
    // "AliESDTZERO,"
    // "TPCVertex,"
    // "SPDVertex,"
    // "PrimaryVertex,"
    // "AliMultiplicity,"
    // "PHOSTrigger,"
    // "EMCALTrigger,"
    // "SPDPileupVertices,"
    // "TrkPileupVertices,"
    "Tracks,"
    "MuonTracks,"
    "PmdTracks,"
    "TrdTracks,"
    "Cascades,"
    "Kinks,"
    "AliRawDataErrorLogs,"
    "AliESDACORDE";

  std::auto_ptr<AliESDEvent> ptrEsd(new AliESDEvent);
  if (!ptrEsd.get()) return -ENOMEM;
  
  ptrEsd->CreateStdContent();

  // remove some of the objects which are not needed
  if (ptrEsd->GetList() && !skipObjects.IsNull()) {
    TObjArray* pTokens=skipObjects.Tokenize(",");
    if (pTokens) {
      const char* id=NULL;
      TIter next(pTokens);
      TObject* pObject=NULL;
      while ((pObject=next())!=NULL) {
	id=pObject->GetName();
	TObject* pObj=ptrEsd->GetList()->FindObject(id);
	if (pObj) {
	  log.LoggingVarargs(kHLTLogInfo, "AliHLTESDEvent", "LoadTemplate" , __FILE__ , __LINE__ , "removing object %s", id);
	  ptrEsd->GetList()->Remove(pObj);
	  delete pObj;
	} else {
	  log.LoggingVarargs(kHLTLogWarning, "AliHLTESDEvent", "LoadTemplate" , __FILE__ , __LINE__ , "failed to remove object '%s' from ESD", id);
	}
      }
      ptrEsd->GetStdContent();
      delete pTokens;
    }
  }

  // add the new objects
  if (ptrEsd->GetList()->FindObject("Tracks")==NULL) {
    TClonesArray* pTracks=new TClonesArray("AliHLTOnlineESDtrack",0);
    pTracks->SetName("Tracks");
    ptrEsd->AddObject(pTracks);
  } else {
    log.LoggingVarargs(kHLTLogWarning, "AliHLTESDEvent", "LoadTemplate" , __FILE__ , __LINE__ , "member 'Tracks' is still in the list, skipping to add customized array");
  }

  fTemplateEsd = ptrEsd.release();

  return 0;
}

void AliHLTESDEvent::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliHLTESDEvent.

  AliHLTLogging log;
  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(AliHLTESDEvent::Class(),this);
    TObject* srcobject=GetList()->FindObject("Tracks");
    if (srcobject) {
      TClonesArray* tclsrc=dynamic_cast<TClonesArray*>(srcobject);
      if (tclsrc && tclsrc->GetClass() && tclsrc->GetClass()==AliHLTOnlineESDtrack::Class()) {
	GetList()->Remove(tclsrc);
	TClonesArray* tcltgt=new TClonesArray("AliHLTESDtrack",0);
	if (tcltgt) {
	  tcltgt->SetName("Tracks");
	  tcltgt->ExpandCreate(tclsrc->GetEntriesFast());
	  for (int i=0; i<tclsrc->GetEntriesFast(); i++) {
	    AliHLTOnlineESDtrack* pHLTTrack=dynamic_cast<AliHLTOnlineESDtrack*>(tclsrc->At(i));
	    if (pHLTTrack==NULL || (*tcltgt)[i]==NULL) continue;
	    AliHLTESDtrack* pESDTrack=dynamic_cast<AliHLTESDtrack*>((*tcltgt)[i]);
	    if (pESDTrack==NULL) continue;
	    *pESDTrack=*pHLTTrack;
	  }
	  AddObject(tcltgt);
	  tclsrc->Clear("C");
	  delete tclsrc;
	  tclsrc=NULL;
	  srcobject=NULL;
	}
      }
    }
  } else {
    if (fTemplateEsd) {
      TObject* esdTracksObject=NULL;
      TObject* onlineTracksObject=NULL;
      fTemplateEsd->Reset();
      TIter nextobject(this->GetList());
      TObject* esdobject=NULL;
      while ((esdobject=nextobject())) {
	TObject* tgtobject=fTemplateEsd->GetList()->FindObject(esdobject->GetName());
	if (!tgtobject) {
	  log.LoggingVarargs(kHLTLogInfo, "AliHLTESDEvent", "Streamer" , __FILE__ , __LINE__ , "skipping object %s", esdobject->GetName());
	  continue;
	}

	if (strcmp(tgtobject->GetName(), "Tracks")==0 && tgtobject->IsA()==TClonesArray::Class()) {
	  // special handling of Tracks array, replace Tracks array with own array
	  TClonesArray* tclsrc=dynamic_cast<TClonesArray*>(esdobject);
	  TClonesArray* tcltgt=dynamic_cast<TClonesArray*>(tgtobject);
	  if (!tclsrc || !tcltgt) continue;
	  tcltgt->Clear("C");

	  if (tcltgt->GetClass()==AliHLTOnlineESDtrack::Class()) {
	    tcltgt->ExpandCreate(tclsrc->GetEntriesFast());
	    for (int i=0; i<tclsrc->GetEntriesFast(); i++) {
	      if (!tclsrc->At(i)) continue;
	      AliESDtrack* pESDTrack=dynamic_cast<AliESDtrack*>(tclsrc->At(i));
	      TObject* tgtobj=(*tcltgt)[i];
	      if (pESDTrack==NULL || tgtobj==NULL) continue;
	      AliHLTOnlineESDtrack* pHLTTrack=dynamic_cast<AliHLTOnlineESDtrack*>(tgtobj);
	      if (!pHLTTrack) continue;
	      *pHLTTrack=*pESDTrack;
	    }
	    esdTracksObject=esdobject;
	    onlineTracksObject=tcltgt;
	    continue;
	  } else if (tcltgt->GetClass()!=AliESDtrack::Class()) {
	    // no handling if not Ali(HLT)ESDtrack
	    continue;
	  }
	}

	// default: just copy the object
	esdobject->Copy(*tgtobject);
      }

      // replace with the optimized objects
      if (esdTracksObject) GetList()->Remove(esdTracksObject);
      if (onlineTracksObject) GetList()->Add(onlineTracksObject);
      
      R__b.WriteClassBuffer(AliHLTESDEvent::Class(),fTemplateEsd);

      // replace with the original objects
      if (onlineTracksObject) GetList()->Remove(onlineTracksObject);
      if (esdTracksObject) GetList()->Add(esdTracksObject);
      esdTracksObject=NULL; onlineTracksObject=NULL;
    } else {
      R__b.WriteClassBuffer(AliHLTESDEvent::Class(),this);
    }
  }
}

void AliHLTESDEvent::Execute(const char *method,  const char *params, Int_t *error)
{
  // handle custom function calls
  if (method && strcmp(method, "LoadTemplate")==0) {
    LoadTemplate(params);
    return;
  }
  AliESDEvent::Execute(method, params, error);
}
