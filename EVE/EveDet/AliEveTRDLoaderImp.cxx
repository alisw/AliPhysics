// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDLoaderImp.h"
#include "AliEveTRDModuleImp.h"

#include <TEveManager.h>

#include "TFile.h"
#include "TTree.h"

#include <TGButton.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTRDrawData.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"

#include "AliTRDv1.h"
#include "AliTRDhit.h"
#include "AliTRDdigitsManager.h"

ClassImp(AliEveTRDLoaderSim)
ClassImp(AliEveTRDLoaderRaw)
ClassImp(AliEveTRDLoaderSimEditor)
//ClassImp(TRDLoaderRawEditor)

///////////////////////////////////////////////////////////
/////////////    AliEveTRDLoaderSim  /////////////////////
///////////////////////////////////////////////////////////


//______________________________________________________________________________
AliEveTRDLoaderSim::AliEveTRDLoaderSim(const Text_t* n, const Text_t* t) :
  AliEveTRDLoader(n, t),
  fRunLoader(0)
{}

//______________________________________________________________________________
AliEveTRDLoaderSim::~AliEveTRDLoaderSim()
{}

//______________________________________________________________________________
Bool_t	AliEveTRDLoaderSim::GoToEvent(int ev)
{
  if(!fChildren.size()){
    AliWarning("Please select first the chamber that you want to monitor from \"Chamber(s) selector\".");
    return kFALSE;
  }
  if(!fLoadHits && !fLoadDigits && !fLoadClusters && !fLoadTracks){
    AliWarning("Please select first the type of data that you want to monitor and then hit the \"Load\" button.");
    return kFALSE;
  }

  fEvent = ev;

  if(!fRunLoader){
    AliError("RunLoader not initialized.");
    return kFALSE;
  }
  fRunLoader->UnloadAll("TRD");
  Unload();

  if(fRunLoader->GetEvent(ev)) return kFALSE;
  TTree *t = 0;
  if(fLoadHits){
    fRunLoader->LoadHits("TRD", "READ");
    t = fRunLoader->GetTreeH("TRD", kFALSE);
    if(!t) return kFALSE;
    fTRD->SetTreeAddress();
    if(!LoadHits(t)) return kFALSE;
  }
  if(fLoadDigits){
    fRunLoader->LoadDigits("TRD", "READ");
    t = fRunLoader->GetTreeD("TRD", kFALSE);
    if(!t) return kFALSE;
    fTRD->SetTreeAddress();
    if(!LoadDigits(t)) return kFALSE;
  }
  if(fLoadClusters){
    fRunLoader->LoadRecPoints("TRD", "READ");
    t = fRunLoader->GetTreeR("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadClusters(t)) return kFALSE;
  }
  if(fLoadTracks){
    fRunLoader->LoadTracks("TRD", "READ");
    t = fRunLoader->GetTreeT("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadTracklets(t)) return kFALSE;
  }

  gEve->Redraw3D();
  return kTRUE;
}


//______________________________________________________________________________
Bool_t	AliEveTRDLoaderSim::LoadHits(TTree *tH)
{
  Info("LoadHits()", "Loading ...");
  if(!fChildren.size()) return kTRUE;

  AliEveTRDChamber *chmb = 0x0;
  AliTRDhit *hit = 0x0;
  Int_t d;
  for(int iTrack=0; iTrack<tH->GetEntries(); iTrack++){
    gAlice->ResetHits();
    if(!tH->GetEvent(iTrack)) continue;
    hit = (AliTRDhit*)fTRD->FirstHit(-1);
    if(!hit) continue;
    d = hit->GetDetector();
    chmb = GetChamber(d);
    while(hit){
      if(d != hit->GetDetector()){
        d = hit->GetDetector();
        chmb = GetChamber(d);
      }
      if(chmb) chmb->AddHit(hit);
      hit = (AliTRDhit*)fTRD->NextHit();
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t	AliEveTRDLoaderSim::Open(const char *filename, const char *dir)
{
  //Info("Open()", "");


  fFilename = filename;
  fDir = dir;
  fDir += "/";

  fRunLoader = AliRunLoader::GetRunLoader();
  if(!fRunLoader) fRunLoader = AliRunLoader::Open(filename,
                                                  AliConfig::GetDefaultEventFolderName(),"read");
  if(!fRunLoader){
    AliError("Couldn't find run loader");
    return kFALSE;
  }
  fRunLoader->SetDirName(fDir);

  gAlice = fRunLoader->GetAliRun();
  if(!gAlice) fRunLoader->LoadgAlice();
  if(!gAlice){
    AliError("Couldn't find gAlice object");
    return kFALSE;
  }
  fTRD = (AliTRDv1*)gAlice->GetDetector("TRD");
  if(!fTRD){
    AliError("Couldn't find TRD");
    return kFALSE;
  }

  return kTRUE;
}



///////////////////////////////////////////////////////////
/////////////   AliEveTRDLoaderRaw    /////////////////////
///////////////////////////////////////////////////////////


//______________________________________________________________________________
AliEveTRDLoaderRaw::AliEveTRDLoaderRaw(const Text_t* n, const Text_t* t) :
  AliEveTRDLoader(n, t),
  fRawDateReader (0),
  fRawRootReader (0),
  fRaw           (0),
  fDataRoot      (kTRUE),
  fEventOld      (-1)
{}

//______________________________________________________________________________
AliEveTRDLoaderRaw::~AliEveTRDLoaderRaw()
{}


//______________________________________________________________________________
Bool_t  AliEveTRDLoaderRaw::Open(const char *filename, const char *dir)
{
  //	Info("Open()", Form("Open %s/%s", dir, filename));
  fFilename = filename;
  fDir = dir;
  fDir += "/";


  if(fRaw) delete fRaw;
  fRaw = new AliTRDrawData();

  if(fDataRoot){
    if(fRawRootReader) delete fRawRootReader;
    fRawRootReader = new AliRawReaderRoot(filename);
  } else {
    if(fRawDateReader) delete fRawDateReader;
    fRawDateReader = new AliRawReaderDate(fDir+fFilename);
  }

  return kTRUE;
}

//______________________________________________________________________________
void AliEveTRDLoaderRaw::SetDataType(TRDDataTypes type)
{
  fDataRoot = (type == kRawRoot) ? kTRUE : kFALSE;
}

//______________________________________________________________________________
Bool_t AliEveTRDLoaderRaw::GoToEvent(int ev)
{
  if(!fChildren.size()){
    AliWarning("Please select first the chamber that you want to monitor from \"Chamber(s) selector\".");
    return kFALSE;
  }

  static const TEveException eH("AliEveTRDLoader::GotoEvent ");
  if(fRawRootReader == 0x0) throw(eH + "data file not opened.");


  if(ev == fEventOld) return kTRUE;
  Bool_t checkEnd;
  if(ev < fEventOld) {
    fRawRootReader->RewindEvents();
    fEventOld = -1;
    checkEnd = kFALSE;
  } else checkEnd = kTRUE;

  do NextEvent(); while(fEventOld != ev && !(checkEnd == kTRUE && fEventOld == 0));
  LoadEvent();
  gEve->Redraw3D();
  //gEve->EnableRedraw();
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliEveTRDLoaderRaw::LoadEvent()
{
  Info("LoadEvent()", "Loading ...");

  static const TEveException eH("AliEveTRDLoader::LoadEvent ");
  if(fRawRootReader == 0x0) throw(eH + "data file not opened.");


  fRawRootReader->Reset();

  AliEveTRDChamber *chmb;
  AliTRDdigitsManager *dm;
  dm = fRaw->Raw2Digits(fRawRootReader);

  for(int idet=0; idet<540; idet++){
    if(!(chmb=GetChamber(idet))) continue;
    chmb->LoadDigits(dm);
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliEveTRDLoaderRaw::NextEvent(Bool_t rewindOnEnd)
{
  static const TEveException eH("AliEveTRDLoader::NextEvent ");
  if(fRawRootReader == 0x0) throw(eH + "data file not opened.");


  if(fRawRootReader->NextEvent() == kTRUE) ++fEventOld;
  else {
    if(fEventOld == -1) throw(eH + "no events available.");
    if(rewindOnEnd) {
      Warning("NextEvent()", Form("Reached end of stream (event=%d), rewinding to first event.", fEventOld));
      fRawRootReader->RewindEvents();
      fRawRootReader->NextEvent();
      fEventOld = 0;
    } else throw(eH + "last event reached.");
  }
}



///////////////////////////////////////////////////////////
//////////// AliEveTRDLoaderSimEditor /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDLoaderSimEditor::AliEveTRDLoaderSimEditor(const TGWindow* p, Int_t width, Int_t height,
                                                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0), fLoadHits(0), fLoadDigits(0), fLoadClusters(0), fLoadTracks(0)
{
  MakeTitle("AliEveTRDLoaderSim");

  // "Data selector" group frame
  TGGroupFrame *fGroupFrame = new TGGroupFrame(this,"Data selector");
  fLoadHits = new TGCheckButton(fGroupFrame,"  Hits");
  fLoadHits->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, "Toggle(=0)");
  fGroupFrame->AddFrame(fLoadHits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fLoadDigits = new TGCheckButton(fGroupFrame,"  Digits");
  fLoadDigits->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, "Toggle(=1)");
  fGroupFrame->AddFrame(fLoadDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fLoadClusters = new TGCheckButton(fGroupFrame,"  Clusters");
  fLoadClusters->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, "Toggle(=2)");
  fGroupFrame->AddFrame(fLoadClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fLoadTracks = new TGCheckButton(fGroupFrame,"  Tracklets ");
  fLoadTracks->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, "Toggle(=3)");
  fGroupFrame->AddFrame(fLoadTracks, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fGroupFrame->SetLayoutManager(new TGVerticalLayout(fGroupFrame));
  //	fGroupFrame->Resize(164,116);
  AddFrame(fGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
}

//______________________________________________________________________________
AliEveTRDLoaderSimEditor::~AliEveTRDLoaderSimEditor()
{}

//______________________________________________________________________________
void AliEveTRDLoaderSimEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveTRDLoaderSim*>(obj);

  Bool_t kFile = kTRUE;
  if(fM->fFilename.CompareTo("") == 0) kFile = kFALSE;

  /*	printf("\thits      %s\n", fM->fLoadHits ? "true" : "false");
    printf("\tdigits    %s\n", fM->fLoadDigits ? "true" : "false");
    printf("\tclusters  %s\n", fM->fLoadClusters ? "true" : "false");
    printf("\ttracklets %s\n", fM->fLoadTracks ? "true" : "false");*/
  fLoadHits->SetEnabled(kFile);
  if(kFile) fLoadHits->SetState(fM->fLoadHits ? kButtonDown : kButtonUp);
  fLoadDigits->SetEnabled(kFile);
  if(kFile) fLoadDigits->SetState(fM->fLoadDigits ? kButtonDown : kButtonUp);
  fLoadClusters->SetEnabled(kFile);
  if(kFile) fLoadClusters->SetState(fM->fLoadClusters ? kButtonDown : kButtonUp);
  fLoadTracks->SetEnabled(kFile);
  if(kFile) fLoadTracks->SetState(fM->fLoadTracks ? kButtonDown : kButtonUp);
}

//______________________________________________________________________________
void AliEveTRDLoaderSimEditor::Toggle(Int_t id)
{
  switch(id){
    case 0:
      fM->fLoadHits = fLoadHits->IsDown() ? kTRUE : kFALSE;
      break;
    case 1:
      fM->fLoadDigits = fLoadDigits->IsDown() ? kTRUE : kFALSE;
      break;
    case 2:
      fM->fLoadClusters = fLoadClusters->IsDown() ? kTRUE : kFALSE;
      break;
    case 3:
      fM->fLoadTracks = fLoadTracks->IsDown() ? kTRUE : kFALSE;
      break;
  }
}

///////////////////////////////////////////////////////////
/////////////   TRDLoaderRawEditor    /////////////////////
///////////////////////////////////////////////////////////

// //________________________________________________________
// TRDLoaderRawEditor::TRDLoaderRawEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options | kVerticalFrame, back)
// {
// 	MakeTitle("AliEveTRDLoaderRaw");
// }
//
// void	TRDLoaderRawEditor::SetModel(TObject* obj)
// {
// 	Info("SetModel()", "");
// 	fM = dynamic_cast<AliEveTRDLoaderRaw*>(obj);
// }
