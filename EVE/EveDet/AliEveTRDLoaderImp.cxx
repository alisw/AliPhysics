// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDLoaderImp.h"
#include "AliEveTRDModuleImp.h"
#include "EveBase/AliEveEventManager.h"

#include <TEveManager.h>

//#include "TFile.h"
#include "TTree.h"

#include <TGButton.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTRDrawData.h"
#include "AliTRDrawStream.h"
#include "AliTRDdigitsManager.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"

ClassImp(AliEveTRDLoaderSim)
ClassImp(AliEveTRDLoaderRaw)
ClassImp(AliEveTRDLoaderSimEditor)
//ClassImp(TRDLoaderRawEditor)

///////////////////////////////////////////////////////////
/////////////    AliEveTRDLoaderSim  /////////////////////
///////////////////////////////////////////////////////////


//______________________________________________________________________________
AliEveTRDLoaderSim::AliEveTRDLoaderSim(const Text_t* n, const Text_t* t) :
  AliEveTRDLoader(n, t)
  ,fRunLoader(0x0)
{
  // Constructor.
  if(gAlice && (fRunLoader = AliEveEventManager::AssertRunLoader())) SetDataLinked();
}

//______________________________________________________________________________
Bool_t	AliEveTRDLoaderSim::GoToEvent(int ev)
{
  // Go to given event.

  if(!fChildren.size()){
    AliWarning("Please select first the chamber that you want to monitor from \"Chamber(s) selector\".");
    return kFALSE;
  }
  if(!fDataType){
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
  if(fDataType&kTRDHits){
    fRunLoader->LoadHits("TRD", "READ");
    t = fRunLoader->GetTreeH("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadHits(t)) return kFALSE;
  }
  if(fDataType&kTRDDigits){
    fRunLoader->LoadDigits("TRD", "READ");
    t = fRunLoader->GetTreeD("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadDigits(t)) return kFALSE;
  }
  if(fDataType&kTRDClusters){
    fRunLoader->LoadRecPoints("TRD", "READ");
    t = fRunLoader->GetTreeR("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadClusters(t)) return kFALSE;
  }
  if(fDataType&kTRDTracklets){
    fRunLoader->LoadTracks("TRD", "READ");
    t = fRunLoader->GetTreeT("TRD", kFALSE);
    if(!t) return kFALSE;
    if(!LoadTracklets(t)) return kFALSE;
  }

  gEve->Redraw3D();
  return kTRUE;
}


//______________________________________________________________________________
Bool_t	AliEveTRDLoaderSim::Open(const char *filename, const char *dir)
{
  // Open file in given dir.

  if(fRunLoader) return kTRUE;
  
  fRunLoader = AliRunLoader::Instance();
  if(!fRunLoader) fRunLoader = AliRunLoader::Open(filename,
         AliConfig::GetDefaultEventFolderName(),"read");
  if(!fRunLoader) return kFALSE;

  gAlice = fRunLoader->GetAliRun();
  if(!gAlice) fRunLoader->LoadgAlice();
  if(!gAlice) return kFALSE;
 
  fFilename = filename;
  fDir = dir;
  fDir += "/";
  fRunLoader->SetDirName(fDir);

  SetDataLinked();
  return kTRUE;
}


///////////////////////////////////////////////////////////
/////////////   AliEveTRDLoaderRaw    /////////////////////
///////////////////////////////////////////////////////////


//______________________________________________________________________________
AliEveTRDLoaderRaw::AliEveTRDLoaderRaw(const Text_t* n, const Text_t* t) :
  AliEveTRDLoader(n, t)
  ,fRawDateReader (0x0)
  ,fRawRootReader (0x0)
  ,fRaw           (0x0)
  ,fEventCnt(-1)
{
  // Constructor.
}

//______________________________________________________________________________
Bool_t  AliEveTRDLoaderRaw::Open(const char *filename, const char *dir)
{
  // Open file in gvenn dir.

  fFilename = filename;
  fDir = dir;
  fDir += "/";

  if(fRaw) delete fRaw;
  fRaw = new AliTRDrawData();

  if(fDataType&kTRDRawRoot){
    if(fRawRootReader) delete fRawRootReader;
    fRawRootReader = new AliRawReaderRoot(filename);
  } else if(fDataType&kTRDRawDate){
    if(fRawDateReader) delete fRawDateReader;
    fRawDateReader = new AliRawReaderDate(fDir+fFilename);
  } else {
    AliError("No data type was set.");
    return kFALSE;
  }
  SetDataLinked();
  return kTRUE;
}


//______________________________________________________________________________
Bool_t AliEveTRDLoaderRaw::GoToEvent(int ev)
{
  // Go to given event.

  //AliInfo(Form("Event %d %d %d", ev, fEvent, fEventCnt));


  if(!fChildren.size()){
    AliWarning("Please select first the chamber that you want to monitor from \"Chamber(s) selector\".");
    return kFALSE;
  }

  static const TEveException kEH("AliEveTRDLoader::GotoEvent ");
  if(fRawRootReader == 0x0) throw(kEH + "data file not opened.");

  fEvent = ev;
  if(ev == fEventCnt) return kTRUE;
  if(ev < fEventCnt) {
    fRawRootReader->RewindEvents();
    fEventCnt = -1;
  }

  Bool_t FOUND = kFALSE;
  while(fRawRootReader->NextEvent()){ 
    fEventCnt++;
    if(fEventCnt == ev){
      FOUND = kTRUE;
      break;
    }  
  }
  if(!FOUND) return kFALSE;
  
  LoadEvent();
  gEve->Redraw3D();

  return kTRUE;
}


//______________________________________________________________________________
Bool_t AliEveTRDLoaderRaw::LoadEvent()
{
  // Load event.
  AliInfo("Loading ...");

  static const TEveException kEH("AliEveTRDLoader::LoadEvent ");
  if(fRawRootReader == 0x0) throw(kEH + "data file not opened.");


  fRawRootReader->Reset();
  fRawRootReader->SelectEquipment(0, 1024, 1041);
  fRawRootReader->Select("TRD");
  
//   AliTRDrawStream::AllowCorruptedData();
//   AliTRDrawStream::DisableStackNumberChecker();
//   AliTRDrawStream::DisableStackLinkNumberChecker();

  AliTRDrawStream *pinput = new AliTRDrawStream(fRawRootReader);
  AliTRDrawStream &input = *pinput;

 // AliInfo(Form("Stream version: %s", input.IsA()->GetName()));

  AliEveTRDChamber *chmb;
  AliTRDdigitsManager *dm = new AliTRDdigitsManager();
  dm->CreateArrays();

  Int_t det    = 0;
  while ((det = input.NextChamber(dm)) >= 0){
    if(!(chmb=GetChamber(det))) continue;
    chmb->LoadDigits(dm);

    dm->RemoveDigits(det);
    dm->RemoveDictionaries(det);
    dm->ClearIndexes(det);
  }


  delete pinput;
  pinput = NULL;


  return kTRUE;
}



///////////////////////////////////////////////////////////
//////////// AliEveTRDLoaderSimEditor /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDLoaderSimEditor::AliEveTRDLoaderSimEditor(const TGWindow* p, Int_t width, Int_t height,
                                                   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back)
  ,fM(0x0)
  ,fCheckedHits(0x0)
  ,fCheckedDigits(0x0)
  ,fCheckedClusters(0x0)
  ,fCheckedTracklets(0x0)
{
  // Constructor.

  MakeTitle("AliEveTRDLoaderSim");

  // "Data selector" group frame
  TGGroupFrame *fGroupFrame = new TGGroupFrame(this,"Data selector");
  fCheckedHits = new TGCheckButton(fGroupFrame,"  Hits");
  fCheckedHits->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, Form("Toggle(=%d)", (Int_t)AliEveTRDLoader::kTRDHits));
  fGroupFrame->AddFrame(fCheckedHits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fCheckedDigits = new TGCheckButton(fGroupFrame,"  Digits");
  fCheckedDigits->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, Form("Toggle(=%d)", (Int_t)AliEveTRDLoader::kTRDDigits));
  fGroupFrame->AddFrame(fCheckedDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fCheckedClusters = new TGCheckButton(fGroupFrame,"  Clusters");
  fCheckedClusters->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, Form("Toggle(=%d)", (Int_t)AliEveTRDLoader::kTRDClusters));
  fGroupFrame->AddFrame(fCheckedClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fCheckedTracklets = new TGCheckButton(fGroupFrame,"  Tracklets ");
  fCheckedTracklets->Connect("Clicked()", "AliEveTRDLoaderSimEditor", this, Form("Toggle(=%d)", (Int_t)AliEveTRDLoader::kTRDTracklets));
  fGroupFrame->AddFrame(fCheckedTracklets, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

  fGroupFrame->SetLayoutManager(new TGVerticalLayout(fGroupFrame));
  //	fGroupFrame->Resize(164,116);
  AddFrame(fGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
}

//______________________________________________________________________________
void AliEveTRDLoaderSimEditor::SetModel(TObject* obj)
{
  // Set model object.

  if(!(fM = dynamic_cast<AliEveTRDLoaderSim*>(obj))) return;

  Bool_t kRL   = (fM->IsDataLinked()) ? kTRUE : kFALSE;

  fCheckedHits->SetEnabled(kRL);
  if(kRL) fCheckedHits->SetState(fM->fDataType&AliEveTRDLoader::kTRDHits ? kButtonDown : kButtonUp);
  fCheckedDigits->SetEnabled(kRL);
  if(kRL) fCheckedDigits->SetState(fM->fDataType&AliEveTRDLoader::kTRDDigits ? kButtonDown : kButtonUp);
  fCheckedClusters->SetEnabled(kRL);
  if(kRL) fCheckedClusters->SetState(fM->fDataType&AliEveTRDLoader::kTRDClusters ? kButtonDown : kButtonUp);
  fCheckedTracklets->SetEnabled(kRL);
  if(kRL) fCheckedTracklets->SetState(fM->fDataType&AliEveTRDLoader::kTRDTracklets ? kButtonDown : kButtonUp);
}

//______________________________________________________________________________
void AliEveTRDLoaderSimEditor::Toggle(Int_t id)
{
  // Toggle given button id.

  switch(id){
  case AliEveTRDLoader::kTRDHits:
    fM->fDataType |= fCheckedHits->IsDown() ? AliEveTRDLoader::kTRDHits : 0;
    break;
  case AliEveTRDLoader::kTRDDigits:
    fM->fDataType |= fCheckedDigits->IsDown() ? AliEveTRDLoader::kTRDDigits : 0;
    break;
  case AliEveTRDLoader::kTRDClusters:
    fM->fDataType |= fCheckedClusters->IsDown() ? AliEveTRDLoader::kTRDClusters : 0;
    break;
  case AliEveTRDLoader::kTRDTracklets:
    fM->fDataType |= fCheckedTracklets->IsDown() ? AliEveTRDLoader::kTRDTracklets : 0;
    break;
  }
}

