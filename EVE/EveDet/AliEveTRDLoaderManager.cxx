// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTRDLoaderManager.h"
#include "AliEveTRDLoader.h"
#include "AliEveTRDLoaderImp.h"

#include <TEveManager.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGListBox.h>
//#include <TGListTree.h>
//#include <TGString.h>
//#include <TGToolTip.h>
#include <TClonesArray.h>

#include "AliLog.h"

ClassImp(AliEveTRDLoaderManager)
ClassImp(AliEveTRDLoaderManagerEditor)

///////////////////////////////////////////////////////////
/////////        AliEveTRDLoaderManager      //////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDLoaderManager::AliEveTRDLoaderManager(const Text_t* n, const Text_t* t) :
  TEveElementList(n, t)
{
  // Constructor.
}

//______________________________________________________________________________
void AliEveTRDLoaderManager::Add(Int_t type, const Text_t *name, const Text_t *title)
{
  // Add something.

  AliEveTRDLoader *trdl = 0x0;
  switch(type){
  case AliEveTRDLoader::kTRDHits:
  case AliEveTRDLoader::kTRDDigits:
  case AliEveTRDLoader::kTRDClusters:
  case AliEveTRDLoader::kTRDTracklets:
    AddElement(trdl = new AliEveTRDLoader(name, title));
    break;
  case AliEveTRDLoader::kTRDRawRoot:
  case AliEveTRDLoader::kTRDRawDate:
    AddElement(trdl = new AliEveTRDLoaderRaw(name, title));
    break;
  default:
    AddElement(trdl = new AliEveTRDLoaderSim(name, title));
    break;
  }
  trdl->SetDataType(type);

  gEve->Redraw3D();
}


//______________________________________________________________________________
void AliEveTRDLoaderManager::Paint(Option_t *option)
{
  // Paint object.

  List_i ichmb = fChildren.begin();
  while(ichmb != fChildren.end()){
    (dynamic_cast<AliEveTRDLoader*>(*ichmb))->Paint(option);
    ichmb++;
  }
}

//______________________________________________________________________________
void AliEveTRDLoaderManager::Remove(Int_t entry)
{
  // Remove something.

  //printf("AliEveTRDLoaderManager::Remove(%d)\n", entry);
  List_i it = fChildren.begin();
  for(int i=0; i<entry; i++) it++;
  gEve->RemoveElement((*it), this);
  fChildren.erase(it);
}

///////////////////////////////////////////////////////////
//////////   AliEveTRDLoaderManagerEditor       ///////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDLoaderManagerEditor::
AliEveTRDLoaderManagerEditor(const TGWindow* p, Int_t width, Int_t height,
			     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0), fSelector(0), fAdd(0), fRemoveButton(0), fGroupFrame(0), fRemove(0)
{
  // Constructor.

  MakeTitle("AliEveTRDLoaderManager");

  // control frame - always there
  TGHorizontalFrame *fHorizontalFrame539 = new TGHorizontalFrame(this, 300, 26, kHorizontalFrame);

  TGLabel *fLabel546 = new TGLabel(fHorizontalFrame539,"Register Loader",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame);
  fLabel546->SetTextJustify(36);
  fHorizontalFrame539->AddFrame(fLabel546, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));

  // combo box
  fSelector = new TGComboBox(fHorizontalFrame539,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
  fSelector->AddEntry("MC (gAlice) ", AliEveTRDLoader::kTRDHits | AliEveTRDLoader::kTRDDigits | AliEveTRDLoader::kTRDClusters);
  fSelector->AddEntry("Hits ", AliEveTRDLoader::kTRDHits);
  fSelector->AddEntry("Digits ", AliEveTRDLoader::kTRDDigits);
  fSelector->AddEntry("Clusters ", AliEveTRDLoader::kTRDClusters);
  fSelector->AddEntry("Tracklets ", AliEveTRDLoader::kTRDTracklets);
  fSelector->AddEntry("Raw (ROOT) ", AliEveTRDLoader::kTRDRawRoot);
  fSelector->AddEntry("Raw (DATE) ", AliEveTRDLoader::kTRDRawDate);
  fSelector->Resize(136,22);
  fHorizontalFrame539->AddFrame(fSelector, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
  //fSelector->SetToolTipText("Select TRD data loader and add it to the list.\nThe loader can be removed by clicking the \"Remove\" button");
  fSelector->Connect("Selected(char*)", "AliEveTRDLoaderManagerEditor", this, "Add(char*)");
  AddFrame(fHorizontalFrame539, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

  fGroupFrame = 0;
  fRemove     = 0;
}


//______________________________________________________________________________
void AliEveTRDLoaderManagerEditor::Add(Char_t *name)
{
  // Slot to add something.
  if(!fGroupFrame){
    // "TRD Loaders" group frame
    fGroupFrame = new TGGroupFrame(this,"TRD Loaders",kVerticalFrame,TGGroupFrame::GetDefaultGC()(),TGGroupFrame::GetDefaultFontStruct());
    fGroupFrame->SetLayoutManager(new TGVerticalLayout(fGroupFrame));
    fGroupFrame->Resize(300,128);
    AddFrame(fGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2,2,2,2));

    fRemove = new TClonesArray("TGTextButton", 3);
  }


  // horizontal frame
  TGHorizontalFrame *fHorizontalFrame = new TGHorizontalFrame(fGroupFrame, 264, 26, kHorizontalFrame);

  TGLabel *fLabel717 = new TGLabel(fHorizontalFrame, name);
  fLabel717->SetTextJustify(36);
  fHorizontalFrame->AddFrame(fLabel717, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));


  Int_t nbutton = fM->fChildren.size();
  fRemoveButton = new((*fRemove)[nbutton]) TGTextButton(fHorizontalFrame, "Remove", nbutton);
  fRemoveButton->SetTextJustify(36);
  fRemoveButton->Resize(53,22);
  fRemoveButton->Connect("Clicked()", "AliEveTRDLoaderManagerEditor", this, Form("Remove(=%d)", nbutton));
  fRemoveButton->SetToolTipText(Form("Remove %s Loader", name));
  fHorizontalFrame->AddFrame(fRemoveButton, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));

  fGroupFrame->AddFrame(fHorizontalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));


  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  char *title[] = {"MC loader", "Single file loader", "Raw data loader"};
  // char *color[] = {"#ff0000", "#0000ff", "#59d454"};
  int id = fSelector->GetSelected(), type;
  switch(id){
  case AliEveTRDLoader::kTRDHits:
  case AliEveTRDLoader::kTRDDigits:
  case AliEveTRDLoader::kTRDClusters:
  case AliEveTRDLoader::kTRDTracklets:
    type = 1;
    break;
  case AliEveTRDLoader::kTRDRawRoot:
  case AliEveTRDLoader::kTRDRawDate:
    type = 2;
    break;
  default:
    type = 0;
    break;
  }
  fM->Add(id, name, title[type]);
}


//______________________________________________________________________________
void AliEveTRDLoaderManagerEditor::Remove(Int_t entry)
{
  // Slot to remove something.

  TIterator *it = fGroupFrame->GetList()->MakeIterator();
  int ientry = 0;
  while(/*TGFrame *f=(TGFrame*)*/it->Next()){
    //printf("%s\n", f->IsA()->GetName());
    if(entry == ientry){
      //fGroupFrame->RemoveFrame(f);
      break;
    }
    ientry++;
  }


  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  //fM->Remove(entry);
}

//______________________________________________________________________________
void AliEveTRDLoaderManagerEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveTRDLoaderManager*>(obj);
}

