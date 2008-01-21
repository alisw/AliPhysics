// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector2DEditor.h"
#include <EveDet/AliEveTPCSector2D.h>

#include <TGButton.h>
#include <TGComboBox.h>
#include <TGLabel.h>


//______________________________________________________________________________
// AliEveTPCSector2DEditor
//

ClassImp(AliEveTPCSector2DEditor)

AliEveTPCSector2DEditor::AliEveTPCSector2DEditor(const TGWindow *p,
				     Int_t width, Int_t height,
				     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fShowMax(0), fAverage(0), fUseTexture(0), fPickEmpty(0), fPickMode(0)
{
  MakeTitle("AliEveTPCSector2D");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fShowMax = new TGCheckButton(f, "ShowMax");
    f->AddFrame(fShowMax, new TGLayoutHints(kLHintsLeft, 3, 16, 1, 0));
    fShowMax->Connect("Toggled(Bool_t)","AliEveTPCSector2DEditor", this, "DoShowMax()");
    fAverage = new TGCheckButton(f, "Average");
    f->AddFrame(fAverage, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 0));
    fAverage->Connect("Toggled(Bool_t)","AliEveTPCSector2DEditor", this, "DoAverage()");
    AddFrame(f);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fUseTexture = new TGCheckButton(f, "UseTexture");
    f->AddFrame(fUseTexture, new TGLayoutHints(kLHintsTop, 3, 9, 1, 0));
    fUseTexture->Connect("Toggled(Bool_t)","AliEveTPCSector2DEditor", this, "DoUseTexture()");
    fPickEmpty = new TGCheckButton(f, "PickEmpty");
    f->AddFrame(fPickEmpty, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fPickEmpty->Connect("Toggled(Bool_t)","AliEveTPCSector2DEditor", this, "DoPickEmpty()");
    AddFrame(f);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "PickMode");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 10, 1, 2));
    fPickMode = new TGComboBox(f);
    fPickMode->AddEntry("Print", 0);
    fPickMode->AddEntry("1D histo", 1);
    fPickMode->AddEntry("2D histo", 2);
    TGListBox* lb = fPickMode->GetListBox();
    lb->Resize(lb->GetWidth(), 3*18);
    fPickMode->Resize(80, 20);
    fPickMode->Connect("Selected(Int_t)", "AliEveTPCSector2DEditor", this, "DoPickMode(Int_t)");
    f->AddFrame(fPickMode, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
    AddFrame(f);
  }
}

AliEveTPCSector2DEditor::~AliEveTPCSector2DEditor()
{}

/******************************************************************************/

void AliEveTPCSector2DEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveTPCSector2D*>(obj);

  fShowMax->SetState(fM->fShowMax ? kButtonDown : kButtonUp);
  SetupAverage();

  fUseTexture->SetState(fM->fUseTexture ? kButtonDown : kButtonUp);
  fPickEmpty->SetState(fM->fPickEmpty ? kButtonDown : kButtonUp);
  fPickMode->Select(fM->fPickMode, kFALSE);
}

/******************************************************************************/

void AliEveTPCSector2DEditor::DoShowMax()
{
  fM->SetShowMax(fShowMax->IsOn());
  SetupAverage();
  Update();
}

void AliEveTPCSector2DEditor::DoAverage()
{
  fM->SetAverage(fAverage->IsOn());
  Update();
}

void AliEveTPCSector2DEditor::SetupAverage()
{
  if(fM->fShowMax) {
    fAverage->SetEnabled(kFALSE);
  } else {
    fAverage->SetEnabled(kTRUE);
    fAverage->SetState(fM->fAverage ? kButtonDown : kButtonUp);
  }
}

/******************************************************************************/

void AliEveTPCSector2DEditor::DoUseTexture()
{
  fM->fUseTexture = fUseTexture->IsOn();
  Update();
}

void AliEveTPCSector2DEditor::DoPickEmpty()
{
  fM->fPickEmpty = fPickEmpty->IsOn();
  // Update();
}

void AliEveTPCSector2DEditor::DoPickMode(Int_t mode)
{
  fM->fPickMode = mode;
  // Update();
}
