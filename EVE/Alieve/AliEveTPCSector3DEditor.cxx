// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector3DEditor.h"
#include <Alieve/AliEveTPCSector3D.h>

#include <TEveGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>


//______________________________________________________________________
// AliEveTPCSector3DEditor
//

ClassImp(AliEveTPCSector3DEditor)

AliEveTPCSector3DEditor::AliEveTPCSector3DEditor(const TGWindow *p,
                                     Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fRnrFrame(0), fDriftVel(0), fPointFrac(0), fPointSize(0)
{
  MakeTitle("AliEveTPCSector3D");

  Int_t labelW = 60;

  fRnrFrame = new TGCheckButton(this, "ShowFrame");
  AddFrame(fRnrFrame, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrFrame->Connect
    ("Toggled(Bool_t)","AliEveTPCSector3DEditor", this, "DoRnrFrame()");

  fDriftVel = new TEveGValuator(this, "Vdrift fac", 110, 0);
  fDriftVel->SetLabelWidth(labelW);
  fDriftVel->SetShowSlider(kFALSE);
  fDriftVel->SetNELength(6);
  fDriftVel->Build();
  fDriftVel->SetLimits(0.1, 10, 1, TGNumberFormat::kNESRealThree);
  fDriftVel->SetToolTip("Drift velocity factor");
  fDriftVel->Connect("ValueSet(Double_t)",
		     "AliEveTPCSector3DEditor", this, "DoDriftVel()");
  AddFrame(fDriftVel, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fPointFrac = new TEveGValuator(this,"Point frac", 200, 0);
  fPointFrac->SetLabelWidth(labelW);
  fPointFrac->SetNELength(4);
  fPointFrac->Build();
  fPointFrac->GetSlider()->SetWidth(101 + 16);
  fPointFrac->SetLimits(0.0, 1.0, 101);
  fPointFrac->SetToolTip("Fraction of signal range displayed as points");
  fPointFrac->Connect("ValueSet(Double_t)",
		      "AliEveTPCSector3DEditor", this, "DoPointFrac()");
  AddFrame(fPointFrac, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fPointSize = new TEveGValuator(this,"Point size", 200, 0);
  fPointSize->SetLabelWidth(labelW);
  fPointSize->SetShowSlider(kFALSE);
  fPointSize->SetNELength(4);
  fPointSize->Build();
  //fPointSize->GetSlider()->SetWidth(101 + 16);
  fPointSize->SetLimits(0.1, 32.0, 1, TGNumberFormat::kNESRealOne);
  fPointSize->SetToolTip("Size of displayed points");
  fPointSize->Connect("ValueSet(Double_t)",
		      "AliEveTPCSector3DEditor", this, "DoPointSize()");
  AddFrame(fPointSize, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
}

AliEveTPCSector3DEditor::~AliEveTPCSector3DEditor()
{}

/**************************************************************************/

void AliEveTPCSector3DEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveTPCSector3D*>(obj);

  fRnrFrame->SetState(fM->fRnrFrame ? kButtonDown : kButtonUp);
  fDriftVel->SetValue(fM->fDriftVel);

  fPointFrac->SetValue(fM->fPointFrac);
  fPointSize->SetValue(fM->fPointSize);
}

/**************************************************************************/

void AliEveTPCSector3DEditor::DoRnrFrame()
{
  fM->SetRnrFrame(fRnrFrame->IsOn());
  Update();
}

void AliEveTPCSector3DEditor::DoDriftVel()
{
  fM->SetDriftVel(fDriftVel->GetValue());
  Update();
}

void AliEveTPCSector3DEditor::DoPointFrac()
{
  fM->SetPointFrac(fPointFrac->GetValue());
  Update();
}

void AliEveTPCSector3DEditor::DoPointSize()
{
  fM->SetPointSize(fPointSize->GetValue());
  Update();
}

