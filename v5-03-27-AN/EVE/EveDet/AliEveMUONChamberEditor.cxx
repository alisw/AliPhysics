// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMUONChamberEditor.h"

#include <EveDet/AliEveMUONChamber.h>

#include <TEveGValuators.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>


//______________________________________________________________________________
// AliEveMUONChamberEditor
//

ClassImp(AliEveMUONChamberEditor)

//______________________________________________________________________________
AliEveMUONChamberEditor::AliEveMUONChamberEditor(const TGWindow *p,
                                     Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fThreshold(0),
  fMaxVal(0),
  fClusterSize(0),
  fHitSize(0)
{
  //
  // constructor
  //

  MakeTitle("AliEveMUONChamber");

  Int_t labelW = 60;

  fThreshold = new TEveGValuator(this, "ADC min", 200, 0);
  fThreshold->SetNELength(4);
  fThreshold->SetLabelWidth(labelW);
  fThreshold->Build();
  fThreshold->GetSlider()->SetWidth(120);
  fThreshold->SetLimits(0,4096);
  fThreshold->Connect("ValueSet(Double_t)",
                      "AliEveMUONChamberEditor", this, "DoThreshold()");
  AddFrame(fThreshold, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fMaxVal = new TEveGValuator(this,"ADC max", 200, 0);
  fMaxVal->SetNELength(4);
  fMaxVal->SetLabelWidth(labelW);
  fMaxVal->Build();
  fMaxVal->GetSlider()->SetWidth(120);
  fMaxVal->SetLimits(0, 4096);
  fMaxVal->Connect("ValueSet(Double_t)",
                   "AliEveMUONChamberEditor", this, "DoMaxVal()");
  AddFrame(fMaxVal, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fClusterSize = new TEveGValuator(this,"Cls size", 200, 0);
  fClusterSize->SetLabelWidth(labelW);
  fClusterSize->SetShowSlider(kFALSE);
  fClusterSize->SetNELength(4);
  fClusterSize->Build();
  fClusterSize->SetLimits(0, 24);
  fClusterSize->SetToolTip("Size of displayed clusters");
  fClusterSize->Connect("ValueSet(Double_t)",
                      "AliEveMUONChamberEditor", this, "DoClusterSize()");
  AddFrame(fClusterSize, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fHitSize = new TEveGValuator(this,"TEveHit size", 200, 0);
  fHitSize->SetLabelWidth(labelW);
  fHitSize->SetShowSlider(kFALSE);
  fHitSize->SetNELength(4);
  fHitSize->Build();
  fHitSize->SetLimits(0, 24);
  fHitSize->SetToolTip("Size of displayed clusters");
  fHitSize->Connect("ValueSet(Double_t)",
                      "AliEveMUONChamberEditor", this, "DoHitSize()");
  AddFrame(fHitSize, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

}

//______________________________________________________________________________
AliEveMUONChamberEditor::~AliEveMUONChamberEditor()
{
  //
  // destructor
  //

}

//______________________________________________________________________________
void AliEveMUONChamberEditor::SetModel(TObject* obj)
{
  //
  // ...
  //

  fM = static_cast<AliEveMUONChamber*>(obj);

  fThreshold->SetValue(fM->fThreshold);
  fMaxVal->SetValue(fM->fMaxVal);
  fClusterSize->SetValue(fM->fClusterSize);
  fHitSize->SetValue(fM->fHitSize);

}

//______________________________________________________________________________
void AliEveMUONChamberEditor::DoThreshold()
{
  //
  // set digit minimum amplitude
  //

  fM->SetThreshold((Short_t) fThreshold->GetValue());
  fThreshold->SetValue(fM->fThreshold);
  Update();

}

//______________________________________________________________________________
void AliEveMUONChamberEditor::DoMaxVal()
{
  //
  // set digit maximum amplitude
  //

  fM->SetMaxVal((Int_t) fMaxVal->GetValue());
  fMaxVal->SetValue(fM->fMaxVal);
  Update();

}

//______________________________________________________________________________
void AliEveMUONChamberEditor::DoClusterSize()
{
  //
  // set the cluster point size
  //

  fM->SetClusterSize((Int_t) fClusterSize->GetValue());
  fClusterSize->SetValue(fM->fClusterSize);
  Update();

}

//______________________________________________________________________________
void AliEveMUONChamberEditor::DoHitSize()
{
  //
  // set the hit point size
  //

  fM->SetHitSize((Int_t) fHitSize->GetValue());
  fHitSize->SetValue(fM->fHitSize);
  Update();

}
