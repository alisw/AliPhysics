// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveV0ListEditor.h"
#include "AliEveV0.h"

#include "TEveGValuators.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include "TGComboBox.h"
#include "TGLabel.h"

//______________________________________________________________________________
// GUI editor for AliEveV0List.
//

ClassImp(AliEveV0ListEditor)

//______________________________________________________________________________
AliEveV0ListEditor::AliEveV0ListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fMinMaxRCut(0),
  fMinMaxDaughterDCA(0),
  fMinMaxPt(0),
  fNegativeSpecies(0),
  fPositiveSpecies(0),
  fNegativeCheckMaxPidProbability(0),
  fPositiveCheckMaxPidProbability(0),
  fNegativeLevelPidProbability(0),
  fPositiveLevelPidProbability(0),
  fMinMaxInvariantMass(0)
{
  // Constructor.

  MakeTitle("AliEveV0List");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::AliEveV0ListEditor", this, "DoXYZZ()");

   fMinMaxRCut = new TEveGDoubleValuator(this,"Radius:", 130, 0);
   fMinMaxRCut->SetNELength(5);
   fMinMaxRCut->SetLabelWidth(74);
   fMinMaxRCut->Build();
   fMinMaxRCut->GetSlider()->SetWidth(200);
   fMinMaxRCut->SetLimits(0, 100, TGNumberFormat::kNESRealOne);
   fMinMaxRCut->Connect("ValueSet()", "AliEveV0ListEditor", this, "DoMinMaxRCut()");
   AddFrame(fMinMaxRCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxDaughterDCA = new TEveGDoubleValuator(this,"DCA:", 130, 0);
   fMinMaxDaughterDCA->SetNELength(5);
   fMinMaxDaughterDCA->SetLabelWidth(74);
   fMinMaxDaughterDCA->Build();
   fMinMaxDaughterDCA->GetSlider()->SetWidth(200);
   fMinMaxDaughterDCA->SetLimits(0, 1, TGNumberFormat::kNESRealTwo);
   fMinMaxDaughterDCA->Connect("ValueSet()", "AliEveV0ListEditor", this, "DoMinMaxDaughterDCA()");
   AddFrame(fMinMaxDaughterDCA, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxPt = new TEveGDoubleValuator(this,"pT:", 80, 0);
   fMinMaxPt->SetNELength(5);
   fMinMaxPt->SetLabelWidth(74);
   fMinMaxPt->Build();
   fMinMaxPt->GetSlider()->SetWidth(200);
   fMinMaxPt->SetLimits(0, 20, TGNumberFormat::kNESRealOne);
   fMinMaxPt->Connect("ValueSet()", "AliEveV0ListEditor", this, "DoMinMaxPt()");
   AddFrame(fMinMaxPt, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   TGHorizontalFrame* fNegFrame = new TGHorizontalFrame(this);
   TGLabel* labNeg = new TGLabel(fNegFrame, "Neg:");
   fNegFrame->AddFrame(labNeg, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 1, 1, 1));
   fNegativeSpecies = new TGComboBox(fNegFrame);
   fNegativeSpecies->AddEntry("e", 11);
   fNegativeSpecies->AddEntry("mu", 13);
   fNegativeSpecies->AddEntry("pi", 211);
   fNegativeSpecies->AddEntry("K", 321);
   fNegativeSpecies->AddEntry("p", 2212);
   TGListBox* lbNeg = fNegativeSpecies->GetListBox();
   lbNeg->Resize(lbNeg->GetWidth(), 2*18);
   fNegativeSpecies->Resize(45, 20);
   fNegativeSpecies->Connect("Selected(Int_t)", "AliEveV0ListEditor", this, "DoSelectNegPid(Int_t)");
   fNegFrame->AddFrame(fNegativeSpecies, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
   AddFrame(fNegFrame);

   fNegativeCheckMaxPidProbability = new TGCheckButton(fNegFrame, "Check");
   fNegFrame->AddFrame(fNegativeCheckMaxPidProbability, new TGLayoutHints(kLHintsLeft, 0, 2, 1, 1));
   fNegativeCheckMaxPidProbability->Connect("Toggled(Bool_t)", "AliEveV0ListEditor", this, "DoCheckNegPid()");

   fNegativeLevelPidProbability = new TGNumberEntry(fNegFrame, 0.5, 3, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax, 0, 1);
   fNegativeLevelPidProbability->Resize(50,20);
   fNegativeLevelPidProbability->Connect("ValueSet(Long_t)", "AliEveV0ListEditor", this, "DoSelectNegProb()");
   fNegFrame->AddFrame(fNegativeLevelPidProbability, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 1, 1));

   TGHorizontalFrame* fPosFrame = new TGHorizontalFrame(this);
   TGLabel* labPos = new TGLabel(fPosFrame, "Pos:");
   fPosFrame->AddFrame(labPos, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 1, 1, 1));
   fPositiveSpecies = new TGComboBox(fPosFrame);
   fPositiveSpecies->AddEntry("e", 11);
   fPositiveSpecies->AddEntry("mu", 13);
   fPositiveSpecies->AddEntry("pi", 211);
   fPositiveSpecies->AddEntry("K", 321);
   fPositiveSpecies->AddEntry("p", 2212);
   TGListBox* lbPos = fPositiveSpecies->GetListBox();
   lbPos->Resize(lbPos->GetWidth(), 2*18);
   fPositiveSpecies->Resize(45, 20);
   fPositiveSpecies->Connect("Selected(Int_t)", "AliEveV0ListEditor", this, "DoSelectPosPid(Int_t)");
   fPosFrame->AddFrame(fPositiveSpecies, new TGLayoutHints(kLHintsTop, 3, 1, 1, 1));
   AddFrame(fPosFrame);

   fPositiveCheckMaxPidProbability = new TGCheckButton(fPosFrame, "Check");
   fPosFrame->AddFrame(fPositiveCheckMaxPidProbability, new TGLayoutHints(kLHintsLeft, 0, 2, 1, 1));
   fPositiveCheckMaxPidProbability->Connect("Toggled(Bool_t)", "AliEveV0ListEditor", this, "DoCheckPosPid()");

   fPositiveLevelPidProbability = new TGNumberEntry(fPosFrame, 0.5, 3, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax, 0, 1);
   fPositiveLevelPidProbability->Resize(50,20);
   fPositiveLevelPidProbability->Connect("ValueSet(Long_t)", "AliEveV0ListEditor", this, "DoSelectPosProb()");
   fPosFrame->AddFrame(fPositiveLevelPidProbability, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 1, 1));

   fMinMaxInvariantMass = new TEveGDoubleValuator(this,"Inv. Mass:", 80, 0);
   fMinMaxInvariantMass->SetNELength(5);
   fMinMaxInvariantMass->SetLabelWidth(74);
   fMinMaxInvariantMass->Build();
   fMinMaxInvariantMass->GetSlider()->SetWidth(200);
   fMinMaxInvariantMass->SetLimits(0, 1.2, TGNumberFormat::kNESRealThree);
   fMinMaxInvariantMass->Connect("ValueSet()", "AliEveV0ListEditor", this, "DoMinMaxInvariantMass()");
   AddFrame(fMinMaxInvariantMass, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0ListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveV0List*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxRCut->SetValues(fM->fMinRCut, fM->fMaxRCut);
  fMinMaxDaughterDCA->SetValues(fM->fMinDaughterDCA, fM->fMaxDaughterDCA);
  fMinMaxPt->SetValues(fM->fMinPt, fM->fMaxPt);
  fMinMaxInvariantMass->SetValues(fM->fMinInvariantMass, fM->fMaxInvariantMass);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveV0ListEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

void AliEveV0ListEditor::DoMinMaxRCut()
{
  fM->FilterByRadius(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}

void AliEveV0ListEditor::DoMinMaxDaughterDCA()
{
  fM->FilterByDaughterDCA(fMinMaxDaughterDCA->GetMin(), fMinMaxDaughterDCA->GetMax());
}

void AliEveV0ListEditor::DoMinMaxPt()
{
  fM->FilterByPt(fMinMaxPt->GetMin(), fMinMaxPt->GetMax());
}

void AliEveV0ListEditor::DoSelectNegPid(Int_t rNegPid)
{
  fM->SetNegCheckedPid(rNegPid);
  Update();
}

void AliEveV0ListEditor::DoCheckNegPid()
{
  Int_t   lNegPid  = fM->GetNegCheckedPid();
  Float_t lNegProb = fM->GetNegCheckedProb();
  if (lNegPid) {
    fM->FilterByCheckedPidMinProb(fNegativeCheckMaxPidProbability->IsOn(),0,lNegPid,lNegProb);
    printf("Selection for negative daughter pid %d prob %.2f \n",lNegPid,lNegProb);
    Update();
  }
}

void AliEveV0ListEditor::DoSelectNegProb()
{
  Float_t rMinNegProb = (Float_t)fNegativeLevelPidProbability->GetNumber();
  fM->SetNegCheckedProb(rMinNegProb);
  Update();
}

void AliEveV0ListEditor::DoSelectPosPid(Int_t rPosPid)
{
  fM->SetPosCheckedPid(rPosPid);
  Update();
}

void AliEveV0ListEditor::DoCheckPosPid()
{
  Int_t   lPosPid  = fM->GetPosCheckedPid();
  Float_t lPosProb = fM->GetPosCheckedProb();
  if (lPosPid) {
    fM->FilterByCheckedPidMinProb(fPositiveCheckMaxPidProbability->IsOn(),1,lPosPid,lPosProb);
    printf("Selection for positive daughter pid %d prob %.2f \n",lPosPid,lPosProb);
    Update();
  }
}

void AliEveV0ListEditor::DoSelectPosProb()
{
  Float_t rMinPosProb = (Float_t)fPositiveLevelPidProbability->GetNumber();
  fM->SetPosCheckedProb(rMinPosProb);
  Update();
}

void AliEveV0ListEditor::DoMinMaxInvariantMass()
{
  Int_t lNegPid = fM->GetNegCheckedPid();
  Int_t lPosPid = fM->GetPosCheckedPid();
  if( lNegPid && lPosPid)
    fM->FilterByInvariantMass(fMinMaxInvariantMass->GetMin(), fMinMaxInvariantMass->GetMax(),lNegPid,lPosPid);
  else 
    fM->FilterByInvariantMass(fMinMaxInvariantMass->GetMin(), fMinMaxInvariantMass->GetMax(),211,-211);
}
