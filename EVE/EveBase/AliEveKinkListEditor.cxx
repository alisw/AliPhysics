// $Id$
// Author: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveKinkListEditor.h"
#include "AliEveKink.h"

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
// GUI editor for AliEveKinkList.
//

ClassImp(AliEveKinkListEditor)

//______________________________________________________________________________
AliEveKinkListEditor::AliEveKinkListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fMinMaxRCut(0),
  fMinMaxKinkAngleCut (0),
  fMinMaxPt(0),
  fMinMaxInvariantMass(0),
  fDaughterSpecies(0),
  fDaughterCheckMaxPidProbability(0),
  fDaughterLevelPidProbability(0)
{
  // Constructor.

  MakeTitle("AliEveKinkList");

   fMinMaxRCut = new TEveGDoubleValuator(this,"Radius:", 130, 0);
   fMinMaxRCut->SetNELength(5);
   fMinMaxRCut->SetLabelWidth(74);
   fMinMaxRCut->Build();
   fMinMaxRCut->GetSlider()->SetWidth(200);
   fMinMaxRCut->SetLimits(0, 250, TGNumberFormat::kNESRealOne);
   fMinMaxRCut->Connect("ValueSet()", "AliEveKinkListEditor", this, "DoMinMaxRCut()");
   AddFrame(fMinMaxRCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
   
   fMinMaxKinkAngleCut = new TEveGDoubleValuator(this,"Kink Angle:", 130, 0);
   fMinMaxKinkAngleCut->SetNELength(5);
   fMinMaxKinkAngleCut->SetLabelWidth(74);
   fMinMaxKinkAngleCut->Build();
   fMinMaxKinkAngleCut->GetSlider()->SetWidth(200);
   fMinMaxKinkAngleCut->SetLimits(0, 280, TGNumberFormat::kNESRealOne);
   fMinMaxKinkAngleCut->Connect("ValueSet()", "AliEveKinkListEditor", this, "DoMinMaxKinkAngleCut()");
   AddFrame(fMinMaxKinkAngleCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));  

   fMinMaxPt = new TEveGDoubleValuator(this,"pT:", 80, 0);
   fMinMaxPt->SetNELength(5);
   fMinMaxPt->SetLabelWidth(74);
   fMinMaxPt->Build();
   fMinMaxPt->GetSlider()->SetWidth(200);
   fMinMaxPt->SetLimits(0, 20, TGNumberFormat::kNESRealOne);
   fMinMaxPt->Connect("ValueSet()", "AliEveKinkListEditor", this, "DoMinMaxPt()");
   AddFrame(fMinMaxPt, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxInvariantMass = new TEveGDoubleValuator(this,"Inv. Mass:", 80, 0);
   fMinMaxInvariantMass->SetNELength(5);
   fMinMaxInvariantMass->SetLabelWidth(74);
   fMinMaxInvariantMass->Build();
   fMinMaxInvariantMass->GetSlider()->SetWidth(200);
   fMinMaxInvariantMass->SetLimits(0, 1.0, TGNumberFormat::kNESRealThree);
   fMinMaxInvariantMass->Connect("ValueSet()", "AliEveKinkListEditor", this, "DoMinMaxInvariantMass()");
   AddFrame(fMinMaxInvariantMass, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
   
   TGHorizontalFrame* fDaugFrame = new TGHorizontalFrame(this);
   TGLabel* labDaug = new TGLabel(fDaugFrame, "Daug:");
   fDaugFrame->AddFrame(labDaug, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 1, 1, 1));
   fDaughterSpecies = new TGComboBox(fDaugFrame);
   fDaughterSpecies->AddEntry("e", 11);
   fDaughterSpecies->AddEntry("mu", 13);
   fDaughterSpecies->AddEntry("pi", 211);
   TGListBox* lbaug = fDaughterSpecies->GetListBox();
   lbaug->Resize(lbaug->GetWidth(), 2*18);
   fDaughterSpecies->Resize(45, 20);
   fDaughterSpecies->Connect("Selected(Int_t)", "AliEveKinkListEditor", this, "DoSelectDaugPid(Int_t)");
   fDaugFrame->AddFrame(fDaughterSpecies, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
   AddFrame(fDaugFrame);

   fDaughterCheckMaxPidProbability = new TGCheckButton(fDaugFrame, "Check");
   fDaugFrame->AddFrame(fDaughterCheckMaxPidProbability, new TGLayoutHints(kLHintsLeft, 0, 2, 1, 1));
   fDaughterCheckMaxPidProbability->Connect("Toggled(Bool_t)", "AliEveKinkListEditor", this, "DoCheckDaugPid()");

   fDaughterLevelPidProbability = new TGNumberEntry(fDaugFrame, 0.5, 3, -1, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax, 0, 1);
   fDaughterLevelPidProbability->Resize(50,20);
   fDaughterLevelPidProbability->Connect("ValueSet(Long_t)", "AliEveKinkListEditor", this, "DoSelectDaugProb()");
   fDaugFrame->AddFrame(fDaughterLevelPidProbability, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 1, 1)); 
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveKinkList*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxRCut->SetValues(fM->fMinRCut, fM->fMaxRCut);
  fMinMaxKinkAngleCut->SetValues(fM->fMinKinkAngle, fM->fMaxKinkAngle);
  fMinMaxPt->SetValues(fM->fMinPt, fM->fMaxPt);
  fMinMaxInvariantMass->SetValues(fM->fMinInvariantMass, fM->fMaxInvariantMass);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________

void AliEveKinkListEditor::DoMinMaxRCut()
{
  fM->FilterByRadius(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}

void AliEveKinkListEditor::DoMinMaxKinkAngleCut()
{
  fM->FilterByKinkAngle(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}

void AliEveKinkListEditor::DoMinMaxPt()
{
  fM->FilterByPt(fMinMaxPt->GetMin(), fMinMaxPt->GetMax());
}

void AliEveKinkListEditor::DoMinMaxInvariantMass()
{
    fM->FilterByInvariantMass(fMinMaxInvariantMass->GetMin(), fMinMaxInvariantMass->GetMax(),13);
}

void AliEveKinkListEditor::DoSelectDaugPid(Int_t rDaugPid)
{
  fM->SetDaugCheckedPid(rDaugPid);
  Update();
}

void AliEveKinkListEditor::DoCheckDaugPid()
{
  Int_t   lDaugPid  = fM->GetDaugCheckedPid();
  Float_t lDaugProb = fM->GetDaugCheckedProb();
  if (lDaugPid) {
    fM->FilterByCheckedPidMinProb(fDaughterCheckMaxPidProbability->IsOn(),lDaugPid,lDaugProb);
    printf("Selection for daughter pid %d prob %.2f \n",lDaugPid,lDaugProb);
    Update();
  }
}

void AliEveKinkListEditor::DoSelectDaugProb()
{
  Float_t rMinDaugProb = (Float_t)fDaughterLevelPidProbability->GetNumber();
  fM->SetDaugCheckedProb(rMinDaugProb);
  Update();
}
