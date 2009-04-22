// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/



#include "TVirtualPad.h"
#include "TColor.h"

#include "TGLabel.h"
#include "TEveGValuators.h"

// Cleanup these includes:

//#include "TGButton.h"
//#include "TGNumberEntry.h"
//#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include <TGComboBox.h>

#include "AliEveCascadeListEditor.h"
#include "AliEveCascade.h"



//______________________________________________________________________________
// GUI editor for AliEveCascadeList.
//

ClassImp(AliEveCascadeListEditor)

//______________________________________________________________________________
AliEveCascadeListEditor::AliEveCascadeListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fCascadeSpecies(0),
  fMinMaxRCut(0),
  fMinMaxDaughterDCA(0),
  fMinMaxPt(0),
  fMinMaxInvariantMass(0)
{
  // Constructor.

  MakeTitle("AliEveCascadeList");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::AliEveCascadeListEditor", this, "DoXYZZ()");

   fMinMaxRCut = new TEveGDoubleValuator(this,"Radius:", 130, 0);
   fMinMaxRCut->SetNELength(5);
   fMinMaxRCut->SetLabelWidth(74);
   fMinMaxRCut->Build();
   fMinMaxRCut->GetSlider()->SetWidth(200);
   fMinMaxRCut->SetLimits(0, 100, TGNumberFormat::kNESRealOne);
   fMinMaxRCut->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxRCut()");
   AddFrame(fMinMaxRCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxDaughterDCA = new TEveGDoubleValuator(this,"DCA:", 130, 0);
   fMinMaxDaughterDCA->SetNELength(5);
   fMinMaxDaughterDCA->SetLabelWidth(74);
   fMinMaxDaughterDCA->Build();
   fMinMaxDaughterDCA->GetSlider()->SetWidth(200);
   fMinMaxDaughterDCA->SetLimits(0, 1, TGNumberFormat::kNESRealTwo);
   fMinMaxDaughterDCA->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxDaughterDCA()");
   AddFrame(fMinMaxDaughterDCA, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxPt = new TEveGDoubleValuator(this,"pT:", 130, 0);
   fMinMaxPt->SetNELength(5);
   fMinMaxPt->SetLabelWidth(74);
   fMinMaxPt->Build();
   fMinMaxPt->GetSlider()->SetWidth(200);
   fMinMaxPt->SetLimits(0, 20, TGNumberFormat::kNESRealOne);
   fMinMaxPt->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxPt()");
   AddFrame(fMinMaxPt, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
   
   
   TGHorizontalFrame* fCascadeFrame = new TGHorizontalFrame(this);
   TGLabel* labPos = new TGLabel(fCascadeFrame, "Cascade:");
   fCascadeFrame->AddFrame(labPos, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 1, 1, 1));
   fCascadeSpecies = new TGComboBox(fCascadeFrame);
   fCascadeSpecies->AddEntry("Xi",     kXiMinus);
   fCascadeSpecies->AddEntry("Omega",  kOmegaMinus);
   TGListBox* lbPos = fCascadeSpecies->GetListBox();
   lbPos->Resize(lbPos->GetWidth(), 2*18);
   fCascadeSpecies->Resize(45, 20);
   fCascadeSpecies->Connect("Selected(Int_t)", "AliEveCascadeListEditor", this, "DoSelectInvMassHyp(Int_t)");
   fCascadeFrame->AddFrame(fCascadeSpecies, new TGLayoutHints(kLHintsTop, 3, 1, 1, 1));
   AddFrame(fCascadeFrame);
   
   
   fMinMaxInvariantMass = new TEveGDoubleValuator(this,"Inv. Mass:", 130, 0);
   fMinMaxInvariantMass->SetNELength(5);
   fMinMaxInvariantMass->SetLabelWidth(74);
   fMinMaxInvariantMass->Build();
   fMinMaxInvariantMass->GetSlider()->SetWidth(200);
   fMinMaxInvariantMass->SetLimits(1, 6, TGNumberFormat::kNESRealThree);
   fMinMaxInvariantMass->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxInvariantMass()");
   AddFrame(fMinMaxInvariantMass, new TGLayoutHints(kLHintsBottom, 1, 1, 1, 1));
   
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveCascadeList*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxRCut->SetValues(fM->fMinRCut, fM->fMaxRCut);
  fMinMaxDaughterDCA->SetValues(fM->fMinDaughterDCA, fM->fMaxDaughterDCA);
  fMinMaxPt->SetValues(fM->fMinPt, fM->fMaxPt);
  fMinMaxInvariantMass->SetValues(fM->fMinInvariantMass,fM->fMaxInvariantMass);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveCascadeListEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

void AliEveCascadeListEditor::DoMinMaxRCut()
{
  // Filter cascade candidates by transverse radius	
	
  fM->FilterByRadius(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}

void AliEveCascadeListEditor::DoMinMaxDaughterDCA()
{
  // Filter cascade candidates by DCA to primary vertex	
	
  fM->FilterByDaughterDCA(fMinMaxDaughterDCA->GetMin(), fMinMaxDaughterDCA->GetMax());
}

void AliEveCascadeListEditor::DoMinMaxPt()
{
   // Filter cascade candidates by transverse momentum
	
  fM->FilterByPt(fMinMaxPt->GetMin(), fMinMaxPt->GetMax());
}



void AliEveCascadeListEditor::DoSelectInvMassHyp(Int_t rInvMassHyp)
{
   // Apply the invariant mass hypothesis according the choice made by user
	
	fM->SetInvMassHyp(rInvMassHyp);
	Update();
}

void AliEveCascadeListEditor::DoMinMaxInvariantMass()
{
    // Filter cascade candidates by invariant mass (under mass hypothesis made by user)
	
	Int_t rInvMassHyp = fM->GetInvMassHyp();
	if(rInvMassHyp)
	fM->FilterByInvariantMass(fMinMaxInvariantMass->GetMin(), fMinMaxInvariantMass->GetMax(), rInvMassHyp);
	else 
	fM->FilterByInvariantMass(fMinMaxInvariantMass->GetMin(), fMinMaxInvariantMass->GetMax(), kXiMinus);
	
}

