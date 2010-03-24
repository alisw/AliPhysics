/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

#include "AliEveTOFStripEditor.h"
#include <EveDet/AliEveTOFStrip.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TEveGValuators.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

//
// AliEveTOFStripEditor class
// Editor for AliEveTOFStrip class
//
// Author A. De Caro (email: decaro@sa.infn.it)
//

ClassImp(AliEveTOFStripEditor)

AliEveTOFStripEditor::AliEveTOFStripEditor(const TGWindow *p, Int_t width, Int_t height,
                                           UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM         (0),
  fThreshold (0),
  fMaxVal    (0)
  // Initialize widget pointers to 0
{
  //ctr

  MakeTitle("AliEveTOFStrip");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "AliEveTOFStripEditor", this, "DoXYZZ()");

  fThreshold = new TEveGValuator(this, "Threshold", 200, 0);
  fThreshold->SetNELength(4);
  fThreshold->SetLabelWidth(60);
  fThreshold->Build();
  fThreshold->GetSlider()->SetWidth(120);
  fThreshold->SetLimits(0,250);
  fThreshold->Connect("ValueSet(Double_t)",
		      "AliEveTOFSectorEditor", this, "DoThreshold()");
  AddFrame(fThreshold, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fMaxVal = new TEveGValuator(this,"MaxVal", 200, 0);
  fMaxVal->SetNELength(4);
  fMaxVal->SetLabelWidth(60);
  fMaxVal->Build();
  fMaxVal->GetSlider()->SetWidth(60);
  fMaxVal->SetLimits(0, 500);
  fMaxVal->Connect("ValueSet(Double_t)",
		   "AliEveTOFSectorEditor", this, "DoMaxVal()");
  AddFrame(fMaxVal, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

}

/******************************************************************************/

void AliEveTOFStripEditor::SetModel(TObject* obj)
{
  // Set object to monitor at visualization level

  fM = dynamic_cast<AliEveTOFStrip*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/******************************************************************************/
void AliEveTOFStripEditor::DoThreshold()
{
  fM->SetThreshold((Short_t) fThreshold->GetValue());
  fThreshold->SetValue(fM->GetThreshold());
  Update();
}

/******************************************************************************/
void AliEveTOFStripEditor::DoMaxVal()
{
  fM->SetMaxVal((Int_t) fMaxVal->GetValue());
  fMaxVal->SetValue(fM->GetMaxVal());
  Update();
}

// Implements callback/slot methods

// void AliEveTOFStripEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
