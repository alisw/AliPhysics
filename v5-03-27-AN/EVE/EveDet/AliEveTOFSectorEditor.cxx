/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

#include "AliEveTOFSectorEditor.h"
#include <EveDet/AliEveTOFSector.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TEveGValuators.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>


//______________________________________________________________________________
// AliEveTOFSectorEditor
//

ClassImp(AliEveTOFSectorEditor)

AliEveTOFSectorEditor::AliEveTOFSectorEditor(const TGWindow *p, Int_t width, Int_t height,
					     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0) ,
  fSectorID  (0), fAutoTrans (0),
  fPlate (0),
  fPlate0(0), fPlate1(0), fPlate2(0), fPlate3(0), fPlate4(0),
  fThreshold (0), fMaxVal(0)
{
  //ctr

  fPlate = new TGCheckButton*[5];
  for (Int_t ii=0; ii<5; ii++) fPlate[ii] = new TGCheckButton;

  //fPriority = 40;
  MakeTitle("AliEveTOFSector");

  fSectorID = new TEveGValuator(this, "SectorID", 110, 0);
  fSectorID->SetLabelWidth(60);
  fSectorID->SetShowSlider(kFALSE);
  fSectorID->SetNELength(4);
  fSectorID->Build();
  fSectorID->SetLimits(0, 17);
  fSectorID->SetToolTip("The 18 Tof Sector's");
  fSectorID->Connect("ValueSet(Double_t)",
		     "AliEveTOFSectorEditor", this, "DoSectorID()");
  // Reuse sectorID for auto-transformation button
  fAutoTrans = new TGCheckButton(fSectorID, "AutoTrans");
  fAutoTrans->SetToolTipText("Automatically set transformation to true position");
  fSectorID->AddFrame(fAutoTrans, new TGLayoutHints(kLHintsLeft, 12, 0, 1, 0));
  fAutoTrans->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoAutoTrans()");
  AddFrame(fSectorID, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "AliEveTOFSectorEditor", this, "DoXYZZ()"); {
  TGHorizontalFrame* f = new TGHorizontalFrame(this);

  Int_t nPlate = 0;
  fPlate0 = new TGCheckButton(f, "Plate0");
  f->AddFrame(fPlate0, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fPlate0->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate0()");
  //fPlate0->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate(Int_t)");

  nPlate = 1;
  fPlate1 = new TGCheckButton(f, "Plate 1");
  f->AddFrame(fPlate1, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fPlate1->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate1()");
  //fPlate1->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate(Int_t)");

  nPlate = 2;
  fPlate2 = new TGCheckButton(f, "Plate 2");
  f->AddFrame(fPlate2, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fPlate2->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate2()");
  //fPlate2->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate(Int_t)");

  nPlate = 3;
  fPlate3 = new TGCheckButton(f, "Plate 3");
  f->AddFrame(fPlate3, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fPlate3->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate3()");
  //fPlate3->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate(Int_t)");

  nPlate = 4;
  fPlate4 = new TGCheckButton(f, "Plate 4");
  f->AddFrame(fPlate4, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fPlate4->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate4()");
  //fPlate4->Connect("Toggled(Bool_t)","AliEveTOFSectorEditor", this, "DoPlate(Int_t)");

  fPlate[0] = fPlate0;
  fPlate[1] = fPlate1;
  fPlate[2] = fPlate2;
  fPlate[3] = fPlate3;
  fPlate[4] = fPlate4;

  AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

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

AliEveTOFSectorEditor::~AliEveTOFSectorEditor()
{
  //dtr

  delete [] fPlate;
}

/******************************************************************************/

void AliEveTOFSectorEditor::SetModel(TObject* obj)
{
  // Set object to monitor at visualization level

  fM = static_cast<AliEveTOFSector*>(obj);

  fSectorID->SetValue(fM->GetSectorID());
  fAutoTrans->SetState(fM->GetAutoTrans()  ? kButtonDown : kButtonUp);

  fPlate0->SetState(fM->GetPlate(0) ? kButtonDown : kButtonUp);
  fPlate1->SetState(fM->GetPlate(1) ? kButtonDown : kButtonUp);
  fPlate2->SetState(fM->GetPlate(2) ? kButtonDown : kButtonUp);
  fPlate3->SetState(fM->GetPlate(3) ? kButtonDown : kButtonUp);
  fPlate4->SetState(fM->GetPlate(4) ? kButtonDown : kButtonUp);
}

/******************************************************************************/
void AliEveTOFSectorEditor::DoSectorID()
{
  fM->SetSectorID((Int_t) fSectorID->GetValue());
  Update();
}

void AliEveTOFSectorEditor::DoAutoTrans()
{
  fM->SetAutoTrans(fAutoTrans->IsOn());
  Update();
}

/******************************************************************************/

void AliEveTOFSectorEditor::DoPlate(Int_t nPlate)
{
  fM->SetPlate(nPlate, fPlate[nPlate]->IsOn());
  Update();
}

void AliEveTOFSectorEditor::DoPlate0()
{
  fM->SetPlate(0, fPlate0->IsOn());
  Update();
}

void AliEveTOFSectorEditor::DoPlate1()
{
  fM->SetPlate(1, fPlate1->IsOn());
  Update();
}

void AliEveTOFSectorEditor::DoPlate2()
{
  fM->SetPlate(2, fPlate2->IsOn());
  Update();
}
void AliEveTOFSectorEditor::DoPlate3()
{
  fM->SetPlate(3, fPlate3->IsOn());
  Update();
}

void AliEveTOFSectorEditor::DoPlate4()
{
  fM->SetPlate(4, fPlate4->IsOn());
  Update();
}


void AliEveTOFSectorEditor::DoThreshold()
{
  fM->SetThreshold((Short_t) fThreshold->GetValue());
  fThreshold->SetValue(fM->GetThreshold());
  Update();
}

void AliEveTOFSectorEditor::DoMaxVal()
{
  fM->SetMaxVal((Int_t) fMaxVal->GetValue());
  fMaxVal->SetValue(fM->GetMaxVal());
  Update();
}

/*
void AliEveTOFSectorEditor::DoTime()
{ 
  fM->SetMinTime((Int_t) fTime->GetMin());
  fM->SetMaxTime((Int_t) fTime->GetMax());
  Update();
}
*/
