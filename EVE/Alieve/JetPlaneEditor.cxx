// $Header$

#include "JetPlaneEditor.h"
#include <Alieve/JetPlane.h>
#include <TEveGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>
#include <TROOT.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
#include <TGFrame.h>
#include <TGTab.h>
using namespace Alieve;

//______________________________________________________________________
// JetPlaneEditor
//

Alieve::JetPlaneEditor::StaticDataWindow* JetPlaneEditor::fgStaticWindow = 0;

ClassImp(JetPlaneEditor)

JetPlaneEditor::JetPlaneEditor(const TGWindow *p, Int_t width, Int_t height,
			       UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fRnrJets(0),
  fRnrTracks(0),
  fEnergyScale(0),
  fEnergyColorScale(0),
  fOneSelection(0),
  fTwoSelection(0),
  fInformationSetup(0)
  // Initialize widget pointers to 0
{
  MakeTitle("JetPlane");
  Int_t labelW = 67;

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::JetPlaneEditor", this, "DoXYZZ()");

  fRnrJets  = new TGCheckButton(this, "Rnr Jets");
  AddFrame(fRnrJets, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  fRnrJets->Connect("Clicked()", "Alieve::JetPlaneEditor", this, "DoRnrJets()");

  fRnrTracks  = new TGCheckButton(this, "Rnr Tracks");
  AddFrame(fRnrTracks, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  fRnrTracks->Connect("Clicked()", "Alieve::JetPlaneEditor", this, "DoRnrTracks()");

  fEnergyScale = new TEveGValuator(this, "Length scale:", 110, 0);
  fEnergyScale->SetLabelWidth(labelW);
  fEnergyScale->SetNELength(6);
  fEnergyScale->Build();
  fEnergyScale->SetLimits(1, 500, 500, TGNumberFormat::kNESRealOne);
  fEnergyScale->SetToolTip("Energy mapped to length of arrow.");
  fEnergyScale->Connect("ValueSet(Double_t)", "Alieve::JetPlaneEditor", this, "DoEnergyScale()");
  AddFrame(fEnergyScale, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fEnergyColorScale = new TEveGValuator(this, "Color scale:", 110, 0);
  fEnergyColorScale->SetLabelWidth(labelW);
  fEnergyColorScale->SetNELength(6);
  fEnergyColorScale->Build();
  fEnergyColorScale->SetLimits(-2, 2, 100, TGNumberFormat::kNESRealOne);
  fEnergyColorScale->SetToolTip("Energy mapped to highest palette color.");
  fEnergyColorScale->Connect("ValueSet(Double_t)", "Alieve::JetPlaneEditor", this, "DoEnergyColorScale()");
  AddFrame(fEnergyColorScale, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fOneSelection = new TGRadioButton(this, "&One TEveTrack/Jet");
  fTwoSelection = new TGRadioButton(this, "&Two TEveTrack/Jet");
  AddFrame(fOneSelection, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  AddFrame(fTwoSelection, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  fOneSelection->Connect("Clicked()", "Alieve::JetPlaneEditor", this, "DoOneSelection()");
  fTwoSelection->Connect("Clicked()", "Alieve::JetPlaneEditor", this, "DoTwoSelection()");

  // fInformationSetup = new TGTextButton(this, "TEveTrack/Jet Print");
  // AddFrame(fInformationSetup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 0, 2, 2));
  // fInformationSetup->Connect("Clicked()", "Alieve::JetPlaneEditor", this, "DoStaticDataWindow()");
}

JetPlaneEditor::~JetPlaneEditor()
{}

/**************************************************************************/

void JetPlaneEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<JetPlane*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
  fRnrJets->SetState(fM->GetRnrJets() ? kButtonDown : kButtonUp);
  fRnrTracks->SetState(fM->GetRnrTracks() ? kButtonDown : kButtonUp);
  fEnergyScale->SetValue(fM->GetEnergyScale());
  fEnergyColorScale->SetValue(fM->GetEnergyColorScale());
  fOneSelection->SetState(fM->GetOneSelection() ? kButtonDown : kButtonUp);
  fTwoSelection->SetState(fM->GetTwoSelection() ? kButtonDown : kButtonUp);
}

/**************************************************************************/

// Implements callback/slot methods

// void JetPlaneEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }

void JetPlaneEditor::DoRnrJets()
{
  fM->SetRnrJets(fRnrJets->IsOn());
  Update();
}

void JetPlaneEditor::DoRnrTracks()
{
  fM->SetRnrTracks(fRnrTracks->IsOn());
  Update();
}

void JetPlaneEditor::DoEnergyColorScale()
{
  fM->SetEnergyColorScale(fEnergyColorScale->GetValue());
  Update();
}

void JetPlaneEditor::DoEnergyScale()
{
  fM->SetEnergyScale(fEnergyScale->GetValue());
  Update();
}

void JetPlaneEditor::DoOneSelection()
{
  fTwoSelection->SetState(kButtonUp);
  fM->SetOneSelection(fOneSelection->IsOn());
  fM->SetTwoSelection(fTwoSelection->IsOn());
  Update();
}

void JetPlaneEditor::DoTwoSelection()
{
  fOneSelection->SetState(kButtonUp);
  fM->SetOneSelection(fOneSelection->IsOn());
  fM->SetTwoSelection(fTwoSelection->IsOn());
  Update();
}

void JetPlaneEditor::DoStaticDataWindow()
{	
  printf("\n Soon available ... \n");
  if (fgStaticWindow == 0)
    fgStaticWindow = new StaticDataWindow(gClient->GetRoot(), this, 400, 200);

  // call fgStaticWindow->ReadValues(); // like setmodel

  // position relative to the parent's window
  fgStaticWindow->MapWindow();
  fgStaticWindow->RaiseWindow();
  fgStaticWindow->CenterOnParent();
}

/**************************************************************************/

ClassImp(JetPlaneEditor::StaticDataWindow)

JetPlaneEditor::StaticDataWindow::StaticDataWindow(const TGWindow *p, const TGWindow *main,
						   UInt_t w, UInt_t h, UInt_t options) :
  TGTransientFrame(p, main, w, h, options),
  fFrame1(0),
  fOkButton(0),
  fCancelButton(0),
  fL1(0),
  fL2(0),
  fL3(0),
  fL5(0),
  fTab(0),
  fChk1(0),fChk2(0),fChk3(0),fChk4(0),fChk5(0)
{
  // Create a dialog window. A dialog window pops up with respect to its
  // "main" window.

  Connect("CloseWindow()", "JetPlaneEditor::StaticDataWindow", this, "DoClose()");
  DontCallClose(); // to avoid double deletions.

  // use hierarchical cleaning
  SetCleanup(kDeepCleanup);

  fFrame1 = new TGHorizontalFrame(this, 60, 20, kFixedWidth);

  fOkButton = new TGTextButton(fFrame1, "&Ok", 1);
  fOkButton->Connect("Clicked()", "JetPlaneEditor::StaticDataWindow", this, "DoOK()");
  fCancelButton = new TGTextButton(fFrame1, "&Cancel", 2);
  fCancelButton->Connect("Clicked()", "JetPlaneEditor::StaticDataWindow", this, "DoCancel()");

  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);

  fFrame1->AddFrame(fOkButton, fL1);
  fFrame1->AddFrame(fCancelButton, fL1);

  fFrame1->Resize(150, fOkButton->GetDefaultHeight());
  AddFrame(fFrame1, fL2);

  // Tabs for one and two track information

  fTab = new TGTab(this, 300, 300);
  fTab->Connect("Selected(Int_t)", "JetPlaneEditor::StaticDataWindow", this, "DoTab(Int_t)");

  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

  TGCompositeFrame *tf = fTab->AddTab("One TEveTrack/Jet");

  //    fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  //    fF1->AddFrame(new TGTextButton(fF1, "&Test button", 0), fL3);
  //    fF1->AddFrame(fTxt1 = new TGTextEntry(fF1, new TGTextBuffer(100)), fL3);
  //    fF1->AddFrame(fTxt2 = new TGTextEntry(fF1, new TGTextBuffer(100)), fL3);
  //    tf->AddFrame(fF1, fL3);
  //    fTxt1->Resize(150, fTxt1->GetDefaultHeight());
  //    fTxt2->Resize(150, fTxt2->GetDefaultHeight());
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2);
  fF2 = new TGCompositeFrame(tf, 60, 60, kVerticalFrame);
  fF2->AddFrame(fChk1 = new TGCheckButton(fF2, "4-Momentum: {pt, px, py, pz} "), fL1);
  fF2->AddFrame(fChk2 = new TGCheckButton(fF2, "4-Momentum: {pt, Phi, Theta}"), fL1);
  fF2->AddFrame(fChk3 = new TGCheckButton(fF2, "Pseudorapidity: Eta"), fL1);
  fF2->AddFrame(fChk4 = new TGCheckButton(fF2, "Energy: E"), fL1);
  fF2->AddFrame(fChk5 = new TGCheckButton(fF2, "Charge and Mass"), fL1);

  tf = fTab->AddTab("Two Tracks/Jets");

  tf->AddFrame(fF2, fL3);

  //    fBtn1->Connect("Clicked()", "TestDialog", this, "HandleButtons()");
  //    fBtn2->Connect("Clicked()", "TestDialog", this, "HandleButtons()");
  //    fChk1->Connect("Clicked()", "TestDialog", this, "HandleButtons()");
  //    fChk2->Connect("Clicked()", "TestDialog", this, "HandleButtons()");
  //    fRad1->Connect("Clicked()", "TestDialog", this, "HandleButtons()");
  //    fRad2->Connect("Clicked()", "TestDialog", this, "HandleButtons()");


  TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);
  AddFrame(fTab, fL5);

  MapSubwindows();
  Resize();

  SetWindowName("TEveTrack/Jet Common Setup");
}

JetPlaneEditor::StaticDataWindow::~StaticDataWindow()
{
  DeleteWindow();
}

void JetPlaneEditor::StaticDataWindow::DoClose()
{
  UnmapWindow();
}

void JetPlaneEditor::StaticDataWindow::DoOK()
{
  // Read data from widgets, copy to static members of JetPlane

  SendCloseMessage();
}

void JetPlaneEditor::StaticDataWindow::DoCancel()
{
  SendCloseMessage();
}

void JetPlaneEditor::StaticDataWindow::DoTab(Int_t /*id*/)
{
  // printf("Tab item %d activated\n", id);
}


