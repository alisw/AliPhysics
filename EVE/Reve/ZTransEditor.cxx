// $Header$

#include "ZTransEditor.h"
#include <Reve/ZTrans.h>
#include <Reve/RGValuators.h>

#include <TVirtualPad.h>

#include <TGButton.h>


using namespace Reve;

//______________________________________________________________________
// ZTransSubEditor
//

ClassImp(ZTransSubEditor)

ZTransSubEditor::ZTransSubEditor(TGWindow* p) :
  TGVerticalFrame(p),
  fTrans  (0)
{
  // --- Top controls

  fTopHorFrame = new TGHorizontalFrame(this);

  fUseTrans  = new TGCheckButton(fTopHorFrame, "UseTrans");
  fTopHorFrame->AddFrame(fUseTrans, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
  fUseTrans->Connect("Toggled(Bool_t)", "Reve::ZTransSubEditor", this, "DoUseTrans()");
  fEditTrans = new TGCheckButton(fTopHorFrame, "EditTrans");
  fTopHorFrame->AddFrame(fEditTrans, new TGLayoutHints(kLHintsLeft, 2,1,0,0));  
  fEditTrans->Connect("Toggled(Bool_t)"," Reve::ZTransSubEditor", this, "DoEditTrans()");

  AddFrame(fTopHorFrame, new TGLayoutHints(kLHintsTop, 0,0,2,1));


  // --- Trans edit part

  fEditTransFrame = new TGVerticalFrame(this);

  fPos = new RGTriVecValuator(fEditTransFrame, "Pos", 160, 20);
  fPos->SetLabelWidth(9);
  fPos->SetNELength(6);
  fPos->Build(kFALSE, "x", "y", "z");
  fPos->SetLimits(-1e5, 1e5, TGNumberFormat::kNESRealThree);
  fEditTransFrame->AddFrame(fPos, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,0,0));

  fRot = new RGTriVecValuator(fEditTransFrame, "Rot", 160, 20);
  fRot->SetLabelWidth(17);
  fRot->SetNELength(5);
  fRot->Build(kFALSE, "Rz", "RY", "Rx");
  fRot->SetLimits(-360, 360, TGNumberFormat::kNESRealOne);
  fEditTransFrame->AddFrame(fRot, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,0,0));

  fScale = new RGTriVecValuator(fEditTransFrame, "Scale", 160, 20);
  fScale->SetLabelWidth(17);
  fScale->SetNELength(5);
  fScale->Build(kFALSE, "Sx", "Sy", "Sz");
  fScale->SetLimits(1e-2, 1e2, TGNumberFormat::kNESRealTwo);
  fEditTransFrame->AddFrame(fScale, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,0,0));

  fPos  ->Connect("ValueSet()", "Reve::ZTransSubEditor", this, "DoTransChanged()");
  fRot  ->Connect("ValueSet()", "Reve::ZTransSubEditor", this, "DoTransChanged()");
  fScale->Connect("ValueSet()", "Reve::ZTransSubEditor", this, "DoTransChanged()");

  {
    TGHorizontalFrame* hframe = new TGHorizontalFrame(fEditTransFrame);

    fAutoUpdate = new TGCheckButton(hframe, "AutoUpdate");
    hframe->AddFrame(fAutoUpdate, new TGLayoutHints(kLHintsLeft, 1,10,1,1));
    fUpdate = new TGTextButton(hframe, "Update");
    hframe->AddFrame(fUpdate, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 5,5,1,1));
    fUpdate->Connect("Clicked()", "Reve::ZTransSubEditor", this, "TransChanged()");

    fEditTransFrame->AddFrame(hframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,0,0));
  }

  AddFrame(fEditTransFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,1,2));
}

/**************************************************************************/

void ZTransSubEditor::SetDataFromTrans(ZTrans* t)
{
  fTrans = t;

  fUseTrans ->SetState(fTrans->fUseTrans  ? kButtonDown : kButtonUp);
  fEditTrans->SetState(fTrans->fEditTrans ? kButtonDown : kButtonUp);
  if (fTrans->fEditTrans)
    fEditTransFrame->MapWindow();
  else
    fEditTransFrame->UnmapWindow();
  ((TGMainFrame*)fEditTransFrame->GetMainFrame())->Layout();

  fPos->SetValues(fTrans->ArrT());
  Float_t a[3];
  fTrans->GetRotAngles(a);
  a[0] *= TMath::RadToDeg();
  a[1] *= TMath::RadToDeg();
  a[2] *= TMath::RadToDeg();
  fRot->SetValues(a);
  Double_t x, y, z;
  fTrans->GetScale(x, y, z);
  fScale->SetValues(x, y, z);
}

void ZTransSubEditor::SetTransFromData()
{
  Double_t v[3];
  fTrans->UnitTrans();
  fRot->GetValues(v);
  fTrans->SetRotByAngles(v[0]*TMath::DegToRad(), v[1]*TMath::DegToRad(), v[2]*TMath::DegToRad());
  fPos->GetValues(v);
  fTrans->SetPos(v);
  fScale->GetValues(v);
  fTrans->Scale(v[0], v[1], v[2]);
}

/**************************************************************************/

void ZTransSubEditor::UseTrans()
{
  Emit("UseTrans()");
}

void ZTransSubEditor::TransChanged()
{
  SetTransFromData();
  Emit("TransChanged()");
}

/**************************************************************************/

void ZTransSubEditor::DoUseTrans()
{
  fTrans->SetUseTrans(fUseTrans->IsOn());
  UseTrans();
}

void ZTransSubEditor::DoEditTrans()
{
  fTrans->SetEditTrans(fEditTrans->IsOn());
  if (fEditTrans->IsOn())
    fEditTransFrame->MapWindow();
  else
    fEditTransFrame->UnmapWindow();
  ((TGMainFrame*)fEditTransFrame->GetMainFrame())->Layout();
}

void ZTransSubEditor::DoTransChanged()
{
  if (fAutoUpdate->IsOn())
    TransChanged();
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// ZTransEditor
//

ClassImp(ZTransEditor)

ZTransEditor::ZTransEditor(const TGWindow *p, Int_t width, Int_t height,
			   UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
  // Initialize widget pointers to 0
{
  MakeTitle("ZTrans");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::ZTransEditor", this, "DoXYZZ()");
}

ZTransEditor::~ZTransEditor()
{}

/**************************************************************************/

void ZTransEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<ZTrans*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void ZTransEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
