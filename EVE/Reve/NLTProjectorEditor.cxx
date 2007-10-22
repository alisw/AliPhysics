// $Header$

#include "NLTProjectorEditor.h"
#include <Reve/NLTProjector.h>

#include <Reve/RGValuators.h>

#include <TColor.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGComboBox.h>
#include <TGLabel.h>

using namespace Reve;

//______________________________________________________________________
// NLTProjectorEditor
//

ClassImp(NLTProjectorEditor)

NLTProjectorEditor::NLTProjectorEditor(const TGWindow *p,
				       Int_t width, Int_t height,
				       UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),

  fType(0),
  fDistortion(0),
  fFixedRadius(0),
  fCurrentDepth(0),

  fAxisColor(0),
  fSIMode(0),
  fSILevel(0)
{
  MakeTitle("Axis");
  {
    TGHorizontalFrame* hf1 = new TGHorizontalFrame(this);

    TGCompositeFrame *labfr = 
      new TGHorizontalFrame(hf1, 60, 15, kFixedSize);
    TGLabel* l = new TGLabel(labfr, "Color");
    labfr->AddFrame(l, new TGLayoutHints(kLHintsLeft|kLHintsBottom));
    hf1->AddFrame(labfr, new TGLayoutHints(kLHintsLeft|kLHintsBottom));


    fAxisColor = new TGColorSelect(hf1, 0, -1);
    hf1->AddFrame(fAxisColor, new TGLayoutHints(kLHintsLeft, 2, 0, 1, 1));
    fAxisColor->Connect
      ("ColorSelected(Pixel_t)",
       "Reve::NLTProjectorEditor", this, "DoAxisColor(Pixel_t)");

    AddFrame(hf1);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "StepMode");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 6, 1, 2));
    fSIMode = new TGComboBox(f, "Position");
    fSIMode->AddEntry("Value", 1);
    fSIMode->AddEntry("Position", 0);
    fSIMode->GetTextEntry()->SetToolTipText("Set tick-marks on equidistant values/screen position.");
    TGListBox* lb = fSIMode->GetListBox();
    lb->Resize(lb->GetWidth(), 2*18);
    fSIMode->Resize(80, 20);
    fSIMode->Connect("Selected(Int_t)", "Reve::NLTProjectorEditor",
		     this, "DoSplitInfoMode(Int_t)");
    f->AddFrame(fSIMode, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
    AddFrame(f);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "SplitLevel");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 8, 1, 2));
    
    fSILevel = new TGNumberEntry(f, 0, 3, -1,TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax, 0, 3);
    fSILevel->GetNumberEntry()->SetToolTipText("Number of tick-marks TMath::Power(2, level).");
    fSILevel->Connect("ValueSet(Long_t)", "Reve::NLTProjectorEditor", this, "DoSplitInfoLevel()");
    f->AddFrame(fSILevel, new TGLayoutHints(kLHintsTop, 1, 1, 1, 2));
    AddFrame(f, new TGLayoutHints(kLHintsTop, 0, 0, 0, 3) );
  }
  MakeTitle("NLTProjection");

  // Create widgets
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "Type");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 31, 1, 2));
    fType = new TGComboBox(f);
    fType->AddEntry("CFishEye", NLTProjection::PT_CFishEye);
    fType->AddEntry("RhoZ",     NLTProjection::PT_RhoZ);
    TGListBox* lb = fType->GetListBox();
    lb->Resize(lb->GetWidth(), 2*18);
    fType->Resize(80, 20);
    fType->Connect("Selected(Int_t)", "Reve::NLTProjectorEditor",
		   this, "DoType(Int_t)");
    f->AddFrame(fType, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
    AddFrame(f);
  }

  Int_t labelW = 60;
  fDistortion = new RGValuator(this, "Distortion:", 90, 0);
  fDistortion->SetNELength(5);
  fDistortion->SetLabelWidth(labelW);
  fDistortion->Build();
  fDistortion->SetLimits(0, 50, 101, TGNumberFormat::kNESRealTwo);
  fDistortion->Connect("ValueSet(Double_t)", "Reve::NLTProjectorEditor",
		       this, "DoDistortion()");
  AddFrame(fDistortion, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));


  fFixedRadius = new RGValuator(this, "FixedR:", 90, 0);
  fFixedRadius->SetNELength(5);
  fFixedRadius->SetLabelWidth(labelW);
  fFixedRadius->Build();
  fFixedRadius->SetLimits(0, 1000, 101, TGNumberFormat::kNESRealOne);
  fFixedRadius->SetToolTip("Radius not scaled by distotion.");
  fFixedRadius->Connect("ValueSet(Double_t)", "Reve::NLTProjectorEditor",
			this, "DoFixedRadius()");
  AddFrame(fFixedRadius, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));


  fCurrentDepth = new RGValuator(this, "CurrentZ:", 90, 0);
  fCurrentDepth->SetNELength(5);
  fCurrentDepth->SetLabelWidth(labelW);
  fCurrentDepth->Build();
  fCurrentDepth->SetLimits(-300, 300, 601, TGNumberFormat::kNESRealOne);
  fCurrentDepth->SetToolTip("Radius not scaled by distotion.");
  fCurrentDepth->Connect("ValueSet(Double_t)", "Reve::NLTProjectorEditor",
			this, "DoFixedRadius()");
  AddFrame(fCurrentDepth, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

NLTProjectorEditor::~NLTProjectorEditor()
{}

/**************************************************************************/

void NLTProjectorEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<NLTProjector*>(obj);

  fAxisColor->SetColor(TColor::Number2Pixel(fM->GetAxisColor()), kFALSE);
  fSIMode->Select(fM->GetSplitInfoMode(), kFALSE);  
  fSILevel->SetNumber(fM->GetSplitInfoLevel());

  fType->Select(fM->GetProjection()->GetType(), kFALSE);  
  fDistortion->SetValue(1000.0f * fM->GetProjection()->GetDistortion());
  fFixedRadius->SetValue(fM->GetProjection()->GetFixedRadius());
  fCurrentDepth->SetValue(fM->GetCurrentDepth());
}

/**************************************************************************/

void NLTProjectorEditor::DoType(Int_t type)
{
  fM->SetProjection((NLTProjection::PType_e)type, 0.001f * fDistortion->GetValue());
  fM->ProjectChildren();
  Update();
}

void NLTProjectorEditor::DoDistortion()
{
  fM->GetProjection()->SetDistortion(0.001f * fDistortion->GetValue());
  fM->UpdateName();
  fM->ProjectChildren();
  Update();
}

void NLTProjectorEditor::DoFixedRadius()
{
  fM->GetProjection()->SetFixedRadius(fFixedRadius->GetValue());
  fM->ProjectChildren();
  Update();
}

void NLTProjectorEditor::DoCurrentDepth()
{
  fM->GetProjection()->SetFixedRadius(fCurrentDepth->GetValue());
  fM->ProjectChildren();
  Update();
}
/**************************************************************************/

void NLTProjectorEditor::DoSplitInfoMode(Int_t type)
{
  fM->SetSplitInfoMode(type);
  Update();
}

void NLTProjectorEditor::DoSplitInfoLevel()
{
  fM->SetSplitInfoLevel((Int_t)fSILevel->GetNumber());
  Update();
}

void NLTProjectorEditor::DoAxisColor(Pixel_t pixel)
{ 
  fM->SetAxisColor(Color_t(TColor::GetColor(pixel)));
  Update();
}
