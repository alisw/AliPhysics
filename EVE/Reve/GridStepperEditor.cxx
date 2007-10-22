// $Header$

#include "GridStepperEditor.h"
#include <Reve/GridStepper.h>
#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGSlider.h>
#include <TGButton.h>
#include <TGNumberEntry.h>

using namespace Reve;

//______________________________________________________________________
// GridStepperSubEditor
//
//

ClassImp(GridStepperSubEditor)

//______________________________________________________________________
GridStepperSubEditor::GridStepperSubEditor(const TGWindow *p) :
  TGVerticalFrame(p),
  fM             (0),

  fNx            (0), 
  fNy            (0),
  fNz            (0),
  fDx            (0), 
  fDy            (0),
  fDz            (0)
{
  Int_t labelW = 15;

  TGHorizontalFrame* VF = new TGHorizontalFrame(this);
    
  {
    TGGroupFrame* f = new TGGroupFrame(VF, "NumRows", kVerticalFrame);
    f->SetWidth(30);
    VF->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

    fNx = new RGValuator(f,"X:", 200, 0);
    fNx->SetNELength(3);
    fNx->SetLabelWidth(labelW);
    fNx->SetShowSlider(kFALSE);
    fNx->Build();
    fNx->SetLimits(1, 15);
    fNx->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoNs()");
    f->AddFrame(fNx, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    fNy = new RGValuator(f,"Y:", 200, 0);
    fNy->SetNELength(3);
    fNy->SetLabelWidth(labelW);
    fNy->SetShowSlider(kFALSE);
    fNy->Build();
    fNy->SetLimits(1, 15);
    fNy->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoNs()");
    f->AddFrame(fNy, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    fNz = new RGValuator(f,"Z:", 200, 0);
    fNz->SetNELength(3);
    fNz->SetLabelWidth(labelW);
    fNz->SetShowSlider(kFALSE);
    fNz->Build();
    fNz->SetLimits(1, 15);
    fNz->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoNs()");
    f->AddFrame(fNz, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    //AddFrame(f, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
  }
  {
    TGGroupFrame* f = new TGGroupFrame(VF, "Step", kVerticalFrame);
    f->SetWidth(130);
    VF->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

    fDx = new RGValuator(f,"X:", 200, 0);
    fDx->SetNELength(5);
    fDx->SetLabelWidth(labelW);
    fDx->SetShowSlider(kFALSE);
    fDx->Build();
    fDx->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
    fDx->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoDs()");
    f->AddFrame(fDx, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    fDy = new RGValuator(f,"Y:", 200, 0);
    fDy->SetNELength(5);
    fDy->SetLabelWidth(labelW);
    fDy->SetShowSlider(kFALSE);
    fDy->Build();
    fDy->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
    fDy->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoDs()");
    f->AddFrame(fDy, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    fDz = new RGValuator(f,"Z:", 200, 0);
    fDz->SetNELength(5);
    fDz->SetLabelWidth(labelW);
    fDz->SetShowSlider(kFALSE);
    fDz->Build();
    fDz->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
    fDz->Connect("ValueSet(Double_t)",
		 "Reve::GridStepperSubEditor", this, "DoDs()");
    f->AddFrame(fDz, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    //AddFrame(f, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
  }
  AddFrame(VF, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
}

//______________________________________________________________________
void GridStepperSubEditor::SetModel(GridStepper* m)
{
  // Set model object.
  fM = m;

  fNx->SetValue(m->Nx);
  fNy->SetValue(m->Ny);
  fNz->SetValue(m->Nz);

  fDx->SetValue(m->Dx);
  fDy->SetValue(m->Dy);
  fDz->SetValue(m->Dz);
}

//______________________________________________________________________
void GridStepperSubEditor::Changed()
{
  // Emit Changed signal.

  Emit("Changed()");
}

//______________________________________________________________________
void GridStepperSubEditor::DoNs()
{
  fM->SetNs((Int_t)fNx->GetValue(), (Int_t)fNy->GetValue(), (Int_t)fNz->GetValue()); 
  Changed();
}

//______________________________________________________________________
void GridStepperSubEditor::DoDs()
{
   // Set some value from some widget
  fM->SetDs(fDx->GetValue(), fDy->GetValue(), fDz->GetValue()); 
  Changed();
}

//______________________________________________________________________
// GridStepperEditor
//
//

ClassImp(GridStepperEditor)

//______________________________________________________________________
GridStepperEditor::GridStepperEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM  (0),
  fSE (0)
{
  // Constructor.

  MakeTitle("GridStepper");

  fSE = new GridStepperSubEditor(this);
  AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
  fSE->Connect("Changed()", "GridStepperEditor", this, "Update()");
}

/**************************************************************************/

//______________________________________________________________________
void GridStepperEditor::SetModel(TObject* obj)
{
  // Set model object.
  fM = dynamic_cast<GridStepper*>(obj);
  fSE->SetModel(fM);
}
