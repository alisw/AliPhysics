// $Header$

#include "ITSModuleStepperEditor.h"
#include <Alieve/ITSModuleStepper.h>
#include <TEveGridStepperEditor.h>
#include <TEveManager.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepperEditor
//

ClassImp(ITSModuleStepperEditor)

ITSModuleStepperEditor::ITSModuleStepperEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fM(0),
  fStepper(0)
{
  MakeTitle("ITSModuleStepper");

  fStepper =  new TEveGridStepperSubEditor(this);
  fStepper->Connect("Changed()", "Alieve::ITSModuleStepperEditor", this, "UpdateStore()");
  AddFrame(fStepper, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));
}

ITSModuleStepperEditor::~ITSModuleStepperEditor()
{}

/**************************************************************************/

void ITSModuleStepperEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<ITSModuleStepper*>(obj);
  fStepper->SetModel(fM->GetStepper());
}

/**************************************************************************/

void ITSModuleStepperEditor::UpdateStore()
{
  fM->Apply();
  Update();
  gEve->Redraw3D(kTRUE);
}
