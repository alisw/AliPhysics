// $Header$

#include "PMDModuleEditor.h"
#include <Alieve/PMDModule.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// PMDModuleEditor
//

ClassImp(PMDModuleEditor)

PMDModuleEditor::PMDModuleEditor(const TGWindow *p, Int_t width, Int_t height,
				 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoLabel(0)
  // Initialize widget pointers to 0
{
  MakeTitle("PMDModule");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::PMDModuleEditor", this, "DoXYZZ()");

  fInfoLabel = new TGLabel(this);
  fInfoLabel->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
					 8, 0, 2, 0));
}

PMDModuleEditor::~PMDModuleEditor()
{}

/**************************************************************************/

void PMDModuleEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<PMDModule*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
  fInfoLabel->SetText(Form("Pads hit: %d", fM->GetNPads()));
}

/**************************************************************************/

// Implements callback/slot methods

// void PMDModuleEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
