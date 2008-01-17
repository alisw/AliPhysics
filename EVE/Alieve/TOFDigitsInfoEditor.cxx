// $Header$

#include "TOFDigitsInfoEditor.h"
#include <Alieve/TOFDigitsInfo.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
using namespace Alieve;

//______________________________________________________________________
// TOFDigitsInfoEditor
//

ClassImp(TOFDigitsInfoEditor)

TOFDigitsInfoEditor::TOFDigitsInfoEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
  // Initialize widget pointers to 0
{
  MakeTitle("TOFDigitsInfo");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::TOFDigitsInfoEditor", this, "DoXYZZ()");
}

TOFDigitsInfoEditor::~TOFDigitsInfoEditor()
{}

/**************************************************************************/

void TOFDigitsInfoEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TOFDigitsInfo*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void TOFDigitsInfoEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
