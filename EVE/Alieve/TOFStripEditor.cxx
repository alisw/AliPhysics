// $Header$

#include "TOFStripEditor.h"
#include <Alieve/TOFStrip.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
using namespace Alieve;

//______________________________________________________________________
// TOFStripEditor
//

ClassImp(TOFStripEditor)

TOFStripEditor::TOFStripEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
  // Initialize widget pointers to 0
{
  MakeTitle("TOFStrip");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::TOFStripEditor", this, "DoXYZZ()");
}

TOFStripEditor::~TOFStripEditor()
{}

/**************************************************************************/

void TOFStripEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TOFStrip*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void TOFStripEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
