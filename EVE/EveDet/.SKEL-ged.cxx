// $Header$

#include "CLASS.h"
#include <Alieve/STEM.h>

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
// CLASS
//

ClassImp(CLASS)

CLASS::CLASS(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
  // Initialize widget pointers to 0
{
  MakeTitle("STEM");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::CLASS", this, "DoXYZZ()");
}

CLASS::~CLASS()
{}

/**************************************************************************/

void CLASS::SetModel(TObject* obj)
{
  fM = dynamic_cast<STEM*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void CLASS::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
