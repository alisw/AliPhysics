// $Header$

#include "CLASS.h"
#include <STEM.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// XXCLASS
//
//

ClassImp(XXCLASS)

//______________________________________________________________________
XXCLASS::XXCLASS(const TGWindow *p) :
  TGVerticalFrame(p),
  fM             (0)
{
  // Constructor.
}

//______________________________________________________________________
void XXCLASS::SetModel(STEM* m)
{
  // Set model object.

  fM = m;
}

//______________________________________________________________________
void XXCLASS::Changed()
{
  // Emit Changed signal.

  Emit("Changed()");
}

//______________________________________________________________________
//void XXCLASS::DoABCD()
//{
//   // Set some value from some widget
//   Changed();
//}


//______________________________________________________________________
// CLASS
//
//

ClassImp(CLASS)

//______________________________________________________________________
CLASS::CLASS(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM  (0),
  fSE (0)
{
  // Constructor.

  MakeTitle("STEM");

  fSE = new XXCLASS(this);
  AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
  fSE->Connect("Changed()", "CLASS", this, "Update()");
}

/**************************************************************************/

//______________________________________________________________________
void CLASS::SetModel(TObject* obj)
{
  // Set model object.
  fM = dynamic_cast<STEM*>(obj);
  fSE->SetModel(fM);
}

/**************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________
// void CLASS::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
