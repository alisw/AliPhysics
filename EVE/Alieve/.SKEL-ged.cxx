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

CLASS::CLASS(const TGWindow *p, Int_t id, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("STEM");

  //!!! create the widgets here ...

  // Register the editor.
  TClass *cl = STEM::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

CLASS::~CLASS()
{}

/**************************************************************************/

void CLASS::SetModel(TVirtualPad* pad, TObject* obj, Int_t event)
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(STEM::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<STEM*>(fModel);

  SetActive();
}

/**************************************************************************/
