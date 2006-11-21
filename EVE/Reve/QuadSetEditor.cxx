// $Header$

#include "QuadSetEditor.h"
#include <Reve/QuadSet.h>

#include <Reve/RGValuators.h>
#include <Reve/ZTransEditor.h>
#include <Reve/RGBAPaletteEditor.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// QuadSetEditor
//

ClassImp(QuadSetEditor)

QuadSetEditor::QuadSetEditor(const TGWindow *p, Int_t width, Int_t height,
			     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fHMTrans   (0),
  fPalette   (0)
  // Initialize widget pointers to 0
{
  MakeTitle("Transformation matrix");

  fHMTrans = new ZTransSubEditor(this);
  fHMTrans->Connect("UseTrans()",     "Reve::QuadSetEditor", this, "Update()");
  fHMTrans->Connect("TransChanged()", "Reve::QuadSetEditor", this, "Update()");
  AddFrame(fHMTrans, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));

  MakeTitle("Palette controls");

  fPalette = new RGBAPaletteSubEditor(this);
  fPalette->Connect("Changed", "Reve::QuadSetEditor", this, "Update()");
  AddFrame(fPalette, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));

  MakeTitle("QuadSet");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::QuadSetEditor", this, "DoXYZZ()");
}

QuadSetEditor::~QuadSetEditor()
{}

/**************************************************************************/

void QuadSetEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<QuadSet*>(obj);

  fHMTrans->SetDataFromTrans(&fM->fHMTrans);

  if (fM->fValueIsColor || fM->fPalette == 0) {
    fPalette->UnmapWindow();
  } else {
    fPalette->SetModel(fM->fPalette);
    fPalette->MapWindow();
  }

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void QuadSetEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }
