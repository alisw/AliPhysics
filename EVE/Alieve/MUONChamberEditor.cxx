#include "MUONChamberEditor.h"

#include <Alieve/MUONChamber.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONChamberEditor
//

ClassImp(MUONChamberEditor)

//______________________________________________________________________
MUONChamberEditor::MUONChamberEditor(const TGWindow *p,
                                     Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
{
  //
  // constructor
  //

  MakeTitle("MUONChamber");

}

//______________________________________________________________________
MUONChamberEditor::~MUONChamberEditor()
{
  //
  // destructor
  //

}

//______________________________________________________________________
void MUONChamberEditor::SetModel(TObject* obj)
{
  //
  // ...
  //

  fM = dynamic_cast<MUONChamber*>(obj);

}
