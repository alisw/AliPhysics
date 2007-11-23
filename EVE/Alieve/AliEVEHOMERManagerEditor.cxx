// $Header$

#include "AliEVEHOMERManagerEditor.h"
#include <Alieve/AliEVEHOMERManager.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// AliEVEHOMERManagerEditor
//

ClassImp(AliEVEHOMERManagerEditor)

AliEVEHOMERManagerEditor::AliEVEHOMERManagerEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fButt(0)
  // Initialize widget pointers to 0
{
  MakeTitle("AliEVEHOMERManager");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Alieve::AliEVEHOMERManagerEditor", this, "DoXYZZ()");
  fButt = new TGTextButton(this, "Ooogadooga");
  AddFrame(fButt); //, new TGLayoutHints(...));
  fButt->Connect("Clicked()", "AliEVEHOMERManagerEditor", this, "DoButt()");

}

AliEVEHOMERManagerEditor::~AliEVEHOMERManagerEditor()
{}

/**************************************************************************/

void AliEVEHOMERManagerEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEVEHOMERManager*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/**************************************************************************/

// Implements callback/slot methods

// void AliEVEHOMERManagerEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }

void AliEVEHOMERManagerEditor::DoButt()
{
  fM->CreateHOMERSourcesList();
}
