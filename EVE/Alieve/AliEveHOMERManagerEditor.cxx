// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERManagerEditor.h"
#include <Alieve/AliEveHOMERManager.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
//______________________________________________________________________________
// AliEveHOMERManagerEditor
//

ClassImp(AliEveHOMERManagerEditor)

AliEveHOMERManagerEditor::AliEveHOMERManagerEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fButt(0)
  // Initialize widget pointers to 0
{
  MakeTitle("AliEveHOMERManager");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "AliEveHOMERManagerEditor", this, "DoXYZZ()");
  fButt = new TGTextButton(this, "Ooogadooga");
  AddFrame(fButt); //, new TGLayoutHints(...));
  fButt->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "DoButt()");

}

AliEveHOMERManagerEditor::~AliEveHOMERManagerEditor()
{}

/******************************************************************************/

void AliEveHOMERManagerEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveHOMERManager*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/******************************************************************************/

// Implements callback/slot methods

// void AliEveHOMERManagerEditor::DoXYZZ()
// {
//   fM->SetXYZZ(fXYZZ->GetValue());
//   Update();
// }

void AliEveHOMERManagerEditor::DoButt()
{
  fM->CreateHOMERSourcesList();
}
