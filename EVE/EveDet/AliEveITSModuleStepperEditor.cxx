// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSModuleStepperEditor.h"
#include <EveDet/AliEveITSModuleStepper.h>
#include <TEveGridStepperEditor.h>
#include <TEveManager.h>

#include <TVirtualPad.h>

//==============================================================================
//==============================================================================
// AliEveITSModuleStepperEditor
//==============================================================================

//______________________________________________________________________________
//
// Editor for AliEveITSModuleStepper.

ClassImp(AliEveITSModuleStepperEditor)

AliEveITSModuleStepperEditor::AliEveITSModuleStepperEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fM(0),
  fStepper(0)
{
  // Constructor.

  MakeTitle("AliEveITSModuleStepper");

  fStepper =  new TEveGridStepperSubEditor(this);
  fStepper->Connect("Changed()", "AliEveITSModuleStepperEditor", this, "UpdateStepper()");
  AddFrame(fStepper, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));
}

/******************************************************************************/

void AliEveITSModuleStepperEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEveITSModuleStepper*>(obj);
  fStepper->SetModel(fM->GetStepper());
}

/******************************************************************************/

void AliEveITSModuleStepperEditor::UpdateStepper()
{
  // Update the model object.

  fM->Apply();
  Update();
  gEve->Redraw3D(kTRUE);
}
