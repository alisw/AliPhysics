/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveVZEROModuleEditor.h"
#include <EveDet/AliEveVZEROModule.h>

#include <AliVZERORawStream.h>

#include <TEveGValuators.h>
#include <TGSlider.h>

//______________________________________________________________________________
//
// Editor for AliEveVZEROModule.

ClassImp(AliEveVZEROModuleEditor)

AliEveVZEROModuleEditor::AliEveVZEROModuleEditor(const TGWindow *p,
                                       Int_t width, Int_t height,
                                       UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fSampleIndex(NULL)
{
  // Constructor.

  MakeTitle("AliEveVZEROModule");

  fSampleIndex = new TEveGValuator(this,"Sample", 200, 0);
  fSampleIndex->SetNELength(4);
  fSampleIndex->SetLabelWidth(60);
  fSampleIndex->Build();
  fSampleIndex->GetSlider()->SetWidth(120);
  fSampleIndex->SetLimits(0, AliVZERORawStream::kNEvOfInt-1, AliVZERORawStream::kNEvOfInt, TGNumberFormat::kNESInteger);
  fSampleIndex->Connect("ValueSet(Double_t)",
		 "AliEveVZEROModuleEditor", this, "DoSampleIndex()");
  AddFrame(fSampleIndex, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  /*
  fSampleIndex = new TGNumberEntry(this,
				   AliVZERORawStream::kNEvOfInt/2, 3, -1,
				   TGNumberFormat::kNESInteger,
				   TGNumberFormat::kNEANonNegative,
				   TGNumberFormat::kNELLimitMinMax,
				   0,AliVZERORawStream::kNEvOfInt);
  AddFrame(fSampleIndex, new TGLayoutHints(kLHintsNormal, 10, 2, 0, 0));
  fSampleIndex->Connect("ValueSet(Double_t)",
			"AliEveVZEROModuleEditor", this, "DoSampleIndex()");
  fSampleIndex->SetText("ADC sample index (between 0 and 21)");
  */
}

/******************************************************************************/

void AliEveVZEROModuleEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveVZEROModule*>(obj);

  fSampleIndex->SetValue(fM->GetSampleIndex());
}

/******************************************************************************/

void AliEveVZEROModuleEditor::DoSampleIndex()
{
  // Slot for SampleIndex.

  fM->SetSampleIndex((Int_t)fSampleIndex->GetValue());
  Update();
}
