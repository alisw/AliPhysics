// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveEventManagerEditor.h"
#include "AliEveEventManager.h"

#include "TVirtualPad.h"
#include "TColor.h"

#include <TEveGValuators.h>
#include <TGButton.h>
#include <TGTextView.h>
#include <TGLabel.h>

//______________________________________________________________________________
// GUI editor for AliEveEventManager.
//

ClassImp(AliEveEventManagerEditor)

//______________________________________________________________________________
AliEveEventManagerEditor::AliEveEventManagerEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fAutoLoad(0),
  fAutoLoadTime(0),
  fNextEvent(0),
  fPrevEvent(0),
  fLastEvent(0),
  fRefresh(0),
  fEventInfo(0)
{
  // Constructor.

  MakeTitle("AliEveEventManager");

  // Create widgets
  {
    fAutoLoadTime = new TEveGValuator(this, "Autoload time:", 110, 0);
    fAutoLoadTime->SetShowSlider(kFALSE);
    fAutoLoadTime->SetNELength(4);
    fAutoLoadTime->Build();
    fAutoLoadTime->SetLimits(0, 1000);
    fAutoLoadTime->SetToolTip("Automatic event loading time in seconds");
    fAutoLoadTime->Connect("ValueSet(Double_t)",
			   "AliEveEventManagerEditor", this, "DoSetAutoLoadTime()");

    fAutoLoad = new TGCheckButton(fAutoLoadTime,"Autoload");
    fAutoLoad->SetToolTipText("Automatic event loading (online)");
    fAutoLoadTime->AddFrame(fAutoLoad, new TGLayoutHints(kLHintsLeft, 12, 0, 2, 0));
    fAutoLoad->Connect("Toggled(Bool_t)",
		       "AliEveEventManagerEditor", this, "DoSetAutoLoad()");
    AddFrame(fAutoLoadTime);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fPrevEvent = new TGTextButton(f, "Previous Event");
    fPrevEvent->SetWidth(100);
    fPrevEvent->ChangeOptions(fPrevEvent->GetOptions() | kFixedWidth);
    f->AddFrame(fPrevEvent, new TGLayoutHints(kLHintsNormal, 0,4,0,0));
    fPrevEvent->Connect("Clicked()",
			"AliEveEventManagerEditor", this, "DoPrevEvent()");
    fNextEvent = new TGTextButton(f, "Next Event");
    fNextEvent->SetWidth(100);
    fNextEvent->ChangeOptions(fNextEvent->GetOptions() | kFixedWidth);
    f->AddFrame(fNextEvent, new TGLayoutHints(kLHintsNormal, 4,0,0,0));
    fNextEvent->Connect("Clicked()",
			"AliEveEventManagerEditor", this, "DoNextEvent()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fLastEvent = new TGTextButton(f, "Last Event");
    fLastEvent->SetWidth(100);
    fLastEvent->ChangeOptions(fLastEvent->GetOptions() | kFixedWidth);
    f->AddFrame(fLastEvent, new TGLayoutHints(kLHintsNormal, 0,4,0,0));
    fLastEvent->Connect("Clicked()",
			"AliEveEventManagerEditor", this, "DoLastEvent()");
    fRefresh = new TGTextButton(f, "Refresh");
    fRefresh->SetWidth(100);
    fRefresh->ChangeOptions(fRefresh->GetOptions() | kFixedWidth);
    f->AddFrame(fRefresh, new TGLayoutHints(kLHintsNormal, 4,0,0,0));
    fRefresh->Connect("Clicked()",
			"AliEveEventManagerEditor", this, "DoRefresh()");
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 8,8,8,0));
  }

  {
    TGVerticalFrame* f = new TGVerticalFrame(this);

    TGLabel *eventInfoLabel = new TGLabel(f, "Event Information:");
    f->AddFrame(eventInfoLabel, new TGLayoutHints(kLHintsNormal, 0,0,6,2));

    fEventInfo = new TGTextView(f, 200, 300);
    f->AddFrame(fEventInfo, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));

    AddFrame(f, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
  }

}

/******************************************************************************/

//______________________________________________________________________________
void AliEveEventManagerEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveEventManager*>(obj);

  // Set values of widgets
  fAutoLoadTime->SetValue(fM->GetAutoLoadTime());
  fAutoLoadTime->SetEnabled(fM->GetAutoLoad());
  fAutoLoad->SetState(fM->GetAutoLoad() ? kButtonDown : kButtonUp);

  fPrevEvent->SetEnabled(!fM->GetIsOnline()); 

  fEventInfo->LoadBuffer(fM->GetEventInfo());
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
void AliEveEventManagerEditor::DoSetAutoLoad()
{
  // Set the auto-load flag
  //
  fM->SetAutoLoad(fAutoLoad->IsOn());
  Update();
}

//______________________________________________________________________________
void AliEveEventManagerEditor::DoSetAutoLoadTime()
{
  // Set the auto-load time in seconds
  //
  fM->SetAutoLoadTime(fAutoLoadTime->GetValue());
  fM->SetAutoLoad(fAutoLoad->IsOn());
  Update();
}

//______________________________________________________________________________
void AliEveEventManagerEditor::DoPrevEvent()
{
  // Load previous event
  fM->PrevEvent();
}

//______________________________________________________________________________
void AliEveEventManagerEditor::DoNextEvent()
{
  // Load next event
  fM->NextEvent();
}

//______________________________________________________________________________
void AliEveEventManagerEditor::DoLastEvent()
{
  // Load last event
  fM->GotoEvent(-1);
}

//______________________________________________________________________________
void AliEveEventManagerEditor::DoRefresh()
{
  // Refresh / reopen data-files
  Int_t ev = fM->GetEventId();
  fM->Close();
  fM->Open();
  fM->GotoEvent(ev);
}
