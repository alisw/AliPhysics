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
  fNextEvent(0),
  fEventInfo(0)
{
  // Constructor.

  MakeTitle("AliEveEventManager");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fNextEvent = new TGTextButton(f, "Next Event");
    fNextEvent->SetWidth(100);
    fNextEvent->ChangeOptions(fNextEvent->GetOptions() | kFixedWidth);
    f->AddFrame(fNextEvent, new TGLayoutHints(kLHintsNormal, 4,0,0,0));
    fNextEvent->Connect("Clicked()",
			"AliEveEventManagerEditor", this, "DoNextEvent()");
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

  fEventInfo->LoadBuffer(fM->GetEventInfoVertical());
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveEventManagerEditor::DoNextEvent()
{
  // Load next event

  fM->NextEvent();
}


//==============================================================================
// AliEveEventManagerWindow
//==============================================================================

//______________________________________________________________________________
//
// Horizontal GUI for AliEveEventManager, to be placed in the
// bottom part of ROOT browser.

ClassImp(AliEveEventManagerWindow)

AliEveEventManagerWindow::AliEveEventManagerWindow() :
  TGMainFrame(gClient->GetRoot(), 400, 100, kVerticalFrame),
  fFirstEvent   (0),
  fPrevEvent    (0),
  fNextEvent    (0),
  fLastEvent    (0),
  fRefresh      (0),
  fEventId      (0),
  fInfoLabel    (0),
  fAutoLoad     (0),
  fAutoLoadTime (0),
  fEventInfo    (0)
{
  // Constructor.

  const TString cls("AliEveEventManagerWindow");
  TGTextButton *b = 0;
  {
    Int_t width = 50;

    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    AddFrame(f, new TGLayoutHints(kLHintsExpandX, 0,0,2,2));

    fFirstEvent = b = MkTxtButton(f, "First", width);
    b->Connect("Clicked()", cls, this, "DoFirstEvent()");
    fPrevEvent = b = MkTxtButton(f, "Prev", width);
    b->Connect("Clicked()", cls, this, "DoPrevEvent()");

    fEventId = new TGNumberEntry(f, 0, 5, -1,TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0, 10000);
    f->AddFrame(fEventId, new TGLayoutHints(kLHintsNormal, 10, 5, 0, 0));
    fEventId->Connect("ValueSet(Long_t)", cls, this, "DoSetEvent()");
    fInfoLabel = new TGLabel(f);
    f->AddFrame(fInfoLabel, new TGLayoutHints(kLHintsNormal, 5, 10, 4, 0));

    fNextEvent = b = MkTxtButton(f, "Next", width);
    b->Connect("Clicked()", cls, this, "DoNextEvent()");
    fLastEvent = b = MkTxtButton(f, "Last", width);
    b->Connect("Clicked()", cls, this, "DoLastEvent()");

    MkLabel(f, "||", 0, 8, 8);

    fRefresh = b = MkTxtButton(f, "Refresh", width + 8);
    b->Connect("Clicked()",cls, this, "DoRefresh()");

    MkLabel(f, "||", 0, 8, 8);

    fAutoLoad = new TGCheckButton(f, "Autoload");
    f->AddFrame(fAutoLoad, new TGLayoutHints(kLHintsLeft, 0, 4, 3, 0));
    fAutoLoad->SetToolTipText("Automatic event loading.");
    fAutoLoad->Connect("Toggled(Bool_t)", cls, this, "DoSetAutoLoad()");

    fAutoLoadTime = new TEveGValuator(f, "Time: ", 110, 0);
    f->AddFrame(fAutoLoadTime);    
    fAutoLoadTime->SetShowSlider(kFALSE);
    fAutoLoadTime->SetNELength(4);
    fAutoLoadTime->Build();
    fAutoLoadTime->SetLimits(0, 1000);
    fAutoLoadTime->SetToolTip("Automatic event loading time in seconds.");
    fAutoLoadTime->Connect("ValueSet(Double_t)", cls, this, "DoSetAutoLoadTime()");
  }

  fEventInfo = new TGTextView(this, 400, 600);
  AddFrame(fEventInfo, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));

  gAliEveEvent->Connect("NewEventLoaded()", cls, this, "Update()");

  SetCleanup(kDeepCleanup);
  Layout();
  MapSubwindows();
  MapWindow();
}

//______________________________________________________________________________
AliEveEventManagerWindow::~AliEveEventManagerWindow()
{
  // Destructor.

  gAliEveEvent->Disconnect("NewEventLoaded()", this);
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoFirstEvent()
{
  // Load previous event
  gAliEveEvent->GotoEvent(0);
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoPrevEvent()
{
  // Load previous event
  gAliEveEvent->PrevEvent();
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoNextEvent()
{
  // Load next event
  gAliEveEvent->NextEvent();
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoLastEvent()
{
  // Load previous event
  gAliEveEvent->GotoEvent(-1);
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoSetEvent()
{
  // Set current event
  gAliEveEvent->GotoEvent((Int_t) fEventId->GetNumber());
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoRefresh()
{
  // Refresh event status.

  Int_t ev = gAliEveEvent->GetEventId();
  gAliEveEvent->Close();
  gAliEveEvent->Open();
  gAliEveEvent->GotoEvent(ev);
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoSetAutoLoad()
{
  // Set the auto-load flag

  gAliEveEvent->SetAutoLoad(fAutoLoad->IsOn());
  Update();
}

//______________________________________________________________________________
void AliEveEventManagerWindow::DoSetAutoLoadTime()
{
  // Set the auto-load time in seconds

  gAliEveEvent->SetAutoLoadTime(fAutoLoadTime->GetValue());
}

//______________________________________________________________________________
void AliEveEventManagerWindow::Update()
{
  // Update current event, number of available events.

  Bool_t autoLoad = gAliEveEvent->GetAutoLoad();
  Bool_t extCtrl  = gAliEveEvent->IsUnderExternalControl();
  Bool_t evNavOn  = !autoLoad && !extCtrl;

  fFirstEvent->SetEnabled(evNavOn);
  fPrevEvent ->SetEnabled(evNavOn);
  fLastEvent ->SetEnabled(evNavOn);
  fNextEvent ->SetEnabled(!autoLoad);
  fRefresh   ->SetEnabled(evNavOn);

  fEventId->SetNumber(gAliEveEvent->GetEventId());
  fEventId->SetState(evNavOn);
  fInfoLabel->SetText(Form("/ %d", gAliEveEvent->GetMaxEventId()));

  // fAutoLoadTime->SetEnabled(gAliEveEvent->GetAutoLoad());
  fAutoLoad->SetState(gAliEveEvent->GetAutoLoad() ? kButtonDown : kButtonUp);
  fAutoLoadTime->SetValue(gAliEveEvent->GetAutoLoadTime());

  fEventInfo->LoadBuffer(gAliEveEvent->GetEventInfoHorizontal());

  Layout();
}

//------------------------------------------------------------------------------
// Protected methods
//------------------------------------------------------------------------------

//______________________________________________________________________________
TGTextButton* AliEveEventManagerWindow::MkTxtButton(TGCompositeFrame* p,
						    const char* txt, Int_t width,
						    Int_t lo, Int_t ro, Int_t to, Int_t bo)
{
  // Create a standard button.
  // If width is not zero, the fixed-width flag is set.

  TGTextButton* b = new TGTextButton(p, txt);
  if (width > 0) {
    b->SetWidth(width);
    b->ChangeOptions(b->GetOptions() | kFixedWidth);
  }
  p->AddFrame(b, new TGLayoutHints(kLHintsNormal, lo,ro,to,bo));
  return b;
}

//______________________________________________________________________________
TGLabel* AliEveEventManagerWindow::MkLabel(TGCompositeFrame* p,
					   const char* txt, Int_t width,
					   Int_t lo, Int_t ro, Int_t to, Int_t bo)
{
  // Create a standard button.
  // If width is not zero, the fixed-width flag is set.

  TGLabel* l = new TGLabel(p, txt);
  if (width > 0) {
    l->SetWidth(width);
    l->ChangeOptions(l->GetOptions() | kFixedWidth);
  }
  p->AddFrame(l, new TGLayoutHints(kLHintsNormal, lo,ro,to,bo));
  return l;
}
