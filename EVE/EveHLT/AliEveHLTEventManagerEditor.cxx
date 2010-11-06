// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHLTEventManagerEditor.h"
#include "AliEveHLTEventManager.h"

#include <TVirtualPad.h>
#include <TColor.h>
#include <TROOT.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>

//______________________________________________________________________________
// AliEveHLTEventManagerEditor
//

ClassImp(AliEveHLTEventManagerEditor)

AliEveHLTEventManagerEditor::AliEveHLTEventManagerEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  
TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fButtonConnect(NULL),
  fButtonWriteToFile(0),
  fButtonNextEvent(0),
  fButtonNavigateBack(0),
  fButtonNavigateFwd(0),
  fButtonPrintScreens(NULL),
  fButtonEventLoopText(0),
  fButtonUpdateEvents(NULL),
  fButtonEventLoop(0),
  fEventLoopStarted(kFALSE),
  fBufferLoopStarted(kTRUE)
{

  MakeTitle("AliEveHLTEventManager");

  fButtonUpdateEvents = new TGTextButton(this, "  Stop buffer loop  ");
  fButtonUpdateEvents->SetWidth(30);
  AddFrame(fButtonUpdateEvents); //, new TGLayoutHints(...));
  fButtonUpdateEvents->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "PollEvents()");

  fButtonNextEvent = new TGTextButton(this, "   NextEvent  ");
  AddFrame(fButtonNextEvent); //, new TGLayoutHints(...));
  fButtonNextEvent->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "NextEvent()");

  fButtonNavigateBack = new TGTextButton(this, "  Navigate Back  ");
  AddFrame(fButtonNavigateBack); //, new TGLayoutHints(...));
  fButtonNavigateBack->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "NavigateBack()");

  fButtonNavigateFwd = new TGTextButton(this, "  Navigate Fwd  ");
  AddFrame(fButtonNavigateFwd); //, new TGLayoutHints(...));
  fButtonNavigateFwd->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "NavigateFwd()");


  fButtonPrintScreens = new TGTextButton(this, "  Save Viewers  ");
  AddFrame(fButtonPrintScreens); //, new TGLayoutHints(...));
  fButtonPrintScreens->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "PrintScreens()");
  
  fButtonWriteToFile = new TGTextButton(this, "  Save to file  ");
  AddFrame(fButtonWriteToFile); //, new TGLayoutHints(...));
  fButtonWriteToFile->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "WriteBlockListToFile()");

  fButtonConnect = new TGTextButton(this, "    Reconnect   ");
  AddFrame(fButtonConnect); //, new TGLayoutHints(...));
  fButtonConnect->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "ConnectToHLT()");

 
 


  // fBoxTriggerSelector = new TGComboBox(this, "Select Trigger");
  // fBoxTriggerSelector->AddEntry("HLT Global Trigger", 0);
  // fBoxTriggerSelector->AddEntry("Barrel multiplicity trigger", 1);
  // fBoxTriggerSelector->AddEntry("PHOS Geometry trigger", 2);
  // fBoxTriggerSelector->AddEntry("No trigger selection", 3);
  // fBoxTriggerSelector->Connect("Selected(Int_t)","AliEveHLTEventManagerEditor", this, "SetTriggerString(int)");
  // fBoxTriggerSelector->SetWidth(150);
  // fBoxTriggerSelector->SetHeight(25);
  // AddFrame(fBoxTriggerSelector);

  fButtonEventLoopText = new TGTextButton(this, "  Loop Events  ");
  AddFrame(fButtonEventLoopText); //, new TGLayoutHints(...));
  fButtonEventLoopText->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "EventLoop()");

  fButtonEventLoop = new TGPictureButton(this, gClient->GetPicture("$ALICE_ROOT/EVE/hlt-macros/HLT-logo.png"));
  AddFrame(fButtonEventLoop); //, new TGLayoutHints(...));
  fButtonEventLoop->Connect("Clicked()", "AliEveHLTEventManagerEditor", this, "EventLoop()");
}

/******************************************************************************/

void AliEveHLTEventManagerEditor::SetModel(TObject* obj) {
  fM = dynamic_cast<AliEveHLTEventManager*>(obj);
}

/******************************************************************************/

void AliEveHLTEventManagerEditor::ConnectToHLT() {
   // Connects to HOMER sources -> to HLT.
  fM->ConnectEventBuffer();
}

void AliEveHLTEventManagerEditor::NextEvent() {
  // call next event from AliEveHOMERManger
  fM->NextEvent();
}

void AliEveHLTEventManagerEditor::WriteBlockListToFile() {
  fM->SaveEveryThing();
}


void AliEveHLTEventManagerEditor::PrintScreens() {
  //Print screens
  fM->PrintScreens();
}

void AliEveHLTEventManagerEditor::NavigateFwd() {
  // navigate forward
  if ( !fEventLoopStarted ) {
    fM->NavigateFwd();
  }
}

void AliEveHLTEventManagerEditor::NavigateBack() {
  // navigate back
  if ( !fEventLoopStarted ) {
    fM->NavigateBack();
  }
}

void AliEveHLTEventManagerEditor::PollEvents() {
  if(fBufferLoopStarted) {
    fButtonUpdateEvents->SetText(  "  Start buffer loop  ");
    fBufferLoopStarted = kFALSE;
  } else {
    fButtonUpdateEvents->SetText("-  Stop buffer loop - ");
    fBufferLoopStarted = kTRUE;
  }
  
  fM->StartBufferMonitor();
}


void AliEveHLTEventManagerEditor::EventLoop() {
  // Start/stop event loop
  if ( !fEventLoopStarted ) {
    fEventLoopStarted = kTRUE;
    fButtonEventLoopText->SetText(" Stop event Loop ");
    fM->StartLoop();
  } else {
    fM->StopLoop();
    fEventLoopStarted = kFALSE;
    fButtonEventLoopText->SetText(" Loop Events ");
  }
}

void AliEveHLTEventManagerEditor::SetTriggerString(int id) {

  if (id < 0 || id > 4) {
    return;
  }
  
  // TString tsa[4] = {"HLTGlobalTrigger", 
  // 		    "BarrelMultiplicityTrigger", 
  // 		    "PHOSgeomTrigger",
  // 		    "ALL"};
   
 
  // fM->SetTriggerString(tsa[id]);
    
}

