// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERManagerEditor.h"
#include "AliEveHOMERManager.h"

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
// AliEveHOMERManagerEditor
//

ClassImp(AliEveHOMERManagerEditor)

AliEveHOMERManagerEditor::AliEveHOMERManagerEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  
TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fButtonConnect(NULL),
  fButtonWriteToFile(0),
  fButtonNextEvent(0),
  fButtonNavigateBack(0),
  fButtonNavigateFwd(0),
  fButtonPrintScreens(NULL),
  fBoxTriggerSelector(0),
//  fBoxEventLoopSpeed(0),
  fButtonEventLoopText(0),
  fButtonEventLoop(0),
 fEventLoopStarted(kFALSE) 
{
  
  MakeTitle("AliEveHOMERManager");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "AliEveHOMERManagerEditor", this, "DoXYZZ()");

  fButtonConnect = new TGTextButton(this, " Reconnect ");
  AddFrame(fButtonConnect); //, new TGLayoutHints(...));
  fButtonConnect->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "ConnectToHLT()");

  fButtonWriteToFile = new TGTextButton(this, " Write to file  ");
  AddFrame(fButtonWriteToFile); //, new TGLayoutHints(...));
  fButtonWriteToFile->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "WriteBlockListToFile()");


  fButtonNextEvent = new TGTextButton(this, "  NextEvent  ");
  AddFrame(fButtonNextEvent); //, new TGLayoutHints(...));
  fButtonNextEvent->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "NextEvent()");

  fButtonNavigateBack = new TGTextButton(this, "  Navigate Back  ");
  AddFrame(fButtonNavigateBack); //, new TGLayoutHints(...));
  fButtonNavigateBack->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "NavigateBack()");

  fButtonNavigateFwd = new TGTextButton(this, "  Navigate Fwd  ");
  AddFrame(fButtonNavigateFwd); //, new TGLayoutHints(...));
  fButtonNavigateFwd->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "NavigateFwd()");


  fButtonPrintScreens = new TGTextButton(this, "  Save Viewers  ");
  AddFrame(fButtonPrintScreens); //, new TGLayoutHints(...));
  fButtonPrintScreens->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "PrintScreens()");


  fBoxTriggerSelector = new TGComboBox(this, "Select Trigger");
  fBoxTriggerSelector->AddEntry("HLT Global Trigger", 0);
  fBoxTriggerSelector->AddEntry("Barrel multiplicity trigger", 1);
  fBoxTriggerSelector->AddEntry("PHOS Geometry trigger", 2);
  fBoxTriggerSelector->AddEntry("No trigger selection", 3);
  fBoxTriggerSelector->Connect("Selected(Int_t)","AliEveHOMERManagerEditor", this, "SetTriggerString(int)");
  fBoxTriggerSelector->SetWidth(150);
  fBoxTriggerSelector->SetHeight(25);
  AddFrame(fBoxTriggerSelector);

  fButtonEventLoopText = new TGTextButton(this, "  Loop Events  ");
  AddFrame(fButtonEventLoopText); //, new TGLayoutHints(...));
  fButtonEventLoopText->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "EventLoop()");

  fButtonEventLoop = new TGPictureButton(this, gClient->GetPicture("$ALICE_ROOT/EVE/hlt-macros/HLT-logo.png"));
  AddFrame(fButtonEventLoop); //, new TGLayoutHints(...));
  fButtonEventLoop->Connect("Clicked()", "AliEveHOMERManagerEditor", this, "EventLoop()");
}

/******************************************************************************/

void AliEveHOMERManagerEditor::SetModel(TObject* obj) {
  fM = dynamic_cast<AliEveHOMERManager*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());
}

/******************************************************************************/

void AliEveHOMERManagerEditor::ConnectToHLT() {
   // Connects to HOMER sources -> to HLT.
  fM->ReConnectHOMER();
}

void AliEveHOMERManagerEditor::NextEvent() {
  // call next event from AliEveHOMERManger
  fM->NextHOMEREvent();
}

void AliEveHOMERManagerEditor::WriteBlockListToFile() {
  gROOT->ProcessLineFast("writeToFile();");

}


void AliEveHOMERManagerEditor::PrintScreens() {
  //Print screens
  fM->PrintScreens();
}

void AliEveHOMERManagerEditor::NavigateFwd() {
  // navigate forward
  
  if ( !fEventLoopStarted ) {
    if ( fM->NavigateEventBufferFwd() == -1 )
      return;

    fM->ProcessEvent();
    //    gROOT->ProcessLineFast("processEvent();");
  }
  return;
}

void AliEveHOMERManagerEditor::NavigateBack() {
  // navigate back

  if ( !fEventLoopStarted ) {
    if ( fM->NavigateEventBufferBack() == -1 )
      return;
    
    fM->ProcessEvent();
    //gROOT->ProcessLineFast("processEvent();");
  }
  return;
}

void AliEveHOMERManagerEditor::EventLoop() {
  // Start/stop event loop
  if ( !fEventLoopStarted ) {
    fEventLoopStarted = kTRUE;
    fButtonEventLoopText->SetText(" Stop Loop ");
    //fM->SetEventLoopStarted(kTRUE);
    fM->StartLoop();
    
    //    gROOT->ProcessLineFast("loopEvent();");
  
  } else {
    
    //gROOT->ProcessLineFast("stopLoopEvent();");
    fM->StopLoop();
    fEventLoopStarted = kFALSE;
    //fM->SetEventLoopStarted(kFALSE);
    fButtonEventLoopText->SetText(" Loop Events ");
  }
}

void AliEveHOMERManagerEditor::SetTriggerString(int id) {

  if (id < 0 || id > 4) {
    return;
  }
  
  TString tsa[4] = {"HLTGlobalTrigger", 
		    "BarrelMultiplicityTrigger", 
		    "PHOSgeomTrigger",
		    "ALL"};
   
 
  fM->SetTriggerString(tsa[id]);
    
}

