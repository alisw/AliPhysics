/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           AliLoginFrame class
//   The class that deals with the login frame of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TGLabel.h"
#include "TGTextBuffer.h"
#include "TGTextEntry.h"
#include "TGridResult.h"
#include "TTimer.h"
#include "TGMsgBox.h"

#include "TGrid.h"

#include "AliAnalysisGUI.h"
#include "AliLoginFrame.h"

ClassImp(AliLoginFrame)

//___________________________________________________________________________
AliLoginFrame::AliLoginFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, UInt_t options) : 
  TGTransientFrame(p, main, w, h, options),
  fLabel1(0), fLabel2(0),
  fTextServer(0), fTextUsername(0),
  fButtonLogIn(0), fButtonCancel(0),
  fVFrame1(0),
  fHFrame1(0), fHFrame2(0), fHFrame3(0) {
  // Constructor.

  SetWindowName("Login");

  // create the Vertical Frame for the text fields and buttons
  fVFrame1 = new TGVerticalFrame(this, 400, 150, kFixedWidth);

  // create the Horizontal Frame 1 
  fHFrame1 = new TGHorizontalFrame(this, 400, 50, kFixedWidth);

  fLabel1 = new TGLabel(fHFrame1, new TGString("Server"));
  fTextServer = new TGTextEntry(fHFrame1, new TGTextBuffer(35));
  fTextServer->SetText("alien://");
  
  fHFrame1->AddFrame(fLabel1, 
		     new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));
  fHFrame1->AddFrame(fTextServer, 
		     new TGLayoutHints(kLHintsRight | kLHintsTop, 5, 5, 5, 5));
  
  fVFrame1->AddFrame(fHFrame1,new TGLayoutHints(kLHintsTop));
  
  // Create the Horizontal for the buttons
  fHFrame3 = new TGHorizontalFrame(this, 400, 50, kFixedWidth);
  
  fButtonLogIn = new TGTextButton(fHFrame3, "Log In", 1);
  fButtonLogIn->Connect("Clicked()", "AliLoginFrame", this, "DoLogIn()");
  
  fButtonCancel = new TGTextButton(fHFrame3, "Cancel", 2);
  fButtonCancel->Connect("Clicked()", "AliLoginFrame", this, "DoCancel()");
  
  fHFrame3->AddFrame(fButtonLogIn,new TGLayoutHints
		     (kLHintsCenterY | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
  fHFrame3->AddFrame(fButtonCancel,new TGLayoutHints
		     (kLHintsCenterY | kLHintsRight | kLHintsExpandX, 2, 2, 2, 2));
  
  fHFrame3->Resize();
  
  fVFrame1->AddFrame(fHFrame3, new TGLayoutHints(kLHintsBottom)); 
  
  AddFrame(fVFrame1, new TGLayoutHints(kLHintsTop)); 
  
  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
AliLoginFrame::~AliLoginFrame() {
  // Destructor.
  
  DeleteWindow();
  SetCleanup(kDeepCleanup);
}

//___________________________________________________________________________
void AliLoginFrame::DoLogIn() {
  // When LogIn Button is pressed.
  
  AliAnalysisGUI * fAnalysisGUI = dynamic_cast<AliAnalysisGUI*> (const_cast<TGWindow*> (GetMain()));
  
  fAnalysisGUI->LogIn(fTextServer->GetText());
  
  if(gGrid) {
    TGridResult* result = gGrid->Command("motd",kFALSE,0);
    int i=0;
    const char* line="";
    TString * msg = new TString();
    while ((line=result->GetKey(i++,""))) {      
      msg->Append(line);
      msg->Append("\n");
    }
    new TGMsgBox(gClient->GetRoot(), this, "Welcome", msg->Data(), 0, kMBOk);
    TTimer::SingleShot(1, "AliLoginFrame", this, "CloseWindow()");
  }
}

//___________________________________________________________________________
void AliLoginFrame::DoCancel() {
  // When Cancel button is pressed.

  TTimer::SingleShot(150, "AliLoginFrame", this, "CloseWindow()");
}
