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
//           AliTagFrame class
//   The class that deals with the event tag tab of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TTimer.h"

#include "AliTagFrame.h"

ClassImp(AliTagFrame)

//___________________________________________________________________________
AliTagFrame::AliTagFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, const char* tagName, Int_t tagId, ETagRangeType range)
  : TGTransientFrame(p, main, w, h, options), fMin(0), fMax(0), fRange(range) {
  //constructor
  SetCleanup(kDeepCleanup);

  SetWindowName(tagName);
/*
   AliTagFrameFunctions tagFunctions [6] = {
     &AliTagFrame::Vx, &AliTagFrame::Vy, &AliTagFrame::Vz,
     &AliTagFrame::Participants, &AliTagFrame::ImpactParameter, &AliTagFrame::PrimaryVertex
   };

   AliTagFrameFunctions tagFunctions [3] = {
     &AliTagFrame::NegMultiplicityRange, &AliTagFrame::VOsRange, &AliTagFrame::NPionRange
   };
*/

  CreateTagName(tagName);
  //   CREATE_TAG_FRAME(this, tagFunctions[tagId])();
  if(fRange == kRangeMinMax){
    CreateTagRange(fVFrame1, fEntry1, "Min");
    CreateTagRange(fVFrame2, fEntry2, "Max");
  }
  else if(fRange == kRangeMin){
    CreateTagRange(fVFrame1, fEntry1, "Min");
  }
  else if(fRange == kRangeMax){
    CreateTagRange(fVFrame1, fEntry1, "Max");
  }
   
  CreateTagButton();

  fButton->Connect("Clicked()", "AliTagFrame", this, "OnClicked()");

  MapSubwindows();
  Resize();
  MapWindow();
  
  gClient->WaitFor(this); 
}

//___________________________________________________________________________
AliTagFrame::~AliTagFrame() {
  // destructor

  //   DeleteWindow();
}

//___________________________________________________________________________
void AliTagFrame::CreateTagName(const char* name) {
  // Creates the Tag Name and insert it at Result Group Frame     

   AddFrame(new TGLabel(this, new TGString(name)), new TGLayoutHints(kLHintsLeft, 5,5,40,5));
}

//___________________________________________________________________________
void AliTagFrame::CreateTagRange(TGVerticalFrame * vFrame, TGNumberEntryField *& entry, const char* name) {
  // Creates the label and entries.
  
  vFrame = new TGVerticalFrame(this);
  AddFrame(vFrame, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));
  
  TGLabel * label2 = new TGLabel(vFrame, new TGString(name));
  vFrame->AddFrame(label2, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5,5,5,5));
  entry = new TGNumberEntryField(vFrame);
  vFrame->AddFrame(entry, new TGLayoutHints(kLHintsBottom, 5,5,5,5));
}

//___________________________________________________________________________
void AliTagFrame::CreateTagButton() {
  // Creates the OK button.
  
  fButton = new TGTextButton(this, "OK", 1);
  AddFrame(fButton, new TGLayoutHints(kLHintsRight, 5,5,35,5));
}

//___________________________________________________________________________
void AliTagFrame::OnClicked() {
  // OnClicked slot.
  
  if(fRange == kRangeMinMax){
    fMin = static_cast<int>(fEntry1->GetNumber());
    fMax = static_cast<int>(fEntry2->GetNumber());
  }
  else if(fRange == kRangeMin){
    fMin = static_cast<int>(fEntry1->GetNumber());
  }
  else if(fRange == kRangeMax){
    fMax = static_cast<int>(fEntry1->GetNumber());
  }
  
  TTimer::SingleShot(1, "AliTagFrame", this, "CloseWindow()");
}

