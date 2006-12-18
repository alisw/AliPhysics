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
//           AliSelectorFrame class
//   The class that deals with the selector tab of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TGFileDialog.h"
#include "TGTextEntry.h"
#include "TGLabel.h"

#include "TObjString.h"

#include "AliTagAnalysisFrame.h"
#include "AliAnalysisGUI.h"

#include "AliSelectorFrame.h"

ClassImp(AliSelectorFrame)

//___________________________________________________________________________
AliSelectorFrame::AliSelectorFrame(const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI *v, AliTagAnalysisFrame* t): 
  TGHorizontalFrame(main, w, h), 
  fVFrame1(0), fVFrame2(0),
  fLabel1(0), fTextSelector(0), fButtonSelect(0), fButtonRun(0),
  fAliAnalysisGUI(v), fTagAnalysisFrame(t) {
  // ctor.
  
  fVFrame1 = new TGVerticalFrame(this, 100, 100);
  AddFrame(fVFrame1, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
  
  fLabel1 = new TGLabel(fVFrame1, new TGString("Selector macro"));
  fVFrame1->AddFrame(fLabel1, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fTextSelector = new TGTextEntry(fVFrame1, new TGTextBuffer(50));
  fVFrame1->AddFrame(fTextSelector, new TGLayoutHints(kLHintsBottom, 5,5,10,5));
  fTextSelector->SetEnabled(false);
  
  fVFrame2 = new TGVerticalFrame(this, 100, 100);
  AddFrame(fVFrame2, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
  
  fButtonSelect = new TGTextButton(fVFrame2, "Select...", 1);
  fVFrame2->AddFrame(fButtonSelect,  new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5,5,5,5));
  fButtonSelect->Connect("Clicked()", "AliSelectorFrame", this, "OnSelect()");
  
  fButtonRun = new TGTextButton(fVFrame2, "Run", 2);
  fVFrame2->AddFrame(fButtonRun,  new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 5,5,5,5));
  
  fButtonRun->Connect("Clicked()", "AliSelectorFrame", this, "OnRun()");
  
  MapWindow();
  Resize();
  MapSubwindows();
}

//___________________________________________________________________________
AliSelectorFrame::AliSelectorFrame(const AliSelectorFrame&):
  fVFrame1(0),
  fVFrame2(0),
  fTextSelector(0),
  fButtonSelect(0),
  fButtonRun(0),
  fAliAnalysisGUI(0),
  fTagAnalysisFrame(0)
{
  //copy constructor
}

//___________________________________________________________________________
void AliSelectorFrame::OnSelect() {
  // When Select button is pressed.

  const char *filetypes[] = { 
    "macro files",    "*.C",
    0,               0 
  };
  
  static TString dir(".");
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  fi.fIniDir    = StrDup(dir);
  
  new TGFileDialog(gClient->GetRoot(), fAliAnalysisGUI, kFDOpen, &fi);
  
  fTextSelector->SetText(fi.fFilename);
}

//___________________________________________________________________________
void AliSelectorFrame::OnRun() {
  // Run the Analysis Selector
  
  TString fname (fTextSelector->GetText());

  TObjArray *a = fname.Tokenize("/");
  TObjString * ostr = (TObjString*)a->At(a->GetEntries()-1);
  
  fTagAnalysisFrame->ProcessSelector(ostr->GetString().Data());
  
  delete a;
}
