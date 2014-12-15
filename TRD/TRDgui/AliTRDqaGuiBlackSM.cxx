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

/* $Id: AliTRDqaGuiBlackSM.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of black (non zero zuppresed) events from TRD. 
// It lets display and browse throu histograms created by the class 
// AliTRDqaBlackEvents.
// The class works in cooperation with AliTRDqaGuiMainBlack.
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiBlackSM.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiBlackSM)

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackSM::AliTRDqaGuiBlackSM() 
  : fIdx(0),
    fIdxType(0),
    fSetRangePed(0),
    fSetRangeNoise(0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0),
    fGSelectType(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 5; i++) {
    fNameList[i] = 0x0;
  }
  for (Int_t i = 0; i < 2; i++) {
    fRangePed[i]   = 0.0;
    fRangeNoise[i] = 0.0;
  }
  for (Int_t i = 0; i < 30; i++) {
    fHistList[i]   = 0x0;
    fCanvasList[i] = 0x0;
  }

  strncpy(fFileName,"",256);

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackSM::AliTRDqaGuiBlackSM(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500),
    fIdx(0),
    fIdxType(0),
    fSetRangePed(0),
    fSetRangeNoise(0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0),
    fGSelectType(0)
{
  //
  // Main constructor
  //
  
  fIdx = 0;
  fIdxType = 0;

  // steering panel 
  
  SetLayoutManager(new TGVerticalLayout(this));

  fGPanel = new TGHorizontalFrame(this);

  // fGLabel = new TGLabel(fGPanel, "Current SM: ");
  fGPrev = new TGTextButton(fGPanel, "Prev SM");
  fGNext = new TGTextButton(fGPanel, "Next SM");

  fGSelect = new TGComboBox(fGPanel);
  for(int i=0; i<18; i++) fGSelect->AddEntry(Form("SM %d", i), i);
  fGSelect->Resize(100, (Int_t)(fGPrev->GetHeight()*1.4));
  fGSelect->Select(fIdx,0);

  const char *textTypes[11] = {
    "pedestals", "noise", "peak-peak", "pedestalDist", "noiseDist", "signal", 
    "entries", "entriesDist", "entriesRM", "errorLocMCM", "errorLocADC"
  };
  
  fGSelectType = new TGComboBox(fGPanel);
  for(int i=0; i<11; i++) fGSelectType->AddEntry(textTypes[i], i);
  fGSelectType->Resize(100, (Int_t)(fGPrev->GetHeight()*1.4));
  fGSelectType->Select(fIdxType, 0);

  //fGPlay = new TGTextButton(fGPanel, "PLAY");

  TGLayoutHints *hint = new TGLayoutHints(kLHintsNormal, 5, 5, 5, 5);
  
  // fGPanel->AddFrame(fGLabel, hint);
  fGPanel->AddFrame(fGPrev, hint);
  fGPanel->AddFrame(fGSelect, hint);
  fGPanel->AddFrame(fGNext, hint);
  //fGPanel->AddFrame(fGPlay, hint);
  fGPanel->AddFrame(fGSelectType, hint);

  AddFrame(fGPanel);

  // panel logic
  fGPrev->Connect("Clicked()", "AliTRDqaGuiBlackSM", this, "PreviusSM()");
  fGNext->Connect("Clicked()", "AliTRDqaGuiBlackSM", this, "NextSM()");
  fGSelect->Connect("Selected(Int_t)", "AliTRDqaGuiBlackSM", this, "SelectSM(Int_t)");
  fGSelectType->Connect("Selected(Int_t)", "AliTRDqaGuiBlackSM", this, "SelectType(Int_t)");
  //fGPlay->Connect("Clicked()", "AliTRDqaGuiBlackSM", this, "Play()");

  // histograms
  /**/
  fGCanvas = new TGCompositeFrame(this);
  fGCanvas->SetLayoutManager(new TGMatrixLayout(fGCanvas,6,5,0,0));

  for(Int_t i=0; i<30; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("pos_%d", i), fGCanvas, 200, 120);
    fGCanvas->AddFrame(fCanvasList[i]);
    fCanvasList[i]->GetCanvas()->SetRightMargin(0.05);
    fCanvasList[i]->GetCanvas()->SetTopMargin(0.05);
  }
  
  for(Int_t i=0; i<30; i++) {
    fHistList[i] = 0;
  }

  AddFrame(fGCanvas);
  /**/
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackSM::SetQAFile(const char *filename) {
  //
  // Set the file with histograms
  //
 
  const char *names[11] = {"ped", "noise", "pp","pedDist", "noiseDist", "signal", 
			  "entries", "entriesDist", "entriesRM", "errorLocMCM", "errorLocADC" };
  const char *opt[11] = {"col", "col", "", "", "", "", "col", "", "col", "col", "col"};
  const Int_t kLogy[11] = {0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0};

  strncpy(fFileName,filename,256);
 
  for(int i=0; i<30; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  TFile *file = new TFile(filename);

  for(int i=0; i<30; i++) {
    Int_t index = i + fIdx * 30;
    const char *nn = Form("%s_%d", names[fIdxType], index);
    //printf("%s\n", nn);
    fHistList[i] = (TH1*)file->Get(nn); //Form("%s_$d", names[fIdxType], index));
    
    Int_t s = i/6;
    Int_t l = i%6;
    Int_t pos = (5-l) * 5 + s;

    fCanvasList[pos]->GetCanvas()->cd();
    fCanvasList[pos]->GetCanvas()->SetLogy(kLogy[fIdxType]);

    if (fHistList[i]) fHistList[i]->Draw(opt[fIdxType]);

    if (fHistList[i] && (fIdxType == 6)) {
      fHistList[i]->SetMinimum(0);
      fHistList[i]->SetMaximum(2);
    }

    if ( fHistList[i] && (fIdxType == 0)  && fSetRangePed) {
      fHistList[i]->SetMinimum(fRangePed[0]);
      fHistList[i]->SetMaximum(fRangePed[1]);
    }
    
    if ( fHistList[i] && (fIdxType == 1) && fSetRangeNoise) {
      fHistList[i]->SetMinimum(fRangeNoise[0]);
      fHistList[i]->SetMaximum(fRangeNoise[1]);
    }

    fCanvasList[pos]->GetCanvas()->Update();
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackSM::SetSM(Int_t idx) {
  //
  // Selects active super module
  //
  
  fIdx = idx; 
  fGSelect->Select(fIdx, 0); 
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackSM::SelectType(Int_t idx) {
  //
  // Selects data type 
  //

  fIdxType = idx;
  fGSelectType->Select(fIdxType, 0);
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////
/*
void AliTRDqaGuiBlackSM::Play() {
  
  SetSM(0);
  for(Int_t i=0; i<18; i++) {
    gSystem->Sleep(2 * 1e3);
    NextSM();
  }
}
*/
//////////////////////////////////////////////////////////////////////////////////

