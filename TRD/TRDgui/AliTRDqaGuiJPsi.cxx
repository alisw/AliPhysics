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

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
//
// S. Radomski 
// Uni-Heidelberg
// April 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiJPsi.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiJPsi)

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiJPsi::AliTRDqaGuiJPsi() 
  : fIdx(0),
    fFileName(0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 6; i++) {
    fNameList[i]   = 0x0;
    fCanvasList[i] = 0x0;
    fHistList[i]   = 0x0;
  }

  //strncpy(fFileName,"",256);

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiJPsi::AliTRDqaGuiJPsi(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500),
    fIdx(0),
    fFileName(0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0)
{
  //
  // Main constructor
  //
  
  fIdx = 0;
  
  // steering panel 
  
  SetLayoutManager(new TGVerticalLayout(this));

  fGPanel = new TGHorizontalFrame(this);

  // fGLabel = new TGLabel(fGPanel, "Current Step: ");
  fGPrev = new TGTextButton(fGPanel, "Prev Step");
  fGNext = new TGTextButton(fGPanel, "Next Step");

  fGSelect = new TGComboBox(fGPanel);
  for(int i=0; i<5; i++) fGSelect->AddEntry(Form("Step %d", i), i);
  fGSelect->Resize(100, fGPrev->GetHeight());
  fGSelect->Select(fIdx);
  
  TGLayoutHints *hint = new TGLayoutHints(kLHintsNormal, 5, 5, 5, 5);
  
  // fGPanel->AddFrame(fGLabel, hint);
  fGPanel->AddFrame(fGPrev, hint);
  fGPanel->AddFrame(fGSelect, hint);
  fGPanel->AddFrame(fGNext, hint);

  AddFrame(fGPanel);

  // panel logic
  fGPrev->Connect("Clicked()", "AliTRDqaGuiJPsi", this, "PreviusStep()");
  fGNext->Connect("Clicked()", "AliTRDqaGuiJPsi", this, "NextStep()");
  fGSelect->Connect("Selected(Int_t", "AliTRDqaGuiJPsi", this, "SelectStep(Int_t)");
 
  // histograms
  /**/
  fGCanvas = new TGCompositeFrame(this);
  fGCanvas->SetLayoutManager(new TGMatrixLayout(fGCanvas,2,3,1,1));

  fNameList[0] = "mass";
  fNameList[1] = "nTracksNeg";
  fNameList[2] = "ptNeg";
  fNameList[3] = "pidNeg";
  //fNameList[3] = "ptAngle";
  fNameList[4] = "nTracksPos";
  fNameList[5] = "ptPos";

  for(Int_t i=0; i<6; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(fNameList[i], fGCanvas, 320, 300);
    fGCanvas->AddFrame(fCanvasList[i]);
    fCanvasList[i]->GetCanvas()->SetRightMargin(0.05);
  }
  
  for(Int_t i=0; i<4; i++) {
    fHistList[i] = 0;
  }

  AddFrame(fGCanvas);
  /**/
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiJPsi::SetQAFile(const char *filename) {
  //
  // Ste file with histograms
  // 

  //strncpy(fFileName,filename,256);
  fFileName = filename;

  for(Int_t i=0; i<6; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  const Int_t logy[6] = {0, 1, 1, 1, 1, 1};
  const Int_t logz[6] = {0, 0, 0, 1, 0, 0};
  const char *opt[6]  = {"", "", "", "", "", ""};

  TFile *file = new TFile(filename);
  
  for(Int_t i=0; i<6; i++) {
    fHistList[i] = (TH1D*)gDirectory->Get(Form("%s_%d", fNameList[i], fIdx));
    if (fHistList[i]) fHistList[i]->SetDirectory(0);
    fCanvasList[i]->GetCanvas()->cd();
    if (fHistList[i]) fHistList[i]->Draw(opt[i]);
    gPad->SetLogy(logy[i]);
    gPad->SetLogz(logz[i]);
    fCanvasList[i]->GetCanvas()->Update();
  }

  file->Close();
  delete file;
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiJPsi::SetStep(Int_t idx) {
  //
  // Sets active supermodule 
  //

  fIdx = idx; 
  fGSelect->Select(fIdx, 0); 
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

