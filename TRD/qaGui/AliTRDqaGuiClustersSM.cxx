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

/* $Id: AliTRDqaGuiClustersSM.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters. It lets display and browse throu histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiClustersSM.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiClustersSM)

const Int_t AliTRDqaGuiClustersSM::fgkLogList[4] = {0,0,0,0};

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiClustersSM::AliTRDqaGuiClustersSM() 
  : fIdx(0),
    fFileName(0x0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0),
    fGPlay(0)
{
  //
  // Default constructor
  //
  
  for (Int_t i = 0; i < 4; i++) {
    fNameList[i]   = 0x0;
    fCanvasList[i] = 0x0;
    fHistList[i]   = 0x0;
  }

  //strncpy(fFileName,"",256);

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiClustersSM::AliTRDqaGuiClustersSM(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500),
    fIdx(0),
    fFileName(0x0),
    fGPanel(0),
    fGCanvas(0),
    fGSelect(0),
    fGPrev(0),
    fGNext(0),
    fGPlay(0)
{
  //
  // Main constructor
  //
  
  fIdx = 0;
  
  // steering panel 
  
  SetLayoutManager(new TGVerticalLayout(this));

  fGPanel = new TGHorizontalFrame(this);

  // fGLabel = new TGLabel(fGPanel, "Current SM: ");
  fGPrev = new TGTextButton(fGPanel, "Prev SM");
  fGNext = new TGTextButton(fGPanel, "Next SM");

  fGSelect = new TGComboBox(fGPanel);
  for(int i=0; i<18; i++) fGSelect->AddEntry(Form("SM %d", i), i);
  fGSelect->Resize(100, (Int_t)(1.4*fGPrev->GetHeight()));
  fGSelect->Select(fIdx);

  fGPlay = new TGTextButton(fGPanel, "PLAY");

   TGLayoutHints *hint = new TGLayoutHints(kLHintsNormal, 5, 5, 5, 5);
  
  // fGPanel->AddFrame(fGLabel, hint);
  fGPanel->AddFrame(fGPrev, hint);
  fGPanel->AddFrame(fGSelect, hint);
  fGPanel->AddFrame(fGNext, hint);
  fGPanel->AddFrame(fGPlay, hint);

  AddFrame(fGPanel);

  // panel logic
  fGPrev->Connect("Clicked()", "AliTRDqaGuiClustersSM", this, "PreviusSM()");
  fGNext->Connect("Clicked()", "AliTRDqaGuiClustersSM", this, "NextSM()");
  fGSelect->Connect("Selected(Int_t", "AliTRDqaGuiClustersSM", this, "SelectSM(Int_t)");
  fGPlay->Connect("Clicked()", "AliTRDqaGuiClustersSM", this, "Play()");

  // histograms
  /**/
  fGCanvas = new TGCompositeFrame(this);
  fGCanvas->SetLayoutManager(new TGMatrixLayout(fGCanvas,2,2,1,1));

  fNameList[0] = "sigTimeShape";
  fNameList[1] = "sigTime";
  fNameList[2] = "totalCharge";
  fNameList[3] = "nCls";

  for(Int_t i=0; i<4; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(fNameList[i], fGCanvas, 480, 300);
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

void AliTRDqaGuiClustersSM::SetQAFile(const char *filename) {
  //
  // Ste file with histograms
  //

  //strncpy(fFileName,filename,256);
  fFileName = filename;

  for(Int_t i=0; i<4; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  TFile *file = new TFile(filename);
  file->cd("TRD/RecPoints");
  
  const char *opt[4] = {"", "colz", "", ""};

  for(int i=0; i<4; i++) {
    fHistList[i] = (TH1D*)gDirectory->Get(Form("qaTRD_recPoints_%s_sm%d", fNameList[i], fIdx));
    fCanvasList[i]->GetCanvas()->cd();
    gPad->SetLogy(fgkLogList[i]);
    if (fHistList[i]) fHistList[i]->Draw(opt[i]);
    fCanvasList[i]->GetCanvas()->Update();
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersSM::SetSM(Int_t idx) {
  //
  // Sets active supermodule 
  //

  fIdx = idx; 
  fGSelect->Select(fIdx, 0); 
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersSM::Play() {
  //
  // Loop throught suermodules
  //
  
  SetSM(0);
  for(Int_t i=0; i<18; i++) {
    gSystem->Sleep(1 * 1000);
    NextSM();
  }
}

//////////////////////////////////////////////////////////////////////////////////

