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

/* $Id: AliTRDqaGuiClustersStack.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters on the stack by stack basis. 
// It lets display and browse throu histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiClustersStack.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include "TPaveText.h"
#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiClustersStack)

const Int_t AliTRDqaGuiClustersStack::fgknSM = 18;
const Int_t AliTRDqaGuiClustersStack::fgknStack = 5;
const Int_t AliTRDqaGuiClustersStack::fgknCh = 6;


AliTRDqaGuiClustersStack::AliTRDqaGuiClustersStack() 
  : fIdxSM (0),
    fIdxStack (0),
    fView (0),
    fFileName (0x0),
    fGPanel (0),
    fGCanvas (0),
    fGSelectSM (0),
    fGSelectStack (0),
    fGSelectView (0),
    fGPrevSM (0),
    fGPrevStack (0),
    fGNextSM (0),
    fGNextStack (0),
    fGPlay (0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 6; i++) {
    fCanvasList[i] = 0x0;
    fHistList[i]   = 0x0;
  }

  //strncpy(fFileName,"",256);

}

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiClustersStack::AliTRDqaGuiClustersStack(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500),
    fIdxSM (0),
    fIdxStack (0),
    fView (0),
    fFileName (0x0),
    fGPanel (0),
    fGCanvas (0),
    fGSelectSM (0),
    fGSelectStack (0),
    fGSelectView (0),
    fGPrevSM (0),
    fGPrevStack (0),
    fGNextSM (0),
    fGNextStack (0),
    fGPlay (0)
{
  //
  // Main constructor
  //
  
  fIdxSM = 0;
  fIdxStack = 0;

  // steering panel 
  
  SetLayoutManager(new TGVerticalLayout(this));

  fGPanel = new TGHorizontalFrame(this);

  // fGLabel = new TGLabel(fGPanel, "Current Stack: ");
  fGPrevSM = new TGTextButton(fGPanel, "Prev SM");
  fGPrevStack = new TGTextButton(fGPanel, "Prev Stack");

  fGNextSM = new TGTextButton(fGPanel, "Next SM");
  fGNextStack = new TGTextButton(fGPanel, "Next Stack");

  fGSelectSM = new TGComboBox(fGPanel);
  for(int i=0; i<fgknSM; i++) fGSelectSM->AddEntry(Form("SM %d", i), i);
  fGSelectSM->Resize(100, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectSM->Select(fIdxSM);

  fGSelectStack = new TGComboBox(fGPanel);
  for(int i=0; i<fgknStack; i++) fGSelectStack->AddEntry(Form("Stack %d", i), i);
  fGSelectStack->Resize(100, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectStack->Select(fIdxStack);

  fGPlay = new TGTextButton(fGPanel, "PLAY");

  fGSelectView = new TGComboBox(fGPanel);
  fGSelectView->AddEntry("amplitude",0);
  fGSelectView->AddEntry("time -- signal MPV", 1);
  fGSelectView->AddEntry("time -- total charge", 2);
  fGSelectView->AddEntry("time -- nClusters", 3);
  fGSelectView->Resize(150, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectView->Select(0);

  TGLayoutHints *hint = new TGLayoutHints(kLHintsNormal, 5, 5, 5, 5);
  
  // fGPanel->AddFrame(fGLabel, hint);
  fGPanel->AddFrame(fGPrevSM, hint);
  fGPanel->AddFrame(fGPrevStack, hint);

  fGPanel->AddFrame(fGSelectSM, hint);
  fGPanel->AddFrame(fGSelectStack, hint);

  fGPanel->AddFrame(fGNextStack, hint);
  fGPanel->AddFrame(fGNextSM, hint);

  fGPanel->AddFrame(fGPlay, hint);

  fGPanel->AddFrame(fGSelectView, hint);
  

  AddFrame(fGPanel);

  // panel logic
  fGPrevStack->Connect("Clicked()", "AliTRDqaGuiClustersStack", this, "PreviusStack()");
  fGNextStack->Connect("Clicked()", "AliTRDqaGuiClustersStack", this, "NextStack()");
  fGPrevSM->Connect("Clicked()", "AliTRDqaGuiClustersStack", this, "PreviusSM()");
  fGNextSM->Connect("Clicked()", "AliTRDqaGuiClustersStack", this, "NextSM()");

  fGSelectSM->Connect("Selected(Int_t)", "AliTRDqaGuiClustersStack", this, "SelectSM(Int_t)");
  fGSelectStack->Connect("Selected(Int_t)", "AliTRDqaGuiClustersStack", this, "SelectStack(Int_t)");
  fGSelectView->Connect("Selected(Int_t)", "AliTRDqaGuiClustersStack", this, "SelectView(Int_t)");
  
  //fGPlay->Connect("Clicked()", "AliTRDqaGuiClustersStack", this, "Play()");

  // histograms
  /**/
  fGCanvas = new TGCompositeFrame(this);
  fGCanvas->SetLayoutManager(new TGMatrixLayout(fGCanvas,2,3,1,1));

  for(Int_t i=0; i<fgknCh; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("L%d",i), fGCanvas, 320, 300);
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

void AliTRDqaGuiClustersStack::SetQAFile(const char *filename) {
  //
  // Sets a file with histograms
  //
  
  //strncpy(fFileName,filename,256);
  fFileName = filename;

  for(Int_t i=0; i<fgknCh; i++) {
    if (fHistList[i]) delete fHistList[i];
    fHistList[i] = 0;
  }
  
  TFile *file = new TFile(filename);
  file->cd("TRD/RecPoints");
  
  // const char *opt[2] = {"colz", ""};
  if (fView == 0) CreateHistAmplitude();
  if (fView == 1) CreateHistTimeMPV();
  if (fView == 2 || fView == 3) CreateHistTimeCharge();
  
  for(Int_t i=0; i<fgknCh; i++) {
    fCanvasList[i]->GetCanvas()->cd();
    if (fHistList[i]) fHistList[i]->Draw(); //opt[i]);
    fCanvasList[i]->GetCanvas()->Update();
  }

  // style
  //  for(Int_t i=0; i<fgknCh; i++) {
  //   TPaveText *title = (TPaveText*)fCanvasList[i]->GetCanvas()->FindObject("title");
  //  title->SetX1NDC(0.7);
  //  title->SetX2NDC(0.95);
  // }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::SetStack(Int_t idxStack) {
  //
  // sets active stack
  // 

  fIdxStack = idxStack; 
  fGSelectSM->Select(fIdxSM, 0); 
  fGSelectStack->Select(fIdxStack, 0);
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::SetSM(Int_t idxSM) {
  //
  // sets active super module
  //

  fIdxSM = idxSM; 
  fGSelectSM->Select(fIdxSM, 0); 
  fGSelectStack->Select(fIdxStack, 0); 
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::SetView(Int_t idxView) {
  //
  // sets active data type
  //
  
  fView = idxView;
  fGSelectView->Select(idxView);
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::CreateHistTimeMPV() {
  //
  // builds histograms of the Most Probable Value distribution 
  // chamber-by-chamber from a 2D distribution
  //
  

  TH2D *fData = (TH2D*)gDirectory->Get(Form("qaTRD_recPoints_sigTime_sm%d", fIdxSM));
  if (!fData) return;

  Int_t nBins = fData->GetYaxis()->GetNbins();
  Double_t min = fData->GetYaxis()->GetXmin();
  Double_t max = fData->GetYaxis()->GetXmax();
  
  for(Int_t i=0; i<fgknCh; i++) {
    
    //printf("I = %d\n", i);

    Int_t det = fgknCh * fIdxStack + i;
    fHistList[i] = new TH1D(Form("det%d",det), Form("Det = %d;time bin; MPV",det), nBins, min, max);
    
    // fill the histograms;
    for(Int_t j=1; j<nBins+2; j++) {

      Double_t c = fHistList[i]->GetBinCenter(j);
      Int_t bin = fData->FindBin(det, c);
      Double_t value = fData->GetBinContent(bin);
      fHistList[i]->SetBinContent(j, value);
      // printf("..j=%d c=%lf bin=%d value=%lf\n", j, c, bin, value);  
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::CreateHistAmplitude() {
  //
  // builds histograms with amplitude
  // 
  
  //Info("createHist", "start");

  TH2D *fData = (TH2D*)gDirectory->Get("qaTRD_recPoints_amp");
  if (!fData) return;

  Int_t nBins = fData->GetYaxis()->GetNbins();
  Double_t min = fData->GetYaxis()->GetXmin();
  Double_t max = fData->GetYaxis()->GetXmax();
  
  for(Int_t i=0; i<fgknCh; i++) {
    
    //printf("I = %d\n", i);
    //if (fHistList[i]) delete fHistList[i];

    Int_t det = fgknCh * fgknStack * fIdxSM + fgknCh * fIdxStack + i;
    fHistList[i] = new TH1D(Form("det%d",det), Form("Det = %d;amplidtude",det), nBins, min, max);
    
    // fill the histograms;
    for(Int_t j=1; j<nBins+2; j++) {
      //printf("..j=%d\n", j);
      Double_t c = fHistList[i]->GetBinCenter(j);
      Int_t bin = fData->FindBin(det, c);
      Double_t value = fData->GetBinContent(bin);
      fHistList[i]->SetBinContent(j, value);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClustersStack::CreateHistTimeCharge() {
  //
  // builds histograms with time-charge distribution
  //

  TH3D *fData = (TH3D*)gDirectory->Get("qaTRD_recPoints_sigTime");
  if (!fData) return;

  Int_t nBins = fData->GetYaxis()->GetNbins();
  //Int_t nBinsCharge = fData->GetZaxis()->GetNbins();
  Double_t min = fData->GetYaxis()->GetXmin();
  Double_t max = fData->GetYaxis()->GetXmax();
  
  for(Int_t i=0; i<fgknCh; i++) {
    
    //printf("I = %d\n", i);
    //if (fHistList[i]) delete fHistList[i];

    Int_t det = fgknCh * fgknStack * fIdxSM + fgknCh * fIdxStack + i;
    const char *yaxis[2] = {"total charge", "number of clusters"};
    fHistList[i] = new TH1D(Form("det%d",det), Form("Det = %d;time bin;%s",det,yaxis[fView-2]), nBins, min, max);
    
    // fill the histograms;
    for(Int_t j=1; j<nBins+1; j++) {
      
      Double_t charge = 0;
      Double_t ncls = 0;
      for(Int_t k=0; k<201; k++) { // needs more robust loop
	Int_t bin = fData->FindBin(det, fHistList[i]->GetBinCenter(j), k);
	Double_t v = fData->GetBinContent(bin);
	charge += k * v;
	ncls += v;
      }
      // if (ncls > 1)
      if (fView == 2)
	fHistList[i]->SetBinContent(j, charge);

      if (fView == 3)
	fHistList[i]->SetBinContent(j, ncls);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
