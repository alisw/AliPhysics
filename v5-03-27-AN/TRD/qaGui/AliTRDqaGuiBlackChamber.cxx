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

/* $Id: AliTRDqaGuiBlackChamber.cxx 23871 2008-02-12 11:48:20Z hristov $ */

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

#include "AliTRDqaGuiBlackChamber.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include "TLine.h"
#include "TPaveText.h"
#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiBlackChamber)

const Int_t AliTRDqaGuiBlackChamber::fgknSM = 18;
const Int_t AliTRDqaGuiBlackChamber::fgknChamber = 30;
//const Int_t AliTRDqaGuiBlackChamber::fgknCh = 6;

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackChamber::AliTRDqaGuiBlackChamber() 
  : fView(0),
    fSetRangePed(0),
    fSetRangeNoise(0),
    fIdxSM(0),
    fIdxChamber(0),
    fFileName(0x0),
    fGPanel(0),
    fGCanvas(0),
    fGCanvasUp(0),
    fGCanvasDown(0),
    fGSelectSM(0),
    fGSelectChamber(0),
    fGSelectView(0),
    fGPrevSM(0),
    fGPrevChamber(0), 
    fGNextSM(0),
    fGNextChamber(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 2; i++) {
    fRangePed[i]   = 0.0;
    fRangeNoise[i] = 0.0;
  }

  for (Int_t j = 0; j < 5; j++) {
    fCanvasList[j] = 0x0;
    fHistList[j]   = 0x0;    
  }

  //strncpy(fFileName,"",256);

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackChamber::AliTRDqaGuiBlackChamber(TGWindow *parent)
  : TGCompositeFrame(parent, 720, 500), 
    fView(0),
    fSetRangePed(0),
    fSetRangeNoise(0),
    fIdxSM(0),
    fIdxChamber(0),
    fFileName(0x0),
    fGPanel(0),
    fGCanvas(0),
    fGCanvasUp(0),
    fGCanvasDown(0),
    fGSelectSM(0),
    fGSelectChamber(0),
    fGSelectView(0),
    fGPrevSM(0),
    fGPrevChamber(0),
    fGNextSM(0),
    fGNextChamber(0)
{
  //
  // main constructor
  // 
  
   // steering panel 
  
  SetLayoutManager(new TGVerticalLayout(this));
  fGPanel = new TGHorizontalFrame(this);

  // fGLabel = new TGLabel(fGPanel, "Current Chamber: ");
  fGPrevSM = new TGTextButton(fGPanel, "Prev SM");
  fGPrevChamber = new TGTextButton(fGPanel, "Prev Chamber");

  fGNextSM = new TGTextButton(fGPanel, "Next SM");
  fGNextChamber = new TGTextButton(fGPanel, "Next Chamber");

  fGSelectSM = new TGComboBox(fGPanel);
  for(int i=0; i<fgknSM; i++) fGSelectSM->AddEntry(Form("SM %d", i), i);
  fGSelectSM->Resize(100, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectSM->Select(fIdxSM);

  fGSelectChamber = new TGComboBox(fGPanel);
  for(int i=0; i<fgknChamber; i++) fGSelectChamber->AddEntry(Form("Chamber %d", i), i);
  fGSelectChamber->Resize(100, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectChamber->Select(fIdxChamber);

  // vew
  fGSelectView = new TGComboBox(fGPanel);
  fGSelectView->AddEntry("pedestals",0);
  fGSelectView->AddEntry("entiries", 1);
  fGSelectView->Resize(150, (Int_t)(fGPrevSM->GetHeight()*1.4));
  fGSelectView->Select(0);


  //fGPlay = new TGTextButton(fGPanel, "PLAY");

  TGLayoutHints *hint = new TGLayoutHints(kLHintsNormal, 5, 5, 5, 5);
  
  // fGPanel->AddFrame(fGLabel, hint);
  fGPanel->AddFrame(fGPrevSM, hint);
  fGPanel->AddFrame(fGPrevChamber, hint);

  fGPanel->AddFrame(fGSelectSM, hint);
  fGPanel->AddFrame(fGSelectChamber, hint);

  fGPanel->AddFrame(fGNextChamber, hint);
  fGPanel->AddFrame(fGNextSM, hint);

  fGPanel->AddFrame(fGSelectView, hint);
  //fGPanel->AddFrame(fGPlay, hint);

  AddFrame(fGPanel);

  // panel logic
  fGPrevChamber->Connect("Clicked()", "AliTRDqaGuiBlackChamber", this, "PreviusChamber()");
  fGNextChamber->Connect("Clicked()", "AliTRDqaGuiBlackChamber", this, "NextChamber()");
  fGPrevSM->Connect("Clicked()", "AliTRDqaGuiBlackChamber", this, "PreviusSM()");
  fGNextSM->Connect("Clicked()", "AliTRDqaGuiBlackChamber", this, "NextSM()");

  fGSelectSM->Connect("Selected(Int_t)", "AliTRDqaGuiBlackChamber", this, "SelectSM(Int_t)");
  fGSelectChamber->Connect("Selected(Int_t)", "AliTRDqaGuiBlackChamber", this, "SelectChamber(Int_t)");
  
  fGSelectView->Connect("Selected(Int_t)", "AliTRDqaGuiBlackChamber", this, "SelectView(Int_t)");

  //fGPlay->Connect("Clicked()", "AliTRDqaGuiBlackChamber", this, "Play()");

  // histograms
  
  // fGCanvas = new TGCompositeFrame(this);
  // fGCanvas->SetLayoutManager(new TGMatrixLayout(fGCanvas,2,2,0,0));
  
  //for(Int_t i=0; i<4; i++) {
  //  fCanvasList[i] = new TRootEmbeddedCanvas(Form("L%d",i), fGCanvas, 480, 300);
  //  fGCanvas->AddFrame(fCanvasList[i]);
  // }


  fGCanvasUp = new TGCompositeFrame(this);
  fGCanvasUp->SetLayoutManager(new TGMatrixLayout(fGCanvasUp,1,2,1,1));

  for(Int_t i=0; i<2; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("L%d",i), fGCanvasUp, 480, 400);
    fGCanvasUp->AddFrame(fCanvasList[i]);    
    fCanvasList[i]->GetCanvas()->SetTopMargin(0.05);
    fCanvasList[i]->GetCanvas()->SetRightMargin(0.15);
    fCanvasList[i]->GetCanvas()->SetBottomMargin(0.05);
  }
  
  
  fGCanvasDown = new TGCompositeFrame(this);
  fGCanvasDown->SetLayoutManager(new TGMatrixLayout(fGCanvasDown, 1,3,1,1));
    
  for(Int_t i=2; i<5; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("L%d",i), fGCanvasDown, 320, 300);
    fGCanvasDown->AddFrame(fCanvasList[i]);   
    fCanvasList[i]->GetCanvas()->SetTopMargin(0.05);
    fCanvasList[i]->GetCanvas()->SetRightMargin(0.05);
  }
  
  for(Int_t i=0; i<5; i++) {
    fHistList[i] = 0;
  }

  AddFrame(fGCanvasUp);
  AddFrame(fGCanvasDown);

}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackChamber::SetQAFile(const char *filename) {
  //
  // sets a file with histograms
  //

  //const char *names[5] = {"ped", "noise", "pedDist", "noiseDist", "signal"};
  const char *names[10] = {
    "ped", "noise", "pedDist", "noiseDist", "signal",
    "entries", "entriesRM", "entriesDist", "", ""
  };
  const char *opt[10] = {"colz", "colz", "", "", "", "colz", "colz", "", "", ""};
  const Int_t kLogy[10] = {0, 0, 1, 1, 1, 0, 0, 1, 1,1};
  
  //strncpy(fFileName,filename,256);
  fFileName = filename;
 
  for(int i=0; i<5; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  TFile *file = new TFile(filename);

  // printf("%d %lf %lf\n", fSetRangePed, fRangePed[0], fRangePed[1]);
  //printf("%d %lf %lf\n", fSetRangeNoise, fRangeNoise[0], fRangeNoise[1]);

  for(Int_t i=0; i<5; i++) {

    Int_t index = fIdxSM * 30 + fIdxChamber;
    const char *nn = Form("%s_%d", names[i+5*fView], index);
    //printf("%s\n", nn);
    fHistList[i] = (TH1*)file->Get(nn); //Form("%s_$d", names[fIdxType], index));
    if (!fHistList[i]) continue;

    if ( (fView == 1) && (i == 0)) {
      fHistList[i]->SetMinimum(0);
      fHistList[i]->SetMaximum(2);
    }

    if ( (fView == 0) && (i == 0)  && fSetRangePed) {
      fHistList[i]->SetMinimum(fRangePed[0]);
      fHistList[i]->SetMaximum(fRangePed[1]);
    }
    
    if ( (fView == 0) && (i == 1) && fSetRangeNoise) {
      fHistList[i]->SetMinimum(fRangeNoise[0]);
      fHistList[i]->SetMaximum(fRangeNoise[1]);
    }

    fCanvasList[i]->GetCanvas()->cd();
    fCanvasList[i]->GetCanvas()->SetLogy(kLogy[i+5*fView]);
    if (fHistList[i]) fHistList[i]->Draw(opt[i+5*fView]);
    //fCanvasList[i]->GetCanvas()->Update();
  }
  
  // mcm lines
  TLine *line;
  for(Int_t i=1; i<8; i++) {

    fCanvasList[0]->GetCanvas()->cd();
    line = new TLine(0, i*18-0.5, 15, i*18-0.5);
    line->SetLineStyle(2);
    if (i!=4) line->SetLineStyle(3);
    line->Draw();

    fCanvasList[1]->GetCanvas()->cd();
    line = new TLine(0, i*18-0.5, 15, i*18-0.5);
    line->SetLineStyle(2);
    if (i!=4) line->SetLineStyle(3);
    line->Draw();    
  }
  
  for(Int_t i=1; i<4; i++) {
    
    fCanvasList[0]->GetCanvas()->cd();
    line = new TLine(i*4-0.5, 0, i*4-0.5, 143);
    line->SetLineStyle(2);
    line->Draw();

    fCanvasList[1]->GetCanvas()->cd();
    line = new TLine(i*4-0.5, 0, i*4-0.5, 143);
    line->SetLineStyle(2);
    line->Draw();  
  }

  


  
  for(Int_t i=0; i<5; i++)
    fCanvasList[i]->GetCanvas()->Update();
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackChamber::SetChamber(Int_t idxChamber) {
  //
  // sets active chamber 
  //

  fIdxChamber = idxChamber; 
  fGSelectSM->Select(fIdxSM, 0); 
  fGSelectChamber->Select(fIdxChamber, 0);
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackChamber::SetSM(Int_t idxSM) {
  //
  // sets active supermodule
  //

  fIdxSM = idxSM; 
  fGSelectSM->Select(fIdxSM, 0); 
  fGSelectChamber->Select(fIdxChamber, 0); 
  SetQAFile(fFileName);
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackChamber::SetView(Int_t idxView) {
  //
  // sets active view
  //
  
  fView = idxView;
  fGSelectView->Select(idxView);
  SetQAFile(fFileName);
}
//////////////////////////////////////////////////////////////////////////////////
