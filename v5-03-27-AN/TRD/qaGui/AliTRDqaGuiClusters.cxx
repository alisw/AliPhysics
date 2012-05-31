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

/* $Id: AliTRDqaGuiClusters.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters on the full detector level.
// It displays histograms created by the AliTRDQADataMakerRec 
// run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiClusters.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiClusters)

const Int_t AliTRDqaGuiClusters::fgkLogList[4] = {0,1,0,1};

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiClusters::AliTRDqaGuiClusters() 
  : TGCompositeFrame() {
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 4; i++) {
    fNameList[i]   = 0x0;
    fHistList[i]   = 0x0;
    fCanvasList[i] = 0x0;
    for (Int_t j = 0; j < 3; j++) {
      fHistRefs[i][j] = 0x0;
    }
  }

}

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiClusters::AliTRDqaGuiClusters(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 400) {
  //
  // main constructor
  //
  
  SetLayoutManager(new TGMatrixLayout(this,2,2,1,1));

  fNameList[0] = "detMap";
  fNameList[1] = "nCls";
  fNameList[2] = "signal";
  fNameList[3] = "ampMPV";

  for(Int_t i=0; i<4; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(fNameList[i], this, 480, 320);
    AddFrame(fCanvasList[i]);
  }
  
  for(Int_t i=0; i<4; i++) {
    fHistList[i] = 0;
  }

}

//////////////////////////////////////////////////////////////////////////////////
void AliTRDqaGuiClusters::SetQAFile(const char *filename) {
  //
  // sets the file with histograms
  //

  for(Int_t i=0; i<4; i++) {
    if (fHistList[i]) delete fHistList[i];
    for(Int_t j=0; j<3; j++) 
      if (fHistRefs[i][j]) delete fHistRefs[i][j];
  }
  
  const char *opt[4] = {"colz", "", "", ""};

  TFile *file = new TFile(filename);
  file->cd("TRD/RecPoints");
  
  for(int i=0; i<4; i++) {

    fCanvasList[i]->GetCanvas()->cd();
    gPad->SetLogy(fgkLogList[i]);

    fHistList[i] = (TH1D*)gDirectory->Get(Form("qaTRD_recPoints_%s", fNameList[i]));
    if (fHistList[i]) fHistList[i]->Draw(opt[i]);
    // if (fgkLogList[i]) fHistList[i]->SetMinimum(0.1);
    
    TH1D *refHist = (TH1D*)gDirectory->Get(Form("qaTRD_recPoints_%s_%s", fNameList[i], "ref"));
    
    if (refHist) {
      BuildColor(i, refHist);
      for(Int_t j=0; j<3; j++) fHistRefs[i][j]->Draw("SAME");
      delete refHist;
    }
    
    fCanvasList[i]->GetCanvas()->Update();
  }
}

//////////////////////////////////////////////////////////////////////////////////
/*
TH2D *AliTRDqaGuiClusters::BuildHisto(TH1D *ref, TH1D *data) {
  
  Int_t nbinsx = ref->GetNbinsX();
  Double_t minx = ref->GetXaxis()->GetXmin();
  Double_t maxx = ref->GetXaxis()->GetXmax();
  Double_t min  = data->GetMinimum(); 
  Double_t max  = data->GetMaximum();
  
  TH2D *pad = new TH2D("pad", "", nbinsx, minx, maxx, 1, 0.7, 0.9);

  // rewriting
  for(Int_t i=0; i<nbinsx; i++) {
    Double_t x = ref->GetBinCenter(i+1);
    pad->Fill(x, 0.8, ref->GetBinContent(i+1));
  }
  
  pad->SetMinimum(-1);
  return pad;
}
*/

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClusters::BuildColor(Int_t i, TH1D *ref) {
  
  const Int_t nHist = 3;
  const Int_t clr[nHist] = {3, 5, 2};

  for(Int_t j=0; j<nHist; j++) {
    fHistRefs[i][j] = (TH1D*)fHistList[i]->Clone(Form("%s_%d", fHistList[i]->GetName(), j));
    fHistRefs[i][j]->SetFillColor(clr[j]);
    fHistRefs[i][j]->SetLineColor(clr[j]);   
  }
  
  for(Int_t k=0; k<ref->GetNbinsX(); k++) {
    Double_t v = ref->GetBinContent(k+1);
    if (v < 0.3) fHistRefs[i][1]->SetBinContent(k+1, 0);
    if (v < 0.7) fHistRefs[i][2]->SetBinContent(k+1, 0);
  }
  
}

//////////////////////////////////////////////////////////////////////////////////
