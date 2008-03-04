/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * withount fee, provided that the abovÅe copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDqaBlackEvents.cxx 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  QA of black events                                                    //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"

#include "AliTRDgeometry.h"
#include "AliTRDrawStreamTB.h"
#include "AliTRDqaBlackEvents.h"

ClassImp(AliTRDqaBlackEvents)


///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents() 
  :TObject() 
  ,fMinNoise(0.5)
  ,fMaxNoise(2) 
  ,fFitType(0)
  ,fnEvents(0)
{
  //
  // Constructor 
  // to create the histograms call Init()
  //
}

///////////////////////////////////////////////////////////////////////////////////////////////////

AliTRDqaBlackEvents::AliTRDqaBlackEvents(const AliTRDqaBlackEvents &qa) 
  :TObject(qa) 
  ,fMinNoise(0.5)
  ,fMaxNoise(2) 
  ,fFitType(0)
  ,fnEvents(0)
{
  //
  // Copy constructor 
  // to create the histograms call Init()
  //
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Init() 
{
  //
  // creates histograms 
  // 

  //TFile *file = new 
  //Info("Init", "Statring");

  fnEvents = 0;

  for(Int_t i=0; i<540; i++) {
    fNPoint[i]  = new TH2D(Form("entries_%d", i), "",  16, -0.5, 15.5, 144, -0.5, 143.5);
    fData[i]    = new TH3D(Form("data_%d", i), "", 16, -0.5, 15.5, 144, -0.5, 143.5, 50, -0.5, 49.5);
    fChPed[i]   = new TH2D(Form("ped_%d", i), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    fChNoise[i] = new TH2D(Form("noise_%d", i), "", 16, -0.5, 15.5, 144, -0.5, 143.5);
    fPed[i]     = new TH1D(Form("pedDist_%d", i), ";pedestals (ADC counts)", 100, 5, 15);
    fNoise[i]   = new TH1D(Form("noiseDist_%d", i), ";noise (ADC counts)", 100, 0, 5); 
    fSignal[i]  = new TH1D(Form("signal_%d", i), "signal (ADC counts)", 100, -0.5, 99.5);
  }

  fOccupancy = new TH1D("occupancy", "", 20, -0.5, 19.5);

  //Info("Init", "Done");
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Reset() 
{
  //
  // Resets the histograms
  //

  for(Int_t i=0; i<540; i++) {
    fData[i]->Reset();
    fChPed[i]->Reset();
    fChNoise[i]->Reset();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDqaBlackEvents::AddEvent(AliTRDrawStreamTB *data) 
{
  //
  // Add an event
  //

  
  Char_t isUsed[540][16][144]; 
  for(Int_t i=0; i<540; i++)
    for(Int_t j=0; j<16; j++)
      for(Int_t k=0; k<144; k++)
  	isUsed[i][j][k] = 0;
 
  Int_t nb = 0;
  while (data->Next()) {

    Int_t det = data->GetDet();
    Int_t row = data->GetRow();
    Int_t col = data->GetCol();
    Int_t *sig = data->GetSignals();
    nb++;

    //printf("det = %d\n", det);
    
    if (det<0 || det>=540) continue;
    isUsed[det][row][col]++;

    // if (!isUsed[det][data->GetRow()][data->GetCol()]) {
    //  isUsed[det][data->GetRow()][data->GetCol()] = 1;
    //  continue;
    // }

    fNPoint[det]->Fill(row, col);
    
    for(Int_t k=0; k<30; k++) { /// to be corrected
      fSignal[det]->Fill(sig[k]);
      fData[det]->Fill(row, col, sig[k]);
    }
  }
  
  for(Int_t i=0; i<540; i++) {
    if (i != 0 && i != 8) continue;
    for(Int_t j=0; j<16; j++)
      for(Int_t k=0; k<144; k++)
  	fOccupancy->Fill(isUsed[i][j][k]);
  }

  
  fnEvents++;
  return nb;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::Process(const char *filename) 
{
  //
  // Process something
  //

  Int_t map[540];
  
  TH1D *hist = new TH1D("fitSignal", "", 50, -0.5, 49.5);
  TF1 *fit = new TF1("fit", "gaus(0)", 0, 20);
  fit->SetParameters(1e3, 10, 1);
    
  for(Int_t i=0; i<540; i++) {
    
    map[i] = 0;
    if (fData[i]->GetSum() < 10) continue;
    map[i] = 1;

    Info("process", "processing chamber %d", i);

    for(Int_t j=0; j<fData[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fData[i]->GetYaxis()->GetNbins(); k++) {
	
	// project the histogramm
	hist->Reset();
	for(Int_t bb=0; bb<50; bb++) {
	  Int_t dataBin = fData[i]->FindBin(j, k, bb);
	  Double_t v = fData[i]->GetBinContent(dataBin);
	  hist->SetBinContent(bb+1, v);
	}

	//TH1D *hist = fData[i]->ProjectionZ(Form("pad_%_%d_%d", i, j, k), j+1, j+1, k+1, k+1);
	
	Int_t bin = fChPed[i]->FindBin(j, k);

	if (hist->GetSum() > 1) {
	  
	  Double_t ped = 0, noise = 0;

	  if (fFitType == 0) {
	    fit->SetParameters(1e3, 10, 1);
	    hist->Fit(fit, "q0", "goff", 0, 20);
	    TF1 *f = hist->GetFunction("fit");
	    ped = TMath::Abs(f->GetParameter(1));
	    noise = TMath::Abs(f->GetParameter(2));
	  } else {
	    ped = hist->GetMean();
	    noise = hist->GetRMS();
	  }

	  fChPed[i]->SetBinContent(bin, ped);
	  fChNoise[i]->SetBinContent(bin, noise);
	  
	  fPed[i]->Fill(ped);
	  fNoise[i]->Fill(noise);

	} else {
	  fChPed[i]->SetBinContent(bin, 0);
	  fChNoise[i]->SetBinContent(bin, 0);
	}
	
	//delete hist;
      }
    }
  }

  // normalize number of entries histos

  
  Int_t max = 0;
  for(Int_t i=0; i<540; i++) { 
    if (!map[i]) continue;
    for(Int_t j=0; j<fNPoint[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fNPoint[i]->GetYaxis()->GetNbins(); k++) {
	Int_t dataBin = fNPoint[i]->FindBin(j, k);
	Double_t v = fNPoint[i]->GetBinContent(dataBin);
	if (v > max) max = (Int_t)v;
      }
    }
  }
  
  for(Int_t i=0; i<540; i++) {
    
    if (!map[i]) continue;
    
    fNPointDist[i] = new TH1D(Form("entriesDist_%d", i), ";number of events", max+2, -0.5, max+1.5);
    
    for(Int_t j=0; j<fNPoint[i]->GetXaxis()->GetNbins(); j++) {
      for(Int_t k=0; k<fNPoint[i]->GetYaxis()->GetNbins(); k++) {
	Int_t dataBin = fNPoint[i]->FindBin(j, k);
	Double_t v = fNPoint[i]->GetBinContent(dataBin);
	//if (v > fnEvents) printf("N = %d V = %lf\n", fnEvents, v);
	fNPointDist[i]->Fill(v); 
      }
    }
    
    fNPoint[i]->Scale(1./fnEvents);
  }
  
  
  TFile *file = new TFile(filename, "UPDATE");
  for(Int_t i=0; i<540; i++) {
    if (!map[i]) continue; 
    fChPed[i]->Write();
    fChNoise[i]->Write();
    fNPoint[i]->Write();
    fNPointDist[i]->Write();
    fPed[i]->Write();
    fNoise[i]->Write();
    fSignal[i]->Write();
  }
  fOccupancy->Write();
  file->Close();
  delete file;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::DrawChamber(const char *filename, Int_t det, Int_t w, Int_t h) 
{
  //
  // Draw raport for one chamber: 
  // pedestal map, noise map, distribution of pedestal and noise
  // 
  // input:
  // name of the file with histograms (created with Process())
  // detector Id (0 - 539)
  // 

  // setup global style
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.05);

  TFile *file = new TFile(filename, "READ");

  TCanvas *c = new TCanvas("blackEvents",Form("blackEvents %d",det), w, h);
  c->SetVertical(kFALSE);
  c->Divide(3,1, 0.01, 0.01);
  c->cd(3);
  
  TPad *mPad = (TPad*) gPad;
  mPad->Divide(1,2,0.01,0.01);
  
  c->cd(1);
  TH2D *h2 = (TH2D*)file->Get(Form("ped_%d",det));
  h2->SetMinimum(5);
  h2->SetMaximum(15);
  h2->SetTitle(";Z direction;#phi direction");
  h2->Draw("colz");
  
  c->cd(2);
  h2 = (TH2D*)file->Get(Form("noise_%d",det));
  h2->SetMinimum(fMinNoise);
  h2->SetMaximum(fMaxNoise);
  h2->SetTitle(";Z direction;#phi direction");
  h2->Draw("colz");
  
  mPad->cd(1);
  //gPad->SetLogy();
  TH1D *h1 = (TH1D*)file->Get(Form("pedDist_%d", det));
  h1->Draw();
  
  mPad->cd(2);
  gPad->SetLogy();
  h1 = (TH1D*)file->Get(Form("noiseDist_%d", det));
  h1->Draw();			 
  
  h1->Fit("gaus");
  TF1 *f = h1->GetFunction("gaus");
  const char *tt = Form("#mu = %.2f #sigma = %0.2f ", f->GetParameter(1),f->GetParameter(2));
  TLatex *ll = new TLatex(2, 100, tt);
  ll->SetTextSize(0.06);
  ll->Draw();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void AliTRDqaBlackEvents::DrawSm(const char *filename, Int_t sm, Int_t w, Int_t h) 
{
  //
  // ????????????
  //
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadBottomMargin(0.05);  
  //gStyle->SetPadLeftMargin(0.02);  
  //gStyle->SetPadRightMargin(0.02);

  TFile *file = new TFile(filename, "READ");

  TCanvas *c = new TCanvas("blackEventsSM",Form("blackEvents SM %d",sm), w, h);
  c->SetVertical(kFALSE);
  c->Divide(5, 6, 0.001, 0.01);
  
  for(Int_t i=0; i<30; i++) {
    
    TH2D *h2 = (TH2D*)file->Get(Form("noise_%d",i+30*sm));
    if (!h2) continue;
    h2->SetMinimum(fMinNoise);
    h2->SetMaximum(fMaxNoise);

    // to be replaced by the official calculation
    Int_t stack = i/6;
    Int_t layer = i%6;
    Int_t index = (5-layer)*5 + stack + 1;
    //printf("%d %d %d %d\n", i, stack, layer, index);
    c->cd(index);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.02);

    h2->Draw("col");
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
