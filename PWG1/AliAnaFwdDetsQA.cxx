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

//------------------------------
// Analysis task for quality-assurance
// of forward detectors ESD
//
// 12/06/2009 cvetan.cheshkov@cern.ch
//------------------------------

#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TParticle.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDMuonTrack.h"
#include "AliESDCaloCluster.h"
#include "AliRun.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliESDVZERO.h"

#include "AliAnaFwdDetsQA.h"

ClassImp(AliAnaFwdDetsQA)

AliAnaFwdDetsQA::AliAnaFwdDetsQA():
AliAnalysisTaskSE("AliAnaFwdDetsQA"),
  fListOfHistos(0),
  fT0vtxRec(0),
  fT0vtxRecGen(0),
  fT0time(0),
  fT0time2(0),
  fT0mult(0),
  fT0vtxRes(0),
  fT0ampl(0),
  fV0a(0),
  fV0c(0),
  fV0multA(0),
  fV0multC(0),
  fV0multAcorr(0),
  fV0multCcorr(0),
  fV0Acorr(0),
  fV0Ccorr(0),
  fV0ampl(0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

AliAnaFwdDetsQA::AliAnaFwdDetsQA(const char* name):
AliAnalysisTaskSE(name),
  fListOfHistos(0),
  fT0vtxRec(0),
  fT0vtxRecGen(0),
  fT0time(0),
  fT0time2(0),
  fT0mult(0),
  fT0vtxRes(0),
  fT0ampl(0),
  fV0a(0),
  fV0c(0),
  fV0multA(0),
  fV0multC(0),
  fV0multAcorr(0),
  fV0multCcorr(0),
  fV0Acorr(0),
  fV0Ccorr(0),
  fV0ampl(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  AliInfo("Constructor AliAnaFwdDetsQA");
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

TH1F * AliAnaFwdDetsQA::CreateHisto(const char* name, const char* title,Int_t nBins, 
					    Double_t xMin, Double_t xMax,
					    const char* xLabel, const char* yLabel)
{
  // helper method which can be used
  // in order to create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  result->SetOption("E");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  result->SetMarkerStyle(kFullCircle);
  return result;
}

TH1F *AliAnaFwdDetsQA::CreateEffHisto(const TH1F* hGen, const TH1F* hRec)
{
  // helper method which can be used
  // in order create an efficiency histogram
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  hEff->GetYaxis()->SetTitle("#epsilon [%]");
  
  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGenEff = hGen->GetBinContent(iBin);
    Double_t nRecEff = hRec->GetBinContent(iBin);
    if (nGenEff > 0) {
      Double_t eff = nRecEff/nGenEff;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt(eff*(1.-eff) / nGenEff);
      if (error < 1e-12) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);			
    }
    else {
      hEff->SetBinContent(iBin, -100.);
      hEff->SetBinError(iBin, 0);
    }
  }
  return hEff;
}


Bool_t AliAnaFwdDetsQA::FitHisto(TH1* histo, Double_t& res, Double_t& resError)
{
  // fit a gaussian to
  // a histogram
  static TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  Double_t maxFitRange = 2;
  
  if (histo->Integral() > 50) {
    Float_t mean = histo->GetMean();
    Float_t rms = histo->GetRMS();
    fitFunc->SetRange(mean - maxFitRange*rms, mean + maxFitRange*rms);
    fitFunc->SetParameters(mean, rms);
    histo->Fit(fitFunc, "QRI0");
    histo->GetFunction("fitFunc")->ResetBit(1<<9);
    res = TMath::Abs(fitFunc->GetParameter(2));
    resError = TMath::Abs(fitFunc->GetParError(2));
    return kTRUE;
  }
  return kFALSE;
}

void AliAnaFwdDetsQA::UserCreateOutputObjects()
{
  // Create histograms
  // Create output container
  AliInfo("AliAnaFwdDetsQA::UserCreateOutputObjects");
  fListOfHistos = new TList();
  
  fT0vtxRec = CreateHisto("hT0vtxRec", "Z vertex reconstructed with T0", 100, -25, 25, "Z_{vtx} [cm]", "");
  fT0time   = CreateHisto("hT0time", "Time0 reconstruction with T0", 5000, 10000, 20000, "t_{0} [ps]", "");
  fT0time2  = CreateHisto("hT0time2", "Time0 reconstruction with T0 (mult > 10)", 5000, 10000, 20000, "t_{0} [ps]", "");
  fT0mult   = CreateHisto("hT0mult", "Total reconstructed multiplicity in T0", 100, -0.5, 99.5, "Multiplicity", "");
  fT0vtxRes = CreateHisto("hT0vtxRes", "T0 Z vertex resolution", 100, -25, 25, "Delta(Z_{vtx}) [cm]", "");
  fT0ampl   = CreateHisto("hT0ampl","T0 multiplicity in single channel (all T0 channels)",400,-0.5,99.5);
  
  fT0vtxRecGen = new TH2F("hT0vtxRecGen", "T0 reconstructed vs generated Z vertex", 100, -25, 25, 100, -25, 25);
  fT0vtxRecGen->GetXaxis()->SetTitle("Generated Z vertex");
  fT0vtxRecGen->GetYaxis()->SetTitle("Reconstructed Z vertex");
  fT0vtxRecGen->SetMarkerStyle(kFullCircle);
  fT0vtxRecGen->SetMarkerSize(0.4);

  fV0a = CreateHisto("hV0a","Number of fired PMTs (V0A)",65,-0.5,64.5);
  fV0c = CreateHisto("hV0c","Number of fired PMTs (V0C)",65,-0.5,64.5);
  fV0multA = CreateHisto("hV0multA","Total reconstructed multiplicity (V0A)",100,0.,1000.);
  fV0multC = CreateHisto("hV0multC","Total reconstructed multiplicity (V0C)",100,0.,1000.);
  fV0ampl  = CreateHisto("hV0ampl","V0 multiplicity in single channel (all V0 channels)",400,-0.5,99.5);

  fV0multAcorr = new TH2F("hV0multAcorr", "Reconstructed vs generated (primaries only) multiplicity (V0A)",100,0.,500.,100,0.,1000.);
  fV0multAcorr->GetXaxis()->SetTitle("# of primaries in V0A acceptance");
  fV0multAcorr->GetYaxis()->SetTitle("Reconstructed mutiplicity in V0A");
  fV0multCcorr = new TH2F("hV0multCcorr", "Reconstructed vs generated (primaries only) multiplicity (V0C)",100,0.,500.,100,0.,1000.);
  fV0multCcorr->GetXaxis()->SetTitle("# of primaries in V0C acceptance");
  fV0multCcorr->GetYaxis()->SetTitle("Reconstructed mutiplicity in V0C");

  fV0Acorr = new TH2F("hV0Acorr", "Number of fired PMTs vs generated (primaries only) multiplicity (V0A)",100,0.,500.,65,-0.5,64.5);
  fV0Acorr->GetXaxis()->SetTitle("# of primaries in V0A acceptance");
  fV0Acorr->GetYaxis()->SetTitle("Number of fired PMTs in V0A");
  fV0Ccorr = new TH2F("hV0Ccorr", "Number of fired PMTs vs generated (primaries only) multiplicity (V0C)",100,0.,500.,65,-0.5,64.5);
  fV0Ccorr->GetXaxis()->SetTitle("# of primaries in V0C acceptance");
  fV0Ccorr->GetYaxis()->SetTitle("Number of fired PMTs in V0C");

  fListOfHistos->Add(fT0vtxRec);
  fListOfHistos->Add(fT0time);
  fListOfHistos->Add(fT0mult);
  fListOfHistos->Add(fT0vtxRecGen);
  fListOfHistos->Add(fT0vtxRes);
  fListOfHistos->Add(fV0a);
  fListOfHistos->Add(fV0c);
  fListOfHistos->Add(fV0multA);
  fListOfHistos->Add(fV0multC);
  fListOfHistos->Add(fV0multAcorr);
  fListOfHistos->Add(fV0multCcorr);
  fListOfHistos->Add(fV0Acorr);
  fListOfHistos->Add(fV0Ccorr);
  fListOfHistos->Add(fT0ampl);
  fListOfHistos->Add(fV0ampl);
  fListOfHistos->Add(fT0time2);
}

void AliAnaFwdDetsQA::UserExec(Option_t */*option*/)
{
  // The analysis code
  // goes here
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }	
	
  //Primary vertex needed
  TArrayF mcvtx(3);  
  mcEvent->GenEventHeader()->PrimaryVertex(mcvtx);

  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    Printf("ERROR: Could not retrieve MC stack");
    return;
  }

  Int_t nV0A = 0;
  Int_t nV0C = 0;
  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {//Check this loop again
    if (!stack->IsPhysicalPrimary(iTracks)) continue;
    AliMCParticle* track = (AliMCParticle*)mcEvent->GetTrack(iTracks);
    TParticle* particle = track->Particle();
    if (!particle) continue;
    if (track->Charge() == 0) continue;
    Double_t eta = particle->Eta();
    if (eta > 2.8 && eta < 5.1) {
      nV0A++;
    }
    if (eta > -3.7 && eta < -1.7) {
      nV0C++;
    }
  }

  // ESD information
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  const AliESDTZERO* esdT0 = esd->GetESDTZERO();
  Double_t t0zvtx = esdT0->GetT0zVertex();
  Double_t t0time = esdT0->GetT0();

  fT0vtxRec->Fill(t0zvtx);
  fT0time->Fill(t0time);
  fT0vtxRecGen->Fill(mcvtx[2],t0zvtx);
  t0zvtx *= -1.0;
  fT0vtxRes->Fill(t0zvtx - mcvtx[2]);

  Double_t t0mult = 0;
  for(Int_t i = 0; i < 24; i++) {
    t0mult += esdT0->GetT0amplitude()[i];
    fT0ampl->Fill(esdT0->GetT0amplitude()[i]);
  }

  fT0mult->Fill(t0mult);
  if (t0mult > 10)
    fT0time2->Fill(t0time);

  AliESDVZERO* esdV0 = esd->GetVZEROData();
  fV0a->Fill(esdV0->GetNbPMV0A());
  fV0c->Fill(esdV0->GetNbPMV0C());
  fV0multA->Fill(esdV0->GetMTotV0A());
  fV0multC->Fill(esdV0->GetMTotV0C());

  fV0multAcorr->Fill((Float_t)nV0A,esdV0->GetMTotV0A());
  fV0multCcorr->Fill((Float_t)nV0C,esdV0->GetMTotV0C());
  fV0Acorr->Fill((Float_t)nV0A,(Float_t)esdV0->GetNbPMV0A());
  fV0Ccorr->Fill((Float_t)nV0C,(Float_t)esdV0->GetNbPMV0C());

  for(Int_t i = 0; i < 64; i++) {
    fV0ampl->Fill(esdV0->GetMultiplicity(i));
  }
  // Post output data.
  PostData(1, fListOfHistos);
}

void AliAnaFwdDetsQA::Terminate(Option_t *)
{
  // Terminate
  // Store the output histos
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }
	
  fT0vtxRec = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  fT0time = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  fT0mult = dynamic_cast<TH1F*>(fListOfHistos->At(2));
  fT0vtxRecGen = dynamic_cast<TH2F*>(fListOfHistos->At(3));
  fT0vtxRes = dynamic_cast<TH1F*>(fListOfHistos->At(4));

  fV0a = dynamic_cast<TH1F*>(fListOfHistos->At(5));
  fV0c = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  fV0multA = dynamic_cast<TH1F*>(fListOfHistos->At(7));
  fV0multC = dynamic_cast<TH1F*>(fListOfHistos->At(8));
  fV0multAcorr = dynamic_cast<TH2F*>(fListOfHistos->At(9));
  fV0multCcorr = dynamic_cast<TH2F*>(fListOfHistos->At(10));
  fV0Acorr = dynamic_cast<TH2F*>(fListOfHistos->At(11));
  fV0Ccorr = dynamic_cast<TH2F*>(fListOfHistos->At(12));

  fT0ampl = dynamic_cast<TH1F*>(fListOfHistos->At(13));
  fV0ampl = dynamic_cast<TH1F*>(fListOfHistos->At(14));
  fT0time2 = dynamic_cast<TH1F*>(fListOfHistos->At(15));


  // draw the histograms if not in batch mode
  if (!gROOT->IsBatch()) {
    new TCanvas;
    fT0vtxRec->DrawCopy("E");
    new TCanvas;
    fT0time->DrawCopy("E");
    new TCanvas;
    fT0mult->DrawCopy("E");
    new TCanvas;
    fT0vtxRes->DrawCopy("E");
    new TCanvas;
    fT0vtxRecGen->DrawCopy();
    new TCanvas;
    fV0a->DrawCopy("E");
    new TCanvas;
    fV0c->DrawCopy("E");
    new TCanvas;
    fV0multA->DrawCopy("E");
    new TCanvas;
    fV0multC->DrawCopy("E");
    new TCanvas;
    fV0multAcorr->DrawCopy();
    new TCanvas;
    fV0multCcorr->DrawCopy();
    new TCanvas;
    fV0Acorr->DrawCopy();
    new TCanvas;
    fV0Ccorr->DrawCopy();
    new TCanvas;
    fT0ampl->DrawCopy("E");
    new TCanvas;
    fV0ampl->DrawCopy("E");
    new TCanvas;
    fT0time2->DrawCopy("E");
  }

  // write the output histograms to a file
  TFile* outputFile = TFile::Open("FwdDetsQA.root", "recreate");
  if (!outputFile || !outputFile->IsOpen())
    {
      Error("AliAnaFwdDetsQA", "opening output file FwdDetsQA.root failed");
      return;
    }
  fT0vtxRec->Write();
  fT0time->Write();
  fT0mult->Write();
  fT0vtxRecGen->Write();
  fT0vtxRes->Write();
  fT0ampl->Write();
  fT0time2->Write();

  fV0a->Write();
  fV0c->Write();
  fV0multA->Write();
  fV0multC->Write();
  fV0multAcorr->Write();
  fV0multCcorr->Write();
  fV0Acorr->Write();
  fV0Ccorr->Write();
  fV0ampl->Write();

  outputFile->Close();
  delete outputFile;

  //delete esd;
  Info("AliAnaFwdDetsQA", "Successfully finished");
}






























































