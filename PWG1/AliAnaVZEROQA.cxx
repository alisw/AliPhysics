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
// of VZERO ESD
//
// 05/12/2009 cvetan.cheshkov@cern.ch
//------------------------------

#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"

#include "AliAnaVZEROQA.h"

ClassImp(AliAnaVZEROQA)

AliAnaVZEROQA::AliAnaVZEROQA():
AliAnalysisTaskSE("AliAnaVZEROQA"),
  fListOfHistos(0),

  fhAdcNoTimeA(0),
  fhAdcWithTimeA(0),
  fhAdcNoTimeC(0),
  fhAdcWithTimeC(0),

  fhAdcPMTNoTime(0),
  fhAdcPMTWithTime(0),
 
  fhTimeA(0),
  fhTimeC(0),

  fhWidthA(0),
  fhWidthC(0),

  fhTimePMT(0),
  fhWidthPMT(0),

  fhAdcWidthA(0),
  fhAdcWidthC(0),

  fhTimeCorr(0),

  fhAdcTimeA(0),
  fhAdcTimeC(0),

  fV0a(0),
  fV0c(0),
  fV0multA(0),
  fV0multC(0),
  fV0ampl(0),

  fhTimePMTCorr(0),
  fhEvents(0),

  fhVtxXYBB(0),
  fhVtxZBB(0),
  fhVtxXYBGA(0),
  fhVtxZBGA(0),
  fhVtxXYBGC(0),
  fhVtxZBGC(0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

AliAnaVZEROQA::AliAnaVZEROQA(const char* name):
AliAnalysisTaskSE(name),
  fListOfHistos(0),
  fhAdcNoTimeA(0),
  fhAdcWithTimeA(0),
  fhAdcNoTimeC(0),
  fhAdcWithTimeC(0),

  fhAdcPMTNoTime(0),
  fhAdcPMTWithTime(0),
 
  fhTimeA(0),
  fhTimeC(0),

  fhWidthA(0),
  fhWidthC(0),

  fhTimePMT(0),
  fhWidthPMT(0),

  fhAdcWidthA(0),
  fhAdcWidthC(0),

  fhTimeCorr(0),

  fhAdcTimeA(0),
  fhAdcTimeC(0),

  fV0a(0),
  fV0c(0),
  fV0multA(0),
  fV0multC(0),
  fV0ampl(0),

  fhTimePMTCorr(0),
  fhEvents(0),

  fhVtxXYBB(0),
  fhVtxZBB(0),
  fhVtxXYBGA(0),
  fhVtxZBGA(0),
  fhVtxXYBGC(0),
  fhVtxZBGC(0)
{
  // Constructor
  AliInfo("Constructor AliAnaVZEROQA");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

TH1F * AliAnaVZEROQA::CreateHisto1D(const char* name, const char* title,Int_t nBins, 
				    Double_t xMin, Double_t xMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  result->SetOption("E");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  result->SetMarkerStyle(kFullCircle);
  return result;
}

TH2F * AliAnaVZEROQA::CreateHisto2D(const char* name, const char* title,Int_t nBinsX, 
				    Double_t xMin, Double_t xMax,
				    Int_t nBinsY,
				    Double_t yMin, Double_t yMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH2F* result = new TH2F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}

void AliAnaVZEROQA::UserCreateOutputObjects()
{
  // Create histograms
  AliInfo("AliAnaVZEROQA::UserCreateOutputObjects");
  // Create output container
  fListOfHistos = new TList();

  fhAdcNoTimeA = CreateHisto1D("hAdcNoTimeA","ADC (no Leading Time) V0A",200,0,200,"ADC charge","Entries");
  fhAdcWithTimeA = CreateHisto1D("hAdcWithTimeA","ADC ( with Leading Time) V0A",200,0,200,"ADC charge","Entries");
  fhAdcNoTimeC = CreateHisto1D("hAdcNoTimeC","ADC (no Leading Time) V0C",200,0,200,"ADC charge","Entries");
  fhAdcWithTimeC = CreateHisto1D("hAdcWithTimeC","ADC ( with Leading Time) V0C",200,0,200,"ADC charge","Entries");
  
  fhAdcPMTNoTime = CreateHisto2D("hadcpmtnotime","ADC vs PMT index (no leading time)",64,-0.5,63.5,200,0,200,"PMT index","ADC charge");
  fhAdcPMTWithTime = CreateHisto2D("hadcpmtwithtime","ADC vs PMT index (with leading time)",64,-0.5,63.5,200,0,200,"PMT index","ADC charge");

  fhTimeA = CreateHisto1D("htimepmtA","Time measured by TDC V0A",400,-100,100,"Leading time (ns)","Entries");
  fhTimeC = CreateHisto1D("htimepmtC","Time measured by TDC V0C",400,-100,100,"Leading time (ns)","Entries");

  fhWidthA = CreateHisto1D("hwidthA","Signal width measured by TDC V0A",200,0,100,"Signal width (ns)","Entries");
  fhWidthC = CreateHisto1D("hwidthC","Signal width measured by TDC V0C",200,0,100,"Signal width (ns)","Entries");

  fhTimePMT = CreateHisto2D("htimepmt","Time measured by TDC vs PMT index",64,-0.5,63.5,200,0,100,"PMT Index","Leading time (ns)");
  fhWidthPMT = CreateHisto2D("hwidthpmt","Time width vs PMT index",64,-0.5,63.5,200,0,100,"PMT Index","Signal width (ns)");

  fhAdcWidthA = CreateHisto2D("hadcwidthA","Time width vs ADC V0A",200,0,200,200,0,100,"ADC charge","Width (ns)");
  fhAdcWidthC = CreateHisto2D("hadcwidthC","Time width vs ADC V0C",200,0,200,200,0,100,"ADC charge","Width (ns)");

  fhTimeCorr = CreateHisto2D("htimecorr","Average time C side vs. A side",200,0,100,200,0,100,"Time V0A (ns)","Time V0C (ns");

  fhAdcTimeA = CreateHisto2D("hAdcTimeA","ADC vs Time V0A",1000,-100,100,200,0,200,"Time (ns)","ADC charge");
  fhAdcTimeC = CreateHisto2D("hAdcTimeC","ADC vs Time V0C",1000,-100,100,200,0,200,"Time (ns)","ADC charge");

  fV0a = CreateHisto1D("hV0a","Number of fired PMTs (V0A)",65,-0.5,64.5);
  fV0c = CreateHisto1D("hV0c","Number of fired PMTs (V0C)",65,-0.5,64.5);
  fV0multA = CreateHisto1D("hV0multA","Total reconstructed multiplicity (V0A)",100,0.,1000.);
  fV0multC = CreateHisto1D("hV0multC","Total reconstructed multiplicity (V0C)",100,0.,1000.);
  fV0ampl  = CreateHisto1D("hV0ampl","V0 multiplicity in single channel (all V0 channels)",400,-0.5,99.5);

  fhTimePMTCorr = CreateHisto2D("htimepmtcorr","Time measured by TDC (corrected for slewing, channels aligned) vs PMT index",64,-0.5,63.5,200,0,100,"PMT Index","Leading time (ns)");
  fhEvents = CreateHisto2D("hEvents","V0C vs V0A (empty,bb,bg)",3,-0.5,2.5,3,-0.5,2.5);

  fhVtxXYBB = CreateHisto2D("fhVtxXYBB","XY SPD vertex (bb)",200,-2,2,200,-2,2);
  fhVtxZBB = CreateHisto1D("fhVtxZBB","Z SPD vertex (bb)",400,-50,50);
  fhVtxXYBGA = CreateHisto2D("fhVtxXYBGA","XY SPD vertex (bga)",200,-2,2,200,-2,2);
  fhVtxZBGA = CreateHisto1D("fhVtxZBGA","Z SPD vertex (bga)",400,-50,50);
  fhVtxXYBGC = CreateHisto2D("fhVtxXYBGC","XY SPD vertex (bgc)",200,-2,2,200,-2,2);
  fhVtxZBGC = CreateHisto1D("fhVtxZBGC","Z SPD vertex (bgc)",400,-50,50);

  fListOfHistos->Add(fhAdcNoTimeA);
  fListOfHistos->Add(fhAdcWithTimeA);
  fListOfHistos->Add(fhAdcNoTimeC);
  fListOfHistos->Add(fhAdcWithTimeC);

  fListOfHistos->Add(fhAdcPMTNoTime);
  fListOfHistos->Add(fhAdcPMTWithTime);
 
  fListOfHistos->Add(fhTimeA);
  fListOfHistos->Add(fhTimeC);

  fListOfHistos->Add(fhWidthA);
  fListOfHistos->Add(fhWidthC);

  fListOfHistos->Add(fhTimePMT);
  fListOfHistos->Add(fhWidthPMT);

  fListOfHistos->Add(fhAdcWidthA);
  fListOfHistos->Add(fhAdcWidthC);

  fListOfHistos->Add(fhTimeCorr);

  fListOfHistos->Add(fhAdcTimeA);
  fListOfHistos->Add(fhAdcTimeC);

  fListOfHistos->Add(fV0a);
  fListOfHistos->Add(fV0c);
  fListOfHistos->Add(fV0multA);
  fListOfHistos->Add(fV0multC);
  fListOfHistos->Add(fV0ampl);

  fListOfHistos->Add(fhTimePMTCorr);
  fListOfHistos->Add(fhEvents);

  fListOfHistos->Add(fhVtxXYBB);
  fListOfHistos->Add(fhVtxZBB);
  fListOfHistos->Add(fhVtxXYBGA);
  fListOfHistos->Add(fhVtxZBGA);
  fListOfHistos->Add(fhVtxXYBGC);
  fListOfHistos->Add(fhVtxZBGC);
}

void AliAnaVZEROQA::UserExec(Option_t */*option*/)
{
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  AliESDVZERO* esdV0 = esd->GetVZEROData();

  Float_t timeA = 0,timeC = 0;
  Int_t ntimeA = 0, ntimeC = 0;
  for (Int_t i=0; i<64; ++i) {
      if (esdV0->GetTime(i) < 1e-6) {
	if (i >= 32) {
	  fhAdcNoTimeA->Fill(esdV0->GetAdc(i));
	}
	else {
	  fhAdcNoTimeC->Fill(esdV0->GetAdc(i));
	}
	fhAdcPMTNoTime->Fill(i,esdV0->GetAdc(i));
      }
      else {
	if (i >= 32) {
	  fhAdcWithTimeA->Fill(esdV0->GetAdc(i));
	}
	else {
	  fhAdcWithTimeC->Fill(esdV0->GetAdc(i));
	}
	fhAdcPMTWithTime->Fill(i,esdV0->GetAdc(i));
      }

      if (i >= 32) {
	fhTimeA->Fill(esdV0->GetTime(i));
	fhWidthA->Fill(esdV0->GetWidth(i));
	fhAdcWidthA->Fill(esdV0->GetAdc(i),esdV0->GetTime(i));
	fhAdcTimeA->Fill(esdV0->GetTime(i),esdV0->GetAdc(i));
      }
      else {
	fhTimeC->Fill(esdV0->GetTime(i));
	fhWidthC->Fill(esdV0->GetWidth(i));
	fhAdcWidthC->Fill(esdV0->GetAdc(i),esdV0->GetTime(i));
	fhAdcTimeC->Fill(esdV0->GetTime(i),esdV0->GetAdc(i));
       }
      fhTimePMT->Fill(i,esdV0->GetTime(i));
      fhWidthPMT->Fill(i,esdV0->GetWidth(i));

      Float_t correctedTime = CorrectLeadingTime(i,esdV0->GetTime(i),esdV0->GetAdc(i));
      fhTimePMTCorr->Fill(i,correctedTime);

      if (esdV0->GetTime(i) > 1e-6) {
	if (i >= 32) {
	  timeA += correctedTime;
	  ntimeA++;
	}
	else {
	  timeC += correctedTime;
	  ntimeC++;
	}
      }
  }

  if (ntimeA > 0) timeA = timeA/ntimeA;
  if (ntimeC > 0) timeC = timeC/ntimeC;

  fhTimeCorr->Fill(timeA,timeC);

  Int_t flaga = 0, flagc = 0;
  TString stra("emptyA"),strc("emptyC");
  if (timeA > 48 && timeA < 62) { flaga = 1; stra = "BBA"; }
  if (timeA > 26  && timeA < 33) { flaga = 2; stra = "BGA"; }
  if (timeC > 49 && timeC < 60) { flagc = 1; strc = "BBC"; }
  if (timeC > 43 && timeC < 48.5) { flagc = 2; strc = "BGC"; }

  fhEvents->Fill(flaga,flagc);

  const AliESDVertex *vtx = esd->GetPrimaryVertexSPD();

  if (flaga <= 1 && flagc <=1) {
    fhVtxXYBB->Fill(vtx->GetXv(),vtx->GetYv());
    fhVtxZBB->Fill(vtx->GetZv());
  }
  else {
    if (flaga == 2) {
      fhVtxXYBGA->Fill(vtx->GetXv(),vtx->GetYv());
      fhVtxZBGA->Fill(vtx->GetZv());
    }
    if (flagc == 2) {
      fhVtxXYBGC->Fill(vtx->GetXv(),vtx->GetYv());
      fhVtxZBGC->Fill(vtx->GetZv());
    }
  }

  fV0a->Fill(esdV0->GetNbPMV0A());
  fV0c->Fill(esdV0->GetNbPMV0C());
  fV0multA->Fill(esdV0->GetMTotV0A());
  fV0multC->Fill(esdV0->GetMTotV0C());
  for(Int_t i = 0; i < 64; i++) {
    fV0ampl->Fill(esdV0->GetMultiplicity(i));
  }

  // Post output data.
  PostData(1, fListOfHistos);
}

void AliAnaVZEROQA::Terminate(Option_t *)
{
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }
	
  fhAdcNoTimeA = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  fhAdcWithTimeA = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  fhAdcNoTimeC = dynamic_cast<TH1F*>(fListOfHistos->At(2));
  fhAdcWithTimeC = dynamic_cast<TH1F*>(fListOfHistos->At(3));

  fhAdcPMTNoTime = dynamic_cast<TH2F*>(fListOfHistos->At(4));
  fhAdcPMTWithTime = dynamic_cast<TH2F*>(fListOfHistos->At(5));

  fhTimeA = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  fhTimeC = dynamic_cast<TH1F*>(fListOfHistos->At(7));

  fhWidthA = dynamic_cast<TH1F*>(fListOfHistos->At(8));
  fhWidthC = dynamic_cast<TH1F*>(fListOfHistos->At(9));

  fhTimePMT = dynamic_cast<TH2F*>(fListOfHistos->At(10));
  fhWidthPMT = dynamic_cast<TH2F*>(fListOfHistos->At(11));

  fhAdcWidthA = dynamic_cast<TH2F*>(fListOfHistos->At(12));
  fhAdcWidthC = dynamic_cast<TH2F*>(fListOfHistos->At(13));

  fhTimeCorr = dynamic_cast<TH2F*>(fListOfHistos->At(14));

  fhAdcTimeA = dynamic_cast<TH2F*>(fListOfHistos->At(15));
  fhAdcTimeC = dynamic_cast<TH2F*>(fListOfHistos->At(16));

  fV0a = dynamic_cast<TH1F*>(fListOfHistos->At(17));
  fV0c = dynamic_cast<TH1F*>(fListOfHistos->At(18));
  fV0multA = dynamic_cast<TH1F*>(fListOfHistos->At(19));
  fV0multC = dynamic_cast<TH1F*>(fListOfHistos->At(20));
  fV0ampl  = dynamic_cast<TH1F*>(fListOfHistos->At(21));

  fhTimePMTCorr = dynamic_cast<TH2F*>(fListOfHistos->At(22));
  fhEvents = dynamic_cast<TH2F*>(fListOfHistos->At(23));

  fhVtxXYBB = dynamic_cast<TH2F*>(fListOfHistos->At(24));
  fhVtxZBB = dynamic_cast<TH1F*>(fListOfHistos->At(25));
  fhVtxXYBGA = dynamic_cast<TH2F*>(fListOfHistos->At(26));
  fhVtxZBGA = dynamic_cast<TH1F*>(fListOfHistos->At(27));
  fhVtxXYBGC = dynamic_cast<TH2F*>(fListOfHistos->At(28));
  fhVtxZBGC = dynamic_cast<TH1F*>(fListOfHistos->At(29));

 // draw the histograms if not in batch mode
  if (!gROOT->IsBatch()) {
    new TCanvas;
    fhTimePMT->DrawCopy();
    new TCanvas;
    fhAdcTimeA->DrawCopy();
    new TCanvas;
    fhAdcTimeC->DrawCopy();
    new TCanvas;
    fhAdcPMTNoTime->DrawCopy();
    new TCanvas;
    fhAdcPMTWithTime->DrawCopy();
    new TCanvas;
    fhTimeCorr->DrawCopy("E");
    new TCanvas;
    fV0ampl->DrawCopy("E");
    new TCanvas;
    fhTimePMTCorr->DrawCopy("colz");
    new TCanvas;
    fhEvents->DrawCopy("colz");
  }

  // write the output histograms to a file
  TFile* outputFile = TFile::Open("VZEROQA.root", "recreate");
  if (!outputFile || !outputFile->IsOpen())
    {
      Error("AliAnaVZEROQA", "opening output file VZEROQA.root failed");
      return;
    }

  fhAdcNoTimeA->Write();
  fhAdcWithTimeA->Write();
  fhAdcNoTimeC->Write();
  fhAdcWithTimeC->Write();

  fhAdcPMTNoTime->Write();
  fhAdcPMTWithTime->Write();
 
  fhTimeA->Write();
  fhTimeC->Write();

  fhWidthA->Write();
  fhWidthC->Write();

  fhTimePMT->Write();
  fhWidthPMT->Write();

  fhAdcWidthA->Write();
  fhAdcWidthC->Write();

  fhTimeCorr->Write();

  fhAdcTimeA->Write();
  fhAdcTimeC->Write();

  fV0a->Write();
  fV0c->Write();
  fV0multA->Write();
  fV0multC->Write();
  fV0ampl->Write();

  fhTimePMTCorr->Write();
  fhEvents->Write();

  fhVtxXYBB->Write();
  fhVtxZBB->Write();
  fhVtxXYBGA->Write();
  fhVtxZBGA->Write();
  fhVtxXYBGC->Write();
  fhVtxZBGC->Write();

  outputFile->Close();
  delete outputFile;

  //delete esd;
  Info("AliAnaVZEROQA", "Successfully finished");
}

Float_t AliAnaVZEROQA::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc)
{
  // Correct for slewing and align the channels

  if (time == 0) return 0;

  // Time offsets between channels
  Float_t timeShift[64] = {30.2914 , 30.0019 , 30.7429 , 30.1997 , 30.1511 , 29.6437 , 30.0609 , 29.5452 , 30.1437 , 30.745 , 30.7537 , 30.446 , 30.2771 , 30.838 , 30.3748 , 30.0635 , 30.1786 , 30.282 , 31.0992 , 30.7491 , 30.624 , 30.9268 , 30.6585 , 30.4895 , 31.5815 , 31.3871 , 31.2032 , 31.5778 , 31.0838 , 31.2259 , 31.2122 , 31.5989 , 28.3792 , 28.8325 , 27.8719 , 28.3475 , 26.9925 , 27.9300 , 28.4223 , 28.4996 , 28.2934 , 28.1281 , 27.209 , 28.5327 , 28.1181 , 28.0888 , 29.5111 , 28.6601 , 29.7705 , 29.6531 , 30.3373 , 30.2345 , 30.5935 , 29.8164 , 30.2235 , 29.6505 , 30.1225 , 31.2045 , 30.8399 , 30.6789 , 30.2784 , 31.7028 , 31.4239 , 30.1814};
  time -= timeShift[i];

  // Slewing correction
  if (adc == 0) return time;

  time += 30.;
  if (adc > 300.) adc = 300.;
  if (adc > 70.) {
    return (time -
	    2.93028e+01 +
	    adc*1.25188e-02 -
	    adc*adc*2.68348e-05);
  }
  else {
    return (time -
	    3.52314e+01 +
	    adc*5.99289e-01 -
	    adc*adc*2.74668e-02 +
	    adc*adc*adc*6.61224e-04 -
	    adc*adc*adc*adc*7.77105e-06 +
	    adc*adc*adc*adc*adc*3.51229e-08);
  }
}





























































