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
  fListOfHistos->SetOwner();

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

  fListOfHistos->Add(fhEvents);

  fListOfHistos->Add(fhVtxXYBB);
  fListOfHistos->Add(fhVtxZBB);
  fListOfHistos->Add(fhVtxXYBGA);
  fListOfHistos->Add(fhVtxZBGA);
  fListOfHistos->Add(fhVtxXYBGC);
  fListOfHistos->Add(fhVtxZBGC);

  PostData(1, fListOfHistos);
}

void AliAnaVZEROQA::UserExec(Option_t */*option*/)
{
  // Fill the QA histograms
  // using ESD data
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
    Printf("ERROR: No ESD event");
    return;
  }
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  if (!esdV0) {
    Printf("ERROR: No ESD VZERO");
    return;
  }

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
	fhAdcWidthA->Fill(esdV0->GetAdc(i),esdV0->GetWidth(i));
	fhAdcTimeA->Fill(esdV0->GetTime(i),esdV0->GetAdc(i));
      }
      else {
	fhTimeC->Fill(esdV0->GetTime(i));
	fhWidthC->Fill(esdV0->GetWidth(i));
	fhAdcWidthC->Fill(esdV0->GetAdc(i),esdV0->GetWidth(i));
	fhAdcTimeC->Fill(esdV0->GetTime(i),esdV0->GetAdc(i));
       }
      fhTimePMT->Fill(i,esdV0->GetTime(i));
      fhWidthPMT->Fill(i,esdV0->GetWidth(i));

  }

  fhTimeCorr->Fill(esdV0->GetV0ATime(),esdV0->GetV0CTime());

  AliESDVZERO::Decision flaga = esdV0->GetV0ADecision();
  AliESDVZERO::Decision flagc = esdV0->GetV0CDecision();

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
  // Store the output histograms
  // from the list
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }
	
  Info("AliAnaVZEROQA", "Successfully finished");
}
