
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliSpectraBothHistoManager class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraBothHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraBothHistoManager)


using namespace AliSpectraNameSpaceBoth;
#include "HistogramNamesBoth.cxx" // generate this automatically running createNames.py 


//const char* kParticleSpecies[] =
 // {
   // "PionPlus",
  //  "KaonPlus",
  //  "ProtonPlus",
  //  "PionMinus",
  //  "KaonMinus",
  //  "ProtonMinus",
 // };


AliSpectraBothHistoManager::AliSpectraBothHistoManager(const char *name,Int_t nrebin): TNamed(name, "AOD Spectra Histo Manager"), fOutputList(0), fNRebin(0)
{
  // ctor
  fNRebin=nrebin;
  fOutputList = new TList;
  fOutputList->SetOwner();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  for (Int_t ihist  = 0; ihist < kNHist ; ihist++)
    {
      if (ihist <= kNPtGenHist) BookPtGenHistogram(kHistNameBoth[ihist]);  // PT histos
      if (ihist > kNPtGenHist && ihist <= kNPtGenAllChHist) BookPtGenAllChHistogram(kHistNameBoth[ihist]);  // PT histos
      if (ihist > kNPtGenAllChHist && ihist <= kNPtRecHist) BookPtRecHistogram(kHistNameBoth[ihist]);  // PT histos
      if (ihist > kNPtRecHist && ihist <= kNPtRecAllChHist) BookPtRecAllChHistogram(kHistNameBoth[ihist]);  // PT histos
      if (ihist > kNPtRecAllChHist && ihist <= kNHistPID) BookPIDHistogram(kHistNameBoth[ihist]);  // PID histos
      if (ihist > kNHistPID && ihist <= kNHistNSig) BookNSigHistogram(kHistNameBoth[ihist]);  // NSigmaSep histos
      if(ihist==kHistGenMulvsRawMul) BookGenMulvsRawMulHistogram(kHistNameBoth[ihist]); 
    }
   
  TH1::AddDirectory(oldStatus);
}

//_______________________________________________________

TH2F* AliSpectraBothHistoManager::BookPtGenHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking pt gen histogram %s", name));
   
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nbinsTempl=52;
   
  TH2F * hist = new TH2F(name,Form("P_{T} distribution (%s)", name),nbinsTempl,templBins,2,-0.5,1.5);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
  hist->GetXaxis()->SetTitle("generated P_{T} (GeV / c)");
  hist->GetYaxis()->SetTitle("IsPhysicalPrimary()");
  hist->SetMarkerStyle(kFullCircle);
  hist->Sumw2();
  fOutputList->Add(hist);
   
  return hist;
}

//_______________________________________________________

TH2F* AliSpectraBothHistoManager::BookPtGenAllChHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking pt gen histogram - no PID %s", name));
   
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.5,6.0,6.5,7,7.5,8,8.5,9,9.5,10};
  Int_t nbinsTempl=62;
   
  TH2F * hist = new TH2F(name,Form("P_{T} distribution (All Ch) (%s)", name),nbinsTempl,templBins,2,-0.5,1.5);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
  hist->GetXaxis()->SetTitle("generated P_{T} (GeV / c)");
  hist->GetYaxis()->SetTitle("IsPhysicalPrimary()");
  hist->SetMarkerStyle(kFullCircle);
  hist->Sumw2();
  fOutputList->Add(hist);
   
  return hist;
}


//_______________________________________________________
TH2F* AliSpectraBothHistoManager::BookPtRecHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking pt rec histogram %s,  rebin:%d", name, fNRebin));
   
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nbinsTempl=52;
   
  TH2F * hist = new TH2F(name,Form("reconstructed P_{T} distribution (%s)", name),nbinsTempl,templBins,3000,-3,3);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
  if(fNRebin!=0)hist->RebinY(fNRebin);
  hist->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  hist->GetYaxis()->SetTitle("DCA xy");
  hist->SetMarkerStyle(kFullCircle);
  hist->Sumw2();
  fOutputList->Add(hist);

  return hist;
}

//_______________________________________________________
TH2F* AliSpectraBothHistoManager::BookPtRecAllChHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking pt rec histogram %s,  rebin:%d", name, fNRebin));
   
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.5,6.0,6.5,7,7.5,8,8.5,9,9.5,10};
  Int_t nbinsTempl=62;
   
  TH2F * hist = new TH2F(name,Form("reconstructed P_{T} distribution (All Ch)  (%s)", name),nbinsTempl,templBins,3000,-3,3);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
  if(fNRebin!=0)hist->RebinY(fNRebin);
  hist->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  hist->GetYaxis()->SetTitle("DCA xy");
  hist->SetMarkerStyle(kFullCircle);
  hist->Sumw2();
  fOutputList->Add(hist);

  return hist;
}

//_____________________________________________________________________________

TH2F* AliSpectraBothHistoManager::BookPIDHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking PID histogram %s, rebin:%d", name, fNRebin));
TString tmp(name);	
  TH2F * hist = new TH2F(name, Form("Particle Identification (%s)", name), 50, 0, 2.5, 200, -1000, 1000);
  if(fNRebin!=0){
    hist->RebinY(fNRebin);
  //  hist->RebinX(fNRebin);
  }
 if(tmp.Contains("Pt"))
	hist->GetXaxis()->SetTitle("P_{t} (GeV / c)");
  else		
  	hist->GetXaxis()->SetTitle("P (GeV / c)");

  hist->GetYaxis()->SetTitle("PID signal");
  //  hist->Sumw2();
  fOutputList->Add(hist);

  return hist;
}

//_____________________________________________________________________________

TH2F* AliSpectraBothHistoManager::BookNSigHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking NSigma histogram %s, rebin:%d", name, fNRebin));
  Int_t nbins=200;
  Float_t miny=-20;
  TString tmp(name);			  
  if(tmp.Contains("TPCTOF"))
  {
	nbins=100;
	miny=0.0;
  }				
  TH2F * hist = new TH2F(name, Form("Particle Identification (%s)", name), 50, 0, 2.5,nbins,miny, 20);
  if(fNRebin!=0){
    hist->RebinY(fNRebin);
    //hist->RebinX(fNRebin);
  }
  if(tmp.Contains("Pt"))
	hist->GetXaxis()->SetTitle("P_{t} (GeV / c)");
  else		
  	hist->GetXaxis()->SetTitle("P (GeV / c)");
  //hist->GetYaxis()->SetTitle("TPC");
  //hist->Sumw2();
  fOutputList->Add(hist);
  
  return hist;
}

//_____________________________________________________________________________

 TH2F*   AliSpectraBothHistoManager::BookGenMulvsRawMulHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking  Gen vs Raw multilicity histogram %s", name));
  
  TH2F * hist = new TH2F(name, Form("Gen vs Raw multilicity (%s)", name), 100, -0.5, 99.5, 100, -0.5, 99.5);
  hist->GetXaxis()->SetTitle("gen mul ");
  hist->GetYaxis()->SetTitle("raw mul ");
  //  hist->Sumw2();
  fOutputList->Add(hist);
  
  return hist;
}
//_____________________________________________________________________________

TH1F* AliSpectraBothHistoManager::GetPtHistogram1D(const char * name,Double_t minDCA,Double_t maxDCA)
{
  //   //return the projection of the TH2 (pt,DCA) in the DCA bin range [firstDCAbin,lastDCAbin]
  //   //if minDCA=-1 && maxDCA=-1 the projection is done using the full DCA range
  TH2F *hist=(TH2F*)fOutputList->FindObject(name);
  TH1F *outhist=0x0;
  Printf("--- Projecting %s on Xaxis[%f,%f]:",name,minDCA,maxDCA);
  if(minDCA==-1 && maxDCA==-1){
    outhist=(TH1F*)hist->ProjectionX("_px",0,-1,"e");
    Printf("Full Range");
  }else {
    Int_t firstbin=hist->GetYaxis()->FindBin(minDCA);
    Int_t lastbin=hist->GetYaxis()->FindBin(maxDCA);
    Printf("firstbin: %d lastbin: %d",firstbin,lastbin);
    outhist=(TH1F*)hist->ProjectionX("_px",firstbin,lastbin,"e");
  }
  Printf("Entries outhist: %.0f   Entries hist: %.0f",outhist->GetEntries(),hist->GetEntries());
  return outhist;
}

//_____________________________________________________________________________

TH1F* AliSpectraBothHistoManager::GetDCAHistogram1D(const char * name,Double_t minPt,Double_t maxPt)
{
  //   //return the projection of the TH2 (pt,DCA) in the DCA bin range [firstDCAbin,lastDCAbin]
  //   //if minPt=-1 && maxPt=-1 the projection is done using the full DCA range
  TH2F *hist=(TH2F*)fOutputList->FindObject(name);
  TH1F *outhist=0x0;
  Printf("--- Projecting %s on Yaxis[%f,%f]:",name,minPt,maxPt);
  if(minPt==-1 && maxPt==-1){
    outhist=(TH1F*)hist->ProjectionY("_py",0,-1,"e");
    Printf("Full Range");
  }else {
    Int_t firstbin=hist->GetXaxis()->FindBin(minPt);
    Int_t lastbin=hist->GetXaxis()->FindBin(maxPt);
    Printf("firstbin: %d lastbin: %d",firstbin,lastbin);
    outhist=(TH1F*)hist->ProjectionY("_py",firstbin,lastbin,"e");
    Printf("GetDCAHistogram1D(%s) BinRange:%d  %d  Pt Range: %f %f",hist->GetName(),firstbin,lastbin,hist->GetXaxis()->GetBinLowEdge(firstbin),hist->GetXaxis()->GetBinLowEdge(firstbin)+hist->GetXaxis()->GetBinWidth(lastbin));
  }
  Printf("Entries outhist: %.0f   Entries hist: %.0f",outhist->GetEntries(),hist->GetEntries());
  return outhist;
}

//_____________________________________________________________________________

Long64_t AliSpectraBothHistoManager::Merge(TCollection* list)
{
  // Merge a list of AliSpectraBothHistoManager objects with this.
  // Returns the number of merged objects (including this).

  //  AliInfo("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraBothHistoManager* entry = dynamic_cast<AliSpectraBothHistoManager*> (obj);
    if (entry == 0) 
      continue;

    TList * hlist = entry->GetOutputList();      
    collections.Add(hlist);
    count++;
  }
  
  fOutputList->Merge(&collections);
  
  delete iter;

  return count+1;
}

//________________________________________________________________________________
TH1* AliSpectraBothHistoManager::GetHistogram1D(UInt_t histoType, UInt_t particleType, UInt_t charge) {
  // GetHistogram using particle ID and histogram type
  Int_t baseId = -1;

  if (particleType == kSpUndefined) {
    AliError ("Trying to get histo for undefined particle");
    return 0;
  }

  switch(histoType) {
  case kHistPtGenTruePrimary:
    baseId = kHistPtGenTruePrimaryPionPlus;
    break;
  case kHistPtRecSigma:
    baseId = kHistPtRecSigmaPionPlus;
    break;
  case kHistPtRecTrue:
    baseId = kHistPtRecTruePionPlus;
    break;
  case kHistPtRecTruePrimary:
    baseId = kHistPtRecTruePrimaryPionPlus;
    break;
  case kHistPtRecPrimary:
    baseId = kHistPtRecPrimaryPionPlus;
    break;
  case kHistPtRecSigmaPrimary:
    baseId = kHistPtRecSigmaPrimaryPionPlus;
    break;
  case kHistPtRecSigmaSecondaryMaterial:
    baseId = kHistPtRecSigmaSecondaryMaterialPionPlus;
    break;
  case kHistPtRecSigmaSecondaryWeakDecay:
    baseId = kHistPtRecSigmaSecondaryWeakDecayPionPlus;
    break;
  case kHistNSigTPC:
    baseId = kHistNSigPionTPC;
    break;
  case kHistNSigTOF:
    baseId = kHistNSigPionTOF;
    break;
  case kHistNSigTPCTOF:
    baseId = kHistNSigPionTPCTOF;
    break;
  default:
    baseId = -1;
  }
  
  if (baseId < 0)
    AliFatal(Form("Wrong histogram type %d", histoType));

  //cout << "T[" << histoType << "] ID["<< baseId <<"] P["<<particleType<<"] C[" << charge 
  //     << " --> ["<< baseId + particleType + 3*(charge) <<"] = " ;

  baseId = baseId + particleType + 3*(charge);

  //cout <<  GetHistogram(baseId)->GetName() << endl;

  return GetHistogram(baseId);
}
//____________________________________________________________________________________________________
TH2* AliSpectraBothHistoManager::GetHistogram2D(UInt_t histoType, UInt_t particleType, UInt_t charge){
  // returns histo based on ids, casting it to TH2*
  return (TH2*) GetHistogram1D(histoType,particleType,charge);


}
