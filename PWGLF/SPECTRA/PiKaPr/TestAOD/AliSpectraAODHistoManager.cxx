
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
//         AliSpectraAODHistoManager class
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
#include "AliSpectraAODHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODHistoManager)


using namespace AliSpectraNameSpace;

AliSpectraAODHistoManager::AliSpectraAODHistoManager(const char *name): TNamed(name, "AOD Spectra Histo Manager"), fOutputList(0)
{
  // ctor
   fOutputList = new TList;
   fOutputList->SetOwner();
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   for (Int_t ihist  = 0; ihist < kNHist ; ihist++)
   {
      if (ihist <= kNPtGenHist) BookPtGenHistogram(kHistName[ihist]);  // PT histos
      if (ihist > kNPtGenHist && ihist <= kNPtRecHist) BookPtRecHistogram(kHistName[ihist]);  // PT histos
      if (ihist > kNPtRecHist && ihist <= kNHistPID) BookPIDHistogram(kHistName[ihist]);  // PID histos
      if (ihist > kNHistPID && ihist <= kNHistNSig) BookNSigHistogram(kHistName[ihist]);  // NSigmaSep histos
      if (ihist > kNHistNSig) BookqVecHistogram(kHistName[ihist]);  // qDistr histos
   }
   //adding quickly o plot to check the centrality distr in two different ways
   TH2F *hist=new TH2F("CentCheck","CentCheck",1000,0,100,1000,0,100);
   hist->SetXTitle("fAOD->GetCentrality()->GetCentralityPercentile(V0M)");
   hist->SetYTitle("fAOD->GetHeader()->GetCentralityP()->GetCentralityPercentileUnchecked(V0M)");
   hist->Sumw2();
   fOutputList->Add(hist);
   
   TH1::AddDirectory(oldStatus);

}

//_______________________________________________________

TH2F* AliSpectraAODHistoManager::BookPtGenHistogram(const char * name)
{
   // Return a pt histogram with predefined binning, set the ID and add it to the output list
   AliInfo(Form("Booking pt histogram %s", name));
   
   //standard histo
   const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
   Int_t nbinsTempl=52;
   
   TH2F * hist = new TH2F(name,Form("P_{T} distribution (%s)", name),nbinsTempl,templBins,2,-0.5,1.5);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
   hist->GetXaxis()->SetTitle("generated P_{T} (GeV / c)");
   hist->GetYaxis()->SetTitle("DCA xy");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);
   
   return hist;
}


//_______________________________________________________
TH2F* AliSpectraAODHistoManager::BookPtRecHistogram(const char * name)
{
   // Return a pt histogram with predefined binning, set the ID and add it to the output list
   AliInfo(Form("Booking pt histogram %s", name));
   
   //standard histo
   const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
   Int_t nbinsTempl=52;
   
   TH2F * hist = new TH2F(name,Form("reconstructed P_{T} distribution (%s)", name),nbinsTempl,templBins,3000,-3,3);//need to be at least 1 becuase the generated are filled with (pt,IsPhysPrim)
   hist->GetXaxis()->SetTitle("P_{T} (GeV / c)");
   hist->GetYaxis()->SetTitle("DCA xy");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);

   return hist;
}

//_____________________________________________________________________________

TH2F* AliSpectraAODHistoManager::BookPIDHistogram(const char * name)
{
   // Return a pt histogram with predefined binning, set the ID and add it to the output list
   AliInfo(Form("Booking pt histogram %s", name));

   TH2F * hist = new TH2F(name, Form("Particle Identification (%s)", name), 100, 0, 1.5, 2000, -1000, 1000);
   hist->GetXaxis()->SetTitle("(GeV / c)");
   hist->GetYaxis()->SetTitle("PID signal");
//  hist->Sumw2();
   fOutputList->Add(hist);

   return hist;
}

//_____________________________________________________________________________

TH2F* AliSpectraAODHistoManager::BookNSigHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking pt histogram %s", name));
  
  TH2F * hist = new TH2F(name, Form("Particle Identification (%s)", name), 200, 0, 1.5, 4000,-40, 40);
  hist->GetXaxis()->SetTitle("P (GeV / c)");
  hist->GetYaxis()->SetTitle("TPC");
  //hist->Sumw2();
  fOutputList->Add(hist);
  
  return hist;
}

//_____________________________________________________________________________

TH2F* AliSpectraAODHistoManager::BookqVecHistogram(const char * name)
{
  // Return a pt histogram with predefined binning, set the ID and add it to the output list
  AliInfo(Form("Booking q Vector histogram %s", name));
  
  TH2F * hist = new TH2F(name, Form("q Vector distribution vs Centrality (%s)", name), 100, 0, 10, 100, 0, 100);
  hist->GetXaxis()->SetTitle("q vector");
  hist->GetYaxis()->SetTitle("Centrality (V0)");
  //  hist->Sumw2();
  fOutputList->Add(hist);
  
  return hist;
}


//_____________________________________________________________________________

TH1F* AliSpectraAODHistoManager::GetPtHistogram1D(const char * name,Double_t minDCA,Double_t maxDCA)
{
  //   //return the projection of the TH2 (pt,DCA) in the DCA bin range [firstDCAbin,lastDCAbin]
  //   //if minDCA=-1 && maxDCA=-1 the projection is done using the full DCA range
  TH2F *hist=(TH2F*)fOutputList->FindObject(name);
  TH1F *outhist=0x0;
  if(minDCA==-1 && maxDCA==-1)outhist=(TH1F*)hist->ProjectionX("",0,-1,"e");
  else {
    Int_t firstbin=hist->GetYaxis()->FindBin(minDCA);
    Int_t lastbin=hist->GetYaxis()->FindBin(maxDCA);
    Printf("firstbin: %d lastbin: %d",firstbin,lastbin);
    outhist=(TH1F*)hist->ProjectionX("",firstbin,lastbin,"e");
  }
  Printf("Entries outhist: %.0f   Entries hist: %.0f",outhist->GetEntries(),hist->GetEntries());
  return outhist;
}

//_____________________________________________________________________________

TH1F* AliSpectraAODHistoManager::GetDCAHistogram1D(const char * name,Double_t minPt,Double_t maxPt)
{
  //   //return the projection of the TH2 (pt,DCA) in the DCA bin range [firstDCAbin,lastDCAbin]
  //   //if minPt=-1 && maxPt=-1 the projection is done using the full DCA range
  TH2F *hist=(TH2F*)fOutputList->FindObject(name);
  TH1F *outhist=0x0;
  if(minPt==-1 && maxPt==-1)outhist=(TH1F*)hist->ProjectionY("",0,-1,"e");
  else {
    Int_t firstbin=hist->GetXaxis()->FindBin(minPt);
    Int_t lastbin=hist->GetXaxis()->FindBin(maxPt);
    outhist=(TH1F*)hist->ProjectionY("",firstbin,lastbin,"e");
  }
  return outhist;
}

//_____________________________________________________________________________

Long64_t AliSpectraAODHistoManager::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODHistoManager objects with this.
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
    AliSpectraAODHistoManager* entry = dynamic_cast<AliSpectraAODHistoManager*> (obj);
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

