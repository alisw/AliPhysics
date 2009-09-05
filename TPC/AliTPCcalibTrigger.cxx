
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

/*
  // Load libraries
  gSystem->Load("libANALYSIS");
    gSystem->Load("libTPCcalib");
 

    .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

  TFile f("CalibObjects.root");
  AliTPCcalibTrigger *calibTrigger = (AliTPCcalibTrigger *)f->Get("TPCCalib")->FindObject("calibTrigger");


*/

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibTrigger.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"

ClassImp(AliTPCcalibTrigger)

AliTPCcalibTrigger::AliTPCcalibTrigger():
  AliTPCcalibBase("calibTrigger","calibTrigger"),
  fHisMap(0)
{

}

AliTPCcalibTrigger::AliTPCcalibTrigger(const char * name, const char * title):
  AliTPCcalibBase(name,title),
  fHisMap(0)
{
  //
  //
  //
  fHisMap = new TMap;
}

Long64_t AliTPCcalibTrigger::Merge(TCollection *li) {
  //
  // Merge histograms
  //
  TIterator* iter = li->MakeIterator();
  AliTPCcalibTrigger* cal = 0;

  while ((cal = (AliTPCcalibTrigger*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibTrigger::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    TMap * addMap=(cal->fHisMap);
    if(!addMap) return 0;
    TIterator* iterator = addMap->MakeIterator();
    iterator->Reset();
    TPair* addPair=0;
    while((addPair=(TPair *)(addMap->FindObject(iterator->Next())))){
      THnSparse* addHist=dynamic_cast<THnSparseF*>(addPair->Value());
      if (!addHist) continue;
      addHist->Print();
      THnSparse* localHist=dynamic_cast<THnSparseF*>(fHisMap->GetValue(addHist->GetName()));
      if(!localHist){
        localHist=MakeHisto(addHist->GetName());
        fHisMap->Add(new TObjString(addHist->GetName()),localHist);
      }
      localHist->Add(addHist);
    }
  }
  return 0;
}



void AliTPCcalibTrigger::Process(AliESDEvent *event){
  //
  //
  //
  if (!event) return;
  const TString &trigger = event->GetFiredTriggerClasses();
  //
  if (!GetHisto(trigger.Data())){
    AddHisto(trigger.Data(),MakeHisto(trigger.Data()));
  }
  if (!GetHisto("all")){
    AddHisto("all",MakeHisto("all"));
  }

  THnSparse *histoAll = GetHisto("all");
  THnSparse *histo = GetHisto(trigger.Data());
  Double_t xcont[8]={0,0,0,0,0,0,0,0};
  
  Int_t ntracks = event->GetNumberOfTracks();
  xcont[0] = ntracks;
  //
  // GetLongest track
  //  
  AliESDtrack * longest=0;
  Int_t nclmax=0;
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (!track) continue;
    if (track->GetTPCNcls()<=nclmax) continue;
    nclmax = track->GetTPCNcls();
    longest= track;
  }
  //
  // get inof of longest track
  /*TString  axisName[8]={
    "ntracks",
    "nclMax",
    "dcaR",
    "dcaZ",
    "alpha",
    "theta",
    "pt",
    "dEdx"
  };
  */
  if (longest){
    Float_t dca[2];
    Double_t pxyz[3];
    longest->GetDZ(0.,0.,0.,event->GetMagneticField(),dca);
    Bool_t status = longest->GetPxPyPz(pxyz);
    xcont[1]=nclmax;
    xcont[2]=dca[0];
    xcont[3]=dca[1];
    xcont[4]=TMath::ATan2(pxyz[1],pxyz[0]);
    xcont[5]=longest->GetParameter()[3];
    xcont[6]=longest->Pt();
    xcont[7]=longest->GetTPCsignal();    
  }
  //
  histoAll->Fill(xcont);
  histo->Fill(xcont);
  
}

THnSparse *AliTPCcalibTrigger::MakeHisto(const char* trigger){
  //
  // Make event/track histograms
  // trigger histo 
  //
  //                 ntracks  nclMax dcaR dcaZ alpha   theta pt dEdx
  Int_t    bins[8] = {50,     20,    20,  20,  18,     25,   25, 25 };
  //Int_t    bins[8] = {50*     20*  25*  25*  18*     25*   25* 25 };
  Double_t xmin[8] = {0.,     0,      0,  0,   -3.14, -1.5,  0, 0};
  Double_t xmax[8] = {50,     160,  150,  250, 3.14,  1.5,   100, 100};
  TString  axisName[8]={
    "ntracks",
    "nclMax",
    "dcaR",
    "dcaZ",
    "alpha",
    "theta",
    "pt",
    "dEdx"
  };
  TString  axisTitle[8]={
    "Number of tracks",
    "N_{cl}",
    "dca_{R} (cm)",
    "dca_{z} (cm)",
    "alpha (mrad)",
    "theta",
    "p_{t} (GeV/c)",
    "dEdx (a.u.)"
  };

  
  THnSparse *sparse = new THnSparseF(Form("his_%s",trigger), Form("his_%s",trigger), 8, bins, xmin, xmax);
  for (Int_t iaxis=0; iaxis<8; iaxis++){
    sparse->GetAxis(iaxis)->SetName(axisName[iaxis]);
    sparse->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
  }
  return sparse;
}

THnSparse * AliTPCcalibTrigger::GetHisto(const char *trigger) { 
  //
  // return histogram for given class
  if (!fHisMap) fHisMap=new TMap;
  return (THnSparse*) fHisMap->GetValue(trigger);
}

void   AliTPCcalibTrigger::AddHisto(const char *trigger, THnSparse *his) { 
  if (!GetHisto(trigger)) {
    TObjString *pstr = new TObjString(trigger);
    fHisMap->Add(pstr,his);
  }
}
