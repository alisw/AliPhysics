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

/* $Id: AliAnalysisTaskPIDqa.cxx 43811 2010-09-23 14:13:31Z wiechula $ */
#include <TList.h>
#include <TVectorD.h>
#include <TObjArray.h>
#include <TH2.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <TChain.h>
#include <TF1.h>
#include <TSpline.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliVTrack.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliITSPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <AliTRDPIDResponse.h>
#include <AliTOFPIDResponse.h>

#include "AliAnalysisTaskPIDqa.h"


ClassImp(AliAnalysisTaskPIDqa)

//______________________________________________________________________________
AliAnalysisTaskPIDqa::AliAnalysisTaskPIDqa():
AliAnalysisTaskSE(),
fPIDResponse(0x0),
fListQA(0x0),
fListQAits(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0)
{
  //
  // Dummy constructor
  //
}

//______________________________________________________________________________
AliAnalysisTaskPIDqa::AliAnalysisTaskPIDqa(const char* name):
AliAnalysisTaskSE(name),
fPIDResponse(0x0),
fListQA(0x0),
fListQAits(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//______________________________________________________________________________
AliAnalysisTaskPIDqa::~AliAnalysisTaskPIDqa()
{
  //
  // Destructor
  //

}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::UserCreateOutputObjects()
{
  //
  // Create the output QA objects
  //

  AliLog::SetClassDebugLevel("AliAnalysisTaskPIDqa",10);

  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliError("PIDResponse object was not created");
  
  //
  fListQA=new TList;
  fListQA->SetOwner();
  
  fListQAits=new TList;
  fListQAits->SetOwner();
  fListQAits->SetName("ITS");
  
  fListQAtpc=new TList;
  fListQAtpc->SetOwner();
  fListQAtpc->SetName("TPC");
  
  fListQAtrd=new TList;
  fListQAtrd->SetOwner();
  fListQAtrd->SetName("TRD");
  
  fListQAtof=new TList;
  fListQAtof->SetOwner();
  fListQAtof->SetName("TOF");
  
  fListQA->Add(fListQAits);
  fListQA->Add(fListQAtpc);
  fListQA->Add(fListQAtrd);
  fListQA->Add(fListQAtof);

  SetupITSqa();
  SetupTPCqa();
  SetupTRDqa();
  SetupTOFqa();

  PostData(1,fListQA);
}


//______________________________________________________________________________
void AliAnalysisTaskPIDqa::UserExec(Option_t */*option*/)
{
  //
  // Setup the PID response functions and fill the QA histograms
  //

  AliVEvent *event=InputEvent();
  if (!event||!fPIDResponse) return;

  
  FillITSqa();
  FillTPCqa();
  FillTOFqa();

  PostData(1,fListQA);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillITSqa()
{
  AliVEvent *event=InputEvent();
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // ITS refit + ITS pid
    if (!( ( (status&0x0004)==0x0004 ) && ( (status&0x0008)==0x0008 ) )) return;
    Double_t mom=track->P();
    
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      TH2 *h=(TH2*)fListQAits->At(ispecie);
      if (!h) continue;
      Double_t nSigma=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType)ispecie);
      h->Fill(mom,nSigma);
    }
    
    TH2 *h=(TH2*)fListQAits->At(AliPID::kSPECIES);
    if (h) {
      Double_t sig=track->GetITSsignal();
      h->Fill(mom,sig);
    }
    
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillTPCqa()
{
  AliVEvent *event=InputEvent();
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);
    
    //
    //basic track cuts
    //
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // TPC refit + ITS refit + TPC pid
    if (!( (status&0x0040)==0x0040) || !( (status&0x0004)==0x0004) ||  !((status&0x0080)==0x0080) ) return;
    
    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }
    
    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;
    
    Double_t mom=track->GetTPCmomentum();
    
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      TH2 *h=(TH2*)fListQAtpc->At(ispecie);
      if (!h) continue;
      Double_t nSigma=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)ispecie);
      h->Fill(mom,nSigma);
    }
    
    TH2 *h=(TH2*)fListQAtpc->At(AliPID::kSPECIES);
    if (h) {
      Double_t sig=track->GetTPCsignal();
      h->Fill(mom,sig);
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillTOFqa()
{
  AliVEvent *event=InputEvent();
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);
    
    //
    //basic track cuts
    //
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // TPC refit + ITS refit +
    // TOF out + TOFpid +
    // kTIME
    if (!((status&0x0040)==0x0040) || !((status&0x0004)==0x0004)
        || !((status&0x2000)==0x2000) || !((status&0x8000)==0x8000)
        || !((status&0x80000000)==0x80000000) ) return;
    
    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }
    
    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;
    
    
    Double_t mom=track->P();
    
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      TH2 *h=(TH2*)fListQAtof->At(ispecie);
      if (!h) continue;
      Double_t nSigma=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)ispecie);
      h->Fill(mom,nSigma);
    }
    
    TH2 *h=(TH2*)fListQAtof->At(AliPID::kSPECIES);
    if (h) {
      Double_t sig=track->GetTOFsignal();
      h->Fill(mom,sig);
    }
    
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupITSqa()
{
  //
  // Create the ITS qa objects
  //
  
  TVectorD *vX=MakeLogBinning(200,.1,30);
  
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaP = new TH2F(Form("hNsigmaP_ITS_%s",AliPID::ParticleName(ispecie)),
                              Form("ITS n#sigma %s vs. p;p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1,vX->GetMatrixArray(),
                              200,-10,10);
    fListQAits->Add(hNsigmaP);
  }
  
  
  TH2F *hSig = new TH2F("hSigP_ITS",
                        "ITS signal vs. p;p [GeV]; ITS signal [arb. units]",
                        vX->GetNrows()-1,vX->GetMatrixArray(),
                        300,0,300);
  
  fListQAits->Add(hSig);
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupTPCqa()
{
  //
  // Create the TPC qa objects
  //
  
  TVectorD *vX=MakeLogBinning(200,.1,30);
  
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaP = new TH2F(Form("hNsigmaP_TPC_%s",AliPID::ParticleName(ispecie)),
                              Form("TPC n#sigma %s vs. p;p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1,vX->GetMatrixArray(),
                              200,-10,10);
    fListQAtpc->Add(hNsigmaP);
  }
  
  
  TH2F *hSig = new TH2F("hSigP_TPC",
                        "TPC signal vs. p;p [GeV]; TPC signal [arb. units]",
                        vX->GetNrows()-1,vX->GetMatrixArray(),
                        300,0,300);
  
  fListQAtpc->Add(hSig);
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupTRDqa()
{
  //
  // Create the TRD qa objects
  //
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupTOFqa()
{
  //
  // Create the TOF qa objects
  //
  
  TVectorD *vX=MakeLogBinning(200,.1,30);
  
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaP = new TH2F(Form("hNsigmaP_TOF_%s",AliPID::ParticleName(ispecie)),
                              Form("TOF n#sigma %s vs. p;p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1,vX->GetMatrixArray(),
                              200,-10,10);
    fListQAtof->Add(hNsigmaP);
  }
  
  
  TH2F *hSig = new TH2F("hSigP_TOF",
                        "TOF signal vs. p;p [GeV]; TOF signal [arb. units]",
                        vX->GetNrows()-1,vX->GetMatrixArray(),
                        300,0,300);
  
  fListQAtof->Add(hSig);
  
}

//______________________________________________________________________________
TVectorD* AliAnalysisTaskPIDqa::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    AliError("For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//______________________________________________________________________________
TVectorD* AliAnalysisTaskPIDqa::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make linear binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliAnalysisTaskPIDqa::MakeArbitraryBinning(const char* bins)
{
  //
  // Make arbitrary binning, bins separated by a ','
  //
  TString limits(bins);
  if (limits.IsNull()){
    AliError("Bin Limit string is empty, cannot add the variable");
    return 0x0;
  }
  
  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    AliError("Need at leas 2 bin limits, cannot add the variable");
    delete arr;
    return 0x0;
  }
  
  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }
  
  delete arr;
  return binLimits;
}
