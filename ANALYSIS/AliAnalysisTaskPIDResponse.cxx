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

/* $Id: AliAnalysisTaskPIDResponse.cxx 43811 2010-09-23 14:13:31Z wiechula $ */
#include <TList.h>
#include <TVectorD.h>
#include <TObjArray.h>
#include <TH2.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <TChain.h>

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

#include "AliAnalysisTaskPIDResponse.h"


ClassImp(AliAnalysisTaskPIDResponse)

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse():
AliAnalysisTaskSE(),
fIsMC(kFALSE),
fTOFTimeZeroTypeUser(-1),
fTOFTimeZeroType(AliPIDResponse::kBest_T0),
fTOFres(100.),
fPIDResponse(0x0),
fListQA(0x0),
fListQAits(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fRecoPass(0),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0)
{
  //
  // Dummy constructor
  //
}

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse(const char* name):
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fTOFTimeZeroTypeUser(-1),
fTOFTimeZeroType(AliPIDResponse::kBest_T0),
fTOFres(100.),
fPIDResponse(0x0),
fListQA(0x0),
fListQAits(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fRecoPass(0),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::~AliAnalysisTaskPIDResponse()
{
  //
  // Destructor
  //

  delete fArrPidResponseMaster;
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::ExecNewRun()
{
  //
  // Things to Execute upon a new run 
  //

  SetRecoInfo();
    
  SetITSParametrisation();

  SetTPCPidResponseMaster();
  SetTPCParametrisation();

  fPIDResponse->GetTOFResponse().SetTimeResolution(fTOFres);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::UserCreateOutputObjects()
{
  //
  // Create the output QA objects
  //

  AliLog::SetClassDebugLevel("AliAnalysisTaskPIDResponse",10);

  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");

  //pid response object
  inputHandler->CreatePIDResponse(fIsMC);
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");
  
  //
  fListQA=new TList;
  fListQAits=new TList;
  fListQAits->SetName("ITS");
  fListQAtpc=new TList;
  fListQAtpc->SetName("TPC");
  fListQAtrd=new TList;
  fListQAtrd->SetName("TRD");
  fListQAtof=new TList;
  fListQAtof->SetName("TOF");
  
  fListQA->Add(fListQAits);
  fListQA->Add(fListQAtpc);
  fListQA->Add(fListQAtrd);
  fListQA->Add(fListQAtof);

  SetupTTSqa();
  SetupTPCqa();
  SetupTRDqa();
  SetupTOFqa();

  PostData(1,fListQA);
}


//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::UserExec(Option_t */*option*/)
{
  //
  // Setup the PID response functions and fill the QA histograms
  //

  AliVEvent *event=InputEvent();
  if (!event) return;
  fRun=event->GetRunNumber();

  if (fRun!=fOldRun){
    ExecNewRun();
    fOldRun=fRun;
  }

  Double_t timeZeroType=fTOFTimeZeroTypeUser;
  if (timeZeroType<0) timeZeroType=fTOFTimeZeroType;
  fPIDResponse->SetTOFResponse(event, (AliPIDResponse::EStartTimeType_t)timeZeroType);
  
  FillITSqa();
  FillTPCqa();
  FillTOFqa();

  PostData(1,fListQA);
}


//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fRecoPass=0;
  fLHCperiod="";
  fMCperiodTPC="";

  fBeamType="";
  
  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;
  
  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  }

  fBeamType="PP";
  
  //find the period by run number (UGLY, but not stored in ESD and AOD... )
  if (fRun>=114737&&fRun<=117223)      { fLHCperiod="LHC10B"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=118503&&fRun<=121040) { fLHCperiod="LHC10C"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=122195&&fRun<=126437) { fLHCperiod="LHC10D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=127719&&fRun<=130850) { fLHCperiod="LHC10E"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=133004&&fRun<=135029) { fLHCperiod="LHC10F"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=135654&&fRun<=136377) { fLHCperiod="LHC10G"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=136851&&fRun<=139517) { fLHCperiod="LHC10H"; fMCperiodTPC="LHC10H8"; fBeamType="PBPB"; }
  else if (fRun>=139699) { fLHCperiod="LHC11A"; fMCperiodTPC="LHC10F6A"; }
  
  //find pass from file name (UGLY, but not stored in ESD... )
  TString fileName(file->GetName());
  if (fileName.Contains("/pass1")) {
    fRecoPass=1;
  } else if (fileName.Contains("/pass2")) {
    fRecoPass=2;
  } else if (fileName.Contains("/pass3")) {
    fRecoPass=3;
  }
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetITSParametrisation()
{
  //
  // Set the ITS parametrisation
  //
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetTPCPidResponseMaster()
{
  //
  // Load the TPC pid response functions from the OADB
  //

  //reset the PID response functions
  delete fArrPidResponseMaster;
  fArrPidResponseMaster=0x0;
  
  TString fileName(Form("%s/COMMON/PID/data/TPCPIDResponse.root", AliAnalysisManager::GetOADBPath()));

  TFile f(fileName.Data());
  if (f.IsOpen() && !f.IsZombie()){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(f.Get("TPCPIDResponse"));
    f.Close();
  }

  if (!fArrPidResponseMaster){
    AliFatal("Could not retrieve the TPC pid response");
    return;
  }
  fArrPidResponseMaster->SetOwner();
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetTPCParametrisation()
{
  //
  // Change BB parametrisation for current run
  //
  
  if (fLHCperiod.IsNull()) {
    AliFatal("No period set, not changing parametrisation");
    return;
  }
  
  //
  // Set default parametrisations for data and MC
  //
  
  //data type
  TString datatype="DATA";
  //in case of mc fRecoPass is per default 1
  if (fIsMC) {
    datatype="MC";
    fRecoPass=1;
  }
  
  //
  //set the PID splines
  //
  TString period=fLHCperiod;
  if (fArrPidResponseMaster){
    TObject *grAll=0x0;
    //for MC don't use period information
//     if (fIsMC) period="[A-Z0-9]*";
    //for MC use MC period information
    if (fIsMC) period=fMCperiodTPC;
//pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
    
    //loop over entries and filter them
    for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp){
      TObject *responseFunction=fArrPidResponseMaster->At(iresp);
      if (responseFunction==0x0) continue;
      TString responseName=responseFunction->GetName();
      
      if (!reg.MatchB(responseName)) continue;
      
      TObjArray *arr=reg.MatchS(responseName);
      TString particleName=arr->At(1)->GetName();
      delete arr;
      if (particleName.IsNull()) continue;
      if (particleName=="ALL") grAll=responseFunction;
      else {
        //find particle id
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
          TString particle=AliPID::ParticleName(ispec);
          particle.ToUpper();
          if ( particle == particleName ){
            //test if there is already a function set. If yes, cleanup
            TObject *old=const_cast<TObject*>(fPIDResponse->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec));
            if (old) delete old;
            fPIDResponse->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
            fPIDResponse->GetTPCResponse().SetUseDatabase(kTRUE);
            AliInfo(Form("Adding graph: %d - %s",ispec,responseFunction->GetName()));
            break;
          }
        }
      }
    }
    
    //set default response function to all particles which don't have a specific one
    if (grAll){
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
        if (!fPIDResponse->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec)){
          fPIDResponse->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,grAll);
          AliInfo(Form("Adding graph: %d - %s",ispec,grAll->GetName()));
        }
      }
    }
  }

  //
  // Setup resolution parametrisation
  //

  //default
  fPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);

  if (fRun>=122195){
    fPIDResponse->GetTPCResponse().SetSigma(2.30176e-02, 5.60422e+02);
  }
//   if ( fBeamType == "PBPB" ){
//     Double_t corrSigma=GetMultiplicityCorrectionSigma(GetTPCMultiplicityBin());
//     fPIDResponse->GetTPCResponse().SetSigma(3.79301e-03*corrSigma, 2.21280e+04);
//   }
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::FillITSqa()
{
  AliVEvent *event=InputEvent();
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // ITS refit + ITS pid
    if (!( (status&0x0004==0x0004) && (status&0x0008==0x0008) )) return;
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
void AliAnalysisTaskPIDResponse::FillTPCqa()
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
    if (!(status&0x0040==0x0040) || !(status&0x0004==0x0004) ||  !(status&0x0080==0x0080) ) return;
    
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
    
    TH2 *hPt=(TH2*)fListQAtpc->At(AliPID::kSPECIES+1);
    if (hPt) {
      Double_t sig=track->GetTPCsignal();
      hPt->Fill(sig);
    }
    
    TH2 *hMom=(TH2*)fListQAtpc->At(AliPID::kSPECIES+2);
    if (hMom) {
      Double_t sig=track->GetTPCmomentum();
      hMom->Fill(sig);
    }
    
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::FillTOFqa()
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
    if (!(status&0x0040==0x0040) || !(status&0x0004==0x0004)
        || !(status&0x2000==0x2000) || !(status&0x8000==0x8000)
        || !(status&0x80000000==0x80000000) ) return;
    
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
void AliAnalysisTaskPIDResponse::SetupTTSqa()
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
void AliAnalysisTaskPIDResponse::SetupTPCqa()
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
  
  TH1F *hPt = new TH1F("hPt","Pt_TPC",100,0,300);
  fListQAtpc->Add(hPt);
  
  TH1F *hMom = new TH1F("hMom","Mom_TPC",100,0,20);
  fListQAtpc->Add(hMom);
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetupTRDqa()
{
  //
  // Create the TRD qa objects
  //
  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetupTOFqa()
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
TVectorD* AliAnalysisTaskPIDResponse::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
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
TVectorD* AliAnalysisTaskPIDResponse::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
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
TVectorD* AliAnalysisTaskPIDResponse::MakeArbitraryBinning(const char* bins)
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
