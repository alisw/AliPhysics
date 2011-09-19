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

#include <AliESDEvent.h>

#include "AliAnalysisTaskPIDqa.h"


ClassImp(AliAnalysisTaskPIDqa)

//______________________________________________________________________________
AliAnalysisTaskPIDqa::AliAnalysisTaskPIDqa():
AliAnalysisTaskSE(),
fPIDResponse(0x0),
fListQA(0x0),
fListQAits(0x0),
fListQAitsSA(0x0),
fListQAitsPureSA(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0),
fListQAemcal(0x0),
fListQAtpctof(0x0)
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
fListQAitsSA(0x0),
fListQAitsPureSA(0x0),
fListQAtpc(0x0),
fListQAtrd(0x0),
fListQAtof(0x0),
fListQAemcal(0x0),
fListQAtpctof(0x0)
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

  fListQAitsSA=new TList;
  fListQAitsSA->SetOwner();
  fListQAitsSA->SetName("ITS_SA");

  fListQAitsPureSA=new TList;
  fListQAitsPureSA->SetOwner();
  fListQAitsPureSA->SetName("ITS_PureSA");

  fListQAtpc=new TList;
  fListQAtpc->SetOwner();
  fListQAtpc->SetName("TPC");
  
  fListQAtrd=new TList;
  fListQAtrd->SetOwner();
  fListQAtrd->SetName("TRD");
  
  fListQAtof=new TList;
  fListQAtof->SetOwner();
  fListQAtof->SetName("TOF");
  
  fListQAemcal=new TList;
  fListQAemcal->SetOwner();
  fListQAemcal->SetName("EMCAL");

  fListQAtpctof=new TList;
  fListQAtpctof->SetOwner();
  fListQAtpctof->SetName("TPC_TOF");

  fListQA->Add(fListQAits);
  fListQA->Add(fListQAitsSA);
  fListQA->Add(fListQAitsPureSA);
  fListQA->Add(fListQAtpc);
  fListQA->Add(fListQAtrd);
  fListQA->Add(fListQAtof);
  fListQA->Add(fListQAemcal);
  fListQA->Add(fListQAtpctof);

  SetupITSqa();
  SetupTPCqa();
  SetupTRDqa();
  SetupTOFqa();
  SetupEMCALqa();
  SetupTPCTOFqa();

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
  FillTRDqa();
  FillTOFqa();
  FillEMCALqa();
  FillTPCTOFqa();

  PostData(1,fListQA);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillITSqa()
{
  //
  // Fill PID qa histograms for the ITS
  //

  AliVEvent *event=InputEvent();
  
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // ITS refit + ITS pid selection
    if (!( ( (status & AliVTrack::kITSrefit)==AliVTrack::kITSrefit ) ||
	   ! ( (status & AliVTrack::kITSpid  )==AliVTrack::kITSpid   ) )) continue;
    Double_t mom=track->P();
    
    TList *theList = 0x0;
    if(( (status & AliVTrack::kTPCin)==AliVTrack::kTPCin )){
      //ITS+TPC tracks
      theList=fListQAits;
    }else{
      if(!( (status & AliVTrack::kITSpureSA)==AliVTrack::kITSpureSA )){ 
	//ITS Standalone tracks
    	theList=fListQAitsSA;
      }else{
	//ITS Pure Standalone tracks
	theList=fListQAitsPureSA;
      }
    }
    
    
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      TH2 *h=(TH2*)theList->At(ispecie);
      if (!h) continue;
      Double_t nSigma=fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType)ispecie);
      h->Fill(mom,nSigma);
    }
    TH2 *h=(TH2*)theList->At(AliPID::kSPECIES);
    if (h) {
      Double_t sig=track->GetITSsignal();
      h->Fill(mom,sig);
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillTPCqa()
{
  //
  // Fill PID qa histograms for the TPC
  //
  
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
    if (!( (status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
        !( (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) ||
        !( (status & AliVTrack::kTPCpid  ) == AliVTrack::kTPCpid  ) ) continue;
    
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
void AliAnalysisTaskPIDqa::FillTRDqa()
{
  //
  // Fill PID qa histograms for the TRD
  //
  AliVEvent *event=InputEvent();
  Int_t ntracks = event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack <  ntracks; itrack++){
    AliVTrack *track = (AliVTrack *)event->GetTrack(itrack);

    //
    //basic track cuts
    //
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // TPC refit + ITS refit + TPC pid + TRD out
    if (!( (status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
        !( (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) ||
        !( (status & AliVTrack::kTPCpid  ) == AliVTrack::kTPCpid  ) ||
        !( (status & AliVTrack::kTRDout  ) == AliVTrack::kTRDout  )) continue;
    
    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }
    
    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;

    Double_t likelihoods[AliPID::kSPECIES];
    if(fPIDResponse->ComputeTRDProbability(track, AliPID::kSPECIES, likelihoods) != AliPIDResponse::kDetPidOk) continue;
    Int_t ntracklets = 0;
    Double_t momentum = -1.;
    for(Int_t itl = 0; itl < 6; itl++)
      if(track->GetTRDmomentum(itl) > 0.){
        ntracklets++;
        if(momentum < 0) momentum = track->GetTRDmomentum(itl);
    } 
    for(Int_t ispecie = 0; ispecie < AliPID::kSPECIES; ispecie++){
      TH2F *hLike = (TH2F *)fListQAtrd->At(ntracklets*AliPID::kSPECIES+ispecie);
      if (hLike) hLike->Fill(momentum,likelihoods[ispecie]);
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillTOFqa()
{
  AliVEvent *event=InputEvent();



  Int_t ntracks=event->GetNumberOfTracks();
  Int_t tracksAtTof = 0;
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);

    //
    //basic track cuts
    //
    ULong_t status=track->GetStatus();
    // TPC refit + ITS refit +
    // TOF out + TOFpid +
    // kTIME
    // (we don't use kTOFmismatch because it depends on TPC....)
    if (!((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
        !((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) ||
        !((status & AliVTrack::kTOFout  ) == AliVTrack::kTOFout  ) ||
        !((status & AliVTrack::kTOFpid  ) == AliVTrack::kTOFpid  ) ||
        !((status & AliVTrack::kTIME    ) == AliVTrack::kTIME    ) ) continue;

    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }

    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;

    tracksAtTof++;

    Double_t mom=track->P();

    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      TH2 *h=(TH2*)fListQAtof->At(ispecie);
      if (!h) continue;
      Double_t nSigma=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)ispecie);
      h->Fill(mom,nSigma);
    }

    TH2 *h=(TH2*)fListQAtof->FindObject("hSigP_TOF");
    if (h) {
      Double_t sig=track->GetTOFsignal();
      h->Fill(mom,sig);
    }

    Double_t mask = (Double_t)fPIDResponse->GetTOFResponse().GetStartTimeMask(mom) + 0.5;
    ((TH1F*)fListQAtof->FindObject("hStartTimeMask_TOF"))->Fill(mask);

    if (mom >= 1.0 && mom <= 2.0 ) {
      Double_t nsigma= fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kKaon);
      if (mask == 0) {
	((TH1F*)fListQAtof->FindObject("hnSigma_TOF_Kaon_T0-Fill"))->Fill(nsigma);
      } else if (mask == 1) {
	((TH1F*)fListQAtof->FindObject("hnSigma_TOF_Kaon_T0-TOF"))->Fill(nsigma);
      } else if ( (mask == 2) || (mask == 4) || (mask == 6) ) {
	((TH1F*)fListQAtof->FindObject("hnSigma_TOF_Kaon_T0-T0"))->Fill(nsigma);
      } else {
	((TH1F*)fListQAtof->FindObject("hnSigma_TOF_Kaon_T0-Best"))->Fill(nsigma);
      }
    }

    Double_t res = (Double_t)fPIDResponse->GetTOFResponse().GetStartTimeRes(mom);
    ((TH1F*)fListQAtof->FindObject("hStartTimeRes_TOF"))->Fill(res);

    AliESDEvent *esd = dynamic_cast<AliESDEvent *>(event);
    if (esd) {
      Double_t startTime = esd->GetT0TOF(0);
      if (startTime < 90000) ((TH1F*)fListQAtof->FindObject("hStartTimeAC_T0"))->Fill(startTime);
      else {
        startTime = esd->GetT0TOF(1);
        if (startTime < 90000) ((TH1F*)fListQAtof->FindObject("hStartTimeA_T0"))->Fill(startTime);
        startTime = esd->GetT0TOF(2);
        if (startTime < 90000) ((TH1F*)fListQAtof->FindObject("hStartTimeC_T0"))->Fill(startTime);
      }
    }
  }
  if (tracksAtTof > 0) {
    ((TH1F* )fListQAtof->FindObject("hnTracksAt_TOF"))->Fill(tracksAtTof);
    Int_t mask = fPIDResponse->GetTOFResponse().GetStartTimeMask(5.);
    if (mask & 0x1) ((TH1F*)fListQAtof->FindObject("hT0MakerEff"))->Fill(tracksAtTof);
  }

}


//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillEMCALqa()
{
  //
  // Fill PID qa histograms for the EMCAL
  //

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
    if (!( (status & AliVTrack::kEMCALmatch) == AliVTrack::kEMCALmatch) ) continue;

    Double_t pt=track->Pt();
   
    //EMCAL nSigma (only for electrons at the moment)
    TH2 *h=(TH2*)fListQAemcal->At(0);
    if (!h) continue;
    Double_t nSigma=fPIDResponse->NumberOfSigmasEMCAL(track, (AliPID::EParticleType)0);
    h->Fill(pt,nSigma);
    
    //EMCAL signal (E/p vs. pT)
    h=(TH2*)fListQAemcal->At(1);
    if (h) {

      Int_t nMatchClus = track->GetEMCALcluster();
      Double_t mom     = track->P();
      Double_t eop     = -1.;

      if(nMatchClus > -1){
    
	AliVCluster *matchedClus = (AliVCluster*)event->GetCaloCluster(nMatchClus);
    
	if(matchedClus){

	  // matched cluster is EMCAL
	  if(matchedClus->IsEMCAL()){
      
	    Double_t fClsE       = matchedClus->E();
	    eop                  = fClsE/mom;

	    h->Fill(pt,eop);
	    
	  }
	}
      }
      else{
	Printf("status status = AliVTrack::kEMCALmatch, BUT no matched cluster!");
      }
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::FillTPCTOFqa()
{
  //
  // Fill PID qa histograms for the TOF
  //   Here also the TPC histograms after TOF selection are filled
  //

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
    if (!((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
        !((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) ||
        !( (status & AliVTrack::kTPCpid  ) == AliVTrack::kTPCpid ) ||
        !((status & AliVTrack::kTOFout  ) == AliVTrack::kTOFout  ) ||
        !((status & AliVTrack::kTOFpid  ) == AliVTrack::kTOFpid  ) ||
        !((status & AliVTrack::kTIME    ) == AliVTrack::kTIME    ) ) continue;

    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }

    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;


    Double_t mom=track->P();
    Double_t momTPC=track->GetTPCmomentum();

    for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
      //TOF nSigma
      Double_t nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)ispecie);
      Double_t nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)ispecie);

      //TPC after TOF cut
      TH2 *h=(TH2*)fListQAtpctof->At(ispecie);
      if (h && TMath::Abs(nSigmaTOF)<3.) h->Fill(momTPC,nSigmaTPC);

      //TOF after TPC cut
      h=(TH2*)fListQAtpctof->At(ispecie+AliPID::kSPECIES);
      if (h && TMath::Abs(nSigmaTPC)<3.) h->Fill(mom,nSigmaTOF);

      //EMCAL after TOF and TPC cut
      h=(TH2*)fListQAtpctof->At(ispecie+2*AliPID::kSPECIES);
      if (h && TMath::Abs(nSigmaTOF)<3. && TMath::Abs(nSigmaTPC)<3. ){

	Int_t nMatchClus = track->GetEMCALcluster();
	Double_t pt      = track->Pt();
	Double_t eop     = -1.;
	
	if(nMatchClus > -1){
	  
	  AliVCluster *matchedClus = (AliVCluster*)event->GetCaloCluster(nMatchClus);
	  
	  if(matchedClus){
	    
	    // matched cluster is EMCAL
	    if(matchedClus->IsEMCAL()){
	      
	      Double_t fClsE       = matchedClus->E();
	      eop                  = fClsE/mom;

	      h->Fill(pt,eop);
 
	      
	    }
	  }
	}
      }
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
  
  //ITS+TPC tracks
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

  //ITS Standalone tracks
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaPSA = new TH2F(Form("hNsigmaP_ITSSA_%s",AliPID::ParticleName(ispecie)),
				Form("ITS n#sigma %s vs. p;p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
				vX->GetNrows()-1,vX->GetMatrixArray(),
				200,-10,10);
    fListQAitsSA->Add(hNsigmaPSA);
  }
  TH2F *hSigSA = new TH2F("hSigP_ITSSA",
			  "ITS signal vs. p;p [GeV]; ITS signal [arb. units]",
			  vX->GetNrows()-1,vX->GetMatrixArray(),
			  300,0,300);
  fListQAitsSA->Add(hSigSA);
  
  //ITS Pure Standalone tracks
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaPPureSA = new TH2F(Form("hNsigmaP_ITSPureSA_%s",AliPID::ParticleName(ispecie)),
				    Form("ITS n#sigma %s vs. p;p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
				    vX->GetNrows()-1,vX->GetMatrixArray(),
				    200,-10,10);
    fListQAitsPureSA->Add(hNsigmaPPureSA);
  }
  TH2F *hSigPureSA = new TH2F("hSigP_ITSPureSA",
			      "ITS signal vs. p;p [GeV]; ITS signal [arb. units]",
			      vX->GetNrows()-1,vX->GetMatrixArray(),
			      300,0,300);
  fListQAitsPureSA->Add(hSigPureSA);
  
  delete vX;  
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

  delete vX;  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupTRDqa()
{
  //
  // Create the TRD qa objects
  //
  TVectorD *vX=MakeLogBinning(200,.1,30);
  for(Int_t itl = 0; itl < 6; ++itl){
    for(Int_t ispecie = 0; ispecie < AliPID::kSPECIES; ispecie++){
      TH2F *hLikeP = new TH2F(Form("hLikeP_TRD_%dtls_%s", itl, AliPID::ParticleName(ispecie)),
                              Form("TRD Likelihood to be %s %s for tracks having %d %s; p (GeV/c); TRD %s Likelihood", ispecie == 0 ? "an" : "a", AliPID::ParticleName(ispecie), itl+1, itl == 0 ? "tracklet" : "tracklets", AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1, vX->GetMatrixArray(),
                              100, 0., 1.);
      fListQAtrd->Add(hLikeP);
    }
  }
  delete vX;
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

  // for Kaons PID we differentiate on Time Zero
  TH1F *hnSigT0Fill = new TH1F("hNsigma_TOF_Kaon_T0-Fill","TOF n#sigma (Kaon) T0-FILL [1-2. GeV/c]",200,-10,10);
  fListQAtof->Add(hnSigT0Fill);
  TH1F *hnSigT0T0 = new TH1F("hNsigma_TOF_Kaon_T0-T0","TOF n#sigma (Kaon) T0-T0 [1-2. GeV/c]",200,-10,10);
  fListQAtof->Add(hnSigT0T0);
  TH1F *hnSigT0TOF = new TH1F("hNsigma_TOF_Kaon_T0-TOF","TOF n#sigma (Kaon) T0-TOF [1.-2. GeV/c]",200,-10,10);
  fListQAtof->Add(hnSigT0TOF);
  TH1F *hnSigT0Best = new TH1F("hNsigma_TOF_Kaon_T0-Best","TOF n#sigma (Kaon) T0-Best [1-2. GeV/c]",200,-10,10);
  fListQAtof->Add(hnSigT0Best);


  TH2F *hSig = new TH2F("hSigP_TOF",
                        "TOF signal vs. p;p [GeV]; TOF signal [arb. units]",
                        vX->GetNrows()-1,vX->GetMatrixArray(),
                        300,0,300);

  delete vX;
  
  fListQAtof->Add(hSig);

  TH1F *hStartTimeMaskTOF = new TH1F("hStartTimeMask_TOF","StartTime mask",8,0,8);
  fListQAtof->Add(hStartTimeMaskTOF);
  TH1F *hStartTimeResTOF = new TH1F("hStartTimeRes_TOF","StartTime resolution [ps]",100,0,500);
  fListQAtof->Add(hStartTimeResTOF);

  TH1F *hnTracksAtTOF = new TH1F("hnTracksAt_TOF","Matched tracks at TOF",20,0,20);
  fListQAtof->Add(hnTracksAtTOF);
  TH1F *hT0MakerEff = new TH1F("hT0MakerEff","Events with T0-TOF vs nTracks",20,0,20);
  fListQAtof->Add(hT0MakerEff);

  // this in principle should stay on a T0 PID QA, but are just the data prepared for TOF use
  TH1F *hStartTimeAT0 = new TH1F("hStartTimeA_T0","StartTime from T0A [ps]",1000,-1000,1000);
  fListQAtof->Add(hStartTimeAT0);
  TH1F *hStartTimeCT0 = new TH1F("hStartTimeC_T0","StartTime from T0C [ps]",1000,-1000,1000);
  fListQAtof->Add(hStartTimeCT0);
  TH1F *hStartTimeACT0 = new TH1F("hStartTimeAC_T0","StartTime from T0AC [ps]",1000,-1000,1000);;
  fListQAtof->Add(hStartTimeACT0);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupEMCALqa()
{
  //
  // Create the EMCAL qa objects
  //

  TVectorD *vX=MakeLogBinning(200,.1,30);
  
  TH2F *hNsigmaPt = new TH2F(Form("hNsigmaPt_EMCAL_%s",AliPID::ParticleName(0)),
			     Form("EMCAL n#sigma %s vs. p_{T};p_{T} [GeV]; n#sigma",AliPID::ParticleName(0)),
			     vX->GetNrows()-1,vX->GetMatrixArray(),
			     200,-10,10);
  fListQAemcal->Add(hNsigmaPt);  
  
  TH2F *hSigPt = new TH2F("hSigPt_EMCAL",
                        "EMCAL signal (E/p) vs. p_{T};p_{T} [GeV]; EMCAL signal (E/p) [arb. units]",
                        vX->GetNrows()-1,vX->GetMatrixArray(),
                        200,0,2);
  fListQAemcal->Add(hSigPt);

  delete vX;  
}

//______________________________________________________________________________
void AliAnalysisTaskPIDqa::SetupTPCTOFqa()
{
  //
  // Create the qa objects for TPC + TOF combination
  //
  
  TVectorD *vX=MakeLogBinning(200,.1,30);

  //TPC signals after TOF cut
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaP = new TH2F(Form("hNsigmaP_TPC_TOF_%s",AliPID::ParticleName(ispecie)),
                              Form("TPC n#sigma %s vs. p (after TOF 3#sigma cut);p_{TPC} [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1,vX->GetMatrixArray(),
                              200,-10,10);
    fListQAtpctof->Add(hNsigmaP);
  }

  //TOF signals after TPC cut
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *hNsigmaP = new TH2F(Form("hNsigmaP_TOF_TPC_%s",AliPID::ParticleName(ispecie)),
                              Form("TOF n#sigma %s vs. p (after TPC n#sigma cut);p [GeV]; n#sigma",AliPID::ParticleName(ispecie)),
                              vX->GetNrows()-1,vX->GetMatrixArray(),
                              200,-10,10);
    fListQAtpctof->Add(hNsigmaP);
  }

  //EMCAL signal after TOF and TPC cut
  for (Int_t ispecie=0; ispecie<AliPID::kSPECIES; ++ispecie){
    TH2F *heopPt = new TH2F(Form("heopPt_TOF_TPC_%s",AliPID::ParticleName(ispecie)),
			    Form("EMCAL signal (E/p) %s vs. p_{T};p_{T} [GeV]; EMCAL signal (E/p) [arb. units]",AliPID::ParticleName(ispecie)),
			    vX->GetNrows()-1,vX->GetMatrixArray(),
			    200,0,2);
    fListQAtpctof->Add(heopPt);
  }

  delete vX;
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

