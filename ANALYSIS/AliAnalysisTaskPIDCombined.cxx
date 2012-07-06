/*************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                      for PID        Analysis                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>

#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>

#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>

#include "AliAnalysisTaskPIDCombined.h"

const char *AliAnalysisTaskPIDCombined::fgkBinMomDesc[AliAnalysisTaskPIDCombined::kPtBins] = {
  " 0 <= p < 0.5 GeV/c",
  " 0.5 <= p < 0.7 GeV/c",
  " 0.7 <= p < 1.0 GeV/c",
  " 1.0 <= p < 1.5 GeV/c",
  " 1.5 <= p < 2.0 GeV/c",
  " p >= 2.0 GeV/c"
};

ClassImp(AliAnalysisTaskPIDCombined)

//_________________________________________________________________________________
AliAnalysisTaskPIDCombined::AliAnalysisTaskPIDCombined() :
  AliAnalysisTaskSE(),
  fHistList(),
  fProbTPCnSigma(),
  fProbTOFnSigma(),
  fProbTPCTOFnSigmaTPC(),
  fProbTPC(),
  fProbTOF(),
  fProbTPCTOF(),
  fPriors(),
  fProbTPCTOFnSigTPCMom(),
  fProbTPCnSigTPCMom(),
  fProbTOFnSigTOFMom(),
  fPriorsUsed(),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fDeDx(NULL),
  fDeDxTuned(NULL)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskPIDCombined::AliAnalysisTaskPIDCombined(const char *name) :
  AliAnalysisTaskSE(name),
  fHistList(),
  fProbTPCnSigma(),
  fProbTOFnSigma(),
  fProbTPCTOFnSigmaTPC(),
  fProbTPC(),
  fProbTOF(),
  fProbTPCTOF(),
  fPriors(),
  fProbTPCTOFnSigTPCMom(),
  fProbTPCnSigTPCMom(),
  fProbTOFnSigTOFMom(),
  fPriorsUsed(),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fDeDx(NULL),
  fDeDxTuned(NULL)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
}

//_________________________________________________________________________________
void AliAnalysisTaskPIDCombined::UserCreateOutputObjects()
{
  //
  // Initialise the framework objects
  //


  // ------- track cuts
  fTrackCuts = new AliESDtrackCuts("fTrackCuts", "Standard");
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetMinNClustersTPC(80);
  fTrackCuts->SetMaxChi2PerClusterTPC(4);
  fTrackCuts->SetMaxDCAToVertexXY(3);
  fTrackCuts->SetMaxDCAToVertexZ(3);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  fTrackFilter->AddCuts(fTrackCuts);



  // ------- setup PIDCombined
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

  // no light nuclei - no need to call it, this is default
  //  fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);


  fHistList.Add(new TH1D("nEvents","Number of Evnets;Selection",2,0,2));

  for (Int_t ispec=0; ispec<5; ++ispec){

 
    fProbTPC[ispec]=new TH2D(Form("prob%s_mom_TPC",AliPID::ParticleName(ispec)),
                                   Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                   100,0.,20.,50,0.,1.);
    fHistList.Add(fProbTPC[ispec]);
    fProbTPCnSigma[ispec]=new TH2D(Form("prob%s_nSigma_TPC",AliPID::ParticleName(ispec)),
                                   Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                   20,-5.,5.,50,0.,1.);
    fHistList.Add(fProbTPCnSigma[ispec]);

    for (Int_t ibin=0;ibin<kPtBins;ibin++) {
      fProbTPCnSigTPCMom[ibin][ispec]=new TH2D(Form("prob%s_nSigma_TPC (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
					       Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
					       20,-5.,5.,50,0.,1.);
      fHistList.Add(fProbTPCnSigTPCMom[ibin][ispec]);
    }



    fProbTOF[ispec]=new TH2D(Form("prob%s_mom_TOF",AliPID::ParticleName(ispec)),
                                   Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                   100,0.,20.,50,0.,1.);
    fHistList.Add(fProbTOF[ispec]);
    fProbTOFnSigma[ispec]=new TH2D(Form("prob%s_nSigma_TOF",AliPID::ParticleName(ispec)),
                                   Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
                                   20,-5.,5.,50,0.,1.);
    fHistList.Add(fProbTOFnSigma[ispec]);
    for (Int_t ibin=0;ibin<kPtBins;ibin++) {
      fProbTOFnSigTOFMom[ibin][ispec]=new TH2D(Form("prob%s_nSigma_TOF (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
					       Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
					       20,-5.,5.,50,0.,1.);
      fHistList.Add(fProbTOFnSigTOFMom[ibin][ispec]);
    }



    fProbTPCTOF[ispec]=new TH2D(Form("prob%s_mom_TPCTOF",AliPID::ParticleName(ispec)),
                                   Form("%s probability vs. momentum;momentum;probability",AliPID::ParticleName(ispec)),
                                   100,0.,20.,50,0.,1.);
    fHistList.Add(fProbTPCTOF[ispec]);
    fProbTPCTOFnSigmaTPC[ispec]=new TH2D(Form("probTPCTOF%s_nSigma_TPC",AliPID::ParticleName(ispec)),
                                   Form("%s TPCTOF probability vs. n#sigmaTPC;n#sigma;probability",AliPID::ParticleName(ispec)),
                                   20,-5.,5.,50,0.,1.);
    fHistList.Add(fProbTPCTOFnSigmaTPC[ispec]);
    for (Int_t ibin=0;ibin<kPtBins;ibin++) {
      fProbTPCTOFnSigTPCMom[ibin][ispec]=new TH2D(Form("probTPCTOF%s_nSigma_TPC (%s)",AliPID::ParticleName(ispec),fgkBinMomDesc[ibin]),
					       Form("%s probability vs. n#sigma;n#sigma;probability",AliPID::ParticleName(ispec)),
					       20,-5.,5.,50,0.,1.);
      fHistList.Add(fProbTPCTOFnSigTPCMom[ibin][ispec]);
    }



    // basic priors
    fPriors[ispec]=new TH1F(Form("%s_priors",AliPID::ParticleName(ispec)),
			    Form("%s priors vs momentum",AliPID::ParticleName(ispec)),
			    100,0.,20.);
    fHistList.Add(fPriors[ispec]);
    switch (ispec) {
    case AliPID::kElectron:
      for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
      break;
    case AliPID::kMuon:
      for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.02);
      break;
    case AliPID::kPion:
      for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.56);
      break;
    case AliPID::kKaon:
      for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
      break;
    case AliPID::kProton:
      for (Int_t ich=1;ich<=100;ich++) fPriors[ispec]->SetBinContent(ich,0.20);
      break;
    default:
      break;
    }
    fPIDCombined->SetPriorDistribution((AliPID::EParticleType)ispec,fPriors[ispec]);

    // priors used
    fPriorsUsed[ispec] = new TH2D(Form("%s_priorsUsed",AliPID::ParticleName(ispec)),
			    Form("%s priors vs transverse momentum;p_{t} (GeV/c);priors",AliPID::ParticleName(ispec)),
				  100,0.,20.,101,0,1.01);      
    fHistList.Add(fPriorsUsed[ispec]);
  }


  fDeDx = new TH2D("hDeDx",";p_{TPC};dE/dx (a.u.)",500,0,5,500,0,500);
  fHistList.Add(fDeDx);
  fDeDxTuned = new TH2D("hDeDxTuned",";p_{TPC};dE/dx (a.u.)",500,0,5,500,0,500);
  fHistList.Add(fDeDxTuned);

  fHistList.SetOwner();
  PostData(1,&fHistList);


}

//_________________________________________________________________________________
void AliAnalysisTaskPIDCombined::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");

  //  Printf(" ---------------------- UserExec PID task ---------------------");
  
  FillHistogram("nEvents",0.);

  AliVEvent *event=InputEvent();
  AliVTrack *track=0x0;
  Int_t ntracks=event->GetNumberOfTracks();

  Double_t probTPC[AliPID::kSPECIES]={0.};
  Double_t probTOF[AliPID::kSPECIES]={0.};
  Double_t probTPCTOF[AliPID::kSPECIES]={0.};
  
  //loop over all tracks
  for (Int_t itrack=0; itrack<ntracks; ++itrack){

    track=(AliVTrack*)event->GetTrack(itrack);

    if ( fTrackFilter->IsSelected(track) ) {

      Double_t mom=track->GetTPCmomentum();
      Int_t ibin=GetMomBin(mom);

      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
     
      if (detUsed != 0) {  // TPC is available
	for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec) {
	  Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)ispec);
	  fProbTPC[ispec]->Fill(mom,probTPC[ispec]);
	  fProbTPCnSigma[ispec]->Fill(nSigmaTPC,probTPC[ispec]);
	  fProbTPCnSigTPCMom[ibin][ispec]->Fill(nSigmaTPC,probTPC[ispec]);
	}

	// compute priors for TPC+TOF, even if we ask just TOF for PID
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
	Double_t priors[5]; 	// check priors used for TOF
	fPIDCombined->GetPriors(track,priors,fPIDResponse,detUsed);
	for(Int_t ispec=0;ispec<5;ispec++) fPriorsUsed[ispec]->Fill(TMath::Abs(track->Pt()),priors[ispec]);

	if (detUsed != 0) {  // TOF is available
	  for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec) {
	    Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)ispec);
	    fProbTOF[ispec]->Fill(mom,probTOF[ispec]);
	    fProbTOFnSigma[ispec]->Fill(nSigmaTOF,probTOF[ispec]);
	    fProbTOFnSigTOFMom[ibin][ispec]->Fill(nSigmaTOF,probTOF[ispec]);
	  }
	}

	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
	if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask() ) {
	  for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec) {
	    Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)ispec);
	    fProbTPCTOF[ispec]->Fill(mom,probTPCTOF[ispec]);
	    fProbTPCTOFnSigmaTPC[ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);
	    fProbTPCTOFnSigTPCMom[ibin][ispec]->Fill(nSigmaTPC,probTPCTOF[ispec]);
	  }
	}

      }

      fPIDResponse->GetTPCsignalTunedOnData(track);

      fDeDx->Fill(mom,track->GetTPCsignal());
      fDeDxTuned->Fill(mom,track->GetTPCsignalTunedOnData());

    }
  }

  PostData(1, &fHistList);
}

//_________________________________________________________________________________
void AliAnalysisTaskPIDCombined::FillHistogram(const char* name, Double_t x, Double_t weight)
{
  //
  // Fill 1D histogram by name
  //
  ((TH1*)fHistList.FindObject(name))->Fill(x,weight);
}

//_________________________________________________________________________________
void AliAnalysisTaskPIDCombined::FillHistogram(const char* name, Double_t x, Double_t y, Double_t weight)
{
  //
  // Fill 2D histogram by name
  //
  ((TH2*)fHistList.FindObject(name))->Fill(x,y,weight);
}


//_________________________________________________________________________________
Int_t AliAnalysisTaskPIDCombined::GetMomBin(Float_t mom)
{
  //
  // Given momentum return histogram to be filled
  //
  if (mom>0. && mom < 0.5) return 0;
  if (mom>=0.5 && mom < 0.7) return 1;
  if (mom>=0.7 && mom < 1.0) return 2;
  if (mom>=1.0 && mom < 1.5) return 3;
  if (mom>=1.5 && mom < 2.0) return 4;
  return kPtBins-1;
}

