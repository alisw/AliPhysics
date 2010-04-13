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
//
// Class for PID QA
// Several studies done on clean samples of electrons, pions and kaons
// coming from V0 PID
// Compatible with both ESDs and AODs
//
// Autors:
//    Matus Kalisky <matus.kalisky@cern.ch>
//    Markus Heide <mheide@uni-muenster.de>
//    Markus Fasel <M.Fasel@gsi.de>
//
#include <TClass.h>
#include <TIterator.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliESDpid.h"
//#include "AliTRDPIDResponseLQ1D.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliExternalTrackParam.h"


#include "AliHFEcollection.h"
#include "AliHFEpidQA.h"
#include "AliHFEV0pid.h"
#include "AliHFEV0pidMC.h"
#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidQA)

//__________________________________________
AliHFEpidQA::AliHFEpidQA():
  fMC(NULL)
  , fV0pid(NULL)
  , fV0pidMC(NULL)
  , fOutput(NULL)
  , fT0(0)
  , fRun(0)
  , fESDpid(NULL)
{
  //
  // Default constructor
  //
  fESDpid = new AliESDpid;
}

//__________________________________________
AliHFEpidQA::~AliHFEpidQA(){
  //
  // Destructor
  //
  if(fV0pid) delete fV0pid;
  if(fV0pidMC) delete fV0pidMC;
  if(fOutput) delete fOutput;
  if(fESDpid) delete fESDpid;
//  if(fTRDpidResponse) delete fTRDpidResponse; 
}

//__________________________________________
void AliHFEpidQA::Init(){
  //
  // Prepare task output
  //

  fV0pid = new AliHFEV0pid;
  if(HasV0pidQA()) fV0pid->InitQA();
  fV0pidMC = new AliHFEV0pidMC();
  fV0pidMC->Init();

  fOutput = new AliHFEcollection("pidQA", "PID QA output");

  // 1st: Histos for purity studies 
  fOutput->CreateTH2F("purityElectron", "Electron Putrity", 2, -0.5, 1.5, 20, 0.1, 10);
  fOutput->BinLogAxis("purityElectron" ,1);
  fOutput->CreateTH2F("purityPionK0", "K0 Pion Putrity", 2, -0.5, 1.5, 20, 0.1, 10);
  fOutput->BinLogAxis("purityPionK0" ,1);
  fOutput->CreateTH2F("purityPionL", "Lambda Pion Putrity", 2, -0.5, 1.5, 20, 0.1, 10);
  fOutput->BinLogAxis("purityPionL" ,1);
  fOutput->CreateTH2F("purityProton", "Proton Putrity", 2, -0.5, 1.5, 20, 0.1, 10);
  fOutput->BinLogAxis("purityProton" ,1);

  // Histograms for TRD Electron Likelihood
  fOutput->CreateTH2F("hTRDelLikeElectron", "TRD Electron Likelihoods for Electrons; p (GeV/c); likelihood", 20, 0.1, 10, 100, 0., 1.);
  fOutput->BinLogAxis("hTRDelLikeElectron", 0);
  fOutput->CreateTH2F("hTRDelLikePionK0", "TRD Electron Likelihoods for K0 Pions; p (GeV/c); likelihood", 20, 0.1, 10, 100, 0., 1.);
  fOutput->BinLogAxis("hTRDelLikePionK0", 0);
  fOutput->CreateTH2F("hTRDelLikePionL", "TRD Electron Likelihoods for Lambda Pions; p (GeV/c); likelihood", 20, 0.1, 10, 100, 0., 1.);
  fOutput->BinLogAxis("hTRDelLikePionL", 0);
  fOutput->CreateTH2F("hTRDelLikeProton", "TRD Electron Likelihoods for Protons; p (GeV/c); likelihood", 20, 0.1, 10, 100, 0., 1.);
  fOutput->BinLogAxis("hTRDelLikeProton", 0);

  // TPC pid response
  fOutput->CreateTH2F("hTPC_dEdx_Electron", "TPC dEdx for conversion electrons; p (GeV/c); dEdx (a.u.)", 20, 0.1, 10, 200, 0, 200);
  fOutput->BinLogAxis("hTPC_dEdx_Electron", 0);
  fOutput->CreateTH2F("hTPC_dEdx_PionK0", "TPC dEdx for K0 pions; p (GeV/c); dEdx (a.u.)", 20, 0.1, 10, 200, 0, 200);
  fOutput->BinLogAxis("hTPC_dEdx_PionK0", 0);
  fOutput->CreateTH2F("hTPC_dEdx_PionL", "TPC dEdx for Lambda pions; p (GeV/c); dEdx (a.u.)", 20, 0.1, 10, 200, 0, 200);
  fOutput->BinLogAxis("hTPC_dEdx_PionL", 0);
  fOutput->CreateTH2F("hTPC_dEdx_Proton", "TPC dEdx for Lambda proton; p (GeV/c); dEdx (a.u.)", 20, 0.1, 10, 200, 0, 200);
  fOutput->BinLogAxis("hTPC_dEdx_Proton", 0);

 fOutput->CreateTH2F("hTPCnSigmaElectron", "TPC number of sigmas for conversion electrons; p (GeV/c); number of sigmas", 20, 0.1, 10, 100, -7, 7);
 fOutput->BinLogAxis("hTPCnSigmaElectron", 0);
 fOutput->CreateTH2F("hTPCnSigmaPionK0", "TPC number of sigmas for K0 pions; p (GeV/c); number of sigmas", 20, 0.1, 10, 100, -7, 7);
 fOutput->BinLogAxis("hTPCnSigmaPionK0", 0);
 fOutput->CreateTH2F("hTPCnSigmaPionL", "TPC number of sigmas for Lambda pions; p (GeV/c); number of sigmas", 20, 0.1, 10, 100, -7, 7);
 fOutput->BinLogAxis("hTPCnSigmaPionL", 0);
 fOutput->CreateTH2F("hTPCnSigmaProton", "TPC number of sigmas for Lambda protons; p (GeV/c); number of sigmas", 20, 0.1, 10, 100, -7, 7);
 fOutput->BinLogAxis("hTPCnSigmaProton", 0);

  fOutput->CreateTH2F("hTPC_PID", "TPC pid all tracks; tpc pid probability; species",100, 0, 1, 5, -0.5, 4.5 );
  fOutput->CreateTH2F("hTPC_PID_p_Electron", "TPC PID for conversion electrons; p (GeV/c); TPC PID", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTPC_PID_p_Electron", 0);
  fOutput->CreateTH2F("hTPC_PID_p_PionK0", "TPC PID for K0 pions; p (GeV/c); TPC PID", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTPC_PID_p_PionK0", 0);
  fOutput->CreateTH2F("hTPC_PID_p_PionL", "TPC PID for Lambda pions; p (GeV/c); TPC PID", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTPC_PID_p_PionL", 0);
  fOutput->CreateTH2F("hTPC_PID_p_Proton", "TPC PID for Lambda protons; p (GeV/c); TPC PID", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTPC_PID_p_Proton", 0);


  // TRD pid response
  fOutput->CreateTH2F("hTRD_Electron_trk", "all Electron candidate tracklets; p (GeV/c); N TRD tracklets", 100, 0.1, 10, 7, -0.5, 6.5);
  fOutput->BinLogAxis("hTRD_Electron_trk", 0);
  fOutput->CreateTH1F("hTRD_Electron_Y_like", "YES - V0 electron eff. for fixed likelihood cut; p (GeV/c); counts", 25, 0.1, 10);
  fOutput->BinLogAxis("hTRD_Electron_Y_like", 0);
  fOutput->CreateTH1F("hTRD_Electron_N_like", "NO - V0 electron eff. for fixed likelihood cut; p (GeV/c); counts", 25, 0.1, 10);
  fOutput->BinLogAxis("hTRD_Electron_N_like", 0);
  fOutput->CreateTH2F("hTRD_El_like_Electron", "V0 electron likelihoods for electrons; p (GeV/c); likelihood", 25, 0.1, 10, 1000, 0., 1.);
  fOutput->BinLogAxis("hTRD_El_like_Electron", 0);
  fOutput->CreateTH2F("hTRD_El_like_Pion", "V0 electron likelihoods for poins; p (GeV/c); likelihood", 25, 0.1, 10, 1000, 0., 1.);
  fOutput->BinLogAxis("hTRD_El_like_Pion", 0);
  fOutput->CreateTH2F("hTRD_El_like_Proton", "V0 electron likelihoods for protons; p (GeV/c); likelihood", 25, 0.1, 10, 1000, 0., 1.);
  fOutput->BinLogAxis("hTRD_El_like_Proton", 0);


  // TOF pid response
  
  fOutput->CreateTH2F("hTOF_PID", "TOF pid all tracks; tof pid probability; species",100, 0, 1,5,  -0.5, 4.5 );

  fOutput->CreateTH2F("hTOF_PID_p_Electron", "TOF PID for gamma converisons; p_T (GeV/c); counts", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTOF_PID_p_Electron", 0);
  fOutput->CreateTH2F("hTOF_PID_p_PionK0", "TOF PID for K0 pions; p_T (GeV/c); counts", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTOF_PID_p_PionK0", 0);
  fOutput->CreateTH2F("hTOF_PID_p_PionL", "TOF PID for Lambda pions; p_T (GeV/c); counts", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTOF_PID_p_PionL", 0);
  fOutput->CreateTH2F("hTOF_PID_p_Proton", "TOF PID for Lambda protons; p_T (GeV/c); counts", 100, 0.1, 10, 5, -0.5, 4.5);
  fOutput->BinLogAxis("hTOF_PID_p_Proton", 0);

  fOutput->CreateTH2F("hTOF_beta_Electron", "TOF beta for gamma conversions; #beta; p (GeV/c)", 120, 0, 1.2, 100, 0.1, 10);
  fOutput->BinLogAxis("hTOF_beta_Electron", 1);
  fOutput->CreateTH2F("hTOF_beta_PionK0", "TOF beta for K0 pions; #beta; p (GeV/c)", 120, 0, 1.2, 100, 0.1, 10);
  fOutput->BinLogAxis("hTOF_beta_PionK0", 1);
  fOutput->CreateTH2F("hTOF_beta_PionL", "TOF beta Lambda pions; #beta; p (GeV/c)", 120, 0, 1.2, 100, 0.1, 10);
  fOutput->BinLogAxis("hTOF_beta_PionL", 1);
  fOutput->CreateTH2F("hTOF_beta_Proton", "TOF beta for Lambda protons; #beta; p (GeV/c)", 120, 0, 1.2, 100, 0.1, 10);
  fOutput->BinLogAxis("hTOF_beta_Proton", 1);
  


  // Prepare TRD PID
/*  if(HasRecalculateTRDpid()){
    fTRDpidResponse = new AliTRDPIDResponseLQ1D;
    fTRDpidResponse->LoadReferences();
  }*/
}

//__________________________________________
void AliHFEpidQA::Process(AliVEvent *inputEvent){
  //
  // Run PID QA
  //

  if(fRun >= 104065 && fRun <= 104892){
    CorrectT0();
  }

  if(!fV0pid){
    AliError("V0pid not available! Forgotten to initialize?");
    return;
  }

  if(fMC) fV0pidMC->SetMCEvent(fMC);

  fV0pid->Process(inputEvent);
  TObjArray *electrons = fV0pid->GetListOfElectrons();
  TObjArray *pionsK0 = fV0pid->GetListOfPionsK0();
  TObjArray *pionsL = fV0pid->GetListOfPionsL();
  TObjArray *protons = fV0pid->GetListOfProtons();

  if(fMC){
    fV0pidMC->Process(electrons, AliHFEV0pid::kRecoElectron);
    fV0pidMC->Process(pionsK0, AliHFEV0pid::kRecoPionK0);
    fV0pidMC->Process(pionsL, AliHFEV0pid::kRecoPionL);
    fV0pidMC->Process(protons, AliHFEV0pid::kRecoProton);
  }

  AliDebug(2, Form("Number of Electrons      : %d", electrons->GetEntries()));
  AliDebug(2, Form("Number of K0 Pions       : %d", pionsK0->GetEntries()));
  AliDebug(2, Form("Number of Lambda Pions   : %d", pionsL->GetEntries()));
  AliDebug(2, Form("Number of Protons        : %d", protons->GetEntries()));
  if(fMC){
    AliDebug(2, "MC Information available. Doing Purity checks...");
    // Calculate the purity of the clean samples using MC 
    MakePurity(electrons, AliHFEV0pid::kRecoElectron);
    MakePurity(pionsK0,  AliHFEV0pid::kRecoPionK0);
    MakePurity(pionsL,  AliHFEV0pid::kRecoPionL);
    MakePurity(protons,  AliHFEV0pid::kRecoProton);
  }
  // Now we can do studies on the PID itself
  // TRD PID: Fill electron Likelihoods for the particle species
  FillTRDelectronLikelihoods(electrons,  AliHFEV0pid::kRecoElectron);
  FillTRDelectronLikelihoods(pionsK0,  AliHFEV0pid::kRecoPionK0);
  FillTRDelectronLikelihoods(pionsL,  AliHFEV0pid::kRecoPionL);
  FillTRDelectronLikelihoods(protons,  AliHFEV0pid::kRecoProton);
  
  FillPIDresponse(electrons, AliHFEV0pid::kRecoElectron);
  FillPIDresponse(pionsK0, AliHFEV0pid::kRecoPionK0);
  FillPIDresponse(pionsL, AliHFEV0pid::kRecoPionL);
  FillPIDresponse(protons, AliHFEV0pid::kRecoProton);

  // Analysis done, flush the containers
  fV0pid->Flush();
}

//__________________________________________
void AliHFEpidQA::MakePurity(TObjArray *tracks, Int_t species){
  //
  // Fill the QA histos for a given species
  //
  if(!fMC) return;
  AliDebug(3, Form("Doing Purity checks for species %d", species));
  Int_t pdg = 0;
  Char_t hname[256];
  switch(species){
    case  AliHFEV0pid::kRecoElectron:
      pdg = TMath::Abs(kElectron);
      sprintf(hname, "purityElectron");
      break;
    case  AliHFEV0pid::kRecoPionK0:
      pdg = TMath::Abs(kPiPlus);
      sprintf(hname, "purityPionK0");
      break;
    case  AliHFEV0pid::kRecoPionL:
      pdg = TMath::Abs(kPiPlus);
      sprintf(hname, "purityPionL");
      break;
    case  AliHFEV0pid::kRecoProton:
      pdg = TMath::Abs(kProton);
      sprintf(hname, "purityProton");
      break;
    default:  // non investigated species
      AliDebug(3, "Species not investigated");
      return;
  }
  AliDebug(3, Form("Number of tracks: %d", tracks->GetEntries()));
  TIterator *trackIter = tracks->MakeIterator();
  AliVParticle *recTrack = NULL, *mcTrack = NULL;
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter->Next()))){
    Int_t label = recTrack->GetLabel();
    AliDebug(4, Form("MC Label %d", label));
    mcTrack =fMC->GetTrack(TMath::Abs(label));
    if(!mcTrack){
      AliDebug(4, "MC track not available");
      continue; // we don't know
    }

    // Get the pdg code
    Int_t trackPdg = 0;
    if(!TString(mcTrack->IsA()->GetName()).CompareTo("AliMCParticle")){
      // case ESD
      AliMCParticle *mcp = dynamic_cast<AliMCParticle *>(mcTrack);
      trackPdg = TMath::Abs(mcp->Particle()->GetPdgCode());
    } else {
      // case AOD
      AliAODMCParticle *aodmcp = dynamic_cast<AliAODMCParticle *>(mcTrack);
      trackPdg = TMath::Abs(aodmcp->GetPdgCode());
    }
    if(trackPdg == pdg)    // Correct identification
      fOutput->Fill(hname, 0., recTrack->Pt());
    else  // Wrong identification
      fOutput->Fill(hname, 1., recTrack->Pt());
  }
  delete trackIter;
}

//__________________________________________
void AliHFEpidQA::FillTRDelectronLikelihoods(TObjArray * const particles, Int_t species){
  //
  // Fill electron Likelihoods for the TRD
  // Required for the calculation of the electron efficiency, 
  // pion and proton efficiency and the thresholds
  //
  Char_t hname[256] = "hTRDelLike";
  switch(species){
    case  AliHFEV0pid::kRecoElectron:
      sprintf(hname, "%sElectron", hname);
      break;
    case  AliHFEV0pid::kRecoPionK0:
      sprintf(hname, "%sPionK0", hname);
      break;
    case  AliHFEV0pid::kRecoPionL:
      sprintf(hname, "%sPionL", hname);
      break;
    case  AliHFEV0pid::kRecoProton:
      sprintf(hname, "%sProton", hname);
      break;
    default:
      AliDebug(2, Form("Species %d not investigated", species));
      return;
  };
  AliVParticle *recTrack = NULL;
  TIterator *trackIter = particles->MakeIterator();
  Double_t quantities[2] = {0., 0.};
  Double_t trdPidProbs[5];
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter->Next()))){
    if(!TString(recTrack->IsA()->GetName()).CompareTo("AliESDtrack")){
      // case ESD
      AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(recTrack);
      if(!esdTrack->GetTRDntracklets()) continue; // require at least 1 tracklet
      // take momentum at the innermost TRD layer
      Double_t p = 0.;
      for(Int_t ily = 0; ily < 6; ily++){
        if((p = esdTrack->GetTRDmomentum(ily)) > 1e-6) break;
      }
      quantities[0] = p;
      if(HasRecalculateTRDpid()) 
        RecalculateTRDpid(esdTrack, trdPidProbs);
      else
        esdTrack->GetTRDpid(trdPidProbs);
      quantities[1] = trdPidProbs[ AliPID::kElectron];
    }
    else{
      AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(recTrack);
      if(!aodTrack->GetDetPid()) continue;
      Float_t *trdMom = aodTrack->GetDetPid()->GetTRDmomentum(), p = 0.;
      for(Int_t ily = 0; ily < 6; ily++){
        if((p = trdMom[ily]) > 1e-6) break;
      }
      quantities[0] = p;
      // case AOD (for the moment lacks)
      if(HasRecalculateTRDpid()){
        RecalculateTRDpid(aodTrack, trdPidProbs); 
        quantities[1] = trdPidProbs[AliPID::kElectron];
      }
      else
        continue;
    }
    fOutput->Fill(hname, quantities[0], quantities[1]);
  }
}
//__________________________________________
void AliHFEpidQA::FillPIDresponse(TObjArray * const particles, Int_t species){
  //
  // Fill the PID response of different detectors to V0 daughter particles
  //
  Char_t hname[256] = "";
  const Char_t *typeName[5] = {"Electron", "PionK0", "PionL", "Kaon", "Proton"};
  const Int_t typePID[5] = {0, 2, 2, 3, 4};
  
  AliHFEpidTRD *pidTRD = new AliHFEpidTRD("TRDpid");
  
  AliVParticle *recTrack = NULL;
  TIterator *trackIter = particles->MakeIterator(); 
  while((recTrack = dynamic_cast<AliVParticle *>(trackIter->Next()))){
    // ESD
    if(!TString(recTrack->IsA()->GetName()).CompareTo("AliESDtrack")){
      // case ESD
      AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(recTrack);
      const AliExternalTrackParam *tpcIn = esdTrack->GetTPCInnerParam();
      if(!tpcIn) continue;
      
      // track kinematics
      Double_t p = tpcIn->P();
      //Double_t pt = tpcIn->Pt();

      // TPC dEdx
      Double_t dEdx = esdTrack->GetTPCsignal();
      sprintf(hname, "hTPC_dEdx_%s", typeName[species]);
      fOutput->Fill(hname, p, dEdx);

      //TPC number of sigmas
      Double_t nsigma = fESDpid->NumberOfSigmasTPC(esdTrack,(AliPID::EParticleType)typePID[species]);
      sprintf(hname, "hTPCnSigma%s",  typeName[species]);
      fOutput->Fill(hname, p, nsigma);

      // TPC PID response
      sprintf(hname, "hTPC_PID_p_%s", typeName[species]);
      Double_t tpcPID[5] = {-1, -1, -1, -1, -1};
      esdTrack->GetTPCpid(tpcPID);
      Int_t ix = 0;
      Double_t tmp = 0.;
      for(Int_t k=0; k<5; ++k){
	if(tpcPID[k] > tmp){
	  ix = k;
	  tmp = tpcPID[k];
	}
	fOutput->Fill("hTPC_PID", tpcPID[k], k);
      }
      if(tpcPID[ix] > 0){
	fOutput->Fill(hname, p, ix);
      }

      // TOF PID response
      sprintf(hname, "hTOF_PID_p_%s", typeName[species]);
      Double_t tofPID[5] = {-1., -1., -1., -1., -1};
      esdTrack->GetTOFpid(tofPID);
      tmp = 0.;
      for(Int_t k=0; k<5; ++k){
	if(tofPID[k] > tmp){
	  ix = k;
	  tmp = tofPID[k];
	}
	if(tofPID[k] > 0)
	  fOutput->Fill("hTOF_PID", tofPID[k], k);
      }
      if(tofPID[ix] > 0){
	fOutput->Fill(hname, p, ix);
      }
      
      //TRD first electron only
      Int_t nTRK = (int)esdTrack->GetTRDntrackletsPID();
      if(AliHFEV0pid::kRecoElectron == species){
	sprintf(hname, "hTRD_%s_trk", typeName[species]);
	fOutput->Fill(hname, p, nTRK);
      }
      Char_t n1[256] = "";
      Char_t n2[256] = "";	
      Double_t pidProbs[AliPID::kSPECIES];
      esdTrack->GetTRDpid(pidProbs);
      Double_t threshold = pidTRD->GetTRDthresholds(0.9, p);
      if(AliHFEV0pid::kRecoElectron == species && 6 == nTRK){
	sprintf(n1, "hTRD_%s_Y_like", typeName[species]);
	sprintf(n2, "hTRD_%s_N_like", typeName[species]);
	if(pidProbs[typePID[0]] > threshold) fOutput->Fill(n1, p);
	else fOutput->Fill(n2, p);
	sprintf(hname, "hTRD_El_like_Electron");
	fOutput->Fill(hname, p, pidProbs[typePID[species]]);
      }
      if( ((AliHFEV0pid::kRecoPionK0 == species) || (AliHFEV0pid::kRecoPionL == species)) && 6 == nTRK ){
	sprintf(hname, "hTRD_El_like_Pion");
	fOutput->Fill(hname, p, pidProbs[typePID[0]]);
      }
      if(AliHFEV0pid::kRecoProton == species && 6 == nTRK){
	sprintf(hname,"hTRD_El_like_Proton");
	fOutput->Fill(hname, p, pidProbs[typePID[0]]);
      }
    

      //TOF beta
      sprintf(hname, "hTOF_beta_%s", typeName[species]);
      Float_t beta = TOFbeta(esdTrack);
      fOutput->Fill(hname, beta, p);
    }
    // AOD - comming soon
    else{
      continue;
    }
  }// .. tracks in TObjArray
  
  if(pidTRD) delete pidTRD;

}

//__________________________________________
TList *AliHFEpidQA::GetOutput(){
  //
  // Getter for Output histograms
  //
  return fOutput->GetList();
}

//__________________________________________
TList *AliHFEpidQA::GetV0pidQA(){
  //
  // Getter for V0 PID QA histograms
  //
  return fV0pid->GetListOfQAhistograms();
}

//__________________________________________
TList *AliHFEpidQA::GetV0pidMC(){
  //
  // Getter for V0 PID QA histograms
  //
  if(fV0pidMC)
    return fV0pidMC->GetListOfQAhistograms();
  return NULL;
}

//__________________________________________
void AliHFEpidQA::RecalculateTRDpid(AliESDtrack * /*track*/, Double_t * /*pidProbs*/) const{
//  fTRDpidResponse->MakePID(track);
//  track->GetTRDpid(pidProbs);
}

//__________________________________________
void AliHFEpidQA::RecalculateTRDpid(AliAODTrack * /*track*/, Double_t * /*pidProbs*/) const{
//  fTRDpidResponse->MakePID(track, pidProbs);
}
//___________________________________________________________________
void AliHFEpidQA::CorrectT0(){
  // temporary solutions for correction the T0 for pass4 & pass5
  // returns corrected T0 for known runs
  // returns 0 if the correction failed
  if(! fRun > 0){
    AliError("Run number not set");
    fT0 = 0.;
    return;
  }
  Bool_t runFound = kFALSE;
  const Int_t corr[31][2] = {{104065, 1771614},
			     {104068, 1771603},
			     {104070, 1771594},
			     {104073, 1771610},
			     {104080, 1771305},
			     {104083, 1771613},
			     {104157, 1771665},
			     {104159, 1771679},
			     {104160, 1771633},
			     {104316, 1764344},
			     {104320, 1764342},
			     {104321, 1764371},
			     {104439, 1771750},
			     {104792, 1771755},
			     {104793, 1771762},
			     {104799, 1771828},
			     {104800, 1771788},
			     {104801, 1771796},
			     {104802, 1771775},
			     {104803, 1771795},
			     {104824, 1771751},
			     {104825, 1771763},
			     {104845, 1771792},
			     {104852, 1771817},
			     {104864, 1771825},
			     {104865, 1771827},
			     {104867, 1771841},
			     {104876, 1771856},
			     {104878, 1771847},
			     {104879, 1771830},
			     {104892, 1771837}};
 
  for(Int_t i=0; i<31; ++i){
    if(fRun == corr[i][0]){
      runFound = kTRUE;
      fT0 = (float)corr[i][1];
      // for the pass4 & pass5
      fT0 -= 37*1024*24.4 - 170.;
      //..
      break;
    }
  }

  if(!runFound){
    TString error = "Setting T0 correction FAILED, no TOF pid available for run: ";
    error += fRun;
    AliError(error);
    fT0 = 0.;
  }
  //cout<<" -D: run: "<<current_run<<" , fT0: "<<fT0<<endl;
  
}
//___________________________________________________________________
Float_t AliHFEpidQA::TOFbeta(AliESDtrack * const track) const {
  // computes the TOF beta
  Double_t l = track->GetIntegratedLength();  // cm
  Double_t t = track->GetTOFsignal();
  Double_t t0 = fT0; // ps
  if(l < 360. || l > 800.) return 0.;
  if(t <= 0.) return 0.;
  if(t0 <= 0.) return 0.;

  t -= t0; // subtract the T0

  l *= 0.01;  // cm ->m
  t *= 1e-12; //ps -> s

  
  Double_t v = l / t;
  Float_t beta = v / TMath::C();

  return beta;
}
 
