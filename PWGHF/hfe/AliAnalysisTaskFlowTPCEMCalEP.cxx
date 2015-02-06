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

// Class for heavy-flavour electron v2 with EMCal triggered events
// Author: Denise Godoy


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskFlowTPCEMCalEP.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliEventplane.h"
#include "AliCentrality.h"

#include "AliSelectNonHFE.h"


ClassImp(AliAnalysisTaskFlowTPCEMCalEP)
//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fVevent(0)
  ,fpidResponse(0)
  ,fMC(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fNonHFE(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fIsMC(kFALSE)
  ,fIsAOD(kFALSE)
  ,fSetMassConstraint(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(1000.)
  ,fInvmassCut(0.05)
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fnonHFEalgorithm("KF")
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,fD0_e(0)
  ,fTot_pi0e(0)
  ,fPhot_pi0e(0)
  ,fPhotBCG_pi0e(0)
  ,fTot_etae(0)
  ,fPhot_etae(0)
  ,fPhotBCG_etae(0)
  ,fInvMass(0)
  ,fInvMassBack(0)
  ,fDCA(0)
  ,fDCABack(0)
  ,fOpAngle(0)
  ,fOpAngleBack(0)
{
  //Named constructor

  for(Int_t k = 0; k < 3; k++) {
    fevPlaneV0[k] = NULL;
    feTPCV2[k] = NULL;
    feV2[k] = NULL;
    fChargPartV2[k] = NULL;
    fMtcPartV2[k] = NULL;

    fPi0Pt[k] = NULL;
    fEtaPt[k] = NULL;
    fInvmassLS[k] = NULL;
    fInvmassULS[k] = NULL;
    fOpeningAngleLS[k] = NULL;
    fOpeningAngleULS[k] = NULL;
  }

  for(Int_t i=0; i<3; i++) {
    for(Int_t j=0; j<8; j++) {
      for(Int_t k=0; k<4; k++) {
        fEoverPsig[i][j][k] = NULL;
	fEoverPuls[i][j][k] = NULL;
	fEoverPls[i][j][k] = NULL;
	fEoverPbcg[i][j][k] = NULL;
      }
    }
  }

  for(Int_t k = 0; k < 6; k++) {
    fDe[k]= NULL;
    fD0e[k]= NULL;
    fDpluse[k]= NULL;
    fDminuse[k]= NULL;
  }

  fPID = new AliHFEpid("hfePid");
  fTrackCuts = new AliESDtrackCuts();
  fNonHFE = new AliSelectNonHFE();

  InitParameters();
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
  ,fESD(0)
  ,fAOD(0)
  ,fVevent(0)
  ,fpidResponse(0)
  ,fMC(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fNonHFE(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fIsMC(kFALSE)
  ,fIsAOD(kFALSE)
  ,fSetMassConstraint(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)       
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(1000.)
  ,fInvmassCut(0.05)	
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fnonHFEalgorithm("KF")
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,fD0_e(0)
  ,fTot_pi0e(0)
  ,fPhot_pi0e(0)
  ,fPhotBCG_pi0e(0)
  ,fTot_etae(0)
  ,fPhot_etae(0)
  ,fPhotBCG_etae(0)
  ,fInvMass(0)
  ,fInvMassBack(0)
  ,fDCA(0)
  ,fDCABack(0)
  ,fOpAngle(0)
  ,fOpAngleBack(0)
{
  
  //Default constructor

  for(Int_t k = 0; k < 3; k++) {
    fevPlaneV0[k] = NULL;
    feTPCV2[k] = NULL;
    feV2[k] = NULL;
    fChargPartV2[k] = NULL;
    fMtcPartV2[k] = NULL;

    fPi0Pt[k] = NULL;
    fEtaPt[k] = NULL;
    fInvmassLS[k] = NULL;
    fInvmassULS[k] = NULL;
    fOpeningAngleLS[k] = NULL;
    fOpeningAngleULS[k] = NULL;
  }

  for(Int_t i=0; i<3; i++) {
    for(Int_t j=0; j<8; j++) {
      for(Int_t k=0; k<4; k++) {
        fEoverPsig[i][j][k] = NULL;
	fEoverPuls[i][j][k] = NULL;
	fEoverPls[i][j][k] = NULL;
	fEoverPbcg[i][j][k] = NULL;
      }
    }
  }

  for(Int_t k = 0; k < 6; k++) {
    fDe[k]= NULL;
    fD0e[k]= NULL;
    fDpluse[k]= NULL;
    fDminuse[k]= NULL;
  }

  fPID = new AliHFEpid("hfePid");
  fTrackCuts = new AliESDtrackCuts();
  fNonHFE = new AliSelectNonHFE();

  InitParameters();


  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskFlowTPCEMCalEP::~AliAnalysisTaskFlowTPCEMCalEP()
{
  //Destructor 
  
  delete fOutputList;
  delete fPID;
  delete fCFM;
  delete fPIDqa;
  delete fTrackCuts;
}
//_________________________________________

void AliAnalysisTaskFlowTPCEMCalEP::UserExec(Option_t*)
{
  //Main loop
  //Called for each event

  // create pointer to event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD){
    printf("ERROR: fESD not available\n");
    return;
  }

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fVevent){
    printf("ERROR: fVEvent not available\n");
    return;
  }

  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    if(fIsAOD)fPID->InitializePID(fAOD->GetRunNumber());
    if(!fIsAOD) fPID->InitializePID(fESD->GetRunNumber());
  }
 
  Bool_t SelColl = kTRUE;
  if(GetCollisionCandidates()==AliVEvent::kAny)
  {
     SelColl = kFALSE;
     TString firedTrigger;
     firedTrigger = fESD->GetFiredTriggerClasses();
     if(firedTrigger.Contains("CVLN_B2-B-NOPF-ALLNOTRD") || firedTrigger.Contains("CVLN_R1-B-NOPF-ALLNOTRD") || firedTrigger.Contains("CSEMI_R1-B-NOPF-ALLNOTRD"))SelColl=kTRUE;
     if(!SelColl)return;
  }

  if(fIsMC)fMC = MCEvent();
  AliStack* stack = NULL;
  if(fIsMC && fMC) stack = fMC->Stack();
 
  Int_t fNOtrks =fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();

  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();

  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(0);

  if(fNOtrks<2) return;

  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }

  fPID->SetPIDResponse(fpidResponse);

  fCFM->SetRecEventInfo(fVevent);

  Float_t cent = -1.;
  Int_t iPt=8, iCent=3, iDeltaphi=4;
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");
  fCent->Fill(cent);
 
  if (cent>=0  && cent<10) iCent=0;
  if (cent>=10 && cent<20) iCent=1;
  if (cent>=20 && cent<40) iCent=2;
  if (cent<0 || cent>=40) return;
 
  //Event planes

  Double_t evPlaneV0A = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0A",fESD,2));
  if(evPlaneV0A > TMath::Pi()) evPlaneV0A = evPlaneV0A - TMath::Pi();

  Double_t evPlaneV0C = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0C",fESD,2));
  if(evPlaneV0C > TMath::Pi()) evPlaneV0C = evPlaneV0C - TMath::Pi();

  Double_t evPlaneV0 = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0",fESD,2));
  if(evPlaneV0 > TMath::Pi()) evPlaneV0 = evPlaneV0 - TMath::Pi();
  fevPlaneV0[iCent]->Fill(evPlaneV0);

  AliEventplane* esdTPCep = fESD->GetEventplane();
  TVector2 *standardQ = 0x0;
  Double_t qx = -999., qy = -999.;
  standardQ = esdTPCep->GetQVector(); 
  if(!standardQ)return;
 
  qx = standardQ->X();
  qy = standardQ->Y();

  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  Float_t evPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;

  //Event plane resolutions

  // --> 2 subevent method (only for TPC EP)

  TVector2 *qsub1a = esdTPCep->GetQsub1();
  TVector2 *qsub2a = esdTPCep->GetQsub2();
  Double_t evPlaneResTPC = -999.;
  if(qsub1a && qsub2a){
    evPlaneResTPC = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }

  fTPCsubEPres->Fill(evPlaneResTPC);

  // --> 3 event method (V0, V0A, and V0C EP)

  Double_t Qx2pos = 0., Qy2pos = 0., Qx2neg = 0., Qy2neg = 0., Qweight = 1;

  for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {

    AliVParticle* vparticle = fVevent->GetTrack(iTracks);
    if (!vparticle){
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    AliESDtrack *trackEP = dynamic_cast<AliESDtrack*>(vparticle);

    if (!trackEP) {
      printf("ERROR: Could not receive track %d\n", iTracks);
     continue;
    }

    if(TMath::Abs(trackEP->Eta())>0.8 || trackEP->Pt() < 0.15 || trackEP->Pt() > 4) continue;

    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, trackEP)) continue;

    if (trackEP->Pt() < 2) Qweight = trackEP->Pt()/2;
    if (trackEP->Pt() >= 2) Qweight = 1;


    if(trackEP->Eta()>0 && trackEP->Eta()<0.8){
      Qx2pos += Qweight*TMath::Cos(2*trackEP->Phi());
      Qy2pos += Qweight*TMath::Sin(2*trackEP->Phi());
    }
    if(trackEP->Eta()<0 && trackEP->Eta()>-0.8){
      Qx2neg += Qweight*TMath::Cos(2*trackEP->Phi());
      Qy2neg += Qweight*TMath::Sin(2*trackEP->Phi());
    }
  }//track loop only for EP 

  Double_t evPlaneTPCneg = TMath::ATan2(Qy2neg, Qx2neg)/2;
  Double_t evPlaneTPCpos = TMath::ATan2(Qy2pos, Qx2pos)/2;

  Double_t evPlaneRes[4]={GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos),
                          GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg),
                          GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg),cent};
  fEPres->Fill(evPlaneRes);

  // Selection of primary pi0 and eta in MC to compute the weight
  if(fIsMC && fMC && stack){
    Int_t nParticles = stack->GetNtrack();
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      int fPDG = particle->GetPdgCode(); 
      double pTMC = particle->Pt();

      Double_t etaMC = particle->Eta();
      if (TMath::Abs(etaMC)>1.2)continue;

      Bool_t isMotherPrimary = IsPi0EtaPrimary(particle,stack);
      Bool_t isFromLMdecay = IsPi0EtaFromLMdecay(particle,stack);
      Bool_t isFromHFdecay = IsPi0EtaFromHFdecay(particle,stack);

      if (isMotherPrimary && !isFromHFdecay && !isFromLMdecay){
        if(fPDG==111) fPi0Pt[iCent]->Fill(pTMC); //pi0
        if(fPDG==221) fEtaPt[iCent]->Fill(pTMC); //eta
     }
    }
  }//MC

  Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};
  Double_t deltaPhiRange[4];
  for(Int_t j=0;j<4;j++){
    deltaPhiRange[j] = j*(TMath::Pi()/4);
  }

  // Track loop 
  for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {

    AliVParticle* vparticle = fVevent->GetTrack(iTracks);
    if (!vparticle){
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    AliESDtrack *track = dynamic_cast<AliESDtrack*>(vparticle);
     
    if(!track) continue;
    
    if (TMath::Abs(track->Eta())>0.7) continue;
 
    fTrackPtBefTrkCuts->Fill(track->Pt());

    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

    if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
      if(track->GetKinkIndex(0) != 0) continue;
    } 

    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

    fTrackPtAftTrkCuts->Fill(track->Pt());

    Double_t clsE=-9.,p=-99.,EovP=-99.,pt=-99.,dEdx=-99.,fTPCnSigma=9.,phi=-9.,m02=-9.,m20=-9.,fEMCalnSigma=9.,dphi=9.,cosdphi=9.;
 
    pt = track->Pt();
    if(pt<1.5) continue;
    fTrkpt->Fill(pt);
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        iPt=i;
        continue;
      }
    }

    Int_t clsId = track->GetEMCALcluster();
    if (clsId>0){
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(clsId);
      if(cluster && cluster->IsEMCAL()){
        clsE = cluster->E();
        m20 = cluster->GetM20();
        m02 = cluster->GetM02();
      }
    }

    p = track->P();
    phi = track->Phi(); 
    dEdx = track->GetTPCsignal();
    EovP = clsE/p;
    fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
    fEMCalnSigma = GetSigmaEMCal(EovP, pt, cent);
    fdEdxBef->Fill(p,dEdx);
    fTPCnsigma->Fill(p,fTPCnSigma);


    dphi = GetDeltaPhi(phi,evPlaneV0);   
    cosdphi = GetCos2DeltaPhi(phi,evPlaneV0);   
    for(Int_t i=0;i<3;i++) {
      if (dphi>=deltaPhiRange[i] && dphi<deltaPhiRange[i+1]){
        iDeltaphi=i;
        continue;
      }
    }

    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t fFlagPhotonicElecBCG = kFALSE;
    Double_t weight = 1.; 

    SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,weight,iCent);
   
    Int_t partPDG = -99;
    Double_t partPt = -99.;
    Bool_t MChijing; 

    if(fIsMC && fMC && stack){
      if(fTPCnSigma < -1 || fTPCnSigma > 3) continue;
      Int_t label = track->GetLabel();
      if(label!=0){
        TParticle *particle = stack->Particle(TMath::Abs(label));	
        if(particle){
          partPDG = particle->GetPdgCode();
          partPt = particle->Pt();
          if (TMath::Abs(partPDG)!=11) continue;

          MChijing = fMC->IsFromBGEvent(TMath::Abs(label));
          int iHijing = 1;
          if(!MChijing) iHijing = 0; // 0 if enhanced sample

          Bool_t pi0Decay = IsElectronFromPi0(particle,stack,weight,cent);
          Bool_t etaDecay = IsElectronFromEta(particle,stack,weight,cent);

          SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,weight,iCent);

          if (pi0Decay){
            fTot_pi0e->Fill(partPt,weight);
            if(fFlagPhotonicElec) fPhot_pi0e->Fill(partPt,weight);
            if(fFlagPhotonicElecBCG) fPhotBCG_pi0e->Fill(partPt,weight);
          }
          if (etaDecay){
            fTot_etae->Fill(partPt,weight);
            if(fFlagPhotonicElec) fPhot_etae->Fill(partPt,weight);
            if(fFlagPhotonicElecBCG) fPhotBCG_etae->Fill(partPt,weight);
          }
        }// end particle
      }// end label
    }//end MC

    if(fTPCnSigma >= 1.5 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,EovP);
    Int_t pidpassed = 1;
    
    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;

    if (m20>0.02 && m02>0.02){
      Double_t corr[7]={iCent,iPt,fTPCnSigma,fEMCalnSigma,m02,dphi,cosdphi};
      fCorr->Fill(corr);
    }
  
    fChargPartV2[iCent]->Fill(iPt,cosdphi); 
    if (clsE>0) fMtcPartV2[iCent]->Fill(iPt,cosdphi); 

    if (pidpassed==0) continue;

    if (fTPCnSigma>=-5 && fTPCnSigma<-3.2) fEoverPbcg[iCent][iPt][iDeltaphi]->Fill(EovP);
    if (fTPCnSigma>=-0.5 && fTPCnSigma<3) feTPCV2[iCent]->Fill(iPt,cosdphi); 
    if (fTPCnSigma>=-0.5 && fTPCnSigma<3 && fEMCalnSigma>-1 && fEMCalnSigma<3) feV2[iCent]->Fill(iPt,cosdphi); 

    fTrkEovPAft->Fill(pt,EovP);
    fdEdxAft->Fill(p,dEdx);

    if(fFlagPhotonicElec){
      fPhotoElecPt->Fill(pt);
    }

    if (!fFlagPhotonicElec) fSemiInclElecPt->Fill(pt);

    if (m20>0.02 && m02>0.02 && m02<0.27 && fTPCnSigma>-1 && fTPCnSigma<3){ 
      fEoverPsig[iCent][iPt][iDeltaphi]->Fill(EovP);
      if (fFlagPhotonicElec) fEoverPuls[iCent][iPt][iDeltaphi]->Fill(EovP);
      if (fFlagPhotonicElecBCG) fEoverPls[iCent][iPt][iDeltaphi]->Fill(EovP);
    }

  }//end of track loop 
  PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::UserCreateOutputObjects()
{
  //--- Check MC
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
    fIsMC = kTRUE;
    printf("+++++ MC Data available");
  }
  //--------Initialize PID
  fPID->SetHasMCData(fIsMC);
  
  if(!fPID->GetNumberOfPIDdetectors()) 
    {
      fPID->AddDetector("TPC", 0);
      fPID->AddDetector("EMCAL", 1);
    }
  
  fPID->SortDetectors(); 
  fPIDqa = new AliHFEpidQAmanager();
  fPIDqa->Initialize(fPID);
  
  //--------Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  fCuts->Initialize(fCFM);
  
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fPIDqa->MakeList("PIDQA"));
  
  fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
  fOutputList->Add(fNoEvents);
  
  fTrkpt = new TH1F("fTrkpt","track pt",100,0,50);
  fOutputList->Add(fTrkpt);
  
  fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",100,0,50);
  fOutputList->Add(fTrackPtBefTrkCuts);
  
  fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",100,0,50);
  fOutputList->Add(fTrackPtAftTrkCuts);
  
  fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",100,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigma);
  
  fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPBef);
  
  fTrkEovPAft = new TH2F("fTrkEovPAft","track E/p after HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPAft);
  
  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxBef);
  
  fdEdxAft = new TH2F("fdEdxAft","track dEdx vs p after HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxAft);
  
  fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",100,0,50);
  fOutputList->Add(fPhotoElecPt);
  
  fSemiInclElecPt = new TH1F("fSemiInclElecPt", "Semi-inclusive electron pt",100,0,50);
  fOutputList->Add(fSemiInclElecPt);
  
  fCent = new TH1F("fCent","Centrality",100,0,100) ;
  fOutputList->Add(fCent);
  
  fTPCsubEPres = new TH1F("fTPCsubEPres","TPC subevent plane resolution",100,-1,1);
  fOutputList->Add(fTPCsubEPres);
  
  Int_t binsv1[4]={100,100,100,90}; // V0-TPCpos, V0-TPCneg, TPCpos-TPCneg, cent
  Double_t xminv1[4]={-1,-1,-1,0};
  Double_t xmaxv1[4]={1,1,1,90};
  fEPres = new THnSparseD ("fEPres","EP resolution",4,binsv1,xminv1,xmaxv1);
  fOutputList->Add(fEPres);
	
  //iCent,iPt,fTPCnSigma,fEMCalnSigma,m02,dphi,cosdphi
  Int_t binsv2[7]={3,8,100,100,100,120,100}; 
  Double_t xminv2[7]={0,0,-5,-5,0,0,-1};
  Double_t xmaxv2[7]={3,8,5,5,2,TMath::TwoPi(),1}; 
  fCorr = new THnSparseD ("fCorr","Correlations",7,binsv2,xminv2,xmaxv2);
  fOutputList->Add(fCorr);
  
  for(Int_t i=0; i<3; i++) {
    fevPlaneV0[i] = new TH1F(Form("fevPlaneV0%d",i),"V0 EP",100,0,TMath::Pi());
    fOutputList->Add(fevPlaneV0[i]);

    feTPCV2[i] = new TH2F(Form("feTPCV2%d",i), "", 8,0,8,100,-1,1);
    fOutputList->Add(feTPCV2[i]);

    feV2[i] = new TH2F(Form("feV2%d",i), "", 8,0,8,100,-1,1);
    fOutputList->Add(feV2[i]);

    fChargPartV2[i] = new TH2F(Form("fChargPartV2%d",i), "", 8,0,8,100,-1,1);
    fOutputList->Add(fChargPartV2[i]);

    fMtcPartV2[i] = new TH2F(Form("fMtcPartV2%d",i), "", 8,0,8,100,-1,1);
    fOutputList->Add(fMtcPartV2[i]);

    fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
    fOutputList->Add(fInvmassLS[i]);

    fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
    fOutputList->Add(fInvmassULS[i]);

    fOpeningAngleLS[i] = new TH2F(Form("fOpeningAngleLS%d",i),"Opening angle for LS pairs",100,0,1,100,0,50);
    fOutputList->Add(fOpeningAngleLS[i]);

    fOpeningAngleULS[i] = new TH2F(Form("fOpeningAngleULS%d",i),"Opening angle for ULS pairs",100,0,1,100,0,50);
    fOutputList->Add(fOpeningAngleULS[i]);

    fPi0Pt[i] = new TH1F(Form("fPi0Pt%d",i), "Pi0 weight",100,0,50);
    fOutputList->Add(fPi0Pt[i]);

    fEtaPt[i] = new TH1F(Form("fEtaPt%d",i), "Eta weight",100,0,50);
    fOutputList->Add(fEtaPt[i]);
  }

  for(Int_t i=0; i<3; i++) {
    for(Int_t j=0; j<8; j++) {
      for(Int_t k=0; k<4; k++) {
        fEoverPsig[i][j][k] = new TH1F(Form("fEoverPsig%d%d%d",i,j,k), "",100,0,3);
        fOutputList->Add(fEoverPsig[i][j][k]);

	fEoverPuls[i][j][k] = new TH1F(Form("fEoverPuls%d%d%d",i,j,k), "",100,0,3);
        fOutputList->Add(fEoverPuls[i][j][k]);

	fEoverPls[i][j][k] = new TH1F(Form("fEoverPls%d%d%d",i,j,k), "",100,0,3);
        fOutputList->Add(fEoverPls[i][j][k]);

	fEoverPbcg[i][j][k] = new TH1F(Form("fEoverPbcg%d%d%d",i,j,k), "",100,0,3);
        fOutputList->Add(fEoverPbcg[i][j][k]);
      }
    }
  }

  fD0_e = new TH2F("fD0_e", "D0 vs e",100,0,50,200,-6.3,6.3);
  fOutputList->Add(fD0_e);
  
  for(Int_t k = 0; k < 6; k++) {
    
    TString De_name = Form("fDe%d",k);
    TString D0e_name = Form("fD0e%d",k);
    TString Dpluse_name = Form("fDpluse%d",k);
    TString Dminuse_name = Form("fDminuse%d",k);
       
    fDe[k] = new TH1F((const char*)De_name,"",100,0,50);
    fD0e[k] = new TH1F((const char*)D0e_name,"",100,0,50);
    fDpluse[k] = new TH1F((const char*)Dpluse_name,"",100,0,50);
    fDminuse[k] = new TH1F((const char*)Dminuse_name,"",100,0,50);
    
    fOutputList->Add(fDe[k]);
    fOutputList->Add(fD0e[k]);
    fOutputList->Add(fDpluse[k]);
    fOutputList->Add(fDminuse[k]);
    
  }

  int nbin_v2 = 8;
  double bin_v2[9] = {2,2.5,3,4,6,8,10,13,18};
  
  fTot_pi0e = new TH1F("fTot_pi0e","fTot_pi0e",nbin_v2,bin_v2);
  fOutputList->Add(fTot_pi0e);
  
  fPhot_pi0e = new TH1F("fPhot_pi0e","fPhot_pi0e",nbin_v2,bin_v2);
  fOutputList->Add(fPhot_pi0e);
  
  fPhotBCG_pi0e = new TH1F("fPhotBCG_pi0e","fPhotBCG_pi0e",nbin_v2,bin_v2);
  fOutputList->Add(fPhotBCG_pi0e);
    
  fTot_etae = new TH1F("fTot_etae","fTot_etae",nbin_v2,bin_v2);
  fOutputList->Add(fTot_etae);
  
  fPhot_etae = new TH1F("fPhot_etae","fPhot_etae",nbin_v2,bin_v2);
  fOutputList->Add(fPhot_etae);
  
  fPhotBCG_etae = new TH1F("fPhotBCG_etae","fPhotBCG_etae",nbin_v2,bin_v2);  
  fOutputList->Add(fPhotBCG_etae);

  fInvMass = new TH1F("fInvMass","",200,0,0.3);
  fOutputList->Add(fInvMass);

  fInvMassBack = new TH1F("fInvMassBack","",200,0,0.3);
  fOutputList->Add(fInvMassBack);

  fDCA = new TH1F("fDCA","",200,0,1);
  fOutputList->Add(fDCA);

  fDCABack = new TH1F("fDCABack","",200,0,1);
  fOutputList->Add(fDCABack);

  fOpAngle = new TH1F("fOpAngle","",200,0,0.5);
  fOutputList->Add(fOpAngle);

  fOpAngleBack = new TH1F("fOpAngleBack","",200,0,0.5);
  fOutputList->Add(fOpAngleBack);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::Terminate(Option_t *)
{
  // Info("Terminate");
	AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetCos2DeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get cos[2(phi-psi_EP)] or cos[2(psi_subEP1 - psi_subEP2)]
  Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB); 
  if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
  Double_t cos2DeltaPhi = TMath::Cos(2*dPhi);
  
  return cos2DeltaPhi;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetDeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get phi-psi_EP
  Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB); 
  if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
  
  return dPhi;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetPi0weight(Double_t mcPi0pT, Float_t cent) const
{
  //Get Pi0 weight
  double weight = 1.0;
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEtaweight(Double_t mcEtapT, Float_t cent) const
{
  //Get eta weight
  double weight = 1.0;
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetSigmaEMCal(Double_t EoverP, Double_t pt, Float_t cent) const
{
  //Get sigma for EMCal PID
  Double_t NumberOfSigmasEMCal = 99.;
  Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};

  if (cent>=0  && cent<10){
    Double_t mean[8]={0.953184,0.957259,0.97798,0.9875,1.03409,1.06257,1.02776,1.04338};
    Double_t sigma[8]={0.130003,0.113493,0.092966,0.0836828,0.101804,0.0893414,0.0950752,0.050427};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (mean[i]-EoverP)/sigma[i];  
        continue;
      }
    }
  }
  if (cent>=10 && cent<20){
    Double_t mean[8]={0.96905,0.952985,0.96871,0.983934,1.00047,0.988736,1.02101,1.04557};
    Double_t sigma[8]={0.0978103,0.103215,0.0958494,0.0797962,0.0719482,0.0672677,0.0754882,0.0461192};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (mean[i]-EoverP)/sigma[i];  
        continue;
      }
    }
  }
  if (cent>=20 && cent<40){
    Double_t mean[8]={0.947362,0.951933,0.959288,0.977004,0.984502,1.02004,1.00489,0.986696};
    Double_t sigma[8]={0.100127,0.0887731,0.0842077,0.0787335,0.0804325,0.0652376,0.0766669,0.0597849};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (mean[i]-EoverP)/sigma[i];  
        continue;
      }
    }
  }



  return NumberOfSigmasEMCal;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsPi0EtaFromHFdecay(TParticle *particle, AliStack* stack) 
{
  // Check if the mother comes from heavy-flavour decays

  Bool_t isHFdecay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  if (TMath::Abs(partPDG)!=111 || TMath::Abs(partPDG)!=221) return isHFdecay; // particle is not pi0 or eta

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = stack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();

    // c decays
    if((TMath::Abs(motherPDG)==411) || (TMath::Abs(motherPDG)==421) || (TMath::Abs(motherPDG)==431) || (TMath::Abs(motherPDG)==4122) || (TMath::Abs(motherPDG)==4132) || (TMath::Abs(motherPDG)==4232) || (TMath::Abs(motherPDG)==43320)) isHFdecay = kTRUE; 

    // b decays
    if((TMath::Abs(motherPDG)==511) || (TMath::Abs(motherPDG)==521) || (TMath::Abs(motherPDG)==531) || (TMath::Abs(motherPDG)==5122) || (TMath::Abs(motherPDG)==5132) || (TMath::Abs(motherPDG)==5232) || (TMath::Abs(motherPDG)==53320)) isHFdecay = kTRUE; 
  }

  return isHFdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsPi0EtaFromLMdecay(TParticle *particle, AliStack* stack) 
{
  // Check if the mother comes from light-meson decays

  Bool_t isLMdecay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  if (TMath::Abs(partPDG)!=111 || TMath::Abs(partPDG)!=221) return isLMdecay; // particle is not pi0 or eta

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = stack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();

    if(motherPDG == 111 || motherPDG == 221 || motherPDG==223 || motherPDG==333 || motherPDG==331 || (TMath::Abs(motherPDG)==113) || (TMath::Abs(motherPDG)==213) || (TMath::Abs(motherPDG)==313) || (TMath::Abs(motherPDG)==323)) isLMdecay = kTRUE;
  }

  return isLMdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsPi0EtaPrimary(TParticle *particle, AliStack* stack) 
{
  // Check if the pi0 or eta are primary

  Bool_t isprimary = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  Bool_t pi0etaprimary = particle->IsPrimary();
  if (pi0etaprimary) isprimary = kTRUE;  
  
  return isprimary;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsElectronFromPi0(TParticle *particle, AliStack* stack, Double_t &weight, Float_t cent) 
{
  // Check if electron comes from primary pi0 not from light-meson and heavy-flavour decays

  Bool_t isPi0Decay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  if (TMath::Abs(partPDG)!=11) return isPi0Decay; // particle is not electron

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = stack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();
    Double_t motherPt = mother->Pt();

    Bool_t isMotherPi0primary = IsPi0EtaPrimary(mother,stack);
    Bool_t isMotherPi0fromHF = IsPi0EtaFromHFdecay(mother,stack);
    Bool_t isMotherPi0fromLM = IsPi0EtaFromLMdecay(mother,stack);

    if (motherPDG==111 && (isMotherPi0primary || (!isMotherPi0fromHF && !isMotherPi0fromLM))){ // pi0 -> e 
      isPi0Decay = kTRUE; 
      weight = GetPi0weight(motherPt,cent);
    }

    Int_t idSecondMother = particle->GetSecondMother(); 
    if (motherPDG==22 && idSecondMother>0){
      TParticle *secondMother = stack->Particle(idSecondMother);
      Int_t secondMotherPDG = secondMother->GetPdgCode();
      Double_t secondMotherPt = secondMother->Pt();
    
      Bool_t isSecondMotherPi0primary = IsPi0EtaPrimary(secondMother,stack);
      Bool_t isSecondMotherPi0fromHF = IsPi0EtaFromHFdecay(secondMother,stack);
      Bool_t isSecondMotherPi0fromLM = IsPi0EtaFromLMdecay(secondMother,stack);

      if (secondMotherPDG==111 && (isSecondMotherPi0primary || (!isSecondMotherPi0fromHF && !isSecondMotherPi0fromLM))){ //pi0 -> gamma -> e 
        isPi0Decay = kTRUE;
        weight = GetPi0weight(secondMotherPt,cent);
      }
    }
  }
  return isPi0Decay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsElectronFromEta(TParticle *particle, AliStack* stack, Double_t &weight, Float_t cent)
{
  // Check if electron comes from primary eta not from light-meson and heavy-flavour decays

  Bool_t isEtaDecay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  if (TMath::Abs(partPDG)!=11) return isEtaDecay; // particle is not electron

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = stack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();
    Double_t motherPt = mother->Pt();

    Bool_t isMotherEtaprimary = IsPi0EtaPrimary(mother,stack);
    Bool_t isMotherEtafromHF = IsPi0EtaFromHFdecay(mother,stack);
    Bool_t isMotherEtafromLM = IsPi0EtaFromLMdecay(mother,stack);
    
    if (motherPDG==221  && (isMotherEtaprimary || (!isMotherEtafromHF && !isMotherEtafromLM))){ //primary eta -> e
      isEtaDecay = kTRUE; 
      weight = GetEtaweight(motherPt,cent);
    }

    Int_t idSecondMother = mother->GetFirstMother();	
    if ((motherPDG==22 || motherPDG==111) && idSecondMother>0){
      TParticle *secondMother = stack->Particle(idSecondMother);
      Int_t secondMotherPDG = secondMother->GetPdgCode();
      Double_t secondMotherPt = secondMother->Pt();

      Bool_t isSecondMotherEtaprimary = IsPi0EtaPrimary(secondMother,stack);
      Bool_t isSecondMotherEtafromHF = IsPi0EtaFromHFdecay(secondMother,stack);
      Bool_t isSecondMotherEtafromLM = IsPi0EtaFromLMdecay(secondMother,stack);

      if (secondMotherPDG==221  && (isSecondMotherEtaprimary || (!isSecondMotherEtafromHF && !isSecondMotherEtafromLM))){ //eta -> pi0/g-> e
        isEtaDecay = kTRUE; 
        weight = GetEtaweight(secondMotherPt,cent);
      }
      Int_t idThirdMother = secondMother->GetFirstMother();
      if (idThirdMother>0){
        TParticle *thirdMother = stack->Particle(idThirdMother);
        Int_t thirdMotherPDG = thirdMother->GetPdgCode();
        Double_t thirdMotherPt = thirdMother->Pt();

        Bool_t isThirdMotherEtaprimary = IsPi0EtaPrimary(thirdMother,stack);
        Bool_t isThirdMotherEtafromHF = IsPi0EtaFromHFdecay(thirdMother,stack);
        Bool_t isThirdMotherEtafromLM = IsPi0EtaFromLMdecay(thirdMother,stack);

        if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221 && (isThirdMotherEtaprimary || (!isThirdMotherEtafromHF && !isThirdMotherEtafromLM))){//p eta->pi0->g-> e 
          isEtaDecay = kTRUE; 
          weight = GetEtaweight(thirdMotherPt,cent);
        }
      }
    }
  }
  return isEtaDecay;
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::SelectPhotonicElectron(Int_t iTracks,AliESDtrack *track,Bool_t &fFlagPhotonicElec, Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent)
{
  //Identify non-heavy flavour electrons using Invariant mass method
  
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.9,0.9);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(4);
  fTrackCuts->SetMinNClustersTPC(80);
  fTrackCuts->SetMaxDCAToVertexZ(3.2);
  fTrackCuts->SetMaxDCAToVertexXY(2.4);
  fTrackCuts->SetDCAToVertex2D(kTRUE);
  
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Bool_t flagPhotonicElec = kFALSE;
  Bool_t flagPhotonicElecBCG = kFALSE;
  
  for(Int_t jTracks = 0; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    
    if(jTracks==iTracks) continue;
    
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }
    
    Double_t pt=-999., ptAsso=-999., nTPCsigmaAsso=-999.;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Double_t openingAngle = -999., mass=999., width = -999;
    Int_t chargeAsso = 0, charge = 0;
    
    nTPCsigmaAsso = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
    pt = track->Pt();
    ptAsso = trackAsso->Pt();
    chargeAsso = trackAsso->Charge();
    charge = track->Charge();
    
    if(ptAsso <0.5) continue;
    if(!fTrackCuts->AcceptTrack(trackAsso)) continue;
    if(TMath::Abs(nTPCsigmaAsso)>3) continue;
    
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
    
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    AliKFParticle ge1(*track, fPDGe1);
    AliKFParticle ge2(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);
    
    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
    
    if(fSetMassConstraint && pVtx) {
      AliKFVertex primV(*pVtx);
      primV += recg;
      primV -= ge1;
      primV -= ge2;
      recg.SetProductionVertex(primV);
      recg.SetMassConstraint(0,0.0001);
    }

    openingAngle = ge1.GetAngle(ge2);

    if(fFlagLS) fOpeningAngleLS[iCent]->Fill(openingAngle,pt);
    if(fFlagULS) fOpeningAngleULS[iCent]->Fill(openingAngle,pt);

    if(openingAngle > fOpeningAngleCut) continue;
    
    recg.GetMass(mass,width);
    
    if(fFlagLS) fInvmassLS[iCent]->Fill(mass,pt,weight);
    if(fFlagULS) fInvmassULS[iCent]->Fill(mass,pt,weight);

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec) flagPhotonicElec = kTRUE;
    if(mass<fInvmassCut && fFlagLS && !flagPhotonicElecBCG) flagPhotonicElecBCG = kTRUE;
    
  }
  fFlagPhotonicElec = flagPhotonicElec;
  fFlagPhotonicElecBCG = flagPhotonicElecBCG;
  
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::InitParameters()
{
  // Init parameters

  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.7,0.7);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts->SetMinNClustersTPC(100);
  fTrackCuts->SetPtRange(0.5,100);

  fNonHFE->SetAODanalysis(kFALSE);
  fNonHFE->SetInvariantMassCut(fInvmassCut);
  fNonHFE->SetOpeningAngleCut(fOpeningAngleCut);
  fNonHFE->SetChi2OverNDFCut(fChi2Cut);
  fNonHFE->SetAlgorithm(fnonHFEalgorithm); //KF or DCA
  if (fnonHFEalgorithm == "DCA") fNonHFE->SetDCACut(fDCAcut);
  fNonHFE->SetTrackCuts(-3.5,3.5,fTrackCuts);

}

