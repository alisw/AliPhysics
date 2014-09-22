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
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.05)
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fminCent(0)
  ,fmaxCent(0)
  ,fnonHFEalgorithm("KF")
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fMCphotoElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fevPlaneV0A(0)
  ,fevPlaneV0C(0)
  ,fevPlaneV0(0)
  ,fevPlaneTPC(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,feTPCV2(0)
  ,feV2(0)
  ,fphoteV2(0)
  ,fChargPartV2(0)
  ,fGammaWeight(0)
  ,fPi0Weight(0)
  ,fEtaWeight(0)
  ,fD0Weight(0)
  ,fDplusWeight(0)
  ,fDminusWeight(0)
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
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)       
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.05)	
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fminCent(0)
  ,fmaxCent(0)
  ,fnonHFEalgorithm("KF")
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fMCphotoElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fevPlaneV0A(0)
  ,fevPlaneV0C(0)
  ,fevPlaneV0(0)
  ,fevPlaneTPC(0)
  ,fTPCsubEPres(0)
  ,fEPres(0)
  ,fCorr(0)
  ,feTPCV2(0)
  ,feV2(0)
  ,fphoteV2(0)
  ,fChargPartV2(0)
  ,fGammaWeight(0)
  ,fPi0Weight(0)
  ,fEtaWeight(0)
  ,fD0Weight(0)
  ,fDplusWeight(0)
  ,fDminusWeight(0)
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
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");
  fCent->Fill(cent);
 
//  cout<<"TEST: "<<fInvmassCut<< "   "<< fOpeningAngleCut<<"   "<<fnonHFEalgorithm<<"   "<<fminCent<<"   "<<fmaxCent<<" here!!!!!!!!!!!" <<endl;

  //Event planes

  Double_t evPlaneV0A = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0A",fESD,2));
  if(evPlaneV0A > TMath::Pi()) evPlaneV0A = evPlaneV0A - TMath::Pi();
  fevPlaneV0A->Fill(evPlaneV0A);

  Double_t evPlaneV0C = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0C",fESD,2));
  if(evPlaneV0C > TMath::Pi()) evPlaneV0C = evPlaneV0C - TMath::Pi();
  fevPlaneV0C->Fill(evPlaneV0C);

  Double_t evPlaneV0 = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0",fESD,2));
  if(evPlaneV0 > TMath::Pi()) evPlaneV0 = evPlaneV0 - TMath::Pi();
  fevPlaneV0->Fill(evPlaneV0);

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
  fevPlaneTPC->Fill(evPlaneTPC);

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

  Double_t Qx2pos = -999., Qy2pos = -999., Qx2neg = -999., Qy2neg = -999., Qweight = 1;

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

  Double_t evPlaneRes[7]={GetCos2DeltaPhi(evPlaneV0A,evPlaneV0C),GetCos2DeltaPhi(evPlaneV0A,evPlaneTPC),
      GetCos2DeltaPhi(evPlaneV0C,evPlaneTPC),GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos),
      GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg),GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg),cent};
  fEPres->Fill(evPlaneRes);

  // MC
  if(fIsMC && fMC && stack){
    Int_t nParticles = stack->GetNtrack();
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      int fPDG = particle->GetPdgCode(); 
      double pTMC = particle->Pt();
      double etaMC = particle->Eta();
      if(fabs(etaMC)>0.7)continue;

      Bool_t MChijing = fMC->IsFromBGEvent(iParticle);
      int iHijing = 1;
      if(!MChijing)iHijing = 0;

      if(fPDG==111)fPi0Weight->Fill(pTMC,iHijing);//pi0
      if(fPDG==221)fEtaWeight->Fill(pTMC,iHijing);//eta
      if(fPDG==421)fD0Weight->Fill(pTMC,iHijing);//D0
      if(fPDG==411)fDplusWeight->Fill(pTMC,iHijing);//D+
      if(fPDG==-411)fDminusWeight->Fill(pTMC,iHijing);//D-
	
      Int_t idMother = particle->GetFirstMother();
      if (idMother>0){
        TParticle *mother = stack->Particle(idMother);
        int motherPDG = mother->GetPdgCode();
        if(fPDG==22 && motherPDG!=111 && motherPDG!=221)fGammaWeight->Fill(pTMC,iHijing);//gamma
      } 
    }
  }

  Double_t ptRange[8] = {2, 2.5, 3, 4, 6, 8, 10, 13};
  Double_t ptDmeson[7] = {2,3,4,6,8,12,16};
  Double_t deltaPhiRange[7];
  for(Int_t j=0;j<7;j++){
    deltaPhiRange[j] = j*(TMath::Pi()/6);
  }

  // Track loop 
  for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {

    AliVParticle* vparticle = fVevent->GetTrack(iTracks);
    if (!vparticle){
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vparticle);
    AliESDtrack *track = dynamic_cast<AliESDtrack*>(vparticle);

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

    Double_t clsE = -999., p = -999., EovP=-999., pt = -999., dEdx=-999., fTPCnSigma=0, phi=-999., 
      wclsE = -999., wEovP = -999., m02= -999., m20= -999.;
 
    pt = track->Pt();
    if(pt<1.5) continue;
    fTrkpt->Fill(pt);

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
    wEovP = wclsE/p;
    fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
    fdEdxBef->Fill(p,dEdx);
    fTPCnsigma->Fill(p,fTPCnSigma);
    
    //Remove electron candidate from the event plane
    Float_t evPlaneCorrTPC = evPlaneTPC;
    if(dEdx>70 && dEdx<90){
      Double_t qX = standardQ->X() - esdTPCep->GetQContributionX(track); 
      Double_t qY = standardQ->Y() - esdTPCep->GetQContributionY(track); 
      TVector2 newQVectorfortrack;
      newQVectorfortrack.Set(qX,qY);
      evPlaneCorrTPC = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 
    }

    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t fFlagPhotonicElecBCG = kFALSE;

    //Non-HFE reconstruction
    fNonHFE->SetPIDresponse(fpidResponse);
    fNonHFE->FindNonHFE(iTracks,vparticle,fVevent);
    Int_t *fUlsPartner = fNonHFE->GetPartnersULS(); // Pointer to the ULS partners index
    Int_t *fLsPartner = fNonHFE->GetPartnersLS(); // Pointer to the LS partners index

    if (fNonHFE->IsULS()) fFlagPhotonicElec=kTRUE;
    if (fNonHFE->IsLS()) fFlagPhotonicElecBCG=kTRUE;
 
    fNonHFE->SetHistAngleBack(fOpAngleBack);
    fNonHFE->SetHistAngle(fOpAngle);
    fNonHFE->SetHistMassBack(fInvMassBack);
    fNonHFE->SetHistMass(fInvMass);
    if (fnonHFEalgorithm == "DCA")fNonHFE->SetHistDCABack(fDCABack);
    if (fnonHFEalgorithm == "DCA")fNonHFE->SetHistDCA(fDCA);


    Double_t corr[11]={phi,fTPCnSigma,cent,pt,EovP,GetDeltaPhi(phi,evPlaneV0),GetCos2DeltaPhi(phi,evPlaneV0), static_cast<Double_t>(fFlagPhotonicElec), 
    static_cast<Double_t>(fFlagPhotonicElecBCG),m02,m20};
    fCorr->Fill(corr);

    Int_t whichFirstMother = 0, whichSecondMother = 0, whichThirdMother = 0; 
    Int_t whichPart = -99;
    Int_t partPDG = -99, motherPDG = -99, secondMotherPDG = -99, thirdMotherPDG = -99;
    Double_t partPt = -99. , motherPt = -99., secondMotherPt = -99.,thirdMotherPt = -99.;
    Double_t weight = 1.; 
    Double_t Dweight = 1.; 
    Bool_t MChijing; 

    Bool_t pi0Decay= kFALSE;
    Bool_t etaDecay= kFALSE;


    Double_t phiD = -999.,phie = -999.,phieRec = -999.,ptD = -999.,pte = -999.,pteRec = -999.;
 
    if(fIsMC && fMC && stack){
      Int_t label = track->GetLabel();
      if(label>0){
        TParticle *particle = stack->Particle(label);
        if(particle){
          partPDG = particle->GetPdgCode();
          partPt = particle->Pt();

          if (TMath::Abs(partPDG)==11) whichPart = 0; //electron
          if (partPDG==22) whichPart = 3; //gamma
          if (partPDG==111) whichPart = 2; //pi0
          if (partPDG==221) whichPart = 1; //eta

          MChijing = fMC->IsFromBGEvent(label);

          int iHijing = 1;
          if(!MChijing) iHijing = 0; // 0 if enhanced sample

          Int_t idMother = particle->GetFirstMother();
          if (idMother>0){
            TParticle *mother = stack->Particle(idMother);
            motherPt = mother->Pt();
            motherPDG = mother->GetPdgCode();

            if (motherPDG==22) whichFirstMother = 3; //gamma
            if (motherPDG==111) whichFirstMother = 2; //pi0
            if (motherPDG==221) whichFirstMother = 1; //eta

            Int_t idSecondMother = particle->GetSecondMother();
            if (idSecondMother>0){
              TParticle *secondMother = stack->Particle(idSecondMother);
              secondMotherPt = secondMother->Pt();
              secondMotherPDG = secondMother->GetPdgCode();
  
              if (secondMotherPDG==111) whichSecondMother = 2; //pi0
              if (secondMotherPDG==221) whichSecondMother = 1; //eta

              Int_t idThirdMother = secondMother->GetFirstMother();
              if (idThirdMother>0){
                TParticle *thirdMother = stack->Particle(idThirdMother);
                thirdMotherPt = thirdMother->Pt();
                thirdMotherPDG = thirdMother->GetPdgCode();

                if (thirdMotherPDG==221) whichThirdMother = 1; //eta
              }//third mother
            }//second mother
  

            if (TMath::Abs(partPDG)==11){

              // D meson decay
              if (motherPDG==421 || TMath::Abs(motherPDG)==411){ // D
                Double_t phi_D = -99., Deltaphi_De=-99.;
                phi_D = mother->Phi();
                Deltaphi_De = phi_D - phi;
  
                fD0_e->Fill(pt,Deltaphi_De);
  
                Dweight= GetDweight(0,motherPt,cent);
                if(iHijing==1) Dweight = 1.0;
                for(Int_t i=0;i<6;i++){
                  if (motherPt>=ptDmeson[i] && motherPt<ptDmeson[i+1]) fDe[i]->Fill(pt,Dweight);
                }
              }
              if (motherPDG==421){ // D0
                Dweight= GetDweight(1,motherPt,cent);
                if(iHijing==1) Dweight = 1.0;
                for(Int_t i=0;i<6;i++){
                  if (motherPt>=ptDmeson[i] && motherPt<ptDmeson[i+1]) fD0e[i]->Fill(pt,Dweight);
                }
              }
              if (motherPDG==411){ // D+
                Dweight= GetDweight(2,motherPt,cent);
                if(iHijing==1) Dweight = 1.0;
                for(Int_t i=0;i<6;i++){
                  if (motherPt>=ptDmeson[i] && motherPt<ptDmeson[i+1]) fDpluse[i]->Fill(pt,Dweight);
                }
              }
              if (motherPDG==-411){ //D-
                Dweight= GetDweight(3,motherPt,cent);
                if(iHijing==1) Dweight = 1.0;
                for(Int_t i=0;i<6;i++){
                  if (motherPt>=ptDmeson[i] && motherPt<ptDmeson[i+1]) fDminuse[i]->Fill(pt,Dweight);
                }
              }
  
              //pi0 decay 
              if (motherPDG==111 && secondMotherPDG!=221){ //not eta -> pi0 -> e
                weight = GetPi0weight(motherPt,cent);
                pi0Decay = kTRUE;
              }
              if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG!=221){ //not eta -> pi0 -> gamma -> e
                weight = GetPi0weight(secondMotherPt,cent);
                pi0Decay = kTRUE;
              }
  
              //eta decay
              if (motherPDG==221){ //eta -> e
                weight = GetEtaweight(motherPt,cent);
                etaDecay = kTRUE;
              }
              if (motherPDG==111 && secondMotherPDG==221){ //eta -> pi0 -> e
                weight = GetEtaweight(secondMotherPt,cent);
                etaDecay = kTRUE;
              }
              if (motherPDG==22 && secondMotherPDG==221){ //eta -> gamma -> e
                weight = GetEtaweight(secondMotherPt,cent);
                etaDecay = kTRUE;
              }
              if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221){ //eta -> pi0 -> gamma -> e
                weight = GetEtaweight(thirdMotherPt,cent);
                etaDecay = kTRUE;
              }
            }// end of electron
  
            if (fTPCnSigma>-1 && fTPCnSigma<3 && EovP>1 && EovP<1.3 && (motherPDG==22 || motherPDG==111 || motherPDG==221)){
              if(iHijing==1) weight = 1.0;
  
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
            }
  
            Double_t mc[15]={EovP,fTPCnSigma,partPt, static_cast<Double_t>(fFlagPhotonicElec),
                static_cast<Double_t>(fFlagPhotonicElecBCG), static_cast<Double_t>(whichPart),cent,pt,
                  static_cast<Double_t>(whichFirstMother), static_cast<Double_t>(whichSecondMother),
            static_cast<Double_t>(whichThirdMother), static_cast<Double_t>(iHijing),motherPt,secondMotherPt,thirdMotherPt};
  
            //fMCphotoElecPt->Fill(mc);
            if (motherPDG==22 || motherPDG==111 || motherPDG==221) fMCphotoElecPt->Fill(mc);// mother = gamma, pi0, eta
          }//end mother          
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
  
    Double_t corrV2[7]={phi,cent,pt,EovP,GetCos2DeltaPhi(phi,evPlaneTPC),GetCos2DeltaPhi(phi,evPlaneV0A),
    GetCos2DeltaPhi(phi,evPlaneV0C)};
    fChargPartV2->Fill(corrV2); 

    if(fTPCnSigma >= -0.5){
      Double_t correctedV2[5]={cent,pt,GetCos2DeltaPhi(phi,evPlaneCorrTPC),GetCos2DeltaPhi(phi,evPlaneV0A),
      GetCos2DeltaPhi(phi,evPlaneV0C)};
      feTPCV2->Fill(correctedV2);
    }

    if(pidpassed==0) continue;

    Double_t correctedV2[5]={cent,pt,GetCos2DeltaPhi(phi,evPlaneCorrTPC),GetCos2DeltaPhi(phi,evPlaneV0A),
     GetCos2DeltaPhi(phi,evPlaneV0C)};

    feV2->Fill(correctedV2);
    fTrkEovPAft->Fill(pt,EovP);
    fdEdxAft->Fill(p,dEdx);

    if(fFlagPhotonicElec){
      fphoteV2->Fill(correctedV2);
      fPhotoElecPt->Fill(pt);
    }

    if(!fFlagPhotonicElec) fSemiInclElecPt->Fill(pt);

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
  
  fevPlaneV0A = new TH1F("fevPlaneV0A","V0A EP",100,0,TMath::Pi());
  fOutputList->Add(fevPlaneV0A);
  
  fevPlaneV0C = new TH1F("fevPlaneV0C","V0C EP",100,0,TMath::Pi());
  fOutputList->Add(fevPlaneV0C);
  
  fevPlaneV0 = new TH1F("fevPlaneV0","V0 EP",100,0,TMath::Pi());
  fOutputList->Add(fevPlaneV0);

  fevPlaneTPC = new TH1F("fevPlaneTPC","TPC EP",100,0,TMath::Pi());
  fOutputList->Add(fevPlaneTPC);
    
  fTPCsubEPres = new TH1F("fTPCsubEPres","TPC subevent plane resolution",100,-1,1);
  fOutputList->Add(fTPCsubEPres);
  
  Int_t binsv1[7]={100,100,100,100,100,100,90}; // V0A-V0C, V0A-TPC, V0C-TPC, V0-TPCpos, V0-TPCneg, TPCpos-TPCneg, cent
  Double_t xminv1[7]={-1,-1,-1,-1,-1,-1,0};
  Double_t xmaxv1[7]={1,1,1,1,1,1,90};
  fEPres = new THnSparseD ("fEPres","EP resolution",7,binsv1,xminv1,xmaxv1);
  fOutputList->Add(fEPres);
	
  //phi,fTPCnSigma,cent,pt,EovP,GetDeltaPhi(phi,evPlaneV0),GetCos2DeltaPhi(phi,evPlaneV0), fFlagPhotonicElec, fFlagPhotonicElecBCG,m02,m20
  Int_t binsv2[11]={100,200,90,100,100,100,100,3,3,100,100}; 
  Double_t xminv2[11]={0,-5,0,0,0,0,-1,-1,-1,0,0};
  Double_t xmaxv2[11]={2*TMath::Pi(),5,90,50,3,TMath::Pi(),1,2,2,1,1}; 
  fCorr = new THnSparseD ("fCorr","Correlations",11,binsv2,xminv2,xmaxv2);
  fOutputList->Add(fCorr);
    
  Int_t binsv3[5]={90,100,100,100,100}; // cent, pt, TPCcos2DeltaPhi, V0Acos2DeltaPhi, V0Ccos2DeltaPhi
  Double_t xminv3[5]={0,0,-1,-1,-1};
  Double_t xmaxv3[5]={90,50,1,1,1}; 
  feV2 = new THnSparseD ("feV2","inclusive electron v2",5,binsv3,xminv3,xmaxv3);
  fOutputList->Add(feV2);
  
  Int_t binsv4[5]={90,100,100,100,100}; // cent, pt, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv4[5]={0,0,-1,-1,-1};
  Double_t xmaxv4[5]={90,50,1,1,1}; 
  fphoteV2 = new THnSparseD ("fphoteV2","photonic electron v2",5,binsv4,xminv4,xmaxv4);
  fOutputList->Add(fphoteV2);
  
  Int_t binsv5[7]={100,90,100,100,100,100,100}; // phi, cent, pt, EovP, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv5[7]={0,0,0,0,-1,-1,-1};
  Double_t xmaxv5[7]={2*TMath::Pi(),90,50,3,1,1,1}; 
  fChargPartV2 = new THnSparseD ("fChargPartV2","Charged particle v2",7,binsv5,xminv5,xmaxv5);
  fOutputList->Add(fChargPartV2);
  
  Int_t binsv6[5]={90,100,100,100,100}; // cent, pt, TPCdeltaPhi, V0AdeltaPhi, V0CdeltaPhi
  Double_t xminv6[5]={0,0,-1,-1,-1};
  Double_t xmaxv6[5]={90,50,1,1,1}; 
  feTPCV2 = new THnSparseD ("feTPCV2","inclusive electron v2 (TPC)",5,binsv6,xminv6,xmaxv6);
  fOutputList->Add(feTPCV2);
  
  //EovP,fTPCnSigma,partPt,fFlagPhotonicElec,fFlagPhotonicElecBCG,whichPart,cent,pt,firstMother,secondMother,thirdMother,iHijing,motherPt,secondMotherPt,thirdMotherPt
  Int_t binsv7[15]={100,100,100,3,3,5,90,100,5,5,5,3,100,100,100}; 
  Double_t xminv7[15]={0,-3.5,0,-1,-1,-1,0,0,-1,-1,-1,-1,0,0,0};
  Double_t xmaxv7[15]={3,3.5,50,2,2,4,90,50,4,4,4,2,50,50,50}; 
  fMCphotoElecPt = new THnSparseD ("fMCphotoElecPt", "pt distribution (MC)",15,binsv7,xminv7,xmaxv7);
  fOutputList->Add(fMCphotoElecPt);

  fGammaWeight = new TH2F("fGammaWeight", "Gamma weight",100,0,50,3,-1,2);
  fOutputList->Add(fGammaWeight);
 
  fPi0Weight = new TH2F("fPi0Weight", "Pi0 weight",100,0,50,3,-1,2);
  fOutputList->Add(fPi0Weight);
  
  fEtaWeight = new TH2F("fEtaWeight", "Eta weight",100,0,50,3,-1,2);
  fOutputList->Add(fEtaWeight);

  fD0Weight = new TH2F("fD0Weight", "D0 weight",100,0,50,3,-1,2);
  fOutputList->Add(fD0Weight);

  fDplusWeight = new TH2F("fDplusWeight", "D+ weight",100,0,50,3,-1,2);
  fOutputList->Add(fDplusWeight);

  fDminusWeight = new TH2F("fDminusWeight", "D- weight",100,0,50,3,-1,2);
  fOutputList->Add(fDminusWeight);

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

  if(mcPi0pT>0.0 && mcPi0pT<5.0){
	  if (cent>20.0 && cent<40.0) weight = (2.877091*mcPi0pT)/(TMath::Power(0.706963+mcPi0pT/3.179309,17.336628)*exp(-mcPi0pT));
		if (cent>30.0 && cent<50.0) weight = (2.392024*mcPi0pT)/(TMath::Power(0.688810+mcPi0pT/3.005145,16.811845)*exp(-mcPi0pT));
  }
  else{
    if (cent>20.0 && cent<40.0)	weight = (0.0004*mcPi0pT)/TMath::Power(-0.176181+mcPi0pT/3.989747,5.629235);
		if (cent>30.0 && cent<50.0) weight = (0.000186*mcPi0pT)/TMath::Power(-0.606279+mcPi0pT/3.158642,4.365540);
  }
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEtaweight(Double_t mcEtapT, Float_t cent) const
{
  //Get eta weight
  double weight = 1.0;

  if (cent>20.0 && cent<40.0) weight = (0.818052*mcEtapT)/(TMath::Power(0.358651+mcEtapT/2.878631,9.494043));
  if (cent>30.0 && cent<50.0) weight = (0.622703*mcEtapT)/(TMath::Power(0.323045+mcEtapT/2.736407,9.180356));
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetDweight(Int_t whichD, Double_t mcDpT, Float_t cent) const
{
  //get D weights
  double weight = 1.0;
    
  if (cent>30.0 && cent<50.0){
    if (whichD == 0) weight = 0.271583*TMath::Landau(mcDpT,3.807103,1.536753,0); // D
    if (whichD == 1) weight = 0.300771*TMath::Landau(mcDpT,3.725771,1.496980,0); // D0
    if (whichD == 2) weight = 0.247280*TMath::Landau(mcDpT,3.746811,1.607551,0); // D+
    if (whichD == 3) weight = 0.249410*TMath::Landau(mcDpT,3.611508,1.632196,0); //D-
  }
  return weight;
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
