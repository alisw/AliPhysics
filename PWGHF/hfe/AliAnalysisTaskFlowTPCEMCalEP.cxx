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

ClassImp(AliAnalysisTaskFlowTPCEMCalEP)

using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fAssPtCut(0.5)
  ,fAssTPCnCut(80)
  ,fAssITSrefitCut(kTRUE)
  ,fESD(0)
  ,fAOD(0)
  ,fVevent(0)
  ,fpidResponse(0)
  ,fMC(0)
  ,fStack(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fAssTrackCuts(0)
  ,fCuts(0)
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
  ,fInvmassCut(0.14)
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fWhichDecay(0)
  ,fPi0EtaWeight(1.)
  ,fCentAftThr(0)
  ,fTrigger(0)
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fCentAftFlt(0)
  ,fTPCsubEPres(0)
  ,fCorr(0)
  ,fElecMC(0)
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
    fevPlaneV0AftThr[k] = NULL;
    feTPCV2[k] = NULL;
    feV2[k] = NULL;
    fChargPartV2[k] = NULL;
    fMtcPartV2[k] = NULL;
    fEPres[k] = NULL;

    fPi0Pt[k] = NULL;
    fEtaPt[k] = NULL;
    fElecPtULSInvmassCut[k] = NULL;
    fElecPtLSInvmassCut[k] = NULL;
    fElecPtInvmassCut[k] = NULL;
    fInclElec[k] = NULL;
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
  fAssTrackCuts = new AliESDtrackCuts();
  
  InitParameters();
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
  ,fAssPtCut(0.5)
  ,fAssTPCnCut(80)
  ,fAssITSrefitCut(kTRUE)
  ,fESD(0)
  ,fAOD(0)
  ,fVevent(0)
  ,fpidResponse(0)
  ,fMC(0)
  ,fStack(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fAssTrackCuts(0)
  ,fCuts(0)
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
  ,fInvmassCut(0.14)	
  ,fChi2Cut(3.5)
  ,fDCAcut(999)
  ,fWhichDecay(0)
  ,fPi0EtaWeight(1.)
  ,fCentAftThr(0)
  ,fTrigger(0)
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fCentAftFlt(0)
  ,fTPCsubEPres(0)
  ,fCorr(0)
  ,fElecMC(0)
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
    fevPlaneV0AftThr[k] = NULL;
    feTPCV2[k] = NULL;
    feV2[k] = NULL;
    fChargPartV2[k] = NULL;
    fMtcPartV2[k] = NULL;
    fEPres[k] = NULL;

    fPi0Pt[k] = NULL;
    fEtaPt[k] = NULL;
    fElecPtULSInvmassCut[k] = NULL;
    fElecPtLSInvmassCut[k] = NULL;
    fElecPtInvmassCut[k] = NULL;
    fInclElec[k] = NULL;
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
  fAssTrackCuts = new AliESDtrackCuts();
  
  InitParameters();

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
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
  delete fAssTrackCuts;
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
  if(fIsMC && fMC) fStack = fMC->Stack();
 
  Int_t fNOtrks =fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  const AliESDVertex *spdVtx = fESD->GetPrimaryVertexSPD();

  Double_t pVtxZ = -999, spdVtxZ= -999;
  pVtxZ = pVtx->GetZ();
  spdVtxZ = spdVtx->GetZ();

  if(TMath::Abs(pVtxZ)>10) return;
  if(TMath::Abs(pVtxZ-spdVtxZ)>0.5) return;
  
  fNoEvents->Fill(0);

  if(fNOtrks<2) return;

  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }

  fPID->SetPIDResponse(fpidResponse);

  fCFM->SetRecEventInfo(fVevent);

  Float_t cent = -1., centTRK = -1.;
  Int_t iPt=8, iCent=3, iDeltaphi=4;
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");
  centTRK = centrality->GetCentralityPercentile("TRK");
  
  if(TMath::Abs(cent-centTRK)>5.) return;
  
  
  fCent->Fill(cent);
 
  if (cent>=0  && cent<10) iCent=0;
  if (cent>=10 && cent<20) iCent=1;
  if (cent>=20 && cent<40) iCent=2;
  if (cent<0 || cent>=40) return;

  // Random rejection of events in MB in 0-10% to flatten the centrality distribution
  
  Bool_t rejectEvent = kFALSE;
  Int_t centBin = fCent->FindBin(cent);
  Double_t centWeight = 1.;
  centWeight =GetCentWeight(centBin);
 
//   rejectEvent = RejectEvent(cent,centBin);
//   if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA && rejectEvent) return;
  
  // Trigger study and selection of V0 low threshold trigger in 10-20%

  fTrigger->Fill((Double_t)0,(Double_t)iCent);
  if (fESD->GetFiredTriggerClasses().Contains("CPBI2EGA-B-NOPF-CENTNOTRD")) fTrigger->Fill((Double_t)1,(Double_t)iCent); // kEMCEGA
  if (fESD->GetFiredTriggerClasses().Contains("CPBI2EJE-B-NOPF-CENTNOTRD")) fTrigger->Fill((Double_t)2,(Double_t)iCent); //kEMCEJE
  if (fESD->GetFiredTriggerClasses().Contains("CVHN-B-NOPF-ALLNOTRD") || fESD->GetFiredTriggerClasses().Contains("CVHN_R2-B-NOPF-ALLNOTRD") 
    || fESD->GetFiredTriggerClasses().Contains("CVLN_R1-B-NOPF-ALLNOTRD") || fESD->GetFiredTriggerClasses().Contains("CVLN-B-NOPF-ALLNOTRD")) fTrigger->Fill((Double_t)3,(Double_t)iCent); //kCentral or kSemiCentral
  if (fESD->GetFiredTriggerClasses().Contains("CPBI1-B-NOPF-ALLNOTRD") || fESD->GetFiredTriggerClasses().Contains("CPBI2_B1-B-NOPF-ALLNOTRD")) fTrigger->Fill((Double_t)4,(Double_t)iCent); //kMB
 
  Bool_t SelColl = kTRUE;
  if(GetCollisionCandidates()==AliVEvent::kAny)
  {
     SelColl = kFALSE;
     TString firedTrigger;
     firedTrigger = fESD->GetFiredTriggerClasses();
     if(firedTrigger.Contains("CVLN_B2-B-NOPF-ALLNOTRD") || firedTrigger.Contains("CVLN_R1-B-NOPF-ALLNOTRD") || firedTrigger.Contains("CSEMI_R1-B-NOPF-ALLNOTRD"))SelColl=kTRUE;
     if(!SelColl)return;
  }
  
  //Event planes

  Double_t evPlaneV0A = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0A",fESD,2));
  if(evPlaneV0A > TMath::Pi()) evPlaneV0A = evPlaneV0A - TMath::Pi();

  Double_t evPlaneV0C = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0C",fESD,2));
  if(evPlaneV0C > TMath::Pi()) evPlaneV0C = evPlaneV0C - TMath::Pi();

  Double_t evPlaneV0 = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0",fESD,2));
  if(evPlaneV0 > TMath::Pi()) evPlaneV0 = evPlaneV0 - TMath::Pi();
  
  Int_t epBin = fevPlaneV0[0]->FindBin(evPlaneV0);
  Double_t EPweight = 1.;
  EPweight = GetEPweight(epBin);

//   Bool_t rejectEventPlane = kFALSE;
//   rejectEventPlane = RejectEventPlane(evPlaneV0,epBin);
//   if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA && rejectEventPlane) return;
 
  Double_t wEvent = EPweight*centWeight;
  
  if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA ){
    fevPlaneV0[iCent]->Fill(evPlaneV0,wEvent);
    fCentAftFlt->Fill(cent,wEvent);
  }
  else fevPlaneV0[iCent]->Fill(evPlaneV0);
  
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

  // 2 subevent method (only for TPC EP)

  TVector2 *qsub1a = esdTPCep->GetQsub1();
  TVector2 *qsub2a = esdTPCep->GetQsub2();
  Double_t evPlaneResTPC = -999.;
  if(qsub1a && qsub2a){
    evPlaneResTPC = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }

  fTPCsubEPres->Fill(evPlaneResTPC);

  // 3 event method (V0, V0A, and V0C EP)

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

  if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA ){
    fEPres[0]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos),wEvent);
    fEPres[1]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg),wEvent);
    fEPres[2]->Fill(cent,GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg),wEvent);
  }
  else{
    fEPres[0]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos));
    fEPres[1]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg));
    fEPres[2]->Fill(cent,GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg));
  }
  
  // Selection of pi0 and eta in MC to compute the weight
  if(fIsMC && fMC && fStack){
    Int_t nParticles = fStack->GetNtrack();
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = fStack->Particle(iParticle);
      int fPDG = particle->GetPdgCode(); 
      double pTMC = particle->Pt();

      Double_t etaMC = particle->Eta();
      if (TMath::Abs(etaMC)>1.2)continue;

      Bool_t isMotherPrimary = IsPrimary(particle);
      Bool_t isFromLMdecay = IsFromLMdecay(particle);
      Bool_t isFromHFdecay = IsFromHFdecay(particle);

      if (isMotherPrimary || (!isFromHFdecay && !isFromLMdecay)){
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

  Bool_t IsSameEvent = kFALSE;
  
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
    if(fIsMC) fEMCalnSigma = GetSigmaEMCalMC(EovP, pt, iCent);
    else fEMCalnSigma = GetSigmaEMCal(EovP, pt, iCent);
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
    Double_t MCweight = 1.;
    Int_t iDecay = 0;

    // checking centrality and event plane distributions for events with electron above the trigger threshold
    if (fTPCnSigma>=-1 && fTPCnSigma<3 && fEMCalnSigma>0 && fEMCalnSigma<3 && pt>=8 && GetCollisionCandidates()==AliVEvent::kEMCEGA && !IsSameEvent){
	fevPlaneV0AftThr[iCent]->Fill(evPlaneV0);
	fCentAftThr->Fill(cent);
        IsSameEvent=kTRUE;
    }
    
    if(fTPCnSigma >= -1 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,EovP);
    
    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) continue;

    // MC part
    Int_t partPDG = -99;
    Double_t partPt = -99.;
    Bool_t MChijing; 

    if(fIsMC && fMC && fStack){
      if(m20<=0.02 && m02<=0.02) continue;
      Int_t label = track->GetLabel();
      if(label!=0){
        TParticle *particle = fStack->Particle(TMath::Abs(label));	
        if(particle){
          partPDG = particle->GetPdgCode();
          partPt = particle->Pt();
          
          MChijing = fMC->IsFromBGEvent(TMath::Abs(label));
          Int_t iHijing = 1;
          if(!MChijing) iHijing = 0; // 0 if enhanced sample

          GetWeightAndDecay(particle,iCent,iDecay,MCweight);

          Double_t corr[8]={cent,pt,fTPCnSigma,fEMCalnSigma,m02,dphi,cosdphi,static_cast<Double_t>(iDecay)};
          fCorr->Fill(corr,MCweight);
            
          if (TMath::Abs(partPDG)!=11) continue;

	      fInclElec[iCent]->Fill(pt,(Double_t)iDecay,MCweight);
	  	  
          SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,MCweight,iCent,iHijing,iDecay,fEMCalnSigma,fTPCnSigma);

        }// end particle
      }// end label
    }//end MC
        
    // data
    else{
      if(m20>0.02 && m02>0.02){ 
	Double_t corr[8]={cent,pt,fTPCnSigma,fEMCalnSigma,m02,dphi,cosdphi,0.};
	if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA ) fCorr->Fill(corr,wEvent);
	else fCorr->Fill(corr);
	
	SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,1,iCent,0,0,fEMCalnSigma,fTPCnSigma);
	fInclElec[iCent]->Fill(pt,0.);
      }
    }
    
    if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA ) fChargPartV2[iCent]->Fill(iPt,cosdphi,wEvent); 
    else fChargPartV2[iCent]->Fill(iPt,cosdphi); 
    
    if (clsE>0){
      if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA )  fMtcPartV2[iCent]->Fill(iPt,cosdphi,wEvent);
      else  fMtcPartV2[iCent]->Fill(iPt,cosdphi); 
    }
    
    if (fTPCnSigma>=-5 && fTPCnSigma<-3.2) fEoverPbcg[iCent][iPt][iDeltaphi]->Fill(EovP);
    if (fTPCnSigma>=-0.5 && fTPCnSigma<3) feTPCV2[iCent]->Fill(iPt,cosdphi); 
    if (fTPCnSigma>=-0.5 && fTPCnSigma<3 && fEMCalnSigma>-1 && fEMCalnSigma<3) feV2[iCent]->Fill(iPt,cosdphi); 

    fTrkEovPAft->Fill(pt,EovP);
    fdEdxAft->Fill(p,dEdx);

    if (m20>0.02 && m02>0.02 && fTPCnSigma>-1 && fTPCnSigma<3){ 
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
  
  /*if(!fPID->GetNumberOfPIDdetectors()) 
    {
      fPID->AddDetector("TPC", 0);
      fPID->AddDetector("EMCAL", 1);
    }*/
  
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
  
  int nbin_v2 = 7;
  double bin_v2[8] = {2,2.5,3,4,6,8,10,13};

  fCentAftThr = new TH1F("fCentAftThr","Centrality for electron pt > 8 GeV/c",100,0,100) ;
  fOutputList->Add(fCentAftThr);
  
  fTrigger = new TH2F("fTrigger","",5,0,5,3,0,3) ;
  fOutputList->Add(fTrigger);
  
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
  
  fCent = new TH1F("fCent","Centrality",100,0,100) ;
  fOutputList->Add(fCent);

  fCentAftFlt = new TH1F("fCentAftFlt","Centrality",100,0,100) ;
  fOutputList->Add(fCentAftFlt);
 
  fTPCsubEPres = new TH1F("fTPCsubEPres","TPC subevent plane resolution",100,-1,1);
  fOutputList->Add(fTPCsubEPres);
  
  //iCent,iPt,fTPCnSigma,fEMCalnSigma,m02,dphi,cosdphi,iDecay
  Int_t binsv2[8]={100,40,100,100,100,120,100,7}; 
  Double_t xminv2[8]={0,0,-5,-5,0,0,-1,0};
  Double_t xmaxv2[8]={100,20,5,5,2,TMath::Pi(),1,7}; 
  fCorr = new THnSparseD ("fCorr","Correlations",8,binsv2,xminv2,xmaxv2);
  fCorr->Sumw2();
  fOutputList->Add(fCorr);
  
  //iCent,pt,mass,fFlagLS,fFlagULS,iHijing,iDecay,,fEMCalnSigma,fTPCnSigma
  Int_t binsv3[9]={3,40,100,40,40,4,7,100,100}; 
  Double_t xminv3[9]={0,0,0,0,0,-1,0,-5,-5};
  Double_t xmaxv3[9]={3,20,0.3,20,20,3,7,5,5}; 
  fElecMC = new THnSparseD ("fElecMC","MC",9,binsv3,xminv3,xmaxv3);
  fElecMC->Sumw2();
  fOutputList->Add(fElecMC);
  
  for(Int_t i=0; i<3; i++) {
    fevPlaneV0[i] = new TH1F(Form("fevPlaneV0%d",i),"V0 EP",100,0,TMath::Pi());
    fevPlaneV0[i]->Sumw2();
    fOutputList->Add(fevPlaneV0[i]);

    fevPlaneV0AftThr[i] = new TH1F(Form("fevPlaneV0AftThr%d",i),"V0 EP for electron pt > 8 GeV/c",100,0,TMath::Pi());
    fevPlaneV0AftThr[i]->Sumw2();
    fOutputList->Add(fevPlaneV0AftThr[i]);
    
    feTPCV2[i] = new TH2F(Form("feTPCV2%d",i), "", 8,0,8,100,-1,1);
    feTPCV2[i]->Sumw2();
    fOutputList->Add(feTPCV2[i]);

    feV2[i] = new TH2F(Form("feV2%d",i), "", 8,0,8,100,-1,1);
    feV2[i]->Sumw2();
    fOutputList->Add(feV2[i]);

    fChargPartV2[i] = new TH2F(Form("fChargPartV2%d",i), "", 8,0,8,100,-1,1);
    fChargPartV2[i]->Sumw2();
    fOutputList->Add(fChargPartV2[i]);

    fMtcPartV2[i] = new TH2F(Form("fMtcPartV2%d",i), "", 8,0,8,100,-1,1);
    fMtcPartV2[i]->Sumw2();
    fOutputList->Add(fMtcPartV2[i]);
    
    fEPres[i] = new TH2F(Form("fEPres%d",i), "", 100,0,100,100,-1,1);
    fEPres[i]->Sumw2();
    fOutputList->Add(fEPres[i]);    

    fElecPtULSInvmassCut[i] = new TH2F(Form("fElecPtULSInvmassCut%d",i), "electron pt, ULS, invariant mass cut",nbin_v2,bin_v2,7,0,7);
    fElecPtULSInvmassCut[i]->Sumw2();
    fOutputList->Add(fElecPtULSInvmassCut[i]);

    fElecPtLSInvmassCut[i] = new TH2F(Form("fElecPtLSInvmassCut%d",i), "electron pt, LS, invariant mass cut",nbin_v2,bin_v2,7,0,7);
    fElecPtLSInvmassCut[i]->Sumw2();
    fOutputList->Add(fElecPtLSInvmassCut[i]);
    
    fElecPtInvmassCut[i] = new TH2F(Form("fElecPtInvmassCut%d",i), "electron pt, invariant mass cut",nbin_v2,bin_v2,7,0,7);
    fElecPtInvmassCut[i]->Sumw2();
    fOutputList->Add(fElecPtInvmassCut[i]);
    
    fInclElec[i] = new TH2F(Form("fInclElec%d",i), "inclusive electron pt", nbin_v2,bin_v2,7,0,7);
    fInclElec[i]->Sumw2();
    fOutputList->Add(fInclElec[i]);
    
    fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
    fInvmassLS[i]->Sumw2();
    fOutputList->Add(fInvmassLS[i]);

    fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
    fInvmassULS[i]->Sumw2();
    fOutputList->Add(fInvmassULS[i]);

    fOpeningAngleLS[i] = new TH2F(Form("fOpeningAngleLS%d",i),"Opening angle for LS pairs",100,0,1,100,0,50);
    fOpeningAngleLS[i]->Sumw2();
    fOutputList->Add(fOpeningAngleLS[i]);

    fOpeningAngleULS[i] = new TH2F(Form("fOpeningAngleULS%d",i),"Opening angle for ULS pairs",100,0,1,100,0,50);
    fOpeningAngleULS[i]->Sumw2();
    fOutputList->Add(fOpeningAngleULS[i]);

    fPi0Pt[i] = new TH1F(Form("fPi0Pt%d",i), "Pi0 weight",100,0,50);
    fPi0Pt[i]->Sumw2();
    fOutputList->Add(fPi0Pt[i]);

    fEtaPt[i] = new TH1F(Form("fEtaPt%d",i), "Eta weight",100,0,50);
    fEtaPt[i]->Sumw2();
    fOutputList->Add(fEtaPt[i]);
  }

  for(Int_t i=0; i<3; i++) {
    for(Int_t j=0; j<8; j++) {
      for(Int_t k=0; k<4; k++) {
        fEoverPsig[i][j][k] = new TH1F(Form("fEoverPsig%d%d%d",i,j,k), "",100,0,3);
	fEoverPsig[i][j][k]->Sumw2();
        fOutputList->Add(fEoverPsig[i][j][k]);

	fEoverPuls[i][j][k] = new TH1F(Form("fEoverPuls%d%d%d",i,j,k), "",100,0,3);
	fEoverPuls[i][j][k]->Sumw2();
        fOutputList->Add(fEoverPuls[i][j][k]);

	fEoverPls[i][j][k] = new TH1F(Form("fEoverPls%d%d%d",i,j,k), "",100,0,3);
	fEoverPls[i][j][k]->Sumw2();
        fOutputList->Add(fEoverPls[i][j][k]);

	fEoverPbcg[i][j][k] = new TH1F(Form("fEoverPbcg%d%d%d",i,j,k), "",100,0,3);
	fEoverPbcg[i][j][k]->Sumw2();
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

  fInvMass = new TH1F("fInvMass","",100,0,0.3);
  fOutputList->Add(fInvMass);

  fInvMassBack = new TH1F("fInvMassBack","",100,0,0.3);
  fOutputList->Add(fInvMassBack);

  fDCA = new TH1F("fDCA","",100,0,1);
  fOutputList->Add(fDCA);

  fDCABack = new TH1F("fDCABack","",100,0,1);
  fOutputList->Add(fDCABack);

  fOpAngle = new TH1F("fOpAngle","",100,0,0.5);
  fOutputList->Add(fOpAngle);

  fOpAngleBack = new TH1F("fOpAngleBack","",100,0,0.5);
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
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetPi0weight(Double_t mcPi0pT, Int_t iCent) const
{
  //Get Pi0 weight
  double weight = 1.0;
    
  if (iCent==0){
     double parLowPt[4] = {0.00386062,0.913163,54.9096,84.0142};
     double parHighPt[4] = {0.02824,0.1246,3.56903,5.60296};
    
    if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
    if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
  }
  if (iCent==1){
     double parLowPt[4] = {0.000197581,0.960398,174.617,269.436};
     double parHighPt[4] = {0.0435973,0.0732613,3.43644,5.56708};
    
    if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
    if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);  
  }
  if (iCent==2){
     double parLowPt[4] = {0.00395183,0.905839,52.3325,78.9736};
     double parHighPt[4] = {0.0639772,0.0954623,3.21665,5.67225};

    if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
    if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
  }
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEtaweight(Double_t mcEtapT, Int_t iCent) const
{
  //Get eta weight
  double weight = 1.0;
    
  if (iCent==0){
     double parLowPt[4] = {0.00218816,0.903496,52.9872,73.6404};
     double parHighPt[4] = {0.0742314,0.296077,3.33914,5.86723};
    
    if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
    if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);  
  }
  if (iCent==1){
     double parLowPt[4] = {0.00218739,0.889904,49.5944,66.8576};
     double parHighPt[4] = {0.124957,0.216147,3.09109,5.76838};

    if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
    if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);  
  }
  if (iCent==2){
     double parLowPt[4] = {0.00326269,0.911628,57.3255,78.6103};
     double parHighPt[4] = {0.134015,0.207723,3.00919,5.83206};

    if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
    if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);  
  }
    
  return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetSigmaEMCal(Double_t EoverP, Double_t pt, Int_t iCent) const
{
  //Get sigma for EMCal PID
  Double_t NumberOfSigmasEMCal = 99.;
  Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};

  if (iCent==0){
     Double_t mean[8]={1.04892,1.04471,1.04397,1.04715,1.04617,1.04147,1.05363,1.04902};
     Double_t sigma[8]={0.157249,0.143196,0.130118,0.118269,0.105798,0.100816,0.0910207,0.0953318};   
     for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  if (iCent==1){
    Double_t mean[8]={1.01201,1.01064,1.01248,1.01728,1.02346,1.02177,1.04038,1.03314};
    Double_t sigma[8]={0.144614,0.126229,0.120568,0.107897,0.0919854,0.0920917,0.0859356,0.085302};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  if (iCent==2){
    Double_t mean[8]={0.975778,0.975963,0.983835,0.988513,0.999726,1.00552,1.01144,1.00319};
    Double_t sigma[8]={0.130389,0.117007,0.10375,0.0971151,0.0893869,0.0873147,0.083138,0.0874688};   
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  return NumberOfSigmasEMCal;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetSigmaEMCalMC(Double_t EoverP, Double_t pt, Int_t iCent) const
{
  //Get sigma for EMCal PID
  Double_t NumberOfSigmasEMCal = 99.;
  Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};

  if (iCent==0){
    Double_t mean[8]={1.01076,1.00735,1.00386,1.00281,1.00114,0.998282,0.995936,0.998286};
    Double_t sigma[8]={0.153704,0.137907,0.127886,0.115947,0.102482,0.0921989,0.0896079,0.0944837};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  if (iCent==1){
    Double_t mean[8]={0.97531,0.973007,0.971888,0.972424,0.97437,0.976057,0.977703,0.984494};
    Double_t sigma[8]={0.132568,0.119308,0.107527,0.099176,0.0873851,0.0779302,0.0779114,0.0834648};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  if (iCent==2){
    Double_t mean[8]={0.954379,0.952449,0.952901,0.955364,0.961415,0.965205,0.968959,0.976448};
    Double_t sigma[8]={0.120315,0.106597,0.0968691,0.0879189,0.0784124,0.0719245,0.0704888,0.080023};
    for(Int_t i=0;i<8;i++) {
      if (pt>=ptRange[i] && pt<ptRange[i+1]){
        NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
        continue;
      }
    }
  }
  return NumberOfSigmasEMCal;
}
//________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::GetWeightAndDecay(TParticle *particle, Int_t iCent, Int_t &decay, Double_t &weight)
{
  //Get pi0/eta weight for MC with enchanced signal and decay channel
  Double_t w = 1.;
  Int_t d = 0;
  Int_t partPDG = particle->GetPdgCode();
    
  if (TMath::Abs(partPDG)==11){
      Int_t idMother = particle->GetFirstMother();
      
      if (idMother>0){
          TParticle *mother = fStack->Particle(idMother);
          Int_t motherPDG = mother->GetPdgCode();
          Double_t motherPt = mother->Pt();
          
          Bool_t isMotherPrimary = IsPrimary(mother);
          Bool_t isMotherFromHF = IsFromHFdecay(mother);
          Bool_t isMotherFromLM = IsFromLMdecay(mother);
          
          if (motherPDG==111 && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // pi0 -> e
              d = 1;
              w = GetPi0weight(motherPt,iCent);
          }
          
          if (motherPDG==221  && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // eta -> e
              d = 2;
              w = GetEtaweight(motherPt,iCent);
          }
          
          //Int_t idSecondMother = particle->GetSecondMother();
          Int_t idSecondMother = mother->GetFirstMother();
          
          if (idSecondMother>0){
              TParticle *secondMother = fStack->Particle(idSecondMother);
              Int_t secondMotherPDG = secondMother->GetPdgCode();
              Double_t secondMotherPt = secondMother->Pt();
              
              Bool_t isSecondMotherPrimary = IsPrimary(secondMother);
              Bool_t isSecondMotherFromHF = IsFromHFdecay(secondMother);
              Bool_t isSecondMotherFromLM = IsFromLMdecay(secondMother);
              
              if (motherPDG==22 && secondMotherPDG==111 && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //pi0 -> g -> e
                  d = 3;
                  w = GetPi0weight(secondMotherPt,iCent);
              }
              
              if (motherPDG==22 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> g -> e
                  d = 4;
                  w = GetEtaweight(secondMotherPt,iCent);
              }
              
              if (motherPDG==111 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> pi0 -> e
                  d = 5;
                  w = GetEtaweight(secondMotherPt,iCent);
              }
              
              Int_t idThirdMother = secondMother->GetFirstMother();
              if (idThirdMother>0){
                  TParticle *thirdMother = fStack->Particle(idThirdMother);
                  Int_t thirdMotherPDG = thirdMother->GetPdgCode();
                  Double_t thirdMotherPt = thirdMother->Pt();
                  
                  Bool_t isThirdMotherPrimary = IsPrimary(thirdMother);
                  Bool_t isThirdMotherFromHF = IsFromHFdecay(thirdMother);
                  Bool_t isThirdMotherFromLM = IsFromLMdecay(thirdMother);
                  
                  if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221 && (isThirdMotherPrimary || (!isThirdMotherFromHF && !isThirdMotherFromLM))){//eta->pi0->g-> e 
                      d = 6;
                      w = GetEtaweight(thirdMotherPt,iCent);
                  }
              }//third mother      
          }//second mother 
      }//mother
  }// if electron
  decay = d;
  weight = w;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetCentWeight(Int_t centbin){
  // Get cebtrality weight for flattening (0-10%)
  Int_t wBin = centbin-1;
  if (wBin<0 || wBin>9) return 1;
  Double_t weightcent[] = {0.996425,0.987564,0.966774,0.966422,0.967739,0.991296,0.983746,0.991528,1.01627,1.15803};
    
  return weightcent[wBin];
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEPweight(Int_t bin) 
{
  //Get event plane weight for flattening (0-10%)
  Int_t wBin = bin-1;
  if (wBin<0 || wBin>99) return 1;
  
  Double_t weightEP[] = {0.982991,0.988171,0.9899237,0.9914497,0.9906325,0.9956888,0.9972689,1.000973,1.002418,1.006948,1.007226,1.008336,1.01335,1.011154,1.018333,
  1.019898,1.026543,1.023092,1.028325,1.026844,1.031437,1.03014,1.031728,1.030307,1.037547,1.03471,1.03722,1.039466,1.037632,1.041682,1.042824,1.037494,1.046057,
  1.046622,1.042124,1.043161,1.040339,1.040997,1.043782,1.039092,1.039026,1.033509,1.035641,1.034528,1.031159,1.029701,1.033969,1.021809,1.02614,1.017396,1.017012,
  1.013525,1.012976,1.007164,1.006868,1.00653,0.9983816,0.9962069,0.9987208,0.9958153,0.9902154,0.9837839,0.9805614,0.9825041,0.9821056,0.9785275,0.9793774,0.9739373,
  0.9722809,0.9728094,0.972367,0.9687113,0.96755,0.9635185,0.9605392,0.9610214,0.9614648,0.9591571,0.9603319,0.9610102,0.9675955,0.9609205,0.9605896,0.9625102,0.9589448,
  0.9624427,0.966783,0.9632197,0.9626284,0.9706073,0.9693101,0.9717702,0.9703041,0.9747158,0.9741852,0.9755416,0.9798203,0.9797912,0.9790047,0.9802287};
  
  return weightEP[wBin];
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::RejectEvent(Double_t cent, Int_t centbin)
{
  // Reject randomly event in 0-10% in order to flatten the centrality distribution in MB events  
  Int_t wBin = centbin-1;
  if (wBin<0 || wBin>9) return kFALSE;
  
  Double_t weight[] = {0.858984,0.853393,0.835121,0.834615,0.835851,0.855657,0.849417,0.856341,0.877473,1};
    
  Double_t centDigits=cent-(Int_t)(cent*100.)/100.;
  
  if(centDigits*100.>weight[wBin]) return kTRUE; 
  
  return kFALSE;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::RejectEventPlane(Double_t EP, Int_t EPbin)
{
  // Reject randomly event plane in 0-10% in order to flatten the EP distribution in MB events  
  Int_t wBin = EPbin-1;
  if (wBin<0 || wBin>99) return kFALSE;
  
  Double_t weight[] = {0.939203,0.944153,0.945827,0.947285,0.946505,0.951336,0.952845,0.956385,0.957765,0.962094,0.962358,0.963419,0.968211,0.966111,0.972971,0.974466,0.980815,0.977518,0.982518,0.981103,0.985491,0.984252,0.985769,0.984412,0.991329,0.988619,0.991017,0.993163,0.99141,0.99528,0.996371,0.991279,0.99946,1,0.995702,0.996693,0.993997,0.994626,0.997287,0.992805,0.992743,0.987471,0.989508,0.988445,0.985226,0.983833,0.987911,0.976292,0.98043,0.972075,0.971709,0.968377,0.967853,0.9623,0.962017,0.961694,0.953908,0.951831,0.954233,0.951456,0.946106,0.939961,0.936882,0.938738,0.938357,0.934939,0.935751,0.930553,0.92897,0.929475,0.929053,0.92556,0.92445,0.920598,0.917752,0.918213,0.918636,0.916431,0.917554,0.918202,0.924494,0.918116,0.9178,0.919635,0.916228,0.91957,0.923717,0.920313,0.919748,0.927371,0.926132,0.928482,0.927082,0.931297,0.93079,0.932086,0.936174,0.936146,0.935395,0.936564};
    
  Double_t centDigits=EP-(Int_t)(EP*100.)/100.;
  
  if(centDigits*100.>weight[wBin]) return kTRUE; 
  
  return kFALSE;
}

//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsFromHFdecay(TParticle *particle) 
{
  // Check if the mother comes from heavy-flavour decays
  Bool_t isHFdecay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = fStack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();

    // c decays
    if((TMath::Abs(motherPDG)==411) || (TMath::Abs(motherPDG)==421) || (TMath::Abs(motherPDG)==431) || (TMath::Abs(motherPDG)==4122) || (TMath::Abs(motherPDG)==4132) || (TMath::Abs(motherPDG)==4232) || (TMath::Abs(motherPDG)==43320)) isHFdecay = kTRUE; 

    // b decays
    if((TMath::Abs(motherPDG)==511) || (TMath::Abs(motherPDG)==521) || (TMath::Abs(motherPDG)==531) || (TMath::Abs(motherPDG)==5122) || (TMath::Abs(motherPDG)==5132) || (TMath::Abs(motherPDG)==5232) || (TMath::Abs(motherPDG)==53320)) isHFdecay = kTRUE; 
  }

  return isHFdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsFromLMdecay(TParticle *particle) 
{
  // Check if the mother comes from light-meson decays
  Bool_t isLMdecay = kFALSE;
  Int_t partPDG = particle->GetPdgCode();

  Int_t idMother = particle->GetFirstMother();
  if (idMother>0){
    TParticle *mother = fStack->Particle(idMother);
    Int_t motherPDG = mother->GetPdgCode();

    if(motherPDG == 111 || motherPDG == 221 || motherPDG==223 || motherPDG==333 || motherPDG==331 || (TMath::Abs(motherPDG)==113) || (TMath::Abs(motherPDG)==213) || (TMath::Abs(motherPDG)==313) || (TMath::Abs(motherPDG)==323)) isLMdecay = kTRUE;
  }

  return isLMdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsPrimary(TParticle *particle) 
{
  // Check if the particle is primary
  Bool_t isprimary = kFALSE;
  if (particle->IsPrimary()) isprimary = kTRUE;  
  
  return isprimary;
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::SelectPhotonicElectron(Int_t iTracks,AliESDtrack *track,Bool_t &fFlagPhotonicElec, Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent, Int_t iHijing, Int_t iDecay, Double_t fEMCalnSigma, Double_t fTPCnSigma)
{
  //Identify non-heavy flavour electrons using Invariant mass method
  fAssTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fAssTrackCuts->SetRequireTPCRefit(kTRUE);
  fAssTrackCuts->SetRequireITSRefit(fAssITSrefitCut);
  fAssTrackCuts->SetEtaRange(-0.9,0.9);
  fAssTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fAssTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fAssTrackCuts->SetMinNClustersTPC(fAssTPCnCut);
  fAssTrackCuts->SetMinNClustersITS(2);
  fAssTrackCuts->SetMaxDCAToVertexZ(3.2);
  fAssTrackCuts->SetMaxDCAToVertexXY(2.4);
  fAssTrackCuts->SetDCAToVertex2D(kTRUE);
  
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
    Int_t chargeAsso = 0, charge = 0, pdgAsso = 0;
    
    nTPCsigmaAsso = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
    pt = track->Pt();
    ptAsso = trackAsso->Pt();
    chargeAsso = trackAsso->Charge();
    charge = track->Charge();
        
    if(ptAsso <fAssPtCut) continue;
    if(!fAssTrackCuts->AcceptTrack(trackAsso)) continue;
    if(TMath::Abs(nTPCsigmaAsso)>3) continue;
    
    /*if(fIsMC && fMC && fStack){
      Int_t labelAsso = trackAsso->GetLabel();
      if(labelAsso!=0){
	TParticle *particleAsso = fStack->Particle(TMath::Abs(labelAsso));	
	if(particleAsso){
	  pdgAsso = particleAsso->GetPdgCode();
	  if(!(TMath::Abs(pdgAsso)==11)) continue;
	}
      }
    }*/
    
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

    if(fFlagLS) fOpeningAngleLS[iCent]->Fill(openingAngle,ptAsso);
    if(fFlagULS) fOpeningAngleULS[iCent]->Fill(openingAngle,ptAsso);

    //if(openingAngle > fOpeningAngleCut) continue;
    
    recg.GetMass(mass,width);
    
    Double_t elecMC[9]={(Double_t)iCent,pt,mass,(Double_t)fFlagLS,(Double_t)fFlagULS,(Double_t)iHijing,(Double_t)iDecay, fEMCalnSigma, fTPCnSigma};
    fElecMC->Fill(elecMC,weight);
    
    if(fFlagLS) fInvmassLS[iCent]->Fill(mass,pt,weight);
    if(fFlagULS) fInvmassULS[iCent]->Fill(mass,pt,weight);
    
    if(mass<fInvmassCut) fElecPtInvmassCut[iCent]->Fill(pt,iDecay,weight);
    if(mass<fInvmassCut && fFlagULS) fElecPtULSInvmassCut[iCent]->Fill(pt,iDecay,weight);
    if(mass<fInvmassCut && fFlagLS) fElecPtLSInvmassCut[iCent]->Fill(pt,iDecay,weight);
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
}
