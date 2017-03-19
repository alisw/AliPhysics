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


///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskHypCrossCheck class
// analysis task for the make some base checks for hypertriton analysis
//
// This task is optimized for ESDs.root
//
// Author:
// S. Trogolo, trogolo@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
//#include <vector>

#include <TArray.h>
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHypCrossCheck.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliVertexerTracks.h"
#include "AliVEvent.h"
#include "AliVTrack.h"


ClassImp(AliAnalysisTaskHypCrossCheck)

using std::vector;

/* $Id$ */
//________________________________________________________________________
AliAnalysisTaskHypCrossCheck::AliAnalysisTaskHypCrossCheck(TString taskname):
  AliAnalysisTaskSE(taskname.Data()),
  fESDevent(0),
  fESDtrackCuts(0x0),
  fESDtrackCutsV0(0x0),
  fPrimaryVertex(0x0),
  fPIDResponse(0x0),
  fStack(0x0),
  fVertexer(0x0),
  fVtx1(0x0),
  fVtx2(0x0),
  fTrkArray(0x0),
  fCustomTPCpid(0),
  fMC(kFALSE),
  fFillTree(kFALSE),
  fFillTGen(kFALSE),
  fRun1PbPb(kFALSE),
  fRun2PbPb(kFALSE),
  fCentrality(0x0),
  fCentralityPercentile(0x0),
  fTriggerConfig(1),
  fEvtEmbedSelection(kTRUE),
  fEvtSpecie(4),
  fRequestITSrefit(kFALSE),
  fRequestITSrefitPion(kFALSE),
  fRequestTPCSigmas(3),
  fRequestTOFPid(kFALSE),
  fRequestTOFSigmas(3),
  fChooseMatter(kTRUE),
  fChooseAntiMatter(kTRUE),
  fMinvSignal(kTRUE),
  fMinvLikeSign(kFALSE),
  fSideBand(kFALSE),
  fDCAPiPVmin(0.4),
  fDCAzHe3PVmax(1.), //999.
  fDCAhe3pi(0.7), //0.2
  fMaxDecayLength(999.),
  fMinDecayLength(0.),
  fDCAPiSVxymax(999.), //0.6
  fDCAPiSVzmax(999.), //0.8
  fDCAHe3SVmax(999.), //0.6
  fMinNormalizedDecL(0.), //0.
  fMaxPtMother(10.), //11.
  fMinPtMother(0.), //0.
  fMaxLifeTime(999.),
  fMinLifeTime(1.), //0.
  fRapidity(0.9), //1.
  fCosPointingAngle(0.99), //0.99
  fAnglehe3pi(TMath::Pi()), //TMath::Pi()/2
  fLowCentrality(0.),
  fHighCentrality(90.),
  fOutput(0x0),
  fHistCount(0x0),
  fHistCentralityClass(0x0),
  fHistCentralityPercentile(0x0),
  fHistTrigger(0x0),
  fHistMultiplicity(0x0),
  fHistZPrimaryVtx(0x0),
  fHistXPrimaryVtx(0x0),
  fHistYPrimaryVtx(0x0),
  fHistChi2perTPCcluster(0x0),
  fHistTrackFlagReco(0x0),
  fHistTPCpid(0x0),
  fHistTPCHe3signal(0x0),
  fHistTPCpionsignal(0x0),
  fHistNsigmaHe3(0x0),
  fHistNsigmaPion(0x0),
  fHistTOFsignal(0x0),
  fHistTOFHe3signal(0x0),
  fHistpionTPCcls(0x0),
  //fHistCorrDCAHe3primary(0x0),
  //fHistCorrDCApiprimary(0x0),
  fHistDCApiprimary(0x0),
  fHistDCAXYpiprimary(0x0),
  fHistDCAZpiprimary(0x0),
  fHistDCAHe3primary(0x0),
  fHistDCAXYHe3primary(0x0),
  fHistDCAZHe3primary(0x0),
  fHistDCAhe3pion(0x0),
  fHistZDecayVtx(0x0),
  fHistXDecayVtx(0x0),
  fHistYDecayVtx(0x0),
  fHistDCAXYhe3vtx(0x0),
  fHistDCAZhe3vtx(0x0),
  fHistDCAhe3vtx(0x0),
  fHistDCAXYpionvtx(0x0),
  fHistDCAZpionvtx(0x0),
  fHistDCApionvtx(0x0),
  fHistDecayLengthH3L(0x0),
  fHistNormalizedDecayL(0x0),
  fHistDeltaPt_Hyper(0x0),
  fHistDeltaPt_Pion(0x0),
  fHistDeltaPt_He3(0x0),
  fHistLifetime(0x0),
  fHistAngle_He3_pion(0x0),
  fHistHyperRapidity(0x0),
  fHistCosPointingAngle(0x0),
  fHistPtPion(0x0),
  fHistPtHelium3(0x0),
  fHistPtHypertriton(0x0),
  fHistMassHypertriton(0x0),
  fHistMassAntiHypertriton(0x0),
  fHistMassVsPt(0x0),
  fHistMassHyperVsPt(0x0),
  fHistMassAntiHyperVsPt(0x0),
  fHistTPCdeusignal_pdg(0x0),
  fHistTPCtrisignal_pdg(0x0),
  fHistTPCHe3signal_pdg(0x0),
  fHistTPCHe4signal_pdg(0x0),
  fHistTPCHe3signal_3LH(0x0),
  fHistTPCpionsignal_3LH(0x0),
  fHistNsigmaHe3_3LH(0x0),
  fHistNsigmaPion_3LH(0x0),
  fHistParticle(0x0),
  fHistParticle_Mass(0x0),
  fHistpionTPCclsMCt(0x0),
  fHisthelium3TPCclsMCt(0x0),
  fHistpTpionMCt(0x0),
  fHistpThe3MCt(0x0),
  fHistMompionMCt(0x0),
  fHistMomHe3MCt(0x0),
  fHistCorrDCAHe3primaryMCt(0x0),
  fHistCorrDCApiprimaryMCt(0x0),
  fHistDCApiprimaryMCt(0x0),
  fHistDCAXYpiprimaryMCt(0x0),
  fHistDCAZpiprimaryMCt(0x0),
  fHistDCAHe3primaryMCt(0x0),
  fHistDCAXYHe3primaryMCt(0x0),
  fHistDCAZHe3primaryMCt(0x0),
  fHistDCAhe3pionMCt(0x0),
  fHistDeltaPt_PionMCt(0x0),
  fHistDeltaPt_He3MCt(0x0),
  fHistZDecayVtxMCt(0x0),
  fHistXDecayVtxMCt(0x0),
  fHistYDecayVtxMCt(0x0),
  fHistDecayLengthH3L_MCt(0x0),
  fHistNormalizedDecayL_MCt(0x0),
  fHistDCAXYhe3vtxMCt(0x0),
  fHistDCAZhe3vtxMCt(0x0),
  fHistDCAhe3vtxMCt(0x0),
  fHistDCAXYpionvtxMCt(0x0),
  fHistDCAZpionvtxMCt(0x0),
  fHistDCApionvtxMCt(0x0),
  fHistDeltaPt_HyperMCt(0x0),
  fHistLifetime_MCt(0x0),
  fHistAngle_He3_pion_MCt(0x0),
  fHistHypertritonMomMCt(0x0),
  fHistPtPionMCt(0x0),
  fHistPtHelium3MCt(0x0),
  fHistPtHypertritonMCt(0x0),
  fHistHyperRapidityMCt(0x0),
  fHistCosPointingAngleMCt(0x0),
  fHistMassHypertritonMCt(0x0),
  fHistMassAntiHypertritonMCt(0x0),
  fHistMassVsPtMCt(0x0),
  fHistMassHyperVsPtMCt(0x0),
  fHistMassAntiHyperVsPtMCt(0x0),
  fHistHypertritonMomGen(0x0),
  fHistHypertritonMomGen_3Body(0x0),
  fHistAntiHypertritonMomGen_3Body(0x0),
  fHistHypertritonMomGen_isPrimary_3Body(0x0),
  fHistHypertritonMomGen_2Body(0x0),
  fHistAntiHypertritonMomGen_2Body(0x0),
  fHistHypertritonMomGen_isPrimary_2Body(0x0),
  fHistHypertritonYGen(0x0),
  fHistHypertritonYGen_3Body(0x0),
  fHistAntiHypertritonYGen_3Body(0x0),
  fHistHypertritonYGen_isPrimary_3Body(0x0),
  fHistHypertritonYGen_2Body(0x0),
  fHistAntiHypertritonYGen_2Body(0x0),
  fHistHypertritonYGen_isPrimary_2Body(0x0),
  fTTree(0x0),
  fTMCtruth(kFALSE),
  fTCentralityPerc(0x0),
  fTchi2NDFhe3(0x0),
  fTMhypoTrkhe3(0x0),
  fTPCclsPIDhe3(0x0),
  fITSrefithe3(0x0),
  fTpTPChe3(0x0),
  fTpXhe3(0x0),
  fTpYhe3(0x0),
  fTpZhe3(0x0),
  fTTPCnsigmahe3(0x0),
  fTDCAXYhe3prvtx(0x0),
  fTDCAZhe3prvtx(0x0),
  fTchi2NDFpion(0x0),
  fTMhypoTrkpion(0x0),
  fTPCclsPIDpion(0x0),
  fITSrefitpion(0x0),
  fTpTPCpion(0x0),
  fTpXpion(0x0),
  fTpYpion(0x0),
  fTpZpion(0x0),
  fTTPCnsigmapion(0x0),
  fTDCAXYpioprvtx(0x0),
  fTDCAZpioprvtx(0x0),
  fTDCAhe3pi(0x0),
  fTDCAXYhevtx(0x0),
  fTDCAZhevtx(0x0),
  fTDCAXYpivtx(0x0),
  fTDCAZpivtx(0x0),
  fTAngle_he3pi(0x0),
  fTphe3_gen_X(0x0),
  fTphe3_gen_Y(0x0),
  fTphe3_gen_Z(0x0),
  fTppio_gen_X(0x0),
  fTppio_gen_Y(0x0),
  fTppio_gen_Z(0x0),
  fTpdgHe3(0x0),
  fTpdgPion(0x0),
  fTmomidHe3(0x0),
  fTmomidPi(0x0),
  fTpdgmomHe3(0x0),
  fTpdgmomPi(0x0),
  fTuniqID_he3(0x0),
  fTuniqID_pion(0x0),
  fTRapidity(0x0),
  fTDecayLength(0x0),
  fTDecayLengthError(0x0),
  fTCosPA(0x0),
  fTInvariantMass(0x0),
  fTGen(0x0),
  fTCentrality_gen(0x0),
  fTHyper_px(0x0),
  fTHyper_py(0x0),
  fTHyper_pz(0x0),
  fTHyper_eta(0x0),
  fTHyper_rapidity(0x0),
  fTHyper_pv_x(0x0),
  fTHyper_pv_y(0x0),
  fTHyper_pv_z(0x0),
  fTHyper_charge(0x0),
  fTHe3_px(0x0),
  fTHe3_py(0x0),
  fTHe3_pz(0x0),
  fTHe3_eta(0x0),
  fTHe3_rapidity(0x0),
  fTHe3_pv_x(0x0),
  fTHe3_pv_y(0x0),
  fTHe3_pv_z(0x0),
  fTHe3_charge(0x0),
  fTPion_px(0x0),
  fTPion_py(0x0),
  fTPion_pz(0x0),
  fTPion_eta(0x0),
  fTPion_rapidity(0x0),
  fTPion_charge(0x0)
{
  //Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());

  //ESD Track cuts
  if(!fESDtrackCuts) fESDtrackCuts = new AliESDtrackCuts();
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(fRequestITSrefit);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);

  //ESD Track cuts V0
  if(!fESDtrackCutsV0) fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsV0->SetMinNClustersTPC(80);
  fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
  fESDtrackCutsV0->SetRequireITSRefit(fRequestITSrefitPion);
  fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
  fESDtrackCutsV0->SetPtRange(0.2,1.5);

}


//________________________________________________________________________
AliAnalysisTaskHypCrossCheck::~AliAnalysisTaskHypCrossCheck(){
  //Destructor
  if(fOutput){
    delete fOutput;
    fOutput = 0;
  }

  if(fTTree) delete fTTree;
  if(fTGen) delete fTGen;

    if(fPIDResponse){
    delete fPIDResponse;
  }

    if(fESDtrackCuts) delete fESDtrackCuts;
    if(fESDtrackCutsV0) delete fESDtrackCutsV0;
    if(fPrimaryVertex) delete fPrimaryVertex;
    if(fVertexer) delete fVertexer;
    if(fStack) delete fStack;
    if(fTrkArray) delete fTrkArray;
    if(fVtx1) delete fVtx1;
    if(fVtx2) delete fVtx2;


} // end of Destructor

//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::PassTriggerSelection(UInt_t PhysSelMask){

  if(!fRun1PbPb && !fRun2PbPb){
    AliWarning("WARNING: No Run(1-2) Pb-Pb period selection =====>>> Please choose one!");
    return kFALSE;
  }

  if(fRun1PbPb && fRun2PbPb){
    AliWarning("WARNING: Run1 AND Run2 Pb-Pb period selected =====>>> Please choose only one!");
    return kFALSE;
  }

  if(fRun1PbPb){ //trigger for 2011 Pb-Pb analysis
    Bool_t isSelectedCentral = (PhysSelMask & AliVEvent::kCentral);
    Bool_t isSelectedSemiCentral = (PhysSelMask & AliVEvent::kSemiCentral);
    Bool_t isSelectedMB = (PhysSelMask & AliVEvent::kMB);

    if(isSelectedMB) { // Minimum Bias
      fHistTrigger->Fill(0);
    }
    if(isSelectedCentral) { // Central
      fHistTrigger->Fill(1);
    }
    if(isSelectedSemiCentral) { // SemiCentral
      fHistTrigger->Fill(2);
    }

    if(fTriggerConfig==1 && !isSelectedMB && !isSelectedCentral && !isSelectedSemiCentral){ //kMB + kCentral + kSemiCentral
      return kFALSE;
    }else if(fTriggerConfig==2 && !isSelectedCentral && !isSelectedSemiCentral){ //kCentral + kSemiCentral
      return kFALSE;
    }else if(fTriggerConfig==3 && !isSelectedSemiCentral){ //kSemiCentral
      return kFALSE;
    }else if(fTriggerConfig == 4 && !isSelectedCentral){ //kCentral
        return kFALSE;
    }
    return kTRUE;
  } // end of fRun1PbPb

  if(fRun2PbPb){ //trigger for 2015 Pb-Pb analysis
    Bool_t isSelectedMB = (PhysSelMask & AliVEvent::kINT7); // V0AND minimum bias trigger
    if(isSelectedMB) { // Minimum Bias
      fHistTrigger->Fill(0);
    }

    if(!isSelectedMB) return kFALSE;
    return kTRUE;
  } // end of fRun2PbPb

  return kFALSE;

}


//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::PassCentralitySelection(){

  if(!fRun1PbPb && !fRun2PbPb){
    AliWarning("WARNING: No Run(1-2) Pb-Pb period selection =====>>> Please choose one!");
    return kFALSE;
  }

  if(fRun1PbPb && fRun2PbPb){
    AliWarning("WARNING: Run1 AND Run2 Pb-Pb period selected =====>>> Please choose only one!");
    return kFALSE;
  }

  if(fRun1PbPb){
    if(fESDevent->GetEventSpecie() == fEvtSpecie){ // Event Specie == 4 == PbPb
      AliWarning("fRun1PbPb: Centrality task");
      AliCentrality *centr=fESDevent->GetCentrality();
      fCentrality = centr->GetCentralityClass10("V0M");
      fCentralityPercentile = centr->GetCentralityPercentile("V0M");
    }
  }
  if(fRun2PbPb){
    AliWarning("fRun2PbPb: Centrality task");
    AliMultSelection *centr=(AliMultSelection*)fESDevent->FindListObject("MultSelection");
    if(!centr){
      AliWarning("AliMultSelection object not found!");
    }
    fCentralityPercentile = centr->GetMultiplicityPercentile("V0M",fEvtEmbedSelection);
  }

  fHistCentralityClass->Fill(fCentrality);
  fHistCentralityPercentile->Fill(fCentralityPercentile);

  if (fCentralityPercentile < fLowCentrality || fCentralityPercentile > fHighCentrality) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::SetConvertedAODVertices(AliESDVertex *ESDvtxp, AliESDVertex *ESDvtxs)const{

  Double_t pos[3], cov[6], chi2perNDF;

  //Conversion of the primary vertex
  ESDvtxp->GetXYZ(pos);
  ESDvtxp->GetCovMatrix(cov);
  chi2perNDF = ESDvtxp->GetChi2toNDF();

  fVtx1->SetPosition(pos[0], pos[1], pos[2]);
  fVtx1->SetCovMatrix(cov);
  fVtx1->SetChi2perNDF(chi2perNDF);
  fVtx1->SetID(ESDvtxp->GetID());
  fVtx1->SetType(AliAODVertex::kPrimary);


  pos[0] = 0.; pos[1] = 0.; pos[2] = 0.;
  cov[0] = 0.; cov[1] = 0.; cov[2] = 0.; cov[3] = 0.; cov[4] = 0.; cov[5] = 0.;
  chi2perNDF = 0.;

  //Conversion of the secondary vertex
  ESDvtxs->GetXYZ(pos);
  ESDvtxs->GetCovMatrix(cov);
  chi2perNDF = ESDvtxs->GetChi2toNDF();

  fVtx2->SetPosition(pos[0], pos[1], pos[2]);
  fVtx2->SetCovMatrix(cov);
  fVtx2->SetChi2perNDF(chi2perNDF);
  fVtx2->SetID(ESDvtxs->GetID());
  fVtx2->SetType(AliAODVertex::kUndef);

}

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::SetCustomTPCpid(Float_t *par, Float_t sigma) {
  if (par == 0x0 && sigma <= 0) {
    fCustomTPCpid.Set(1);
  } else {
    fCustomTPCpid.Set(6);
    for (int i = 0; i < 5; ++i)
      fCustomTPCpid.AddAt(par[i],i);
    fCustomTPCpid.AddAt(sigma, 5);
  }
}

//_________________________________________________________________________
Float_t AliAnalysisTaskHypCrossCheck::GetBetheAlephTPCnsigmas(AliESDtrack* t, AliPID::EParticleType part) {
    const Float_t massZ = AliPID::ParticleMassZ(part);
    const Float_t p = t->GetInnerParam()->GetP() / massZ;
    const Float_t pCharge = AliPID::ParticleCharge(part);
    const Float_t r = pCharge * pCharge * AliExternalTrackParam::BetheBlochAleph(p, fCustomTPCpid[0], fCustomTPCpid[1],
        fCustomTPCpid[2], fCustomTPCpid[3],
        fCustomTPCpid[4]);
    return (t->GetTPCsignal() - r) / (fCustomTPCpid[5] * r);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::HasTOF(AliESDtrack *trk, float &beta_tof){
    // TOF signal
    Bool_t TOFon = kFALSE;
    Bool_t isTOFout = trk->GetStatus() & AliESDtrack::kTOFout;
    Bool_t isTOFtime = trk->GetStatus() & AliESDtrack::kTIME;
    const float TOFlength = trk->GetIntegratedLength();
    Bool_t isTOFreached = Bool_t(isTOFout & isTOFtime) && TOFlength > 350.;
    if(!isTOFreached){
        TOFon = kFALSE;
    }
    else{
         const float time = trk->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(trk->P());

         if (time < 0 || time < TOFlength/2.99792457999999984e-02) {
             TOFon = kFALSE;
         } else {
             beta_tof = TOFlength / (2.99792457999999984e-02 * time);
             //const float gamma = 1/TMath::Sqrt(1-beta*beta);
             //const float mass = p/TMath::Sqrt(gamma*gamma - 1);
             TOFon = kTRUE;
         }
    }
    return TOFon;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::PassPIDSelection(AliESDtrack *trk, int specie, Bool_t isTOFin){
    //PID selection
    bool tofPID = kTRUE, tpcPID = kTRUE;
    //TPC-pid
    float const nsigmaTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::EParticleType (specie)));
    if(nsigmaTPC>fRequestTPCSigmas) tpcPID = kFALSE;
    else tpcPID = kTRUE;
    //TOF-pid
    if(fRequestTOFPid){
        if(isTOFin){
            float const nsigmaTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::EParticleType (specie)));
            if(nsigmaTOF>fRequestTOFSigmas) tofPID = kFALSE;
            else tofPID = kTRUE;
        }
    }
    return tpcPID && tofPID;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::CheckPrimaryDistribution(AliESDtrack *trk, Int_t &pdgParticle){

  Int_t lb = trk->GetLabel();
  TParticle *part = fStack->Particle(TMath::Abs(lb));
  if(part->GetUniqueID() != kPDecay) return kFALSE;
  pdgParticle = part->GetPdgCode();
  if(TMath::Abs(pdgParticle) != 211 && TMath::Abs(pdgParticle) != 1000020030) return kFALSE;
  TParticle *mom = fStack->Particle(TMath::Abs(part->GetFirstMother()));
  if(TMath::Abs(mom->GetPdgCode()) != 1010010030) return kFALSE;
  if(TMath::Abs(fStack->Particle(mom->GetLastDaughter()-1)->GetPdgCode()) == 2212) return kFALSE;
  if(TMath::Abs(pdgParticle) == 211 && (mom->GetLastDaughter() != lb) ) return kFALSE;
  if(TMath::Abs(pdgParticle) == 1000020030 && ((mom->GetLastDaughter()-1) != lb) ) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypCrossCheck::DaughtersFromSameMother(int lpion, int lhelium3, Double_t &pThyp){

    if(!fStack->IsSecondaryFromWeakDecay(TMath::Abs(lpion))) return kFALSE;
    if(!fStack->IsSecondaryFromWeakDecay(TMath::Abs(lhelium3))) return kFALSE;

    TParticle *pPi = fStack->Particle(TMath::Abs(lpion));
    TParticle *pHe3 = fStack->Particle(TMath::Abs(lhelium3));

    Int_t labelM_pio = pPi->GetFirstMother();
    Int_t labelM_he3 = pHe3->GetFirstMother();

    TParticle *tparticleMotherPi = fStack->Particle(TMath::Abs(labelM_pio));
    TParticle *tparticleMotherHe3 = fStack->Particle(TMath::Abs(labelM_he3));

    if(labelM_he3 != labelM_pio) return kFALSE;

    if(TMath::Abs(tparticleMotherHe3->GetPdgCode()) != 1010010030) return kFALSE;

    //if(pHe3->GetUniqueID() != kPDecay || pPi->GetUniqueID() != kPDecay) return kFALSE;


    if((pHe3->GetPdgCode() == 1000020030 && pPi->GetPdgCode() == -211) || (pHe3->GetPdgCode() == -1000020030 && pPi->GetPdgCode() == 211)){
      pThyp = tparticleMotherHe3->Pt();
      return kTRUE;
    }
    else return kFALSE;

}

//________________________________________________________________________
/*void AliAnalysisTaskHypCrossCheck::SetParametersAtPrimaryVertex(AliESDtrack *trk, AliESDtrack *extpar)const{

  const Double_t *par=extpar->GetParameter();
  const Double_t *cov=extpar->GetCovariance();
  Double_t alpha=extpar->GetAlpha();
  Double_t x=extpar->GetX();
  trk->Set(x,alpha,par,cov);
  return;
}*/

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::CheckGenerated(){

  int pdgCodePart = 0.;
  int pdgHe3, pdgPi = 0.;

  TParticle *candidateDau = 0x0;
  TParticle *dauHe3 = 0x0;
  TParticle *dauPi = 0x0;
  Bool_t isMotherPrimary = kFALSE;
  Bool_t isThreeBody, isTwoBody = kFALSE;

  for (Int_t iMC=0; iMC<fStack->GetNtrack(); iMC++){
      TParticle *p0 = fStack->Particle(iMC);

      isMotherPrimary = kFALSE;
      isThreeBody = kFALSE;
      isTwoBody = kFALSE;

      if(!p0) continue;

      pdgCodePart = p0->GetPdgCode();
      isMotherPrimary = fStack->IsPhysicalPrimary(iMC);


      if(pdgCodePart == 11) fHistParticle->Fill(0); // e-
      if(pdgCodePart == -11) fHistParticle->Fill(1); // e+
      if(pdgCodePart == 211) fHistParticle->Fill(2); // pi+
      if(pdgCodePart == -211) fHistParticle->Fill(3); // pi-
      if(pdgCodePart == 321) fHistParticle->Fill(4); // k+
      if(pdgCodePart == -321) fHistParticle->Fill(5); // k-
      if(pdgCodePart == 2212) fHistParticle->Fill(6); // p
      if(pdgCodePart == -2212) fHistParticle->Fill(7); // pbar
      if(pdgCodePart == 1000010020) fHistParticle->Fill(8); // d
      if(pdgCodePart == -1000010020) fHistParticle->Fill(9); // dbar
      if(pdgCodePart == 1000010030) fHistParticle->Fill(10); // t
      if(pdgCodePart == -1000010030) fHistParticle->Fill(11); // tbar
      if(pdgCodePart == 1000020030) fHistParticle->Fill(12); // He3
      if(pdgCodePart == -1000020030) fHistParticle->Fill(13); // He3bar
      if(pdgCodePart == 1000020040) fHistParticle->Fill(14); // He4
      if(pdgCodePart == -1000020040) fHistParticle->Fill(15); // He4bar
      if(pdgCodePart == 3122) fHistParticle->Fill(20); // Lambda
      if(pdgCodePart == 1020010020) fHistParticle->Fill(21); //Xi0-proton
      if(pdgCodePart == -1020010020) fHistParticle->Fill(22); //anti Xi0-proton

      if(TMath::Abs(pdgCodePart)==1000010020 || TMath::Abs(pdgCodePart)==1000010030 || TMath::Abs(pdgCodePart)==1000020030 || TMath::Abs(pdgCodePart)==1000020040
            || TMath::Abs(pdgCodePart)==1010010030 || TMath::Abs(pdgCodePart)==1010010040 || TMath::Abs(pdgCodePart)==1010020040) fHistParticle_Mass->Fill(p0->GetMass());

      if(pdgCodePart == 1010010030 || pdgCodePart == -1010010030){ // is Hypertriton?
          fHistHypertritonMomGen->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
          fHistHypertritonYGen->Fill(p0->Y());
          candidateDau = fStack->Particle(p0->GetLastDaughter()-1);

          if(TMath::Abs(candidateDau->GetPdgCode()) == 2212 && candidateDau->GetUniqueID() == kPDecay) {
              isThreeBody = kTRUE;
              if(pdgCodePart == 1010010030){
                fHistParticle->Fill(16); //3LH in three body
                fHistHypertritonMomGen_3Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                fHistHypertritonYGen_3Body->Fill(p0->Y());
              }
              if(pdgCodePart == -1010010030) {
                fHistParticle->Fill(17); //anti-3LH in three body
                fHistAntiHypertritonMomGen_3Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                fHistAntiHypertritonYGen_3Body->Fill(p0->Y());
              }
              if(isMotherPrimary) {
                  fHistHypertritonMomGen_isPrimary_3Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                  fHistHypertritonYGen_isPrimary_3Body->Fill(p0->Y());
              }
          } // end of isThreeBody

          if(TMath::Abs(candidateDau->GetPdgCode()) == 1000020030 && candidateDau->GetUniqueID() == kPDecay) {
              isTwoBody = kTRUE;
              if(pdgCodePart == 1010010030) {
                fHistParticle->Fill(18); //3LH in two body
                fHistHypertritonMomGen_2Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                fHistHypertritonYGen_2Body->Fill(p0->Y());
              }
              if(pdgCodePart == -1010010030) {
                fHistParticle->Fill(19); //anti-3LH in two body
                fHistAntiHypertritonMomGen_2Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                fHistAntiHypertritonYGen_2Body->Fill(p0->Y());
              }
              if(isMotherPrimary) {
                  fHistHypertritonMomGen_isPrimary_2Body->Fill(TMath::Sqrt((p0->Px()*p0->Px()) + (p0->Py()*p0->Py())));
                  fHistHypertritonYGen_isPrimary_2Body->Fill(p0->Y());
              }
          } // end of isTwoBody

          if(isTwoBody){
              dauHe3 = fStack->Particle(p0->GetLastDaughter()-1);
              pdgHe3 = dauHe3->GetPdgCode();
              if(TMath::Abs(pdgHe3) != 1000020030 || dauHe3->GetUniqueID() != kPDecay) {
                  printf("FATAL ERROR: helium-3 not correctly identified!\n");
              }
              dauPi = fStack->Particle(p0->GetLastDaughter());
              pdgPi = dauPi->GetPdgCode();
              if(TMath::Abs(pdgPi) != 211 || dauPi->GetUniqueID() != kPDecay) {
                  printf("FATAL ERROR: pion not correctly identified!\n");
              }

              //if(fChooseMatter && pdgCodePart == -1010010030) continue;
              //if(fChooseAntiMatter && pdgCodePart == 1010010030) continue;

              if(fFillTGen){ //filling Generated TTree
                  fTCentrality_gen = fCentralityPercentile;
                  fTHyper_px = p0->Px();
                  fTHyper_py = p0->Py();
                  fTHyper_pz = p0->Pz();
                  fTHyper_eta = p0->Eta();
                  fTHyper_rapidity = p0->Y();
                  fTHyper_pv_x = p0->Vx();
                  fTHyper_pv_y = p0->Vy();
                  fTHyper_pv_z = p0->Vz();
                  if(pdgCodePart == 1010010030) fTHyper_charge = 1;
                  else if(pdgCodePart == -1010010030) fTHyper_charge = -1;
                  fTHe3_px = dauHe3->Px();
                  fTHe3_py = dauHe3->Py();
                  fTHe3_pz = dauHe3->Pz();
                  fTHe3_eta = dauHe3->Eta();
                  fTHe3_rapidity = dauHe3->Y();
                  fTHe3_pv_x = dauHe3->Vx();
                  fTHe3_pv_y = dauHe3->Vy();
                  fTHe3_pv_z = dauHe3->Vz();
                  if(pdgHe3 == 1000020030) fTHe3_charge = 1;
                  else if(pdgHe3 == -1000020030) fTHe3_charge = -1;
                  fTPion_px = dauPi->Px();
                  fTPion_py = dauPi->Py();
                  fTPion_pz = dauPi->Pz();
                  fTPion_eta = dauPi->Eta();
                  fTPion_rapidity = dauPi->Y();
                  if(pdgPi == 211) fTPion_charge = 1;
                  else if(pdgPi == -211) fTPion_charge = -1;

                  fTGen->Fill();
                  PostData(3,fTGen);
              } // end of fill fTGen
          } // end of isTwoBody info search
      } // end of if for Hypertriton ID
  } // end of loop on fStack MC track
}

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::CombineTwoTracks(Bool_t isMatter, TArrayI arrHe3, TArrayI arrPi, Int_t charge){
  //Method to combine tracks
  //Here implemented topological and kinematical cuts
  //Total charge predefined with the correct type of array assigned to this method

  //Define variables
  //===Combining tracks loop===//
  Double_t bz = fESDevent->GetMagneticField();
  fVertexer->SetFieldkG(bz);

  Double_t dlh[3] = {0,0,0}; //array for the coordinates of the decay length
  Double_t dca_he3pi, angle_he3pi = 0.;
  Float_t dcapvhe[2] = {0.,0.};
  Float_t dcapvpi[2] = {0.,0.};
  Double_t dcahe3[2] = {0.,0.}; // dca between the candidate d,p,pi and the candidate decay vertex
  Double_t dcapi[2] = {0.,0.}; // dcad[0]= transverse plane coordinate; dcad[1]= z coordinate
  Double_t decayLengthH3L, normalizedDecayL, rapidity, pointingAngleH, ctau= 0.;
  Double_t lPi, lHe3 = 0;
  Double_t decVt[3] = {0.,0.,0.};
  Bool_t brotherHood = kFALSE;
  Int_t labelM_he3, labelM_pio = 0;
  TLorentzVector posHe3, negPi; //Lorentz vector of helium-3 and pion in the LAB
  AliESDtrack *trackHe3 = 0x0;
  AliESDtrack *trackNPi = 0x0;

  Double_t xthiss(0.0);
  Double_t xpp(0.0);

  AliESDVertex *decayVtx = 0x0;

  TLorentzVector Hypertriton;
  TVector3 h1, he1, pi1;
  Double_t pTotHyper = 0.;
  Double_t TransvMomGen = 0.;

  TParticle *tparticleHe3 = 0x0;
  TParticle *tparticlePi = 0x0;

  // -------------------------------------------------------
  // Loop for Invariant Mass
  // -------------------------------------------------------

    for(Int_t m=0; m<arrHe3.GetSize(); m++){ // candidate helium-3 loop

      trackHe3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(arrHe3[m]));

      for(Int_t s=0; s<arrPi.GetSize(); s++ ){ // candidate pion loop cpion.size()

          fTrkArray->Clear();
          Hypertriton.Clear();
          posHe3.Clear();
          negPi.Clear();
          h1.Clear();
          he1.Clear();
          pi1.Clear();

          trackNPi = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(arrPi[s]));
          brotherHood = kFALSE;
          TransvMomGen = 0.;

          if(fMC){
            lHe3 = trackHe3->GetLabel();
            tparticleHe3 = fStack->Particle(TMath::Abs(lHe3));
            lPi = trackNPi->GetLabel();
            tparticlePi = fStack->Particle(TMath::Abs(lPi));
            brotherHood = DaughtersFromSameMother(lPi, lHe3,TransvMomGen);
          }

          //if(brotherHood) continue;

          if(trackNPi->GetID() == trackHe3->GetID()) continue;

          trackHe3->GetImpactParameters(dcapvhe[0],dcapvhe[1]);
          trackNPi->GetImpactParameters(dcapvpi[0],dcapvpi[1]);

          dca_he3pi = trackNPi->GetDCA(trackHe3,bz,xthiss,xpp);
          fHistDCAhe3pion->Fill(dca_he3pi);
          if(dca_he3pi > fDCAhe3pi) continue;

          if(fMC){
            fHistDeltaPt_Pion->Fill(trackNPi->Pt(),tparticlePi->Pt()-trackNPi->Pt());
            fHistDeltaPt_He3->Fill(charge*trackHe3->Pt(),tparticleHe3->Pt()-charge*trackHe3->Pt());
          }
          fTrkArray->AddAt(trackHe3,0);
          fTrkArray->AddAt(trackNPi,1);


          fVertexer->SetVtxStart(fPrimaryVertex);
          decayVtx = (AliESDVertex*)fVertexer->VertexForSelectedESDTracks(fTrkArray);

          SetConvertedAODVertices(fPrimaryVertex,decayVtx);

          fHistZDecayVtx->Fill(decayVtx->GetZ());
          fHistXDecayVtx->Fill(decayVtx->GetX());
          fHistYDecayVtx->Fill(decayVtx->GetY());
          decVt[0] = decayVtx->GetX();
          decVt[1] = decayVtx->GetY();
          decVt[2] = decayVtx->GetZ();

          dlh[0]=fESDevent->GetPrimaryVertex()->GetX() - decayVtx->GetX();
          dlh[1]=fESDevent->GetPrimaryVertex()->GetY() - decayVtx->GetY();
          dlh[2]=fESDevent->GetPrimaryVertex()->GetZ() - decayVtx->GetZ();

          decayLengthH3L = TMath::Sqrt((dlh[0]*dlh[0]) + (dlh[1]*dlh[1]) + (dlh[2]*dlh[2]));
          normalizedDecayL = fVtx2->DistanceToVertex(fVtx1)/fVtx2->ErrorDistanceToVertex(fVtx1);

          fHistDecayLengthH3L->Fill(decayLengthH3L);
          fHistNormalizedDecayL->Fill(normalizedDecayL);

          trackHe3->PropagateToDCA(decayVtx, bz, 10,dcahe3);
          fHistDCAXYhe3vtx->Fill(dcahe3[0]);
          fHistDCAZhe3vtx->Fill(dcahe3[1]);
          fHistDCAhe3vtx->Fill(TMath::Sqrt( (dcahe3[0]*dcahe3[0])+(dcahe3[1]*dcahe3[1]) ));

          trackNPi->PropagateToDCA(decayVtx, bz, 10,dcapi);
          fHistDCAXYpionvtx->Fill(dcapi[0]);
          fHistDCAZpionvtx->Fill(dcapi[1]);
          fHistDCApionvtx->Fill(TMath::Sqrt( (dcapi[0]*dcapi[0])+(dcapi[1]*dcapi[1]) ));

          //trackHe3->PropagateToDCA(fPrimaryVertex,bz,100);
          //trackNPi->PropagateToDCA(fPrimaryVertex,bz,100);

          delete decayVtx;

          //if(decayLengthH3L > fMaxDecayLength || decayLengthH3L < fMinDecayLength) continue;
          //if(TMath::Sqrt((dcahe3[0]*dcahe3[0])+(dcahe3[1]*dcahe3[1])) > fDCAHe3SVmax) continue;
          //if(TMath::Abs(dcapi[0]) > fDCAPiSVxymax) continue;
          //if(TMath::Abs(dcapi[1]) > fDCAPiSVzmax) continue;

          //if(normalizedDecayL < fMinNormalizedDecL) continue;

          posHe3.SetXYZM(charge*trackHe3->Px(),charge*trackHe3->Py(),charge*trackHe3->Pz(),2.80923);
          negPi.SetXYZM(trackNPi->Px(),trackNPi->Py(),trackNPi->Pz(),0.13957);

          Hypertriton=posHe3+negPi;

          if(brotherHood){
            fHistDeltaPt_Hyper->Fill(Hypertriton.Pt(),TransvMomGen-Hypertriton.Pt());
          }
          pTotHyper = Hypertriton.P();
          if(Hypertriton.Pt() < fMinPtMother || Hypertriton.Pt() > fMaxPtMother) continue;

          if(fSideBand == kTRUE && (Hypertriton.M() < 3.08 || Hypertriton.M() > 3.18)) continue;
          ctau = (Hypertriton.M()*decayLengthH3L)/pTotHyper;
          fHistLifetime->Fill(ctau);

          if(ctau < fMinLifeTime || ctau > fMaxLifeTime) continue;

          //Angular correlation

          he1.SetXYZ(charge *trackHe3->Px(),charge*trackHe3->Py(),charge*trackHe3->Pz());
          pi1.SetXYZ(trackNPi->Px(),trackNPi->Py(),trackNPi->Pz());

          angle_he3pi = he1.Angle(pi1);
          fHistAngle_He3_pion->Fill(angle_he3pi);

          //if(angle_he3pi > fAnglehe3pi) continue;


          fHistPtHypertriton->Fill(Hypertriton.Pt());
          fHistPtHelium3->Fill(trackHe3->Pt());
          fHistPtPion->Fill(trackNPi->Pt());

          rapidity = Hypertriton.Rapidity();
          fHistHyperRapidity->Fill(rapidity);

          if(TMath::Abs(rapidity) > fRapidity) continue;


          h1.SetXYZ(-dlh[0],-dlh[1],-dlh[2]);
          pointingAngleH = Hypertriton.Angle(h1);
          fHistCosPointingAngle->Fill(TMath::Cos(pointingAngleH));

          if(fMC){
            if(brotherHood){
            fHistDCAhe3pionMCt->Fill(dca_he3pi);
            fHistDeltaPt_PionMCt->Fill(trackNPi->Pt(),tparticlePi->Pt()-trackNPi->Pt());
            fHistDeltaPt_He3MCt->Fill(trackHe3->Pt(),tparticleHe3->Pt()-trackHe3->Pt());
            fHistDCAXYhe3vtxMCt->Fill(dcahe3[0]);
            fHistDCAZhe3vtxMCt->Fill(dcahe3[1]);
            fHistDCAhe3vtxMCt->Fill(TMath::Sqrt( (dcahe3[0]*dcahe3[0])+(dcahe3[1]*dcahe3[1]) ));
            fHistDCAXYpionvtxMCt->Fill(dcapi[0]);
            fHistDCAZpionvtxMCt->Fill(dcapi[1]);
            fHistDCApionvtxMCt->Fill(TMath::Sqrt( (dcapi[0]*dcapi[0])+(dcapi[1]*dcapi[1]) ));
            fHistZDecayVtxMCt->Fill(decVt[2]);
            fHistXDecayVtxMCt->Fill(decVt[0]);
            fHistYDecayVtxMCt->Fill(decVt[1]);
            fHistDecayLengthH3L_MCt->Fill(decayLengthH3L);
            fHistNormalizedDecayL_MCt->Fill(normalizedDecayL);
            fHistDeltaPt_HyperMCt->Fill(Hypertriton.Pt(),TransvMomGen-Hypertriton.Pt());
            fHistLifetime_MCt->Fill(ctau);
            fHistAngle_He3_pion_MCt->Fill(angle_he3pi);
            fHistHypertritonMomMCt->Fill(pTotHyper);
            fHistPtPionMCt->Fill(trackNPi->Pt());
            fHistPtHelium3MCt->Fill(trackHe3->Pt());
            fHistPtHypertritonMCt->Fill(Hypertriton.Pt());
            fHistHyperRapidityMCt->Fill(rapidity);
            fHistCosPointingAngleMCt->Fill(TMath::Cos(pointingAngleH));
            fHistMassVsPtMCt->Fill(Hypertriton.Pt(),Hypertriton.M());
            if(isMatter) {
              fHistMassHypertritonMCt->Fill(Hypertriton.M());
              fHistMassHyperVsPtMCt->Fill(Hypertriton.Pt(),Hypertriton.M());
            }
            if(!isMatter){
              fHistMassAntiHypertritonMCt->Fill(Hypertriton.M());
              fHistMassAntiHyperVsPtMCt->Fill(Hypertriton.Pt(),Hypertriton.M());
            }
          }
        } //end of if(fMC)

          if (TMath::Cos(pointingAngleH) < fCosPointingAngle) continue;

          fHistMassVsPt->Fill(Hypertriton.Pt(),Hypertriton.M());

          if(isMatter){
            fHistMassHypertriton->Fill(Hypertriton.M());
            fHistMassHyperVsPt->Fill(Hypertriton.Pt(),Hypertriton.M());
          }
          if(!isMatter){
            fHistMassAntiHypertriton->Fill(Hypertriton.M());
            fHistMassAntiHyperVsPt->Fill(Hypertriton.Pt(),Hypertriton.M());
          }



	  if(fFillTree){
	    fTMCtruth = brotherHood;
	    fTCentralityPerc = fCentralityPercentile;
	    //helium-3
	    fTchi2NDFhe3 = trackHe3->GetTPCchi2()/trackHe3->GetTPCclusters(0);
	    fTMhypoTrkhe3 = trackHe3->GetPIDForTracking();
	    fTPCclsPIDhe3 = trackHe3->GetTPCsignalN();
      if(trackHe3->GetStatus()&AliVTrack::kITSrefit) fITSrefithe3 = kTRUE;
      else fITSrefithe3 = kFALSE;
	    fTpTPChe3 = trackHe3->GetInnerParam()->GetP();
	    fTpXhe3 = trackHe3->Px();
	    fTpYhe3 = trackHe3->Py();
	    fTpZhe3 = trackHe3->Pz();
	    fTTPCnsigmahe3 = fPIDResponse->NumberOfSigmasTPC(trackHe3,AliPID::kHe3);
      fTDCAXYhe3prvtx = dcapvhe[0];
      fTDCAZhe3prvtx = dcapvhe[1];
	    //pion
	    fTchi2NDFpion = trackNPi->GetTPCchi2()/trackNPi->GetTPCclusters(0);
	    fTMhypoTrkpion = trackNPi->GetPIDForTracking();
	    fTPCclsPIDpion = trackNPi->GetTPCsignalN();
      if(trackNPi->GetStatus()&AliVTrack::kITSrefit) fITSrefitpion = kTRUE;
      fITSrefitpion = kFALSE;
	    fTpTPCpion = trackNPi->GetInnerParam()->GetP();
	    fTpXpion = trackNPi->Px();
	    fTpYpion = trackNPi->Py();
	    fTpZpion = trackNPi->Pz();
	    fTTPCnsigmapion = fPIDResponse->NumberOfSigmasTPC(trackNPi,AliPID::kPion);
      fTDCAXYpioprvtx = dcapvpi[0];
      fTDCAZpioprvtx = dcapvpi[1];
	    //triplets
	    fTDCAhe3pi = dca_he3pi;
	    fTDCAXYhevtx = dcahe3[0];
	    fTDCAZhevtx = dcahe3[1];
	    fTDCAXYpivtx = dcapi[0];
	    fTDCAZpivtx = dcapi[1];
	    fTAngle_he3pi = he1.Angle(pi1);
	    fTRapidity = Hypertriton.Rapidity();
	    fTDecayLength = decayLengthH3L;
	    fTDecayLengthError = fVtx2->ErrorDistanceToVertex(fVtx1);
	    fTCosPA = TMath::Cos(pointingAngleH);

	    if(fMC){
	      labelM_he3 = tparticleHe3->GetFirstMother();
	      labelM_pio = tparticlePi->GetFirstMother();
	      TParticle *tpmomHe = fStack->Particle(TMath::Abs(labelM_he3));
	      TParticle *tpmomPi = fStack->Particle(TMath::Abs(labelM_pio));
	      fTphe3_gen_X = tparticleHe3->Px();
	      fTphe3_gen_Y = tparticleHe3->Py();
	      fTphe3_gen_Z = tparticleHe3->Pz();
	      fTppio_gen_X = tparticlePi->Px();
	      fTppio_gen_Y = tparticlePi->Py();
	      fTppio_gen_Z = tparticlePi->Pz();
	      fTpdgHe3 = tparticleHe3->GetPdgCode();
	      fTpdgPion = tparticlePi->GetPdgCode();
	      fTmomidHe3 = tparticleHe3->GetFirstMother();
	      fTmomidPi = tparticlePi->GetFirstMother();
	      fTpdgmomHe3 = tpmomHe->GetPdgCode();
	      fTpdgmomPi = tpmomPi->GetPdgCode();
	      fTuniqID_he3 = tparticleHe3->GetUniqueID();
	      fTuniqID_pion = tparticlePi->GetUniqueID();
	    }

	    if(!isMatter) fTInvariantMass = -Hypertriton.M();
	    else fTInvariantMass = Hypertriton.M();

	    fTTree->Fill();
	    PostData(2,fTTree);
	  } //end of Fill Tree
      } // end of candidate pion loop
    } // end of candidate helium-3 loop


/*
  Method to set and restore track parameters at Primary Vertex
  - void SetParametersAtPrimaryVertex(...)
  - TClonesArray *tracksAtPrimaryVertex = new TClonesArray("AliESDtrack",3);
  - tracksAtPrimaryVertex->RemoveAt(0); // 0 for deuteron, 1 for helium-3, 2 for pion
  - new ((*tracksAtPrimaryVertex)[1]) AliESDtrack(*trackHe3);
    new ((*tracksAtPrimaryVertex)[2]) AliESDtrack(*trackNPi);
    to add a copy of a track, not the direct track we are using
    possible variation: instead of adding a copy of a track to restore parameter,
                        propagate the copy of the track to the SV and then delete the copyright
                        so, no need to ri-propagate to primary vertex
  - SetParametersAtPrimaryVertex(trackHe3, (AliESDtrack*)tracksAtPrimaryVertex->At(1));
    SetParametersAtPrimaryVertex(trackNPi, (AliESDtrack*)tracksAtPrimaryVertex->At(2));
  - delete tracksAtPrimaryVertex;
*/

  return;
}

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::UserCreateOutputObjects(){
  // Create a TList with Histograms
  // Called once
  //printf("**************************************************\n");
  //printf("AliAnalysisTaskHypCrossCheck::CreateOutputObjects()\n");
  //printf("**************************************************\n");

  fVertexer = new AliVertexerTracks();
  fTrkArray = new TObjArray(2);
  fVtx1 = new AliAODVertex();
  fVtx2 = new AliAODVertex();

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("clistHypertriton");

  fHistCount = new TH1F("fHistCount","Counter histogram",4,-0.5,3.5);
  fHistCount->GetXaxis()->SetBinLabel(1,"Reco Event");
  fHistCount->GetXaxis()->SetBinLabel(2,"Trigger sel");
  fHistCount->GetXaxis()->SetBinLabel(3,"Centrality sel");
  fHistCount->GetXaxis()->SetBinLabel(4,"Primary vtx sel");

  fHistCentralityClass = new TH1F("fHistCentralityClass","Centrality Class; centrality class; entries",11,-0.5,10.5);
  fHistCentralityPercentile = new TH1F("fHistCentralityPercentile","Centrality; centrality percentile; entries",101,-0.5,100.5);

  fHistTrigger = new TH1F("fHistTrigger","Trigger statistics",4,-0.5,3.5);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"kMB");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");

  fHistMultiplicity = new TH1F("fHistMultiplicity ", "Multiplicity; multiplicity; entries", 100, 0, 20000);

  fHistZPrimaryVtx = new TH1F("fHistZPrimaryVtx","primary vertex - Z coordinate; z_{primary vtx} (cm); entries",40,-10.,10.);

  fHistXPrimaryVtx = new TH1F("fHistXPrimaryVtx","primary vertex - X coordinate; x_{primary vtx} (cm); entries",800,-1.,1.);

  fHistYPrimaryVtx = new TH1F("fHistYPrimaryVtx","primary vertex - Y coordinate; y_{primary vtx} (cm); entries",800,-1.,1.);

  fHistChi2perTPCcluster = new TH1F("fHistChi2perTPCcluster","#Chi^{2}/TPC cluster distribution; #Chi^{2}/TPCcls; entries",44,-0.5,10.5);

  fHistTrackFlagReco = new TH1F("fHistTrackFlagReco","Check track flag",6,-0.5,5.5);
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(1,"kTPCrefit");
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(2,"kITSin");
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(3,"kITSout");
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(4,"kITSrefit");
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(5,"kITSin & !kITSrefit");
  fHistTrackFlagReco->GetXaxis()->SetBinLabel(6,"kITSout & !kITSrefit");


  //TPC
  fHistTPCpid = new TH2F("fHistTPCpid", "dE/dx for all tracks; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);

  //Hypertriton prongs
  fHistTPCHe3signal = new TH2F("fHistTPCHe3signal", "dE/dx after ^{3}He PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistTPCpionsignal= new TH2F("fHistTPCpionsignal", "dE/dx after  #pi^{-} PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistNsigmaHe3 = new TH2F("fHistNsigmaHe3","n#sigma candidate ^{3}He; p_{TPC}(GeV/c); n#sigma",100,0,10,160,-8,8);
  fHistNsigmaPion = new TH2F("fHistNsigmaPion","n#sigma candidate #pi; p_{TPC}(GeV/c); n#sigma",100,0,10,160,-8,8);

  //TOF

  fHistTOFsignal = new TH2F("fHistTOFsignal","TOF signal; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFHe3signal = new TH2F("fHistTOFHe3signal","#beta vs TPCmom - ^{3}He; p_{TPC} (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistpionTPCcls = new TH1F("fHistpionTPCcls","#pi^{-} TPC clusters; TPC clusters; entries",201,-0.5,200.5);
  //fHistCorrDCAHe3primary = new TH2F("fHistCorrDCAHe3primary","DCA_{PV,xy} vs DCA_{PV,z} - helium-3; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
  //fHistCorrDCApiprimary = new TH2F("fHistCorrDCApiprimary","DCA_{PV,xy} vs DCA_{PV,z} - pion; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
  fHistDCApiprimary = new TH1F("fHistDCApiprimary","DCA pion-primary vertex; DCA (cm); entries",800,0.f,20.f);
  fHistDCAXYpiprimary = new TH1F("fHistDCAXYpiprimary","DCAxy pion-primary vertex; DCA_{xy} (cm); entries",800,-10.f,10.f);
  fHistDCAZpiprimary = new TH1F("fHistDCAZpiprimary","DCAz pion-primary vertex; DCA_{z} (cm); entries",800,-10.f,10.f);

  fHistDCAHe3primary = new TH1F("fHistDCAHe3primary","DCA ^{3}He-primary vertex; DCA (cm); entries",800,0.f,20.f);
  fHistDCAXYHe3primary = new TH1F("fHistDCAXYHe3primary","DCAxy ^{3}He-primary vertex; DCA_{xy} (cm); entries",800,-10.f,10.f);
  fHistDCAZHe3primary = new TH1F("fHistDCAZHe3primary","DCAz ^{3}He-primary vertex; DCA_{z} (cm); entries",800,-10.f,10.f);

  //DCA prongs
  fHistDCAhe3pion = new TH1F("fHistDCAhe3pion","DCA ^{3}He-#pi tracks;DCA (cm);entries",550,-0.5,5.0);

  //Decay histo
  fHistZDecayVtx = new TH1F("fHistZDecayVtx","decay vertex - z coordinate; z_{decay vtx} (cm); entries",800,-20.f,20.f);
  fHistXDecayVtx = new TH1F("fHistXDecayVtx","decay vertex - x coordinate; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistYDecayVtx = new TH1F("fHistYDecayVtx","decay vertex - y coordinate; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistDCApionvtx = new TH1F("fHistDCApionvtx","DCA candidate #pi^{-}-decay vertex; DCA (cm); entries",1500,0.,15.);
  fHistDCAXYpionvtx = new TH1F("fHistDCAXYpionvtx","DCA candidate #pi^{-}-decay vertex - xy coordinate; DCA_{xy} (cm); entries",1000,-5.,5.);
  fHistDCAZpionvtx = new TH1F("fHistDCAZpionvtx","DCA candidate #pi^{-}-decay vertex - z coordinate; DCA_{z} (cm); entries",2000,-10.,10.);
  fHistDCAhe3vtx = new TH1F("fHistDCAhe3vtx","DCA candidate ^{3}He-decay vertex; DCA (cm); entries",1000,0.,15.);
  fHistDCAXYhe3vtx = new TH1F("fHistDCAXYhe3vtx","DCA candidate ^{3}He-decay vertex - xy coordinate; DCA_{xy} (cm); entries",1000,-5.,5.);
  fHistDCAZhe3vtx = new TH1F("fHistDCAZhe3vtx","DCA candidate ^{3}He-decay vertex - z coordinate; DCA_{z} (cm); entries",2000,-10.,10.);
  fHistDecayLengthH3L = new TH1F("fHistDecayLengthH3L","decay length ^{3}H_{#Lambda}; decay length (cm); entries",400,0.,400.);
  fHistNormalizedDecayL = new TH1F("fHistNormalizedDecayL","normalized decay length; decL/#sigma_{dL}; entries",400,0.,100.);
  fHistDeltaPt_Hyper = new TH2F("fHistDeltaPt_Hyper","delta p_{T} ^{3}H_{#Lambda};p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
  fHistDeltaPt_Pion = new TH2F("fHistDeltaPt_Pion","delta p_{T} #pi; p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
  fHistDeltaPt_He3 = new TH2F("fHistDeltaPt_He3","delta p_{T} ^{3}He; p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
  fHistLifetime = new TH1F("fHistLifetime","ct ^{3}H_{#Lambda}; ct(cm); entries",400,0.,40.);

  fHistAngle_He3_pion = new TH1F("fHistAngle_He3_pion","Angle between ^{3}He and #pi; #alpha_{He3_#pi} (rad); entries/(0.03 rad)",100,0.,TMath::Pi());

  fHistHyperRapidity = new TH1F("fHistHyperRapidity","rapidity distribution of ^{3}H_{#Lambda}; y; entries",400,-2.f,2.f);

  fHistCosPointingAngle= new TH1F("fHistCosPointingAngle", "Cosine of pointing angle distribution; cos(#theta_{pointing}); entries", 220, -1.1, 1.1);

  fHistPtPion = new TH1F("fHistPtPion","#pi p_{T}; p_{T} (GeV/c); entries",1000.,0.,10.);
  fHistPtHelium3 = new TH1F("fHistPtHelium3","^{3}He p_{T}; p_{T} (GeV/c); entries",1000.,0.,10.);
  fHistPtHypertriton = new TH1F("fHistPtHypertriton","candidate Hypertriton - p_{T} distribution; p_{T} (GeV/c); entries",130.,0.,13.);



  if(fMinvSignal){
      fHistMassHypertriton = new TH1F("fHistMassHypertriton", "Invariant mass distribution ^{3}He+#pi^{-};invariant mass ^{3}He+#pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton = new TH1F("fHistMassAntiHypertriton", "Invariant mass distribution ^{3}#bar{He}+#pi^{+};invariant mass #bar{d} + #bar{p} + #pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassVsPt = new TH2F("fHistMassVsPt","^{3}H_{#Lambda} and #bar{^{3}H_{#Lambda}} mass vs p_{T} ; p_{T}(GeV/c); ^{3}He+#pi (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
      fHistMassHyperVsPt = new TH2F("fHistMassHyperVsPt","^{3}H_{#Lambda} invariant mass vs p_{T} ; p_{T}(GeV/c); ^{3}He+#pi^{-} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
      fHistMassAntiHyperVsPt = new TH2F("fHistMassAntiHyperVsPt","#bar{^{3}H_{#Lambda}} invariant mass vs p_{T} ; p_{T}(GeV/c); ^{3}#bar{He}+#pi^{+} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
    }
  if(fMinvLikeSign){
      fHistMassHypertriton = new TH1F("fHistMassHypertriton_LS", "Invariant mass distribution - Like Sign;invariant mass ^{3}He+#pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton = new TH1F("fHistMassAntiHypertriton_LS", "Invariant mass distribution - Like Sign;invariant mass ^{3}#bar{He}+#pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassVsPt = new TH2F("fHistMassVsPt","^{3}H_{#Lambda} and #bar{^{3}H_{#Lambda}} mass vs p_{T} ; p_{T}(GeV/c); ^{3}He+#pi (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
      fHistMassHyperVsPt = new TH2F("fHistMassHyperVsPt","^{3}H_{#Lambda} invariant mass vs p_{T} ; p_{T}(GeV/c); ^{3}He+#pi^{-} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
      fHistMassAntiHyperVsPt = new TH2F("fHistMassAntiHyperVsPt","#bar{^{3}H_{#Lambda}} invariant mass vs p_{T} ; p_{T}(GeV/c); ^{3}#bar{He}+#pi^{-} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
    }


  if(fMC){

    fHistTPCdeusignal_pdg = new TH2F("fHistTPCdeusignal_pdg", "dE/dx after d PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
    fHistTPCtrisignal_pdg = new TH2F("fHistTPCtrisignal_pdg", "dE/dx after t PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
    fHistTPCHe3signal_pdg = new TH2F("fHistTPCHe3signal_pdg", "dE/dx after ^{3}He PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
    fHistTPCHe4signal_pdg = new TH2F("fHistTPCHe4signal_pdg", "dE/dx after ^{4}He PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);

    fHistTPCHe3signal_3LH = new TH2F("fHistTPCHe3signal_3LH", "dE/dx after ^{3}He from ^{3}H_{#Lambda} PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);

    fHistTPCpionsignal_3LH= new TH2F("fHistTPCpionsignal_3LH", "dE/dx after  #pi^{-} from ^{3}H_{#Lambda} PID; p/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);

    fHistNsigmaHe3_3LH = new TH2F("fHistNsigmaHe3_3LH","n#sigma candidate p from ^{3}H_{#Lambda}; p_{TPC}(GeV/c); n#sigma",100,0,10,160,-8,8);
    fHistNsigmaPion_3LH = new TH2F("fHistNsigmaPion_3LH","n#sigma candidate #pi from ^{3}H_{#Lambda}; p_{TPC}(GeV/c); n#sigma",100,0,10,160,-8,8);

    fHistParticle = new TH1F("fHistParticle","Check particle candidate",23,-0.5,22.5);
    fHistParticle->GetXaxis()->SetBinLabel(1,"electron");
    fHistParticle->GetXaxis()->SetBinLabel(2,"positron");
    fHistParticle->GetXaxis()->SetBinLabel(3,"#pi^{+}");
    fHistParticle->GetXaxis()->SetBinLabel(4,"#pi^{-}");
    fHistParticle->GetXaxis()->SetBinLabel(5,"K^{+}");
    fHistParticle->GetXaxis()->SetBinLabel(6,"K^{-}");
    fHistParticle->GetXaxis()->SetBinLabel(7,"proton");
    fHistParticle->GetXaxis()->SetBinLabel(8,"anti-proton");
    fHistParticle->GetXaxis()->SetBinLabel(9,"deuteron");
    fHistParticle->GetXaxis()->SetBinLabel(10,"anti-deuteron");
    fHistParticle->GetXaxis()->SetBinLabel(11,"triton");
    fHistParticle->GetXaxis()->SetBinLabel(12,"anti-triton");
    fHistParticle->GetXaxis()->SetBinLabel(13,"He3");
    fHistParticle->GetXaxis()->SetBinLabel(14,"anti-He3");
    fHistParticle->GetXaxis()->SetBinLabel(15,"#alpha");
    fHistParticle->GetXaxis()->SetBinLabel(16,"anti-#alpha");
    fHistParticle->GetXaxis()->SetBinLabel(17,"^{3}H_{#Lambda}");
    fHistParticle->GetXaxis()->SetBinLabel(18,"anti-^{3}H_{#Lambda}");
    fHistParticle->GetXaxis()->SetBinLabel(19,"^{3}H_{#Lambda} - 2 body");
    fHistParticle->GetXaxis()->SetBinLabel(20,"anti-^{3}H_{#Lambda} - 2 body");
    fHistParticle->GetXaxis()->SetBinLabel(21,"#Lambda - 2 body");
    fHistParticle->GetXaxis()->SetBinLabel(22,"Xi0-proton");
    fHistParticle->GetXaxis()->SetBinLabel(23,"Anti-Xi0-proton");

    fHistParticle_Mass = new TH1F("fHistParticle_Mass","Mass of generated particles; mass (GeV/c^{2})",4000,1.,5.);

    fHistpionTPCclsMCt = new TH1F("fHistpionTPCclsMCt","#pi^{-} TPC clusters - MCtruth; TPC clusters; entries",201,-0.5,200.5);
    fHisthelium3TPCclsMCt = new TH1F("fHisthelium3TPCclsMCt","^{3}He TPC clusters - MCtruth; TPC clusters; entries",201,-0.5,200.5);
    fHistpTpionMCt = new TH1F("fHistpTpionMCt","pion p_{T} distribution; p_{T} (GeV/c);entries",1000,0.,10.);
    fHistpThe3MCt = new TH1F("fHistpThe3MCt","^{3}He p_{T} distribution; p_{T} (GeV/c);entries",1000,0.,10.);
    fHistMompionMCt = new TH1F("fHistMompionMCt","pion p distribution; p (GeV/c);entries",1000,0.,10.);
    fHistMomHe3MCt = new TH1F("fHistMomHe3MCt","^{3}He p distribution; p (GeV/c);entries",1000,0.,10.);
    fHistCorrDCAHe3primaryMCt = new TH2F("fHistCorrDCAHe3primaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - ^{3}He MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
    fHistCorrDCApiprimaryMCt = new TH2F("fHistCorrDCApiprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - pion MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
    fHistDCApiprimaryMCt = new TH1F("fHistDCApiprimaryMCt","DCA pion-primary vertex; DCA (cm); entries",800,0.f,20.f);
    fHistDCAXYpiprimaryMCt = new TH1F("fHistDCAXYpiprimaryMCt","DCAxy pion-primary vertex; DCA_{xy} (cm); entries",800,-10.f,10.f);
    fHistDCAZpiprimaryMCt = new TH1F("fHistDCAZpiprimaryMCt","DCAz pion-primary vertex; DCA_{z} (cm); entries",800,-10.f,10.f);

    fHistDCAHe3primaryMCt = new TH1F("fHistDCAHe3primaryMCt","DCA ^{3}He-primary vertex; DCA (cm); entries",800,0.f,20.f);
    fHistDCAXYHe3primaryMCt = new TH1F("fHistDCAXYHe3primaryMCt","DCAxy ^{3}He-primary vertex; DCA_{xy} (cm); entries",800,-10.f,10.f);
    fHistDCAZHe3primaryMCt = new TH1F("fHistDCAZHe3primaryMCt","DCAz ^{3}He-primary vertex; DCA_{z} (cm); entries",800,-10.f,10.f);

    fHistDCAhe3pionMCt = new TH1F("fHistDCAhe3pionMCt","DCA ^{3}He-#pi tracks - MCtruth; DCA (cm);entries",550,-0.5,5.0);
    fHistDeltaPt_PionMCt = new TH2F("fHistDeltaPt_PionMCt","delta p_{T} #pi - MCtruth; p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
    fHistDeltaPt_He3MCt = new TH2F("fHistDeltaPt_He3MCt","delta p_{T} ^{3}He - MCtruth; p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
    fHistZDecayVtxMCt = new TH1F("fHistZDecayVtxMCt","decay vertex - z coordinate - MCtruth; z_{decay vtx} (cm); entries",800,-20.f,20.f);
    fHistXDecayVtxMCt = new TH1F("fHistXDecayVtxMCt","decay vertex - x coordinate - MCtruth; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
    fHistYDecayVtxMCt = new TH1F("fHistYDecayVtxMCt","decay vertex - y coordinate - MCtruth; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
    fHistDecayLengthH3L_MCt = new TH1F("fHistDecayLengthH3L_MCt","decay length ^{3}H_{#Lambda} - MCtruth; decay length (cm); entries",400,0.,400.);
    fHistNormalizedDecayL_MCt = new TH1F("fHistNormalizedDecayL_MCt","normalized decay length - MCtruth; decL/#sigma_{dL}; entries",400,0.,100.);
    fHistDCAXYhe3vtxMCt = new TH1F("fHistDCAXYhe3vtxMCt","DCA candidate p-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",1000,-5.,5.);
    fHistDCAZhe3vtxMCt = new TH1F("fHistDCAZhe3vtxMCt","DCA candidate p-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",2000,-10.,10.);
    fHistDCAhe3vtxMCt = new TH1F("fHistDCAhe3vtxMCt","DCA candidate p-decay vertex; DCA (cm); entries",1500,0.f,15.f);
    fHistDCAXYpionvtxMCt = new TH1F("fHistDCAXYpionvtxMCt","DCA candidate #pi^{-}-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",1000,-5.,5.);
    fHistDCAZpionvtxMCt = new TH1F("fHistDCAZpionvtxMCt","DCA candidate #pi^{-}-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",2000,-10.,10.);
    fHistDCApionvtxMCt = new TH1F("fHistDCApionvtxMCt","DCA candidate #pi-decay vertex; DCA (cm); entries",1500,0.f,15.f);
    fHistDeltaPt_HyperMCt = new TH2F("fHistDeltaPt_HyperMCt","delta p_{T} ^{3}H_{#Lambda} - MCtruth;p^{rec}_{T} (GeV/c);p^{MC}_{T} - p^{rec}_{T} (GeV/c)",130,0.,13.,200,-1.,1.);
    fHistLifetime_MCt = new TH1F("fHistLifetime_MCt","ct ^{3}H_{#Lambda} - MCtruth; ct(cm); entries",400,0.,40.);
    fHistAngle_He3_pion_MCt = new TH1F("fHistAngle_He3_pion_MCt","Angle between ^{3}He and #pi - MCtruth; #alpha_{3He_#pi} (rad) - MCtruth; entries/(0.03 rad)",100,0.,TMath::Pi());
    fHistHypertritonMomMCt = new TH1F("fHistHypertritonMomMCt","^{3}H_{#Lambda} momentum - MCtruth; p_{^{3}H_{#Lambda}} (GeV/c); entries/0.01",1100.,0.,11.);
    fHistPtPionMCt = new TH1F("fHistPtPionMCt","#pi p_{T}; p_{T} (GeV/c); entries",1000.,0.,10.);
    fHistPtHelium3MCt = new TH1F("fHistPtHelium3MCt","^{3}He p_{T}; p_{T} (GeV/c); entries",1000.,0.,10.);
    fHistPtHypertritonMCt = new TH1F("fHistPtHypertritonMCt","^{3}H_{#Lambda} p_{T} - MCtruth; p_{T} (GeV/c); entries",130,0.,13.);

    fHistHyperRapidityMCt =  new TH1F("fHistHyperRapidityMCt","rapidity distribution of ^{3}H_{#Lambda} - MC truth; y; entries",400,-2.f,2.f);
    fHistCosPointingAngleMCt = new TH1F("fHistCosPointingAngleMCt","Cosine of pointing angle distribution - MCtruth; cos(#theta_{pointing}); entries", 220, -1.1, 1.1);
    fHistMassHypertritonMCt = new TH1F("fHistMassHypertritonMCt","^{3}H_{#Lambda} invariant mass - MCtruth; invariant mass ^{3}He+#pi^{-} (GeV/c^{2}); entries", 500,2.9,3.4);
    fHistMassAntiHypertritonMCt = new TH1F("fHistMassAntiHypertritonMCt","#bar{^{3}H_{#Lambda}} invariant mass - MCtruth; invariant mass ^{3}#bar{He}+#pi^{+} (GeV/c^{2}); entries", 500,2.9,3.4);
    fHistMassVsPtMCt = new TH2F("fHistMassVsPtMCt","^{3}H_{#Lambda} and #bar{^{3}H_{#Lambda}} mass vs p_{T} - MCtruth;p_{T}(GeV/c); ^{3}He+#pi (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
    fHistMassHyperVsPtMCt = new TH2F("fHistMassHyperVsPtMCt","^{3}H_{#Lambda} invariant mass vs p_{T} - MCtruth;p_{T}(GeV/c); ^{3}He+#pi^{-} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);
    fHistMassAntiHyperVsPtMCt = new TH2F("fHistMassAntiHyperVsPtMCt","#bar{^{3}H_{#Lambda}} invariant mass vs p_{T} - MCtruth;p_{T}(GeV/c); ^{3}#bar{He}+#pi^{+} (GeV/c^{2})",100,0.,10.,500,2.9,3.4);


    fHistHypertritonMomGen = new TH1F("fHistHypertritonMomGen","Hypertriton p_{T}; p_{T} (GeV/c); entries",100,0.,10.);
    fHistHypertritonMomGen_3Body =  new TH1F("fHistHypertritonMomGen_3Body","Hypertriton p_{T} - 3b; p_{T} (GeV/c); entries",100,0.,10.);
    fHistAntiHypertritonMomGen_3Body =  new TH1F("fHistAntiHypertritonMomGen_3Body","anti-Hypertriton p_{T} - 3b; p_{T} (GeV/c); entries",100,0.,10.);
    fHistHypertritonMomGen_isPrimary_3Body = new TH1F("fHistHypertritonMomGen_isPrimary_3Body","Hypertriton p_{T} - 3b - isPrimary; p_{T} (GeV/c); entries",100,0.,10.);
    fHistHypertritonMomGen_2Body = new TH1F("fHistHypertritonMomGen_2Body","Hypertriton p_{T} - 2b; p_{T} (GeV/c); entries",100,0.,10.);
    fHistAntiHypertritonMomGen_2Body = new TH1F("fHistAntiHypertritonMomGen_2Body","anti-Hypertriton p_{T} - 2b; p_{T} (GeV/c); entries",100,0.,10.);
    fHistHypertritonMomGen_isPrimary_2Body = new TH1F("fHistHypertritonMomGen_isPrimary_2Body","Hypertriton p_{T} - 2b - isPrimary; p_{T} (GeV/c); entries",100,0.,10.);

    fHistHypertritonYGen = new TH1F("fHistHypertritonYGen","Hypertriton Y",200,-1.f,1.f);
    fHistHypertritonYGen_3Body = new TH1F("fHistHypertritonYGen_3Body","Hypertriton Y - 3b",200,-1.f,1.f);
    fHistAntiHypertritonYGen_3Body = new TH1F("fHistAntiHypertritonYGen_3Body","anti-Hypertriton Y - 3b",200,-1.f,1.f);
    fHistHypertritonYGen_isPrimary_3Body = new TH1F("fHistHypertritonYGen_isPrimary_3Body","Hypertriton Y - 3b - isPrimary",200,-1.f,1.f);
    fHistHypertritonYGen_2Body = new TH1F("fHistHypertritonYGen_2Body","Hypertriton Y - 2b",200,-1.f,1.f);
    fHistAntiHypertritonYGen_2Body = new TH1F("fHistAntiHypertritonYGen_2Body","anti-Hypertriton Y - 2b",200,-1.f,1.f);
    fHistHypertritonYGen_isPrimary_2Body = new TH1F("fHistHypertritonYGen_isPrimary_2Body","Hypertriton Y - 2b - isPrimary",200,-1.f,1.f);

  }

  fOutput->Add(fHistCount);
  fOutput->Add(fHistCentralityClass);
  fOutput->Add(fHistCentralityPercentile);
  fOutput->Add(fHistTrigger);
  fOutput->Add(fHistMultiplicity);
  fOutput->Add(fHistZPrimaryVtx);
  fOutput->Add(fHistXPrimaryVtx);
  fOutput->Add(fHistYPrimaryVtx);
  fOutput->Add(fHistChi2perTPCcluster);
  fOutput->Add(fHistTrackFlagReco);
  fOutput->Add(fHistTPCpid);
  fOutput->Add(fHistTPCHe3signal);
  fOutput->Add(fHistTPCpionsignal);
  fOutput->Add(fHistNsigmaHe3);
  fOutput->Add(fHistNsigmaPion);
  fOutput->Add(fHistTOFsignal);
  fOutput->Add(fHistTOFHe3signal);
  fOutput->Add(fHistpionTPCcls);
  //fOutput->Add(fHistCorrDCAHe3primary);
  //fOutput->Add(fHistCorrDCApiprimary);
  fOutput->Add(fHistDCApiprimary);
  fOutput->Add(fHistDCAXYpiprimary);
  fOutput->Add(fHistDCAZpiprimary);
  fOutput->Add(fHistDCAHe3primary);
  fOutput->Add(fHistDCAXYHe3primary);
  fOutput->Add(fHistDCAZHe3primary);
  fOutput->Add(fHistDCAhe3pion);
  fOutput->Add(fHistZDecayVtx);
  fOutput->Add(fHistXDecayVtx);
  fOutput->Add(fHistYDecayVtx);
  fOutput->Add(fHistDCAXYhe3vtx);
  fOutput->Add(fHistDCAZhe3vtx);
  fOutput->Add(fHistDCAhe3vtx);
  fOutput->Add(fHistDCAXYpionvtx);
  fOutput->Add(fHistDCAZpionvtx);
  fOutput->Add(fHistDCApionvtx);
  fOutput->Add(fHistDecayLengthH3L);
  fOutput->Add(fHistNormalizedDecayL);
  fOutput->Add(fHistDeltaPt_Hyper);
  fOutput->Add(fHistDeltaPt_Pion);
  fOutput->Add(fHistDeltaPt_He3);
  fOutput->Add(fHistLifetime);
  fOutput->Add(fHistAngle_He3_pion);
  fOutput->Add(fHistHyperRapidity);
  fOutput->Add(fHistCosPointingAngle);
  fOutput->Add(fHistPtPion);
  fOutput->Add(fHistPtHelium3);
  fOutput->Add(fHistPtHypertriton);
  fOutput->Add(fHistMassHypertriton);
  fOutput->Add(fHistMassAntiHypertriton);
  fOutput->Add(fHistMassVsPt);
  fOutput->Add(fHistMassHyperVsPt);
  fOutput->Add(fHistMassAntiHyperVsPt);

  if(fMC){
    fOutput->Add(fHistTPCdeusignal_pdg);
    fOutput->Add(fHistTPCtrisignal_pdg);
    fOutput->Add(fHistTPCHe3signal_pdg);
    fOutput->Add(fHistTPCHe4signal_pdg);
    fOutput->Add(fHistTPCHe3signal_3LH);
    fOutput->Add(fHistTPCpionsignal_3LH);
    fOutput->Add(fHistNsigmaHe3_3LH);
    fOutput->Add(fHistNsigmaPion_3LH);
    fOutput->Add(fHistParticle);
    fOutput->Add(fHistParticle_Mass);
    fOutput->Add(fHistpionTPCclsMCt);
    fOutput->Add(fHisthelium3TPCclsMCt);
    fOutput->Add(fHistpTpionMCt);
    fOutput->Add(fHistpThe3MCt);
    fOutput->Add(fHistMompionMCt);
    fOutput->Add(fHistMomHe3MCt);
    fOutput->Add(fHistCorrDCAHe3primaryMCt);
    fOutput->Add(fHistCorrDCApiprimaryMCt);
    fOutput->Add(fHistDCApiprimaryMCt);
    fOutput->Add(fHistDCAXYpiprimaryMCt);
    fOutput->Add(fHistDCAZpiprimaryMCt);
    fOutput->Add(fHistDCAHe3primaryMCt);
    fOutput->Add(fHistDCAXYHe3primaryMCt);
    fOutput->Add(fHistDCAZHe3primaryMCt);
    fOutput->Add(fHistDCAhe3pionMCt);
    fOutput->Add(fHistDeltaPt_PionMCt);
    fOutput->Add(fHistDeltaPt_He3MCt);
    fOutput->Add(fHistZDecayVtxMCt);
    fOutput->Add(fHistXDecayVtxMCt);
    fOutput->Add(fHistYDecayVtxMCt);
    fOutput->Add(fHistDecayLengthH3L_MCt);
    fOutput->Add(fHistNormalizedDecayL_MCt);
    fOutput->Add(fHistDCAXYhe3vtxMCt);
    fOutput->Add(fHistDCAZhe3vtxMCt);
    fOutput->Add(fHistDCAhe3vtxMCt);
    fOutput->Add(fHistDCAXYpionvtxMCt);
    fOutput->Add(fHistDCAZpionvtxMCt);
    fOutput->Add(fHistDCApionvtxMCt);
    fOutput->Add(fHistDeltaPt_HyperMCt);
    fOutput->Add(fHistLifetime_MCt);
    fOutput->Add(fHistAngle_He3_pion_MCt);
    fOutput->Add(fHistHypertritonMomMCt);
    fOutput->Add(fHistPtPionMCt);
    fOutput->Add(fHistPtHelium3MCt);
    fOutput->Add(fHistPtHypertritonMCt);
    fOutput->Add(fHistHyperRapidityMCt);
    fOutput->Add(fHistCosPointingAngleMCt);
    fOutput->Add(fHistMassHypertritonMCt);
    fOutput->Add(fHistMassAntiHypertritonMCt);
    fOutput->Add(fHistMassVsPtMCt);
    fOutput->Add(fHistMassHyperVsPtMCt);
    fOutput->Add(fHistMassAntiHyperVsPtMCt);
    fOutput->Add(fHistHypertritonMomGen);
    fOutput->Add(fHistHypertritonMomGen_3Body);
    fOutput->Add(fHistAntiHypertritonMomGen_3Body);
    fOutput->Add(fHistHypertritonMomGen_isPrimary_3Body);
    fOutput->Add(fHistHypertritonMomGen_2Body);
    fOutput->Add(fHistAntiHypertritonMomGen_2Body);
    fOutput->Add(fHistHypertritonMomGen_isPrimary_2Body);
    fOutput->Add(fHistHypertritonYGen);
    fOutput->Add(fHistHypertritonYGen_3Body);
    fOutput->Add(fHistAntiHypertritonYGen_3Body);
    fOutput->Add(fHistHypertritonYGen_isPrimary_3Body);
    fOutput->Add(fHistHypertritonYGen_2Body);
    fOutput->Add(fHistAntiHypertritonYGen_2Body);
    fOutput->Add(fHistHypertritonYGen_isPrimary_2Body);
  }

  // Post output data.
  PostData(1,fOutput);
  //printf("**************************************************\n");
  //printf("end of fOutput\n");
  //printf("**************************************************\n");

  if(fFillTree){
    OpenFile(2);
    fTTree = new TTree("hypertriton","hypertriton candidates");
    fTTree->Branch("MCtruth",&fTMCtruth,"MCtruth/O");
    fTTree->Branch("CentralityPerc",&fTCentralityPerc,"CentralityPerc/F");
    fTTree->Branch("Chi2NDFhe3",&fTchi2NDFhe3,"Chi2NDFhe3/F");
    fTTree->Branch("MassHypoTrkHe3",&fTMhypoTrkhe3,"MassHypoTrkHe3/s");
    fTTree->Branch("TPCclsPIDhe3",&fTPCclsPIDhe3,"TPCclsPIDhe3/s");
    fTTree->Branch("ITSrefithe3",&fITSrefithe3,"ITSrefithe3/O");
    fTTree->Branch("pTPChe3",&fTpTPChe3,"pTPChe3/F");
    fTTree->Branch("phe3_x",&fTpXhe3,"phe3_x/F");
    fTTree->Branch("phe3_y",&fTpYhe3,"phe3_y/F");
    fTTree->Branch("phe3_z",&fTpZhe3,"phe3_z/F");
    fTTree->Branch("TPCnsigmahe3",&fTTPCnsigmahe3,"TPCnsigmahe3/F");
    fTTree->Branch("DCAxyhe3prim",&fTDCAXYhe3prvtx,"DCAxyhe3prim/F");
    fTTree->Branch("DCAzhe3prim",&fTDCAZhe3prvtx,"DCAzhe3prim/F");
    fTTree->Branch("Chi2NDFpion",&fTchi2NDFpion,"Chi2NDFpion/F");
    fTTree->Branch("MassHypoTrkPion",&fTMhypoTrkpion,"MassHypoTrkPion/s");
    fTTree->Branch("TPCclsPIDpion",&fTPCclsPIDpion,"TPCclsPIDpion/s");
    fTTree->Branch("ITSrefitpion",&fITSrefitpion,"ITSrefitpion/O");
    fTTree->Branch("pTPCpion",&fTpTPCpion,"pTPCpion/F");
    fTTree->Branch("ppion_x",&fTpXpion,"ppion_x/F");
    fTTree->Branch("ppion_y",&fTpYpion,"ppion_y/F");
    fTTree->Branch("ppion_z",&fTpZpion,"ppion_z/F");
    fTTree->Branch("TPCnsigmapion",&fTTPCnsigmapion,"TPCnsigmapion/F");
    fTTree->Branch("DCAxypioprim",&fTDCAXYpioprvtx,"DCAxypioprim/F");
    fTTree->Branch("DCAzpioprim",&fTDCAZpioprvtx,"DCAzpioprim/F");
    fTTree->Branch("DCAhe3pi",&fTDCAhe3pi,"DCAhe3pi/F");
    fTTree->Branch("DCAXYhe3vtx",&fTDCAXYhevtx,"DCAXYhe3vtx/F");
    fTTree->Branch("DCAZhe3vtx",&fTDCAZhevtx,"DCAZhe3vtx/F");
    fTTree->Branch("DCAXYpionvtx",&fTDCAXYpivtx,"DCAXYpionvtx/F");
    fTTree->Branch("DCAZpionvtx",&fTDCAZpivtx,"DCAZpionvtx/F");
    fTTree->Branch("Angle_he3pi",&fTAngle_he3pi,"Angle_he3pi/F");
    fTTree->Branch("phe3Gen_x",&fTphe3_gen_X,"phe3Gen_x/F");
    fTTree->Branch("phe3Gen_y",&fTphe3_gen_Y,"phe3Gen_y/F");
    fTTree->Branch("phe3Gen_z",&fTphe3_gen_Z,"phe3Gen_z/F");
    fTTree->Branch("ppioGen_x",&fTppio_gen_X,"ppioGen_x/F");
    fTTree->Branch("ppioGen_y",&fTppio_gen_Y,"ppioGen_y/F");
    fTTree->Branch("ppioGen_z",&fTppio_gen_Z,"ppioGen_z/F");
    fTTree->Branch("pdgCandHe3",&fTpdgHe3,"pdgCandHe3/I");
    fTTree->Branch("pdgCandPion",&fTpdgPion,"pdgCandPion/I");
    fTTree->Branch("momId_he3",&fTmomidHe3,"momId_he3/I");
    fTTree->Branch("momId_pion",&fTmomidPi,"momId_pion/I");
    fTTree->Branch("pdgMomHe3",&fTpdgmomHe3,"pdgMomHe3/F");
    fTTree->Branch("pdgMomPion",&fTpdgmomPi,"pdgMomPion/F");
    fTTree->Branch("uniqueID_he3",&fTuniqID_he3,"uniqueID_he3/I");
    fTTree->Branch("uniqueID_pion",&fTuniqID_pion,"uniqueID_pion/I");
    fTTree->Branch("Rapidity",&fTRapidity,"Rapidity/F");
    fTTree->Branch("DecayLength",&fTDecayLength,"DecayLength/F");
    fTTree->Branch("DecayLengthError",&fTDecayLengthError,"DecayLengthError/F");
    fTTree->Branch("CosPA",&fTCosPA,"CosPA/F");
    fTTree->Branch("InvariantMass",&fTInvariantMass,"InvariantMass/F");

    fTTree->SetAutoSave(150000000);
    PostData(2,fTTree);
  } else {
    fTTree = new TTree();
  }

  if(fFillTGen){
    OpenFile(3);
    fTGen = new TTree("hypertriton_2","hypertriton simulated");
    fTGen->Branch("Centrality",&fTCentrality_gen,"Centrality/F");
    fTGen->Branch("HyperPx",&fTHyper_px,"HyperPx/F");
    fTGen->Branch("HyperPy",&fTHyper_py,"HyperPy/F");
    fTGen->Branch("HyperPz",&fTHyper_pz,"HyperPz/F");
    fTGen->Branch("HyperEta",&fTHyper_eta,"HyperEta/F");
    fTGen->Branch("HyperY",&fTHyper_rapidity,"HyperY/F");
    fTGen->Branch("HyperVx",&fTHyper_pv_x,"HyperVx/F");
    fTGen->Branch("HyperVy",&fTHyper_pv_y,"HyperVy/F");
    fTGen->Branch("HyperVz",&fTHyper_pv_z,"HyperVz/F");
    fTGen->Branch("HyperCharge",&fTHyper_charge,"HyperCharge/F");
    fTGen->Branch("He3Px",&fTHe3_px,"He3Px/F");
    fTGen->Branch("He3Py",&fTHe3_py,"He3Py/F");
    fTGen->Branch("He3Pz",&fTHe3_pz,"He3Pz/F");
    fTGen->Branch("He3Eta",&fTHe3_eta,"He3Eta/F");
    fTGen->Branch("He3Y",&fTHe3_rapidity,"He3Y/F");
    fTGen->Branch("He3Vx",&fTHe3_pv_x,"He3Vx/F");
    fTGen->Branch("He3Vy",&fTHe3_pv_y,"He3Vy/F");
    fTGen->Branch("He3Vz",&fTHe3_pv_z,"He3Vz/F");
    fTGen->Branch("He3Charge",&fTHe3_charge,"He3Charge/F");
    fTGen->Branch("PionPx",&fTPion_px,"PionPx/F");
    fTGen->Branch("PionPy",&fTPion_py,"PionPy/F");
    fTGen->Branch("PionPz",&fTPion_pz,"PionPz/F");
    fTGen->Branch("PionEta",&fTPion_eta,"PionEta/F");
    fTGen->Branch("PionY",&fTPion_rapidity,"PionY/F");
    fTGen->Branch("PionCharge",&fTPion_charge,"PionCharge/F");


    fTGen->SetAutoSave(150000000);
    PostData(3,fTGen);
  } else {
    fTGen = new TTree();
  }


  //printf("**************************************************\n");
  //printf("end of CreateOutputObjects\n");
  //printf("**************************************************\n");

} // end of UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::UserExec(Option_t *){
  // Main loop
  // Called for each event
  //----------------------------------------------------------------------------
  // Mass definition
  Double_t pionMass            =     0.13957; // pion mass = TDatabasePDG::Instance()->GetParticle(211)->Mass() GeV/c2
  Double_t helium3Mass         =     2.80923; // helium-3 mass
  Double_t hypertritonMass     =     2.991106;

  //define PDGCodes
  Long_t pdgPionPlus           =              211;
  Long_t pdgPionMinus          =             -211;
  Long_t pdgHelium3            =       1000020030;
  Long_t pdgAntiHelium3        =      -1000020030;
  Long_t pdgHypertriton        =       1010010030;
  Long_t pdgAntiHypertriton    =      -1010010030;
  //----------------------------------------------------------------------------


  //==========Define variables==========//
  //===PID loop===//
  Bool_t useTOF = kFALSE;
  Int_t ntracks,label = 0;
  Double_t chi2PerClusterTPC, nClustersTPC=0.;
  Double_t p, pOverZ, pT, ptot = 0.;
  Float_t beta = 0.;
  Float_t dca_prim= 0.;
  Float_t dcaprim[2] = {0.,0.}; //array for the dca from primary vertex coordinates (xy,z)
  Float_t dcaprimc[3] = {0.,0.,0.}; // covariance matrix elements of the TPC only impact parameter
  AliESDtrack *track = 0x0;


  //==========ESD event==========
  fESDevent=(AliESDEvent*)InputEvent();
  if(!fESDevent){
    printf("AliAnalysisTaskHypCrossCheck::Exec(): bad ESD\n");
    PostData(1,fOutput);
    return;
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler*)mgr->GetInputEventHandler();

  fHistCount->Fill(0); // number of reco events opened

 //==========MC info==========
  if(fMC){//MC info and sample selection
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (mgr->GetMCtruthEventHandler());
    if (!eventHandler) {
      printf("ERROR: Could not retrieve MC event handler");
      PostData(1,fOutput);
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      printf("ERROR: Could not retrieve MC event");
      PostData(1,fOutput);
      return;
    }
    fStack = mcEvent->Stack();
    if (!fStack) {
      printf("ERROR: fStack not available\n");
      PostData(1,fOutput);
      return;
    }
  } // end of MC info

  //==========Trigger class==========
  UInt_t maskPhysSel = handl->IsEventSelected();
  if(!PassTriggerSelection(maskPhysSel)){
    PostData(1,fOutput);
    return;
  }

  fHistCount->Fill(1); // number of events passing the Trigger Selection
  //==========Centrality==========
  if(!PassCentralitySelection()) {
    PostData(1,fOutput);
    return; //0 bis 80 %
  }
  fHistCount->Fill(2); // number of reco events passing Centrality Selection


  //==========Multiplicity==========
  Int_t refMultTpc = AliESDtrackCuts::GetReferenceMultiplicity(fESDevent, kTRUE);
  fHistMultiplicity->Fill(refMultTpc);

  //==========Primary Vertex==========
  fPrimaryVertex = (AliESDVertex*)fESDevent->GetPrimaryVertexTracks();
  if (fPrimaryVertex->GetNContributors()<1) {
      // SPD vertex
    fPrimaryVertex = (AliESDVertex*)fESDevent->GetPrimaryVertexSPD();
    if(fPrimaryVertex->GetNContributors()<1) {
      PostData(1,fOutput);
      return;
    }
  }

  if (TMath::Abs(fPrimaryVertex->GetZ()) > 10) {
    PostData(1,fOutput);
    return;
  }

  fHistCount->Fill(3); // number of reco events passing Primary Vertex Selection

  fHistZPrimaryVtx->Fill(fPrimaryVertex->GetZ());
  fHistXPrimaryVtx->Fill(fPrimaryVertex->GetX());
  fHistYPrimaryVtx->Fill(fPrimaryVertex->GetY());

  //********TREE GENERATO*********
  //---------------- MC Generation ----------------
  if(fMC){
    CheckGenerated();
  } // end of if(fMC)


  // -------------------------------------------------------
  // Loop for PID on ESD tracks
  // -------------------------------------------------------

  fPIDResponse = handl->GetPIDResponse();

  ntracks = fESDevent->GetNumberOfTracks();

  TArrayI chelium3(ntracks);
  UInt_t nHe3TPC = 0;
  TArrayI cantihelium3(ntracks);
  UInt_t nAntiHe3TPC = 0;
  TArrayI cpionplus(ntracks);
  UInt_t nPioPlusTPC = 0;
  TArrayI cpionminus(ntracks);
  UInt_t nPioMinusTPC = 0;

  Bool_t positive = kFALSE;
  Bool_t negative = kFALSE;
  Int_t  checkPdg = 0.;
  Bool_t isRealHyp2 = kFALSE;
  //vector <Float_t> cmassd;
  //vector <Float_t> cmassp;


  for(Int_t i=0; i < ntracks; i++) {
    track = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(i));
    beta = 0.;
    useTOF = kFALSE;
    positive = kFALSE;
    negative = kFALSE;
    isRealHyp2 = kFALSE;
    ULong_t status = track->GetStatus();
    // Chi2/TPCcls

    nClustersTPC = track->GetTPCclusters(0);
    chi2PerClusterTPC = track->GetTPCchi2()/nClustersTPC;
    fHistChi2perTPCcluster->Fill(chi2PerClusterTPC);

    if(fMC){if(track->GetLabel()<0) continue;}
    if(track->GetID()<0) continue;

    if(!fESDtrackCuts->AcceptTrack(track)) continue;
    //if((status&AliVTrack::kITSin) && !(status&AliVTrack::kITSrefit)) continue;
    if(!track->GetInnerParam()) continue;

    if(track->GetTPCsignalN()<80) continue;

    p = track->P(); //track->GetTPCmomentum()
    ptot = track->GetInnerParam()->GetP();
    pOverZ = ptot*track->GetSign();
    pT = track->Pt();
    if(track->GetSign() > 0) positive = kTRUE;
    else negative = kTRUE;

    fHistTPCpid->Fill(pOverZ, track->GetTPCsignal());
    /*if(fRequestTOFPid){
        useTOF = HasTOF(track, beta);
        if(useTOF)fHistTOFsignal->Fill(p,beta);
    }*/
    //Filling PID histo
    label = track->GetLabel();

    fHistTrackFlagReco->Fill(0);
    if(status&AliVTrack::kITSin) fHistTrackFlagReco->Fill(1);
    if(status&AliVTrack::kITSout) fHistTrackFlagReco->Fill(2);
    if(status&AliVTrack::kITSrefit) fHistTrackFlagReco->Fill(3);
    if((status&AliVTrack::kITSout) && !(status&AliVTrack::kITSrefit)) fHistTrackFlagReco->Fill(5);
    if((status&AliVTrack::kITSin) && !(status&AliVTrack::kITSrefit)){
      fHistTrackFlagReco->Fill(4);
      //continue;
    }
    track->GetImpactParameters(dcaprim,dcaprimc);
    dca_prim = TMath::Sqrt((dcaprim[0]*dcaprim[0])+(dcaprim[1]*dcaprim[1]));
    if(fMC) isRealHyp2 = CheckPrimaryDistribution(track,checkPdg);

    if(isRealHyp2){
      if(TMath::Abs(checkPdg)==pdgPionPlus){
        fHistpTpionMCt->Fill(track->Pt());
        fHistMompionMCt->Fill(track->P());
        fHistpionTPCclsMCt->Fill(track->GetTPCclusters(0));
        fHistNsigmaPion_3LH->Fill(pOverZ,fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
        fHistTPCpionsignal_3LH->Fill(pOverZ, track->GetTPCsignal());
      }
      if(TMath::Abs(checkPdg)==pdgHelium3){
        fHistpThe3MCt->Fill(track->Pt());
        fHistMomHe3MCt->Fill(track->P());
        fHisthelium3TPCclsMCt->Fill(track->GetTPCclusters(0));
        fHistNsigmaHe3_3LH->Fill(pOverZ,fPIDResponse->NumberOfSigmasTPC(track,AliPID::kHe3));
        fHistTPCHe3signal_3LH->Fill(pOverZ, track->GetTPCsignal());
      }
    }
    if(fMC){
        TParticle *tparticleDaughter = fStack->Particle(TMath::Abs(label));
        if(tparticleDaughter->GetPdgCode() == 1000010020 || tparticleDaughter->GetPdgCode() == -1000010020) fHistTPCdeusignal_pdg->Fill(pOverZ, track->GetTPCsignal());
        if(tparticleDaughter->GetPdgCode() == 1000010030 || tparticleDaughter->GetPdgCode() == -1000010030) fHistTPCtrisignal_pdg->Fill(pOverZ, track->GetTPCsignal());
        if(tparticleDaughter->GetPdgCode() == 1000020030 || tparticleDaughter->GetPdgCode() == -1000020030) fHistTPCHe3signal_pdg->Fill(pOverZ, track->GetTPCsignal());
        if(tparticleDaughter->GetPdgCode() == 1000020040 || tparticleDaughter->GetPdgCode() == -1000020040) fHistTPCHe4signal_pdg->Fill(pOverZ, track->GetTPCsignal());


        //if(pT > 7) continue;
        if(pT >= 0. && TMath::Abs(dcaprim[1]) <= fDCAzHe3PVmax){
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kHe3)) <= 3 ) { // p  TMath::Abs(tparticleDaughter->GetPdgCode()) == 1000020030  || tparticleDaughter->GetPdgCode() == -2212

              fHistNsigmaHe3->Fill(pOverZ,fPIDResponse->NumberOfSigmasTPC(track,AliPID::kHe3));
              fHistTPCHe3signal->Fill(pOverZ, track->GetTPCsignal());

              fHistDCAHe3primary->Fill(dca_prim);
              fHistDCAXYHe3primary->Fill(dcaprim[0]);
              fHistDCAZHe3primary->Fill(dcaprim[1]);

              if(isRealHyp2 && TMath::Abs(checkPdg)==pdgHelium3){
                fHistDCAHe3primaryMCt->Fill(dca_prim);
                fHistDCAXYHe3primaryMCt->Fill(dcaprim[0]);
                fHistDCAZHe3primaryMCt->Fill(dcaprim[1]);
              }
              if(useTOF){
                fHistTOFHe3signal->Fill(p,beta);
                //fHistTOFpromass->Fill(mass);
              }

          if(positive) chelium3[nHe3TPC++] = i;
          if(negative) cantihelium3[nAntiHe3TPC++] = i;
          //chelium3.push_back(track);
          //cmassp.push_back(mass); // da spostare nel if length > 350.?
          }
        }

        if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)) <= 3) { // pi+ TMath::Abs(tparticleDaughter->GetPdgCode()) == 211 tparticleDaughter->GetPdgCode() == 211 ||

            fHistNsigmaPion->Fill(pOverZ,fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
            fHistTPCpionsignal->Fill(pOverZ, track->GetTPCsignal());

            fHistDCApiprimary->Fill(dca_prim);
            fHistDCAXYpiprimary->Fill(dcaprim[0]);
            fHistDCAZpiprimary->Fill(dcaprim[1]);
            fHistpionTPCcls->Fill(track->GetTPCclusters(0));

            if(isRealHyp2 && TMath::Abs(checkPdg)==pdgPionPlus){
              fHistDCApiprimaryMCt->Fill(dca_prim);
              fHistDCAXYpiprimaryMCt->Fill(dcaprim[0]);
              fHistDCAZpiprimaryMCt->Fill(dcaprim[1]);
            }

          if(!fESDtrackCutsV0->AcceptTrack(track)) continue;
          if(dca_prim < fDCAPiPVmin) continue;

          if(positive) cpionplus[nPioPlusTPC++] = i;
          if(negative) cpionminus[nPioMinusTPC++] = i;

        }
    } // end of MC PID
    else{
      //AliPID::EParticleType parType = AliPID::kHe3;
      if(TMath::Abs(GetBetheAlephTPCnsigmas(track,AliPID::kHe3)) <= fRequestTPCSigmas) { // helium-3 PassPIDSelection(track,AliPID::kHe3, useTOF)
        fHistTPCHe3signal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
        fHistDCAHe3primary->Fill(dca_prim);
        fHistDCAXYHe3primary->Fill(dcaprim[0]);
        fHistDCAZHe3primary->Fill(dcaprim[1]);
        if(useTOF){
          fHistTOFHe3signal->Fill(p,beta);
          //fHistTOFpromass->Fill(mass);
          //cmassp.push_back(mass);
          }
        if(positive) chelium3[nHe3TPC++] = i;
        if(negative) cantihelium3[nAntiHe3TPC++] = i;
      }

        if(!fESDtrackCutsV0->AcceptTrack(track)) continue;
        if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)) <= 3) { //pion^+
            fHistTPCpionsignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
            fHistDCApiprimary->Fill(dca_prim);
            fHistDCAXYpiprimary->Fill(dcaprim[0]);
            fHistDCAZpiprimary->Fill(dcaprim[1]);
            fHistpionTPCcls->Fill(track->GetTPCclusters(0));
            if(dca_prim < fDCAPiPVmin) continue;
            if(positive) cpionplus[nPioPlusTPC++] = i;
            if(negative) cpionminus[nPioMinusTPC++] = i;
        }
    } // end of Data PID (or MC PID done as on Data)
} // end of PID loop

chelium3.Set(nHe3TPC);
cantihelium3.Set(nAntiHe3TPC);
cpionplus.Set(nPioPlusTPC);
cpionminus.Set(nPioMinusTPC);

if(fMinvSignal){
  //Hypertriton Invariant Mass
  if(fChooseMatter) CombineTwoTracks(kTRUE, chelium3, cpionminus,2);

  //Anti-Hypertriton Invariant Mass
  if(fChooseAntiMatter) CombineTwoTracks(kFALSE, cantihelium3, cpionplus,2);
}

if(fMinvLikeSign){
  //Hypertriton Invariant Mass - LS
  if(fChooseMatter) CombineTwoTracks(kTRUE, chelium3, cpionplus,2);

  //Anti-Hypertriton Invariant Mass - LS
  if(fChooseAntiMatter) CombineTwoTracks(kFALSE, cantihelium3, cpionminus,2);
}

chelium3.Reset();
cantihelium3.Reset();
cpionplus.Reset();
cpionminus.Reset();

  //cmassd.clear();
  //cmassp.clear();

  // Post output data.
  PostData(1,fOutput);
  if(fFillTree) PostData(2,fTTree);
  if(fFillTGen) PostData(3,fTGen);


} // end of UserExec

//________________________________________________________________________
void AliAnalysisTaskHypCrossCheck::Terminate(Option_t *){
  // Merge output
  // Called once at the end of the query

  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  printf("end of Terminate");
  return;

} // end of Terminate
