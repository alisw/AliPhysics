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
// AliAnalysisTaskHypertriton3 class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
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
#include "AliAnalysisTaskHypertriton3.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliVertexerTracks.h"
#include "AliVEvent.h"
#include "AliVTrack.h"


ClassImp(AliAnalysisTaskHypertriton3)

using std::cout;
using std::endl;
using std::vector;

/* $Id$ */
//________________________________________________________________________
AliAnalysisTaskHypertriton3::AliAnalysisTaskHypertriton3(TString taskname):
  AliAnalysisTaskSE(taskname.Data()),
  fESDevent(0),
  fESDtrackCuts(0x0),
  fESDtrackCutsV0(0x0),
  fPrimaryVertex(0x0),
  fPIDResponse(0x0),
  fVertexer(0x0),
  fVtx1(0x0),
  fVtx2(0x0),
  fTrkArray(0x0),
  fMC(kFALSE),
  fFillTree(kFALSE),
  fRun1PbPb(kTRUE),
  fRun2PbPb(kFALSE),
  fCutMassUp(3.1),
  fRequireMassRange(kFALSE),
  fEvtSpecie(4),
  fEvtEmbedSelection(kFALSE),
  fCentrality(0x0),
  fCentralityPercentile(0x0),
  fCentralityClass(0x0),
  fTriggerConfig(1),
  fRequireITSclusters(kFALSE),
  fMinITSclustersN(0),
  fRequireITSrefit(kFALSE),
  fRequireITSrefitPion(kFALSE),
  fRequireITSin(kFALSE),
  fPionTPCSigmas(3),
  fProtonTPCSigmas(3),
  fDeuteronTPCSigmas(3),
  fRequireTOFPid(kFALSE),
  fRequestTOFSigmas(3),
  fChooseMatter(kTRUE),
  fChooseAntiMatter(kTRUE),
  fMinvSignal(kTRUE),
  fMinvLikeSign(kFALSE),
  fSideBand(kFALSE),
  fTriangularDCAtracks(kFALSE),
  fMinPtDeuteron(0),
  fMaxPtDeuteron(10.),
  fMinPtProton(0),
  fMaxPtProton(10.),
  fDCAPiPVmin(0.1),
  fDCAzPPVmax(999.),
  fDCAzDPVmax(999.),
  fCosPointingAngle(0.998),
  fMaxDecayLength(15.),
  fMinDecayLength(0.),
  fMinNormalizedDecL(0.),
  fMaxLifeTime(9999.),
  fMinLifeTime(0.),
  fRapidity(0.5),
  fMaxPtMother(10.),
  fMinPtMother(2.),
  fDCAPiSVxymax(0.6),
  fDCAPiSVzmax(0.8),
  fDCAProSVmax(0.7),
  fDCADeuSVmax(0.6),
  fDCAdp(0.2),
  fDCApip(0.5),
  fDCAdpi(0.5),
  fAngledp(TMath::Pi()),
  fAngledpi(TMath::Pi()),
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
  fHistITSclusters(0x0),
  fHistTPCpid(0x0),
  fHistTPCdeusignal(0x0),
  fHistTPCprosignal(0x0),
  fHistTPCpionsignal(0x0),
  fHistTOFsignal(0x0),
  fHistTOFdeusignal(0x0),
  fHistTOFprosignal(0x0),
  fHistTOFpionsignal(0x0),
  //fHistTOFdeumass(0x0),
  //fHistTOFpromass(0x0),
  fHistpionTPCcls(0x0),
  //fHistCorrDCAdprimary(0x0),
  //fHistCorrDCApprimary(0x0),
  //fHistCorrDCApiprimary(0x0),
  fHistDCApiprimary(0x0),
  fHistDCAXYpiprimary(0x0),
  fHistDCAZpiprimary(0x0),
  fHistDCApprimary(0x0),
  fHistDCAXYpprimary(0x0),
  fHistDCAZpprimary(0x0),
  fHistDCAdprimary(0x0),
  fHistDCAXYdprimary(0x0),
  fHistDCAZdprimary(0x0),
  fHistDCAdeupro(0x0),
  fHistDCApiondeu(0x0),
  fHistDCApionpro(0x0),
  fHistDCAdpdpi(0x0),
  fHistDCApdppi(0x0),
  fHistDCApidpip(0x0),
  fHistZDecayVtx(0x0),
  fHistXDecayVtx(0x0),
  fHistYDecayVtx(0x0),
  fHistDCApionvtx(0x0),
  fHistDCAXYpionvtx(0x0),
  fHistDCAZpionvtx(0x0),
  fHistSigmaDcaXYpionvtx(0x0),
  fHistSigmaDcaZpionvtx(0x0),
  fHistSigmaDcapionvtx(0x0),
  fHistDCAprovtx(0x0),
  fHistDCAXYprovtx(0x0),
  fHistDCAZprovtx(0x0),
  fHistSigmaDcaXYprovtx(0x0),
  fHistSigmaDcaZprovtx(0x0),
  fHistSigmaDcaprovtx(0x0),
  fHistDCAdeuvtx(0x0),
  fHistDCAXYdeuvtx(0x0),
  fHistDCAZdeuvtx(0x0),
  fHistSigmaDcaXYdeuvtx(0x0),
  fHistSigmaDcaZdeuvtx(0x0),
  fHistSigmaDcadeuvtx(0x0),
  fHistDecayLengthH3L(0x0),
  fHistNormalizedDecayL(0x0),
  fHistLifetime(0x0),
  fHistAngle_deu_pro(0x0),
  fHistAngle_deu_pion(0x0),
  fHistAngle_pro_pion(0x0),
  fHistAngleCorr_dp_dpi(0x0),
  fHistAngleCorr_dp_ppi(0x0),
  fHistAngleCorr_ppi_dpi(0x0),
  fHistHyperRapidity(0x0),
  fHistCosPointingAngle(0x0),
  fHistPtPion(0x0),
  fHistPtDeuteron(0x0),
  fHistPtProton(0x0),
  fHistPtHypertriton(0x0),
  fHistMassHypertriton(0x0),
  fHistMassAntiHypertriton(0x0),
  fHistMassHypertriton_Cent(0x0),
  fHistMassAntiHypertriton_Cent(0x0),
  fHistMassHypertriton_SemiCent(0x0),
  fHistMassAntiHypertriton_SemiCent(0x0),
  fHistMassHyp_Lifetime(0x0),
  fHistMassHyp_Lifetime_M(0x0),
  fHistMassHyp_Lifetime_A(0x0),
  fHistParticle(0x0),
  fHistpionTPCclsMCt(0x0),
  fHistpTpionMCt(0x0),
  fHistpTproMCt(0x0),
  fHistpTdeuMCt(0x0),
  fHistMompionMCt(0x0),
  fHistMomproMCt(0x0),
  fHistMomdeuMCt(0x0),
  fHistCorrDCAdprimaryMCt(0x0),
  fHistCorrDCApprimaryMCt(0x0),
  fHistCorrDCApiprimaryMCt(0x0),
  fHistDCApiprimaryMCt(0x0),
  fHistDCApprimaryMCt(0x0),
  fHistDCAdprimaryMCt(0x0),
  fHistDCAdeuproMCt(0x0),
  fHistDCApiondeuMCt(0x0),
  fHistDCApionproMCt(0x0),
  fHistZDecayVtxMCt(0x0),
  fHistXDecayVtxMCt(0x0),
  fHistYDecayVtxMCt(0x0),
  fHistDCAXYdeuvtxMCt(0x0),
  fHistDCAZdeuvtxMCt(0x0),
  fHistDCAXYprovtxMCt(0x0),
  fHistDCAZprovtxMCt(0x0),
  fHistDCAXYpionvtxMCt(0x0),
  fHistDCAZpionvtxMCt(0x0),
  fHistNormalizedDecayL_MCt(0x0),
  fHistLifetime_MCt(0x0),
  fHistAngle_deu_pro_MCt(0x0),
  fHistAngle_deu_pion_MCt(0x0),
  fHistAngle_pro_pion_MCt(0x0),
  fHistAngleCorr_dp_dpi_MCt(0x0),
  fHistAngleCorr_dp_ppi_MCt(0x0),
  fHistAngleCorr_ppi_dpi_MCt(0x0),
  fHistHypertritonMomMCt(0x0),
  fHistHyperRapidityMCt(0x0),
  fHistMassHypertritonMCt(0x0),
  fHistMassAntiHypertritonMCt(0x0),
  fTTree(0x0),
  fTMCtruth(kFALSE),
  fTCentralityPerc(0x0),
  fTchi2NDFdeu(0x0),
  fTPCclsdeu(0x0),
  fITSdclsmap(0x0),
  fTpTdeu(0x0),
  fTpdeu(0x0),
  fTTPCnsigmadeu(0x0),
  fTTOFnsigmadeu(0x0),
  fTchi2NDFpro(0x0),
  fTPCclspro(0x0),
  fITSpclsmap(0x0),
  fTpTpro(0x0),
  fTppro(0x0),
  fTTPCnsigmapro(0x0),
  fTTOFnsigmapro(0x0),
  fTchi2NDFpion(0x0),
  fTPCclspion(0x0),
  fITSpiclsmap(0x0),
  fTpTpion(0x0),
  fTppion(0x0),
  fTTPCnsigmapion(0x0),
  fTTOFnsigmapion(0x0),
  fTDCAdp(0x0),
  fTDCAdpi(0x0),
  fTDCAppi(0x0),
  fTDCAXYdvtx(0x0),
  fTDCAZdvtx(0x0),
  fTDCAXYpvtx(0x0),
  fTDCAZpvtx(0x0),
  fTDCAXYpivtx(0x0),
  fTDCAZpivtx(0x0),
  fTAngle_dp(0x0),
  fTAngle_dpi(0x0),
  fTAngle_ppi(0x0),
  fTpdeu_gen_X(0x0),
  fTpdeu_gen_Y(0x0),
  fTpdeu_gen_Z(0x0),
  fTppro_gen_X(0x0),
  fTppro_gen_Y(0x0),
  fTppro_gen_Z(0x0),
  fTppio_gen_X(0x0),
  fTppio_gen_Y(0x0),
  fTppio_gen_Z(0x0),
  fTpdgDeu(0x0),
  fTpdgPro(0x0),
  fTpdgPion(0x0),
  fTmomidD(0x0),
  fTmomidP(0x0),
  fTmomidPi(0x0),
  fTpdgmomD(0x0),
  fTpdgmomP(0x0),
  fTpdgmomPi(0x0),
  fTuniqID_deu(0x0),
  fTuniqID_pro(0x0),
  fTuniqID_pion(0x0),
  fTRapidity(0x0),
  fTDecayLengthProper(0x0),
  fTDecayLengthNorm(0x0),
  fTCosPA(0x0),
  fTTransverseMom(0x0),
  fTInvariantMass(0x0)
{
  //Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

  //ESD Track cuts
  if(!fESDtrackCuts) fESDtrackCuts = new AliESDtrackCuts();
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(fRequireITSrefit);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);

  //ESD Track cuts V0
  if(!fESDtrackCutsV0) fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
  fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsV0->SetMinNClustersTPC(100);
  fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
  fESDtrackCutsV0->SetRequireITSRefit(fRequireITSrefitPion);
  fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
  fESDtrackCutsV0->SetPtRange(0.2,1.2);

}


//________________________________________________________________________
AliAnalysisTaskHypertriton3::~AliAnalysisTaskHypertriton3(){
  //Destructor
  if(fOutput){
    delete fOutput;
    fOutput = 0;
  }

  if(fTTree) delete fTTree;

    if(fPIDResponse){
    delete fPIDResponse;
  }

    if(fESDtrackCuts) delete fESDtrackCuts;
    if(fESDtrackCutsV0) delete fESDtrackCutsV0;
    if(fPrimaryVertex) delete fPrimaryVertex;
    if(fVertexer) delete fVertexer;
    if(fTrkArray) delete fTrkArray;
    if(fVtx1) delete fVtx1;
    if(fVtx2) delete fVtx2;


} // end of Destructor

//________________________________________________________________________
Double_t AliAnalysisTaskHypertriton3::GetDCAcut(Int_t part, Double_t dca)const{
  Double_t cut;
  if(part == 5){ // p-d vs pi-d
    cut = fDCAdpi*(1-(dca/fDCAdp));
    return cut;
  }
  else if(part == 4){ // p-d vs p-pi
    cut = fDCApip*(1-(dca/fDCAdp));
    return cut;
  }
  return -1;
}


//________________________________________________________________________
void AliAnalysisTaskHypertriton3::SetConvertedAODVertices(AliESDVertex *ESDvtxp, AliESDVertex *ESDvtxs)const{

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
Bool_t AliAnalysisTaskHypertriton3::PassTriggerSelection(UInt_t PhysSelMask){

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
Bool_t AliAnalysisTaskHypertriton3::PassCentralitySelection(){

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
    if(fESDevent->GetEventSpecie() == fEvtSpecie){ // Event Specie == 4 == PbPb
      AliWarning("fRun2PbPb: Centrality task");
      AliMultSelection *centr=(AliMultSelection*)fESDevent->FindListObject("MultSelection");
      if(!centr){
        AliWarning("AliMultSelection object not found!");
      }
      fCentralityPercentile = centr->GetMultiplicityPercentile("V0M",fEvtEmbedSelection);
    }
  }

  fHistCentralityClass->Fill(fCentrality);
  fHistCentralityPercentile->Fill(fCentralityPercentile);

  if (fCentralityPercentile < fLowCentrality || fCentralityPercentile > fHighCentrality) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypertriton3::PassTrackSelection(AliESDtrack *trk){

  if(!fESDtrackCuts->AcceptTrack(trk)) return kFALSE;
  if(!trk->GetInnerParam()) return kFALSE;
  if(trk->GetTPCsignalN()<80) return kFALSE;

  fHistTrackFlagReco->Fill(0);
  ULong_t status = trk->GetStatus();
  if(status&AliVTrack::kITSin) fHistTrackFlagReco->Fill(1);
  if(status&AliVTrack::kITSout) fHistTrackFlagReco->Fill(2);
  if(status&AliVTrack::kITSrefit) fHistTrackFlagReco->Fill(3);
  if((status&AliVTrack::kITSout) && !(status&AliVTrack::kITSrefit)) fHistTrackFlagReco->Fill(5);
  if((status&AliVTrack::kITSin) && !(status&AliVTrack::kITSrefit)){
    fHistTrackFlagReco->Fill(4);
    if(fRequireITSin) return kFALSE;
  }

  if(fRequireITSclusters && trk->GetITSNcls()<fMinITSclustersN) return kFALSE;
  fHistITSclusters->Fill(trk->GetITSNcls());

  return kTRUE;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskHypertriton3::HasTOF(AliESDtrack *trk, float &beta_tof){
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
Bool_t AliAnalysisTaskHypertriton3::PassPIDSelection(AliESDtrack *trk, int specie, Bool_t isTOFin, Float_t nsigma_cut){
    //PID selection
    bool tofPID = kTRUE, tpcPID = kTRUE;
    //TPC-pid
    float const nsigmaTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::EParticleType (specie)));
    if(nsigmaTPC > nsigma_cut) tpcPID = kFALSE;
    else tpcPID = kTRUE;
    //TOF-pid
    if(fRequireTOFPid){
        if(isTOFin){
            float const nsigmaTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::EParticleType (specie)));
            if(nsigmaTOF>fRequestTOFSigmas) tofPID = kFALSE;
            else tofPID = kTRUE;
        }
    }
    return tpcPID && tofPID;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHypertriton3::ComputeSigma(Double_t dc[2], Double_t dc_cov[3]){
  Double_t d = (dc[0]*dc[0])+(dc[1]*dc[1]);
  Double_t t_xy = (dc[0]*dc[0])*dc_cov[0];
  Double_t t_z = (dc[1]*dc[1])*dc_cov[1];
  Double_t t_cov = dc[0]*dc[1]*dc_cov[2];

  Double_t sgm_d = TMath::Sqrt((t_xy + t_z + t_cov)/d);

  return sgm_d;
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::CombineThreeTracks(Bool_t isMatter, TArrayI arrD, TArrayI arrP, TArrayI arrPi, Bool_t cent0, Bool_t cent1){
//Method to combine tracks
//Here implemented topological and kinematical cuts
//Total charge predefined with the correct type of array assigned to this method

//Define variables
//===Combining tracks loop===//
Double_t bz = fESDevent->GetMagneticField();
fVertexer->SetFieldkG(bz);

Double_t dlh[3] = {0,0,0}; //array for the coordinates of the decay length
Double_t dca_dp, dca_dpi, dca_ppi, angle_dp, angle_dpi, angle_ppi = 0.;
Double_t dcad[2] = {0.,0.}; // dca between the candidate d,p,pi
Double_t dcap[2] = {0.,0.}; // and the candidate decay vertex
Double_t dcapi[2] = {0.,0.}; // dcad[0]= transverse plane coordinate; dcad[1]= z coordinate
Double_t dcapi_cov[3] = {0.,0.,0.};
Double_t dcap_cov[3] = {0.,0.,0.};
Double_t dcad_cov[3] = {0.,0.,0.};
Double_t decayLengthH3L, normalizedDecayL, rapidity, pointingAngleH, ctau= 0.;
Double_t lD, lP, lPi = 0;
Double_t decVt[3] = {0.,0.,0.};
Bool_t brotherHood = kFALSE;
TLorentzVector posD, posP, negPi; //Lorentz vector of deuteron, proton and pion in the LAB
AliESDtrack *trackD = 0x0;
AliESDtrack *trackP = 0x0;
AliESDtrack *trackNPi = 0x0;

Double_t xthiss(0.0);
Double_t xpp(0.0);
Float_t piprim[2] = {0.,0.};
Float_t piprimc[3] = {0.,0.,0.};
Float_t nsd, nsp, nspi = 0.;
Float_t nsd_t, nsp_t, nspi_t, b_t = 0.;
AliESDVertex *decayVtx = 0x0;

TLorentzVector Hypertriton;
TVector3 h1, d1, p1, pi1;
Double_t pTotHyper = 0.;

TParticle *tparticleD = 0x0;
TParticle *tparticleP = 0x0;
TParticle *tparticlePi = 0x0;


// -------------------------------------------------------
// Loop for Invariant Mass
// -------------------------------------------------------


for(Int_t j=0; j<arrD.GetSize(); j++){ // candidate deuteron loop cdeuteron.size()

  trackD = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(arrD[j]));


  for(Int_t m=0; m<arrP.GetSize(); m++){ // candidate proton loop cproton.size()

    trackP = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(arrP[m]));

    if(trackD->GetID() == trackP->GetID()) continue;

    dca_dp = trackD->GetDCA(trackP,bz,xthiss,xpp);

    fHistDCAdeupro->Fill(dca_dp);

    if(dca_dp > fDCAdp) continue;


    for(Int_t s=0; s<arrPi.GetSize(); s++ ){ // candidate pion loop cpion.size()

      fTrkArray->Clear();
      Hypertriton.Clear();
      posD.Clear();
      posP.Clear();
      negPi.Clear();
      h1.Clear();
      d1.Clear();
      p1.Clear();
      pi1.Clear();

      trackNPi = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(arrPi[s]));
      brotherHood = kFALSE;


      if(trackNPi->GetID() == trackP->GetID()) continue;
      if(trackNPi->GetID() == trackD->GetID()) continue;


      dca_dpi = trackNPi->GetDCA(trackD,bz,xthiss,xpp);
      dca_ppi = trackNPi->GetDCA(trackP,bz,xthiss,xpp);


      fHistDCAdpdpi->Fill(dca_dp,dca_dpi);
      fHistDCApdppi->Fill(dca_dp,dca_ppi);
      fHistDCApidpip->Fill(dca_ppi,dca_dpi);

      if(fTriangularDCAtracks){
        if(dca_dpi > GetDCAcut(5,dca_dp)) continue;
        if(dca_ppi > GetDCAcut(4,dca_dp)) continue;
      } else{
        if(dca_dpi > fDCAdpi) continue;
        if(dca_ppi > fDCApip) continue;
      }

      fHistDCApiondeu->Fill(dca_dpi);
      fHistDCApionpro->Fill(dca_ppi);

      fTrkArray->AddAt(trackD,0);
      fTrkArray->AddAt(trackP,1);
      fTrkArray->AddAt(trackNPi,2);

      fVertexer->SetVtxStart(fPrimaryVertex);
      decayVtx = (AliESDVertex*)fVertexer->VertexForSelectedESDTracks(fTrkArray);

      SetConvertedAODVertices(fPrimaryVertex,decayVtx);

      fHistZDecayVtx->Fill(decayVtx->GetZ());
      fHistXDecayVtx->Fill(decayVtx->GetX());
      fHistYDecayVtx->Fill(decayVtx->GetY());
      decVt[0] = decayVtx->GetX();
      decVt[1] = decayVtx->GetY();
      decVt[2] = decayVtx->GetZ();

      dlh[0] = decayVtx->GetX() - fESDevent->GetPrimaryVertex()->GetX();
      dlh[1] = decayVtx->GetY() - fESDevent->GetPrimaryVertex()->GetY();
      dlh[2] = decayVtx->GetZ() - fESDevent->GetPrimaryVertex()->GetZ();

      decayLengthH3L = TMath::Sqrt((dlh[0]*dlh[0]) + (dlh[1]*dlh[1]) + (dlh[2]*dlh[2]));

      normalizedDecayL = fVtx2->DistanceToVertex(fVtx1)/fVtx2->ErrorDistanceToVertex(fVtx1);

      fHistDecayLengthH3L->Fill(decayLengthH3L);
      fHistNormalizedDecayL->Fill(normalizedDecayL);

      if(normalizedDecayL < fMinNormalizedDecL) {
        delete decayVtx;
        continue;
      }

      AliExternalTrackParam trkPi(*trackNPi);
      trkPi.PropagateToDCA(decayVtx, bz, 10,dcapi,dcapi_cov);
      fHistDCAXYpionvtx->Fill(dcapi[0]);
      fHistSigmaDcaXYpionvtx->Fill(TMath::Sqrt(dcapi_cov[0]));
      fHistDCAZpionvtx->Fill(dcapi[1]);
      fHistSigmaDcaZpionvtx->Fill(TMath::Sqrt(dcapi_cov[1]));
      fHistDCApionvtx->Fill(TMath::Sqrt((dcapi[0]*dcapi[0]) + (dcapi[1]*dcapi[1])));
      fHistSigmaDcapionvtx->Fill(ComputeSigma(dcapi,dcapi_cov));

      if(TMath::Abs(dcapi[0]) > fDCAPiSVxymax || TMath::Abs(dcapi[1]) > fDCAPiSVzmax) {
        delete decayVtx;
        continue;
      }

      AliExternalTrackParam trkP(*trackP);
      trkP.PropagateToDCA(decayVtx, bz, 10,dcap,dcapi_cov);
      fHistDCAXYprovtx->Fill(dcap[0]);
      fHistSigmaDcaXYprovtx->Fill(TMath::Sqrt(dcapi_cov[0]));
      fHistDCAZprovtx->Fill(dcap[1]);
      fHistSigmaDcaZprovtx->Fill(TMath::Sqrt(dcapi_cov[1]));
      fHistDCAprovtx->Fill(TMath::Sqrt((dcap[0]*dcap[0]) + (dcap[1]*dcap[1])));
      fHistSigmaDcaprovtx->Fill(ComputeSigma(dcap,dcapi_cov));

      if(TMath::Sqrt((dcap[0]*dcap[0])+(dcap[1]*dcap[1])) > fDCAProSVmax) {
        delete decayVtx;
        continue;
      }

      AliExternalTrackParam trkD(*trackD);
      trkD.PropagateToDCA(decayVtx, bz, 10,dcad,dcapi_cov);
      fHistDCAXYdeuvtx->Fill(dcad[0]);
      fHistSigmaDcaXYdeuvtx->Fill(TMath::Sqrt(dcapi_cov[0]));
      fHistDCAZdeuvtx->Fill(dcad[1]);
      fHistSigmaDcaZdeuvtx->Fill(TMath::Sqrt(dcapi_cov[1]));
      fHistDCAdeuvtx->Fill(TMath::Sqrt((dcad[0]*dcad[0]) + (dcad[1]*dcad[1])));
      fHistSigmaDcadeuvtx->Fill(ComputeSigma(dcad,dcapi_cov));

      if(TMath::Sqrt((dcad[0]*dcad[0])+(dcad[1]*dcad[1])) > fDCADeuSVmax) {
        delete decayVtx;
        continue;
      }

      delete decayVtx;

      posD.SetXYZM(trkD.Px(),trkD.Py(),trkD.Pz(),1.87561);
      posP.SetXYZM(trkP.Px(),trkP.Py(),trkP.Pz(),0.93827);
      negPi.SetXYZM(trkPi.Px(),trkPi.Py(),trkPi.Pz(),0.13957);


      Hypertriton=posD+posP+negPi;
      if(fRequireMassRange && Hypertriton.M()>fCutMassUp) continue;


      if(decayLengthH3L > fMaxDecayLength || decayLengthH3L < fMinDecayLength) continue;


      pTotHyper = Hypertriton.P();
      fHistPtHypertriton->Fill(Hypertriton.Pt());
      fHistPtProton->Fill(trackP->Pt());
      fHistPtDeuteron->Fill(trackD->Pt());
      fHistPtPion->Fill(trackNPi->Pt());

      if(Hypertriton.Pt() < fMinPtMother || Hypertriton.Pt() > fMaxPtMother) continue;

      h1.SetXYZ(dlh[0],dlh[1],dlh[2]);
      pointingAngleH = Hypertriton.Angle(h1);
      fHistCosPointingAngle->Fill(TMath::Cos(pointingAngleH));
      if(TMath::Cos(pointingAngleH) < fCosPointingAngle) continue;

      if(fSideBand == kTRUE && (Hypertriton.M() < 3.08 || Hypertriton.M() > 3.18)) continue;
      ctau = (Hypertriton.M()*decayLengthH3L)/pTotHyper;
      fHistLifetime->Fill(ctau);

      if(ctau < fMinLifeTime || ctau > fMaxLifeTime) continue;

      rapidity = Hypertriton.Rapidity();
      fHistHyperRapidity->Fill(rapidity);
      if(TMath::Abs(rapidity) > fRapidity) continue;

      //Angular correlation

      d1.SetXYZ(trkD.Px(),trkD.Py(),trkD.Pz());
      p1.SetXYZ(trkP.Px(),trkP.Py(),trkP.Pz());
      pi1.SetXYZ(trkPi.Px(),trkPi.Py(),trkPi.Pz());

      angle_dp = d1.Angle(p1);
      angle_dpi = d1.Angle(pi1);
      angle_ppi = p1.Angle(pi1);
      fHistAngle_deu_pro->Fill(angle_dp);
      fHistAngle_deu_pion->Fill(angle_dpi);
      fHistAngle_pro_pion->Fill(angle_ppi);

      fHistAngleCorr_dp_dpi->Fill(angle_dp,angle_dpi);
      fHistAngleCorr_dp_ppi->Fill(angle_dp,angle_ppi);
      fHistAngleCorr_ppi_dpi->Fill(angle_ppi,angle_dpi);

      if(angle_dp > fAngledp) continue;
      if(angle_dpi > fAngledpi) continue;

      fHistMassHyp_Lifetime->Fill(Hypertriton.M(),ctau);
      if(isMatter)	{ //
          fHistMassHypertriton->Fill(Hypertriton.M());
          fHistMassHyp_Lifetime_M->Fill(Hypertriton.M(),ctau);
          if(cent0) fHistMassHypertriton_Cent->Fill(Hypertriton.M());
          if(cent1) fHistMassHypertriton_SemiCent->Fill(Hypertriton.M());
       }
      if(!isMatter)	{ //
          fHistMassAntiHypertriton->Fill(Hypertriton.M());
          fHistMassHyp_Lifetime_A->Fill(Hypertriton.M(),ctau);
          if(cent0) fHistMassAntiHypertriton_Cent->Fill(Hypertriton.M());
          if(cent1) fHistMassAntiHypertriton_SemiCent->Fill(Hypertriton.M());
      }


/*  if(fMC){
  lD = trackD->GetLabel();
  lP = trackP->GetLabel();
  lPi = trackNPi->GetLabel();

  tparticleD = stack->Particle(TMath::Abs(lD));
  tparticleP = stack->Particle(TMath::Abs(lP));
  tparticlePi = stack->Particle(TMath::Abs(lPi));

  if((tparticleD->GetPdgCode() == 1000010020 && tparticleP->GetPdgCode() == 2212 && tparticlePi->GetPdgCode() == -211) || (tparticleD->GetPdgCode() == -1000010020 && tparticleP->GetPdgCode() == -2212 && tparticlePi->GetPdgCode() == 211)){
        labelM_deu = tparticleD->GetFirstMother();
        labelM_pro = tparticleP->GetFirstMother();
        labelM_pio = tparticlePi->GetFirstMother();

        TParticle *tparticleMotherD = stack->Particle(TMath::Abs(labelM_deu));
        TParticle *tparticleMotherP = stack->Particle(TMath::Abs(labelM_pro));
        TParticle *tparticleMotherPi = stack->Particle(TMath::Abs(labelM_pio));

    if(labelM_deu == labelM_pro && labelM_pro == labelM_pio){
          if(TMath::Abs(tparticleMotherP->GetPdgCode()) == 1010010030){
              if(tparticleD->GetUniqueID() == kPDecay && tparticleP->GetUniqueID() == kPDecay && tparticlePi->GetUniqueID() == kPDecay){

                  brotherHood = kTRUE;

                  fHistpionTPCclsMCt->Fill(trackNPi->GetTPCclusters(0));
                  fHistpTpionMCt->Fill(trackNPi->Pt());
                  fHistpTproMCt->Fill(trackP->Pt());
                  fHistpTdeuMCt->Fill(trackD->Pt());
                  fHistMompionMCt->Fill(trackNPi->P());
                  fHistMomproMCt->Fill(trackP->P());
                  fHistMomdeuMCt->Fill(trackD->P());
                  fHistDCAdeuproMCt->Fill(dca_dp);
                  fHistDCApiondeuMCt->Fill(dca_dpi);
                  fHistDCApionproMCt->Fill(dca_ppi);
                  fHistCorrDCApiprimaryMCt->Fill(piprim[0],piprim[1]);
                  fHistCorrDCApprimaryMCt->Fill(pprim[0],pprim[1]);
                  fHistCorrDCAdprimaryMCt->Fill(dprim[0],dprim[1]);
                  fHistDCApiprimaryMCt->Fill(dcapiprim);
                  fHistDCApprimaryMCt->Fill(dcapprim);
                  fHistDCAdprimaryMCt->Fill(dcadprim);
                  fHistDCAXYdeuvtxMCt->Fill(dcad[0]);
                  fHistDCAZdeuvtxMCt->Fill(dcad[1]);
                  fHistDCAXYprovtxMCt->Fill(dcap[0]);
                  fHistDCAZprovtxMCt->Fill(dcap[1]);
                  fHistDCAXYpionvtxMCt->Fill(dcapi[0]);
                  fHistDCAZpionvtxMCt->Fill(dcapi[1]);
                  fHistZDecayVtxMCt->Fill(decVt[2]);
                  fHistXDecayVtxMCt->Fill(decVt[0]);
                  fHistYDecayVtxMCt->Fill(decVt[1]);
                  fHistNormalizedDecayL_MCt->Fill(normalizedDecayL);
                  fHistLifetime_MCt->Fill(ctau);
                  fHistAngle_deu_pro_MCt->Fill(d1.Angle(p1));
                  fHistAngle_deu_pion_MCt->Fill(d1.Angle(pi1));
                  fHistAngle_pro_pion_MCt->Fill(p1.Angle(pi1));
                  fHistAngleCorr_dp_dpi_MCt->Fill(d1.Angle(p1),d1.Angle(pi1));
                  fHistAngleCorr_dp_ppi_MCt->Fill(d1.Angle(p1),p1.Angle(pi1));
                  fHistAngleCorr_ppi_dpi_MCt->Fill(p1.Angle(pi1),d1.Angle(pi1));
                  fHistHypertritonMomMCt->Fill(pTotHyper);
                  fHistHyperRapidityMCt->Fill(rapidity);
                  if(charge_pi<0) fHistMassHypertritonMCt->Fill(Hypertriton.M());
                  if(charge_pi>0) fHistMassAntiHypertritonMCt->Fill(Hypertriton.M());
          } //end of check kPDecay
        } //end of check pdgCodeHypertriton
      } //end of check labelM_deu, labelM_pro, labelM_pio
  } //end of check pdgCode daughters
} //end of if(MC)
*/

if(fFillTree){
  trackNPi->GetImpactParameters(piprim,piprimc);
  nsd = fPIDResponse->NumberOfSigmasTPC(trackD,AliPID::kDeuteron);
  nsp = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kProton);
  nspi = fPIDResponse->NumberOfSigmasTPC(trackNPi,AliPID::kPion);
  fTCentralityPerc = fCentralityClass;
  if(fRequireTOFPid){
    nsd_t = HasTOF(trackD, b_t) ? fPIDResponse->NumberOfSigmasTOF(trackD,AliPID::kDeuteron) : -999;
    nsp_t = HasTOF(trackP, b_t) ? fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kProton) : -999;
    nspi_t = HasTOF(trackNPi, b_t) ? fPIDResponse->NumberOfSigmasTOF(trackNPi,AliPID::kPion) : -999;

    if(nsd_t < 0)fTTOFnsigmadeu = TMath::Floor(nsd_t/0.25);
    else fTTOFnsigmadeu = TMath::Ceil(nsd_t/0.25);
    if(nsp_t < 0)fTTOFnsigmapro = TMath::Floor(nsp_t/0.25);
    else fTTOFnsigmapro = TMath::Ceil(nsp_t/0.25);
    if(nspi_t < 0)fTTOFnsigmapion = TMath::Floor(nspi_t/0.25);
    else fTTOFnsigmapion = TMath::Ceil(nspi_t/0.25);
  }
  //deuteron
  fTchi2NDFdeu = TMath::Floor(trackD->GetTPCchi2()/(0.5*trackD->GetTPCclusters(0)));
  fTPCclsdeu = trackD->GetTPCclusters(0);
  fITSdclsmap = trackD->GetITSClusterMap();
  fTpTdeu = TMath::Floor(trackD->Pt()/0.000107692);
  fTpdeu = TMath::Floor(trackD->P()/0.000107692);
  if(nsd < 0)fTTPCnsigmadeu = TMath::Floor(nsd/0.25);
  else fTTPCnsigmadeu = TMath::Ceil(nsd/0.25);
  //proton
  fTchi2NDFpro = TMath::Floor(trackP->GetTPCchi2()/(0.5*trackP->GetTPCclusters(0)));
  fTPCclspro = trackP->GetTPCclusters(0);
  fITSpclsmap = trackP->GetITSClusterMap();
  fTpTpro = TMath::Floor(trackP->Pt()/0.000107692);
  fTppro = TMath::Floor(trackP->P()/0.000107692);
  if(nsp < 0)fTTPCnsigmapro = TMath::Floor(nsp/0.25);
  else fTTPCnsigmapro = TMath::Ceil(nsp/0.25);
  //pion
  fTchi2NDFpion = TMath::Floor(trackNPi->GetTPCchi2()/(0.5*trackNPi->GetTPCclusters(0)));
  fTPCclspion = trackNPi->GetTPCclusters(0);
  fITSpiclsmap = trackNPi->GetITSClusterMap();
  fTpTpion = TMath::Floor(trackNPi->Pt()/0.000107692);
  fTppion = TMath::Floor(trackNPi->P()/0.000107692);
  if(nspi < 0)fTTPCnsigmapion = TMath::Floor(nspi/0.25);
  else fTTPCnsigmapion = TMath::Ceil(nspi/0.25);
  if(TMath::Abs(piprim[0])<819)  fTDCAXYpioprvtx = TMath::Ceil(piprim[0]/0.025);
  else fTDCAXYpioprvtx = 32765;
  if(TMath::Abs(piprim[1])<819)  fTDCAZpioprvtx = TMath::Ceil(piprim[1]/0.025);
  else fTDCAZpioprvtx = 32765;

  //triplets
  fTDCAdp = TMath::Ceil(dca_dp/0.001);
  fTDCAdpi = TMath::Ceil(dca_dpi/0.001);
  fTDCAppi = TMath::Ceil(dca_ppi/0.001);
  if(dcad[0]<0)fTDCAXYdvtx = TMath::Floor(dcad[0]/0.001);
  else fTDCAXYdvtx = TMath::Ceil(dcad[0]/0.001);

  if(dcad[1]<0)fTDCAZdvtx = TMath::Floor(dcad[1]/0.001);
  else fTDCAZdvtx = TMath::Ceil(dcad[1]/0.001);

  if(dcap[0]<0)fTDCAXYpvtx = TMath::Floor(dcap[0]/0.001);
  else fTDCAXYpvtx = TMath::Ceil(dcap[0]/0.001);

  if(dcap[1]<0)fTDCAZpvtx = TMath::Floor(dcap[1]/0.001);
  else fTDCAZpvtx = TMath::Ceil(dcap[1]/0.001);

  if(dcapi[0]<0)fTDCAXYpivtx = TMath::Floor(dcapi[0]/0.001);
  else fTDCAXYpivtx = TMath::Ceil(dcapi[0]/0.001);

  if(dcapi[1]<0)fTDCAZpivtx = TMath::Floor(dcapi[1]/0.001);
  else fTDCAZpivtx = TMath::Ceil(dcapi[1]/0.001);

  fTAngle_dp = TMath::Ceil(angle_dp);
  fTAngle_dpi = TMath::Ceil(angle_dpi);
  fTAngle_ppi = TMath::Ceil(angle_ppi);

  fTRapidity = TMath::Ceil(TMath::Abs(rapidity)/0.1);
  fTDecayLengthProper = ctau;
  fTDecayLengthNorm = TMath::Floor(normalizedDecayL);
  fTCosPA = TMath::Cos(pointingAngleH);
  fTTransverseMom = TMath::Floor(Hypertriton.Pt()/0.001);

  /*if(fMC){
      fTMCtruth = brotherHood;
      labelM_deu = tparticleD->GetFirstMother();
      labelM_pro = tparticleP->GetFirstMother();
      labelM_pio = tparticlePi->GetFirstMother();
      TParticle *tpmomD = stack->Particle(TMath::Abs(labelM_deu));
      TParticle *tpmomP = stack->Particle(TMath::Abs(labelM_pro));
      TParticle *tpmomPi = stack->Particle(TMath::Abs(labelM_pio));
      fTpdeu_gen_X = tparticleD->Px();
      fTpdeu_gen_Y = tparticleD->Py();
      fTpdeu_gen_Z = tparticleD->Pz();
      fTppro_gen_X = tparticleP->Px();
      fTppro_gen_Y = tparticleP->Py();
      fTppro_gen_Z = tparticleP->Pz();
      fTppio_gen_X = tparticlePi->Px();
      fTppio_gen_Y = tparticlePi->Py();
      fTppio_gen_Z = tparticlePi->Pz();
      fTpdgDeu = tparticleD->GetPdgCode();
      fTpdgPro = tparticleP->GetPdgCode();
      fTpdgPion = tparticlePi->GetPdgCode();
      fTmomidD = tparticleD->GetFirstMother();
      fTmomidP = tparticleP->GetFirstMother();
      fTmomidPi = tparticlePi->GetFirstMother();
      fTpdgmomD = tpmomD->GetPdgCode();
      fTpdgmomP = tpmomP->GetPdgCode();
      fTpdgmomPi = tpmomPi->GetPdgCode();
      fTuniqID_deu = tparticleD->GetUniqueID();
      fTuniqID_pro = tparticleP->GetUniqueID();
      fTuniqID_pion = tparticlePi->GetUniqueID();
  }
*/
  if(isMatter) fTInvariantMass = -Hypertriton.M();
  else fTInvariantMass = Hypertriton.M();

  fTTree->Fill();
  PostData(2,fTTree);
     } //end of Fill Tree
    } // end of candidate pion loop
  } // end of candidate proton loop
}// end of candidate deuteron loop

}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserCreateOutputObjects(){
  // Create a TList with Histograms
  // Called once
  //printf("**************************************************\n");
  //printf("AliAnalysisTaskHypertriton3::CreateOutputObjects()\n");
  //printf("**************************************************\n");

  fVertexer = new AliVertexerTracks();
  fTrkArray = new TObjArray(3);
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

  fHistITSclusters = new TH1F("fHistITSclusters","fHistITSclusters; nITScls; counts",7,-0.5,6.5);

  //TPC
  fHistTPCpid = new TH2F("fHistTPCpid", "dE/dx for all tracks; p_{TPC} (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistTPCpid->SetOption("scat");
  fHistTPCpid->SetMarkerStyle(kFullCircle);

  //Hypertriton prongs
  fHistTPCdeusignal = new TH2F("fHistTPCdeusignal", "dE/dx after d(#bar{d}) PID; p_{TPC}/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistTPCdeusignal->SetOption("scat");
  fHistTPCdeusignal->SetMarkerStyle(kFullCircle);

  fHistTPCprosignal = new TH2F("fHistTPCprosignal", "dE/dx after p(#bar{p}) PID; p_{TPC}/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistTPCprosignal->SetOption("scat");
  fHistTPCprosignal->SetMarkerStyle(kFullCircle);

  fHistTPCpionsignal= new TH2F("fHistTPCpionsignal", "dE/dx after  #pi^{+}(#pi^{-}) PID; p_{TPC}/Z (GeV/c); TPC signal", 1000, -10., 10, 300, 0, 2100);
  fHistTPCpionsignal->SetOption("scat");
  fHistTPCpionsignal->SetMarkerStyle(kFullCircle);


  //TOF

  fHistTOFsignal = new TH2F("fHistTOFsignal","TOF signal; p (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFdeusignal = new TH2F("fHistTOFdeusignal","#beta vs p - deuteron; p (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFprosignal = new TH2F("fHistTOFprosignal","#beta vs p - proton; p (GeV/c); #beta",400,0.,4.,400,0.,1.1);

  fHistTOFpionsignal = new TH2F("fHistTOFpionsignal","#beta vs p - pion; p (GeV/c); #beta",400,0.,4,400,0.,1.1);

  //fHistTOFdeumass = new TH1F("fHistTOFdeumass","deuteron mass distribution - TOF; mass (GeV/c^{2}); entries",400,0.8,2.8);

  //fHistTOFpromass = new TH1F("fHistTOFpromass","proton mass distribution - TOF; mass (GeV/c^{2}); entries",200,0.5,1.5);


  fHistpionTPCcls = new TH1F("fHistpionTPCcls","#pi^{-} TPC clusters; TPC clusters; entries",201,-0.5,200.5);
  //fHistCorrDCAdprimary = new TH2F("fHistCorrDCAdprimary","DCA_{PV,xy} vs DCA_{PV,z} - deuteron; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
  //fHistCorrDCApprimary = new TH2F("fHistCorrDCApprimary","DCA_{PV,xy} vs DCA_{PV,z} - proton; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
  //fHistCorrDCApiprimary = new TH2F("fHistCorrDCApiprimary","DCA_{PV,xy} vs DCA_{PV,z} - pion; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
  fHistDCApiprimary = new TH1F("fHistDCApiprimary","DCA pion-primary vertex; DCA (cm); entries",800,0.f,20.f);
  fHistDCAXYpiprimary = new TH1F("fHistDCAXYpiprimary","DCAxy pion-primary vertex; DCA (cm); entries",800,-10.f,10.f);
  fHistDCAZpiprimary = new TH1F("fHistDCAZpiprimary","DCAz pion-primary vertex; DCA (cm); entries",800,-10.f,10.f);

  fHistDCApprimary = new TH1F("fHistDCApprimary","DCA proton-primary vertex; DCA (cm); entries",800,0.f,20.f);
  fHistDCAXYpprimary = new TH1F("fHistDCAXYpprimary","DCAxy proton-primary vertex; DCA (cm); entries",800,-10.f,10.f);
  fHistDCAZpprimary = new TH1F("fHistDCAZpprimary","DCAz proton-primary vertex; DCA (cm); entries",800,-10.f,10.f);

  fHistDCAdprimary = new TH1F("fHistDCAdprimary","DCA deuteron-primary vertex; DCA (cm); entries",800,0.f,20.f);
  fHistDCAXYdprimary = new TH1F("fHistDCAXYdprimary","DCAxy deuteron-primary vertex; DCA (cm); entries",800,-10.f,10.f);
  fHistDCAZdprimary = new TH1F("fHistDCAZdprimary","DCAz deuteron-primary vertex; DCA (cm); entries",800,-10.f,10.f);

  //DCA prongs
  fHistDCAdeupro = new TH1F("fHistDCAdeupro","DCA d-p tracks;d-p DCA (cm);entries",550,-0.5,5.0);
  fHistDCApiondeu = new TH1F("fHistDCApiondeu","DCA #pi^{-}-d tracks; #pi^{-}-d DCA (cm); entries",550,-0.5,5.0);
  fHistDCApionpro = new TH1F("fHistDCApionpro","DCA #pi^{-}-p tracks; #pi^{-}-p DCA (cm); entries",550,-0.5,5.0);
  fHistDCAdpdpi = new TH2F("fHistDCAdpdpi","DCA deu-pro vs DCA deu-pion;DCA_{dp} (cm); DCA_{d#pi} (cm)",100,0.,1.,100,0.,1.);
  fHistDCApdppi = new TH2F("fHistDCApdppi","DCA deu-pro vs DCA pion-pro; DCA_{dp} (cm); DCA_{#pip} (cm)",100,0.,1.,100,0.,1.);
  fHistDCApidpip = new TH2F("fHistDCApidpip","DCA pion-pro vs DCA pion-deu;DCA_{#pip} (cm); DCA_{#pid} (cm)",100,0.,1.,100,0.,1.);

  //Decay histo
  fHistZDecayVtx = new TH1F("fHistZDecayVtx","decay vertex - z coordinate; z_{decay vtx} (cm); entries",800,-20.f,20.f);
  fHistXDecayVtx = new TH1F("fHistXDecayVtx","decay vertex - x coordinate; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistYDecayVtx = new TH1F("fHistYDecayVtx","decay vertex - y coordinate; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
  fHistDCApionvtx = new TH1F("fHistDCApionvtx","DCA candidate #pi^{-}-decay vertex; DCA (cm); entries",1500,0.,15.);
  fHistDCAXYpionvtx = new TH1F("fHistDCAXYpionvtx","DCA candidate #pi^{-}-decay vertex - xy coordinate; DCA_{xy} (cm); entries",1000,-5.,5.);
  fHistDCAZpionvtx = new TH1F("fHistDCAZpionvtx","DCA candidate #pi^{-}-decay vertex - z coordinate; DCA_{z} (cm); entries",2000,-10.,10.);
  fHistSigmaDcaXYpionvtx = new TH1F("fHistSigmaDcaXYpionvtx","#sigma_{DCA, xy} - candiate #pi-decay vtx; #sigma_{DCA, xy} (cm)",202,-0.05,5);
  fHistSigmaDcaZpionvtx = new TH1F("fHistSigmaDcaZpionvtx","#sigma_{DCA, z} - candiate #pi-decay vtx; #sigma_{DCA, z} (cm)",202,-0.05,5);
  fHistSigmaDcapionvtx = new TH1F("fHistSigmaDcapionvtx","#sigma_{DCA} - candidate #pi-decay vtx; #sigma_{DCA} (cm)",202,-0.05,5);
  fHistDCAprovtx = new TH1F("fHistDCAprovtx","DCA candidate p-decay vertex; DCA (cm); entries",1000,0.,15.);
  fHistDCAXYprovtx = new TH1F("fHistDCAXYprovtx","DCA candidate p-decay vertex - xy coordinate; DCA_{xy} (cm); entries",1000,-5.,5.);
  fHistDCAZprovtx = new TH1F("fHistDCAZprovtx","DCA candidate p-decay vertex - z coordinate; DCA_{z} (cm); entries",2000,-10.,10.);
  fHistSigmaDcaXYprovtx = new TH1F("fHistSigmaDcaXYprovtx","#sigma_{DCA, xy} - candiate p-decay vtx; #sigma_{DCA, xy} (cm)",202,-0.05,5);
  fHistSigmaDcaZprovtx = new TH1F("fHistSigmaDcaZprovtx","#sigma_{DCA, z} - candiate p-decay vtx; #sigma_{DCA, z} (cm)",202,-0.05,5);
  fHistSigmaDcaprovtx = new TH1F("fHistSigmaDcaprovtx","#sigma_{DCA} - candidate p-decay vtx; #sigma_{DCA} (cm)",202,-0.05,5);
  fHistDCAdeuvtx = new TH1F("fHistDCAdeuvtx","DCA candidate d-decay vertex; DCA (cm); entries",1500,0.,15.);
  fHistDCAXYdeuvtx = new TH1F("fHistDCAXYdeuvtx","DCA candidate d-decay vertex - xy coordinate; DCA_{xy} (cm); entries",1000,-5.,5.);
  fHistDCAZdeuvtx = new TH1F("fHistDCAZdeuvtx","DCA candidate d-decay vertex - z coordinate; DCA_{z} (cm); entries",2000,-10.,10.);
  fHistSigmaDcaXYdeuvtx = new TH1F("fHistSigmaDcaXYdeuvtx","#sigma_{DCA, xy} - candiate d-decay vtx; #sigma_{DCA, xy} (cm)",202,-0.05,5);
  fHistSigmaDcaZdeuvtx = new TH1F("fHistSigmaDcaZdeuvtx","#sigma_{DCA, z} - candiate d-decay vtx; #sigma_{DCA, z} (cm)",202,-0.05,5);
  fHistSigmaDcadeuvtx = new TH1F("fHistSigmaDcadeuvtx","#sigma_{DCA} - candidate d-decay vtx; #sigma_{DCA} (cm)",202,-0.05,5);
  fHistDecayLengthH3L = new TH1F("fHistDecayLengthH3L","decay length ^{3}H_{#Lambda}; decay length (cm); entries",400,0.,400.);
  fHistNormalizedDecayL = new TH1F("fHistNormalizedDecayL","normalized decay length; decL/#sigma_{dL}; entries",400,0.,100.);
  fHistLifetime = new TH1F("fHistLifetime","ct ^{3}H_{#Lambda}; ct(cm); entries",400,0.,40.);

  fHistAngle_deu_pro = new TH1F("fHistAngle_deu_pro","Angle between d and p; #alpha_{d_p} (rad); entries/(0.03 rad)",100,0.,TMath::Pi());
  fHistAngle_deu_pion = new TH1F("fHistAngle_deu_pion","Angle between d and #pi; #beta_{d_#pi} (rad); entries/(0.03 rad)",100,0.,TMath::Pi());
  fHistAngle_pro_pion = new TH1F("fHistAngle_pro_pion","Angle between p and #pi; #gamma_{p_#pi} (rad);entries/(0.03 rad)",100,0.,TMath::Pi());
  fHistAngleCorr_dp_dpi = new TH2F("fHistAngleCorr_dp_dpi","Correlation: #alpha_{d_p} vs #beta_{d_#pi};#alpha_{d_p};#beta_{d_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());
  fHistAngleCorr_dp_ppi = new TH2F("fHistAngleCorr_dp_ppi","Correlation: #alpha_{d_p} vs #gamma_{p_#pi};#alpha_{d_p};#gamma_{p_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());
  fHistAngleCorr_ppi_dpi = new TH2F("fHistAngleCorr_ppi_dpi","Correlation: #gamma_{p_#pi} vs #beta_{d_#pi};#gamma_{p_#pi};#beta_{d_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());
  fHistHyperRapidity = new TH1F("fHistHyperRapidity","rapidity distribution of ^{3}H_{#Lambda}; y; entries",400,-2.f,2.f);

  fHistCosPointingAngle= new TH1F("fHistCosPointingAngle", "Cos pointing angle distribution; cos point angle; entries", 220, -1.1, 1.1);
  fHistPtPion = new TH1F("fHistPtPion","#pi p_{T}; p_{T} (GeV/c); entries",100.,0.,10.);
  fHistPtDeuteron = new TH1F("fHistPtDeuteron","d p_{T}; p_{T} (GeV/c); entries",100.,0.,10.);
  fHistPtProton = new TH1F("fHistPtProton","p p_{T}; p_{T} (GeV/c); entries",100.,0.,10.);
  fHistPtHypertriton = new TH1F("fHistPtHypertriton","candidate Hypertriton - p_{T} distribution; p_{T} (GeV/c); entries",100.,0.,10.);
  if(fMinvSignal){
      fHistMassHypertriton = new TH1F("fHistMassHypertriton", "Invariant mass distribution d+p+#pi^{-};invariant mass d+p+#pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton = new TH1F("fHistMassAntiHypertriton", "Invariant mass distribution #bar{d} + #bar{p} + #pi^{+};invariant mass #bar{d} + #bar{p} + #pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassHypertriton_Cent = new TH1F("fHistMassHypertriton_Cent","Invariant mass distribution d+p+#pi^{-};invariant mass d+p+#pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton_Cent = new TH1F("fHistMassAntiHypertriton_Cent", "Invariant mass distribution #bar{d} + #bar{p} + #pi^{+};invariant mass #bar{d} + #bar{p} + #pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassHypertriton_SemiCent = new TH1F("fHistMassHypertriton_SemiCent","Invariant mass distribution d+p+#pi^{-};invariant mass d+p+#pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton_SemiCent = new TH1F("fHistMassAntiHypertriton_SemiCent", "Invariant mass distribution #bar{d} + #bar{p} + #pi^{+};invariant mass #bar{d} + #bar{p} + #pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassHyp_Lifetime = new TH2F("fHistMassHyp_Lifetime","Inv mass vs c#it{t}; mass (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);
      fHistMassHyp_Lifetime_M = new TH2F("fHistMassHyp_Lifetime_M","Inv mass vs c#it{t}; d+p+#pi^{-} (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);
      fHistMassHyp_Lifetime_A = new TH2F("fHistMassHyp_Lifetime_A","Inv mass vs c#it{t}; #bar{d} + #bar{p} + #pi^{+} (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);
  }
  if(fMinvLikeSign){
      fHistMassHypertriton = new TH1F("fHistMassHypertriton_LS", "Invariant mass distribution - Like Sign;invariant mass d+p+#pi^{+} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassAntiHypertriton = new TH1F("fHistMassAntiHypertriton_LS", "Invariant mass distribution - Like Sign;invariant mass #bar{d} + #bar{p} + #pi^{-} (GeV/c^{2}); entries ", 500, 2.9, 3.4);
      fHistMassHypertriton_Cent = new TH1F("fHistMassHypertriton_LS_Cent","Invarian mass spectrum - Like Sign; inv mass d+p+#pi^{+} (GeV/c^{2});entries",500,2.9,3.4);
      fHistMassAntiHypertriton_Cent = new TH1F("fHistMassAntiHypertriton_LS_Cent","Invarian mass spectrum - Like Sign; inv mass #bar{d}+#bar{p}+#pi^{-} (GeV/c^{2});entries",500,2.9,3.4);
      fHistMassHypertriton_SemiCent = new TH1F("fHistMassHypertriton_LS_SemiCent","Invarian mass spectrum - Like Sign; inv mass d+p+#pi^{+} (GeV/c^{2});entries",500,2.9,3.4);
      fHistMassAntiHypertriton_SemiCent = new TH1F("fHistMassAntiHypertriton_LS_SemiCent","Invarian mass spectrum - Like Sign; inv mass #bar{d}+#bar{p}+#pi^{-} (GeV/c^{2});entries",500,2.9,3.4);
      fHistMassHyp_Lifetime = new TH2F("fHistMassHyp_Lifetime_LS","Inv mass vs c#it{t}; mass (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);
      fHistMassHyp_Lifetime_M = new TH2F("fHistMassHyp_Lifetime_LS_M","Inv mass vs c#it{t}; d+p+#pi^{+} (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);
      fHistMassHyp_Lifetime_A = new TH2F("fHistMassHyp_Lifetime_LS_A","Inv mass vs c#it{t}; #bar{d} + #bar{p} + #pi^{-} (GeV/#it{c}^{2}); c#it{t} (cm)",500,2.9,3.4,80,0.,40.);

  }

  if(fMC){
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

    fHistpionTPCclsMCt = new TH1F("fHistpionTPCclsMCt","#pi^{-} TPC clusters - MCtruth; TPC clusters; entries",201,-0.5,200.5);
    fHistpTpionMCt = new TH1F("fHistpTpionMCt","pion p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistpTproMCt = new TH1F("fHistpTproMCt","proton p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistpTdeuMCt = new TH1F("fHistpTdeuMCt","deuteron p_{T} distribution; p_{T} (GeV/c);entries",800,0.,8.);
    fHistMompionMCt = new TH1F("fHistMompionMCt","pion p distribution; p (GeV/c);entries",1000,0.,10.);
    fHistMomproMCt = new TH1F("fHistMomproMCt","proton p distribution; p (GeV/c);entries",1000,0.,10.);
    fHistMomdeuMCt = new TH1F("fHistMomdeuMCt","deuteron p distribution; p (GeV/c);entries",1000,0.,10.);
    fHistCorrDCAdprimaryMCt = new TH2F("fHistCorrDCAdprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - deuteron MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
    fHistCorrDCApprimaryMCt = new TH2F("fHistCorrDCApprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - proton MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
    fHistCorrDCApiprimaryMCt = new TH2F("fHistCorrDCApiprimaryMCt","DCA_{PV,xy} vs DCA_{PV,z} - pion MCtruth; DCA_{xy} (cm); DCA_{z} (cm)",320,-20.f,20.f,320,-20.f,20.f);
    fHistDCApiprimaryMCt = new TH1F("fHistDCApiprimaryMCt","DCA pion-primary vertex; DCA (cm); entries",3200,0.f,80.f);
    fHistDCApprimaryMCt = new TH1F("fHistDCApprimaryMCt","DCA proton-primary vertex; DCA (cm); entries",3200,0.f,80.f);
    fHistDCAdprimaryMCt = new TH1F("fHistDCAdprimaryMCt","DCA deuteron-primary vertex; DCA (cm); entries",3200,0.f,80.f);

    fHistDCAdeuproMCt = new TH1F("fHistDCAdeuproMCt","DCA d-p tracks - MCtruth;d-p DCA (cm);entries",250,-1.0,1.0);
    fHistDCApiondeuMCt = new TH1F("fHistDCApiondeuMCt","DCA #pi^{-}-d tracks - MCtruth;#pi^{-}-d DCA (cm);entries",250,-1.0,1.0);
    fHistDCApionproMCt = new TH1F("fHistDCApionproMCt","DCA #pi^{-}-p tracks - MCtruth;#pi^{-}-p DCA (cm);entries",250,-1.0,1.0);

    fHistZDecayVtxMCt = new TH1F("fHistZDecayVtxMCt","decay vertex - z coordinate - MCtruth; z_{decay vtx} (cm); entries",800,-20.f,20.f);
    fHistXDecayVtxMCt = new TH1F("fHistXDecayVtxMCt","decay vertex - x coordinate - MCtruth; x_{decay vtx} (cm); entries",8000,-20.f,20.f);
    fHistYDecayVtxMCt = new TH1F("fHistYDecayVtxMCt","decay vertex - y coordinate - MCtruth; y_{decay vtx} (cm); entries",8000,-20.f,20.f);
    fHistDCAXYdeuvtxMCt = new TH1F("fHistDCAXYdeuvtxMCt","DCA candidate d-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZdeuvtxMCt = new TH1F("fHistDCAZdeuvtxMCt","DCA candidate d-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);
    fHistDCAXYprovtxMCt = new TH1F("fHistDCAXYprovtxMCt","DCA candidate p-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZprovtxMCt = new TH1F("fHistDCAZprovtxMCt","DCA candidate p-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);
    fHistDCAXYpionvtxMCt = new TH1F("fHistDCAXYpionvtxMCt","DCA candidate #pi^{-}-decay vertex - xy coordinate - MCtruth; DCA_{xy} (cm); entries",200,-5.,5.);
    fHistDCAZpionvtxMCt = new TH1F("fHistDCAZpionvtxMCt","DCA candidate #pi^{-}-decay vertex - z coordinate - MCtruth; DCA_{z} (cm); entries",200,-10.,10.);

    fHistNormalizedDecayL_MCt = new TH1F("fHistNormalizedDecayL_MCt","normalized decay length - MCtruth; decL/#sigma_{dL}; entries",400,0.,100.);
    fHistLifetime_MCt = new TH1F("fHistLifetime_MCt","ct ^{3}H_{#Lambda} - MCtruth; ct(cm); entries",400,0.,400.);

    fHistAngle_deu_pro_MCt = new TH1F("fHistAngle_deu_pro_MCt","Angle between d and p - MCtruth; #alpha_{d_p} (rad) - MCtruth; entries/(0.03 rad)",100,0.,TMath::Pi());
    fHistAngle_deu_pion_MCt = new TH1F("fHistAngle_deu_pion_MCt","Angle between d and #pi - MCtruth; #beta_{d_#pi} (rad); entries/(0.03 rad)",100,0.,TMath::Pi());
    fHistAngle_pro_pion_MCt = new TH1F("fHistAngle_pro_pion_MCt","Angle between p and #pi - MCtruth; #gamma_{p_#pi} (rad);entries/(0.03 rad)",100,0.,TMath::Pi());

    fHistAngleCorr_dp_dpi_MCt = new TH2F("fHistAngleCorr_dp_dpi_MCt","Correlation: #alpha_{d_p} vs #beta_{d_#pi} - MCtruth;#alpha_{d_p};#beta_{d_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());
    fHistAngleCorr_dp_ppi_MCt = new TH2F("fHistAngleCorr_dp_ppi_MCt","Correlation: #alpha_{d_p} vs #gamma_{p_#pi} - MCtruth;#alpha_{d_p};#gamma_{p_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());
    fHistAngleCorr_ppi_dpi_MCt = new TH2F("fHistAngleCorr_ppi_dpi_MCt","Correlation: #gamma_{p_#pi} vs #beta_{d_#pi} - MCtruth;#gamma_{p_#pi};#beta_{d_#pi}",100,0.,TMath::Pi(),100,0.,TMath::Pi());

    fHistHypertritonMomMCt = new TH1F("fHistHypertritonMomMCt","^{3}H_{#Lambda} momentum - MCtruth; p_{^{3}H_{#Lambda}} (GeV/c); entries/0.01",1100.,0.,11.);

    fHistHyperRapidityMCt =  new TH1F("fHistHyperRapidityMCt","rapidity distribution of ^{3}H_{#Lambda} - MC truth; y; entries",400,-2.f,2.f);

    fHistMassHypertritonMCt = new TH1F("fHistMassHypertritonMCt","^{3}H_{#Lambda} invariant mass - MCtruth; invariant mass d+p+#pi^{-} (GeV/c^{2}); entries", 400,2.9,3.1);
    fHistMassAntiHypertritonMCt = new TH1F("fHistMassAntiHypertritonMCt","#bar{^{3}H_{#Lambda}} invariant mass - MCtruth; invariant mass #bar{d}+#bar{p}+#pi^{+} (GeV/c^{2}); entries", 400,2.9,3.1);
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
  fOutput->Add(fHistITSclusters);
  fOutput->Add(fHistTPCpid);
  fOutput->Add(fHistTPCdeusignal);
  fOutput->Add(fHistTPCprosignal);
  fOutput->Add(fHistTPCpionsignal);
  fOutput->Add(fHistTOFsignal);
  fOutput->Add(fHistTOFdeusignal);
  fOutput->Add(fHistTOFprosignal);
  fOutput->Add(fHistTOFpionsignal);
  //fOutput->Add(fHistTOFdeumass);
  //fOutput->Add(fHistTOFpromass);
  fOutput->Add(fHistpionTPCcls);
  //fOutput->Add(fHistCorrDCAdprimary);
  //fOutput->Add(fHistCorrDCApprimary);
  //fOutput->Add(fHistCorrDCApiprimary);
  fOutput->Add(fHistDCApiprimary);
  fOutput->Add(fHistDCAXYpiprimary);
  fOutput->Add(fHistDCAZpiprimary);
  fOutput->Add(fHistDCApprimary);
  fOutput->Add(fHistDCAXYpprimary);
  fOutput->Add(fHistDCAZpprimary);
  fOutput->Add(fHistDCAdprimary);
  fOutput->Add(fHistDCAXYdprimary);
  fOutput->Add(fHistDCAZdprimary);
  fOutput->Add(fHistDCAdeupro);
  fOutput->Add(fHistDCApiondeu);
  fOutput->Add(fHistDCApionpro);
  fOutput->Add(fHistDCAdpdpi);
  fOutput->Add(fHistDCApdppi);
  fOutput->Add(fHistDCApidpip);
  fOutput->Add(fHistZDecayVtx);
  fOutput->Add(fHistXDecayVtx);
  fOutput->Add(fHistYDecayVtx);
  fOutput->Add(fHistDecayLengthH3L);
  fOutput->Add(fHistNormalizedDecayL);
  fOutput->Add(fHistDCAXYpionvtx);
  fOutput->Add(fHistDCAZpionvtx);
  fOutput->Add(fHistDCApionvtx);
  fOutput->Add(fHistSigmaDcaXYpionvtx);
  fOutput->Add(fHistSigmaDcaZpionvtx);
  fOutput->Add(fHistSigmaDcapionvtx);
  fOutput->Add(fHistDCAXYprovtx);
  fOutput->Add(fHistDCAZprovtx);
  fOutput->Add(fHistDCAprovtx);
  fOutput->Add(fHistSigmaDcaXYprovtx);
  fOutput->Add(fHistSigmaDcaZprovtx);
  fOutput->Add(fHistSigmaDcaprovtx);
  fOutput->Add(fHistDCAXYdeuvtx);
  fOutput->Add(fHistDCAZdeuvtx);
  fOutput->Add(fHistDCAdeuvtx);
  fOutput->Add(fHistSigmaDcaXYdeuvtx);
  fOutput->Add(fHistSigmaDcaZdeuvtx);
  fOutput->Add(fHistSigmaDcadeuvtx);
  fOutput->Add(fHistPtPion);
  fOutput->Add(fHistPtDeuteron);
  fOutput->Add(fHistPtProton);
  fOutput->Add(fHistPtHypertriton);
  fOutput->Add(fHistCosPointingAngle);
  fOutput->Add(fHistLifetime);
  fOutput->Add(fHistHyperRapidity);
  fOutput->Add(fHistAngle_deu_pro);
  fOutput->Add(fHistAngle_deu_pion);
  fOutput->Add(fHistAngle_pro_pion);
  fOutput->Add(fHistAngleCorr_dp_dpi);
  fOutput->Add(fHistAngleCorr_dp_ppi);
  fOutput->Add(fHistAngleCorr_ppi_dpi);
  fOutput->Add(fHistMassHypertriton);
  fOutput->Add(fHistMassAntiHypertriton);
  fOutput->Add(fHistMassHypertriton_Cent);
  fOutput->Add(fHistMassAntiHypertriton_Cent);
  fOutput->Add(fHistMassHypertriton_SemiCent);
  fOutput->Add(fHistMassAntiHypertriton_SemiCent);
  fOutput->Add(fHistMassHyp_Lifetime);
  fOutput->Add(fHistMassHyp_Lifetime_M);
  fOutput->Add(fHistMassHyp_Lifetime_A);


  if(fMC){
    fOutput->Add(fHistParticle);
    fOutput->Add(fHistpionTPCclsMCt);
    fOutput->Add(fHistpTpionMCt);
    fOutput->Add(fHistpTproMCt);
    fOutput->Add(fHistpTdeuMCt);
    fOutput->Add(fHistMompionMCt);
    fOutput->Add(fHistMomproMCt);
    fOutput->Add(fHistMomdeuMCt);
    fOutput->Add(fHistCorrDCAdprimaryMCt);
    fOutput->Add(fHistCorrDCApprimaryMCt);
    fOutput->Add(fHistCorrDCApiprimaryMCt);
    fOutput->Add(fHistDCApiprimaryMCt);
    fOutput->Add(fHistDCApprimaryMCt);
    fOutput->Add(fHistDCAdprimaryMCt);
    fOutput->Add(fHistDCAdeuproMCt);
    fOutput->Add(fHistDCApiondeuMCt);
    fOutput->Add(fHistDCApionproMCt);
    fOutput->Add(fHistZDecayVtxMCt);
    fOutput->Add(fHistXDecayVtxMCt);
    fOutput->Add(fHistYDecayVtxMCt);
    fOutput->Add(fHistDCAXYdeuvtxMCt);
    fOutput->Add(fHistDCAZdeuvtxMCt);
    fOutput->Add(fHistDCAXYprovtxMCt);
    fOutput->Add(fHistDCAZprovtxMCt);
    fOutput->Add(fHistDCAXYpionvtxMCt);
    fOutput->Add(fHistDCAZpionvtxMCt);
    fOutput->Add(fHistNormalizedDecayL_MCt);
    fOutput->Add(fHistLifetime_MCt);
    fOutput->Add(fHistAngle_deu_pro_MCt);
    fOutput->Add(fHistAngle_deu_pion_MCt);
    fOutput->Add(fHistAngle_pro_pion_MCt);
    fOutput->Add(fHistAngleCorr_dp_dpi_MCt);
    fOutput->Add(fHistAngleCorr_dp_ppi_MCt);
    fOutput->Add(fHistAngleCorr_ppi_dpi_MCt);
    fOutput->Add(fHistHypertritonMomMCt);
    fOutput->Add(fHistHyperRapidityMCt);
    fOutput->Add(fHistMassHypertritonMCt);
    fOutput->Add(fHistMassAntiHypertritonMCt);
  }

  // Post output data.
  PostData(1,fOutput);
  //printf("**************************************************\n");
  //printf("end of fOutput\n");
  //printf("**************************************************\n");

  if(fFillTree){
  OpenFile(2);
  fTTree = new TTree("hypertriton","hypertriton candidates");
  fTTree->Branch("CentralityPerc",&fTCentralityPerc,"CentralityPerc/b");
  fTTree->Branch("Chi2NDFdeu",&fTchi2NDFdeu,"Chi2NDFdeu/b");
  fTTree->Branch("TPCclsdeu",&fTPCclsdeu,"TPCclsdeu/s");
  fTTree->Branch("ITSdclsmap",&fITSdclsmap,"ITSdclsmap/b");
  fTTree->Branch("pTdeu",&fTpTdeu,"pTdeu/s");
  fTTree->Branch("pdeu",&fTpdeu,"pdeu/s");
  fTTree->Branch("TPCnsigmadeu",&fTTPCnsigmadeu,"TPCnsigmadeu/B");
  fTTree->Branch("Chi2NDFpro",&fTchi2NDFpro,"Chi2NDFpro/b");
  fTTree->Branch("TPCclspro",&fTPCclspro,"TPCclspro/s");
  fTTree->Branch("ITSpclsmap",&fITSpclsmap,"ITSpclsmap/b");
  fTTree->Branch("pTpro",&fTpTpro,"pTpro/s");
  fTTree->Branch("ppro",&fTppro,"ppro/s");
  fTTree->Branch("TPCnsigmapro",&fTTPCnsigmapro,"TPCnsigmapro/B");
  fTTree->Branch("Chi2NDFpion",&fTchi2NDFpion,"Chi2NDFpion/b");
  fTTree->Branch("TPCclspion",&fTPCclspion,"TPCclspion/s");
  fTTree->Branch("ITSpiclsmap",&fITSpiclsmap,"ITSpiclsmap/b");
  fTTree->Branch("pTpion",&fTpTpion,"pTpion/s");
  fTTree->Branch("ppion",&fTppion,"ppion/s");
  fTTree->Branch("TPCnsigmapion",&fTTPCnsigmapion,"TPCnsigmapion/B");
  fTTree->Branch("DCAxypioprim",&fTDCAXYpioprvtx,"DCAxypioprim/S");
  fTTree->Branch("DCAzpioprim",&fTDCAZpioprvtx,"DCAzpioprim/S");
  fTTree->Branch("DCAdp",&fTDCAdp,"DCAdp/s");
  fTTree->Branch("DCAdpi",&fTDCAdpi,"DCAdpi/s");
  fTTree->Branch("DCAppi",&fTDCAppi,"DCAppi/s");
  fTTree->Branch("DCAXYdeuvtx",&fTDCAXYdvtx,"DCAXYdeuvtx/S");
  fTTree->Branch("DCAZdeuvtx",&fTDCAZdvtx,"DCAZdeuvtx/S");
  fTTree->Branch("DCAXYprovtx",&fTDCAXYpvtx,"DCAXYprovtx/S");
  fTTree->Branch("DCAZprovtx",&fTDCAZpvtx,"DCAZprovtx/S");
  fTTree->Branch("DCAXYpionvtx",&fTDCAXYpivtx,"DCAXYpionvtx/S");
  fTTree->Branch("DCAZpionvtx",&fTDCAZpivtx,"DCAZpionvtx/S");
  fTTree->Branch("Angle_dp",&fTAngle_dp,"Angle_dp/s");
  fTTree->Branch("Angle_dpi",&fTAngle_dpi,"Angle_dpi/s");
  fTTree->Branch("Angle_ppi",&fTAngle_ppi,"Angle_ppi/s");
  fTTree->Branch("Rapidity",&fTRapidity,"Rapidity/b");
  fTTree->Branch("DecayLengthProper",&fTDecayLengthProper,"DecayLengthProper/F");
  fTTree->Branch("DecayLengthNorm",&fTDecayLengthNorm,"DecayLengthNorm/s");
  fTTree->Branch("CosPA",&fTCosPA,"CosPA/F");
  fTTree->Branch("TransverseMom",&fTTransverseMom,"TransverseMom/s");
  fTTree->Branch("InvariantMass",&fTInvariantMass,"InvariantMass/F");
  /*fTTree->Branch("TPCclsPIDdeu",&fTPCclsPIDdeu,"TPCclsPIDdeu/s");
  //fTTree->Branch("DCAxydeuprim",&fTDCAXYdeuprvtx,"DCAxydeuprim/F");
  //fTTree->Branch("DCAzdeuprim",&fTDCAZdeuprvtx,"DCAzdeuprim/F");
  //fTTree->Branch("TPCclsPIDpro",&fTPCclsPIDpro,"TPCclsPIDpro/s");
  //fTTree->Branch("DCAxyproprim",&fTDCAXYproprvtx,"DCAxyproprim/F");
  //fTTree->Branch("DCAzproprim",&fTDCAZproprvtx,"DCAzproprim/F");
  //fTTree->Branch("TPCclsPIDpion",&fTPCclsPIDpion,"TPCclsPIDpion/s");*/
  if(fRequireTOFPid){
    fTTree->Branch("TOFnsigmadeu",&fTTOFnsigmadeu,"TOFnsigmadeu/S");
    fTTree->Branch("TOFnsigmapro",&fTTOFnsigmapro,"TOFnsigmapro/S");
    fTTree->Branch("TOFnsigmapion",&fTTOFnsigmapion,"TOFnsigmapion/S");
  }

  if(fMC){
    fTTree->Branch("MCtruth",&fTMCtruth,"MCtruth/F");
    fTTree->Branch("pdeuGen_x",&fTpdeu_gen_X,"pdeuGen_x/F");
    fTTree->Branch("pdeuGen_y",&fTpdeu_gen_Y,"pdeuGen_y/F");
    fTTree->Branch("pdeuGen_z",&fTpdeu_gen_Z,"pdeuGen_z/F");
    fTTree->Branch("pproGen_x",&fTppro_gen_X,"pproGen_x/F");
    fTTree->Branch("pproGen_y",&fTppro_gen_Y,"pproGen_y/F");
    fTTree->Branch("pproGen_z",&fTppro_gen_Z,"pproGen_z/F");
    fTTree->Branch("ppioGen_x",&fTppio_gen_X,"ppioGen_x/F");
    fTTree->Branch("ppioGen_y",&fTppio_gen_Y,"ppioGen_y/F");
    fTTree->Branch("ppioGen_z",&fTppio_gen_Z,"ppioGen_z/F");
    fTTree->Branch("pdgCandDeu",&fTpdgDeu,"pdgCandDeu/I");
    fTTree->Branch("pdgCandPro",&fTpdgPro,"pdgCandPro/I");
    fTTree->Branch("pdgCandPion",&fTpdgPion,"pdgCandPion/I");
    fTTree->Branch("momId_deu",&fTmomidD,"momId_deu/I");
    fTTree->Branch("momId_pro",&fTmomidP,"momId_pro/I");
    fTTree->Branch("momId_pion",&fTmomidPi,"momId_pion/I");
    fTTree->Branch("pdgMomDeu",&fTpdgmomD,"pdgMomDeu/F");
    fTTree->Branch("pdgMomPro",&fTpdgmomP,"pdgMomPro/F");
    fTTree->Branch("pdgMomPion",&fTpdgmomPi,"pdgMomPion/F");
    fTTree->Branch("uniqueID_deu",&fTuniqID_deu,"uniqueID_deu/I");
    fTTree->Branch("uniqueID_pro",&fTuniqID_pro,"uniqueID_pro/I");
    fTTree->Branch("uniqueID_pion",&fTuniqID_pion,"uniqueID_pion/I");
  }

  fTTree->SetAutoSave(150000000);
  PostData(2,fTTree);
  } else {
    fTTree = new TTree();
  }


  //printf("**************************************************\n");
  //printf("end of CreateOutputObjects\n");
  //printf("**************************************************\n");

} // end of UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::UserExec(Option_t *){
  // Main loop
  // Called for each event
  //----------------------------------------------------------------------------
  // Mass definition
  Double_t pionMass            =     0.13957; // pion mass = TDatabasePDG::Instance()->GetParticle(211)->Mass() GeV/c2
  Double_t protonMass          =     0.93827; // proton mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass() GeV/c2
  Double_t deuteronMass        =     1.87561; // deuteron mass = AliPID::ParticleMass(AliPID::kDeuteron) GeV/c2
  Double_t hypertritonMass     =     2.991106;

  //define PDGCodes
  Long_t pdgPionPlus           =              211;
  Long_t pdgPionMinus          =             -211;
  Long_t pdgProton             =             2212;
  Long_t pdgAntiProton         =            -2212;
  Long_t pdgDeuteron           =       1000010020;
  Long_t pdgAntiDeuteron       =      -1000010020;
  Long_t pdgHypertriton        =       1010010030;
  Long_t pdgAntiHypertriton    =      -1010010030;
  //----------------------------------------------------------------------------


  //==========Define variables==========//
  //===PID loop===//
  Bool_t useTOF = kFALSE;
  Int_t ntracks,label = 0 ;
  Int_t labelM_deu, labelM_pro, labelM_pio = 0;
  Double_t chi2PerClusterTPC, nClustersTPC=0.;
  Double_t p, pOverZ, pT = 0.;
  Float_t beta = 0.;
  Float_t dca_prim= 0.;
  Float_t dcaprim[2] = {0.,0.}; //array for the dca from primary vertex coordinates (xy,z) - deuteron
  Float_t dcaprimc[3] = {0.,0.,0.}; // covariance matrix elements of the TPC only impact parameter
  AliESDtrack *track = 0x0;


  //==========ESD event==========
  fESDevent=(AliESDEvent*)InputEvent();
  if(!fESDevent){
    printf("AliAnalysisTaskHypertriton3::Exec(): bad ESD\n");
    PostData(1,fOutput);
    return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler*)mgr->GetInputEventHandler();

  fHistCount->Fill(0); // number of reco events opened

 //==========MC info==========
  AliStack *stack=0x0;
  AliMCEvent *mcEvent = 0x0;
  if(fMC){//MC info and sample selection
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (mgr->GetMCtruthEventHandler());
    if (!eventHandler) {
      printf("ERROR: Could not retrieve MC event handler");
      PostData(1,fOutput);
      return;
    }
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      printf("ERROR: Could not retrieve MC event");
      PostData(1,fOutput);
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      printf("ERROR: stack not available\n");
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
  //printf("====================\n");
  if (!PassCentralitySelection()) {
    PostData(1,fOutput);
    return; //0 bis 80 %
  }
  Bool_t isCent = kFALSE;
  Bool_t isSemiCent = kFALSE;

  if(fCentralityPercentile >= 0. && fCentralityPercentile < 10.) isCent = kTRUE;
  if(fCentralityPercentile >= 10. && fCentralityPercentile < 50.) isSemiCent = kTRUE;

  fCentralityClass = TMath::Floor(fCentralityPercentile/5);

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


  // -------------------------------------------------------
  // Loop for PID on ESD tracks
  // -------------------------------------------------------

  fPIDResponse = handl->GetPIDResponse();

  ntracks = fESDevent->GetNumberOfTracks();

  TArrayI cdeuteron(ntracks);
  UInt_t nDeuTPC = 0;
  TArrayI cantideuteron(ntracks);
  UInt_t nAntiDeuTPC = 0;
  TArrayI cproton(ntracks);
  UInt_t nProTPC = 0;
  TArrayI cantiproton(ntracks);
  UInt_t nAntiProTPC = 0;
  TArrayI cpionplus(ntracks);
  UInt_t nPioPlusTPC = 0;
  TArrayI cpionminus(ntracks);
  UInt_t nPioMinusTPC = 0;

  Bool_t positive = kFALSE;
  Bool_t negative = kFALSE;


  for(Int_t i=0; i < ntracks; i++) {
    track = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(i));
    beta = 0.;
    useTOF = kFALSE;
    positive = kFALSE;
    negative = kFALSE;

    if(fMC){if(track->GetLabel()<0) continue;}
    if(track->GetID()<0) continue;

    // Chi2/TPCcls
    nClustersTPC = track->GetTPCclusters(0);
    chi2PerClusterTPC = track->GetTPCchi2()/nClustersTPC;
    fHistChi2perTPCcluster->Fill(chi2PerClusterTPC);

    if(!PassTrackSelection(track)) continue;

    p = track->P(); //track->GetTPCmomentum()
    pOverZ = p*track->GetSign();
    pT = track->Pt();
    if(track->GetSign() > 0) positive = kTRUE;
    else negative = kTRUE;
    //if(p<0.2) continue;

    //Filling PID histo
    fHistTPCpid->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
    if(fRequireTOFPid){
        useTOF = HasTOF(track, beta);
        if(useTOF)fHistTOFsignal->Fill(p,beta);
    }

    track->GetImpactParameters(dcaprim,dcaprimc);
    dca_prim = TMath::Sqrt((dcaprim[0]*dcaprim[0])+(dcaprim[1]*dcaprim[1]));

    if(fMC){
      TParticle *tparticleDaughter = stack->Particle(TMath::Abs(track->GetLabel()));
        if(pT > fMaxPtDeuteron) continue;
        if(pT >= fMinPtDeuteron){
            if(tparticleDaughter->GetPdgCode() == 1000010020 || tparticleDaughter->GetPdgCode() == -1000010020){ //d
                fHistTPCdeusignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
                if(useTOF){
                    fHistTOFdeusignal->Fill(p,beta);
                    //fHistTOFdeumass->Fill(mass);
                }
                if(positive) cdeuteron[nDeuTPC++] = i;
                if(negative) cantideuteron[nAntiDeuTPC++] = i;
            }
        }
        if(pT >= fMinPtProton && pT <= fMaxPtProton){
          if (tparticleDaughter->GetPdgCode() == 2212 || tparticleDaughter->GetPdgCode() == -2212) { // p
              fHistTPCprosignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
              if(useTOF){
                  fHistTOFprosignal->Fill(p,beta);
                  //fHistTOFpromass->Fill(mass);
              }
              if(positive) cproton[nProTPC++] = i;
              if(negative) cantiproton[nAntiProTPC++] = i;
          }
        }
        if(!fESDtrackCutsV0->AcceptTrack(track)) continue;
        if(tparticleDaughter->GetPdgCode() == 211 || tparticleDaughter->GetPdgCode() == -211) { // pi+
            fHistTPCpionsignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
            if(positive) cpionplus[nPioPlusTPC++] = i;
            if(negative) cpionminus[nPioMinusTPC++] = i;
        }
    } // end of MC PID
    else {
    if(pT > fMaxPtDeuteron) continue;
    if(pT >= fMinPtDeuteron && TMath::Abs(dcaprim[1])<= fDCAzDPVmax){
        if(PassPIDSelection(track, AliPID::kDeuteron, useTOF,fDeuteronTPCSigmas)) { //deuteron
          fHistTPCdeusignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
          fHistDCAdprimary->Fill(dca_prim);
          fHistDCAXYdprimary->Fill(dcaprim[0]);
          fHistDCAZdprimary->Fill(dcaprim[1]);
          if(useTOF){
              fHistTOFdeusignal->Fill(p,beta);
              //fHistTOFdeumass->Fill(mass);
          }
          if(positive) cdeuteron[nDeuTPC++] = i;
          if(negative) cantideuteron[nAntiDeuTPC++] = i;
       }
     }
    if(pT >= fMinPtProton && pT <= fMaxPtProton){
      if(PassPIDSelection(track,AliPID::kProton, useTOF,fProtonTPCSigmas)) { // proton
          fHistTPCprosignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
          fHistDCApprimary->Fill(dca_prim);
          fHistDCAXYpprimary->Fill(dcaprim[0]);
          fHistDCAZpprimary->Fill(dcaprim[1]);
          if(useTOF){
              fHistTOFprosignal->Fill(p,beta);
              //fHistTOFpromass->Fill(mass);
            }
            if(positive) cproton[nProTPC++] = i;
            if(negative) cantiproton[nAntiProTPC++] = i;
      }
    }
      if(!fESDtrackCutsV0->AcceptTrack(track)) continue;

      if(PassPIDSelection(track,AliPID::kPion, useTOF,fPionTPCSigmas)) { //pion^+
          fHistTPCpionsignal->Fill(track->GetTPCmomentum()*track->GetSign(), track->GetTPCsignal());
          fHistDCApiprimary->Fill(dca_prim);
          fHistDCAXYpiprimary->Fill(dcaprim[0]);
          fHistDCAZpiprimary->Fill(dcaprim[1]);
          fHistpionTPCcls->Fill(track->GetTPCclusters(0));
          if(useTOF){
            fHistTOFpionsignal->Fill(p,beta);
          }
          if(dca_prim < fDCAPiPVmin) continue;
          if(positive) cpionplus[nPioPlusTPC++] = i;
          if(negative) cpionminus[nPioMinusTPC++] = i;
      }
    }
  } // end of PID loop

  cdeuteron.Set(nDeuTPC);
  cantideuteron.Set(nAntiDeuTPC);
  cproton.Set(nProTPC);
  cantiproton.Set(nAntiProTPC);
  cpionplus.Set(nPioPlusTPC);
  cpionminus.Set(nPioMinusTPC);

  if(fMinvSignal){
    //Hypertriton Invariant Mass
    if(fChooseMatter) CombineThreeTracks(kTRUE,cdeuteron, cproton, cpionminus,isCent, isSemiCent);

    //Anti-Hypertriton Invariant Mass
    if(fChooseAntiMatter) CombineThreeTracks(kFALSE,cantideuteron, cantiproton, cpionplus, isCent, isSemiCent);
  }

  if(fMinvLikeSign){
    //Hypertriton Invariant Mass - LS
    CombineThreeTracks(kTRUE,cdeuteron, cproton, cpionplus, isCent, isSemiCent);

    //Anti-Hypertriton Invariant Mass - LS
    CombineThreeTracks(kFALSE,cantideuteron, cantiproton, cpionminus, isCent, isSemiCent);
  }

  if(fMC){
    for(Int_t iMC=0; iMC<stack->GetNtrack(); iMC++){ // check MC stack content
      TParticle *pstack = stack->Particle(iMC);
      if(pstack->GetPdgCode() == 11) fHistParticle->Fill(0); // e-
      if(pstack->GetPdgCode() == -11) fHistParticle->Fill(1); // e+
      if(pstack->GetPdgCode() == 211) fHistParticle->Fill(2); // pi+
      if(pstack->GetPdgCode() == -211) fHistParticle->Fill(3); // pi-
      if(pstack->GetPdgCode() == 321) fHistParticle->Fill(4); // k+
      if(pstack->GetPdgCode() == -321) fHistParticle->Fill(5); // k-
      if(pstack->GetPdgCode() == 2212) fHistParticle->Fill(6); // p
      if(pstack->GetPdgCode() == -2212) fHistParticle->Fill(7); // pbar
      if(pstack->GetPdgCode() == 1000010020) fHistParticle->Fill(8); // d
      if(pstack->GetPdgCode() == -1000010020) fHistParticle->Fill(9); // dbar
      if(pstack->GetPdgCode() == 1000010030) fHistParticle->Fill(10); // t
      if(pstack->GetPdgCode() == -1000010030) fHistParticle->Fill(11); // tbar
      if(pstack->GetPdgCode() == 1000020030) fHistParticle->Fill(12); // He3
      if(pstack->GetPdgCode() == -1000020030) fHistParticle->Fill(13); // He3bar
      if(pstack->GetPdgCode() == 1000020040) fHistParticle->Fill(14); // He4
      if(pstack->GetPdgCode() == -1000020040) fHistParticle->Fill(15); // He4bar
      if(pstack->GetPdgCode() == 1010010030 && pstack->GetNDaughters() == 3) fHistParticle->Fill(16); // H3L
      if(pstack->GetPdgCode() == -1010010030 && pstack->GetNDaughters() == 3) fHistParticle->Fill(17); // H3Lbar
      if(pstack->GetPdgCode() == 1010010030 && pstack->GetNDaughters() == 2) fHistParticle->Fill(18); // H3L
      if(pstack->GetPdgCode() == -1010010030 && pstack->GetNDaughters() == 2) fHistParticle->Fill(19); // H3Lbar
      if(pstack->GetPdgCode() == 3122 && pstack->GetNDaughters() == 2) fHistParticle->Fill(20); // Lambda
    } //end of MC content truth check
  }

  cdeuteron.Reset();
  cantideuteron.Reset();
  cproton.Reset();
  cantiproton.Reset();
  cpionplus.Reset();
  cpionminus.Reset();

  // Post output data.
  PostData(1,fOutput);
  if(fFillTree) PostData(2,fTTree);



} // end of UserExec

//________________________________________________________________________
void AliAnalysisTaskHypertriton3::Terminate(Option_t *){
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
