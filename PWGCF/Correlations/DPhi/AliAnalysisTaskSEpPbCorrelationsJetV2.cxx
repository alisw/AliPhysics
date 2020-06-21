/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/


/**************************************************************************************************
* AliAnalysisTaskSEpPbCorrelationsJetV2:
* This task is developed from AliAnalysisTaskSEpPbCorrelationsYS.cxx, aimed to calculate the TPC pairs 
* - FMD correlation, which can obtain the V2 of jet particles.
**************************************************************************************************/


#include "AliAnalysisManager.h"
#include "TGrid.h"
#include "AliLog.h"
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TH3.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMap.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TTree.h>
#include <TProfile.h>
#include <TObjectTable.h>
#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliCFContainer.h"
#include "AliGenEventHeader.h"
#include "AliTHn.h"

#include "AliAODEvent.h"
#include "AliESDAD.h"
#include "AliESDEvent.h"
#include "AliVAD.h"
#include "AliAODForwardMult.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliESDVertex.h"
//#include "AliAODPid.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAODcascade.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliCentrality.h"
#include "AliEventPoolManager.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "Riostream.h"

//#include "AliFlowEventSimple.h"
///#include "AliFlowVector.h"/
//#include "AliFlowTrackSimple.h"
//#include "AliAnalysisTaskMSP.h"
//#include "AliFlowAnalysisWithMSP.h"
#include "AliMultSelection.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliForwardCorrectionManager.h"
//#include "AliForwardUtil.h"

#include "AliAnalysisTaskSEpPbCorrelationsJetV2.h"
#include "AliAnalysisTaskSEpPbCorrelationsYS.h"
ClassImp(AliAnalysisTaskSEpPbCorrelationsJetV2)
ClassImp(AliAssociatedTrackYS)
ClassImp(AliAssociatedTPCPairs)
ClassImp(AliMixTrackYS)
ClassImp(AliAssociatedVZEROYS)

AliAnalysisTaskSEpPbCorrelationsJetV2::AliAnalysisTaskSEpPbCorrelationsJetV2()
    : AliAnalysisTaskSE(),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fptdiff(kFALSE),
      fmakehole(kFALSE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fTPCTPClist(0),
      //fPIDResponse(0),
      ffilterbit(5),
      fnoClusters(70),
      fCutChargedDCAzMax(5),
      fCutChargedDCAxyMax(15),
      fPtMin(0.2),
      fPtMax(3.0),
      fCenMin(0.),
      fCenMax(10.),
      fEtaMax(0.8),
      fEtaMaxExtra(0.),
      fEtaMinExtra(-0.8),
      fEtaMaxV0(0.8),
      fEtaMinV0(0.),
      fdcaDaughtersToPrimVtx(0.06),
      fdcaBetweenDaughters(1.0),
      fRadiMin(0.5),
      fRadiMax(100),
      fcutcTauLam(30),
      fcutcTauK0(20),
      fcosMinK0s(0.97),
      fcosMinLambda(0.995),
      fMaxnSigmaTPCV0(5),
      hv0dcharge(0),
      fclustermin(70),
      fratiocluster(0.8),
      fEtaMaxDaughter(0.8),
      fEtaMinDaughter(0.),
      fHistMass_K0s(0),
      fHistMass_K0s_MC(0),
      fHistMass_Lambda(0),
      fHistMass_ALambda(0),
      fHistMass_ALambda_MC(0),
      fHistMassXiMinus(0),
      fHistMassXiPlus(0),
      fHistMassOmegaMinus(0),
      fHistMassOmegaPlus(0),
      fHistMass_bumpcorr(0),
      fHist_V0QA(0),
      fHist_CascadeQA(0),
      fHistMass_Lambda_MC(0),
      fEventCuts(0),
      fUtils(),
      fEvent(0),
      mcEvent(0),
      lPrimaryBestVtx(0),
      fPrimaryZVtx(0),
      fvzero(0),
      fPoolMgr(0),
      fPoolMgr1(0),
      poolmin(0),
      poolmax(0),
      fPoolMaxNEvents(2000),
      fPoolMinNTracks(50000),
      fMinEventsToMix(5),
      fNzVtxBins(0),
      fNCentBins(0),
      fMaxnSigmaTPCTOF(3.),
      fh2_pt_trig_asso(0),
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentrality_beforeFMDMulcut(0),
      fHistCentzvertex(0),
      fHistCentV0vsTracklets(0),
      fHistCentV0vsTrackletsbefore(0),
      mixedDist(0),
      mixedDist2(0),
      //fHistLeadQA(0),      
      fHistPIDQA(0),
      fhistmcprim(0),
      fhmcprimvzeta(0),
      frefetaa(0),
      frefetac(0),
      frefvz(0),
      fhmcprimpdgcode(0),
      fh2_FMD_acceptance_prim(0),
      fh2_FMD_eta_phi_prim(0),
      fh2_FMD_acceptance(0),
      fh2_ITS_acceptance(0),
      fh2_SPD_multcorr(0),
      fh2_SPDV0_multcorr(0),
      fh2_SPDtrack_multcorr(0),
      fhtrackletsdphi(0),
      fh2_FMD_eta_phi(0),
      fHist_NeventRun(0),
      fHist_V0AMultRun(0),
      fHist_V0CMultRun(0),
      fHist_FMDAMultRun(0),
      fHist_FMDCMultRun(0),
      fhistfmdphiacc(0),  
      fhFMDmultchannel(0),
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fFMDV0same(0),
      fFMDV0same_post(0),
      fFMDV0Asame(0),
      fFMDV0Asame_post(0),
      fFMDV0Csame(0),
      fFMDV0Csame_post(0),
      fHist_vzeromult(0),
      fHist_vzeromultEqweighted(0),
      fHist2dmult(0),
      fHistVZERO(0),
      fHist_Stat(0),
      fHist_V0Stat(0),
      fHistPhiDTPCNSig(0),
      fHistPhiDTOFNSig(0),
      fHistPhiDTPCTOFNSig(0),
      fHistMass_PhiMeson(0),
      fHistMass_PhiMeson_MIX(0),
      fHist_PhiQA(0),
      fHistTriggerTrack(0),
      fHistReconstTrack(0),
      fHistTriggerTrackMix(0),
      fHistReconstTrackMix(0),
      fHistQna(0),
      fHistQnc(0),
      fHistQn(0),
      fHistQna_VZERO(0),
      fHistQnc_VZERO(0),
      fHistQn_VZERO(0),
      fHistVn(0),
      SP_TPCATPCC(0),
      SP_TPCATPCC_default(0),
      SP_V0AV0C_default(0),
      SP_V0ATPC_default(0),
      SP_V0CTPC_default(0),
      fHist_V0AV0C(0),
      fHist_V0ATPC(0),
      fHist_V0CTPC(0),
      SP_uTPCA(0),
      SP_uTPCC(0){
  for (Int_t iBin = 0; iBin < 100; iBin++) {
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  for (Int_t i = 0; i < 3; i++) {
    tPrimaryVtxPosition[i] = 0;
  }
  for (Int_t i = 0; i < 6; i++) {
    fHistPosNsig[i] = 0;
    fHistNegNsig[i] = 0;
    fHistPosNsigQA[i] = 0;
    fHist_AP[i] = 0;
  }
  for(Int_t i=0;i<3;i++){
    fh3NegNsig[i]=0;
    fh3PosNsig[i]=0;
  }
  for (Int_t i = 0; i < 6; i++) {
    fHistNsig[i]=0;
    fHistNsigcorr[i]=0;
  }
  for (Int_t i = 0; i < 4; i++) {
    fHistQAQB[i]=0;
    fHistQAQB_VZERO[i]=0;
    fHistCorrQna[i]=0;
    fHistCorrQnc[i]=0;
  }
  for (Int_t i = 0; i < 8; i++) {
    SP_uTPC_PP[i]=0;
    SP_uTPC[i]=0;
    SP_uTPC1[i]=0;
    SP_uTPC2[i]=0;
    SP_uTPC3[i]=0;
    SP_uVZEROA_PP[i]=0;
    SP_uVZEROA[i]=0;
    SP_uVZEROA1[i]=0;
    SP_uVZEROA2[i]=0;
    SP_uVZEROA3[i]=0;
    SP_uVZEROC_PP[i]=0;
    SP_uVZEROC[i]=0;
    SP_uVZEROC1[i]=0;
    SP_uVZEROC2[i]=0;
    SP_uVZEROC3[i]=0;
  }
  for(Int_t i=0;i<4;i++){
    fhrefetaFMD[i]=0;
    fhrefphiFMD[i]=0;
  }
  for(Int_t i=0;i<10;i++){
    fhcorr[i]=0;
  }
  for(Int_t i=0;i<31;i++){
    fhFMDmult_runbyrun_cside[i]=0;
  }
  for(Int_t i=0;i<65;i++){
    fhFMDmult_runbyrun_aside[i]=0;
  }
     
}
AliAnalysisTaskSEpPbCorrelationsJetV2::AliAnalysisTaskSEpPbCorrelationsJetV2(const char *name)
    : AliAnalysisTaskSE(name),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fptdiff(kFALSE),
      fmakehole(kFALSE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fTPCTPClist(0),
      //fPIDResponse(0),
      ffilterbit(5),
      fnoClusters(70),
      fCutChargedDCAzMax(5),
      fCutChargedDCAxyMax(15),
      fPtMin(0.2),
      fPtMax(3.0),
      fCenMin(0.),
      fCenMax(10.),
      fEtaMax(0.8),
      fEtaMaxExtra(0.),
      fEtaMinExtra(-0.8),
      fEtaMaxV0(0.8),
      fEtaMinV0(0.),
      fdcaDaughtersToPrimVtx(0.06),
      fdcaBetweenDaughters(1.0),
      fRadiMin(0.5),
      fRadiMax(100),
      fcutcTauLam(30),
      fcutcTauK0(20),
      fcosMinK0s(0.97),
      fcosMinLambda(0.995),
      fMaxnSigmaTPCV0(5),
      hv0dcharge(0),
      fclustermin(70),
      fratiocluster(0.8),
      fEtaMaxDaughter(0.8),
      fEtaMinDaughter(0.),
      fHistMass_K0s(0),
      fHistMass_K0s_MC(0),
      fHistMass_Lambda(0),
      fHistMass_ALambda(0),
      fHistMass_ALambda_MC(0),
      fHistMassXiMinus(0),
      fHistMassXiPlus(0),
      fHistMassOmegaMinus(0),
      fHistMassOmegaPlus(0),
      fHistMass_bumpcorr(0),
      fHist_V0QA(0),
      fHist_CascadeQA(0),
      fHistMass_Lambda_MC(0),
      fEventCuts(0),
      fUtils(),
      fEvent(0),
      mcEvent(0),
      lPrimaryBestVtx(0),
      fPrimaryZVtx(0),
      fvzero(0),
      fPoolMgr(0),
      fPoolMgr1(0),
      poolmin(0),
      poolmax(0),
      fPoolMaxNEvents(2000),
      fPoolMinNTracks(50000),
      fMinEventsToMix(5),
      fNzVtxBins(0),
      fNCentBins(0),
      fMaxnSigmaTPCTOF(3.),
      fh2_pt_trig_asso(0),
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentrality_beforeFMDMulcut(0),
      fHistCentzvertex(0),
      fHistCentV0vsTracklets(0),
      fHistCentV0vsTrackletsbefore(0),
      mixedDist(0),
      mixedDist2(0),
      //fHistLeadQA(0),      
      fHistPIDQA(0),
      fhistmcprim(0),
      fhmcprimvzeta(0),
      frefetaa(0),
      frefetac(0),
      frefvz(0),
      fhmcprimpdgcode(0),
      fh2_FMD_acceptance_prim(0),
      fh2_FMD_eta_phi_prim(0),
      fh2_FMD_acceptance(0),
      fh2_ITS_acceptance(0),
      fh2_SPD_multcorr(0),
      fh2_SPDV0_multcorr(0),
      fh2_SPDtrack_multcorr(0),
      fhtrackletsdphi(0),
      fh2_FMD_eta_phi(0),
      fHist_NeventRun(0),
      fHist_V0AMultRun(0),
      fHist_V0CMultRun(0),
      fHist_FMDAMultRun(0),
      fHist_FMDCMultRun(0),
      fhistfmdphiacc(0),  
      fhFMDmultchannel(0),
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fFMDV0same(0),
      fFMDV0same_post(0),
      fFMDV0Asame(0),
      fFMDV0Asame_post(0),
      fFMDV0Csame(0),
      fFMDV0Csame_post(0),
      fHist_vzeromult(0),
      fHist_vzeromultEqweighted(0),
      fHist2dmult(0),
      fHistVZERO(0),
      fHist_Stat(0),
      fHist_V0Stat(0),
      fHistPhiDTPCNSig(0),
      fHistPhiDTOFNSig(0),
      fHistPhiDTPCTOFNSig(0),
      fHistMass_PhiMeson(0),
      fHistMass_PhiMeson_MIX(0),
      fHist_PhiQA(0),
      fHistTriggerTrack(0),
      fHistReconstTrack(0),
      fHistTriggerTrackMix(0),
      fHistReconstTrackMix(0),
      fHistQna(0),
      fHistQnc(0),
      fHistQn(0),
      fHistQna_VZERO(0),
      fHistQnc_VZERO(0),
      fHistQn_VZERO(0),
      fHistVn(0),
      SP_TPCATPCC(0),
      SP_TPCATPCC_default(0),
      SP_V0AV0C_default(0),
      SP_V0ATPC_default(0),
      SP_V0CTPC_default(0),
      fHist_V0AV0C(0),
      fHist_V0ATPC(0),
      fHist_V0CTPC(0),
      SP_uTPCA(0),
      SP_uTPCC(0){
        for (Int_t iBin = 0; iBin < 100; iBin++) {
          fZvtxBins[iBin] = 0.;
          fCentBins[iBin] = 0.;
        }
        for (Int_t i = 0; i < 3; i++) {
          tPrimaryVtxPosition[i] = 0;
        }
	
	for (Int_t i = 0; i < 6; i++) {
	  fHistPosNsig[i] = 0;
	  fHistNegNsig[i] = 0;
	  fHistPosNsigQA[i] = 0 ;
	  fHist_AP[i] = 0;
	}
	for(Int_t i=0;i<3;i++){
	  fh3NegNsig[i]=0;
	  fh3PosNsig[i]=0;
	}
	
        for (Int_t i = 0; i < 6; i++) {
          fHistNsig[i] = 0;
          fHistNsigcorr[i]=0;
        }
	
        for (Int_t i = 0; i < 4; i++) {
          fHistQAQB[i]=0;
          fHistQAQB_VZERO[i]=0;
          fHistCorrQna[i]=0;
          fHistCorrQnc[i]=0;
        }
	
        for (Int_t i = 0; i < 8; i++) {
          SP_uTPC_PP[i]=0;
          SP_uTPC[i]=0;
          SP_uTPC1[i]=0;
          SP_uTPC2[i]=0;
          SP_uTPC3[i]=0;
          SP_uVZEROA_PP[i]=0;
          SP_uVZEROA[i]=0;
          SP_uVZEROA1[i]=0;
          SP_uVZEROA2[i]=0;
          SP_uVZEROA3[i]=0;
          SP_uVZEROC_PP[i]=0;
          SP_uVZEROC[i]=0;
          SP_uVZEROC1[i]=0;
          SP_uVZEROC2[i]=0;
          SP_uVZEROC3[i]=0;
        }
        for(Int_t i=0;i<4;i++){
          fhrefetaFMD[i]=0;
          fhrefphiFMD[i]=0;
        }
        for(Int_t i=0;i<10;i++){
          fhcorr[i]=0;
        }
	for(Int_t i=0;i<31;i++){
	  fhFMDmult_runbyrun_cside[i]=0;
	}
	
	for(Int_t i=0;i<65;i++){
	  fhFMDmult_runbyrun_aside[i]=0;
	}

	  
        DefineOutput(1, TList::Class());
        DefineOutput(2, TList::Class());
        DefineOutput(3, TList::Class());
      }

AliAnalysisTaskSEpPbCorrelationsJetV2::~AliAnalysisTaskSEpPbCorrelationsJetV2()
{
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
    fOutputList = 0x0;
  }
  
  if (fOutputList1 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList1;
    fOutputList1 = 0x0;
  }
  
  if (fOutputList2 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList2;
    fOutputList2 = 0x0;
  }
/*  
  if (fPIDResponse) {
    delete fPIDResponse;
    fPIDResponse = 0;
  }
*/
}
void AliAnalysisTaskSEpPbCorrelationsJetV2::UserCreateOutputObjects() {
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fOutputList->SetName("global");
  DefineGeneralOutput();
  PostData(1, fOutputList);

  fOutputList1 = new TList();
  fOutputList1->SetOwner(kTRUE);
  fOutputList1->SetName("anahistos");
  DefineCorrOutput();

  //DefineVZEROOutput();
  PostData(2, fOutputList1);

  fOutputList2 = new TList();
  fOutputList2->SetOwner(kTRUE);
  fOutputList2->SetName("QA");
  DefinedQAHistos();
  PostData(3, fOutputList2);

  fEventCuts.AddQAplotsToList(fOutputList);

  frefetac=new TH1F("frefetac","frefetac",30,-3.4,-1.9);
  fOutputList2->Add(frefetac);
  frefetaa=new TH1F("frefetaa","frefetaa",62,1.9,5.0);
  fOutputList2->Add(frefetaa);
  
  frefvz=new TH1F("frefvz","z-vertex",10,-10,10);
  fOutputList2->Add(frefvz);
 
   fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
   if (!fPoolMgr)
   return;
   fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);


   fPoolMgr1 = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
   if (!fPoolMgr1)
   return;
   fPoolMgr1->SetTargetValues(fPoolMinNTracks, 0.1, 5);
 }
 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefineGeneralOutput() {

   fHist_Stat = new TH1F("fHist_Stat", "Stat Histogram", 14, -0.5, 13.5);
   fHist_Stat->GetXaxis()->SetBinLabel(1, "All Events");
   fHist_Stat->GetXaxis()->SetBinLabel(2, "Analyzed Events");
   fHist_Stat->GetXaxis()->SetBinLabel(3, "MultSelection OK");
   fHist_Stat->GetXaxis()->SetBinLabel(4, "Vertex OK");
   fHist_Stat->GetXaxis()->SetBinLabel(5, "Centrality OK");
   fHist_Stat->GetXaxis()->SetBinLabel(6, "FMD OK");
   fHist_Stat->GetXaxis()->SetBinLabel(7, "FMD multi cut");
   fHist_Stat->GetXaxis()->SetBinLabel(8, "--");
   fHist_Stat->GetXaxis()->SetBinLabel(9, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(10, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(11, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(12, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(13, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(14, "-");
   fOutputList->Add(fHist_Stat);

   fHist_V0Stat = new TH1F("fHist_V0Stat", "Stat Histogram", 16, -0.5, 15.5);
   fHist_V0Stat->GetXaxis()->SetBinLabel(1, "all");
   fHist_V0Stat->GetXaxis()->SetBinLabel(2, "On-Fly");
   fHist_V0Stat->GetXaxis()->SetBinLabel(3, "Off-Fly");
   fHist_V0Stat->GetXaxis()->SetBinLabel(4, "V0 pseudorapidity");
   fHist_V0Stat->GetXaxis()->SetBinLabel(5, "DCA Dau. tracks to PV");
   fHist_V0Stat->GetXaxis()->SetBinLabel(6, "DCA dauthers");
   fHist_V0Stat->GetXaxis()->SetBinLabel(7, "Fiducial volume");
   fHist_V0Stat->GetXaxis()->SetBinLabel(8, "Pass IsAcceptedV0");
   fHist_V0Stat->GetXaxis()->SetBinLabel(9, "track cut");
   fHist_V0Stat->GetXaxis()->SetBinLabel(10, "charge");
   fHist_V0Stat->GetXaxis()->SetBinLabel(11, "PID for K0s");
   fHist_V0Stat->GetXaxis()->SetBinLabel(12, "ctau for k0s");
   fHist_V0Stat->GetXaxis()->SetBinLabel(13, "AP cut for K0s");
   fHist_V0Stat->GetXaxis()->SetBinLabel(14, "-");
   fHist_V0Stat->GetXaxis()->SetBinLabel(15, "-");
   fHist_V0Stat->GetXaxis()->SetBinLabel(16, "-");
   fOutputList->Add(fHist_V0Stat);


   fHistzvertex = new TH1F("fHistzvertex", ";VZ;count", 60, -15, 15);
   fOutputList->Add(fHistzvertex);

   Double_t fmaxcent;
   if(fCentType=="Manual"){
     fmaxcent=200;
   }else{
     if (fcollisiontype=="HMPP") fmaxcent=1.;
     else fmaxcent=100.;
   }
   fHistCentrality = new TH1F("fHistCentrality", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality);

   fHistCentrality_beforecut = new TH1F("fHistCentrality_beforecut", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality_beforecut);
 
   fHistCentrality_beforeFMDMulcut = new TH1F("fHistCentrality_beforeFMDMulcut", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality_beforeFMDMulcut);

 
   fHistCentzvertex = new TH2F("fHistCentzvertex", "Cent;VZ;count", 100,0, fmaxcent, 60, -15, 15);
   fOutputList->Add(fHistCentzvertex);

   fHistCentV0vsTracklets=new TH2F("fHistCentV0vsTracklets","fHistCentV0vsTracklets",100,0,100,1000,0,10000);
   fOutputList->Add(fHistCentV0vsTracklets);
   
   fHistCentV0vsTrackletsbefore=new TH2F("fHistCentV0vsTrackletsbefore","fHistCentV0vsTrackletsbefore",100,0,100,1000,0,10000);
   fOutputList->Add(fHistCentV0vsTrackletsbefore);
   
   
   TTree *settingsTree = new TTree("UEAnalysisSettings", "Analysis Settings in UE estimation");
   settingsTree->Branch("fZVertex", &fZVertex, "fZVertex/D");
   settingsTree->Branch("fEtaMax", &fEtaMax, "fEtaMax/D");
   settingsTree->Branch("fPtMin", &fPtMin, "fPtMin/D");
   settingsTree->Branch("fMaxnSigmaTPCTOF", &fMaxnSigmaTPCTOF, "fMaxnSigmaTPCTOF/D");

   settingsTree->Fill();
//   fOutputList->Add(settingsTree);
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefinedQAHistos() {
   Int_t ncentmax;
   if(fCentType=="Manual") ncentmax=200;
   else ncentmax=100;
      
   mixedDist=new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, ncentmax, 200, 0, fPoolMinNTracks*1.5 );
   mixedDist2=new TH2F("mixedDist2", ";centrality;events;events", 101, 0,ncentmax, 100, 0, 1000) ;
   fOutputList2->Add(mixedDist);
   fOutputList2->Add(mixedDist2);


   const Int_t ipidBin[4] = {11, 40, 72, 15};
   Double_t binning_pt_lead[12] = {0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
                                   2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
   Double_t binning_eta[41] = {-1.,   -0.95, -0.9,  -0.85, -0.8,  -0.75, -0.7,
                               -0.65, -0.6,  -0.55, -0.5,  -0.45, -0.4,  -0.35,
                               -0.3,  -0.25, -0.2,  -0.15, -0.1,  -0.05, 0.,
                               0.05,  0.1,   0.15,  0.2,   0.25,  0.3,   0.35,
                               0.4,   0.45,  0.5,   0.55,  0.6,   0.65,  0.7,
                               0.75,  0.8,   0.85,  0.9,   0.95,  1.0};
   Double_t binning_dphi[73] = {
       -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464,
       -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865,
       -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266,
       0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332,
       0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931,
       1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530,
       1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129,
       2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727,
       2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326,
       3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925,
       3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524,
       4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123,
       4.712389};
   Double_t binning_cent[16] = {0.,  1.,  2.,  3.,  4.,  5.,  10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1};
     Double_t binning_zvx[11] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
   if(fasso=="PID" && fQA){
   fHistPIDQA = new AliTHn("fHistPIDQA", "fHistPIDQA", 3, 4, ipidBin);
   fHistPIDQA->SetBinLimits(0, binning_pt_lead);
   fHistPIDQA->SetBinLimits(1, binning_eta);
   fHistPIDQA->SetBinLimits(2, binning_dphi);
   fHistPIDQA->SetBinLimits(3, binning_cent);
   fHistPIDQA->SetVarTitle(0, "pt");
   fHistPIDQA->SetVarTitle(1, "eta");
   fHistPIDQA->SetVarTitle(2, "phi");
   fHistPIDQA->SetVarTitle(3, "centrality");
   fOutputList1->Add(fHistPIDQA);
   }

/*
   fHistLeadQA = new AliTHn("fHistLeadQA", "fHistLeadQA", 1, 4, ipidBin);
   fHistLeadQA->SetBinLimits(0, binning_pt_lead);
   fHistLeadQA->SetBinLimits(1, binning_eta);
   fHistLeadQA->SetBinLimits(2, binning_dphi);
   if(fCentType=="Manual")fHistLeadQA->SetBinLimits(3,0,300);
   else fHistLeadQA->SetBinLimits(3, binning_cent);
   fHistLeadQA->SetVarTitle(0, "pt");
   fHistLeadQA->SetVarTitle(1, "eta");
   fHistLeadQA->SetVarTitle(2, "phi");
   fHistLeadQA->SetVarTitle(3, "centrality");
   fOutputList1->Add(fHistLeadQA);
*/

   const Int_t imcprimbin[4]={11,55,30,10};
   Double_t binning_pt_mcprim[12] = {0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
   if(!fDataType){
     fhistmcprim=new AliTHn("fhistmcprim","fhistmcprim",1,4,imcprimbin);
     fhistmcprim->SetBinLimits(0,binning_pt_mcprim);
     fhistmcprim->SetBinLimits(1,-5.5,5.5);
     fhistmcprim->SetBinLimits(2,0.,2*TMath::Pi());
     fhistmcprim->SetBinLimits(3,0.,100.);
     fhistmcprim->SetVarTitle(0,"pt");
     fhistmcprim->SetVarTitle(1,"eta");
     fhistmcprim->SetVarTitle(2,"phi");
     fhistmcprim->SetVarTitle(3,"centrality");
     //fOutputList2->Add(fhistmcprim);
     fhmcprimvzeta=new TH2D("fhmcprimvzeta","fhmcprimvzeta",200,-4,6,20,-10,10);
     //fOutputList2->Add(fhmcprimvzeta);
     fhmcprimpdgcode=new TH1D("fhmcprimpdgcode","fhmcprimpdgcode",4000,-0.5,3999.5);
     //fOutputList2->Add(fhmcprimpdgcode);
     fh2_FMD_acceptance_prim=new TH2D("fh2_FMD_acceptance_prim","fh2_FMD_acceptance_prim",200,-4,6,200,-10,10);
     //fOutputList2->Add(fh2_FMD_acceptance_prim);
     fh2_FMD_eta_phi_prim=new TH2D("fh2_FMD_eta_phi_prim","fh2_FMD_eta_phi_prim",200,-4,6,20,0,2*TMath::Pi());
     //fOutputList2->Add(fh2_FMD_eta_phi_prim);

     
     
     for(Int_t i=0;i<4;i++){
       fhrefetaFMD[i]=new TH1D(Form("fhrefetaFMD_%d",i),Form("fhrefetaFMD_%d",i),200,-4,6);
       fhrefphiFMD[i]=new TH1D(Form("fhrefphiFMD_%d",i),Form("fhrefphiFMD_%d",i),100,0,2*TMath::Pi());
      // fOutputList2->Add(fhrefetaFMD[i]);
      // fOutputList2->Add(fhrefphiFMD[i]);
     }
   }

	 fHist_NeventRun=new TH1F("fHist_NeventRun","fHist_NeventRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_NeventRun);

	 fHist_V0AMultRun=new TH1F("fHist_V0AMultRun","fHist_V0AMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_V0AMultRun);
	 fHist_V0CMultRun=new TH1F("fHist_V0CMultRun","fHist_V0CMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_V0CMultRun);
	 fHist_FMDAMultRun=new TH1F("fHist_FMDAMultRun","fHist_FMDAMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_FMDAMultRun);
	 fHist_FMDCMultRun=new TH1F("fHist_FMDCMultRun","fHist_FMDCMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_FMDCMultRun);



	 Float_t nmutFoward;
	 if(fcollisiontype=="PbPb") nmutFoward=2000.;
	 else nmutFoward=1000.;
	 
	 fFMDV0 = new TH2F("FMDV0", "FMD vs V0 pre cut;FMD;V0;",2000, 0, 2*nmutFoward, 2000, 0, 2*nmutFoward);
//	 fOutputList2->Add(fFMDV0);
	 fFMDV0_post=new TH2F("FMDV0_post", "FMD vs V0 post cut;FMD;V0;",2000, 0, 2*nmutFoward, 2000, 0, 2*nmutFoward);
//	 fOutputList2->Add(fFMDV0_post);
	 fFMDV0A = new TH2F("FMDV0A", "FMD vs V0A;FMD;V0A;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0A);
	 fFMDV0A_post = new TH2F("FMDV0A_post", "FMD vs V0A post cut;FMD;V0A;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0A_post);
	 fFMDV0C = new TH2F("FMDV0C", "FMD vs V0C;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0C);
	 fFMDV0C_post = new TH2F("FMDV0C_post", "FMD vs V0C post cut;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0C_post);
	 
	 fFMDV0same = new TH2F("FMDV0same", "FMD vs V0 pre cut;FMD;V0;",2000, 0, 2000, 2*nmutFoward, 0, 2*nmutFoward);
//	 fOutputList2->Add(fFMDV0same);
	 fFMDV0same_post=new TH2F("FMDV0same_post", "FMD vs V0 post cut;FMD;V0;",2000, 0, 2000, 2*nmutFoward, 0, 2*nmutFoward);
//	 fOutputList2->Add(fFMDV0same_post);
	 fFMDV0Asame = new TH2F("FMDV0Asame", "FMD vs V0A;FMD;V0A;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0Asame);
	 fFMDV0Asame_post = new TH2F("FMDV0Asame_post", "FMD vs V0A post cut;FMD;V0A;",2000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0Asame_post);
	 fFMDV0Csame = new TH2F("FMDV0Csame", "FMD vs V0C;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0Csame);
	 fFMDV0Csame_post = new TH2F("FMDV0Csame_post", "FMD vs V0C post cut;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
//	 fOutputList2->Add(fFMDV0Csame_post);
	 
	 fh2_ITS_acceptance=new TH2D("fh2_ITS_acceptance","fh2_ITS_acceptance",200,-10,10,200,-4,6);
//	 fOutputList2->Add(fh2_ITS_acceptance);
	 
	 fh2_SPD_multcorr=new TH2F("fh2_SPD_multcorr","fh2_SPD_multcorr",400,0,400,2000,0,2000);
//	 fOutputList2->Add(fh2_SPD_multcorr);
	 
	 fh2_SPDV0_multcorr=new TH2F("fh2_SPDV0_multcorr","fh2_SPDV0_multcorr",400,0,400,2000,0,2000);
//	 fOutputList2->Add(fh2_SPDV0_multcorr);
	 
	 fh2_SPDtrack_multcorr=new TH2F("fh2_SPDtrack_multcorr","fh2_SPDtrack_multcorr",400,0,400,400,0,400);
//	 fOutputList2->Add(fh2_SPDtrack_multcorr);
	 
	 fhtrackletsdphi=new TH1F("fhtrackletsdphi","dphi tracklets",100,-100,100);
//	 fOutputList2->Add(fhtrackletsdphi);
	 
	 
	 fh2_FMD_acceptance=new TH2D("fh2_FMD_acceptance","fh2_FMD_acceptance",200,-4,6,200,-10,10);
//	 fOutputList2->Add(fh2_FMD_acceptance);
	 
	 
	 fh2_FMD_eta_phi=new TH2D("fh2_FMD_eta_phi","fh2_FMD_eta_phi",200,-4,6,20,0,2*TMath::Pi());
//	 fOutputList2->Add(fh2_FMD_eta_phi);
	 
	 fhistfmdphiacc=new TH2D("fhistfmdphiacc","fhistfmdphiacc",200,-4,6,20,0,ncentmax);
//	 fOutputList2->Add(fhistfmdphiacc);
	 
	 fhFMDmultchannel=new TH2F("fhFMDmultchannel","fhFMDmultchannel",200,-4,6,100,0,100);
//	 fOutputList2->Add(fhFMDmultchannel);
	 
	 for(Int_t i=0;i<31;i++){
	   fhFMDmult_runbyrun_cside[i]=new TH2D(Form("fhFMDmult_runbyrun_cside_%d",i),Form("fhFMDmult_runbyrun_%d",i),200,-0.5,199.5,100,0,100);
	   // old     fOutputList2->Add(fhFMDmult_runbyrun_cside[i]);
	 }
	 for(Int_t i=0;i<65;i++){
	   fhFMDmult_runbyrun_aside[i]=new TH2D(Form("fhFMDmult_runbyrun_aside_%d",i),Form("fhFMDmult_runbyrun_%d",i),200,-0.5,199.5,100,0,100);
	   // old     fOutputList2->Add(fhFMDmult_runbyrun_aside[i]);
	 }
	 
	 const Int_t ifmdbin[4]={200,20,20,20};
	 fhistfmd=new AliTHn("fhistfmd","fhistfmd",1,4,ifmdbin);
	 fhistfmd->SetBinLimits(0,-4.,6.);
	 fhistfmd->SetBinLimits(1,0.,2*TMath::Pi());
	 if(fCentType=="Manual")fhistfmd->SetBinLimits(2,0.,200.);
	 else fhistfmd->SetBinLimits(2,0.,100.);
	 fhistfmd->SetBinLimits(3,-10.,10.);
	 fhistfmd->SetVarTitle(0,"eta");
	 fhistfmd->SetVarTitle(1,"phi");
	 fhistfmd->SetVarTitle(2,"centrality");
	 fhistfmd->SetVarTitle(3,"vzy");
// old	 fOutputList2->Add(fhistfmd);
	 
	 const Int_t iitsbin[4]={200,20,20,40};
	 const	 Double_t MinITS[4]={-4,0.,0.,-10.};
	 const	 Double_t MaxITS[4]={6,2*TMath::Pi(),100.,10.};
	 fhistits=new THnSparseF("fhistits","fhistits",4,iitsbin,MinITS,MaxITS);
	 // old  fOutputList2->Add(fhistits);
	 
	 const Double_t binning_etafmd[51]={
					    -3.4,-3.3,-3.2,-3.1,-3.0,
					    -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,
					    -1.9,-1.8,-1.7,
					    1.7,1.8,1.9,
					    2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
					    3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
					    4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
	 
	 const Int_t binfmdsec[5]={160,50,15,20,10};
	 fhSecFMD= new AliTHn("fhSecFMD","fhSecFMD",2,5,binfmdsec);
	 fhSecFMD->SetBinLimits(0,-4.025,3.975);
	 fhSecFMD->SetBinLimits(1,binning_etafmd);
	 fhSecFMD->SetBinLimits(2,binning_cent);
	 fhSecFMD->SetBinLimits(3,-0.55*TMath::Pi(),1.45*TMath::Pi());
	 fhSecFMD->SetBinLimits(4,binning_zvx);
	 fhSecFMD->SetVarTitle(0,"#Delta#eta");
	 fhSecFMD->SetVarTitle(1,"FMD Eta");
	 fhSecFMD->SetVarTitle(2,"centrality");
	 fhSecFMD->SetVarTitle(3,"#Delta#phi");
	 fhSecFMD->SetVarTitle(4,"z vertex");
// old	 fOutputList2->Add(fhSecFMD);
	 
	 if(fasso=="PID"){
	   for(Int_t i=0;i<6;i++){
	     fHistNsig[i]=new TH2D(Form("fHistNsig_%d",i),Form("HistNsig_%d",i), 160, 0., 8., 600, -30., 30);
//	     fOutputList2->Add(fHistNsig[i]);
	     fHistNsigcorr[i]=new TH2D(Form("fHistNsigcorr_%d",i),"fHistNsigcorr",500,-10,10,500,-10,10);
//	     fOutputList2->Add(fHistNsigcorr[i]);
	   }
	 }
	 
	 if(fasso=="Phi"){
	   fHistPhiDTPCNSig = new TH2D("fHistPhiDTPCNSig", "fHistPhiDTPCNSig", 150, 0.,15., 200, -10., 10);
	   fOutputList2->Add(fHistPhiDTPCNSig);
	   fHistPhiDTOFNSig = new TH2D("fHistPhiDTOFNSig", "fHistPhiDTOFNSig", 150, 0.,15., 200, -10., 10);
	   fOutputList2->Add(fHistPhiDTOFNSig);
	   fHistPhiDTPCTOFNSig = new TH2D("fHistPhiDTPCTOFNSig", "fHistPhiDTPCTOFNSig",150, 0., 15., 200, -10., 10);
	   fOutputList2->Add(fHistPhiDTPCTOFNSig);
	 }
	 Int_t nBins = 400;
	 Double_t mphiMin = 1.02 - 0.1;
	 Double_t mphiMax = 1.02 + 0.1;
	 
	 Int_t nCentralityBins = 20;
	 Double_t centBins1[16] = {0.,  1.,  2.,  3.,  4.,  5.,  10., 20.,
				   30., 40., 50., 60., 70., 80., 90., 100.0};
	 const Double_t *centralityBins = centBins1;

   Int_t nPtBinsV0 = 150;
   const Double_t PtBinsV0[12] = {0,   0.5, 0.75, 1.0, 1.5, 2.0,
                                  2.5, 3.0, 3.5,  4.0, 8.0, 15.0};

   Double_t mk0sMin = 0.5 - 0.1;
   Double_t mk0sMax = 0.5 + 0.1;
   Int_t netabins=4;
   const Int_t spBins[3] = {nBins, nPtBinsV0, nCentralityBins};
   const Int_t spBinsV0[4] = {nBins, nPtBinsV0, nCentralityBins,netabins};
   const Int_t spBinsBump[3] = {500, nPtBinsV0, nCentralityBins};
   // v0
   const Double_t spMink0s[4] = {mk0sMin, PtBinsV0[0], centralityBins[0],-0.8};
   const Double_t spMaxk0s[4] = {mk0sMax, PtBinsV0[11], centralityBins[15],0.8};
   Double_t mlambdaMin = 1.15 - 0.1;
   Double_t mlambdaMax = 1.15 + 0.1;
   const Double_t spMinLambda[4] = {mlambdaMin, PtBinsV0[0], centralityBins[0],-0.8};
   const Double_t spMaxLambda[4] = {mlambdaMax, PtBinsV0[11],centralityBins[15],0.8};
   const Double_t spMinBump[3] = {0, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxBump[3] = {2.5, PtBinsV0[11], centralityBins[15]};
   if(fasso=="V0"){
     fHistMass_K0s = new THnSparseF("fHistMass_K0s", "mass for K0s", 4, spBinsV0, spMink0s, spMaxk0s);
     fOutputList2->Add(fHistMass_K0s);
     fHistMass_K0s_MC = new THnSparseF("fHistMass_K0s_MC", "mass for K0s", 4, spBinsV0, spMink0s, spMaxk0s);
     fOutputList2->Add(fHistMass_K0s_MC);

     fHistMass_Lambda = new THnSparseF("fHistMass_Lambda", "mass for Lambda", 4,spBinsV0, spMinLambda, spMaxLambda);
     fOutputList2->Add(fHistMass_Lambda);
     fHistMass_Lambda_MC = new THnSparseF("fHistMass_Lambda_MC", "MC mass for Lambda", 4,spBinsV0, spMinLambda, spMaxLambda);
     fOutputList2->Add(fHistMass_Lambda_MC);

     fHistMass_ALambda = new THnSparseF("fHistMass_ALambda", "mass for Anti Lambda", 4, spBinsV0,spMinLambda, spMaxLambda);
     fOutputList2->Add(fHistMass_ALambda);
     fHistMass_ALambda_MC = new THnSparseF("fHistMass_ALambda_MC", "mass for Anti Lambda", 4, spBinsV0,spMinLambda, spMaxLambda);
     fOutputList2->Add(fHistMass_ALambda_MC);

     fHistMass_bumpcorr =new TH2D("fHistMass_bumpcorr", "mass for Lambda bump correlation", 400,mlambdaMin, mlambdaMax, 1000, 0, 1);
     fOutputList2->Add(fHistMass_bumpcorr);

     const Int_t spBinsQAV0[4] = {40, 72, 7, 20};
     const Double_t spMinV0QA[4] = {-1., 0, -0.5, 0.};
     const Double_t spMaxV0QA[4] = {1., TMath::TwoPi(), 6.5, 100.0};
     fHist_V0QA = new THnSparseF("fHist_V0QA", "QA for V0 particle", 4, spBinsQAV0, spMinV0QA, spMaxV0QA);
     fOutputList2->Add(fHist_V0QA);

     hv0dcharge = new TH1D("hv0dcharge", "hv0dcharge", 3, -0.5, 2.5);
     fOutputList2->Add(hv0dcharge);
     for (Int_t i = 0; i < 6; i++) {
       fHistPosNsig[i] = new TH2D(Form("fHistPosNsig_%d", i), "fHistPosNsig", 160, 0., 8., 600, -30., 30);
       fHistNegNsig[i] = new TH2D(Form("fHistNegNsig_%d", i), "fHistNegNsig", 160, 0., 8., 600, -30., 30);
       fOutputList2->Add(fHistPosNsig[i]);
       fOutputList2->Add(fHistNegNsig[i]);
       fHistPosNsigQA[i] = new TH2D(Form("fHistPosNsigQA_%d", i), "fHistPosNsigQA",160, 0., 8., 600, -30., 30);
       fOutputList2->Add(fHistPosNsigQA[i]);
     }
     for(Int_t i=0;i<3;i++){
       fh3NegNsig[i]=new TH3D(Form("fh3NegNsig_%d",i),Form("fh3NegNsig_%d",i),40,0,8,200,-10,10,200,-10,10);
       fh3PosNsig[i]=new TH3D(Form("fh3PosNsig_%d",i),Form("fh3PosNsig_%d",i),40,0,8,200,-10,10,200,-10,10);
       fOutputList2->Add(fh3NegNsig[i]);
       fOutputList2->Add(fh3PosNsig[i]);
     }
     for (Int_t i = 0; i < 6; i++) {
       fHist_AP[i] = new TH2D(Form("fHist_AP_%d", i), Form("fHist_AP_%d", i), 200, -1, 1, 200, 0, 0.4);
       fOutputList2->Add(fHist_AP[i]);
     }

   }
   if(fasso=="Cascade"){
     // QA Plot for Cascade
     Double_t mxiMin = 1.3 - 0.1;
     Double_t mxiMax = 1.3 + 0.1;
     Double_t momegaMin = 1.65 - 0.1;
     Double_t momegaMax = 1.65 + 0.1;
     const Double_t spMinXi[3] = {mxiMin, PtBinsV0[0], centralityBins[0]};
     const Double_t spMaxXi[3] = {mxiMax, PtBinsV0[11], centralityBins[15]};
     const Double_t spMinOmega[3] = {momegaMin, PtBinsV0[0], centralityBins[0]};
     const Double_t spMaxOmega[3] = {momegaMax, PtBinsV0[11], centralityBins[15]};
     fHistMassXiMinus = new THnSparseF("fHistMassXiMinus", "mass for Xi-", 3,
     spBins, spMinXi, spMaxXi);
     fHistMassXiPlus = new THnSparseF("fHistMassXiPlus", "mass for Xi+", 3, spBins,
     spMinXi, spMaxXi);
     fHistMassOmegaMinus = new THnSparseF("fHistMassOmegaMinus", "mass for Omega-",
     3, spBins, spMinOmega, spMaxOmega);
     fHistMassOmegaPlus = new THnSparseF("fHistMassOmegaPlus", "mass for Omega+",
     3, spBins, spMinOmega, spMaxOmega);
     fOutputList2->Add(fHistMassXiMinus);
     fOutputList2->Add(fHistMassXiPlus);
     fOutputList2->Add(fHistMassOmegaMinus);
     fOutputList2->Add(fHistMassOmegaPlus);

     const Int_t spBinsQACasc[5] = {20, 40, 72, 7, 20};
     const Double_t spMinCascQA[5] = {0, -1., 0, -0.5, 0.};
     const Double_t spMaxCascQA[5] = {10, 1., TMath::TwoPi(), 6.5, 100.0};

     fHist_CascadeQA = new THnSparseF("fHist_CascadeQA", "QA for Cascade particle", 5, spBinsQACasc, spMinCascQA, spMaxCascQA);
     fOutputList2->Add(fHist_CascadeQA);
   }
   if(fasso=="Phi"){
     // QA Plot for Phi meson
     const Double_t spMinPhi[3] = {mphiMin, PtBinsV0[0], centralityBins[0]};
     const Double_t spMaxPhi[3] = {mphiMax, PtBinsV0[11], centralityBins[15]};
     // Phimeson
     fHistMass_PhiMeson =
     new THnSparseF("fHistMass_PhiMeson", "mass for phi meson", 3, spBins,
     spMinPhi, spMaxPhi);
     fOutputList2->Add(fHistMass_PhiMeson);

     fHistMass_PhiMeson_MIX = new THnSparseF("fHistMass_PhiMeson_MIX",
     "mass for phi meson of mixed events",
     3, spBins, spMinPhi, spMaxPhi);
     fOutputList2->Add(fHistMass_PhiMeson_MIX);

     const Int_t spBinsQA[4] = {40, 72, 7, 20};
     const Double_t spMinPhiQA[4] = {-1., 0, -0.5, 0.};
     const Double_t spMaxPhiQA[4] = {1., TMath::TwoPi(), 6.5, 100.0};
     fHist_PhiQA = new THnSparseF("fHist_PhiQA", "QA for Phimeson", 4, spBinsQA, spMinPhiQA, spMaxPhiQA);
     fOutputList2->Add(fHist_PhiQA);
   }
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefineCorrOutput() {

   Double_t binning_pt_assoc[12] = {0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
                                    2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
   Double_t binning_pt_lead[12] = {0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
                                   2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
   Double_t binning_cent[12] = {0., 5.,  10., 20.,
                                30., 40., 50., 60., 70., 80., 90., 100.1};
   Double_t binning_cent_HMPP[12] = {0., 0.01,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,1.};

   Double_t binning_deta[49] = {
       -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5,
       -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
       -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5,
       0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,
       1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4};

   Double_t binning_dphi[73] = {
       -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464,
       -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865,
       -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266,
       0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332,
       0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931,
       1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530,
       1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129,
       2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727,
       2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326,
       3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925,
       3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524,
       4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123,
       4.712389};

   Int_t nCFStepstrig=1;
/*
   if(fasso=="hadron")  nCFStepstrig=1;
   else  if (fasso == "V0" || fasso == "Phi")    nCFStepstrig = 7;
   else  if (fasso == "Cascade")    nCFStepstrig = 6;
   else  if(fasso=="PID")    nCFStepstrig=3;
*/
   Double_t binning_eta_vzero[11]={-3.7,-3.2,-2.7,-2.2,-1.7,0.,2.8,3.4,3.9,4.5,5.1};
   Double_t binning_phi_vzero[9]={0.,0.7853,1.5707,2.3561,3.1415,3.9269,4.7123,5.4977,6.2831};
   //Bins for FMD
 
  const Double_t binning_etafmd[2]={1.7,3.3};

   const Double_t binning_etafmdc[2]={-3.3, -1.7};

   Double_t binning_zvx[11] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
     
     Double_t binning_pt_lead_trig[] = {0.5, 1., 1.5, 2., 3., 5.};

     //Double_t binning_cent_trig[8] = {0., 5.,  10., 20., 40., 60., 70., 100.1};
     Double_t binning_cent_trig[] = {0., 10., 60., 100.1};
     Int_t nbinning_cent_trig = sizeof(binning_cent_trig) /sizeof(Double_t) - 1;
     Double_t binning_cent_trig_PbPb[9] = {0., 10.,  20., 30., 40., 50.,60.,70.,80.};
     Double_t binning_mult_trig[10]={0,20,40,60,80,100,120,140,160,200};

     const Int_t nEvtVarsFMD = 3;
     Int_t ntpcpt = sizeof(binning_pt_lead_trig)/sizeof(Double_t) - 1;
     
     const Int_t iEvtBinFMD[] = {ntpcpt,10,12};
     
     fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsFMD, iEvtBinFMD);
     
     fHistTriggerTrack->SetBinLimits(0, binning_pt_lead_trig);
     fHistTriggerTrack->SetBinLimits(1, -10.,10.);
     fHistTriggerTrack->SetBinLimits(2, 0.,12.);
     fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrack->SetVarTitle(1, "zvertex");
     fHistTriggerTrack->SetVarTitle(2, "Random Number");
 
     const Int_t nEvtVars = 2;
     const Int_t iEvtBin[2] = {ntpcpt, 12};
    
     fHistTriggerTrackMix = new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVars, iEvtBin);
     fHistTriggerTrackMix->SetBinLimits(0, binning_pt_lead_trig);
     fHistTriggerTrackMix->SetBinLimits(1, 0.,12.);
     fHistTriggerTrackMix->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrackMix->SetVarTitle(1, "Random Number");
     
   fOutputList1->Add(fHistTriggerTrack);
   fOutputList1->Add(fHistTriggerTrackMix);
   
   const Int_t nTrackVars = 5;
   const Int_t iTrackBin[5] = {48, 11, 11, 15, 72};
   //////////////////////////////////////////
   //Containers two particle correlation
   //////////////////////////////////////////

   Int_t nCFSteps = 1;
/*
   if(fasso=="hadron")    nCFSteps=1;
   else if (fasso == "V0" || fasso == "Phi")    nCFSteps = 7;
   else  if (fasso == "Cascade")    nCFSteps = 6;
   else    if(fasso=="PID")    nCFSteps=3;
*/
   Double_t binning_dphi_vzero[9]={-1.178097,-0.392699,0.392699,1.178097,1.963495,2.748893,3.534291,4.319689,5.105088};
   if(fAnaMode=="TPCTPC") {
     Double_t binning_deta_tpctpc[33] = {-1.6, -1.5,
					 -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
					 -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5,
					 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,
					 1.6};
     const Double_t binning_cent_tpctpc[8]={0.,5.,10.,20.,40.,60.,70.,100.1};
     const Int_t iTrackBin_TPCTPC[6] = {32, 11, 11, 7, 72, 10};
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 6, iTrackBin_TPCTPC);
     fHistReconstTrack->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrack->SetBinLimits(1, binning_pt_assoc);
     fHistReconstTrack->SetBinLimits(2, binning_pt_lead);
     fHistReconstTrack->SetBinLimits(3, binning_cent_tpctpc);
     fHistReconstTrack->SetBinLimits(4, binning_dphi);
     fHistReconstTrack->SetBinLimits(5, -10,10);
     fHistReconstTrack->SetVarTitle(0, "#Delta#eta");
     fHistReconstTrack->SetVarTitle(1, "p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(2, "leading p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(3, "centrality");
     fHistReconstTrack->SetVarTitle(4, "#Delta#phi");
     fHistReconstTrack->SetVarTitle(5, "Vz");
     fHistReconstTrackMix =  new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 6, iTrackBin_TPCTPC);
     fHistReconstTrackMix->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrackMix->SetBinLimits(1, binning_pt_assoc);
     fHistReconstTrackMix->SetBinLimits(2, binning_pt_lead);
     fHistReconstTrackMix->SetBinLimits(3, binning_cent_tpctpc);
     fHistReconstTrackMix->SetBinLimits(4, binning_dphi);
     fHistReconstTrackMix->SetBinLimits(5, -10,10);
     fHistReconstTrackMix->SetVarTitle(0, "#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1, "p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2, "leading p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(3, "centrality");
     fHistReconstTrackMix->SetVarTitle(4, "#Delta#phi");
     fHistReconstTrackMix->SetVarTitle(5, "Vz");
   }else if (fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
 /* 
    Double_t binning_deta_tpctpc[] = {-1.6, -1.5,
                                         -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
                                         -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5,
                                         0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,
                                         1.6};
 */
/*
    Double_t binning_deta_tpctpc[] = {-1.6, -1.4, -1.2, -1.0, -0.8,
                                        -0.6, -0.4, -0.2, 0,  0.2, 0.4,
                                         0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
 */
    Double_t binning_deta_tpctpc[] = {-1.6, -1.3, -1.1, -0.9, -0.7 ,-0.5, -0.3, -0.1,  0.1,  0.3,  0.5, 0.7, 0.9, 1.1, 1.3, 1.6};
 

     const Double_t binning_etafmd_tpcfmda[]={
     1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,
     4.1,4.3,4.5,4.7};

/*
     Double_t binning_detaFMDTPC[]={-5.6,-5.4,-5.2,-5.0,
				      -4.8,-4.6,-4.4,-4.2,-4.,
				      -3.8,-3.6,-3.4,-3.2,-3.,
				      -2.8,-2.6,-2.4,-2.2,-2.,
				      -1.8,-1.6,-1.4,-1.2};
*/
   Double_t binning_detaFMDTPC[]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1.,-0.8};

/*
   Double_t binning_detaFMDCTPC[]={ 0.9,1.1,1.3,1.5,1.7,1.9,
				       2.1,2.3,2.5,2.7,2.9,
				       3.1,3.3,3.5,3.7,3.9,4.1};
*/

   Double_t binning_detaFMDCTPC[]={ 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};



/*
     Double_t binning_dphi_reduce[37] = {
       -1.570796,  -1.396263,  -1.221730,
       -1.047198,   -0.872665, -0.698132,
       -0.523599,  -0.349066, -0.174533, 
       0.0,       0.174533,   0.349066, 
       0.523599,  0.698132,   0.872665, 
       1.047198,  1.221730,   1.396263,  
       1.570796,  1.745329,   1.919862,  
       2.094395,  2.268928,   2.443461,  
       2.617994,  2.792527,   2.967060,  
       3.141593,  3.316126,   3.490659,  
       3.665191,  3.839724,   4.014257,  
       4.188790,  4.363323,   4.537856,  
       4.712389};
*/


     Double_t binning_dphi_reduce[] = {
       -1.570796,  -1.221730,
       -0.872665, -0.523599,
       -0.174533,  0.174533,  0.523599, 
       0.872665,  1.221730,   1.570796,  
       1.919862,  2.268928,   2.617994,  
       2.967060,  3.316126,  3.665191,  
       4.014257,  4.363323,   4.712389};
    Int_t nbinning_dphi_reduce = sizeof(binning_dphi_reduce)/sizeof(Double_t) - 1;
    
 //   Double_t binning_dphi_reduce_TPCTPC[] = {-0.5*TMath::Pi(), -0.3*TMath::Pi(), -0.1*TMath::Pi(), 0.1*TMath::Pi(), 0.3*TMath::Pi(), 0.5*TMath::Pi(),1.5*TMath::Pi()};
 
/*   
    Double_t binning_dphi_reduce_TPCTPC[] = {
       -1.570796,  -1.221730, 
       -0.872665, -0.523599, 
       -0.174533,  0.174533,  0.523599, 
       0.872665,  1.221730,   1.570796,  
       1.919862,  2.268928,   2.617994,  
       2.967060,  3.316126,  3.665191,  
       4.014257,  4.363323,   4.712389};
*/

/*
    Double_t Example[] = {
      -1.570796,  -1.221730,
      -0.872665, -0.523599, 
      -0.174533,  0.174533,  0.523599, 
      0.872665, 1.221730, 1.570796};
*/
/*
    Double_t Example[] = { 
     -1.570796,  -1.396263,  -1.221730,
     -1.047198,   -0.872665, -0.698132,
     -0.523599,  -0.349066, -0.174533,
     0.0,       0.174533,   0.349066,
     0.523599,  0.698132,   0.872665,
     1.047198,  1.221730,   1.396263,
     1.570796, 4.712389};
*/

    //Double_t binning_dphi_reduce_TPCTPC_FMDC[] = {-0.5*TMath::Pi(), -0.4*TMath::Pi(), -0.3*TMath::Pi(), -0.2*TMath::Pi(), -0.1*TMath::Pi(), 0.1*TMath::Pi(), 0.2*TMath::Pi(), 0.3*TMath::Pi(), 0.4*TMath::Pi(), 0.5*TMath::Pi(),1.5*TMath::Pi()  };

    //Double_t binning_dphi_reduce_TPCTPC_FMDA[] = { -0.5*TMath::Pi(), -0.3*TMath::Pi(), -0.1*TMath::Pi(), 0.1*TMath::Pi(), 0.3*TMath::Pi(), 0.5*TMath::Pi(),1.5*TMath::Pi()};

/*
    Double_t binning_dphi_reduce_TPCTPC[] = {
      -1.570796,  -1.396263,  -1.221730,
       -1.047198,   -0.872665, -0.698132,
       -0.523599,  -0.349066, -0.174533, 
       0.0,       0.174533,   0.349066, 
       0.523599,  0.698132,   0.872665, 
       1.047198,  1.221730,   1.396263,  
       1.570796,  1.745329,   1.919862,  
       2.094395,  2.268928,   2.443461,  
       2.617994,  2.792527,   2.967060,  
       3.141593,  3.316126,   3.490659,  
       3.665191,  3.839724,   4.014257,  
       4.188790,  4.363323,   4.537856,  
       4.712389};
*/ 

	
     Int_t ndetatpcfmd;
     Int_t nbinning_dphi_reduce_TPCTPC;            
     if(fAnaMode=="TPCFMD") {
       ndetatpcfmd=sizeof(binning_detaFMDTPC)/sizeof(Double_t) - 1;
       //nbinning_dphi_reduce_TPCTPC = sizeof(binning_dphi_reduce_TPCTPC_FMDA)/sizeof(Double_t) - 1;
     }else{
       ndetatpcfmd = sizeof(binning_detaFMDCTPC)/sizeof(Double_t) - 1;
       //nbinning_dphi_reduce_TPCTPC = sizeof(binning_dphi_reduce_TPCTPC_FMDC)/sizeof(Double_t) - 1;
     }
     
     Double_t binning_pt_fmdtpc_asso[] = {0.5, 1., 1.5, 2., 3., 5.};
     Double_t binning_pt_fmdtpc[]      = {0.5, 1., 1.5, 2., 3., 5.};
   

     Int_t ndetatpctpc = sizeof(binning_deta_tpctpc)/sizeof(Double_t)-1;            
     Int_t ntpcpt_asso = sizeof(binning_pt_fmdtpc_asso)/sizeof(Double_t)-1;
     Int_t ntpcpt      = sizeof(binning_pt_fmdtpc)/sizeof(Double_t)-1;

     fh2_pt_trig_asso = new TH2D("fh2_pt_trig_asso","",ntpcpt,binning_pt_fmdtpc,ntpcpt_asso,binning_pt_fmdtpc_asso);

     //const Double_t binning_cent_fmdfmd[8]={0.,5.,10.,20.,40.,60.,70.,100.1};
     const Double_t binning_cent_fmdfmd[]={0.,10., 60.,100.1};
     const Double_t binning_cent_fmdfmd_PbPb[9] = {0., 10.,  20., 30., 40., 50.,60.,70.,80.};
     
     const Double_t binning_cent_fmdfmd_HMPP[9]={0.,0.01,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

     Double_t binning_mult[10]={0,20,40,60,80,100,120,140,160,200};

     Int_t ncentbin=sizeof(binning_cent_fmdfmd)/sizeof(Double_t) - 1;

     Int_t nCFStepstpcfmd=1;
     //const Int_t iTrackBin_tpcfmd[8]={ndetatpctpc, ntpcpt_asso, ntpcpt, ncentbin, nbinning_dphi_reduce_TPCTPC, 10, ndetatpcfmd, nbinning_dphi_reduce};
     const Int_t iTrackBin_tpcfmd[7]={ndetatpctpc, ntpcpt_asso, ntpcpt, 10, ndetatpcfmd, nbinning_dphi_reduce, 12};

     //fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFStepstpcfmd, 8, iTrackBin_tpcfmd);
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFStepstpcfmd, 7, iTrackBin_tpcfmd);
     
     //fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, 8,iTrackBin_tpcfmd);
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, 7, iTrackBin_tpcfmd);

     if(fAnaMode=="TPCFMD") {
       fHistReconstTrack->SetBinLimits(4,binning_detaFMDTPC);
//       fHistReconstTrack->SetBinLimits(3,binning_dphi_reduce_TPCTPC_FMDA);
       //Mixed Events
       fHistReconstTrackMix->SetBinLimits(4,binning_detaFMDTPC);
//       fHistReconstTrackMix->SetBinLimits(3,binning_dphi_reduce_TPCTPC_FMDA);
     } else if(fAnaMode=="TPCFMDC") { 
       fHistReconstTrack->SetBinLimits(4,binning_detaFMDCTPC);
//       fHistReconstTrack->SetBinLimits(3,binning_dphi_reduce_TPCTPC_FMDC);
       //Mixed Events
       fHistReconstTrackMix->SetBinLimits(4,binning_detaFMDCTPC);
       //fHistReconstTrackMix->SetBinLimits(3,binning_dphi_reduce_TPCTPC_FMDC);
     }
     
     fHistReconstTrack->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrack->SetBinLimits(1, binning_pt_fmdtpc_asso);
     fHistReconstTrack->SetBinLimits(2, binning_pt_fmdtpc);
     
     fHistReconstTrack->SetBinLimits(3,-10.,10.);
     fHistReconstTrack->SetBinLimits(5,binning_dphi_reduce);
     fHistReconstTrack->SetBinLimits(6, 0., 12.);
     
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
     fHistReconstTrack->SetVarTitle(1,"p_{T}_{TPC-asso} GeV/c");
     fHistReconstTrack->SetVarTitle(2,"p_{T}_{TPC-trig} GeV/c");
//     fHistReconstTrack->SetVarTitle(3,"#Delta#phi_{TPC-TPC}");
     fHistReconstTrack->SetVarTitle(3,"z vertex");
     fHistReconstTrack->SetVarTitle(4,"#Delta#eta_{TPC-FMD}");
     fHistReconstTrack->SetVarTitle(5,"#Delta#phi_{TPC-FMD}");
     fHistReconstTrack->SetVarTitle(6,"Random Number");

     fHistReconstTrackMix->SetBinLimits(0,binning_deta_tpctpc);
     fHistReconstTrackMix->SetBinLimits(1,binning_pt_fmdtpc_asso);
     fHistReconstTrackMix->SetBinLimits(2,binning_pt_fmdtpc);
     

     fHistReconstTrackMix->SetBinLimits(3,-10.,10.);
     fHistReconstTrackMix->SetBinLimits(5,binning_dphi_reduce);
     fHistReconstTrackMix->SetBinLimits(6, 0.,12.);

     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
     fHistReconstTrackMix->SetVarTitle(1,"p_{T}_{TPC-asso} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2,"p_{T}_{TPC-trig} GeV/c");
//     fHistReconstTrackMix->SetVarTitle(3,"#Delta#phi_{TPC-TPC}");
     fHistReconstTrackMix->SetVarTitle(3,"z vertex");
     fHistReconstTrackMix->SetVarTitle(4,"#Delta#eta_{TPC-FMD}");
     fHistReconstTrackMix->SetVarTitle(5,"#Delta#phi_{TPC-FMD}");
     fHistReconstTrackMix->SetVarTitle(6,"Random Number");
   }
   fOutputList1->Add(fHistReconstTrack);
   fOutputList1->Add(fHistReconstTrackMix);
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::UserExec(Option_t *) {

   DumpTObjTable("Start analysis");
   AliAnalysisManager *mgr        = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inEvMain = (AliInputEventHandler *)(mgr->GetInputEventHandler());
   
   if (!inEvMain)    return;
   
//   fPIDResponse = inEvMain->GetPIDResponse();
//   if (!fPIDResponse)    return;

   if(!fDataType){
     AliMCEventHandler* mctruth = (AliMCEventHandler*)(mgr->GetMCtruthEventHandler());
     if(!mctruth)       return;
     mcEvent=mctruth->MCEvent();//AliMCEvent
   }

   fEvent = dynamic_cast<AliAODEvent *>(inEvMain->GetEvent());
   if (!fEvent) {
     AliWarning("ERROR: fEvent not available \n");
     return;
   }

   
   fHist_Stat->Fill(0);
   if(fcollisiontype=="pPb" || fcollisiontype=="PP" || fcollisiontype=="PbPb"){
     //   if(fcollisiontype=="pPb"){
   if (!fEventCuts.AcceptEvent(fEvent)) {
     PostData(1, fOutputList);
     PostData(2, fOutputList1);
     PostData(3, fOutputList2);
     return;
   }
   }else if(fcollisiontype=="HMPP"){
     UInt_t maskIsSelected = inEvMain->IsEventSelected();
     Bool_t isSelected     = kFALSE;
     if(fCentType=="V0M")  isSelected = ((maskIsSelected & AliVEvent::kHighMultV0)== AliVEvent::kHighMultV0);
     //	 else if(fCentType=="Manual") isSelected = ((maskIsSelected & AliVEvent::kHighMultSPD)== AliVEvent::kHighMultSPD);
     else  isSelected = ((maskIsSelected & AliVEvent::kHighMultSPD)== AliVEvent::kHighMultSPD);
     if (!isSelected) {
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     }
   }else if (fcollisiontype.Contains("MBPP")){
     UInt_t maskIsSelected = inEvMain->IsEventSelected();
     Bool_t isSelected     = kFALSE;
     isSelected = ((maskIsSelected & AliVEvent::kINT7)== AliVEvent::kINT7);//Both for data and 
     
     if (!isSelected) {
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     }
   }
   
   /*
	 //Vertex
	 Bool_t IsGoodVtx=kFALSE;
	 Bool_t IsValidVtx=kFALSE;
	 const AliVVertex* spdVtx   = fEvent->GetPrimaryVertexSPD() ;
	 if( spdVtx ){
	   if( spdVtx->GetNContributors() > 0.5 ) IsValidVtx = kTRUE;
m	   auto fZ = spdVtx->GetZ();
	   auto  zbin = binZ.FindBin(fZ) -1;
	   if( fabs(fZ) < 10 && !(zbin < 0 )) IsGoodVtx = kTRUE;
	 }
	 
	 if(!IsGoodVtx || !IsValidVtx){
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
	 }
	 */
	 //Pile up
   if(fcollisiontype.Contains("MBPP") || fcollisiontype.Contains("HMPP")){
     Bool_t IsNotPileup = kTRUE;
     if( !fEvent->IsPileupFromSPDInMultBins() ) IsNotPileup = kTRUE;
     if(!IsNotPileup){
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     }
     fHist_Stat->Fill(8);
   }

	 
   
   fHist_Stat->Fill(1);      

   AliMultSelection *multSelection =    (AliMultSelection *)fEvent->FindListObject("MultSelection");
   if(!multSelection) return;
   fHist_Stat->Fill(2);
   
   //Pileu rejection by MultSelection
   if(fcollisiontype.Contains("HMPP") || fcollisiontype.Contains("MBPP")){
     Bool_t IsSelectedFromAliMultSelection=kFALSE;
     if( multSelection->GetThisEventIsNotPileup() &&
	 multSelection->GetThisEventIsNotPileupInMultBins() &&
	 multSelection->GetThisEventHasNoInconsistentVertices() &&
	 multSelection->GetThisEventPassesTrackletVsCluster() ){
       IsSelectedFromAliMultSelection = kTRUE;
     }
     if(!IsSelectedFromAliMultSelection) {
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     }
     fHist_Stat->Fill(9);
   }

   
   lPrimaryBestVtx = fEvent->GetPrimaryVertex();
   if(fcollisiontype.Contains("HMPP") || fcollisiontype.Contains("MBPP")){
	 Int_t nTracksPrim = lPrimaryBestVtx->GetNContributors();
	 if (nTracksPrim < 0.5) {
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
	 }
	 fHist_Stat->Fill(10);
   }
   
   if ((TMath::Abs(lPrimaryBestVtx->GetZ())) >= fZVertex)   {
     PostData(1, fOutputList);
     PostData(2, fOutputList1);
     PostData(3, fOutputList2);
     return;
   }
   
   tPrimaryVtxPosition[0] = lPrimaryBestVtx->GetX();
   tPrimaryVtxPosition[1] = lPrimaryBestVtx->GetY();
   tPrimaryVtxPosition[2] = lPrimaryBestVtx->GetZ();

   fPrimaryZVtx = lPrimaryBestVtx->GetZ();
   fHist_Stat->Fill(3);

   bSign = 0.;
   bSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

   AliVMultiplicity *tracklets = ((AliAODEvent*)fEvent)->GetTracklets();
   if (!tracklets) return;
   Int_t nTracklets = tracklets->GetNumberOfTracklets();
   
   // Multiplicity Object
   //if(fcollisiontype=="pPb" || fcollisiontype=="PP"){
     if(fCentType=="Manual"){
       //       AliVMultiplicity *tracklets = ((AliAODEvent*)fEvent)->GetTracklets();
       //    if (!tracklets) return;
       //       Int_t nTracklets = tracklets->GetNumberOfTracklets();
       Int_t nITScluster= tracklets->GetNumberOfITSClusters(0)+tracklets->GetNumberOfITSClusters(1);
       Int_t nTracks = fEvent->GetNumberOfTracks();
       fh2_SPD_multcorr->Fill(nTracklets,nITScluster);
       fh2_SPDtrack_multcorr->Fill(nTracklets,nTracks);
       lCentrality=nTracklets;

     }else{
       lCentrality = multSelection->GetMultiplicityPercentile(fCentType);
       Int_t qual = multSelection->GetEvSelCode();
       if (qual == 199)  lCentrality = -999;
       if (lCentrality < 0. || lCentrality > 100. - 0.0000001)   return;
     }

     Double_t *CentBins = fCentBins;
     poolmin = CentBins[0];
     poolmax = CentBins[fNCentBins];
     fHist_Stat->Fill(4);

   fHistCentV0vsTrackletsbefore->Fill(lCentrality,nTracklets);
   
   
   fHistCentrality_beforecut->Fill(lCentrality);
   
   //   fHist_Stat->Fill(5);
   DumpTObjTable("After event selection");
   MakeAna();
      
   PostData(1, fOutputList);
   PostData(2, fOutputList1);
   PostData(3, fOutputList2);
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::Terminate(Option_t *) {
   //  AliInfo(Form("Number of Correlation
   DumpTObjTable("End of the analysis");
   Printf("Entries======================%d",fNEntries);
   if (fPoolMgr)    delete fPoolMgr;   // PoolMgr->ClearPools();
   if (fPoolMgr1)    delete fPoolMgr1; // fPoolMgr1->ClearPools();
 }

void AliAnalysisTaskSEpPbCorrelationsJetV2::MakeAna() {

   DumpTObjTable("start correlation analysis");
   TObjArray *selectedTracksLeading = new TObjArray;
   selectedTracksLeading->SetOwner(kTRUE);
   TObjArray *selectedTracksAssociated = new TObjArray;
   selectedTracksAssociated->SetOwner(kTRUE);
   
   fvzero = fEvent->GetVZEROData();

   Float_t nFMD_fwd_hits=0;
   Float_t nFMD_bwd_hits=0;
   Float_t nFMD_fwdV0acc_hits=0;
   Float_t nFMD_bwdV0acc_hits=0;
   AliAODForwardMult*aodForward=static_cast<AliAODForwardMult*>(fEvent->FindListObject("Forward"));

   // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
   Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
   const TH2D& d2Ndetadphi = aodForward->GetHistogram();
   TH1*hphiacceptance=aodForward->GetPhiAcceptance();
   Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
   Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
   Double_t pt = 0;

   if(fAnaMode=="TPCFMD"||fAnaMode=="TPCFMDC"||fAnaMode=="ITSFMD"||fAnaMode=="ITSFMDC"||fAnaMode.Contains("FMDFMD")){
     for (Int_t iEta = 1; iEta <= nEta; iEta++) {
       Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
       if (!valid) {
	 continue;
       }
       
       Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
       Float_t phiacc=hphiacceptance->GetBinContent(iEta);
       fhistfmdphiacc->Fill(eta,lCentrality,phiacc);
       for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
	 // Bin content is most likely number of particles!
	 Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	 
	 Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
	 fh2_FMD_acceptance->Fill(eta,tPrimaryVtxPosition[2],mostProbableN);
	 //Float_t corrfactor=fhcorr[ivzbin-1]->GetBinContent(iEta,iPhi);
	 
	 if (mostProbableN > 0) {
	   if(eta>0){
	     nFMD_fwd_hits+=mostProbableN;
	     if(2.8<eta && eta<5.03) nFMD_fwdV0acc_hits+=mostProbableN;
	   }else{
	     nFMD_bwd_hits+=mostProbableN;
	     if(-3.4<eta && eta<-2.01) nFMD_bwdV0acc_hits+=mostProbableN;
	   }
	 }

          
	 
	 if(fmakehole){
	   if((eta>-2.9 && eta<-2.7) && (5*2*TMath::Pi()/20.<phi && 7*2*TMath::Pi()/20.>phi)) continue;
	   if((eta>-2.7 && eta<-2.5) && (1*2*TMath::Pi()/20.<phi && 2*2*TMath::Pi()/20.>phi)) continue;
	   if((eta>-2.1 && eta<-1.9) && (17*2*TMath::Pi()/20.<phi && 20*2*TMath::Pi()/20.>phi)) continue;
	 }

	 if (mostProbableN > 0) {
	   //if(eta>0){
	   if(eta < 4.8 && eta > 2){
	     Int_t nfmdetabin1=frefetaa->FindBin(eta);
	     Double_t runbin1=ConvertRunNumber(fEvent->GetRunNumber());
	     //	     fhFMDmult_runbyrun_aside[nfmdetabin1-1]->Fill(runbin1,mostProbableN);
	     
	     if(fAnaMode=="TPCFMD" || fAnaMode=="ITSFMD") selectedTracksAssociated->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));			
	     else if(fAnaMode=="FMDFMD") selectedTracksLeading->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));	
	     else if(fAnaMode=="FMDFMD_Ctrig")	       selectedTracksAssociated->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));			
	     
	   //}else if(eta<0){
	   }else if( eta < -1.8 && eta > -3.2){
	     Int_t nfmdetabin=frefetac->FindBin(eta);
	     Double_t runbin=ConvertRunNumber(fEvent->GetRunNumber());
	     //	     fhFMDmult_runbyrun_cside[nfmdetabin-1]->Fill(runbin,mostProbableN);
	     
	     if(fAnaMode=="TPCFMDC" || fAnaMode=="ITSFMDC" ||fAnaMode=="FMDFMD") selectedTracksAssociated->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));
	     else if(fAnaMode=="FMDFMD_Ctrig") selectedTracksLeading->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));	
	   }
	   fhFMDmultchannel->Fill(eta,mostProbableN);
	   
	   Double_t cont[4]={eta,phi,lCentrality,fPrimaryZVtx};
	   fhistfmd->Fill(cont,0,mostProbableN);
	   fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	 }
       }
     }
     delete hphiacceptance;
     
     //events cuts
     
     if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0){
       selectedTracksLeading->Clear();
       delete selectedTracksLeading;
       selectedTracksAssociated->Clear();
       delete selectedTracksAssociated;
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     } //events cuts
      
     
     fHist_Stat->Fill(5);
     fHistCentrality_beforeFMDMulcut->Fill(lCentrality);  
    
     DumpTObjTable("End of fill fmd tracks");
     
     Float_t nV0A_hits = fvzero->GetMTotV0A();
     Float_t nV0C_hits = fvzero->GetMTotV0C();
     
     fFMDV0->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
     fFMDV0A->Fill(nFMD_fwd_hits, nV0A_hits);
     fFMDV0C->Fill(nFMD_bwd_hits, nV0C_hits);
     
     fFMDV0same->Fill(nFMD_bwdV0acc_hits + nFMD_fwdV0acc_hits, nV0C_hits + nV0A_hits);
     fFMDV0Asame->Fill(nFMD_fwdV0acc_hits, nV0A_hits);
     fFMDV0Csame->Fill(nFMD_bwdV0acc_hits, nV0C_hits);
     
     if(fCentType=="Manual")fh2_SPDV0_multcorr->Fill(lCentrality,nV0C_hits+nV0A_hits);
       if(fFMDcut){
       Double_t FMDcutapar0=0.;
       Double_t FMDcutapar1=0.;
       Double_t FMDcutcpar0=0.;
       Double_t FMDcutcpar1=0.;
       switch(fFMDcutmode){
       case 1:
	 FMDcutapar0=1.3;
	 FMDcutapar1=200;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=200;
       	 break;
       case 2:
	 FMDcutapar0=1.3;
	 FMDcutapar1=600;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=600;
	 break;
       case 3:
	 FMDcutapar0=1.5;
	 FMDcutapar1=100;
	 FMDcutcpar0=2.3;
	 FMDcutcpar1=100;
	 break;
       case 4:
	 FMDcutapar0=1.5;
	 FMDcutapar1=300;
	 FMDcutcpar0=2.3;
	 FMDcutcpar1=300;
	 break;
       case 5:
	 FMDcutapar0=1.3;
	 FMDcutapar1=400;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=400;
	 break;
       default: break;
       }

       if(fcollisiontype=="PbPb") {
	 if (nV0A_hits + nV0C_hits < 1.5*(nFMD_fwdV0acc_hits + nFMD_bwdV0acc_hits) - 20) {
	   selectedTracksLeading->Clear();
	   delete selectedTracksLeading;
	   selectedTracksAssociated->Clear();
	   delete selectedTracksAssociated;
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
	 }
       }else{
	 if((nV0A_hits<(FMDcutapar0*nFMD_fwd_hits-FMDcutapar1)) || (nV0C_hits<(FMDcutcpar0*nFMD_bwd_hits-FMDcutcpar1)) ){
	   selectedTracksLeading->Clear();
	   delete selectedTracksLeading;
	   selectedTracksAssociated->Clear();
	   delete selectedTracksAssociated;
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
       }
     }
     
     fFMDV0_post->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
     fFMDV0A_post->Fill(nFMD_fwd_hits, nV0A_hits);
     fFMDV0C_post->Fill(nFMD_bwd_hits, nV0C_hits);
     fFMDV0same_post->Fill(nFMD_bwdV0acc_hits + nFMD_fwdV0acc_hits, nV0C_hits + nV0A_hits);
     fFMDV0Asame_post->Fill(nFMD_fwdV0acc_hits, nV0A_hits);
     fFMDV0Csame_post->Fill(nFMD_bwdV0acc_hits, nV0C_hits);
     }
   }
   fHist_Stat->Fill(6);

   
   if(fcollisiontype=="PbPb"){
     if(!NotSPDClusterVsTrackletBG()) {
       	   selectedTracksLeading->Clear();
	   delete selectedTracksLeading;
	   selectedTracksAssociated->Clear();
	   delete selectedTracksAssociated;
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
     }
   }
   fHist_Stat->Fill(7);
   fHistCentrality->Fill(lCentrality);

   //Centrality Selection  
   if(lCentrality>fCenMax || lCentrality<fCenMin) 
   {
    selectedTracksLeading->Clear();
    delete selectedTracksLeading;
    selectedTracksAssociated->Clear();
    delete selectedTracksAssociated;
    PostData(1, fOutputList);
    PostData(2, fOutputList1);
    PostData(3, fOutputList2);
    return;
   }


   fHistzvertex->Fill(tPrimaryVtxPosition[2]);
   fHistCentzvertex->Fill(lCentrality, tPrimaryVtxPosition[2]);
   
   AliVMultiplicity *tracklets = ((AliAODEvent*)fEvent)->GetTracklets();
   if (!tracklets) return;
   Int_t nTracklets = tracklets->GetNumberOfTracklets();
   fHistCentV0vsTracklets->Fill(lCentrality,nTracklets);
   
   DumpTObjTable("End of FMD vs V0 cuts");
   
if(fAnaMode=="TPCTPC"){
//  if(fasso=="hadron") 
    selectedTracksAssociated=GetAcceptedTracksLeading(fEvent,kFALSE,selectedTracksAssociated);
//  else if (fasso == "Phi")    selectedTracksAssociated = GetAcceptedTracksAssociated(fEvent);
//  else if (fasso == "V0")    selectedTracksAssociated = GetAcceptedV0Tracks(fEvent);
//  else if (fasso == "PID")    selectedTracksAssociated = GetAcceptedTracksPID(fEvent);
//  else if (fasso == "Cascade")    selectedTracksAssociated = GetAcceptedCascadeTracks(fEvent);
 }
// Leading Particle
 if(fAnaMode=="TPCFMD" || fAnaMode=="TPCTPC" || fAnaMode=="TPCFMDC"){
   selectedTracksLeading=GetAcceptedTracksLeading(fEvent,kTRUE,selectedTracksLeading);
 }
 
 DumpTObjTable("End of TPC/ITS track fill");

 TObjArray *selectedTracksAssociated_TPC = new TObjArray; selectedTracksAssociated_TPC->SetOwner(kTRUE);
 selectedTracksAssociated_TPC = GetAcceptedTracksLeading(fEvent,kFALSE,selectedTracksAssociated_TPC);

 TObjArray *selected_TPC_Pairs = new TObjArray; selected_TPC_Pairs->SetOwner(kTRUE);

 FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksAssociated,selectedTracksAssociated_TPC,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0, selected_TPC_Pairs);
 FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedTracksAssociated,selectedTracksAssociated_TPC,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0,selected_TPC_Pairs);
 DumpTObjTable("End of fill  Correlation");

 selected_TPC_Pairs->Clear();
 delete selected_TPC_Pairs;
 selectedTracksLeading->Clear();
 delete selectedTracksLeading;
 selectedTracksAssociated->Clear();
 delete selectedTracksAssociated;
 selectedTracksAssociated_TPC->Clear();
 delete selectedTracksAssociated_TPC;
 
 //DumpTObjTable("after delete TObjects");
 //fNEntries++;
}

TObjArray* AliAnalysisTaskSEpPbCorrelationsJetV2::GetFMDhitsYS(Bool_t Aside){
  TObjArray *tracks1 = new TObjArray;
    tracks1->SetOwner(kTRUE);
    AliAODForwardMult* aodForward =static_cast<AliAODForwardMult*>(fEvent->FindListObject("Forward"));
    // Shape of d2Ndetadphi: 200, -4, 6, q20, 0, 2pi
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
    //AliAnalysisTaskValidation::Tracks ret_vector;
    // FMD has no pt resolution!
    Float_t pt = 0;
    for (Int_t iEta = 1; iEta <= nEta; iEta++) {
      Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
      if (!valid) {
	// No data expected for this eta
	continue;
        }
      Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
      for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
	// Bin content is most likely number of particles!
	Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
	if (mostProbableN > 0) {
	  Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	  //ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, pt, mostProbableN));
	  if(Aside){
	    if(eta<0) continue;
	  } else{
	    if(eta>0) continue;
	  }
	  tracks1->Add(new AliAssociatedVZEROYS(mostProbableN,eta,phi,0,0,0));
	  Double_t cont[3]={eta,phi,lCentrality};
	  fhistfmd->Fill(cont,0,mostProbableN);
	  fh2_FMD_acceptance->Fill(eta,tPrimaryVtxPosition[2]);
	  fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	}
      }
    }    
   return tracks1;
  }



TObjArray *AliAnalysisTaskSEpPbCorrelationsJetV2::GetAcceptedTracksLeading(AliAODEvent *fAOD,Bool_t leading,TObjArray*tracks) {
  //TObjArray *tracks = new TObjArray;
  //tracks->SetOwner(kTRUE);
  Int_t nTracks = fAOD->GetNumberOfTracks();
  Double_t pidqa[4];
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (!aodTrack)      continue;
    if (!IsAcceptedTrack(aodTrack))      continue;
    //    if (aodTrack->Eta()<fEtaMinExtra || aodTrack->Eta()>fEtaMaxExtra) continue;
    if (aodTrack->Charge() == 0)      continue;
    if(leading){
     pidqa[0]=aodTrack->Pt();
     pidqa[1]=aodTrack->Eta();
     pidqa[2]=RangePhi(aodTrack->Phi());
     pidqa[3]=lCentrality;
     //fHistLeadQA->Fill(pidqa,0);
    }
    Int_t SpAsso=0;
    tracks->Add(new AliAssociatedTrackYS(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, SpAsso, 1));
  }
  return tracks;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsJetV2::IsAcceptedTrack(const AliAODTrack *aodTrack) {
  if (!aodTrack)
    return kFALSE;
  //  if(!aodTrack->TestFilterMask(BIT(5))) return kFALSE; // standard cut with
  //  tight DCA cut
  if(fcollisiontype=="PbPb"){
    if (!aodTrack->TestFilterMask(ffilterbit))  return kFALSE; 
    Float_t nCrossedRowsTPC =aodTrack->GetTPCClusterInfo(2,1);
    if (nCrossedRowsTPC < fnoClusters) return kFALSE;
    if(fCutChargedDCAzMax > 0. ||fCutChargedDCAxyMax > 0.){
      Double_t dTrackXYZ[3]  = {0.,0.,0.};
      Double_t dVertexXYZ[3] = {0.,0.,0.};
      Double_t dDCAXYZ[3]    = {0.,0.,0.};
      aodTrack->GetXYZ(dTrackXYZ);
      lPrimaryBestVtx->GetXYZ(dVertexXYZ);
      for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }
      
      if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) >fCutChargedDCAzMax) return kFALSE;
      if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
    }
  }else{
    if (!aodTrack->TestFilterMask(BIT(ffilterbit)))  return kFALSE; 
  }
   if (aodTrack->Pt() < fPtMin)    return kFALSE;
   if(aodTrack->Pt()>fPtMax) return kFALSE;
  //if(!fptdiff) {if (aodTrack->Pt() > fPtMax) return kFALSE;}
  if (aodTrack->Pt() < 0.5) return kFALSE;
  
  if (TMath::Abs(aodTrack->Eta()) > fEtaMax)
    return kFALSE;
  return kTRUE;
}

void AliAnalysisTaskSEpPbCorrelationsJetV2::FillCorrelationTracks( Double_t centrality, TObjArray *triggerArray, TObjArray *selectedTrackArray, TObjArray *selectedTrackArray_TPC, AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut, Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,Float_t bSign, Int_t step, TObjArray *selected_TPC_Pairs)
{
 twoTrackEfficiencyCut=kFALSE;
 twoTrackEfficiencyCutValue=0;
 fTwoTrackCutMinRadius=0;
 bSign=0;
 step=1;//default
 
 if (!triggerHist || !associateHist)    return;
 Double_t binscontTrig[3] = {0.};
 Double_t binscont[7] = {0.}; 


 for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
 {  
  AliAssociatedTrackYS *trigger = (AliAssociatedTrackYS *)triggerArray->At(i);
  if (!trigger)    continue;
  Double_t triggerPt = trigger->Pt();
  Double_t triggerEta = trigger->Eta();
  Double_t triggerPhi = trigger->Phi();
  Int_t trigID = trigger->GetID();
  binscontTrig[0] = triggerPt;
  binscontTrig[1] = fPrimaryZVtx;
  binscontTrig[2] = rand()%12 + 0.5;
 // triggerHist->Fill(binscontTrig, 0);
  for (Int_t j = 0; j < selectedTrackArray_TPC->GetEntriesFast(); j++) {
    AliAssociatedTrackYS *associate_TPC =   (AliAssociatedTrackYS*)selectedTrackArray_TPC->At(j);
    if (!associate_TPC)        continue;
    if (triggerPt < associate_TPC->Pt())          continue;
    if (trigID == associate_TPC->GetID())          continue;
    if((trigger->Charge())*(associate_TPC->Charge())<0) continue;
    triggerHist->Fill(binscontTrig, 0);
    Double_t dTPC_Pairs_Eta = triggerEta-associate_TPC->Eta();
    Double_t dTPC_Pairs_phi = RangePhi(triggerPhi-associate_TPC->Phi());
    selected_TPC_Pairs->Add(new AliAssociatedTPCPairs(trigger->Charge(), triggerEta, triggerPhi, triggerPt, associate_TPC->Pt(), trigger->GetID(), -999, -999, 0, 1, dTPC_Pairs_Eta,dTPC_Pairs_phi));
    }
  }

  for(Int_t j2 = 0;j2 < selected_TPC_Pairs->GetEntriesFast(); j2++)
  {
   AliAssociatedTPCPairs *associate_TPC_Pair = (AliAssociatedTPCPairs*)selected_TPC_Pairs->At(j2);
   if(!associate_TPC_Pair) continue;
   for(Int_t k = 0; k < selectedTrackArray->GetEntriesFast(); k++)
   {
    binscont[0] = associate_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
    binscont[1] = associate_TPC_Pair->Pt_Asso(); // associate TPC pt  
    binscont[2] = associate_TPC_Pair->Pt();; // trigger TPC pt

    //cout<<fh2_pt_trig_asso->GetXaxis()->FindBin(binscont[2])-1<<"ddddddddddd"<<fh2_pt_trig_asso->GetYaxis()->FindBin(binscont[1])-1<<endl;

    Double_t sigma = ((TF2*)fTPCTPClist->FindObject(Form("f2pc_%d_%d",fh2_pt_trig_asso->GetXaxis()->FindBin(binscont[2])-1,fh2_pt_trig_asso->GetYaxis()->FindBin(binscont[1])-1)))->GetParameter(4);
    
    if((associate_TPC_Pair->Getdphi_pairs()<-1*3*sigma)||(associate_TPC_Pair->Getdphi_pairs()>1*3*sigma)) continue;

    //binscont[3] = associate_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2
    binscont[3] = fPrimaryZVtx; //Vz
    AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(k);
    if(!associate) continue;
    Double_t dFMD_Eta = associate->Eta();
    //if(TMath::Abs(dFMD_Eta)<1.7 || TMath::Abs(dFMD_Eta)>4.7) continue;
    binscont[4] = associate_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
    binscont[5] = RangePhi(associate_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
    binscont[6] = rand()%12 + 0.5; ; // Random
    associateHist->Fill(binscont, 0, (Double_t)associate->Multiplicity()); 
   }
  }
}

void AliAnalysisTaskSEpPbCorrelationsJetV2::FillCorrelationTracksMixing(Double_t centrality, Double_t pvxMix, Double_t poolmax, Double_t poolmin,    TObjArray *triggerArray, TObjArray *selectedTrackArray, TObjArray *selectedTrackArray_TPC, AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut,    Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius, Float_t bSign, Int_t step, TObjArray *selected_TPC_Pairs)
{
  Bool_t twoTrackEfficiencyCut_1= twoTrackEfficiencyCut;
  Double_t   twoTrackEfficiencyCutValue_1=  twoTrackEfficiencyCutValue;
  Double_t fTwoTrackCutMinRadius1=fTwoTrackCutMinRadius;
  Float_t bSign1=bSign;
  Int_t step1=step;

  Double_t poolmax1=poolmax;
  Double_t poolmin1=poolmin;

  Double_t counterMix = 0;
  AliEventPool *pool = fPoolMgr->GetEventPool(centrality, pvxMix);
  if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality,
                  pvxMix));
  }

 if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
  Double_t binscontTrig[2];

      mixedDist ->Fill(centrality, pool->NTracksInPool());
      mixedDist2->Fill(centrality, pool->GetCurrentNEvents());      
      Int_t nMix = pool->GetCurrentNEvents();
      Double_t binscont[7];
      for(Int_t j2 = 0; j2 < selected_TPC_Pairs->GetEntriesFast(); j2++)
      {
       AliAssociatedTPCPairs *associate_TPC_Pair = (AliAssociatedTPCPairs*)selected_TPC_Pairs->At(j2);
       if(!associate_TPC_Pair) continue;

       binscontTrig[0] = associate_TPC_Pair->Pt();
       binscontTrig[1] = rand()%12 + 0.5; 
    //   triggerHist->Fill(binscontTrig, 0);

       binscont[0] = associate_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
       binscont[1] = associate_TPC_Pair->Pt_Asso(); // associate TPC pt  
       binscont[2] = associate_TPC_Pair->Pt(); // trigger TPC pt
//       binscont[3] = centrality; // centrality

       Double_t sigma = ((TF2*)fTPCTPClist->FindObject(Form("f2pc_%d_%d",fh2_pt_trig_asso->GetXaxis()->FindBin(binscont[2])-1,fh2_pt_trig_asso->GetYaxis()->FindBin(binscont[1])-1)))->GetParameter(4);

       if((associate_TPC_Pair->Getdphi_pairs()<-1*3*sigma)||(associate_TPC_Pair->Getdphi_pairs()>1*3*sigma)) continue;

       //binscont[3] = associate_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2
       binscont[3] = pvxMix; //Vz
       for (Int_t jMix = 0; jMix < nMix; jMix++) {
        TObjArray *mixEvents = pool->GetEvent(jMix);
        for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
        {
         AliAssociatedTrackYS* associate=(AliAssociatedTrackYS*)  mixEvents->At(k);
         if(!associate)continue;
         Double_t dFMD_Eta = associate->Eta();
         //if(TMath::Abs(dFMD_Eta)<1.7 || TMath::Abs(dFMD_Eta)>4.7) continue;
         binscont[4] = associate_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
         binscont[5] = RangePhi(associate_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
         binscont[6] = rand()%12 + 0.5; ; // Random
         associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix);
        }
       }
     }
  }  
  

  TObjArray* tracksClone=CloneTrack(selectedTrackArray);
  pool->UpdatePool(tracksClone);
}


TObjArray* AliAnalysisTaskSEpPbCorrelationsJetV2::CloneTrack(TObjArray*selectedTrackArray){
  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliAssociatedTrackYS *particle =  (AliAssociatedTrackYS *)selectedTrackArray->At(i);
    tracksClone->Add(new AliAssociatedTrackYS(particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
					      particle->GetID(), particle->GetIDFirstDaughter(),
					      particle->GetIDSecondDaughter(), particle->WhichCandidate(),
					      particle->Multiplicity()));
  }
  
  return tracksClone;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi_FMD(Double_t DPhi) {
  //if (DPhi < (-TMath::Pi() / 2 -0.0001))   DPhi += 2 * TMath::Pi();
  //if (DPhi > (3 * TMath::Pi() / 2-0.0001)) DPhi -= 2*TMath::Pi();
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi2(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < -1.178097)    DPhi += 2 * TMath::Pi();
  return DPhi;
}



void AliAnalysisTaskSEpPbCorrelationsJetV2::DumpTObjTable(const char* note)
{
  if(note) {
    //    printf("TObjectTable::%s",note);
  }
  //  gObjectTable->Print();
}

 Int_t AliAnalysisTaskSEpPbCorrelationsJetV2::ConvertRunNumber(Int_t run){
   if( (265308<run && run< 265526) || (267160<run && run<267167)){
   switch(run){
   case  265309 : return 0;
   case  265332 : return 1;
   case  265334 : return 2;
   case  265336 : return 3;
   case  265338 : return 4;
   case  265339 : return 5;
  case  265342 : return 6;
  case  265343 : return 7;
  case  265344 : return 8;
  case  265377 : return 9;
  case  265378 : return 10;
  case  265381 : return 11;
  case  265383 : return 12;
  case  265384 : return 13;
  case  265385 : return 14;
  case  265387 : return 15;
  case  265388 : return 16;
  case  265419 : return 17;
  case  265420 : return 18;
  case  265421 : return 19;
  case  265422 : return 20;
  case  265424 : return 21;
  case  265425 : return 22;
  case  265426 : return 23;
  case  265427 : return 24;
  case  265435 : return 25;
  case  265499 : return 26;
  case  265500 : return 27;
  case  265501 : return 28;
  case  265521 : return 29;
  case  265525 : return 30;
  case  267166 : return 31; //16t
  case  267165 : return 32; //16t
  case  267164 : return 33; //16t
  case  267163 : return 34; //16t
  case  267161 : return 35; //16t
  default : return 199;
  }
 } else if(run> 280281 && run<281962){
     switch(run){
     case 281961  : return 0;
       //     case 281959 : return 1;
     case 281956 : return 1;
     case 281953:return 2;
     case 281940:return 3;
     case 281939:return 4;
     case 281932:return 5;
     case 281931:return 6;
     case 281928:return 7;
     case 281920:return 8;
     case 281918:return 9;
     case 281915:return 10;//
     case 281895:return 11;//
     case 281894:return 12;//
     case 281892:return 13;//
     case 281633:return 14;//
     case 281583:return 15;
       //     case 281581:return 12;
       //     case 281580:return 13;
     case 281574:return 16;
     case 281569:return 17;
     case 281568:return 18;
     case 281562:return 19;
     case 281557:return 20;
     case 281511:return 21;
     case 281509:return 22;
     case 281477:return 23;
     case 281475:return 24;
     case 281450:return 25;
     case 281449:return 26;
     case 281444:return 27;
     case 281443:return 28;
     case 281441:return 29;
     case 281415:return 30;
     case 281321:return 31;
     case 281301:return 32;
     case 281277:return 33;
     case 281275:return 34;
     case 281273:return 35;
     case 281271:return 36;
     case 281243:return 37;
     case 281242:return 38;
     case 281241:return 39;
     case 281240:return 40;
     case 281213:return 41;
     case 281212:return 42;
     case 281191:return 43;
     case 281190:return 44;
     case 281189:return 45;
     case 281181:return 46;
     case 281180:return 47;//
     case 281179:return 48;
     case 281081:return 49;
     case 281080:return 50;
     case 281062:return 51;
     case 281061:return 52;
     case 281060:return 53;
     case 280999:return 54;
     case 280998:return 55;
     case 280997:return 56;
     case 280994:return 57;
     case 280990:return 58;
     case 280947:return 59;
     case 280940:return 60;
     case 280936:return 61;
     case 280897:return 62;
     case 280890:return 63;//
     case 280881:return 64;//
     case 280880:return 65;
     case 280856:return 66;
     case 280849:return 67;
     case 280848:return 68;
     case 280847:return 69;
     case 280845:return 70;//
     case 280844:return 71;
     case 280842:return 72;
     case 280793:return 73;
     case 280792:return 74;
     case 280787:return 75;
     case 280786:return 76;
     case 280768:return 77;
     case 280767:return 78;
     case 280766:return 79;
     case 280765:return 80;
     case 280764:return 81;
     case 280763:return 82;
     case 280762:return 83;
     case 280761:return 84;
     case 280757:return 85;
     case 280756:return 86;
     case 280755:return 87;
     case 280754:return 88;
     case 280753:return 89;
     case 280729:return 90;
     case 280706:return 91;
     case 280705:return 92;
     case 280681:return 93;
     case 280679:return 94;
       //     case 280676:return 88;
       //     case 280673:return 89;
     case 280671:return 95;
       //     case 280650:return 91;
       //     case 280648:return 92;
     case 280647:return 96;
     case 280645:return 97;
     case 280639:return 98;
     case 280637:return 99;
     case 280636:return 100;
     case 280634:return 101;
     case 280613:return 102;
     case 280583:return 103;
     case 280581:return 104;
     case 280576:return 105;//
     case 280575:return 106;//
     case 280574:return 107;
     case 280551:return 108;
     case 280550:return 109;
     case 280547:return 110;
     case 280546:return 111;
     case 280519:return 112;
     case 280518:return 113;
     case 280499:return 114;
     case 280448:return 115;
     case 280447:return 116;
     case 280446:return 117;
     case 280445:return 118;
     case 280443:return 119;
     case 280419:return 120;
     case 280415:return 121;
     case 280413:return 122;//
     case 280406:return 123;
     case 280405:return 124;
     case 280403:return 125;
     case 280375:return 126;
     case 280374:return 127;
       //     case 280352:return 122;
     case 280351:return 128;
     case 280350:return 129;
     case 280349:return 130;
     case 280348:return 131;
     case 280312:return 132;
     case 280310:return 133;
     case 280290:return 134;
     case 280286:return 135;
     case 280285:return 136;
     case 280284:return 137;
     case 280283:return 138;
     case 280282:return 139;
     default : return 199;
     }
   }     else if (run>274689  && run<286509){
     switch(run){
       //LHC17k
     case 276508:return 0;
     case 276507:return 1;
     case 276506:return 2;
     case 276462:return 3;
     case 276439:return 4;
     case 276438:return 5;
     case 276437:return 6;
     case 276435:return 7;
     case 276351:return 8;
     case 276348:return 9;       
     case 276312:return 10;
     case 276307:return 11;
     case 276302:return 12;
     case 276297:return 13;
     case 276294:return 14;
     case 276292:return 15;
     case 276291:return 16;
     case 276290:return 17;
     case 276259:return 18;
     case 276257:return 19;
     case 276230:return 20;
     case 276205:return 21;
     case 276178:return 22;
     case 276170:return 23;
     case 276104:return 24;
     case 276102:return 25;
     case 276099:return 26;
     case 276098:return 27;
     case 276097:return 28;
     case 276045:return 29;
     case 276041:return 30;
     case 276040:return 31;
     case 276020:return 32;
     case 276019:return 33;
     case 276017:return 34;
     case 276013:return	35;
     case 276012:return 36;
     case 275925:return 37;
     case 275924:return 38;
     case 275847:return 39;
     case 275664:return 40;
     case 275661:return 41;
     case 275657:return 42;
     case 275650:return 43;
     case 275648:return 44;
     case 275647:return 45;
     case 275624:return 46;
     case 275623:return 47;       
     case 275622:return 48;
     case 275621:return 49;
     case 275617:return 50;
     case 275612:return 51;
     case 275559:return 52;
     case 275558:return 53;
     case 275515:return 54;
     case 275472:return 55;
     case 275471:return 56;
     case 275467:return 57;
     case 275459:return 58;
     case 275457:return 59;
     case 275456:return 60;
     case 275453:return 61;
     case 275452:return 62;
     case 275448:return 63;
     case 275443:return 64;
     case 275406:return 65;
     case 275404:return 66;
     case 275401:return 67;
     case 275395:return 68;
     case 275394:return 69;
     case 275372:return 70;
     case 275369:return 71;
     case 275361:return 72;
     case 275360:return 73;
     case 275333:return 74;
     case 275332:return 75;
     case 275328:return 76;
     case 275326:return 77;
     case 275324:return 78;
     case 275322:return 79;
     case 275314:return 80;
     case 275283:return 81;
     case 275247:return 82;
     case 275246:return 83;
     case 275245:return 84;
     case 275239:return 85;
     case 275188:return 86;
     case 275184:return 87;
     case 275180:return 88;
     case 275177:return 89;
     case 275174:return 90;
     case 275173:return 91;
     case 275151:return 92;
     case 275150:return 93;
     case 275149:return 94;
     case 275076:return 95;
     case 275075:return 96;
     case 275073:return 97;
     case 275068:return 98;
     case 275067:return 99;
     case 274979:return 100;
     case 274978:return 101;
     case 274886:return 102;
     case 274882:return 103;
     case 274878:return 104;
     case 274877:return 105;
     case 274822:return 106;
     case 274821:return 107;
     case 274817:return 108;
     case 274815:return 109;
     case 274811:return 110;
     case 274807:return 111;
     case 274806:return 112;
     case 274803:return 113;
     case 274802:return 114;
     case 274801:return 115;
     case 274708:return 116;
     case 274690:return 117;
     default:return 199;
     }
   }     else if (run> 271867 && run<273104){
     switch(run){
       //LHC17h
     case 273103:return 0;
     case 273100:return 1;
     case 273099:return 2;
     case 272949:return 3;
     case 272947:return 4;
     case 272939:return 5;
     case 272935:return 6;
     case 272934:return 7;
     case 272933:return 8;
     case 272932:return 9;
     case 272905:return 10;
     case 272903:return 11;
     case 272871:return 12;
     case 272870:return 13;
     case 272836:return 14;
     case 272833:return 15;
     case 272829:return 16;
     case 272828:return 17;
     case 272784:return 18;
     case 272782:return 19;
     case 272764:return 20;
     case 272763:return 21;
     case 272760:return 22;
     case 272749:return 23;
     case 272747:return 24;
     case 272746:return 25;
     case 272712:return 26;
     case 272692:return 27;
     case 272610:return 28;
     case 272608:return 29;
     case 272607:return 30;
     case 272585:return 31;
     case 272577:return 32;
     case 272575:return 33;
     case 272574:return 34;
     case 272521:return 35;
     case 272468:return 36;
     case 272463:return 37;
     case 272462:return 38;
     case 272461:return 39;
     case 272413:return 40;
     case 272411:return 41;
     case 272400:return 42;
     case 272399:return 43;
     case 272395:return 44;
     case 272394:return 45;
     case 272389:return 46;
     case 272388:return 47;
     case 272360:return 48;
     case 272359:return 49;
     case 272335:return 50;
     case 272194:return 51;
     case 272156:return 52;
     case 272155:return 53;
     case 272154:return 54;
     case 272153:return 55;
     case 272152:return 56;
     case 272123:return 57;
     case 272101:return 58;
     case 272100:return 59;
     case 272041:return 60;
     case 272040:return 61;
     case 272039:return 62;
     case 272038:return 63;
     case 272036:return 64;
     case 271886:return 65;
     case 271881:return 66;
     case 271880:return 67;
     case 271874:return 68;
     case 271873:return 69;
     case 271871:return 70;
     case 271870:return 71;
     case 271868:return 72;
     default:return 199;
     }
   }     else if (run>273590  && run<27443){
     switch(run){
       //LHC17k
     case 274442:return 0;
     case 274390:return 1;
     case 274389:return 2;
     case 274388:return 3;
     case 274387:return 4;
     case 274386:return 5;
     case 274385:return 6;
     case 274364:return 7;
     case 274363:return 8;
     case 274360:return 9;
     case 274352:return 10;
     case 274351:return 11;
     case 274329:return 12;
     case 274283:return 13;
     case 274281:return 14;
     case 274280:return 15;
     case 274278:return 16;
     case 274276:return 17;
     case 274271:return 18;
     case 274270:return 19;
     case 274269:return 20;
     case 274268:return 21;
     case 274266:return 22;
     case 274264:return 23;
     case 274263:return 24;
     case 274259:return 25;
     case 274258:return 26;
     case 274232:return 27;
     case 274212:return 28;
     case 274174:return 29;
     case 274148:return 30;
     case 274147:return 31;
     case 274125:return 32;
     case 274094:return 33;
     case 274092:return 34;
     case 274064:return 35;
     case 274058:return 36;
     case 273986:return 37;
     case 273985:return 38;
     case 273946:return 39;
     case 273943:return 40;
     case 273942:return 41;
     case 273918:return 42;
     case 273889:return 43;
     case 273887:return 44;
     case 273886:return 45;
     case 273885:return 46;
     case 273825:return 47;
     case 273824:return 48;
     case 273719:return 49;
     case 273711:return 50;
     case 273709:return 51;
     case 273695:return 52;
     case 273690:return 53;
     case 273689:return 54;
     case 273687:return 55;
     case 273654:return 56;
     case 273653:return 57;
     case 273593:return 58;
     case 273592:return 59;
     case 273591:return 60;
     default :return 199;
       }
   }     else if (run>274652 && run<274672){
     switch(run){
       //LHC17j
     case 274671:return 0;
     case 274669:return 1;
     case 274667:return 2;
     case 274657:return 3;
     case 274653:return 4;
     default:return 199;
     }
   } else{
     return 199;
   }

}


