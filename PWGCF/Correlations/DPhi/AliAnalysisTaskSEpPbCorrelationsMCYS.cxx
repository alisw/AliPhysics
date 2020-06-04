/**************************************************************************************************
*      Leading Charged Track+V0 Correlations.(Works for Real  Data)                               *
*      Yuko Sekiguchi * Center for Nuclear Study(CNS) , University of Tokyo                       *
*      Email:y_sekiguchi@cns.s.u-tokyo.ac.jp                                                      *
**************************************************************************************************/
#include "AliAnalysisManager.h"
#include "TGrid.h"
#include "AliLog.h"
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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

#include "AliAnalysisTaskSEpPbCorrelationsMCYS.h"

ClassImp(AliAnalysisTaskSEpPbCorrelationsMCYS)
ClassImp(AliAssociatedTrackYSMC)
ClassImp(AliMixTrackYSMC)
ClassImp(AliAssociatedVZEROYSMC)

AliAnalysisTaskSEpPbCorrelationsMCYS::AliAnalysisTaskSEpPbCorrelationsMCYS()
    : AliAnalysisTaskSE(),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      fcentcalib(kFALSE),
      frun2(kTRUE),
      fQA(kFALSE),
      fMCclosure(kFALSE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fptdiff(kFALSE),
      fextractsec(kTRUE),
      fmakehole(kFALSE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      ffillcorrelation(kTRUE),
      fprim(kTRUE),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fPIDResponse(0),
      ffilterbit(5),
      fPtMin(0.2),
      fPtMax(3.0),
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
      fUtils(0x0),
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
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentzvertex(0),
      mixedDist(0),
      mixedDist2(0),
      fHistLeadQA(0),      
      fHistPIDQA(0),
      fhistmcprim(0),
      fhistmcprimfinal(0),     
      fNTrackCorrMC(0),
      fhmcprimvzeta(0),
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
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fV0Amultprim(0),
      fV0Amultmodi(0),
      fh2_V0A(0),
      fh2_V0A_all(0),      
      fh2_V0C(0),
      fh2_V0A_comp(0),
      fh2_V0A_comp_prim(0),
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

  for(Int_t i=0;i<10;i++){
    fhcorreffi[i]=0;
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
}
AliAnalysisTaskSEpPbCorrelationsMCYS::AliAnalysisTaskSEpPbCorrelationsMCYS(const char *name)
    : AliAnalysisTaskSE(name),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      fcentcalib(kFALSE),
      frun2(kTRUE),
      fQA(kFALSE),
      fMCclosure(kFALSE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fptdiff(kFALSE),
      fextractsec(kTRUE),
      fmakehole(kFALSE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      ffillcorrelation(kTRUE),
      fprim(kTRUE),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fPIDResponse(0),
      ffilterbit(5),
      fPtMin(0.2),
      fPtMax(3.0),
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
      fUtils(0x0),
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
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentzvertex(0),
      mixedDist(0),
      mixedDist2(0),
      fHistLeadQA(0),      
      fHistPIDQA(0),
      fhistmcprim(0),
      fhistmcprimfinal(0),
      fNTrackCorrMC(0),
      fhmcprimvzeta(0),
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
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fV0Amultprim(0),
      fV0Amultmodi(0),
      fh2_V0A(0),
      fh2_V0A_all(0),
      fh2_V0C(0),
      fh2_V0A_comp(0),
      fh2_V0A_comp_prim(0),
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
	

	for(Int_t i=0;i<10;i++){
	  fhcorreffi[i]=0;
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
        DefineOutput(1, TList::Class());
        DefineOutput(2, TList::Class());
        DefineOutput(3, TList::Class());
      }

AliAnalysisTaskSEpPbCorrelationsMCYS::~AliAnalysisTaskSEpPbCorrelationsMCYS()
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
  
  if (fPIDResponse) {
    delete fPIDResponse;
    fPIDResponse = 0;
  }
}
void AliAnalysisTaskSEpPbCorrelationsMCYS::UserCreateOutputObjects() {
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fOutputList->SetName("global");
  DefineGeneralOutput();
  PostData(1, fOutputList);

  fOutputList1 = new TList();
  fOutputList1->SetOwner(kTRUE);
  fOutputList1->SetName("anahistos");
  DefineCorrOutput();

  DefineVZEROOutput();
  PostData(2, fOutputList1);

  fOutputList2 = new TList();
  fOutputList2->SetOwner(kTRUE);
  fOutputList2->SetName("QA");
  DefinedQAHistos();
  PostData(3, fOutputList2);

  fEventCuts.AddQAplotsToList(fOutputList);


  frefvz=new TH1F("frefvz","z-vertex",10,-10,10);
  fOutputList2->Add(frefvz);

  if(fcentcalib){
  TGrid::Connect("alien://");
  TFile*file=TFile::Open("alien:///alice/cern.ch/user/y/ysekiguc/fcalibration_centrality_AMPT_modireco.root");
  //  TFile*file=TFile::Open("alien:///alice/cern.ch/user/y/ysekiguc/fcalibration_centrality_AMPT_prim.root");
  //  TFile*file=TFile::Open("/home/yuko/work/local_alicework/MCESDanalysis/draw_result/correction.root");
  
  if(!file) AliError("No correction factor");
  fhcorr[0]=(TH1D*)file->Get("hcent");
  fOutputList->Add(fhcorr[0]);
  }

  
  if(fMCclosure){
  TGrid::Connect("alien://");
  TFile*file1=TFile::Open("alien:///alice/cern.ch/user/y/ysekiguc/corrections/fcorrection_efficiency.root");
  if(!file1) AliError("No correction factor");
  for(Int_t i=0;i<10;i++){
    fhcorreffi[i]=(TH3D*)file1->Get(Form("effi_%d",i));
    // fOutputList2->Add(fhcorr[i]);
  }
  }
  //for(Int_t i=0;i<10;i++){
  //fhcorr[i]=(TH2D*)file->Get(Form("fRefetaphiclone_%d",i));
  // fOutputList2->Add(fhcorr[i]);
  // }

 
   fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
   if (!fPoolMgr)
   return;
   fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);


   fPoolMgr1 = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
   if (!fPoolMgr1)
   return;
   fPoolMgr1->SetTargetValues(fPoolMinNTracks, 0.1, 5);
 }
 void AliAnalysisTaskSEpPbCorrelationsMCYS::DefineGeneralOutput() {

   fHist_Stat = new TH1F("fHist_Stat", "Stat Histogram", 14, -0.5, 13.5);
   fHist_Stat->GetXaxis()->SetBinLabel(1, "All Events");
   fHist_Stat->GetXaxis()->SetBinLabel(2, "Analyzed Events");
   fHist_Stat->GetXaxis()->SetBinLabel(3, "MultSelection OK");
   fHist_Stat->GetXaxis()->SetBinLabel(4, "Vertex OK");
   fHist_Stat->GetXaxis()->SetBinLabel(5, "Centrality OK");
   fHist_Stat->GetXaxis()->SetBinLabel(6, "FMD OK");
   fHist_Stat->GetXaxis()->SetBinLabel(7, "FMD multi cut");
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
   fOutputList->Add(fHist_V0Stat);

   fHistzvertex = new TH1F("fHistzvertex", ";VZ;count", 60, -15, 15);
   fOutputList->Add(fHistzvertex);

   Double_t fmaxcent;
   if (fcollisiontype=="HMPP") fmaxcent=1.;
   else fmaxcent=100.;
   
   fHistCentrality = new TH1F("fHistCentrality", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality);

   fHistCentrality_beforecut = new TH1F("fHistCentrality_beforecut", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality_beforecut);
  
   fHistCentzvertex = new TH2F("fHistCentzvertex", "Cent;VZ;count", 100,0, fmaxcent, 60, -15, 15);
   fOutputList->Add(fHistCentzvertex);

   TTree *settingsTree = new TTree("UEAnalysisSettings", "Analysis Settings in UE estimation");
   settingsTree->Branch("fZVertex", &fZVertex, "fZVertex/D");
   settingsTree->Branch("fEtaMax", &fEtaMax, "fEtaMax/D");
   settingsTree->Branch("fPtMin", &fPtMin, "fPtMin/D");
   settingsTree->Branch("fMaxnSigmaTPCTOF", &fMaxnSigmaTPCTOF, "fMaxnSigmaTPCTOF/D");

   //  settingsTree->Branch("fanamode",&fAnaMode,"fAnaMode/B");
   //  settingsTree->Branch("fanalysisasso",&fanalysisasso,"fanalysisasso/I");
   // settingsTree->Branch("fanalysiscent",&fanalysiscent,"fanalysiscent/I");
   settingsTree->Fill();
   fOutputList->Add(settingsTree);
 }
 void AliAnalysisTaskSEpPbCorrelationsMCYS::DefineVZEROOutput() {

   const Int_t nVZEROBins[3] = {10, 8, 15};
   Double_t binning_eta_vzero[11] = {-3.7, -3.2, -2.7, -2.2, -1.7, 0., 2.8,  3.4,  3.9,  4.5,  5.1};
   Double_t binning_phi_vzero[9] = {0., 0.7853, 1.5707, 2.3561, 3.1415, 3.9269, 4.7123, 5.4977, 6.2831};
   Double_t binning_cent[16] = {0.,  0.1, 1.,  2.,  3.,  5.,  10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1};


   fHist_vzeromult = new TH2F("fHist_vzeromult", "fHist_vzeromult", 64, -0.5, 63.5, 500, 0, 500);
   fOutputList1->Add(fHist_vzeromult);
   fHist_vzeromultEqweighted =  new TH2F("fHist_vzeromultEqweighted", "fHist_vzeromultEqweighted", 64, -0.5, 63.5, 500, 0, 500);
   fOutputList1->Add(fHist_vzeromultEqweighted);
   fHist2dmult = new TH3F("fHist2dmult", "fHist2dmult", 64, -0.5, 63.5, 500, 0, 500, 500, 0, 500);
   fOutputList1->Add(fHist2dmult);
   
   
   fHistVZERO = new AliTHn("fHistVZERO", "fHistVZERO", 1, 3, nVZEROBins);
   fHistVZERO->SetBinLimits(0, binning_eta_vzero);
   fHistVZERO->SetBinLimits(1, binning_phi_vzero);
   fHistVZERO->SetBinLimits(2, binning_cent);
   fOutputList1->Add(fHistVZERO);


 }
 void AliAnalysisTaskSEpPbCorrelationsMCYS::DefinedQAHistos() {

   mixedDist=new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, 101, 200, 0, fPoolMinNTracks*1.5 );
   mixedDist2=new TH2F("mixedDist2", ";centrality;events;events", 101, 0, 101, 100, -0.5, 99.5) ;
   fOutputList2->Add(mixedDist);
   fOutputList2->Add(mixedDist2);
   
   const Int_t ipidBin[5] = {12, 40, 72, 15, 10};
   Int_t ipidBin_effi[]={100,20,60,20,10};
   Double_t binning_pt_lead[13] = {0,0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
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
   Double_t binning_cent[16] = {0., 0.1,  1.,  2.,  3.,   5.,  10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1};
   Double_t binning_zvx[11] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
   if(fasso=="PID" && fQA){
   fHistPIDQA = new AliTHn("fHistPIDQA", "fHistPIDQA", 3, 5, ipidBin);
   fHistPIDQA->SetBinLimits(0, binning_pt_lead);
   fHistPIDQA->SetBinLimits(1, binning_eta);
   fHistPIDQA->SetBinLimits(2, binning_dphi);
   fHistPIDQA->SetBinLimits(3, binning_cent);
   fHistPIDQA->SetBinLimits(4, -10.,10.);
   fHistPIDQA->SetVarTitle(0, "pt");
   fHistPIDQA->SetVarTitle(1, "eta");
   fHistPIDQA->SetVarTitle(2, "phi");
   fHistPIDQA->SetVarTitle(3, "centrality");
   fHistPIDQA->SetVarTitle(4, "vz");
   fOutputList1->Add(fHistPIDQA);
   }

   fHistLeadQA = new AliTHn("fHistLeadQA", "fHistLeadQA", 1, 5, ipidBin_effi);
   fHistLeadQA->SetBinLimits(0,0,5);
   fHistLeadQA->SetBinLimits(1,-1.,1);
   fHistLeadQA->SetBinLimits(2,0.,2*TMath::Pi());
   fHistLeadQA->SetBinLimits(3,0.,100);
   fHistLeadQA->SetBinLimits(4,-10.,10.);
   fHistLeadQA->SetVarTitle(0, "pt");
   fHistLeadQA->SetVarTitle(1, "eta");
   fHistLeadQA->SetVarTitle(2, "phi");
   fHistLeadQA->SetVarTitle(3, "centrality");
   fHistLeadQA->SetVarTitle(4, "vz");
   if(fQA) fOutputList1->Add(fHistLeadQA);



   if(!fDataType){
     fhistmcprim=new AliTHn("fhistmcprim","fhistmcprim",1,5,ipidBin_effi);
     //     fhistmcprim->SetBinLimits(0,binning_pt_lead);
     fhistmcprim->SetBinLimits(0,0,5);
     fhistmcprim->SetBinLimits(1,-1.,1);
     fhistmcprim->SetBinLimits(2,0.,2*TMath::Pi());
     fhistmcprim->SetBinLimits(3,0.,100);
     fhistmcprim->SetBinLimits(4,-10.,10.);
     fhistmcprim->SetVarTitle(0,"pt");
     fhistmcprim->SetVarTitle(1,"eta");
     fhistmcprim->SetVarTitle(2,"phi");
     fhistmcprim->SetVarTitle(3,"centrality");
     fhistmcprim->SetVarTitle(4,"vz");
     if(fQA)     fOutputList2->Add(fhistmcprim);


     const Int_t ipidBinfmd[5] = {1, 200, 20, 15,10};
     fhistmcprimfinal=new AliTHn("fhistmcprimfinal","fhistmcprimfinal",1,5,ipidBinfmd);
     //     fhistmcprimfinal->SetBinLimits(0,binning_pt_lead);
     fhistmcprimfinal->SetBinLimits(0,0,3000);
     fhistmcprimfinal->SetBinLimits(1,-4,6);
     fhistmcprimfinal->SetBinLimits(2,0.,2*TMath::Pi());
     fhistmcprimfinal->SetBinLimits(3,binning_cent);
     fhistmcprimfinal->SetBinLimits(4,-10.,10.);
     fhistmcprimfinal->SetVarTitle(0,"pt");
     fhistmcprimfinal->SetVarTitle(1,"eta");
     fhistmcprimfinal->SetVarTitle(2,"phi");
     fhistmcprimfinal->SetVarTitle(3,"centrality");
     fhistmcprimfinal->SetVarTitle(4,"vz");
     fOutputList2->Add(fhistmcprimfinal);

     if(fcollisiontype=="pp"||fcollisiontype=="HMPP"||fcollisiontype=="MBPP") fNTrackCorrMC=new TH2D("fNTrackCorrMC","fNTrackCorrMC",150,0,150,150,0,150);
     else if(fcollisiontype=="pPb") fNTrackCorrMC=new TH2D("fNTrackCorrMC","fNTrackCorrMC",200,0,200,200,0,200);
     else if(fcollisiontype=="PbPb") fNTrackCorrMC=new TH2D("fNTrackCorrMC","fNTrackCorrMC",1000,0,3500,1000,0,3500);
     fOutputList2->Add(fNTrackCorrMC);
     
     fhmcprimvzeta=new TH2D("fhmcprimvzeta","fhmcprimvzeta",200,-4,6,20,-10,10);
     fOutputList2->Add(fhmcprimvzeta);
     fhmcprimpdgcode=new TH1D("fhmcprimpdgcode","fhmcprimpdgcode",4000,-0.5,3999.5);
     fOutputList2->Add(fhmcprimpdgcode);
     fh2_FMD_acceptance_prim=new TH2D("fh2_FMD_acceptance_prim","fh2_FMD_acceptance_prim",200,-4,6,200,-10,10);
     fOutputList2->Add(fh2_FMD_acceptance_prim);
     fh2_FMD_eta_phi_prim=new TH2D("fh2_FMD_eta_phi_prim","fh2_FMD_eta_phi_prim",200,-4,6,20,0,2*TMath::Pi());
     fOutputList2->Add(fh2_FMD_eta_phi_prim);

     for(Int_t i=0;i<4;i++){
       fhrefetaFMD[i]=new TH1D(Form("fhrefetaFMD_%d",i),Form("fhrefetaFMD_%d",i),200,-4,6);
       fhrefphiFMD[i]=new TH1D(Form("fhrefphiFMD_%d",i),Form("fhrefphiFMD_%d",i),100,0,2*TMath::Pi());
       fOutputList2->Add(fhrefetaFMD[i]);
       fOutputList2->Add(fhrefphiFMD[i]);
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
   

   //   if(fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC" || fAnaMode=="ITSFMD" || fAnaMode=="ITSFMDC" || fAnaMode=="FMDFMD" || fAnaMode=="SECA"|| fAnaMode=="SECC"){
   fFMDV0 = new TH2F("FMDV0", "FMD vs V0 pre cut;FMD;V0;",2000, 0, 2000, 2000, 0, 2000);
   fOutputList2->Add(fFMDV0);
   fFMDV0_post=new TH2F("FMDV0_post", "FMD vs V0 post cut;FMD;V0;",2000, 0, 2000, 2000, 0, 2000);
   fOutputList2->Add(fFMDV0_post);
   fFMDV0A = new TH2F("FMDV0A", "FMD vs V0A;FMD;V0A;",1000, 0, 1000, 1000, 0, 1000);
   fOutputList2->Add(fFMDV0A);
   fFMDV0A_post = new TH2F("FMDV0A_post", "FMD vs V0A post cut;FMD;V0A;",1000, 0, 1000, 1000, 0, 1000);
   fOutputList2->Add(fFMDV0A_post);
   fFMDV0C = new TH2F("FMDV0C", "FMD vs V0C;FMD;V0C;",1000, 0, 1000, 1000, 0, 1000);
   fOutputList2->Add(fFMDV0C);
   fFMDV0C_post = new TH2F("FMDV0C_post", "FMD vs V0C post cut;FMD;V0C;",1000, 0, 1000, 1000, 0, 1000);
   fOutputList2->Add(fFMDV0C_post);
   
   fV0Amultprim = new TH1F("fV0Amultprim", " V0mult",1000,-0.5,999.5);
   fOutputList2->Add(fV0Amultprim);
   fV0Amultmodi = new TH1F("fV0Amultmodi", " V0multmodi",1000,-0.5,999.5);
   fOutputList2->Add(fV0Amultmodi);
   
   fh2_V0A = new TH2F("fh2_V0A", " V0 vs primary tracks in V0",250, 0, 1000, 250, 0, 1000);
   fOutputList2->Add(fh2_V0A);
   fh2_V0A_all = new TH2F("fh2_V0A_all", " V0 vs primary tracks in V0",250, 0, 1000, 250, 0, 1000);
   fOutputList2->Add(fh2_V0A_all);
   fh2_V0C=new TH2F("fh2_V0C", " V0  vs primary tracks in V0",250, 0, 1000, 250, 0, 1000);
   fOutputList2->Add(fh2_V0C);
   
   fh2_V0A_comp = new TH2F("fh2_V0A_comp", "V0 modi reco vs V0 raw",250, 0, 1000, 250, 0, 1000);
   fOutputList2->Add(fh2_V0A_comp);
   fh2_V0A_comp_prim = new TH2F("fh2_V0A_comp_prim", "V0 modi reco vs all prim",250, 0, 1000, 250, 0, 1000);
   fOutputList2->Add(fh2_V0A_comp_prim);
   
   
   fh2_ITS_acceptance=new TH2D("fh2_ITS_acceptance","fh2_ITS_acceptance",200,-10,10,200,-4,6);
   fOutputList2->Add(fh2_ITS_acceptance);
   
   fh2_SPD_multcorr=new TH2F("fh2_SPD_multcorr","fh2_SPD_multcorr",400,0,400,2000,0,2000);
   fOutputList2->Add(fh2_SPD_multcorr);
   
   fh2_SPDV0_multcorr=new TH2F("fh2_SPDV0_multcorr","fh2_SPDV0_multcorr",400,0,400,2000,0,2000);
   fOutputList2->Add(fh2_SPDV0_multcorr);
   
   fh2_SPDtrack_multcorr=new TH2F("fh2_SPDtrack_multcorr","fh2_SPDtrack_multcorr",400,0,400,400,0,400);
   fOutputList2->Add(fh2_SPDtrack_multcorr);
   
   fhtrackletsdphi=new TH1F("fhtrackletsdphi","dphi tracklets",100,-100,100);
   fOutputList2->Add(fhtrackletsdphi);
   
   fh2_FMD_acceptance=new TH2D("fh2_FMD_acceptance","fh2_FMD_acceptance",200,-4,6,200,-10,10);
   fOutputList2->Add(fh2_FMD_acceptance);
   
   
   fh2_FMD_eta_phi=new TH2D("fh2_FMD_eta_phi","fh2_FMD_eta_phi",200,-4,6,20,0,2*TMath::Pi());
   fOutputList2->Add(fh2_FMD_eta_phi);
   
   fhistfmdphiacc=new TH2D("fhistfmdphiacc","fhistfmdphiacc",200,-4,6,15,binning_cent);
   fOutputList2->Add(fhistfmdphiacc);
   
   const Int_t ifmdbin[4]={200,20,15,20};
   fhistfmd=new AliTHn("fhistfmd","fhistfmd",1,4,ifmdbin);
   fhistfmd->SetBinLimits(0,-4.,6.);
   fhistfmd->SetBinLimits(1,0.,2*TMath::Pi());
   fhistfmd->SetBinLimits(2,binning_cent);
   //   fhistfmd->SetBinLimits(3,-10.,10.);
   fhistfmd->SetBinLimits(3,-10,10);
   fhistfmd->SetVarTitle(0,"eta");
   fhistfmd->SetVarTitle(1,"phi");
   fhistfmd->SetVarTitle(2,"centrality");
   fhistfmd->SetVarTitle(3,"vzy");
   fOutputList2->Add(fhistfmd);
   
   const Int_t iitsbin[4]={200,20,20,40};
   const	 Double_t MinITS[4]={-4,0.,0.,-10.};
   const	 Double_t MaxITS[4]={6,2*TMath::Pi(),100.,10.};
   fhistits=new THnSparseF("fhistits","fhistits",4,iitsbin,MinITS,MaxITS);
   fOutputList2->Add(fhistits);
     
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
   fOutputList2->Add(fhSecFMD);
   
   if(fasso=="PID"){
     for(Int_t i=0;i<6;i++){
       fHistNsig[i]=new TH2D(Form("fHistNsig_%d",i),Form("HistNsig_%d",i), 160, 0., 8., 600, -30., 30);
       fOutputList2->Add(fHistNsig[i]);
       fHistNsigcorr[i]=new TH2D(Form("fHistNsigcorr_%d",i),"fHistNsigcorr",500,-10,10,500,-10,10);
       fOutputList2->Add(fHistNsigcorr[i]);
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

void AliAnalysisTaskSEpPbCorrelationsMCYS::DefineCorrOutput() {
  
  Double_t binning_pt_assoc[13] = {0., 0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
				   2.0, 2.5, 3.0,  3.5, 5.0,  8.};
  Double_t binning_pt_lead[13] = {0., 0.2, 0.5, 0.75, 1.0, 1.25, 1.5,
				  2.0, 2.5, 3.0,  3.5, 5.0, 8.0};
  Double_t binning_cent[12] = {0., 5.,  10., 20.,
			       30., 40., 50., 60., 70., 80., 90., 100.1};
  //Double_t binning_cent_HMPP[12] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,1.1};                

  
  Double_t binning_deta[49] = {-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5,
			       -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
			       -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5,
			       0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,
			       1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4};
  
  Double_t binning_dphi[73] = {-1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464,
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
  const Int_t nEvtVars = 2;
  const Int_t iEvtBin[2] = {11, 11};

  //centrality bin
  Double_t binning_cent_HMPP[12] = {0., 0.01,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,1.};
  const Double_t binning_cent_fmdfmd_PbPb[10] = {0., 5., 10.,  20., 30., 40., 50.,60.,70.,80};
  const  Double_t binning_cent_MBPP[8]={0.,0.1,1.,10.,20.,40.,60.,100.1};

  Double_t binning_cent_trig[8] = {0., 5.,  10., 20., 40., 60.,70.,100.1};
  
  Int_t ncentbin;
  if(fCentType=="Manual") {
    ncentbin=9;
  }else {
    if(fcollisiontype=="HMPP") ncentbin=11;
    else if(fcollisiontype=="PbPb") ncentbin=9;
    else if (fcollisiontype=="MBPP") ncentbin=7;
    else ncentbin=7;
  }
     
  Int_t nCFStepstrig=1;
  if(fasso=="hadron")  nCFStepstrig=1;
   else  if (fasso == "V0" || fasso == "Phi")    nCFStepstrig = 7;
   else  if (fasso == "Cascade")    nCFStepstrig = 6;
   else  if(fasso=="PID")    nCFStepstrig=3;

   Double_t binning_eta_vzero[11]={-3.7,-3.2,-2.7,-2.2,-1.7,0.,2.8,3.4,3.9,4.5,5.1};
   Double_t binning_phi_vzero[9]={0.,0.7853,1.5707,2.3561,3.1415,3.9269,4.7123,5.4977,6.2831};

   //Bins for FMD
   const Double_t binning_etafmd[33]={
     1.7,1.8,1.9,
     2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
     3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
     4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
   const Double_t binning_etafmdc[18]={
     -3.4,-3.3,-3.2,-3.1,-3.0,
     -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,
     -1.9,-1.8,-1.7};
   Double_t binning_zvx[11] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
   Double_t binning_mult_trig[10]={0,20,40,60,80,100,120,140,160,200};
   if(fAnaMode=="FMDFMD"){
     const Int_t nEvtVarsV0Leading=3;
     const Int_t iEvtBinV0Leading[3]={11,32,ncentbin};
     fHistTriggerTrack= new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsV0Leading, iEvtBinV0Leading);
     if(fcollisiontype=="HMPP")  fHistTriggerTrack->SetBinLimits(0,binning_cent_HMPP);
     else if(fcollisiontype=="PbPb")  fHistTriggerTrack->SetBinLimits(1,binning_cent_fmdfmd_PbPb);
     else if(fcollisiontype=="MBPP")  fHistTriggerTrack->SetBinLimits(1,binning_cent_MBPP);
     else fHistTriggerTrack->SetBinLimits(0,binning_cent_trig);
     fHistTriggerTrack->SetBinLimits(1,binning_etafmd);
     fHistTriggerTrack->SetBinLimits(2,-10.,10.);
     fHistTriggerTrack->SetVarTitle(0,"centrality");
     fHistTriggerTrack->SetVarTitle(1,"eta");
     fHistTriggerTrack->SetVarTitle(2,"z vertex");
     const Int_t nEvtVarsV0Leadingmix=2;
     const Int_t iEvtBinV0Leadingmix[2]={11,32};
     fHistTriggerTrackMix= new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVarsV0Leadingmix, iEvtBinV0Leadingmix);
     fHistTriggerTrackMix->SetBinLimits(0,binning_cent);
     if(fAnaMode=="FMDFMD")    fHistTriggerTrackMix->SetBinLimits(1,binning_etafmd);
     fHistTriggerTrackMix->SetVarTitle(0,"centrality");
     fHistTriggerTrackMix->SetVarTitle(1,"eta");
   }else if(fAnaMode=="TPCFMD" ||fAnaMode=="TPCFMDC"||fAnaMode=="ITSFMD" || fAnaMode=="ITSFMDC"){
     //     Double_t binning_cent_trig[9] = {0., 5.,  10., 20.				      , 40., 60.,70.,80.,100.1};
     const Int_t nEvtVarsFMD = 4;
     Int_t netabin;
     if(fAnaMode=="TPCFMD"||fAnaMode=="TPCFMDC")  netabin=4;
     else    netabin=18;
     const Int_t iEvtBinFMD[4] = {12,ncentbin,10,netabin};
     Double_t binning_eta_tpcfmd[5]={-0.8,-0.4,-0.,0.4,0.8};
     Double_t binning_eta_itsfmd[19]={-1.7, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.,0.2,  0.4,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6,  1.7};
     fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsFMD, iEvtBinFMD);
     fHistTriggerTrack->SetBinLimits(0, binning_pt_lead);
     if(fCentType=="Manual"){
       fHistTriggerTrack->SetBinLimits(1,binning_mult_trig);
     }else{
       if(fcollisiontype=="HMPP")  fHistTriggerTrack->SetBinLimits(1,binning_cent_HMPP);
       else if(fcollisiontype=="PbPb")  fHistTriggerTrack->SetBinLimits(1,binning_cent_fmdfmd_PbPb);
       else if(fcollisiontype=="MBPP")  fHistTriggerTrack->SetBinLimits(1,binning_cent_MBPP);
       else fHistTriggerTrack->SetBinLimits(1,binning_cent_trig);
     }
     fHistTriggerTrack->SetBinLimits(2, -10.,10.);
     if(fAnaMode=="TPCFMD"||fAnaMode=="TPCFMDC")	   fHistTriggerTrack->SetBinLimits(3, binning_eta_tpcfmd);
     else fHistTriggerTrack->SetBinLimits(3, binning_eta_itsfmd); 
     
     fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrack->SetVarTitle(1, "centrality");
     fHistTriggerTrack->SetVarTitle(2, "zvertex");
     fHistTriggerTrack->SetVarTitle(3, "TPC/Eta eta");
     
     fHistTriggerTrackMix = new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVars, iEvtBin);
     fHistTriggerTrackMix->SetBinLimits(0, binning_pt_lead);
     fHistTriggerTrackMix->SetBinLimits(1, binning_cent_trig);
     fHistTriggerTrackMix->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrackMix->SetVarTitle(1, "centrality");
     
   }else{
     const Int_t nEvtVars_tpctpc = 3;
     const Int_t iEvtBin_tpctpc[3] = {11, 11,20};
     
     fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVars_tpctpc, iEvtBin_tpctpc);
     fHistTriggerTrack->SetBinLimits(0, binning_pt_lead);
     fHistTriggerTrack->SetBinLimits(1, binning_cent);
     fHistTriggerTrack->SetBinLimits(2, -10.,10.);
     
     fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrack->SetVarTitle(1, "centrality");
     fHistTriggerTrack->SetVarTitle(2, "vz(cm)");
     
     
     fHistTriggerTrackMix = new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVars, iEvtBin);
     fHistTriggerTrackMix->SetBinLimits(0, binning_pt_lead);
     fHistTriggerTrackMix->SetBinLimits(1, binning_cent);
     fHistTriggerTrackMix->SetVarTitle(0, "leading p_{T} GeV/c");
     fHistTriggerTrackMix->SetVarTitle(1, "centrality");
   }
   fOutputList1->Add(fHistTriggerTrack);
   fOutputList1->Add(fHistTriggerTrackMix);
   
   const Int_t nTrackVars = 5;
   const Int_t iTrackBin[5] = {48, 11, 11, 15, 72};
   //////////////////////////////////////////
   //Containers two particle correlation
   //////////////////////////////////////////

   Int_t nCFSteps = 1;
   if(fasso=="hadron")    nCFSteps=1;
   else if (fasso == "V0" || fasso == "Phi")    nCFSteps = 7;
   else  if (fasso == "Cascade")    nCFSteps = 6;
   else    if(fasso=="PID")    nCFSteps=3;
   Double_t binning_dphi_vzero[9]={-1.178097,-0.392699,0.392699,1.178097,1.963495,2.748893,3.534291,4.319689,5.105088};
   if(fAnaMode=="TPCTPC") {
     const Double_t binning_cent_tpctpc[7]={0.,5.,10.,20.,40.,60.,100.1};
     const Int_t iTrackBin_TPCTPC[6] = {48, 11, 11, 6, 72, 10};
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 6, iTrackBin_TPCTPC);
     fHistReconstTrack->SetBinLimits(0, binning_deta);
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
     fHistReconstTrackMix->SetBinLimits(0, binning_deta);
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
   }else if(fAnaMode=="TPCV0A"||fAnaMode=="TPCV0C"){
     const Int_t iTrackBin_VZEROA[5]={66,11,10,11,72};
     const Int_t iTrackBin_VZEROC[5]={62,11,10,11,72};
     Double_t binning_detaVZEROATPC[67]={-5.6,-5.55,-5.5,-5.45,-5.4,-5.35,-5.3,-5.25,-5.2,-5.15,-5.1,-5.05,-5.0,-4.95, -4.9,-4.85, -4.8, -4.75, -4.7, -4.65,-4.6,-4.55,-4.5,-4.45,-4.4,-4.35,-4.3, -4.25,-4.2,-4.15,-4.1,-4.05,-4.0,-3.95,-3.9,-3.85,-3.8,-3.75,-3.7,-3.65,-3.6,-3.55,-3.5,-3.45,-3.4,-3.35,-3.3,-3.25,-3.2,-3.15,-3.1,-3.05,-3.0,-2.95,-2.9,-2.85,-2.8,-2.75,-2.7,-2.65,-2.6,-2.55,-2.5,-2.45,-2.4,-2.35,-2.3};
     Double_t binning_detaVZEROCTPC[63]={1.15, 1.2, 1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2,2.25,2.3, 2.35, 2.4, 2.45, 2.5, 2.55,2.6, 2.65, 2.7, 2.75, 2.8, 2.85,2.9, 2.95, 3.0, 3.05, 3.1,3.15, 3.2, 3.25, 3.3,3.35, 3.4, 3.45,3.5,3.55, 3.6,3.65, 3.7, 3.75,3.8, 3.85, 3.9,3.95, 4.0,4.05, 4.1,4.15, 4.2, 4.25};
     if(fAnaMode=="TPCV0A"){
       fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars, iTrackBin_VZEROA);
       fHistReconstTrack->SetBinLimits(0,binning_detaVZEROATPC);
     }else{
       fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars, iTrackBin_VZEROC);
       fHistReconstTrack->SetBinLimits(0,binning_detaVZEROCTPC);
     }
     fHistReconstTrack->SetBinLimits(1,binning_pt_lead);
     fHistReconstTrack->SetBinLimits(2,binning_eta_vzero);
     fHistReconstTrack->SetBinLimits(3,binning_cent);
     fHistReconstTrack->SetBinLimits(4,binning_dphi);
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrack->SetVarTitle(1,"p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(2,"Vzero Eta");
     fHistReconstTrack->SetVarTitle(3,"centrality");
     fHistReconstTrack->SetVarTitle(4,"#Delta#phi");

     if(fAnaMode=="TPCV0A"){
       fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars,iTrackBin_VZEROA);
       fHistReconstTrackMix->SetBinLimits(0,binning_detaVZEROATPC);
     }else{
       fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars,iTrackBin_VZEROC);
       fHistReconstTrackMix->SetBinLimits(0,binning_detaVZEROCTPC);
     }
     fHistReconstTrackMix->SetBinLimits(1,binning_pt_lead);
     fHistReconstTrackMix->SetBinLimits(2,binning_eta_vzero);
     fHistReconstTrackMix->SetBinLimits(3,binning_cent);
     fHistReconstTrackMix->SetBinLimits(4,binning_dphi);
     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1,"p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2,"Vzero Eta");
     fHistReconstTrackMix->SetVarTitle(3,"centrality");
     fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");

   }else if (fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
     
     Double_t binning_detaFMDTPC[49]={
				      -5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.0,
				      -4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,
				      -3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,
				      -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,
				      -1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,
				      -0.9};
     Double_t binning_detaFMDCTPC[34]={
				       0.9,
				       1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
				       2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
				       3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
				       4.0,4.1,4.2};

     
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
	
     Int_t ndetatpcfmd;
     Int_t nfmdbin;
     if(fAnaMode=="TPCFMD") {
       ndetatpcfmd=48;
       nfmdbin=32;
     }else{// if(fAnaMode=="TPCFMDC") {
       ndetatpcfmd=33;
       nfmdbin=17;
     }
     Double_t binning_pt_fmdtpc[5] = {0.2, 0.5, 1.0, 3.0, 8.0};
     Int_t ntpcpt;
     if(fptdiff) ntpcpt=4;
     else ntpcpt=1;
 	    
     Int_t nbineta=72;
     
     //     const Double_t binning_cent_fmdfmd[9]={0.,5.,10.,20.,40.,60.,70,80.,100.1};
     //     const Double_t binning_cent_fmdfmd_HMPP[9]={0.,0.01,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
     //     const Double_t binning_cent_fmdfmd[8]={0.,5.,10.,20.,40.,60.,70.,80.,100.1};
	 //     const Int_t iTrackBin_tpcfmd[7]={ndetatpcfmd,1,nfmdbin,8,72,10,4};
     const Int_t iTrackBin_tpcfmd[7]={ndetatpcfmd,ntpcpt,nfmdbin,ncentbin,nbineta,10,4};
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 7, iTrackBin_tpcfmd);
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 7,iTrackBin_tpcfmd);
     if(fAnaMode=="TPCFMD") {
       fHistReconstTrack->SetBinLimits(0,binning_detaFMDTPC);
       fHistReconstTrack->SetBinLimits(2,binning_etafmd);
       //Mixed Events
       fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDTPC);
       fHistReconstTrackMix->SetBinLimits(2,binning_etafmd);
     } else if(fAnaMode=="TPCFMDC") { 
       fHistReconstTrack->SetBinLimits(0,binning_detaFMDCTPC);
       fHistReconstTrack->SetBinLimits(2,binning_etafmdc);
       //Mixed Events
       fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDCTPC);
       fHistReconstTrackMix->SetBinLimits(2,binning_etafmdc);
     }
     
     if(fptdiff)fHistReconstTrack->SetBinLimits(1,binning_pt_fmdtpc);
     else fHistReconstTrack->SetBinLimits(1,fPtMin,fPtMax);
					  
     if(fcollisiontype=="HMPP") fHistReconstTrack->SetBinLimits(3,binning_cent_HMPP);
     else if(fcollisiontype=="PbPb") fHistReconstTrack->SetBinLimits(3,binning_cent_fmdfmd_PbPb);
     else if(fcollisiontype=="MBPP") fHistReconstTrack->SetBinLimits(3,binning_cent_MBPP);
     else fHistReconstTrack->SetBinLimits(3,binning_cent_trig);
     //     fHistReconstTrack->SetBinLimits(4,binning_dphi_reduce);
     fHistReconstTrack->SetBinLimits(4,binning_dphi);
     fHistReconstTrack->SetBinLimits(5,-10.,10.);
     fHistReconstTrack->SetBinLimits(6,-0.8,0.8);
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrack->SetVarTitle(1,"p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(2,"FMD Eta");
     fHistReconstTrack->SetVarTitle(3,"centrality");
     fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
     fHistReconstTrack->SetVarTitle(5,"z vertex");
     fHistReconstTrack->SetVarTitle(6,"TPC eta");
     
     if(fptdiff) fHistReconstTrackMix->SetBinLimits(1,binning_pt_fmdtpc);
     else  fHistReconstTrackMix->SetBinLimits(1,fPtMin,fPtMax);
     
     if(fcollisiontype=="HMPP") fHistReconstTrackMix->SetBinLimits(3,binning_cent_HMPP);
     else if(fcollisiontype=="PbPb") fHistReconstTrackMix->SetBinLimits(3,binning_cent_fmdfmd_PbPb);
     else if(fcollisiontype=="MBPP") fHistReconstTrackMix->SetBinLimits(3,binning_cent_MBPP);
     else fHistReconstTrackMix->SetBinLimits(3,binning_cent_trig);
     //     fHistReconstTrackMix->SetBinLimits(4,binning_dphi_reduce);
     fHistReconstTrackMix->SetBinLimits(4,binning_dphi);
     fHistReconstTrackMix->SetBinLimits(5,-10.,10.);
     fHistReconstTrackMix->SetBinLimits(6,-0.8,0.8);
     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1,"p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2,"FMD Eta");
     fHistReconstTrackMix->SetVarTitle(3,"centrality");
     fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
     fHistReconstTrackMix->SetVarTitle(5,"z vertex");
     fHistReconstTrackMix->SetVarTitle(6,"TPC eta");	
     
   }else if(fAnaMode=="FMDFMD"){
     const Int_t nTrackVars_fmdfmd = 6;
     //     const Double_t binning_cent_fmdfmd[9]={0.,5.,10.,20.,40.,60.,70.,80.,100.1};
     //	 const Double_t binning_cent_fmdfmd_HMPP[9]={0.,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
     //const Double_t binning_cent_fmdfmd_HMPP[9]={0.,0.01,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

     const Int_t iTrackBin_fmdfmd[6]={49,17,32,ncentbin,20,10};

     fHistReconstTrack= new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
     fHistReconstTrack->SetBinLimits(0,3.425,8.325);
     //    fHistReconstTrack->SetBinLimits(0,3.525,8.325);
     fHistReconstTrack->SetBinLimits(1,binning_etafmdc);
     fHistReconstTrack->SetBinLimits(2,binning_etafmd);
     //     if(fcollisiontype=="HMPP")  fHistReconstTrack->SetBinLimits(3,binning_cent_fmdfmd_HMPP);
     //	 else fHistReconstTrack->SetBinLimits(3,binning_cent_fmdfmd);
     if(fcollisiontype=="HMPP") fHistReconstTrack->SetBinLimits(3,binning_cent_HMPP);
     else if(fcollisiontype=="PbPb") fHistReconstTrack->SetBinLimits(3,binning_cent_fmdfmd_PbPb);
     else if(fcollisiontype=="MBPP") fHistReconstTrack->SetBinLimits(3,binning_cent_MBPP);
     else fHistReconstTrack->SetBinLimits(3,binning_cent_trig);

     if(!fprim)   fHistReconstTrack->SetBinLimits(4,-0.55*TMath::Pi(),1.45*TMath::Pi());
     else  fHistReconstTrack->SetBinLimits(4,-0.5*TMath::Pi()-0.0001,1.5*TMath::Pi()-0.0001);
     fHistReconstTrack->SetBinLimits(5,-10.,10.);
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrack->SetVarTitle(1,"FMD(Asso) Eta");
     fHistReconstTrack->SetVarTitle(2,"FMD(Trigger) Eta");
     fHistReconstTrack->SetVarTitle(3,"centrality");
     fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
     fHistReconstTrack->SetVarTitle(5,"z vertex");
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
     fHistReconstTrackMix->SetBinLimits(0,3.425,8.325);
     //fHistReconstTrackMix->SetBinLimits(0,3.525,8.325);
     fHistReconstTrackMix->SetBinLimits(1,binning_etafmdc);
     fHistReconstTrackMix->SetBinLimits(2,binning_etafmd);
     
     if(fcollisiontype=="HMPP") fHistReconstTrackMix->SetBinLimits(3,binning_cent_HMPP);
     else if(fcollisiontype=="PbPb") fHistReconstTrackMix->SetBinLimits(3,binning_cent_fmdfmd_PbPb);
     else if(fcollisiontype=="MBPP") fHistReconstTrackMix->SetBinLimits(3,binning_cent_MBPP);
     else fHistReconstTrackMix->SetBinLimits(3,binning_cent_trig);
     
     //     if(fcollisiontype=="HMPP")fHistReconstTrackMix->SetBinLimits(3,binning_cent_fmdfmd_HMPP);
     //     else fHistReconstTrackMix->SetBinLimits(3,binning_cent_fmdfmd);



     //    fHistReconstTrackMix->SetBinLimits(4,-0.551*TMath::Pi(),1.449*TMath::Pi());   
     if(!fprim)fHistReconstTrackMix->SetBinLimits(4,-0.55*TMath::Pi(),1.45*TMath::Pi());
     else fHistReconstTrackMix->SetBinLimits(4,-0.5*TMath::Pi()-0.0001,1.5*TMath::Pi()-0.0001);
     fHistReconstTrackMix->SetBinLimits(5,-10,10.);
     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1,"FMD(Asso) Eta");
     fHistReconstTrackMix->SetVarTitle(2,"FMD(Trigger) Eta");
     fHistReconstTrackMix->SetVarTitle(3,"centrality");
     fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
     fHistReconstTrackMix->SetVarTitle(5,"z vertex");
   }else if(fAnaMode=="SECA" || fAnaMode=="SECC"){
     const Int_t nTrackVars_fmdfmd = 5;
     const Double_t binning_cent_fmdfmd[12]={0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.1};
     Int_t binsec;
     if(fAnaMode=="SECA") binsec=64;
     else  binsec=34;

     const Int_t iTrackBin_fmdfmd[5]={130,binsec,11,20,10};
     fHistReconstTrack= new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
     fHistReconstTrack->SetBinLimits(0,-3.275,3.225);
     if(fAnaMode=="SECA") fHistReconstTrack->SetBinLimits(1,1.7,4.9);
     else if (fAnaMode=="SECC") fHistReconstTrack->SetBinLimits(1,-3.4,-1.7);
     fHistReconstTrack->SetBinLimits(2,binning_cent_fmdfmd);
     fHistReconstTrack->SetBinLimits(3,-0.55*TMath::Pi(),1.45*TMath::Pi());
     fHistReconstTrack->SetBinLimits(4,binning_zvx);
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrack->SetVarTitle(1,"FMD(Trigger) Eta");
     fHistReconstTrack->SetVarTitle(2,"centrality");
     fHistReconstTrack->SetVarTitle(3,"#Delta#phi");
     fHistReconstTrack->SetVarTitle(4,"z vertex");
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
     fHistReconstTrackMix->SetBinLimits(0,-3.275,3.225);
     if(fAnaMode=="SECA") fHistReconstTrackMix->SetBinLimits(1,1.7,4.9);
     else if(fAnaMode=="SECC") fHistReconstTrackMix->SetBinLimits(1,-3.4,-1.7);
     fHistReconstTrackMix->SetBinLimits(2,binning_cent_fmdfmd);
     fHistReconstTrackMix->SetBinLimits(3,-0.55*TMath::Pi(),1.45*TMath::Pi());
     fHistReconstTrackMix->SetBinLimits(4,binning_zvx);
     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1,"FMD(Trigger) Eta");
     fHistReconstTrackMix->SetVarTitle(2,"centrality");
     fHistReconstTrackMix->SetVarTitle(3,"#Delta#phi");
     fHistReconstTrackMix->SetVarTitle(4,"z vertex");
   }else if(fAnaMode=="V0AV0C"){
     const Int_t nTrackVars_v0av0c = 4;
     const Int_t iTrackBin_VZEROAVZEROC[4]={10,10,11,8};
     fHistReconstTrack= new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars_v0av0c,iTrackBin_VZEROAVZEROC);
     fHistReconstTrack->SetBinLimits(0,binning_eta_vzero);
     fHistReconstTrack->SetBinLimits(1,binning_eta_vzero);
     fHistReconstTrack->SetBinLimits(2,binning_cent);
     fHistReconstTrack->SetBinLimits(3,binning_dphi_vzero);
     fHistReconstTrack->SetVarTitle(0,"Vzero(Asso) Eta");
     fHistReconstTrack->SetVarTitle(1,"Vzero(Trigger) Eta");
     fHistReconstTrack->SetVarTitle(2,"centrality");
     fHistReconstTrack->SetVarTitle(3,"#Delta#phi");
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars_v0av0c,iTrackBin_VZEROAVZEROC);
     fHistReconstTrackMix->SetBinLimits(0,binning_eta_vzero);
     fHistReconstTrackMix->SetBinLimits(1,binning_eta_vzero);
     fHistReconstTrackMix->SetBinLimits(2,binning_cent);
     fHistReconstTrackMix->SetBinLimits(3,binning_dphi_vzero);
     fHistReconstTrackMix->SetVarTitle(0,"Vzero(Asso) Eta");
     fHistReconstTrackMix->SetVarTitle(1,"Vzero(Trigger) Eta");
     fHistReconstTrackMix->SetVarTitle(2,"centrality");
     fHistReconstTrackMix->SetVarTitle(3,"#Delta#phi");
   }else if(fAnaMode=="ITSFMD" || fAnaMode=="ITSFMDC"){
       Double_t binning_detaFMDITS[66]={
		 -6.6,-6.5,-6.4,-6.3,-6.2,-6.1,-6.0,
		 -5.9,-5.8,-5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.0,
		 -4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,
		 -3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,
		 -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,
		 -1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,
		 -0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1};
	   Double_t binning_detaFMDCITS[50]={
		 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
		 1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
		 2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
		 3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
		 4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
		 5.0};
	   Double_t binning_dphi_itsfmd[37] = {
		 -1.570796,  -1.396263,  -1.221730, 
		 -1.047198,  -0.872665,  -0.698132, 
		 -0.523599,  -0.349066,  -0.174533,  
		 0.0,         0.174533,   0.349066,  
		 0.523599,    0.698132,   0.872665, 
		 1.047198,    1.221730,   1.396263, 
		 1.570796,    1.745329,   1.919862,  
		 2.094395,    2.268928,   2.443461,  
		 2.617994,    2.792527,   2.967060,  
		 3.141593,    3.316126,   3.490659,  
		 3.665191,    3.839724,   4.014257,  
		 4.188790,    4.363323,   4.537856,  
		 4.712389};
	   Double_t binning_cent_its[8]={0.,5.,10.,20.,40.,60.,80.,100.};
	   Int_t nbinitsfmddeltaeta;
	   Int_t nbinetafmd;
	   if(fAnaMode=="ITSFMD"){
		 nbinitsfmddeltaeta=65;
		 nbinetafmd=32;
	   }else{// if(fAnaMode=="ITSFMDC"){
		 nbinitsfmddeltaeta=49;
		 nbinetafmd=17;
	   }
       const Int_t iTrackBin_tpcfmd[6]={nbinitsfmddeltaeta,7,nbinetafmd,36,20,6};
	   Double_t binning_zvx_fmdits[2]={-10,10};
	   //	   Double_t binning_itseta[19]={-1.7, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.,0.2,  0.4,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6,  1.7};
	   Double_t binning_itseta[7]={-1.7,-1.2,-0.6,0.,0.6,1.2,1.7};

	   fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 6, iTrackBin_tpcfmd);
	   if(fAnaMode=="ITSFMD")  {
	     fHistReconstTrack->SetBinLimits(0,binning_detaFMDITS);
	     fHistReconstTrack->SetBinLimits(2,binning_etafmd);
	   }else if(fAnaMode=="ITSFMDC") {
	     fHistReconstTrack->SetBinLimits(0,binning_detaFMDCITS);
	     fHistReconstTrack->SetBinLimits(2,binning_etafmdc);
	 }
	   fHistReconstTrack->SetBinLimits(1,binning_cent_its);
	   fHistReconstTrack->SetBinLimits(3,binning_dphi_itsfmd);
	   fHistReconstTrack->SetBinLimits(4,-10.,10.);
	   fHistReconstTrack->SetBinLimits(5,binning_itseta);
	   fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
	   fHistReconstTrack->SetVarTitle(2,"FMD Eta");
	   fHistReconstTrack->SetVarTitle(1,"centrality");
	   fHistReconstTrack->SetVarTitle(3,"#Delta#phi");
	   fHistReconstTrack->SetVarTitle(4,"z vertex");
	   fHistReconstTrack->SetVarTitle(5,"ITS Eta");
	   
	   fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 6,iTrackBin_tpcfmd);
	   if(fAnaMode=="ITSFMD") {
	     fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDITS);
	     fHistReconstTrackMix->SetBinLimits(2,binning_etafmd);
	   }else if(fAnaMode=="ITSFMDC") {
	     fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDCITS);
	     fHistReconstTrackMix->SetBinLimits(2,binning_etafmdc);
	   }
	   fHistReconstTrackMix->SetBinLimits(1,binning_cent_its);
	   fHistReconstTrackMix->SetBinLimits(3,binning_dphi_itsfmd);
	   fHistReconstTrackMix->SetBinLimits(4,-10.,10.);
	   fHistReconstTrackMix->SetBinLimits(5,binning_itseta);
	   fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
	   fHistReconstTrackMix->SetVarTitle(1,"centrality");
	   fHistReconstTrackMix->SetVarTitle(2,"FMD Eta");
	   fHistReconstTrackMix->SetVarTitle(3,"#Delta#phi");
	   fHistReconstTrackMix->SetVarTitle(4,"z vertex");
	   fHistReconstTrackMix->SetVarTitle(5,"ITS Eta");
   }
   fOutputList1->Add(fHistReconstTrack);
   fOutputList1->Add(fHistReconstTrackMix);
   /*
 if(fAnaMode=="SP"){
   fHistQna=new TH2D("fHistQna","fHistQna",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQna);

   fHistQnc=new TH2D("fHistQnc","fHistQnc",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQnc);


   fHistQn=new TH2D("fHistQn","fHistQn",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQn);

   fHistQna_VZERO=new TH2D("fHistQna_VZERO","fHistQna_VZERO",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQna_VZERO);

   fHistQnc_VZERO=new TH2D("fHistQnc_VZERO","fHistQnc_VZERO",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQnc_VZERO);

   fHistQn_VZERO=new TH2D("fHistQn_VZERO","fHistQn_VZERO",200,0.,10.,10,0,100);
   fOutputList1->Add(fHistQn_VZERO);

   fHistVn=new TH1D("fHistVn","fHistVn",200,-1.,1.);
   fOutputList1->Add(fHistVn);

   for(Int_t i=0;i<4;i++){
     fHistQAQB[i]=new TH1D(Form("fHistQAQB_%d",i),Form("fHistQAQB_%d",i),400,-2.,2.);
     fOutputList1->Add(fHistQAQB[i]);
     fHistQAQB_VZERO[i]=new TH1D(Form("fHistQAQB_VZERO_%d",i),Form("fHistQAQB_VZERO_%d",i),400,-2.,2.);
     fOutputList1->Add(fHistQAQB_VZERO[i]);
     fHistCorrQna[i]=new TH2D(Form("fHistCorrQna_%d",i),"fHistCorrQna",200,0.,10.,200,0,10);
     fOutputList1->Add(fHistCorrQna[i]);
     fHistCorrQnc[i]=new TH2D(Form("fHistCorrQnc_%d",i),"fHistCorrQnc",200,0.,10.,200,0,10);
     fOutputList1->Add(fHistCorrQnc[i]);
   }

   Double_t binning_cent_QAQC[5]={0.,20.,40.,60.,100.};
   SP_TPCATPCC = new TProfile("SP_TPCATPCC","QAQC",4,binning_cent_QAQC,-3,+3,"s");
   fOutputList1->Add(SP_TPCATPCC);
   SP_TPCATPCC_default = new TProfile("SP_TPCATPCC_default","QAQC",4,binning_cent_QAQC,-3,+3);
   fOutputList1->Add(SP_TPCATPCC_default);

   SP_V0AV0C_default = new TProfile("SP_V0AV0C_default","QAQC",4,binning_cent_QAQC,-3,+3);
   fOutputList1->Add(SP_V0AV0C_default);
   SP_V0ATPC_default = new TProfile("SP_V0ATPC_default","QAQC",4,binning_cent_QAQC,-3,+3);
   fOutputList1->Add(SP_V0ATPC_default);
   SP_V0CTPC_default = new TProfile("SP_V0CTPC_default","QAQC",4,binning_cent_QAQC,-3,+3);
   fOutputList1->Add(SP_V0CTPC_default);

   fHist_V0AV0C = new TH1F("fHist_V0AV0C","QAQC",200,-1,1);
   fOutputList1->Add(fHist_V0AV0C);
   fHist_V0ATPC= new TH1F("fHist_V0ATPC","QAQC",200,-1,1);
   fOutputList1->Add(fHist_V0ATPC);
   fHist_V0CTPC = new TH1F("fHist_V0CTPC","QAQC",200,-1,1);
   fOutputList1->Add(fHist_V0CTPC);


   SP_uTPCA = new TProfile("SP_uTPCA","u x Q_{TPCA}",11,binning_pt_assoc,-3,+3);
   fOutputList1->Add(SP_uTPCA);

   SP_uTPCC = new TProfile("SP_uTPCC","u x Q_{TPCC}",11,binning_pt_assoc,-3,+3);
   fOutputList1->Add(SP_uTPCC);
   Int_t nbin_uTPC=11;
   Double_t binning_pt_assoc_uTPC[12] = {0.3, 0.5, 0.75, 1.0, 1.25, 1.5,2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
   Int_t     nbin_uTPCPhi=4;
   Double_t binning_pt_assoc_uTPCPhi[5] = {0., 0.5, 2.0, 4.0, 8.0};
   Int_t nbin_uTPCCas=3;
   Double_t binning_pt_assoc_uTPCCas[4] = {0., 1., 4., 8.};
   for(Int_t i=0;i<8;i++){
     SP_uVZEROA_PP[i] = new TProfile(Form("SP_uVZEROA_PP_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROA_PP[i]);
     SP_uVZEROA[i] = new TProfile(Form("SP_uVZEROA_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROA[i]);
     SP_uVZEROA1[i] = new TProfile(Form("SP_uVZEROA1_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROA1[i]);
     SP_uVZEROA2[i] = new TProfile(Form("SP_uVZEROA2_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROA2[i]);
     SP_uVZEROA3[i] = new TProfile(Form("SP_uVZEROA3_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROA3[i]);
     SP_uVZEROC_PP[i] = new TProfile(Form("SP_uVZEROC_PP_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROC_PP[i]);
     SP_uVZEROC[i] = new TProfile(Form("SP_uVZEROC_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROC[i]);
     SP_uVZEROC1[i] = new TProfile(Form("SP_uVZEROC1_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROC1[i]);
     SP_uVZEROC2[i] = new TProfile(Form("SP_uVZEROC2_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROC2[i]);
     SP_uVZEROC3[i] = new TProfile(Form("SP_uVZEROC3_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
     fOutputList1->Add(SP_uVZEROC3[i]);
   }
   if(fasso=="PID" || fasso=="hadron" || fasso=="V0"){
     for(Int_t i=0;i<8;i++){
       SP_uTPC_PP[i] = new TProfile(Form("SP_uTPC_PP_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
       fOutputList1->Add(SP_uTPC_PP[i]);
       SP_uTPC[i] = new TProfile(Form("SP_uTPC_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
       fOutputList1->Add(SP_uTPC[i]);
       SP_uTPC1[i] = new TProfile(Form("SP_uTPC1_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
       fOutputList1->Add(SP_uTPC1[i]);
       SP_uTPC2[i] = new TProfile(Form("SP_uTPC2_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
       fOutputList1->Add(SP_uTPC2[i]);
       SP_uTPC3[i] = new TProfile(Form("SP_uTPC3_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPC,-3,+3);
       fOutputList1->Add(SP_uTPC3[i]);
   }
   }else if(fasso=="Phi"){
     for(Int_t i=0;i<8;i++){
       SP_uTPC_PP[i] = new TProfile(Form("SP_uTPC_PP_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPCPhi,-3,+3);
       fOutputList1->Add(SP_uTPC_PP[i]);
       SP_uTPC[i] = new TProfile(Form("SP_uTPC_%d",i),"u x Q_{TPC}",nbin_uTPCPhi,binning_pt_assoc_uTPCPhi,-3,+3);
       fOutputList1->Add(SP_uTPC[i]);
       SP_uTPC1[i] = new TProfile(Form("SP_uTPC1_%d",i),"u x Q_{TPC}",nbin_uTPCPhi,binning_pt_assoc_uTPCPhi,-3,+3);
       fOutputList1->Add(SP_uTPC1[i]);
       SP_uTPC2[i] = new TProfile(Form("SP_uTPC2_%d",i),"u x Q_{TPC}",nbin_uTPCPhi,binning_pt_assoc_uTPCPhi,-3,+3);
       fOutputList1->Add(SP_uTPC2[i]);
       SP_uTPC3[i] = new TProfile(Form("SP_uTPC3_%d",i),"u x Q_{TPC}",nbin_uTPCPhi,binning_pt_assoc_uTPCPhi,-3,+3);
       fOutputList1->Add(SP_uTPC3[i]);
     }
   }else if(fasso=="Cascade"){
      for(Int_t i=0;i<8;i++){
       SP_uTPC_PP[i] = new TProfile(Form("SP_uTPC_PP_%d",i),"u x Q_{TPC}",nbin_uTPC,binning_pt_assoc_uTPCCas,-3,+3);
       fOutputList1->Add(SP_uTPC_PP[i]);
       SP_uTPC[i] = new TProfile(Form("SP_uTPC_%d",i),"u x Q_{TPC}",nbin_uTPCCas,binning_pt_assoc_uTPCCas,-3,+3);
       fOutputList1->Add(SP_uTPC[i]);
       SP_uTPC1[i] = new TProfile(Form("SP_uTPC1_%d",i),"u x Q_{TPC}",nbin_uTPCCas,binning_pt_assoc_uTPCCas,-3,+3);
       fOutputList1->Add(SP_uTPC1[i]);
       SP_uTPC2[i] = new TProfile(Form("SP_uTPC2_%d",i),"u x Q_{TPC}",nbin_uTPCCas,binning_pt_assoc_uTPCCas,-3,+3);
       fOutputList1->Add(SP_uTPC2[i]);
       SP_uTPC3[i] = new TProfile(Form("SP_uTPC3_%d",i),"u x Q_{TPC}",nbin_uTPCCas,binning_pt_assoc_uTPCCas,-3,+3);
       fOutputList1->Add(SP_uTPC3[i]);
     }
   }
   }
   */


 }

 void AliAnalysisTaskSEpPbCorrelationsMCYS::UserExec(Option_t *) {
   DumpTObjTable("Start analysis");

   AliAnalysisManager *mgr        = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inEvMain = (AliInputEventHandler *)(mgr->GetInputEventHandler());
   if (!inEvMain)    return;
   
   fPIDResponse = inEvMain->GetPIDResponse();
   if (!fPIDResponse)    return;

   fEvent = dynamic_cast<AliAODEvent *>(inEvMain->GetEvent());
   if (!fEvent) {
     AliWarning("ERROR: fEvent not available \n");
     return;
   }

   
   
   fHist_Stat->Fill(0);
   if(fcollisiontype=="pPb"|| fcollisiontype=="PbPb" || fcollisiontype=="PP"){
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
	 else  isSelected = ((maskIsSelected & AliVEvent::kHighMultSPD)== AliVEvent::kHighMultSPD);
	 //	 isSelected = ((maskIsSelected & AliVEvent::kHighMultV0)== AliVEvent::kHighMultV0);//Both for data and 
	 //	 isSelected = ((maskIsSelected & AliVEvent::kINT7)== AliVEvent::kINT7);//Both for data and 
	 //	 cout<<"Entry=="<<fNEntries<<endl;
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
   AliMultSelection*multSelection;
   if(!fcentcalib){
     multSelection= (AliMultSelection *)fEvent->FindListObject("MultSelection");
     //     AliMultSelectionTask::SetPreferSuperCalib(kTRUE);//recommended by Ionut
     if(!multSelection) return;
     fHist_Stat->Fill(2);
   }
   
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

   //   if(fcollisiontype.Contains("HMPP")AliMultSelectionTask::SetHighMultQABinning(kTRUE);

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

   // Multiplicity Object
   //if(fcollisiontype=="pPb" || fcollisiontype=="PP"){
   fvzero = fEvent->GetVZEROData();

   if(!fcentcalib){
     lCentrality = multSelection->GetMultiplicityPercentile(fCentType);
     Int_t qual = multSelection->GetEvSelCode();
     if (qual == 199)  lCentrality = -999;
   } else{
     Float_t sum = 0., max = 0.;
     for(Int_t i = 32; i < 64; ++i)
       {      sum +=fvzero->GetMultiplicity(i);
	 if (fvzero->GetMultiplicity(i) > max) max = fvzero->GetMultiplicity(i);
       }
     sum -= max;
     //     fV0Amultmodi->Fill(sum);

     Int_t nbinmult= fhcorr[0]->GetXaxis()->FindBin(sum);
     lCentrality=fhcorr[0]->GetBinContent(nbinmult);
   }      
     Float_t sum = 0., max = 0.;
     for(Int_t i = 32; i < 64; ++i)
       {      sum +=fvzero->GetMultiplicity(i);
	 if (fvzero->GetMultiplicity(i) > max) max = fvzero->GetMultiplicity(i);
       }
     sum -= max;
     //     fV0Amultmodi->Fill(sum);
     
     Float_t nV0A_hits = fvzero->GetMTotV0A();
     fh2_V0A_comp->Fill(sum,nV0A_hits);     

     /*
     Float_t v0amult=fvzero->GetMTotV0A();
     Int_t nbinmult= fhcorr[0]->GetXaxis()->FindBin(v0amult);
     lCentrality=fhcorr[0]->GetBinContent(nbinmult);
     */
     
     AliAODMCHeader* aodMCheader=(AliAODMCHeader*)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
     TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
     if(!mcArray){
        Printf("No MC particle branch found");
        return;
      }
     Int_t nMCAllTracks = mcArray->GetEntriesFast();
      Int_t ntrackv0aprimary=0;
      Int_t ntrackv0aprimaryall=0;
      for (Int_t i = 0; i < nMCAllTracks; i++){
	AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
	if (!mcTrack) {
	  Error("ReadEventAODMC", "Could not receive particle %d", i);
	  continue;
	}
	Bool_t TrIsPrim=mcTrack->IsPhysicalPrimary();
	Float_t mcTrackEta = mcTrack->Eta();
	Bool_t TrCharge=mcTrack->Charge()!=0;
	if(!TrCharge)        continue;
	if(mcTrackEta>2.8 && mcTrackEta<5.1) ntrackv0aprimaryall++;
	if(!TrIsPrim)	     continue;
	if(mcTrackEta>2.8 && mcTrackEta<5.1) ntrackv0aprimary++;
      }
      
      fV0Amultprim->Fill(ntrackv0aprimary);
      fh2_V0A_comp_prim->Fill(ntrackv0aprimaryall,sum);     
      
      if (lCentrality < 0. || lCentrality > 100. - 0.0000001) {
	return;
      }
      Double_t *CentBins = fCentBins;
      poolmin = CentBins[0];
      poolmax = CentBins[fNCentBins];
      fHist_Stat->Fill(4);
      
      fHistCentrality_beforecut->Fill(lCentrality);
      
      //   fHist_Stat->Fill(5);
      DumpTObjTable("After event selection");


      MakeAna();
      
      PostData(1, fOutputList);
      PostData(2, fOutputList1);
      PostData(3, fOutputList2);
 }

void AliAnalysisTaskSEpPbCorrelationsMCYS::Terminate(Option_t *) {
  //  AliInfo(Form("Number of Correlation
  DumpTObjTable("End of the analysis");
  Printf("Entries======================%d",fNEntries);
  if (fPoolMgr)    delete fPoolMgr;   // PoolMgr->ClearPools();
  if (fPoolMgr1)    delete fPoolMgr1; // fPoolMgr1->ClearPools();
}

void AliAnalysisTaskSEpPbCorrelationsMCYS::MakeAna() {
  
  DumpTObjTable("start correlation analysis");
  TObjArray *selectedTracksLeading = new TObjArray;
  selectedTracksLeading->SetOwner(kTRUE);
  TObjArray *selectedTracksAssociated = new TObjArray;
  selectedTracksAssociated->SetOwner(kTRUE);
  
  TObjArray* selectedTracksMC1=new TObjArray;
  selectedTracksMC1->SetOwner(kTRUE);
  TObjArray* selectedTracksMC2=new TObjArray;
  selectedTracksMC2->SetOwner(kTRUE);
  
  Double_t eta_min;
  Double_t eta_max;
  Double_t eta_ave;
  Double_t phi_vzero;
  Double_t mult_vzero;
  Double_t vzeroqa[3];
  Double_t mult_vzero_eq;
  Float_t nV0A_hits_fmdacc=0;
  Float_t nV0C_hits_fmdacc=0;
  
  for (Int_t imod = 0; imod < 64; imod++) {
    eta_min = fvzero->GetVZEROEtaMin(imod);
    eta_max = fvzero->GetVZEROEtaMax(imod);
    phi_vzero = fvzero->GetVZEROAvgPhi(imod);
    mult_vzero = fvzero->GetMultiplicity(imod);
    mult_vzero_eq = fEvent->GetVZEROEqMultiplicity(imod);
    eta_ave = (eta_min + eta_max) / 2.;
    if(eta_ave>2.8 && eta_ave<5.03) nV0A_hits_fmdacc+=mult_vzero_eq;
    else if(eta_ave>-3.4 && eta_ave<-2.01) nV0C_hits_fmdacc+=mult_vzero_eq;
    
    fHist_vzeromult->Fill(imod, mult_vzero);
    fHist_vzeromultEqweighted->Fill(imod, mult_vzero_eq);
    fHist2dmult->Fill(imod, mult_vzero_eq, mult_vzero);
    
    vzeroqa[0] = eta_ave;
    vzeroqa[1] = phi_vzero;
    vzeroqa[2] = lCentrality;
    fHistVZERO->Fill(vzeroqa, 0, (Double_t)mult_vzero_eq);
    if(imod>31) {
      if(fAnaMode=="TPCV0A") selectedTracksAssociated->Add(new AliAssociatedTrackYSMC(-999,eta_ave,phi_vzero,-999,-999,-999,-999,-999,mult_vzero_eq));
      if(fAnaMode=="V0AV0C")selectedTracksLeading->Add(new AliAssociatedTrackYSMC(-999,eta_ave,phi_vzero,-999,-999,-999,-999,-999,mult_vzero_eq));
    }else if(imod<32) {
      if(fAnaMode=="TPCV0C" ||fAnaMode=="V0AV0C")selectedTracksAssociated->Add(new AliAssociatedTrackYSMC(-999,eta_ave,phi_vzero,-999,-999,-999,-999,-999,mult_vzero_eq));
    }
  }
   Float_t nFMD_fwd_hits=0;
   Float_t nFMD_bwd_hits=0;
   Float_t nFMD_fwdV0acc_hits=0;
   Float_t nFMD_bwdV0acc_hits=0;
   
   Float_t nV0A_hits = fvzero->GetMTotV0A();
   Float_t nV0C_hits = fvzero->GetMTotV0C();


   
   AliAODForwardMult*aodForward=static_cast<AliAODForwardMult*>(fEvent->FindListObject("Forward"));
   /*
     if(!aodForward) {
     selectedTracksLeading->Clear();
     delete selectedTracksLeading;
     selectedTracksAssociated->Clear();
     delete selectedTracksAssociated;
     PostData(1, fOutputList);
     PostData(2, fOutputList1);
     PostData(3, fOutputList2);
     return;
     }*/
   if(aodForward){
     // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
     Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
     const TH2D& d2Ndetadphi = aodForward->GetHistogram();
     TH1*hphiacceptance=aodForward->GetPhiAcceptance();
     Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
     Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
     Double_t pt = 0;
     for (Int_t iEta = 1; iEta <= nEta; iEta++) {
	 Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
	 if (!valid) {
	   continue;
	 }
	 
	 Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
	 Float_t phiacc=hphiacceptance->GetBinContent(iEta);
	 for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
	   // Bin content is most likely number of particles!
	   Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	 
	   Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
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
	     fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	   }
	 }
       }
     
       //delete hphiacceptance;
       
       if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0){
	 selectedTracksLeading->Clear();
	 delete selectedTracksLeading;
	 selectedTracksAssociated->Clear();
	 delete selectedTracksAssociated;
	 PostData(1, fOutputList);
	 PostData(2, fOutputList1);
	 PostData(3, fOutputList2);
	 return;
       }
       
       fHist_Stat->Fill(5);
       
       DumpTObjTable("End of fill fmd tracks");
       
       fFMDV0->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
       fFMDV0A->Fill(nFMD_fwd_hits, nV0A_hits);
       fFMDV0C->Fill(nFMD_bwd_hits, nV0C_hits);
       
       /*
	 fHist_NeventRun->Fill(ConvertRunNumber(fEvent->GetRunNumber()));
	 fHist_V0AMultRun->Fill(ConvertRunNumber(fEvent->GetRunNumber()),nV0A_hits);
	 fHist_V0CMultRun->Fill(ConvertRunNumber(fEvent->GetRunNumber()),nV0C_hits);
	 fHist_FMDAMultRun->Fill(ConvertRunNumber(fEvent->GetRunNumber()),nFMD_fwd_hits);
	 fHist_FMDCMultRun->Fill(ConvertRunNumber(fEvent->GetRunNumber()),nFMD_bwd_hits);
       */
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
	 case 6:
	   FMDcutapar0=1.75;
	   FMDcutapar1=150;
	   FMDcutcpar0=1.4;
	   FMDcutcpar1=120;
	   break;
	 case 7:
	   FMDcutapar0=1.64755;
	   FMDcutapar1=119.602;
	   FMDcutcpar0=2.73426;
	   FMDcutcpar1=150.31;
	   break;
	 case 8:
	   FMDcutapar0=1.64755;
	   FMDcutapar1=159.47;
	   FMDcutcpar0=2.73426;
	   FMDcutcpar1=200.413;
	   break;
	 case 9:
	   FMDcutapar0=1.2031;
	   FMDcutapar1=73.123;
	   FMDcutcpar0=2.25453;
	   FMDcutcpar1=104.941;
	   break;
	 case 10:
	   FMDcutapar0=1.2031;
	   FMDcutapar1=97.4973;
	   FMDcutcpar0=2.25453;
	   FMDcutcpar1=139.921;
	   break;
	 case 12:
	   FMDcutapar0=0;
	   FMDcutapar1=0;
	   FMDcutcpar0=0;
	   FMDcutcpar1=0;
	 default: break;
	 }
	 
	 if(fcollisiontype=="PbPb") {
	   if ((nV0A_hits_fmdacc + nV0C_hits_fmdacc) < 1.5*(nFMD_fwdV0acc_hits + nFMD_bwdV0acc_hits) - 20) {
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
     }
     
       fFMDV0_post->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
       fFMDV0A_post->Fill(nFMD_fwd_hits, nV0A_hits);
       fFMDV0C_post->Fill(nFMD_bwd_hits, nV0C_hits);
       
       
       for (Int_t iEta = 1; iEta <= nEta; iEta++) {
	 Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
	 if (!valid) {
	   continue;
	 }
	 Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
	 Float_t phiacc=hphiacceptance->GetBinContent(iEta);
	 fhistfmdphiacc->Fill(eta,lCentrality,phiacc);
	 for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
	   Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	   Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
	   Double_t cont[4]={eta,phi,lCentrality,fPrimaryZVtx};
	   fhistfmd->Fill(cont,0,mostProbableN);
	   if (mostProbableN > 0) {
	     if(eta>0){
	       if(fAnaMode=="TPCFMD" || fAnaMode=="ITSFMD") selectedTracksAssociated->Add(new AliAssociatedTrackYSMC(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));			
	       if(fAnaMode=="FMDFMD") selectedTracksLeading->Add(new AliAssociatedTrackYSMC(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));	
	       
	     }else if(eta<0){
	       if(fAnaMode=="TPCFMDC" || fAnaMode=="ITSFMDC" ||fAnaMode=="FMDFMD") selectedTracksAssociated->Add(new AliAssociatedTrackYSMC(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));
	     }
	   }
	 }
       }
       
       delete hphiacceptance;
       
     
   }
   fHist_Stat->Fill(6);
   fHistCentrality->Fill(lCentrality);
   fHistzvertex->Fill(tPrimaryVtxPosition[2]);
   fHistCentzvertex->Fill(lCentrality, tPrimaryVtxPosition[2]);


   

   
   DumpTObjTable("End of FMD vs V0 cuts");
   

if(fAnaMode=="TPCTPC"){
  if(fasso=="hadron") selectedTracksAssociated=GetAcceptedTracksLeading(fEvent,kFALSE,selectedTracksAssociated);
  else if (fasso == "Phi")    selectedTracksAssociated = GetAcceptedTracksAssociated(fEvent);
  else if (fasso == "V0")    selectedTracksAssociated = GetAcceptedV0Tracks(fEvent);
  else if (fasso == "PID")    selectedTracksAssociated = GetAcceptedTracksPID(fEvent);
  else if (fasso == "Cascade")    selectedTracksAssociated = GetAcceptedCascadeTracks(fEvent);
 }
// Leading Particle
 if(fAnaMode=="TPCFMD" || fAnaMode=="TPCTPC" || fAnaMode=="TPCFMDC"){
   selectedTracksLeading=GetAcceptedTracksLeading(fEvent,kTRUE,selectedTracksLeading);
 }
 
 DumpTObjTable("End of TPC/ITS track fill");
 Int_t pdgcode=0;
 Double_t conmcprim[5];
 Bool_t TrIsPrim=kFALSE;
 Bool_t TrIsSecondMate=kFALSE;
 Bool_t TrIsSecondWeak=kFALSE;
 Bool_t TrIsOthers=kFALSE;
 Double_t mcTrackEta=-999.;
 Double_t mcTrackPt=-999.;
 Double_t mcTrackPhi=-999.;
 Bool_t TrCharge=kFALSE;
 
  
 AliAODMCHeader* aodMCheader=(AliAODMCHeader*)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
 TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
 if(!mcArray){
   Printf("No MC particle branch found");
   return;
 }
 
 Int_t nMCAllTracks = mcArray->GetEntriesFast();
 Int_t nMCtrackssamecut=0;
 Int_t ntrackv0aprimary=0;
 Int_t ntrackv0aall=0;
 Int_t ntrackv0cprimary=0;

 for (Int_t i = 0; i < nMCAllTracks; i++){
   AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
   if (!mcTrack) {
     Error("ReadEventAODMC", "Could not receive particle %d", i);
     continue;
   }
   
   TrIsPrim=mcTrack->IsPhysicalPrimary();
   TrIsSecondMate=mcTrack->IsSecondaryFromMaterial();
   TrIsSecondWeak=mcTrack->IsSecondaryFromWeakDecay();
   TrIsOthers=!TrIsPrim && !TrIsSecondMate && !TrIsSecondWeak;
   
   mcTrackEta = mcTrack->Eta();
   mcTrackPt  = mcTrack->Pt();
   mcTrackPhi = mcTrack->Phi();
   
   
   TrCharge=mcTrack->Charge()!=0;
   pdgcode=TMath::Abs(mcTrack->PdgCode());
   if(!TrCharge)        continue;
   if(mcTrackEta>2.8 && mcTrackEta<5.1) ntrackv0aall++;
   conmcprim[0]=mcTrackPt;
   conmcprim[1]=mcTrackEta;
   conmcprim[2]=mcTrackPhi;
   conmcprim[3]=lCentrality;
   conmcprim[4]=fPrimaryZVtx;
   
   if(TrIsPrim){
     if(fQA) if(mcTrackEta>-1. && mcTrackEta<1.) fhistmcprim->Fill(conmcprim,0);//primay charged partilce distribution(no mother particle)
     if(mcTrackEta>-0.8&& mcTrackEta<0.8) {
       if(mcTrackPt>fPtMin && mcTrackPt<fPtMax) {
	 nMCtrackssamecut++;
       }
     }
     fhmcprimvzeta->Fill(mcTrackEta,mcTrack->Zv());
     fhmcprimpdgcode->Fill(pdgcode);
     fh2_FMD_eta_phi_prim->Fill(mcTrackEta,mcTrackPhi);
     
     fhistmcprimfinal->Fill(conmcprim,0);
     Double_t mcTrackEta1=mcTrackEta;
     
     if(mcTrackEta>2.8 && mcTrackEta<5.1) ntrackv0aprimary++;
     if(mcTrackEta>-3.7 && mcTrackEta<-1.7) ntrackv0cprimary++;

     if(fAnaMode=="TPCTPC") {
       if(mcTrackEta<-0.8 || mcTrackEta>0.8) continue;
       if(mcTrackPt<fPtMin || mcTrackPt>fPtMax) continue;
       selectedTracksMC1->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));         
       selectedTracksMC2->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));  			
     }else if(fAnaMode=="TPCFMD"||fAnaMode=="TPCFMDC"){
       if(mcTrackEta>-0.8 && mcTrackEta<0.8){
	 if(mcTrackPt<fPtMin || mcTrackPt>fPtMax) continue;
	 selectedTracksMC1->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
       }else{
	 if(fAnaMode=="TPCFMD"){
	   if(mcTrackEta>1.7 && mcTrackEta<4.9)  selectedTracksMC2->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));  
	 }else if(fAnaMode=="TPCFMDC" ||fAnaMode=="FMDFMD"){
	   if(mcTrackEta>-3.4  && mcTrackEta<-1.7) selectedTracksMC2->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
	 }
       }
     }else if(fAnaMode=="FMDFMD"){
       if(mcTrackEta>1.7 && mcTrackEta<4.9) selectedTracksMC1->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
       if(mcTrackEta>-3.4  && mcTrackEta<-1.7) selectedTracksMC2->Add(new AliAssociatedTrackYSMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
     }
   }			
 }	  
 

 Int_t nrecotrack=selectedTracksLeading->GetEntriesFast();
 // cout<<nrecotrack<<endl;
 //   return;              

 fNTrackCorrMC->Fill(nrecotrack,nMCtrackssamecut); 

 fh2_V0A_all->Fill(ntrackv0aall,nV0A_hits);
 fh2_V0A->Fill(ntrackv0aprimary,nV0A_hits);
 fh2_V0C->Fill(ntrackv0cprimary,nV0A_hits);



 if(ffillcorrelation){
   if(!fprim){
   if(fextractsec){
     FillCorrelationTracks(lCentrality,selectedTracksMC1,selectedTracksAssociated,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
     FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksMC1,selectedTracksAssociated,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
   }else{
     FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksAssociated,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
     FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedTracksAssociated,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
   }
     DumpTObjTable("End of fill  Correlation");
  }else{
    FillCorrelationTracks(lCentrality,selectedTracksMC1,selectedTracksMC2,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksMC1,selectedTracksMC2,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
    DumpTObjTable("End of fill  Correlation");
  }

 }
  
  selectedTracksLeading->Clear();
  delete selectedTracksLeading;
  selectedTracksAssociated->Clear();
  delete selectedTracksAssociated;
  selectedTracksMC1->Clear();
  delete selectedTracksMC1;
  selectedTracksMC2->Clear();
  delete selectedTracksMC2;
  
  DumpTObjTable("after delete TObjects");
  fNEntries++;
 }

TObjArray* AliAnalysisTaskSEpPbCorrelationsMCYS::GetFMDhitsYS(Bool_t Aside){
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
	  tracks1->Add(new AliAssociatedVZEROYSMC(mostProbableN,eta,phi,0,0,0));
	  Double_t cont[3]={eta,phi,lCentrality};
	  fhistfmd->Fill(cont,0,mostProbableN);
	  fh2_FMD_acceptance->Fill(eta,tPrimaryVtxPosition[2]);
	  fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	}
      }
    }
    


    /*
      TObjArray *tracks1 = new TObjArray;
      tracks1->SetOwner(kTRUE);
      Int_t i=0;
    for (auto const &track: ret_vector) {
           //cout<<i<<" "<<track.eta<<" "<<track.phi<<" "<<track.weight<<endl;
           if(Aside==kTRUE){
            if(track.eta>0) tracks1->Add(new AliAssociatedVZEROYSMC(track.eta,track.phi,track.weight,0,0,0));
            } else{
            if(track.eta<0) tracks1->Add(new AliAssociatedVZEROYSMC(track.eta,track.phi,track.weight,0,0,0));
          }
           i++;
    }
    */
        return tracks1;
    //  std::random_device rd;
    //    std::default_random_engine engine{rd()};
    //    std::shuffle(std::begin(ret_vector), std::end(ret_vector), engine);

    //    return ret_vector;
  }


void AliAnalysisTaskSEpPbCorrelationsMCYS::CalculateSP(){
/*
  //Scalar Product
  Int_t fHarmonic=2;
  Double_t fQTPCCCos,fQTPCCSin,fQTPCC;
  Double_t fQTPCACos,fQTPCASin,fQTPCA;
  Double_t fQTPCCos,fQTPCSin,fQTPC;
  Double_t fQTPCCCosArray[3],fQTPCCSinArray[3];
  Double_t fQTPCACosArray[3],fQTPCASinArray[3];
  fQTPCCCos=0;
  fQTPCCSin=0;
  fQTPCC=0;
  fQTPCACos=0;
  fQTPCASin=0;
  fQTPCA=0;
  fQTPCCos=0;
  fQTPCSin=0;
  fQTPC=0;
  for(Int_t i=0;i<selectedTracksLeading->GetEntriesFast();i++){
    AliAssociatedTrackYSMC* trackRP=(AliAssociatedTrackYSMC*)selectedTracksLeading->At(i);
    if(!trackRP) continue;
    Double_t phiRP=trackRP->Phi();
    Double_t ptRP=trackRP->Pt();
    Double_t etaRP=trackRP->Eta();
    Double_t w=1.;
    fQTPCCos+=w*TMath::Cos(fHarmonic*phiRP);
    fQTPCSin+=w*TMath::Sin(fHarmonic*phiRP);
    fQTPC+=w;
    if(etaRP<0.4 && etaRP>-0.4) continue;
    if(etaRP<0){
      fQTPCCCos+=w*TMath::Cos(fHarmonic*phiRP);
      fQTPCCSin+=w*TMath::Sin(fHarmonic*phiRP);
      fQTPCC+=w;
    }else{
      fQTPCACos+=w*TMath::Coss(fHarmonic*phiRP);
      fQTPCASin+=w*TMath::Sin(fHarmonic*phiRP);
      fQTPCA+=w;hi
    }
    if(etaRP<0) {
      fQTPCCCosArray[0]=fQTPCCCos;
      fQTPCCSinArray[0]=fQTPCCSin;
    }
    if(etaRP>0) {
      fQTPCACosArray[0]=fQTPCACos;
      fQTPCASinArray[0]=fQTPCASin;
    }
    if(etaRP<-0.4) {
      fQTPCCCosArray[1]=fQTPCfCCos;
      fQTPCCSinArray[1]=fQTPCCSin;
    }
    if(etaRP>0.4) {
      fQTPCACosArray[1]=fQTPCACos;
      fQTPCASinArray[1]=fQTPCASin;
    }

  }

  Double_t tpc_qmcos=fQTPCCos/fQTPC;
  Double_t tpc_qmsin=fQTPCSin/fQTPC;
  Double_t tpcc_qmcos=fQTPCCCos/fQTPCC;
  Double_t tpcc_qmsin=fQTPCCSin/fQTPCC;
  Double_t tpca_qmcos=fQTPCACos/fQTPCA;
  Double_t tpca_qmsin=fQTPCASin/fQTPCA;
  Double_t qna=TMath::Sqrt((fQTPCACos*fQTPCACos+fQTPCASin*fQTPCASin)/fQTPCA);
  if(fQTPCA>1)fHistQna->Fill(qna,lCentrality);
  Double_t qnc=TMath::Sqrt((fQTPCCCos*fQTPCCCos+fQTPCCSin*fQTPCCSin)/fQTPCC);
  if(fQTPCC>1)fHistQnc->Fill(qnc,lCentrality);
  Double_t qn=TMath::Sqrt((fQTPCCos*fQTPCCos+fQTPCSin*fQTPCSin)/fQTPC);
  if(fQTPC>1)fHistQn->Fill(qn,lCentrality);
  //  cout<<fQTPCCos<<" "<<fQTPCSin<<" "<<fQTPC<<" "<<qn<<endl;
  Double_t qaqc=tpca_qmcos*tpcc_qmcos+tpca_qmsin*tpcc_qmsin;//Q_a*Q_b/(M_a*M_b)
  //Double_t qaqcsqrt=TMath::Sqrt(tpca_qmcos*tpcc_qmcos+tpca_qmsin*tpcc_qmsin);//Q_a*Q_b/(M_a*M_b)
  if(fQTPCC>0 && fQTPCA>0) {
  fHistVn->Fill(qaqc);
  if(lCentrality>=0. && lCentrality<20.) fHistQAQB[0]->Fill(qaqc);
  if(lCentrality>=20. && lCentrality<40.) fHistQAQB[1]->Fill(qaqc);
  if(lCentrality>=40. && lCentrality<60.) fHistQAQB[2]->Fill(qaqc);
  if(lCentrality>=60. && lCentrality<100.) fHistQAQB[3]->Fill(qaqc);
  SP_TPCATPCC->Fill(lCentrality,qaqc);
  SP_TPCATPCC_default->Fill(lCentrality,qaqc);
  }

  //VZERO
  fvzero = fEvent->GetVZEROData();
  Double_t eta_min;
  Double_t eta_max;
  Double_t eta_ave;
  Double_t phi_vzero;
  Double_t mult_vzero;
  Double_t vzeroqa[3];
  Double_t mult_vzero_eq;
  Double_t fQVZEROCCos,fQVZEROCSin,fQVZEROC;
  Double_t fQVZEROACos,fQVZEROASin,fQVZEROA;
  Double_t fQVZEROCos,fQVZEROSin,fQVZERO;
  fQVZEROCCos=0;
  fQVZEROCSin=0;
  fQVZEROC=0;
  fQVZEROACos=0;
  fQVZEROASin=0;
  fQVZEROA=0;
  fQVZEROCos=0;
  fQVZEROSin=0;
  fQVZERO=0;
  for (Int_t imod = 0; imod < 64; imod++) {
    eta_min = fvzero->GetVZEROEtaMin(imod);
    eta_max = fvzero->GetVZEROEtaMax(imod);
    phi_vzero = fvzero->GetVZEROAvgPhi(imod);
    mult_vzero = fvzero->GetMultiplicity(imod);
    mult_vzero_eq = fEvent->GetVZEROEqMultiplicity(imod);
    eta_ave = (eta_min + eta_max) / 2.;
    fHist_vzeromult->Fill(imod, mult_vzero);
    fHist_vzeromultEqweighted->Fill(imod, mult_vzero_eq);
    fHist2dmult->Fill(imod, mult_vzero_eq, mult_vzero);
    vzeroqa[0] = eta_ave;
    vzeroqa[1] = phi_vzero;
    vzeroqa[2] = lCentrality;
    if (fQA)   fHistVZERO->Fill(vzeroqa, 0, (Double_t)mult_vzero_eq);
    Double_t phiRPVZERO=phi_vzero;
    Double_t etaRPVZERO=eta_ave;
    Double_t w=mult_vzero;
    fQVZEROCos+=w*TMath::Cos(fHarmonic*phiRPVZERO);
    fQVZEROSin+=w*TMath::Sin(fHarmonic*phiRPVZERO);
    fQVZERO+=w;
    if(etaRPVZERO<0){
      fQVZEROCCos+=w*TMath::Cos(fHarmonic*phiRPVZERO);
      fQVZEROCSin+=w*TMath::Sin(fHarmonic*phiRPVZERO);
      fQVZEROC+=w;
    }else{
      fQVZEROACos+=w*TMath::Cos(fHarmonic*phiRPVZERO);
      fQVZEROASin+=w*TMath::Sin(fHarmonic*phiRPVZERO);
      fQVZEROA+=w;
    }

    if(imod>31) selectedTrackV0A->Add(new AliAssociatedVZEROYSMC(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
    if(imod<32) selectedTrackV0C->Add(new AliAssociatedVZEROYSMC(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
    }
  Double_t vzeroc_qmcos=fQVZEROCCos/fQVZEROC;
  Double_t vzeroc_qmsin=fQVZEROCSin/fQVZEROC;
  Double_t vzeroa_qmcos=fQVZEROACos/fQVZEROA;
  Double_t vzeroa_qmsin=fQVZEROASin/fQVZEROA;
  Double_t qna_vzero=TMath::Sqrt((fQVZEROACos*fQVZEROACos+fQVZEROASin*fQVZEROASin)/fQVZEROA);
  if(fQVZEROA>1) fHistQna_VZERO->Fill(qna_vzero,lCentrality);
  Double_t qnc_vzero=TMath::Sqrt((fQVZEROCCos*fQVZEROCCos+fQVZEROCSin*fQVZEROCSin)/fQVZEROC);
  if(fQVZEROC>1) fHistQnc_VZERO->Fill(qnc_vzero,lCentrality);
  Double_t qn_vzero=TMath::Sqrt((fQVZEROCos*fQVZEROCos+fQVZEROSin*fQVZEROSin)/fQVZERO);
  if(fQVZERO>1)  fHistQn_VZERO->Fill(qn_vzero,lCentrality);
  Double_t qaqc_vzero=vzeroa_qmcos*vzeroc_qmcos+vzeroa_qmsin*vzeroc_qmsin;//Q_a*Q_b/(M_a*M_b)
  Double_t qaqc_vzeroatpc=vzeroa_qmcos*tpc_qmcos+vzeroa_qmsin*tpc_qmsin;//Q_a*Q_b/(M_a*M_b)
  Double_t qaqc_vzeroctpc=vzeroc_qmcos*tpc_qmcos+vzeroc_qmsin*tpc_qmsin;//Q_a*Q_b/(M_a*M_b)
  if(fQVZEROC>1 && fQTPCC>1) {
    if(lCentrality>=0. && lCentrality<20.)    fHistCorrQnc[0]->Fill(qnc,qnc_vzero);
    if(lCentrality>=20. && lCentrality<40.) fHistCorrQnc[1]->Fill(qnc,qnc_vzero);
    if(lCentrality>=40. && lCentrality<60. ) fHistCorrQnc[2]->Fill(qnc,qnc_vzero);
    if(lCentrality>=60. && lCentrality<100. ) fHistCorrQnc[3]->Fill(qnc,qnc_vzero);

  }
  if(fQVZEROA>1 && fQTPCA>1) {
    if(lCentrality>=0. && lCentrality<20.)    fHistCorrQna[0]->Fill(qna,qna_vzero);
    if(lCentrality>=20. && lCentrality<40.) fHistCorrQna[1]->Fill(qna,qna_vzero);
    if(lCentrality>=40. && lCentrality<60. ) fHistCorrQna[2]->Fill(qna,qna_vzero);
    if(lCentrality>=60. && lCentrality<100. ) fHistCorrQna[3]->Fill(qna,qna_vzero);
  }
  if(fQVZEROC>0 && fQVZEROA>0) {
    if(lCentrality>=0. && lCentrality<20.) fHistQAQB_VZERO[0]->Fill(qaqc_vzero);
    if(lCentrality>=20. && lCentrality<40.) fHistQAQB_VZERO[1]->Fill(qaqc_vzero);
    if(lCentrality>=40. && lCentrality<60.) fHistQAQB_VZERO[2]->Fill(qaqc_vzero);
    if(lCentrality>=60. && lCentrality<100.) fHistQAQB_VZERO[3]->Fill(qaqc_vzero);
    fHist_V0AV0C->Fill(qaqc_vzero);
    SP_V0AV0C_default->Fill(lCentrality,qaqc_vzero,1);
  }
  if(fQVZEROC>0 && fQTPC>0){
    fHist_V0CTPC->Fill(qaqc_vzeroctpc);
    SP_V0CTPC_default->Fill(lCentrality,qaqc_vzeroctpc,1);
  }
  if(fQVZEROA>0 && fQTPC>0){
    fHist_V0ATPC->Fill(qaqc_vzeroatpc);
    SP_V0ATPC_default->Fill(lCentrality,qaqc_vzeroatpc,1);
  }


  //Calculate uQ
  Double_t uQ=0;
  Double_t uQ_vzeroa=0;
  Double_t uQ_vzeroc=0;
  for(Int_t i=0;i<selectedTracksAssociated->GetEntriesFast();i++){
    AliAssociatedTrackYSMC* trackPOI=(AliAssociatedTrackYSMC*)selectedTracksAssociated->At(i);
    Double_t phi=trackPOI->Phi();
    Double_t eta=trackPOI->Eta();
    Double_t pt=trackPOI->Pt();
    Int_t SpAssoc=trackPOI->WhichCandidate();
    Double_t cosn = TMath::Cos(fHarmonic*phi);
    Double_t sinn = TMath::Sin(fHarmonic*phi);
    //VZERO-TPC-VZERO
    uQ_vzeroa=cosn*vzeroa_qmcos+sinn*vzeroa_qmsin;
    uQ_vzeroc=cosn*vzeroc_qmcos+sinn*vzeroc_qmsin;
    Double_t w=1;
    if(fQVZEROA>0){
      SP_uVZEROA_PP[SpAssoc]->Fill(pt,uQ_vzeroa,1);
      if(lCentrality>=0. && lCentrality<20.) SP_uVZEROA[SpAssoc]->Fill(pt,uQ_vzeroa,w);
      if(lCentrality>=20. && lCentrality<40.) SP_uVZEROA1[SpAssoc]->Fill(pt,uQ_vzeroa,w);
      if(lCentrality>=40. && lCentrality<60.) SP_uVZEROA2[SpAssoc]->Fill(pt,uQ_vzeroa,w);
      if(lCentrality>=60. && lCentrality<100.) SP_uVZEROA3[SpAssoc]->Fill(pt,uQ_vzeroa,w);
    }
    if(fQVZEROC>0){
      SP_uVZEROC_PP[SpAssoc]->Fill(pt,uQ_vzeroa,1);
      if(lCentrality>=0. && lCentrality<20.) SP_uVZEROC[SpAssoc]->Fill(pt,uQ_vzeroc,w);
      if(lCentrality>=20. && lCentrality<40.) SP_uVZEROC1[SpAssoc]->Fill(pt,uQ_vzeroc,w);
      if(lCentrality>=40. && lCentrality<60.) SP_uVZEROC2[SpAssoc]->Fill(pt,uQ_vzeroc,w);
      if(lCentrality>=60. && lCentrality<100.) SP_uVZEROC3[SpAssoc]->Fill(pt,uQ_vzeroc,w);
    }
    //TPC-TPC
    if(fQTPCC<1 || fQTPCA<1) continue;
    if( eta<0.4 && eta>-0.4 ) continue;
    if(eta<0){
      uQ=(cosn*tpca_qmcos+sinn*tpca_qmsin);  //u x Q/M_a
      SP_uTPCA->Fill(pt,uQ,1);
    }else{
      uQ=(cosn*tpcc_qmcos+sinn*tpcc_qmsin);
      SP_uTPCC->Fill(pt,uQ,1);
    }
    SP_uTPC_PP[SpAssoc]->Fill(pt,uQ,1);
    if(lCentrality>=0. && lCentrality<20.) SP_uTPC[SpAssoc]->Fill(pt,uQ,1);
    if(lCentrality>=20. && lCentrality<40.) SP_uTPC1[SpAssoc]->Fill(pt,uQ,1);
    if(lCentrality>=40. && lCentrality<60.) SP_uTPC2[SpAssoc]->Fill(pt,uQ,1);
    if(lCentrality>=60. && lCentrality<100.) SP_uTPC3[SpAssoc]->Fill(pt,uQ,1);
  }

*/




}


TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedTracksLeading(AliAODEvent *fAOD,Bool_t leading,TObjArray*tracks) {
  //TObjArray *tracks = new TObjArray;
  //tracks->SetOwner(kTRUE);
  Int_t nTracks = fAOD->GetNumberOfTracks();
  Double_t pidqa[5];
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (!aodTrack)      continue;
    if (!IsAcceptedTrack(aodTrack))      continue;
    if (aodTrack->Charge() == 0)      continue;
    if (aodTrack->Pt() > fPtMax) continue;
    
    if(leading){
     pidqa[0]=aodTrack->Pt();
     pidqa[1]=aodTrack->Eta();
     pidqa[2]=aodTrack->Phi();
     pidqa[3]=lCentrality;
     pidqa[4]=fPrimaryZVtx;
     if(!fDataType){
       Int_t myTrackLabel = TMath::Abs(aodTrack->GetLabel());
       TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
       AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(myTrackLabel);
       if(fQA) if(mcTrack->IsPhysicalPrimary()) fHistLeadQA->Fill(pidqa,0);
     }else {
        fHistLeadQA->Fill(pidqa,0);
     }
    }

    Int_t SpAsso=0;
    tracks->Add(new AliAssociatedTrackYSMC(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, SpAsso, 1));
  }
  return tracks;
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedTracksPID(AliAODEvent *fAOD) {
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  Int_t nTracks = fAOD->GetNumberOfTracks();
  Double_t pidqa[4];
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
    Int_t SpPID=-999;
    if (!aodTrack)  continue;
    if (!IsAcceptedTrack(aodTrack))    continue;
    if (aodTrack->Charge() == 0)      continue;
    Double_t nSigmaKaonTPC = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kKaon);
    Double_t nSigmaPionTPC = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kPion);
    Double_t nSigmaProtonTPC = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kProton);
    Double_t nSigmaKaonTOF = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kKaon);
    Double_t nSigmaPionTOF = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kPion);
    Double_t nSigmaProtonTOF = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kProton);
    fHistNsig[0]->Fill(aodTrack->Pt(),nSigmaPionTPC);
    fHistNsig[1]->Fill(aodTrack->Pt(),nSigmaKaonTPC);
    fHistNsig[2]->Fill(aodTrack->Pt(),nSigmaProtonTPC);
    fHistNsig[3]->Fill(aodTrack->Pt(),nSigmaPionTOF);
    fHistNsig[4]->Fill(aodTrack->Pt(),nSigmaKaonTOF);
    fHistNsig[5]->Fill(aodTrack->Pt(),nSigmaProtonTOF);
    fHistNsigcorr[3]->Fill(nSigmaPionTPC,nSigmaPionTOF);
    fHistNsigcorr[4]->Fill(nSigmaKaonTPC,nSigmaKaonTOF);
    fHistNsigcorr[5]->Fill(nSigmaProtonTPC,nSigmaProtonTOF);

    Double_t d2nsigmakaon =  nSigmaKaonTPC * nSigmaKaonTPC + nSigmaKaonTOF * nSigmaKaonTOF;
    Double_t d2nsigmapion =  nSigmaPionTPC * nSigmaPionTPC + nSigmaPionTOF * nSigmaPionTOF;
    Double_t d2nsigmaproton =  nSigmaProtonTPC * nSigmaProtonTPC + nSigmaProtonTOF * nSigmaProtonTOF;

    Bool_t fPIDTOF = kTRUE;
    if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fAOD->GetTrack(i)) == 0)      fPIDTOF = kFALSE;
    else      fPIDTOF = kTRUE;

    Double_t nSigmaKaonTOFTPC;
    Double_t nSigmaPionTOFTPC;
    Double_t nSigmaProtonTOFTPC;

    if (fPIDTOF && aodTrack->Pt() > 0.5) {
      nSigmaKaonTOFTPC = TMath::Sqrt(d2nsigmakaon);
      nSigmaPionTOFTPC = TMath::Sqrt(d2nsigmapion);
      nSigmaProtonTOFTPC = TMath::Sqrt(d2nsigmaproton);
    } else {
      nSigmaKaonTOFTPC = TMath::Abs(nSigmaKaonTPC);
      nSigmaPionTOFTPC = TMath::Abs(nSigmaPionTPC);
      nSigmaProtonTOFTPC = TMath::Abs(nSigmaProtonTPC);
    }
    if ((nSigmaKaonTOFTPC < fMaxnSigmaTPCTOF) &&  (nSigmaKaonTOFTPC < nSigmaPionTOFTPC) &&   (nSigmaKaonTOFTPC < nSigmaProtonTOFTPC)) SpPID = 0;
    if ((nSigmaPionTOFTPC < fMaxnSigmaTPCTOF) &&  (nSigmaPionTOFTPC < nSigmaKaonTOFTPC) &&   (nSigmaPionTOFTPC < nSigmaProtonTOFTPC)) SpPID = 1;
    if ((nSigmaProtonTOFTPC < fMaxnSigmaTPCTOF) &&  (nSigmaProtonTOFTPC < nSigmaKaonTOFTPC) &&   (nSigmaProtonTOFTPC < nSigmaPionTOFTPC)) SpPID = 2;

    pidqa[0]=aodTrack->Pt();
    pidqa[1]=aodTrack->Eta();
    pidqa[2]=RangePhi(aodTrack->Phi());
    pidqa[3]=lCentrality;
    if(SpPID<0) continue;
    if(fQA) fHistPIDQA->Fill(pidqa, SpPID);
    if(SpPID==1) fHistNsigcorr[0]->Fill(nSigmaPionTPC,nSigmaPionTOF);
    if(SpPID==0) fHistNsigcorr[1]->Fill(nSigmaKaonTPC,nSigmaKaonTOF);
    if(SpPID==2) fHistNsigcorr[2]->Fill(nSigmaProtonTPC,nSigmaProtonTOF);

    tracks->Add(new AliAssociatedTrackYSMC(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, SpPID, 1));
  }
  return tracks;
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedTracksAssociated(AliAODEvent *fAODEvent) {
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  TObjArray *dtrack = new TObjArray;
  dtrack->SetOwner(kTRUE);

  Int_t nTracks = fAODEvent->GetNumberOfTracks();
  for (Int_t i = 0; i < nTracks - 1; i++) {
    AliAODTrack *aodTrack1 = dynamic_cast<AliAODTrack *>(fAODEvent->GetTrack(i));
    if (!aodTrack1)  continue;
    Double_t nSigmaKaonTPC_phi1 =        fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kKaon);
    Double_t nSigmaKaonTOF_phi1 =        fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kKaon);
    Double_t nSigmaPionTPC_phi1 =        fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kPion);
    Double_t nSigmaPionTOF_phi1 =        fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kPion);
    Double_t nSigmaProtonTPC_phi1 =        fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kProton);
    Double_t nSigmaProtonTOF_phi1 =        fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kProton);
    Double_t d2sigmaphi1kaontpctof = nSigmaKaonTPC_phi1 * nSigmaKaonTPC_phi1 + nSigmaKaonTOF_phi1 * nSigmaKaonTOF_phi1;
    Double_t d2sigmaphi1piontpctof = nSigmaPionTPC_phi1 * nSigmaPionTPC_phi1 + nSigmaPionTOF_phi1 * nSigmaPionTOF_phi1;
    Double_t d2sigmaphi1protontpctof = nSigmaProtonTPC_phi1 * nSigmaProtonTPC_phi1 + nSigmaProtonTOF_phi1 * nSigmaProtonTOF_phi1;
    Bool_t fPIDTOF_phi1;

    if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, aodTrack1) == 0)
      fPIDTOF_phi1 = kFALSE;
    else
      fPIDTOF_phi1 = kTRUE;

    Double_t nSigmaKaonTOFTPC_phi1;
    Double_t nSigmaPionTOFTPC_phi1;
    Double_t nSigmaProtonTOFTPC_phi1;
    if (fPIDTOF_phi1 && aodTrack1->Pt() > 0.5) {
      nSigmaKaonTOFTPC_phi1 = TMath::Sqrt(d2sigmaphi1kaontpctof);
      nSigmaPionTOFTPC_phi1 = TMath::Sqrt(d2sigmaphi1piontpctof);
      nSigmaProtonTOFTPC_phi1 = TMath::Sqrt(d2sigmaphi1protontpctof);
    } else {
      nSigmaKaonTOFTPC_phi1 = TMath::Abs(nSigmaKaonTPC_phi1);
      nSigmaPionTOFTPC_phi1 = TMath::Abs(nSigmaPionTPC_phi1);
      nSigmaProtonTOFTPC_phi1 = TMath::Abs(nSigmaProtonTPC_phi1);
    }
    if (!IsAcceptedPhiDaughterTrack(aodTrack1))
      continue;

    fHistPhiDTPCNSig->Fill(aodTrack1->Pt(), nSigmaKaonTPC_phi1);
    fHistPhiDTOFNSig->Fill(aodTrack1->Pt(), nSigmaKaonTOF_phi1);
    fHistPhiDTPCTOFNSig->Fill(aodTrack1->Pt(),
    TMath::Sqrt(d2sigmaphi1kaontpctof));

    Bool_t isKaon1 = kFALSE;
    // if(nSigmaKaonTOFTPC_phi1<5.0 &&
    // nSigmaKaonTOFTPC_phi1<nSigmaPionTOFTPC_phi1 &&
    // nSigmaKaonTOFTPC_phi1<nSigmaProtonTOFTPC_phi1) isKaon1=kTRUE;
    if (TMath::Abs(nSigmaKaonTPC_phi1) < 5.0 &&
    TMath::Abs(nSigmaKaonTPC_phi1) < TMath::Abs(nSigmaPionTPC_phi1) &&
    TMath::Abs(nSigmaKaonTPC_phi1) < TMath::Abs(nSigmaProtonTPC_phi1))
    isKaon1 = kTRUE;
    if (!isKaon1)
    continue;

    dtrack->Add(new AliMixTrackYSMC(
      aodTrack1->Charge(), aodTrack1->Eta(), aodTrack1->Phi(),
      aodTrack1->Pt(), aodTrack1->Px(), aodTrack1->Py(), aodTrack1->Pz()));

      for (Int_t j = i + 1; j < nTracks; j++) {
        AliAODTrack *aodTrack2 =
        dynamic_cast<AliAODTrack *>(fAODEvent->GetTrack(j));
        if (!aodTrack2)
        continue;
        Double_t nSigmaKaonTPC_phi2 =
        fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kKaon);
        Double_t nSigmaKaonTOF_phi2 =
        fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kKaon);
        Double_t nSigmaPionTPC_phi2 =
        fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kPion);
        Double_t nSigmaPionTOF_phi2 =
        fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kPion);
        Double_t nSigmaProtonTPC_phi2 =
        fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kProton);
        Double_t nSigmaProtonTOF_phi2 =
        fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kProton);
        Double_t d2sigmaphi2kaontpctof = nSigmaKaonTPC_phi2 * nSigmaKaonTPC_phi2 +
        nSigmaKaonTOF_phi2 * nSigmaKaonTOF_phi2;
        Double_t d2sigmaphi2piontpctof = nSigmaPionTPC_phi2 * nSigmaPionTPC_phi2 +
        nSigmaPionTOF_phi2 * nSigmaPionTOF_phi2;
        Double_t d2sigmaphi2protontpctof =
        nSigmaProtonTPC_phi2 * nSigmaProtonTPC_phi2 +
        nSigmaProtonTOF_phi2 * nSigmaProtonTOF_phi2;
        Bool_t fPIDTOF_phi2;

        if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, aodTrack2) == 0)
        fPIDTOF_phi2 = kFALSE;
        else
        fPIDTOF_phi2 = kTRUE;

        Double_t nSigmaKaonTOFTPC_phi2;
        Double_t nSigmaPionTOFTPC_phi2;
        Double_t nSigmaProtonTOFTPC_phi2;
        if (fPIDTOF_phi2 && aodTrack2->Pt() > 0.5) {
          nSigmaKaonTOFTPC_phi2 = TMath::Sqrt(d2sigmaphi2kaontpctof);
          nSigmaPionTOFTPC_phi2 = TMath::Sqrt(d2sigmaphi2piontpctof);
          nSigmaProtonTOFTPC_phi2 = TMath::Sqrt(d2sigmaphi2protontpctof);
        } else {
          nSigmaKaonTOFTPC_phi2 = TMath::Abs(nSigmaKaonTPC_phi2);
          nSigmaPionTOFTPC_phi2 = TMath::Abs(nSigmaPionTPC_phi2);
          nSigmaProtonTOFTPC_phi2 = TMath::Abs(nSigmaProtonTPC_phi2);
        }

      Bool_t isKaon2 = kFALSE;
      // if(nSigmaKaonTOFTPC_phi2<5.0 &&
      // nSigmaKaonTOFTPC_phi2<nSigmaPionTOFTPC_phi2 &&
      // nSigmaKaonTOFTPC_phi2<nSigmaProtonTOFTPC_phi2) isKaon2=kTRUE;
      if (TMath::Abs(nSigmaKaonTPC_phi2) < 5.0 &&
          TMath::Abs(nSigmaKaonTPC_phi2) < TMath::Abs(nSigmaPionTPC_phi2) &&
          TMath::Abs(nSigmaKaonTPC_phi2) < TMath::Abs(nSigmaProtonTPC_phi2))
        isKaon2 = kTRUE;
      if (!isKaon2)        continue;

      if (!IsAcceptedPhiDaughterTrack(aodTrack2))        continue;
      if (aodTrack1->GetID() == aodTrack2->GetID())        continue;
      if (aodTrack1->Charge() == aodTrack2->Charge())        continue;

      Double_t mass_phi;
      Double_t px_phi = aodTrack1->Px() + aodTrack2->Px();
      Double_t py_phi = aodTrack1->Py() + aodTrack2->Py();
      Double_t pz_phi = aodTrack1->Pz() + aodTrack2->Pz();
      Double_t p_phi =
          TMath::Sqrt(px_phi * px_phi + py_phi * py_phi + pz_phi * pz_phi);
      Double_t phi_phi = atan2(py_phi, px_phi);
      if (phi_phi < 0.)
        phi_phi += 2 * TMath::Pi();
      Double_t px1 = aodTrack1->Px();
      Double_t py1 = aodTrack1->Py();
      Double_t pz1 = aodTrack1->Pz();
      Double_t E1 =
          TMath::Sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + 0.493677 * 0.493677);
      Double_t px2 = aodTrack2->Px();
      Double_t py2 = aodTrack2->Py();
      Double_t pz2 = aodTrack2->Pz();
      Double_t E2 =
          TMath::Sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + 0.493677 * 0.493677);
      mass_phi =
          TMath::Sqrt((E1 + E2) * (E1 + E2) - (px1 + px2) * (px1 + px2) -
                      (py1 + py2) * (py1 + py2) - (pz1 + pz2) * (pz1 + pz2));
      Double_t pt_phi =
          TMath::Sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2));

      // if(pt_phi<0.5) continue;

      Double_t eta_phi = 0.5 * log((p_phi + pz_phi) / (p_phi - pz_phi));

      Double_t PhiQA[3] = {mass_phi, pt_phi, lCentrality};
      fHistMass_PhiMeson->Fill(PhiQA);
      if (TMath::Abs(eta_phi) > 0.8)
        continue;

      Int_t SpPhi = -999.;

      const Double_t fPhiMass =
          TDatabasePDG::Instance()->GetParticle(333)->Mass();
      if (TMath::Abs(mass_phi - fPhiMass) < 0.006)
        SpPhi = 0; // same for step
      if (mass_phi > 0.99 && mass_phi < fPhiMass - 0.006)
        SpPhi = 1; // same as step
      if (mass_phi > fPhiMass + 0.006 && mass_phi < fPhiMass + 0.018)
        SpPhi = 2; // same as step
      if (mass_phi > fPhiMass + 0.018 && mass_phi < fPhiMass + 0.030)
        SpPhi = 3; // same as step
      if (mass_phi > fPhiMass + 0.030 && mass_phi < fPhiMass + 0.042)
        SpPhi = 4; // same as step
      if (mass_phi > fPhiMass + 0.042 && mass_phi < fPhiMass + 0.054)
        SpPhi = 5; // same as step
      if (mass_phi > fPhiMass + 0.054 && mass_phi < fPhiMass + 0.066)
        SpPhi = 6; // same as step

      if(SpPhi<0) continue;
      tracks->Add(new AliAssociatedTrackYSMC(0, eta_phi, phi_phi, pt_phi, -999,
                                           aodTrack1->GetID(),
                                           aodTrack2->GetID(), SpPhi, 1));

      Double_t spPhiQA[4] = {eta_phi, phi_phi, (Float_t)SpPhi, lCentrality};
      fHist_PhiQA->Fill(spPhiQA);
    }
  }

  // Mixed Event
  Double_t pvxMix = fPrimaryZVtx;
  Double_t counterMix = 0;
  AliEventPool *pool = fPoolMgr1->GetEventPool(lCentrality, pvxMix);
  if (!pool)
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCentrality,pvxMix));
  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||    pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix = 0; jMix < nMix; jMix++) {
      TObjArray *mixEvents = pool->GetEvent(jMix);
      for (Int_t i = 0; i < dtrack->GetEntriesFast(); i++) {
        AliMixTrackYSMC *dtrack1 = (AliMixTrackYSMC *)dtrack->At(i);
        if (!dtrack1)
          continue;
        Double_t pdx1 = dtrack1->Px();
        Double_t pdy1 = dtrack1->Py();
        Double_t pdz1 = dtrack1->Pz();
        Double_t Ed1 = TMath::Sqrt(pdx1 * pdx1 + pdy1 * pdy1 + pdz1 * pdz1 +
                                   0.493677 * 0.493677);
        counterMix++;
        for (Int_t j = 0; j < mixEvents->GetEntriesFast(); j++) {
          AliMixTrackYSMC *dtrack2 = (AliMixTrackYSMC *)mixEvents->At(j);
          if (!dtrack2)
            continue;
          if (dtrack1->Charge() == dtrack2->Charge())
            continue;
          Double_t pdx2 = dtrack2->Px();
          Double_t pdy2 = dtrack2->Py();
          Double_t pdz2 = dtrack2->Pz();
          Double_t Ed2 = TMath::Sqrt(pdx2 * pdx2 + pdy2 * pdy2 + pdz2 * pdz2 +
                                     0.493677 * 0.493677);
          Double_t mass_phi_mix = TMath::Sqrt(
              (Ed1 + Ed2) * (Ed1 + Ed2) - (pdx1 + pdx2) * (pdx1 + pdx2) -
              (pdy1 + pdy2) * (pdy1 + pdy2) - (pdz1 + pdz2) * (pdz1 + pdz2));
          Double_t pt_phi_mix = TMath::Sqrt((pdx1 + pdx2) * (pdx1 + pdx2) +
                                            (pdy1 + pdy2) * (pdy1 + pdy2));
          // if(pt_phi_mix<0.5) continue;
          Double_t px_phi_mix = pdx1 + pdx2;
          Double_t py_phi_mix = pdy1 + pdy2;
          Double_t pz_phi_mix = pdz1 + pdz2;
          Double_t p_phi_mix =
              TMath::Sqrt(px_phi_mix * px_phi_mix + py_phi_mix * py_phi_mix +
                          pz_phi_mix * pz_phi_mix);
          Double_t phi_phi_mix = atan2(py_phi_mix, px_phi_mix);
          if (phi_phi_mix < 0.)
            phi_phi_mix += 2 * TMath::Pi();
          Double_t eta_phi_mix =
              0.5 * log((p_phi_mix + pz_phi_mix) / (p_phi_mix - pz_phi_mix));
          if (TMath::Abs(eta_phi_mix) > 0.8)
            continue;
          Double_t PhiQAMix[3] = {mass_phi_mix, pt_phi_mix, lCentrality};
          fHistMass_PhiMeson_MIX->Fill(PhiQAMix, 1. / (Double_t)nMix);
        }
      }
    }
  }

  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i = 0; i < dtrack->GetEntriesFast(); i++) {
    AliMixTrackYSMC *particle = (AliMixTrackYSMC *)dtrack->At(i);
    tracksClone->Add(new AliMixTrackYSMC(
        particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
        particle->Px(), particle->Py(), particle->Pz()));
  }
  
  pool->UpdatePool(tracksClone);

  return tracks;
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedV0Tracks(const AliAODEvent *fAODEvent) {
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  // V0Particles
  Int_t nv0s = fAODEvent->GetNumberOfV0s();
  for (Int_t iv0 = 0; iv0 < nv0s; iv0++) {
    AliAODv0 *aodv0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(iv0));
    if (!aodv0) {
      AliError(Form("ERROR: Could not retrieve aodv0 %d", iv0));
      continue;
    }
    fHist_V0Stat->Fill(0);
    if(!IsAcceptedV0(aodv0)) continue;
    fHist_V0Stat->Fill(7);
    // c*taup;ppp;
    const Double_t kLambdaMass =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    const Double_t kK0Mass = TDatabasePDG::Instance()->GetParticle(311)->Mass();
    // Double_t cutcTauLam  = fcutctau*7.89;//c*mean life time=2.632*2.9979 cm
    // Double_t cutcTauK0   = fcutctau*2.68; //c*mean life time=0.8954*2.9979 cm
    // life time < 20(K0s), 30(lambda)
    fcutcTauLam=3.8*7.89;
    fcutcTauK0=5*2.68;
    Bool_t cutK0ctau = IsAcceptedDecayLength(aodv0,kK0Mass,fcutcTauK0);
    Bool_t cutLambdactau =IsAcceptedDecayLength(aodv0,kLambdaMass,fcutcTauLam);
    Bool_t cutAntiLambdactau = IsAcceptedDecayLength(aodv0,kLambdaMass,fcutcTauLam);

    // cpa
    Double_t cpa = aodv0->CosPointingAngle(fAODEvent->GetPrimaryVertex());
    fcosMinK0s=0.97;
    fcosMinLambda=0.99;
    Bool_t cpaK0s = cpa > fcosMinK0s;
    Bool_t cpaLambda = cpa > fcosMinLambda;

    const AliAODTrack* myTrackPos=0;
    const AliAODTrack* myTrackNeg=0;



    
    
    AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack *>(aodv0->GetDaughter(0)); // The first dauther track, which should be positive
    AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack *>(aodv0->GetDaughter(1)); // The second dauther track, which should be negative

    if (!myTrackPosTest || !myTrackNegTest) {
      Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
      continue;
    }

    if (!IsAcceptedDaughterTrack(myTrackPosTest) || !IsAcceptedDaughterTrack(myTrackNegTest))   continue;

    Double_t alpha=aodv0->AlphaV0();
    Double_t qt = aodv0->PtArmV0();
    fHist_V0Stat->Fill(8);
    if ((myTrackPosTest->Charge() == 1) && (myTrackNegTest->Charge() == -1)) {
      hv0dcharge->Fill(0);
      myTrackPos = myTrackPosTest;
      myTrackNeg = myTrackNegTest;
    }

    if ((myTrackPosTest->Charge() == -1) && (myTrackNegTest->Charge() == 1)) {
      hv0dcharge->Fill(1);
      myTrackPos = myTrackNegTest;
      myTrackNeg =  myTrackPosTest;
    }


    if (myTrackPosTest->Charge() == myTrackNegTest->Charge()) {
      hv0dcharge->Fill(2);
      continue;
    }

    fHist_AP[0]->Fill(alpha,qt);
    fHist_V0Stat->Fill(9);
    // PID Cut
    // Positive tracks
    Double_t nSigmaPosPionTPC   =  fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
    Double_t nSigmaPosPionTOF   =  fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kPion);
    Double_t nSigmaPosKaonTPC   =  fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kKaon);
    Double_t nSigmaPosKaonTOF   =  fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kKaon);
    Double_t nSigmaPosProtonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
    Double_t nSigmaPosProtonTOF =  fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kProton);
    // negative tracks
    Double_t nSigmaNegPionTPC   =  fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion);
    Double_t nSigmaNegPionTOF   =  fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kPion);
    Double_t nSigmaNegKaonTPC   =  fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kKaon);
    Double_t nSigmaNegKaonTOF   =  fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kKaon);
    Double_t nSigmaNegProtonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton);
    Double_t nSigmaNegProtonTOF =  fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kProton);

    Bool_t bpPion = TMath::Abs(nSigmaPosPionTPC) <= fMaxnSigmaTPCV0;
    Bool_t bpProton = TMath::Abs(nSigmaPosProtonTPC) <= fMaxnSigmaTPCV0;
    Bool_t bnPion = TMath::Abs(nSigmaNegPionTPC) <= fMaxnSigmaTPCV0;
    Bool_t bnProton = TMath::Abs(nSigmaNegProtonTPC) <= fMaxnSigmaTPCV0;

    Bool_t bpPion_tof = TMath::Abs(nSigmaPosPionTPC) <= 4.;
    Bool_t bpProton_tof = TMath::Abs(nSigmaPosProtonTPC) <=4;
    Bool_t bnPion_tof = TMath::Abs(nSigmaNegPionTOF) <= 4;
    Bool_t bnProton_tof = TMath::Abs(nSigmaNegProtonTOF) <=4;

    Bool_t fTOFV0=kTRUE;

    if(fTOFV0 && myTrackPos->Pt()>0.3 && fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackPos) != 0){
      bpPion=bpPion && bpPion_tof;
      bpProton=bpProton && bpProton_tof;
    }

    if(fTOFV0 && myTrackNeg->Pt()>0.3 && fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackNeg) != 0){
      bnPion=bnPion && bnPion_tof;
      bnProton=bnProton && bnProton_tof;
    }

    // Arme ntros podoranski cutOD
    Bool_t k0APcut = (aodv0->PtArmV0() > (TMath::Abs(0.2 * aodv0->AlphaV0())));
    Bool_t kGammaconvcut = !(TMath::Power(aodv0->AlphaV0() / 0.95, 2) + TMath::Power(aodv0->PtArmV0() / 0.05, 2) < 1);

    if (k0APcut) fHist_AP[4]->Fill(aodv0->AlphaV0(), aodv0->PtArmV0());
    if (kGammaconvcut)  fHist_AP[5]->Fill(aodv0->AlphaV0(), aodv0->PtArmV0());

    Bool_t fPIDV0=kFALSE;

    //if(fPIDV0 && fTOFV0) if ((myTrackPos->Pt()>0.3 && fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackPos) == 0)  || (myTrackNeg->Pt()>0.3 && fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackNeg) == 0)) continue;

    fHistPosNsig[0]->Fill(myTrackPos->Pt(), nSigmaPosPionTPC);
    fHistPosNsig[1]->Fill(myTrackPos->Pt(), nSigmaPosProtonTPC);
    fHistPosNsig[2]->Fill(myTrackPos->Pt(), nSigmaPosKaonTPC);
    fHistPosNsig[3]->Fill(myTrackPos->Pt(), nSigmaPosPionTOF);
    fHistPosNsig[4]->Fill(myTrackPos->Pt(), nSigmaPosProtonTOF);
    fHistPosNsig[5]->Fill(myTrackPos->Pt(), nSigmaPosKaonTOF);

    fHistNegNsig[0]->Fill(myTrackNeg->Pt(), nSigmaNegPionTPC);
    fHistNegNsig[1]->Fill(myTrackNeg->Pt(), nSigmaNegProtonTPC);
    fHistNegNsig[2]->Fill(myTrackNeg->Pt(), nSigmaNegKaonTPC);
    fHistNegNsig[3]->Fill(myTrackNeg->Pt(), nSigmaNegPionTOF);
    fHistNegNsig[4]->Fill(myTrackNeg->Pt(), nSigmaNegProtonTOF);
    fHistNegNsig[5]->Fill(myTrackNeg->Pt(), nSigmaNegKaonTOF);

    /*
    if (fQA) {
    if (!kGammaconvcut)   fHistPosNsigQA[0]->Fill(myTrackPos->Pt(),nSigmaPosPionTPC); // Nsigma of Pion if the AP diagram indicate gamma conversion
    }
    */
    Bool_t cutK0Pid = (bpPion && bnPion);
    Bool_t cutLambdaPid = (bpProton && bnPion);
    Bool_t cutAntiLambdaPid = (bpPion && bnProton);

    // Mass cut
    Double_t lInvMassLambda = aodv0->MassLambda();
    Double_t lInvMassK0 = aodv0->MassK0Short();
    Double_t lInvMassAntiLambda = aodv0->MassAntiLambda();

    const Double_t kProtonMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    const Double_t kPionMass =  TDatabasePDG::Instance()->GetParticle(211)->Mass();

    Double_t mispidmass= -999.;
    if (fQA) {
      // QA to evaluate bump
      Double_t px1 = myTrackPos->Px();
      Double_t py1 = myTrackPos->Py();
      Double_t pz1 = myTrackPos->Pz();
      Double_t px2 = myTrackNeg->Px();
      Double_t py2 = myTrackNeg->Py();
      Double_t pz2 = myTrackNeg->Pz();
      Double_t px_v0 = px1 + px2;
      Double_t py_v0 = py1 + py2;
      Double_t pz_v0 = pz1 + pz2;
      Double_t p2_v0 = px_v0 * px_v0 + py_v0 * py_v0 + pz_v0 * pz_v0;
      Double_t E_pos_proton = TMath::Sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + kProtonMass * kProtonMass);
      Double_t E_neg_proton = TMath::Sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + kProtonMass * kProtonMass);
      Double_t E_neg_pion = TMath::Sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + kPionMass * kPionMass);
      Double_t E_pos_pion = TMath::Sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + kPionMass * kPionMass);
      mispidmass = TMath::Sqrt((E_pos_pion + E_neg_pion) * (E_pos_pion + E_neg_pion) - p2_v0);
    }

    Bool_t mispid=TMath::Abs(mispidmass-kK0Mass)<0.01;

    Bool_t cutMassLambda = ((lInvMassLambda > 1.05) && (lInvMassLambda < 1.25));
    Bool_t cutMassAntiLambda = ((lInvMassAntiLambda > 1.05) && (lInvMassAntiLambda < 1.25));
    Bool_t cutMassK0 = (lInvMassK0 > 0.4) && (lInvMassK0 < 0.6);

    if(cutK0Pid){
      fHist_V0Stat->Fill(10);
      if(cutK0ctau){
        fHist_V0Stat->Fill(11);
        if(k0APcut&&kGammaconvcut)  fHist_V0Stat->Fill(12);
      }
    }
    if(cutLambdaPid){
      fHist_V0Stat->Fill(13);
      if(cutLambdactau){
        fHist_V0Stat->Fill(14);
        if(!k0APcut&&kGammaconvcut) fHist_V0Stat->Fill(15);
      }
    }
    /*
    Bool_t IDK0s = cutK0ctau && cpaK0s && k0APcut && (!cutMassLambda) && (!cutMassAntiLambda);//&& kGammaconvcut;
    Bool_t IDLa =  cutLambdactau && cpaLambda && !k0APcut && (!cutMassK0); //&& kGammaconvcut;
    Bool_t IDALa =  cutAntiLambdactau && cpaLambda && !k0APcut && (!cutMassK0);// && kGammaconvcut;
    */
    /*
    if(fPIDV0){
      IDK0s=IDK0s && cutK0Pid ;
      IDLa= IDLa && cutLambdaPid && !mispid;
      IDALa= IDALa && cutAntiLambdaPid && !mispid;
    }
    */
    Bool_t IDK0s = cutK0ctau && cpaK0s && k0APcut;
    Bool_t IDLa =  cutLambdactau && cpaLambda && !k0APcut;
    Bool_t IDALa =  cutAntiLambdactau && cpaLambda && !k0APcut;


    if (IDK0s)    fHist_AP[1]->Fill(aodv0->AlphaV0(), aodv0->PtArmV0());
    if (IDLa  || IDALa)    fHist_AP[2]->Fill(aodv0->AlphaV0(), aodv0->PtArmV0());

    Double_t K0sMassQA[4] = {lInvMassK0, aodv0->Pt(), lCentrality,aodv0->Eta()};
      if (IDK0s)      fHistMass_K0s->Fill(K0sMassQA);
      Double_t LambdaMassQA[4] = {lInvMassLambda, aodv0->Pt(), lCentrality,aodv0->Eta()};
      if (IDLa)       fHistMass_Lambda->Fill(LambdaMassQA);
      Double_t ALambdaMassQA[4] = {lInvMassAntiLambda, aodv0->Pt(), lCentrality,aodv0->Eta()};
      if (IDALa)      fHistMass_ALambda->Fill(ALambdaMassQA);

      if (cutLambdaPid && cutLambdactau && cpaLambda && !k0APcut && kGammaconvcut && mispidmass > 0) fHistMass_bumpcorr->Fill(lInvMassLambda, mispidmass);
      Int_t SpV0 = -999;
      if ((IDK0s && (TMath::Abs(lInvMassK0 - kK0Mass) < 0.01)))
      SpV0 = 0; // same for step
      if (IDK0s && (lInvMassK0 > 0.40) && (lInvMassK0 < 0.44))
      SpV0 = 1; // same as step
      if (IDK0s && (lInvMassK0 > 0.44) && (lInvMassK0 < kK0Mass - 0.01))
      SpV0 = 2; // same as step
      if (IDK0s && (lInvMassK0 > kK0Mass + 0.01) && (lInvMassK0 < 0.56))
      SpV0 = 3; // same as step
      if (IDK0s && (lInvMassK0 > 0.56) && (lInvMassK0 < 0.60))
      SpV0 = 4; // same as step
      if ((IDLa && TMath::Abs(lInvMassLambda - kLambdaMass) < 0.005) || (IDALa && TMath::Abs(lInvMassAntiLambda - kLambdaMass) < 0.005))
      SpV0 = 5; // same as step
      if ((IDLa  && (lInvMassLambda > kLambdaMass + 0.005) && (lInvMassLambda < 1.25)) || (IDALa && (lInvMassAntiLambda > kLambdaMass + 0.005) && (lInvMassAntiLambda < 1.25)))
      SpV0 = 6; // same as step
      if (SpV0 < 0)          continue;
        Double_t spV0QA[4] = {aodv0->Eta(), aodv0->Phi(), (Float_t)SpV0,  lCentrality};
        fHist_V0QA->Fill(spV0QA);
        if (fQA) {
          if (IDK0s){
            fHistPosNsigQA[0]->Fill(myTrackPos->Pt(), nSigmaPosPionTPC);
            fHistPosNsigQA[1]->Fill(myTrackNeg->Pt(), nSigmaNegPionTPC);
	    fh3NegNsig[0]->Fill(myTrackNeg->Pt(),nSigmaNegPionTPC,nSigmaNegPionTOF);
	    fh3PosNsig[0]->Fill(myTrackPos->Pt(),nSigmaPosPionTPC,nSigmaPosPionTOF);
          }
          if (IDLa){
            fHistPosNsigQA[2]->Fill(myTrackPos->Pt(), nSigmaPosProtonTPC);
	    fHistPosNsigQA[3]->Fill(myTrackNeg->Pt(), nSigmaNegPionTPC);
	    fh3NegNsig[1]->Fill(myTrackNeg->Pt(),nSigmaNegPionTPC,nSigmaNegPionTOF);
	    fh3PosNsig[1]->Fill(myTrackPos->Pt(),nSigmaPosProtonTPC,nSigmaPosProtonTOF);
	  }
          if (IDALa){
            fHistPosNsigQA[4]->Fill(myTrackPos->Pt(), nSigmaPosPionTPC);
            fHistPosNsigQA[5]->Fill(myTrackNeg->Pt(), nSigmaNegProtonTPC);
	    fh3NegNsig[2]->Fill(myTrackNeg->Pt(),nSigmaNegProtonTPC,nSigmaNegProtonTOF);
	    fh3PosNsig[2]->Fill(myTrackPos->Pt(),nSigmaPosPionTPC,nSigmaPosPionTOF);
          }
        }
        tracks->Add(new AliAssociatedTrackYSMC(aodv0->Charge(), aodv0->Eta(), aodv0->Phi(), aodv0->Pt(),aodv0->GetID(), myTrackPos->GetID(), myTrackNeg->GetID(), SpV0, 1));

      if(!fDataType){
        TClonesArray *mcArray1 = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
        if(!mcArray1){
          Printf("No MC particle branch found");
          continue;
        }
        Int_t MotherOfMotherPdg =0;
        Int_t myTrackPosLabel = TMath::Abs(myTrackPos->GetLabel());
        Int_t myTrackNegLabel = TMath::Abs(myTrackNeg->GetLabel());
        AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray1->At(myTrackPosLabel);
        if (!mcPosTrack) continue;         //check if positive daughter track is MC particle
        AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray1->At(myTrackNegLabel);
        if (!mcNegTrack) continue;         //check if negative daughter track is MC particle
        Int_t PosTrackPdg = mcPosTrack->GetPdgCode();
        Int_t NegTrackPdg = mcNegTrack->GetPdgCode();
        Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
        Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();
        if ((myTrackPosMotherLabel<0)||(myTrackNegMotherLabel<0)) continue;    // check if dauter tracks have mother track
        if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;              // require the same mother particle for positive and negative daughter tracks

        AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray1->At(myTrackPosMotherLabel);
        if (!mcPosMother) continue;
        Int_t MotherPdg = mcPosMother->GetPdgCode();
        Int_t MotherOfMother = mcPosMother->GetMother();
        Bool_t IsK0FromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==310)&&(PosTrackPdg==211)&&(NegTrackPdg==-211);            //Moter track is K0 short, Pos is Pion, Neg is Pion
        Bool_t IsLambdaFromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==3122)&&(PosTrackPdg==2212)&&(NegTrackPdg==-211);      //Moter track is Lambda, Pos is Proton, Neg is Pion
        Bool_t IsAntiLambdaFromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==-3122)&&(PosTrackPdg==211)&&(NegTrackPdg==-2212); //Moter track is Anti-Lambda, Pos is pion, Neg is Proton
        if(IsK0FromMC)          fHistMass_K0s_MC->Fill(K0sMassQA);
        if(IsLambdaFromMC)      fHistMass_Lambda_MC->Fill(LambdaMassQA);
        if(IsAntiLambdaFromMC)  fHistMass_ALambda_MC->Fill(ALambdaMassQA);
      }
    }
  return tracks;
}



TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedCascadeTracks(AliAODEvent *fAODEvent) {
  // To select Cascade Particle
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  Int_t nCascades = fAODEvent->GetNumberOfCascades();
  Double_t lInvMassXiMinus=-999.;
  Double_t lInvMassXiPlus=-999.;
  Double_t lInvMassOmegaMinus=-999.;
  Double_t lInvMassOmegaPlus=-999.;
  for (Int_t icasc = 0; icasc < nCascades; icasc++) {
    AliAODcascade *casc = fAODEvent->GetCascade(icasc);
    if (!casc)      continue;
    const AliAODTrack *myTrackCascPos=0;
    const AliAODTrack *myTrackCascNeg=0;
    //    const AliAODTrack*myTrackCascBach;
    AliAODTrack *myTrackCascPosTest = dynamic_cast<AliAODTrack *>(casc->GetDaughter(0));
    AliAODTrack *myTrackCascNegTest = dynamic_cast<AliAODTrack *>(casc->GetDaughter(1));
    const AliAODTrack *myTrackCascBach = dynamic_cast<AliAODTrack *>(casc->GetDecayVertexXi()->GetDaughter(0));
    // myTrackCascBach=myTrackCascBachTest;
    if ((myTrackCascPosTest->Charge() == 1) && (myTrackCascNegTest->Charge() == -1)) {
      myTrackCascPos = myTrackCascPosTest;
      myTrackCascNeg = myTrackCascNegTest;
    }
    if ((myTrackCascPosTest->Charge() == -1) && (myTrackCascNegTest->Charge() == 1)) {
      myTrackCascPos = myTrackCascNegTest;
      myTrackCascNeg = myTrackCascPosTest;
    }
    if (!myTrackCascPos || !myTrackCascNeg || !myTrackCascBach) {
      AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
      continue;
    }
    UInt_t lIdxPosCasc = (UInt_t)TMath::Abs(myTrackCascPos->GetID());
    UInt_t lIdxNegCasc = (UInt_t)TMath::Abs(myTrackCascNeg->GetID());
    UInt_t lCascBachIdx = (UInt_t)TMath::Abs(myTrackCascBach->GetID());
    if (lCascBachIdx == lIdxNegCasc) {
      AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
      continue;
    }
    if (lCascBachIdx == lIdxPosCasc) {
      AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
      continue;
    }
//    if (!IsAcceptedCascade(casc)) continue;
    Double_t lInvMassxi = casc->MassXi();
    Double_t lInvMassomega = casc->MassOmega();
    Short_t lChargeXi = casc->ChargeXi();

    if (!IsAcceptedDaughterTrack(myTrackCascPos) || !IsAcceptedDaughterTrack(myTrackCascNeg) || !IsAcceptedDaughterTrack(myTrackCascBach)) continue;
    Double_t lCasx = casc->DecayVertexXiX();
    Double_t lCasy = casc->DecayVertexXiY();
    Double_t lCasz = casc->DecayVertexXiZ();
    const Double_t kximass =   TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    const Double_t komegamass =  TDatabasePDG::Instance()->GetParticle(3334)->Mass();

    Bool_t topoXi=IsAcceptedCascade(casc);
    Bool_t topoOmega=IsAcceptedCascadeOmega(casc);
    if(!topoXi && !topoOmega) continue;

    Double_t cutctauxi = 6 * 4.91;
    Double_t cutctauomega = 6 * 2.461;
    Double_t lCasDecayLength = TMath::Sqrt(TMath::Power(lCasx - tPrimaryVtxPosition[0], 2) + TMath::Power(lCasy - tPrimaryVtxPosition[1], 2) + TMath::Power(lCasz - tPrimaryVtxPosition[2], 2));
    Double_t lPCas = TMath::Sqrt((casc->MomXiX()) * (casc->MomXiX()) +(casc->MomXiY()) * (casc->MomXiY()) +(casc->MomXiZ()) * (casc->MomXiZ()));
    Double_t lctauXi = (lCasDecayLength * kximass) / lPCas;
    Double_t lctauOmega = (lCasDecayLength * komegamass) / lPCas;
    // Bool_t cutXictau=(lctauXi < cutctauxi);
    // Bool_t cutOmegactau=(lctauOmega<cutctauomega);
    Bool_t cutXictau = kTRUE;
    Bool_t cutOmegactau = kTRUE;



    if (lChargeXi < 0)      lInvMassXiMinus = lInvMassxi;
    if (lChargeXi > 0)      lInvMassXiPlus = lInvMassxi;
    if (lChargeXi < 0)      lInvMassOmegaMinus = lInvMassomega;
    if (lChargeXi > 0)      lInvMassOmegaPlus = lInvMassomega;

    Float_t nSigmaCascBachKaonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascBach, AliPID::kKaon);
    Float_t nSigmaCascBachPionTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascBach, AliPID::kPion);
    Float_t nSigmaCascBachProtonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascBach, AliPID::kProton);
    Float_t nSigmaCascBachKaonTOF =  fPIDResponse->NumberOfSigmasTOF(myTrackCascBach, AliPID::kKaon);
    Float_t nSigmaCascBachPionTOF =  fPIDResponse->NumberOfSigmasTOF(myTrackCascBach, AliPID::kPion);
    Float_t nSigmaCascBachProtonTOF =  fPIDResponse->NumberOfSigmasTOF(myTrackCascBach, AliPID::kProton);

    Float_t nSigmaCascPosKaonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascPos, AliPID::kKaon);
    Float_t nSigmaCascPosPionTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascPos, AliPID::kPion);
    Float_t nSigmaCascPosProtonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascPos, AliPID::kProton);
    Float_t nSigmaCascPosKaonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascPos, AliPID::kKaon);
    Float_t nSigmaCascPosPionTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascPos, AliPID::kPion);
    Float_t nSigmaCascPosProtonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascPos, AliPID::kProton);

    Float_t nSigmaCascNegKaonTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascNeg, AliPID::kKaon);
    Float_t nSigmaCascNegPionTPC =  fPIDResponse->NumberOfSigmasTPC(myTrackCascNeg, AliPID::kPion);
    Float_t nSigmaCascNegProtonTPC = fPIDResponse->NumberOfSigmasTPC(myTrackCascNeg, AliPID::kProton);
    Float_t nSigmaCascNegKaonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascNeg, AliPID::kKaon);
    Float_t nSigmaCascNegPionTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascNeg, AliPID::kPion);
    Float_t nSigmaCascNegProtonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackCascNeg, AliPID::kProton);

    Float_t d2sigmacascbachkaontpctof = nSigmaCascBachKaonTPC * nSigmaCascBachKaonTPC +  nSigmaCascBachKaonTOF * nSigmaCascBachKaonTOF;
    Float_t d2sigmacascbachpiontpctof = nSigmaCascBachPionTPC * nSigmaCascBachPionTPC +  nSigmaCascBachPionTOF * nSigmaCascBachPionTOF;
    Float_t d2sigmacascbachprotontpctof = nSigmaCascBachProtonTPC * nSigmaCascBachProtonTPC + nSigmaCascBachProtonTOF * nSigmaCascBachProtonTOF;

    Float_t d2sigmacascposkaontpctof = nSigmaCascPosKaonTPC * nSigmaCascPosKaonTPC + nSigmaCascPosKaonTOF * nSigmaCascPosKaonTOF;
    Float_t d2sigmacascpospiontpctof = nSigmaCascPosPionTPC * nSigmaCascPosPionTPC + nSigmaCascPosPionTOF * nSigmaCascPosPionTOF;
    Float_t d2sigmacascposprotontpctof = nSigmaCascPosProtonTPC * nSigmaCascPosProtonTPC + nSigmaCascPosProtonTOF * nSigmaCascPosProtonTOF;

    Float_t d2sigmacascnegkaontpctof = nSigmaCascNegKaonTPC * nSigmaCascNegKaonTPC + nSigmaCascNegKaonTOF * nSigmaCascNegKaonTPC;
    Float_t d2sigmacascnegpiontpdtof = nSigmaCascNegPionTPC * nSigmaCascNegPionTPC + nSigmaCascNegPionTOF * nSigmaCascNegPionTOF;
    Float_t d2sigmacascnegprotontpctof = nSigmaCascNegProtonTPC * nSigmaCascNegProtonTPC + nSigmaCascNegProtonTOF * nSigmaCascNegProtonTOF;

    Bool_t fPIDTOF_cascbach;
    Bool_t fPIDTOF_cascpos;
    Bool_t fPIDTOF_cascneg;
    if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackCascBach) == 0)      fPIDTOF_cascbach = kFALSE;
    else      fPIDTOF_cascbach = kTRUE;
    if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackCascPos) == 0)     fPIDTOF_cascpos = kFALSE;
    else      fPIDTOF_cascpos = kTRUE;
    if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, myTrackCascNeg) == 0)     fPIDTOF_cascneg = kFALSE;
    else      fPIDTOF_cascneg = kTRUE;

    Float_t nSigmaCascBachKaonTOFTPC;
    Float_t nSigmaCascBachPionTOFTPC;
    Float_t nSigmaCascBachProtonTOFTPC;

    Float_t nSigmaCascPosKaonTOFTPC;
    Float_t nSigmaCascPosPionTOFTPC;
    Float_t nSigmaCascPosProtonTOFTPC;

    Float_t nSigmaCascNegKaonTOFTPC;
    Float_t nSigmaCascNegPionTOFTPC;
    Float_t nSigmaCascNegProtonTOFTPC;
    if (fPIDTOF_cascbach && myTrackCascBach->Pt() > 0.5) {
      nSigmaCascBachKaonTOFTPC = TMath::Sqrt(d2sigmacascbachkaontpctof);
      nSigmaCascBachPionTOFTPC = TMath::Sqrt(d2sigmacascbachpiontpctof);
      nSigmaCascBachProtonTOFTPC = TMath::Sqrt(d2sigmacascbachprotontpctof);
    } else {
      nSigmaCascBachKaonTOFTPC = TMath::Abs(nSigmaCascBachKaonTPC);
      nSigmaCascBachPionTOFTPC = TMath::Abs(nSigmaCascBachPionTPC);
      nSigmaCascBachProtonTOFTPC = TMath::Abs(nSigmaCascBachProtonTPC);
    }

    if (fPIDTOF_cascpos && myTrackCascPos->Pt() > 0.5) {
      nSigmaCascPosKaonTOFTPC = TMath::Sqrt(d2sigmacascbachkaontpctof);
      nSigmaCascPosPionTOFTPC = TMath::Sqrt(d2sigmacascbachpiontpctof);
      nSigmaCascPosProtonTOFTPC = TMath::Sqrt(d2sigmacascbachprotontpctof);
    } else {
      nSigmaCascPosKaonTOFTPC = TMath::Abs(nSigmaCascPosKaonTPC);
      nSigmaCascPosPionTOFTPC = TMath::Abs(nSigmaCascPosPionTPC);
      nSigmaCascPosProtonTOFTPC = TMath::Abs(nSigmaCascPosProtonTPC);
    }

    if (fPIDTOF_cascpos && myTrackCascNeg->Pt() > 0.5) {
      nSigmaCascNegKaonTOFTPC = TMath::Sqrt(d2sigmacascbachkaontpctof);
      nSigmaCascNegPionTOFTPC = TMath::Sqrt(d2sigmacascbachpiontpctof);
      nSigmaCascNegProtonTOFTPC = TMath::Sqrt(d2sigmacascbachprotontpctof);
    } else {
      nSigmaCascNegKaonTOFTPC = TMath::Abs(nSigmaCascNegKaonTPC);
      nSigmaCascNegPionTOFTPC = TMath::Abs(nSigmaCascNegPionTPC);
      nSigmaCascNegProtonTOFTPC = TMath::Abs(nSigmaCascNegProtonTPC);
    }

    Int_t SpPID_cascbach = 999; // 1:kaon 2:pion 3:proton
    Int_t SpPID_cascpos = 999;
    Int_t SpPID_cascneg = 999;

    /*
    if( (nSigmaCascBachKaonTOFTPC<5.0) &&
    (nSigmaCascBachKaonTOFTPC<nSigmaCascBachPionTOFTPC) &&
    (nSigmaCascBachKaonTOFTPC<nSigmaCascBachProtonTOFTPC) ) SpPID_cascbach=1;
    if( (nSigmaCascBachPionTOFTPC<5.0) &&
    (nSigmaCascBachPionTOFTPC<nSigmaCascBachKaonTOFTPC) &&
    (nSigmaCascBachPionTOFTPC<nSigmaCascBachProtonTOFTPC) ) SpPID_cascbach=2;

    if( (nSigmaCascPosPionTOFTPC<5.0) &&
    (nSigmaCascPosPionTOFTPC<nSigmaCascPosKaonTOFTPC) &&
    (nSigmaCascPosPionTOFTPC<nSigmaCascPosProtonTOFTPC) ) SpPID_cascpos=2; if(
    (nSigmaCascPosProtonTOFTPC<5.0) &&
    (nSigmaCascPosPionTOFTPC<nSigmaCascPosKaonTOFTPC) &&
    (nSigmaCascPosProtonTOFTPC<nSigmaCascPosPionTOFTPC) ) SpPID_cascpos=3;

    if( (nSigmaCascNegPionTOFTPC<5.0) &&
    (nSigmaCascNegPionTOFTPC<nSigmaCascNegKaonTOFTPC) &&
    (nSigmaCascNegPionTOFTPC<nSigmaCascNegProtonTOFTPC) ) SpPID_cascneg=2; if(
    (nSigmaCascNegProtonTOFTPC<5.0) &&
    (nSigmaCascNegPionTOFTPC<nSigmaCascNegKaonTOFTPC) &&
    (nSigmaCascNegProtonTOFTPC<nSigmaCascNegPionTOFTPC) ) SpPID_cascneg=3;
    */

    if (TMath::Abs(nSigmaCascBachKaonTPC) < 5.0)  SpPID_cascbach = 1;
    if (TMath::Abs(nSigmaCascBachPionTPC) < 5.0)      SpPID_cascbach = 2;
    if (TMath::Abs(nSigmaCascPosPionTPC) < 5.0)      SpPID_cascpos = 2;
    if (TMath::Abs(nSigmaCascPosProtonTPC) < 5.0)      SpPID_cascpos = 3;
    if (TMath::Abs(nSigmaCascNegPionTPC) < 5.0)      SpPID_cascneg = 2;
    if (TMath::Abs(nSigmaCascNegProtonTPC) < 5.0)      SpPID_cascneg = 3;

    Bool_t cutXiMinusPID =    (SpPID_cascbach == 2 && SpPID_cascneg == 2 && SpPID_cascpos == 3);
    Bool_t cutXiPlusPID =        (SpPID_cascbach == 2 && SpPID_cascneg == 3 && SpPID_cascpos == 2);
    Bool_t cutOmegaMinusPID = (SpPID_cascbach == 1 && SpPID_cascneg == 2 && SpPID_cascpos == 3);
    Bool_t cutOmegaPlusPID =        (SpPID_cascbach == 1 && SpPID_cascneg == 3 && SpPID_cascpos == 2);

    Bool_t cutMassXi = ((lInvMassxi > 1.2) && (lInvMassxi < 1.4));
    Bool_t cutMassOmega = ((lInvMassomega > 1.55) && (lInvMassomega < 1.75));



    Bool_t cutXiMinussc = cutXictau && cutXiMinusPID && topoXi && (!cutMassOmega);
    Bool_t cutXiPlussc = cutXictau &&  cutXiPlusPID && topoXi && (!cutMassOmega);
    Bool_t cutOmegaMinussc = cutOmegactau && cutOmegaMinusPID && topoOmega && (!cutMassXi);
    Bool_t cutOmegaPlussc = cutOmegactau && cutOmegaPlusPID && topoOmega && (!cutMassXi);

    Double_t lPtCas = TMath::Sqrt((casc->MomXiX()) * (casc->MomXiX()) +  (casc->MomXiY()) * (casc->MomXiY()));

    Double_t spXiMinus[3] = {lInvMassXiMinus, lPtCas, lCentrality};
    Double_t spXiPlus[3] = {lInvMassXiPlus, lPtCas, lCentrality};
    Double_t spOmegaMinus[3] = {lInvMassOmegaMinus, lPtCas, lCentrality};
    Double_t spOmegaPlus[3] = {lInvMassOmegaPlus, lPtCas, lCentrality};

    if (cutXiMinussc && cutMassXi)            fHistMassXiMinus->Fill(spXiMinus);
    if (cutXiPlussc && cutMassXi)             fHistMassXiPlus->Fill(spXiPlus);
    if (cutOmegaMinussc && cutMassOmega)      fHistMassOmegaMinus->Fill(spOmegaMinus);
    if (cutOmegaPlussc && cutMassOmega)       fHistMassOmegaPlus->Fill(spOmegaPlus);

//    if (cutOmegaMinussc)   cout << lInvMassOmegaMinus << endl;
    Double_t etaCas =  0.5 * TMath::Log((lPCas + casc->MomXiZ()) / (lPCas - casc->MomXiZ()));
    if (TMath::Abs(etaCas) > 0.8) continue;
    Double_t phiCas = TMath::ATan2(casc->MomXiY(), casc->MomXiX());
    if (phiCas < 0.)  phiCas += 2 * TMath::Pi();

    Int_t SpCas = -999;
    Bool_t IDXi = (cutXiMinussc || cutXiPlussc) && cutMassXi;
    Bool_t IDOmega = (cutOmegaMinussc || cutOmegaPlussc) && cutMassOmega;
    Bool_t XiSignal = ((lInvMassxi > 1.314) && (lInvMassxi < 1.33));
    Bool_t XiSignalBg1 = ((lInvMassxi < 1.314) && (lInvMassxi > 1.2));
    Bool_t XiSignalBg2 = ((lInvMassxi > 1.33) && (lInvMassxi < 1.4));
    Bool_t OmegaSignal = ((lInvMassomega > 1.667) && (lInvMassomega < 1.678));
    Bool_t OmegaSignalBg1 = ((lInvMassomega < 1.667) && (lInvMassomega > 1.55));
    Bool_t OmegaSignalBg2 = ((lInvMassomega > 1.678) && (lInvMassomega < 1.75));
    if (IDXi && XiSignal)         SpCas = 0;
    if (IDXi && XiSignalBg1)      SpCas = 1;
    if (IDXi && XiSignalBg2)      SpCas = 2;
    if (IDOmega && OmegaSignal)      SpCas = 3;
    if (IDOmega && OmegaSignalBg1)      SpCas = 4;
    if (IDOmega && OmegaSignalBg2)      SpCas = 5;

    if(SpCas<0) continue;
    Double_t spCascadeQA[5] = {casc->Pt(),casc->Eta(), casc->Phi(), (Float_t)SpCas,  lCentrality};
    fHist_CascadeQA->Fill(spCascadeQA);
    tracks->Add(new AliAssociatedTrackYSMC(1., etaCas, phiCas, lPtCas, myTrackCascPos->GetID(),myTrackCascNeg->GetID(), myTrackCascBach->GetID(), SpCas, 1));
  }
  return tracks;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedCascade(
    const AliAODcascade *casc) {
  if (!casc)
    return kFALSE;
  /*
  AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
  AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
  AliAODTrack *btrack = (AliAODTrack*)
  (casc->GetDecayVertexXi()->GetDaughter(0)); if(!ptrack||!ntrack||!btrack)
  return kFALSE;
  */
  const Double_t mLPDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Double_t mxiPDG = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Double_t momegaPDG = TDatabasePDG::Instance()->GetParticle(3334)->Mass();

  Double_t massLambda=-1.;
  Double_t massAntiLambda=-1.;
  if(casc->ChargeXi()<0){
    massLambda = casc->MassLambda();
    if(TMath::Abs(massLambda-mLPDG)>0.00775) return kFALSE;
  }else{
    massAntiLambda = casc->MassAntiLambda();
    if(TMath::Abs(massAntiLambda-mLPDG)>0.00775) return kFALSE;
  }
  Double_t lPosXi[3];
  lPosXi[0] = casc->DecayVertexXiX();
  lPosXi[1] = casc->DecayVertexXiY();
  lPosXi[2] = casc->DecayVertexXiZ();
  Double_t decayvertXi = TMath::Sqrt(lPosXi[0] * lPosXi[0] + lPosXi[1] * lPosXi[1]);
  Double_t lPosV0[3];
  lPosV0[0] = casc->DecayVertexV0X();
  lPosV0[1] = casc->DecayVertexV0Y();
  lPosV0[2] = casc->DecayVertexV0Z();
  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0] * lPosV0[0] + lPosV0[1] * lPosV0[1]);
  if (decayvertV0 < 1.55)
    return kFALSE; // Decay vertex V0
  if (decayvertXi < 0.29)
    return kFALSE; // Decay Vertex Xi
  Double_t lDcaXiDaughters = casc->DcaXiDaughters();
  Double_t lDcaV0Daughters = casc->DcaV0Daughters();
  if (lDcaXiDaughters > 1.83)
    return kFALSE; // DCA between Xi Daughters
  if (lDcaV0Daughters > 1.33)
    return kFALSE; // DCA between V0 Daughters
  Double_t lDcaBachToPrimVertex = casc->DcaBachToPrimVertex();
  Double_t lDcaV0ToPrimVertex = casc->DcaV0ToPrimVertex();
  Double_t lDcaPosToPrimVertex = casc->DcaPosToPrimVertex();
  Double_t lDcaNegToPrimVertex = casc->DcaNegToPrimVertex();
  if (lDcaBachToPrimVertex < 0.0146)
    return kFALSE; // DCA between Bach track and Primary vertex
  if (lDcaV0ToPrimVertex < 0.02)
    return kFALSE; // DCA between V0 deday vertex and primary vertex
  if (lDcaPosToPrimVertex < 0.061)
    return kFALSE; // DCA between Pos track and Primary vertex
  if (lDcaNegToPrimVertex < 0.061)
    return kFALSE; // DCA between Neg track and Primary vertex
  Double_t lXiCosineOfPointingAngle = casc->CosPointingAngleXi(tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], tPrimaryVtxPosition[2]);
  Double_t lV0CosineOfPointingAngleXi = casc->CosPointingAngle(lPosXi);

  if (lXiCosineOfPointingAngle < 0.9813)
    return kFALSE; // Pointing angle of Xi
    if (lV0CosineOfPointingAngleXi < 0.97)
  return kFALSE; /// Pointing angle of V0

  return kTRUE;
}


Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedCascadeOmega(const AliAODcascade *casc) {
  if (!casc)    return kFALSE;
  const Double_t mLPDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Double_t mxiPDG = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Double_t momegaPDG = TDatabasePDG::Instance()->GetParticle(3334)->Mass();

  Double_t massLambda=-1.;
  Double_t massAntiLambda=-1.;
  if(casc->ChargeXi()<0){
    massLambda = casc->MassLambda();
    if(TMath::Abs(massLambda-mLPDG)>0.00775) return kFALSE;
  }else{
    massAntiLambda = casc->MassAntiLambda();
    if(TMath::Abs(massAntiLambda-mLPDG)>0.00775) return kFALSE;
  }
  Double_t lPosXi[3];
  lPosXi[0] = casc->DecayVertexXiX();
  lPosXi[1] = casc->DecayVertexXiY();
  lPosXi[2] = casc->DecayVertexXiZ();
  Double_t decayvertXi = TMath::Sqrt(lPosXi[0] * lPosXi[0] + lPosXi[1] * lPosXi[1]);
  Double_t lPosV0[3];
  lPosV0[0] = casc->DecayVertexV0X();
  lPosV0[1] = casc->DecayVertexV0Y();
  lPosV0[2] = casc->DecayVertexV0Z();
  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0] * lPosV0[0] + lPosV0[1] * lPosV0[1]);
  if (decayvertV0 < 1.918)
    return kFALSE; // Decay vertex V0
  if (decayvertXi < 0.59)
    return kFALSE; // Decay Vertex Xi
  Double_t lDcaXiDaughters = casc->DcaXiDaughters();
  Double_t lDcaV0Daughters = casc->DcaV0Daughters();
  if (lDcaXiDaughters > 1.091)
    return kFALSE; // DCA between Xi Daughters
  if (lDcaV0Daughters > 1.206)
    return kFALSE; // DCA between V0 Daughters
  Double_t lDcaBachToPrimVertex = casc->DcaBachToPrimVertex();
  Double_t lDcaV0ToPrimVertex = casc->DcaV0ToPrimVertex();
  Double_t lDcaPosToPrimVertex = casc->DcaPosToPrimVertex();
  Double_t lDcaNegToPrimVertex = casc->DcaNegToPrimVertex();
  if (lDcaBachToPrimVertex < 0.041)
    return kFALSE; // DCA between Bach track and Primary vertex
  if (lDcaV0ToPrimVertex < 0.068)
    return kFALSE; // DCA between V0 deday vertex and primary vertex
  if (lDcaPosToPrimVertex < 0.061)
    return kFALSE; // DCA between Pos track and Primary vertex
  if (lDcaNegToPrimVertex < 0.061)
    return kFALSE; // DCA between Neg track and Primary vertex
  Double_t lXiCosineOfPointingAngle = casc->CosPointingAngleXi(tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], tPrimaryVtxPosition[2]);
  Double_t lV0CosineOfPointingAngleXi = casc->CosPointingAngle(lPosXi);

  if (lXiCosineOfPointingAngle < 0.9811)
    return kFALSE; // Pointing angle of Xi
    if (lV0CosineOfPointingAngleXi < 0.933)
  return kFALSE; /// Pointing angle of V0

  return kTRUE;
}


Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedPhiDaughterTrack(const AliAODTrack *aodTrack) {
  if (!aodTrack->TestFilterMask(BIT(ffilterbit)))    return kFALSE;
  if (aodTrack->Pt() < 0.15)    return kFALSE;
  if (TMath::Abs(aodTrack->Eta()) > 0.8)    return kFALSE;
  if (aodTrack->Charge() == 0)    return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedDaughterTrack(const AliAODTrack *itrack) {//apply for V0s and Cascade
  if (TMath::Abs(itrack->Eta()) > 0.8)  return kFALSE;
  //  if(itrack->Eta()>fEtaMaxDaughter || itrack->Eta()<fEtaMinDaughter) return kFALSE;
  if (!itrack->IsOn(AliAODTrack::kTPCrefit))  return kFALSE;
  Float_t nCrossedRowsTPC = itrack->GetTPCClusterInfo(2, 1);
  if (nCrossedRowsTPC < fclustermin)  return kFALSE;
  Int_t findable = itrack->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC / findable < fratiocluster) return kFALSE;
  if (itrack->GetKinkIndex(0) > 0) return kFALSE;
  //  if(itrack->Pt()<0.15) return kFALSE;
  if(itrack->Pt()<0.1) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedTrack(const AliAODTrack *aodTrack) {
  if (!aodTrack)
    return kFALSE;
  //  if(!aodTrack->TestFilterMask(BIT(5))) return kFALSE; // standard cut with
  //  tight DCA cut
  //  if (!aodTrack->TestFilterMask(BIT(ffilterbit)))
  if (!aodTrack->TestFilterBit(ffilterbit)) return kFALSE;
    /*
  if (!aodTrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC =aodTrack->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=aodTrack->GetTPCNclsF();
  if (findable <= 0)return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;s
  */
  if (aodTrack->Pt() < fPtMin)    return kFALSE;

  //  if (aodTrack->Pt() > 3.) return kFALSE;
  if (TMath::Abs(aodTrack->Eta()) > fEtaMax)
    return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedV0(const AliAODv0 *aodV0) {
  /*
    Pseudo rapidity V0 < 0.8
    daughter tracks DCA to P >0.06cm
    DCA between daughters <1 sigma
    V0 2D decay radius from 0.5cm to 200cm
  */
  if (!aodV0)  return kFALSE;
  if(aodV0->GetOnFlyStatus()) fHist_V0Stat->Fill(1); //Onfly
  else fHist_V0Stat->Fill(2); //Onfly
    if (fOnfly) {
      if (!aodV0->GetOnFlyStatus()) return kFALSE; // onfly reconstraction only
    } else {
    if (aodV0->GetOnFlyStatus())  return kFALSE; // ofline reconstraction only
    }
  Double_t leta_v0 = aodV0->PseudoRapV0();
  if (TMath::Abs(leta_v0) > 0.8)  return kFALSE; // default 0.8
  //  if (leta_v0 < fEtaMinV0 || fEtaMaxV0<leta_v0)  return kFALSE; // default 0.8
  else fHist_V0Stat->Fill(3);

  // DCA to primary vertex
  Float_t xyn = aodV0->DcaNegToPrimVertex();
  Float_t xyp = aodV0->DcaPosToPrimVertex();
  if ((TMath::Abs(xyn) > fdcaDaughtersToPrimVtx)  || (TMath::Abs(xyp) > fdcaDaughtersToPrimVtx))  fHist_V0Stat->Fill(4);
  if (TMath::Abs(xyn) < fdcaDaughtersToPrimVtx)    return kFALSE; // default  0.06
  if (TMath::Abs(xyp) < fdcaDaughtersToPrimVtx)    return kFALSE; // default  0.06
  // DCA between dauther tracks
  Double_t dca = aodV0->DcaV0Daughters();
  if (dca > fdcaBetweenDaughters)  return kFALSE; // default 1.0
  else 	 fHist_V0Stat->Fill(5);
  // Fiducial volume
  Double_t xyz[3];
  aodV0->GetSecondaryVtx(xyz);
  Double_t r2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
  if (r2 > fRadiMin * fRadiMin && r2 < fRadiMax * fRadiMax)  	 fHist_V0Stat->Fill(6);
  if (r2 < fRadiMin * fRadiMin)    return kFALSE; // default 0.5 cm
  if (r2 > fRadiMax * fRadiMax)    return kFALSE; // default 100cm


  return kTRUE;
}
Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedDecayLength(const AliAODv0*aodv0,Double_t mass,Double_t maxctau){
  Double_t lDVx = aodv0->GetSecVtxX();
  Double_t lDVy = aodv0->GetSecVtxY();
  Double_t lDVz = aodv0->GetSecVtxZ();
  Double_t lV0DecayLength =TMath::Sqrt(TMath::Power(lDVx - tPrimaryVtxPosition[0], 2) +  TMath::Power(lDVy - tPrimaryVtxPosition[1], 2) + TMath::Power(lDVz - tPrimaryVtxPosition[2], 2));
  Double_t lPV0 = TMath::Sqrt((aodv0->Pt()) * (aodv0->Pt()) +(aodv0->Pz()) * (aodv0->Pz()));
  Double_t lcTau = (lV0DecayLength * mass) / lPV0;
  if(lcTau>maxctau) return kFALSE;
  return kTRUE;
}

void AliAnalysisTaskSEpPbCorrelationsMCYS::FillCorrelationTracks( Double_t centrality, TObjArray *triggerArray, TObjArray *selectedTrackArray,AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut, Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,Float_t bSign, Int_t step) {
  twoTrackEfficiencyCut=kFALSE;
  twoTrackEfficiencyCutValue=0;
  fTwoTrackCutMinRadius=0;
  bSign=0;
  step=1;//default
  if (!triggerHist || !associateHist)    return;
  if(fAnaMode=="TPCTPC"){
      Double_t binscontTrig[3];
      Double_t binscont[6];
      for (Int_t i = 0; i < triggerArray->GetEntriesFast(); i++) {
	AliAssociatedTrackYSMC *trigger = (AliAssociatedTrackYSMC *)triggerArray->At(i);
	if (!trigger)    continue;
	Float_t triggerPt = trigger->Pt();
	Float_t triggerEta = trigger->Eta();
	Float_t triggerPhi = trigger->Phi();
	Int_t trigFirstID = trigger->GetIDFirstDaughter();
	Int_t trigSecondID = trigger->GetIDSecondDaughter();
	Int_t trigID = trigger->GetID();
	
	Float_t efficiency=999;
	Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
	if(fMCclosure && !fprim && !fextractsec){
	Int_t iPt=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(triggerPt);
	Int_t iEta=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(triggerEta);
	Int_t iPhi=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(triggerPhi);
	efficiency=fhcorreffi[ivzbin-1]->GetBinContent(iPt,iEta,iPhi);
	}else efficiency=1.;
	if(efficiency==0.) continue;
	binscontTrig[0] = triggerPt;
	binscontTrig[1] = centrality;
	binscontTrig[2] = fPrimaryZVtx;
	triggerHist->Fill(binscontTrig, 0, 1./efficiency);
	for (Int_t j = 0; j < selectedTrackArray->GetEntriesFast(); j++) {
	  AliAssociatedTrackYSMC *associate =   (AliAssociatedTrackYSMC*)selectedTrackArray->At(j);
	  if (!associate)        continue;
	  Int_t AssoFirstID = associate->GetIDFirstDaughter();
	  Int_t AssoSecondID = associate->GetIDSecondDaughter();
	  if (fasso == "V0" || fasso == "Phi"){
	    if (trigID == AssoFirstID || trigID == AssoSecondID){
	      continue;
	    }
	  }
	  if (fasso == "hadron" || fasso=="PID") {
	    if (triggerPt <= associate->Pt())          continue;
	    if (trigID == associate->GetID())          continue;
	  }
	  if (fasso == "Cascade")  if (trigID == associate->GetID() || trigID == AssoFirstID ||  trigID == AssoSecondID)          continue;
	  binscont[0] = triggerEta - associate->Eta();
	  binscont[1] = associate->Pt();
	  binscont[2] = triggerPt;
	  binscont[3] = centrality;
	  binscont[4] = RangePhi(triggerPhi - associate->Phi());
	  binscont[5] = fPrimaryZVtx;
	  Int_t SpAsso = associate->WhichCandidate();
	  Float_t efficiency1=999.;
	  if(fMCclosure && !fprim && !fextractsec){
	  Int_t iPt1=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(associate->Pt());
	  Int_t iEta1=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(associate->Eta());
	  Int_t iPhi1=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(associate->Phi());
	  efficiency1=fhcorreffi[ivzbin-1]->GetBinContent(iPt1,iEta1,iPhi1);
	  }else  efficiency1=1.;
	  if(efficiency1==0.) continue;
	  

	  if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" ||  (fasso == "PID")) {
	    if (SpAsso < 0)          continue;
	    associateHist->Fill(binscont, SpAsso);
	  }else if(fasso=="hadron"){
	    associateHist->Fill(binscont, 0, 1./(efficiency*efficiency1));
	  }
	}
      }
  }else if (fAnaMode=="TPCV0A" || fAnaMode=="TPCV0C" || fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
    Double_t binscontTrig[4];
    Double_t binscont[7];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
      if(!trigger)continue;
      Float_t  triggerEta  = trigger->Eta();
      Float_t  triggerPhi  = trigger->Phi();
      Float_t  triggerPt  = trigger->Pt();
      binscontTrig[0]=triggerPt;
      binscontTrig[1]=centrality;
      binscontTrig[2]=fPrimaryZVtx;
      binscontTrig[3]=triggerEta;
      Float_t efficiency=999;
      Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
      if(fMCclosure && !fprim && !fextractsec){
	Int_t iPt=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(triggerPt);
	Int_t iEta=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(triggerEta);
	Int_t iPhi=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(triggerPhi);
	efficiency=fhcorreffi[ivzbin-1]->GetBinContent(iPt,iEta,iPhi);
      }else efficiency=1.;
      if(efficiency==0.) continue;
      
      Int_t SpAsso= trigger->WhichCandidate();
      triggerHist->Fill(binscontTrig,SpAsso,1./efficiency);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) selectedTrackArray->At(j);
        if(!associate)continue;
        Float_t associatemultiplicity=associate->Multiplicity();
        Float_t assophi=associate->Phi();
        Float_t assoeta=associate->Eta();
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=triggerPt;
        binscont[2]=associate->Eta();
        binscont[3]=centrality;
        binscont[4]=RangePhi(triggerPhi-associate->Phi());
        binscont[5]=fPrimaryZVtx;
	binscont[6]=triggerEta;
	
        if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" ||  (fasso == "PID")) {
          if (SpAsso < 0)          continue;
          associateHist->Fill(binscont, SpAsso,(Double_t)associate->Multiplicity()/efficiency);
        }else if(fasso=="hadron"){
          associateHist->Fill(binscont, 0, (Double_t)associate->Multiplicity()/efficiency);
        }
      }
    }
  }else if (fAnaMode=="ITSFMD" || fAnaMode=="ITSFMDC"){
    Double_t binscontTrig[4];
    Double_t binscont[6];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
      if(!trigger)continue;
      Float_t  triggerEta  = trigger->Eta();
      Float_t  triggerPhi  = trigger->Phi();
      Float_t  triggerPt  = 0.5;
      binscontTrig[0]=triggerPt;
      binscontTrig[1]=centrality;
      binscontTrig[2]=fPrimaryZVtx;
      binscontTrig[3]=triggerEta;
      Int_t SpAsso= trigger->WhichCandidate();
      triggerHist->Fill(binscontTrig,SpAsso);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) selectedTrackArray->At(j);
        if(!associate)continue;
        Float_t associatemultiplicity=associate->Multiplicity();
        Float_t assophi=associate->Phi();
        Float_t assoeta=associate->Eta();
        binscont[0]=triggerEta-associate->Eta();
	binscont[1]=centrality;	
	binscont[2]=associate->Eta();
	binscont[3]=RangePhi(triggerPhi-associate->Phi());
        binscont[4]=fPrimaryZVtx;
        binscont[5]=triggerEta;
	associateHist->Fill(binscont, 0, (Double_t)associate->Multiplicity());
      }
    }
 } else if(fAnaMode=="FMDFMD"){
    Double_t binscontTrig[3];
    Double_t binscont[6];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
      if(!trigger)continue;
      Float_t  triggerMultiplicity= trigger->Multiplicity();
      Float_t  triggerEta  = trigger->Eta();
      Float_t  triggerPhi  = trigger->Phi();
      binscontTrig[0]=centrality;
      binscontTrig[1]=triggerEta;
      binscontTrig[2]=fPrimaryZVtx;
      triggerHist->Fill(binscontTrig,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) selectedTrackArray->At(j);
        if(!associate)continue;
        Float_t assophi=associate->Phi();
        Float_t assoeta=associate->Eta();
        Float_t mult=associate->Multiplicity()*triggerMultiplicity;
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=associate->Eta();
        binscont[2]=triggerEta;
        binscont[3]=centrality;
	binscont[4]=RangePhi_FMD(triggerPhi-associate->Phi());
        binscont[5]=fPrimaryZVtx;
	
        Float_t dphivzero=triggerPhi-associate->Phi();
        if(triggerPhi==assophi && triggerEta==assoeta) continue;
        associateHist->Fill(binscont,0,(Double_t)mult);
      }
    }
  }else if(fAnaMode=="SECA"||fAnaMode=="SECC"){
    Double_t binscontTrig[3];
    Double_t binscont[5];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
      if(!trigger)continue;
      Float_t  triggerMultiplicity= trigger->Multiplicity();
      Float_t  triggerEta  = trigger->Eta();
      Float_t  triggerPhi  = trigger->Phi();
      binscontTrig[0]=centrality;
      binscontTrig[1]=triggerEta;
      binscontTrig[2]=fPrimaryZVtx;
      triggerHist->Fill(binscontTrig,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) selectedTrackArray->At(j);
        if(!associate)continue;
        Float_t assophi=associate->Phi();
        Float_t assoeta=associate->Eta();
        Float_t mult=associate->Multiplicity()*triggerMultiplicity;
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=triggerEta;
        binscont[2]=centrality;
        binscont[3]=RangePhi_FMD(triggerPhi-associate->Phi());
        binscont[4]=fPrimaryZVtx;
        //        if(triggerEta>2.7 && triggerEta<2.8)cout<<triggerEta<<" "<<associate->Eta()<<" "<<binscont[0]<<endl;
        Float_t dphivzero=triggerPhi-associate->Phi();
        if(triggerPhi==assophi && triggerEta==assoeta) continue;
        associateHist->Fill(binscont,0,(Double_t)mult);
      }
    }
  }else if(fAnaMode=="V0AV0C"){
    Double_t  binscont1[4];
    Double_t  binscontTrig1[3];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t  triggerMultiplicity= trigger->Multiplicity();
      Double_t  triggerEta  = trigger->Eta();
      Double_t  triggerPhi  = trigger->Phi();
      binscontTrig1[0]=centrality;
      binscontTrig1[1]=triggerEta;
      binscontTrig1[2]=triggerPhi;
      triggerHist->Fill(binscontTrig1,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) selectedTrackArray->At(j);
        if(!associate)continue;
        Double_t associatemultiplicity=associate->Multiplicity();
        Double_t assophi=associate->Phi();
        Double_t assoeta=associate->Eta();
        binscont1[0]=associate->Eta();
        binscont1[1]=triggerEta;
        binscont1[2]=centrality;
        binscont1[3]=RangePhi2(triggerPhi-associate->Phi());
        if(triggerPhi==assophi && triggerEta==assoeta) continue;
        Double_t dphivzero=triggerPhi-associate->Phi();
        associateHist->Fill(binscont1,0,(Double_t)triggerMultiplicity*associatemultiplicity);
      }
    }
  }


}

void AliAnalysisTaskSEpPbCorrelationsMCYS::FillCorrelationTracksMixing(Double_t centrality, Double_t pvxMix, Double_t poolmax, Double_t poolmin,    TObjArray *triggerArray, TObjArray *selectedTrackArray, AliTHn *triggerHist,    AliTHn *associateHist, Bool_t twoTrackEfficiencyCut,    Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,    Float_t bSign, Int_t step) {
 
  Bool_t twoTrackEfficiencyCut_1= twoTrackEfficiencyCut;
  Double_t   twoTrackEfficiencyCutValue_1=  twoTrackEfficiencyCutValue;
  Double_t fTwoTrackCutMinRadius1=fTwoTrackCutMinRadius;
  Float_t bSign1=bSign;
  Int_t step1=step;
  
  Double_t poolmax1=poolmax;
  Double_t poolmin1=poolmin;
  /*
  if (!triggerHist || !associateHist){
    return;
  }
  */
  Double_t counterMix = 0;
  AliEventPool *pool = fPoolMgr->GetEventPool(centrality, pvxMix);
  if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality,
		  pvxMix));
  }

  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    mixedDist ->Fill(centrality, pool->NTracksInPool());
    mixedDist2->Fill(centrality, pool->GetCurrentNEvents());
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix = 0; jMix < nMix; jMix++) {  
      TObjArray *mixEvents = pool->GetEvent(jMix);
      if(fAnaMode=="TPCTPC"){
        Double_t binscontTrig[2];
        Double_t binscont[6];
	for (Int_t i = 0; i < triggerArray->GetEntriesFast(); i++) {
          AliAssociatedTrackYSMC *trig = (AliAssociatedTrackYSMC *)triggerArray->At(i);
          if (!trig)          continue;
          Double_t triggerPhi = trig->Phi();
          Double_t triggerEta = trig->Eta();
          Double_t triggerPt = trig->Pt();
          counterMix++;
          binscontTrig[0] = triggerPt;
          binscontTrig[1] = centrality;
	  Float_t efficiency=999;
	  Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
	  if(fMCclosure && !fprim && !fextractsec){
	    Int_t iPt=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(triggerPt);
	    Int_t iEta=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(triggerEta);
	    Int_t iPhi=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(triggerPhi);
	    efficiency=fhcorreffi[ivzbin-1]->GetBinContent(iPt,iEta,iPhi);
	  }else efficiency=1.;
	  if(efficiency==0.) continue;
      
          triggerHist->Fill(binscontTrig, 0);
          for (Int_t j = 0; j < mixEvents->GetEntriesFast(); j++) {
            AliAssociatedTrackYSMC *associate =  (AliAssociatedTrackYSMC *)mixEvents->At(j);
            if (!associate) continue;
	    binscont[0] = triggerEta - associate->Eta();
	    binscont[1] = associate->Pt();
	    binscont[2] = triggerPt;
	    binscont[3] = centrality;
	    binscont[4] = RangePhi(triggerPhi - associate->Phi());
	    binscont[5] = pvxMix;
	    Float_t efficiency1=999.;
	    if(fMCclosure && !fprim && !fextractsec){
	      Int_t iPt1=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(associate->Pt());
	      Int_t iEta1=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(associate->Eta());
	      Int_t iPhi1=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(associate->Phi());
	      efficiency1=fhcorreffi[ivzbin-1]->GetBinContent(iPt1,iEta1,iPhi1);
	    }else  efficiency1=1.;
	    if(efficiency1==0.) continue;
	    
	    Int_t SpAsso = associate->WhichCandidate();
	    if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" || (fasso == "PID")) {
	      if (SpAsso < 0)   continue;
	      associateHist->Fill(binscont, SpAsso, 1. / (Double_t)nMix/(efficiency*efficiency1));
	    }else if(fasso=="hadron"){
	      associateHist->Fill(binscont, 0,1./(Double_t)nMix/(efficiency*efficiency1));
	    }
	  }
	}
      }else if(fAnaMode=="TPCV0A" || fAnaMode=="TPCV0C" || fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
        Double_t binscontTrig[2];
        Double_t binscont[7];
        for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
          AliAssociatedTrackYSMC* trigger =(AliAssociatedTrackYSMC*) triggerArray->At(i);
          if(!trigger)continue;
          Double_t triggerPt   = trigger->Pt();
          Double_t triggerEta  = trigger->Eta();
          Double_t triggerPhi  = trigger->Phi();
	  Int_t SpAsso=trigger->WhichCandidate();
          counterMix++;
          binscontTrig[0]=triggerPt;
          binscontTrig[1]=centrality;
          triggerHist->Fill(binscontTrig,SpAsso);
	  Float_t efficiency=999;
	  Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
	  if(fMCclosure && !fprim && !fextractsec){
	    Int_t iPt=fhcorreffi[ivzbin-1]->GetXaxis()->FindBin(triggerPt);
	    Int_t iEta=fhcorreffi[ivzbin-1]->GetYaxis()->FindBin(triggerEta);
	    Int_t iPhi=fhcorreffi[ivzbin-1]->GetZaxis()->FindBin(triggerPhi);
	    efficiency=fhcorreffi[ivzbin-1]->GetBinContent(iPt,iEta,iPhi);
	  }else efficiency=1.;
	  if(efficiency==0.) continue;
      
	  
          for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
            AliAssociatedTrackYSMC* associate=(AliAssociatedTrackYSMC*)  mixEvents->At(j);
            
            if(!associate)continue;
            binscont[0]=triggerEta-associate->Eta();
            binscont[1]=triggerPt;
            binscont[2]=associate->Eta();
            binscont[3]=centrality;
            binscont[4]=RangePhi(triggerPhi-associate->Phi());
            binscont[5]=pvxMix;
	    binscont[6]=triggerEta;

	    if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" || (fasso == "PID")) {
              if (SpAsso < 0)   continue;
              associateHist->Fill(binscont, SpAsso, (Double_t)associate->Multiplicity()/(Double_t)nMix/efficiency);
            }else if(fasso=="hadron"){
              associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix/efficiency);
            }
          }
        }
      }else if(fAnaMode=="ITSFMD" || fAnaMode=="ITSFMDC"){
        Double_t binscontTrig[2];
        Double_t binscont[6];
        for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
          AliAssociatedTrackYSMC* trigger =(AliAssociatedTrackYSMC*) triggerArray->At(i);
          if(!trigger)continue;
          Double_t triggerPt   = 0.5;
          Double_t triggerEta  = trigger->Eta();
          Double_t triggerPhi  = trigger->Phi();
		  Int_t SpAsso=trigger->WhichCandidate();
          counterMix++;
          binscontTrig[0]=triggerPt;
          binscontTrig[1]=centrality;
          triggerHist->Fill(binscontTrig,SpAsso);
          for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
            AliAssociatedTrackYSMC* associate=(AliAssociatedTrackYSMC*)  mixEvents->At(j);
	    if(!associate)continue;
            binscont[0]=triggerEta-associate->Eta();
            binscont[1]=centrality;
            binscont[2]=associate->Eta();
            binscont[3]=RangePhi(triggerPhi-associate->Phi());
            binscont[4]=pvxMix;
            binscont[5]=triggerEta;
	    associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix);
	  }
        }
      }else if(fAnaMode=="FMDFMD"){
        Double_t binscontTrig[2];
        Double_t binscont[6];
        for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
          AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
        if(!trigger)continue;
        Float_t triggerMultiplicity= trigger->Multiplicity();
        Float_t  triggerEta  = trigger->Eta();
        Float_t  triggerPhi  = trigger->Phi();
        counterMix++;
        binscontTrig[0]=centrality;
        binscontTrig[1]=triggerEta;
	//       triggerHist->Fill(binscontTrig,step,(Double_t)triggerMultiplicity);
        for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
          AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) mixEvents->At(j);
          if(!associate)continue;
          //          Double_t associatemultiplicity=associate->Multiplicity();
          Float_t assophi=associate->Phi();
          Float_t assoeta=associate->Eta();
          Float_t mult=triggerMultiplicity*associate->Multiplicity();
          binscont[0]=triggerEta-associate->Eta();
          binscont[1]=associate->Eta();
          binscont[2]=triggerEta;
          binscont[3]=centrality;
	  binscont[4]=RangePhi_FMD(triggerPhi-associate->Phi());
	  binscont[5]=pvxMix;
          // if(triggerPhi==assophi && triggerEta==assoeta) continue;
          //  associateHist->Fill(binscont,0,(Double_t)triggerMultiplicity*associatemultiplicity/(Double_t)nMix);
          associateHist->Fill(binscont,0,mult/(Double_t)nMix);
        }
        }
      }else if(fAnaMode=="SECA"||fAnaMode=="SECC"){
        Double_t binscontTrig[2];
        Double_t binscont[5];
        for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
          AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
	  if(!trigger)continue;
	  Float_t triggerMultiplicity= trigger->Multiplicity();
	  Float_t  triggerEta  = trigger->Eta();
	  Float_t  triggerPhi  = trigger->Phi();
	  counterMix++;
	  binscontTrig[0]=centrality;
	  binscontTrig[1]=triggerEta;
	  triggerHist->Fill(binscontTrig,step,(Double_t)triggerMultiplicity);
	  for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
	    AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) mixEvents->At(j);
	    if(!associate)continue;
	    //          Double_t associatemultiplicity=associate->Multiplicity();
	    Float_t assophi=associate->Phi();
	    Float_t assoeta=associate->Eta();
	    Float_t mult=triggerMultiplicity*associate->Multiplicity();
	    binscont[0]=triggerEta-associate->Eta();
	    binscont[1]=triggerEta;
	    binscont[2]=centrality;
	    binscont[3]=RangePhi_FMD(triggerPhi-associate->Phi());
	    binscont[4]=pvxMix;
	    //     if(triggerPhi==assophi && triggerEta==assoeta) continue;
	    //  associateHist->Fill(binscont,0,(Double_t)triggerMultiplicity*associatemultiplicity/(Double_t)nMix);
	    associateHist->Fill(binscont,0,mult/(Double_t)nMix);
	  }
        }
      }else if (fAnaMode=="V0AV0C"){
	Double_t  binscont1[4];
	Double_t  binscontTrig1[3];
	for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
	  AliAssociatedTrackYSMC* trigger = (AliAssociatedTrackYSMC*) triggerArray->At(i);
	  if(!trigger)continue;
	  Double_t  triggerMultiplicity= trigger->Multiplicity();
	  Double_t  triggerEta  = trigger->Eta();
	  Double_t  triggerPhi  = trigger->Phi();
	  counterMix++;
	  binscontTrig1[0]=centrality;
	  binscontTrig1[1]=triggerEta;
	  binscontTrig1[2]=triggerPhi;
	  triggerHist->Fill(binscontTrig1,step,(Double_t)triggerMultiplicity);
	  for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
	    AliAssociatedTrackYSMC* associate = (AliAssociatedTrackYSMC*) mixEvents->At(j);
	    if(!associate)continue;
	    Double_t associatemultiplicity=associate->Multiplicity();
	    Double_t assophi=associate->Phi();
	    Double_t assoeta=associate->Eta();
	    binscont1[0]=associate->Eta();
	    binscont1[1]=triggerEta;
	    binscont1[2]=centrality;
	    binscont1[3]=RangePhi2(triggerPhi-associate->Phi());
	    if(triggerPhi==assophi && triggerEta==assoeta) continue;
	    associateHist->Fill(binscont1,0,(Double_t)triggerMultiplicity*associatemultiplicity/(Double_t)nMix);
	  }
	}
      }
    }
  }
  
  TObjArray* tracksClone=CloneTrack(selectedTrackArray);
  pool->UpdatePool(tracksClone);
  
}
TObjArray* AliAnalysisTaskSEpPbCorrelationsMCYS::CloneTrack(TObjArray*selectedTrackArray){
  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliAssociatedTrackYSMC *particle =  (AliAssociatedTrackYSMC *)selectedTrackArray->At(i);
    tracksClone->Add(new AliAssociatedTrackYSMC(particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
					      particle->GetID(), particle->GetIDFirstDaughter(),
					      particle->GetIDSecondDaughter(), particle->WhichCandidate(),
					      particle->Multiplicity()));
  }
  
  return tracksClone;
}

Double_t AliAnalysisTaskSEpPbCorrelationsMCYS::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsMCYS::RangePhi_FMD(Double_t DPhi) {
  //if (DPhi < (-TMath::Pi() / 2 -0.0001))   DPhi += 2 * TMath::Pi();
  //if (DPhi > (3 * TMath::Pi() / 2-0.0001)) DPhi -= 2*TMath::Pi();
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsMCYS::RangePhi2(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < -1.178097)    DPhi += 2 * TMath::Pi();
  return DPhi;
}



void AliAnalysisTaskSEpPbCorrelationsMCYS::DumpTObjTable(const char* note)
{
  if(note) {
    //    printf("TObjectTable::%s",note);
  }
  //  gObjectTable->Print();
}

 Int_t AliAnalysisTaskSEpPbCorrelationsMCYS::ConvertRunNumber(Int_t run){
   if( (265309<run && run< 265525) || (267161<run && run<267166)){
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
     case 281959 : return 1;
     case 281953:return 2;
     case 281940:return 3;
     case 281939:return 4;
     case 281932:return 5;
     case 281931:return 6;
     case 281928:return 7;
     case 281920:return 8;
     case 281918:return 9;
     case 281633:return 10;
     case 281583:return 11;
     case 281581:return 12;
     case 281580:return 13;
     case 281574:return 14;
     case 281569:return 15;
     case 281568:return 16;
     case 281562:return 17;
     case 281557:return 18;
     case 281511:return 19;
     case 281509:return 20;
     case 281477:return 21;
     case 281475:return 22;
     case 281450:return 23;
     case 281449:return 24;
     case 281444:return 25;
     case 281443:return 26;
     case 281441:return 27;
     case 281415:return 28;
     case 281321:return 29;
     case 281301:return 30;
     case 281277:return 31;
     case 281275:return 32;
     case 281273:return 33;
     case 281271:return 34;
     case 281243:return 35;
     case 281242:return 36;
     case 281241:return 37;
     case 281240:return 38;
     case 281213:return 39;
     case 281212:return 40;
     case 281191:return 41;
     case 281190:return 42;
     case 281189:return 43;
     case 281181:return 44;
     case 281179:return 45;
     case 281081:return 46;
     case 281080:return 47;
     case 281062:return 48;
     case 281061:return 49;
     case 281060:return 50;
     case 280999:return 51;
     case 280998:return 52;
     case 280997:return 53;
     case 280994:return 54;
     case 280990:return 55;
     case 280947:return 56;
     case 280940:return 57;
     case 280936:return 58;
     case 280897:return 59;
     case 280880:return 60;
     case 280856:return 61;
     case 280849:return 62;
     case 280848:return 63;
     case 280847:return 64;
     case 280844:return 65;
     case 280842:return 66;
     case 280793:return 67;
     case 280792:return 68;
     case 280787:return 69;
     case 280786:return 70;
     case 280768:return 71;
     case 280767:return 72;
     case 280766:return 72;
     case 280765:return 73;
     case 280764:return 74;
     case 280763:return 75;
     case 280762:return 76;
     case 280761:return 77;
     case 280757:return 78;
     case 280756:return 79;
     case 280755:return 80;
     case 280754:return 81;
     case 280753:return 82;
     case 280729:return 83;
     case 280706:return 84;
     case 280705:return 85;
     case 280681:return 86;
     case 280679:return 87;
     case 280676:return 88;
     case 280673:return 89;
     case 280671:return 90;
     case 280650:return 91;
     case 280648:return 92;
     case 280647:return 93;
     case 280645:return 94;
     case 280639:return 95;
     case 280637:return 96;
     case 280636:return 97;
     case 280634:return 98;
     case 280613:return 99;
     case 280583:return 100;
     case 280581:return 101;
     case 280574:return 102;
     case 280551:return 103;
     case 280550:return 104;
     case 280547:return 105;
     case 280546:return 106;
     case 280519:return 107;
     case 280518:return 108;
     case 280499:return 109;
     case 280448:return 110;
     case 280447:return 111;
     case 280446:return 112;
     case 280445:return 113;
     case 280443:return 114;
     case 280419:return 115;
     case 280415:return 116;
     case 280406:return 117;
     case 280405:return 118;
     case 280403:return 119;
     case 280375:return 120;
     case 280374:return 121;
     case 280352:return 122;
     case 280351:return 123;
     case 280350:return 124;
     case 280349:return 125;
     case 280348:return 126;
     case 280312:return 127;
     case 280310:return 128;
     case 280290:return 129;
     case 280286:return 130;
     case 280285:return 131;
     case 280284:return 132;
     case 280282:return 133;
     default : return 199;
     }
   } else{
     return 199;
   }

}


