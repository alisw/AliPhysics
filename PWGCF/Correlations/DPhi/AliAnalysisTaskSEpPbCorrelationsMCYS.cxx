/**************************************************************************************************
*      Leading Charged Track+V0 Correlations.(Works for Real  Data)                               *
*      Yuko Sekiguchi * Center for Nuclear Study(CNS) , University of Tokyo                       *
*      Email:y_sekiguchi@cns.s.u-tokyo.ac.jp                                                      *
**************************************************************************************************/
#include "AliAnalysisManager.h"
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
#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliESDtrackCuts.h"
#include "AliCFContainer.h"
#include "AliGenEventHeader.h"
#include "AliTHn.h"

#include "AliAODEvent.h"
#include "AliESDAD.h"
#include "AliESDEvent.h"
#include "AliVAD.h"
#include "AliVEvent.h"
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

#include "AliAODForwardMult.h"
//#include "AliForwardUtil.h"
#include "AliAnalysisTaskSEpPbCorrelationsMCYS.h"

ClassImp(AliAnalysisTaskSEpPbCorrelationsMCYS)
ClassImp(AliAssociatedTrackYS)
ClassImp(AliMixTrackYS)
ClassImp(AliAssociatedVZEROYS)

AliAnalysisTaskSEpPbCorrelationsMCYS::AliAnalysisTaskSEpPbCorrelationsMCYS()
    : AliAnalysisTaskSE(),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fIsAOD(kTRUE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fPIDResponse(0),
      multSelection(0),
      ffilterbit(5),
      fPtMin(0.3),
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
      fNEntries(0),
      inputEvent(0),
      fEvent(0),
      fESD(0),
      mctruth(0),
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
      fHistCentrality_aftercut(0),
      fHistLeadQA(0),
      fHistPIDQA(0),
      fhistmcprim(0),
      fhistmcprimfinal(0),
      fhmcprimvzeta(0),
      fhmcprimpdgcode(0),
      fh2_FMD_acceptance_prim(0),
      fh2_FMD_eta_phi_prim(0),
      fh2_FMD_acceptance(0),
      fh2_FMD_eta_phi(0),
      fhistfmd(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
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
      fHistReconstTrackMix(0)
      {
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
		
		
		for(Int_t i=0;i<4;i++){
		  fhrefetaFMD[i]=0;
		  fhrefphiFMD[i]=0;
		}
		
		
		for(Int_t i=0;i<4;i++){
		  fhistitsdeltaphi[i]=0;
		  fhistitsrefdeltaphi[i]=0;
		  fhistitsrefdeltaphiaftercut[i]=0;
		  fhistitsrefdeltaetaaftercut[i]=0;
		}
		
	  }

AliAnalysisTaskSEpPbCorrelationsMCYS::AliAnalysisTaskSEpPbCorrelationsMCYS(const char *name)
    : AliAnalysisTaskSE(name),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fIsAOD(kTRUE),
      fOnfly(kFALSE),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fPID(kFALSE),
      fCentType("ZNA"),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      fPIDResponse(0),
      multSelection(0),
      ffilterbit(5),
      fPtMin(0.3),
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
      fNEntries(0),
      inputEvent(0),
      fEvent(0),
      fESD(0),
      mctruth(0),
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
      fHistCentrality_aftercut(0),
      fHistLeadQA(0),
      fHistPIDQA(0),
      fhistmcprim(0),
      fhistmcprimfinal(0),
      fhmcprimvzeta(0),
      fhmcprimpdgcode(0),
      fh2_FMD_acceptance_prim(0),
      fh2_FMD_eta_phi_prim(0),
      fh2_FMD_acceptance(0),
      fh2_FMD_eta_phi(0),
      fhistfmd(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
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
      fHistReconstTrackMix(0)
      {
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
		for(Int_t i=0;i<4;i++){
          fhrefetaFMD[i]=0;
          fhrefphiFMD[i]=0;
        }
		for(Int_t i=0;i<4;i++){
		  fhistitsdeltaphi[i]=0;
		  fhistitsrefdeltaphi[i]=0;
		  fhistitsrefdeltaphiaftercut[i]=0;
		  fhistitsrefdeltaetaaftercut[i]=0;
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

  fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
  if (!fPoolMgr)  return;
  fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

  fPoolMgr1 = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
  if (!fPoolMgr1)  return;
  fPoolMgr1->SetTargetValues(fPoolMinNTracks, 0.1, 5);
}

void AliAnalysisTaskSEpPbCorrelationsMCYS::DefineGeneralOutput() {

  fHist_Stat = new TH1F("fHist_Stat", "Stat Histogram", 11, -0.5, 10.5);
  fHist_Stat->GetXaxis()->SetBinLabel(1, "All Events");
  fHist_Stat->GetXaxis()->SetBinLabel(2, "Analyzed Events");
  fHist_Stat->GetXaxis()->SetBinLabel(3, "Trigger OK");
  fHist_Stat->GetXaxis()->SetBinLabel(4, "Vertex OK");
  fHist_Stat->GetXaxis()->SetBinLabel(5, "Centrality OK");
  fHist_Stat->GetXaxis()->SetBinLabel(6, "Pile-Up rejected");
  fHist_Stat->GetXaxis()->SetBinLabel(7, "SPD vetex OK");
  fHist_Stat->GetXaxis()->SetBinLabel(8, "FMD multi cut");
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

  fHistCentrality = new TH1F("fHistCentrality", ";centrality;count", 20, 0, 100);
  fOutputList->Add(fHistCentrality);
  fHistCentrality_aftercut = new TH1F("fHistCentrality_aftercut", ";centrality;count", 20, 0, 100);
  fOutputList->Add(fHistCentrality_aftercut);
  TTree *settingsTree = new TTree("UEAnalysisSettings", "Analysis Settings in UE estimation");
  settingsTree->Branch("fZVertex", &fZVertex, "fZVertex/D");
  settingsTree->Branch("fEtaMax", &fEtaMax, "fEtaMax/D");
  settingsTree->Branch("fPtMin", &fPtMin, "fPtMin/D");
  settingsTree->Branch("fMaxnSigmaTPCTOF", &fMaxnSigmaTPCTOF, "fMaxnSigmaTPCTOF/D");

  // V0 Particle
  settingsTree->Branch("fEtaMinV0", &fEtaMinV0, "fEtaMinV0/D");
  settingsTree->Branch("fdcaDaughtersToPrimVtx", &fdcaDaughtersToPrimVtx, "fdcaDaughtersToPrimVtx/D");
  settingsTree->Branch("fdcaBetweenDaughters", &fdcaBetweenDaughters, "fdcaBetweenDaughters/D");
  settingsTree->Branch("fRadiMin", &fRadiMin, "fRadiMin/D");
  settingsTree->Branch("fRadiMax", &fRadiMax, "fRadiMax/D");
  settingsTree->Branch("fcutcTauK0", &fcutcTauK0, "fcutcTauK0");
  settingsTree->Branch("fcutcTauLam", &fcutcTauLam, "fcutcTauLam");
  settingsTree->Branch("fcosMinK0s", &fcosMinK0s, "fcosMinK0s");
  settingsTree->Branch("fcosMinLambda", &fcosMinLambda, "fcosMinLambda");
  settingsTree->Branch("fMaxnSigmaTPCV0", &fMaxnSigmaTPCV0, "fMaxnSigmaTPCV0");
  // Phi
  settingsTree->Branch("ffilterbit", &ffilterbit, "ffilterbit/I");

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
  Double_t binning_cent[16] = {0.,  1.,  2.,  3.,  4.,  5.,  10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1};

  if(fAnaMode=="TPCV0A"||fAnaMode=="TPCV0C"||fAnaMode=="V0AV0C"){
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

}
void AliAnalysisTaskSEpPbCorrelationsMCYS::DefinedQAHistos() {

  const Int_t ipidBin[5] = {11, 40, 72, 15,20};
  Double_t binning_pt_lead[12] = {0.3, 0.5, 0.75, 1.0, 1.25, 1.5,
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

  fHistLeadQA = new AliTHn("fHistLeadQA", "fHistLeadQA", 1, 5, ipidBin);
  fHistLeadQA->SetBinLimits(0, binning_pt_lead);
  fHistLeadQA->SetBinLimits(1, binning_eta);
  fHistLeadQA->SetBinLimits(2, binning_dphi);
  fHistLeadQA->SetBinLimits(3, binning_cent);
  fHistLeadQA->SetBinLimits(4, -10.,10.);
  fHistLeadQA->SetVarTitle(0, "pt");
  fHistLeadQA->SetVarTitle(1, "eta");
  fHistLeadQA->SetVarTitle(2, "phi");
  fHistLeadQA->SetVarTitle(3, "centrality");
  fHistLeadQA->SetVarTitle(4, "vz");
  fOutputList1->Add(fHistLeadQA);


  const Int_t imcprimbin[5]={12,200,20,11,10};
  Double_t centBins_mcprim[12] = {0.,  5.,  10., 20.,
                            30., 40., 50., 60., 70., 80., 90., 100.0};
    
  Double_t binning_pt_mcprim[13] = {0.,0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
  if(!fDataType){
	
	for(Int_t i=0;i<4;i++){
	fhistitsdeltaphi[i]=new TH1F(Form("fhistitsdeltaphi_%d",i),"fhistitsdeltaphi",100,0,100);
	fOutputList2->Add(fhistitsdeltaphi[i]);
	fhistitsrefdeltaphi[i]=new TH2F(Form("fhistitsrefdeltaphi_%d",i),"fhistitsrefdeltaphi",300,0,3,200,-100,100);
	fOutputList2->Add(fhistitsrefdeltaphi[i]);
	fhistitsrefdeltaphicorr[i]=new TH2F(Form("fhistitsrefdeltaphicorr_%d",i),"fhistitsrefdeltaphicorr",300,0,3,200,-100,100);
	fOutputList2->Add(fhistitsrefdeltaphicorr[i]);

	fhistitsrefdeltaphiaftercut[i]=new TH2F(Form("fhistitsrefdeltaphiaftercut_%d",i),"fhistitsrefdeltaphiaftercut",300,0,3,200,-100,100);
	fOutputList2->Add(fhistitsrefdeltaphiaftercut[i]);
	fhistitsrefdeltaphicorraftercut[i]=new TH2F(Form("fhistitsrefdeltaphicorraftercut_%d",i),"fhistitsrefdeltaphicorraftercut",300,0,3,200,-100,100);
	fOutputList2->Add(fhistitsrefdeltaphicorraftercut[i]);

	fhistitsrefdeltaetaaftercut[i]=new TH2F(Form("fhistitsrefdeltaetaaftercut_%d",i),"fhistitsrefdeltaetaaftercut",300,0,3,200,-0.1,0.1);
	fOutputList2->Add(fhistitsrefdeltaetaaftercut[i]);
	}


    fhistmcprim=new AliTHn("fhistmcprim","fhistmcprim",1,5,imcprimbin);
    fhistmcprim->SetBinLimits(0,binning_pt_mcprim);
    fhistmcprim->SetBinLimits(1,-4.,6.);
    fhistmcprim->SetBinLimits(2,0.,2*TMath::Pi());
    fhistmcprim->SetBinLimits(3,centBins_mcprim);
    fhistmcprim->SetBinLimits(4,-10.,10.);
    fhistmcprim->SetVarTitle(0,"pt");
    fhistmcprim->SetVarTitle(1,"eta");
    fhistmcprim->SetVarTitle(2,"phi");
    fhistmcprim->SetVarTitle(3,"centrality");
    fhistmcprim->SetVarTitle(4,"vz");
    fOutputList2->Add(fhistmcprim);

    fhistmcprimfinal=new AliTHn("fhistmcprimfinal","fhistmcprimfinal",1,5,imcprimbin);
    fhistmcprimfinal->SetBinLimits(0,binning_pt_mcprim);
    fhistmcprimfinal->SetBinLimits(1,-4.,6.);
    fhistmcprimfinal->SetBinLimits(2,0.,2*TMath::Pi());
    fhistmcprimfinal->SetBinLimits(3,centBins_mcprim);
    fhistmcprimfinal->SetBinLimits(4,-10.,10.);
    fhistmcprimfinal->SetVarTitle(0,"pt");
    fhistmcprimfinal->SetVarTitle(1,"eta");
    fhistmcprimfinal->SetVarTitle(2,"phi");
    fhistmcprimfinal->SetVarTitle(3,"centrality");
    fhistmcprimfinal->SetVarTitle(4,"vz");
    fOutputList2->Add(fhistmcprimfinal);

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

  const Int_t ifmdbin[4]={200,20,11,20};
  if(fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC" || fAnaMode=="FMDFMD"){
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

    fh2_FMD_acceptance=new TH2D("fh2_FMD_acceptance","fh2_FMD_acceptance",200,-4,6,200,-10,10);
    fOutputList2->Add(fh2_FMD_acceptance);
    fh2_FMD_eta_phi=new TH2D("fh2_FMD_eta_phi","fh2_FMD_eta_phi",200,-4,6,20,0,2*TMath::Pi());
    fOutputList2->Add(fh2_FMD_eta_phi);

    fhistfmd=new AliTHn("fhistfmd","fhistfmd",1,4,ifmdbin);
    fhistfmd->SetBinLimits(0,-4.,6.);
    fhistfmd->SetBinLimits(1,0.,2*TMath::Pi());
    fhistfmd->SetBinLimits(2,centBins_mcprim);
    fhistfmd->SetBinLimits(3,-10.,10.);
    fhistfmd->SetVarTitle(0,"eta");
    fhistfmd->SetVarTitle(1,"phi");
    fhistfmd->SetVarTitle(2,"centrality");
    fhistfmd->SetVarTitle(3,"zvertex");
    fOutputList2->Add(fhistfmd);
  }


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

  Double_t binning_pt_assoc[12] = {0.3, 0.5, 0.75, 1.0, 1.25, 1.5,
                                   2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
  Double_t binning_pt_lead[13] = {0.,0.3, 0.5, 0.75, 1.0, 1.25, 1.5,
                                  2.0, 2.5, 3.0,  3.5, 4.0,  8.0};
  Double_t binning_cent[16] = {0.,  1.,  2.,  3.,  4.,  5.,  10., 20.,
                               30., 40., 50., 60., 70., 80., 90., 100.1};
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
  const Int_t nEvtVars = 2;
  const Int_t iEvtBin[2] = {12, 15};

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

  if(fAnaMode=="V0AV0C"){
    const Int_t nEvtVarsV0Leading=3;
    const Int_t iEvtBinV0Leading[3]={15,10,8};
    fHistTriggerTrack= new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsV0Leading, iEvtBinV0Leading);
    fHistTriggerTrack->SetBinLimits(0,binning_cent);
    fHistTriggerTrack->SetBinLimits(1,binning_eta_vzero);
    fHistTriggerTrack->SetBinLimits(2,binning_phi_vzero);
    fHistTriggerTrack->SetVarTitle(0,"centrality");
    fHistTriggerTrack->SetVarTitle(1,"eta");
    fHistTriggerTrack->SetVarTitle(2,"phi");
    fHistTriggerTrackMix= new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVarsV0Leading, iEvtBinV0Leading);
    fHistTriggerTrackMix->SetBinLimits(0,binning_cent);
    fHistTriggerTrackMix->SetBinLimits(1,binning_eta_vzero);
    fHistTriggerTrackMix->SetBinLimits(2,binning_phi_vzero);
    fHistTriggerTrackMix->SetVarTitle(0,"centrality");
    fHistTriggerTrackMix->SetVarTitle(1,"eta");
    fHistTriggerTrackMix->SetVarTitle(2,"phi");
  }else if(fAnaMode=="FMDFMD"){
    const Int_t nEvtVarsV0Leading=2;
    const Int_t iEvtBinV0Leading[2]={15,32};
    fHistTriggerTrack= new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsV0Leading, iEvtBinV0Leading);
    fHistTriggerTrack->SetBinLimits(0,binning_cent);
    fHistTriggerTrack->SetBinLimits(1,binning_etafmd);
    fHistTriggerTrack->SetVarTitle(0,"centrality");
    fHistTriggerTrack->SetVarTitle(1,"eta");
    fHistTriggerTrackMix= new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVarsV0Leading, iEvtBinV0Leading);
    fHistTriggerTrackMix->SetBinLimits(0,binning_cent);
    fHistTriggerTrackMix->SetBinLimits(1,binning_etafmd);
    fHistTriggerTrackMix->SetVarTitle(0,"centrality");
    fHistTriggerTrackMix->SetVarTitle(1,"eta");
  }else{
    fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVars, iEvtBin);
    fHistTriggerTrack->SetBinLimits(0, binning_pt_lead);
    fHistTriggerTrack->SetBinLimits(1, binning_cent);
    fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
    fHistTriggerTrack->SetVarTitle(1, "centrality");

    fHistTriggerTrackMix = new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVars, iEvtBin);
    fHistTriggerTrackMix->SetBinLimits(0, binning_pt_lead);
    fHistTriggerTrackMix->SetBinLimits(1, binning_cent);
    fHistTriggerTrackMix->SetVarTitle(0, "leading p_{T} GeV/c");
    fHistTriggerTrackMix->SetVarTitle(1, "centrality");
  }
  fOutputList1->Add(fHistTriggerTrack);
  fOutputList1->Add(fHistTriggerTrackMix);

  const Int_t nTrackVars = 5;
  const Int_t iTrackBin[5] = {48, 11, 12, 15, 72};

  //////////////////////////////////////////
  //Containers two particle correlation
  //////////////////////////////////////////
  Int_t nCFSteps = 1;
  if(fasso=="hadron")    nCFSteps=1;
  else if (fasso == "V0" || fasso == "Phi")    nCFSteps = 7;
  else  if (fasso == "Cascade")    nCFSteps = 6;
  else    if(fasso=="PID")    nCFSteps=3;
  Double_t binning_dphi_vzero[9]={-1.178097,-0.392699,0.392699,1.178097,1.963495,2.748893,3.534291,4.319689,5.105088};
  if(fAnaMode=="TPCTPC"){
    fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars, iTrackBin);
    fHistReconstTrack->SetBinLimits(0, binning_deta);
    fHistReconstTrack->SetBinLimits(1, binning_pt_assoc);
    fHistReconstTrack->SetBinLimits(2, binning_pt_lead);
    fHistReconstTrack->SetBinLimits(3, binning_cent);
    fHistReconstTrack->SetBinLimits(4, binning_dphi);
    fHistReconstTrack->SetVarTitle(0, "#Delta#eta");
    fHistReconstTrack->SetVarTitle(1, "p_{T} GeV/c");
    fHistReconstTrack->SetVarTitle(2, "leading p_{T} GeV/c");
    fHistReconstTrack->SetVarTitle(3, "centrality");
    fHistReconstTrack->SetVarTitle(4, "#Delta#phi");
    fHistReconstTrackMix =  new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars, iTrackBin);
    fHistReconstTrackMix->SetBinLimits(0, binning_deta);
    fHistReconstTrackMix->SetBinLimits(1, binning_pt_assoc);
    fHistReconstTrackMix->SetBinLimits(2, binning_pt_lead);
    fHistReconstTrackMix->SetBinLimits(3, binning_cent);
    fHistReconstTrackMix->SetBinLimits(4, binning_dphi);
    fHistReconstTrackMix->SetVarTitle(0, "#Delta#eta");
    fHistReconstTrackMix->SetVarTitle(1, "p_{T} GeV/c");
    fHistReconstTrackMix->SetVarTitle(2, "leading p_{T} GeV/c");
    fHistReconstTrackMix->SetVarTitle(3, "centrality");
    fHistReconstTrackMix->SetVarTitle(4, "#Delta#phi");
  }else if(fAnaMode=="TPCV0A"||fAnaMode=="TPCV0C"){
    const Int_t iTrackBin_VZEROA[5]={66,12,10,15,72};
    const Int_t iTrackBin_VZEROC[5]={62,12,10,15,72};
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
  }else if (fAnaMode=="TPCFMD"){
/*
      Double_t binning_detaFMDTPC[97]={-5.7,-5.65,-5.6,-5.55,-5.5,-5.45,-5.4,-5.35,-5.3,-5.25,-5.2,-5.15,-5.1,-5.05,-5.0,-4.95, -4.9,-4.85, -4.8, -4.75, -4.7, -4.65,-4.6,-4.55,-4.5,-4.45,-4.4,-4.35,-4.3,
      -4.25,-4.2,-4.15,-4.1,-4.05,-4.0,-3.95,-3.9,-3.85,-3.8,-3.75,-3.7,-3.65,-3.6,-3.55,-3.5,-3.45,-3.4,-3.35,-3.3,-3.25,-3.2,-3.15,-3.1,-3.05,-3.0,-2.95,-2.9,-2.85,-2.8,-2.75,-2.7,-2.65,
      -2.6,-2.55,-2.5,-2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,
      -0.95,-0.9};
      */
      Double_t binning_detaFMDTPC[49]={
      -5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.0,
      -4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,
      -3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,
      -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,
      -1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,
      -0.9};

	  const Int_t iTrackBin_tpcfmd[5]={48,12,32,15,72};//deta,pt,etafmd,centrality,dphi
	  fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 5, iTrackBin_tpcfmd);
	  fHistReconstTrack->SetBinLimits(0,binning_detaFMDTPC);
	  fHistReconstTrack->SetBinLimits(1,binning_pt_lead);
	  fHistReconstTrack->SetBinLimits(2,binning_etafmd);
	  fHistReconstTrack->SetBinLimits(3,binning_cent);
	  fHistReconstTrack->SetBinLimits(4,binning_dphi);
	  fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
	  fHistReconstTrack->SetVarTitle(1,"p_{T} GeV/c");
	  fHistReconstTrack->SetVarTitle(2,"FMD Eta");
	  fHistReconstTrack->SetVarTitle(3,"centrality");
	  fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
	  fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 5,iTrackBin_tpcfmd);
	  fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDTPC);
	  fHistReconstTrackMix->SetBinLimits(1,binning_pt_lead);
	  fHistReconstTrackMix->SetBinLimits(2,binning_etafmd);
	  fHistReconstTrackMix->SetBinLimits(3,binning_cent);
	  fHistReconstTrackMix->SetBinLimits(4,binning_dphi);
	  fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
	  fHistReconstTrackMix->SetVarTitle(1,"p_{T} GeV/c");
	  fHistReconstTrackMix->SetVarTitle(2,"FMD Eta");
	  fHistReconstTrackMix->SetVarTitle(3,"centrality");
	  fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
  }else if (fAnaMode=="TPCFMDC"){
    /*
    Double_t binning_detaFMDCTPC[67]={0.9,0.95,
    1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,
    2.0,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,
    3.0,3.05,3.1,3.15,3.2,3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,
    4.0,4.05,4.1,4.15,4.2};
    */
    Double_t binning_detaFMDCTPC[34]={
      0.9,
      1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
      2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
      3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
      4.0,4.1,4.2};

    const Int_t iTrackBin_tpcfmdc[5]={33,12,17,15,72};
    fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 5, iTrackBin_tpcfmdc);
    fHistReconstTrack->SetBinLimits(0,binning_detaFMDCTPC);
    fHistReconstTrack->SetBinLimits(1,binning_pt_lead);
    fHistReconstTrack->SetBinLimits(2,binning_etafmdc);
    fHistReconstTrack->SetBinLimits(3,binning_cent);
    fHistReconstTrack->SetBinLimits(4,binning_dphi);
    fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrack->SetVarTitle(1,"p_{T} GeV/c");
    fHistReconstTrack->SetVarTitle(2,"FMD Eta");
    fHistReconstTrack->SetVarTitle(3,"centrality");
    fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
    fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 5,iTrackBin_tpcfmdc);
    fHistReconstTrackMix->SetBinLimits(0,binning_detaFMDCTPC);
    fHistReconstTrackMix->SetBinLimits(1,binning_pt_lead);
    fHistReconstTrackMix->SetBinLimits(2,binning_etafmdc);
    fHistReconstTrackMix->SetBinLimits(3,binning_cent);
    fHistReconstTrackMix->SetBinLimits(4,binning_dphi);
    fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrackMix->SetVarTitle(1,"p_{T} GeV/c");
    fHistReconstTrackMix->SetVarTitle(2,"FMD Eta");
    fHistReconstTrackMix->SetVarTitle(3,"centrality");
    fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
  }else if(fAnaMode=="FMDFMD"){
    const Int_t nTrackVars_fmdfmd = 5;
    const Int_t iTrackBin_fmdfmd[5]={50,17,32,15,20};
    //const Double_t binning_detafmdfmd={};
    fHistReconstTrack= new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
    //fHistReconstTrack->SetBinLimits(1,binning_detafmdfmd);
    fHistReconstTrack->SetBinLimits(0,3.425,8.425);
    fHistReconstTrack->SetBinLimits(1,binning_etafmdc);
    fHistReconstTrack->SetBinLimits(2,binning_etafmd);
    fHistReconstTrack->SetBinLimits(3,binning_cent);
    fHistReconstTrack->SetBinLimits(4,-0.55*TMath::Pi(),1.45*TMath::Pi());
    fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrack->SetVarTitle(1,"FMD(Asso) Eta");
    fHistReconstTrack->SetVarTitle(2,"FMD(Trigger) Eta");
    fHistReconstTrack->SetVarTitle(3,"centrality");
    fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
    fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars_fmdfmd,iTrackBin_fmdfmd);
    fHistReconstTrackMix->SetBinLimits(0,3.425,8.425);
    fHistReconstTrackMix->SetBinLimits(1,binning_etafmdc);
    fHistReconstTrackMix->SetBinLimits(2,binning_etafmd);
    fHistReconstTrackMix->SetBinLimits(3,binning_cent);
    fHistReconstTrackMix->SetBinLimits(4,-0.55*TMath::Pi(),1.45*TMath::Pi());
    fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrackMix->SetVarTitle(1,"FMD(Asso) Eta");
    fHistReconstTrackMix->SetVarTitle(2,"FMD(Trigger) Eta");
    fHistReconstTrackMix->SetVarTitle(3,"centrality");
    fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
  }else if(fAnaMode=="V0AV0C"){
    const Int_t nTrackVars_v0av0c = 4;
    const Int_t iTrackBin_VZEROAVZEROC[4]={10,10,15,8};
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
  }
  
  fOutputList1->Add(fHistReconstTrack);
  fOutputList1->Add(fHistReconstTrackMix);
}
void AliAnalysisTaskSEpPbCorrelationsMCYS::UserExec(Option_t *) {
  AliAnalysisManager *mgr        = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inEvMain = (AliInputEventHandler *)(mgr->GetInputEventHandler());
  if (!inEvMain)    return;
  fPIDResponse = inEvMain->GetPIDResponse();
  if (!fPIDResponse)    return;
  if(!fDataType){
    if(!fIsAOD){
      mctruth = (AliMCEventHandler*)(mgr->GetMCtruthEventHandler());
      if(!mctruth)       return;
      if (!mctruth->InitOk()) return;
      if (!mctruth->TreeK()) return;
      mctruth->CreateLabelMap();
      mcEvent=mctruth->MCEvent();//AliMCEvent
      mcEvent->InitEvent();
      mcEvent->PreReadAll();
    }else{
      //      mcEvent=((AliAODInputHandler*)inEvMain)->MCEvent();
    }
  }
  if(fIsAOD){
    fEvent = dynamic_cast<AliAODEvent *>(inEvMain->GetEvent());
    if (!fEvent) {
      AliWarning("ERROR: fEvent not available \n");
      return;
    }
  }else{
    fESD=dynamic_cast<AliESDEvent*>(inEvMain->GetEvent());
    if (!fESD) {
      AliWarning("ERROR: fESD not available \n");
      return;
    }
  }


  inputEvent=fEvent;
  if(!inputEvent) inputEvent=fESD;

  fHist_Stat->Fill(0);
  multSelection =    (AliMultSelection *)inputEvent->FindListObject("MultSelection");

  if(!multSelection) return;
  fHist_Stat->Fill(1);

  UInt_t maskIsSelected = inEvMain->IsEventSelected();
  Bool_t isSelected     = kFALSE;
  isSelected = ((maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7);//Both for data and MC
  if (!isSelected)  return;
  fHist_Stat->Fill(2);

  // Primary Vllllertex
  lPrimaryBestVtx = inputEvent->GetPrimaryVertex();
  if (!lPrimaryBestVtx)    return;
  Int_t nTracksPrim = lPrimaryBestVtx->GetNContributors();
  if (nTracksPrim < 1)    return;

  if ((TMath::Abs(lPrimaryBestVtx->GetZ())) >= fZVertex)    return;
  tPrimaryVtxPosition[0] = lPrimaryBestVtx->GetX();
  tPrimaryVtxPosition[1] = lPrimaryBestVtx->GetY();
  tPrimaryVtxPosition[2] = lPrimaryBestVtx->GetZ();
  fHistzvertex->Fill(tPrimaryVtxPosition[2]);
  fPrimaryZVtx = lPrimaryBestVtx->GetZ();
  fHist_Stat->Fill(3);

  bSign = 0.;
  bSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

  // Multiplicity Object
  if(fcollisiontype=="pPb"){
    if(frun2){
      //AliMultSelection *multSelection =    (AliMultSelection *)fEvent->FindListObject("MultSelection");
      lCentrality = multSelection->GetMultiplicityPercentile(fCentType);
      Int_t qual = multSelection->GetEvSelCode();
      if (qual == 199)  lCentrality = -999;
    } else{
      AliCentrality *centobj = 0;
      //centobj = fEvent->GetCentrality();
      centobj = inputEvent->GetCentrality();
      lCentrality = centobj->GetCentralityPercentile(fCentType);
      if(!centobj) lCentrality=-1.;
    }

    if (lCentrality < 0. || lCentrality > 100. - 0.0000001)   return;
    Double_t *CentBins = fCentBins;
    poolmin = CentBins[0];
    poolmax = CentBins[fNCentBins];
    fHist_Stat->Fill(4);
  }else{
    //    AliMultSelection *multSelection =    (AliMultSelection *)fEvent->FindListObject("MultSelection");
    // if (!multSelection) {
	///  AliWarning ("AliMultSelection could not be found in the aod event list of objects");
    //  }
  }

  AliAnalysisUtils *fUtils=new AliAnalysisUtils;
  //if(fcollisiontype=="pp") if(fUtils->IsPileUpSPD(fEvent)) return;
  if(fUtils->IsPileUpMV(inputEvent)) return;
  //if(fUtils->IsPileUpSPD(fEvent)) return;
  //  if(fEvent->IsPileupFromSPD(5,0.8,3.,2.,5.)) return;
  fHist_Stat->Fill(5);
  // SPD vertex selection
  const AliVVertex* vtxSPD = dynamic_cast<const AliVVertex*>(inputEvent->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return;
  fHist_Stat->Fill(6);
  //if(fNEntries==2) return;
  if(fcollisiontype=="pPb") fHistCentrality->Fill(lCentrality);


  MakeAna();

  
  AliEventplane*evtPlane=fEvent->GetEventplane();
  Double_t qx=0;
  Double_t qy=0;
  Double_t evtPlanePhi=-999;
  if(evtPlane){
    evtPlanePhi=evtPlane->CalculateVZEROEventPlane(fEvent,10,2,qx,qy);
    // cout<<evtPlanePhi<<endl;
  }
  
  
  fNEntries++;
  
  PostData(1, fOutputList);
  PostData(2, fOutputList1);
  PostData(3, fOutputList2);
}

void AliAnalysisTaskSEpPbCorrelationsMCYS::Terminate(Option_t *) {
  //  AliInfo(Form("Number of Correlation
  Printf("Entries======================%f",fNEntries);
  if (fPoolMgr)    delete fPoolMgr;   // PoolMgr->ClearPools();
  if (fPoolMgr1)    delete fPoolMgr1; // fPoolMgr1->ClearPools();
}

void AliAnalysisTaskSEpPbCorrelationsMCYS::MakeAna() {
  TObjArray *selectedTracksLeading = new TObjArray;
  selectedTracksLeading->SetOwner(kTRUE);
  TObjArray *selectedTracksAssociated = new TObjArray;
  selectedTracksAssociated->SetOwner(kTRUE);

  TObjArray *selectedTrackV0A=new TObjArray;
  selectedTrackV0A->SetOwner(kTRUE);
  TObjArray *selectedTrackV0C=new TObjArray;
  selectedTrackV0C->SetOwner(kTRUE);

  TObjArray *selectedFMDArray1=new TObjArray;
  selectedFMDArray1->SetOwner(kTRUE);
  TObjArray *selectedFMDArray2=new TObjArray;
  selectedFMDArray2->SetOwner(kTRUE);
  
  TObjArray* selectedTracksMC=new TObjArray;
  selectedTracksMC->SetOwner(kTRUE);

  // Leading Particle

  fvzero = fEvent->GetVZEROData();
  if(fAnaMode=="TPCV0A"||fAnaMode=="TPCV0C"||fAnaMode=="V0AV0C"){
  Double_t eta_min;
  Double_t eta_max;
  Double_t eta_ave;
  Double_t phi_vzero;
  Double_t mult_vzero;
  Double_t vzeroqa[3];
  Double_t mult_vzero_eq;
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
	//  if(fAnaMode=="TPCV0A" || fAnaMode=="TPCV0C" || fAnaMode=="V0AV0C"){
	//    if(imod>31) selectedTrackV0A->Add(new Aliassociatedvzeroys(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
	if(imod>31) selectedTrackV0A->Add(new AliAssociatedTrackYS(0,eta_ave,phi_vzero,0,0,0,0,0,mult_vzero_eq));
	// if(imod<32) selectedTrackV0C->Add(new AliAssociatedVZEROYS(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
   	if(imod<32) selectedTrackV0C->Add(new AliAssociatedTrackYS(0,eta_ave,phi_vzero,0,0,0,0,0,mult_vzero_eq));
    //}
  }
}


if(fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC" || fAnaMode=="FMDFMD"){
   selectedFMDArray1=GetFMDhitsYS(kTRUE);//A-side
    selectedFMDArray2=GetFMDhitsYS(kFALSE);//C-side
    Float_t nV0A_hits = fvzero->GetMTotV0A();
    Float_t nV0C_hits = fvzero->GetMTotV0C();
    Float_t nFMD_fwd_hits=0;
    Float_t nFMD_bwd_hits=0;

    for(Int_t i=0;i<selectedFMDArray1->GetEntriesFast();i++){
	  //      AliAssociatedVZEROYS* FMDA=(AliAssociatedVZEROYS*)selectedFMDArray1->At(i);
      AliAssociatedTrackYS* FMDA=(AliAssociatedTrackYS*)selectedFMDArray1->At(i);
      if(!FMDA) continue;
      nFMD_fwd_hits+=FMDA->Multiplicity();
    }
    for(Int_t i=0;i<selectedFMDArray2->GetEntriesFast();i++){
      AliAssociatedTrackYS* FMDC=(AliAssociatedTrackYS*)selectedFMDArray2->At(i);
      if(!FMDC) continue;
      nFMD_bwd_hits+=FMDC->Multiplicity();
    }


    fFMDV0->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
    fFMDV0A->Fill(nFMD_fwd_hits, nV0A_hits);
    fFMDV0C->Fill(nFMD_bwd_hits, nV0C_hits);

    if (nV0A_hits + nV0C_hits < 1.5*(nFMD_fwd_hits + nFMD_bwd_hits) - 20) return; //events cuts

    fFMDV0_post->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
    fFMDV0A_post->Fill(nFMD_fwd_hits, nV0A_hits);
    fFMDV0C_post->Fill(nFMD_bwd_hits, nV0C_hits);
  }



  fHist_Stat->Fill(7);
  if(fcollisiontype=="pPb") fHistCentrality_aftercut->Fill(lCentrality);


  //finish event selection, everything should be filled after this line
  if(fIsAOD){
  selectedTracksLeading = GetAcceptedTracksLeading(fEvent,kTRUE);
  if (fasso == "Phi")    selectedTracksAssociated = GetAcceptedTracksAssociated(fEvent);
  if (fasso == "V0")    selectedTracksAssociated = GetAcceptedV0Tracks(fEvent);
  if (fasso == "PID")    selectedTracksAssociated = GetAcceptedTracksPID(fEvent);
  if (fasso == "Cascade")    selectedTracksAssociated = GetAcceptedCascadeTracks(fEvent);
  if (fAnaMode=="TPCTPC" && fasso == "hadron")    selectedTracksAssociated = GetAcceptedTracksLeading(fEvent,kFALSE);
  }

  // AliAODTracklets *tracklets = ((AliAODEvent*)fEvent)->GetTracklets();
  AliVMultiplicity *tracklets = ((AliAODEvent*)fEvent)->GetTracklets();
  if (!tracklets) return;
  Int_t nTracklets = tracklets->GetNumberOfTracklets();
  
  for (Int_t i = 0; i < nTracklets; i++) {
	Double_t dphi  = tracklets->GetDeltaPhi(i);
	/*
	  if (TMath::Abs(dphi) * 1000 > 5) {
	  continue;
	  }
	*/
	
	Double_t theta = tracklets->GetTheta(i);
	Double_t etaits   = -TMath::Log(TMath::Tan(theta/2));
	// Drop everything outside of -1.7 < eta 1.7 to avoid overlas with the FMD
	/*
	if (etaits < -1.7 || etaits > 1.7) {
	  continue;
	}
	*/
	Double_t phiits   = tracklets->GetPhi(i);
	Double_t phiitscorr = phiits+ tracklets->GetDeltaPhi(i)*37./39.;//correction dphi*39./34. (Dphi in rad)	
	if (phiitscorr<0) phiitscorr+=TMath::TwoPi();
	if (phiitscorr>TMath::TwoPi()) phiitscorr-=TMath::TwoPi();

	//	   fh2_ITS_acceptance->Fill(etaits,tPrimaryVtxPosition[2]);
	Double_t itsqa[3]={etaits,phiits,lCentrality};
	// fhistits->Fil0l(itsqa);
	
   	Int_t label1 = tracklets->GetLabel(i,0);
	Int_t label2 = tracklets->GetLabel(i,1);
	//	if(label1!=label2)cout<<i<<" "<<label1<<" "<<label2<<endl;

	AliVParticle* particle = fMCEvent->GetTrack(label1);
	if (!particle) 	  continue;
   
	Short_t charge  = particle->Charge();
	Float_t ptMC    = particle->Pt();
	Float_t etaMC   = particle->Eta();
	Float_t phiMC   = particle->Phi();
	Float_t pdg     = particle->PdgCode();
	Bool_t primary  = particle->InheritsFrom("AliAODMCParticle") ? ((AliAODMCParticle*) particle)->IsPhysicalPrimary() : fMCEvent->IsPhysicalPrimary(label1);

	Bool_t itsdphicut=TMath::Abs(dphi) * 1000 <5.;
		
   	fhistitsdeltaphi[0]->Fill(abs(dphi*1000));
	fhistitsrefdeltaphi[0]->Fill(ptMC,1000.*(phiits-phiMC));
	fhistitsrefdeltaphicorr[0]->Fill(ptMC,1000.*(phiitscorr-phiMC));
	if(itsdphicut) 	{
	  fhistitsrefdeltaphiaftercut[0]->Fill(ptMC,1000.*(phiits-phiMC));
	  fhistitsrefdeltaphicorraftercut[0]->Fill(ptMC,1000.*(phiitscorr-phiMC));
	  fhistitsrefdeltaetaaftercut[0]->Fill(ptMC,etaits-etaMC);
	}
	
	if(label1!=label2) {
	  fhistitsdeltaphi[1]->Fill(abs(dphi*1000));
	  fhistitsrefdeltaphi[1]->Fill(ptMC,1000.*(phiits-phiMC));
	  fhistitsrefdeltaphicorr[2]->Fill(ptMC,1000.*(phiitscorr-phiMC));
	  if(itsdphicut) {
		fhistitsrefdeltaphiaftercut[1]->Fill(ptMC,1000.*(phiits-phiMC));
		fhistitsrefdeltaphicorraftercut[1]->Fill(ptMC,1000.*(phiitscorr-phiMC));
		fhistitsrefdeltaetaaftercut[1]->Fill(ptMC,etaits-etaMC);
	  }
	}
	
	if(label1==label2){
	  if(primary){
		fhistitsdeltaphi[2]->Fill(abs(dphi*1000));
		fhistitsrefdeltaphi[2]->Fill(ptMC,1000.*(phiits-phiMC));
		fhistitsrefdeltaphicorr[2]->Fill(ptMC,1000.*(phiitscorr-phiMC));
		if(itsdphicut) {
		  fhistitsrefdeltaphiaftercut[2]->Fill(ptMC,1000.*(phiits-phiMC));
		  fhistitsrefdeltaphicorraftercut[2]->Fill(ptMC,1000.*(phiitscorr-phiMC));
		}
		}
	  if(!primary){
		fhistitsrefdeltaphi[3]->Fill(ptMC,1000.*(phiits-phiMC));
		fhistitsdeltaphi[3]->Fill(abs(dphi*1000));
		fhistitsrefdeltaphicorr[3]->Fill(ptMC,1000.*(phiitscorr-phiMC));
		if(itsdphicut) 	{
		  fhistitsrefdeltaphiaftercut[3]->Fill(ptMC,1000.*(phiits-phiMC));
		  fhistitsrefdeltaphicorraftercut[3]->Fill(ptMC,1000.*(phiitscorr-phiMC));
		}
	  }
	}
	//selectedTracksLeading->Add(new AliAssociatedTrackYS(-999, etaits, phiits, -999, 0, -999,-999, 0, 1));
  }
  //MC
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
  if(!fDataType){
    if(fIsAOD){
      AliAODMCHeader* aodMCheader=(AliAODMCHeader*)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
      TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
      if(!mcArray){
        Printf("No MC particle branch found");
        return;
      }

      Int_t nMCAllTracks = mcArray->GetEntriesFast();
	  Int_t hoge1=0;
	  for (Int_t i = 0; i < nMCAllTracks; i++)
      {
        AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
        if (!mcTrack) {
          Error("ReadEventAODMC", "Could not receive particle %d", i);
          continue;
        }
        
        TrIsPrim=mcTrack->IsPhysicalPrimary();
        TrIsSecondMate=mcTrack->IsSecondaryFromMaterial();
        TrIsSecondWeak=mcTrack->IsSecondaryFromWeakDecay();
        TrIsOthers=!TrIsPrim && !TrIsSecondMate && !TrIsSecondWeak;

		/*
        while(!mcTrack->IsPhysicalPrimary()){
         if(mcTrack->GetMother()<0) break;
         mcTrack=(AliAODMCParticle*)mcArray->At(((AliAODMCParticle * )mcTrack)->GetMother());
         if(!mcTrack) break;
        }
		*/
        
		mcTrackEta = mcTrack->Eta();
        mcTrackPt  = mcTrack->Pt();
        mcTrackPhi = mcTrack->Phi();
        TrCharge=mcTrack->Charge()!=0;
		pdgcode=TMath::Abs(mcTrack->PdgCode());
        if(!TrCharge)        continue;
        
        Int_t daughterlabel1=mcTrack->GetFirstDaughter();
        if(daughterlabel1>0){
          AliAODMCParticle*d0=(AliAODMCParticle*)mcArray->At(daughterlabel1);
        }
        Int_t motherlabel1 = mcTrack->GetMother();
        if(motherlabel1>0){
            AliAODMCParticle*m0=(AliAODMCParticle*)mcArray->At(motherlabel1);
        }
        conmcprim[0]=mcTrackPt;
        conmcprim[1]=mcTrackEta;
        conmcprim[2]=mcTrackPhi;
        conmcprim[3]=lCentrality;
        conmcprim[4]=fPrimaryZVtx;
	 	
		if(TrIsPrim){
          fhistmcprim->Fill(conmcprim,0);//primay charged partilce distribution(no mother particle)
          fhmcprimvzeta->Fill(mcTrackEta,mcTrack->Zv());
          fhmcprimpdgcode->Fill(pdgcode);
          fh2_FMD_eta_phi_prim->Fill(mcTrackEta,mcTrackPhi);
		  //		  fhistmcprimfinal->Fill(conmcprim,0);
		  selectedTracksMC->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),-999,-999,-999,0, 1));        
		}
        
      
      }
	  //	  cout<<nMCAllTracks<<" "<<hoge1<<endl;
        //  cout<<fNEntries<<" "<<"number of primary"<<hoge1<<endl;
    }else{
      /*
      /////
      int np=mcEvent->GetNumberOfTracks();
      int nprim=mcEvent->GetNumberOfPrimaries();
      int j=0;
      int error=0;
      for(int i=0;i<np;i++){
	if(mctruth->IsParticleSelected(i)){
	  
	  if(mctruth->GetNewLabel(i)!=j){
          error++;
	  }
	  
      j++;
    }
  }
//  cout<<"# of IsSelected Particles="<<j<<" : number of label error="<<error<<" nprim="<<nprim<<" np="<<np<<endl;

      Int_t hoge2=0;
      for (Int_t ip = 0; ip < np; ip++){
	AliMCParticle* mcpart = (AliMCParticle*) mcEvent->GetTrack(ip);
	TParticle* part = mcpart->Particle();
	Float_t xv = part->Vx();
	Float_t yv = part->Vy();
	Float_t zv = part->Vz();
	Float_t rv = TMath::Sqrt(xv * xv + yv * yv);

	Bool_t write = kFALSE;
	if (ip < nprim) {
	  // Select the primary event
	  
	  write = kTRUE;
	}else if (part->GetUniqueID() == kPDecay) {
	  // Particles from decay
	  // Check that the decay chain ends at a primary particle
	  AliMCParticle* mother = mcpart;
	  Int_t imo = mcpart->GetMother();
	  while((imo >= nprim) && (mother->Particle()->GetUniqueID() == kPDecay)) {
	    mother =  (AliMCParticle*) mcEvent->GetTrack(imo);
	    imo =  mother->GetMother();
	  }
	  // Select according to pseudorapidity and production point of primary ancestor
	  if (imo < nprim){write = kTRUE;
	    
	  }
	  // if(!Select(((AliMCParticle*) mcE->GetTrack(imo))->Particle(), rv, zv))write = kFALSE; // selection on eta and \
	  phi of the mother
	} else if (part->GetUniqueID() == kPPair) {
	  // Now look for pair production
	  Int_t imo = mcpart->GetMother();
	  if (imo < nprim) {
	    // Select, if the gamma is a primary
	    write = kTRUE;
	    
	  } else {
	    // Check if the gamma comes from the decay chain of a primary particle
	    AliMCParticle* mother =  (AliMCParticle*) mcEvent->GetTrack(imo);
	    imo = mother->GetMother();
	    while((imo >= nprim) && (mother->Particle()->GetUniqueID() == kPDecay)) {
	      mother =   (AliMCParticle*) mcEvent->GetTrack(imo);
	      imo =  mother->GetMother();
	    }
	    // Select according to pseudorapidity and production point
        //if (imo < nprim && Select(mother->Particle(), rv, zv))
	    if (imo < nprim){
            hoge2++;
          Float_t eta=((TParticle*)mother->Particle())->Eta();
          Bool_t select=kFALSE;
          if (TMath::Abs(eta) < 2.5 && rv < 170) select=kTRUE;
          if(eta > -4.2 && eta < -2.3 && zv > -500.) select=kTRUE;

          if(select)write = kTRUE;
        }
      }
    }
    if (write) {
     //if(mcH)mcH->SelectParticle(ip);
    // hoge2++;
   }
 }
      */
 //cout<<fNEntries<<" "<<"number of setparticles=="<<hoge2<<endl;
      //  AliStack*rstack=mcEvent->Stack();

      //cout<<fNEntries<<" "<<mcEvent->GetNumberOfTracks()<<" "<<rstack->GetNtrack()<<endl;
  //  cout<<fNEntries<<" "<<mcEvent->GetNumberOfTracks()<<" "<<mcEvent->GetNumberOfPrimaries()<<" "<<rstack->GetNtrack()<<" "<<((TTree*)(rstack->TreeK()))->GetEntries()<<endl;
      /*
    Double_t  nTracks_vzero=0;
    Int_t nTrCount=0;
    Int_t nTrRefs=0;
    Int_t nV0A = 0;
    Int_t nV0C = 0;
    Double_t conmcprim[4];
    Int_t hoge=0;
    for(Int_t iTracks=0;iTracks<mcEvent->GetNumberOfTracks();iTracks++){
        AliMCParticle* track = (AliMCParticle*)mcEvent->GetTrack(iTracks);
        TParticle *rParticle=mcEvent->Particle(iTracks);
      if (!track) {
        Error("ReadEventMC", "Could not receive particle %d", iTracks);
        continue;
      }
      pdgcode = TMath::Abs(track->PdgCode());
      TrIsPrim=mcEvent->IsPhysicalPrimary(iTracks);
      TrIsSecondMate=mcEvent->IsSecondaryFromMaterial(iTracks);
      TrIsSecondWeak=mcEvent->IsSecondaryFromWeakDecay(iTracks);
      TrIsOthers=!TrIsPrim && !TrIsSecondMate && !TrIsSecondWeak;

      //if(mctruth->IsParticleSelected(iTracks)) hoge++;
      //cout<<fNEntries<<" "<<iTracks<<" "<<pdgcode<<" "<<<<" "<<TrIsPrim<<" "<<TrIsSecondMate<<" "<<TrIsSecondWeak<<endl;
    //  if(TrIsPrim){
    //  cout<<fNEntries<<" "<<iTracks<<" "<<track->PdgCode()<<endl;
      //rstack->DumpPart(iTracks);
    //  }


      mcTrackEta = track->Eta();
      mcTrackPt  = track->Pt();
      mcTrackPhi = track->Phi();
      TrCharge=track->Charge()!=0;
      pdgcode = TMath::Abs(track->PdgCode());
      if(!TrCharge)        continue;
      Int_t daughterindex1=track->GetFirstDaughter();
      if(daughterindex1>0){
        AliMCParticle* dparticle=(AliMCParticle*)mcEvent->GetTrack(daughterindex1);
      }
      Int_t motherIndex1 = track->GetMother();
      if(motherIndex1>0){
          TParticle *mParticle1=mcEvent->Particle(motherIndex1);
      }
      conmcprim[0]=mcTrackPt;
      conmcprim[1]=mcTrackEta;
      conmcprim[2]=mcTrackPhi;
      conmcprim[3]=lCentrality;
      conmcprim[4]=fPrimaryZVtx;
      if(TrIsPrim && motherIndex1<0) {
         fhistmcprim->Fill(conmcprim,0);//primay charged partilce distribution(no mother particle)
         fhmcprimvzeta->Fill(mcTrackEta,track->Zv());
         fhmcprimpdgcode->Fill(pdgcode);
         fh2_FMD_eta_phi_prim->Fill(mcTrackEta,mcTrackPhi);
       }

       //rstack-> DumpLoadedStack();//Dump particle(Too heavy, don't forget comment out before job submit)
       //rstack->DumpPart(iTracks);

     }
*/
  //   cout<<fNEntries<<" "<<"number of physicsprimary==="<<hoge<<endl;
      // cout<<fNEntries<<" "<<mcEvent->GetNumberOfTracks()<<" "<<hoge<<" "<<mcEvent->GetNumberOfTracks()-hoge<<endl;
/*
      Int_t nref=0;
      Int_t nfmdhits=0;
      for (Int_t iref = 0; iref < track->GetNumberOfTrackReferences(); iref++) {
         AliTrackReference *ref = track->GetTrackReference(iref);
         if (ref->DetectorId() != 12) continue;// Select FMD
         nfmdhits++;
         Float_t r = TMath::Sqrt(ref->X()*ref->X()+ref->Y()*ref->Y());
         Float_t z=ref->Z();
      //   Float_t zv=TMath::Abs(z-rParticle->Vz());
         Float_t zv=TMath::Abs(z-tPrimaryVtxPosition[2]);
         if(z<0) zv=-1.*zv;
         Double_t theta=TMath::ATan2(r,zv);
         Double_t etav= -1.*TMath::Log(TMath::Tan(theta/2.));
         Double_t phiv = TMath::ATan2(ref->Y(),ref->X());
         if(phiv<0) phiv=phiv+2*TMath::Pi();
         if(TrIsPrim) {
           fh2_FMD_acceptance_prim->Fill(etav,rParticle->Vz());
           //fh2_FMD_eta_phi_prim->Fill(etav,phiv);
         }
         if(TrIsPrim && motherIndex1<0) {
           fhrefetaFMD[0]->Fill(etav);
           fhrefphiFMD[0]->Fill(phiv);
         } else if(TrIsSecondMate) {
           fhrefetaFMD[1]->Fill(etav);
           fhrefphiFMD[1]->Fill(phiv);
         } else if(TrIsSecondWeak) {
           fhrefetaFMD[2]->Fill(etav);
           fhrefphiFMD[2]->Fill(phiv);
         }else if(TrIsOthers) {
           fhrefetaFMD[3]->Fill(etav);
           fhrefphiFMD[3]->Fill(phiv);
         }
       }
     }
*/
   }

 }

  
  FillCorrelationTracks(lCentrality,selectedTracksMC,selectedTracksMC,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
  FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksMC,selectedTracksMC,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  selectedTracksMC->Clear();
  delete selectedTracksMC;
  
  return;

  if(fAnaMode=="TPCTPC"){
    FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksAssociated,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedTracksAssociated,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }else if(fAnaMode=="TPCV0A"){
    FillCorrelationTracks(lCentrality,selectedTracksLeading ,selectedTrackV0A,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedTrackV0A,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }else if(fAnaMode=="TPCV0C"){
    FillCorrelationTracks(lCentrality,selectedTracksLeading ,selectedTrackV0C,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading ,selectedTrackV0C,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }  else if(fAnaMode=="V0AV0C"){
    FillCorrelationTracks(lCentrality,selectedTrackV0A,selectedTrackV0C,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTrackV0A,selectedTrackV0C,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }else if(fAnaMode=="TPCFMD"){
    FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedFMDArray1,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedFMDArray1,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }else if(fAnaMode=="TPCFMDC"){
    FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedFMDArray2,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedFMDArray2,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }else if(fAnaMode=="FMDFMD"){
    FillCorrelationTracks(lCentrality,selectedFMDArray1,selectedFMDArray2,fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0);
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedFMDArray1,selectedFMDArray2,fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
  }

  selectedFMDArray1->Clear();
  delete selectedFMDArray1;
  selectedFMDArray2->Clear();
  delete selectedFMDArray2;
  selectedTracksLeading->Clear();
  delete selectedTracksLeading;
  selectedTracksAssociated->Clear();
  delete selectedTracksAssociated;
  selectedTrackV0A->Clear();
  delete selectedTrackV0A;
  selectedTrackV0C->Clear();
  delete selectedTrackV0C;

}

TObjArray* AliAnalysisTaskSEpPbCorrelationsMCYS::GetFMDhitsYS(Bool_t Aside){
    TObjArray *tracks1 = new TObjArray;
    tracks1->SetOwner(kTRUE);

    AliAODForwardMult* aodForward =static_cast<AliAODForwardMult*>(fEvent->FindListObject("ForwardMC"));
	//cout<<aodForward->IsSecondaryCorrected()<<endl;
	//cout<<aodForward->GetCentrality()<<" "<<lCentrality<<endl;
	//cout<<aodForward->IsAcceptanceCorrected()<<endl;
	//cout<<aodForward->IsVertexBiasCorrected() <<endl;
	//cout<<aodForward->IsMergingEfficiencyCorrected()<<endl;

    //if(!aodForward) cout<<"no aodforadmult"<<endl;
	// Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
      const TH2D& d2Ndetadphi = aodForward->GetHistogram();

      Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
      Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
      //AliAnalysisTaskValidation::Tracks ret_vector;
      // FMD has no pt resolution!
      Double_t pt = 0;

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
          //  fh2_FMD_acceptance->Fill(eta,tPrimaryVtxPosition[2]);
        ///    fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);

            if(Aside){
              if(eta<0) continue;
             } else{
              if(eta>0) continue;
            }

			//            tracks1->Add(new AliAssociatedVZEROYS(mostProbableN,eta,phi,0,0,0));
            tracks1->Add(new AliAssociatedTrackYS(0,eta,phi,0,0,0,0,0,mostProbableN));
            Double_t cont[4]={eta,phi,lCentrality,tPrimaryVtxPosition[2]};
            //            cout<<eta<<" "<<phi<<" "<<lCentrality<<" "<<tPrimaryVtxPosition[2]<<endl;

			//			if( abs(tPrimaryVtxPosition[2])<1) {
            fhistfmd->Fill(cont,0,mostProbableN);
            fh2_FMD_acceptance->Fill(eta,tPrimaryVtxPosition[2],mostProbableN);
            fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
			//}

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
            if(track.eta>0) tracks1->Add(new AliAssociatedVZEROYS(track.eta,track.phi,track.weight,0,0,0));
            } else{
            if(track.eta<0) tracks1->Add(new AliAssociatedVZEROYS(track.eta,track.phi,track.weight,0,0,0));
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

TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedTracksLeading(AliAODEvent *event,Bool_t leading) {
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  Int_t nTracks = event->GetNumberOfTracks();
  Double_t pidqa[5];
  Double_t  conmcprim[5];
  TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack *fTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(i));
    //AliVTrack *fTrack = dynamic_cast<AliVTrack *>(event->GetTrack(i));
    if (!fTrack)      continue;
  //  if (!IsAcceptedTrackLeading(fTrack))      continue;
    if (!IsAcceptedTrack(fTrack))      continue;
    //if(!fIsAOD)if(!fTrackCuts->AcceptTrack((AliESDTrack)*fTrack)) continue;
    if (fTrack->Charge() == 0)      continue;
    if(leading){
     pidqa[0]=fTrack->Pt();
     pidqa[1]=fTrack->Eta();
     pidqa[2]=RangePhi(fTrack->Phi());
     pidqa[3]=lCentrality;
     pidqa[4]=fPrimaryZVtx;
     fHistLeadQA->Fill(pidqa,0);
	 /*
	 if(!fDataType) {
       Int_t a=0;
       Int_t label=abs(fTrack->GetLabel());
       AliAODMCParticle*mcTrack=(AliAODMCParticle *)mcArray->At(label);
       if(!mcTrack){
         continue;
       }
       while(!mcTrack->IsPhysicalPrimary()){
         if(mcTrack->GetMother()<0) break;
         mcTrack=(AliAODMCParticle*)mcArray->At(((AliAODMCParticle * )mcTrack)->GetMother());
         if(!mcTrack) break;
       }
       Double_t mcTrackEta = mcTrack->Eta();
       Double_t mcTrackPt  = mcTrack->Pt();
       Double_t mcTrackPhi = mcTrack->Phi();
       Int_t  pdgcode=TMath::Abs(mcTrack->PdgCode());      
      

       conmcprim[0]=mcTrackPt;
       conmcprim[1]=mcTrackEta;
       conmcprim[2]=mcTrackPhi;
       conmcprim[3]=lCentrality;
       conmcprim[4]=fPrimaryZVtx;
       
     }*/
    }
    
    Int_t SpAsso=0;
    tracks->Add(new AliAssociatedTrackYS(fTrack->Charge(), fTrack->Eta(), fTrack->Phi(), fTrack->Pt(), fTrack->GetID(), -999, -999, SpAsso, 1));
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

    tracks->Add(new AliAssociatedTrackYS(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, SpPID, 1));
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

    dtrack->Add(new AliMixTrackYS(
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
      tracks->Add(new AliAssociatedTrackYS(0, eta_phi, phi_phi, pt_phi, -999,
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
        AliMixTrackYS *dtrack1 = (AliMixTrackYS *)dtrack->At(i);
        if (!dtrack1)
          continue;
        Double_t pdx1 = dtrack1->Px();
        Double_t pdy1 = dtrack1->Py();
        Double_t pdz1 = dtrack1->Pz();
        Double_t Ed1 = TMath::Sqrt(pdx1 * pdx1 + pdy1 * pdy1 + pdz1 * pdz1 +
                                   0.493677 * 0.493677);
        counterMix++;
        for (Int_t j = 0; j < mixEvents->GetEntriesFast(); j++) {
          AliMixTrackYS *dtrack2 = (AliMixTrackYS *)mixEvents->At(j);
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
    AliMixTrackYS *particle = (AliMixTrackYS *)dtrack->At(i);
    tracksClone->Add(new AliMixTrackYS(
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

    const AliAODTrack* myTrackPos;
    const AliAODTrack* myTrackNeg;

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
        tracks->Add(new AliAssociatedTrackYS(aodv0->Charge(), aodv0->Eta(), aodv0->Phi(), aodv0->Pt(),aodv0->GetID(), myTrackPos->GetID(), myTrackNeg->GetID(), SpV0, 1));

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



TObjArray *AliAnalysisTaskSEpPbCorrelationsMCYS::GetAcceptedCascadeTracks(const AliAODEvent *fAODEvent) {
  // To select Cascade Particle
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  Int_t nCascades = fAODEvent->GetNumberOfCascades();
  Double_t lInvMassXiMinus;
  Double_t lInvMassXiPlus;
  Double_t lInvMassOmegaMinus;
  Double_t lInvMassOmegaPlus;
  for (Int_t icasc = 0; icasc < nCascades; icasc++) {
    const AliAODcascade *casc = fAODEvent->GetCascade(icasc);
    if (!casc)      continue;
    const AliAODTrack *myTrackCascPos;
    const AliAODTrack *myTrackCascNeg;
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
    tracks->Add(new AliAssociatedTrackYS(1., etaCas, phiCas, lPtCas, myTrackCascPos->GetID(),myTrackCascNeg->GetID(), myTrackCascBach->GetID(), SpCas, 1));
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


Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedTrackLeading(const AliVTrack *fTrack) {
  if (!fTrack)  return kFALSE;
  if(fIsAOD){
    if (!((AliAODTrack*)fTrack)->TestFilterMask(BIT(5)))  return kFALSE;
  }
    /*
    if (!fTrack->IsOn(AliVTrack::kTPCrefit)) return kFALSE;
    Float_t nCrossedRowsTPC =fTrack->GetTPCClusterInfo(2,1);
    if (nCrossedRowsTPC < 70) return kFALSE;
    Int_t findable=fTrack->GetTPCNclsF();
    if (findable <= 0)return kFALSE;
    if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  */
  if (fTrack->Pt() < fPtMin) return kFALSE;
  if (TMath::Abs(fTrack->Eta()) > fEtaMax)  return kFALSE;
  return kTRUE;
}


Bool_t AliAnalysisTaskSEpPbCorrelationsMCYS::IsAcceptedTrack(const AliAODTrack *aodTrack) {
  if (!aodTrack)
    return kFALSE;
  //  if(!aodTrack->TestFilterMask(BIT(5))) return kFALSE; // standard cut with
  //  tight DCA cut
  if (!aodTrack->TestFilterMask(BIT(5)))
    return kFALSE; // only tpc cut
  /*
  if (!aodTrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC =aodTrack->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=aodTrack->GetTPCNclsF();
  if (findable <= 0)return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  */
  if (aodTrack->Pt() < fPtMin)
    return kFALSE;
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
/*
void AliAnalysisTaskSEpPbCorrelationsMCYS::FillCorrelationTracksCentralForward( Double_t centrality, TObjArray *triggerArray, TObjArray *selectedTrackArray,AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut, Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,Float_t bSign, Int_t step) {

   Double_t  binscont[5];
    Double_t  binscontTrig[2];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++)     {
      AssociatedTrack* trigger =(AssociatedTrack*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t triggerPt   = trigger->Pt();
      Double_t triggerEta  = trigger->Eta();
      Double_t triggerPhi  = trigger->Phi();
      binscontTrig[0]=triggerPt;
      binscontTrig[1]=centrality;

      triggerHist->Fill(binscontTrig,step);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AssociatedVZERO* associate=(AssociatedVZERO*) selectedTrackArray->At(j);
        if(!associate)continue;
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=triggerPt;
        binscont[2]=associate->Eta();
        binscont[3]=centrality;
        binscont[4]=RangePhi(triggerPhi-associate->Phi());
        associateHist->Fill(binscont,step,(Double_t)associate->Multiplicity());
      }
    }



}
*/
void AliAnalysisTaskSEpPbCorrelationsMCYS::FillCorrelationTracks( Double_t centrality, TObjArray *triggerArray, TObjArray *selectedTrackArray,AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut, Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,Float_t bSign, Int_t step) {
  if (!triggerHist || !associateHist)    return;
  Double_t binscontTrig[2];
  Double_t binscont[5];

  if(fAnaMode=="TPCTPC"){
    for (Int_t i = 0; i < triggerArray->GetEntriesFast(); i++) {
      AliAssociatedTrackYS *trigger = (AliAssociatedTrackYS *)triggerArray->At(i);
      if (!trigger)    continue;
      Double_t triggerPt = trigger->Pt();
      Double_t triggerEta = trigger->Eta();
      Double_t triggerPhi = trigger->Phi();
      Int_t trigFirstID = trigger->GetIDFirstDaughter();
      Int_t trigSecondID = trigger->GetIDSecondDaughter();
      Int_t trigID = trigger->GetID();
      binscontTrig[0] = triggerPt;
      binscontTrig[1] = centrality;
      triggerHist->Fill(binscontTrig, 0);
      for (Int_t j = 0; j < selectedTrackArray->GetEntriesFast(); j++) {
        AliAssociatedTrackYS *associate =   (AliAssociatedTrackYS*)selectedTrackArray->At(j);
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
        Int_t SpAsso = associate->WhichCandidate();
        if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" ||  (fasso == "PID")) {
          if (SpAsso < 0)          continue;
          associateHist->Fill(binscont, SpAsso);
        }else if(fasso=="hadron"){
          associateHist->Fill(binscont, 0);
        }
      }
    }
  }else if (fAnaMode=="TPCV0A" || fAnaMode=="TPCV0C" || fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t  triggerEta  = trigger->Eta();
      Double_t  triggerPhi  = trigger->Phi();
      Double_t  triggerPt  = trigger->Pt();
      binscontTrig[1]=centrality;
      binscontTrig[0]=triggerPt;
      Int_t SpAsso= trigger->WhichCandidate();
      triggerHist->Fill(binscontTrig,SpAsso);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(j);
        if(!associate)continue;
        Double_t associatemultiplicity=associate->Multiplicity();
        Double_t assophi=associate->Phi();
        Double_t assoeta=associate->Eta();
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=triggerPt;
        binscont[2]=associate->Eta();
        binscont[3]=centrality;
        binscont[4]=RangePhi(triggerPhi-associate->Phi());
        if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" ||  (fasso == "PID")) {
          if (SpAsso < 0)          continue;
          associateHist->Fill(binscont, SpAsso,(Double_t)associate->Multiplicity());
        }else{
          associateHist->Fill(binscont, 0, (Double_t)associate->Multiplicity());
        }
      }
    }
  }else if(fAnaMode=="FMDFMD"){
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t  triggerMultiplicity= trigger->Multiplicity();
      Double_t  triggerEta  = trigger->Eta();
      Double_t  triggerPhi  = trigger->Phi();
      binscontTrig[0]=centrality;
      binscontTrig[1]=triggerEta;
      triggerHist->Fill(binscontTrig,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(j);
        if(!associate)continue;
        Double_t associatemultiplicity=associate->Multiplicity();
        Double_t assophi=associate->Phi();
        Double_t assoeta=associate->Eta();
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=associate->Eta();
        binscont[2]=triggerEta;
        binscont[3]=centrality;
        binscont[4]=RangePhi_FMD(triggerPhi-associate->Phi());
        //triggerEta-associate->Eta();
        Double_t dphivzero=triggerPhi-associate->Phi();
        if(triggerPhi==assophi && triggerEta==assoeta) continue;
        associateHist->Fill(binscont,0,(Double_t)triggerMultiplicity*associatemultiplicity);
      }
    }
  }else if(fAnaMode=="V0AV0C"){
    Double_t  binscont1[4];
    Double_t  binscontTrig1[3];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t  triggerMultiplicity= trigger->Multiplicity();
      Double_t  triggerEta  = trigger->Eta();
      Double_t  triggerPhi  = trigger->Phi();
      binscontTrig1[0]=centrality;
      binscontTrig1[1]=triggerEta;
      binscontTrig1[2]=triggerPhi;
      triggerHist->Fill(binscontTrig1,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(j);
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
  if (!triggerHist || !associateHist) return;
  Double_t counterMix = 0;
  AliEventPool *pool = fPoolMgr->GetEventPool(centrality, pvxMix);
  if (!pool)    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality,
  pvxMix));
  //    cout<<pool->NTracksInPool()<<" "<<fPoolMinNTracks<<" "<<pool->GetCurrentNEvents()<<" "<<fMinEventsToMix<<endl;
  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix = 0; jMix < nMix; jMix++) {   TObjArray *mixEvents = pool->GetEvent(jMix);
      Double_t binscontTrig[2];
      Double_t binscont[5];
      if(fAnaMode=="TPCTPC"){
        for (Int_t i = 0; i < triggerArray->GetEntriesFast(); i++) {
          AliAssociatedTrackYS *trig = (AliAssociatedTrackYS *)triggerArray->At(i);
          if (!trig)          continue;
          Double_t triggerPhi = trig->Phi();
          Double_t triggerEta = trig->Eta();
          Double_t triggerPt = trig->Pt();
          counterMix++;
          binscontTrig[0] = triggerPt;
          binscontTrig[1] = centrality;
          triggerHist->Fill(binscontTrig, step);
          for (Int_t j = 0; j < mixEvents->GetEntriesFast(); j++) {
            AliAssociatedTrackYS *associate =  (AliAssociatedTrackYS *)mixEvents->At(j);
            if (!associate) continue;
            /*
            if(triggerEta<=0){
            if(associate->Eta()<=0) continue;
          }else if(triggerEta>0){
          if(associate->Eta()>0) continue;
        }
        */
        binscont[0] = triggerEta - associate->Eta();
        binscont[1] = associate->Pt();
        binscont[2] = triggerPt;
        binscont[3] = centrality;
        binscont[4] = RangePhi(triggerPhi - associate->Phi());
        Int_t SpAsso = associate->WhichCandidate();
        if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" || (fasso == "PID")) {
          if (SpAsso < 0)   continue;
          associateHist->Fill(binscont, SpAsso, 1. / (Double_t)nMix);
        }else if(fasso=="hadron"){
          associateHist->Fill(binscont, 0,1./(Double_t)nMix);
        }
      }
    }
  }else if(fAnaMode=="TPCV0A" || fAnaMode=="TPCV0C" || fAnaMode=="TPCFMD" || fAnaMode=="TPCFMDC"){
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYS* trigger =(AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t triggerPt   = trigger->Pt();
      Double_t triggerEta  = trigger->Eta();
      Double_t triggerPhi  = trigger->Phi();
      Int_t SpAsso=trigger->WhichCandidate();
      counterMix++;
      binscontTrig[0]=triggerPt;
      binscontTrig[1]=centrality;
      triggerHist->Fill(binscontTrig,SpAsso);
	  for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
        AliAssociatedTrackYS* associate=(AliAssociatedTrackYS*)  mixEvents->At(j);
        if(!associate)continue;
        binscont[0]=triggerEta-associate->Eta();
        binscont[1]=triggerPt;
        binscont[2]=associate->Eta();
        binscont[3]=centrality;
        binscont[4]=RangePhi(triggerPhi-associate->Phi());
		//		cout<<centrality<<endl;
        if (fasso == "V0" || fasso == "Phi" || fasso == "Cascade" || (fasso == "PID")) {
          if (SpAsso < 0)   continue;
          associateHist->Fill(binscont, SpAsso, (Double_t)associate->Multiplicity()/(Double_t)nMix);
        }else if(fasso=="hadron"){
          associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix);
        }
      }
    }
    }else if(fAnaMode=="FMDFMD"){
      Double_t  binscont1[4];
      Double_t  binscontTrig1[3];
      for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
        AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
        if(!trigger)continue;
        Double_t  triggerMultiplicity= trigger->Multiplicity();
        Double_t  triggerEta  = trigger->Eta();
        Double_t  triggerPhi  = trigger->Phi();
        counterMix++;
        binscontTrig[0]=centrality;
        binscontTrig[1]=triggerEta;
        triggerHist->Fill(binscontTrig,step,(Double_t)triggerMultiplicity);
        for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
          AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) mixEvents->At(j);
          if(!associate)continue;
          Double_t associatemultiplicity=associate->Multiplicity();
          Double_t assophi=associate->Phi();
          Double_t assoeta=associate->Eta();
          binscont[0]=triggerEta-associate->Eta();
          binscont[1]=associate->Eta();
          binscont[2]=triggerEta;
          binscont[3]=centrality;
          binscont[4]=RangePhi_FMD(triggerPhi-associate->Phi());
          if(triggerPhi==assophi && triggerEta==assoeta) continue;
          associateHist->Fill(binscont,0,(Double_t)triggerMultiplicity*associatemultiplicity/(Double_t)nMix);
        }
      }
    }else if (fAnaMode=="V0AV0C"){
      Double_t  binscont1[4];
      Double_t  binscontTrig1[3];
      for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
        AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
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
          AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) mixEvents->At(j);
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

  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  //  cout<<selectedTrackArray->GetEntriesFast()<<endl;
  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliAssociatedTrackYS *particle =  (AliAssociatedTrackYS *)selectedTrackArray->At(i);
    tracksClone->Add(new AliAssociatedTrackYS(
      particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
      particle->GetID(), particle->GetIDFirstDaughter(),
      particle->GetIDSecondDaughter(), particle->WhichCandidate(),
      particle->Multiplicity()));
    }
   pool->UpdatePool(tracksClone);
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
