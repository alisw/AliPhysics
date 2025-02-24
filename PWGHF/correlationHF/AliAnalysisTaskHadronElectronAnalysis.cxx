

#include "AliAnalysisTaskHadronElectronAnalysis.h"
#include "AliKFParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"

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

#include "TGeoGlobalMagField.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"


#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskHadronElectronAnalysis;
ClassImp(AliAnalysisTaskHadronElectronAnalysis);


AliAnalysisTaskHadronElectronAnalysis::AliAnalysisTaskHadronElectronAnalysis() :
    AliAnalysisTaskSE(),
    fVevent{0},
    fMCarray{0},
    fMCHeader{0},
    fNTotMCpart{0},
    fNembMCpi0{0},
    fNembMCeta{0},
    fNpureMC{0},
fSprsPi0EtaWeightCal{0},
fIsFrmEmbPi0{kFALSE},
fIsFrmEmbEta{kFALSE},
ftype{-1},
fWeightPi0{1},
fWeightEta{1},
fWeight{1},
fRealInclsElecPt{0},
fNonHFeTrkPt{0},
fEtaWeight{0},
fPi0Weight{0},
fNonHFeEmbTrkPt{0},
fPi0eEmbWeightTrkPt{0},
fNonHFeEmbWeightTrkPt{0},
fRecoNonHFeTrkPt{0},
fEtaeEmbWeightTrkPt{0},
fRecoNonHFeEmbTrkPt{0},
fRecoPi0eEmbWeightTrkPt{0},
fRecoNonHFeEmbWeightTrkPt{0},
fRecoEtaeEmbWeightTrkPt{0},
fIsNew{kFALSE},
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHHCont{0},
    fDphiHHContEff{0},
    fDphiHElec{0},
    fDphiHElecTriggered{0},
    fDphiHUSElec{0},
    fDphiHLSElec{0},
fDphiHUSElecEff{0},
fDphiHLSElecEff{0},
    fDphiHElecEff{0},
    fDphiHElecTriggeredEff{0},
    fDphiHH{0},
    fDphiHHEff{0},
    fDphiHElecMixed{0},
    fDphiHULSElecMixed{0},
    fDphiHLSElecMixed{0},
    fDphiHHMixed{0},
fTriggEffVsPt{0},
fAssoEffVsPt{0},
    fCorPoolMgr{0},
    fAssociatedEff{0},
    fTriggerEff{0},
fTriggerEffElec{0},
fAssociatedEffElec{0},
fTriggerEffH{0},
fAssociatedEffH{0},
fTrigPhiDist{0},
 fAssHPhiDist{0},
 fAssElecPhiDist{0},
fmultPercentile{0},
fmultPercentileCuts{0},
    fDedx{0},
    fDedxCuts{0},
    fDedxTOF{0},
    fDedxTOFCuts{0},
    fBetaDedx{0},
    fBeta{0},
    fNSigma{0},
    fNSigmaCuts{0},
f2DPPtHist{0},
fNSigmaAssoNoCuts{0},
fNSigmaAssoNoCuts_woTOFwSPDKAny{0},
fNSigmaAssoNoCuts_woTOFwSPDKBoth{0},
fNSigmaAssoNoCuts_wTOFwSPDKAny{0},
fNSigmaAssoNoCuts_wTOFwSPDKBoth{0},
fNSigmaAssoCuts{0},
fNSigmaAssoCuts_woTOFwSPDKAny{0},
fNSigmaAssoCuts_woTOFwSPDKBoth{0},
fNSigmaAssoCuts_wTOFwSPDKAny{0},
fNSigmaAssoCuts_wTOFwSPDKBoth{0},
fNSigmaAssoCutsPt{0},
fNSigmaAssoCutsPt_woTOFwSPDKAny{0},
fNSigmaAssoCutsPt_woTOFwSPDKBoth{0},
fNSigmaAssoCutsPt_wTOFwSPDKAny{0},
fNSigmaAssoCutsPt_wTOFwSPDKBoth{0},
f2DHEpTDist{0},
f2DHHpTDist{0},
fAllHadronsPt{0},
fAssHadronsPt{0},
fTriggerPt{0},
fElecPtBeforeCuts{0},
fElecPtAfterCuts{0},
fNumberTriggersPerEvent{0},
fNumberElecsPerEvent{0},
fTOFHits{0},
fHEassNoTrig{0},
fHEnoAssTrig{0},
fHHassNoTrig{0},
 fHHnoAssTrig{0},
fTPCCrossedRows{0},
fUnlikesignMassDist{0},
fLikesignMassDist{0},
fHEDeltaEtaDeltaPhi{0},
fEtaCut{0},
fInvariantMassCut{0},
fPartnerElectronPtCut{0},
fTPCNSigmaElecMin{0},
fTPCNSigmaElecMax{0},
fPartnerTPCNSigmaElecMin{0},
fPartnerTPCNSigmaElecMax{0},
fAssocTrkBitElec{0},
fAssocTrkBitH{0},
fPartnerTrkBit{0},
fTrigTrkBit{0},
fUseSPDKAny{0},
fUseSPDKBoth{0},
fRunOnMC{0},
fUseTOFCut{0}
{
    // DO NOT TOUCH
    MULT_LOW = 0;
    MULT_HIGH = 100;
    CENT_ESTIMATOR = "V0A";
    EFF_FILE_PATH = "eff_out.root";
    }
//_____________________________________________________________________________
AliAnalysisTaskHadronElectronAnalysis::AliAnalysisTaskHadronElectronAnalysis(const char *name) :
    AliAnalysisTaskSE(name),
    fVevent{0},
    fMCarray{0},
    fMCHeader{0},
   fNTotMCpart{0},
    fNembMCpi0{0},
    fNembMCeta{0},
    fNpureMC{0},
fSprsPi0EtaWeightCal{0},
fIsFrmEmbPi0{kFALSE},
fIsFrmEmbEta{kFALSE},
ftype{-1},
fWeightPi0{1},
fWeightEta{1},
fWeight{1},
fRealInclsElecPt{0},
fNonHFeTrkPt{0},
fEtaWeight{0},
fPi0Weight{0},
fNonHFeEmbTrkPt{0},
fPi0eEmbWeightTrkPt{0},
fNonHFeEmbWeightTrkPt{0},
fRecoNonHFeTrkPt{0},
fEtaeEmbWeightTrkPt{0},
fRecoNonHFeEmbTrkPt{0},
fRecoPi0eEmbWeightTrkPt{0},
fRecoNonHFeEmbWeightTrkPt{0},
fRecoEtaeEmbWeightTrkPt{0},
fIsNew{kFALSE},
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHHCont{0},
    fDphiHHContEff{0},
    fDphiHElec{0},
    fDphiHElecTriggered{0},
    fDphiHUSElec{0},
    fDphiHLSElec{0},
fDphiHUSElecEff{0},
fDphiHLSElecEff{0},
    fDphiHElecEff{0},
    fDphiHElecTriggeredEff{0},
    fDphiHH{0},
    fDphiHHEff{0},
    fDphiHElecMixed{0},
    fDphiHULSElecMixed{0},
    fDphiHLSElecMixed{0},
    fDphiHHMixed{0},
fTriggEffVsPt{0},
fAssoEffVsPt{0},
    fCorPoolMgr{0},
    fAssociatedEff{0},
    fTriggerEff{0},
fTriggerEffElec{0},
fAssociatedEffElec{0},
fTriggerEffH{0},
fAssociatedEffH{0},
fTrigPhiDist{0},
 fAssHPhiDist{0},
 fAssElecPhiDist{0},
fmultPercentile{0},
fmultPercentileCuts{0},
    fDedx{0},
    fDedxCuts{0},
    fDedxTOF{0},
    fDedxTOFCuts{0},
    fBetaDedx{0},
    fBeta{0},
    fNSigma{0},
    fNSigmaCuts{0},
f2DPPtHist{0},
fNSigmaAssoNoCuts{0},
fNSigmaAssoNoCuts_woTOFwSPDKAny{0},
fNSigmaAssoNoCuts_woTOFwSPDKBoth{0},
fNSigmaAssoNoCuts_wTOFwSPDKAny{0},
fNSigmaAssoNoCuts_wTOFwSPDKBoth{0},
fNSigmaAssoCuts{0},
fNSigmaAssoCuts_woTOFwSPDKAny{0},
fNSigmaAssoCuts_woTOFwSPDKBoth{0},
fNSigmaAssoCuts_wTOFwSPDKAny{0},
fNSigmaAssoCuts_wTOFwSPDKBoth{0},
fNSigmaAssoCutsPt{0},
fNSigmaAssoCutsPt_woTOFwSPDKAny{0},
fNSigmaAssoCutsPt_woTOFwSPDKBoth{0},
fNSigmaAssoCutsPt_wTOFwSPDKAny{0},
fNSigmaAssoCutsPt_wTOFwSPDKBoth{0},
f2DHEpTDist{0},
f2DHHpTDist{0},
fAllHadronsPt{0},
fAssHadronsPt{0},
fTriggerPt{0},
fElecPtBeforeCuts{0},
fElecPtAfterCuts{0},
fNumberTriggersPerEvent{0},
fNumberElecsPerEvent{0},
fTOFHits{0},
fTPCCrossedRows{0},
fHEassNoTrig{0},
fHEnoAssTrig{0},
fHHassNoTrig{0},
 fHHnoAssTrig{0},
fUnlikesignMassDist{0},
fLikesignMassDist{0},
fHEDeltaEtaDeltaPhi{0},
fEtaCut{0},
fInvariantMassCut{0},
fPartnerElectronPtCut{0},
fTPCNSigmaElecMin{0},
fTPCNSigmaElecMax{0},
fPartnerTPCNSigmaElecMin{0},
fPartnerTPCNSigmaElecMax{0},
fAssocTrkBitElec{0},
fAssocTrkBitH{0},
fPartnerTrkBit{0},
fTrigTrkBit{0},
fUseSPDKAny{0},
fUseSPDKBoth{0},
fRunOnMC{0},
fUseTOFCut{0}
{
  // DO NOT TOUCH
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    MULT_LOW = 0;
    MULT_HIGH =100;
    CENT_ESTIMATOR = "V0A";
    EFF_FILE_PATH = "eff_out.root";
}
//_____________________________________________________________________________
AliAnalysisTaskHadronElectronAnalysis::~AliAnalysisTaskHadronElectronAnalysis()
{

    if(fOutputList) delete fOutputList;
    if(fTriggerEff) delete fTriggerEff;
    if(fAssociatedEff) delete fAssociatedEff;

}

void AliAnalysisTaskHadronElectronAnalysis::SetEtaCut(float etaCut) {
  if (TMath::Abs(etaCut) > 4.0) {
    AliFatal("Eta cut set above for, please check input parameter ETA_CUT!");
  }
  fEtaCut = etaCut;
}

void AliAnalysisTaskHadronElectronAnalysis::SetPartnerElectronPtCut(double partnerElectronPtCut) {
 
  fPartnerElectronPtCut = partnerElectronPtCut;
}

void AliAnalysisTaskHadronElectronAnalysis::SetTPCNSigmaElecMin(double TPCNSigmaElecMin){
  fTPCNSigmaElecMin = TPCNSigmaElecMin;
}

void AliAnalysisTaskHadronElectronAnalysis::SetTPCNSigmaElecMax(double TPCNSigmaElecMax){
  fTPCNSigmaElecMax = TPCNSigmaElecMax;
}

void AliAnalysisTaskHadronElectronAnalysis::SetPartnerTPCNSigmaElecMin(double partnerTPCNSigmaElecMin){
  fPartnerTPCNSigmaElecMin = partnerTPCNSigmaElecMin;
}

void AliAnalysisTaskHadronElectronAnalysis::SetPartnerTPCNSigmaElecMax(double partnerTPCNSigmaElecMax){
  fPartnerTPCNSigmaElecMax = partnerTPCNSigmaElecMax;
}

void AliAnalysisTaskHadronElectronAnalysis::SetAssocTrkBitElec(double assocTrkBitElec){
  fAssocTrkBitElec = assocTrkBitElec;
}

void AliAnalysisTaskHadronElectronAnalysis::SetAssocTrkBitH(double assocTrkBitH){
  fAssocTrkBitH = assocTrkBitH;
}

void AliAnalysisTaskHadronElectronAnalysis::SetPartnerTrkBit(double partnerTrkBit){
  fPartnerTrkBit = partnerTrkBit;
}

void AliAnalysisTaskHadronElectronAnalysis::SetTrigTrkBit(double trigTrkBit){
  fTrigTrkBit = trigTrkBit;
}

void AliAnalysisTaskHadronElectronAnalysis::SetUseSPDKAny(bool useSPDKAny){
    fUseSPDKAny = useSPDKAny;
}

void AliAnalysisTaskHadronElectronAnalysis::SetUseSPDKBoth(bool useSPDKBoth){
  fUseSPDKBoth = useSPDKBoth;
}

void AliAnalysisTaskHadronElectronAnalysis::SetRunOnMC(bool runOnMC){
  fRunOnMC = runOnMC;
}

void AliAnalysisTaskHadronElectronAnalysis::SetUseTOFCut(bool useTOFCut){
  fUseTOFCut = useTOFCut;
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //dedx plot
    
    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {MULT_LOW, MULT_HIGH};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    //Distribution axes are: Pt, Phi, Eta, zVtx, multiplicity percentile
    int dist_bins[5] = {100, 16, 20, 10, 100};
    double dist_mins[5] = {0, 0, -1, -10, 0};
    double dist_maxes[5] = {15, 6.28, 1, 10, 100};

    fLooseDist = new THnSparseF("fLooseDist", "All Hadron Distribution", 5, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fLooseDist);

    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 5, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fTriggerDist);

    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 5, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fAssociatedHDist);

    fElectronDist = new THnSparseF("fElectronDist", "Electron Distribution", 5, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fElectronDist);
    
//Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx, multiplicity percentile
   // int hl_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    int hl_cor_bins[7] = {11, 9, 16, 20, 100, 10, 100};//from 32 to 16
    double hl_cor_mins[7] = {1.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10, 0};
    double hl_cor_maxes[7] = {12.0, 10, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10, 100};

    fDphiHHCont = new THnSparseF("fDphiHHCont", "Hadron-Hadron Contamination Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHHCont);

    fDphiHHContEff = new THnSparseF("fDphiHHContEff", "Efficiency Corrected Hadron-Hadron Contamination Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHHContEff);

    fDphiHElec = new THnSparseF("fDphiHElec", "Hadron-Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHElec);

    fDphiHElecTriggered = new THnSparseF("fDphiHElecTriggered", "Hadron-Electron Correlation Histogram With Triggered Event", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHElecTriggered);

    fDphiHUSElec = new THnSparseF("fDphiHUSElec", "Hadron-US Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHUSElec);
    
    fDphiHLSElec = new THnSparseF("fDphiHLSElec", "Hadron-LS Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHLSElec);
    
    fDphiHUSElecEff = new THnSparseF(" fDphiHUSElecEff", "Efficiency Corrected Hadron-ULS Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add( fDphiHUSElecEff);
    
    fDphiHLSElecEff = new THnSparseF("fDphiHLSElecEff", "Efficiency Corrected Hadron-LS Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHLSElecEff);
    
    fDphiHElecEff = new THnSparseF("fDphiHElecEff", "Efficiency-corrected Hadron-Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHElecEff);

    fDphiHElecTriggeredEff = new THnSparseF("fDphiHElecTriggeredEff", "Efficiency-corrected Hadron-Electron Correlation Histogram with Trigered Event", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHElecTriggeredEff);

    fDphiHElecMixed = new THnSparseF("fDphiHElecMixed", "Mixed Hadron-Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHElecMixed);

    fDphiHULSElecMixed= new THnSparseF("fDphiHULSElecMixed", "Mixed Hadron-ULS Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHULSElecMixed);
    
    fDphiHLSElecMixed = new THnSparseF("fDphiHLSElecMixed", "Mixed Hadron-LS Electron Correlation Histogram", 7, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fOutputList->Add(fDphiHLSElecMixed);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx, mulitplicity percentile
    int hh_cor_bins[6] = {11, 9, 16, 20, 10, 100};
    double hh_cor_mins[6] = {1, 1, -1.0*TMath::Pi()/2.0, -2.0, -10, 0};
    double hh_cor_maxes[6] = {12, 10, 3.0*TMath::Pi()/2.0, 2.0, 10, 100};

    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 6, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiHH);

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 6, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiHHEff);

    fDphiTriggerTrigger = new THnSparseF("fDphiTriggerTrigger", "Trigger-Trigger Correlation Histogram", 6, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiTriggerTrigger);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 6, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiHHMixed);

    fDphiTriggerTriggerMixed = new THnSparseF("fDphiTriggerTriggerMixed", "MixedTrigger-Trigger Correlation Histogram", 6, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fOutputList->Add(fDphiTriggerTriggerMixed);

//Mixed Event Trigger and Associated vs pt
    fTriggerEffElec= new TH2D("TriggerEffElec", "TriggerEffElec;p_{T} (GeV/c);counts", 100, 0, 10, 1000, 0, 100);
    fOutputList->Add(fTriggerEffElec);
    
    fAssociatedEffElec= new TH2D("AssociatedEffElec", "AssociatedEffElec;p_{T} (GeV/c);counts", 100, 0, 10, 1000, 0, 100);
    fOutputList->Add( fAssociatedEffElec);
    
    fTriggerEffH= new TH2D("TriggerEffH", "TriggerEffH;p_{T} (GeV/c);counts", 100, 0, 10, 1000, 0, 100);
    fOutputList->Add(fTriggerEffH);
    
    fAssociatedEffH= new TH2D("AssociatedEffH", "AssociatedEffH;p_{T} (GeV/c);counts", 100, 0, 10, 1000,0, 100);
    fOutputList->Add(fAssociatedEffH);
    
    
    
    //phi distributions
    fTrigPhiDist= new TH1D("TriggerPhi", "TriggerPhi", 100, 0, 10);
    fOutputList->Add(fTrigPhiDist);
    
    fAssHPhiDist = new TH1D("AssociatedHPhi", "AssociatedHPhi", 100, 0, 10);
    fOutputList->Add(fAssHPhiDist);
    
    fAssElecPhiDist = new TH1D("AssociatedElecPhi", "AssociatedElecPhi", 100, 0, 10);
    fOutputList->Add(fAssElecPhiDist);
    
    //Multiplicity Percentile Distribution
    fmultPercentile= new TH1D("MultPercentile", "MultPercentile", 100, 0, 100);
    fOutputList->Add( fmultPercentile);
    
    fmultPercentileCuts= new TH1D("MultPercentileCuts", "MultPercentileCuts", 100, 0, 100);
    fOutputList->Add( fmultPercentileCuts);
    
    //axes are pion 0: dca pion 1: dca electron
    //              2: pt pion  3: pt electron
    //              4: pt L     5: mass L
    int dca_bins[6] = {100, 100, 20, 20, 20, 50};
    double dca_mins[6] = {-2.4, -2.4, 0, 0, 0, 1.08};
    double dca_maxes[6] = {2.4, 2.4, 10, 10, 10, 1.16};

    
    //fDedx
    
    fDedx= new TH2D("dedx", "dedx", 1000, -10, 10, 12000, 0, 120);
    fOutputList->Add(fDedx);
    
    fDedxCuts= new TH2D("dedxCuts", "dedxCuts", 100,-10, 10, 1200, 0, 120);
    fOutputList->Add(fDedxCuts);
    
    fDedxTOF= new TH2D("dedxTOF", "dedxTOF", 100, -10, 10, 1200, 0, 120);
    fOutputList->Add(fDedxTOF);
    
    fDedxTOFCuts= new TH2D("dedxTOFCuts", "dedxTOFCuts", 100, -10, 10, 120, 0, 120);
    fOutputList->Add(fDedxTOFCuts);
    
    fBetaDedx= new TH2D("betaDedx", "betaDedx", 1200, 0, 120, 1200, 0, 120);
    fOutputList->Add(fBetaDedx);
    
    fBeta= new TH2D("beta", "beta", 100, 0, 10, 500000, 0, 2);
    fOutputList->Add(fBeta);
   
    fNSigma= new TH2D("NSigma", "NSigma", 1000, 0, 10, 12000, -15 , 10);
    fOutputList->Add(fNSigma);
    
    fNSigmaCuts= new TH2D("NSigmaCuts", "NSigmaCuts", 100, 0, 10, 1200, -10, 10);
    fOutputList->Add(fNSigmaCuts);
    
    f2DPPtHist= new TH2D("2DPPtHist", "2DPPtHist", 100, 0, 10, 100, 0, 10);
    fOutputList->Add(f2DPPtHist);

    fNSigmaAssoNoCuts= new TH2D("NSigmaAssoNoCuts", "NSigmaAssoNOCuts", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoNoCuts);

    fNSigmaAssoNoCuts_woTOFwSPDKAny= new TH2D("NSigmaAssoNoCuts_woTOFwSPDKany", "NSigmaAssoNOCuts_woTOFwSPDKany", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoNoCuts_woTOFwSPDKAny);

    fNSigmaAssoNoCuts_woTOFwSPDKBoth= new TH2D("NSigmaAssoNoCuts_woTOFwSPDKBoth", "NSigmaAssoNoCuts_woTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoNoCuts_woTOFwSPDKBoth);

    fNSigmaAssoNoCuts_wTOFwSPDKAny= new TH2D("NSigmaAssoNoCuts_wTOFwSPDKany", "NSigmaAssoNoCuts_wTOFwSPDKany", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoNoCuts_wTOFwSPDKAny);

    fNSigmaAssoNoCuts_wTOFwSPDKBoth= new TH2D("NSigmaAssoNoCuts_wTOFwSPDKBoth", "NSigmaAssoNoCuts_wTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoNoCuts_wTOFwSPDKBoth);

    fNSigmaAssoCuts= new TH2D("NSigmaAssoCuts", "NSigmaAssoCuts", 100, 0, 10, 1200, -15, 10);//change later to 2500 or something. 
    fOutputList->Add(fNSigmaAssoCuts);

    fNSigmaAssoCuts_woTOFwSPDKAny= new TH2D("NSigmaAssoCuts_woTOFwSPDKAny", "NSigmaAssoCuts_woTOFwSPDKAny", 100, 0, 10, 1200, -15, 10);//change later to 2500 or something. 
    fOutputList->Add(fNSigmaAssoCuts_woTOFwSPDKAny);
    
     fNSigmaAssoCuts_woTOFwSPDKBoth= new TH2D("NSigmaAssoCuts_woTOFwSPDKBoth", "NSigmaAssoCuts_woTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);//change later to 2500 or something. 
    fOutputList->Add(fNSigmaAssoCuts_woTOFwSPDKBoth);

     fNSigmaAssoCuts_wTOFwSPDKAny= new TH2D("NSigmaAssoCuts_wTOFwSPDKAny", "NSigmaAssoCuts_wTOFwSPDKAny", 100, 0, 10, 1200, -15, 10);//change later to 2500 or something. 
    fOutputList->Add(fNSigmaAssoCuts_wTOFwSPDKAny);
 
    fNSigmaAssoCuts_wTOFwSPDKBoth= new TH2D("NSigmaAssoCuts_wTOFwSPDKBoth", "NSigmaAssoCuts_wTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);//change later to 2500 or something. 
    fOutputList->Add(fNSigmaAssoCuts_wTOFwSPDKBoth);

    fNSigmaAssoCutsPt= new TH2D("NSigmaAssoCutsPt", "NSigmaAssoCutsPt", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoCutsPt);

    fNSigmaAssoCutsPt_woTOFwSPDKAny= new TH2D("NSigmaAssoCutsPt_woTOFwSPDAny", "NSigmaAssoCutsPt_woTOFwSPDKAny", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoCutsPt_woTOFwSPDKAny);

    fNSigmaAssoCutsPt_woTOFwSPDKBoth= new TH2D("NSigmaAssoCutsPt_woTOFwSPDKBoth", "NSigmaAssoCutsPt_woTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoCutsPt_woTOFwSPDKBoth);

    fNSigmaAssoCutsPt_wTOFwSPDKAny= new TH2D("NSigmaAssoCutsPt_wTOFwSPDKAny", "NSigmaAssoCutsPt_wTOFwSPDKAny", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoCutsPt_wTOFwSPDKAny);

    fNSigmaAssoCutsPt_wTOFwSPDKBoth= new TH2D("NSigmaAssoCutsPt_wTOFwSPDKBoth", "NSigmaAssoCutsPt_wTOFwSPDKBoth", 100, 0, 10, 1200, -15, 10);
    fOutputList->Add(fNSigmaAssoCutsPt_wTOFwSPDKBoth);

    f2DHEpTDist= new TH2D("f2DHEpTDist", "f2DHEpTDist", 10, 0, 10,10,0,10);
    fOutputList->Add(f2DHEpTDist);
    
    f2DHHpTDist= new TH2D("f2DHHpTDist", "f2DHHpTDist", 10, 0, 10,10,0,10);
    fOutputList->Add(f2DHHpTDist);
    
    fTriggEffVsPt= new TH2D(" fTriggEffVsPt", " fTriggEffVsPt", 100, 0, 100,10,0,10);
    fOutputList->Add( fTriggEffVsPt);
    
    fAssoEffVsPt= new TH2D("fAssoEffVsPt", "fAssoEffVsPt", 100, 0, 100,10,0,10);
    fOutputList->Add(fAssoEffVsPt);
    
    fAllHadronsPt=new TH1D("AllHadronsPt", "AllHadronsPt", 100, 0, 10);
    fOutputList->Add(fAllHadronsPt);
    
    fAssHadronsPt=new TH1D("AssHadronsPt", "AssHadronsPt", 100, 0, 10);
    fOutputList->Add(fAssHadronsPt);
    
    fTriggerPt=new TH1D("TriggerPt", "TriggerPt", 100, 0, 10);
    fOutputList->Add(fTriggerPt);
    
    fElecPtBeforeCuts= new TH1D("ElecPtBeforeCuts", "ElecPtBeforeCuts", 100, 0, 10);
    fOutputList->Add(fElecPtBeforeCuts);
    
    fElecPtAfterCuts= new TH1D("ElecPtAfterCuts", "ElecPtAfterCuts", 100, 0, 10);
    fOutputList->Add(fElecPtAfterCuts);
    
    fNumberTriggersPerEvent= new TH1D("NumberTriggersPerEvent", "NumberTriggersPerEvent", 10, 0, 10);
    fOutputList->Add(fNumberTriggersPerEvent);
    
    fNumberElecsPerEvent= new TH1D("NumberElecsPerEvent", "NumberElecsPerEvent", 10, 0, 10);
    fOutputList->Add(fNumberElecsPerEvent);
    
    fTOFHits = new TH1D("fTOFHits", "TOF hit or no", 2, 0, 2);
    fOutputList->Add(fTOFHits);
    
    fTPCCrossedRows = new TH1D("fTPCCrossedRows", "TPC Crossed Rows", 150, 0, 150);
    fOutputList->Add(fTPCCrossedRows);
    
    fHEassNoTrig = new TH1D("fHEassNoTrig", "fHEassNoTrig", 2, 0, 2);
    fOutputList->Add(fHEassNoTrig);
    
    fHEnoAssTrig = new TH1D("fHEnoAssTrig", "fHEnoAssTrig", 2, 0, 2);
    fOutputList->Add(fHEnoAssTrig);
    
    fHHassNoTrig = new TH1D("fHHassNoTrig", "fHHassNoTrig", 2, 0, 2);
    fOutputList->Add(fHHassNoTrig);
    
    fHHnoAssTrig = new TH1D("fHHnoAssTrig", "fHHnoAssTrig", 2, 0, 2);
    fOutputList->Add(fHHnoAssTrig);
    
    fUnlikesignMassDist = new TH1D("fUnlikesignMassDist", "fUnlikesignMassDist", 1000, 0, 1);
    fOutputList->Add(fUnlikesignMassDist);
    
    fLikesignMassDist = new TH1D("fLikesignMassDist", "fLikesignMassDist", 1000, 0, 1);
    fOutputList->Add(fLikesignMassDist);
    
    fHEDeltaEtaDeltaPhi= new TH2D("HEDeltaEtaDeltaPhi", "HEDeltaEtaDeltaPhi", 1000, -10, 10, 1000, -10, 10);
    fOutputList->Add(fHEDeltaEtaDeltaPhi);
    
    
    //For calculating the Pi0 eta weight
    Int_t binw[5] =     {250,30,2,10}; //pT, PDG, EnhancedSigOrNot, pi0etaType.
        Double_t xminWt[5] = {0,0,0,-1};
        Double_t xmaxWt[5] = {50,3,2,9};

        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;EnhanceSigOrNot;pi0etaType;SPDntrCorr;",4,binw,xminWt,xmaxWt);
        fSprsPi0EtaWeightCal->GetAxis(0)->SetName("pT");
        fSprsPi0EtaWeightCal->GetAxis(1)->SetName("PDG");
        fSprsPi0EtaWeightCal->GetAxis(2)->SetName("EnhancedSigOrNot");
        fSprsPi0EtaWeightCal->GetAxis(3)->SetName("pi0etaType");
        fSprsPi0EtaWeightCal->Sumw2();
        fOutputList->Add(fSprsPi0EtaWeightCal);
    
    fRealInclsElecPt = new TH1F("fRealInclsElecPt","p_{T} distribution of MC tagged inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
    fOutputList->Add(fRealInclsElecPt);
    
    fNonHFeTrkPt = new TH1F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
       fNonHFeTrkPt->Sumw2();
       fOutputList->Add(fNonHFeTrkPt);
    
    Double_t pi = TMath::Pi();
      fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
      fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    fNonHFeEmbTrkPt = new TH1F("fNonHFeEmbTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
       fNonHFeEmbTrkPt->Sumw2();
       fOutputList->Add(fNonHFeEmbTrkPt);
    
    fPi0eEmbWeightTrkPt = new TH1F("fPi0eEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        fPi0eEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fPi0eEmbWeightTrkPt);
    
    fNonHFeEmbWeightTrkPt = new TH1F("fNonHFeEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,50);
    fNonHFeEmbWeightTrkPt->Sumw2();
    fOutputList->Add(fNonHFeEmbWeightTrkPt);
    
    fRecoNonHFeTrkPt = new TH1F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        fRecoNonHFeTrkPt->Sumw2();
        fOutputList->Add(fRecoNonHFeTrkPt);

    fEtaeEmbWeightTrkPt = new TH1F("fEtaeEmbWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
    fEtaeEmbWeightTrkPt->Sumw2();
    fOutputList->Add(fEtaeEmbWeightTrkPt);
    
    fRecoNonHFeEmbTrkPt = new TH1F("fRecoNonHFeEmbTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
     fRecoNonHFeEmbTrkPt->Sumw2();
     fOutputList->Add(fRecoNonHFeEmbTrkPt);
    
    fRecoPi0eEmbWeightTrkPt = new TH1F("fRecoPi0eEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
    fRecoPi0eEmbWeightTrkPt->Sumw2();
    fOutputList->Add(fRecoPi0eEmbWeightTrkPt);
    
    fRecoNonHFeEmbWeightTrkPt = new TH1F("fRecoNonHFeEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
       fRecoNonHFeEmbWeightTrkPt->Sumw2();
       fOutputList->Add(fRecoNonHFeEmbWeightTrkPt);
    
    fRecoEtaeEmbWeightTrkPt = new TH1F("fRecoEtaeEmbWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
      fRecoEtaeEmbWeightTrkPt->Sumw2();
      fOutputList->Add(fRecoEtaeEmbWeightTrkPt);
    
    //add variables to the output list so they show up in the file
    TNamed *etaCut = new TNamed(Form("etaCut = %f", fEtaCut), Form("%f",fEtaCut));
    fOutputList->Add(etaCut);

    TNamed *invariantMassCut = new TNamed(Form("invariantMassCut = %f", fInvariantMassCut), Form("%f",fInvariantMassCut));
    fOutputList->Add(invariantMassCut);

    TNamed *partnerElectronPtCut = new TNamed(Form("partnerElectronPtCut = %f", fPartnerElectronPtCut), Form("%f",fPartnerElectronPtCut));
    fOutputList->Add(partnerElectronPtCut);

    TNamed *TPCNSigmaElecMin = new TNamed(Form("TPCNSigmaElecMin = %f", fTPCNSigmaElecMin), Form("%f", fTPCNSigmaElecMin));
    fOutputList->Add(TPCNSigmaElecMin);

    TNamed *TPCNSigmaElecMax = new TNamed(Form("TPCNSigmaElecMax = %f", fTPCNSigmaElecMax), Form("%f", fTPCNSigmaElecMax));
    fOutputList->Add(TPCNSigmaElecMax);

    TNamed *partnerTPCNSigmaElecMin = new TNamed(Form("partnerTPCNSigmaElecMin = %f", fPartnerTPCNSigmaElecMin),Form("%f", fPartnerTPCNSigmaElecMin));
    fOutputList->Add(partnerTPCNSigmaElecMin);

    TNamed *partnerTPCNSigmaElecMax = new TNamed(Form("partnerTPCNSigmaElecMax = %f", fPartnerTPCNSigmaElecMax),Form("%f", fPartnerTPCNSigmaElecMax));
    fOutputList->Add(partnerTPCNSigmaElecMax);

    TNamed *assocTrkBitElec = new TNamed(Form("assocTrkBitElec = %f", fAssocTrkBitElec), Form("%f", fAssocTrkBitElec));
    fOutputList->Add(assocTrkBitElec);

    TNamed *assocTrkBitH = new TNamed(Form("assocTrkBitH = %f", fAssocTrkBitH), Form("%f", fAssocTrkBitH));
    fOutputList->Add(assocTrkBitH);

    TNamed *partnerTrkBit = new TNamed(Form("partnerTrkBit = %f", fPartnerTrkBit), Form("%f", fPartnerTrkBit));
    fOutputList->Add(partnerTrkBit);

    TNamed *trigTrkBit = new TNamed(Form("trigTrkBit = %f", fTrigTrkBit), Form("%f", fTrigTrkBit));
    fOutputList->Add(trigTrkBit);

    TNamed *useSPDKAny = new TNamed(Form("useSPDKAny = %d", fUseSPDKAny),Form("%d", fUseSPDKAny));
    fOutputList->Add(useSPDKAny);
    
    TNamed *useSPDKBoth = new TNamed(Form("useSPDKBoth = %d", fUseSPDKBoth),Form("%d", fUseSPDKBoth));
    fOutputList->Add(useSPDKBoth);

    TNamed *runOnMC = new TNamed(Form("runOnMC = %d", fRunOnMC), Form("%d", fRunOnMC));
    fOutputList->Add(runOnMC);

    TNamed *useTOFCut = new TNamed(Form("useTOFCut = %d", fUseTOFCut), Form("%d", fUseTOFCut));
    fOutputList->Add(useTOFCut);

    PostData(1, fOutputList);

}

//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, double multPercentile)
{

    double dist_points[5]; //Pt, Phi, Eta, zVtx, multiplicity percentile
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        dist_points[4] = multPercentile;
        fDist->Fill(dist_points);
    }

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::FillSingleParticleDist(std::vector<AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist, double multPercentile)
{

    double dist_points[5]; //Pt, Phi, Eta, zVtx, multiplicity percentile
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = zVtx;
        dist_points[4] = multPercentile;
        fDist->Fill(dist_points);
    }

}
//_____________________________________________________________________________
AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer AliAnalysisTaskHadronElectronAnalysis::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer mom;
    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}
//_____________________________________________________________________________
AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer AliAnalysisTaskHadronElectronAnalysis::RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle)
{

    AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer mom;
    // Rotating track1
    TVector3 track1Vector(track1->Px(), track1->Py(), track1->Pz());
    track1Vector.RotateZ(angle);
    mom.particle.SetPx(track1Vector(0) + track2->Px());
    mom.particle.SetPy(track1Vector(1) + track2->Py());
    mom.particle.SetPz(track1Vector(2) + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::LoadEfficiencies() {

    TFile* effFile = TFile::Open(EFF_FILE_PATH);

    if(!effFile) {
        AliFatal("NULL INPUT FILE WHEN LOADING EFFICIENCIES, EXITING");
    }
    
    fAssociatedEff = (TH1D*) effFile->Get("fAssociatedEff")->Clone("fAssociatedEffClone");
    if(!fAssociatedEff) {
        AliFatal("UNABLE TO FIND ASSOCIATED EFF, EXITING");
    }

    fTriggerEff = (TH1D*) effFile->Get("fTriggerEff")->Clone("fTriggerEffClone");
    if(!fTriggerEff) {
        AliFatal("UNABLE TO FIND TRIGGER EFF, EXITING");
    }

}
//_____________________________________________________________________________
AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer AliAnalysisTaskHadronElectronAnalysis::FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskHadronElectronAnalysis::AliMotherContainer mom;
    // Flipping track1
    mom.particle.SetPx(-track1->Px() + track2->Px());
    mom.particle.SetPy(-track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::MakeSameHElecCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> elec_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff)
{

    double dphi_point[7];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();
        
       

        for(int i = 0; i < (int)elec_list.size(); i++) {

            if(elec_list[i]->GetID() == trigger->GetID()){
                continue;
            }
            dphi_point[1] = elec_list[i]->Pt();
            dphi_point[2] = trigger->Phi() - elec_list[i]->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - elec_list[i]->Eta();
            dphi_point[4] = elec_list[i]->M();
            dphi_point[5] = zVtx;
            dphi_point[6] = multPercentile;
            
            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5)
                               && (elec_list[i]->Pt() < 10 && elec_list[i]->Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int associatedBin = fAssociatedEff->FindBin(elec_list[i]->Pt());
                double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;
                double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point, totalScale);
                
                fTriggEffVsPt->Fill(trigEff,trigger->Pt());
                fAssoEffVsPt->Fill(associatedEff,elec_list[i]->Pt());
                

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff)
{

    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_h_list.size(); i++) {
            auto associate = associated_h_list[i];
           
            if(associate->GetID() == trigger->GetID()){
             continue;//NOT 'return' because this method is being called for the entire trigger/associated lists from a given event
            }
           
            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;
            dphi_point[5] = multPercentile;
            
            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5)
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));
            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;
                double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff)
{

    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = j+1; i < (int)trigger_list.size(); i++) {
            auto associate = trigger_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;
            dphi_point[5] = multPercentile;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5)
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));

            if(eff && in_pt_range) {
                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int associatedBin = fTriggerEff->FindBin(associate->Pt());
                double associatedEff = fTriggerEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;
                double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point, totalScale);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::MakeMixedHElecCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> corrElec_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff)
{

    double dphi_point[7];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)corrElec_list.size(); j++) {
                auto electron = corrElec_list[j];

                dphi_point[1] = electron->Pt();
                dphi_point[2] = trigger->Phi() - electron->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - electron->Eta();
                dphi_point[4] = electron->M();
                dphi_point[5] = zVtx;
                dphi_point[6] = multPercentile;
                
                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5)
                                && (electron->Pt() < 10 && electron->Pt() > 0.5));
                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    fTriggerEffElec->Fill(trigger->Pt(),trigEff);
                    
                    int associatedBin = fAssociatedEff->FindBin(electron->Pt());
                    double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                    double associatedScale = 1.0/associatedEff;
                    fAssociatedEffElec->Fill(electron->Pt(),associatedEff);
                    
                    double totalScale = triggerScale*associatedScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, double multPercentile, bool eff)
{

    double dphi_point[6];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_h_list.size(); j++) {
                auto associate = associated_h_list[j];

                dphi_point[1] = associate->Pt();
                dphi_point[2] = trigger->Phi() - associate->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - associate->Eta();
                dphi_point[4] = zVtx;
                dphi_point[5] = multPercentile;

                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5)
                                && (associate->Pt() < 10 && associate->Pt() > 0.5));

                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    fTriggerEffH->Fill(trigger->Pt(),trigEff);
                    
                    int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                    double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                    double associatedScale = 1.0/associatedEff;
                    fAssociatedEffH->Fill(associate->Pt(),associatedEff);
                    
                    double totalScale = triggerScale*associatedScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }

}

//_____________________________________________________________________________
bool AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsHadrons(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);

    pass = pass && track->TestFilterMask(fAssocTrkBitH);

    return pass;
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsElectronsWoTOFWSPDKBoth(AliAODTrack *track, std::string analysisType){
  Bool_t pass = kTRUE;
if(fRunOnMC == kFALSE){//remove this cut when running on MC
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);
    }
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fAssocTrkBitElec);//
   
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
 //SPDKBoth Cut
 if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return;
 if (pass){
  if(analysisType == "he"){
   fNSigmaAssoNoCuts_woTOFwSPDKBoth->Fill(track->P(), TPCNSigmaElec);
    pass = pass && (TPCNSigmaElec >= fTPCNSigmaElecMin && TPCNSigmaElec <= fTPCNSigmaElecMax);
    if (pass){
     fNSigmaAssoCuts_woTOFwSPDKBoth->Fill(track->P(),TPCNSigmaElec);
     fNSigmaAssoCutsPt_woTOFwSPDKBoth->Fill(track->Pt(),TPCNSigmaElec);
    }
  }
}
}
void AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsElectronsWoTOFWSPDKAny(AliAODTrack *track, std::string analysisType){
    Bool_t pass = kTRUE;
 if(fRunOnMC == kFALSE){//remove this cut when running on MC
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);
    }
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fAssocTrkBitElec);//
  //std::cout << "Made it to line 1179" << std::endl;  
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
 //SPDKAny Cut
  //std::cout << "Made it to line 1182" << std::endl; 
 if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return;
 // std::cout << "Made it to line 1184" << std::endl; 
 if (pass){
  if(analysisType == "he"){
   fNSigmaAssoNoCuts_woTOFwSPDKAny->Fill(track->P(), TPCNSigmaElec);
  //  std::cout << "Made it to line 1188" << std::endl; 
    pass = pass && (TPCNSigmaElec >= fTPCNSigmaElecMin && TPCNSigmaElec <= fTPCNSigmaElecMax);
    //    std::cout << "Made it to line 1190" << std::endl; 
    if (pass){
     // std::cout << "Made it to line 1192" << std::endl; 
     fNSigmaAssoCuts_woTOFwSPDKAny->Fill(track->P(),TPCNSigmaElec);
     //std::cout << "Made it to line 1194" << std::endl; 
     fNSigmaAssoCutsPt_woTOFwSPDKAny->Fill(track->Pt(),TPCNSigmaElec);
    // std::cout << "Made it to line 1196" << std::endl; 
    }
  }
} 
}
void AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsElectronsWTOFWSPDKBoth(AliAODTrack *track, std::string analysisType){
    Bool_t pass = kTRUE;
if(fRunOnMC == kFALSE){//remove this cut when running on MC
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);
    }
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fAssocTrkBitElec);//
   
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
 //TOF Cut
    double TOFNSigmaElec = -999;
     TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
     pass = pass && (TMath::Abs(TOFNSigmaElec) <= 3);
     
 //SPDKBoth Cut
 if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return;
 if (pass){
  if(analysisType == "he"){
   fNSigmaAssoNoCuts_wTOFwSPDKBoth->Fill(track->P(), TPCNSigmaElec);
    pass = pass && (TPCNSigmaElec >= fTPCNSigmaElecMin && TPCNSigmaElec <= fTPCNSigmaElecMax);
    if (pass){
     fNSigmaAssoNoCuts_wTOFwSPDKBoth->Fill(track->P(),TPCNSigmaElec);
     fNSigmaAssoCutsPt_wTOFwSPDKBoth->Fill(track->Pt(),TPCNSigmaElec);
    }
  }
}
}
void AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsElectronsWTOFWSPDKAny(AliAODTrack *track, std::string analysisType){
    Bool_t pass = kTRUE;
 if(fRunOnMC == kFALSE){//remove this cut when running on MC
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);
    }
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fAssocTrkBitElec);//
   
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
 //TOF Cut
    double TOFNSigmaElec = -999;
     TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
     pass = pass && (TMath::Abs(TOFNSigmaElec) <= 3);
      
 //SPDKAny Cut
 if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return;
 if (pass){
  if(analysisType == "he"){
   fNSigmaAssoNoCuts_wTOFwSPDKAny->Fill(track->P(), TPCNSigmaElec);
    pass = pass && (TPCNSigmaElec >= fTPCNSigmaElecMin && TPCNSigmaElec <= fTPCNSigmaElecMax);
    if (pass){
     fNSigmaAssoCuts_wTOFwSPDKAny->Fill(track->P(),TPCNSigmaElec);
     fNSigmaAssoCutsPt_wTOFwSPDKAny->Fill(track->Pt(),TPCNSigmaElec);
    }
  }
} 
}

bool AliAnalysisTaskHadronElectronAnalysis::PassAssociatedCutsElectrons(AliAODTrack *track, std::string analysisType){

    Bool_t pass = kTRUE;
    
    // if(PassTPCCuts(track) && (TMath::Abs(TOFNSigmaElec) <= 2 || TOFNSigmaElec == 1000)) {
    // if(track->Charge() == 1){//charge is 0 if its negative?
     //working with both electrons and positrons because I'm not sorting by charge
    
    if(fRunOnMC == kFALSE){//remove this cut when running on MC
    pass = pass && (track->Pt() >= 2 && track->Pt() <= 4);
    }
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fAssocTrkBitElec);//
   //std::cout << "Made it to line 1258" << std::endl; 
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
     //For Pions, +- half a sigma around the maximum
    //pass = pass && (TPCNSigmaElec >= -5.5 &&TPCNSigmaElec >= -4.5);
  //std::cout << "Made it to line 1262" << std::endl; 
        //SPD Cut
    if(fUseSPDKAny == kTRUE){
    if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return kFALSE;
    }
  // std::cout << "Made it to line 1267" << std::endl;  
    if(fUseSPDKBoth == kTRUE){
    if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return kFALSE; //standard to use &, unless you are using emcal.
    }
     //  std::cout << "Made it to line 1271" << std::endl; 
    //TOF Cut
    double TOFNSigmaElec = -999;
     TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    if(fUseTOFCut == kTRUE){
     pass = pass && (TMath::Abs(TOFNSigmaElec) <= 3);
    }
//std::cout << "Made it to line 1278" << std::endl; 
    //For All Hadrons, -4 -> -10
    if(analysisType == "hhShape"){
     pass = pass && (TPCNSigmaElec <= -4 && TPCNSigmaElec >= -10);
    }
   // std::cout << "Made it to line 1283" << std::endl; 
    //For electrons
    if(analysisType == "he"){
   
     // std::cout << "The value of pass is " << pass << std::endl;
     // std::cout << "The value of TPCNSigmaElec is " << TPCNSigmaElec << std::endl;

    if (pass) fNSigmaAssoNoCuts->Fill(track->P(), TPCNSigmaElec);
   // std::cout << "Made it to line 1288" << std::endl;  
    pass = pass && (TPCNSigmaElec >= fTPCNSigmaElecMin && TPCNSigmaElec <= fTPCNSigmaElecMax);// cut for when running on electrons
     //   std::cout << "Made it to line 1290" << std::endl; 
    if (pass){
     //W/ TOF cut and w/ SPDKBoth
     fNSigmaAssoCuts->Fill(track->P(),TPCNSigmaElec);
     fNSigmaAssoCutsPt->Fill(track->Pt(),TPCNSigmaElec);
    }
    }
  //std::cout << "Made it to line 1297" << std::endl; 

    return pass;
}
//_____________________________________________________________________________
bool AliAnalysisTaskHadronElectronAnalysis::PassPartnerCutsElectrons(AliAODTrack *track){

    Bool_t pass = kTRUE;
    pass = pass && (track->Pt() >= fPartnerElectronPtCut);
    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
    pass = pass && track->TestFilterMask(fPartnerTrkBit);//
    // pass = pass && track->TestBit(TRIG_TRK_BIT);//very important to use TestBit. second generation filter mask thingy
    
    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    pass = pass && (TPCNSigmaElec >= fPartnerTPCNSigmaElecMin && TPCNSigmaElec <= fPartnerTPCNSigmaElecMax);
    
//std::cout << "DID WE PASS? : " << pass << std::endl;
    return pass;
}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskHadronElectronAnalysis::PassTriggerCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= fEtaCut);
  //  pass = pass && (track->Pt() >= 1 &&track->Pt() <= 12);
    pass = pass && (track->Pt() >= 4 &&track->Pt() <= 8);// now the tirgger list contains nothign but actual triggers and you don't have to cut offline. but now the number of triggers per event is done online, so I had to make this change.

    pass = pass && track->TestBit(fTrigTrkBit);

    return pass;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskHadronElectronAnalysis::PassTPCCuts(AliAODTrack *track){//cuts in PassAssociatedCutsElectron are stronger than this so its not really doing anything.
    Bool_t pass = kTRUE;
//not sure about the greater than equal to
    fpidResponse = fInputHandler->GetPIDResponse();

    double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    double TOFNSigmaElec=1000;
    TOFNSigmaElec=fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    
   //pass = pass && (TPCNSigmaElec >= -0.90);###
    pass = pass && (TPCNSigmaElec >= -2);
  
    return pass;
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::Make2DPtDistributionHistogram(std::vector<AliAODTrack*> trigger_list,std::vector<AliAODTrack*> corrElec_list, std::vector<AliAODTrack*> associated_h_list){
    //2D HE pt dist
    for(int trigg_track = 0; trigg_track < trigger_list.size();trigg_track++) {
        for(int ass_track = 0; ass_track < corrElec_list.size(); ass_track++) {
           f2DHEpTDist->Fill(trigger_list[trigg_track]->Pt(),corrElec_list[ass_track]->Pt());
            fHEassNoTrig->Fill(0);
            fHEnoAssTrig->Fill(0);
        }
    }
    //Now for exceptions, ie if there are no triggers or no associated particles
    if(trigger_list.size()==0 && corrElec_list.size()!=0){
        for(int ass_track = 0; ass_track < corrElec_list.size(); ass_track++) {
            f2DHEpTDist->Fill(0.0,corrElec_list[ass_track]->Pt());
            fHEassNoTrig->Fill(1);
        }
    }
    
    if(trigger_list.size()!=0 && corrElec_list.size()==0){
        for(int trigg_track = 0; trigg_track < trigger_list.size();trigg_track++) {
            f2DHEpTDist->Fill(trigger_list[trigg_track]->Pt(),0);//this line is messing things up
            fHEnoAssTrig->Fill(1);
        }
    }
    //2D HH pt dist
    for(int trigg_track = 0; trigg_track < trigger_list.size();trigg_track++) {
        for(int ass_track = 0; ass_track < associated_h_list.size(); ass_track++) {
           f2DHHpTDist->Fill(trigger_list[trigg_track]->Pt(),associated_h_list[ass_track]->Pt());
        
            fHHassNoTrig->Fill(0);
            fHHnoAssTrig->Fill(0);
        }
    }
    //Now for exceptions, ie if there are no triggers or no associated particles
    if(trigger_list.size()==0 && associated_h_list.size()!=0){
        for(int ass_track = 0; ass_track < associated_h_list.size(); ass_track++) {
            f2DHHpTDist->Fill(0.0,associated_h_list[ass_track]->Pt());
            fHHassNoTrig->Fill(1);
        }
    }
    
    if(trigger_list.size()!=0 && associated_h_list.size()==0){
        for(int trigg_track = 0; trigg_track < trigger_list.size();trigg_track++) {
            f2DHHpTDist->Fill(trigger_list[trigg_track]->Pt(),0);
            fHHnoAssTrig->Fill(1);
        }
    }
}

//_____________________________________________________________________________

std::vector<AliAODTrack*> AliAnalysisTaskHadronElectronAnalysis::GetUnlikeSignVector(std::vector<AliAODTrack*> corrElec_list,std::vector<AliAODTrack*> partnerElec_list) {
    std::vector<AliAODTrack*> unlikeSign;
    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t flagPhotonicElec = kFALSE;//not sure if this line belongs here
    Bool_t EffiDenom = fRunOnMC;
    Bool_t EffiNumTag = fRunOnMC;
    for(int i = 0; i < (int) corrElec_list.size(); i++) {
//cout<<"i US: "<< i<<" \n";
        //Monte Carlo Stuff
        if(fRunOnMC==true){
           
            fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
            fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
            if(!fMCarray){
                AliError("MC Array not found.");
                continue;
            }
            if(!fMCHeader){
                printf("Header is empty!\n");
              
            }
            //set the pi0 and eta weight
            fPi0Weight->SetParameters( 9.77213e+02,-1.74570e-01, 2.36061e-03,1.55257e+00,4.60931e+00);
            fEtaWeight->SetParameters( 7.32301e+02,-1.28601e-01,1.07473e-03,1.42236e+00,4.65662e+00);
        //Getting number of particles produced by the MC generator
            if(fIsNew) GetNMCPartProducedpPbNew();
             else GetNMCPartProducedpPbOld();
            /////Remove in bunch pileup events in MC//// only need for newer MC sets
          //  if(AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCHeader)) return;
            
              //Pi0 and Eta weight cal//
          if(fIsNew)  GetPi0EtaWeightpPbNew(fSprsPi0EtaWeightCal);
          else GetPi0EtaWeightpPbOld(fSprsPi0EtaWeightCal);
            //tagging efficiency
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            //Bool_t EffiDenom=MC;
           if(fMCHeader){
           if(fIsNew)  EffiDenom = GetNonHFEEffiDenompPbNew(corrElec_list[i]);
           else EffiDenom = GetNonHFEEffiDenompPbOld(corrElec_list[i]);
             //numerator is calculated in the unlikesign method
            }
        }//end MC
     // cout<<"corrElec_list.size() US "<<corrElec_list.size()<<"\n";
        for(int j = 0; j < (int) partnerElec_list.size(); j++) {
         //  cout<<"j US: "<< j<<" \n";
         //   cout<<"UNLIKE SIGN ELECTRON LIST SIZE IN METHOD: "<< unlikeSign.size()<<" \n";
            if(corrElec_list[i]->GetID()==partnerElec_list[j]->GetID()) continue;
            if(corrElec_list[i]->Charge() != partnerElec_list[j]->Charge()){//1 if positive, 0 if negative
                fVevent = dynamic_cast<AliVEvent*>(InputEvent());
                AliKFParticle::SetField(fVevent->GetMagneticField());
               
                Int_t chargeAsso = partnerElec_list[j]->Charge();
                Int_t charge = corrElec_list[i]->Charge();
                
                Int_t fPDGe1 = 11;
                Int_t fPDGe2 = 11;
                if(charge>0) fPDGe1 = -11;
                    if(chargeAsso>0) fPDGe2 = -11;
                
                AliKFParticle ge1 = AliKFParticle(*corrElec_list[i], fPDGe1);
                AliKFParticle ge2 = AliKFParticle(*partnerElec_list[j], fPDGe2);
                AliKFParticle recg(ge1, ge2);

                Double_t mass=-999., width = -999;
                Int_t MassCorrect;
                MassCorrect = recg.GetMass(mass,width);
                fUnlikesignMassDist->Fill(mass);
        if(mass<fInvariantMassCut){ //if they have invariant mass less than 140 MeV usually
                   
            unlikeSign.push_back(corrElec_list[i]);//add the first electron to the list of unlikesign electrons
         //if MC=true statement here
            if (fRunOnMC==true){
            flagPhotonicElec = kTRUE;
            fFlagPhotonicElec = flagPhotonicElec;
            
            if(fFlagPhotonicElec){
            EffiNumTag = GetNonHFEEffiRecoTag(corrElec_list[i]);//corrElec_list[i] instead of track because its baisically the same thing
            }
            }
          
                //cout<<"UNLIKE SIGN ELECTRON LIST SIZE IN METHOD: "<< unlikeSign.size()<<" \n";
              }//end if statement
              
              
            }//end if statement
        }//end second for loop
    }//end first for loop
    return unlikeSign;
}

std::vector<AliAODTrack*> AliAnalysisTaskHadronElectronAnalysis::GetLikeSignVector(std::vector<AliAODTrack*> corrElec_list,std::vector<AliAODTrack*> partnerElec_list) {
    
    std::vector<AliAODTrack*> likeSign;
    
    for(int i = 0; i < (int) corrElec_list.size(); i++) {
       // cout<<"i LS: "<< i<<" \n";
      //  cout<<"corrElec_list.size() LS "<<corrElec_list.size()<<"\n";
        for(int j = 0; j < (int) partnerElec_list.size(); j++) {
       //     cout<<"j LS: "<< j<<" \n";
            
            
            if(corrElec_list[i]->GetID()== partnerElec_list[j]->GetID()) continue;
            
            if(corrElec_list[i]->Charge() == partnerElec_list[j]->Charge()){
               // std::cout << corrElec_list[i]->GetID() << " : " << partnerElec_list[j]->GetID() << std::endl;
                //std::cout << "\t" << corrElec_list[i]->Pt() << " : " << partnerElec_list[j]->Pt() << std::endl;
//1 if positive, 0 if negative
                fVevent = dynamic_cast<AliVEvent*>(InputEvent());
                AliKFParticle::SetField(fVevent->GetMagneticField());
                
                Int_t chargeAsso = partnerElec_list[j]->Charge();
                Int_t charge = corrElec_list[i]->Charge();
                
                Int_t fPDGe1 = 11;
                Int_t fPDGe2 = 11;
                if(charge>0) fPDGe1 = -11;
                    if(chargeAsso>0) fPDGe2 = -11;
                
                AliKFParticle ge1 = AliKFParticle(*corrElec_list[i], fPDGe1);
                AliKFParticle ge2 = AliKFParticle(*partnerElec_list[j], fPDGe2);
                AliKFParticle recg(ge1, ge2);

                Double_t mass=-999., width = -999;
                Int_t MassCorrect;
                MassCorrect = recg.GetMass(mass,width);
                fLikesignMassDist->Fill(mass);
        if(mass<fInvariantMassCut){ //if they have invariant mass less than 140 MeV
                   
                    likeSign.push_back(corrElec_list[i]);//add the first electron to the list of likesign electrons
                  // cout<<"LIKE SIGN ELECTRON LIST SIZE IN METHOD: " << likeSign.size()<<" \n";
             }//end if statement
            }//end if statement
        }//end second for loop
    }//end first for loop
    
    
    return likeSign;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskHadronElectronAnalysis::GetNMCPartProducedpPbNew()
{
 //Get number of MC particles produced by generators.

  TList *lh = fMCHeader->GetCocktailHeaders();
  fNTotMCpart = 0;
  fNembMCpi0 = 0;
  fNembMCeta = 0;
  fNpureMC = 0;

  TString MCgen;
  TString embpi0("pi");
  TString embeta("eta");
  TString mb("Hijing");

  if(!lh){
    AliError("no MC header");
    return (0);
  }

for(int igene=0; igene<lh->GetEntries(); igene++)
  {
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
    if(!gh) continue;

    MCgen =  gh->GetName();

    //Order of the generators added- cele or bele, pi0, eta, MB
    if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
    if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
    if(MCgen.Contains(mb))    fNpureMC = fNTotMCpart;

    fNTotMCpart += gh->NProduced();
  }
  return kTRUE;
}
//_________________________
Bool_t AliAnalysisTaskHadronElectronAnalysis::GetNMCPartProducedpPbOld()
{
  //Get number of MC particles produced by generators.
//hijing is at end of the list as opposed to the beginning at NEW,
    //run locally in test mode and use cout, how the mc particles are ordered.because it could be hijing first and embedded pi0 and eta

    TList *lh = fMCHeader->GetCocktailHeaders();
  
   // cout<<"MADE IT HERE!!!!!"<<endl;
    fNTotMCpart = 0;
  fNembMCpi0 = 0;
  fNembMCeta = 0;
  fNpureMC = 0;

      TString MCgen;
  TString embpi0("pi");
  TString embeta("eta");

  if(!lh){
    AliError("no MC header");
    return (0);
  }

  for(int igene=0; igene<lh->GetEntries(); igene++)
  {
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
    if(!gh) continue;

    MCgen =  gh->GetName();

    //Order of the generators added- MB, pi0, eta, cele or bele
    if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING

    if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
    if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;

    fNTotMCpart += gh->NProduced();
  }
  return kTRUE;
}
//_______________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::GetPi0EtaWeightpPbNew(THnSparse *SparseWeight)
{
//Get pi0 and eta information for weight calculation

Double_t fvalue[4] = {-999,-999,-999,-999};

for(int imc=0; imc< fNTotMCpart; imc++)
{

AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(imc);

if(TMath::Abs(AODMCtrack->Eta()) > fEtaCut) continue;

//-------Get PDG-------//
Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;

Double_t fPartPDGid = -999;
if (TrackPDG == 111) fPartPDGid = 0.2;
if (TrackPDG == 221) fPartPDGid = 1.2;
if (TrackPDG == 22) fPartPDGid = 2.2;

Double_t fTrkPt = AODMCtrack->Pt();

//-------Check if the particle is from Enhanced signal or not
Bool_t fFromEnhance = kEnhance;
if(imc > fNpureMC) fFromEnhance = kMB;

//------Get type of the particle
Int_t fType = GetPi0EtaType(AODMCtrack);

fvalue[0] = fTrkPt;
fvalue[1] = fPartPDGid;
fvalue[2] = fFromEnhance;
fvalue[3] = fType;

SparseWeight->Fill(fvalue);
}
}
//_______________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::GetPi0EtaWeightpPbOld(THnSparse *SparseWeight)
{
  //Get pi0 and eta information for weight calculation

  Double_t fvalue[4] = {-999,-999,-999,-999};

  for(int imc=0; imc< fNTotMCpart; imc++)
  {

    AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(imc);

    if(TMath::Abs(AODMCtrack->Eta()) > fEtaCut) continue;

    //-------Get PDG
    Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
    if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;

    Double_t fPartPDGid = -999;
    if (TrackPDG == 111) fPartPDGid = 0.2;
    if (TrackPDG == 221) fPartPDGid = 1.2;
    if (TrackPDG == 22) fPartPDGid = 2.2;

    Double_t fTrkPt = AODMCtrack->Pt();

    //-------Check if the particle is from Enhanced signal or not
    Bool_t fFromEnhance = kMB;
    if(imc >= fNpureMC)fFromEnhance = kEnhance;

    //------Get type of the particle
    Int_t fType = GetPi0EtaType(AODMCtrack);

    fvalue[0] = fTrkPt;
    fvalue[1] = fPartPDGid;
    fvalue[2] = fFromEnhance;
    fvalue[3] = fType;

    SparseWeight->Fill(fvalue);
  }
}
//_____________________________________________
Int_t AliAnalysisTaskHadronElectronAnalysis::GetPi0EtaType(AliAODMCParticle *part)
{
  // Return the type of particle

  // IsPrimary
  Bool_t primMC = part->IsPrimary();
  if(!primMC) return kNotIsPrimary;

  // Mother
  Int_t motherlabel = part->GetMother();
  if(motherlabel<0) return kNoMother;

  else {
    AliAODMCParticle *mother = (AliAODMCParticle*)fMCarray->At(motherlabel);
    Int_t motherpdg = TMath::Abs(mother->GetPdgCode());

    if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;

    if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
    if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
    return kNoFeedDown;
  }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskHadronElectronAnalysis::GetNonHFEEffiDenompPbNew(AliVTrack *track)
{
  //Calculate Non-HFE efficiency demoninator

  fIsFrmEmbPi0 = kFALSE, fIsFrmEmbEta = kFALSE;
  ftype = -1, fWeightPi0 = 1.0, fWeightEta = 1.0, fWeight=1.0;
  Bool_t fFromMB = kTRUE;

  Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
  Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
  Double_t MomPt =-999.0;

  AliAODMCParticle *MCPart = 0;
  AliAODMCParticle *MCPartMom = 0;
  AliAODMCParticle *MCPartGMom = 0;
  AliAODMCParticle *MCPartGGMom = 0;
  AliAODMCParticle *MCPartGGGMom = 0;

  Double_t TrkPt = track->Pt();
  Int_t iTrklabel = TMath::Abs(track->GetLabel());
  if(iTrklabel == 0) return kFALSE;

  MCPart = (AliAODMCParticle*)fMCarray->At(iTrklabel);
  if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
  fRealInclsElecPt->Fill(TrkPt);

  Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
  if(!fNonHFE) return kFALSE;
  fNonHFeTrkPt->Fill(TrkPt);

  MCPartMom = (AliAODMCParticle*)fMCarray->At(iMCmom);
  iMCgmom = MCPartMom->GetMother();
  if(iMCgmom > 0){
    MCPartGMom = (AliAODMCParticle*)fMCarray->At(iMCgmom);
    GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());

    iMCggmom = MCPartGMom->GetMother();
    if(iMCggmom > 0){
      MCPartGGMom = (AliAODMCParticle*)fMCarray->At(iMCggmom);
      GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());

      iMCgggmom = MCPartGGMom->GetMother();
      if(iMCgggmom > 0){
        MCPartGGGMom = (AliAODMCParticle*)fMCarray->At(iMCgggmom);
        GGGMomPDG = TMath::Abs(MCPartGGGMom->GetPdgCode());
      }
    }
  }

  //cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e
  if(MomPDG == 221){
    if(iMCmom >= fNembMCeta && iMCmom < fNpureMC) { //from eta event
      fIsFrmEmbEta = kTRUE; //eta->e
      fWeightEta = fEtaWeight->Eval(MCPartMom->Pt());
    }
  }

  if(MomPDG == 111) {
    if(iMCmom >= fNembMCpi0 && iMCmom < fNembMCeta){ //from pi0 event
      fIsFrmEmbPi0 = kTRUE; //pi0 -> e
      fWeightPi0 = fPi0Weight->Eval(MCPartMom->Pt());
    }

    if(GMomPDG == 221){
      if(iMCgmom >= fNembMCeta && iMCgmom < fNpureMC) { //from eta event
        fIsFrmEmbEta = kTRUE; //eta->pi0-> e
        fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
      }
    }
  }

  if(MomPDG == 22){
    if(GMomPDG == 221){
      if(iMCgmom >= fNembMCeta && iMCgmom < fNpureMC) { //from eta event
        fIsFrmEmbEta = kTRUE; //eta->gamma-> e
        fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
      }
    }

    if(GMomPDG == 111){
      if(iMCgmom >= fNembMCpi0 && iMCgmom < fNembMCeta) { //from pi0 event
        fIsFrmEmbPi0 = kTRUE; //pi0-> gamma-> e
        fWeightPi0 = fPi0Weight->Eval(MCPartGMom->Pt());
      }

      if(GGMomPDG == 221){
        if(iMCggmom >= fNembMCeta && iMCggmom < fNpureMC) { //from eta event
          fIsFrmEmbEta = kTRUE; //eta->pi0->gamma-> e
          fWeightEta = fEtaWeight->Eval(MCPartGGMom->Pt());
        }
      }
    }
  }

  if(fIsFrmEmbPi0 || fIsFrmEmbEta){
    fNonHFeEmbTrkPt->Fill(TrkPt);

    if(fIsFrmEmbPi0) {
      fWeight = fWeightPi0;
      fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
      fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
    }
    if(fIsFrmEmbEta){
      fWeight = fWeightEta;
      fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
      fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskHadronElectronAnalysis::GetNonHFEEffiDenompPbOld(AliVTrack *track)
{
  //Calculate Non-HFE efficiency demoninator

  fIsFrmEmbPi0 = kFALSE, fIsFrmEmbEta = kFALSE;
  ftype = -1, fWeightPi0 = 1.0, fWeightEta = 1.0, fWeight=1.0;
  Bool_t fFromMB = kTRUE;

  Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
  Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
  Double_t MomPt =-999.0;

  AliAODMCParticle *MCPart = 0;
  AliAODMCParticle *MCPartMom = 0;
  AliAODMCParticle *MCPartGMom = 0;
  AliAODMCParticle *MCPartGGMom = 0;
  AliAODMCParticle *MCPartGGGMom = 0;

  Double_t TrkPt = track->Pt();
  Int_t iTrklabel = TMath::Abs(track->GetLabel());
  if(iTrklabel == 0) return kFALSE;

  MCPart = (AliAODMCParticle*)fMCarray->At(iTrklabel);
  if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
  fRealInclsElecPt->Fill(TrkPt);

  Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
  if(!fNonHFE) return kFALSE;
  fNonHFeTrkPt->Fill(TrkPt);

  MCPartMom = (AliAODMCParticle*)fMCarray->At(iMCmom);
  iMCgmom = MCPartMom->GetMother();
  if(iMCgmom > 0){
    MCPartGMom = (AliAODMCParticle*)fMCarray->At(iMCgmom);
    GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());

    iMCggmom = MCPartGMom->GetMother();
    if(iMCggmom > 0){
      MCPartGGMom = (AliAODMCParticle*)fMCarray->At(iMCggmom);
      GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());

      iMCgggmom = MCPartGGMom->GetMother();
      if(iMCgggmom > 0){
        MCPartGGGMom = (AliAODMCParticle*)fMCarray->At(iMCgggmom);
        GGGMomPDG = TMath::Abs(MCPartGGGMom->GetPdgCode());
      }
    }
  }

  //cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e
  if(MomPDG == 221){
    if(iMCmom >= fNembMCeta && iMCmom < fNTotMCpart) { //from eta event
      fIsFrmEmbEta = kTRUE; //eta->e
      fWeightEta = fEtaWeight->Eval(MCPartMom->Pt());
    }
  }

  if(MomPDG == 111) {
    if(iMCmom >= fNembMCpi0 && iMCmom < fNembMCeta){ //from pi0 event
      fIsFrmEmbPi0 = kTRUE; //pi0 -> e
      fWeightPi0 = fPi0Weight->Eval(MCPartMom->Pt());
    }

    if(GMomPDG == 221){
      if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
        fIsFrmEmbEta = kTRUE; //eta->pi0-> e
        fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
      }
    }
  }

  if(MomPDG == 22){
    if(GMomPDG == 221){
      if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
        fIsFrmEmbEta = kTRUE; //eta->gamma-> e
        fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
      }
    }

    if(GMomPDG == 111){
      if(iMCgmom >= fNembMCpi0 && iMCgmom < fNembMCeta) { //from pi0 event
        fIsFrmEmbPi0 = kTRUE; //pi0-> gamma-> e
        fWeightPi0 = fPi0Weight->Eval(MCPartGMom->Pt());
      }

      if(GGMomPDG == 221){
        if(iMCggmom >= fNembMCeta && iMCggmom < fNTotMCpart) { //from eta event
          fIsFrmEmbEta = kTRUE; //eta->pi0->gamma-> e
          fWeightEta = fEtaWeight->Eval(MCPartGGMom->Pt());
        }
      }
    }
  }

  if(fIsFrmEmbPi0 || fIsFrmEmbEta){
    fNonHFeEmbTrkPt->Fill(TrkPt);

    if(fIsFrmEmbPi0) {
      fWeight = fWeightPi0;
      fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
      fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
    }
    if(fIsFrmEmbEta){
      fWeight = fWeightEta;
      fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
      fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
    }
  }

  return kTRUE;
}
//_________________________________________


//does this only work when the MC is running as well? because the if statements would never be true otherwise
Bool_t AliAnalysisTaskHadronElectronAnalysis::GetNonHFEEffiRecoTag(AliVTrack *track)
{
  Double_t TrkPt = track->Pt();

  fRecoNonHFeTrkPt->Fill(TrkPt);
  if(fIsFrmEmbPi0 || fIsFrmEmbEta){
    fRecoNonHFeEmbTrkPt->Fill(TrkPt);

    if(fIsFrmEmbPi0) {
      fRecoPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
      fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
    }
    if(fIsFrmEmbEta){
      fRecoEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
      fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
    }
  }

  return kTRUE;
}
//___________________________________________
Bool_t  AliAnalysisTaskHadronElectronAnalysis::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
  //Is electron from pi0, eta and gamma

  iMCmom = MCPart->GetMother();
  AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCarray->At(iMCmom);
  MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
  MomPt = MCPartMom->Pt();

  if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
    if(iMCmom >= fNpureMC)fFromMB = kFALSE;
    type = GetPi0EtaType(MCPartMom);
    return kTRUE;
  }
  else return kFALSE;
}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::UserExec(Option_t*)
{

 // std::cout << "Made it to the User Exec" << std::endl;
  
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
    }


    fpidResponse = fInputHandler->GetPIDResponse();
   
    //Event cuts
    TString cent_estimator = CENT_ESTIMATOR;
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;
    
   

    if(multPercentile < MULT_LOW || multPercentile > MULT_HIGH) return;
  //  cout<<"MUltiplicity Percentile: "<<multPercentile<< endl;
    fmultPercentile->Fill(multPercentile);
    
    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    fmultPercentileCuts->Fill(multPercentile);
    
    int numTracks = fAOD->GetNumberOfTracks();

       
    std::vector<AliAODTrack*> elec_list;
    std::vector<AliAODTrack*> triggered_elec_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> all_hadron_list;
   // std::vector<AliAODTrack*> k_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;
        
        //List for comparison with cuts/filter bits
        all_hadron_list.push_back(track);

       fAllHadronsPt->Fill(track->Pt());
        
       bool isTriggeredEvent = false;
        
        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            //get phi distributions for trigger h
            fTrigPhiDist->Fill(track->Phi());
  
            fTriggerPt->Fill(track->Pt());
            
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);

            isTriggeredEvent = true;
          }//if we pass trigger cuts then add it to the list.

        
        if(PassAssociatedCutsHadrons(track)) {
            associated_h_list.push_back(track);
            fAssHadronsPt->Fill(track->Pt());
            //get phi distributions for associated h
             fAssHPhiDist ->Fill(track->Phi());
        }


            //electrons that pass your TPC/TOF cuts
          //Calculate Beta//
            double time = track->GetTOFsignal() - fpidResponse->GetTOFResponse().GetStartTime(track->P());
            double length = track->GetIntegratedLength();
          //  double v = //length/time
           double beta;
           if(time>0){
            beta = (length/time)/0.0299792;//v/(speed of light in cm/ps)
           }
           else{
            beta=-1;
             }
            fBeta->Fill(track->P(),beta);
            
            double TPCNSigmaElec = -999;
            double TOFNSigmaElec = -999;


            TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            fDedx->Fill(track->P(),track->GetTPCsignal());
            fBetaDedx->Fill(track->GetTPCsignal(),beta);
            f2DPPtHist->Fill(track->P(),track->Pt());
        
           if( TMath::Abs(TOFNSigmaElec) <= 3) {
           fNSigma->Fill(track->P(),TPCNSigmaElec);
             }
            
            if(PassTPCCuts(track)){
                
                fDedxCuts->Fill(track->P(),track->GetTPCsignal());
                fNSigmaCuts->Fill(track->P(),TPCNSigmaElec);
            }
            
            if(TOFNSigmaElec != -999){
                
                fDedxTOF->Fill(track->P(),track->GetTOFsignal());
            }
            if(TMath::Abs(TOFNSigmaElec) <= 2 ){
                
                fDedxTOFCuts->Fill(track->P(),track->GetTOFsignal());
            }
          elec_list.push_back(track);

          if(isTriggeredEvent) triggered_elec_list.push_back(track);

    }//end track loop

    //Making list of possible electrons
    std::vector<AliAODTrack*> corrElec_list;
    //Making list of possible electrons from triggered events
    std::vector<AliAODTrack*> triggered_corrElec_list;
    //Making list of possible electrons with looser cuts
    std::vector<AliAODTrack*> partnerElec_list;


    std::string analysisType[2] = {"he", "hhShape"};
      for(int analysisTypeInt = 0; analysisTypeInt<2;analysisTypeInt++){

        for(int j = 0; j < (int) elec_list.size(); j++) {
       // find the electron for the Hadron Electron correlation
            if(analysisType[analysisTypeInt]== "he"){
            fElecPtBeforeCuts->Fill(elec_list[j]->Pt());//still placing daughter cuts on it.
            }
              //std::cout << "Made it to the new methods, but not running yet" << std::endl;
            //some methods purely to make the necessary histograms for nsigma comparisons
            PassAssociatedCutsElectronsWoTOFWSPDKBoth(elec_list[j], analysisType[analysisTypeInt]);
            //  std::cout << "Made it past WoTOFWSPDKBoth" << std::endl;
            PassAssociatedCutsElectronsWoTOFWSPDKAny(elec_list[j], analysisType[analysisTypeInt]);
          // std::cout << "Made it past WoTOFWSPDKAny" << std::endl;
            PassAssociatedCutsElectronsWTOFWSPDKBoth(elec_list[j], analysisType[analysisTypeInt]);
          //  std::cout << "Made it past WTOFWSPDKBoth" << std::endl;
            PassAssociatedCutsElectronsWTOFWSPDKAny(elec_list[j], analysisType[analysisTypeInt]);
          //  std::cout << "Made it past WTOFWSPDKAny" << std::endl;

            if(PassAssociatedCutsElectrons(elec_list[j], analysisType[analysisTypeInt])){
                //        std::cout << "Made it past PassAssociatedCutsElectrons" << std::endl;  
                if(fRunOnMC==kFALSE && elec_list[j]->Pt() >= 2 && elec_list[j]->Pt() <= 4){
               corrElec_list.push_back(elec_list[j]);
                    
                    if(analysisType[analysisTypeInt]== "he"){
                    //get phi distributions for associated elec
                     fAssElecPhiDist ->Fill(elec_list[j]->Phi());
                    }
                }
             // std::cout << "Made it to line 2159" << std::endl; 
                if(fRunOnMC==kTRUE){//fill the list without the pt cuts if running over MC
                    corrElec_list.push_back(elec_list[j]);
                }
                if(analysisType[analysisTypeInt]== "he"){
                fElecPtAfterCuts->Fill(elec_list[j]->Pt());
               
               double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(elec_list[j], AliPID::kElectron);
                double TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(elec_list[j], AliPID::kElectron);
                 
                
                }
            }
        }
  if(analysisType[analysisTypeInt]== "he"){
    for(int j = 0; j < (int) elec_list.size(); j++) {
   // find the electron for the partner electron list

        if(PassPartnerCutsElectrons(elec_list[j])){
           partnerElec_list.push_back(elec_list[j]);
        }
    }
    
    for(int j = 0; j < (int) triggered_elec_list.size(); j++) {
   // find the electron for only triggered events
         if(PassAssociatedCutsElectrons(triggered_elec_list[j], analysisType[analysisTypeInt])){
                
           if(fRunOnMC==kFALSE && triggered_elec_list[j]->Pt() >= 2 && triggered_elec_list[j]->Pt() <= 4){
               triggered_corrElec_list.push_back(triggered_elec_list[j]);
                
            }  
          }

    }
// std::cout << "Made it to line 2193" << std::endl;
    
    //Find unlikesign and likesign
    std::vector<AliAODTrack*> USElec_list = GetUnlikeSignVector(corrElec_list,partnerElec_list);//unlikesign electron list
    std::vector<AliAODTrack*> LSElec_list = GetLikeSignVector(corrElec_list,partnerElec_list);//likesign electron list

    if (corrElec_list.size() > 0) {
     // std::cout << corrElec_list.size() << " is the size of corr elec\n";
    }
    // if (partnerElec_list.size() > 0) {
    //   std::cout << partnerElec_list.size() << " is the size of part elec\n";
    // }
    
    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDist, multPercentile);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist, multPercentile);
    FillSingleParticleDist(all_hadron_list, primZ, fLooseDist, multPercentile);

    // Filling all of our correlation histograms
    MakeSameHElecCorrelations(trigger_list, corrElec_list, fDphiHElec, primZ, multPercentile, false);
    MakeSameHElecCorrelations(trigger_list, corrElec_list, fDphiHElecEff, primZ, multPercentile, true);
    MakeSameHElecCorrelations(trigger_list, triggered_corrElec_list, fDphiHElecTriggered, primZ, multPercentile, false);
    MakeSameHElecCorrelations(trigger_list, triggered_corrElec_list, fDphiHElecTriggeredEff, primZ, multPercentile, true);
    MakeSameHElecCorrelations(trigger_list, USElec_list, fDphiHUSElec, primZ, multPercentile, false);
    MakeSameHElecCorrelations(trigger_list, LSElec_list, fDphiHLSElec, primZ, multPercentile, false);
    MakeSameHElecCorrelations(trigger_list, USElec_list, fDphiHUSElecEff, primZ, multPercentile, true);
    MakeSameHElecCorrelations(trigger_list, LSElec_list, fDphiHLSElecEff, primZ,multPercentile, true);
    
    
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ, multPercentile, false);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, multPercentile, true);
    MakeSameTriggerTriggerCorrelations(trigger_list, fDphiTriggerTrigger, multPercentile, primZ);
  

  

       
  if(corrElec_list.size() > 0 || associated_h_list.size() > 0) {
    // std::cout << "Made it to mixed event section" << std::endl;
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHElecCorrelations(fCorPool, corrElec_list, fDphiHElecMixed, primZ, multPercentile);
                MakeMixedHElecCorrelations(fCorPool, USElec_list, fDphiHULSElecMixed, primZ, multPercentile);
                MakeMixedHElecCorrelations(fCorPool, LSElec_list, fDphiHLSElecMixed, primZ, multPercentile);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ, multPercentile);
                MakeMixedHHCorrelations(fCorPool, trigger_list, fDphiTriggerTriggerMixed, primZ, multPercentile);
                
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }
      }//end if analysisType == "he"
      if(analysisType[analysisTypeInt]== "hhShape"){
      MakeSameHElecCorrelations(trigger_list, corrElec_list, fDphiHHCont, primZ, multPercentile, false);
      MakeSameHElecCorrelations(trigger_list, corrElec_list, fDphiHHContEff, primZ, multPercentile, true);
    }
     }//end for loop for hadron contamination shape and h-e


   
//trigger per event counter
    fNumberTriggersPerEvent->Fill(trigger_list.size());
    fNumberElecsPerEvent->Fill(corrElec_list.size());
   
   // Make2DPtDistributionHistogram(trigger_list,corrElec_list, associated_h_list);//method throws segmentation violation error

    PostData(1, fOutputList);//writes to output list

}
//_____________________________________________________________________________
void AliAnalysisTaskHadronElectronAnalysis::Terminate(Option_t *option)
{
}




