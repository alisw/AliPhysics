/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskKaon2PC:
// Description: Analysis task to calculate Two-Particle Angular Correlation Functions of Neutral and Charged Kaons
// Author: Anjaly Sasikumar Menon
// (anjaly.sasikumar.menon@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TMarker.h"
#include "TObjArray.h"
#include "TObject.h"
#include "AliLog.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TString.h"
#include <vector>
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliEventPoolManager.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliAODPid.h"
#include "AliMCEvent.h"
#include "AliMultiInputEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisAlien.h"
#include "AliEventCuts.h"
#include "AliAODVZERO.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliVParticle.h"
#include <AliHeader.h>
#include "TLorentzVector.h"
#include "AliAnalysisTaskKaon2PC.h"

class AliAnalysisTaskKaon2PC;       // This analysis class
using namespace std; // std namespace: so you can do things like 'cout'
ClassImp(AliAnalysisTaskKaon2PC)    // classimp: necessary for root
//ClassImp(AliV0XiParticle)

const Double_t Pi = TMath::Pi();

AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC() : AliAnalysisTaskSE(),  
fAOD(0),
fmcEvent(0), 
fOutputList(0), 
fPIDResponse(0),
fPidpTDependentMethod(kTRUE),
fRejectEventPileUp(kTRUE),
fMinBias(kTRUE),
fCentral(kFALSE),
fSemiCentral(kFALSE),
fData(kFALSE),
fMCTruth(kTRUE),
fMCReconstructed(kTRUE),
fMCArray(NULL),
PVx(0), PVy(0), PVz(0),
//track cuts 
fLpTCut(0.4), 
fUpTCut(0.8), 
fEtaCut(0.8), 
fSigCut(3.0),
fBit(96),
fPVzCut(8),
fCentMin(0),
fCentMax(10),
//V0 cuts
fDecayLv0Cut(8.05), 
fLpTv0Cut(0.4), 
fUpTv0Cut(0.8), 
fEtav0Cut(0.8), 
fDcaPosToPrimVtxv0Cut(0.1), 
fDcaNegToPrimVtxv0Cut(0.1), 
fEtaPosv0Cut(0.8), 
fEtaNegv0Cut(0.8), 
fCosPACut(0.99), 
fSigPosv0Cut(3.0), 
fSigNegv0Cut(3.0),
fnumOfTPCcrossedRows(70),
fTPCrowsFindableRatio(0.8),
//PID histograms
fEnergy(0),
fEnergyCuts(0),
fPID(0), 
fPIDKaon(0), 
fPIDK(0), 
fPIDKpiCut(0),
fnsigmakaon(0),
fNsigmaKaon(0),
fNsigmaTOFK(0),
fNsigmaTOFKaon(0),
fNsigmaTPCTOFK(0),
fHistTOFKch(0),
//track and event observables
fVtx(0), 
fClusters(0),
fHistPVz(0),
fHistNEvents(0),
fMCEvents(0),
fTracksCounter(0),
fMCEvents_pileup(0),
fHistNV0(0),  
fHistEta(0), 
fHistDEta(0), 
fHistPhi(0), 
fHistDPhi(0), 
fHistMult(0), 
fHistCent(0),
fHistCent_mcgen(0),
fHistTPCTracksVsClusters(0),
//single particle histograms
fHistMK0(0), 
fHistMK0Cuts(0), 
fHistKChPt(0),
fHistKChPtfullRange(0),
fHistKChPtMix(0),  
fHistK0Pt(0),
fHistK0PtfullRange(0),
fHistKPlusPt(0),
fHistKPlusPtfullRange(0),
fHistKMinusPt(0),
fHistKMinusPtfullRange(0),
fHistKChPhi(0), 
fHistK0Phii(0), 
fHistKpPhi(0), 
fHistKnPhi(0),
fHistPPionPhi(0), 
fHistNPionPhi(0),
//single particle 2D histograms
fHistK0Phi(0),
fHistK0Eta(0),  
fHistChPhi(0), 
fHistChEta(0), 
fHistChRap(0),
fHistPosPhi(0),
fHistPosEta(0), 
fHistPosRap(0), 
fHistNegPhi(0), 
fHistNegEta(0), 
fHistNegRap(0),
//single particle 3D histograms
fHistK0PhiEta(0),
fHistPosPhiEta(0),
fHistNegPhiEta(0),    
//correlations
fHistCFPhi(0), 
fHistCFEta(0), 
fHistKChKChPhi(0),
fHistKPosKNegPhi(0),    
fHistCF(0),
fHistCFz(0),
fHistKChKCh(0),
fHistKPosKNeg(0),
fHistKPosKNegz(0),
fHistKPosKPos(0),
fHistKNegKNeg(0), 
//mixing
fSelectedKCh(0),
fSelectedK0s(0),
fSelectedKpos(0),
fSelectedKneg(0),
fPoolMgr(0x0),
fPoolMaxNEvents(1000),
fPoolMinNTracks(5000),
fMinEventsToMix(10),
fNzVtxBins(20),
fNCentBins(15),
fKpKnCorr(kTRUE),
fK0KchCorr(kTRUE),
fKpKpCorr(kTRUE),
fKnKnCorr(kTRUE),
/*
fNOfSamples(1.0),
fSampleIndex(0.0),
*/
hPt(0),
hPt_kPos(0),
fHistCF_Bg(0),
fHistCF_KpKn_Bg(0),
fHistCF_KpKp_Bg(0),
fHistCF_KnKn_Bg(0),
fHistCF_Bgz(0),
fHistCF_KpKn_Bgz(0),
//MC Truth
fMCSelectedK0s(0),
fMCSelectedKCh(0),
fMCSelectedKpos(0),
fMCSelectedKneg(0),
fMCK0Pt(0),
fMCK0PtfullRange(0),
fMCKPlusPt(0),
fMCKPlusPtfullRange(0),
fMCKMinusPt(0),
fMCKMinusPtfullRange(0),
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fMCKch(0),
fHistK0KchMC(0),
fHistKpKnMC(0),
fHistKpKpMC(0),
fHistKnKnMC(0),
fHistK0KchMC_Bg(0),
fHistKpKnMC_Bg(0),
fHistKpKpMC_Bg(0),
fHistKnKnMC_Bg(0),
fHistGenMultiplicity(0),
//eventcuts
fEventCuts(0)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC(const char* name) : AliAnalysisTaskSE(name),
fAOD(0),
fmcEvent(0), 
fOutputList(0), 
fPIDResponse(0),
fPidpTDependentMethod(kTRUE),
fRejectEventPileUp(kTRUE),
fMinBias(kTRUE),
fCentral(kFALSE),
fSemiCentral(kFALSE),
fData(kFALSE),
fMCTruth(kTRUE),
fMCReconstructed(kTRUE),
fMCArray(NULL),
PVx(0), PVy(0), PVz(0),
//track cuts 
fLpTCut(0.4), 
fUpTCut(0.8), 
fEtaCut(0.8), 
fSigCut(3.0),
fBit(96),
fPVzCut(8),
fCentMin(0),
fCentMax(10),
//V0 cuts
fDecayLv0Cut(8.05), 
fLpTv0Cut(0.4), 
fUpTv0Cut(0.8), 
fEtav0Cut(0.8), 
fDcaPosToPrimVtxv0Cut(0.1), 
fDcaNegToPrimVtxv0Cut(0.1), 
fEtaPosv0Cut(0.8), 
fEtaNegv0Cut(0.8), 
fCosPACut(0.99), 
fSigPosv0Cut(3.0), 
fSigNegv0Cut(3.0),
fnumOfTPCcrossedRows(70),
fTPCrowsFindableRatio(0.8),
//PID histograms
fEnergy(0),
fEnergyCuts(0),
fPID(0), 
fPIDKaon(0), 
fPIDK(0), 
fPIDKpiCut(0),
fnsigmakaon(0),
fNsigmaKaon(0),
fNsigmaTOFK(0),
fNsigmaTOFKaon(0),
fNsigmaTPCTOFK(0),
fHistTOFKch(0),
//track and event observables
fVtx(0), 
fClusters(0),
fHistPVz(0),
fHistNEvents(0),
fMCEvents(0),
fTracksCounter(0),
fMCEvents_pileup(0),
fHistNV0(0),  
fHistEta(0), 
fHistDEta(0), 
fHistPhi(0), 
fHistDPhi(0), 
fHistMult(0), 
fHistCent(0),
fHistCent_mcgen(0),
fHistTPCTracksVsClusters(0),
//single particle 1D histograms
fHistMK0(0), 
fHistMK0Cuts(0), 
fHistKChPt(0),
fHistKChPtfullRange(0),
fHistKChPtMix(0),  
fHistK0Pt(0),
fHistK0PtfullRange(0),
fHistKPlusPt(0),
fHistKPlusPtfullRange(0),
fHistKMinusPt(0),
fHistKMinusPtfullRange(0), 
fHistKChPhi(0), 
fHistK0Phii(0), 
fHistKpPhi(0), 
fHistKnPhi(0),
fHistPPionPhi(0), 
fHistNPionPhi(0),
//single particle 2D histograms
fHistK0Phi(0),
fHistK0Eta(0),  
fHistChPhi(0), 
fHistChEta(0), 
fHistChRap(0),
fHistPosPhi(0),
fHistPosEta(0), 
fHistPosRap(0), 
fHistNegPhi(0), 
fHistNegEta(0), 
fHistNegRap(0),
//single particle 3D histograms
fHistK0PhiEta(0),
fHistPosPhiEta(0),
fHistNegPhiEta(0),    
//correlations
fHistCFPhi(0), 
fHistCFEta(0), 
fHistKChKChPhi(0),
fHistKPosKNegPhi(0),    
fHistCF(0),
fHistCFz(0),
fHistKChKCh(0),
fHistKPosKNeg(0),
fHistKPosKNegz(0),
fHistKPosKPos(0),
fHistKNegKNeg(0),  
//mixing
fSelectedKCh(0),
fSelectedK0s(0),
fSelectedKpos(0),
fSelectedKneg(0),
fPoolMgr(0x0),
fPoolMaxNEvents(1000),
fPoolMinNTracks(5000),
fMinEventsToMix(10),
fNzVtxBins(20),
fNCentBins(15),
fKpKnCorr(kTRUE),
fK0KchCorr(kTRUE),
fKpKpCorr(kTRUE),
fKnKnCorr(kTRUE),
/*
fNOfSamples(1.0),
fSampleIndex(0.0),
*/
hPt(0),
hPt_kPos(0),
fHistCF_Bg(0),
fHistCF_KpKn_Bg(0),
fHistCF_KpKp_Bg(0),
fHistCF_KnKn_Bg(0),
fHistCF_Bgz(0),
fHistCF_KpKn_Bgz(0),
//MC Truth
fMCSelectedK0s(0),
fMCSelectedKCh(0),
fMCSelectedKpos(0),
fMCSelectedKneg(0),
fMCK0Pt(0),
fMCK0PtfullRange(0),
fMCKPlusPt(0),
fMCKPlusPtfullRange(0),
fMCKMinusPt(0),
fMCKMinusPtfullRange(0),
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fMCKch(0),
fHistK0KchMC(0),
fHistKpKnMC(0),
fHistKpKpMC(0),
fHistKnKnMC(0),
fHistK0KchMC_Bg(0),
fHistKpKnMC_Bg(0),
fHistKpKpMC_Bg(0),
fHistKnKnMC_Bg(0),
fHistGenMultiplicity(0),
//eventcuts
fEventCuts(0)


{
    // constructor
    DefineInput(0, TChain::Class());    
    DefineOutput(1, TList::Class());    
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PC::~AliAnalysisTaskKaon2PC()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
    if(fSelectedK0s) delete fSelectedK0s;
    if(fSelectedKCh) delete fSelectedKCh;
    if(fSelectedKpos) delete fSelectedKpos;
    if(fSelectedKneg) delete fSelectedKneg;
    if(fMCSelectedK0s) delete fMCSelectedK0s;
    if(fMCSelectedKCh) delete fMCSelectedKCh;
    if(fMCSelectedKpos) delete fMCSelectedKpos;
    if(fMCSelectedKneg) delete fMCSelectedKneg;
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file.

    //fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0}; // 10 bins
    fzVtxBins = {-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}; // 20 bins
    //fCentBins = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100};
    fCentBins = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 }; // 15 bins
    //fsampleBins = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();     //to get pid response object
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler) fPIDResponse = inputHandler->GetPIDResponse();
    }

    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE); 

    //========== event mixing -> poolmanager initialization ======== 

    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);      

    //PID histograms
    fEnergy = new TH1F("fEnergy", "particle yield vs energy loss",800,0,200);
    fEnergyCuts = new TH1F("fEnergyCuts", "Energyloss distribution after pion, electron rejection",500,30,110);

    fPID= new TH2F("fPID","Particle Identification",500,0.2,1.5,500,0.0,1000.0);
    fPID->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPID->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPID->SetOption("colz");

    fPIDKaon= new TH2F("fPIDKaon","Kaons with |nSigma|<3",500,0.2,1.5,500,0.0,1000.0);
    fPIDKaon->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDKaon->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKaon->SetOption("colz");

    fPIDK = new TH2F("fPIDK","Kaons with |nSigma|<2",500,0.2,1.5,500,0.0,1000.0);
    fPIDK->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDK->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDK->SetOption("colz");

    fPIDKpiCut= new TH2F("fPIDKpiCut","Kaon PID with Pion Rejection",500,0.2,1.5,500,0.0,1000.0);
    fPIDKpiCut->GetXaxis()->SetTitle("p_{} [GeV/c]");
    fPIDKpiCut->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKpiCut->SetOption("colz");

    fnsigmakaon = new TH2F("fnsigmakaon","n#sigma of Kaon Vs Transverse Momentum",500,0.2,2,700,-20.0,20.0);
    fnsigmakaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fnsigmakaon->GetYaxis()->SetTitle("n#sigma of Kaon");
    fnsigmakaon->SetOption("colz");

    fNsigmaKaon = new TH2F("fNsigmaKaon","|n#sigma|<2 Vs Pt of Kaons",700,0.2,2,700,-20.0,20.0);
    fNsigmaKaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaKaon->GetYaxis()->SetTitle("n#sigma of Kaon");
    fNsigmaKaon->SetOption("colz");

    fNsigmaTOFK = new TH2F("fNsigmaTOFK","",700,0.2,2,700,-60.0,60.0);
    fNsigmaTOFK->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaTOFK->GetYaxis()->SetTitle("n#sigma_{TOF} of Kaon");
    fNsigmaTOFK->SetOption("colz");

    fNsigmaTOFKaon = new TH2F("fNsigmaTOFKaon","",700,0.2,2,700,-10.0,10.0);
    fNsigmaTOFKaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF} of Kaon");
    fNsigmaTOFKaon->SetOption("colz");

    fNsigmaTPCTOFK = new TH2F("fNsigmaTPCTOFK","",700,0.2,2,700,-20.0,20.0);
    fNsigmaTPCTOFK->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaTPCTOFK->GetYaxis()->SetTitle("n#sigma_{TOF} of Kaon");
    fNsigmaTPCTOFK->SetOption("colz");

    fHistTOFKch = new TH2F("fHistTOFKch", "", 500, 0.2,2,700,0,2);

    //track and event observables
    fVtx = new TH1F("fVtx", "PV_{z} distribution of Tracks", 100, -13, 13);
    fClusters = new TH1F("fClusters", "TPCClusters distribution", 500, 1, 170);
    fHistPVz = new TH1F("fHistPVz", "PVz Distribution", 20, -10, 10);
    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 6, 0, 6);
    fMCEvents = new TH1F("fMCEvents", "fMCEvents", 6, 0, 6);
    fTracksCounter = new TH1F("fTracksCounter", "fTracksCounter", 4, 0, 4);
    fMCEvents_pileup = new TH1F("fMCEvents_pileup", "fMCEvents_pileup", 1, 0, 1);
    fHistNV0 = new TH1F("fHistNV0","Number of V0s",100, 0, 5000);
    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -5, 5);
    fHistDEta = new TH1F("fHistDEta", "fHistDEta", 100, -10, 10);
    fHistPhi = new TH1F("fHistPhi", "Phi Distribution", 100, 0, 7);
    fHistDPhi = new TH1F("fHistDPhi", "fHistDPhi", 100, 0, 10);
    fHistMult = new TH1F("fHistMult", "Number of tracks", 100, 0, 100);
    fHistCent = new TH1F("fHistCent", "CentV0M", 100, 0, 100);
    fHistCent_mcgen = new TH1F("fHistCent_mcgen", "fHistCent_mcgen", 100, 0, 100);
    
    fHistTPCTracksVsClusters =  new TH2F("fHistTPCTracksVsClusters","fHistTPCTracksVsClusters",400,0,400,50000,0,50000);
    fHistTPCTracksVsClusters->Sumw2();
    fHistTPCTracksVsClusters->GetXaxis()->SetTitle("TPC tracks");
    fHistTPCTracksVsClusters->GetYaxis()->SetTitle("TPC Clusters");

    //single particle histograms
    fHistMK0=new TH1F("fHistMK0", "Invariant Mass Distribution of Neutral Kaons", 100, 0.4, 0.6);
    fHistMK0Cuts=new TH1F("fHistMK0Cuts", "Invariant Masss Distribution of Neutral Kaons After cuts", 100, 0.4, 0.6);

    fHistKChPt = new TH1F("fHistKChPt", "p_{T} distribution of all Charged Kaon Tracks", 100, fLpTCut, fUpTCut);
    fHistKChPt->SetOption("HIST E p");
    fHistKChPtfullRange = new TH1F("fHistKChPtfullRange", "p_{T} distribution of all Charged Kaon Tracks", 100, 0, 2);
    fHistKChPtfullRange->SetOption("HIST E p");

    fHistKPlusPt = new TH1F("fHistKPlusPt", "", 100, fLpTCut, fUpTCut);
    fHistKPlusPt->SetOption("HIST E p");
    fHistKPlusPtfullRange = new TH1F("fHistKPlusPtfullRange", "", 100, 0, 2);
    fHistKPlusPtfullRange->SetOption("HIST E p");

    fHistKMinusPt = new TH1F("fHistKMinusPt", "", 100, fLpTCut, fUpTCut);
    fHistKMinusPt->SetOption("HIST E p");
    fHistKMinusPtfullRange = new TH1F("fHistKMinusPtfullRange", "", 100, 0, 2);
    fHistKMinusPtfullRange->SetOption("HIST E p");

    fHistK0Pt = new TH1F("fHistK0Pt", "p_{T} distribution of all Neutral Kaon Tracks", 100, fLpTv0Cut, fUpTv0Cut);
    fHistK0Pt->SetOption("HIST E p");
    fHistK0PtfullRange = new TH1F("fHistK0PtfullRange", "", 100, 0, 2);
    fHistK0PtfullRange->SetOption("HIST E p");

    fHistKChPtMix= new TH1F("fHistKChPtMix", "p_{T} distribution of all Charged Kaon Tracks", 100, 0, 2);
    fHistKChPtMix->SetOption("HIST E p");

    fHistKChPhi = new TH1F("fHistKChPhi", "#phi distribution of all Charged Kaon Tracks", 100, 0, 2*Pi);
    fHistKChPhi->SetOption("HIST E p");
    fHistK0Phi = new TH1F("fHistK0Phi", "#phi distribution of all Neutral Kaon Tracks", 100, 0, 2*Pi);
    fHistK0Phi->SetOption("HIST E p");

    fHistKpPhi = new TH1F("fHistKpPhi", "#phi distribution of all K^{+} Tracks", 100, 0, 2*Pi);
    fHistKpPhi->SetOption("HIST E p");
    fHistKnPhi = new TH1F("fHistKnPhi", "#phi distribution of all K^{-} Tracks", 100, 0, 2*Pi);
    fHistKnPhi->SetOption("HIST E p");

    fHistPPionPhi = new TH1F("fHistPPionPhi","Number of +ve piondaughters Vs Track Phi", 100,0,2*Pi);
    fHistNPionPhi = new TH1F("fHistNPionPhi","Number of -ve piondaughters Vs Track Phi", 100,0,2*Pi);
    
    fHistK0Phii = new TH1F("fHistK0Phii", "Number of K0s Vs Track Phi; Centrality", 16,0,2*Pi);
    fHistK0Phii->GetXaxis()->SetTitle("V0 Phi");
    
    fHistK0Eta = new TH1F("fHistK0Eta", "Number of K0s Vs Track Eta", 16,-0.8, 0.8);
    fHistK0Eta->GetXaxis()->SetTitle("V0 Eta");

    fHistChPhi = new TH1F("fHistChPhi", "Number of charged particles Vs Track Phi", 16,0,2*Pi);
    fHistChPhi->GetXaxis()->SetTitle("Track Phi (in radians)");

    fHistChEta = new TH1F("fHistChEta", "Number of charged particles Vs Track Eta", 16,-0.8, 0.8);
    fHistChEta->GetXaxis()->SetTitle("Track Eta");

    fHistChRap = new TH1F("fHistChRap", "Number of charged particles Vs Rapidity", 16,-0.8, 0.8);
    fHistChRap->GetXaxis()->SetTitle("Rapidity");

    fHistPosPhi = new TH1F("fHistPosPhi","Number of +ve Kaons Vs Track Phi", 16,0,2*Pi);
    fHistPosPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    
    fHistPosEta = new TH1F("fHistPosEta", "", 16,-0.8, 0.8);
    fHistPosEta->GetXaxis()->SetTitle("Track Eta");
    
    fHistPosRap = new TH1F("fHistPosRap", "", 16,-0.8, 0.8);
    fHistPosRap->GetXaxis()->SetTitle("Rapidity");
    
    fHistNegPhi = new TH1F("fHistNegPhi", "Number of -ve Kaons Vs Track Phi", 16,0,2*Pi);
    fHistNegPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    
    fHistNegEta = new TH1F("fHistNegEta", "", 16,-0.8, 0.8);
    fHistNegEta->GetXaxis()->SetTitle("Track Eta");

    fHistNegRap = new TH1F("fHistNegRap", "", 16,-0.8, 0.8);
    fHistNegRap->GetXaxis()->SetTitle("Rapidity");

    fHistPosPhiEta = new TH2F("fHistPosPhiEta", "",60,0,2*Pi, 60,-0.8, 0.8);
    fHistPosPhiEta->GetXaxis()->SetTitle("Track Eta");
    fHistPosPhiEta->SetOption("SURF1");

    fHistNegPhiEta = new TH2F("fHistNegPhiEta", "",60,0,2*Pi, 60,-0.8, 0.8);
    fHistNegPhiEta->GetXaxis()->SetTitle("Track Eta");
    fHistNegPhiEta->SetOption("SURF1");

    fHistK0PhiEta = new TH2F("fHistK0PhiEta", "", 60,0,2*Pi,60,-0.8, 0.8);
    fHistK0PhiEta->GetXaxis()->SetTitle("V0 Phi (in radians)");
    fHistK0PhiEta->SetOption("SURF1");

    //correlations
    fHistCFPhi = new TH2F("fHistCFPhi","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi and #Delta#eta",72,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCFPhi->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhi->GetYaxis()->SetTitle("#Delta#eta");
    fHistCFPhi->SetOption("SURF1");
    
    fHistCFEta = new TH2F("fHistCFEta","Number of pairs of K^{+/-} and K^{0} Vs #Delta#eta",32,-1.6, 1.6, 100,0,100);
    fHistCFEta->GetXaxis()->SetTitle("#Delta#eta ");
    fHistCFEta->GetYaxis()->SetTitle("Centrality in %");
    fHistCFEta->SetOption("SURF1");

    fHistKPosKNegPhi = new TH2F("fHistKPosKNegPhi", "C(#Delta#phi) of K^{+}-K{-} ",32,-0.5*Pi,1.5*Pi,100,0,100);
    fHistKPosKNegPhi->SetOption("SURF1");
    fHistKChKChPhi = new TH2F("fHistKChKChPhi", "C(#Delta#phi) of Kch-Kch ",32,-0.5*Pi,1.5*Pi,100,0,100);
    fHistKChKChPhi->SetOption("SURF1"); 

    fHistCF = new TH2F("fHistCF","Number of pairs of K^{+/-} and K^{0}",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistCF->GetXaxis()->SetTitle("#Delta#phi ");
    fHistCF->GetYaxis()->SetTitle("#Delta#eta");
    fHistCF->SetOption("SURF1");

    fHistCFz = new TH3F("fHistCFz","Number of pairs of K^{+/-} and K^{0}: PVz",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6, 20,-10,10);
    fHistCFz->GetXaxis()->SetTitle("#Delta#phi ");
    fHistCFz->GetYaxis()->SetTitle("#Delta#eta");

    fHistKPosKNeg = new TH2F("fHistKPosKNeg","K^{+}-K^{-} Correlation",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKPosKNeg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKNeg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKPosKNeg->SetOption("SURF1");

    fHistKPosKNegz = new TH3F("fHistKPosKNegz","K^{+}-K^{-} Correlation: PVz",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6,20,-10,10);
    fHistKPosKNegz->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKNegz->GetYaxis()->SetTitle("#Delta#eta");

    fHistKChKCh = new TH2F("fHistKChKCh","Kch-Kch Correlation",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKChKCh->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKChKCh->GetYaxis()->SetTitle("#Delta#eta");
    fHistKChKCh->SetOption("SURF1");

    fHistKPosKPos = new TH2F("fHistKPosKPos","K^{+}-K^{+} Correlation",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),32,-1.6, 1.6 );
    fHistKPosKPos->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKPos->GetYaxis()->SetTitle("#Delta#eta");
    fHistKPosKPos->SetOption("SURF1");

    fHistKNegKNeg = new TH2F("fHistKNegKNeg","K^{-}-K^{-} Correlation",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),32,-1.6, 1.6 );
    fHistKNegKNeg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKNegKNeg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKNegKNeg->SetOption("SURF1");

    //mixing histograms
    //const Int_t sizeOfSamples = (Int_t) fNOfSamples;

    fHistCF_Bg = new TH2F("fHistCF_Bg","Background for CF of K^{+/-} and K^{0}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_Bg->SetOption("SURF1");

    fHistCF_KpKn_Bg = new TH2F("fHistCF_KpKn_Bg","Background for CF of K^{+} and K^{-}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_KpKn_Bg->SetOption("SURF1");

    fHistCF_KpKp_Bg = new TH2F("fHistCF_KpKp_Bg","Background for CF of K^{+} and K^{+}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_KpKp_Bg->SetOption("SURF1");

    fHistCF_KnKn_Bg = new TH2F("fHistCF_KnKn_Bg","Background for CF of K^{-} and K^{-}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_KnKn_Bg->SetOption("SURF1");

    fHistCF_Bgz = new TH3F("fHistCF_Bgz","Background for CF of K^{+/-} and K^{0}: pVz",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6,20,-10,10);
    fHistCF_Bgz->SetOption("SURF1");

    fHistCF_KpKn_Bgz = new TH3F("fHistCF_KpKn_Bgz","Background for CF of K^{+} and K^{-}: PVz",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6,20,-10,10);
    fHistCF_KpKn_Bgz->SetOption("SURF1");

    hPt = new TH1F("hPt", "Track pT distribution", 100, 0, 2);
    hPt_kPos = new TH1F("hPt_kPos", "Track pT distribution", 100, 0, 2);

    // MC Truth histograms

    //+++++++++++++++++++++ MC ++++++++++++++++++++++++++

    Int_t bins[4] = {100,100,100,100};
    Double_t min[4] = {fLpTCut,0,-1.0,0.4};
    Double_t max[4] = {fUpTCut,2*Pi,1.0,0.6};

    fMCK0Pt = new TH1F("fMCK0Pt", "", 100, fLpTv0Cut, fUpTv0Cut);
    fMCK0Pt->SetOption("HIST E p");
    fMCK0PtfullRange = new TH1F("fMCK0PtfullRange", "", 100, 0, 2);
    fMCK0PtfullRange->SetOption("HIST E p");

    fMCKPlusPt = new TH1F("fMCKPlusPt", "", 100, fLpTCut, fUpTCut);
    fMCKPlusPt->SetOption("HIST E p");
    fMCKPlusPtfullRange = new TH1F("fMCKPlusPtfullRange", "", 100, 0, 2);
    fMCKPlusPtfullRange->SetOption("HIST E p");

    fMCKMinusPt = new TH1F("fMCKMinusPt", "", 100, fLpTCut, fUpTCut);
    fMCKMinusPt->SetOption("HIST E p");
    fMCKMinusPtfullRange = new TH1F("fMCKMinusPtfullRange", "", 100, 0, 2);
    fMCKMinusPtfullRange->SetOption("HIST E p");
    
    fMCK0 = new THnSparseF("fMCK0","fMCK0",4,bins,min,max);
    fMCK0->GetAxis(0)->SetTitle("p_{T} of K^{0}_{S}");
    fMCK0->GetAxis(1)->SetTitle("#phi");
    fMCK0->GetAxis(2)->SetTitle("#eta");
    fMCK0->GetAxis(3)->SetTitle("mass");

    fMCKpos = new THnSparseF("fMCKpos","fMCKpos",4,bins,min,max);
    fMCKpos->GetAxis(0)->SetTitle("p_{T} of K^{+}");
    fMCKpos->GetAxis(1)->SetTitle("#phi of K^{+}");
    fMCKpos->GetAxis(2)->SetTitle("#eta of K^{+}");
    fMCKpos->GetAxis(3)->SetTitle("mass of K^{+}");

    fMCKneg = new THnSparseF("fMCKneg","fMCKneg",4,bins,min,max);
    fMCKneg->GetAxis(0)->SetTitle("p_{T} of K^{-}");
    fMCKneg->GetAxis(1)->SetTitle("#phi of K^{-}");
    fMCKneg->GetAxis(2)->SetTitle("#eta of K^{-}");
    fMCKneg->GetAxis(3)->SetTitle("mass of K^{-}");

    fMCKch = new THnSparseF("fMCKch","fMCKch",4,bins,min,max);
    fMCKch->GetAxis(0)->SetTitle("p_{T} of K^{ch}");
    fMCKch->GetAxis(1)->SetTitle("#phi of K^{ch}");
    fMCKch->GetAxis(2)->SetTitle("#eta of K^{ch}");
    fMCKch->GetAxis(3)->SetTitle("mass of K^{ch}");

    fHistK0KchMC = new TH2F("fHistK0KchMC","K^{0}-K^{ch} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistK0KchMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistK0KchMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistK0KchMC->SetOption("SURF1");

    fHistKpKnMC = new TH2F("fHistKpKnMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKnMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKnMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKnMC->SetOption("SURF1");

    fHistKpKpMC = new TH2F("fHistKpKpMC","K^{+}-K^{+} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKpMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKpMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKpMC->SetOption("SURF1");

    fHistKnKnMC = new TH2F("fHistKnKnMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKnKnMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKnKnMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistKnKnMC->SetOption("SURF1");

    fHistK0KchMC_Bg = new TH2F("fHistK0KchMC_Bg","K^{0}-K^{ch} Background Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistK0KchMC_Bg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistK0KchMC_Bg->GetYaxis()->SetTitle("#Delta#eta");
    fHistK0KchMC_Bg->SetOption("SURF1");

    fHistKpKnMC_Bg = new TH2F("fHistKpKnMC_Bg","K^{+}-K^{-} Correlation Background for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKnMC_Bg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKnMC_Bg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKnMC_Bg->SetOption("SURF1");

    fHistKpKpMC_Bg = new TH2F("fHistKpKpMC_Bg","K^{+}-K^{+} Correlation Background for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKpMC_Bg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKpMC_Bg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKpMC_Bg->SetOption("SURF1");

    fHistKnKnMC_Bg = new TH2F("fHistKnKnMC_Bg","K^{-}-K^{-} Correlation Background for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKnKnMC_Bg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKnKnMC_Bg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKnKnMC_Bg->SetOption("SURF1");

    fHistGenMultiplicity = new TH1D ("fHistGenMultiplicity","fHistGenMultiplicity",500,0,500);


    fOutputList->Add(fEnergy);
    fOutputList->Add(fEnergyCuts);
    fOutputList->Add(fPID);
    fOutputList->Add(fPIDKaon);
    fOutputList->Add(fPIDK);
    fOutputList->Add(fPIDKpiCut);
    fOutputList->Add(fnsigmakaon);
    fOutputList->Add(fNsigmaKaon);
    fOutputList->Add(fNsigmaTOFK);
    fOutputList->Add(fNsigmaTOFKaon);
    fOutputList->Add(fNsigmaTPCTOFK);
    fOutputList->Add(fHistTOFKch);

    fOutputList->Add(fVtx);
    fOutputList->Add(fClusters);
    fOutputList->Add(fHistPVz);
    fOutputList->Add(fHistNEvents);
    fOutputList->Add(fMCEvents);
    fOutputList->Add(fTracksCounter);
    fOutputList->Add(fMCEvents_pileup);
    fOutputList->Add(fHistNV0);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistDEta);
    fOutputList->Add(fHistPhi);
    fOutputList->Add(fHistDPhi);
    fOutputList->Add(fHistMult);
    fOutputList->Add(fHistCent);
    fOutputList->Add(fHistCent_mcgen);
    fOutputList->Add(fHistTPCTracksVsClusters);

    fOutputList->Add(fHistMK0);
    fOutputList->Add(fHistMK0Cuts);
    fOutputList->Add(fHistKChPt);
    fOutputList->Add(fHistKChPtfullRange);
    fOutputList->Add(fHistKChPtMix);
    fOutputList->Add(fHistK0Pt);
    fOutputList->Add(fHistK0PtfullRange);
    fOutputList->Add(fHistKPlusPt);
    fOutputList->Add(fHistKPlusPtfullRange);
    fOutputList->Add(fHistKMinusPt);
    fOutputList->Add(fHistKMinusPtfullRange);
    fOutputList->Add(fHistKChPhi);
    fOutputList->Add(fHistK0Phii);
    fOutputList->Add(fHistKpPhi);
    fOutputList->Add(fHistKnPhi);  
    fOutputList->Add(fHistPPionPhi);
    fOutputList->Add(fHistNPionPhi);

    fOutputList->Add(fHistK0Phi);
    fOutputList->Add(fHistK0Eta);
    fOutputList->Add(fHistChPhi);
    fOutputList->Add(fHistChEta);
    fOutputList->Add(fHistChRap);
    fOutputList->Add(fHistPosPhi);
    fOutputList->Add(fHistPosEta);
    fOutputList->Add(fHistPosRap);
    fOutputList->Add(fHistNegPhi);
    fOutputList->Add(fHistNegEta);
    fOutputList->Add(fHistNegRap);

    fOutputList->Add(fHistK0PhiEta);
    fOutputList->Add(fHistPosPhiEta);
    fOutputList->Add(fHistNegPhiEta);
    
    fOutputList->Add(fHistCFPhi);
    fOutputList->Add(fHistCFEta);
    fOutputList->Add(fHistKChKChPhi);
    fOutputList->Add(fHistKPosKNegPhi);
    fOutputList->Add(fHistCF);
    fOutputList->Add(fHistCFz);
    fOutputList->Add(fHistKChKCh);
    fOutputList->Add(fHistKPosKNeg);
    fOutputList->Add(fHistKPosKPos);
    fOutputList->Add(fHistKNegKNeg);
    fOutputList->Add(fHistKPosKNegz);
    fOutputList->Add(hPt);
    fOutputList->Add(hPt_kPos);

    fOutputList->Add(fHistCF_Bg);
    fOutputList->Add(fHistCF_KpKn_Bg);
    fOutputList->Add(fHistCF_KpKp_Bg);
    fOutputList->Add(fHistCF_KnKn_Bg);
    fOutputList->Add(fHistCF_Bgz);
    fOutputList->Add(fHistCF_KpKn_Bgz);

    fOutputList->Add(fMCK0);
    fOutputList->Add(fMCKpos);
    fOutputList->Add(fMCKneg);
    fOutputList->Add(fMCKch);
    fOutputList->Add(fMCK0Pt);
    fOutputList->Add(fMCK0PtfullRange);
    fOutputList->Add(fMCKPlusPt);
    fOutputList->Add(fMCKPlusPtfullRange);
    fOutputList->Add(fMCKMinusPt);
    fOutputList->Add(fMCKMinusPtfullRange);
    fOutputList->Add(fHistK0KchMC);
    fOutputList->Add(fHistKpKnMC);
    fOutputList->Add(fHistKpKpMC);
    fOutputList->Add(fHistKnKnMC);
    fOutputList->Add(fHistK0KchMC_Bg);
    fOutputList->Add(fHistKpKnMC_Bg);
    fOutputList->Add(fHistKpKpMC_Bg);
    fOutputList->Add(fHistKnKnMC_Bg);
    fOutputList->Add(fHistGenMultiplicity);
    
    fEventCuts.AddQAplotsToList(fOutputList);
    PostData(1, fOutputList);
    
}

void AliAnalysisTaskKaon2PC::SetTrackCuts(Double_t c1, Double_t c2, Double_t c3, Double_t c4) {
    
    fLpTCut = c1;
    fUpTCut = c2;
    fEtaCut = c3;
    fSigCut = c4;
}

void AliAnalysisTaskKaon2PC::SetV0TrackCuts(Double_t c5, Double_t c6, Double_t c7, Double_t c8, Double_t c9, Double_t c10, Double_t c11, Double_t c12, Double_t c13, Double_t c14, Double_t c15) {
    
    fDecayLv0Cut = c5;
    fLpTv0Cut = c6;
    fUpTv0Cut = c7;
    fEtav0Cut = c8;
    fDcaPosToPrimVtxv0Cut = c9;
    fDcaNegToPrimVtxv0Cut = c10;
    fEtaPosv0Cut = c11;
    fEtaNegv0Cut = c12;
    fCosPACut = c13;
    fSigPosv0Cut = c14;
    fSigNegv0Cut = c15;
}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskKaon2PC::AcceptTrack(const AliAODTrack *Trk) {

    if (!Trk->TestFilterBit(96)) return kFALSE;
    if (Trk->Charge() == 0) return kFALSE;         //excluding neutral particles
    if (fabs(Trk->Eta()) > fEtaCut) return kFALSE; // eta cut

    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kKaon);
    Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kPion);
    Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kElectron);
    Double_t nSigmaproton =TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kProton)) ;
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kKaon);
    Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kElectron);
    Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kPion);

    //if (fabs(nSigmaelectron) < 2.0) return kFALSE;  // excluding electrons via TPC
    //if (fabs(nSigmapion) < 2.0) return kFALSE;      // excluding pions via TPC
    //if (fabs(nSigmaproton) < 2.0) return kFALSE;      // excluding pions via TPC
    //if (fabs(nSigmakaon < 3.0) && (nSigmapion < 3.0)) return kFALSE;
    //if (fabs(nSigmakaon < 3.0) && (nSigmaproton < 3.0)) return kFALSE;
    //if (fabs(nSigmakaon < 3.0) && (nSigmaelectron < 3.0)) return kFALSE;

    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskKaon2PC::IsKaonNSigma3(float mom, float nSigmakaon, float nSigmaTOFkaon)
{
  if (mom > 0.5) {

    if (TMath::Hypot( nSigmakaon, nSigmaTOFkaon ) < fSigCut)
      return kTRUE;
  }
  else {
    if (TMath::Abs(nSigmakaon) < fSigCut)
      return kTRUE;
  }

  return kFALSE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskKaon2PC::AcceptPosTrack(const AliAODTrack *Trk) {
    if(Trk->Charge() < 0) return kFALSE;        // excluding negative tracks

    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskKaon2PC::AcceptNegTrack(const AliAODTrack *Trk) {
    if(Trk->Charge() > 0) return kFALSE;        // excluding positive tracks

    return kTRUE;
}

//_____________________________________________________________________________

Double_t AliAnalysisTaskKaon2PC::Beta(const AliAODTrack *track)
{
    Double_t startTime     = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());         //in ps
    Double_t stoptime      = track->GetTOFsignal();
    Double_t c             = TMath::C()*1.E-9;                                                              // m/ns
    Double_t length        = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kKaon)*1E-3*c;
    stoptime -= startTime;      
    Double_t scaleStopTime = stoptime*1E-3;          
    scaleStopTime          = scaleStopTime*c;
    return length/scaleStopTime;
}

//___________________________  v0 selection with no pT cut _________________________

Bool_t AliAnalysisTaskKaon2PC::AcceptV0(const AliAODv0 *v0, Double_t *vertex) {

    if (v0->GetOnFlyStatus()) return kFALSE;
    Double_t length = v0->DecayLengthV0(vertex);
    if (length > fDecayLv0Cut) return kFALSE;
    Double_t dcaDau = v0->DcaV0Daughters();
    if (TMath::Abs(dcaDau)  > 0.8) return kFALSE;
    Double_t pT = v0->Pt();
    //Double_t pT=TMath::Sqrt(v0->Pt2V0());
    if (fabs(v0->Eta()) > fEtav0Cut) return kFALSE;
    Double_t DCAtoPV = v0->DcaV0ToPrimVertex();
    if (DCAtoPV > 0.1 ) return kFALSE;
    Double_t dcaPosToPV = v0->DcaPosToPrimVertex();
    if (TMath::Abs(dcaPosToPV) < fDcaPosToPrimVtxv0Cut && fDcaPosToPrimVtxv0Cut != -999) return kFALSE;
    Double_t dcaNegToPV = v0->DcaNegToPrimVertex();
    if (TMath::Abs(dcaNegToPV) < fDcaNegToPrimVtxv0Cut && fDcaNegToPrimVtxv0Cut != -999) return kFALSE;
    Double_t etaPos = v0->PseudoRapPos();
    if (fabs(etaPos) > fEtaPosv0Cut && fEtaPosv0Cut != -999) return kFALSE;
    Double_t etaNeg = v0->PseudoRapNeg();
    if (fabs(etaNeg) > fEtaNegv0Cut && fEtaNegv0Cut != -999) return kFALSE;
    Double_t cosPA= v0->CosPointingAngle(vertex);
    if (cosPA < fCosPACut && fCosPACut != -999) return kFALSE;
    Double_t armpt = v0->PtArmV0();
    Double_t alpha = v0->AlphaV0();
    if (armpt < 0.2*fabs(alpha)) return kFALSE;
    AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0);
    AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1);
    Double_t nSigmaPionPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
    if (fabs(nSigmaPionPos) > fSigPosv0Cut) return kFALSE;
    Double_t nSigmaPionNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
    if (fabs(nSigmaPionNeg) > fSigNegv0Cut) return kFALSE;
    Double_t pos_clusters = pTrack->GetTPCNcls();
    Double_t neg_clusters = nTrack->GetTPCNcls();
    if (pos_clusters < 70) return kFALSE;
    if (neg_clusters < 70) return kFALSE;
    Double_t massLam = v0->MassLambda();
    Double_t massALam = v0->MassAntiLambda();
    if (massLam > 1.11 && massLam < 1.12) return kFALSE;
    if (massALam > 1.11 && massALam < 1.12) return kFALSE;

    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskKaon2PC::IsMyGoodDaughterTrack(const AliAODTrack *t) {
	// TPC refit
 	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < fnumOfTPCcrossedRows) return kFALSE;
	Int_t findable=t->GetTPCNclsF();

    if (findable <= 0) return kFALSE;
    if (nCrossedRowsTPC/findable < fTPCrowsFindableRatio) return kFALSE;


    if (TMath::Abs(t->Eta())>=fEtaCut) return kFALSE;

	return kTRUE;

}
//_____________________________________________________________________________

// MC Truth Booleans

Bool_t AliAnalysisTaskKaon2PC::SelectK0TracksMC(AliMCParticle *mcTrack ) {
    if (mcTrack->Pt() <= fLpTCut || mcTrack->Pt() >= fUpTCut) return kFALSE; // mc pt cut
    if (fabs(mcTrack->Eta()) > fEtaCut) return kFALSE; // mc eta cut
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectK0 =  mcPartPdg==310&& (isPhysPrim); 
    if (SelectK0) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKPosTracksMC(AliMCParticle *mcTrack ) {
    if (mcTrack->Pt() <= fLpTCut || mcTrack->Pt() >= fUpTCut) return kFALSE; // mc pt cut
    if (fabs(mcTrack->Eta()) > fEtaCut) return kFALSE; // mc eta cut
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKpos =  mcPartPdg==321&& (isPhysPrim); 
    if (SelectKpos) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKNegTracksMC(AliMCParticle *mcTrack ) {
    if (mcTrack->Pt() <= fLpTCut || mcTrack->Pt() >= fUpTCut) return kFALSE; // mc pt cut
    if (fabs(mcTrack->Eta()) > fEtaCut) return kFALSE; // mc eta cut
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKneg =  mcPartPdg==-321&& (isPhysPrim); 
    if (SelectKneg) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKchTracksMC(AliMCParticle *mcTrack ) {
    if (mcTrack->Pt() <= fLpTCut || mcTrack->Pt() >= fUpTCut) return kFALSE; // mc pt cut
    if (fabs(mcTrack->Eta()) > fEtaCut) return kFALSE; // mc eta cut
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKch =  (mcPartPdg==321 || mcPartPdg==-321) &&(isPhysPrim); 
    if (SelectKch) return kFALSE;
    return kTRUE;
}


//_____________________________________________________________________________

void AliAnalysisTaskKaon2PC::RunData() {

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;

    fHistNEvents->Fill(0.5);

    //  data trigger selection
	Bool_t isSelected = kFALSE;

    if(fMinBias) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    if(fCentral) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
    if(fSemiCentral)isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
	if (!isSelected) return;    

    fHistNEvents->Fill(1.5);                               

    //==================== Pile Up==========================

    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (fRejectEventPileUp){
        if (!fEventCuts.AcceptEvent(fAOD)) return;
    }

    fHistNEvents->Fill(2.5);

    //primary Vertex    
    const AliVVertex* primVertex = fEventCuts.GetPrimaryVertex(); 
    if (!primVertex) return;
    
    Double_t PVx = primVertex->GetX();
    Double_t PVy = primVertex->GetY();
    Double_t PVz = primVertex->GetZ();

    Double_t vertex[3] = { -100.0, -100.0, -100.0 }; //?
    primVertex->GetXYZ(vertex);  

    //zvertex cut
    if ( ( TMath::Abs(PVz) ) >= fPVzCut) return ;
    fVtx->Fill(PVz);

    fHistNEvents->Fill(3.5);

    //Multiplicity selection
    //AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    //if(!MultSelection) return;
    
    //centrality
    //double CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    Double_t CentV0M = fEventCuts.GetCentrality(); //centrality
    cout << "centrality of event is " << CentV0M << endl;
    if ((CentV0M < fCentMin)||(CentV0M > fCentMax)) return;

    fHistNEvents->Fill(4.5);

//======== PID loop (no pT cut) ===========

    Int_t iTracks = (fAOD->GetNumberOfTracks());
    cout << "Number of tracks in an event is " << iTracks << endl;
    Double_t MultTPC = fAOD->GetNumberOfTPCTracks();

    Int_t nv0s = (fAOD->GetNumberOfV0s());
    Double_t tpcClusters = fAOD->GetNumberOfTPCClusters();
    
    fHistTPCTracksVsClusters->Fill(MultTPC,tpcClusters);
    //fHistV0AmplitudeVsPVposition->Fill(PVz,VZEROmultiplicity);           
    
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));       
        Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        Double_t nSigmaproton =TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) ;
        Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        Double_t nClustersTPC = track->GetTPCNcls();
        Double_t kchPt = track->Pt();
        Double_t charge = track->Charge();
        if(!track) continue;                               // if we failed, skip this track
        if(abs(track->Eta())>0.8) continue;                // eta cut
        if (!track->TestFilterBit(96)) continue;           // filterbit selection
        if (track->Charge() == 0) continue;  
        fClusters->Fill(nClustersTPC);
        fEnergy->Fill(track->GetTPCsignal());
        fPID->Fill(track->Pt(),track->GetTPCsignal());
        if (fabs(nSigmaelectron) < 2.0) continue;          // excluding electrons via TPC
        if (fabs(nSigmapion) < 2.0) continue;              // excluding pions via TPC
        if (fabs(nSigmaproton) < 2.0) continue;              // excluding protons via TPC
        if (fabs(nSigmakaon)<2.0) {
            fPIDK->Fill(track->Pt(),track->GetTPCsignal());
            fNsigmaKaon->Fill(track->Pt(), nSigmakaon);
            fEnergyCuts->Fill(track->GetTPCsignal()); }
        if ((fabs(nSigmakaon)<2.0) && (0.0 < kchPt< 2.0)) {hPt->Fill(kchPt);}
        if ((fabs(nSigmakaon)<2.0) && (0.0 < kchPt< 2.0) && (charge > 0)) {hPt_kPos->Fill(kchPt);}
        if (fabs(nSigmakaon) < 20.0 ) {fnsigmakaon->Fill(track->Pt(), nSigmakaon);}
        fNsigmaTOFK->Fill(track->Pt(), nSigmaTOFkaon);
        if (fabs(nSigmaTOFkaon)<3.0) {fNsigmaTOFKaon->Fill(track->Pt(), nSigmaTOFkaon);}
        if ((fabs(nSigmakaon)<3.0) && (fabs(nSigmaelectron) > 3.0 ) && (fabs(nSigmapion) > 3.0 ) &&
            (fabs(nSigmaTOFelectron) > 3.0 ) && (fabs(nSigmaTOFpion) > 3.0 ) && (fabs(nSigmaTOFkaon)<3.0)){
            fPIDKpiCut->Fill(track->Pt(),track->GetTPCsignal()); }
         
    }   //end of PID loop

//======== v0 loop (no invariant mass cut) =========

    fHistNV0->Fill(nv0s); 
    for(Int_t i = 0; i < nv0s; i++)  {
        AliAODv0 *v0=fAOD->GetV0(i);                    
        if(!v0) {
        cout<<"No V0 "<<endl;
        continue;
        }
        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        fHistMK0->Fill(v0->MassK0Short());
        if(!AcceptV0(v0, vertex)) continue;
        Double_t pT = v0->Pt();
        if (pT <= fLpTv0Cut || pT >= fUpTv0Cut) continue;
        fHistMK0Cuts->Fill(v0->MassK0Short());
    }

//======== Charged Kaon Track Selection  ==========
fSelectedKCh = new TObjArray;
fSelectedKCh->SetOwner(kTRUE);

fSelectedKpos = new TObjArray;
fSelectedKpos->SetOwner(kTRUE);

fSelectedKneg = new TObjArray;
fSelectedKneg->SetOwner(kTRUE);

if(!fPidpTDependentMethod){


for(Int_t i=0; i < iTracks; i++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  
    if(!track) continue;                                                 
    if (!AcceptTrack(track)) continue;
    Double_t TOFsignal = track->GetTOFsignal();
    Float_t beta          = 0.0;
    beta                  = Beta(track);
    Int_t chargetrack = track->Charge();
    Double_t trackPhi = track->Phi();
    Double_t trackEta = track->Eta();
    Double_t trackPt = track->Pt();

    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    
    if (!(fabs(nSigmakaon)<fSigCut)) continue;
    if (!(fabs(nSigmaTOFkaon)<fSigCut)) continue;

    fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());
    fNsigmaTPCTOFK->Fill(track->Pt(), nSigmaTOFkaon);
    fHistTOFKch->Fill(trackPt, beta);
    fHistKChPtfullRange->Fill(trackPt);

    if (chargetrack > 0) {fHistKPlusPtfullRange->Fill(trackPt);}
    if (chargetrack < 0) {fHistKMinusPtfullRange->Fill(trackPt);}

    if (!(trackPt <= fUpTCut)) continue;
    if (!(trackPt >= fLpTCut)) continue;

    fSelectedKCh->Add((AliAODTrack*)track);

    //fill single particle charged kaon histograms
    fHistKChPhi->Fill(trackPhi);
    fHistKChPt->Fill(trackPt);
    fHistChPhi->Fill(trackPhi);
    fHistChEta->Fill(trackEta);
    fHistChRap->Fill(track->Y());

    if (chargetrack > 0) {
        fTracksCounter->Fill(0.5);
        fSelectedKpos->Add((AliAODTrack*)track);
        fHistKPlusPt->Fill(trackPt);
        fHistKpPhi->Fill(trackPhi);
        fHistPosPhi->Fill(trackPhi);          
        fHistPosEta->Fill(trackEta);          
        fHistPosPhiEta->Fill(trackPhi,trackEta);
        fHistPosRap->Fill(track->Y(0.493));
        }
    if (chargetrack < 0) {
        fTracksCounter->Fill(1.5);
        fSelectedKneg->Add((AliAODTrack*)track);
        fHistKMinusPt->Fill(trackPt);
        fHistKnPhi->Fill(trackPhi);
        fHistNegPhi->Fill(trackPhi);          
        fHistNegEta->Fill(trackEta);          
        fHistNegPhiEta->Fill(trackPhi,trackEta);
        fHistNegRap->Fill(track->Y(0.493));
        }
    fHistPhi->Fill(track->Phi());
    fHistEta->Fill(track->Eta());
}

}

Bool_t isKaonNsigma = kFALSE;

if(fPidpTDependentMethod){

//cout << "entering pid pt dependent" << endl;

for(Int_t i=0; i < iTracks; i++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  
    if(!track) continue;                                                 
    if (!AcceptTrack(track)) continue;

    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

    isKaonNsigma = (IsKaonNSigma3(track->Pt(),nSigmakaon, nSigmaTOFkaon));

    if (!isKaonNsigma) continue;

    Double_t TOFsignal = track->GetTOFsignal();
    Float_t beta          = 0.0;
    beta                  = Beta(track);
    Int_t chargetrack = track->Charge();
    Double_t trackPhi = track->Phi();
    Double_t trackEta = track->Eta();
    Double_t trackPt = track->Pt();

    fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());
    fNsigmaTPCTOFK->Fill(track->Pt(), nSigmaTOFkaon);
    fHistTOFKch->Fill(trackPt, beta);
    fHistKChPtfullRange->Fill(trackPt);

    if (chargetrack > 0) {fHistKPlusPtfullRange->Fill(trackPt);}
    if (chargetrack < 0) {fHistKMinusPtfullRange->Fill(trackPt);}

    if (!(trackPt <= fUpTCut)) continue;
    if (!(trackPt >= fLpTCut)) continue;

    fSelectedKCh->Add((AliAODTrack*)track);

    //fill single particle charged kaon histograms
    fHistKChPhi->Fill(trackPhi);
    fHistKChPt->Fill(trackPt);
    fHistChPhi->Fill(trackPhi);
    fHistChEta->Fill(trackEta);
    fHistChRap->Fill(track->Y());

    if (chargetrack > 0) {
        fTracksCounter->Fill(0.5);
        fSelectedKpos->Add((AliAODTrack*)track);
        fHistKPlusPt->Fill(trackPt);
        fHistKpPhi->Fill(trackPhi);
        fHistPosPhi->Fill(trackPhi);          
        fHistPosEta->Fill(trackEta);          
        fHistPosPhiEta->Fill(trackPhi,trackEta);
        fHistPosRap->Fill(track->Y(0.493));
        }
    if (chargetrack < 0) {
        fTracksCounter->Fill(1.5);
        fSelectedKneg->Add((AliAODTrack*)track);
        fHistKMinusPt->Fill(trackPt);
        fHistKnPhi->Fill(trackPhi);
        fHistNegPhi->Fill(trackPhi);          
        fHistNegEta->Fill(trackEta);          
        fHistNegPhiEta->Fill(trackPhi,trackEta);
        fHistNegRap->Fill(track->Y(0.493));
        }
    fHistPhi->Fill(track->Phi());
    fHistEta->Fill(track->Eta());

        }

    }

Int_t nSelectedKCh = fSelectedKCh->GetEntries();
//cout << " nSelectedKCh is " << nSelectedKCh << endl;

Int_t nSelectedKpos = fSelectedKpos->GetEntries();
//cout << " nSelectedKpos is " << nSelectedKpos << endl;

Int_t nSelectedKneg = fSelectedKneg->GetEntries();
//cout << " nSelectedKneg is " << nSelectedKneg << endl;

//======== Neutral Kaon Selection ==========

fSelectedK0s = new TObjArray;
fSelectedK0s->SetOwner(kTRUE);

for(Int_t j=0; j < nv0s; j++) {
    AliAODv0 *v0=fAOD->GetV0(j);

    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

    // Track cuts for daughter tracks
   	if ( !(IsMyGoodDaughterTrack(pTrack)) || !(IsMyGoodDaughterTrack(nTrack)) ) continue;

    if(pTrack->Charge()==0 || nTrack->Charge()==0) continue;

    Double_t V0Phi = v0->Phi();
    Double_t V0Eta = v0->Eta();
    Double_t V0Pt = v0->Pt();

    if(!v0) continue;
    if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
    if(!AcceptV0(v0, vertex)) continue;
    fHistK0PtfullRange->Fill(V0Pt);

    Double_t pT = v0->Pt();
    if (pT <= fLpTv0Cut || pT >= fUpTv0Cut) continue;

    fTracksCounter->Fill(2.5);
    fSelectedK0s->Add(v0);
    //fill single particle neutral kaon histograms
    fHistK0Pt->Fill(V0Pt);
    fHistK0Phii->Fill(V0Phi);
    fHistK0PhiEta->Fill(V0Phi,V0Eta);          
    fHistK0Phi->Fill(V0Phi); 
    fHistK0Eta->Fill(V0Eta);                  
    fHistPPionPhi->Fill(pTrack->Phi());
    fHistNPionPhi->Fill(nTrack->Phi());
}
Int_t nSelectedK0s = fSelectedK0s->GetEntries();
//cout << " nSelectedK0s is " << nSelectedK0s << endl;
    
//======== KchKch Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKCh->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(i));  
        if(!track1) continue;                            
        if (!AcceptTrack(track1)) continue;

        Int_t chargetrack = track1->Charge();
        Double_t track1Phi = track1->Phi();
        Double_t track1Eta = track1->Eta();

        for(Int_t j(i+1); j < fSelectedKCh->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(j));  
            if(!track2) continue;                            
            if (!AcceptTrack(track2)) continue;

            Double_t track2Phi = track2->Phi();
            Double_t track2Eta = track2->Eta();

            Double_t DEta = fabs(track1Eta - track2Eta);
            Double_t DPhi = fabs(track1Phi - track2Phi);

            FillDPhiHist(DPhi,fHistKChKChPhi,CentV0M);
            Fill2DHist(DPhi,DEta,fHistKChKCh);                                  
        }     
        
    }   
    
//======== K+K- Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKpos->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(i));  
        if(!track1) continue;                            
        if (!AcceptTrack(track1)) continue;
        if (!AcceptPosTrack(track1)) continue;

        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();

        for(Int_t j(0); j < fSelectedKneg->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(i));  
            if(!track2) continue;                            
            if (!AcceptTrack(track2)) continue;
            if (!AcceptNegTrack(track2)) continue;

            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();

            Double_t DPhiPN = fabs(track1PosPhi - track2NegPhi);
            Double_t DEtaPN = fabs(track1PosEta - track2NegEta);

            FillDPhiHist(DPhiPN,fHistKPosKNegPhi,CentV0M);
            Fill2DHist(DPhiPN,DEtaPN,fHistKPosKNeg);
            Fill3DHist(DPhiPN,DEtaPN,PVz,fHistKPosKNegz); 

        }
    }

//======== K+K+ Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKpos->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(i));  
        if(!track1) continue;
        if (!AcceptTrack(track1)) continue;                           
        if (!AcceptPosTrack(track1)) continue;

        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();

        for(Int_t j(i+1); j < fSelectedKpos->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(j));  
            if(!track2) continue;
            if (!AcceptTrack(track2)) continue;                            
            if (!AcceptPosTrack(track2)) continue;

            Double_t track2PosPhi = track2->Phi();
            Double_t track2PosEta = track2->Eta();

            Double_t DPhiPP = fabs(track1PosPhi - track2PosPhi);
            Double_t DEtaPP = fabs(track1PosEta - track2PosEta);

            Fill2DHist(DPhiPP,DEtaPP,fHistKPosKPos);
        }
    }

//======== K-K- Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKneg->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(i));  
        if(!track1) continue;
        if (!AcceptTrack(track1)) continue;                           
        if (!AcceptNegTrack(track1)) continue;

        Double_t track1NegPhi = track1->Phi();
        Double_t track1NegEta = track1->Eta();

        for(Int_t j(i+1); j < fSelectedKneg->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(j));  
            if(!track2) continue;
            if (!AcceptTrack(track2)) continue;                            
            if (!AcceptNegTrack(track2)) continue;

            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();

            Double_t DPhiNN = fabs(track1NegPhi - track2NegPhi);
            Double_t DEtaNN = fabs(track1NegEta - track2NegEta);

            Fill2DHist(DPhiNN,DEtaNN,fHistKNegKNeg);
        }
    }

//======== K0s-Kch Correlation Analysis ==========
    for(Int_t i = 0; i < fSelectedK0s->GetEntriesFast(); i++){
        AliAODv0 * v0 = (AliAODv0*)fSelectedK0s->At(i);
        if(!v0) continue;

        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        Double_t V0Pt = v0->Pt();

        AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
        AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

        Int_t pTrkID = pTrack->GetID()  >= 0 ? pTrack->GetID() : -1-pTrack->GetID();
        Int_t nTrkID = nTrack->GetID() >= 0 ? nTrack->GetID() : -1-nTrack->GetID();

        for (Int_t iTrk=0; iTrk < fSelectedKCh->GetEntriesFast(); iTrk++) {                                       
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(iTrk));  
        if(!track) continue;

        Int_t chargetrack = track->Charge();
        Double_t trackPhi = track->Phi();
        Double_t trackEta = track->Eta();
        Double_t trackPt = track->Pt();

        Int_t trID = track->GetID() >= 0 ? track->GetID() : -1-track->GetID();

        //remove auto correlations
        if( pTrkID == trID || nTrkID == trID ) continue;
        // cout<<"Hello! We are moving auto-correlations here!!"<<endl;                                                 

        Double_t deltaEta = fabs(trackEta-V0Eta);
        Double_t deltaPhi = fabs(trackPhi-V0Phi);

        fHistDEta->Fill(deltaEta);
        fHistDEta->Fill(-deltaEta);

        //FillDPhiHist(deltaPhi,fHistCFPhi,CentV0M);

        if (deltaPhi > Pi) deltaPhi = Pi-(deltaPhi-Pi);

        fHistCFPhi->Fill(deltaPhi, deltaEta);
        fHistCFPhi->Fill(deltaPhi, -deltaEta);
        if (deltaPhi < 0.5*Pi ) {
        fHistCFPhi->Fill(-deltaPhi, deltaEta);
        fHistCFPhi->Fill(-deltaPhi,-deltaEta); 
        }
        else {
        fHistCFPhi->Fill(2*Pi-(deltaPhi), deltaEta);
        fHistCFPhi->Fill(2*Pi-(deltaPhi),-deltaEta);
        }

        Fill2DHist(deltaPhi,deltaEta,fHistCF);
        Fill3DHist(deltaPhi,deltaEta,PVz,fHistCFz); 

        if (deltaPhi > Pi) deltaPhi = Pi-(deltaPhi-Pi);
        fHistDPhi->Fill(deltaPhi);                     
        if (deltaPhi < 0.5*Pi) {fHistDPhi->Fill(-deltaPhi);}
        else {fHistDPhi->Fill(2*Pi-deltaPhi);}

        if (0 < deltaPhi < 0.5*Pi) {    
            fHistCFEta->Fill(deltaEta,CentV0M);
            fHistCFEta->Fill(-deltaEta,CentV0M);  
        }                                  

        }  // iTrk kch loop ends  
    }      // ik0s loop ends
//=====================================================
fHistPVz->Fill(PVz);
fHistNEvents->Fill(5.5);
fHistMult->Fill(iTracks);
fHistCent->Fill(CentV0M); 

//================== mixing ============================ (You can create a seperate function FillCorrelationsMixed() later)

AliEventPool *pool = fPoolMgr->GetEventPool(CentV0M, PVz);
if(pool) {cout << "Good news....!!!!!!! Pool found.. " << endl; }
if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", CentV0M,PVz)); return; }

if (fK0KchCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedK0s->GetEntries(); iTrig++){
            AliVParticle* K0Trig = dynamic_cast<AliVParticle*>(fSelectedK0s->At(iTrig));
            if(!K0Trig) continue;

            for (Int_t iAss(0); iAss < bgTracks->GetEntries(); iAss++){
            AliVParticle* KChAssoc = dynamic_cast<AliVParticle*> (bgTracks->At(iAss));
            //AliAODTrack* KChAssoc = dynamic_cast<AliAODTrack*> (bgTracks->At(iAss));
            if (!KChAssoc) continue;
            if ((jMix==0) && (iTrig==0)) {fHistKChPtMix->Fill(KChAssoc->Pt());}

            Double_t DPhiMix = fabs(K0Trig->Phi() - KChAssoc->Phi());
            Double_t DEtaMix = fabs(K0Trig->Eta() - KChAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_Bg);
            Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_Bgz);

            } //end of kch loop end

        } // end of k0 trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKpKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks2 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonTrig = dynamic_cast<AliVParticle*>(fSelectedKpos->At(iTrig));
            if(!KaonTrig) continue;

            for (Int_t iAss(0); iAss < bgTracks2->GetEntries(); iAss++){

            AliVParticle* KaonAssoc = dynamic_cast<AliVParticle*> (bgTracks2->At(iAss));
            if(!KaonAssoc) continue;
            if( KaonAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KpKn_Bg);
            Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_KpKn_Bgz);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKpKpCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks3 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonPosTrig = dynamic_cast<AliVParticle*>(fSelectedKpos->At(iTrig));
            if(!KaonPosTrig) continue;

            for (Int_t iAss(0); iAss < bgTracks3->GetEntries(); iAss++){

            AliVParticle* KaonPosAssoc = dynamic_cast<AliVParticle*> (bgTracks3->At(iAss));
            if(!KaonPosAssoc) continue;
            if( KaonPosAssoc->Charge() < 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonPosTrig->Phi() - KaonPosAssoc->Phi());
            Double_t DEtaMix = fabs(KaonPosTrig->Eta() - KaonPosAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KpKp_Bg);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKnKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks4 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKneg->GetEntries(); iTrig++){
            AliVParticle* KaonNegTrig = dynamic_cast<AliVParticle*>(fSelectedKneg->At(iTrig));
            if(!KaonNegTrig) continue;

            for (Int_t iAss(0); iAss < bgTracks4->GetEntries(); iAss++){

            AliVParticle* KaonNegAssoc = dynamic_cast<AliVParticle*> (bgTracks4->At(iAss));
            if(!KaonNegAssoc) continue;
            if( KaonNegAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonNegTrig->Phi() - KaonNegAssoc->Phi());
            Double_t DEtaMix = fabs(KaonNegTrig->Eta() - KaonNegAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KnKn_Bg);

            } //end of k- loop end

        } // end of k- trigger loop 

    } //end of mixing event loop 

}//end of pool

}

    TObjArray* tracksClone = (TObjArray*) fSelectedKCh->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);   


} //end of RunData function

//=================== MC Truth Correlation function ========

void AliAnalysisTaskKaon2PC::RunMCTruth() {

fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
if(!fAOD) return;

fMCEvents->Fill(0.5);

//  data trigger selection
Bool_t isSelected = kFALSE;

if(fMinBias) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
if(fCentral) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
if(fSemiCentral)isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
if (!isSelected) return; 

fMCEvents->Fill(1.5);

//==================== Pile Up==========================

fEventCuts.fUseITSTPCCluCorrelationCut = true;
if (fRejectEventPileUp){
    if (!fEventCuts.AcceptEvent(fAOD)) return;

}

fMCEvents->Fill(2.5);

Double_t CentV0M = fEventCuts.GetCentrality(); //centrality for MC
if ((CentV0M < fCentMin)||(CentV0M > fCentMax)) return;
cout << "centrality values for MC gen are" << CentV0M << endl;

Int_t nAcceptedParticles =0;
AliMCParticle *mcTrack = 0x0;

fmcEvent  = dynamic_cast<AliMCEvent*> (MCEvent());
if(!fmcEvent){
    Printf("No MC particle branch found");
    return;
    }

AliAODMCHeader *mcHeader = 0;
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
        return;
    }

// Float_t vzMC = mcHeader->GetVtxZ();
// cout << "pvz values from mcheader" << vzMC << endl;
// if (TMath::Abs(vzMC) >= fPVzCut) return;

//AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(fmcEvent,"Hijing") return; 
//AliAnalysisUtils::IsPileupInGeneratedEvent(fmcEvent,"Hijing") return;

AliVVertex * mcVertex = (AliVVertex*)fmcEvent->GetPrimaryVertex();
Double_t vzMC = mcVertex->GetZ();
if (TMath::Abs(vzMC) >= fPVzCut) return;
cout << "pvz values from truth are" << vzMC << endl;

fMCEvents->Fill(3.5);

//retreive MC particles from event 
fMCArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
if(!fMCArray){
    Printf("No MC particle branch found");
    return;
}

//Int_t nMCTracks = fmcEvent->GetNumberOfTracks(); // MC Truth, Total number of tracks per event

Int_t nMCAllTracks =fMCArray->GetEntriesFast();
cout << "mc tracks without" << nMCAllTracks << endl;
// new tracks array - without injected signal
TObjArray * mcTracks = new TObjArray;
AliVTrack *genTrackMix = 0x0;

for (Int_t i = 0; i < nMCAllTracks; i++){
    AliAODMCParticle *mcTrack = (AliAODMCParticle*)fMCArray->At(i);
    if (!mcTrack) {
        Error("ReadEventAODMC", "Could not receive particle %d", i);
        continue;
    }

    mcTracks->Add(mcTrack);
    }

Int_t nMCTracks = mcTracks->GetEntriesFast();

for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!mcTrack) {
        Error("ReadEventAODMC", "Could not receive particle %d", i);
        continue;
    }

    Double_t trackPseudorap = mcTrack->Eta();
    if( mcTrack->IsPhysicalPrimary()&&mcTrack->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8))) {
                nAcceptedParticles += 1;
            }
}

cout << "number of accepted particles from MC tracks is"<< nAcceptedParticles << endl;

fMCEvents->Fill(4.5);
fHistGenMultiplicity->Fill(nAcceptedParticles);

fHistCent_mcgen->Fill(CentV0M);

AliMCParticle *mcMotherParticle = 0x0;
AliMCParticle* daughter0 = 0x0;
AliMCParticle* daughter1 = 0x0;
Bool_t SelectK0;
Bool_t SelectKpos;
Bool_t SelectKneg;
Bool_t SelectKch;

fMCSelectedK0s = new TObjArray;
fMCSelectedK0s->SetOwner(kTRUE);

fMCSelectedKCh = new TObjArray;
fMCSelectedKCh->SetOwner(kTRUE);

fMCSelectedKpos = new TObjArray;
fMCSelectedKpos->SetOwner(kTRUE);

fMCSelectedKneg = new TObjArray;
fMCSelectedKneg->SetOwner(kTRUE);

fMCEvents->Fill(5.5); 

for (Int_t i = 0; i < nMCTracks; i++){
	mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);

    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fmcEvent)) continue; 

    if (!mcTrack) continue;

    if(! (fabs(mcTrack->Eta()) < fEtaCut)) continue; 
        
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();

    SelectK0 =  mcPartPdg==310&& (isPhysPrim);  // (8010 total. When we include only primary, 7546)-310 for neutral
    SelectKpos = mcPartPdg==321&& (isPhysPrim); // 321 is the code for positive kaons
    SelectKneg = mcPartPdg==-321&& (isPhysPrim); // 321 is the code for positive kaons
    SelectKch = (mcPartPdg==321 || mcPartPdg==-321) && (isPhysPrim);

    Double_t TrackPt = mcTrack->Pt();
    Double_t TrackPhi = mcTrack->Phi();
    Double_t TrackEta = mcTrack->Eta();
    Double_t TrackMass = mcTrack->M();

    
    if(SelectK0) fMCK0PtfullRange->Fill(TrackPt);
    if(SelectKpos) fMCKPlusPtfullRange->Fill(TrackPt);
    if(SelectKneg) fMCKMinusPtfullRange->Fill(TrackPt);

    if (!(TrackPt <= fUpTCut)) continue;
    if (!(TrackPt >= fLpTCut)) continue;
    
    Double_t KaonVariables[4]= {TrackPt, TrackPhi, TrackEta, TrackMass};
    if(SelectK0) {
        fMCK0->Fill(KaonVariables);
        fMCK0Pt->Fill(TrackPt);
    }
    if(SelectKpos) {
        fMCKpos->Fill(KaonVariables);
        fMCKPlusPt->Fill(TrackPt);
    }
    
    if(SelectKneg) {
        fMCKneg->Fill(KaonVariables);
        fMCKMinusPt->Fill(TrackPt);
    }
    if(SelectKch) fMCKch->Fill(KaonVariables);

    if(SelectK0) fMCSelectedK0s->Add(mcTrack);
    if(SelectKch) fMCSelectedKCh->Add(mcTrack);
    if(SelectKpos) fMCSelectedKpos->Add(mcTrack);
    if(SelectKneg) fMCSelectedKneg->Add(mcTrack);

	Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
	Bool_t TrCharge = (mcTrack->Charge())!=0;
    Short_t cha;
    if (mcTrack->Charge()>0) cha=1.;
    else if (mcTrack->Charge()<0) cha= -1.;
    else cha =0;

}
//Int_t nmck0 =  fMCSelectedK0s->GetEntries();
//cout << "nmck0 is" << nmck0 << endl;

//============ MC Generated K0-Kch Correlation ============

for (Int_t i = 0; i < fMCSelectedK0s->GetEntries(); i++){
    AliVParticle* mcTrack1 = dynamic_cast<AliVParticle*>(fMCSelectedK0s->At(i));
    //if (!SelectK0TracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    //if(! (mcTrack1->IsPhysicalPrimary()&&mcTrack1->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = 0; j < fMCSelectedKCh->GetEntries(); j++){
        AliVParticle* mcTrack2 = dynamic_cast<AliVParticle*>(fMCSelectedKCh->At(j));
        //if (!SelectKchTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();

        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistK0KchMC);
    }

}

//============ MC Generated K+K- Correlation ============

for (Int_t i = 0; i < fMCSelectedKpos->GetEntries(); i++){
    AliVParticle* mcTrack1 = dynamic_cast<AliVParticle*>(fMCSelectedKpos->At(i));
    if (!mcTrack1) continue;
    //if (!SelectKPosTracksMC(mcTrack1)) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = 0; j < fMCSelectedKneg->GetEntries(); j++){
        AliVParticle* mcTrack2 = dynamic_cast<AliVParticle*>(fMCSelectedKneg->At(j));
        if (!mcTrack2) continue;
        //if (!SelectKNegTracksMC(mcTrack2)) continue;
        Double_t trackPseudorap = mcTrack2->Eta();
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistKpKnMC);
    }

}

//============ MC Generated K+K+ Correlation ============

for (Int_t i = 0; i < fMCSelectedKpos->GetEntries(); i++){
    AliVParticle* mcTrack1 = dynamic_cast<AliVParticle*>(fMCSelectedKpos->At(i));
    //if (!SelectKPosTracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < fMCSelectedKpos->GetEntries(); j++){
        AliVParticle* mcTrack2 = dynamic_cast<AliVParticle*>(fMCSelectedKpos->At(j));
        //if (!SelectKPosTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistKpKpMC);
    }

}

//============ MC Generated K-K- Correlation ============

for (Int_t i = 0; i < fMCSelectedKneg->GetEntries(); i++){
    AliVParticle* mcTrack1 = dynamic_cast<AliVParticle*>(fMCSelectedKneg->At(i));
    //if (!SelectKNegTracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < fMCSelectedKneg->GetEntries(); j++){
        AliVParticle* mcTrack2 = dynamic_cast<AliVParticle*>(fMCSelectedKneg->At(j));
        //if (!SelectKNegTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistKnKnMC);
    }

}


//================== MC Gen mixing ============================ 

AliEventPool *pool = fPoolMgr->GetEventPool(CentV0M, vzMC);
if(pool) {cout << "Good news....!!!!!!! Pool found.. " << endl; }
if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", CentV0M,vzMC)); return; }

if (fK0KchCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fMCSelectedK0s->GetEntries(); iTrig++){
            AliVParticle* K0Trig = dynamic_cast<AliVParticle*>(fMCSelectedK0s->At(iTrig));
            if(!K0Trig) continue;

            for (Int_t iAss(0); iAss < bgTracks->GetEntries(); iAss++){
            AliVParticle* KChAssoc = dynamic_cast<AliVParticle*> (bgTracks->At(iAss));
            //AliAODTrack* KChAssoc = dynamic_cast<AliAODTrack*> (bgTracks->At(iAss));
            if (!KChAssoc) continue;

            Double_t DPhiMix = fabs(K0Trig->Phi() - KChAssoc->Phi());
            Double_t DEtaMix = fabs(K0Trig->Eta() - KChAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistK0KchMC_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_Bgz);

            } //end of kch MC Gen loop end

        } // end of k0 MC Gen trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKpKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks2 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fMCSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonTrig = dynamic_cast<AliVParticle*>(fMCSelectedKpos->At(iTrig));
            if(!KaonTrig) continue;

            for (Int_t iAss(0); iAss < bgTracks2->GetEntries(); iAss++){
            AliVParticle* KaonAssoc = dynamic_cast<AliVParticle*> (bgTracks2->At(iAss));
            if(!KaonAssoc) continue;
            if( KaonAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistKpKnMC_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_KpKn_Bgz);

            } //end of k- loop end

        } // end of k+ trigger loop []

    } //end of mixing event loop 

}//end of pool

}

if (fKpKpCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks2 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fMCSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonTrig = dynamic_cast<AliVParticle*>(fMCSelectedKpos->At(iTrig));
            if(!KaonTrig) continue;
            if( KaonTrig->Charge() < 0.0 ) continue;

            for (Int_t iAss(iTrig+1); iAss < bgTracks2->GetEntries(); iAss++){
            AliVParticle* KaonAssoc = dynamic_cast<AliVParticle*> (bgTracks2->At(iAss));
            if(!KaonAssoc) continue;
            if( KaonAssoc->Charge() < 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistKpKpMC_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_KpKn_Bgz);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKnKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks2 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fMCSelectedKneg->GetEntries(); iTrig++){
            AliVParticle* KaonTrig = dynamic_cast<AliVParticle*>(fMCSelectedKneg->At(iTrig));
            if(!KaonTrig) continue;
            if( KaonTrig->Charge() > 0.0 ) continue;

            for (Int_t iAss(iTrig+1); iAss < bgTracks2->GetEntries(); iAss++){
            AliVParticle* KaonAssoc = dynamic_cast<AliVParticle*> (bgTracks2->At(iAss));
            if(!KaonAssoc) continue;
            if( KaonAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistKnKnMC_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_KpKn_Bgz);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}


TObjArray* tracksClone = (TObjArray*) fMCSelectedKCh->Clone();
tracksClone->SetOwner(kTRUE);
pool->UpdatePool(tracksClone); 

}  // end of MCTruth function

//=================== MC Reconstructed Correlation function ========

void AliAnalysisTaskKaon2PC::RunMCReconstructed() {

    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;

    fHistNEvents->Fill(0.5);

    //  data trigger selection
	Bool_t isSelected = kFALSE;

    if(fMinBias) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    if(fCentral) isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
    if(fSemiCentral)isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
	if (!isSelected) return;    

    fHistNEvents->Fill(1.5);                               

    //==================== Pile Up==========================

    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (fRejectEventPileUp){
        if (!fEventCuts.AcceptEvent(fAOD)) return;
    }

    fHistNEvents->Fill(2.5);

    //primary Vertex    
    const AliVVertex* primVertex = fEventCuts.GetPrimaryVertex(); 
    if (!primVertex) return;
    
    Double_t PVx = primVertex->GetX();
    Double_t PVy = primVertex->GetY();
    Double_t PVz = primVertex->GetZ();

    Double_t vertex[3] = { -100.0, -100.0, -100.0 }; //?
    primVertex->GetXYZ(vertex);  

    //zvertex cut
    if ( ( TMath::Abs(PVz) ) >= fPVzCut) return ;
    fVtx->Fill(PVz);

    fHistNEvents->Fill(3.5);

    //Multiplicity selection
    //AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    //if(!MultSelection) return;
    
    //centrality
    //double CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    Double_t CentV0M = fEventCuts.GetCentrality(); //centrality
    cout << "centrality of event is " << CentV0M << endl;
    if ((CentV0M < fCentMin)||(CentV0M > fCentMax)) return;

    fHistNEvents->Fill(4.5);

//======== PID loop (no pT cut) ===========

    Int_t iTracks = (fAOD->GetNumberOfTracks());
    cout << "Number of tracks in an event is " << iTracks << endl;
    Double_t MultTPC = fAOD->GetNumberOfTPCTracks();

    Int_t nv0s = (fAOD->GetNumberOfV0s());
    Double_t tpcClusters = fAOD->GetNumberOfTPCClusters();
    
    fHistTPCTracksVsClusters->Fill(MultTPC,tpcClusters);
    //fHistV0AmplitudeVsPVposition->Fill(PVz,VZEROmultiplicity);           
    
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));       
        Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        Double_t nSigmaproton =TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) ;
        Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        Double_t nClustersTPC = track->GetTPCNcls();
        Double_t kchPt = track->Pt();
        Double_t charge = track->Charge();
        if(!track) continue;                               // if we failed, skip this track
        if(abs(track->Eta())>0.8) continue;                // eta cut
        if (!track->TestFilterBit(96)) continue;           // filterbit selection
        if (track->Charge() == 0) continue;  
        fClusters->Fill(nClustersTPC);
        fEnergy->Fill(track->GetTPCsignal());
        fPID->Fill(track->Pt(),track->GetTPCsignal());
        if (fabs(nSigmaelectron) < 2.0) continue;          // excluding electrons via TPC
        if (fabs(nSigmapion) < 2.0) continue;              // excluding pions via TPC
        if (fabs(nSigmaproton) < 2.0) continue;              // excluding protons via TPC
        if (fabs(nSigmakaon)<2.0) {
            fPIDK->Fill(track->Pt(),track->GetTPCsignal());
            fNsigmaKaon->Fill(track->Pt(), nSigmakaon);
            fEnergyCuts->Fill(track->GetTPCsignal()); }
        if ((fabs(nSigmakaon)<2.0) && (0.0 < kchPt< 2.0)) {hPt->Fill(kchPt);}
        if ((fabs(nSigmakaon)<2.0) && (0.0 < kchPt< 2.0) && (charge > 0)) {hPt_kPos->Fill(kchPt);}
        if (fabs(nSigmakaon) < 20.0 ) {fnsigmakaon->Fill(track->Pt(), nSigmakaon);}
        fNsigmaTOFK->Fill(track->Pt(), nSigmaTOFkaon);
        if (fabs(nSigmaTOFkaon)<3.0) {fNsigmaTOFKaon->Fill(track->Pt(), nSigmaTOFkaon);}
        if ((fabs(nSigmakaon)<3.0) && (fabs(nSigmaelectron) > 3.0 ) && (fabs(nSigmapion) > 3.0 ) &&
            (fabs(nSigmaTOFelectron) > 3.0 ) && (fabs(nSigmaTOFpion) > 3.0 ) && (fabs(nSigmaTOFkaon)<3.0)){
            fPIDKpiCut->Fill(track->Pt(),track->GetTPCsignal()); }
         
    }   //end of PID loop

//======== v0 loop (no invariant mass cut) =========

    fHistNV0->Fill(nv0s); 
    for(Int_t i = 0; i < nv0s; i++)  {
        AliAODv0 *v0=fAOD->GetV0(i);                    
        if(!v0) {
        cout<<"No V0 "<<endl;
        continue;
        }
        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        fHistMK0->Fill(v0->MassK0Short());
        if(!AcceptV0(v0, vertex)) continue;
        Double_t pT = v0->Pt();
        if (pT <= fLpTv0Cut || pT >= fUpTv0Cut) continue;
        fHistMK0Cuts->Fill(v0->MassK0Short());
    }

//======== Charged Kaon Track Selection  ==========
fSelectedKCh = new TObjArray;
fSelectedKCh->SetOwner(kTRUE);

fSelectedKpos = new TObjArray;
fSelectedKpos->SetOwner(kTRUE);

fSelectedKneg = new TObjArray;
fSelectedKneg->SetOwner(kTRUE);

if(!fPidpTDependentMethod){


for(Int_t i=0; i < iTracks; i++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  
    if(!track) continue;                                                 
    if (!AcceptTrack(track)) continue;
    Double_t TOFsignal = track->GetTOFsignal();
    Float_t beta          = 0.0;
    beta                  = Beta(track);
    Int_t chargetrack = track->Charge();
    Double_t trackPhi = track->Phi();
    Double_t trackEta = track->Eta();
    Double_t trackPt = track->Pt();

    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    
    if (!(fabs(nSigmakaon)<fSigCut)) continue;
    if (!(fabs(nSigmaTOFkaon)<fSigCut)) continue;

    fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());
    fNsigmaTPCTOFK->Fill(track->Pt(), nSigmaTOFkaon);
    fHistTOFKch->Fill(trackPt, beta);
    fHistKChPtfullRange->Fill(trackPt);

    if (chargetrack > 0) {fHistKPlusPtfullRange->Fill(trackPt);}
    if (chargetrack < 0) {fHistKMinusPtfullRange->Fill(trackPt);}

    if (!(trackPt <= fUpTCut)) continue;
    if (!(trackPt >= fLpTCut)) continue;

    fSelectedKCh->Add((AliAODTrack*)track);

    //fill single particle charged kaon histograms
    fHistKChPhi->Fill(trackPhi);
    fHistKChPt->Fill(trackPt);
    fHistChPhi->Fill(trackPhi);
    fHistChEta->Fill(trackEta);
    fHistChRap->Fill(track->Y());

    if (chargetrack > 0) {
        fTracksCounter->Fill(0.5);
        fSelectedKpos->Add((AliAODTrack*)track);
        fHistKPlusPt->Fill(trackPt);
        fHistKpPhi->Fill(trackPhi);
        fHistPosPhi->Fill(trackPhi);          
        fHistPosEta->Fill(trackEta);          
        fHistPosPhiEta->Fill(trackPhi,trackEta);
        fHistPosRap->Fill(track->Y(0.493));
        }
    if (chargetrack < 0) {
        fTracksCounter->Fill(1.5);
        fSelectedKneg->Add((AliAODTrack*)track);
        fHistKMinusPt->Fill(trackPt);
        fHistKnPhi->Fill(trackPhi);
        fHistNegPhi->Fill(trackPhi);          
        fHistNegEta->Fill(trackEta);          
        fHistNegPhiEta->Fill(trackPhi,trackEta);
        fHistNegRap->Fill(track->Y(0.493));
        }
    fHistPhi->Fill(track->Phi());
    fHistEta->Fill(track->Eta());
}

}

Bool_t isKaonNsigma = kFALSE;

if(fPidpTDependentMethod){

//cout << "entering pid pt dependent" << endl;

for(Int_t i=0; i < iTracks; i++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  
    if(!track) continue;                                                 
    if (!AcceptTrack(track)) continue;

    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

    isKaonNsigma = (IsKaonNSigma3(track->Pt(),nSigmakaon, nSigmaTOFkaon));

    if (!isKaonNsigma) continue;

    Double_t TOFsignal = track->GetTOFsignal();
    Float_t beta          = 0.0;
    beta                  = Beta(track);
    Int_t chargetrack = track->Charge();
    Double_t trackPhi = track->Phi();
    Double_t trackEta = track->Eta();
    Double_t trackPt = track->Pt();
    fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());
    fNsigmaTPCTOFK->Fill(track->Pt(), nSigmaTOFkaon);
    fHistTOFKch->Fill(trackPt, beta);
    fHistKChPtfullRange->Fill(trackPt);

    if (chargetrack > 0) {fHistKPlusPtfullRange->Fill(trackPt);}
    if (chargetrack < 0) {fHistKMinusPtfullRange->Fill(trackPt);}

    if (!(trackPt <= fUpTCut)) continue;
    if (!(trackPt >= fLpTCut)) continue;

    fSelectedKCh->Add((AliAODTrack*)track);

    //fill single particle charged kaon histograms
    fHistKChPhi->Fill(trackPhi);
    fHistKChPt->Fill(trackPt);
    fHistChPhi->Fill(trackPhi);
    fHistChEta->Fill(trackEta);
    fHistChRap->Fill(track->Y());

    if (chargetrack > 0) {
        fTracksCounter->Fill(0.5);
        fSelectedKpos->Add((AliAODTrack*)track);
        fHistKPlusPt->Fill(trackPt);
        fHistKpPhi->Fill(trackPhi);
        fHistPosPhi->Fill(trackPhi);          
        fHistPosEta->Fill(trackEta);          
        fHistPosPhiEta->Fill(trackPhi,trackEta);
        fHistPosRap->Fill(track->Y(0.493));
        }
    if (chargetrack < 0) {
        fTracksCounter->Fill(1.5);
        fSelectedKneg->Add((AliAODTrack*)track);
        fHistKMinusPt->Fill(trackPt);
        fHistKnPhi->Fill(trackPhi);
        fHistNegPhi->Fill(trackPhi);          
        fHistNegEta->Fill(trackEta);          
        fHistNegPhiEta->Fill(trackPhi,trackEta);
        fHistNegRap->Fill(track->Y(0.493));
        }
    fHistPhi->Fill(track->Phi());
    fHistEta->Fill(track->Eta());

        }

    }

Int_t nSelectedKCh = fSelectedKCh->GetEntries();
//cout << " nSelectedKCh is " << nSelectedKCh << endl;

Int_t nSelectedKpos = fSelectedKpos->GetEntries();
//cout << " nSelectedKpos is " << nSelectedKpos << endl;

Int_t nSelectedKneg = fSelectedKneg->GetEntries();
//cout << " nSelectedKneg is " << nSelectedKneg << endl;

//======== Neutral Kaon Selection ==========

fSelectedK0s = new TObjArray;
fSelectedK0s->SetOwner(kTRUE);

for(Int_t j=0; j < nv0s; j++) {
    AliAODv0 *v0=fAOD->GetV0(j);

    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

    // Track cuts for daughter tracks
   	if ( !(IsMyGoodDaughterTrack(pTrack)) || !(IsMyGoodDaughterTrack(nTrack)) ) continue;

    if(pTrack->Charge()==0 || nTrack->Charge()==0) continue;

    Double_t V0Phi = v0->Phi();
    Double_t V0Eta = v0->Eta();
    Double_t V0Pt = v0->Pt();

    if(!v0) continue;
    if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
    if(!AcceptV0(v0, vertex)) continue;
    fHistK0PtfullRange->Fill(V0Pt);

    Double_t pT = v0->Pt();
    if (pT <= fLpTv0Cut || pT >= fUpTv0Cut) continue;

    fTracksCounter->Fill(2.5);
    fSelectedK0s->Add(v0);
    //fill single particle neutral kaon histograms
    fHistK0Pt->Fill(V0Pt);
    fHistK0Phii->Fill(V0Phi);
    fHistK0PhiEta->Fill(V0Phi,V0Eta);          
    fHistK0Phi->Fill(V0Phi); 
    fHistK0Eta->Fill(V0Eta);                  
    fHistPPionPhi->Fill(pTrack->Phi());
    fHistNPionPhi->Fill(nTrack->Phi());
}
Int_t nSelectedK0s = fSelectedK0s->GetEntries();
//cout << " nSelectedK0s is " << nSelectedK0s << endl;

//======== KchKch Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKCh->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(i));  
        if(!track1) continue;                            
        if (!AcceptTrack(track1)) continue;

        Int_t chargetrack = track1->Charge();
        Double_t track1Phi = track1->Phi();
        Double_t track1Eta = track1->Eta();

        for(Int_t j(i+1); j < fSelectedKCh->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(j));  
            if(!track2) continue;                            
            if (!AcceptTrack(track2)) continue;

            Double_t track2Phi = track2->Phi();
            Double_t track2Eta = track2->Eta();

            Double_t DEta = fabs(track1Eta - track2Eta);
            Double_t DPhi = fabs(track1Phi - track2Phi);

            FillDPhiHist(DPhi,fHistKChKChPhi,CentV0M);
            Fill2DHist(DPhi,DEta,fHistKChKCh);                                  
        }     
        
    }   
    
//======== K+K- Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKpos->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(i));  
        if(!track1) continue;                            
        if (!AcceptTrack(track1)) continue;
        if (!AcceptPosTrack(track1)) continue;

        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();

        for(Int_t j(0); j < fSelectedKneg->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(j));  
            if(!track2) continue;                            
            if (!AcceptTrack(track2)) continue;
            if (!AcceptNegTrack(track2)) continue;

            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();

            Int_t labelTrig=track1-> GetLabel();
            //cout << "label of track 1 neg is" << labelTrig << endl;
            Int_t labelAssoc=track2-> GetLabel();
            //cout << "label of track 2 neg is" << labelAssoc << endl; 
            //if(labelTrig==labelAssoc) continue;

            Double_t DPhiPN = fabs(track1PosPhi - track2NegPhi);
            Double_t DEtaPN = fabs(track1PosEta - track2NegEta);

            FillDPhiHist(DPhiPN,fHistKPosKNegPhi,CentV0M);
            Fill2DHist(DPhiPN,DEtaPN,fHistKPosKNeg);
            Fill3DHist(DPhiPN,DEtaPN,PVz,fHistKPosKNegz); 

        }
    }
//======== K+K+ Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKpos->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(i));  
        if(!track1) continue;
        if (!AcceptTrack(track1)) continue;                           
        if (!AcceptPosTrack(track1)) continue;

        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();

        for(Int_t j(0); j < fSelectedKpos->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKpos->At(j)); 

            if(!track2) continue;
            if (!AcceptTrack(track2)) continue;                            
            if (!AcceptPosTrack(track2)) continue;

            Int_t labelTrig=track1->GetLabel();
            //cout << "label of track 1 is" << labelTrig << endl;
            Int_t labelAssoc=track2->GetLabel();
            //cout << "label of track 2 is" << labelAssoc << endl; 
            if(labelTrig==labelAssoc) continue;

            Double_t track2PosPhi = track2->Phi();
            Double_t track2PosEta = track2->Eta();

            Double_t DPhiPP = fabs(track1PosPhi - track2PosPhi);
            Double_t DEtaPP = fabs(track1PosEta - track2PosEta);

            Fill2DHist(DPhiPP,DEtaPP,fHistKPosKPos);
        }
    }

//======== K-K- Correlation Loop ==========

    for(Int_t i(0); i < fSelectedKneg->GetEntriesFast(); i++) {                 
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(i));  
        if(!track1) continue;
        if (!AcceptTrack(track1)) continue;                           
        if (!AcceptNegTrack(track1)) continue;

        Double_t track1NegPhi = track1->Phi();
        Double_t track1NegEta = track1->Eta();

        for(Int_t j(0); j < fSelectedKneg->GetEntriesFast(); j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fSelectedKneg->At(j));  
            if(!track2) continue;
            if (!AcceptTrack(track2)) continue;                            
            if (!AcceptNegTrack(track2)) continue;

            Int_t labelTrig=track1-> GetLabel();
            //cout << "label of track 1 neg is" << labelTrig << endl;
            Int_t labelAssoc=track2-> GetLabel();
            //cout << "label of track 2 neg is" << labelAssoc << endl; 
            if(labelTrig==labelAssoc) continue;

            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();

            double DPhiNN = fabs(track1NegPhi - track2NegPhi);
            double DEtaNN = fabs(track1NegEta - track2NegEta);
            //cout << "dphi is" << DPhiNN << endl;

            Fill2DHist(DPhiNN,DEtaNN,fHistKNegKNeg);
        }
    }

    //======== K0s-Kch Correlation Analysis ==========
    for(Int_t i = 0; i < fSelectedK0s->GetEntriesFast(); i++){
        AliAODv0 * v0 = (AliAODv0*)fSelectedK0s->At(i);
        if(!v0) continue;

        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        Double_t V0Pt = v0->Pt();

        AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
        AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

        Int_t pTrkID = pTrack->GetID()  >= 0 ? pTrack->GetID() : -1-pTrack->GetID();
        Int_t nTrkID = nTrack->GetID() >= 0 ? nTrack->GetID() : -1-nTrack->GetID();

        for (Int_t iTrk=0; iTrk < fSelectedKCh->GetEntriesFast(); iTrk++) {                                       
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fSelectedKCh->At(iTrk));  
        if(!track) continue;

        Int_t chargetrack = track->Charge();
        Double_t trackPhi = track->Phi();
        Double_t trackEta = track->Eta();
        Double_t trackPt = track->Pt();

        Int_t trID = track->GetID() >= 0 ? track->GetID() : -1-track->GetID();

        //remove auto correlations
        if( pTrkID == trID || nTrkID == trID ) continue;
        // cout<<"Hello! We are moving auto-correlations here!!"<<endl;                                                 

        Double_t deltaEta = fabs(trackEta-V0Eta);
        Double_t deltaPhi = fabs(trackPhi-V0Phi);

        fHistDEta->Fill(deltaEta);
        fHistDEta->Fill(-deltaEta);

        //FillDPhiHist(deltaPhi,fHistCFPhi,CentV0M);

        if (deltaPhi > Pi) deltaPhi = Pi-(deltaPhi-Pi);

        fHistCFPhi->Fill(deltaPhi, deltaEta);
        fHistCFPhi->Fill(deltaPhi, -deltaEta);
        if (deltaPhi < 0.5*Pi ) {
        fHistCFPhi->Fill(-deltaPhi, deltaEta);
        fHistCFPhi->Fill(-deltaPhi,-deltaEta); 
        }
        else {
        fHistCFPhi->Fill(2*Pi-(deltaPhi), deltaEta);
        fHistCFPhi->Fill(2*Pi-(deltaPhi),-deltaEta);
        }

        Fill2DHist(deltaPhi,deltaEta,fHistCF);
        Fill3DHist(deltaPhi,deltaEta,PVz,fHistCFz); 

        if (deltaPhi > Pi) deltaPhi = Pi-(deltaPhi-Pi);
        fHistDPhi->Fill(deltaPhi);                     
        if (deltaPhi < 0.5*Pi) {fHistDPhi->Fill(-deltaPhi);}
        else {fHistDPhi->Fill(2*Pi-deltaPhi);}

        if (0 < deltaPhi < 0.5*Pi) {    
            fHistCFEta->Fill(deltaEta,CentV0M);
            fHistCFEta->Fill(-deltaEta,CentV0M);  
        }                                  

        }  // iTrk kch loop ends  
    }      // ik0s loop ends
//=====================================================
fHistPVz->Fill(PVz);
fHistMult->Fill(iTracks);
fHistCent->Fill(CentV0M);
fHistNEvents->Fill(5.5); 

//================== mixing ============================ 

AliEventPool *pool = fPoolMgr->GetEventPool(CentV0M, PVz);
if(pool) {cout << "Good news....!!!!!!! Pool found.. " << endl; }
if(!pool) { AliError(Form("No pool found for centrality = %f, zVtx = %f", CentV0M,PVz)); return; }

if (fK0KchCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedK0s->GetEntries(); iTrig++){
            AliVParticle* K0Trig = dynamic_cast<AliVParticle*>(fSelectedK0s->At(iTrig));
            if(!K0Trig) continue;

            for (Int_t iAss(0); iAss < bgTracks->GetEntries(); iAss++){
            AliVParticle* KChAssoc = dynamic_cast<AliVParticle*> (bgTracks->At(iAss));
            if (!KChAssoc) continue;
            if ((jMix==0) && (iTrig==0)) {fHistKChPtMix->Fill(KChAssoc->Pt());}

            Double_t DPhiMix = fabs(K0Trig->Phi() - KChAssoc->Phi());
            Double_t DEtaMix = fabs(K0Trig->Eta() - KChAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_Bgz);

            } //end of kch loop end

        } // end of k0 trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKpKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks2 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonTrig = dynamic_cast<AliVParticle*>(fSelectedKpos->At(iTrig));
            if(!KaonTrig) continue;

            for (Int_t iAss(0); iAss < bgTracks2->GetEntries(); iAss++){

            AliVParticle* KaonAssoc = dynamic_cast<AliVParticle*> (bgTracks2->At(iAss));
            if(!KaonAssoc) continue;
            if( KaonAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KpKn_Bg);
            //Fill3DHist(DPhiMix,DEtaMix,PVz,fHistCF_KpKn_Bgz);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKpKpCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks3 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKpos->GetEntries(); iTrig++){
            AliVParticle* KaonPosTrig = dynamic_cast<AliVParticle*>(fSelectedKpos->At(iTrig));
            if(!KaonPosTrig) continue;

            for (Int_t iAss(iTrig+1); iAss < bgTracks3->GetEntries(); iAss++){

            AliVParticle* KaonPosAssoc = dynamic_cast<AliVParticle*> (bgTracks3->At(iAss));
            if(!KaonPosAssoc) continue;
            if( KaonPosAssoc->Charge() < 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonPosTrig->Phi() - KaonPosAssoc->Phi());
            Double_t DEtaMix = fabs(KaonPosTrig->Eta() - KaonPosAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KpKp_Bg);

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

if (fKnKnCorr){

if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() >= fMinEventsToMix) {

Int_t nMix = pool->GetCurrentNEvents();

    for (Int_t jMix=0; jMix< nMix; jMix++){
        TObjArray* bgTracks4 = pool->GetEvent(jMix);  //bgTracks are from the mixed events

        for(Int_t iTrig(0); iTrig < fSelectedKneg->GetEntries(); iTrig++){
            AliVParticle* KaonNegTrig = dynamic_cast<AliVParticle*>(fSelectedKneg->At(iTrig));
            if(!KaonNegTrig) continue;

            for (Int_t iAss(iTrig+1); iAss < bgTracks4->GetEntries(); iAss++){

            AliVParticle* KaonNegAssoc = dynamic_cast<AliVParticle*> (bgTracks4->At(iAss));
            if(!KaonNegAssoc) continue;
            if( KaonNegAssoc->Charge() > 0.0 ) continue;

            Double_t DPhiMix = fabs(KaonNegTrig->Phi() - KaonNegAssoc->Phi());
            Double_t DEtaMix = fabs(KaonNegTrig->Eta() - KaonNegAssoc->Eta());

            Fill2DHist(DPhiMix,DEtaMix,fHistCF_KnKn_Bg);

            } //end of k- loop end

        } // end of k- trigger loop 

    } //end of mixing event loop 

}//end of pool

}

    TObjArray* tracksClone = (TObjArray*) fSelectedKCh->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);   

} //end of MC RunConstructed function

//_____________________________________________________________________________

void AliAnalysisTaskKaon2PC::UserExec(Option_t *)
{
    if (fData) { RunData(); }
    if (fMCTruth) { RunMCTruth(); }
    if (fMCReconstructed) {RunMCReconstructed(); }

    //deleting TObjArrays
    // fSelectedK0s->Clear();
    // fSelectedKpos->Clear();
    // fSelectedKneg->Clear();
    // fSelectedKCh->Clear();
    // fMCSelectedK0s->Clear();
    // fMCSelectedKCh->Clear();
    // fMCSelectedKneg->Clear();
    // fMCSelectedKpos->Clear();

    PostData(1, fOutputList);
}

//====================  Filling functions  =========================

void AliAnalysisTaskKaon2PC::Fill2DHist(Double_t DPhi, Double_t DEta, TH2F* hist){
    DPhi = fabs(DPhi);
    DEta = fabs(DEta);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,DEta);
    hist->Fill(DPhi,-DEta);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi, DEta);
        hist->Fill(-DPhi,-DEta); 
    }
    else {
        hist->Fill(2*Pi-(DPhi), DEta);
        hist->Fill(2*Pi-(DPhi),-DEta);
    }
}

void AliAnalysisTaskKaon2PC::Fill3DHist(Double_t DPhi, Double_t DEta, Double_t PVz, TH3F* hist){
    DPhi = fabs(DPhi);
    DEta = fabs(DEta);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,DEta, PVz);
    hist->Fill(DPhi,-DEta, PVz);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi, DEta, PVz);
        hist->Fill(-DPhi,-DEta, PVz); 
    }
    else {
        hist->Fill(2*Pi-(DPhi), DEta, PVz);
        hist->Fill(2*Pi-(DPhi),-DEta, PVz);
    }
}

void AliAnalysisTaskKaon2PC::FillDPhiHist(Double_t DPhi, TH2F* hist, Double_t fWeight=1){
    DPhi = fabs(DPhi);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,fWeight);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi,fWeight);
    }
    else {
        hist->Fill(2*Pi-(DPhi),fWeight);
    }
}

void AliAnalysisTaskKaon2PC::Fill2DHistMCTruth(Double_t DPhi, Double_t DEta, TH2F* hist){
    DPhi = fabs(DPhi);
    DEta = fabs(DEta);
    if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
    hist->Fill(DPhi,DEta);
    hist->Fill(DPhi,-DEta);
    if (DPhi < 0.5*Pi ) {
        hist->Fill(-DPhi, DEta);
        hist->Fill(-DPhi,-DEta); 
    }
    else {
        hist->Fill(2*Pi-(DPhi), DEta);
        hist->Fill(2*Pi-(DPhi),-DEta);
    }
}

//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::Terminate(Option_t *)
{
    if(fPoolMgr) delete fPoolMgr;
}
//_____________________________________________________________________________

