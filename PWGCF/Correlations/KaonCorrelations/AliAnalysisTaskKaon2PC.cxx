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
#include "AliVParticle.h"
#include "TLorentzVector.h"
#include "AliAnalysisTaskKaon2PC.h"


class AliAnalysisTaskKaon2PC;       // This analysis class
using namespace std; // std namespace: so you can do things like 'cout'
ClassImp(AliAnalysisTaskKaon2PC)    // classimp: necessary for root
const Double_t Pi = TMath::Pi();

AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC() : AliAnalysisTaskSE(),  
fAOD(0),
fmcEvent(0), 
fOutputList(0), 
fPIDResponse(0),
fRejectEventPileUp(kTRUE),
fAnalysisMC(kFALSE),
PVx(0), PVy(0), PVz(0),
//track cuts 
fLpTCut(0.4), 
fUpTCut(0.8), 
fEtaCut(0.8), 
fSigCut(2.0),
fBit(96),
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
fSigPosv0Cut(2.0), 
fSigNegv0Cut(2.0),
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
fHistNEvents(0),
fHistNV0(0),  
fHistEta(0), 
fHistDEta(0), 
fHistPhi(0), 
fHistDPhi(0), 
fHistMult(0), 
fHistCent(0),
fHistTPCTracksVsClusters(0),
//single particle histograms
fHistMK0(0), 
fHistMK0Cuts(0), 
fHistKChPt(0), 
fHistK0Pt(0), 
fHistKChPhi(0), 
fHistK0Phi(0), 
fHistKpPhi(0), 
fHistKnPhi(0),
fHistPPionPhi(0), 
fHistNPionPhi(0),
//single particle 2D histograms
f2DHistK0Phi(0),
f2DHistK0Eta(0),  
f2DHistChPhi(0), 
f2DHistChEta(0), 
f2DHistChRap(0),
f2DHistPosPhi(0),
f2DHistPosEta(0), 
f2DHistPosRap(0), 
f2DHistNegPhi(0), 
f2DHistNegEta(0), 
f2DHistNegRap(0),
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
fHistKChKCh(0),
fHistKPosKNeg(0), 
//mixing
fSelectedKCh(0),
fSelectedK0s(0),
fSelectedKpos(0),
fSelectedKneg(0),
fPoolMgr(0x0),
fPoolMaxNEvents(1000),
fPoolMinNTracks(5000),
fMinEventsToMix(10),
fNzVtxBins(10),
fNCentBins(15),
fKpKnCorr(kTRUE),
fK0KchCorr(kTRUE),
/*
fNOfSamples(1.0),
fSampleIndex(0.0),
*/
hPt(0),
hPt_kPos(0),
fHistCF_Bg(0),
fHistCF_KpKn_Bg(0),
//MC Truth
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fHistKpKnMC(0),
fHistK0KchMC(0),
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
fRejectEventPileUp(kTRUE),
fAnalysisMC(kFALSE),
PVx(0), PVy(0), PVz(0),
//track cuts 
fLpTCut(0.4), 
fUpTCut(0.8), 
fEtaCut(0.8), 
fSigCut(2.0),
fBit(96),
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
fSigPosv0Cut(2.0), 
fSigNegv0Cut(2.0),
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
fHistNEvents(0),
fHistNV0(0),  
fHistEta(0), 
fHistDEta(0), 
fHistPhi(0), 
fHistDPhi(0), 
fHistMult(0), 
fHistCent(0),
fHistTPCTracksVsClusters(0),
//single particle 1D histograms
fHistMK0(0), 
fHistMK0Cuts(0), 
fHistKChPt(0), 
fHistK0Pt(0), 
fHistKChPhi(0), 
fHistK0Phi(0), 
fHistKpPhi(0), 
fHistKnPhi(0),
fHistPPionPhi(0), 
fHistNPionPhi(0),
//single particle 2D histograms
f2DHistK0Phi(0),
f2DHistK0Eta(0),  
f2DHistChPhi(0), 
f2DHistChEta(0), 
f2DHistChRap(0),
f2DHistPosPhi(0),
f2DHistPosEta(0), 
f2DHistPosRap(0), 
f2DHistNegPhi(0), 
f2DHistNegEta(0), 
f2DHistNegRap(0),
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
fHistKChKCh(0),
fHistKPosKNeg(0), 
//mixing
fSelectedKCh(0),
fSelectedK0s(0),
fSelectedKpos(0),
fSelectedKneg(0),
fPoolMgr(0x0),
fPoolMaxNEvents(1000),
fPoolMinNTracks(5000),
fMinEventsToMix(10),
fNzVtxBins(10),
fNCentBins(15),
fKpKnCorr(kTRUE),
fK0KchCorr(kTRUE),
/*
fNOfSamples(1.0),
fSampleIndex(0.0),
*/
hPt(0),
hPt_kPos(0),
fHistCF_Bg(0),
fHistCF_KpKn_Bg(0),
//MC Truth
fMCK0(0),
fMCKpos(0),
fMCKneg(0),
fHistKpKnMC(0),
fHistK0KchMC(0),
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
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file.

    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0}; // 10 bins
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
    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 1, 0, 1);
    fHistNV0 = new TH1F("fHistNV0","Number of V0s",100, 0, 5000);
    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -5, 5);
    fHistDEta = new TH1F("fHistDEta", "fHistDEta", 100, -10, 10);
    fHistPhi = new TH1F("fHistPhi", "Phi Distribution", 100, 0, 7);
    fHistDPhi = new TH1F("fHistDPhi", "fHistDPhi", 100, 0, 10);
    fHistMult = new TH1F("fHistMult", "Number of tracks", 100, 0, 100);
    fHistCent = new TH1F("fHistCent", "CentV0M", 100, 0, 100);
    
    fHistTPCTracksVsClusters =  new TH2F("fHistTPCTracksVsClusters","fHistTPCTracksVsClusters",400,0,400,50000,0,50000);
    fHistTPCTracksVsClusters->Sumw2();
    fHistTPCTracksVsClusters->GetXaxis()->SetTitle("TPC tracks");
    fHistTPCTracksVsClusters->GetYaxis()->SetTitle("TPC Clusters");

    //single particle histograms
    fHistMK0=new TH1F("fHistMK0", "Invariant Mass Distribution of Neutral Kaons", 100, 0.4, 0.6);
    fHistMK0Cuts=new TH1F("fHistMK0Cuts", "Invariant Masss Distribution of Neutral Kaons After cuts", 100, 0.4, 0.6);

    fHistKChPt = new TH1F("fHistKChPt", "p_{T} distribution of all Charged Kaon Tracks", 100, 0, 2);
    fHistKChPt->SetOption("HIST E p");
    fHistK0Pt = new TH1F("fHistK0Pt", "p_{T} distribution of all Neutral Kaon Tracks", 100, 0, 2);
    fHistK0Pt->SetOption("HIST E p");

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
    
    f2DHistK0Phi = new TH2F("f2DHistK0Phi", "Number of K0s Vs Track Phi; Centrality", 16,0,2*Pi,100,0,100);
    f2DHistK0Phi->GetXaxis()->SetTitle("V0 Phi");
    f2DHistK0Phi->SetOption("colz");
    
    f2DHistK0Eta = new TH2F("f2DHistK0Eta", "Number of K0s Vs Track Eta; Centrality", 16,-0.8, 0.8,100,0,100);
    f2DHistK0Eta->GetXaxis()->SetTitle("V0 Eta");
    f2DHistK0Eta->SetOption("colz");

    f2DHistChPhi = new TH2F("f2DHistChPhi", "Number of charged particles Vs Track Phi; Centrality", 16,0,2*Pi,100,0,100);
    f2DHistChPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    f2DHistChPhi->SetOption("colz");

    f2DHistChEta = new TH2F("f2DHistChEta", "Number of charged particles Vs Track Eta; Centrality", 16,-0.8, 0.8 ,100,0,100);
    f2DHistChEta->GetXaxis()->SetTitle("Track Eta");
    f2DHistChEta->SetOption("colz");

    f2DHistChRap = new TH2F("f2DHistChRap", "Number of charged particles Vs Rapidity; Centrality", 16,-0.8, 0.8 ,100,0,100);
    f2DHistChRap->GetXaxis()->SetTitle("Rapidity");
    f2DHistChRap->SetOption("colz");

    f2DHistPosPhi = new TH2F("f2DHistPosPhi","Number of +ve Kaons Vs Track Phi", 16,0,2*Pi,100,0,100);
    f2DHistPosPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    f2DHistPosPhi->SetOption("colz");
    
    f2DHistPosEta = new TH2F("f2DHistPosEta", "", 16,-0.8, 0.8,100,0,100);
    f2DHistPosEta->GetXaxis()->SetTitle("Track Eta");
    f2DHistPosEta->SetOption("colz");
    
    f2DHistPosRap = new TH2F("f2DHistPosRap", "", 16,-0.8, 0.8,100,0,100);
    f2DHistPosRap->GetXaxis()->SetTitle("Rapidity");
    f2DHistPosRap->SetOption("colz");
    
    f2DHistNegPhi = new TH2F("f2DHistNegPhi", "Number of -ve Kaons Vs Track Phi", 16,0,2*Pi,100,0,100);
    f2DHistNegPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    f2DHistNegPhi->SetOption("colz");
    
    f2DHistNegEta = new TH2F("f2DHistNegEta", "", 16,-0.8, 0.8,100,0,100);
    f2DHistNegEta->GetXaxis()->SetTitle("Track Eta");
    f2DHistNegEta->SetOption("colz");

    f2DHistNegRap = new TH2F("f2DHistNegRap", "", 16,-0.8, 0.8,100,0,100);
    f2DHistNegRap->GetXaxis()->SetTitle("Rapidity");
    f2DHistNegRap->SetOption("colz");

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
    fHistCFPhi = new TH2F("fHistCFPhi","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi and #Delta#eta",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
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

    fHistCF = new TH2F("fHistCF","Number of pairs of K^{+/-} and K^{0}",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6 );
    fHistCF->GetXaxis()->SetTitle("#Delta#phi ");
    fHistCF->GetYaxis()->SetTitle("#Delta#eta");
    fHistCF->SetOption("SURF1");

    fHistKPosKNeg = new TH2F("fHistKPosKNeg","K^{+}-K^{-} Correlation",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6 );
    fHistKPosKNeg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKNeg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKPosKNeg->SetOption("SURF1");

    fHistKChKCh = new TH2F("fHistKChKCh","Kch-Kch Correlation",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKChKCh->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKChKCh->GetYaxis()->SetTitle("#Delta#eta");
    fHistKChKCh->SetOption("SURF1");

    //mixing histograms
    //const Int_t sizeOfSamples = (Int_t) fNOfSamples;

    fHistCF_Bg = new TH2F("fHistCF_Bg","Background for CF of K^{+/-} and K^{0}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_Bg->SetOption("SURF1");

    fHistCF_KpKn_Bg = new TH2F("fHistCF_KpKn_Bg","Background for CF of K^{+} and K^{-}",32,-0.5*Pi,1.5*Pi,32,-1.6,1.6);
    fHistCF_KpKn_Bg->SetOption("SURF1");

    hPt = new TH1F("hPt", "Track pT distribution", 100, 0, 2);
    hPt_kPos = new TH1F("hPt_kPos", "Track pT distribution", 100, 0, 2);

    // MC Truth histograms

    //+++++++++++++++++++++ MC ++++++++++++++++++++++++++

    Int_t bins[4] = {100,32,32,100};
    Double_t min[4] = {0.2,0,-0.8,0.4};
    Double_t max[4] = {1.0,2*Pi,0.8,0.6};
    
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

    fHistKpKnMC = new TH2F("fHistKpKnMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistKpKnMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKpKnMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistKpKnMC->SetOption("SURF1");

    fHistK0KchMC = new TH2F("fHistK0KchMC","K^{+}-K^{-} Correlation for MC Truth",32,-0.5*Pi,1.5*Pi,32,-1.6, 1.6);
    fHistK0KchMC->GetXaxis()->SetTitle("#Delta#phi ");
    fHistK0KchMC->GetYaxis()->SetTitle("#Delta#eta");
    fHistK0KchMC->SetOption("SURF1");

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
    fOutputList->Add(fHistNEvents);
    fOutputList->Add(fHistNV0);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistDEta);
    fOutputList->Add(fHistPhi);
    fOutputList->Add(fHistDPhi);
    fOutputList->Add(fHistMult);
    fOutputList->Add(fHistCent);
    fOutputList->Add(fHistTPCTracksVsClusters);

    fOutputList->Add(fHistMK0);
    fOutputList->Add(fHistMK0Cuts);
    fOutputList->Add(fHistKChPt);
    fOutputList->Add(fHistK0Pt);
    fOutputList->Add(fHistKChPhi);
    fOutputList->Add(fHistK0Phi);
    fOutputList->Add(fHistKpPhi);
    fOutputList->Add(fHistKnPhi);  
    fOutputList->Add(fHistPPionPhi);
    fOutputList->Add(fHistNPionPhi);

    fOutputList->Add(f2DHistK0Phi);
    fOutputList->Add(f2DHistK0Eta);
    fOutputList->Add(f2DHistChPhi);
    fOutputList->Add(f2DHistChEta);
    fOutputList->Add(f2DHistChRap);
    fOutputList->Add(f2DHistPosPhi);
    fOutputList->Add(f2DHistPosEta);
    fOutputList->Add(f2DHistPosRap);
    fOutputList->Add(f2DHistNegPhi);
    fOutputList->Add(f2DHistNegEta);
    fOutputList->Add(f2DHistNegRap);

    fOutputList->Add(fHistK0PhiEta);
    fOutputList->Add(fHistPosPhiEta);
    fOutputList->Add(fHistNegPhiEta);
    
    fOutputList->Add(fHistCFPhi);
    fOutputList->Add(fHistCFEta);
    fOutputList->Add(fHistKChKChPhi);
    fOutputList->Add(fHistKPosKNegPhi);
    fOutputList->Add(fHistCF);
    fOutputList->Add(fHistKChKCh);
    fOutputList->Add(fHistKPosKNeg);
    fOutputList->Add(hPt);
    fOutputList->Add(hPt_kPos);

    fOutputList->Add(fHistCF_Bg);
    fOutputList->Add(fHistCF_KpKn_Bg);

    fOutputList->Add(fMCK0);
    fOutputList->Add(fMCKpos);
    fOutputList->Add(fMCKneg);
    fOutputList->Add(fHistKpKnMC);
    fOutputList->Add(fHistK0KchMC);
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
    if (Trk->Pt() <= fLpTCut || Trk->Pt() >= fUpTCut) return kFALSE; // pt cut
    if (fabs(Trk->Eta()) > fEtaCut) return kFALSE; // eta cut
    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kKaon);
    Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kPion);
    Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kElectron);
    Double_t nSigmaproton =TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kProton)) ;
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kKaon);
    Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kElectron);
    Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kPion);
    if (fabs(nSigmaelectron) < 2.0) return kFALSE;  // excluding electrons via TPC
    if (fabs(nSigmapion) < 2.0) return kFALSE;      // excluding pions via TPC
    if (fabs(nSigmaproton) < 2.0) return kFALSE;      // excluding pions via TPC
    //if (fabs(nSigmakaon < 3.0) && (nSigmapion < 3.0)) return kFALSE;
    //if (fabs(nSigmakaon < 3.0) && (nSigmaproton < 3.0)) return kFALSE;
    //if (fabs(nSigmakaon < 3.0) && (nSigmaelectron < 3.0)) return kFALSE;
    if (fabs(nSigmakaon) > fSigCut) return kFALSE;
    //if (fabs(nSigmaTOFkaon) > 3.0) return kFALSE;

    return kTRUE;
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

// MC Truth Booleans

Bool_t AliAnalysisTaskKaon2PC::SelectK0TracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectK0 =  mcPartPdg==310&& (isPhysPrim); 
    if (SelectK0) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKPosTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKpos =  mcPartPdg==321&& (isPhysPrim); 
    if (SelectKpos) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKNegTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKneg =  mcPartPdg==-321&& (isPhysPrim); 
    if (SelectKneg) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::SelectKchTracksMC(AliMCParticle *mcTrack ) {
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();
    Bool_t SelectKch =  mcPartPdg==321&& mcPartPdg==-321&& (isPhysPrim); 
    if (SelectKch) return kFALSE;
    return kTRUE;
}


//_____________________________________________________________________________

void AliAnalysisTaskKaon2PC::RunData() {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;                                   
    
    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) return;

    //==================== Pile Up==========================

    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (fRejectEventPileUp){
        if (!fEventCuts.AcceptEvent(fAOD)) return;
    }

    //primary Vertex    
    const AliVVertex* primVertex = fEventCuts.GetPrimaryVertex(); 
    if (!primVertex) return;
    
    Double_t PVx = primVertex->GetX();
    Double_t PVy = primVertex->GetY();
    Double_t PVz = primVertex->GetZ();

    Double_t vertex[3] = { -100.0, -100.0, -100.0 }; //?
    primVertex->GetXYZ(vertex);  

    //zvertex cut
    if ( ( TMath::Abs(PVz) ) >= 10.0) return ;
    fVtx->Fill(PVz);

    //Multiplicity selection
    //AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    //if(!MultSelection) return;
    
    //centrality
    //double CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    Double_t CentV0M = fEventCuts.GetCentrality(); //centrality
    if ((CentV0M < fCentMin)||(CentV0M > fCentMax)) return;

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
        if (fabs(nSigmakaon)<3.0) {fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());}
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
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    //cout << "TOFsignal is" << TOFsignal << endl;
    fNsigmaTPCTOFK->Fill(trackPt, nSigmaTOFkaon);
    fHistTOFKch->Fill(trackPt, beta);
    fSelectedKCh->Add((AliAODTrack*)track);

    //fill single particle charged kaon histograms
    fHistKChPhi->Fill(trackPhi);
    fHistKChPt->Fill(trackPt);
    f2DHistChPhi->Fill(trackPhi,CentV0M);
    f2DHistChEta->Fill(trackEta,CentV0M);
    f2DHistChRap->Fill(track->Y(),CentV0M);

    if (chargetrack > 0) {
        fSelectedKpos->Add((AliAODTrack*)track);
        fHistKpPhi->Fill(trackPhi);
        f2DHistPosPhi->Fill(trackPhi,CentV0M);          
        f2DHistPosEta->Fill(trackEta,CentV0M);          
        fHistPosPhiEta->Fill(trackPhi,trackEta);
        f2DHistPosRap->Fill(track->Y(0.493),CentV0M);
        }
    if (chargetrack < 0) {
        fSelectedKneg->Add((AliAODTrack*)track);
        fHistKnPhi->Fill(trackPhi);
        f2DHistNegPhi->Fill(trackPhi,CentV0M);          
        f2DHistNegEta->Fill(trackEta,CentV0M);          
        fHistNegPhiEta->Fill(trackPhi,trackEta);
        f2DHistNegRap->Fill(track->Y(0.493),CentV0M);
        }
    fHistPhi->Fill(track->Phi());
    fHistEta->Fill(track->Eta());
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

    if(pTrack->Charge()==0 || nTrack->Charge()==0) continue;

    Double_t V0Phi = v0->Phi();
    Double_t V0Eta = v0->Eta();
    Double_t V0Pt = v0->Pt();

    if(!v0) continue;
    if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
    if(!AcceptV0(v0, vertex)) continue;
    fHistK0Pt->Fill(V0Pt);
    Double_t pT = v0->Pt();
    if (pT <= fLpTv0Cut || pT >= fUpTv0Cut) continue;
    fSelectedK0s->Add(v0);
    //fill single particle neutral kaon histograms
    fHistK0Phi->Fill(V0Phi);
    fHistK0PhiEta->Fill(V0Phi,V0Eta);          
    f2DHistK0Phi->Fill(V0Phi,CentV0M); 
    f2DHistK0Eta->Fill(V0Eta,CentV0M);                  
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

fHistNEvents->Fill(1);
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
            if(!KChAssoc) continue;

            Double_t DPhiMix = fabs(K0Trig->Phi() - KChAssoc->Phi());
            Double_t DEtaMix = fabs(K0Trig->Eta() - KChAssoc->Eta());

            if (DPhiMix > Pi) DPhiMix = Pi-(DPhiMix-Pi);

            fHistCF_Bg->Fill(DPhiMix, DEtaMix);
            fHistCF_Bg->Fill(DPhiMix, -DEtaMix);

            if (DPhiMix < 0.5*Pi ) {
            fHistCF_Bg->Fill(-DPhiMix, DEtaMix);
            fHistCF_Bg->Fill(-DPhiMix,-DEtaMix); 
            }
            else {
            fHistCF_Bg->Fill(2*Pi-(DPhiMix), DEtaMix);
            fHistCF_Bg->Fill(2*Pi-(DPhiMix),-DEtaMix);
            }

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
            //if( KaonAssoc->Charge() != KaonTrig->Charge() ) continue;

            Double_t DPhiMix = fabs(KaonTrig->Phi() - KaonAssoc->Phi());
            Double_t DEtaMix = fabs(KaonTrig->Eta() - KaonAssoc->Eta());

            if (DPhiMix > Pi) DPhiMix = Pi-(DPhiMix-Pi);

            fHistCF_KpKn_Bg->Fill(DPhiMix, DEtaMix);
            fHistCF_KpKn_Bg->Fill(DPhiMix, -DEtaMix);

            if (DPhiMix < 0.5*Pi ) {
            fHistCF_KpKn_Bg->Fill(-DPhiMix, DEtaMix);
            fHistCF_KpKn_Bg->Fill(-DPhiMix,-DEtaMix); 
            }
            else {
            fHistCF_KpKn_Bg->Fill(2*Pi-(DPhiMix), DEtaMix);
            fHistCF_KpKn_Bg->Fill(2*Pi-(DPhiMix),-DEtaMix);
            }

            } //end of k- loop end

        } // end of k+ trigger loop 

    } //end of mixing event loop 

}//end of pool

}

    TObjArray* tracksClone = (TObjArray*) fSelectedKCh->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);   


} //end of RunData function

//=================== MC Truth Correlation function ========

void AliAnalysisTaskKaon2PC::RunMC() {

Int_t nAcceptedParticles =0;
AliMCParticle *mcTrack = 0x0;
fmcEvent  = dynamic_cast<AliMCEvent*> (MCEvent());
if(!fmcEvent){
    Printf("No MC particle branch found");
    return;
    }

AliVVertex * mcVertex = (AliVVertex*)fmcEvent->GetPrimaryVertex();
fPV[2] = mcVertex->GetZ();
if (TMath::Abs(fPV[2])>=10) return;

Int_t nMCTracks = fmcEvent->GetNumberOfTracks(); // MC Truth, Total number of tracks per event
AliVTrack *genTrackMix = 0x0;

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

//cout << "number of accepted particles from MC tracks is"<< nAcceptedParticles << endl;
fHistGenMultiplicity->Fill(nAcceptedParticles);

AliMCParticle *mcMotherParticle = 0x0;
AliMCParticle* daughter0 = 0x0;
AliMCParticle* daughter1 = 0x0;
Bool_t SelectK0;
Bool_t SelectKpos;
Bool_t SelectKneg;

for (Int_t i = 0; i < nMCTracks; i++){
	mcTrack = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!mcTrack) continue;
        
    Int_t mcPartPdg = mcTrack->PdgCode();
    Bool_t isPhysPrim = mcTrack->IsPhysicalPrimary();

    SelectK0 =  mcPartPdg==310&& (isPhysPrim);  // (8010 total. When we include only primary, 7546)-310 for neutral
    SelectKpos = mcPartPdg==321&& (isPhysPrim); // 321 is the code for positive kaons
    SelectKneg = mcPartPdg==-321&& (isPhysPrim); // 321 is the code for positive kaons

    Double_t TrackPt = mcTrack->Pt();
    Double_t TrackPhi = mcTrack->Phi();
    Double_t TrackEta = mcTrack->Eta();
    Double_t TrackMass = mcTrack->M();
    
    Double_t KaonVariables[4]= {TrackPt, TrackPhi, TrackEta, TrackMass};
    if(SelectK0) fMCK0->Fill(KaonVariables);
    if(SelectKpos) fMCKpos->Fill(KaonVariables);
    if(SelectKneg) fMCKneg->Fill(KaonVariables);
   
	Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
	Bool_t TrCharge = (mcTrack->Charge())!=0;
    Short_t cha;
    if (mcTrack->Charge()>0) cha=1.;
    else if (mcTrack->Charge()<0) cha= -1.;
    else cha =0;

}


for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack1 = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!SelectKPosTracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t trackPseudorap = mcTrack1->Eta();
    if(! (mcTrack1->IsPhysicalPrimary()&&mcTrack1->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < nMCTracks; j++){
        AliMCParticle *mcTrack2 = (AliMCParticle*)fmcEvent->GetTrack(j);
        if (!SelectKNegTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t trackPseudorap = mcTrack2->Eta();
        if(! (mcTrack2->IsPhysicalPrimary()&&mcTrack2->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistKpKnMC);
    }

}

for (Int_t i = 0; i < nMCTracks; i++){
    AliMCParticle *mcTrack1 = (AliMCParticle*)fmcEvent->GetTrack(i);
    if (!SelectK0TracksMC(mcTrack1)) continue;
    if (!mcTrack1) continue;
    Double_t trackPseudorap = mcTrack1->Eta();
    if(! (mcTrack1->IsPhysicalPrimary()&&mcTrack1->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
    Double_t phi1 = mcTrack1->Phi();
    Double_t eta1 = mcTrack1->Eta();
    
    for (Int_t j = i+1; j < nMCTracks; j++){
        AliMCParticle *mcTrack2 = (AliMCParticle*)fmcEvent->GetTrack(j);
        if (!SelectKchTracksMC(mcTrack2)) continue;
        if (!mcTrack2) continue;
        Double_t trackPseudorap = mcTrack2->Eta();
        if(! (mcTrack2->IsPhysicalPrimary()&&mcTrack2->Charge()!=0&&((trackPseudorap>-0.8&&trackPseudorap<0.8)))) continue;
        Double_t phi2 = mcTrack2->Phi();
        Double_t eta2 = mcTrack2->Eta();
        Double_t DEta = fabs(eta1 - eta2);
        Double_t DPhi = fabs(phi1 - phi2);
        Fill2DHistMCTruth(DPhi,DEta,fHistK0KchMC);
    }

}

}

//_____________________________________________________________________________

void AliAnalysisTaskKaon2PC::UserExec(Option_t *)
{
    if (!fAnalysisMC) { RunData(); }

    // deleting TObjArrays
    //fSelectedK0s->Clear();
    //fSelectedKCh->Clear();

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

