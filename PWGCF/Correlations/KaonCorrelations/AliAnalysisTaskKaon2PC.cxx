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
#include "TRandom3.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliMultiInputEventHandler.h"
#include "AliAnalysisAlien.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskKaon2PC.h"


class AliAnalysisTaskKaon2PC;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskKaon2PC) // classimp: necessary for root

Double_t Pi = TMath::Pi();
AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fPIDResponse(0), PVz(0), fLpTCut(0.2), fUpTCut(0.5), fEtaCut(0.8), fSigCut(2.0), 
    fDecayLv0Cut(8.05), fLpTv0Cut(0.2), fUpTv0Cut(0.5), fEtav0Cut(0.8), fDcaPosToPrimVtxv0Cut(0.1), 
    fDcaNegToPrimVtxv0Cut(0.1), fEtaPosv0Cut(0.8), fEtaNegv0Cut(0.8), fCosPACut(0.99), fSigPosv0Cut(2.0), 
    fSigNegv0Cut(2.0), fHistMK0(0), fHistMK0Cuts(0), fHistNV0(0), fHistPt(0), fVtx(0), fClusters(0), fPID(0), fPIDKa(0), fPIDKaon(0), fPIDK(0), fPIDKeCut(0), fPIDKpiCut(0), 
    fHistK0PhiEta(0), fHistK0Phi(0), fHistK0Eta(0), fHistChPhi(0), fHistChEta(0), fHistChRap(0), fHistPosPhi(0), fHistPosEta(0), 
    fHistPosPhiEta(0), fHistPosRap(0), fHistNegPhi(0), fHistNegEta(0), fHistNegPhiEta(0), fHistNegRap(0), fnsigmakaon(0), fNsigmaKaon(0),
    fNsigmaTOFK(0), fNsigmaTOFKaon(0),  fNsigmaTPCTOFK(0), fHistNEvents(0), fHistEta(0), fHistDEta(0), fHistPhi(0), fHistDPhi(0), fHistMult(0), 
    fHistCent(0), fHistSigCent(0), fHistInvCent(0), fHistCFPhi(0), fHistCFEta(0), fHistCF(0),  
    fHistKChCh(0), fHistKPosKPos(0), fHistKPosKNeg(0), fHistKNegKNeg(0), fHistK0K0(0), fHistPPionPhi(0), fHistNPionPhi(0), fEventCuts(0)


{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fPIDResponse(0), PVz(0), fLpTCut(0.2), fUpTCut(0.5), fEtaCut(0.8), fSigCut(2.0), 
    fDecayLv0Cut(8.05), fLpTv0Cut(0.2), fUpTv0Cut(0.5), fEtav0Cut(0.8), fDcaPosToPrimVtxv0Cut(0.1), 
    fDcaNegToPrimVtxv0Cut(0.1), fEtaPosv0Cut(0.8), fEtaNegv0Cut(0.8), fCosPACut(0.99), fSigPosv0Cut(2.0), 
    fSigNegv0Cut(2.0), fHistMK0(0), fHistMK0Cuts(0), fHistNV0(0), fHistPt(0), fVtx(0), fClusters(0), fPID(0), fPIDKa(0), fPIDKaon(0), fPIDK(0), fPIDKeCut(0), fPIDKpiCut(0), 
    fHistK0PhiEta(0), fHistK0Phi(0), fHistK0Eta(0), fHistChPhi(0), fHistChEta(0), fHistChRap(0), fHistPosPhi(0), fHistPosEta(0), 
    fHistPosPhiEta(0), fHistPosRap(0), fHistNegPhi(0), fHistNegEta(0), fHistNegPhiEta(0), fHistNegRap(0), fnsigmakaon(0), fNsigmaKaon(0),
    fNsigmaTOFK(0), fNsigmaTOFKaon(0),  fNsigmaTPCTOFK(0), fHistNEvents(0), fHistEta(0), fHistDEta(0), fHistPhi(0), fHistDPhi(0), fHistMult(0), 
    fHistCent(0), fHistSigCent(0), fHistInvCent(0), fHistCFPhi(0), fHistCFEta(0), fHistCF(0),  
    fHistKChCh(0), fHistKPosKPos(0), fHistKPosKNeg(0), fHistKNegKNeg(0), fHistK0K0(0), fHistPPionPhi(0), fHistNPionPhi(0), fEventCuts(0)


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
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file.


    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       
    
    // example of a histogram
    fHistMK0=new TH1F("fHistMK0", "Invariant Mass Distribution of Neutral Kaons", 100, 0.4, 0.6);
    fHistMK0Cuts=new TH1F("fHistMK0Cuts", "Invariant Mass Distribution of Neutral Kaons", 100, 0.4, 0.6);

    fHistNV0 = new TH1F("fHistNV0","Number of V0s",100, 0, 5000);

    fHistPt = new TH1F("fHistPt", "p_{T} distribution of all Charged Kaon Tracks", 100, 0, 1);
    fHistPt->SetOption("HIST E p");

    fVtx = new TH1F("fVtx", "PV_{z} distribution of Tracks", 100, -15, 15);

    fClusters = new TH1F("fClusters", "TPCClusters distribution", 500, 1, 170);

    fPID= new TH2F("fPID","Particle Identification",800,0.2,1.0,800,0.0,700.0);
    fPID->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPID->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPID->SetOption("colz");
    
    fPIDKa= new TH2F("fPIDKa","Kaons with |nSigma|<5",800,0.2,1.0,800,0.0,700.0);
    fPIDKa->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDKa->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKa->SetOption("colz");

    fPIDKaon= new TH2F("fPIDKaon","Kaons with |nSigma|<3",800,0.2,1.0,800,0.0,700.0);
    fPIDKaon->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDKaon->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKaon->SetOption("colz");

    fPIDK = new TH2F("fPIDK","Kaons with |nSigma|<2",800,0.2,1.0,800,0.0,700.0);
    fPIDK->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDK->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDK->SetOption("colz");

    fPIDKeCut = new TH2F("fPIDKeCut","Kaon PID with Electron Rejection",800,0.2,2.0,800,0.0,700.0);
    fPIDKeCut->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fPIDKeCut->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKeCut->SetOption("colz");

    fPIDKpiCut= new TH2F("fPIDKpiCut","Kaon PID with Pion Rejection",800,0.2,2.0,800,0.0,700.0);
    fPIDKpiCut->GetXaxis()->SetTitle("p_{} [GeV/c]");
    fPIDKpiCut->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKpiCut->SetOption("colz");
    
    fHistK0PhiEta = new TH3F("fHistK0PhiEta", "Number of K0s Vs Track Phi and Track Eta; Centrality", 60,0,2*Pi,60,-0.8, 0.8,100,0,100);
    fHistK0PhiEta->GetXaxis()->SetTitle("V0 Phi (in radians)");
    fHistK0PhiEta->SetOption("SURF1");

    fHistK0Phi = new TH2F("fHistK0Phi", "Number of K0s Vs Track Phi; Centrality", 16,0,2*Pi,100,0,100);
    fHistK0Phi->GetXaxis()->SetTitle("V0 Phi");
    fHistK0Phi->SetOption("colz");
    fHistK0Phi->GetYaxis()->SetRangeUser(0,6000);
    
    fHistK0Eta = new TH2F("fHistK0Eta", "Number of K0s Vs Track Eta; Centrality", 16,-0.8, 0.8,100,0,100);
    fHistK0Eta->GetXaxis()->SetTitle("V0 Eta");
    fHistK0Eta->SetOption("colz");

    fHistChPhi = new TH2F("fHistChPhi", "Number of charged particles Vs Track Phi; Centrality", 16,0,2*Pi,100,0,100);
    fHistChPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistChPhi->SetOption("SURF1");

    fHistChEta = new TH2F("fHistChEta", "Number of charged particles Vs Track Eta; Centrality", 16,-0.8, 0.8 ,100,0,100);
    fHistChEta->GetXaxis()->SetTitle("Track Eta");
    fHistChEta->SetOption("SURF1");

    fHistChRap = new TH2F("fHistChRap", "Number of charged particles Vs Rapidity; Centrality", 16,-0.8, 0.8 ,100,0,100);
    fHistChRap->GetXaxis()->SetTitle("Rapidity");
    fHistChRap->SetOption("colz");

    fHistPosPhi = new TH2F("fHistPosPhi","Number of +ve Kaons Vs Track Phi", 16,0,2*Pi,100,0,100);
    fHistPosPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistPosPhi->SetOption("colz");
    
    fHistPosEta = new TH2F("fHistPosEta", "", 16,-0.8, 0.8,100,0,100);
    fHistPosEta->GetXaxis()->SetTitle("Track Eta");
    fHistPosEta->SetOption("colz");

    fHistPosPhiEta = new TH3F("fHistPosPhiEta", "",120,0,2*Pi, 120,-0.8, 0.8,100,0,100);
    fHistPosPhiEta->GetXaxis()->SetTitle("Track Eta");
    fHistPosPhiEta->SetOption("SURF1");
    
    fHistPosRap = new TH2F("fHistPosRap", "", 16,-0.8, 0.8,100,0,100);
    fHistPosRap->GetXaxis()->SetTitle("Rapidity");
    fHistPosRap->SetOption("colz");
    
    fHistNegPhi = new TH2F("fHistNegPhi", "Number of -ve Kaons Vs Track Phi", 16,0,2*Pi,100,0,100);
    fHistNegPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistNegPhi->SetOption("colz");
    
    fHistNegEta = new TH2F("fHistNegEta", "", 16,-0.8, 0.8,100,0,100);
    fHistNegEta->GetXaxis()->SetTitle("Track Eta");
    fHistNegEta->SetOption("colz");

    fHistNegPhiEta = new TH3F("fHistNegPhiEta", "",120,0,2*Pi, 120,-0.8, 0.8,100,0,100);
    fHistNegPhiEta->GetXaxis()->SetTitle("Track Eta");
    fHistNegPhiEta->SetOption("SURF1");

    fHistNegRap = new TH2F("fHistNegRap", "", 16,-0.8, 0.8,100,0,100);
    fHistNegRap->GetXaxis()->SetTitle("Rapidity");
    fHistNegRap->SetOption("colz");

    fnsigmakaon = new TH2F("fnsigmakaon","NSigma of Kaon Vs Transverse Momentum",700,0.1,2,700,-20.0,20.0);
    fnsigmakaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fnsigmakaon->GetYaxis()->SetTitle("NSigma of Kaon");
    fnsigmakaon->SetOption("colz");

    fNsigmaKaon = new TH2F("fNsigmaKaon","|NSigma|<2 Vs Pt of Kaons",700,0.1,2,700,-20.0,20.0);
    fNsigmaKaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaKaon->GetYaxis()->SetTitle("NSigma of Kaon");
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

    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 1, 0, 1);
    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -5, 5);
    fHistDEta = new TH1F("fHistDEta", "fHistDEta", 100, -10, 10);
    fHistPhi = new TH1F("fHistPhi", "Phi Distribution", 100, 0, 7);
    fHistDPhi = new TH1F("fHistDPhi", "fHistDPhi", 100, 0, 10);
    fHistMult = new TH1F("fHistMult", "Number of tracks", 100, 0, 100);
    fHistCent = new TH1F("fHistCent", "CentV0M", 100, 0, 100);

    fHistSigCent = new TH1F("fHistSigCent","", 100,-10,10);
    //fHistSigCent->SetOption("colz");

    fHistInvCent = new TH2F("fHistInvCent", "", 100,0.4,0.6,100,0,100);
    fHistInvCent->SetOption("colz");

    fHistCFPhi = new TH2F("fHistCFPhi","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi",32,-0.5*Pi,1.5*Pi, 100,0,100);
    fHistCFPhi->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhi->GetYaxis()->SetTitle("Centrality in %");
    fHistCFPhi->SetOption("SURF1");
    
    fHistCFEta = new TH2F("fHistCFEta","Number of pairs of K^{+/-} and K^{0} Vs #Delta#eta",32,-1.6, 1.6, 100,0,100);
    fHistCFEta->GetXaxis()->SetTitle("#Delta#eta ");
    fHistCFEta->GetYaxis()->SetTitle("Centrality in %");
    fHistCFEta->SetOption("SURF1");

    fHistCF = new TH3F("fHistCF","Number of pairs of K^{+/-} and K^{0}",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistCF->GetXaxis()->SetTitle("#Delta#phi ");
    fHistCF->GetYaxis()->SetTitle("#Delta#eta");
    fHistCF->SetOption("SURF1");

    fHistKChCh = new TH3F("fHistKChCh","Kch-Kch Correlation",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistKChCh->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKChCh->GetYaxis()->SetTitle("#Delta#eta");
    fHistKChCh->SetOption("SURF1");

    fHistKPosKPos = new TH3F("fHistKPosKPos","K^{+}-K^{+} Correlation",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistKPosKPos->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKPos->GetYaxis()->SetTitle("#Delta#eta");
    fHistKPosKPos->SetOption("SURF1");

    fHistKPosKNeg = new TH3F("fHistKPosKNeg","K^{+}-K^{-} Correlation",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistKPosKNeg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKPosKNeg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKPosKNeg->SetOption("SURF1");

    fHistKNegKNeg = new TH3F("fHistKNegKNeg","K^{-}-K^{-} Correlation",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistKNegKNeg->GetXaxis()->SetTitle("#Delta#phi ");
    fHistKNegKNeg->GetYaxis()->SetTitle("#Delta#eta");
    fHistKNegKNeg->SetOption("SURF1");

    fHistK0K0 = new TH3F("fHistK0K0","K0-K0 Correlation",120,-0.5*Pi,1.5*Pi,120,-1.6, 1.6,100,0,100 );
    fHistK0K0->GetXaxis()->SetTitle("#Delta#phi ");
    fHistK0K0->GetYaxis()->SetTitle("#Delta#eta");
    fHistK0K0->SetOption("SURF1");

    fHistPPionPhi = new TH2F("fHistPPionPhi","Number of +ve piondaughters Vs Track Phi", 16,0,2*Pi,100,0,100);
    fHistNPionPhi = new TH2F("fHistNPionPhi","Number of -ve piondaughters Vs Track Phi", 16,0,2*Pi,100,0,100);


    fOutputList->Add(fHistMK0);
    fOutputList->Add(fHistMK0Cuts);
    fOutputList->Add(fHistNV0);
    fOutputList->Add(fHistPt);
    fOutputList->Add(fVtx);
    fOutputList->Add(fClusters);
    fOutputList->Add(fPID);
    fOutputList->Add(fPIDKa);
    fOutputList->Add(fPIDKaon);
    fOutputList->Add(fPIDK);
    fOutputList->Add(fPIDKeCut);
    fOutputList->Add(fPIDKpiCut);
    fOutputList->Add(fHistK0PhiEta);
    fOutputList->Add(fHistK0Phi);
    fOutputList->Add(fHistK0Eta);
    fOutputList->Add(fHistChPhi);
    fOutputList->Add(fHistChEta);
    fOutputList->Add(fHistChRap);
    fOutputList->Add(fHistPosPhi);
    fOutputList->Add(fHistPosEta);
    fOutputList->Add(fHistPosPhiEta);
    fOutputList->Add(fHistNegPhiEta);
    fOutputList->Add(fHistPosRap);
    fOutputList->Add(fHistNegPhi);
    fOutputList->Add(fHistNegEta);
    fOutputList->Add(fHistNegRap);
    fOutputList->Add(fnsigmakaon);
    fOutputList->Add(fNsigmaKaon);
    fOutputList->Add(fNsigmaTOFK);
    fOutputList->Add(fNsigmaTOFKaon);
    fOutputList->Add(fNsigmaTPCTOFK);
    fOutputList->Add(fHistNEvents);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistDEta);
    fOutputList->Add(fHistPhi);
    fOutputList->Add(fHistDPhi);
    fOutputList->Add(fHistMult);
    fOutputList->Add(fHistCent);
    fOutputList->Add(fHistSigCent);
    fOutputList->Add(fHistInvCent);
    fOutputList->Add(fHistCFPhi);
    fOutputList->Add(fHistCFEta);
    fOutputList->Add(fHistCF);
    fOutputList->Add(fHistKChCh);
    fOutputList->Add(fHistKPosKPos);
    fOutputList->Add(fHistKPosKNeg);
    fOutputList->Add(fHistKNegKNeg);
    fOutputList->Add(fHistK0K0);
    fOutputList->Add(fHistPPionPhi);
    fOutputList->Add(fHistNPionPhi);
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

Bool_t AliAnalysisTaskKaon2PC::AcceptTrack(const AliAODTrack *Trk) {
    if (!Trk->TestFilterBit(272)) return kFALSE;
    if (Trk->Charge() == 0) return kFALSE;         //excluding neutral particles
    if (Trk->Pt() <= fLpTCut || Trk->Pt() >= fUpTCut) return kFALSE; // pt cut
    if (fabs(Trk->Eta()) > fEtaCut) return kFALSE; // eta cut
    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kKaon);
    //cout << "nsigma kaon values are" << nSigmakaon << endl;
    Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kPion);
    //cout << "nsigma pion values are" << nSigmapion << endl;
    Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kElectron);
    Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kKaon);
    Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kElectron);
    Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(Trk, AliPID::kPion);
    if (fabs(nSigmakaon) > fSigCut) return kFALSE;
    if (fabs(nSigmaelectron) < 2.0) return kFALSE;  // excluding electrons via TPC
    if (fabs(nSigmapion) < 2.0) return kFALSE;      // excluding pions via TPC
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::AcceptPosTrack(const AliAODTrack *Trk) {
    if(Trk->Charge() < 0) return kFALSE;        // excluding negative tracks
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::AcceptNegTrack(const AliAODTrack *Trk) {
    if(Trk->Charge() > 0) return kFALSE;        // excluding positive tracks
    return kTRUE;
}

Bool_t AliAnalysisTaskKaon2PC::AcceptV0(const AliAODv0 *v0, Double_t *vertex) {

    Double_t length = v0->DecayLengthV0(vertex);
    if (length > fDecayLv0Cut) return kFALSE;
    Double_t pT = v0->Pt();
    if (pT < fLpTv0Cut || pT > fUpTv0Cut) return kFALSE;
    if (fabs(v0->Eta()) > fEtav0Cut) return kFALSE;
    Double_t dcaPosToPV = v0->DcaPosToPrimVertex();
    if (dcaPosToPV < fDcaPosToPrimVtxv0Cut) return kFALSE;
    Double_t dcaNegToPV = v0->DcaNegToPrimVertex();
    if (dcaNegToPV < fDcaNegToPrimVtxv0Cut) return kFALSE;
    Double_t etaPos = v0->PseudoRapPos();
    if (fabs(etaPos) > fEtaPosv0Cut) return kFALSE;
    Double_t etaNeg = v0->PseudoRapNeg();
    if (fabs(etaNeg) > fEtaNegv0Cut) return kFALSE;
    Double_t cosPA= v0->CosPointingAngle(vertex);
    if (cosPA < fCosPACut) return kFALSE;
    Double_t armpt = v0->PtArmV0();
    Double_t alpha = v0->AlphaV0();
    if (armpt < 0.2*fabs(alpha)) return kFALSE;
    AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0);
    AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1);
    // GetPosID()
    //if (!ptrack) continue;
    //if (!ntrack) continue;
    //cout << pTrack << endl;
    Double_t nSigmaPionPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
    if (fabs(nSigmaPionPos) > fSigPosv0Cut) return kFALSE;
    Double_t nSigmaPionNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
    if (fabs(nSigmaPionNeg) > fSigNegv0Cut) return kFALSE;
    Double_t pos_clusters = pTrack->GetTPCNcls();
    Double_t neg_clusters = nTrack->GetTPCNcls();
    if (pos_clusters < 70) return kFALSE;
    if (neg_clusters < 70) return kFALSE;
    //if (pTrack == track->TrackId() || nTrack == track->TrackId()) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserExec(Option_t *)
{
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();     //to get pid response object
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler) fPIDResponse = inputHandler->GetPIDResponse();
    }

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
    
    //fRunNumber = fAOD->GetRunNumber();
    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) return;
    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (!fEventCuts.AcceptEvent(fAOD)) return;
    //cout << "Number of tracks is"<< iTracks << endl;
    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    //making a cut in pvz -10 to 10cm
    //cout << "primary vertices are " << PrimaryVertex << endl;

    Double_t vertex[3] = { -100.0, -100.0, -100.0 };
    const AliAODVertex *vertexAOD = fAOD->GetPrimaryVertex();
    vertexAOD->GetXYZ(vertex);  //explaination??
    if(!vertexAOD) return;

    PVz = vertexAOD->GetZ();
    if(fabs(PVz)>10.0) return;
    fVtx->Fill(PVz);
    
    //Multiplicity selection
    //AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    //if(!MultSelection) return;
    
    //double CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    float CentV0M = fEventCuts.GetCentrality(); //centrality
    //const PrimaryVertex = fEventCuts.GetPrimaryVertex(); //primary vertex
    //if (CentV0M>10) return;
    //cout << "centrality are " << CentV0M << endl;
    
    //********************************* PID Loop ********************************************************
    
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Double_t nSigmaelectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        Double_t nSigmaTOFkaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        Double_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        Double_t nClustersTPC = track->GetTPCNcls();
        fHistSigCent->Fill(nSigmakaon);
        if(!track) continue;                            // if we failed, skip this track
        if(abs(track->Eta())>0.8) continue;               // eta cut
        if (!track->TestFilterBit(272)) continue;       // filterbit selection
        //cout << "nsigma kaon values are" << nSigmakaon << endl;
        //if (track->Pt() <= fLpTCut || track->Pt() >= fUpTCut) continue; // pt cut (avoided for viewing purposes)
        fClusters->Fill(nClustersTPC);
        fPID->Fill(track->Pt(),track->GetTPCsignal());
        Int_t chargetrack = track->Charge();
        if (fabs(nSigmakaon)<5.0) {fPIDKa->Fill(track->Pt(),track->GetTPCsignal());}
        if (fabs(nSigmakaon)<3.0) {fPIDKaon->Fill(track->Pt(),track->GetTPCsignal());}
        if (fabs(nSigmakaon)<2.0) {
            fPIDK->Fill(track->Pt(),track->GetTPCsignal());
            fNsigmaKaon->Fill(track->Pt(), nSigmakaon);}
        if (fabs(nSigmakaon) < 20.0 ) {fnsigmakaon->Fill(track->Pt(), nSigmakaon);}
        if ((nSigmapion > 3.0) && (nSigmakaon < 3.0) && (nSigmaelectron > 3.0) ) { fPIDKeCut->Fill(track->Pt(),track->GetTPCsignal());  }
        if ((nSigmapion > 3.0) && (nSigmakaon < 3.0)) { fPIDKpiCut->Fill(track->Pt(),track->GetTPCsignal()); }
        fNsigmaTOFK->Fill(track->Pt(), nSigmaTOFkaon);
        if (fabs(nSigmaTOFkaon)<3.0) {fNsigmaTOFKaon->Fill(track->Pt(), nSigmaTOFkaon);}
        if ((fabs(nSigmakaon)<3.0) && (fabs(nSigmaelectron) > 3.0 ) && (fabs(nSigmapion) > 3.0 ) &&
            (fabs(nSigmaTOFelectron) > 3.0 ) && (fabs(nSigmaTOFpion) > 3.0 ) && (fabs(nSigmaTOFkaon)<3.0)  )
            {fNsigmaTPCTOFK->Fill(track->Pt(), nSigmaTOFkaon);}
         
    }   //end of PID loop

    //********************************     V0 Loop  *****************************************************

    Int_t nv0s(fAOD->GetNumberOfV0s());                 //number of decay or v0 vertices in the event
    cout << "Number of v0 vertices are"<< nv0s << endl;
    fHistNV0->Fill(nv0s);
    for(Int_t i = 0; i < nv0s; i++)  {
        AliAODv0 *v0=fAOD->GetV0(i);                    // pointer to reconstructed v0
        if(!v0) {
        cout<<"No V0 "<<endl;
        continue;
        }
        AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0);
        AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1);
        fHistMK0->Fill(v0->MassK0Short());
        if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
        if(!AcceptV0(v0, vertex)) continue;
        fHistMK0Cuts->Fill(v0->MassK0Short());
        fHistInvCent->Fill(CentV0M, v0->MassK0Short());
        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        fHistK0PhiEta->Fill(V0Phi,V0Eta,CentV0M);                  // Yield of neutral kaons in Phi
        fHistK0Phi->Fill(V0Phi,CentV0M); 
        fHistK0Eta->Fill(V0Eta,CentV0M);                  // Yield of neutral kaons in Eta
        fHistPPionPhi->Fill(pTrack->Phi(),CentV0M);
        fHistNPionPhi->Fill(pTrack->Phi(),CentV0M);
    }

    //***************************** Charged-Charged Correlation Loop *****************************************************
    

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  // get a track (type AliAODTrack) from the event
        if(!track1) continue;                            // if we failed, skip this track
        if (!AcceptTrack(track1)) continue;
        Int_t chargetrack = track1->Charge();
        Double_t track1Phi = track1->Phi();
        Double_t track1Eta = track1->Eta();
        for(Int_t j(i+1); j < iTracks; j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));  // get a track (type AliAODTrack) from the event
            if(!track2) continue;                            // if we failed, skip this track
            if (!AcceptTrack(track2)) continue;
            Double_t track2Phi = track2->Phi();
            Double_t track2Eta = track2->Eta();
            Double_t DEta = fabs(track1Eta - track2Eta);
            Double_t DPhi = fabs(track1Phi - track2Phi);
            if (DPhi > Pi) DPhi = Pi-(DPhi-Pi);
            fHistKChCh->Fill(DPhi,DEta, CentV0M);
            fHistKChCh->Fill(DPhi,-DEta, CentV0M);
            if (DPhi < (0.5*Pi)) {
                fHistKChCh->Fill(-DPhi,DEta,CentV0M);
                fHistKChCh->Fill(-DPhi,-DEta,CentV0M); 
                    }
            else {
                fHistKChCh->Fill(2*Pi-(DPhi),DEta,CentV0M);
                fHistKChCh->Fill(2*Pi-(DPhi),-DEta,CentV0M);}                                     
         }     // end of track2 for loop
        
     }   // end of track1 for loop
    
    //*************************** KPos-KPos Correlation Loop ******************************************************

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  // get a track (type AliAODTrack) from the event
        if(!track1) continue;
        if (!AcceptTrack(track1)) continue;                            // if we failed, skip this track
        if (!AcceptPosTrack(track1)) continue;
        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();
        for(Int_t j(i+1); j < iTracks; j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));  // get a track (type AliAODTrack) from the event
            if(!track2) continue;
            if (!AcceptTrack(track2)) continue;                            // if we failed, skip this track
            if (!AcceptPosTrack(track2)) continue;
            Double_t track2PosPhi = track2->Phi();
            Double_t track2PosEta = track2->Eta();
            Double_t DEtaPP = fabs(track1PosEta - track2PosEta);
            Double_t DPhiPP = fabs(track1PosPhi - track2PosPhi);
            if (DPhiPP > Pi) DPhiPP = Pi-(DPhiPP-Pi);
            fHistKPosKPos->Fill(DPhiPP,DEtaPP,CentV0M);
            fHistKPosKPos->Fill(DPhiPP,-DEtaPP,CentV0M);
            if (DPhiPP < (0.5*Pi)) {
                fHistKPosKPos->Fill(-DPhiPP,DEtaPP,CentV0M);
                fHistKPosKPos->Fill(-DPhiPP,-DEtaPP,CentV0M); 
                }
            else {
                fHistKPosKPos->Fill(2*Pi-(DPhiPP),DEtaPP,CentV0M);
                fHistKPosKPos->Fill(2*Pi-(DPhiPP),-DEtaPP,CentV0M);}
        }
    }

    //*************************** KPos-KNeg Correlation Loop ******************************************************

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  // get a track (type AliAODTrack) from the event
        if(!track1) continue;                            // if we failed, skip this track
        if (!AcceptTrack(track1)) continue;
        if (!AcceptPosTrack(track1)) continue;
        Double_t track1PosPhi = track1->Phi();
        Double_t track1PosEta = track1->Eta();
        for(Int_t j(i+1); j < iTracks; j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));  // get a track (type AliAODTrack) from the event
            if(!track2) continue;                            // if we failed, skip this track
            if (!AcceptTrack(track2)) continue;
            if (!AcceptNegTrack(track2)) continue;
            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();
            Double_t DEtaPN = fabs(track1PosEta - track2NegEta);
            Double_t DPhiPN = fabs(track1PosPhi - track2NegPhi);
            if (DPhiPN > Pi) DPhiPN = Pi-(DPhiPN-Pi);
            fHistKPosKNeg->Fill(DPhiPN,DEtaPN,CentV0M);
            fHistKPosKNeg->Fill(DPhiPN,-DEtaPN,CentV0M);
            if (DPhiPN < (0.5*Pi)) {
                fHistKPosKNeg->Fill(-DPhiPN,DEtaPN,CentV0M);
                fHistKPosKNeg->Fill(-DPhiPN,-DEtaPN,CentV0M); 
                }
            else {
                fHistKPosKNeg->Fill(2*Pi-(DPhiPN),DEtaPN,CentV0M);
                fHistKPosKNeg->Fill(2*Pi-(DPhiPN),-DEtaPN,CentV0M);}
        }
    }

    //*************************** KNeg-KNeg Correlation Loop ******************************************************

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  // get a track (type AliAODTrack) from the event
        if(!track1) continue;                            // if we failed, skip this track
        if (!AcceptTrack(track1)) continue;
        if (!AcceptNegTrack(track1)) continue;
        Double_t track1NegPhi = track1->Phi();
        Double_t track1NegEta = track1->Eta();
        for(Int_t j(i+1); j < iTracks; j++) {
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));  // get a track (type AliAODTrack) from the event
            if(!track2) continue;                            // if we failed, skip this track
            if (!AcceptTrack(track2)) continue;
            if (!AcceptNegTrack(track2)) continue;
            Double_t track2NegPhi = track2->Phi();
            Double_t track2NegEta = track2->Eta();
            Double_t DEtaNN = fabs(track1NegEta - track2NegEta);
            Double_t DPhiNN = fabs(track1NegPhi - track2NegPhi); 
            if (DPhiNN > Pi) DPhiNN = Pi-(DPhiNN-Pi);
            fHistKNegKNeg->Fill(DPhiNN,DEtaNN,CentV0M);
            fHistKNegKNeg->Fill(DPhiNN,-DEtaNN,CentV0M);
            if (DPhiNN < (0.5*Pi)) {
                fHistKNegKNeg->Fill(-DPhiNN,DEtaNN,CentV0M);
                fHistKNegKNeg->Fill(-DPhiNN,-DEtaNN,CentV0M); 
                }
            else {
                fHistKNegKNeg->Fill(2*Pi-(DPhiNN),DEtaNN,CentV0M);
                fHistKNegKNeg->Fill(2*Pi-(DPhiNN),-DEtaNN,CentV0M);}                               
         }     // end of track2 for loop
        
    }   // end of track1 for loop
    //*****************************  Neutral-Neutral Correlation Loop ******************************************************
    
    for(Int_t i(0); i < nv0s; i++) {
        AliAODv0 *v01=fAOD->GetV0(i);
        if(!v01) continue;
        if(v01->MassK0Short() < 0.49 || v01->MassK0Short() > 0.51) continue;
        if(!AcceptV0(v01, vertex)) continue;
        for(Int_t j(i+1); j < nv0s; j++) {
            AliAODv0 *v02=fAOD->GetV0(j);
            if(!v02) continue;
            if(v02->MassK0Short() < 0.49 || v02->MassK0Short() > 0.51) continue;
            if(!AcceptV0(v02, vertex)) continue;
            Int_t PosIDTrack1 = v01->GetPosID();
            Int_t NegIDTrack1 = v01->GetNegID();
            Int_t PosIDTrack2 = v02->GetPosID();
            Int_t NegIDTrack2 = v02->GetNegID();
            if (PosIDTrack1 == PosIDTrack2) continue;
            if (NegIDTrack1 == NegIDTrack2) continue;
            Double_t V01Phi = v01->Phi();
            Double_t V01Eta = v01->Eta();
            Double_t V02Phi = v02->Phi();
            Double_t V02Eta = v02->Eta();
            Double_t DK0Phi = fabs(V01Phi-V02Phi);
            Double_t DK0Eta = fabs(V01Eta- V02Eta);
            if (DK0Phi > Pi) DK0Phi = Pi-(DK0Phi-Pi);
            fHistK0K0->Fill(DK0Phi,DK0Eta,CentV0M);
            fHistK0K0->Fill(DK0Phi,-DK0Eta,CentV0M);
            if (DK0Phi < (0.5*Pi)) {
                fHistK0K0->Fill(-DK0Phi,DK0Eta,CentV0M);
                fHistK0K0->Fill(-DK0Phi,-DK0Eta,CentV0M); 
                    }
            else {
                fHistK0K0->Fill(2*Pi-(DK0Phi),DK0Eta,CentV0M);
                fHistK0K0->Fill(2*Pi-(DK0Phi),-DK0Eta,CentV0M);}

        }
    }

    //******************************** Charged-Neutral Correlation Loop  *****************************************************

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));  // get a track (type AliAODTrack) from the event
        if(!track) continue;                            // if we failed, skip this track
        if (!AcceptTrack(track)) continue;
        Int_t chargetrack = track->Charge();
        fHistChPhi->Fill(track->Phi(),CentV0M);
        fHistChEta->Fill(track->Eta(),CentV0M);
        fHistChRap->Fill(track->Y(),CentV0M);
        Double_t trackPhi = track->Phi();
        Double_t trackEta = track->Eta();
        if (chargetrack > 0) {
            fHistPosPhi->Fill(track->Phi(),CentV0M);          //Yield vs Phi of +ve charged particles histo
            fHistPosEta->Fill(track->Eta(),CentV0M);          //Yield vs Eta of +ve charged particles histo
            fHistPosPhiEta->Fill(track->Phi(),track->Eta(),CentV0M);
            fHistPosRap->Fill(track->Y(0.493),CentV0M);
            }
        if (chargetrack < 0) {
            fHistNegPhi->Fill(track->Phi(),CentV0M);          //Yield vs Phi of -ve charged particles histo
            fHistNegEta->Fill(track->Eta(),CentV0M);          //Yield vs Eta of -ve charged particles histo
            fHistNegPhiEta->Fill(track->Phi(),track->Eta(),CentV0M);
            fHistNegRap->Fill(track->Y(0.493),CentV0M);
            }
        fHistPt->Fill(track->Pt());
        fHistPhi->Fill(track->Phi());
        fHistEta->Fill(track->Eta());
        for(Int_t j(0); j < nv0s; j++) {
            AliAODv0 *v0=fAOD->GetV0(j);
            if(!v0) continue;
            if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
            if(!AcceptV0(v0, vertex)) continue;
            Int_t KChID = track->GetID();
            Int_t PosIDTrack = v0->GetPosID();
            Int_t NegIDTrack = v0->GetNegID();
            if (KChID == PosIDTrack) continue;
            if (KChID == NegIDTrack) continue;
            Double_t V0Phi = v0->Phi();
            Double_t V0Eta = v0->Eta();
            Double_t deltaEta = fabs(trackEta-V0Eta);
            Double_t deltaPhi = fabs(trackPhi-V0Phi);
            fHistDEta->Fill(deltaEta);
            fHistDEta->Fill(-deltaEta);
            if (deltaPhi > Pi) deltaPhi = Pi-(deltaPhi-Pi);
            fHistDPhi->Fill(deltaPhi);
            fHistCFPhi->Fill(deltaPhi,CentV0M);
            fHistCF->Fill(deltaPhi,deltaEta,CentV0M);
            fHistCF->Fill(deltaPhi,-deltaEta,CentV0M);                           
            if (deltaPhi < 0.5*Pi) {
                fHistDPhi->Fill(-deltaPhi);
                fHistCFPhi->Fill(-deltaPhi,CentV0M);
                fHistCF->Fill(-deltaPhi,deltaEta,CentV0M);
                fHistCF->Fill(-deltaPhi,-deltaEta,CentV0M); 
                }
            else {
                fHistDPhi->Fill(2*Pi-deltaPhi);
                fHistCFPhi->Fill(2*Pi-deltaPhi,CentV0M);
                fHistCF->Fill(2*Pi-(deltaPhi),deltaEta,CentV0M);
                fHistCF->Fill(2*Pi-(deltaPhi),-deltaEta,CentV0M);
                }
            if (0 < deltaPhi < 0.5*Pi) {    
                fHistCFEta->Fill(deltaEta,CentV0M);
                fHistCFEta->Fill(-deltaEta,CentV0M);  // filling DeltaEta 
            }                                  
        }     // end of v0 track for loop
    }   // end of track for loop
    

fHistNEvents->Fill(0);
fHistMult->Fill(iTracks);
fHistCent->Fill(CentV0M);                          //Nevents vs centrality
PostData(1, fOutputList);
}

//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
