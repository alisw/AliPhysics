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


#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TMarker.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskKaon2PC.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliMultiInputEventHandler.h"
#include "AliAnalysisAlien.h"
#include "AliEventCuts.h"

class AliAnalysisTaskKaon2PC;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskKaon2PC) // classimp: necessary for root

AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC() : AliAnalysisTaskSE(),
    fAOD(0), PVz(0), fOutputList(0), fPIDResponse(0), fLpTCut(0.4), fUpTCut(0.8), fEtaCut(0.8), fSigCut(2.0), fDecayLv0Cut(8.05), fLpTv0Cut(0.2), fUpTv0Cut(1.0), fEtav0Cut(0.8), fDcaPosToPrimVtxv0Cut(0.1), fDcaNegToPrimVtxv0Cut(0.1), fEtaPosv0Cut(0.8), fEtaNegv0Cut(0.8), fCosPACut(0.99), fSigPosv0Cut(1.0), fSigNegv0Cut(1.0), fHistMK0(0), fHistPt(0),  fPID(0), fPIDKa(0), fPIDKaon(0), fPIDK(0), fHistK0Phi(0), fHistK0Eta(0), fHistChPhi(0), fHistChEta(0), fHistChRap(0), fHistPosPhi(0), fHistPosEta(0), fHistPosRap(0), fHistNegPhi(0), fHistNegEta(0),  fHistNegRap(0), fnsigmakaon(0), fNsigmaKaon(0), fHistNEvents(0), fHistEta(0), fHistDEta(0), fHistPhi(0), fHistDPhi(0), fHistMult(0), fHistCent(0), fHistCFPhi(0), fHistCFPhiCuts(0), fHistCFPhiLCuts(0), fHistCFEta(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), PVz(0), fOutputList(0), fPIDResponse(0), fLpTCut(0.4), fUpTCut(0.8), fEtaCut(0.8), fSigCut(2.0), fDecayLv0Cut(8.05), fLpTv0Cut(0.2), fUpTv0Cut(1.0), fEtav0Cut(0.8), fDcaPosToPrimVtxv0Cut(0.1), fDcaNegToPrimVtxv0Cut(0.1), fEtaPosv0Cut(0.8), fEtaNegv0Cut(0.8), fCosPACut(0.99), fSigPosv0Cut(1.0), fSigNegv0Cut(1.0), fHistMK0(0), fHistPt(0),  fPID(0), fPIDKa(0), fPIDKaon(0), fPIDK(0), fHistK0Phi(0), fHistK0Eta(0), fHistChPhi(0), fHistChEta(0), fHistChRap(0), fHistPosPhi(0), fHistPosEta(0), fHistPosRap(0), fHistNegPhi(0), fHistNegEta(0),  fHistNegRap(0), fnsigmakaon(0), fNsigmaKaon(0), fHistNEvents(0), fHistEta(0), fHistDEta(0), fHistPhi(0), fHistDPhi(0), fHistMult(0), fHistCent(0), fHistCFPhi(0), fHistCFPhiCuts(0), fHistCFPhiLCuts(0), fHistCFEta(0)
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
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       
    
    // example of a histogram
    
    
    
    fHistMK0=new TH1F("fHistMK0", "Invariant Mass Distribution of Neutral Kaons", 100, 0.4, 0.6);
    
    fHistPt = new TH1F("fHistPt", "Transverse Momentum distribution", 100, 0, 10);
    
    fPID= new TH2F("fPID","Particle Identification",800,0.2,1.0,800,0.0,700.0);
    fPID->GetXaxis()->SetTitle("p [GeV/c]");
    fPID->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPID->SetOption("colz");
    
    fPIDKa= new TH2F("fPIDKa","Kaons with |nSigma|<5",800,0.2,1.0,800,0.0,700.0);
    fPIDKa->GetXaxis()->SetTitle("p [GeV/c]");
    fPIDKa->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKa->SetOption("colz");

    fPIDKaon= new TH2F("fPIDKaon","Kaons with |nSigma|<3",800,0.2,1.0,800,0.0,700.0);
    fPIDKaon->GetXaxis()->SetTitle("p [GeV/c]");
    fPIDKaon->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDKaon->SetOption("colz");

    fPIDK = new TH2F("fPIDK","Kaons with |nSigma|<2",800,0.2,1.0,800,0.0,700.0);
    fPIDK->GetXaxis()->SetTitle("p [GeV/c]");
    fPIDK->GetYaxis()->SetTitle("TPC dE/dx [arb. units]");
    fPIDK->SetOption("colz");
    
    fHistK0Phi = new TH2F("fHistK0Phi", "Number of K0s Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistK0Phi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistK0Phi->SetOption("colz");
    
    fHistK0Eta = new TH2F("fHistK0Eta", "Number of K0s Vs Track Eta; Centrality", 17,-0.8, 0.8,100,0,100);
    fHistK0Eta->GetXaxis()->SetTitle("Track Eta");
    fHistK0Eta->SetOption("colz");

    fHistChPhi = new TH2F("fHistChPhi", "Number of charged particles Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistChPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistChPhi->SetOption("colz");

    fHistChEta = new TH2F("fHistChEta", "Number of charged particles Vs Track Eta; Centrality", 17,-0.8, 0.8 ,100,0,100);
    fHistChEta->GetXaxis()->SetTitle("Track Eta");
    fHistChEta->SetOption("colz");

    fHistChRap = new TH2F("fHistChRap", "Number of charged particles Vs Rapidity; Centrality", 17,-0.8, 0.8 ,100,0,100);
    fHistChRap->GetXaxis()->SetTitle("Rapidity");
    fHistChRap->SetOption("colz");

    fHistPosPhi = new TH2F("fHistPosPhi", "Number of +ve Kaons Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistPosPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistPosPhi->SetOption("colz");
    
    fHistPosEta = new TH2F("fHistPosEta", "", 17,-0.8, 0.8,100,0,100);
    fHistPosEta->GetXaxis()->SetTitle("Track Eta");
    fHistPosEta->SetOption("colz");
    
    fHistPosRap = new TH2F("fHistPosRap", "", 17,-0.8, 0.8,100,0,100);
    fHistPosRap->GetXaxis()->SetTitle("Rapidity");
    fHistPosRap->SetOption("colz");
    
    fHistNegPhi = new TH2F("fHistNegPhi", "Number of -ve Kaons Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistNegPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistNegPhi->SetOption("colz");

    fHistNegEta = new TH2F("fHistNegEta", "", 17,-0.8, 0.8,100,0,100);
    fHistNegEta->GetXaxis()->SetTitle("Track Eta");
    fHistNegEta->SetOption("colz");

    fHistNegRap = new TH2F("fHistNegRap", "", 17,-0.8, 0.8,100,0,100);
    fHistNegRap->GetXaxis()->SetTitle("Rapidity");
    fHistNegRap->SetOption("colz");

    fNsigmaKaon = new TH2F("fNsigmaKaon","|NSigma|<2 Vs Pt of Kaons",700,0.1,2,700,-20.0,20.0);
    fNsigmaKaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fNsigmaKaon->GetYaxis()->SetTitle("NSigma of Kaon");
    fNsigmaKaon->SetOption("colz");

    fnsigmakaon = new TH2F("fnsigmakaon","NSigma of Kaon Vs Transverse Momentum",700,0.1,2,700,-20.0,20.0);
    fnsigmakaon->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    fnsigmakaon->GetYaxis()->SetTitle("NSigma of Kaon");
    fnsigmakaon->SetOption("colz");

    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 1, 0, 1);
    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -10, 10);
    fHistDEta = new TH1F("fHistDEta", "fHistDEta", 100, -10, 10);
    fHistPhi = new TH1F("fHistPhi", "Phi Distribution", 100, 0, 10);
    fHistDPhi = new TH1F("fHistDPhi", "fHistDPhi", 100, 0, 10);
    fHistMult = new TH1F("fHistMult", "Number of tracks", 100, 0, 100);
    fHistCent = new TH1F("fHistCent", "CentV0M", 100, 0, 100);
    
    fHistCFPhi = new TH2F("fHistCFPhi","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhi->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhi->GetYaxis()->SetTitle("Centrality in %");
    fHistCFPhi->SetOption("colz");
        
    fHistCFPhiCuts = new TH2F("fHistCFPhiCuts","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi with |#Delta#eta>0.5|",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhiCuts->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhiCuts->GetYaxis()->SetTitle("Centrality in %");
    fHistCFPhiCuts->SetOption("colz");
        
    fHistCFPhiLCuts = new TH2F("fHistCFPhiLCuts","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi with |#Delta#eta<0.5|",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhiLCuts->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhiLCuts->GetYaxis()->SetTitle("Centrality in %");
    fHistCFPhiLCuts->SetOption("colz");
    
    fHistCFEta = new TH2F("fHistCFEta","Number of pairs of K^{+/-} and K^{0} Vs #Delta#eta",32,-1.6, 1.6, 100,0,100);
    fHistCFEta->GetXaxis()->SetTitle("#Delta#eta ");
    fHistCFEta->GetYaxis()->SetTitle("Centrality in %");
    fHistCFEta->SetOption("colz");

    fOutputList->Add(fHistMK0);
    fOutputList->Add(fHistPt);
    fOutputList->Add(fPID);
    fOutputList->Add(fPIDKa);
    fOutputList->Add(fPIDKaon);
    fOutputList->Add(fPIDK);
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
    fOutputList->Add(fnsigmakaon);
    fOutputList->Add(fNsigmaKaon);
    fOutputList->Add(fHistNEvents);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistDEta);
    fOutputList->Add(fHistPhi);
    fOutputList->Add(fHistDPhi);
    fOutputList->Add(fHistMult);
    fOutputList->Add(fHistCent);
    fOutputList->Add(fHistCFPhi);
    fOutputList->Add(fHistCFPhiCuts);
    fOutputList->Add(fHistCFPhiLCuts);
    fOutputList->Add(fHistCFEta);
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
    if (!Trk->TestFilterBit(768)) return kFALSE;
    if(Trk->Charge() == 0) return kFALSE;         //excluding neutral particles
    if (Trk->Pt() <= fLpTCut || Trk->Pt() >= fUpTCut) return kFALSE; // pt cut
    if (fabs(Trk->Eta()) > fEtaCut) return kFALSE; // eta cut
    Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kKaon);
    Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kPion);
    Double_t nSigmaproton = fPIDResponse->NumberOfSigmasTPC(Trk, AliPID::kProton);
    if (fabs(nSigmakaon) > fSigCut) return kFALSE;
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
    Double_t nSigmaPionPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
    if (fabs(nSigmaPionPos) > fSigPosv0Cut) return kFALSE;
    Double_t nSigmaPionNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
    if (fabs(nSigmaPionNeg) > fSigNegv0Cut) return kFALSE;
    return kTRUE;
}




//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    //cout << "Number of tracks is"<< iTracks << endl;
    //fRunNumber = fAOD->GetRunNumber();

    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (!fEventCuts.AcceptEvent(fAOD)) return;

    //making a cut in pvz -10 to 10cm
    const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
    if(!PrimaryVertex) return;
    //cout << "primary vertices are " << PrimaryVertex << endl;
    PVz = PrimaryVertex->GetZ();
    if(fabs(PVz)>10.0) return;

    Double_t vertex[3] = { -100.0, -100.0, -100.0 };
    const AliAODVertex *vertexAOD = fAOD->GetPrimaryVertex();
    vertexAOD->GetXYZ(vertex);  //explaination??

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();     //to get pid response object
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
            }
    
    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) return;

    //Multiplicity selection
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    double CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality

    //Loop for PID with |eta|<0.8
    for(Int_t i(0); i < iTracks; i++) {                 
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        if(!track) continue;                            // if we failed, skip this track
        if(abs(track->Eta())>0.8) continue;               // eta cut
        if (!track->TestFilterBit(768)) continue;       // filterbit selection
        fPID->Fill(track->P(),track->GetTPCsignal());
        Int_t chargetrack = track->Charge();
        if (fabs(nSigmakaon)<5.0) {fPIDKa->Fill(track->P(),track->GetTPCsignal());}
        if (fabs(nSigmakaon)<3.0) {fPIDKaon->Fill(track->P(),track->GetTPCsignal());}
        if (fabs(nSigmakaon)<2.0) {
            fPIDK->Fill(track->P(),track->GetTPCsignal());
            fNsigmaKaon->Fill(track->Pt(), nSigmakaon);}
        if (fabs(nSigmakaon) < 20.0 ) {fnsigmakaon->Fill(track->Pt(), nSigmakaon);} 
        }   //end of PID loop 

    Int_t nv0s(fAOD->GetNumberOfV0s());                 //number of decay or v0 vertices in the event
    for(Int_t i = 0; i < nv0s; i++)  {
	    AliAODv0 *v0=fAOD->GetV0(i);                    // pointer to reconstructed v0          
	    if(!v0) { 
	    cout<<"No V0 "<<endl;
	    continue; 
	    }
        fHistMK0->Fill(v0->MassK0Short());
        if(!AcceptV0(v0, vertex)) continue;
        if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
        Double_t V0Phi = v0->Phi();
        Double_t V0Eta = v0->Eta();
        fHistK0Phi->Fill(V0Phi, CentV0M);                  // Yield of neutral kaons in Phi
        fHistK0Eta->Fill(V0Eta, CentV0M);                  // Yield of neutral kaons in Eta                 
    }
    

    for(Int_t i(0); i < iTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        Double_t nSigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Double_t nSigmapion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
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
            fHistPosRap->Fill(track->Y(0.493),CentV0M);
            }
        if (chargetrack <0) {
            fHistNegPhi->Fill(track->Phi(),CentV0M);          //Yield vs Phi of -ve charged particles histo
            fHistNegEta->Fill(track->Eta(),CentV0M);          //Yield vs Eta of -ve charged particles histo
            fHistNegRap->Fill(track->Y(0.493),CentV0M);
            }
            
        fHistPt->Fill(track->Pt());
        fHistEta->Fill(track->Eta());
        fHistPhi->Fill(track->Phi());
        for(Int_t j(0); j < nv0s; j++) {
                    AliAODv0 *v0=fAOD->GetV0(j);
                    if(!v0) continue;
                    if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
                    if(!AcceptV0(v0, vertex)) continue;
                    Double_t V0Phi = v0->Phi();
                    Double_t V0Eta = v0->Eta();
                    Double_t deltaEta = fabs(trackEta- V0Eta);
                    Double_t deltaPhi = trackPhi-V0Phi;
                    fHistDEta->Fill(deltaEta);
                    fHistDPhi->Fill(deltaPhi);
                    if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
                    if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
                    fHistCFPhi->Fill(deltaPhi,CentV0M);                           // filling deltaphi
                    if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhi->Fill(-deltaPhi,CentV0M);
                    else fHistCFPhi->Fill(2*TMath::Pi()-deltaPhi,CentV0M);
                    fHistCFEta->Fill(deltaEta,CentV0M);
                    fHistCFEta->Fill(-deltaEta,CentV0M);                           // filling deltaeta
                    if(deltaEta > 0.5) {
                        if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
                        if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
                        fHistCFPhiCuts->Fill(deltaPhi,CentV0M);
                        if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhiCuts->Fill(-deltaPhi,CentV0M);
                        else fHistCFPhiCuts->Fill(2*TMath::Pi()-deltaPhi,CentV0M);
                        }
                    if(deltaEta < 0.5) {
                        if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
                        if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
                        fHistCFPhiLCuts->Fill(deltaPhi,CentV0M);
                        if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhiLCuts->Fill(-deltaPhi,CentV0M);
                        else fHistCFPhiLCuts->Fill(2*TMath::Pi()-deltaPhi,CentV0M);
                        }
                    }                                   // end of v0 track for loop
        
    }
                                                       // end of track for loop
    fHistNEvents->Fill(0);
    fHistMult->Fill(iTracks);
    fHistCent->Fill(CentV0M);                          //nevents vs centrality
    PostData(1, fOutputList);                           
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
