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



//#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TNtuple.h"
//#include <iostream>
#include "TString.h"
//#include "TPDGCode.h"
//#include "THnSparse.h"
//#include "TRandom3.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
//#include "AliMCEvent.h"
//#include "AliMCVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
//#include "AliExternalTrackParam.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODTrack.h"
#include "AliPID.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
//#include "AliAODRecoDecay.h"
#include "AliInputEventHandler.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisTaskKaon2PC.h"
#include "AliEventCuts.h"


//class AliAnalysisTaskMyTask;    // your analysis class
//using namespace std;              // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskKaon2PC)  // classimp: necessary for root

AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutput(0),
    fPrimaryVtx(0),
    fPIDResponse(0),
    fEventCuts(0),
    fHistCFPhi(0),
    fHistCFPhiCuts(0),
    fHistCFPhiLCuts(0),
    fHistK0Phi(0),
    fHistKPosPhi(0),
    fHistKNegPhi(0),
    fHistCFEta(0),
    fHistKPosEta(0),
    fHistKNegEta(0),
    fHistK0Eta(0),
    fHistV0M(0),
    fLpTCut(0.4),
    fUpTCut(0.8),
    fEtaCut(0.8),
    fSigCut(3.0),
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
    fSigNegv0Cut(3.0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________


AliAnalysisTaskKaon2PC::AliAnalysisTaskKaon2PC(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutput(0),
    fPrimaryVtx(0),
    fPIDResponse(0),
    fEventCuts(0),
    fHistCFPhi(0),
    fHistCFPhiCuts(0),
    fHistCFPhiLCuts(0),
    fHistK0Phi(0),
    fHistKPosPhi(0),
    fHistKNegPhi(0),
    fHistCFEta(0),
    fHistKPosEta(0),
    fHistKNegEta(0),
    fHistK0Eta(0),
    fHistV0M(0),
    fLpTCut(0.4),
    fUpTCut(0.8),
    fEtaCut(0.8),
    fSigCut(3.0),
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
    fSigNegv0Cut(3.0)
{
    //constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of analysis: here it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskKaon2PC::~AliAnalysisTaskKaon2PC()
{
    // destructor
    if(fOutput) {
        delete fOutput;                 //at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    
    fOutput = new TList();          // this is a list that contain all of your histograms
                                    // at the end of the analysis, the contents of this list are written to output file
    fOutput->SetOwner(kTRUE);

    fEventCuts = new AliEventCuts();
    fEventCuts->fUseITSTPCCluCorrelationCut = true;

    // create histograms
    
    fHistCFPhi = new TH2F("fHistCFPhi","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhi->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhi->GetYaxis()->SetTitle("Centrality");
    fHistCFPhi->SetOption("colz");
    
    fHistCFPhiCuts = new TH2F("fHistCFPhiCuts","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi with |#Delta#eta>0.5|",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhiCuts->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhiCuts->GetYaxis()->SetTitle("Centrality");
    fHistCFPhiCuts->SetOption("colz");
    
    fHistCFPhiLCuts = new TH2F("fHistCFPhiLCuts","Number of pairs of K^{+/-} and K^{0} Vs #Delta#Phi with |#Delta#eta<0.5|",32,-0.5*TMath::Pi(),1.5*TMath::Pi(),100,0,100);
    fHistCFPhiLCuts->GetXaxis()->SetTitle("#Delta#Phi (in radians)");
    fHistCFPhiLCuts->GetYaxis()->SetTitle("Centrality");
    fHistCFPhiLCuts->SetOption("colz");

    fHistKPosPhi = new TH2F("fHistKPosPhi", "Number of K^{+} Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistKPosPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistKPosPhi->GetYaxis()->SetTitle("Centrality");
    fHistKPosPhi->SetOption("colz");
    
    fHistKNegPhi = new TH2F("fHistKNegPhi", "Number of K^{-} Vs Track Phi; Centrality", 32,0,2*TMath::Pi(),100,0,100);
    fHistKNegPhi->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistKNegPhi->GetYaxis()->SetTitle("Centrality");
    fHistKNegPhi->SetOption("colz");

    fHistK0Phi = new TH2F("fHistK0Phi", "Number of K^{0}s Vs V0 Phi; Centrality",32,0,2*TMath::Pi(),100,0,100);
    fHistK0Phi->GetXaxis()->SetTitle("V0 Phi (in radians)");
    fHistK0Phi->GetYaxis()->SetTitle("Centrality");
    fHistK0Phi->SetOption("colz");
    
    fHistCFEta = new TH2F("fHistCFEta","Number of pairs of K^{+/-} and K^{0} Vs #Delta#eta",32,-1.6, 1.6, 100,0,100);
    fHistCFEta->GetXaxis()->SetTitle("#Delta#eta (in radians)");
    fHistCFEta->GetYaxis()->SetTitle("Centrality");
    fHistCFEta->SetOption("colz");
    
    fHistKPosEta = new TH2F("fHistKPosEta", "Number of K^{+} Vs Track Eta; Centrality", 16,-0.8, 0.8,100,0,100);
    fHistKPosEta->GetYaxis()->SetTitle("Centrality");
    fHistKPosEta->GetXaxis()->SetTitle("Track Eta (in radians)");
    fHistKPosEta->SetOption("colz");

    fHistKNegEta = new TH2F("fHistKNegEta", "Number of K^{-} Vs Track Eta; Centrality", 16,-0.8, 0.8,100,0,100);
    fHistKNegEta->GetYaxis()->SetTitle("Centrality");
    fHistKNegEta->GetXaxis()->SetTitle("Track Eta (in radians)");
    fHistKNegEta->SetOption("colz");

    fHistK0Eta = new TH2F("fHistK0Eta", "Number of K^{0}s Vs V0 Eta; Centrality",16,-0.8,0.8,100,0,100);
    fHistK0Eta->GetXaxis()->SetTitle("V0 Eta (in radians)");
    fHistK0Eta->GetYaxis()->SetTitle("Centrality");
    fHistK0Eta->SetOption("colz");

    fHistV0M = new TH1F("fHistV0M", "Centrality Distribution", 100,0,100);
    fHistV0M->GetXaxis()->SetTitle("Centrality");
    fHistV0M->GetYaxis()->SetTitle("Frequency");
    fHistV0M->SetOption("colz");
    
    /*
    fHistV0MEta = new TH1F("fHistV0MEta", "Eta Distribution", 100,0,100);
    fHistV0MEta->GetXaxis()->SetTitle("Centrality");
    fHistV0MEta->GetYaxis()->SetTitle("Frequency");
    fHistV0MEta->SetOption("colz");
    */
    
    fOutput->Add(fHistCFPhi);
    fOutput->Add(fHistCFPhiCuts);
    fOutput->Add(fHistCFPhiLCuts);
    fOutput->Add(fHistK0Phi);
    fOutput->Add(fHistKPosPhi);
    fOutput->Add(fHistKNegPhi);
    fOutput->Add(fHistCFEta);
    fOutput->Add(fHistKPosEta);
    fOutput->Add(fHistKNegEta);
    fOutput->Add(fHistK0Eta);
    fOutput->Add(fHistV0M);
    
    PostData(1, fOutput);               // postdata will notify the analysis manager of changes / updates to the
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}

void AliAnalysisTaskKaon2PC::SetTrackCuts(Double_t c1, Double_t c2, Double_t c3, Double_t c4){
    
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

Bool_t AliAnalysisTaskKaon2PC::AcceptTrack(const AliAODTrack *track) {
    if (!track->TestFilterBit(1)) return kFALSE;
    if (track->Pt() < fLpTCut || track->Pt() > fUpTCut) return kFALSE;
    if (fabs(track->Eta()) > fEtaCut) return kFALSE;
    Double_t nsigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    if (fabs(nsigmakaon) > fSigCut) return kFALSE;
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
    
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // have access to the current event.
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
    }
    if (!fEventCuts->AcceptEvent(fAOD)) return;
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();     //to get pid response object
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
            }
    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) return;
    
    Int_t nv0s(fAOD->GetNumberOfV0s());
    Int_t iTracks(fAOD->GetNumberOfTracks());
    
    Double_t vertex[3] = { -100.0, -100.0, -100.0 };
    const AliAODVertex *vertexAOD = fAOD->GetPrimaryVertex();
    vertexAOD->GetXYZ(vertex);
    
    Float_t  centralityV0M = -100.;
    AliMultSelection* MultSelection =0x0;
    MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");

    if (MultSelection) {
        centralityV0M = MultSelection->GetMultiplicityPercentile("V0M");
    } else {
        AliInfo("Didn't find MultSelection!");
        //centralityV0M = 999.;
    }
    
    Double_t PVz = fAOD->GetPrimaryVertex()->GetZ();
    if(fabs(PVz) > 10) return;
        
    if (fAOD) {
    fHistV0M->Fill(centralityV0M);
              }
    
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!AcceptTrack(track)) continue;
        Double_t trackPhi = track->Phi();
        Double_t trackEta = track->Eta();
        Int_t chtrack = track->Charge();
        if (chtrack > 0) {
            fHistKPosPhi->Fill(trackPhi,centralityV0M);
            fHistKPosEta->Fill(trackEta,centralityV0M);
        }
        if (chtrack < 0) {
            fHistKNegPhi->Fill(trackPhi,centralityV0M);
            fHistKNegEta->Fill(trackEta,centralityV0M);
        }
        for(Int_t j(0); j < nv0s; j++) {
            AliAODv0 *v0=fAOD->GetV0(j);
            if(!v0) continue;
            if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
            if(!AcceptV0(v0, vertex)) continue;
            Double_t V0Phi = v0->Phi();
            Double_t V0Eta = v0->Eta();
            fHistK0Phi->Fill(V0Phi, centralityV0M);
            fHistK0Eta->Fill(V0Eta, centralityV0M);
            Double_t deltaEta = fabs(trackEta- V0Eta);
            Double_t deltaPhi = trackPhi-V0Phi;
            if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
            if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
            fHistCFPhi->Fill(deltaPhi,centralityV0M);
            if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhi->Fill(-deltaPhi,centralityV0M);
            else fHistCFPhi->Fill(2*TMath::Pi()-deltaPhi,centralityV0M);
            fHistCFEta->Fill(deltaEta,centralityV0M);
            fHistCFEta->Fill(-deltaEta,centralityV0M);
            if(deltaEta > 0.5) {
                if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
                if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
                fHistCFPhiCuts->Fill(deltaPhi,centralityV0M);
                if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhiCuts->Fill(-deltaPhi,centralityV0M);
                else fHistCFPhiCuts->Fill(2*TMath::Pi()-deltaPhi,centralityV0M);
                }
            if(deltaEta < 0.5) {
                if (deltaPhi < 0) deltaPhi = V0Phi-trackPhi;
                if (deltaPhi > TMath::Pi()) deltaPhi = TMath::Pi()-(deltaPhi-TMath::Pi());
                fHistCFPhiLCuts->Fill(deltaPhi,centralityV0M);
                if (deltaPhi < (0.5*TMath::Pi())) fHistCFPhiLCuts->Fill(-deltaPhi,centralityV0M);
                else fHistCFPhiLCuts->Fill(2*TMath::Pi()-deltaPhi,centralityV0M);
                }
            }
        }
    
    PostData(1, fOutput);
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::Terminate(Option_t *)
{
    
}
//_____________________________________________________________________________


