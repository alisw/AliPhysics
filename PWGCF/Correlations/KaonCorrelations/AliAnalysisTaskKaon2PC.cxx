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
    fHistChKaons(0),
    fHistEtaCuts(0),
    fHistPosKaons(0),
    fHistNegKaons(0),
    fHistV0s(0),
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
    fHistChKaons(0),
    fHistEtaCuts(0),
    fHistPosKaons(0),
    fHistNegKaons(0),
    fHistV0s(0),
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
    
    fHistChKaons = new TH2F("fHistChKaons","Number of pairs of K^{+/-} and K^{0} Vs DeltaPhi and Centrality",30,-1*(TMath::Pi()/2),3*TMath::Pi()/2,100,0,100);
    fHistChKaons->GetXaxis()->SetTitle("Delta Phi (in radians)");
    fHistChKaons->GetYaxis()->SetTitle("Centrality");
    fHistChKaons->SetOption("colz");
    
    fHistEtaCuts = new TH2F("fHistEtaCuts","Number of pairs of K^{+/-} and K^{0} Vs DeltaPhi and Centrality with EtaCuts",30,0,2*TMath::Pi(),100,0,100);
    fHistEtaCuts->GetXaxis()->SetTitle("Delta Phi (in radians)");
    fHistEtaCuts->GetYaxis()->SetTitle("Centrality");
    fHistEtaCuts->SetOption("colz");
    
    fHistPosKaons = new TH2F("fHistPosKaons", "Number of K^{+} Vs Track Phi; Centrality", 30, 0, 2*TMath::Pi(),100,0,100);
    fHistPosKaons->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistPosKaons->GetYaxis()->SetTitle("Centrality");
    fHistPosKaons->SetOption("colz");
    
    fHistNegKaons = new TH2F("fHistNegKaons", "Number of K^{-} Vs Track Phi; Centrality", 30, 0, 2*TMath::Pi(),100,0,100);
    fHistNegKaons->GetXaxis()->SetTitle("Track Phi (in radians)");
    fHistNegKaons->GetYaxis()->SetTitle("Centrality");
    fHistNegKaons->SetOption("colz");
    
    fHistV0s = new TH2F("fHistV0s", "Number of V0s Vs V0 Phi; Centrality", 30, 0, 2*TMath::Pi(),100,0,100);
    fHistV0s->GetXaxis()->SetTitle("V0 Phi (in radians)");
    fHistV0s->GetYaxis()->SetTitle("Centrality");
    fHistV0s->SetOption("colz");

    fOutput->Add(fHistChKaons);
    fOutput->Add(fHistEtaCuts);
    fOutput->Add(fHistPosKaons);
    fOutput->Add(fHistNegKaons);
    fOutput->Add(fHistV0s);
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
        centralityV0M = 999.;
    }
    
    Double_t PVz = fAOD->GetPrimaryVertex()->GetZ();
    if(fabs(PVz) > 10 ) return;
    
    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!AcceptTrack(track)) continue;
        Double_t trackPhi = track->Phi();
        Int_t chtrack = track->Charge();
        if (chtrack > 0) fHistPosKaons->Fill(fabs(trackPhi),centralityV0M);
        fHistPosKaons-> Fill(fabs(fabs(trackPhi)-2*TMath::Pi()),centralityV0M);
        if (chtrack < 0) fHistNegKaons->Fill(fabs(trackPhi),centralityV0M);
        fHistNegKaons-> Fill(fabs(fabs(trackPhi)-2*TMath::Pi()),centralityV0M);
        
        for(Int_t j(0); j < nv0s; j++) {
            AliAODv0 *v0=fAOD->GetV0(j);
            if(!v0) continue;
            if(v0->MassK0Short() < 0.49 || v0->MassK0Short() > 0.51) continue;
            if(!AcceptV0(v0, vertex)) continue;
            Double_t V0Phi = v0->Phi();
            Double_t deltaPhi = fabs(trackPhi-V0Phi);
            Double_t deltaEta = (fabs(track->Eta()) - fabs(v0->Eta()));
            fHistChKaons->Fill(deltaPhi,centralityV0M);
            fHistChKaons->Fill(-deltaPhi,centralityV0M);
            if(deltaEta > 0.5) fHistEtaCuts->Fill(fabs(deltaPhi),centralityV0M);
            fHistEtaCuts-> Fill(fabs(fabs(deltaPhi)-2*TMath::Pi()),centralityV0M);
            fHistV0s->Fill(fabs(V0Phi),centralityV0M);
            fHistV0s->Fill(fabs(fabs(V0Phi)-2*TMath::Pi()),centralityV0M);
                                      }
                                    }
    PostData(1, fOutput);
}
//_____________________________________________________________________________
void AliAnalysisTaskKaon2PC::Terminate(Option_t *)
{
    
}
//_____________________________________________________________________________


