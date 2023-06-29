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

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliSPDppAnalysisTaskSim.h"
#include "AliAODMCParticle.h"

#include "TH3D.h"
#include "TTree.h"
#include "TString.h"

#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"
#include <AliHeader.h>
#include "TCanvas.h"

class AliSPDppAnalysisTaskSim;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliSPDppAnalysisTaskSim) // classimp: necessary for root

AliSPDppAnalysisTaskSim::AliSPDppAnalysisTaskSim() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), MultDist05(0), fMCEvent(0), responseMatrix(0), RecMultDist05(0), fEventCuts(0), fAODV0(0), eventcount1(0), fUseINT1(true), MultDist05inelgr0(0), responseMatrixinelgr0(0), RecMultDist05inelgr0(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliSPDppAnalysisTaskSim::AliSPDppAnalysisTaskSim(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), MultDist05(0), fMCEvent(0), responseMatrix(0), RecMultDist05(0), fEventCuts(0), fAODV0(0), eventcount1(0), fUseINT1(true), MultDist05inelgr0(0), responseMatrixinelgr0(0), RecMultDist05inelgr0(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliSPDppAnalysisTaskSim::~AliSPDppAnalysisTaskSim()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliSPDppAnalysisTaskSim::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    fQAList = new TList();
    fQAList -> SetOwner();
    
    if (fUseINT1) { fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT1,true); }
    
    
    // example of a histogram
    RecMultDist05 = new TH1F("RecMultDist05", "RecMultDist05", 51, -.5, 50.5);
    MultDist05 = new TH1F("MultDist05", "MultDist05", 51, -.5, 50.5);       // create your histogra
    responseMatrix = new TH2F("responseMatrix", "Response Matrix;Generated Multiplicity;Reconstructed Multiplicity", 51, -0.5, 50.5, 51, -0.5, 50.5);
    
    RecMultDist05inelgr0 = new TH1F("RecMultDist05inelgr0", "RecMultDist05inelgr0", 51, -.5, 50.5);
    MultDist05inelgr0 = new TH1F("MultDist05inelgr0", "MultDist05inelgr0", 51, -.5, 50.5);       // create your histogra
    responseMatrixinelgr0 = new TH2F("responseMatrixinelgr0", "Response Matrix;Generated Multiplicity;Reconstructed Multiplicity", 51, -0.5, 50.5, 51, -0.5, 50.5);
    
    fOutputList->Add(responseMatrix);
    fOutputList->Add(RecMultDist05);
    fOutputList->Add(MultDist05);          // don't forget to add it to the list! the list will be written to file, so if you want
    fOutputList->Add(responseMatrixinelgr0);
    fOutputList->Add(RecMultDist05inelgr0);
    fOutputList->Add(MultDist05inelgr0);
                                        // your histogram in the output file, add it to the list!
    fEventCuts.AddQAplotsToList(fOutputList);
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliSPDppAnalysisTaskSim::UserExec(Option_t *)
{
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
        // check if there actually is an event:
    if(!fAOD) { return; }
    
    
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if(fMCEvent) ProcessMCParticles();
    if(!fMCEvent) ProcessData();
    
    if (!fEventCuts.AcceptEvent(fAOD)) { return; }
    
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}

void AliSPDppAnalysisTaskSim::ProcessMCParticles(){
    
//    fAODV0 = fAOD->GetVZEROData();
//    float fV0Amult = fAODV0->GetMTotV0A(); //returns total multiplicity in V0A
//    float fV0Cmult = fAODV0->GetMTotV0C(); //returns total multiplicity in V0C
    
    const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    if(!spdVtx) return;
    Double_t number = spdVtx->GetNContributors();
    float spdVtxZ = spdVtx->GetZ();
    //if(TMath::Abs(spdVtxZ) > 10) return;
    
    //if (fV0Amult == 0 && fV0Cmult == 0 && number == 0) return;

    TClonesArray *stack = nullptr;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) { return; }

//    AliAnalysisUtils util;
//    util.SetMinPlpContribMV(5);
//    util.SetMaxPlpChi2MV(5);
//    util.SetMinWDistMV(15);
//    util.SetCheckPlpFromDifferentBCMV(kFALSE);
    
    fMultiplicity = fAOD -> GetMultiplicity();
    Int_t nTracklets = fMultiplicity->GetNumberOfTracklets();
    
    //tracklet loop
    Int_t reconstructedMultiplicity = 0;
    Int_t reconstructedMultiplicity1 = 0;
    for (auto it = 0; it<nTracklets; it++) {
        
        Double_t eta = fMultiplicity->GetEta(it);
        if (TMath::Abs(eta) < .5){ //removes tracklets outside |eta|<.5
            reconstructedMultiplicity++;
        }
        if (TMath::Abs(eta) < 1){
            reconstructedMultiplicity1++;
        }
    }
    

    
    /*Int_t reconstructedMultiplicity = 0;
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if (!track) continue;
        //if (!track->TestFilterBit(1)) continue;
        int label = TMath::Abs(track->GetLabel());
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
        
        Double_t vz1 = mcTrack->Zv();
        //if (TMath::Abs(vz1) > 10.0) { continue; }
        
    
        if (!mcTrack) continue;
        if (!mcTrack->IsPhysicalPrimary()) continue;
        if(TMath::Abs(mcTrack->Charge()) < 0.1)continue;
        
        Double_t Pt = mcTrack->Pt();
        if(Pt<.15)continue;
        if(Pt>5)continue;
        
        Double_t eta = mcTrack->Eta();
        if (TMath::Abs(eta) < 0.5) {
            reconstructedMultiplicity++;
        }
    }*/
    
    
    
    
    
    int nMCTracks;
        if (!stack) nMCTracks = 0;
        else nMCTracks = stack->GetEntries();

        AliAODMCHeader *mcHeader = 0;
        mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
        if(!mcHeader) {
            printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
            return;
        }
    
    TClonesArray* AODMCTrackArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
        if (!AODMCTrackArray){
            Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
            this->Dump();
            return;
        }

    
    TArrayF MC_Vtx_true_XYZ(3);
    fMCEvent->GenEventHeader()->PrimaryVertex(MC_Vtx_true_XYZ);
    Double_t VertexZ = MC_Vtx_true_XYZ[2];
    if(TMath::Abs(VertexZ) > 10)return;
    
    
    Int_t TrueMCParticles = 0;
    Int_t TrueMCParticles1 = 0;
    for (Int_t mcTrack = 0;  mcTrack < (fMCEvent->GetNumberOfTracks()); mcTrack++){
        
        AliMCParticle *particle = (AliMCParticle*)fMCEvent->GetTrack(mcTrack);
        if (!particle) continue;
        
        //if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(mcTrack, mcHeader, stack)) { continue; }
        
        
        Double_t vz = particle->Zv();
        //if (TMath::Abs(vz) > 10.0) { continue; }
        
        Bool_t TrackIsPrim = particle->IsPhysicalPrimary();
        if (!TrackIsPrim)continue;
        
        Double_t pt = particle->Pt();
        //if(pt<.15)continue;
        //if(pt>5)continue;
        //cout<<pt<<endl;
        Double_t trueMCeta = particle->Eta();
        Double_t trueMCcharge = particle->Charge();
        if(trueMCcharge == 0)continue;
        
        if(TMath::Abs(trueMCeta) < 0.5){
            TrueMCParticles++;
        }
        if(TMath::Abs(trueMCeta) < 1){
            TrueMCParticles1++;
        }
    }

    
    
    // Generated Mult vs Reconstructed Mult
    // Fill the response matrix for the first half of the events
    if (eventcount1 % 2 == 0) {
        responseMatrix->Fill(reconstructedMultiplicity, TrueMCParticles);
        if(TrueMCParticles1 != 0){
            responseMatrixinelgr0->Fill(reconstructedMultiplicity, TrueMCParticles);
        }
        //cout<<eventcount1<<endl;
    }
    // Fill the new multiplicity matrix for the other half of the events
    else {
        RecMultDist05->Fill(reconstructedMultiplicity);
        if(reconstructedMultiplicity1 != 0){
            RecMultDist05inelgr0->Fill(reconstructedMultiplicity);
        }
        MultDist05->Fill(TrueMCParticles);
        if(TrueMCParticles1 != 0){
            MultDist05inelgr0->Fill(TrueMCParticles);
        }
        
    }
    eventcount1++;
    
    
    
    
}

void AliSPDppAnalysisTaskSim::ProcessData(){
    
    
    
}


    
    
    
//_____________________________________________________________________________
void AliSPDppAnalysisTaskSim::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
