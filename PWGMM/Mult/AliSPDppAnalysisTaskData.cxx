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
#include "AliSPDppAnalysisTaskData.h"

#include "TTree.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliHeader.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDVZERO.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPWG0Helper.h"
#include "AliEventCuts.h"

class AliSPDppAnalysisTaskData;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliSPDppAnalysisTaskData) // classimp: necessary for root

AliSPDppAnalysisTaskData::AliSPDppAnalysisTaskData() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), MultDist05(0), fEventCuts(), fAODV0(0), fUseINT1(true), MultDist05Inelgr0(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliSPDppAnalysisTaskData::AliSPDppAnalysisTaskData(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), MultDist05(0), fEventCuts(), fAODV0(0), fUseINT1(true), MultDist05Inelgr0(0)
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
AliSPDppAnalysisTaskData::~AliSPDppAnalysisTaskData()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliSPDppAnalysisTaskData::UserCreateOutputObjects()
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
    MultDist05 = new TH1F("MultDist05", "MultDist05", 51, -.5, 50.5);       // create your histogra
    MultDist05Inelgr0 = new TH1F("MultDist05Inelgr0", "MultDist05Inelgr0", 51, -.5, 50.5);
    fOutputList->Add(MultDist05);          // don't forget to add it to the list! the list will be written to file, so if you want
    fOutputList->Add(MultDist05Inelgr0);
                                        // your histogram in the output file, add it to the list!
    fEventCuts.AddQAplotsToList(fOutputList);
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliSPDppAnalysisTaskData::UserExec(Option_t *)
{
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file

    
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it
    
    
    //event cuts
    if (!fEventCuts.AcceptEvent(fAOD)) {
            return;
        }

    
//    AliAnalysisUtils util;
//    util.SetMinPlpContribMV(5);
//    util.SetMaxPlpChi2MV(5);
//    util.SetMinWDistMV(15);
//    util.SetCheckPlpFromDifferentBCMV(kFALSE);
    
  
//    fAODV0 = fAOD->GetVZEROData();
//    float fV0Amult = fAODV0->GetMTotV0A(); //returns total multiplicity in V0A
//    float fV0Cmult = fAODV0->GetMTotV0C(); //returns total multiplicity in V0C
    
    const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    if(!spdVtx) return;
    Double_t number = spdVtx->GetNContributors();
    float spdVtxZ = spdVtx->GetZ();
    //if(TMath::Abs(spdVtxZ) > 10) return;
    
    //if (fV0Amult == 0 && fV0Cmult == 0 && number == 0) return;
    
    //Gets tracklets
    fMultiplicity = fAOD -> GetMultiplicity();
    Int_t nTracklets = fMultiplicity->GetNumberOfTracklets();
    
    //tracklet loop
    Int_t Tracklets05 = 0;
    Int_t Tracklets1 = 0;
    for (auto it = 0; it<nTracklets; it++) {
        
        Double_t eta = fMultiplicity->GetEta(it);
        if (TMath::Abs(eta) < .5){ //removes tracklets outside |eta|<.5
            Tracklets05++;
        }
        if (TMath::Abs(eta) < 1){
            Tracklets1++;
        }
    }
    
    if (Tracklets1 != 0){
        MultDist05Inelgr0->Fill(Tracklets05);
    }
    
    
    /*Int_t iTracks(fAOD->GetNumberOfTracks());
    Int_t tracks1 = 0;
    for(Int_t a(0); a < iTracks; a++) {                 // loop ove rall these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(a));         // get a track (type AliAODTrack) from the event
        if(!track || !track->TestFilterBit(128)) continue;                            // if we failed, skip this track
        
        //if (std::abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3 )continue; //removes everything but pions
        if(track->Eta() > .5)continue;                  //cuts on ETA range right
        if(track->Eta() < -.5)continue;                 //cuts on ETA range right
                          //cuts on ETA range left
        tracks1++;
     
        //cout<<x_1<<endl;
        
    }
    
    MultDist05->Fill(tracks1);*/
    
    MultDist05->Fill(Tracklets05); //fills multiplicity distribution for |eta|<.5
    

    
    
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}
//_____________________________________________________________________________
void AliSPDppAnalysisTaskData::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
