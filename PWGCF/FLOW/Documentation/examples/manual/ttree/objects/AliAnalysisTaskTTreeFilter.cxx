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

/* analysis task which extracts some kinematic info from AliVEvents in the 
 * aliroot analysis framework and stores then in a ttree
 * see macros in PWGCF/FLOW/Documentation/examples/manual/ttree/macros
 * for usage info
 * author: redmer alexander bertens (rbertens@cern.ch)
 * based on code received from Alexandru Dorbin
 */


#include "AliAnalysisTaskTTreeFilter.h"

// ROOT includes
#include <TTree.h>
#include <TMath.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TString.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliLog.h>
#include <AliVEvent.h>
#include <AliVVertex.h>
#include <AliCentrality.h>
#include <AliVTrack.h>

// local includes
#include "AliFlowTTreeEvent.h"
#include "AliFlowTTreeTrack.h"

ClassImp(AliAnalysisTaskTTreeFilter)

//_____________________________________________________________________________
AliAnalysisTaskTTreeFilter::AliAnalysisTaskTTreeFilter():
    AliAnalysisTaskSE(),
    fEvent(0x0),
    fTrackArray(0x0)
{
    // default constructor for root I/O
}
//______________________________________________________________________________
AliAnalysisTaskTTreeFilter::AliAnalysisTaskTTreeFilter(const char *name):
    AliAnalysisTaskSE(name),
    fEvent(0x0),
    fTrackArray(0x0)
{
    // constructor
    DefineOutput(1, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTTreeFilter::~AliAnalysisTaskTTreeFilter()
{
    // destructor
    if (fTree) {
        delete fTree;
        fTree = 0x0;
    }
}
//______________________________________________________________________________
void AliAnalysisTaskTTreeFilter::UserCreateOutputObjects()
{ 
    // check for manager
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler) AliFatal(" > no input detected - aborting <");

    // open file at slot 1 for large output to avoid buffer overflow
    AliAnalysisTask::OpenFile(1);

    // create the ttree
    fTree = new TTree("tree", "Event data");

    // init the custom event 
    fEvent = new AliFlowTTreeEvent();
    // add the event branch to the tree
    fTree->Branch("event", &fEvent);
  
    // init the track tclonesarray
    fTrackArray = new TClonesArray("AliFlowTTreeTrack", 1000);
    // add clones array as branch via bronch
    fTree->Bronch("track", "TClonesArray", &fTrackArray);

    // Post output data.
    PostData(1, fTree);
}
//______________________________________________________________________________
void AliAnalysisTaskTTreeFilter::UserExec(Option_t *) 
{
    // parse accepted input event
    if(ParseEvent(dynamic_cast<AliVEvent*>(InputEvent()))) PostData(1, fTree);
    else return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTTreeFilter::ParseEvent(AliVEvent* event) 
{
    // parse the input event
    if(!PassesCuts(event)) return kFALSE;

    // store some event info
    fEvent->SetRun(event->GetRunNumber());
    fEvent->SetV0M(event->GetCentrality()->GetCentralityPercentile("V0M"));
    fEvent->SetTRK(event->GetCentrality()->GetCentralityPercentile("TRK"));
    fEvent->SetZvtx(event->GetPrimaryVertex()->GetZ());
  
    // parse the tracks
    ParseTracks(event);

    // write the tree and perform cleanup
    PushToTTree();
    
    // jay !
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskTTreeFilter::ParseTracks(AliVEvent* event)
{
    // parse tracks
  
    for(Int_t i(0), acceptedTracks(0); i < event->GetNumberOfTracks(); i++) {
        // track loop
        AliVTrack* track(static_cast<AliVTrack*>(event->GetTrack(i)));
        if(!PassesCuts(track)) continue;

        // push accepted track to tree
        AliFlowTTreeTrack* acceptedTrack = new((*fTrackArray)[acceptedTracks]) AliFlowTTreeTrack();
        acceptedTracks++;
        // add info
        acceptedTrack->SetPt(track->Pt());
        acceptedTrack->SetEta(track->Eta());
        acceptedTrack->SetPhi(track->Phi());
        acceptedTrack->SetCharge(track->Charge());
    }
    return;
}
//________________________________________________________________________
void AliAnalysisTaskTTreeFilter::PushToTTree()
{
    // push info to tree and do cleanup for next iteration
    fTree->Fill();
    fTrackArray->Clear();
}
//________________________________________________________________________
Bool_t AliAnalysisTaskTTreeFilter::PassesCuts(AliVEvent* event)
{
    // event cuts would go here
    if(!event) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskTTreeFilter::PassesCuts(AliVTrack* track)
{
    // track cuts would go here
    if(!track) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskTTreeFilter::Terminate(Option_t *)
{ 
    // terminate
}
