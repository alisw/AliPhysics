/*************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

// author: I.Lofnes, ingrid.mckibben.lofnes@cern.ch
// Analysis task for determining the vdm-scan stability
//

#include <iostream>
using namespace std;

#include <TChain.h>
#include <TH1D.h>

#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliTriggerAnalysis.h>
#include <AliESDVZERO.h>

#include "AliAnalysisTaskVdmStability.h"
ClassImp(AliAnalysisTaskVdmStability)

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability():
AliAnalysisTaskSE(),
fOutputList(),
fRunNumber(0),
fVtxZ(-999.),
fnVtxCont(0),
fIsGoodZ(kFALSE),
fSelectPhysics(kFALSE),
fIsV0ANDfired(kFALSE),
fIsT0fired(kFALSE),
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0)
{
    // ROOT IO constructor, don't allocate memory here!
}

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability(const char* taskname):
AliAnalysisTaskSE(taskname),
fOutputList(),
fRunNumber(0),
fVtxZ(-999.),
fnVtxCont(0),
fIsGoodZ(kFALSE),
fSelectPhysics(kFALSE),
fIsV0ANDfired(kFALSE),
fIsT0fired(kFALSE),
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());    //TList of event statistics
    DefineOutput(2, TTree::Class());    //event tree information
    
}

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::~AliAnalysisTaskVdmStability()
{
    if(fEventStatV0)       { delete fEventStatV0;       fEventStatV0=0; }
    if(fEventStatT0)       { delete fEventStatT0;       fEventStatT0=0; }
    if (fEventTree)         { delete fEventTree;        fEventTree=0; }
}

//defíne output
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::UserCreateOutputObjects()
{
    // set ownership of list
    fOutputList.SetOwner();
    // ---| skip list initialisation if already done |----------------------------
    if (!fOutputList.IsEmpty()) return;
    
    //Initialize output tree
    //fEventTree = std::unique_ptr<TTree>(new TTree("event", "event"));
    fEventTree = new TTree("event", "event");
    AddEventTreeVariables(fEventTree);
    
    
    // create our histo and add it to the list
    Int_t nbins=kNbinsEvent;
    if (!fEventStatV0){
        fEventStatV0=new TH1D("hEventStatV0","Event statistics V0",nbins,0,nbins);
        fEventStatV0->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
        fEventStatV0->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
        fEventStatV0->GetXaxis()->SetBinLabel(3,"V0 timing");
        fEventStatV0->GetXaxis()->SetBinLabel(4,"PU");
        fEventStatV0->GetXaxis()->SetBinLabel(5,"good z");
        fEventStatV0->GetXaxis()->SetBinLabel(6,"V0 timing + PU");
        fEventStatV0->GetXaxis()->SetBinLabel(7,"V0 timing + good z");
        fEventStatV0->GetXaxis()->SetBinLabel(8,"V0 timing + PU + good z");
        
    }
    fOutputList.Add(fEventStatV0);
    
    if (!fEventStatT0){
        fEventStatT0=new TH1D("hEventStatT0","Event statistics T0",nbins,0,nbins);
        fEventStatT0->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
        fEventStatT0->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
        fEventStatT0->GetXaxis()->SetBinLabel(3,"V0 timing");
        fEventStatT0->GetXaxis()->SetBinLabel(4,"PU");
        fEventStatT0->GetXaxis()->SetBinLabel(5,"good z");
        fEventStatT0->GetXaxis()->SetBinLabel(6,"V0 timing + PU");
        fEventStatT0->GetXaxis()->SetBinLabel(7,"V0 timing + good z");
        fEventStatT0->GetXaxis()->SetBinLabel(8,"V0 timing + PU + good z");
        
    }
    fOutputList.Add(fEventStatT0);
    
    // add the list to our output file
    PostData(1,&fOutputList);
    PostData(2,fEventTree);
}

//Event loop
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::UserExec(Option_t *)
{
    //Get analysis manager
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
    Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
    
    if (isAOD) {
        cout << "Error: AOD not inititalized yet " << endl;
        return;
    }
    
    //Get input handler
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if (!inputHandler) return;
    
    //get event
    fEvent = static_cast<AliESDEvent*>(InputEvent());
    
    //Check if V0 event
    TString ftc = fEvent->GetFiredTriggerClasses();
    if (ftc.Contains("CINT7-B"))
        fIsV0ANDfired = kTRUE;
    else fIsV0ANDfired = kFALSE;
    if (!fIsV0ANDfired) return;
    
    //check if T0 trigger was fired
    fIsT0fired = fEvent->GetHeader()->IsTriggerInputFired("0TVX");
    
    //Check physics selection
    fSelectPhysics = kTRUE;
    if( inputHandler){
        if((isESD && inputHandler->GetEventSelection())){
            fSelectPhysics = inputHandler->IsEventSelected();
        }
    }
    
    
    //Set remaining event variables
    SetEventVariables();
    
    
    //Fill event statistics histograms
    //-----------------------------------------------------
    //Trigger information
    fEventStatV0->Fill(kAllEvents);
    if (fIsT0fired)
        fEventStatT0->Fill(kAllEvents);
    //Physics selection
    if (fSelectPhysics) {
        fEventStatV0->Fill(kSelectedEvents);
        if (fIsT0fired) fEventStatT0->Fill(kSelectedEvents);
    }
    //V0 timing
    if ((fSelectPhysics) && (fGoodTime)) {
        fEventStatV0->Fill(kV0TimeEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0TimeEvents);
    }
    //pileup
    if ((fSelectPhysics) && !(fPileupEvent)) {
        fEventStatV0->Fill(kPileupEvents);
        if (fIsT0fired) fEventStatT0->Fill(kPileupEvents);
    }
    //Z-vertex
    if ((fSelectPhysics) && (fIsGoodZ)) {
        fEventStatV0->Fill(kGoodZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kGoodZEvents);
    }
    //Timing + PU
    if ((fSelectPhysics) && (fGoodTime) && !(fPileupEvent) ) {
        fEventStatV0->Fill(kV0andPUEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andPUEvents);
    }
    //Timing + z
    if ((fSelectPhysics) && (fGoodTime) && (fIsGoodZ) ) {
        fEventStatV0->Fill(kV0andZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andZEvents);
    }
    //Timing + PU + z
    if ((fSelectPhysics) && (fGoodTime) && !(fPileupEvent) && (fIsGoodZ) ) {
        fEventStatV0->Fill(kV0andPUandZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andPUandZEvents);
    }
    
    //Fill tree with event information
    fEventTree->Fill();
    
    PostData(1,&fOutputList);
    PostData(2,fEventTree);
    
}

// Set branch names for event tree variables
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::AddEventTreeVariables(TTree* &tree)
{
    tree->Branch("RunNumber", &fRunNumber);
    tree->Branch("VtxZ", &fVtxZ);
    tree->Branch("nVtxCont", &fnVtxCont);
    //tree->Branch("goodZvertex",&fIsGoodZ);
    tree->Branch("PhysSelected",&fSelectPhysics);
    tree->Branch("V0Trigger",&fIsV0ANDfired);
    tree->Branch("T0Trigger",&fIsT0fired);
    tree->Branch("PileupEvent",&fPileupEvent);
    tree->Branch("timeV0A",&fTV0A);
    tree->Branch("timeV0C",&fTV0C);
    tree->Branch("goodTime",&fGoodTime);
}


// Set the event variables
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::SetEventVariables() {
    
    //set run number
    fRunNumber = static_cast<UInt_t >(InputEvent()->GetRunNumber());
    
    //Set z-vertex information
    const AliESDVertex* trackVtx    = fEvent->GetPrimaryVertexTracks();
    fVtxZ = trackVtx->GetZ();
    fnVtxCont = trackVtx->GetNContributors();
    if (!(fVtxZ < -10 || fVtxZ > 10) && fnVtxCont > 1)
        fIsGoodZ = kTRUE;
    else fIsGoodZ = kFALSE;
    
    //Set V0 timing
    AliESDVZERO * esdV0 = fEvent->GetVZEROData();
    fTV0A = esdV0->GetV0ATime();
    fTV0C = esdV0->GetV0CTime();
    Float_t tV0sum = fTV0A + fTV0C;
    Float_t tV0diff = fTV0A - fTV0C;
    
    if (((4 < tV0diff) && (tV0diff < 12) ) && ((12 < tV0sum)&& (tV0sum < 18)))
        fGoodTime = kTRUE;
    else fGoodTime = kFALSE;
    
    //set pileup
    if (fEvent->IsPileupFromSPD(3,0.8,3.,2.,5.))
        fPileupEvent = kTRUE;
    else fPileupEvent = kFALSE;
    
}


//Terminate
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::Terminate(Option_t *)
{
    
}

