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

/* 
 * AliAnaysisTaskJet
 * author: Redmer Alexander Bertens
 * rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 
 *
 * jet correlation task - correlates jets to the vzero reaction plane using
 * the scalar product method 
 *
 * this task should be run AFTER 
 * $ALICE_ROOT/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskRhoVnModulation()
 *
 * */

#include <TChain.h>
#include <TH1F.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliFlowEventSimple.h>
#include <AliFlowTrack.h>
#include <AliFlowTrackCuts.h>
#include <AliFlowEventCuts.h>
#include <AliFlowCommonConstants.h>
#include <TString.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliEmcalJet.h>

#include "AliAnalysisTaskJetFlow.h"

class AliAnalysisTaskJetFlow;

using namespace std;

ClassImp(AliAnalysisTaskJetFlow)

AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow() : AliAnalysisTaskSE(), 
    fDebug(-1), fJetsName(0), fOutputList(0), fDataType(kESD), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fCutsRP(0), fCutsPOI(0), fCutsNull(0), fCutsEvent(0), fFlowEvent(0), fHistAnalysisSummary(0)
{ /* default constructor */ }
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow(const char* name) : AliAnalysisTaskSE(name),
    fDebug(-1), fJetsName(0), fOutputList(0), fDataType(kESD), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fCutsRP(0), fCutsPOI(0), fCutsNull(0), fCutsEvent(0), fFlowEvent(0), fHistAnalysisSummary(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, AliFlowEventSimple::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::~AliAnalysisTaskJetFlow()
{
    // destructor
    if(fOutputList)     delete fOutputList;
    if(fFlowEvent)      delete fFlowEvent;
    if(fCutsEvent)      delete fCutsEvent;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::LocalInit()
{
    // executed once
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
    else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
    fCutsEvent = new AliFlowEventCuts();
    fCutsEvent->SetRefMultMethod(AliESDtrackCuts::kTrackletsITSTPC); 
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::UserCreateOutputObjects()
{
    // create output objects
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    // histograms
    fHistAnalysisSummary = new TH1F("fHistAnalysisSummary", "fHistAnalysisSummary", 4, -.5, 3.5);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(1, "fDataType");
    fHistAnalysisSummary->SetBinContent(1, (int)fDataType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fCentralityMin");
    fHistAnalysisSummary->SetBinContent(2, fCentralityMin);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fCentralityMax");
    fHistAnalysisSummary->SetBinContent(3, fCentralityMax);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "pt bias");
    fOutputList->Add(fHistAnalysisSummary);
    PostData(1, fOutputList);
    // create the flow event and configure the static cc object
    fFlowEvent = new AliFlowEvent(1000);
    PostData(2, fFlowEvent);
    AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
    cc->SetPtMax(fCCMaxPt+fPtBump);
    cc->SetNbinsPt(fCCBinsInPt);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::UserExec(Option_t *)
{
    // user exec
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    
    if(!InputEvent() || !fCutsNull || !fCutsRP) return; // coverity (and sanity)
    // get the jet array, which is added as an extension to the AliVEvent by the jetfinder
    TClonesArray* jets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName.Data()));
    Int_t nAcceptedJets(0);
    if(jets) {
        Int_t iJets = jets->GetEntriesFast();
        if(iJets <= 0) {
            if(fDebug>0) printf(" > Retrieved empty AliVEvent extension, aborting ! < \n ");
            return;
        }
        // prepare the track selection for RP's and the flow event
        fCutsRP->SetEvent(InputEvent(), MCEvent());
        fCutsNull->SetEvent(InputEvent(), MCEvent());
        fFlowEvent->ClearFast();
        fFlowEvent->Fill(fCutsRP, fCutsNull);
        fFlowEvent->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
        // loop over jets and inject them as POI's
        for(Int_t i(0); i < iJets; i++) {
            AliEmcalJet* jet = static_cast<AliEmcalJet*>(jets->At(i));
            if(jet) {
                if(jet->PtSub() + fPtBump <= .15) {
                    fHistAnalysisSummary->SetBinContent(4, 1);
                    continue;
                }
                nAcceptedJets++;
                AliFlowTrack* flowTrack = new AliFlowTrack(jet);
                flowTrack->SetPt(jet->PtSub() + fPtBump);
                flowTrack->SetForPOISelection(kTRUE);
                flowTrack->SetForRPSelection(kFALSE);
                fFlowEvent->InsertTrack(flowTrack);
            }
        }
    }
    else if(fDebug > 0 ) printf(" Failed to find TClones Jet array while using name %s \n ", fJetsName.Data());
    if(nAcceptedJets < 1) {
        if(fDebug > 0) printf(" > No accepted jets in event ! < " );
        return;
    }
    PostData(1, fOutputList);
    PostData(2, fFlowEvent);
    if(fDebug>0) {
        printf(" Event summary \n");
        printf(" > number of POI's (jets) %i ", fFlowEvent->NumberOfTracks() - fFlowEvent->GetNumberOfRPs());
        printf(" > number of RP's %i \n", fFlowEvent->GetNumberOfRPs());
    }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetFlow::PassesCuts(AliVEvent* event)
{
    // event cuts
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!event) return kFALSE;
    if(TMath::Abs(InputEvent()->GetPrimaryVertex()->GetZ()) > 10.) return kFALSE;
    // aod and esd specific checks
    switch (fDataType) {
       case kESD: {
            AliESDEvent* esdEvent = static_cast<AliESDEvent*>(InputEvent());
            if( (!esdEvent) || (TMath::Abs(esdEvent->GetPrimaryVertexSPD()->GetZ() - esdEvent->GetPrimaryVertex()->GetZ()) > .5) ) return kFALSE; 
       } break;
       case kAOD: {
            AliAODEvent* aodEvent = static_cast<AliAODEvent*>(InputEvent());
            if( (!aodEvent) || (TMath::Abs(aodEvent->GetPrimaryVertexSPD()->GetZ() - aodEvent->GetPrimaryVertex()->GetZ()) > .5) ) return kFALSE; 
       } break;
       default: break;
    }
    Float_t cent(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    return (cent <= fCentralityMin || cent > fCentralityMax || TMath::Abs(cent-InputEvent()->GetCentrality()->GetCentralityPercentile("TRK")) > 5.) ? kFALSE : kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::Terminate(Option_t *)
{
    // terminate
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
}
//_____________________________________________________________________________
