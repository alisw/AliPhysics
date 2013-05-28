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
 * AliAnaysisTaskJetFlow
 * author: Redmer Alexander Bertens
 * rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 
 *
 * Interface task between EMCal jet framework and flow package
 *
 * This task expects POI's in a TClonesArray  (e.g. from running it after 
 * $ALICE_ROOT/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskRhoVnModulation() )
 * and connects them into an AliFlowEvent which is filled with either VZERO tracks or 
 * TPC trakcs. 
 * 
 * POI's can be supplied as AliEmcalJets or as AliVParticles, note the different behavior
 * with respect to the Pt value: for AliEmcalJets subtracted Pt is used by default, for VParticles
 * Pt is taken directly.
 *
 * */

// root includes
#include <TChain.h>
#include <TH1F.h>
#include <TArrayD.h>
#include <TProfile.h>
#include <TString.h>
#include <TList.h>
// aliroot includes
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliEmcalJet.h>
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliFlowEventSimple.h>
#include <AliFlowEvent.h>
#include <AliFlowTrack.h>
#include <AliFlowTrackCuts.h>
#include <AliFlowEventCuts.h>
#include <AliFlowCommonConstants.h>
// local includes
#include "AliAnalysisTaskJetFlow.h"

class AliAnalysisTaskJetFlow;

using namespace std;

ClassImp(AliAnalysisTaskJetFlow)

AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow() : AliAnalysisTaskSE(), 
    fDebug(-1), fJetsName(0), fOutputList(0), fDataType(kESD), fVParticleAnalysis(kFALSE), fDoTestFlowAnalysis(kFALSE), fInitialized(kFALSE), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPOIPtMin(0.15), fPOIPtMax(150), fPtBins(0), fCutsRP_TPC(0), fCutsRP_VZERO(0), fCutsPOI(0), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fHistAnalysisSummary(0), fCentralitySelection(0), fv2VZEROA(0), fv2VZEROC(0)
{ /* default constructor for ROOT IO */ }
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow(const char* name) : AliAnalysisTaskSE(name),
    fDebug(-1), fJetsName(0), fOutputList(0), fDataType(kESD), fVParticleAnalysis(kFALSE), fDoTestFlowAnalysis(kFALSE), fInitialized(kFALSE), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPOIPtMin(0.15), fPOIPtMax(150), fPtBins(0), fCutsRP_TPC(0), fCutsRP_VZERO(0), fCutsPOI(0), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fHistAnalysisSummary(0), fCentralitySelection(0), fv2VZEROA(0), fv2VZEROC(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, AliFlowEventSimple::Class());
    DefineOutput(3, AliFlowEventSimple::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::~AliAnalysisTaskJetFlow()
{
    // destructor
    if(fOutputList)             delete fOutputList;
    if(fFlowEvent_TPC)          delete fFlowEvent_TPC;
    if(fFlowEvent_VZERO)        delete fFlowEvent_VZERO;
    if(fCutsEvent)              delete fCutsEvent;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::LocalInit()
{
    // executed once
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
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
    if(fDoTestFlowAnalysis) { // set up the binning for the test flow analysis
        if((!fPtBins) && fVParticleAnalysis) {
            Double_t pt[51];
            for(Int_t i(0); i < 51; i++) pt[i] = .1*i;
            fPtBins = new TArrayD(51, pt);  // assume they will be charged particles
        } else if (!fPtBins) {
            Double_t pt[] = {0., 10., 20., 30., 40., 50., 80., 110., 140., 170., 200.};
            fPtBins = new TArrayD(sizeof(pt)/sizeof(pt[0]), pt);     // assuming jets
        }
        Double_t bounds[fPtBins->GetSize()];
        for(Int_t i(0); i < fPtBins->GetSize(); i++) bounds[i] = fPtBins->At(i);
        fv2VZEROA = new TProfile("v2_EP_VZEROA", "v2_EP_VZEROA", fPtBins->GetSize()-1, bounds);
        fv2VZEROC = new TProfile("v2_EP_VZEROC", "v2_EP_VZEROC", fPtBins->GetSize()-1, bounds);
        fv2VZEROA->GetXaxis()->SetTitle("Pt [GeV/c]");
        fv2VZEROA->GetYaxis()->SetTitle("v_{2}^{obs}");
        fv2VZEROC->GetXaxis()->SetTitle("Pt [GeV/c]");
        fv2VZEROC->GetYaxis()->SetTitle("v_{2}^{obs}");
        fOutputList->Add(fv2VZEROA);
        fOutputList->Add(fv2VZEROC);
    }
    // qa
    fCentralitySelection = new TH1F("fCentralitySelection", "fCentralitySelection", 100, 0, 100);
    fOutputList->Add(fCentralitySelection);
    PostData(1, fOutputList);
    // create the flow event and configure the static cc object

    (fVParticleAnalysis) ? fFlowEvent_VZERO = new AliFlowEvent(10000) : fFlowEvent_VZERO = new AliFlowEvent(100);
    PostData(2, fFlowEvent_VZERO);
    fFlowEvent_TPC = new AliFlowEvent(10000);
    PostData(3, fFlowEvent_TPC);
    AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
    cc->SetPtMax(fCCMaxPt+fPtBump);
    cc->SetNbinsPt(fCCBinsInPt);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::UserExec(Option_t *)
{
    // user exec
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!InputEvent() || !fCutsNull || !fCutsRP_TPC || !fCutsRP_VZERO) return; // coverity (and sanity)
    if(!fInitialized) { 
        if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
        else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
        fInitialized = kTRUE;
    }
    if(!PassesCuts(InputEvent())) return;               // check the event cuts
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
        fCutsRP_VZERO->SetEvent(InputEvent(), MCEvent());
        fCutsRP_TPC->SetEvent(InputEvent(), MCEvent());
        fCutsNull->SetEvent(InputEvent(), MCEvent());
        fFlowEvent_VZERO->ClearFast();
        fFlowEvent_TPC->ClearFast();
        // the event is filled with rp's only, poi's will be added manually
        fFlowEvent_VZERO->Fill(fCutsRP_VZERO, fCutsNull);
        fFlowEvent_TPC->Fill(fCutsRP_TPC, fCutsNull);
        fFlowEvent_VZERO->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
        fFlowEvent_TPC->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
        // loop over jets and inject them as POI's
        if(fVParticleAnalysis) {
            for(Int_t i(0); i < iJets; i++) {
                AliVParticle* jet = static_cast<AliVParticle*>(jets->At(i));
                if(jet) {
                    if(jet->Pt() + fPtBump <= fPOIPtMin || jet->Pt() > fPOIPtMax) {
                        fHistAnalysisSummary->SetBinContent(4, 1);
                        continue;
                    }
                    nAcceptedJets++;
                    AliFlowTrack* flowTrack = new AliFlowTrack(jet);
                    flowTrack->SetPt(jet->Pt() + fPtBump);
                    flowTrack->SetForPOISelection(kTRUE);
                    flowTrack->SetForRPSelection(kFALSE);
                    fFlowEvent_TPC->InsertTrack(flowTrack);
                    fFlowEvent_VZERO->InsertTrack(flowTrack);
                }
            }
        } else {
            for(Int_t i(0); i < iJets; i++) {
                AliEmcalJet* jet = static_cast<AliEmcalJet*>(jets->At(i));
                if(jet) {
                    if(jet->PtSub() + fPtBump <= fPOIPtMin || jet->PtSub() > fPOIPtMax) {
                        fHistAnalysisSummary->SetBinContent(4, 1);
                        continue;
                    }
                    nAcceptedJets++;
                    AliFlowTrack* flowTrack = new AliFlowTrack(jet);
                    flowTrack->SetPt(jet->PtSub() + fPtBump);
                    flowTrack->SetForPOISelection(kTRUE);
                    flowTrack->SetForRPSelection(kFALSE);
                    fFlowEvent_TPC->InsertTrack(flowTrack);
                    fFlowEvent_VZERO->InsertTrack(flowTrack);
                }
            }
        }
    }
    else if(fDebug > 0 ) printf(" Failed to find TClones Jet array while using name %s \n ", fJetsName.Data());
    if(nAcceptedJets < 1) {
        if(fDebug > 0) printf(" > No accepted jets in event ! < " );
        return;
    }
    fFlowEvent_TPC->TagSubeventsInEta(-10, 0, 0, 10);
    fFlowEvent_VZERO->TagSubeventsInEta(-10, 0, 0, 10);
    if(fDoTestFlowAnalysis) DoTestFlowAnalysis();
    fCentralitySelection->Fill(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    PostData(1, fOutputList);
    PostData(2, fFlowEvent_VZERO);
    PostData(3, fFlowEvent_TPC);
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
void AliAnalysisTaskJetFlow::DoTestFlowAnalysis()
{
    // get a crude estimate of v2 based on the event plane method FIXME not tested !!!
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t _a(0), _b(0), _c(0), _d(0);        // dummmy's
    Double_t Q2a(InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, _a, _b));
    Double_t Q2c(InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, _c, _d));
    TProfile* a = (TProfile*)fv2VZEROA->Clone("temp_a");
    TProfile* c = (TProfile*)fv2VZEROC->Clone("temp_c");
    if(!(a||c)) return; // coverity
    TClonesArray* jets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName.Data()));
    if(jets) {
        Int_t iJets = jets->GetEntriesFast();
        if(fVParticleAnalysis) {
            for(Int_t i(0); i < iJets; i++) {
                AliVParticle* jet = static_cast<AliVParticle*>(jets->At(i));
                if(jet && jet->Pt() + fPtBump >= fPOIPtMin && jet->Pt() < fPOIPtMax) {
                    a->Fill(jet->Pt(), TMath::Cos(2.*(jet->Phi()-Q2a)));
                    c->Fill(jet->Pt(), TMath::Cos(2.*(jet->Phi()-Q2c)));
                }
            }
        } else {
            for(Int_t i(0); i < iJets; i++) {
                AliEmcalJet* jet = static_cast<AliEmcalJet*>(jets->At(i));
                if(jet && jet->PtSub() + fPtBump >= fPOIPtMin && jet->PtSub() < fPOIPtMax) {
                    a->Fill(jet->PtSub(), TMath::Cos(2.*(jet->Phi()-Q2a)));
                    c->Fill(jet->PtSub(), TMath::Cos(2.*(jet->Phi()-Q2c)));
                }
            }
        }
        for(Int_t i(0); i < fv2VZEROA->GetXaxis()->GetNbins(); i++) {
            fv2VZEROA->Fill(fPtBins->At(i)+(fPtBins->At(i)+fPtBins->At(1+i))/2., a->GetBinContent(i+1));
            fv2VZEROC->Fill(fPtBins->At(i)+(fPtBins->At(i)+fPtBins->At(1+i))/2., c->GetBinContent(i+1));
        }
        delete a;
        delete c;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::Terminate(Option_t *)
{
    // terminate
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
}
//_____________________________________________________________________________
