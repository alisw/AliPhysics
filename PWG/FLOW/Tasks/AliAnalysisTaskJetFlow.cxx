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
    fDebug(-1), fExplicitOutlierCut(-1), fJetsName(0), fTracksName(0), fOutputList(0), fDataType(kESD), fVParticleAnalysis(kFALSE), fMinimizeDiffBins(kTRUE), fDoTestFlowAnalysis(kTRUE), fDoMultWeight(kTRUE), fInitialized(kFALSE), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPOIPtMin(0.15), fPOIPtMax(150), fPtBins(0), fCutsRP_TPC(0), fCutsRP_VZERO(0), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fHistAnalysisSummary(0), fCentralitySelection(0), fVZEROAEP(0), fVZEROCEP(0), fv2VZEROA(0), fv2VZEROC(0), fTempA(0), fTempC(0)
{ /* default constructor for ROOT IO */ }
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow(
        const char* name,
        AliFlowTrackCuts* rpCutsTPC,
        AliFlowTrackCuts* rpCutsVZERO,
        TString jetName,
        TString picoName  ) : AliAnalysisTaskSE(name),
    fDebug(-1), fExplicitOutlierCut(-1), fJetsName(jetName), fTracksName(picoName), fOutputList(0), fDataType(kESD), fVParticleAnalysis(kFALSE), fMinimizeDiffBins(kTRUE), fDoTestFlowAnalysis(kTRUE), fDoMultWeight(kTRUE), fInitialized(kFALSE), fPtBump(0), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPOIPtMin(0.15), fPOIPtMax(150), fPtBins(0), fCutsRP_TPC(rpCutsTPC), fCutsRP_VZERO(rpCutsVZERO), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fHistAnalysisSummary(0), fCentralitySelection(0), fVZEROAEP(0), fVZEROCEP(0), fv2VZEROA(0), fv2VZEROC(0), fTempA(0), fTempC(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    if(fCutsRP_VZERO || fCutsRP_TPC)    DefineOutput(2, AliFlowEventSimple::Class());
    if(fCutsRP_VZERO && fCutsRP_TPC)    DefineOutput(3, AliFlowEventSimple::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::~AliAnalysisTaskJetFlow()
{
    // destructor
    if(fOutputList)             delete fOutputList;
    if(fFlowEvent_TPC)          delete fFlowEvent_TPC;
    if(fFlowEvent_VZERO)        delete fFlowEvent_VZERO;
    if(fCutsEvent)              delete fCutsEvent;
    if(fCutsRP_VZERO)           delete fCutsRP_VZERO;
    if(fCutsRP_TPC)             delete fCutsRP_TPC;
    if(fCutsNull)               delete fCutsNull;
    if(fPtBins)                 delete fPtBins;
    if(fTempA)                  delete fTempA;
    if(fTempC)                  delete fTempC;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::LocalInit()
{
    // executed once
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fCutsRP_VZERO || fCutsRP_TPC) {
        fCutsEvent = new AliFlowEventCuts();
        fCutsEvent->SetRefMultMethod(AliESDtrackCuts::kTrackletsITSTPC); 
        fCutsNull = new AliFlowTrackCuts("CutsNull");
        fCutsNull->SetParamType(AliFlowTrackCuts::kGlobal);
        fCutsNull->SetEtaRange(+1, -1);
        fCutsNull->SetPtRange(+1, -1);
    }
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
        fv2VZEROA = new TProfile("v2_EP_VZEROA", "v2_EP_VZEROA", fPtBins->GetSize()-1, fPtBins->GetArray());
        fv2VZEROC = new TProfile("v2_EP_VZEROC", "v2_EP_VZEROC", fPtBins->GetSize()-1, fPtBins->GetArray());
        fv2VZEROA->GetXaxis()->SetTitle("Pt [GeV/c]");
        fv2VZEROA->GetYaxis()->SetTitle("v_{2}^{obs}");
        fv2VZEROC->GetXaxis()->SetTitle("Pt [GeV/c]");
        fv2VZEROC->GetYaxis()->SetTitle("v_{2}^{obs}");
        fOutputList->Add(fv2VZEROA);
        fOutputList->Add(fv2VZEROC);
        fVZEROAEP = new TH1F("V0A_EP", "V0A_EP", 100, -TMath::Pi()/2., TMath::Pi()/2.);
        fVZEROCEP = new TH1F("V0C_EP", "V0C_EP", 100, -TMath::Pi()/2., TMath::Pi()/2.);
        fOutputList->Add(fVZEROAEP);
        fOutputList->Add(fVZEROCEP);
        // bookkeeping histo's, will not be stored
        fTempA = (TProfile*)fv2VZEROA->Clone("temp_a");
        fTempC = (TProfile*)fv2VZEROC->Clone("temp_c");
    }
    // qa
    fCentralitySelection = new TH1F("fCentralitySelection", "fCentralitySelection", 100, 0, 100);
    fOutputList->Add(fCentralitySelection);
    PostData(1, fOutputList);
    // create the flow event and configure the static cc object
    Bool_t slotTwoFilled(kFALSE);
    if(fCutsRP_VZERO) {
        (fVParticleAnalysis) ? fFlowEvent_VZERO = new AliFlowEvent(10000) : fFlowEvent_VZERO = new AliFlowEvent(100);
        PostData(2, fFlowEvent_VZERO);
        slotTwoFilled = kTRUE;
    }
    if(fCutsRP_TPC) {
        fFlowEvent_TPC = new AliFlowEvent(10000);
        (slotTwoFilled) ? PostData(3, fFlowEvent_TPC) : PostData(2, fFlowEvent_TPC);
    }
    if(fFlowEvent_VZERO || fFlowEvent_TPC) {
        AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
        cc->SetPtMax(fCCMaxPt+fPtBump);
        cc->SetNbinsPt(fCCBinsInPt);
        if(fMinimizeDiffBins) { // minimize differential bins to reduce the risk of numerical instability
            cc->SetNbinsMult(1);    // only reduces output size
            cc->SetNbinsPhi(1);     // only reduces output size
            cc->SetNbinsEta(1);     // reduces instability
            cc->SetNbinsQ(1);       // only reduces output size
            cc->SetNbinsMass(1);    // reduces instability but should be one in any case ...
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::UserExec(Option_t *)
{
    // user exec
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!InputEvent()) return; // coverity (and sanity)
    if(!fInitialized) { 
        if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
        else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
        fInitialized = kTRUE;
    }
    if(!PassesCuts(InputEvent())) return;               // check the event cuts
    // get the jet array, which is added as an extension to the AliVEvent by the jetfinder
    TClonesArray* pois(0x0);
    (fVParticleAnalysis) ? pois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName.Data())) : pois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName.Data()));
    Int_t nAcceptedJets(0);
    if(pois) {
        Int_t iPois = pois->GetEntriesFast();
        if(iPois <= 0) {
            if(fDebug>0) printf(" > Retrieved empty AliVEvent extension, aborting ! < \n ");
            return;
        }
        // prepare the track selection for RP's and the flow event
        if(fFlowEvent_VZERO) {
            fCutsNull->SetEvent(InputEvent(), MCEvent());
            fCutsRP_VZERO->SetEvent(InputEvent(), MCEvent());
            fFlowEvent_VZERO->ClearFast();
            fFlowEvent_VZERO->Fill(fCutsRP_VZERO, fCutsNull);
            fFlowEvent_VZERO->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
        }
        if(fFlowEvent_TPC) {
            fCutsNull->SetEvent(InputEvent(), MCEvent());
            fCutsRP_TPC->SetEvent(InputEvent(), MCEvent());
            fFlowEvent_TPC->ClearFast();
            fFlowEvent_TPC->Fill(fCutsRP_TPC, fCutsNull);
            fFlowEvent_TPC->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
        }
        // the event is filled with rp's only, poi's will be added manually
        // loop over jets and inject them as POI's
        if(fFlowEvent_VZERO || fFlowEvent_TPC) {
            if(fVParticleAnalysis) {
                for(Int_t i(0); i < iPois; i++) {
                    AliVParticle* poi = static_cast<AliVParticle*>(pois->At(i));
                    if(poi) {
                        if(poi->Pt() + fPtBump <= fPOIPtMin || poi->Pt() > fPOIPtMax) {
                            fHistAnalysisSummary->SetBinContent(4, 1);
                            continue;
                        }
                        nAcceptedJets++;
                        // AliFlowTracks are created on the stack 
                        AliFlowTrack flowTrack = AliFlowTrack(poi);
                        flowTrack.SetPt(poi->Pt() + fPtBump);
                        flowTrack.SetForPOISelection(kTRUE);
                        flowTrack.SetForRPSelection(kFALSE);
                        if(fFlowEvent_TPC)      fFlowEvent_TPC->InsertTrack(&flowTrack);
                        if(fFlowEvent_VZERO)    fFlowEvent_VZERO->InsertTrack(&flowTrack);
                    }
                }
            } else {
                for(Int_t i(0); i < iPois; i++) {
                    AliEmcalJet* poi = static_cast<AliEmcalJet*>(pois->At(i));
                    if(poi) {
                        if(poi->PtSub() + fPtBump <= fPOIPtMin || poi->PtSub() > fPOIPtMax) {
                            fHistAnalysisSummary->SetBinContent(4, 1);
                            continue;
                        }
                        nAcceptedJets++;
                        AliFlowTrack flowTrack = AliFlowTrack(poi);
                        flowTrack.SetPt(poi->PtSub() + fPtBump);
                        flowTrack.SetForPOISelection(kTRUE);
                        flowTrack.SetForRPSelection(kFALSE);
                        if(fFlowEvent_TPC)      fFlowEvent_TPC->InsertTrack(&flowTrack);
                        if(fFlowEvent_VZERO)    fFlowEvent_VZERO->InsertTrack(&flowTrack);
                    }
                }
            }
        }
        else if(fDebug > 0 ) printf(" Failed to find TClones Jet array while using name %s \n ", fJetsName.Data());
        if((fFlowEvent_VZERO || fFlowEvent_TPC) && nAcceptedJets < 1) {
            if(fDebug > 0) printf(" > No accepted jets in event ! < " );
            return;
        }
        if(fFlowEvent_TPC)      fFlowEvent_TPC->TagSubeventsInEta(-10, 0, 0, 10);
        if(fFlowEvent_VZERO)    fFlowEvent_VZERO->TagSubeventsInEta(-10, 0, 0, 10);
    }
    if(fDoTestFlowAnalysis) DoTestFlowAnalysis();
    fCentralitySelection->Fill(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    PostData(1, fOutputList);
    Bool_t slotTwoFilled(kFALSE);
    if(fFlowEvent_VZERO) {
        PostData(2, fFlowEvent_VZERO);
        slotTwoFilled = kTRUE;
    }
    if(fFlowEvent_TPC) {
        (slotTwoFilled) ? PostData(3, fFlowEvent_TPC) : PostData(2, fFlowEvent_TPC);
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
    if(fExplicitOutlierCut == 2010 || fExplicitOutlierCut == 2011) {
       if(!PassesCuts(fExplicitOutlierCut)) return kFALSE;
    }
    return (cent <= fCentralityMin || cent > fCentralityMax || TMath::Abs(cent-InputEvent()->GetCentrality()->GetCentralityPercentile("TRK")) > 5.) ? kFALSE : kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetFlow::PassesCuts(Int_t year) 
{
    // additional centrality cut based on relation between tpc and global multiplicity
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    AliAODEvent* event(dynamic_cast<AliAODEvent*>(InputEvent()));
    if(!event) return kFALSE;
    Int_t multTPC(0), multGlob(0), nTracks(InputEvent()->GetNumberOfTracks());
    for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) { 
        AliAODTrack* track = event->GetTrack(iTracks);
        if(!track) continue;
        if (!track || track->Pt() < .2 || track->Pt() > 5.0 || TMath::Abs(track->Eta()) > .8 || track->GetTPCNcls() < 70 || !track->GetDetPid() || track->GetDetPid()->GetTPCsignal() < 10.0)  continue;  // general quality cut
        if (track->TestFilterBit(1) && track->Chi2perNDF() > 0.2) multTPC++;
        if (!track->TestFilterBit(16) || track->Chi2perNDF() < 0.1) continue;
        Double_t b[2] = {-99., -99.};
        Double_t bCov[3] = {-99., -99., -99.};
        if (track->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) multGlob++;
    }
    if(year == 2010 && multTPC > (-40.3+1.22*multGlob) && multTPC < (32.1+1.59*multGlob)) return kTRUE;
    if(year == 2011  && multTPC > (-36.73 + 1.48*multGlob) && multTPC < (62.87 + 1.78*multGlob)) return kTRUE;
    return kFALSE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::DoTestFlowAnalysis()
{
    // get a crude estimate of v2 based on the event plane method FIXME not tested !!!
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t _a(0), _b(0), _c(0), _d(0);        // dummmy's
    Double_t Q2a(InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, _a, _b));
    Double_t Q2c(InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, _c, _d));
    fTempA->Reset();            fTempC->Reset();        // clear the containers for a new iteration
    fVZEROAEP->Fill(Q2a);       fVZEROCEP->Fill(Q2c);
    TClonesArray* pois(0x0);
    (fVParticleAnalysis) ? pois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName.Data())) : pois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName.Data()));
    if(pois) {
        Int_t iPois = pois->GetEntriesFast();
        if(fVParticleAnalysis) {
            for(Int_t i(0); i < iPois; i++) {
                AliVParticle* poi = static_cast<AliVParticle*>(pois->At(i));
                if(poi && poi->Pt() + fPtBump >= fPOIPtMin && poi->Pt() < fPOIPtMax) {
                    if(!fDoMultWeight) {
                        fTempA->Fill(poi->Pt(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2a), 2)));
                        fTempC->Fill(poi->Pt(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2c), 2)));
                    } else {
                        fv2VZEROA->Fill(poi->Pt(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2a), 2)));
                        fv2VZEROC->Fill(poi->Pt(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2c), 2)));
                    }
                }
            }
        } else {
            for(Int_t i(0); i < iPois; i++) {
                AliEmcalJet* poi = static_cast<AliEmcalJet*>(pois->At(i));
                if(poi && poi->PtSub() + fPtBump >= fPOIPtMin && poi->PtSub() < fPOIPtMax) {
                    if(!fDoMultWeight) {
                        fTempA->Fill(poi->PtSub(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2a), 2)));
                        fTempC->Fill(poi->PtSub(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2c), 2)));
                    } else {
                        fv2VZEROA->Fill(poi->PtSub(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2a), 2)));
                        fv2VZEROC->Fill(poi->PtSub(), TMath::Cos(2.*PhaseShift((poi->Phi()-Q2c), 2)));
                    }
                }
            }
        }
        if(!fDoMultWeight) {
            for(Int_t i(0); i < fv2VZEROA->GetXaxis()->GetNbins(); i++) {
                fv2VZEROA->Fill(fPtBins->At(i)+(fPtBins->At(i)+fPtBins->At(1+i))/2., fTempA->GetBinContent(i+1));
                fv2VZEROC->Fill(fPtBins->At(i)+(fPtBins->At(i)+fPtBins->At(1+i))/2., fTempC->GetBinContent(i+1));
            }
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::Terminate(Option_t *)
{
    // terminate
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
}
//_____________________________________________________________________________
