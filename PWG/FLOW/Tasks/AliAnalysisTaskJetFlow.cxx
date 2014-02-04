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
 * POI's can be supplied as AliEmcalJets or as AliVTracks, note the different behavior
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
#include <TClonesArray.h>
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
// EMCAL jet framework includes
#include <AliRhoParameter.h>
#include <AliLocalRhoParameter.h>
#include <AliPicoTrack.h>
#include <AliAnalysisTaskRhoVnModulation.h>

class AliAnalysisTaskJetFlow;

using namespace std;

ClassImp(AliAnalysisTaskJetFlow)

AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow() : AliAnalysisTaskSE(), 
    fDebug(-1), fJetsName(0), fJetRadius(0.3), fTracksName(0), fLocalRhoName(0), fPois(0x0), fRPs(0x0), fLocalRho(0x0), fOutputList(0), fDataType(kESD), fVParticleAnalysis(kFALSE), fMinimizeDiffBins(kTRUE), fDoVZEROFlowAnalysis(kTRUE), fDoGappedQC2Analysis(kTRUE), fDoQC2FlowAnalysis(kTRUE), fDoQC4FlowAnalysis(kFALSE), fDoQCFPAnalysis(kFALSE), fDoSPFPAnalysis(kFALSE), fDoMultWeight(kTRUE), fDoPtWeight(0), fInitialized(kFALSE), fUsePtWeight(kFALSE), fCCMinPt(1), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPtBins(0), fCutsRP_VZERO(0), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fRhoVn(0), fHistAnalysisSummary(0), fCentralitySelection(0), fVZEROAEP(0), fVZEROCEP(0), fv2VZEROA(0), fv2VZEROC(0), fRefCumulants(0), fDiffCumlantsV2(0), fDiffCumlantsV3(0), fQC2v2(0), fQC2v3(0), fTempA(0), fTempC(0)
{ /* default constructor for ROOT IO */ }
//_____________________________________________________________________________
AliAnalysisTaskJetFlow::AliAnalysisTaskJetFlow(
        const char* name,
        AliAnalysisTaskRhoVnModulation* rhoTask, 
        Bool_t VPart,
        Bool_t VZEROEP,
        Bool_t GQC2, 
        Bool_t QC2,
        Bool_t QC4,
        Bool_t FlowPackageSP,
        Bool_t FlowPackageQC  
        ) : AliAnalysisTaskSE(name),
    fDebug(-1), fJetsName(0), fJetRadius(0.3), fTracksName(0), fLocalRhoName(0), fPois(0x0), fRPs(0x0), fLocalRho(0x0), fOutputList(0), fDataType(kESD), fVParticleAnalysis(VPart), fMinimizeDiffBins(kTRUE), fDoVZEROFlowAnalysis(VZEROEP), fDoGappedQC2Analysis(GQC2), fDoQC2FlowAnalysis(QC2), fDoQC4FlowAnalysis(QC4), fDoQCFPAnalysis(FlowPackageQC), fDoSPFPAnalysis(FlowPackageSP), fDoMultWeight(kTRUE), fDoPtWeight(0), fInitialized(kFALSE), fUsePtWeight(kFALSE), fCCMinPt(1), fCCMaxPt(150), fCCBinsInPt(50), fCentralityMin(-1), fCentralityMax(-1), fPtBins(0), fCutsRP_VZERO(0x0), fCutsNull(0), fCutsEvent(0), fFlowEvent_TPC(0), fFlowEvent_VZERO(0), fRhoVn(rhoTask), fHistAnalysisSummary(0), fCentralitySelection(0), fVZEROAEP(0), fVZEROCEP(0), fv2VZEROA(0), fv2VZEROC(0), fRefCumulants(0), fDiffCumlantsV2(0), fDiffCumlantsV3(0), fQC2v2(0), fQC2v3(0), fTempA(0), fTempC(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    fJetsName = rhoTask->GetJetsName();
    fTracksName = rhoTask->GetTracksName(); 
    if(GQC2 && QC2) {
        printf(" > Warning, QC2 and gapped QC2 method are both called <\n   will only run gapped QC2 !");
        fDoQC2FlowAnalysis = kFALSE;
    }
    if(FlowPackageSP || FlowPackageQC)    DefineOutput(2, AliFlowEventSimple::Class());
    if(FlowPackageSP && FlowPackageQC)    DefineOutput(3, AliFlowEventSimple::Class());
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
    if(fCutsNull)               delete fCutsNull;
    if(fPtBins)                 delete fPtBins;
    if(fTempA)                  delete fTempA;
    if(fTempC)                  delete fTempC;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::UserCreateOutputObjects()
{
    // create output objects
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // create the cut objects
    if(fDoSPFPAnalysis || fDoQCFPAnalysis) {
        fCutsEvent = new AliFlowEventCuts();
        fCutsEvent->SetRefMultMethod(AliESDtrackCuts::kTrackletsITSTPC); 
        fCutsNull = new AliFlowTrackCuts("CutsNull");
        fCutsNull->SetParamType(AliFlowTrackCuts::kGlobal);
        fCutsNull->SetEtaRange(+1, -1);
        fCutsNull->SetPtRange(+1, -1);
        fCutsRP_VZERO = new AliFlowTrackCuts();
        fCutsRP_VZERO = fCutsRP_VZERO->GetStandardVZEROOnlyTrackCuts();
        if(fDoSPFPAnalysis) {
            (fVParticleAnalysis) ? fFlowEvent_VZERO = new AliFlowEvent(10000) : fFlowEvent_VZERO = new AliFlowEvent(100);
        }
        if(fDoQCFPAnalysis) fFlowEvent_TPC = new AliFlowEvent(10000);
    }
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
    fOutputList->Add(fHistAnalysisSummary);
    if(fVParticleAnalysis && ! fPtBins) { // FIXME inherits from flow package now
        Double_t pt[fCCBinsInPt+1];
        Double_t width = (fCCMaxPt - fCCMinPt ) / (float)fCCBinsInPt;
        for(Int_t i(0); i < fCCBinsInPt+1; i++) pt[i] = fCCMinPt+width*i;
        fPtBins = new TArrayD(fCCBinsInPt+1, pt);  // assume they will be charged particles
    } else if (!fPtBins) {
        Double_t pt[] = {0., 10., 20., 30., 40., 50., 80., 110., 140., 170., 200.};
        fPtBins = new TArrayD(sizeof(pt)/sizeof(pt[0]), pt);     // assuming jets
    }
    if(fDoVZEROFlowAnalysis) {
        fv2VZEROA = new TProfile("Differential v_{2}^{obs} VZEROA", "Differential v_{2}^{obs} VZEROA", fPtBins->GetSize()-1, fPtBins->GetArray());
        fv2VZEROC = new TProfile("Differential v_{2}^{obs} VZEROC", "Differential v_{2}^{obs} VZEROC", fPtBins->GetSize()-1, fPtBins->GetArray());
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
    if(fDoQC2FlowAnalysis || fDoGappedQC2Analysis) {
        fRefCumulants = new TProfile("Reference cumulants", "Reference cumulants", 2, -0.5, 1.5);
        fRefCumulants->GetXaxis()->SetBinLabel(1, "c_{2}[2]");
        fRefCumulants->GetXaxis()->SetBinLabel(2, "c_{3}[2]");
        fOutputList->Add(fRefCumulants);
        fDiffCumlantsV2 = new TProfile("Differential cumulants v2", "Differential cumulants v2", fPtBins->GetSize()-1, fPtBins->GetArray());
        fDiffCumlantsV2->GetXaxis()->SetTitle("Pt [GeV/c]");
        fDiffCumlantsV2->GetYaxis()->SetTitle("v_{2}[2]");
        fOutputList->Add(fDiffCumlantsV2);
        fDiffCumlantsV3 = new TProfile("Differential cumulants v3", "Differential cumulants v3", fPtBins->GetSize()-1, fPtBins->GetArray());
        fDiffCumlantsV3->GetXaxis()->SetTitle("Pt [GeV/c]");
        fDiffCumlantsV3->GetYaxis()->SetTitle("v_{3}[2]");
        fOutputList->Add(fDiffCumlantsV3);
        fQC2v2 = new TH1F("Differential v_{2}[2]", "Differential v_{2}[2]", fPtBins->GetSize()-1, fPtBins->GetArray());
        fQC2v2->Sumw2();
        fQC2v3 = new TH1F("Differential v_{3}[2]", "Differential v_{3}[2]", fPtBins->GetSize()-1, fPtBins->GetArray());
        fQC2v3->Sumw2();
        fOutputList->Add(fQC2v2); 
        fOutputList->Add(fQC2v3);
    }
    // qa
    fCentralitySelection = new TH1F("fCentralitySelection", "fCentralitySelection", 100, 0, 100);
    fOutputList->Add(fCentralitySelection);
    PostData(1, fOutputList);
    // create the flow event and configure the static cc object
    Bool_t slotTwoFilled(kFALSE);
    if(fFlowEvent_VZERO) {
        PostData(2, fFlowEvent_VZERO);
        slotTwoFilled = kTRUE;
    }
    if(fFlowEvent_TPC) {
        (slotTwoFilled) ? PostData(3, fFlowEvent_TPC) : PostData(2, fFlowEvent_TPC);
    }
    if(fFlowEvent_VZERO || fFlowEvent_TPC) {
        AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
        cc->SetPtMin(fCCMinPt);
        cc->SetPtMax(fCCMaxPt);
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
    // see if the analysis is initialized 
    if(!fInitialized) { 
        if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
        else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
        (fVParticleAnalysis) ? fPois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName.Data())) : fPois = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName.Data()));
        fRPs = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName.Data()));
        if(!fPois || !fRPs) return; // couldn't get expected input data
        fInitialized = kTRUE;
        fLocalRho = InputEvent()->FindListObject(fLocalRhoName.Data());
        if(!fLocalRho && !fVParticleAnalysis) {
            AliFatal(Form("Couldn't find %s, aborting!", fLocalRhoName.Data()));
        }
    }
    if(!PassesCuts()) return; // event quality cuts and centrality determination
    // execute the requested flow methods
    if(fDoVZEROFlowAnalysis)            DoVZEROFlowAnalysis();
    if(fDoGappedQC2Analysis)            DoGappedQC2Analysis();
    if(fDoQC2FlowAnalysis)              DoQC2FlowAnalysis();
    if(fDoQC4FlowAnalysis)              DoQC4FlowAnalysis();
    Bool_t post(0x0);   // post only when analysis succeeded
    if(fFlowEvent_TPC || fFlowEvent_VZERO) DoFlowPackageFlowAnalysis();
    // push the output for the different analyses to file
    PostData(1, fOutputList);
    if(fFlowEvent_VZERO) {
        PostData(2, fFlowEvent_VZERO);
        post = kTRUE;
    }
    if(fFlowEvent_TPC) { 
        (post) ? PostData(3, fFlowEvent_TPC) : PostData(2, fFlowEvent_TPC);
    }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetFlow::PassesCuts() 
{
    // passes event cuts
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fRhoVn->PassesCuts(InputEvent())) return kFALSE;          // event quality cuts
    if(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M") < fCentralityMin || InputEvent()->GetCentrality()->GetCentralityPercentile("V0M") > fCentralityMax) return kFALSE;  // cutting on centrality quality is done in event quality cuts
    fCentralitySelection->Fill(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::DoVZEROFlowAnalysis()
{
    // flow with the VZERO event plane method
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t vzero[2][2];
    fRhoVn->CalculateEventPlaneVZERO(vzero);
    Double_t Q2a(vzero[0][0]), Q2c(vzero[1][0]);        // just for readability
    fTempA->Reset();            fTempC->Reset();        // clear the containers for a new iteration
    fVZEROAEP->Fill(Q2a);       fVZEROCEP->Fill(Q2c);
    Int_t iPois(fPois->GetEntriesFast());
    if(fVParticleAnalysis) {
        for(Int_t i(0); i < iPois; i++) {
            AliVTrack* poi = static_cast<AliVTrack*>(fPois->At(i));
            if(fRhoVn->PassesCuts(poi)) {
                if(!fDoMultWeight) {
                    fTempA->Fill(poi->Pt(), TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2a), 2)));
                    fTempC->Fill(poi->Pt(), TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2c), 2)));
                } else {
                    fv2VZEROA->Fill(poi->Pt(), TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2a), 2)));
                    fv2VZEROC->Fill(poi->Pt(), TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2c), 2)));
                }
            }
        }
    } else {
        for(Int_t i(0); i < iPois; i++) {
            AliEmcalJet* poi = static_cast<AliEmcalJet*>(fPois->At(i));
            if(fRhoVn->PassesCuts(poi)) {
                AliLocalRhoParameter* localRho = static_cast<AliLocalRhoParameter*>(fLocalRho);
                if(!localRho) break;
                Double_t rho(localRho->GetLocalVal(poi->Phi(), fJetRadius, localRho->GetVal()));
                Double_t pt(poi->Pt() - poi->Area() * rho);
                if(!fDoMultWeight) {
                    fTempA->Fill(pt, TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2a), 2)));
                    fTempC->Fill(pt, TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2c), 2)));
                } else {
                    fv2VZEROA->Fill(pt, TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2a), 2)));
                    fv2VZEROC->Fill(pt, TMath::Cos(2.*fRhoVn->PhaseShift((poi->Phi()-Q2c), 2)));
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
 //_____________________________________________________________________________
void AliAnalysisTaskJetFlow::DoGappedQC2Analysis()
{
    // do q-cumulant analysis with eta gaps (avoiding autocorrelation of rps and jet constituents)
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // first step, get lhs tpc rp's
    fRhoVn->SetTrackEtaLimits(-0.9, -0.7);
    // get LHS rp's multiplicity and q-vector
    Double_t LHSreQn(0), LHSimQn(0), LHSmQ(0);
    (fDoPtWeight) ? fRhoVn->QCnQnk(2, 1, LHSreQn, LHSimQn) : fRhoVn->QCnQnk(2, 0, LHSreQn, LHSimQn);
    (fDoPtWeight) ? LHSmQ = fRhoVn->QCnM11() : LHSmQ = fRhoVn->QCnM();
    // get the RHS rp's multiplicity and q-vector
    fRhoVn->SetTrackEtaLimits(0.7, 0.9);
    Double_t RHSreQn(0), RHSimQn(0), RHSmQ(0);
    (fDoPtWeight) ? fRhoVn->QCnQnk(2, 1, RHSreQn, RHSimQn) : fRhoVn->QCnQnk(2, 0, RHSreQn, RHSimQn);
    (fDoPtWeight) ? RHSmQ = fRhoVn->QCnM11() : RHSmQ = fRhoVn->QCnM();
    // differential flow vectors
    Double_t repn[fPtBins->GetSize()-1];      // real part of q-vector of all poi's
    Double_t impn[fPtBins->GetSize()-1];      // im part of q-vector of all poi's
    Double_t mp[fPtBins->GetSize()-1];        // poi multiplicity
    Double_t reqn[fPtBins->GetSize()-1];      // real part of q-vectors of poi's labeled as rp
    Double_t imqn[fPtBins->GetSize()-1];      // im part of q-vectors of poi's labeled as rp
    Double_t mq[fPtBins->GetSize()-1];        // multiplicity of poi's labeled as rp
    for(Int_t i(0); i < fPtBins->GetSize(); i++) {
        repn[i] = 0;
        impn[i] = 0;
        mp[i] = 0;
        reqn[i] = 0;
        imqn[i] = 0;
        mq[i] = 0;
    }
    // calculate differential q-vectors and fill the profile with cumulants
    (fVParticleAnalysis) ? fRhoVn->SetTrackEtaLimits(-0.7, 0.7) : fRhoVn->SetLocalJetMinMaxEta(fJetRadius+.2);       // avoid overlap in poi and rp region 
    QCnDifferentialFlowVectors(repn, impn, mp, reqn, imqn, mq, 2);
    // do the calculation
    if(RHSmQ*LHSmQ < 1) {
        fRhoVn->SetTrackEtaLimits(-0.9, 0.9);
        fRhoVn->SetLocalJetMinMaxEta(fJetRadius);
        return;
    }
    fRefCumulants->Fill(0., (LHSreQn*RHSreQn+LHSimQn*RHSimQn)/(RHSmQ*LHSmQ), RHSmQ*LHSmQ);
    for(Int_t i(0); i < fPtBins->GetSize(); i++) {
        if(LHSmQ*mp[i] < 1. ) continue;        // avoid division by zero
        Double_t atPt(fPtBins->At(i)+0.5*(fPtBins->At(i+1)-fPtBins->At(i)));      // pt value
        Double_t diffC((repn[i]*LHSreQn+impn[i]*LHSimQn)/(LHSmQ*mp[i]));
        Double_t eventW(mp[i]*LHSmQ);
        fDiffCumlantsV2->Fill(atPt, diffC, eventW);
    }
    // last step: roll back the eta cuts of the EmcalTask
    fRhoVn->SetTrackEtaLimits(-0.9, 0.9);
    fRhoVn->SetLocalJetMinMaxEta(fJetRadius);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::DoQC2FlowAnalysis()
{
    // flow analysis with the qc2 method
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // reference flow is taken from the pico track selection and evaluated in 
    // the AliAnalysisTaskiRhoVnModulation class
    for(Int_t harm(2); harm < 4; harm++) {        // loop over harmonics
        Double_t reQn(0), imQn(0), mQ(0);         // total q-vector
        // get the total q-vectors and reference cumulants
        (fDoPtWeight) ? fRhoVn->QCnQnk(harm, 1, reQn, imQn) : fRhoVn->QCnQnk(harm, 0, reQn, imQn);
        (fDoPtWeight) ? mQ = fRhoVn->QCnM11() : mQ = fRhoVn->QCnM();
        if(mQ < 2) continue;    // avoid division by zero
        fRefCumulants->Fill((double)(harm-2), ((reQn*reQn+imQn*imQn)-mQ)/(mQ*(mQ-1)), mQ*(mQ-1));
        // differential flow vectors
        Double_t repn[fPtBins->GetSize()-1];      // real part of q-vector of all poi's
        Double_t impn[fPtBins->GetSize()-1];      // im part of q-vector of all poi's
        Double_t mp[fPtBins->GetSize()-1];        // poi multiplicity
        Double_t reqn[fPtBins->GetSize()-1];      // real part of q-vectors of poi's labeled as rp
        Double_t imqn[fPtBins->GetSize()-1];      // im part of q-vectors of poi's labeled as rp
        Double_t mq[fPtBins->GetSize()-1];        // multiplicity of poi's labeled as rp
        for(Int_t i(0); i < fPtBins->GetSize(); i++) {
            repn[i] = 0;
            impn[i] = 0;
            mp[i] = 0;
            reqn[i] = 0;
            imqn[i] = 0;
            mq[i] = 0;
        }
        // calculate differential q-vectors and fill the profile with cumulants
        QCnDifferentialFlowVectors(repn, impn, mp, reqn, imqn, mq, harm);
        // FIXME differential evnet weights
        for(Int_t i(0); i < fPtBins->GetSize(); i++) {
            if(mp[i]*mQ - mq[i] <= 0 ) continue;        // avoid division by zero
            Double_t atPt(fPtBins->At(i)+0.5*(fPtBins->At(i+1)-fPtBins->At(i)));      // pt value
            Double_t diffC(((repn[i]*reQn+impn[i]*imQn)-mq[i])/(mp[i]*mQ-mq[i]));
            Double_t eventW(mp[i]*mQ-mq[i]);
            (harm == 2 ) ? fDiffCumlantsV2->Fill(atPt, diffC, eventW) : fDiffCumlantsV3->Fill(atPt, diffC, eventW);
        }
    } 
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::DoQC4FlowAnalysis()
{
    // flow analysis with the qc4 method - see comments at qc2 FIXME not implemented yet
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskJetFlow::DoFlowPackageFlowAnalysis() 
{
    // name's a bit misleading: this function does anlaysis using FlowAnalysis classes in /PWG/FLOW/Base/
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    // get the jet array, which is added as an extension to the AliVEvent by the jetfinder
    Int_t nAcceptedJets(0), iPois(fPois->GetEntriesFast());
    if(iPois <= 0) return kFALSE;
    if(fFlowEvent_VZERO) {
        fCutsNull->SetEvent(InputEvent(), MCEvent());
        fCutsRP_VZERO->SetEvent(InputEvent(), MCEvent());
        fFlowEvent_VZERO->ClearFast();
        // the event is filled with rp's only, poi's will be added manually
        fFlowEvent_VZERO->Fill(fCutsRP_VZERO, fCutsNull);
        fFlowEvent_VZERO->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
    }
    if(fFlowEvent_TPC) {
        fCutsNull->SetEvent(InputEvent(), MCEvent());
        // in this case, both poi's and rp's will be added manually
        fFlowEvent_TPC->ClearFast();
        fFlowEvent_TPC->SetReferenceMultiplicity(fCutsEvent->RefMult(InputEvent(), MCEvent()));
    }
    // loop over jets and inject them as POI's
    if(fVParticleAnalysis) {
        for(Int_t i(0); i < iPois; i++) {
            AliVTrack* poi = static_cast<AliVTrack*>(fPois->At(i));
            if(fRhoVn->PassesCuts(poi)) {
                nAcceptedJets++;
                // AliFlowTracks are created on the stack 
                AliFlowTrack flowTrack = AliFlowTrack(poi);
                flowTrack.SetForPOISelection(kTRUE);
                if(fFlowEvent_TPC) {
                    fFlowEvent_TPC->SetNumberOfRPs(fFlowEvent_TPC->GetNumberOfRPs()+1);
                    fFlowEvent_TPC->SetNumberOfPOIs(fFlowEvent_TPC->GetNumberOfPOIs()+1);
                    flowTrack.SetForRPSelection(kTRUE);
                    fFlowEvent_TPC->InsertTrack(&flowTrack);
                }
                if(fFlowEvent_VZERO) {
                    flowTrack.SetForRPSelection(kFALSE);
                    fFlowEvent_VZERO->InsertTrack(&flowTrack);
                    fFlowEvent_VZERO->SetNumberOfPOIs(fFlowEvent_VZERO->GetNumberOfPOIs()+1);
                }
            }
        }
    } else {
        // add the jets as pois
        for(Int_t i(0); i < iPois; i++) {
            AliEmcalJet* poi = static_cast<AliEmcalJet*>(fPois->At(i));
            if(fRhoVn->PassesCuts(poi)) {
                nAcceptedJets++;
                AliFlowTrack flowTrack = AliFlowTrack(poi);
                AliLocalRhoParameter* localRho = static_cast<AliLocalRhoParameter*>(fLocalRho);
                if(!localRho) break;
                Double_t rho(localRho->GetLocalVal(poi->Phi(), fJetRadius, localRho->GetVal()));
                flowTrack.SetPt(poi->Pt() - poi->Area() * rho);
                flowTrack.SetForPOISelection(kTRUE);
                flowTrack.SetForRPSelection(kFALSE);
                if(fFlowEvent_TPC) {
                    fFlowEvent_TPC->InsertTrack(&flowTrack);
                    fFlowEvent_TPC->SetNumberOfPOIs(fFlowEvent_TPC->GetNumberOfPOIs()+1);
                }
                if(fFlowEvent_VZERO) {
                    fFlowEvent_VZERO->InsertTrack(&flowTrack);
                    fFlowEvent_VZERO->SetNumberOfPOIs(fFlowEvent_VZERO->GetNumberOfPOIs()+1);
                }
            }
        }
        // then add the reference section only for the TPC reference case
        for(Int_t i(0); i < fRPs->GetEntriesFast(); i++) {
            AliVTrack* rp = static_cast<AliVTrack*>(fRPs->At(i));
            if(fRhoVn->PassesCuts(rp) && rp->Pt() >= .15 && rp->Pt() <= 5.) {
                AliFlowTrack flowTrack = AliFlowTrack(rp);
                flowTrack.SetForPOISelection(kFALSE);
                flowTrack.SetForRPSelection(kTRUE);
                if(fFlowEvent_TPC) {
                    fFlowEvent_TPC->SetNumberOfRPs(fFlowEvent_TPC->GetNumberOfRPs()+1);
                    fFlowEvent_TPC->InsertTrack(&flowTrack);
                }
            }
        }
    }
    if(fFlowEvent_TPC)      fFlowEvent_TPC->TagSubeventsInEta(-10, -1, 1, 10);
    if(fFlowEvent_VZERO)    fFlowEvent_VZERO->TagSubeventsInEta(-10, -1, 1, 10);
    return (nAcceptedJets < 1) ? kFALSE : kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlow::QCnDifferentialFlowVectors(Double_t* repn, Double_t* impn, Double_t *mp, Double_t *reqn, Double_t *imqn, Double_t* mq, Int_t n) 
{
    // get (for now) unweighted differential flow vectors
    // FIXME move (part of) this code to AliAnalysisTaskRhoVnModulation
    Int_t iPois(fPois->GetEntriesFast());
    if(fVParticleAnalysis) {
        for(Int_t i(0); i < iPois; i++) {
            for(Int_t ptBin(0); ptBin < fPtBins->GetSize()-1; ptBin++) {
                AliVTrack* poi = static_cast<AliVTrack*>(fPois->At(i));
                if(fRhoVn->PassesCuts(poi)) {       // inherit cuts from mother task
                    if(poi->Pt() >= fPtBins->At(ptBin) && poi->Pt() < fPtBins->At(ptBin+1)) {
                            // fill the flow vectors assuming that all poi's are in the rp selection (true by design)
                            repn[ptBin]+=TMath::Cos(((Double_t)n)*poi->Phi());
                            impn[ptBin]+=TMath::Sin(((Double_t)n)*poi->Phi());
                            mp[ptBin]++;
                            reqn[ptBin]+=TMath::Cos(((Double_t)n)*poi->Phi());
                            imqn[ptBin]+=TMath::Sin(((Double_t)n)*poi->Phi());
                            mq[ptBin]++;
                    }
                }
            }
        }
    } else {
        for(Int_t i(0); i < iPois; i++) {
            for(Int_t ptBin(0); ptBin < fPtBins->GetSize()-1; ptBin++) {
                AliEmcalJet* poi = static_cast<AliEmcalJet*>(fPois->At(i));
                AliLocalRhoParameter* localRho = static_cast<AliLocalRhoParameter*>(fLocalRho);
                if(!localRho) break;
                Double_t rho(localRho->GetLocalVal(poi->Phi(), fJetRadius, localRho->GetVal()));
                Double_t pt(poi->Pt() - poi->Area() * rho);
                if(fRhoVn->PassesCuts(poi)) {    
                    if(pt >= fPtBins->At(ptBin) && pt < fPtBins->At(ptBin+1)) {    
                            // fill the flow vectors assuming that all poi's are in the rp selection (true by design)  
                            repn[ptBin]+=TMath::Cos(((Double_t)n)*poi->Phi());
                            impn[ptBin]+=TMath::Sin(((Double_t)n)*poi->Phi());
                            mp[ptBin]++;        // qn isn't filled, no overlap between poi's and rp's
                    }
                }
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
