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
 * analysis task for jet flow preparation
 *
 * this task is part of the emcal jet framework and should be run in the emcaljet train
 * the following extensions to an accepted AliVEvent are expected:
 *      - (anti-kt) jets
 *      - background estimate rho
 *      - pico tracks
 *      aod's and esd's are handled transparently
 * the task will attempt to estimate a phi-dependent background density rho 
 * by fitting vn harmonics
 *
 * author: Redmer Alexander Bertens, Utrecht Univeristy, Utrecht, Netherlands
 * rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 
 */

#include <TStyle.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliVVertex.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>

#include <AliPicoTrack.h>
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>

#include "AliAnalysisTaskRhoVnModulation.h"


class AliAnalysisTaskRhoVnModulation;
using namespace std;

ClassImp(AliAnalysisTaskRhoVnModulation)

AliAnalysisTaskRhoVnModulation::AliAnalysisTaskRhoVnModulation() : AliAnalysisTaskEmcalJet("AliAnalysisTaskRhoVnModulation", kTRUE), 
    fDebug(0), fInitialized(0), fFillQAHistograms(kTRUE), fFitModulationType(kNoFit), fDetectorType(kTPC), fFitModulationOptions("Q"), fRunModeType(kGrid), fDataType(kESD), fRandom(0), fMappedRunNumber(0), fInCentralitySelection(-1), fFitModulation(0), fNameJetClones(0), fNamePicoTrackClones(0), fNameRho(0), fAbsVertexZ(10), fHistCentrality(0), fHistVertexz(0), fHistRunnumbersPhi(0), fHistRunnumbersEta(0), fMinDisanceRCtoLJ(0), fRandomConeRadius(0.4), fOutputList(0), fOutputListGood(0), fOutputListBad(0), fHistAnalysisSummary(0), fHistSwap(0), fProfVn(0), fHistPsi2(0), fHistPsi2Spread(0), fHistPsiVZEROA(0), fHistPsiVZEROC(0), fHistPsiTPC(0), 
   fHistRhoVsMult(0), fHistRhoVsCent(0), fHistRhoAVsMult(0), fHistRhoAVsCent(0) {
    for(Int_t i(0); i < 10; i++) {
        fHistPicoTrackPt[i] = 0;
        fHistPicoTrackPhi[i] = 0;
        fHistPicoTrackEta[i] = 0;
        fHistPicoCat1[i] = 0;
        fHistPicoCat2[i] = 0;
        fHistPicoCat3[i] = 0;
        /* fHistClusterPt[i] = 0; */
        /* fHistClusterPhi[i] = 0; */
        /* fHistClusterEta[i] = 0; */ 
        /* fHistClusterCorrPt[i] = 0; */
        /* fHistClusterCorrPhi[i] = 0; */
        /* fHistClusterCorrEta[i] = 0; */
        fHistRhoPackage[i] = 0;
        fHistRho[i] = 0;
        fHistRCPhiEta[i] = 0;
        fHistRhoVsRCPt[i] = 0;
        fHistRCPt[i] = 0;
        fHistDeltaPtRC[i] = 0;
        fHistRCPhiEtaExLJ[i] = 0;
        fHistRhoVsRCPtExLJ[i] = 0;
        fHistRCPtExLJ[i] = 0;
        fHistDeltaPtRCExLJ[i] = 0;
        fHistRCPhiEtaRand[i] = 0;
        fHistRhoVsRCPtRand[i] = 0;
        fHistRCPtRand[i] = 0;
        fHistDeltaPtRCRand[i] = 0;
        fHistJetPtRaw[i] = 0;
        fHistJetPt[i] = 0;
        fHistJetEtaPhi[i] = 0;
        fHistJetPtArea[i] = 0;
        fHistJetPtConstituents[i] = 0;
        fHistJetPtInPlaneVZEROA[i] = 0;
        fHistJetPtOutPlaneVZEROA[i] = 0;
        fHistJetPtMidPlaneVZEROA[i] = 0; 
        fHistJetPtInPlaneVZEROC[i] = 0;
        fHistJetPtOutPlaneVZEROC[i] = 0;
        fHistJetPtMidPlaneVZEROC[i] = 0;
        fHistJetPtInPlaneTPC[i] = 0;
        fHistJetPtOutPlaneTPC[i] = 0;
        fHistJetPtMidPlaneTPC[i] = 0;
        fHistJetPsiTPCPt[i] = 0;
        fHistJetPsiVZEROAPt[i] = 0;
        fHistJetPsiVZEROCPt[i] = 0;
        fHistDeltaPhiVZEROA[i] = 0;
        fHistDeltaPhiVZEROC[i] = 0;
        fHistDeltaPhiTPC[i] = 0;
    }
    // default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskRhoVnModulation::AliAnalysisTaskRhoVnModulation(const char* name, runModeType type) : AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(0), fInitialized(0), fFillQAHistograms(kTRUE), fFitModulationType(kNoFit), fDetectorType(kTPC), fFitModulationOptions("Q"), fRunModeType(type), fDataType(kESD), fRandom(0), fMappedRunNumber(0), fInCentralitySelection(-1), fFitModulation(0), fNameJetClones(0), fNamePicoTrackClones(0), fNameRho(0), fAbsVertexZ(10), fHistCentrality(0), fHistVertexz(0), fHistRunnumbersPhi(0), fHistRunnumbersEta(0), fMinDisanceRCtoLJ(0), fRandomConeRadius(0.4), fOutputList(0), fOutputListGood(0), fOutputListBad(0), fHistAnalysisSummary(0), fHistSwap(0), fProfVn(0), fHistPsi2(0), fHistPsi2Spread(0), fHistPsiVZEROA(0), fHistPsiVZEROC(0), fHistPsiTPC(0), 
   fHistRhoVsMult(0), fHistRhoVsCent(0), fHistRhoAVsMult(0), fHistRhoAVsCent(0) {
    for(Int_t i(0); i < 10; i++) {
        fHistPicoTrackPt[i] = 0;
        fHistPicoTrackPhi[i] = 0;
        fHistPicoTrackEta[i] = 0;
        fHistPicoCat1[i] = 0;
        fHistPicoCat2[i] = 0;
        fHistPicoCat3[i] = 0;
        /* fHistClusterPt[i] = 0; */
        /* fHistClusterPhi[i] = 0; */
        /* fHistClusterEta[i] = 0; */
        /* fHistClusterCorrPt[i] = 0; */
        /* fHistClusterCorrPhi[i] = 0; */
        /* fHistClusterCorrEta[i] = 0; */
        fHistRhoPackage[i] = 0;
        fHistRho[i] = 0;
        fHistRCPhiEta[i] = 0;
        fHistRhoVsRCPt[i] = 0;
        fHistRCPt[i] = 0;
        fHistDeltaPtRC[i] = 0;
        fHistRCPhiEtaExLJ[i] = 0;
        fHistRhoVsRCPtExLJ[i] = 0;
        fHistRCPtExLJ[i] = 0;
        fHistDeltaPtRCExLJ[i] = 0;
        fHistRCPhiEtaRand[i] = 0;
        fHistRhoVsRCPtRand[i] = 0;
        fHistRCPtRand[i] = 0;
        fHistDeltaPtRCRand[i] = 0;
        fHistJetPtRaw[i] = 0;
        fHistJetPt[i] = 0;
        fHistJetEtaPhi[i] = 0;
        fHistJetPtArea[i] = 0;
        fHistJetPtConstituents[i] = 0;
        fHistJetPtInPlaneVZEROA[i] = 0;
        fHistJetPtOutPlaneVZEROA[i] = 0;
        fHistJetPtMidPlaneVZEROA[i] = 0; 
        fHistJetPtInPlaneVZEROC[i] = 0;
        fHistJetPtOutPlaneVZEROC[i] = 0;
        fHistJetPtMidPlaneVZEROC[i] = 0;
        fHistJetPtInPlaneTPC[i] = 0;
        fHistJetPtOutPlaneTPC[i] = 0;
        fHistJetPtMidPlaneTPC[i] = 0;
        fHistJetPsiTPCPt[i] = 0;
        fHistJetPsiVZEROAPt[i] = 0;
        fHistJetPsiVZEROCPt[i] = 0;
        fHistDeltaPhiVZEROA[i] = 0;
        fHistDeltaPhiVZEROC[i] = 0;
        fHistDeltaPhiTPC[i] = 0;
    }
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    switch (fRunModeType) {
        case kLocal : {
            gStyle->SetOptFit(1);
            DefineOutput(2, TList::Class());
            DefineOutput(3, TList::Class());
        } break;
        default: fDebug = -1;   // suppress debug info explicitely when not running locally
    }
}
//_____________________________________________________________________________
AliAnalysisTaskRhoVnModulation::~AliAnalysisTaskRhoVnModulation()
{
    // destructor
    if(fOutputList)     delete fOutputList;
    if(fOutputListGood) delete fOutputListGood;
    if(fOutputListBad)  delete fOutputListBad;
    if(fFitModulation)  delete fFitModulation;
    if(fHistSwap)       delete fHistSwap;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskRhoVnModulation::InitializeAnalysis() 
{
    // initialize the anaysis
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(fMinDisanceRCtoLJ==0) fMinDisanceRCtoLJ = .5*fJetRadius;
    if(dynamic_cast<AliAODEvent*>(InputEvent())) fDataType = kAOD; // determine the datatype
    else if(dynamic_cast<AliESDEvent*>(InputEvent())) fDataType = kESD;
    if(!fRandom) fRandom = new TRandom3(0);  // get a randomized if one hasn't been user-supplied
    switch (fFitModulationType)  {
        case kNoFit : { SetModulationFit(new TF1("fix_kNoFit", "[0]", 0, TMath::TwoPi())); } break;
        case kV2 : {
            SetModulationFit(new TF1("fit_kV2", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
            fFitModulation->SetParameter(0, 0.);        // normalization
            fFitModulation->SetParameter(3, 0.2);       // v2
            fFitModulation->FixParameter(1, 1.);        // constant
            fFitModulation->FixParameter(2, 2.);        // constant
        } break;
        case kV3: {
            SetModulationFit(new TF1("fit_kV3", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi()));
            fFitModulation->SetParameter(0, 0.);        // normalization
            fFitModulation->SetParameter(3, 0.2);       // v3
            fFitModulation->FixParameter(1, 1.);        // constant
            fFitModulation->FixParameter(2, 3.);        // constant
        } break;
        case kCombined : {
             SetModulationFit(new TF1("fit_kCombined", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))"));
             fFitModulation->SetParameter(0, 0.);       // normalization
             fFitModulation->SetParameter(3, 0.2);      // v2
             fFitModulation->FixParameter(1, 1.);       // constant
             fFitModulation->FixParameter(2, 2.);       // constant
             fFitModulation->FixParameter(5, 3.);       // constant
             fFitModulation->SetParameter(7, 0.2);      // v3
        } break;
        default: break;
    }
    switch (fRunModeType) {
        case kGrid : { fFitModulationOptions += "N0"; } break;
        default : break;
    }
    return kTRUE;
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskRhoVnModulation::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c)
{
    // book a TH1F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fOutputList) return 0x0;
    if(c!=-1) name = Form("%s_%i", name, c);
    TH1F* histogram = new TH1F(name, name, bins, min, max);
    histogram->GetXaxis()->SetTitle(x);
    histogram->GetYaxis()->SetTitle("[counts]");
    histogram->Sumw2();
    fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskRhoVnModulation::BookTH2F(const char* name, const char* x, const char*y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c)
{
    // book a TH2F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fOutputList) return 0x0;
    if(c!=-1) name = Form("%s_%i", name, c);
    TH2F* histogram = new TH2F(name, name, binsx, minx, maxx, binsy, miny, maxy);
    histogram->GetXaxis()->SetTitle(x);
    histogram->GetYaxis()->SetTitle(y);
    histogram->Sumw2();
    fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::UserCreateOutputObjects()
{
    // create output objects
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    // global QA
    fHistCentrality =           BookTH1F("fHistCentrality", "centrality \%", 102, -2, 100);
    fHistVertexz =              BookTH1F("fHistVertexz", "vertex z (cm)", 100, -12, 12);

    // pico track kinematics
    for(Int_t i(0); i < 10; i++) { 
        fHistPicoTrackPt[i] =          BookTH1F("fHistPicoTrackPt", "p_{t} [GeV/c]", 100, 0, 50, i);
        fHistPicoTrackPhi[i] =         BookTH1F("fHistPicoTrackPhi", "#phi", 100, 0, TMath::TwoPi(), i);
        fHistPicoTrackEta[i] =         BookTH1F("fHistPicoTrackEta", "#eta", 100, -1, 1, i);
        if(fFillQAHistograms) {
            fHistPicoCat1[i] =             BookTH2F("fHistPicoCat1", "#eta", "#phi", 100, -1, 1, 100, 0, TMath::TwoPi(), i);
            fHistPicoCat2[i] =             BookTH2F("fHistPicoCat2", "#eta", "#phi", 100, -1, 1, 100, 0, TMath::TwoPi(), i);
            fHistPicoCat3[i] =             BookTH2F("fHistPicoCat3", "#eta", "#phi", 100, -1, 1, 100, 0, TMath::TwoPi(), i);
        }
        // emcal kinematics
        /* fHistClusterPt[i] =            BookTH1F("fHistClusterPt", "p_{t} [GeV/c]", 100, 0, 100, i); */
        /* fHistClusterPhi[i] =           BookTH1F("fHistClusterPhi", "#phi", 100, 0, TMath::TwoPi(), i); */
        /* fHistClusterEta[i] =           BookTH1F("fHistClusterEta", "#eta", 100, -5, 5); */

        // emcal kinematics after hadronic correction
        /* fHistClusterCorrPt[i] =        BookTH1F("fHistClusterCorrPt", "p_{t} [GeV/c]", 100, 0, 100, i); */
        /* fHistClusterCorrPhi[i] =       BookTH1F("fHistClusterCorrPhi", "#phi", 100, 0, TMath::TwoPi(), i); */
        /* fHistClusterCorrEta[i] =       BookTH1F("fHistClusterCorrEta", "#eta", 100, -5, 5, i); */
    }

    // event plane estimates and quality
    fHistPsi2 =                 new TProfile("fHistPsi2", "fHistPsi2", 3, 0, 3);
    fHistPsi2->Sumw2();
    fHistPsi2Spread =           new TProfile("fHistPsi2Spread", "fHistPsi2Spread", 3, 0, 3);
    fHistPsi2Spread->Sumw2();
    fHistPsi2->GetXaxis()->SetBinLabel(1, "<#Psi_{2, VZEROA}>");
    fHistPsi2->GetXaxis()->SetBinLabel(2, "<#Psi_{2, VZEROC}>");
    fHistPsi2->GetXaxis()->SetBinLabel(3, "<#Psi_{2, TPC}>");
    fHistPsi2Spread->GetXaxis()->SetBinLabel(1, "<#Psi_{2, VZEROA} - #Psi_{2, VZEROC}>");
    fHistPsi2Spread->GetXaxis()->SetBinLabel(2, "<#Psi_{2, VZEROC} - #Psi_{2, TPC}>");
    fHistPsi2Spread->GetXaxis()->SetBinLabel(3, "<#Psi_{2, VZEROC} - #Psi_{2, TPC}>");
    fOutputList->Add(fHistPsi2);
    fOutputList->Add(fHistPsi2Spread);
    fHistPsiVZEROA =            BookTH1F("fHistPsiVZEROA", "#Psi_{VZEROA}", 100, -.5*TMath::Pi(), .5*TMath::Pi());
    fHistPsiVZEROC =            BookTH1F("fHistPsiVZEROC", "#Psi_{VZEROC}", 100, -.5*TMath::Pi(), .5*TMath::Pi());
    fHistPsiTPC =               BookTH1F("fHistPsiTPC", "#Psi_{TPC}", 100, -.5*TMath::Pi(), .5*TMath::Pi());

    // background
    for(Int_t i(0); i < 10; i ++) {
        fHistRhoPackage[i] =           BookTH1F("fHistRhoPackage",  "#rho [GeV/c]", 100, 0, 150, i);
        fHistRho[i] =                  BookTH1F("fHistRho", "#rho [GeV/c]", 100, 0, 150, i);
    }
    fHistRhoVsMult =            BookTH2F("fHistRhoVsMult", "multiplicity", "#rho [GeV/c]", 100, 0, 4000, 100, 0, 250);
    fHistRhoVsCent =            BookTH2F("fHistRhoVsCent", "centrality", "#rho [GeV/c]", 100, 0, 100, 100, 0, 250);
    fHistRhoAVsMult =           BookTH2F("fHistRhoAVsMult", "multiplicity", "#rho * A (jet) [GeV/c]", 100, 0, 4000, 100, 0, 50);
    fHistRhoAVsCent =           BookTH2F("fHistRhoAVsCent", "centrality", "#rho * A (jet) [GeV/c]", 100, 0, 100, 100, 0, 50);

    // delta pt distributions
    for(Int_t i(0); i < 10; i ++) {
        fHistRCPhiEta[i] =             BookTH2F("fHistRCPhiEta", "#phi (RC)", "#eta (RC)", 100, 0, TMath::TwoPi(), 100, -1, 1, i);
        fHistRhoVsRCPt[i] =            BookTH2F("fHistRhoVsRCPt", "p_{t} (RC) [GeV/c]", "#rho * A (RC) [GeV/c]", 100, 0, 300, 100, 0, 250, i);
        fHistRCPt[i] =                 BookTH1F("fHistRCPt", "p_{t} (RC) [GeV/c]", 130, -20, 150, i);
        fHistDeltaPtRC[i] =            BookTH1F("fHistDeltaPtRC", "#delta p_{t} [GeV/c]", 180, -50, 150, i);
        fHistRCPhiEtaExLJ[i] =         BookTH2F("fHistRCPhiEtaExLJ", "#phi (RC)", "#eta (RC)", 100, 0, TMath::TwoPi(), 100, -1, 1, i);
        fHistRhoVsRCPtExLJ[i] =        BookTH2F("fHistRhoVsRCPtExLJ", "p_{t} (RC) [GeV/c]", "#rho * A (RC) [GeV/c]", 100, 0, 300, 100, 0, 250, i);
        fHistRCPtExLJ[i] =             BookTH1F("fHistRCPtExLJ", "p_{t} (RC) [GeV/c]", 130, -20, 150, i);
        fHistDeltaPtRCExLJ[i] =        BookTH1F("fHistDeltaPtRCExLJ", "#delta p_{t} [GeV/c]", 180, -50, 150, i);
        fHistRCPhiEtaRand[i] =         BookTH2F("fHistRCPhiEtaRand", "#phi (RC)", "#eta (RC)", 100, 0, TMath::TwoPi(), 100, -1, 1, i);
        fHistRhoVsRCPtRand[i] =        BookTH2F("fHistRhoVsRCPtRand", "p_{t} (RC) [GeV/c]", "#rho * A (RC) [GeV/c]", 100, 0, 300, 100, 0, 250, i);
        fHistRCPtRand[i] =             BookTH1F("fHistRCPtRand", "p_{t} (RC) [GeV/c]", 130, -20, 150, i);
        fHistDeltaPtRCRand[i] =        BookTH1F("fHistDeltaPtRCRand", "#delta p_{t} [GeV/c]", 180, -50, 150, i);

        // jet histograms (after kinematic cuts)
        fHistJetPtRaw[i] =             BookTH1F("fHistJetPtRaw", "p_{t} RAW [GeV/c]", 200, -50, 150, i);
        fHistJetPt[i] =                BookTH1F("fHistJetPt", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetEtaPhi[i] =            BookTH2F("fHistJetEtaPhi", "#eta", "#phi", 100, -1, 1, 100, 0, TMath::TwoPi(), i);
        fHistJetPtArea[i] =            BookTH2F("fHistJetPtArea", "p_{t} [GeV/c]", "Area", 200, -50, 150, 60, 0, 0.3, i);
        fHistJetPtConstituents[i] =    BookTH2F("fHistJetPtConstituents", "p_{t} [GeV/c]", "Area", 200, -50, 150, 60, 0, 150, i);

        // in plane and out of plane spectra
        fHistJetPtInPlaneVZEROA[i] =   BookTH1F("fHistJetPtInPlaneVZEROA", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtOutPlaneVZEROA[i] =  BookTH1F("fHistJetPtOutPlaneVZEROA", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtMidPlaneVZEROA[i] =  BookTH1F("fHistJetPtMidPlaneVZEROA", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtInPlaneVZEROC[i] =   BookTH1F("fHistJetPtInPlaneVZEROC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtOutPlaneVZEROC[i] =  BookTH1F("fHistJetPtOutPlaneVZEROC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtMidPlaneVZEROC[i] =  BookTH1F("fHistJetPtMidPlaneVZEROC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtInPlaneTPC[i] =      BookTH1F("fHistJetPtInPlaneTPC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtOutPlaneTPC[i] =     BookTH1F("fHistJetPtOutPlaneTPC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPtMidPlaneTPC[i] =     BookTH1F("fHistJetPtMidPlaneTPC", "p_{t} [GeV/c]", 200, -50, 150, i);
        fHistJetPsiTPCPt[i] =          BookTH2F("fHistJetPsiTPCPt", "#phi_{jet} - #Psi_{TPC}", "p_{t} [GeV/c]", 100, 0., TMath::TwoPi(), 100, -50, 100, i);
        fHistJetPsiVZEROAPt[i] =       BookTH2F("fHistJetPsiVZEROAPt", "#phi_{jet} - #Psi_{VZEROA}", "p_{t} [GeV/c]", 100, 0., TMath::TwoPi(), 100, -50, 100, i);
        fHistJetPsiVZEROCPt[i] =       BookTH2F("fHistJetPsiVZEROCPt", "#phi_{jet} - #Psi_{VZEROC}", "p_{t} [GeV/c]", 100, 0., TMath::TwoPi(), 100, -50, 100, i);

        // phi minus psi
        fHistDeltaPhiVZEROA[i] =       BookTH1F("fHistDeltaPhiVZEROA", "#phi_{jet} - #Psi_{VZEROA}", 100, 0, TMath::TwoPi(), i);
        fHistDeltaPhiVZEROC[i] =       BookTH1F("fHistDeltaPhiVZEROC", "#phi_{jet} - #Psi_{VZEROC}", 100, 0, TMath::TwoPi(), i);
        fHistDeltaPhiTPC[i] =          BookTH1F("fHistDeltaPhiTPC", "#phi_{jet} - #Psi_{TPC}", 100, 0, TMath::TwoPi(), i);
    }

    // analysis summary histrogram, saves all relevant analysis settigns
    fHistAnalysisSummary = BookTH1F("fHistAnalysisSummary", "flag", 37, -0.5, 37.5);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(1, "fjetRadius"); 
    fHistAnalysisSummary->SetBinContent(1, fJetRadius);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(2, "fPtBiasJetTrack");
    fHistAnalysisSummary->SetBinContent(2, fPtBiasJetTrack);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(3, "fPtBiasJetClus");
    fHistAnalysisSummary->SetBinContent(3, fPtBiasJetClus);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(4, "fJetPtCut");
    fHistAnalysisSummary->SetBinContent(4, fJetPtCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(5, "fJetAreaCut");
    fHistAnalysisSummary->SetBinContent(5, fJetAreaCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(6, "fPercAreaCut");
    fHistAnalysisSummary->SetBinContent(6, fPercAreaCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(7, "fAreaEmcCut");
    fHistAnalysisSummary->SetBinContent(7, fAreaEmcCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(8, "fJetMinEta");
    fHistAnalysisSummary->SetBinContent(8, fJetMinEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(9, "fJetMaxEta");
    fHistAnalysisSummary->SetBinContent(9, fJetMaxEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(10, "fJetMinPhi");
    fHistAnalysisSummary->SetBinContent(10, fJetMinPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(11, "fJetMaxPhi");
    fHistAnalysisSummary->SetBinContent(11, fJetMaxPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(12, "fMaxClusterPt");
    fHistAnalysisSummary->SetBinContent(12, fMaxClusterPt);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(13, "fMaxTrackPt");
    fHistAnalysisSummary->SetBinContent(13, fMaxTrackPt);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(14, "fLeadingHadronType");
    fHistAnalysisSummary->SetBinContent(14, fLeadingHadronType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(15, "fAnaType");
    fHistAnalysisSummary->SetBinContent(15, fAnaType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(16, "fForceBeamType");
    fHistAnalysisSummary->SetBinContent(16, fForceBeamType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(17, "fMinCent");
    fHistAnalysisSummary->SetBinContent(17, fMinCent);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(18, "fMaxCent");
    fHistAnalysisSummary->SetBinContent(18, fMaxCent);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(19, "fMinVz");
    fHistAnalysisSummary->SetBinContent(19, fMinVz);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(20, "fMaxVz");
    fHistAnalysisSummary->SetBinContent(20, fMaxVz);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(21, "fOffTrigger");
    fHistAnalysisSummary->SetBinContent(21, fOffTrigger);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(22, "fClusPtCut");
    fHistAnalysisSummary->SetBinContent(22, fClusPtCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(23, "fTrackPtCut");
    fHistAnalysisSummary->SetBinContent(23, fTrackPtCut);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(24, "fTrackMinEta");
    fHistAnalysisSummary->SetBinContent(24, fTrackMinEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(25, "fTrackMaxEta");
    fHistAnalysisSummary->SetBinContent(25, fTrackMaxEta);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(26, "fTrackMinPhi");
    fHistAnalysisSummary->SetBinContent(26, fTrackMinPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(27, "fTrackMaxPhi");
    fHistAnalysisSummary->SetBinContent(27, fTrackMaxPhi);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(28, "fClusTimeCutLow");
    fHistAnalysisSummary->SetBinContent(28, fClusTimeCutLow);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(29, "fClusTimeCutUp");
    fHistAnalysisSummary->SetBinContent(29, fClusTimeCutUp);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(30, "fMinPtTrackInEmcal");
    fHistAnalysisSummary->SetBinContent(30, fMinPtTrackInEmcal);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(31, "fEventPlaneVsEmcal");
    fHistAnalysisSummary->SetBinContent(31, fEventPlaneVsEmcal);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(32, "fMinEventPlane");
    fHistAnalysisSummary->SetBinContent(32, fMaxEventPlane);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(33, "fRandomConeRadius");
    fHistAnalysisSummary->SetBinContent(33, fRandomConeRadius);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(34, "fitModulationType");
    fHistAnalysisSummary->SetBinContent(34, (int)fFitModulationType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(35, "runModeType");
    fHistAnalysisSummary->SetBinContent(35, (int)fRunModeType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(36, "data type");
    fHistAnalysisSummary->SetBinContent(36, (int)fDataType);
    fHistAnalysisSummary->GetXaxis()->SetBinLabel(37, "iterator");
    fHistAnalysisSummary->SetBinContent(37, 1.);

    if(fFillQAHistograms) {
        fHistRunnumbersEta = new TH2F("fHistRunnumbersEta", "fHistRunnumbersEta", 100, -.5, 99.5, 100, -1.1, 1.1);
        fHistRunnumbersEta->Sumw2();
        fOutputList->Add(fHistRunnumbersEta);
        fHistRunnumbersPhi = new TH2F("fHistRunnumbersPhi", "fHistRunnumbersPhi", 100, -.5, 99.5, 100, -0.2, TMath::TwoPi()+0.2);
        fHistRunnumbersPhi->Sumw2();
        fOutputList->Add(fHistRunnumbersPhi);
    }

    fHistSwap = new TH1F("fHistSwap", "fHistSwap", 20, 0, TMath::TwoPi());
    fHistSwap->Sumw2();
    fProfVn = new TProfile("fProfVn", "fProfVn", 2, -0.5, 1.5);
    fProfVn->GetXaxis()->SetBinLabel(1, "v_{2}(EBYE)");
    fProfVn->GetXaxis()->SetBinLabel(2, "v_{2}(EBYE)");

    fOutputList->Add(fProfVn);
    PostData(1, fOutputList);

    switch (fRunModeType) {
        case kLocal : {
            fOutputListGood = new TList();
            fOutputListGood->SetOwner(kTRUE);
            fOutputListBad = new TList();
            fOutputListBad->SetOwner(kTRUE);
            PostData(2, fOutputListGood);
            PostData(3, fOutputListBad);
        } break;
        default: break;
    }


}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskRhoVnModulation::Run()
{
    // user exec: execute once for each event
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!fInitialized) fInitialized = InitializeAnalysis();
    // reject the event if expected data is missing
    if(!PassesCuts(InputEvent())) return kFALSE;
    if(!(fTracks||fJets||fRho)) return kFALSE;
    if(!fCaloClusters && fDebug > 0) printf(" > warning: couldn't retreive calo clusters! < \n");

    // [0][0] psi2a     [1,0]   psi2c
    // [0][1] psi3a     [1,1]   psi3c
    Double_t vzero[2][2];
    CalculateEventPlaneVZERO(vzero);
    // [0] psi2         [1] psi3
    Double_t tpc[2];
    CalculateEventPlaneTPC(tpc);
    
    // arrays which will hold the fit parameters
    Double_t fitParameters[] = {0,0,0,0,0,0,0,0,0};
    Double_t psi2(-1), psi3(-1);
    switch (fDetectorType) {    // determine the detector type for the rho fit
        case kTPC :     { psi2 = tpc[0];        psi3 = tpc[1]; }
        case kVZEROA :  { psi2 = vzero[0][0];   psi3 = vzero[0][1]; }
        case kVZEROC :  { psi2 = vzero[1][0];   psi3 = vzero[1][1]; }
        default : break;
    }

    switch (fFitModulationType) { // do the fits
        case kNoFit : { fFitModulation->FixParameter(0, RhoVal()); } break;
        case kV2 : {
            CorrectRho(fitParameters, psi2, psi3);
            fProfVn->Fill((double)0, fFitModulation->GetParameter(3));
        } break;
        case kV3 : {
            CorrectRho(fitParameters, psi2, psi3);
            fProfVn->Fill((double)1, fFitModulation->GetParameter(3));
        } break;
        case kCombined : {
            CorrectRho(fitParameters, psi2, psi3);
            fProfVn->Fill((double)0, fFitModulation->GetParameter(3));
            fProfVn->Fill((double)1, fFitModulation->GetParameter(7));
        } break;
        default : break;
    }
    // fill a number of histograms 
    FillHistogramsAfterSubtraction(vzero, tpc);

    // send the output to the connected output container
    PostData(1, fOutputList);
    switch (fRunModeType) {
        case kLocal : {
            PostData(2, fOutputListGood);
            PostData(3, fOutputListBad);
        } break;
        default: break;
    }
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::CalculateEventPlaneVZERO(Double_t vzero[2][2]) const 
{
    // grab the UNCALIBRATED vzero event plane
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Double_t qxa2(0), qya2(0), qxc2(0), qyc2(0);    // for psi2
    Double_t qxa3(0), qya3(0), qxc3(0), qyc3(0);    // for psi3
    for(Int_t iVZERO(0); iVZERO < 64; iVZERO++) {
        Double_t phi(TMath::PiOver4()*(.5+iVZERO%8)), /* eta(0), */ weight(InputEvent()->GetVZEROEqMultiplicity(iVZERO));
//        (iVZERO<32) ? eta = -3.45+.5*(iVZERO/8) : eta = 4.8-.6*((iVZERO/8)-4);
        if(iVZERO<32) {
            qxa2 += weight*TMath::Cos(2.*phi);
            qya2 += weight*TMath::Sin(2.*phi);
            qxa3 += weight*TMath::Cos(3.*phi);
            qya3 += weight*TMath::Sin(3.*phi);
        }
        else {
            qxc2 += weight*TMath::Cos(2.*phi);
            qyc2 += weight*TMath::Sin(2.*phi);
            qxc3 += weight*TMath::Cos(3.*phi);
            qyc3 += weight*TMath::Sin(3.*phi);
       }
    }
    vzero[0][0] = .5*TMath::ATan2(qya2, qxa2);
    vzero[1][0] = .5*TMath::ATan2(qyc2, qxc2);
    vzero[0][1] = (1./3.)*TMath::ATan2(qya3, qxa3);
    vzero[1][1] = (1./3.)*TMath::ATan2(qyc3, qxc3);
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::CalculateEventPlaneTPC(Double_t* tpc) const
{
   // grab the TPC event plane
   if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
   Double_t qx2(0), qy2(0);     // for psi2
   Double_t qx3(0), qy3(0);     // for psi3
   if(fTracks) {
       Int_t iTracks(fTracks->GetEntriesFast());
       for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
           AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
           if(!PassesCuts(track)) continue;
           qx2+= TMath::Cos(2.*track->Phi());
           qy2+= TMath::Sin(2.*track->Phi());
           qx3+= TMath::Cos(3.*track->Phi());
           qy3+= TMath::Sin(3.*track->Phi());
       }
   }
   tpc[0] = .5*TMath::ATan2(qy2, qx2);
   tpc[1] = (1./3.)*TMath::ATan2(qy3, qx3);
} 
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::CalculateRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, 
        AliEmcalJet* jet, Bool_t randomize) const
{
    // get a random cone
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    pt = 0; eta = 0; phi = 0;
    Float_t etaJet(999), phiJet(999), dJet(999);        // no jet: same as jet very far away
    if(jet) { // if a leading jet is given, use its kinematic properties
        etaJet = jet->Eta();
        phiJet = jet->Phi();
    }
    // force the random cones to at least be within detector acceptance
    Float_t minPhi(fJetMinPhi), maxPhi(fJetMaxPhi);
    if(maxPhi > TMath::TwoPi()) maxPhi = TMath::TwoPi();
    if(minPhi < 0 ) minPhi = 0;
    Float_t diffRcRJR(TMath::Abs(fRandomConeRadius-fJetRadius));
    // construct a random cone and see if it's far away enough from the leading jet
    Int_t attempts(1000);
    while(kTRUE) {
        attempts--;
        eta = gRandom->Uniform(fJetMinEta+diffRcRJR, fJetMaxEta-diffRcRJR);
        phi = gRandom->Uniform(minPhi, maxPhi);

        dJet = TMath::Sqrt((etaJet-eta)*(etaJet-eta)+(phiJet-phi)*(phiJet-phi));
        if(dJet > fMinDisanceRCtoLJ) break;
        else if (attempts == 0) {
            printf(" > No random cone after 1000 tries, giving up ... !\n");
            return;
        }
    }
    if(fTracks) {
        Int_t iTracks(fTracks->GetEntriesFast());
        for(Int_t i(0); i < iTracks; i++) {
            AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
            if(!PassesCuts(track)) continue;
            Float_t etaTrack(track->Eta()), phiTrack(track->Phi()), ptTrack(track->Pt());
            // if requested, randomize eta and phi to destroy any correlated fluctuations
            if(randomize) {
                etaTrack = gRandom->Uniform(fTrackMinEta, fTrackMaxEta);
                phiTrack = gRandom->Uniform(minPhi, maxPhi);
            }
            // get distance from cone
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi + TMath::TwoPi())) phiTrack+=TMath::TwoPi();
            if(TMath::Abs(phiTrack-phi) > TMath::Abs(phiTrack - phi - TMath::TwoPi())) phiTrack-=TMath::TwoPi();
            if(TMath::Sqrt(TMath::Abs((etaTrack-eta)*(etaTrack-eta)+(phiTrack-phi)*(phiTrack-phi))) <= fRandomConeRadius) pt+=ptTrack;
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::CorrectRho(Double_t* params, Double_t psi2, Double_t psi3) const
{
    // get rho' -> rho(phi)
    // the fit is constrained based on the switch of fFitModulationType which is called 
    // in the initialization of this class. 
    // after fitting, an array of fit parameters is set which can be re-created at
    // any point from the TF1 pointer fFitModulation
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    TString detector("");
    switch (fDetectorType) {
        case kTPC : detector+="TPC";
            break;
        case kVZEROA : detector+="VZEROA";
            break;
        case kVZEROC : detector+="VZEROC";
            break;
        default: break;
    }
    Int_t iTracks(fTracks->GetEntriesFast());
    fHistSwap->Reset();     // clear the histogram
    for(Int_t i(0); i < iTracks; i++) {
            AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
            if(!PassesCuts(track) || track->Pt() > 5 || track->Pt() < 0.15) continue;
            fHistSwap->Fill(track->Phi(), track->Pt());
    }
    fFitModulation->SetParameter(0, RhoVal());
    switch (fFitModulationType) {
        case kNoFit : { fFitModulation->FixParameter(0, RhoVal() ); }
        case kV2 : { fFitModulation->FixParameter(4, psi2); } break;
        case kV3 : { fFitModulation->FixParameter(4, psi3); } break;
        case kCombined : { 
            fFitModulation->FixParameter(4, psi2); 
            fFitModulation->FixParameter(6, psi3);
        } break;
        default : break;
    }
    fHistSwap->Fit(fFitModulation, fFitModulationOptions.Data(), "", 0, TMath::TwoPi());
    for(Int_t i(0); i < fFitModulation->GetNpar(); i++) params[i] = fFitModulation->GetParameter(i);
    // for LOCAL didactic purposes, save the  best and the worst fits
    // this routine can produce a lot of output histograms and will not work on GRID 
    // since the output will become unmergeable
    switch (fRunModeType) {
        case kGrid : break;
        case kLocal : {
            static Int_t didacticCounterBest(0);
            static Int_t didacticCounterWorst(0);
            static Double_t bestFitP(.05);      // threshold for significance
            static Double_t worstFitP(.05);
            Double_t p(ChiSquare(fFitModulation->GetNDF(), fFitModulation->GetChisquare()));
            if(p > bestFitP || p > 0.12) {
                TProfile* didacticProfile = (TProfile*)fHistSwap->Clone(Form("Fit_%i_p_%.3f_cen_%i_%s", didacticCounterBest, p, fInCentralitySelection, detector.Data()));
                fOutputListGood->Add(didacticProfile);
                didacticCounterBest++;
                bestFitP = p;
             }
             else if(p < worstFitP) { 
                TProfile* didacticProfile = (TProfile*)fHistSwap->Clone(Form("Fit_%i_p_%.3f_cen_%i_%s", didacticCounterWorst, p, fInCentralitySelection, detector.Data() ));
                fOutputListBad->Add(didacticProfile);
                didacticCounterWorst++;
                worstFitP = p;
             }
         } break;
         default : break;
    }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskRhoVnModulation::PassesCuts(AliVEvent* event)
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
    fCent = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
    if(fCent <= 0 || fCent >= 100 || TMath::Abs(fCent-InputEvent()->GetCentrality()->GetCentralityPercentile("TRK")) > 5.) return kFALSE;
    fInCentralitySelection = TMath::FloorNint(fCent/10.);
    if(fFillQAHistograms) FillQAHistograms(event);
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskRhoVnModulation::PassesCuts(const AliVCluster* cluster) const
{
    // cluster cuts
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!cluster) return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskRhoVnModulation::PassesCuts(const AliEmcalJet* jet) const
{
    // jet cuts
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    if(!jet) return kFALSE;
    if(TMath::Abs(jet->Eta()) > .5) return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillHistogramsAfterSubtraction(Double_t vzero[2][2], Double_t* tpc) const
{
    // fill histograms 
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    FillTrackHistograms();
    /* FillClusterHistograms(); */
    FillJetHistograms(vzero, tpc); 
    /* FillCorrectedClusterHistograms(); */
    FillEventPlaneHistograms(vzero, tpc);
    FillRhoHistograms();
    FillDeltaPtHistograms();
    FillDeltaPhiHistograms(vzero, tpc);
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillTrackHistograms() const
{
    // fill track histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t iTracks(fTracks->GetEntriesFast());
    for(Int_t i(0); i < iTracks; i++) {
        AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(i));
        if(!PassesCuts(track)) continue;
        fHistPicoTrackPt[fInCentralitySelection]->Fill(track->Pt());
        fHistPicoTrackEta[fInCentralitySelection]->Fill(track->Eta());
        fHistPicoTrackPhi[fInCentralitySelection]->Fill(track->Phi());
        if(fFillQAHistograms) FillQAHistograms(track);
    }
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillClusterHistograms() const
{
    // fill cluster histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    /* Int_t iClusters(fCaloClusters->GetEntriesFast());
    for(Int_t i(0); i < iClusters; i++) {
        AliVCluster* cluster = static_cast<AliVCluster*>(fCaloClusters->At(iClusters));
        if (!PassesCuts(cluster)) continue;
        TLorentzVector clusterLorentzVector;
        cluster->GetMomentum(clusterLorentzVector, const_cast<Double_t*>(fVertex));
        fHistClusterPt[fInCentralitySelection]->Fill(clusterLorentzVector.Pt());
        fHistClusterEta[fInCentralitySelection]->Fill(clusterLorentzVector.Eta());
        fHistClusterPhi[fInCentralitySelection]->Fill(clusterLorentzVector.Phi());
    }
    return; */
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillCorrectedClusterHistograms() const
{
    // fill clusters after hadronic correction FIXME implement
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillEventPlaneHistograms(Double_t vzero[2][2], Double_t* tpc) const
{
    // fill event plane histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistPsi2->Fill(0.5, vzero[0][0]);
    fHistPsi2->Fill(1.5, vzero[1][0]);
    fHistPsi2->Fill(2.5, tpc[0]);
    fHistPsiVZEROA->Fill(vzero[0][0]);
    fHistPsiVZEROC->Fill(vzero[1][0]);
    fHistPsiTPC->Fill(tpc[0]);
    fHistPsi2Spread->Fill(0.5, vzero[0][0]-vzero[1][0]);
    fHistPsi2Spread->Fill(1.5, vzero[0][0]-tpc[0]);
    fHistPsi2Spread->Fill(2.5, vzero[1][0]-tpc[0]);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillRhoHistograms() const
{
    // fill rho histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fHistRhoPackage[fInCentralitySelection]->Fill(RhoVal());    // save the rho estimate from the emcal jet package
    // get multiplicity FIXME inefficient
    Int_t iTracks(fTracks->GetEntriesFast()), mult(0), iJets(fJets->GetEntriesFast());
    for(Int_t i(0); i < iTracks; i ++) { if(PassesCuts(static_cast<AliVTrack*>(fTracks->At(i)))) mult++; }
    Double_t rho(RhoVal(TMath::Pi(), TMath::Pi(), fRho->GetVal()));
    fHistRho[fInCentralitySelection]->Fill(rho);
    fHistRhoVsMult->Fill(mult, rho);
    fHistRhoVsCent->Fill(fCent, rho);
    for(Int_t i(0); i < iJets; i++) {
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
        if(!PassesCuts(jet)) continue;
        fHistRhoAVsMult->Fill(mult, rho * jet->Area());
        fHistRhoAVsCent->Fill(fCent, rho * jet->Area());
    }
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillDeltaPtHistograms() const
{
    // fill delta pt histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    static Int_t sJets[9999] = {-1};
    GetSortedArray(sJets, fJets);
//    if(sJets[0] <= 0) return;
    Int_t i(0), maxCones(20);
    AliEmcalJet* leadingJet(0x0);
    do { // get the leading jet 
        leadingJet = static_cast<AliEmcalJet*>(fJets->At(sJets[i]));
        i++;
    }
    while (!PassesCuts(leadingJet)&&i<fJets->GetEntriesFast()); 
    if(!leadingJet && fDebug > 0) printf(" > failed to retrieve leading jet ! < \n");
    const Float_t areaRC = fRandomConeRadius*fRandomConeRadius*TMath::Pi();
    // we're retrieved the leading jet, now get a random cone
    for(i = 0; i < maxCones; i++) {
       Float_t pt(0), eta(0), phi(0);
       // get a random cone without constraints on leading jet position
       CalculateRandomCone(pt, eta, phi, 0x0);
       if(pt > 0) {
           fHistRCPhiEta[fInCentralitySelection]->Fill(phi, eta);
           fHistRhoVsRCPt[fInCentralitySelection]->Fill(pt, RhoVal(phi, fJetRadius, fRho->GetVal())*areaRC);
           fHistRCPt[fInCentralitySelection]->Fill(pt);
           fHistDeltaPtRC[fInCentralitySelection]->Fill(pt - areaRC*RhoVal(phi, fJetRadius, fRho->GetVal()));
       }
       // get a random cone excluding leading jet area
       CalculateRandomCone(pt, eta, phi, leadingJet);
       if(pt > 0) {
           fHistRCPhiEtaExLJ[fInCentralitySelection]->Fill(phi, eta);
           fHistRhoVsRCPtExLJ[fInCentralitySelection]->Fill(pt, RhoVal(phi, fJetRadius, fRho->GetVal())*areaRC);
           fHistRCPtExLJ[fInCentralitySelection]->Fill(pt);
           fHistDeltaPtRCExLJ[fInCentralitySelection]->Fill(pt - areaRC*RhoVal(phi, fJetRadius, fRho->GetVal()));
       }
       // get a random cone in an event with randomized phi and eta
       CalculateRandomCone(pt, eta, phi, 0x0, kTRUE);
       if( pt > 0) {
           fHistRCPhiEtaRand[fInCentralitySelection]->Fill(phi, eta);
           fHistRhoVsRCPtRand[fInCentralitySelection]->Fill(pt, RhoVal(phi, fJetRadius, fRho->GetVal())*areaRC);
           fHistRCPtRand[fInCentralitySelection]->Fill(pt);
           fHistDeltaPtRCRand[fInCentralitySelection]->Fill(pt - areaRC*RhoVal(phi, fJetRadius, fRho->GetVal()));
       }
    } 
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillJetHistograms(Double_t vzero[2][2], Double_t* tpc) const
{
    // fill jet histograms
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    Int_t iJets(fJets->GetEntriesFast());
    for(Int_t i(0); i < iJets; i++) {
        AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
        if(!PassesCuts(jet)) continue;
        Double_t pt(jet->Pt()), area(jet->Area()), eta(jet->Eta()), phi(jet->Phi());
        Double_t rho(RhoVal(phi, fJetRadius, fRho->GetVal()));
        fHistJetPtRaw[fInCentralitySelection]->Fill(pt);
        fHistJetPt[fInCentralitySelection]->Fill(pt-area*rho);
        fHistJetEtaPhi[fInCentralitySelection]->Fill(eta, phi);
        fHistJetPtArea[fInCentralitySelection]->Fill(pt-area*rho, area);
        Double_t dPhiA = PhaseShift(phi-vzero[0][0]);
        Double_t dPhiC = PhaseShift(phi-vzero[1][0]);
        Double_t dPhiTPC = PhaseShift(phi-tpc[0]);
        Double_t PiE = TMath::PiOver4()/2.;
        fHistJetPsiTPCPt[fInCentralitySelection]->Fill(dPhiTPC, pt-area*rho);
        fHistJetPsiVZEROAPt[fInCentralitySelection]->Fill(dPhiA, pt-area*rho);
        fHistJetPsiVZEROCPt[fInCentralitySelection]->Fill(dPhiC, pt-area*rho);
        if(dPhiA > 15.*PiE && dPhiA < PiE) fHistJetPtInPlaneVZEROA[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiA > 7.*PiE && dPhiA < 9.*PiE) fHistJetPtInPlaneVZEROA[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiA > 3.*PiE && dPhiA < 5.*PiE) fHistJetPtOutPlaneVZEROA[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiA > 11.*PiE && dPhiA < 13.*PiE) fHistJetPtOutPlaneVZEROA[fInCentralitySelection]->Fill(pt-area*rho);
        else fHistJetPtMidPlaneVZEROA[fInCentralitySelection]->Fill(pt-area*rho);
        if(dPhiC > 15.*PiE && dPhiC < PiE) fHistJetPtInPlaneVZEROC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiC > 7.*PiE && dPhiC < 9.*PiE) fHistJetPtInPlaneVZEROC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiC > 3.*PiE && dPhiC < 5.*PiE) fHistJetPtOutPlaneVZEROC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiC > 11.*PiE && dPhiC < 13.*PiE) fHistJetPtOutPlaneVZEROC[fInCentralitySelection]->Fill(pt-area*rho);
        else fHistJetPtMidPlaneVZEROC[fInCentralitySelection]->Fill(pt-area*rho);
        if(dPhiTPC > 15.*PiE && dPhiTPC < PiE) fHistJetPtInPlaneTPC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiTPC > 7.*PiE && dPhiTPC < 9.*PiE) fHistJetPtInPlaneTPC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiTPC > 3.*PiE && dPhiTPC < 5.*PiE) fHistJetPtOutPlaneTPC[fInCentralitySelection]->Fill(pt-area*rho);
        else if(dPhiTPC > 11.*PiE && dPhiTPC < 13.*PiE) fHistJetPtOutPlaneTPC[fInCentralitySelection]->Fill(pt-area*rho);
        else fHistJetPtMidPlaneTPC[fInCentralitySelection]->Fill(pt-area*rho);
        // last but not least, set the subtracted pt
        jet->SetPtSub(jet->PtSub(rho));
        fHistJetPtConstituents[fInCentralitySelection]->Fill(jet->PtSub(), jet->Nch());
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillDeltaPhiHistograms(Double_t vzero[2][2], Double_t* tpc) const
{
   // fill phi minus psi histograms
   if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
   if(fTracks) {
       Int_t iTracks(fTracks->GetEntriesFast());
       for(Int_t iTPC(0); iTPC < iTracks; iTPC++) {
           AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(iTPC));
           if(!PassesCuts(track)) continue;
           fHistDeltaPhiVZEROA[fInCentralitySelection]->Fill(PhaseShift(track->Phi()-vzero[0][0]));
           fHistDeltaPhiVZEROC[fInCentralitySelection]->Fill(PhaseShift(track->Phi()-vzero[1][0]));
           fHistDeltaPhiTPC[fInCentralitySelection]->Fill(PhaseShift(track->Phi()-tpc[0]));
       }
   }
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillQAHistograms(AliVTrack* vtrack) const
{
    // fill qa histograms for pico tracks
    if(!vtrack) return;
    AliPicoTrack* track = static_cast<AliPicoTrack*>(vtrack);
    fHistRunnumbersPhi->Fill(fMappedRunNumber, track->Phi());
    fHistRunnumbersEta->Fill(fMappedRunNumber, track->Eta());
    if(!track) return;
    Int_t type((int)(track->GetTrackType()));
    switch (type) {
        case 0:
           fHistPicoCat1[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        case 1:
           fHistPicoCat2[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        case 2:
           fHistPicoCat3[fInCentralitySelection]->Fill(track->Eta(), track->Phi()); 
           break;
        default: break;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::FillQAHistograms(AliVEvent* vevent) 
{
    // fill qa histograms for events
    if(!vevent) return;
    fHistVertexz->Fill(vevent->GetPrimaryVertex()->GetZ());
    fHistCentrality->Fill(fCent);
    Int_t runNumber(InputEvent()->GetRunNumber());
    Int_t runs[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, 169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
    for(fMappedRunNumber = 0; fMappedRunNumber < 65; fMappedRunNumber++) {
        if(runs[fMappedRunNumber]==runNumber) break;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskRhoVnModulation::Terminate(Option_t *)
{
    // terminate
    switch (fRunModeType) {
        case kLocal : {
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
        if(fFillQAHistograms) {
            Int_t runs[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, 169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
            for(Int_t i(0); i < 64; i++) { 
                fHistRunnumbersPhi->GetXaxis()->SetBinLabel(i+1, Form("%i", runs[i]));
                fHistRunnumbersEta->GetXaxis()->SetBinLabel(i+1, Form("%i", runs[i]));
            }
            fHistRunnumbersPhi->GetXaxis()->SetBinLabel(65, "undetermined");
            fHistRunnumbersEta->GetXaxis()->SetBinLabel(65, "undetermined");
        }
        AliAnalysisTaskRhoVnModulation::Dump();
        for(Int_t i(0); i < fHistAnalysisSummary->GetXaxis()->GetNbins(); i++) printf( " > flag: %s \t content %.2f \n", fHistAnalysisSummary->GetXaxis()->GetBinLabel(1+i), fHistAnalysisSummary->GetBinContent(1+i));
        } break;
        default : break;
    }
}
//_____________________________________________________________________________
