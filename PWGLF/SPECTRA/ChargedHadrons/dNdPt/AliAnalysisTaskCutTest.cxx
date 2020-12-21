#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TTree.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AlidNdPtHelper.h"

#include "AliAnalysisTaskCutTest.h"

ClassImp(AliAnalysisTaskCutTest)

    //________________________________________________________________________
    AliAnalysisTaskCutTest::AliAnalysisTaskCutTest(const char* name)
    : AliAnalysisTaskSE(name), fEsdTrackCuts(0), fAccCuts(0), fEventCuts(0),
      fESD(0), fMC(0), fUseCentrality(kFALSE), fCentralityMin(0),
      fCentralityMax(100), fDeadZoneWidth(2), fMaxZ(220), fTrackHistRefitTPC(0),
      fTrackHistChi2TPC(0), fTrackHistRowsTPC(0), fTrackHistFindableTPC(0),
      fTrackHistSharedTPC(0), fTrackHistHitsITS(0), fTrackHistRefitITS(0),
      fTrackHistChi2ITS(0), fTrackHistHitsSPD(0), fTrackHistDCAz(0),
      fTrackHistDCAxy(0), fTrackHistChi2TPCITS(0), fTrackHistGeoLengthTPC(0),
      fTrackHistGeoNcrTPC(0), fTrackHistGeoNclTPC(0), fTrackHistNclTPC(0),
      fEventHist(0), fTrackHist(0), fTrackHist2(0), fPtHist(0), fPtHist2(0),
      fMCTrackHist(0), fMCTrackHist2(0), fUseMCInfo(0), fOutputContainer(0) {
    // Constructor

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 writes into a TH1 container
    DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
void AliAnalysisTaskCutTest::ConnectInputData(Option_t*) {
    // Connect ESD or AOD here
    // Called once

    TTree* tree = dynamic_cast<TTree*>(GetInputData(0));
    if (!tree) {
        Printf("ERROR: Could not read chain from input slot 0");
    } else {
        // Disable all branches and enable only the needed ones
        // The next two lines are different when data produced as AliESDEvent is
        // read
        // tree->SetBranchStatus("*", kFALSE);
        // tree->SetBranchStatus("fTracks.*", kTRUE);
        // tree->SetBranchStatus("Tracks.*", kTRUE);

        AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(
            AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

        if (!esdH) {
            Printf("ERROR: Could not get ESDInputHandler");
        } else
            fESD = dynamic_cast<AliESDEvent*>(esdH->GetEvent());
    }
}

//________________________________________________________________________
void AliAnalysisTaskCutTest::UserCreateOutputObjects() {

    fPtHist = new TH2D("fPtHist", "1pt:sigma1pt", 200, -10., 10., 1000, 0., 1.);
    fPtHist->GetXaxis()->SetTitle("1/pt");
    fPtHist->GetYaxis()->SetTitle("#sigma(1/pt)");
    fPtHist->Sumw2();

    fPtHist2 =
        new TH2D("fPtHist2", "1pt:sigma1pt", 200, -1., 1., 10000, 0., 0.1);
    fPtHist2->GetXaxis()->SetTitle("1/pt");
    fPtHist2->GetYaxis()->SetTitle("#sigma(1/pt)");
    fPtHist2->Sumw2();

    const Int_t ptBin = 76;
    Double_t ptMin = 0;
    Double_t ptMax = 200;
    Double_t binsPt[77] = {
        0.0,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35,  0.4,   0.45,  0.5,
        0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,   0.95,  1.0,   1.1,
        1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,   2.0,   2.2,   2.4,
        2.6,  2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,   4.5,   5.0,   5.5,
        6.0,  6.5,  7.0,  8.0,  9.0,  10.0, 11.0, 12.0,  13.0,  14.0,  15.0,
        16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,  32.0,  34.0,  36.0,
        40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 200.0};

    const Int_t etaBin = 8;
    Double_t etaMin = -0.8;
    Double_t etaMax = 0.8;

    Int_t phiBin = 1;
    Double_t phiMin = 0;
    Double_t phiMax = 2. * TMath::Pi();

    const Int_t multBin = 1024 * 8;
    Double_t multMin = -0.5;
    Double_t multMax = multBin - 0.5;

    Int_t PtNbins = ptBin;

    Int_t EtaNbins = 8;
    Double_t binsEta[9] = {-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};

    Int_t MultNbins = 16;
    Double_t binsMult[17] = {-0.5,  0.5,    1.5,    2.5,    5.5,   10.5,
                             20.5,  30.5,   50.5,   100.5,  150.5, 250.5,
                             500.5, 1000.5, 2000.5, 4000.5, 8000.5};
    Int_t ZvNbins = 4;
    Double_t binsZv[5] = {
        -10., -5., 0., 5., 10.,
    };
    Int_t PhiNbins = phiBin;
    Double_t* binsPhi = new Double_t[PhiNbins + 1];
    for (Int_t i = 0; i <= PhiNbins; i++) {
        binsPhi[i] = (i * 2. * TMath::Pi()) / ((Double_t)PhiNbins);
    }
    Int_t Phi2Nbins = 144;
    Double_t* binsPhi2 = new Double_t[Phi2Nbins + 1];
    for (Int_t i = 0; i <= Phi2Nbins; i++) {
        binsPhi2[i] = (i * 2. * TMath::Pi()) / ((Double_t)Phi2Nbins);
    }
    Int_t SigmaPtNbins = 500;
    Double_t* binsSigmaPt = new Double_t[SigmaPtNbins + 1];
    for (Int_t i = 0; i <= SigmaPtNbins; i++) {
        binsSigmaPt[i] = ((Double_t)i) / ((Double_t)SigmaPtNbins);
    }
    Int_t DeltaPtNbins = 500;
    Double_t* binsDeltaPt = new Double_t[DeltaPtNbins + 1];
    for (Int_t i = 0; i <= DeltaPtNbins; i++) {
        binsDeltaPt[i] = -1. + 2. * ((Double_t)i) / ((Double_t)DeltaPtNbins);
    }
    Int_t ChargeNbins = 2;
    Double_t binsCharge[3] = {-1.5, 0, 1.5};
    Int_t DeltaSigmaPtNbins = 1000;
    Double_t* binsDeltaSigmaPt = new Double_t[DeltaSigmaPtNbins + 1];
    for (Int_t i = 0; i <= DeltaSigmaPtNbins; i++) {
        binsDeltaSigmaPt[i] =
            -50. + 100. * ((Double_t)i) / ((Double_t)DeltaSigmaPtNbins);
    }

    Int_t MultEventNbins = 1024 * 8;
    Double_t binsMultEvent[MultEventNbins + 1];
    for (int i = 0; i <= MultEventNbins; i++) {
        binsMultEvent[i] = -0.5 + i;
    }

    Int_t binTrackHistRefitTPC[5] = {etaBin, phiBin, ptBin, multBin, 2};
    Double_t minTrackHistRefitTPC[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistRefitTPC[5] = {etaMax, phiMax, ptMax, multMax, 1.5};

    Int_t binTrackHistChi2TPC[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistChi2TPC[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistChi2TPC[5] = {etaMax, phiMax, ptMax, multMax, 50.};

    Int_t binTrackHistRowsTPC[5] = {etaBin, phiBin, ptBin, multBin, 161};
    Double_t minTrackHistRowsTPC[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistRowsTPC[5] = {etaMax, phiMax, ptMax, multMax, 160.5};

    Int_t binTrackHistFindableTPC[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistFindableTPC[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistFindableTPC[5] = {etaMax, phiMax, ptMax, multMax, 5.};

    Int_t binTrackHistSharedTPC[5] = {etaBin, phiBin, ptBin, multBin, 200};
    Double_t minTrackHistSharedTPC[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistSharedTPC[5] = {etaMax, phiMax, ptMax, multMax, 2.};

    Int_t binTrackHistHitsITS[5] = {etaBin, phiBin, ptBin, multBin, 8};
    Double_t minTrackHistHitsITS[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistHitsITS[5] = {etaMax, phiMax, ptMax, multMax, 7.5};

    Int_t binTrackHistRefitITS[5] = {etaBin, phiBin, ptBin, multBin, 2};
    Double_t minTrackHistRefitITS[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistRefitITS[5] = {etaMax, phiMax, ptMax, multMax, 1.5};

    Int_t binTrackHistChi2ITS[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistChi2ITS[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistChi2ITS[5] = {etaMax, phiMax, ptMax, multMax, 100.};

    Int_t binTrackHistHitsSPD[5] = {etaBin, phiBin, ptBin, multBin, 3};
    Double_t minTrackHistHitsSPD[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistHitsSPD[5] = {etaMax, phiMax, ptMax, multMax, 2.5};

    Int_t binTrackHistDCAz[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistDCAz[5] = {etaMin, phiMin, ptMin, multMin, 0};
    Double_t maxTrackHistDCAz[5] = {etaMax, phiMax, ptMax, multMax, 50.};

    Int_t binTrackHistDCAxy[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistDCAxy[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistDCAxy[5] = {etaMax, phiMax, ptMax, multMax, 50.};

    Int_t binTrackHistChi2TPCITS[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistChi2TPCITS[5] = {etaMin, phiMin, ptMin, multMin, 0.};
    Double_t maxTrackHistChi2TPCITS[5] = {etaMax, phiMax, ptMax, multMax, 100.};

    Int_t binTrackHistGeoLengthTPC[5] = {etaBin, phiBin, ptBin, multBin, 201};
    Double_t minTrackHistGeoLengthTPC[5] = {etaMin, phiMin, ptMin, multMin,
                                            -0.5};
    Double_t maxTrackHistGeoLengthTPC[5] = {etaMax, phiMax, ptMax, multMax,
                                            200.5};

    Int_t binTrackHistGeoNcrTPC[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistGeoNcrTPC[5] = {etaMin, phiMin, ptMin, multMin, 0};
    Double_t maxTrackHistGeoNcrTPC[5] = {etaMax, phiMax, ptMax, multMax, 5};

    Int_t binTrackHistGeoNclTPC[5] = {etaBin, phiBin, ptBin, multBin, 500};
    Double_t minTrackHistGeoNclTPC[5] = {etaMin, phiMin, ptMin, multMin, 0};
    Double_t maxTrackHistGeoNclTPC[5] = {etaMax, phiMax, ptMax, multMax, 5};

    Int_t binTrackHistNclTPC[5] = {etaBin, phiBin, ptBin, multBin, 161};
    Double_t minTrackHistNclTPC[5] = {etaMin, phiMin, ptMin, multMin, -0.5};
    Double_t maxTrackHistNclTPC[5] = {etaMax, phiMax, ptMax, multMax, 160.5};

    fTrackHistRefitTPC = new THnSparseD(
        "fTrackHistRefitTPC", "eta:phi:pt:mult:chi2c", 5, binTrackHistRefitTPC,
        minTrackHistRefitTPC, maxTrackHistRefitTPC);
    fTrackHistRefitTPC->SetBinEdges(2, binsPt);
    fTrackHistRefitTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistRefitTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistRefitTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistRefitTPC->GetAxis(4)->SetTitle("RefitTPC");
    //   fTrackHistRefitTPC->Sumw2();

    fTrackHistChi2TPC = new THnSparseD(
        "fTrackHistChi2TPC", "eta:phi:pt:mult:chi2c", 5, binTrackHistChi2TPC,
        minTrackHistChi2TPC, maxTrackHistChi2TPC);
    fTrackHistChi2TPC->SetBinEdges(2, binsPt);
    fTrackHistChi2TPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistChi2TPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistChi2TPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistChi2TPC->GetAxis(4)->SetTitle("Chi2TPC");
    //   fTrackHistChi2TPC->Sumw2();

    fTrackHistRowsTPC = new THnSparseD(
        "fTrackHistRowsTPC", "eta:phi:pt:mult:chi2c", 5, binTrackHistRowsTPC,
        minTrackHistRowsTPC, maxTrackHistRowsTPC);
    fTrackHistRowsTPC->SetBinEdges(2, binsPt);
    fTrackHistRowsTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistRowsTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistRowsTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistRowsTPC->GetAxis(4)->SetTitle("RowsTPC");
    //   fTrackHistRowsTPC->Sumw2();

    fTrackHistFindableTPC =
        new THnSparseD("fTrackHistFindableTPC", "eta:phi:pt:mult:chi2c", 5,
                       binTrackHistFindableTPC, minTrackHistFindableTPC,
                       maxTrackHistFindableTPC);
    fTrackHistFindableTPC->SetBinEdges(2, binsPt);
    fTrackHistFindableTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistFindableTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistFindableTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistFindableTPC->GetAxis(4)->SetTitle("FindableTPC");
    //   fTrackHistFindableTPC->Sumw2();

    fTrackHistSharedTPC = new THnSparseD(
        "fTrackHistSharedTPC", "eta:phi:pt:mult:chi2c", 5,
        binTrackHistSharedTPC, minTrackHistSharedTPC, maxTrackHistSharedTPC);
    fTrackHistSharedTPC->SetBinEdges(2, binsPt);
    fTrackHistSharedTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistSharedTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistSharedTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistSharedTPC->GetAxis(4)->SetTitle("SharedTPC");
    //   fTrackHistSharedTPC->Sumw2();

    fTrackHistHitsITS = new THnSparseD(
        "fTrackHistHitsITS", "eta:phi:pt:mult:chi2c", 5, binTrackHistHitsITS,
        minTrackHistHitsITS, maxTrackHistHitsITS);
    fTrackHistHitsITS->SetBinEdges(2, binsPt);
    fTrackHistHitsITS->GetAxis(0)->SetTitle("#eta");
    fTrackHistHitsITS->GetAxis(1)->SetTitle("#phi");
    fTrackHistHitsITS->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistHitsITS->GetAxis(4)->SetTitle("HitsITS");
    //   fTrackHistHitsITS->Sumw2();

    fTrackHistRefitITS = new THnSparseD(
        "fTrackHistRefitITS", "eta:phi:pt:mult:chi2c", 5, binTrackHistRefitITS,
        minTrackHistRefitITS, maxTrackHistRefitITS);
    fTrackHistRefitITS->SetBinEdges(2, binsPt);
    fTrackHistRefitITS->GetAxis(0)->SetTitle("#eta");
    fTrackHistRefitITS->GetAxis(1)->SetTitle("#phi");
    fTrackHistRefitITS->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistRefitITS->GetAxis(4)->SetTitle("RefitITS");
    //   fTrackHistRefitITS->Sumw2();

    fTrackHistChi2ITS = new THnSparseD(
        "fTrackHistChi2ITS", "eta:phi:pt:mult:chi2c", 5, binTrackHistChi2ITS,
        minTrackHistChi2ITS, maxTrackHistChi2ITS);
    fTrackHistChi2ITS->SetBinEdges(2, binsPt);
    fTrackHistChi2ITS->GetAxis(0)->SetTitle("#eta");
    fTrackHistChi2ITS->GetAxis(1)->SetTitle("#phi");
    fTrackHistChi2ITS->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistChi2ITS->GetAxis(4)->SetTitle("Chi2ITS");
    //   fTrackHistChi2ITS->Sumw2();

    fTrackHistHitsSPD = new THnSparseD(
        "fTrackHistHitsSPD", "eta:phi:pt:mult:chi2c", 5, binTrackHistHitsSPD,
        minTrackHistHitsSPD, maxTrackHistHitsSPD);
    fTrackHistHitsSPD->SetBinEdges(2, binsPt);
    fTrackHistHitsSPD->GetAxis(0)->SetTitle("#eta");
    fTrackHistHitsSPD->GetAxis(1)->SetTitle("#phi");
    fTrackHistHitsSPD->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistHitsSPD->GetAxis(4)->SetTitle("HitsSPD");
    //   fTrackHistHitsSPD->Sumw2();

    fTrackHistDCAz =
        new THnSparseD("fTrackHistDCAz", "eta:phi:pt:mult:chi2c", 5,
                       binTrackHistDCAz, minTrackHistDCAz, maxTrackHistDCAz);
    fTrackHistDCAz->SetBinEdges(2, binsPt);
    fTrackHistDCAz->GetAxis(0)->SetTitle("#eta");
    fTrackHistDCAz->GetAxis(1)->SetTitle("#phi");
    fTrackHistDCAz->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistDCAz->GetAxis(4)->SetTitle("DCAz");
    //   fTrackHistDCAz->Sumw2();

    fTrackHistDCAxy =
        new THnSparseD("fTrackHistDCAxy", "eta:phi:pt:mult:chi2c", 5,
                       binTrackHistDCAxy, minTrackHistDCAxy, maxTrackHistDCAxy);
    fTrackHistDCAxy->SetBinEdges(2, binsPt);
    fTrackHistDCAxy->GetAxis(0)->SetTitle("#eta");
    fTrackHistDCAxy->GetAxis(1)->SetTitle("#phi");
    fTrackHistDCAxy->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistDCAxy->GetAxis(4)->SetTitle("DCAxy");
    //   fTrackHistDCAxy->Sumw2();

    fTrackHistChi2TPCITS = new THnSparseD(
        "fTrackHistChi2TPCITS", "eta:phi:pt:mult:chi2c", 5,
        binTrackHistChi2TPCITS, minTrackHistChi2TPCITS, maxTrackHistChi2TPCITS);
    fTrackHistChi2TPCITS->SetBinEdges(2, binsPt);
    fTrackHistChi2TPCITS->GetAxis(0)->SetTitle("#eta");
    fTrackHistChi2TPCITS->GetAxis(1)->SetTitle("#phi");
    fTrackHistChi2TPCITS->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistChi2TPCITS->GetAxis(4)->SetTitle("Chi2TPCITS");
    //   fTrackHistChi2TPCITS->Sumw2();

    fTrackHistGeoLengthTPC =
        new THnSparseD("fTrackHistGeoLengthTPC", "eta:phi:pt:mult:geo", 5,
                       binTrackHistGeoLengthTPC, minTrackHistGeoLengthTPC,
                       maxTrackHistGeoLengthTPC);
    fTrackHistGeoLengthTPC->SetBinEdges(2, binsPt);
    fTrackHistGeoLengthTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistGeoLengthTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistGeoLengthTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistGeoLengthTPC->GetAxis(4)->SetTitle("GeoLengthTPC");
    //   fTrackHistGeoLengthTPC->Sumw2();

    fTrackHistGeoNcrTPC = new THnSparseD(
        "fTrackHistGeoNcrTPC", "eta:phi:pt:mult:geo", 5, binTrackHistGeoNcrTPC,
        minTrackHistGeoNcrTPC, maxTrackHistGeoNcrTPC);
    fTrackHistGeoNcrTPC->SetBinEdges(2, binsPt);
    fTrackHistGeoNcrTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistGeoNcrTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistGeoNcrTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistGeoNcrTPC->GetAxis(4)->SetTitle("GeoNcrTPC");
    //   fTrackHistGeoNcrTPC->Sumw2();

    fTrackHistGeoNclTPC = new THnSparseD(
        "fTrackHistGeoNclTPC", "eta:phi:pt:mult:geo", 5, binTrackHistGeoNclTPC,
        minTrackHistGeoNclTPC, maxTrackHistGeoNclTPC);
    fTrackHistGeoNclTPC->SetBinEdges(2, binsPt);
    fTrackHistGeoNclTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistGeoNclTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistGeoNclTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistGeoNclTPC->GetAxis(4)->SetTitle("GeoNclTPC");
    //   fTrackHistGeoNclTPC->Sumw2();

    fTrackHistNclTPC = new THnSparseD("fTrackHistNclTPC", "eta:phi:pt:mult:geo",
                                      5, binTrackHistNclTPC, minTrackHistNclTPC,
                                      maxTrackHistNclTPC);
    fTrackHistNclTPC->SetBinEdges(2, binsPt);
    fTrackHistNclTPC->GetAxis(0)->SetTitle("#eta");
    fTrackHistNclTPC->GetAxis(1)->SetTitle("#phi");
    fTrackHistNclTPC->GetAxis(2)->SetTitle("p_{T}");
    fTrackHistNclTPC->GetAxis(4)->SetTitle("NclTPC");
    //   fTrackHistNclTPC->Sumw2();

    Int_t binsEventHist[4] = {ZvNbins, MultEventNbins, MultEventNbins,
                              MultEventNbins};
    fEventHist = new THnSparseD("fEventHist", "zV:multMB:multRec:multTrue", 4,
                                binsEventHist);
    fEventHist->SetBinEdges(0, binsZv);
    fEventHist->SetBinEdges(1, binsMultEvent);
    fEventHist->SetBinEdges(2, binsMultEvent);
    fEventHist->SetBinEdges(3, binsMultEvent);
    fEventHist->GetAxis(0)->SetTitle("Zv (cm)");
    fEventHist->GetAxis(1)->SetTitle("multiplicity MB");
    fEventHist->GetAxis(2)->SetTitle("multiplicity Rec");
    fEventHist->GetAxis(3)->SetTitle("multiplicity True (MC)");
    //   fEventHist->Sumw2();

    Int_t binsTrackHist2[5] = {ZvNbins, PtNbins, EtaNbins, PhiNbins,
                               ChargeNbins};
    fTrackHist2 = new THnSparseD("fTrackHist2", "zV:multMB:multRec:multTrue", 5,
                                 binsTrackHist2);
    fTrackHist2->SetBinEdges(0, binsZv);
    fTrackHist2->SetBinEdges(1, binsPt);
    fTrackHist2->SetBinEdges(2, binsEta);
    fTrackHist2->SetBinEdges(3, binsPhi);
    fTrackHist2->SetBinEdges(4, binsCharge);
    fTrackHist2->GetAxis(0)->SetTitle("Zv (cm)");
    fTrackHist2->GetAxis(1)->SetTitle("pT");
    fTrackHist2->GetAxis(2)->SetTitle("#eta");
    fTrackHist2->GetAxis(3)->SetTitle("#phi");
    fTrackHist2->GetAxis(4)->SetTitle("charge");
    //   fTrackHist2->Sumw2();

    Int_t binsTrackHist[8] = {ZvNbins,   PtNbins,     EtaNbins,
                              PhiNbins,  ChargeNbins, SigmaPtNbins,
                              MultNbins, DeltaPtNbins};
    fTrackHist = new THnSparseD("fTrackHist",
                                "zV:pt:eta:phi:charge:sigmaPt:multMB:deltaPt",
                                8, binsTrackHist);
    fTrackHist->SetBinEdges(0, binsZv);
    fTrackHist->SetBinEdges(1, binsPt);
    fTrackHist->SetBinEdges(2, binsEta);
    fTrackHist->SetBinEdges(3, binsPhi);
    fTrackHist->SetBinEdges(4, binsCharge);
    fTrackHist->SetBinEdges(5, binsSigmaPt);
    fTrackHist->SetBinEdges(6, binsMult);
    fTrackHist->SetBinEdges(7, binsDeltaPt);
    fTrackHist->GetAxis(0)->SetTitle("Zv (cm)");
    fTrackHist->GetAxis(1)->SetTitle("pT");
    fTrackHist->GetAxis(2)->SetTitle("#eta");
    fTrackHist->GetAxis(3)->SetTitle("#phi");
    fTrackHist->GetAxis(4)->SetTitle("charge");
    fTrackHist->GetAxis(5)->SetTitle("sigmaPt");
    fTrackHist->GetAxis(6)->SetTitle("multiplicity MB");
    fTrackHist->GetAxis(7)->SetTitle("deltaPt");
    //   fTrackHist->Sumw2();

    Int_t binsMCTrackHist[8] = {ZvNbins,   PtNbins,     EtaNbins,
                                PhiNbins,  ChargeNbins, SigmaPtNbins,
                                MultNbins, DeltaPtNbins};
    fMCTrackHist = new THnSparseD(
        "fMCTrackHist", "mczV:mcpt:mceta:mcphi:charge:sigmaPt:multMB:deltaPt",
        8, binsMCTrackHist);
    fMCTrackHist->SetBinEdges(0, binsZv);
    fMCTrackHist->SetBinEdges(1, binsPt);
    fMCTrackHist->SetBinEdges(2, binsEta);
    fMCTrackHist->SetBinEdges(3, binsPhi);
    fMCTrackHist->SetBinEdges(4, binsCharge);
    fMCTrackHist->SetBinEdges(5, binsSigmaPt);
    fMCTrackHist->SetBinEdges(6, binsMult);
    fMCTrackHist->SetBinEdges(7, binsDeltaPt);
    fMCTrackHist->GetAxis(0)->SetTitle("Zv (cm)");
    fMCTrackHist->GetAxis(1)->SetTitle("pT");
    fMCTrackHist->GetAxis(2)->SetTitle("#eta");
    fMCTrackHist->GetAxis(3)->SetTitle("#phi");
    fMCTrackHist->GetAxis(4)->SetTitle("charge");
    fMCTrackHist->GetAxis(5)->SetTitle("sigmaPt");
    fMCTrackHist->GetAxis(6)->SetTitle("multiplicity MB");
    fMCTrackHist->GetAxis(7)->SetTitle("deltaPt");
    //   fMCTrackHist->Sumw2();

    Int_t binsMCTrackHist2[8] = {PtNbins, PtNbins, DeltaSigmaPtNbins};
    fMCTrackHist2 = new THnSparseD(
        "fMCTrackHist2", "pt:mcpt:deltaPtOverSigmaPt", 8, binsMCTrackHist2);
    fMCTrackHist2->SetBinEdges(0, binsPt);
    fMCTrackHist2->SetBinEdges(1, binsPt);
    fMCTrackHist2->SetBinEdges(2, binsDeltaSigmaPt);
    fMCTrackHist2->GetAxis(0)->SetTitle("p_{T}");
    fMCTrackHist2->GetAxis(1)->SetTitle("p_{T,MC}");
    fMCTrackHist2->GetAxis(2)->SetTitle("#Delta(1/pT)/#sigma(1/pT)");
    //   fMCTrackHist2->Sumw2();

    fOutputContainer = new TObjArray(13);
    fOutputContainer->SetName(GetName());
    fOutputContainer->SetOwner(kTRUE);
    fOutputContainer->Add(fTrackHistRefitTPC);
    fOutputContainer->Add(fTrackHistChi2TPC);
    fOutputContainer->Add(fTrackHistRowsTPC);
    fOutputContainer->Add(fTrackHistFindableTPC);
    fOutputContainer->Add(fTrackHistSharedTPC);
    fOutputContainer->Add(fTrackHistHitsITS);
    fOutputContainer->Add(fTrackHistRefitITS);
    fOutputContainer->Add(fTrackHistChi2ITS);
    fOutputContainer->Add(fTrackHistHitsSPD);
    fOutputContainer->Add(fTrackHistDCAz);
    fOutputContainer->Add(fTrackHistDCAxy);
    fOutputContainer->Add(fTrackHistChi2TPCITS);
    fOutputContainer->Add(fTrackHistGeoLengthTPC);
    fOutputContainer->Add(fTrackHistGeoNcrTPC);
    fOutputContainer->Add(fTrackHistGeoNclTPC);
    fOutputContainer->Add(fTrackHistNclTPC);
    fOutputContainer->Add(fEventHist);
    fOutputContainer->Add(fTrackHist);
    fOutputContainer->Add(fTrackHist2);
    fOutputContainer->Add(fPtHist);
    fOutputContainer->Add(fPtHist2);
    fOutputContainer->Add(fMCTrackHist);
    fOutputContainer->Add(fMCTrackHist2);

    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskCutTest::UserExec(Option_t*) {
    // Main loop
    // Called for each event

    if (!fESD) {
        Printf("ERROR: fESD not available");
        return;
    }

    // get reconstructed vertex
    const AliESDVertex* vtxESD = 0;
    vtxESD = fESD->GetPrimaryVertexTracks();

    if (!vtxESD) {
        Printf("ERROR: vtxESD not avaiable");
        return;
    }

    if (!fEventCuts->AcceptEvent(fESD)) {
        return;
    }

    if (fUseMCInfo) {
        fMC = MCEvent();
        if (!fMC) {
            Printf("ERROR: MC event not available");
            return;
        }
    }

    Int_t multMBTracks = 0;
    if (vtxESD->GetStatus()) {
        multMBTracks = vtxESD->GetNContributors();
    } else {
        return;
    }

    if (fUseCentrality) {
        AliCentrality* centrality = fESD->GetCentrality();
        Double_t centralityF = centrality->GetCentralityPercentile("V0M");
        if (centralityF > fCentralityMax)
            return;
        if (centralityF < fCentralityMin)
            return;
    }

    // here starts the cut variable stuff
    // track loop
    Int_t nTracks = 0;

    for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
        AliESDtrack* track = fESD->GetTrack(iTracks);
        if (!track) {
            Printf("ERROR: Could not receive track %d", iTracks);
            continue;
        }

        Double_t eta = track->Eta();
        Double_t phi = track->Phi();
        Double_t pt = track->Pt();
        Double_t mult = fESD->GetNumberOfTracks();

        Bool_t acceptRefitTPC;
        Bool_t acceptChi2TPC;
        Bool_t acceptRowsTPC;
        Bool_t acceptFindableTPC;
        Bool_t acceptSharedTPC;
        Bool_t acceptRefitITS;
        Bool_t acceptChi2ITS;
        Bool_t acceptHitsSPD;
        Bool_t acceptDCAz;
        Bool_t acceptDCAxy;
        Bool_t acceptChi2TPCITS;
        Bool_t acceptGeoLengthTPC;
        Bool_t acceptGeoNclTPC;
        Bool_t acceptGeoNcrTPC;
        Double_t vRefitTPC;
        Double_t vChi2TPC;
        Double_t vRowsTPC;
        Double_t vFindableTPC;
        Double_t vSharedTPC;
        Double_t vHitsITS;
        Double_t vRefitITS;
        Double_t vChi2ITS;
        Double_t vHitsSPD;
        Double_t vDCAz;
        Double_t vDCAxy;
        Double_t vChi2TPCITS;
        Double_t vGeoLengthTPC;
        Double_t vGeoNclTPC;
        Double_t vGeoNcrTPC;
        Double_t vNclTPC;

        UInt_t status = track->GetStatus();
        Double_t nClustersTPC = track->GetTPCclusters(0);
        Double_t nClustersTPCShared = track->GetTPCnclsS();
        Double_t nClusterTPCfindable = track->GetTPCNclsF();
        Double_t nCrossedRowsTPC = track->GetTPCCrossedRows();

        Float_t b[2];
        Float_t bCov[3];
        track->GetImpactParameters(b, bCov);

        // fill values
        vRefitTPC = ((status & AliESDtrack::kTPCrefit) != 0);
        vChi2TPC = (nClustersTPC > 0) ? track->GetTPCchi2() / nClustersTPC : -1;
        vRowsTPC = nCrossedRowsTPC;
        vFindableTPC =
            (nClusterTPCfindable > 0) ? vRowsTPC / nClusterTPCfindable : 1;
        vSharedTPC = nClustersTPCShared / nClustersTPC;
        vHitsITS = track->GetITSclusters(0);
        vRefitITS = ((status & AliESDtrack::kITSrefit) != 0);
        vChi2ITS = (vHitsITS > 0) ? (track->GetITSchi2() / vHitsITS) : -1;
        vHitsSPD = track->HasPointOnITSLayer(0) + track->HasPointOnITSLayer(1);
        vDCAz = TMath::Abs(b[1]);
        vDCAxy = TMath::Abs(b[0]) / (0.0182 + 0.0350 / TMath::Power(pt, 1.01));
        vChi2TPCITS = -2;   // done later
        vGeoLengthTPC = -2; // done later
        vGeoNcrTPC = -2;    // done later
        vGeoNclTPC = -2;    // done later
        vNclTPC = nClustersTPC;

        Int_t nCut = 0;

        acceptRefitTPC = (vRefitTPC > 0);
        acceptChi2TPC = (vChi2TPC <= 4);
        acceptRowsTPC = (vRowsTPC >= 120);
        acceptFindableTPC = (vFindableTPC >= 0.8);
        acceptSharedTPC = (vSharedTPC <= 0.4);
        acceptRefitITS = (vRefitITS > 0);
        acceptChi2ITS = (vChi2ITS <= 36);
        acceptHitsSPD = (vHitsSPD > 0);
        acceptDCAz = (vDCAz <= 2);
        acceptDCAxy = (vDCAxy <= 7);

        if (!acceptRefitTPC)
            nCut++;
        if (!acceptChi2TPC)
            nCut++;
        if (!acceptRowsTPC)
            nCut++;
        if (!acceptFindableTPC)
            nCut++;
        if (!acceptSharedTPC)
            nCut++;
        if (!acceptRefitITS)
            nCut++;
        if (!acceptChi2ITS)
            nCut++;
        if (!acceptHitsSPD)
            nCut++;
        if (!acceptDCAz)
            nCut++;
        if (!acceptDCAxy)
            nCut++;

        if (nCut > 1)
            continue; // if the track is removed by more than one cut already
                      // don't calculate everything.

        //
        vGeoLengthTPC = track->GetLengthInActiveZone(
            1, fDeadZoneWidth, fMaxZ, track->GetESDEvent()->GetMagneticField());
        acceptGeoLengthTPC =
            (vGeoLengthTPC >= (130.0 - TMath::Power(1. / pt, -1.5)));

        vGeoNcrTPC = nCrossedRowsTPC / vGeoLengthTPC;
        acceptGeoNcrTPC = (vGeoNcrTPC > 0.8);

        vGeoNclTPC = nClustersTPC / vGeoLengthTPC;
        acceptGeoNclTPC = (vGeoNclTPC > 0.7);

        if (!acceptGeoLengthTPC)
            nCut++;
        if (!acceptGeoNcrTPC)
            nCut++;
        if (!acceptGeoNclTPC)
            nCut++;

        if (nCut > 1)
            continue; // if the track is removed by more than one cut already
                      // don't calculate everything.

        // calculate chi2 tpc its contrained
        const AliESDVertex* vertex = 0;
        vertex = fESD->GetPrimaryVertexTracks();
        if ((!vertex || !vertex->GetStatus())) {
            vertex = fESD->GetPrimaryVertexSPD();
        }
        if (vertex->GetStatus()) {
            vChi2TPCITS = track->GetChi2TPCConstrainedVsGlobal(vertex);
        }

        acceptChi2TPCITS = (vChi2TPCITS <= 36);
        if (!acceptChi2TPCITS)
            nCut++;

        if (nCut > 1)
            continue;

        Double_t vTrack[5] = {eta, phi, pt, mult, -1};

        vTrack[4] = vRefitTPC;
        if (nCut == 0 || !acceptRefitTPC)
            fTrackHistRefitTPC->Fill(vTrack);
        vTrack[4] = vChi2TPC;
        if (nCut == 0 || !acceptChi2TPC)
            fTrackHistChi2TPC->Fill(vTrack);
        vTrack[4] = vRowsTPC;
        if (nCut == 0 || !acceptRowsTPC)
            fTrackHistRowsTPC->Fill(vTrack);
        vTrack[4] = vFindableTPC;
        if (nCut == 0 || !acceptFindableTPC)
            fTrackHistFindableTPC->Fill(vTrack);
        vTrack[4] = vSharedTPC;
        if (nCut == 0 || !acceptSharedTPC)
            fTrackHistSharedTPC->Fill(vTrack);
        vTrack[4] = vHitsITS;
        if (nCut == 0 || !acceptRefitITS)
            fTrackHistHitsITS->Fill(vTrack);
        vTrack[4] = vRefitITS;
        if (nCut == 0 || !acceptRefitITS)
            fTrackHistRefitITS->Fill(vTrack);
        vTrack[4] = vChi2ITS;
        if (nCut == 0 || !acceptChi2ITS)
            fTrackHistChi2ITS->Fill(vTrack);
        vTrack[4] = vHitsSPD;
        if (nCut == 0 || !acceptHitsSPD)
            fTrackHistHitsSPD->Fill(vTrack);
        vTrack[4] = vDCAz;
        if (nCut == 0 || !acceptDCAz)
            fTrackHistDCAz->Fill(vTrack);
        vTrack[4] = vDCAxy;
        if (nCut == 0 || !acceptDCAxy)
            fTrackHistDCAxy->Fill(vTrack);
        vTrack[4] = vChi2TPCITS;
        if (nCut == 0 || !acceptChi2TPCITS)
            fTrackHistChi2TPCITS->Fill(vTrack);
        vTrack[4] = vGeoLengthTPC;
        if (nCut == 0 || !acceptGeoLengthTPC)
            fTrackHistGeoLengthTPC->Fill(vTrack);
        vTrack[4] = vGeoNcrTPC;
        if (nCut == 0 || !acceptGeoNcrTPC)
            fTrackHistGeoNcrTPC->Fill(vTrack);
        vTrack[4] = vGeoNclTPC;
        if (nCut == 0 || !acceptGeoNclTPC)
            fTrackHistGeoNclTPC->Fill(vTrack);
        vTrack[4] = vNclTPC;
        fTrackHistNclTPC->Fill(vTrack);

        nTracks++;

    } // track loop

    TObjArray* allChargedTracks = 0;
    // Int_t multAll=0, multAcc=0, multRec=0;
    Int_t multAll = 0, multRec = 0;
    Int_t *labelsAll = 0, *labelsAcc = 0, *labelsRec = 0;

    allChargedTracks =
        AlidNdPtHelper::GetAllChargedTracks(fESD, AlidNdPtHelper::kTPCITS);
    Int_t entries = allChargedTracks->GetEntries();

    // calculate mult of reconstructed tracks
    Int_t multRecTemp = 0;
    for (Int_t i = 0; i < entries; ++i) {
        AliESDtrack* track = (AliESDtrack*)allChargedTracks->At(i);
        if (!track)
            continue;
        if (track->Charge() == 0)
            continue;
        if (fEsdTrackCuts->AcceptTrack(track)) {
            if (fAccCuts->AcceptTrack(track))
                multRecTemp++;
        }
    }

    AliHeader* header = 0;
    AliGenEventHeader* genHeader = 0;
    AliStack* stack = 0;
    TArrayF vtxMC(3);
    AliPWG0Helper::MCProcessType evtType = AliPWG0Helper::kInvalidProcess;

    Int_t multMCTrueTracks = 0;
    if (IsUseMCInfo()) {
        //
        if (!fMC) {
            AliDebug(AliLog::kError, "mcEvent not available");
            return;
        }
        // get MC event header
        header = fMC->Header();
        if (!header) {
            AliDebug(AliLog::kError, "Header not available");
            return;
        }
        // MC particle stack
        stack = fMC->Stack();
        if (!stack) {
            AliDebug(AliLog::kError, "Stack not available");
            return;
        }
        // get event type (ND=0x1, DD=0x2, SD=0x4)
        evtType = AlidNdPtHelper::GetEventProcessTypePA(header);
        AliDebug(AliLog::kDebug + 1, Form("Found process type %d", evtType));

        // get MC vertex
        genHeader = header->GenEventHeader();
        if (!genHeader) {
            AliDebug(AliLog::kError,
                     "Could not retrieve genHeader from Header");
            return;
        }
        genHeader->PrimaryVertex(vtxMC);

        // multipliticy of all MC primary tracks
        // in Zv, pt and eta ranges)
        multMCTrueTracks =
            AlidNdPtHelper::GetMCTrueTrackMult(fMC, fEventCuts, fAccCuts);
    } // end bUseMC

    // FillEventHistogram
    Double_t valuesEvent[4] = {vtxESD->GetZ(),
                               static_cast<Double_t>(multMBTracks),
                               static_cast<Double_t>(multRecTemp),
                               static_cast<Double_t>(multMCTrueTracks)};
    fEventHist->Fill(valuesEvent);

    labelsAll = new Int_t[entries];
    labelsAcc = new Int_t[entries];
    labelsRec = new Int_t[entries];

    for (Int_t i = 0; i < entries; ++i) {
        AliESDtrack* track = (AliESDtrack*)allChargedTracks->At(i);
        if (!track)
            continue;
        if (track->Charge() == 0)
            continue;

        labelsAll[multAll] = TMath::Abs(track->GetLabel());
        multAll++;

        // esd track selection
        // check track cuts

        if (fEsdTrackCuts->AcceptTrack(track) && fAccCuts->AcceptTrack(track)) {
            FillHistograms(track, stack, vtxESD->GetZ(),
                           AlidNdPtHelper::kRecTracks, multMBTracks,
                           multRecTemp, multMCTrueTracks);
            labelsRec[multRec] = TMath::Abs(track->GetLabel());
            multRec++;
        }
    }

    // Post output data.
    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskCutTest::Terminate(Option_t*) {
    // Draw result to the screen
    // Called once at the end of the query
    //   Printf("TERMINATE");
}

//_____________________________________________________________________________
void AliAnalysisTaskCutTest::FillHistograms(
    AliESDtrack* const esdTrack, AliStack* const stack, const Double_t zv,
    AlidNdPtHelper::TrackObject trackObj, Int_t multMB, Int_t multRecTrack,
    Int_t multTrueMC) {
    if (trackObj != AlidNdPtHelper::kRecTracks)
        return;
    //
    // Fill ESD track and MC histograms
    //
    if (!esdTrack)
        return;

    Double_t q = esdTrack->Charge();
    if (TMath::Abs(q) < 0.001)
        return;

    Double_t pt = esdTrack->Pt();
    // Float_t qpt = esdTrack->Pt() * q;
    Double_t eta = esdTrack->Eta();
    Double_t phi = esdTrack->Phi();

    Float_t dca[2], bCov[3];
    esdTrack->GetImpactParameters(dca, bCov);

    Int_t nClust = esdTrack->GetTPCclusters(0);
    Double_t chi2PerCluster = 0.;
    if (nClust > 0.)
        chi2PerCluster = esdTrack->GetTPCchi2() / Float_t(nClust);

    Double_t sigmaPt = pt * TMath::Sqrt(esdTrack->GetSigma1Pt2());

    fPtHist->Fill(esdTrack->GetSigned1Pt(),
                  TMath::Sqrt(esdTrack->GetSigma1Pt2()));
    if (TMath::Abs(esdTrack->GetSigned1Pt()) < 1.) {
        fPtHist2->Fill(esdTrack->GetSigned1Pt(),
                       TMath::Sqrt(esdTrack->GetSigma1Pt2()));
    }

    //
    // Fill rec vs MC information
    //
    Double_t deltaPt = 0;

    if (stack) {
        Int_t label = TMath::Abs(esdTrack->GetLabel());
        // if(label == 0) return;

        if (label > stack->GetNtrack())
            return;
        TParticle* particle = stack->Particle(label);
        if (!particle)
            return;

        // Bool_t prim = stack->IsPhysicalPrimary(label);
        // Int_t pid = AlidNdPtHelper::ConvertPdgToPid(particle);

        Int_t motherPdg = -1;
        TParticle* mother = 0;

        // TParticle* prim_mother =
        // AlidNdPtHelper::FindPrimaryMother(stack,label);
        Int_t motherLabel = particle->GetMother(0);
        if (motherLabel > 0)
            mother = stack->Particle(motherLabel);
        if (mother)
            motherPdg = TMath::Abs(
                mother->GetPdgCode()); // take abs for visualisation only
        // Int_t mech = particle->GetUniqueID(); // production mechanism

        if (!particle->GetPDG())
            return;
        Double_t gq = particle->GetPDG()->Charge() / 3.0; // Charge units |e|/3
        if (TMath::Abs(gq) < 0.001)
            return;
        Double_t gpt = particle->Pt();
        Double_t geta = particle->Eta();
        Double_t deltaoversigma = -42;
        // Float_t qgpt = particle->Pt() * gq;
        // Float_t gphi = particle->Phi();

        // printf("pt %f, gpt %f \n",pt,gpt);
        if (gpt) {
            // Double_t dpt = 0;
            // dpt = (pt-gpt)/gpt;
            // changed 24-07-2013 was: deltaPt = gpt*((1./pt)-(1./gpt));
            deltaPt = pt * ((1. / gpt) - (1. / pt));
            deltaoversigma = (TMath::Sqrt(esdTrack->GetSigma1Pt2()) > 0)
                                 ? ((1. / gpt) - (1. / pt)) /
                                       TMath::Sqrt(esdTrack->GetSigma1Pt2())
                                 : -42;
        }

        //   Double_t deta = (eta-geta);
        //   fMCTrackHist->Fill(valuesMC);
        Double_t valuesMC2[3] = {pt, gpt, deltaoversigma};
        fMCTrackHist2->Fill(valuesMC2);

    } // end if(stack)

    Double_t values[8] = {
        zv, pt, eta, phi, q, sigmaPt, static_cast<Double_t>(multMB), deltaPt};
    Double_t values2[5] = {zv, pt, eta, phi, q};
    fTrackHist2->Fill(values2);
    if (TMath::Abs(eta) > 0.8)
        return;
    if (TMath::Abs(zv) > 10.0)
        return;
    fTrackHist->Fill(values);
}
