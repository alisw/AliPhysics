#include "AliAnalysisTaskCmeEse.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TRandom.h>
#include <TSpline.h>
#include <TGrid.h>
#include <THnSparse.h>

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "AliAODVZERO.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"

// STL includes
#include <iostream>
#include <ctime>
#include <sys/time.h>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCmeEse)
    //_____________________________________________________________________________
    AliAnalysisTaskCmeEse::AliAnalysisTaskCmeEse():
        AliAnalysisTaskSE(),
        fAOD(0),
        fRun(-1),
        fMultV0(0),
        fQx2mV0A(0),
        fQy2mV0A(0),
        fQx2sV0A(0),
        fQy2sV0A(0),
        fQx2mV0C(0),
        fQy2mV0C(0),
        fQx2sV0C(0),
        fQy2sV0C(0),
        fFitRecLow(0),
        fFitRecHigh(0),
        fRecEff(0),
        fVtxCut(10.0),  
        fFilterbit(128),
        fEtaCut(0.8),
        fNoClus(70),
        fMinPt(0.2),
        fMaxPt(20.0),
        fLHC10h(kTRUE),
        fPileUp(kTRUE),
        fCent(0),
        fQAV0(kFALSE),
        fTrkQA(kFALSE),
        fListOfObjects(0),
        fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0), fMultCorBeforeCuts(0), fMultCorAfterCuts(0), fPercqc2(0), fAllQA(0),
        fQxavsV0Bef(0), fQyavsV0Bef(0), fQxcvsV0Bef(0), fQycvsV0Bef(0),
        fQxavsVtxZBef(0), fQyavsVtxZBef(0), fQxcvsVtxZBef(0), fQycvsVtxZBef(0),
        fQxavsV0Aft(0), fQyavsV0Aft(0), fQxcvsV0Aft(0), fQycvsV0Aft(0),
        fQxavsVtxZAft(0), fQyavsVtxZAft(0), fQxcvsVtxZAft(0), fQycvsVtxZAft(0),
        fV2A(0), fV0AV0Cv2(0), fV0ATPCv2(0), fV0CTPCv2(0),
        fCMESameQPos(0), fCMESameQNeg(0), fCMEOppQ(0), fCMESameQPosCos(0), fCMESameQPosSin(0), fCMESameQPosv1(0), fCMESameQNegCos(0), fCMESameQNegSin(0), fCMESameQNegv1(0), fCMEOppQCos(0), fCMEOppQSin(0), fCMEOppQv1(0)
{

    for (Int_t i = 0; i < 90; i++){

        fSplQ2c[i] = 0;

        if (i < 10){

            fV2Aqc2[i] = 0;
            fV0AV0Cv2qc2[i] = 0;
            fV0ATPCv2qc2[i] = 0;
            fV0CTPCv2qc2[i] = 0;
            fCMESameQPosqc2[i] = 0;
            fCMESameQNegqc2[i] = 0;
            fCMEOppQqc2[i] = 0;
            fCMESameQPosCosqc2[i] = 0;
            fCMESameQPosSinqc2[i] = 0;
            fCMESameQPosv1qc2[i] = 0;
            fCMESameQNegCosqc2[i] = 0;
            fCMESameQNegSinqc2[i] = 0;
            fCMESameQNegv1qc2[i] = 0;
            fCMEOppQCosqc2[i] = 0;
            fCMEOppQSinqc2[i] = 0;
            fCMEOppQv1qc2[i] = 0;
        }

        if (i < 9){

            fV2Pt[i] = 0;

            fCMESameQPosEta[i] = 0;
            fCMESameQPosPtDif[i] = 0;
            fCMESameQPosPtSum[i] = 0;

            fCMESameQNegEta[i] = 0;
            fCMESameQNegPtDif[i] = 0;
            fCMESameQNegPtSum[i] = 0;

            fCMEOppQEta[i] = 0;
            fCMEOppQPtDif[i] = 0;
            fCMEOppQPtSum[i] = 0;

            fPsiA[i] = 0;
            fPsiC[i] = 0;
            fPsiAvsPsiC[i] = 0;


            for (Int_t j = 0; j < 10; j++){

                fV2Ptqc2[i][j] = 0;

                fCMESameQPosEtaqc2[i][j] = 0;
                fCMESameQPosPtDifqc2[i][j] = 0;
                fCMESameQPosPtSumqc2[i][j] = 0;

                fCMESameQNegEtaqc2[i][j] = 0;
                fCMESameQNegPtDifqc2[i][j] = 0;
                fCMESameQNegPtSumqc2[i][j] = 0;

                fCMEOppQEtaqc2[i][j] = 0;
                fCMEOppQPtDifqc2[i][j] = 0;
                fCMEOppQPtSumqc2[i][j] = 0;
            }

        }

    }

}

//______________________________________________________________________________
AliAnalysisTaskCmeEse::AliAnalysisTaskCmeEse(const char *name):
    AliAnalysisTaskSE(name),
    fAOD(0),
    fRun(-1),
    fMultV0(0),
    fQx2mV0A(0),
    fQy2mV0A(0),
    fQx2sV0A(0),
    fQy2sV0A(0),
    fQx2mV0C(0),
    fQy2mV0C(0),
    fQx2sV0C(0),
    fQy2sV0C(0),
    fFitRecLow(0),
    fFitRecHigh(0),
    fRecEff(0),
    fVtxCut(10.0),  
    fFilterbit(128),
    fEtaCut(0.8),
    fNoClus(70),
    fMinPt(0.2),
    fMaxPt(20.0),
    fLHC10h(kTRUE),
    fPileUp(kTRUE),
    fCent(0),
    fQAV0(kFALSE),
    fTrkQA(kFALSE),
    fListOfObjects(0),
    fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0), fMultCorBeforeCuts(0), fMultCorAfterCuts(0), fPercqc2(0), fAllQA(0),
    fQxavsV0Bef(0), fQyavsV0Bef(0), fQxcvsV0Bef(0), fQycvsV0Bef(0),
    fQxavsVtxZBef(0), fQyavsVtxZBef(0), fQxcvsVtxZBef(0), fQycvsVtxZBef(0),
    fQxavsV0Aft(0), fQyavsV0Aft(0), fQxcvsV0Aft(0), fQycvsV0Aft(0),
    fQxavsVtxZAft(0), fQyavsVtxZAft(0), fQxcvsVtxZAft(0), fQycvsVtxZAft(0),
    fV2A(0), fV0AV0Cv2(0), fV0ATPCv2(0), fV0CTPCv2(0),
    fCMESameQPos(0), fCMESameQNeg(0), fCMEOppQ(0), fCMESameQPosCos(0), fCMESameQPosSin(0), fCMESameQPosv1(0), fCMESameQNegCos(0), fCMESameQNegSin(0), fCMESameQNegv1(0), fCMEOppQCos(0), fCMEOppQSin(0), fCMEOppQv1(0)
{

    for (Int_t i = 0; i < 90; i++){

        fSplQ2c[i] = 0;

        if (i < 10){

            fV2Aqc2[i] = 0;
            fV0AV0Cv2qc2[i] = 0;
            fV0ATPCv2qc2[i] = 0;
            fV0CTPCv2qc2[i] = 0;
            fCMESameQPosqc2[i] = 0;
            fCMESameQNegqc2[i] = 0;
            fCMEOppQqc2[i] = 0;
            fCMESameQPosCosqc2[i] = 0;
            fCMESameQPosSinqc2[i] = 0;
            fCMESameQPosv1qc2[i] = 0;
            fCMESameQNegCosqc2[i] = 0;
            fCMESameQNegSinqc2[i] = 0;
            fCMESameQNegv1qc2[i] = 0;
            fCMEOppQCosqc2[i] = 0;
            fCMEOppQSinqc2[i] = 0;
            fCMEOppQv1qc2[i] = 0;
        }

        if (i < 9){

            fV2Pt[i] = 0;

            fCMESameQPosEta[i] = 0;
            fCMESameQPosPtDif[i] = 0;
            fCMESameQPosPtSum[i] = 0;

            fCMESameQNegEta[i] = 0;
            fCMESameQNegPtDif[i] = 0;
            fCMESameQNegPtSum[i] = 0;

            fCMEOppQEta[i] = 0;
            fCMEOppQPtDif[i] = 0;
            fCMEOppQPtSum[i] = 0;

            fPsiA[i] = 0;
            fPsiC[i] = 0;
            fPsiAvsPsiC[i] = 0;


            for (Int_t j = 0; j < 10; j++){

                fV2Ptqc2[i][j] = 0;

                fCMESameQPosEtaqc2[i][j] = 0;
                fCMESameQPosPtDifqc2[i][j] = 0;
                fCMESameQPosPtSumqc2[i][j] = 0;

                fCMESameQNegEtaqc2[i][j] = 0;
                fCMESameQNegPtDifqc2[i][j] = 0;
                fCMESameQNegPtSumqc2[i][j] = 0;

                fCMEOppQEtaqc2[i][j] = 0;
                fCMEOppQPtDifqc2[i][j] = 0;
                fCMEOppQPtSumqc2[i][j] = 0;
            }

        }

    }

    // Output slot #1 writes into a TTree
    DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliAnalysisTaskCmeEse::~AliAnalysisTaskCmeEse()
{
    // Destructor
    if (fListOfObjects) 
        delete fListOfObjects;

}

//______________________________________________________________________________
void AliAnalysisTaskCmeEse::UserCreateOutputObjects()
{ 

    if (!gGrid) {
        TGrid::Connect("alien://");
    }


    TFile* fRecEf = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/hist_recEff_pbpb_run1.root");
    if (!fRecEf){
        printf("Rec eff file cannot be opened \n");
        return;
    }
    fFitRecLow = (TF1*)fRecEf->Get("fitEff1_0");
    fFitRecHigh = (TF1*)fRecEf->Get("fitEff2_0");



    TFile* fSpl = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibSpV0CRun1.root");
    if (!fSpl){
        printf("Spline file cannot be opened \n");
        return;
    }
    for (Int_t isp = 0; isp < 90; isp++)
        fSplQ2c[isp] = (TSpline3*)fSpl->Get(Form("hqc2Int_%d", isp));



    timeval a;
    gettimeofday(&a, 0);
    Int_t randomSeed = a.tv_usec;
    gRandom->SetSeed(randomSeed);



    OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner();


    // Histograms
    fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
    fListOfObjects->Add(fVtx);

    fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
    fListOfObjects->Add(fVtxBeforeCuts);

    fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
    fListOfObjects->Add(fVtxAfterCuts);

    fMultCorBeforeCuts = new TH2F("fMultCorBeforeCuts", "TPC vs Global multiplicity (Before cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fListOfObjects->Add(fMultCorBeforeCuts);

    fMultCorAfterCuts = new TH2F("fMultCorAfterCuts", "TPC vs Global multiplicity (After cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fListOfObjects->Add(fMultCorAfterCuts);


    fPercqc2 = new TH2D("fPercqc2", "; centrality percentile; q_{2} percentile", 100, 0, 100, 100, 0, 100);
    fListOfObjects->Add(fPercqc2);



    const Int_t nCenB = 9;
    Float_t cenBins[nCenB+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80.};

    const Int_t nPtB = 23;
    Double_t ptBins[nPtB+1] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0};


    if (fTrkQA){

        Int_t    binsAll[5] = {     nCenB, nPtB,    3,   18,             72};
        Double_t xminAll[5] = {      -0.5,    0, -1.5, -0.9,              0};
        Double_t xmaxAll[5] = { nCenB-0.5,   5.,  1.5,  0.9, 2.*TMath::Pi()};

        fAllQA = new THnSparseF("fAllQA", "cent:pT:q:eta:phi", 5, binsAll, xminAll, xmaxAll);
        fAllQA->SetBinEdges(1, ptBins);
        fAllQA->GetAxis(0)->SetTitle("centrality");
        fAllQA->GetAxis(1)->SetTitle("p_{T} (Gev/c");
        fAllQA->GetAxis(2)->SetTitle("q");
        fAllQA->GetAxis(3)->SetTitle("#eta");
        fAllQA->GetAxis(4)->SetTitle("#varphi");
        fListOfObjects->Add(fAllQA);

    }



    if (fQAV0){

        fQxavsV0Bef = new TH2D("fQxavsV0Bef", " ; centrality (V0); Q_{x} (V0A)", 100, 0, 100, 1400, -1400., 1400.);
        fListOfObjects->Add(fQxavsV0Bef);

        fQyavsV0Bef = new TH2D("fQyavsV0Bef", " ; centrality (V0); Q_{y} (V0A)", 100, 0, 100, 1400, -1400., 1400.);
        fListOfObjects->Add(fQyavsV0Bef);

        fQxcvsV0Bef = new TH2D("fQxcvsV0Bef", " ; centrality (V0); Q_{x} (V0C)", 100, 0, 100, 1400, -1400., 1400.);
        fListOfObjects->Add(fQxcvsV0Bef);

        fQycvsV0Bef = new TH2D("fQycvsV0Bef", " ; centrality (V0); Q_{y} (V0C)", 100, 0, 100, 1400, -1400., 1400.);
        fListOfObjects->Add(fQycvsV0Bef);


        fQxavsVtxZBef = new TH2D("fQxavsVtxZBef", " ; vertexZ; Q_{x} (V0A)", 20, -10., 10., 1400, -1400., 1400.);
        fListOfObjects->Add(fQxavsVtxZBef);

        fQyavsVtxZBef = new TH2D("fQyavsVtxZBef", " ; vertexZ; Q_{y} (V0A)", 20, -10., 10., 1400, -1400., 1400.);
        fListOfObjects->Add(fQyavsVtxZBef);

        fQxcvsVtxZBef = new TH2D("fQxcvsVtxZBef", " ; vertexZ; Q_{x} (V0C)", 20, -10., 10., 1400, -1400., 1400.);
        fListOfObjects->Add(fQxcvsVtxZBef);

        fQycvsVtxZBef = new TH2D("fQycvsVtxZBef", " ; vertexZ; Q_{y} (V0C)", 20, -10., 10., 1400, -1400., 1400.);
        fListOfObjects->Add(fQycvsVtxZBef);


        fQxavsV0Aft = new TH2D("fQxavsV0Aft", " ; centrality (V0); Q_{x} (V0A)", 100, 0, 100, 200, -10., 10.);
        fListOfObjects->Add(fQxavsV0Aft);

        fQyavsV0Aft = new TH2D("fQyavsV0Aft", " ; centrality (V0); Q_{y} (V0A)", 100, 0, 100, 200, -10., 10.);
        fListOfObjects->Add(fQyavsV0Aft);

        fQxcvsV0Aft = new TH2D("fQxcvsV0Aft", " ; centrality (V0); Q_{x} (V0C)", 100, 0, 100, 200, -10., 10.);
        fListOfObjects->Add(fQxcvsV0Aft);

        fQycvsV0Aft = new TH2D("fQycvsV0Aft", " ; centrality (V0); Q_{y} (V0C)", 100, 0, 100, 200, -10., 10.);
        fListOfObjects->Add(fQycvsV0Aft);


        fQxavsVtxZAft = new TH2D("fQxavsVtxZAft", " ; vertexZ; Q_{x} (V0A)", 20, -10., 10., 200, -10., 10);
        fListOfObjects->Add(fQxavsVtxZAft);

        fQyavsVtxZAft = new TH2D("fQyavsVtxZAft", " ; vertexZ; Q_{y} (V0A)", 20, -10., 10., 200, -10., 10);
        fListOfObjects->Add(fQyavsVtxZAft);

        fQxcvsVtxZAft = new TH2D("fQxcvsVtxZAft", " ; vertexZ; Q_{x} (V0C)", 20, -10., 10., 200, -10., 10);
        fListOfObjects->Add(fQxcvsVtxZAft);

        fQycvsVtxZAft = new TH2D("fQycvsVtxZAft", " ; vertexZ; Q_{y} (V0C)", 20, -10., 10., 200, -10., 10);
        fListOfObjects->Add(fQycvsVtxZAft);

    }



    fV2A = new TProfile("fV2A", "; centrality percentile; v_{2}{EP V0A}", nCenB, cenBins);
    fListOfObjects->Add(fV2A);

    fV0AV0Cv2 = new TProfile("fV0AV0Cv2", "; centrality percentile; V0A-V0C correlations", nCenB, cenBins);
    fListOfObjects->Add(fV0AV0Cv2);

    fV0ATPCv2 = new TProfile("fV0ATPCv2", "; centrality percentile; V0A-TPC correlations", nCenB, cenBins);
    fListOfObjects->Add(fV0ATPCv2);

    fV0CTPCv2 = new TProfile("fV0CTPCv2", "; centrality percentile; V0C-TPC correlations", nCenB, cenBins);
    fListOfObjects->Add(fV0CTPCv2);


    fCMESameQPos = new TProfile("fCMESameQPos", "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQPos);

    fCMESameQNeg = new TProfile("fCMESameQNeg", "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQNeg);

    fCMEOppQ = new TProfile("fCMEOppQ", "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMEOppQ);


    fCMESameQPosCos = new TProfile("fCMESameQPosCos", "; centrality percentile; #LT cos(#varphi_1)cos(#varphi_2) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQPosCos);

    fCMESameQPosSin = new TProfile("fCMESameQPosSin", "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQPosSin);

    fCMESameQPosv1 = new TProfile("fCMESameQPosv1", "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQPosv1);


    fCMESameQNegCos = new TProfile("fCMESameQNegCos", "; centrality percentile; #LT cos(#varphi_{1})cos(#varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQNegCos);

    fCMESameQNegSin = new TProfile("fCMESameQNegSin", "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQNegSin);

    fCMESameQNegv1 = new TProfile("fCMESameQNegv1", "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMESameQNegv1);


    fCMEOppQCos = new TProfile("fCMEOppQCos", "; centrality percentile; #LT cos(#varphi_{1})cos(#varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMEOppQCos);

    fCMEOppQSin = new TProfile("fCMEOppQSin", "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMEOppQSin);

    fCMEOppQv1 = new TProfile("fCMEOppQv1", "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
    fListOfObjects->Add(fCMEOppQv1);


    for (Int_t i = 0 ; i < 10; i++){

        fV2Aqc2[i] = new TProfile(Form("fV2Aqc2_%d", i), "; centrality percentile; v_{2}{EP V0A}", nCenB, cenBins);
        fListOfObjects->Add(fV2Aqc2[i]);

        fV0AV0Cv2qc2[i] = new TProfile(Form("fV0AV0Cv2qc2_%d", i), "; centrality percentile; V0A-V0C correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0AV0Cv2qc2[i]);

        fV0ATPCv2qc2[i] = new TProfile(Form("fV0ATPCv2qc2_%d", i), "; centrality percentile; V0A-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0ATPCv2qc2[i]);

        fV0CTPCv2qc2[i] = new TProfile(Form("fV0CTPCv2qc2_%d", i), "; centrality percentile; V0C-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0CTPCv2qc2[i]);



        fCMESameQPosqc2[i] = new TProfile(Form("fCMESameQPosqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQPosqc2[i]);

        fCMESameQNegqc2[i] = new TProfile(Form("fCMESameQNegqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQNegqc2[i]);

        fCMEOppQqc2[i] = new TProfile(Form("fCMEOppQqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMEOppQqc2[i]);


        fCMESameQPosCosqc2[i] = new TProfile(Form("fCMESameQPosCosqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1})cos(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQPosCosqc2[i]);

        fCMESameQPosSinqc2[i] = new TProfile(Form("fCMESameQPosSinqc2_%d", i), "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQPosSinqc2[i]);

        fCMESameQPosv1qc2[i] = new TProfile(Form("fCMESameQPosv1qc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQPosv1qc2[i]);


        fCMESameQNegCosqc2[i] = new TProfile(Form("fCMESameQNegCosqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1})cos(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQNegCosqc2[i]);

        fCMESameQNegSinqc2[i] = new TProfile(Form("fCMESameQNegSinqc2_%d", i), "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQNegSinqc2[i]);

        fCMESameQNegv1qc2[i] = new TProfile(Form("fCMESameQNegv1qc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMESameQNegv1qc2[i]);


        fCMEOppQCosqc2[i] = new TProfile(Form("fCMEOppQCosqc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1})cos(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMEOppQCosqc2[i]);

        fCMEOppQSinqc2[i] = new TProfile(Form("fCMEOppQSinqc2_%d", i), "; centrality percentile; #LT sin(#varphi_{1})sin(#varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMEOppQSinqc2[i]);

        fCMEOppQv1qc2[i] = new TProfile(Form("fCMEOppQv1qc2_%d", i), "; centrality percentile; #LT cos(#varphi_{1} - #varphi_{2}) #GT", nCenB, cenBins);
        fListOfObjects->Add(fCMEOppQv1qc2[i]);



        if (i < 9){

            fV2Pt[i] = new TProfile(Form("fV2Pt_%d", i), "; p_{T} (GeV/c); v_{2}", nPtB, ptBins);
            fListOfObjects->Add(fV2Pt[i]);

            fCMESameQPosEta[i] = new TProfile(Form("fCMESameQPosEta_%d", i), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
            fListOfObjects->Add(fCMESameQPosEta[i]);

            fCMESameQPosPtDif[i] = new TProfile(Form("fCMESameQPosPtDif_%d", i), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMESameQPosPtDif[i]);

            fCMESameQPosPtSum[i] = new TProfile(Form("fCMESameQPosPtSum_%d", i), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMESameQPosPtSum[i]);


            fCMESameQNegEta[i] = new TProfile(Form("fCMESameQNegEta_%d", i), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
            fListOfObjects->Add(fCMESameQNegEta[i]);

            fCMESameQNegPtDif[i] = new TProfile(Form("fCMESameQNegPtDif_%d", i), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMESameQNegPtDif[i]);

            fCMESameQNegPtSum[i] = new TProfile(Form("fCMESameQNegPtSum_%d", i), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMESameQNegPtSum[i]);


            fCMEOppQEta[i] = new TProfile(Form("fCMEOppQEta_%d", i), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
            fListOfObjects->Add(fCMEOppQEta[i]);

            fCMEOppQPtDif[i] = new TProfile(Form("fCMEOppQPtDif_%d", i), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMEOppQPtDif[i]);

            fCMEOppQPtSum[i] = new TProfile(Form("fCMEOppQPtSum_%d", i), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
            fListOfObjects->Add(fCMEOppQPtSum[i]);


            if (fQAV0){

                fPsiA[i] = new TH1D(Form("fPsiA_%d", i), "; #Psi_{n} (V0A); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiA[i]);

                fPsiC[i] = new TH1D(Form("fPsiC_%d", i), "; #Psi_{n} (V0C); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiC[i]);

                fPsiAvsPsiC[i] = new TH2D(Form("fPsiAvsPsiC_%d", i), "; #Psi_{n} (V0A); #Psi_{n} (V0C)", 36, -1.8, 1.8, 36, -1.8, 1.8);
                fListOfObjects->Add(fPsiAvsPsiC[i]);

            }



            for (Int_t j = 0; j < 10; j++){

                fV2Ptqc2[i][j] = new TProfile(Form("fV2Ptqc2_%d_%d", i, j), "; p_{T} (GeV/c); v_{2}", nPtB, ptBins);
                fListOfObjects->Add(fV2Ptqc2[i][j]);

                fCMESameQPosEtaqc2[i][j] = new TProfile(Form("fCMESameQPosEtaqc2_%d_%d", i, j), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
                fListOfObjects->Add(fCMESameQPosEtaqc2[i][j]);

                fCMESameQPosPtDifqc2[i][j] = new TProfile(Form("fCMESameQPosPtDifqc2_%d_%d", i, j), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMESameQPosPtDifqc2[i][j]);

                fCMESameQPosPtSumqc2[i][j] = new TProfile(Form("fCMESameQPosPtSumqc2_%d_%d", i, j), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMESameQPosPtSumqc2[i][j]);


                fCMESameQNegEtaqc2[i][j] = new TProfile(Form("fCMESameQNegEtaqc2_%d_%d", i, j), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
                fListOfObjects->Add(fCMESameQNegEtaqc2[i][j]);

                fCMESameQNegPtDifqc2[i][j] = new TProfile(Form("fCMESameQNegPtDifqc2_%d_%d", i, j), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMESameQNegPtDifqc2[i][j]);

                fCMESameQNegPtSumqc2[i][j] = new TProfile(Form("fCMESameQNegPtSumqc2_%d_%d", i, j), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMESameQNegPtSumqc2[i][j]);


                fCMEOppQEtaqc2[i][j] = new TProfile(Form("fCMEOppQEtaqc2_%d_%d", i, j), "; #Delta#eta = |#eta_{1} - #eta_{2}|; #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 8, 0, 1.6);
                fListOfObjects->Add(fCMEOppQEtaqc2[i][j]);

                fCMEOppQPtDifqc2[i][j] = new TProfile(Form("fCMEOppQPtDifqc2_%d_%d", i, j), "; |p_{T}^{1} - p_{T}^{2}| (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMEOppQPtDifqc2[i][j]);

                fCMEOppQPtSumqc2[i][j] = new TProfile(Form("fCMEOppQPtSumqc2_%d_%d", i, j), "; (p_{T}^{1} + p_{T}^{2})/2 (GeV/c); #LT cos(#varphi_{1} + #varphi_{2} - 2*#Psi_{2}) #GT", 7, 0, 2.8);
                fListOfObjects->Add(fCMEOppQPtSumqc2[i][j]);

            }

        }

    }


    // Post output data.
    PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskCmeEse::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD || !fAOD->GetHeader()){
        Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }

    Int_t run = fAOD->GetRunNumber();
    if(run != fRun){
        // Load the calibrations run dependent
        OpenInfoCalbration(run);
        fRun = run;
    } 

    Float_t zvtx = GetVertex(fAOD);

    if(zvtx< -990){
        fVtx->Fill(0);
    } else {
        fVtx->Fill(1);
        fVtxBeforeCuts->Fill(zvtx);
        if (TMath::Abs(zvtx) < fVtxCut) {
            fMultCorBeforeCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));

            if (fPileUp){

                if (plpMV(fAOD))
                    return;

                if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0)
                    return;

                if ((fLHC10h) && ((fFilterbit == 128) || (fFilterbit == 1)) && ((Float_t(GetTPCMult(fAOD)) < (-40.3+1.22*GetGlobalMult(fAOD))) || (Float_t(GetTPCMult(fAOD)) > (32.1+1.59*GetGlobalMult(fAOD))) ) )
                    return;

                if ((!fLHC10h) && ((fFilterbit == 128) || (fFilterbit == 1)) && ((Float_t(GetTPCMult(fAOD)) < (-36.73+1.48*GetGlobalMult(fAOD))) || (Float_t(GetTPCMult(fAOD)) > (62.87+1.78*GetGlobalMult(fAOD))) ) )
                    return;

            }

            fMultCorAfterCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));
            Analyze(fAOD, zvtx);
            fVtxAfterCuts->Fill(zvtx);

        }

    }

    // Post output data.
    PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskCmeEse::Analyze(AliAODEvent* aod, Float_t vtxZ)
{  


    AliCentrality* centrality = ((AliAODHeader*)aod->GetHeader())->GetCentralityP();
    if (!centrality){
        //Printf("No centrality - abort \n");
        return;
    }
    Double_t centTrk = centrality->GetCentralityPercentile("TRK");
    Double_t centV0  = centrality->GetCentralityPercentile("V0M");
    Double_t centSPD = centrality->GetCentralityPercentile("CL1");

    if (TMath::Abs(centV0 - centTrk) > 7.5 || centV0 >= 80 || centV0 <= 0)
        return;


    Int_t iCentV0 = (Int_t)centV0;

    if (fCent == 1)
        centV0 = centTrk;
    else if (fCent == 2)
        centV0 = centSPD;

    Int_t iCentSPD = (Int_t)centSPD;
    if (iCentSPD >= 90)
        return;

    Short_t centrCode = -10;    
    if ((centV0 > 0) &&  (centV0 < 5.))
        centrCode = 0;
    else if ((centV0 >= 5.) && (centV0 < 10.))
        centrCode = 1;
    else if ((centV0 >= 10.) && (centV0 < 20.))
        centrCode = 2;
    else if ((centV0 >= 20.) && (centV0 < 30.))
        centrCode = 3;
    else if ((centV0 >= 30.) && (centV0 < 40.))
        centrCode = 4;
    else if ((centV0 >= 40.) && (centV0 < 50.))
        centrCode = 5;
    else if ((centV0 >= 50.) && (centV0 < 60.))
        centrCode = 6;
    else if ((centV0 >= 60.) && (centV0 < 70.))
        centrCode = 7;
    else if ((centV0 >= 70.) && (centV0 < 80.))
        centrCode = 8;

    if (centrCode < 0)
        return;


    //V0 info
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;    
    Double_t sumMa = 0, sumMc = 0;

    AliAODVZERO* aodV0 = aod->GetVZEROData();

    for (Int_t iV0 = 0; iV0 < 64; iV0++) {

        Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);

        Float_t multv0 = aodV0->GetMultiplicity(iV0);

        if (iV0 < 32){

            Double_t multCorC = -10;

            if (iV0 < 8)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
            else if (iV0 >= 8 && iV0 < 16)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
            else if (iV0 >= 16 && iV0 < 24)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
            else if (iV0 >= 24 && iV0 < 32)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);

            if (multCorC < 0){
                //cout<<"Problem with multiplicity in V0C"<<endl;
                continue;
            }

            Qxc2 += TMath::Cos(2.*phiV0) * multCorC;
            Qyc2 += TMath::Sin(2.*phiV0) * multCorC;

            sumMc = sumMc + multCorC;

        } else {

            Double_t multCorA = -10;

            if (iV0 >= 32 && iV0 < 40)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
            else if (iV0 >= 40 && iV0 < 48)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
            else if (iV0 >= 48 && iV0 < 56)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
            else if (iV0 >= 56 && iV0 < 64)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);

            if (multCorA < 0){
                //cout<<"Problem with multiplicity in V0A"<<endl;
                continue;
            }

            Qxa2 += TMath::Cos(2.*phiV0) * multCorA;
            Qya2 += TMath::Sin(2.*phiV0) * multCorA;

            sumMa = sumMa + multCorA;

        }
    }


    if (sumMa <= 0 || sumMc <= 0)
        return;


    Double_t Qxa2Cor = (Qxa2 - fQx2mV0A->GetBinContent(iCentSPD+1))/fQx2sV0A->GetBinContent(iCentSPD+1);
    Double_t Qya2Cor = (Qya2 - fQy2mV0A->GetBinContent(iCentSPD+1))/fQy2sV0A->GetBinContent(iCentSPD+1);
    Double_t Qxc2Cor = (Qxc2 - fQx2mV0C->GetBinContent(iCentSPD+1))/fQx2sV0C->GetBinContent(iCentSPD+1);
    Double_t Qyc2Cor = (Qyc2 - fQy2mV0C->GetBinContent(iCentSPD+1))/fQy2sV0C->GetBinContent(iCentSPD+1);

    Double_t evPlAngV0A = TMath::ATan2(Qya2Cor, Qxa2Cor)/2.;
    Double_t evPlAngV0C = TMath::ATan2(Qyc2Cor, Qxc2Cor)/2.;


    if (fQAV0){

        fPsiA[centrCode]->Fill(evPlAngV0A);
        fPsiC[centrCode]->Fill(evPlAngV0C);
        fPsiAvsPsiC[centrCode]->Fill(evPlAngV0A, evPlAngV0C);

        fQxavsV0Bef->Fill(centV0, Qxa2);
        fQyavsV0Bef->Fill(centV0, Qya2);
        fQxcvsV0Bef->Fill(centV0, Qxc2);
        fQycvsV0Bef->Fill(centV0, Qyc2);

        fQxavsVtxZBef->Fill(vtxZ, Qxa2);
        fQyavsVtxZBef->Fill(vtxZ, Qya2);
        fQxcvsVtxZBef->Fill(vtxZ, Qxc2);
        fQycvsVtxZBef->Fill(vtxZ, Qyc2);

        fQxavsV0Aft->Fill(centV0, Qxa2Cor);
        fQyavsV0Aft->Fill(centV0, Qya2Cor);
        fQxcvsV0Aft->Fill(centV0, Qxc2Cor);
        fQycvsV0Aft->Fill(centV0, Qyc2Cor);

        fQxavsVtxZAft->Fill(vtxZ, Qxa2Cor);
        fQyavsVtxZAft->Fill(vtxZ, Qya2Cor);
        fQxcvsVtxZAft->Fill(vtxZ, Qxc2Cor);
        fQycvsVtxZAft->Fill(vtxZ, Qyc2Cor);

    }


    Double_t Qxc2Corese = Qxc2 - fQx2mV0C->GetBinContent(iCentSPD+1);
    Double_t Qyc2Corese = Qyc2 - fQy2mV0C->GetBinContent(iCentSPD+1);

    Double_t qc2 = TMath::Sqrt((Qxc2Corese*Qxc2Corese + Qyc2Corese*Qyc2Corese)/sumMc);

    Double_t percqc2 = 100.*fSplQ2c[iCentV0]->Eval(qc2);
    Short_t percCqc2 = GetPercCode(percqc2);
    if (percCqc2 < 0){
        //Printf("Problem with percentile: negative code! \n");
        return;
    }

    fPercqc2->Fill(centV0, percqc2);


    const Int_t nAODTracks = aod->GetNumberOfTracks();

    AliAODTrack* trkV[nAODTracks];
    Double_t Qxt2 = 0, Qyt2 = 0;
    Int_t sumMt = 0;


    for (Int_t it1 = 0; it1 < nAODTracks; it1++) {

        AliAODTrack* aodTrk1 = (AliAODTrack*)aod->GetTrack(it1);

        if (!aodTrk1){
            delete aodTrk1;
            continue;
        }

        if (!(aodTrk1->TestFilterBit(fFilterbit)))
            continue;

        if ((TMath::Abs(aodTrk1->Eta()) > fEtaCut) || (aodTrk1->GetTPCNcls() < fNoClus) || (aodTrk1->Pt() < fMinPt) || (aodTrk1->Pt() > fMaxPt))
            continue;


        if (fRecEff){

            Double_t eff = -1.;
            if (aodTrk1->Pt() < 3.5)
                eff = fFitRecLow->Eval(aodTrk1->Pt());
            else
                eff = fFitRecHigh->Eval(aodTrk1->Pt());

            if (eff < 0){
                cout<<"Problem with reconstruction efficiency"<<endl;
                continue;
            }

            Double_t rnd = gRandom->Rndm();
            Double_t minRec = fFitRecLow->Eval(0.2);

            if (rnd > minRec/eff){
                //cout<<"Track rejected: "<<it1<<"  from "<<nAODTracks<<endl;
                continue;
            }

        }


        trkV[sumMt] = aodTrk1;

        Qxt2 += TMath::Cos(2.*aodTrk1->Phi());
        Qyt2 += TMath::Sin(2.*aodTrk1->Phi());

        sumMt++;

        Double_t v2a = TMath::Cos(2.*(aodTrk1->Phi() - evPlAngV0A));
        fV2A->Fill(centV0, v2a);
        fV2Pt[centrCode]->Fill(aodTrk1->Pt(), v2a);
        fV2Aqc2[percCqc2]->Fill(centV0, v2a);
        fV2Ptqc2[centrCode][percCqc2]->Fill(aodTrk1->Pt(), v2a);

        if (fTrkQA){

            Double_t allQA[5] = {Double_t(centrCode),
                aodTrk1->Pt(),
                Double_t(aodTrk1->Charge()),
                aodTrk1->Eta(),
                aodTrk1->Phi()};
            fAllQA->Fill(allQA);
        }

    }


    Double_t evPlAngTPC = TMath::ATan2(Qyt2, Qxt2)/2.;


    fV0AV0Cv2->Fill(centV0, TMath::Cos(2.*(evPlAngV0A - evPlAngV0C)));
    fV0ATPCv2->Fill(centV0, TMath::Cos(2.*(evPlAngV0A - evPlAngTPC)));
    fV0CTPCv2->Fill(centV0, TMath::Cos(2.*(evPlAngV0C - evPlAngTPC)));

    fV0AV0Cv2qc2[percCqc2]->Fill(centV0, TMath::Cos(2.*(evPlAngV0A - evPlAngV0C)));
    fV0ATPCv2qc2[percCqc2]->Fill(centV0, TMath::Cos(2.*(evPlAngV0A - evPlAngTPC)));
    fV0CTPCv2qc2[percCqc2]->Fill(centV0, TMath::Cos(2.*(evPlAngV0C - evPlAngTPC)));



    for (Int_t it = 0; it < sumMt-1; it++){

        for (Int_t jt = it+1; jt < sumMt; jt++){

            if (trkV[it]->Charge() > 0 && trkV[jt]->Charge() > 0){

                fCMESameQPos->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQPosqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

                fCMESameQPosCos->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));
                fCMESameQPosCosqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));

                fCMESameQPosSin->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));
                fCMESameQPosSinqc2[percCqc2]->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));

                fCMESameQPosv1->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));
                fCMESameQPosv1qc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));


                fCMESameQPosEta[centrCode]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQPosPtDif[centrCode]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQPosPtSum[centrCode]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));


                fCMESameQPosEtaqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQPosPtDifqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQPosPtSumqc2[centrCode][percCqc2]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

            }


            if (trkV[it]->Charge() < 0 && trkV[jt]->Charge() < 0){

                fCMESameQNeg->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQNegqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

                fCMESameQNegCos->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));
                fCMESameQNegCosqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));

                fCMESameQNegSin->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));
                fCMESameQNegSinqc2[percCqc2]->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));

                fCMESameQNegv1->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));
                fCMESameQNegv1qc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));


                fCMESameQNegEta[centrCode]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQNegPtDif[centrCode]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQNegPtSum[centrCode]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));


                fCMESameQNegEtaqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQNegPtDifqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMESameQNegPtSumqc2[centrCode][percCqc2]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

            }



            if ((trkV[it]->Charge() < 0 && trkV[jt]->Charge() > 0) || (trkV[it]->Charge() > 0 && trkV[jt]->Charge() < 0)){

                fCMEOppQ->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMEOppQqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

                fCMEOppQCos->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));
                fCMEOppQCosqc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - evPlAngV0A)*TMath::Cos(trkV[jt]->Phi() - evPlAngV0A));

                fCMEOppQSin->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));
                fCMEOppQSinqc2[percCqc2]->Fill(centV0, TMath::Sin(trkV[it]->Phi() - evPlAngV0A)*TMath::Sin(trkV[jt]->Phi() - evPlAngV0A));

                fCMEOppQv1->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));
                fCMEOppQv1qc2[percCqc2]->Fill(centV0, TMath::Cos(trkV[it]->Phi() - trkV[jt]->Phi()));


                fCMEOppQEta[centrCode]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMEOppQPtDif[centrCode]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMEOppQPtSum[centrCode]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));


                fCMEOppQEtaqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Eta() - trkV[jt]->Eta()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMEOppQPtDifqc2[centrCode][percCqc2]->Fill(TMath::Abs(trkV[it]->Pt() - trkV[jt]->Pt()), TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));
                fCMEOppQPtSumqc2[centrCode][percCqc2]->Fill((trkV[it]->Pt() + trkV[jt]->Pt())/2., TMath::Cos(trkV[it]->Phi() + trkV[jt]->Phi() - 2.*evPlAngV0A));

            }

        }
    }

}

//_____________________________________________________________________________
void AliAnalysisTaskCmeEse::OpenInfoCalbration(Int_t run)
{

    if (!gGrid) {
        TGrid::Connect("alien://");
    }

    TFile* foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0Run1.root");

    if(!foadb){
        printf("OADB V0 calibration file cannot be opened\n");
        return;
    }

    AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
    if(!cont){
        printf("OADB object hMultV0BefCorr is not available in the file\n");
        return;
    }
    if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
    }
    fMultV0 = ((TH1D*) cont->GetObject(run));



    AliOADBContainer* contQx2am = (AliOADBContainer*) foadb->Get("fqxa2m");
    if(!contQx2am){
        printf("OADB object fqxa2m is not available in the file\n");
        return;
    }
    if(!(contQx2am->GetObject(run))){
        printf("OADB object fqxa2m is not available for run %i\n", run);
        return;
    }
    fQx2mV0A = ((TH1D*) contQx2am->GetObject(run));


    AliOADBContainer* contQy2am = (AliOADBContainer*) foadb->Get("fqya2m");
    if(!contQy2am){
        printf("OADB object fqya2m is not available in the file\n");
        return;
    }
    if(!(contQy2am->GetObject(run))){
        printf("OADB object fqya2m is not available for run %i\n", run);
        return;
    }
    fQy2mV0A = ((TH1D*) contQy2am->GetObject(run));


    AliOADBContainer* contQx2as = (AliOADBContainer*) foadb->Get("fqxa2s");
    if(!contQx2as){
        printf("OADB object fqxa2s is not available in the file\n");
        return;
    }
    if(!(contQx2as->GetObject(run))){
        printf("OADB object fqxa2s is not available for run %i\n", run);
        return;
    }
    fQx2sV0A = ((TH1D*) contQx2as->GetObject(run));


    AliOADBContainer* contQy2as = (AliOADBContainer*) foadb->Get("fqya2s");
    if(!contQy2as){
        printf("OADB object fqya2s is not available in the file\n");
        return;
    }
    if(!(contQy2as->GetObject(run))){
        printf("OADB object fqya2s is not available for run %i\n", run);
        return;
    }
    fQy2sV0A = ((TH1D*) contQy2as->GetObject(run));



    AliOADBContainer* contQx2cm = (AliOADBContainer*) foadb->Get("fqxc2m");
    if(!contQx2cm){
        printf("OADB object fqxc2m is not available in the file\n");
        return;
    }
    if(!(contQx2cm->GetObject(run))){
        printf("OADB object fqxc2m is not available for run %i\n", run);
        return;
    }
    fQx2mV0C = ((TH1D*) contQx2cm->GetObject(run));


    AliOADBContainer* contQy2cm = (AliOADBContainer*) foadb->Get("fqyc2m");
    if(!contQy2cm){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
    if(!(contQy2cm->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
    }
    fQy2mV0C = ((TH1D*) contQy2cm->GetObject(run));



    AliOADBContainer* contQx2cs = (AliOADBContainer*) foadb->Get("fqxc2s");
    if(!contQx2cs){
        printf("OADB object fqxc2s is not available in the file\n");
        return;
    }
    if(!(contQx2cs->GetObject(run))){
        printf("OADB object fqxc2s is not available for run %i\n", run);
        return;
    }
    fQx2sV0C = ((TH1D*) contQx2cs->GetObject(run));


    AliOADBContainer* contQy2cs = (AliOADBContainer*) foadb->Get("fqyc2s");
    if(!contQy2cs){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
    if(!(contQy2cs->GetObject(run))){
        printf("OADB object fqyc2s is not available for run %i\n", run);
        return;
    }
    fQy2sV0C = ((TH1D*) contQy2cs->GetObject(run));

}


//____________________________________________________________________
Int_t AliAnalysisTaskCmeEse::GetTPCMult(AliVEvent* ev) const
{

    Int_t multTPC = 0;

    AliAODEvent* aod = (AliAODEvent*)ev;

    const Int_t nGoodTracks = aod->GetNumberOfTracks();

    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {

        AliAODTrack* trackAOD = (AliAODTrack*)aod->GetTrack(iTracks);

        if (!trackAOD){
            delete trackAOD;
            continue;
        }

        if (!(trackAOD->TestFilterBit(1)))
            continue;

        if ((trackAOD->Pt() < fMinPt) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > fEtaCut) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0))
            continue;

        multTPC++;

    }

    return multTPC;

}

//____________________________________________________________________
Int_t AliAnalysisTaskCmeEse::GetGlobalMult(AliVEvent* ev) const
{

    Int_t multGlobal = 0;

    AliAODEvent* aod = (AliAODEvent*)ev;

    const Int_t nGoodTracks = aod->GetNumberOfTracks();

    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {

        AliAODTrack* trackAOD = (AliAODTrack*)aod->GetTrack(iTracks);

        if (!trackAOD){
            delete trackAOD;
            continue;
        }

        if (!(trackAOD->TestFilterBit(16)))
            continue;

        if ((trackAOD->Pt() < fMinPt) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > fEtaCut) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) )
            continue;


        // clone for constraining
        Double_t b[2] = { -99., -99.};
        Double_t bCov[3] = { -99., -99., -99.};

        AliAODTrack* trackAODC = new AliAODTrack(*trackAOD);
        if (!trackAODC) {
            AliWarning("Clone of AOD track failed.");
            delete trackAODC;
            continue;
        }

        if (!trackAODC->PropagateToDCA(aod->GetPrimaryVertex(), aod->GetMagneticField(), 100., b, bCov)){
            delete trackAODC;
            continue;
        } else {
            delete trackAODC;
        }


        if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3))
            continue;

        multGlobal++;

    }//track loop

    return multGlobal;

}


//_____________________________________________________________________________
Float_t AliAnalysisTaskCmeEse::GetVertex(AliAODEvent* aod) const
{

    Float_t vtxz = -999;

    const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
    if (!trkVtx || trkVtx->GetNContributors()<=0) 
        return vtxz;
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks")) 
        return vtxz;
    const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
    if (!spdVtx || spdVtx->GetNContributors()<=0) 
        return vtxz;
    TString vtxTyp = spdVtx->GetTitle();
    Double_t cov[6]={0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) 
        return vtxz;
    if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) 
        return vtxz;

    vtxz = trkVtx->GetZ();

    return vtxz;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskCmeEse::GetPercCode(Double_t perc) const
{

    Short_t percCode = -1;

    if ((perc >= 0) && (perc <= 10.0))
        percCode = 0;
    else if ((perc > 10.0) && (perc <= 20.0))
        percCode = 1;
    else if ((perc > 20.0) && (perc <= 30.0))
        percCode = 2;
    else if ((perc > 30.0) && (perc <= 40.0))
        percCode = 3;
    else if ((perc > 40.0) && (perc <= 50.0))
        percCode = 4;
    else if ((perc > 50.0) && (perc <= 60.0))
        percCode = 5;
    else if ((perc > 60.0) && (perc <= 70.0))
        percCode = 6;
    else if ((perc > 70.0) && (perc <= 80.0))
        percCode = 7;
    else if ((perc > 80.0) && (perc <= 90.0))
        percCode = 8;
    else if (perc > 90.0)
        percCode = 9;  

    return percCode;

}

//_____________________________________________________________________________
Double_t AliAnalysisTaskCmeEse::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{

    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
        return 0;
    }
    static TMatrixDSym vVb(3);
    double dist = -1;
    double dx = v0->GetX()-v1->GetX();
    double dy = v0->GetY()-v1->GetY();
    double dz = v0->GetZ()-v1->GetZ();
    double cov0[6],cov1[6];
    v0->GetCovarianceMatrix(cov0);
    v1->GetCovarianceMatrix(cov1);
    vVb(0,0) = cov0[0]+cov1[0];
    vVb(1,1) = cov0[2]+cov1[2];
    vVb(2,2) = cov0[5]+cov1[5];
    vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
    vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
    vVb.InvertFast();
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
        +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;

}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskCmeEse::plpMV(const AliVEvent *event)
{
    // check for multi-vertexer pile-up
    const AliAODEvent *aod = (const AliAODEvent*)event;
    const AliESDEvent *esd = (const AliESDEvent*)event;
    //
    const int    kMinPlpContrib = 5;
    const double kMaxPlpChi2 = 5.0;
    const double kMinWDist = 15;
    //
    if (!aod && !esd) {
        printf("Event is neither of AOD nor ESD\n");
        exit(1);
    }
    //
    const AliVVertex* vtPrm = 0;
    const AliVVertex* vtPlp = 0;
    int nPlp = 0;
    //
    if (aod) {
        if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
        vtPrm = aod->GetPrimaryVertex();
        if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
    }
    else {
        if ( !(nPlp=esd->GetNumberOfPileupVerticesTracks())) return kFALSE;
        vtPrm = esd->GetPrimaryVertexTracks();
        if (((AliESDVertex*)vtPrm)->GetStatus()!=1) return kTRUE; // there are pile-up vertices but no primary
    }

    //int bcPrim = vtPrm->GetBC();
    //
    for (int ipl=0;ipl<nPlp;ipl++) {
        vtPlp = aod ? (const AliVVertex*)aod->GetPileupVertexTracks(ipl) : (const AliVVertex*)esd->GetPileupVertexTracks(ipl);
        //
        if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
        if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
        //  int bcPlp = vtPlp->GetBC();
        //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
        //
        double wDst = GetWDist(vtPrm,vtPlp);
        if (wDst<kMinWDist) continue;
        //
        return kTRUE; // pile-up: well separated vertices
    }
    //
    return kFALSE;
    //
}

//_____________________________________________________________________________
void AliAnalysisTaskCmeEse::Terminate(Option_t *)
{ 
    // Terminate loop
    Printf("Terminate()");
}
