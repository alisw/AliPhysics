#include "AliAnalysisTaskPiKpK0Lamba.h"

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
#include <TVector3.h>
#include <THnSparse.h>
#include <TProfile2D.h>
#include <TGrid.h>


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
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliAODv0.h"
#include "AliMultSelection.h"
#include "AliAODcascade.h"
#include "AliAODTracklets.h"

// STL includes
#include <iostream>
#include <ctime>
#include <sys/time.h>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPiKpK0Lamba)
    //_____________________________________________________________________________
    AliAnalysisTaskPiKpK0Lamba::AliAnalysisTaskPiKpK0Lamba():
        AliAnalysisTaskSE(),
        fAOD(0),
        fPIDResponse(0),
        fRun(-1),
        fMultV0(0),
        fQxnmV0A(0),
        fQynmV0A(0),
        fQxnsV0A(0),
        fQynsV0A(0),
        fQxnmV0C(0),
        fQynmV0C(0),
        fQxnsV0C(0),
        fQynsV0C(0),
        fLowCut(0),
        fHighCut(0),
        fVtxCut(10.0),  
        fFilterbit(128),
        fEtaCut(0.8),
        fNoClus(70),
        fMinPt(0.2),
        fMaxPt(20.0),
        fPileUp(kTRUE),
        fPileUpTOF(kFALSE),
        fCent(0),
        fNHarm(2.),
        fNsigCut(3.),
        fQA(kFALSE),
        fPhiCut(kFALSE),
        fPhiCutLow(0),
        fPhiCutHigh(0),
        fMultTOFLowCut(0),
        fMultTOFHighCut(0),
        fMinPiCut(1.),
        fMaxPiCut(7.),
        fMinPCut(-25.),
        fMaxPCut(-14.),
        fQPos(kFALSE),
        fQNeg(kFALSE),
        fExclPID(kFALSE),
        fNoClusPid(70),
        fQAV0(kFALSE),
        fMultCentLowCut(0),
        fQAOutl(0),
        fBin1Cent(0),
        fRemChV0A(0),
        fEtaRange(2),
        fRemPhiReg(kFALSE),
        fCutMultESDdif(15000.),
        fCrsRowsFrcShCls(kFALSE),
        flagPsi42A(kTRUE),
        fNPtBins(17),
        fNcrFind(0.8),
        fDCADghtPV(0.1),
        fMaxDCADght(0.5),
        fCosPA(0.998),
        fMinRad(5.),
        fMaxRad(100.),
        fArmPodCut(kFALSE),
        fMinPtDght(kFALSE),
        fListOfObjects(0),
        fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0),
        fMultvsCentBef(0), fMultvsCentAft(0),
        fCentrBef(0), fCentrAft(0),
        fCenCL0vsV0MBef(0), fCenCL1vsV0MBef(0), fCenCL0vsCL1Bef(0),
        fCenCL0vsV0MAft(0), fCenCL1vsV0MAft(0), fCenCL0vsCL1Aft(0),
        fMultFBvsMultFBTOFBef(0), fMultFBvsMultFBTOFAft(0),
        fMultESDDifvsMultTPCBef(0), fMultESDDifvsMultTPCAft(0),
        fSPclsvsSPDtrksBef(0), fSPclsvsSPDtrksAft(0),
        fMultV0vsMultTPCoutBef(0), fMultV0vsMultTPCoutAft(0),
        fPidQA(0), fAllQA(0), fV0QA(0), fPidQAB1C(0), fAllQAB1C(0),
        fQxavsV0Bef(0), fQyavsV0Bef(0), fQxcvsV0Bef(0), fQycvsV0Bef(0),
        fQxavsVtxZBef(0), fQyavsVtxZBef(0), fQxcvsVtxZBef(0), fQycvsVtxZBef(0),
        fQxavsV0Aft(0), fQyavsV0Aft(0), fQxcvsV0Aft(0), fQycvsV0Aft(0),
        fQxavsVtxZAft(0), fQyavsVtxZAft(0), fQxcvsVtxZAft(0), fQycvsVtxZAft(0),
        fV0AV0Cvn(0), fV0ATPCvn(0), fV0CTPCvn(0),
        fV0AV0CvnB1C(0), fV0ATPCvnB1C(0), fV0CTPCvnB1C(0),
        fV0AV0Cvnsq(0), fV0ATPCvnsq(0), fV0CTPCvnsq(0),
        fV0AV0CvnB1Csq(0), fV0ATPCvnB1Csq(0), fV0CTPCvnB1Csq(0)
{

    for (Int_t j = 0; j < fNPtBins+1; j++)
        fPtBins[j] = 0;

    for (Int_t i = 0; i < 90; i++){

        if (i < 10){

            fVnAllA[i] = 0;
            fVnPiA[i] = 0;
            fVnKA[i] = 0;
            fVnAntiPA[i] = 0;
            fVnPihighPtA[i] = 0;
            fVnPhighPtA[i] = 0;

            fVnAllC[i] = 0;
            fVnPiC[i] = 0;
            fVnKC[i] = 0;
            fVnAntiPC[i] = 0;
            fVnPihighPtC[i] = 0;
            fVnPhighPtC[i] = 0;


            fPsiA[i] = 0;
            fPsiC[i] = 0;
            fPsiAvsPsiC[i] = 0;

            fSinTrkCosV0A[i] = 0;
            fCosTrkSinV0A[i] = 0;

            fSinTrkCosV0C[i] = 0;
            fCosTrkSinV0C[i] = 0;

            fSinTrkSinV0A[i] = 0;
            fCosTrkCosV0A[i] = 0;

            fSinTrkSinV0C[i] = 0;
            fCosTrkCosV0C[i] = 0;

            for (Int_t ii = 0; ii < fNPtBins; ii++){

                fInvMassK0[ii][i] = 0;
                fVnK0A[ii][i] = 0;
                fVnK0C[ii][i] = 0;
                fInvMassL[ii][i] = 0;
                fVnLA[ii][i] = 0;
                fVnLC[ii][i] = 0;

                if (i == 0){
                    fInvMassK0B1C[ii][i] = 0;
                    fVnK0AB1C[ii][i] = 0;
                    fVnK0CB1C[ii][i] = 0;
                    fInvMassLB1C[ii][i] = 0;
                    fVnLAB1C[ii][i] = 0;
                    fVnLCB1C[ii][i] = 0;
                }

            }

        }


        fVnAllAB1C[i] = 0;
        fVnPiAB1C[i] = 0;
        fVnKAB1C[i] = 0;
        fVnAntiPAB1C[i] = 0;
        fVnPihighPtAB1C[i] = 0;
        fVnPhighPtAB1C[i] = 0;

        fVnAllCB1C[i] = 0;
        fVnPiCB1C[i] = 0;
        fVnKCB1C[i] = 0;
        fVnAntiPCB1C[i] = 0;
        fVnPihighPtCB1C[i] = 0;
        fVnPhighPtCB1C[i] = 0;


        fPsiAB1C[i] = 0;
        fPsiCB1C[i] = 0;
        fPsiAvsPsiCB1C[i] = 0;

        fSinTrkCosV0AB1C[i] = 0;
        fCosTrkSinV0AB1C[i] = 0;

        fSinTrkCosV0CB1C[i] = 0;
        fCosTrkSinV0CB1C[i] = 0;

        fSinTrkSinV0AB1C[i] = 0;
        fCosTrkCosV0AB1C[i] = 0;

        fSinTrkSinV0CB1C[i] = 0;
        fCosTrkCosV0CB1C[i] = 0;

    }

}

//______________________________________________________________________________
AliAnalysisTaskPiKpK0Lamba::AliAnalysisTaskPiKpK0Lamba(const char *name):
    AliAnalysisTaskSE(name),
    fAOD(0),
    fPIDResponse(0),
    fRun(-1),
    fMultV0(0),
    fQxnmV0A(0),
    fQynmV0A(0),
    fQxnsV0A(0),
    fQynsV0A(0),
    fQxnmV0C(0),
    fQynmV0C(0),
    fQxnsV0C(0),
    fQynsV0C(0),
    fLowCut(0),
    fHighCut(0),
    fVtxCut(10.0),
    fFilterbit(128),
    fEtaCut(0.8),
    fNoClus(70),
    fMinPt(0.2),
    fMaxPt(20.0),
    fPileUp(kTRUE),
    fPileUpTOF(kFALSE),
    fCent(0),
    fNHarm(2.),
    fNsigCut(3.),
    fQA(kFALSE),
    fPhiCut(kFALSE),
    fPhiCutLow(0),
    fPhiCutHigh(0),
    fMultTOFLowCut(0),
    fMultTOFHighCut(0),
    fMinPiCut(1.),
    fMaxPiCut(7.),
    fMinPCut(-25.),
    fMaxPCut(-14.),
    fQPos(kFALSE),
    fQNeg(kFALSE),
    fExclPID(kFALSE),
    fNoClusPid(70),
    fQAV0(kFALSE),
    fMultCentLowCut(0),
    fQAOutl(0),
    fBin1Cent(0),
    fRemChV0A(0),
    fEtaRange(2),
    fRemPhiReg(kFALSE),
    fCutMultESDdif(15000.),
    fCrsRowsFrcShCls(kFALSE),
    flagPsi42A(kTRUE),
    fNPtBins(17),
    fNcrFind(0.8),
    fDCADghtPV(0.1),
    fMaxDCADght(0.5),
    fCosPA(0.998),
    fMinRad(5.),
    fMaxRad(100.),
    fArmPodCut(kFALSE),
    fMinPtDght(kFALSE),
    fListOfObjects(0),
    fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0),
    fMultvsCentBef(0), fMultvsCentAft(0),
    fCentrBef(0), fCentrAft(0),
    fCenCL0vsV0MBef(0), fCenCL1vsV0MBef(0), fCenCL0vsCL1Bef(0),
    fCenCL0vsV0MAft(0), fCenCL1vsV0MAft(0), fCenCL0vsCL1Aft(0),
    fMultFBvsMultFBTOFBef(0), fMultFBvsMultFBTOFAft(0),
    fMultESDDifvsMultTPCBef(0), fMultESDDifvsMultTPCAft(0),
    fSPclsvsSPDtrksBef(0), fSPclsvsSPDtrksAft(0),
    fMultV0vsMultTPCoutBef(0), fMultV0vsMultTPCoutAft(0),
    fPidQA(0), fAllQA(0), fV0QA(0), fPidQAB1C(0), fAllQAB1C(0),
    fQxavsV0Bef(0), fQyavsV0Bef(0), fQxcvsV0Bef(0), fQycvsV0Bef(0),
    fQxavsVtxZBef(0), fQyavsVtxZBef(0), fQxcvsVtxZBef(0), fQycvsVtxZBef(0),
    fQxavsV0Aft(0), fQyavsV0Aft(0), fQxcvsV0Aft(0), fQycvsV0Aft(0),
    fQxavsVtxZAft(0), fQyavsVtxZAft(0), fQxcvsVtxZAft(0), fQycvsVtxZAft(0),
    fV0AV0Cvn(0), fV0ATPCvn(0), fV0CTPCvn(0),
    fV0AV0CvnB1C(0), fV0ATPCvnB1C(0), fV0CTPCvnB1C(0),
    fV0AV0Cvnsq(0), fV0ATPCvnsq(0), fV0CTPCvnsq(0),
    fV0AV0CvnB1Csq(0), fV0ATPCvnB1Csq(0), fV0CTPCvnB1Csq(0)
{

    for (Int_t j = 0; j < fNPtBins+1; j++)
        fPtBins[j] = 0;

    for (Int_t i = 0; i < 90; i++){

        if (i < 10){

            fVnAllA[i] = 0;
            fVnPiA[i] = 0;
            fVnKA[i] = 0;
            fVnAntiPA[i] = 0;
            fVnPihighPtA[i] = 0;
            fVnPhighPtA[i] = 0;

            fVnAllC[i] = 0;
            fVnPiC[i] = 0;
            fVnKC[i] = 0;
            fVnAntiPC[i] = 0;
            fVnPihighPtC[i] = 0;
            fVnPhighPtC[i] = 0;


            fPsiA[i] = 0;
            fPsiC[i] = 0;
            fPsiAvsPsiC[i] = 0;

            fSinTrkCosV0A[i] = 0;
            fCosTrkSinV0A[i] = 0;

            fSinTrkCosV0C[i] = 0;
            fCosTrkSinV0C[i] = 0;

            fSinTrkSinV0A[i] = 0;
            fCosTrkCosV0A[i] = 0;

            fSinTrkSinV0C[i] = 0;
            fCosTrkCosV0C[i] = 0;

            for (Int_t ii = 0; ii < fNPtBins; ii++){

                fInvMassK0[ii][i] = 0;
                fVnK0A[ii][i] = 0;
                fVnK0C[ii][i] = 0;
                fInvMassL[ii][i] = 0;
                fVnLA[ii][i] = 0;
                fVnLC[ii][i] = 0;

                if (i == 0){
                    fInvMassK0B1C[ii][i] = 0;
                    fVnK0AB1C[ii][i] = 0;
                    fVnK0CB1C[ii][i] = 0;
                    fInvMassLB1C[ii][i] = 0;
                    fVnLAB1C[ii][i] = 0;
                    fVnLCB1C[ii][i] = 0;
                }

            }

        }


        fVnAllAB1C[i] = 0;
        fVnPiAB1C[i] = 0;
        fVnKAB1C[i] = 0;
        fVnAntiPAB1C[i] = 0;
        fVnPihighPtAB1C[i] = 0;
        fVnPhighPtAB1C[i] = 0;

        fVnAllCB1C[i] = 0;
        fVnPiCB1C[i] = 0;
        fVnKCB1C[i] = 0;
        fVnAntiPCB1C[i] = 0;
        fVnPihighPtCB1C[i] = 0;
        fVnPhighPtCB1C[i] = 0;


        fPsiAB1C[i] = 0;
        fPsiCB1C[i] = 0;
        fPsiAvsPsiCB1C[i] = 0;

        fSinTrkCosV0AB1C[i] = 0;
        fCosTrkSinV0AB1C[i] = 0;

        fSinTrkCosV0CB1C[i] = 0;
        fCosTrkSinV0CB1C[i] = 0;

        fSinTrkSinV0AB1C[i] = 0;
        fCosTrkCosV0AB1C[i] = 0;

        fSinTrkSinV0CB1C[i] = 0;
        fCosTrkCosV0CB1C[i] = 0;

    }

    // Output slot #1 writes into a TTree
    DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliAnalysisTaskPiKpK0Lamba::~AliAnalysisTaskPiKpK0Lamba()
{
    // Destructor
    if (fListOfObjects) 
        delete fListOfObjects;

}

//______________________________________________________________________________
void AliAnalysisTaskPiKpK0Lamba::UserCreateOutputObjects()
{ 


    fLowCut = new TF1("fLowCut", "[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fHighCut = new TF1("fHighCut", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);


    //fPhiCutLow = new TF1("fPhiCutLow",  "-0.01/x+pi/18.0-0.015", 0, 100);
    fPhiCutLow = new TF1("fPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);

    //fPhiCutHigh = new TF1("fPhiCutHigh", "0.55/x/x+pi/18.0+0.03", 0, 100);
    fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);


    fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);

    fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);


    fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);


    fLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
    fHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);

    fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);

    fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);




    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
        AliFatal("Input handler needed");

    //pid response object
    fPIDResponse=inputHandler->GetPIDResponse();
    if (!fPIDResponse)
        AliError("PIDResponse object was not created");


    OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner();


    const Int_t nCenB = 10;
    Float_t cenBins[nCenB+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};

    const Int_t nPtB = 36;
    Double_t ptBins[nPtB+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0, 30.0, 40.0, 50.0};

    const Int_t nPtBL = 20;
    Double_t ptBinsL[nPtBL+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};

    const Int_t nPtBH = 15;
    Double_t ptBinsH[nPtBH+1] = {3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 12.0, 16.0, 20.0, 30.0, 40.0, 50.0};


    if (fQAOutl){

        fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
        fListOfObjects->Add(fVtx);

        fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
        fListOfObjects->Add(fVtxBeforeCuts);

        fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
        fListOfObjects->Add(fVtxAfterCuts);


        fMultvsCentBef = new TH2F("fMultvsCentBef", "; centrality; multiplicity", 100, 0, 100, 1000, -0.5, 3999.5);
        fListOfObjects->Add(fMultvsCentBef);

        fMultvsCentAft = new TH2F("fMultvsCentAft", "; centrality; multiplicity", 100, 0, 100, 1000, -0.5, 3999.5);
        fListOfObjects->Add(fMultvsCentAft);


        fCentrBef = new TH1F("fCentrBef", "; centrality; Counts", 100, 0, 100);
        fListOfObjects->Add(fCentrBef);

        fCentrAft = new TH1F("fCentrAft", "; centrality; Counts", 100, 0, 100);
        fListOfObjects->Add(fCentrAft);


        fCenCL0vsV0MBef = new TH2F("fCenCL0vsV0MBef", "; centrality V0M; centrality CL0", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL0vsV0MBef);

        fCenCL1vsV0MBef = new TH2F("fCenCL1vsV0MBef", "; centrality V0M; centrality CL1", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL1vsV0MBef);

        fCenCL0vsCL1Bef = new TH2F("fCenCL0vsCL1Bef", "; centrality CL1; centrality CL0", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL0vsCL1Bef);


        fCenCL0vsV0MAft = new TH2F("fCenCL0vsV0MAft", "; centrality V0M; centrality CL0", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL0vsV0MAft);

        fCenCL1vsV0MAft = new TH2F("fCenCL1vsV0MAft", "; centrality V0M; centrality CL1", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL1vsV0MAft);

        fCenCL0vsCL1Aft = new TH2F("fCenCL0vsCL1Aft", "; centrality CL1; centrality CL0", 100, 0, 100, 100, 0, 100);
        fListOfObjects->Add(fCenCL0vsCL1Aft);



        fMultFBvsMultFBTOFBef = new TH2I("fMultFBvsMultFBTOFBef", "; multiplicity (FB32); multiplicity (FB32 + TOF)", 1000, -0.5, 3999.5, 1000, -0.5, 1999.5);
        fListOfObjects->Add(fMultFBvsMultFBTOFBef);

        fMultFBvsMultFBTOFAft = new TH2I("fMultFBvsMultFBTOFAft", "; multiplicity (FB32); multiplicity (FB32 + TOF)", 1000, -0.5, 3999.5, 1000, -0.5, 1999.5);
        fListOfObjects->Add(fMultFBvsMultFBTOFAft);


        fMultESDDifvsMultTPCBef = new TH2D("fMultESDDifvsMultTPCBef", "; multiplicity (TPC); multiplicity (ESD - FB128*3.38)", 1000, 0, 7000, 1000, -1000, 29000);
        fListOfObjects->Add(fMultESDDifvsMultTPCBef);

        fMultESDDifvsMultTPCAft = new TH2D("fMultESDDifvsMultTPCAft", "; multiplicity (TPC); multiplicity (ESD - FB128*3.38)", 1000, 0, 7000, 1000, -1000, 29000);
        fListOfObjects->Add(fMultESDDifvsMultTPCAft);


        fSPclsvsSPDtrksBef = new TH2I("fSPclsvsSPDtrksBef", "; SPD N_{tracklets}; SPD N_{clusters}", 1000, -0.5, 6999.5, 1000, -0.5, 24999.5);
        fListOfObjects->Add(fSPclsvsSPDtrksBef);

        fSPclsvsSPDtrksAft = new TH2I("fSPclsvsSPDtrksAft", "; SPD N_{tracklets}; SPD N_{clusters}", 1000, -0.5, 6999.5, 1000, -0.5, 24999.5);
        fListOfObjects->Add(fSPclsvsSPDtrksAft);



        fMultV0vsMultTPCoutBef = new TH2F("fMultV0vsMultTPCoutBef", "; TPC out; V0 multiplicity", 1000, 0, 30000, 1000, 0, 40000);
        fListOfObjects->Add(fMultV0vsMultTPCoutBef);

        fMultV0vsMultTPCoutAft = new TH2F("fMultV0vsMultTPCoutAft", "; TPC out; V0 multiplicity", 1000, 0, 30000, 1000, 0, 40000);
        fListOfObjects->Add(fMultV0vsMultTPCoutAft);

    }


    if (fQA){

        if (!fBin1Cent){

            Int_t    binsPid[6] = {   10, 125,   65,    3,    3, 50};
            Double_t xminPid[6] = { -0.5,   0, -0.1, -1.5, -0.5,  0};
            Double_t xmaxPid[6] = {  9.5, 250,  1.2,  1.5,  2.5, 5.};


            Int_t    binsAll[5] = {   10,   nPtB,    3,   18,             72};
            Double_t xminAll[5] = { -0.5, fMinPt, -1.5, -0.9,              0};
            Double_t xmaxAll[5] = {  9.5, fMaxPt,  1.5,  0.9, 2.*TMath::Pi()};


            fPidQA = new THnSparseF("fPidQA", "cent:dEdx:beta:q:pidFl:p", 6, binsPid, xminPid, xmaxPid);
            //fPidQA->Sumw2();
            fPidQA->GetAxis(0)->SetTitle("centrality");
            fPidQA->GetAxis(1)->SetTitle("dE/dx");
            fPidQA->GetAxis(2)->SetTitle("#beta (TOF)");
            fPidQA->GetAxis(3)->SetTitle("q");
            fPidQA->GetAxis(4)->SetTitle("pid flag");
            fPidQA->GetAxis(5)->SetTitle("p (Gev/c");
            fListOfObjects->Add(fPidQA);


            fAllQA = new THnSparseF("fAllQA", "cent:pT:q:eta:phi", 5, binsAll, xminAll, xmaxAll);
            //fAllQA->Sumw2();
            fAllQA->SetBinEdges(1, ptBins);
            fAllQA->GetAxis(0)->SetTitle("centrality");
            fAllQA->GetAxis(1)->SetTitle("p_{T} (Gev/c");
            fAllQA->GetAxis(2)->SetTitle("q");
            fAllQA->GetAxis(3)->SetTitle("#eta");
            fAllQA->GetAxis(4)->SetTitle("#varphi");
            fListOfObjects->Add(fAllQA);



            Int_t    bins[8] = {   10, fNPtBins, 30, 30,   30,  30,   30, 2};
            Double_t xmin[8] = { -0.5,   fMinPt,  0,  0, 0.98,   0,    0, 0};
            Double_t xmax[8] = {  9.5,   fMaxPt, 3., 3.,   1., 0.6, 100., 2};

            fV0QA = new THnSparseF("fV0QA", "cent:pT:DCATrkPos:DCATrkNeg:cosPA:DCAdght:Rad:V0fl", 8, bins, xmin, xmax);
            fV0QA->SetBinEdges(1, fPtBins);
            fV0QA->GetAxis(0)->SetTitle("centrality");
            fV0QA->GetAxis(1)->SetTitle("p_{T} (Gev/c");
            fV0QA->GetAxis(2)->SetTitle("DCATrkPosPV");
            fV0QA->GetAxis(3)->SetTitle("DCATrkNegPV");
            fV0QA->GetAxis(4)->SetTitle("cosPA");
            fV0QA->GetAxis(5)->SetTitle("DCAdght");
            fV0QA->GetAxis(6)->SetTitle("Radius");
            fV0QA->GetAxis(7)->SetTitle("V0 flag");
            fListOfObjects->Add(fV0QA);

        } else {

            Int_t    binsPidB1C[6] = { 90, 125,   65,    3,    3, 50};
            Double_t xminPidB1C[6] = {  0,   0, -0.1, -1.5, -0.5,  0};
            Double_t xmaxPidB1C[6] = { 90, 250,  1.2,  1.5,  2.5, 5.};


            Int_t    binsAllB1C[5] = { 90,   nPtB,    3,   18,             72};
            Double_t xminAllB1C[5] = {  0, fMinPt, -1.5, -0.9,              0};
            Double_t xmaxAllB1C[5] = { 90, fMaxPt,  1.5,  0.9, 2.*TMath::Pi()};


            fPidQAB1C = new THnSparseF("fPidQAB1C", "cent:dEdx:beta:q:pidFl:p", 6, binsPidB1C, xminPidB1C, xmaxPidB1C);
            //fPidQA->Sumw2();
            fPidQAB1C->GetAxis(0)->SetTitle("centrality");
            fPidQAB1C->GetAxis(1)->SetTitle("dE/dx");
            fPidQAB1C->GetAxis(2)->SetTitle("#beta (TOF)");
            fPidQAB1C->GetAxis(3)->SetTitle("q");
            fPidQAB1C->GetAxis(4)->SetTitle("pid flag");
            fPidQAB1C->GetAxis(5)->SetTitle("p (Gev/c");
            fListOfObjects->Add(fPidQAB1C);


            fAllQAB1C = new THnSparseF("fAllQAB1C", "cent:pT:q:eta:phi", 5, binsAllB1C, xminAllB1C, xmaxAllB1C);
            //fAllQA->Sumw2();
            fAllQAB1C->SetBinEdges(1, ptBins);
            fAllQAB1C->GetAxis(0)->SetTitle("centrality");
            fAllQAB1C->GetAxis(1)->SetTitle("p_{T} (Gev/c");
            fAllQAB1C->GetAxis(2)->SetTitle("q");
            fAllQAB1C->GetAxis(3)->SetTitle("#eta");
            fAllQAB1C->GetAxis(4)->SetTitle("#varphi");
            fListOfObjects->Add(fAllQAB1C);

        }

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


    if (!fBin1Cent){

        fV0AV0Cvn = new TProfile("fV0AV0Cvn", "; centrality percentile; V0A-V0C correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0AV0Cvn);

        fV0ATPCvn = new TProfile("fV0ATPCvn", "; centrality percentile; V0A-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0ATPCvn);

        fV0CTPCvn = new TProfile("fV0CTPCvn", "; centrality percentile; V0C-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0CTPCvn);


        fV0AV0Cvnsq = new TProfile("fV0AV0Cvnsq", "; centrality percentile; V0A-V0C correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0AV0Cvnsq);

        fV0ATPCvnsq = new TProfile("fV0ATPCvnsq", "; centrality percentile; V0A-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0ATPCvnsq);

        fV0CTPCvnsq = new TProfile("fV0CTPCvnsq", "; centrality percentile; V0C-TPC correlations", nCenB, cenBins);
        fListOfObjects->Add(fV0CTPCvnsq);



        for (Int_t ic = 0; ic < nCenB; ic++){

            fVnAllA[ic] = new TProfile(Form("fVnAllA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
            fListOfObjects->Add(fVnAllA[ic]);

            fVnPiA[ic] = new TProfile(Form("fVnPiA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnPiA[ic]);

            fVnKA[ic] = new TProfile(Form("fVnKA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnKA[ic]);

            fVnAntiPA[ic] = new TProfile(Form("fVnAntiPA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnAntiPA[ic]);

            fVnPihighPtA[ic] = new TProfile(Form("fVnPihighPtA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPihighPtA[ic]);

            fVnPhighPtA[ic] = new TProfile(Form("fVnPhighPtA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPhighPtA[ic]);



            fVnAllC[ic] = new TProfile(Form("fVnAllC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
            fListOfObjects->Add(fVnAllC[ic]);

            fVnPiC[ic] = new TProfile(Form("fVnPiC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnPiC[ic]);

            fVnKC[ic] = new TProfile(Form("fVnKC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnKC[ic]);

            fVnAntiPC[ic] = new TProfile(Form("fVnAntiPC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnAntiPC[ic]);

            fVnPihighPtC[ic] = new TProfile(Form("fVnPihighPtC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPihighPtC[ic]);

            fVnPhighPtC[ic] = new TProfile(Form("fVnPhighPtC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPhighPtC[ic]);


            if (fQAV0){

                fPsiA[ic] = new TH1D(Form("fPsiA_%d", ic), "; #Psi_{n} (V0A); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiA[ic]);

                fPsiC[ic] = new TH1D(Form("fPsiC_%d", ic), "; #Psi_{n} (V0C); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiC[ic]);

                fPsiAvsPsiC[ic] = new TH2D(Form("fPsiAvsPsiC_%d", ic), "; #Psi_{n} (V0A); #Psi_{n} (V0C)", 36, -1.8, 1.8, 36, -1.8, 1.8);
                fListOfObjects->Add(fPsiAvsPsiC[ic]);


                fSinTrkCosV0A[ic] = new TProfile(Form("fSinTrkCosV0A_%d", ic), "; p_{T} (GeV/c); #LT sin*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkCosV0A[ic]);

                fCosTrkSinV0A[ic] = new TProfile(Form("fCosTrkSinV0A_%d", ic), "; p_{T} (GeV/c); #LT cos*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkSinV0A[ic]);


                fSinTrkCosV0C[ic] = new TProfile(Form("fSinTrkCosV0C_%d", ic), "; p_{T} (GeV/c); #LT sin*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkCosV0C[ic]);

                fCosTrkSinV0C[ic] = new TProfile(Form("fCosTrkSinV0C_%d", ic), "; p_{T} (GeV/c); #LT cos*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkSinV0C[ic]);



                fSinTrkSinV0A[ic] = new TProfile(Form("fSinTrkSinV0A_%d", ic), "; p_{T} (GeV/c); #LT sin*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkSinV0A[ic]);

                fCosTrkCosV0A[ic] = new TProfile(Form("fCosTrkCosV0A_%d", ic), "; p_{T} (GeV/c); #LT cos*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkCosV0A[ic]);


                fSinTrkSinV0C[ic] = new TProfile(Form("fSinTrkSinV0C_%d", ic), "; p_{T} (GeV/c); #LT sin*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkSinV0C[ic]);

                fCosTrkCosV0C[ic] = new TProfile(Form("fCosTrkCosV0C_%d", ic), "; p_{T} (GeV/c); #LT cos*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkCosV0C[ic]);

            }



            for (Int_t ib = 0; ib < fNPtBins; ib++){

                fInvMassK0[ib][ic] = new TH1D(Form("fInvMassK0_%d_%d", ib, ic), "; M_{K^{0}}; Counts", 40, 0.4, 0.6);
                fListOfObjects->Add(fInvMassK0[ib][ic]);

                fVnK0A[ib][ic] = new TProfile(Form("fVnK0A_%d_%d", ib, ic), "; M_{K^{0}}; v_{n}", 20, 0.4, 0.6);
                fListOfObjects->Add(fVnK0A[ib][ic]);

                fVnK0C[ib][ic] = new TProfile(Form("fVnK0C_%d_%d", ib, ic), "; M_{K^{0}}; v_{n}", 20, 0.4, 0.6);
                fListOfObjects->Add(fVnK0C[ib][ic]);


                fInvMassL[ib][ic] = new TH1D(Form("fInvMassL_%d_%d", ib, ic), "; M_{#lambda}; Counts", 40, 1.07, 1.17);
                fListOfObjects->Add(fInvMassL[ib][ic]);

                fVnLA[ib][ic] = new TProfile(Form("fVnLA_%d_%d", ib, ic), "; M_{#lambda}; v_{n}", 20, 1.07, 1.17);
                fListOfObjects->Add(fVnLA[ib][ic]);

                fVnLC[ib][ic] = new TProfile(Form("fVnLC_%d_%d", ib, ic), "; M_{#lambda}; v_{n}", 20, 1.07, 1.17);
                fListOfObjects->Add(fVnLC[ib][ic]);


                if (ic == 0){

                    fInvMassK0B1C[ib][ic] = new TH1D(Form("fInvMassK0B1C_%d_%d", ib, ic), "; M_{K^{0}}; Counts", 40, 0.4, 0.6);
                    fListOfObjects->Add(fInvMassK0B1C[ib][ic]);

                    fVnK0AB1C[ib][ic] = new TProfile(Form("fVnK0AB1C_%d_%d", ib, ic), "; M_{K^{0}}; v_{n}", 20, 0.4, 0.6);
                    fListOfObjects->Add(fVnK0AB1C[ib][ic]);

                    fVnK0CB1C[ib][ic] = new TProfile(Form("fVnK0CB1C_%d_%d", ib, ic), "; M_{K^{0}}; v_{n}", 20, 0.4, 0.6);
                    fListOfObjects->Add(fVnK0CB1C[ib][ic]);


                    fInvMassLB1C[ib][ic] = new TH1D(Form("fInvMassLB1C_%d_%d", ib, ic), "; M_{#lambda}; Counts", 40, 1.07, 1.17);
                    fListOfObjects->Add(fInvMassLB1C[ib][ic]);

                    fVnLAB1C[ib][ic] = new TProfile(Form("fVnLAB1C_%d_%d", ib, ic), "; M_{#lambda}; v_{n}", 20, 1.07, 1.17);
                    fListOfObjects->Add(fVnLAB1C[ib][ic]);

                    fVnLCB1C[ib][ic] = new TProfile(Form("fVnLCB1C_%d_%d", ib, ic), "; M_{#lambda}; v_{n}", 20, 1.07, 1.17);
                    fListOfObjects->Add(fVnLCB1C[ib][ic]);

                }

            }

        }

    } else {

        fV0AV0CvnB1C = new TProfile("fV0AV0CvnB1C", "; centrality percentile; V0A-V0C correlations", 90, 0, 90);
        fListOfObjects->Add(fV0AV0CvnB1C);

        fV0ATPCvnB1C = new TProfile("fV0ATPCvnB1C", "; centrality percentile; V0A-TPC correlations", 90, 0, 90);
        fListOfObjects->Add(fV0ATPCvnB1C);

        fV0CTPCvnB1C = new TProfile("fV0CTPCvnB1C", "; centrality percentile; V0C-TPC correlations", 90, 0, 90);
        fListOfObjects->Add(fV0CTPCvnB1C);


        fV0AV0CvnB1Csq = new TProfile("fV0AV0CvnB1Csq", "; centrality percentile; V0A-V0C correlations", 90, 0, 90);
        fListOfObjects->Add(fV0AV0CvnB1Csq);

        fV0ATPCvnB1Csq = new TProfile("fV0ATPCvnB1Csq", "; centrality percentile; V0A-TPC correlations", 90, 0, 90);
        fListOfObjects->Add(fV0ATPCvnB1Csq);

        fV0CTPCvnB1Csq = new TProfile("fV0CTPCvnB1Csq", "; centrality percentile; V0C-TPC correlations", 90, 0, 90);
        fListOfObjects->Add(fV0CTPCvnB1Csq);


        for (Int_t ica = 0; ica < 90; ica++){

            fVnAllAB1C[ica] = new TProfile(Form("fVnAllAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
            fListOfObjects->Add(fVnAllAB1C[ica]);

            fVnPiAB1C[ica] = new TProfile(Form("fVnPiAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnPiAB1C[ica]);

            fVnKAB1C[ica] = new TProfile(Form("fVnKAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnKAB1C[ica]);

            fVnAntiPAB1C[ica] = new TProfile(Form("fVnAntiPAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnAntiPAB1C[ica]);

            fVnPihighPtAB1C[ica] = new TProfile(Form("fVnPihighPtAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPihighPtAB1C[ica]);

            fVnPhighPtAB1C[ica] = new TProfile(Form("fVnPhighPtAB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPhighPtAB1C[ica]);



            fVnAllCB1C[ica] = new TProfile(Form("fVnAllCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
            fListOfObjects->Add(fVnAllCB1C[ica]);

            fVnPiCB1C[ica] = new TProfile(Form("fVnPiCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnPiCB1C[ica]);

            fVnKCB1C[ica] = new TProfile(Form("fVnKCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnKCB1C[ica]);

            fVnAntiPCB1C[ica] = new TProfile(Form("fVnAntiPCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBL, ptBinsL);
            fListOfObjects->Add(fVnAntiPCB1C[ica]);

            fVnPihighPtCB1C[ica] = new TProfile(Form("fVnPihighPtCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPihighPtCB1C[ica]);

            fVnPhighPtCB1C[ica] = new TProfile(Form("fVnPhighPtCB1C_%d", ica), "; p_{T} (GeV/c); v_{n}", nPtBH, ptBinsH);
            fListOfObjects->Add(fVnPhighPtCB1C[ica]);



            if (fQAV0){

                fPsiAB1C[ica] = new TH1D(Form("fPsiAB1C_%d", ica), "; #Psi_{n} (V0A); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiAB1C[ica]);

                fPsiCB1C[ica] = new TH1D(Form("fPsiCB1C_%d", ica), "; #Psi_{n} (V0C); Counts", 80, -2., 2.);
                fListOfObjects->Add(fPsiCB1C[ica]);

                fPsiAvsPsiCB1C[ica] = new TH2D(Form("fPsiAvsPsiCB1C_%d", ica), "; #Psi_{n} (V0A); #Psi_{n} (V0C)", 36, -1.8, 1.8, 36, -1.8, 1.8);
                fListOfObjects->Add(fPsiAvsPsiCB1C[ica]);


                fSinTrkCosV0AB1C[ica] = new TProfile(Form("fSinTrkCosV0AB1C_%d", ica), "; p_{T} (GeV/c); #LT sin*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkCosV0AB1C[ica]);

                fCosTrkSinV0AB1C[ica] = new TProfile(Form("fCosTrkSinV0AB1C_%d", ica), "; p_{T} (GeV/c); #LT cos*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkSinV0AB1C[ica]);


                fSinTrkCosV0CB1C[ica] = new TProfile(Form("fSinTrkCosV0CB1C_%d", ica), "; p_{T} (GeV/c); #LT sin*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkCosV0CB1C[ica]);

                fCosTrkSinV0CB1C[ica] = new TProfile(Form("fCosTrkSinV0CB1C_%d", ica), "; p_{T} (GeV/c); #LT cos*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkSinV0CB1C[ica]);



                fSinTrkSinV0AB1C[ica] = new TProfile(Form("fSinTrkSinV0AB1C_%d", ica), "; p_{T} (GeV/c); #LT sin*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkSinV0AB1C[ica]);

                fCosTrkCosV0AB1C[ica] = new TProfile(Form("fCosTrkCosV0AB1C_%d", ica), "; p_{T} (GeV/c); #LT cos*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkCosV0AB1C[ica]);


                fSinTrkSinV0CB1C[ica] = new TProfile(Form("fSinTrkSinV0CB1C_%d", ica), "; p_{T} (GeV/c); #LT sin*sin #GT", nPtB, ptBins);
                fListOfObjects->Add(fSinTrkSinV0CB1C[ica]);

                fCosTrkCosV0CB1C[ica] = new TProfile(Form("fCosTrkCosV0CB1C_%d", ica), "; p_{T} (GeV/c); #LT cos*cos #GT", nPtB, ptBins);
                fListOfObjects->Add(fCosTrkCosV0CB1C[ica]);

            }

        }

    }

    // Post output data.
    PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskPiKpK0Lamba::UserExec(Option_t *) 
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

        if (fQAOutl)
            fVtx->Fill(0);

    } else {

        if (fQAOutl){
            fVtx->Fill(1);
            fVtxBeforeCuts->Fill(zvtx);
        }

        if (TMath::Abs(zvtx) < fVtxCut)
            Analyze(fAOD, zvtx);

    }

    // Post output data.
    PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskPiKpK0Lamba::Analyze(AliAODEvent* aod, Float_t vtxZ)
{  

    //Centrality
    Float_t v0Centr    = -100.;
    Float_t cl1Centr   = -100.;
    Float_t cl0Centr   = -100.;

    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) aod->FindListObject("MultSelection");
    if( !MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return;
    } else {
        v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
        cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
        cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
    }

    if (v0Centr >= 90. || v0Centr < 0)
        return;


    Int_t nITSClsLy0 = aod->GetNumberOfITSClusters(0);
    Int_t nITSClsLy1 = aod->GetNumberOfITSClusters(1);
    Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

    AliAODTracklets* aodTrkl = (AliAODTracklets*)aod->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

    //cout<<nITSClsLy0<<"    "<<nITSClsLy1<<"    "<<nITSCls<<"   "<<nITSTrkls<<endl;


    const Int_t nTracks = aod->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

    Int_t multTrk = 0;
    Int_t multTrkBefC = 0;
    Int_t multTrkTOFBefC = 0;
    Int_t multTPC = 0;
    Int_t multTPCout = 0;

    for (Int_t it = 0; it < nTracks; it++) {

        AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);

        if (!aodTrk){
            delete aodTrk;
            continue;
        }

        if (aodTrk->GetFlags()&AliESDtrack::kTPCout)
            multTPCout++;

        if (aodTrk->TestFilterBit(32)){

            multTrkBefC++;

            if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
                multTrkTOFBefC++;

            if ((TMath::Abs(aodTrk->Eta()) < fEtaCut) && (aodTrk->GetTPCNcls() >= fNoClus) && (aodTrk->Pt() >= fMinPt) && (aodTrk->Pt() < fMaxPt))
                multTrk++;

        }

        if (aodTrk->TestFilterBit(128))
            multTPC++;

    }


    Float_t multTPCn = multTPC;
    Float_t multEsdn = multEsd;
    Float_t multESDTPCDif = multEsdn - multTPCn*3.38;

    AliAODVZERO* aodV0 = aod->GetVZEROData();
    Float_t multV0a = aodV0->GetMTotV0A();
    Float_t multV0c = aodV0->GetMTotV0C();
    Float_t multV0Tot = multV0a + multV0c;


    if (fQAOutl){

        fCentrBef->Fill(v0Centr);

        fMultFBvsMultFBTOFBef->Fill(multTrkBefC, multTrkTOFBefC);
        fMultESDDifvsMultTPCBef->Fill(multTPCn, multESDTPCDif);

        fMultvsCentBef->Fill(v0Centr, multTrk);

        fCenCL0vsV0MBef->Fill(v0Centr, cl0Centr);
        fCenCL1vsV0MBef->Fill(v0Centr, cl1Centr);
        fCenCL0vsCL1Bef->Fill(cl1Centr, cl0Centr);

        fSPclsvsSPDtrksBef->Fill(nITSTrkls, nITSCls);

        fMultV0vsMultTPCoutBef->Fill(multTPCout, multV0Tot);

    }


    //new vertex selection
    const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
    const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)
        return; // one of vertices is missing

    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);

    double dz = vtTrc->GetZ() - vtSPD->GetZ();

    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;

    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
        return; // bad vertexing


    // vertex cut from old selection
    TString vtxTyp = vtSPD->GetTitle();
    Double_t zRes = TMath::Sqrt(covSPD[5]);
    if ((vtxTyp.Contains("vertexer:Z")) && (zRes>0.25) && (vtSPD->GetNContributors() < 20))
        return;




    //pileup cuts
    if (fPileUp){

        if (cl0Centr < fLowCut->Eval(v0Centr))
            return;

        if (cl0Centr > fHighCut->Eval(v0Centr))
            return;


        if (multESDTPCDif > fCutMultESDdif)
            return;

        if (Float_t(multTrkTOFBefC) < fMultTOFLowCut->Eval(Float_t(multTrkBefC)))
            return;

        if (Float_t(multTrkTOFBefC) > fMultTOFHighCut->Eval(Float_t(multTrkBefC)))
            return;

        if (plpMV(aod))
            return;

        Short_t isPileup = aod->IsPileupFromSPD(3);
        if (isPileup != 0)
            return;

        if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0)
            return;

        //new function for 2015 to remove incomplete events
        if (aod->IsIncompleteDAQ())
            return;


        //new cut to remove outliers
        if (Float_t(multTrk) < fMultCentLowCut->Eval(v0Centr))
            return;

    }



    if (fCent == 1)
        v0Centr = cl0Centr;
    else if (fCent == 2)
        v0Centr = cl1Centr;


    Short_t centrCode = -10;
    if ((v0Centr >= 0) && (v0Centr < 5.))
        centrCode = 0;
    else if ((v0Centr >= 5.) && (v0Centr < 10.))
        centrCode = 1;
    else if ((v0Centr >= 10.) && (v0Centr < 20.))
        centrCode = 2;
    else if ((v0Centr >= 20.) && (v0Centr < 30.))
        centrCode = 3;
    else if ((v0Centr >= 30.) && (v0Centr < 40.))
        centrCode = 4;
    else if ((v0Centr >= 40.) && (v0Centr < 50.))
        centrCode = 5;
    else if ((v0Centr >= 50.) && (v0Centr < 60.))
        centrCode = 6;
    else if ((v0Centr >= 60.) && (v0Centr < 70.))
        centrCode = 7;
    else if ((v0Centr >= 70.) && (v0Centr < 80.))
        centrCode = 8;
    else if ((v0Centr >= 80.) && (v0Centr < 90.))
        centrCode = 9;

    if (centrCode < 0)
        return;


    Int_t iCentV0 = Int_t(v0Centr);
    if (iCentV0 >= 90)
        return;


    Int_t iCentSPD = Int_t(cl1Centr);
    if (iCentSPD >= 90)
        return;



    //V0 info
    Double_t Qxan = 0, Qyan = 0;
    Double_t Qxcn = 0, Qycn = 0;
    Double_t sumMa = 0, sumMc = 0;



    for (Int_t iV0 = 0; iV0 < 64; iV0++) {

        if (fRemChV0A){

            if (iV0 == 46)
                continue;

        }


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
                cout<<"Problem with multiplicity in V0C"<<endl;
                continue;
            }

            Qxcn += TMath::Cos(fNHarm*phiV0) * multCorC;
            Qycn += TMath::Sin(fNHarm*phiV0) * multCorC;

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
                cout<<"Problem with multiplicity in V0A"<<endl;
                continue;
            }

            Qxan += TMath::Cos(fNHarm*phiV0) * multCorA;
            Qyan += TMath::Sin(fNHarm*phiV0) * multCorA;

            sumMa = sumMa + multCorA;

        }

    }

    if (sumMa < 0 || sumMc < 0)
        return;


    Double_t QxanCor = Qxan;
    Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCentSPD+1))/fQynsV0A->GetBinContent(iCentSPD+1);
    Double_t QxcnCor = Qxcn;
    Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCentSPD+1))/fQynsV0C->GetBinContent(iCentSPD+1);

    if (fNHarm != 4.){
        QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCentSPD+1))/fQxnsV0A->GetBinContent(iCentSPD+1);
        QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCentSPD+1))/fQxnsV0C->GetBinContent(iCentSPD+1);
    }




    if (fQAOutl){

        fVtxAfterCuts->Fill(vtxZ);
        fCentrAft->Fill(v0Centr);

        fMultFBvsMultFBTOFAft->Fill(multTrkBefC, multTrkTOFBefC);
        fMultESDDifvsMultTPCAft->Fill(multTPCn, multESDTPCDif);

        fMultvsCentAft->Fill(v0Centr, multTrk);

        fCenCL0vsV0MAft->Fill(v0Centr, cl0Centr);
        fCenCL1vsV0MAft->Fill(v0Centr, cl1Centr);
        fCenCL0vsCL1Aft->Fill(cl1Centr, cl0Centr);

        fSPclsvsSPDtrksAft->Fill(nITSTrkls, nITSCls);

        fMultV0vsMultTPCoutAft->Fill(multTPCout, multV0Tot);
    }


    Double_t psina = TMath::ATan2(QyanCor, QxanCor)/fNHarm;
    Double_t psinc = TMath::ATan2(QycnCor, QxcnCor)/fNHarm;

    if (fQAV0){

        if (!fBin1Cent){
            fPsiA[centrCode]->Fill(psina);
            fPsiC[centrCode]->Fill(psinc);
            fPsiAvsPsiC[centrCode]->Fill(psina, psinc);
        } else {
            fPsiAB1C[iCentV0]->Fill(psina);
            fPsiCB1C[iCentV0]->Fill(psinc);
            fPsiAvsPsiCB1C[iCentV0]->Fill(psina, psinc);
        }

        fQxavsV0Bef->Fill(v0Centr, Qxan);
        fQyavsV0Bef->Fill(v0Centr, Qyan);
        fQxcvsV0Bef->Fill(v0Centr, Qxcn);
        fQycvsV0Bef->Fill(v0Centr, Qycn);

        fQxavsVtxZBef->Fill(vtxZ, Qxan);
        fQyavsVtxZBef->Fill(vtxZ, Qyan);
        fQxcvsVtxZBef->Fill(vtxZ, Qxcn);
        fQycvsVtxZBef->Fill(vtxZ, Qycn);

        fQxavsV0Aft->Fill(v0Centr, QxanCor);
        fQyavsV0Aft->Fill(v0Centr, QyanCor);
        fQxcvsV0Aft->Fill(v0Centr, QxcnCor);
        fQycvsV0Aft->Fill(v0Centr, QycnCor);

        fQxavsVtxZAft->Fill(vtxZ, QxanCor);
        fQyavsVtxZAft->Fill(vtxZ, QyanCor);
        fQxcvsVtxZAft->Fill(vtxZ, QxcnCor);
        fQycvsVtxZAft->Fill(vtxZ, QycnCor);

    }



    Double_t Qxtn = 0, Qytn = 0;
    Int_t sumMt = 0;

    for (Int_t it1 = 0; it1 < nTracks; it1++) {

        AliAODTrack* aodTrk1 = (AliAODTrack*)aod->GetTrack(it1);

        if (!aodTrk1){
            delete aodTrk1;
            continue;
        }

        if (aodTrk1->TestFilterBit(768) && TMath::Abs(aodTrk1->Eta()) < fEtaCut && aodTrk1->GetTPCNcls() >= 70 && aodTrk1->Pt() >= fMinPt && aodTrk1->Pt() < 5.){

            Qxtn += TMath::Cos(fNHarm*aodTrk1->Phi());
            Qytn += TMath::Sin(fNHarm*aodTrk1->Phi());
            sumMt++;

        }


        if (!(aodTrk1->TestFilterBit(fFilterbit)))
            continue;

        if ((TMath::Abs(aodTrk1->Eta()) >= fEtaCut) || (aodTrk1->GetTPCNcls() < fNoClus) || (aodTrk1->Pt() < fMinPt) || (aodTrk1->Pt() >= fMaxPt))
            continue;


        if (fPhiCut){
            Double_t phimod = aodTrk1->Phi();
            if(aod->GetMagneticField() < 0)   // for negative polarity field
                phimod = TMath::TwoPi() - phimod;
            if(aodTrk1->Charge() < 0) // for negative charge
                phimod = TMath::TwoPi() - phimod;
            if(phimod < 0)
                cout << "Warning!!!!! phi < 0: " << phimod << endl;

            phimod += TMath::Pi()/18.0; // to center gap in the middle
            phimod = fmod(phimod, TMath::Pi()/9.0);
            if(phimod < fPhiCutHigh->Eval(aodTrk1->Pt()) && phimod > fPhiCutLow->Eval(aodTrk1->Pt()))
                continue; // reject track            
        }


        if (fRemPhiReg){
            if (aodTrk1->Phi() > 1.7 && aodTrk1->Phi() < 2.7)
                continue;
        }


        if (fPileUpTOF){
            if(TMath::Abs(aodTrk1->GetTOFsignalDz())>10.)
                continue;
            if(aodTrk1->GetTOFsignal() < 12000.)
                continue;
            if(aodTrk1->GetTOFsignal() > 25000.)
                continue;
        }


        if (fQPos){
            if (aodTrk1->Charge() < 0)
                continue;
        }


        if (fQNeg){
            if (aodTrk1->Charge() > 0)
                continue;
        }


        if (fEtaRange == 0){
            if (aodTrk1->Eta() >= 0)
                continue;
        }

        if (fEtaRange == 1){
            if (aodTrk1->Eta() < 0)
                continue;
        }


        if (fCrsRowsFrcShCls){

            Float_t nrowscr = aodTrk1->GetTPCNCrossedRows();

            Float_t clsFind = aodTrk1->GetTPCNclsF();
            if (clsFind <= 0)
                continue;

            Float_t ratCrRowsClsFind = nrowscr/clsFind;

            if (nrowscr < 120)
                continue;

            if (ratCrRowsClsFind < 0.9)
                continue;

        }


        Double_t vnSPV0A = TMath::Cos(fNHarm*aodTrk1->Phi())*QxanCor + TMath::Sin(fNHarm*aodTrk1->Phi())*QyanCor;
        Double_t vnSPV0C = TMath::Cos(fNHarm*aodTrk1->Phi())*QxcnCor + TMath::Sin(fNHarm*aodTrk1->Phi())*QycnCor;


        if (fNHarm == 2. && flagPsi42A){
            vnSPV0A = TMath::Cos(4.*aodTrk1->Phi())*(QxanCor*QxanCor - QyanCor*QyanCor) + TMath::Sin(4.*aodTrk1->Phi())*2.*QxanCor*QyanCor;
            vnSPV0C = TMath::Cos(4.*aodTrk1->Phi())*(QxcnCor*QxcnCor - QycnCor*QycnCor) + TMath::Sin(4.*aodTrk1->Phi())*2.*QxcnCor*QycnCor;
        }


        if (!fBin1Cent){

            fVnAllA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
            fVnAllC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);

        } else {

            fVnAllAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
            fVnAllCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);

        }


        if (fQAV0){

            Double_t sinTRKcosV0A = TMath::Sin(fNHarm*aodTrk1->Phi())*QxanCor;
            Double_t cosTRKsinV0A = TMath::Cos(fNHarm*aodTrk1->Phi())*QyanCor;

            Double_t sinTRKcosV0C = TMath::Sin(fNHarm*aodTrk1->Phi())*QxcnCor;
            Double_t cosTRKsinV0C = TMath::Cos(fNHarm*aodTrk1->Phi())*QycnCor;

            Double_t sinTRKsinV0A = TMath::Sin(fNHarm*aodTrk1->Phi())*QyanCor;
            Double_t cosTRKcosV0A = TMath::Cos(fNHarm*aodTrk1->Phi())*QxanCor;

            Double_t sinTRKsinV0C = TMath::Sin(fNHarm*aodTrk1->Phi())*QycnCor;
            Double_t cosTRKcosV0C = TMath::Cos(fNHarm*aodTrk1->Phi())*QxcnCor;

            if (!fBin1Cent){

                fSinTrkCosV0A[centrCode]->Fill(aodTrk1->Pt(), sinTRKcosV0A);
                fCosTrkSinV0A[centrCode]->Fill(aodTrk1->Pt(), cosTRKsinV0A);

                fSinTrkCosV0C[centrCode]->Fill(aodTrk1->Pt(), sinTRKcosV0C);
                fCosTrkSinV0C[centrCode]->Fill(aodTrk1->Pt(), cosTRKsinV0C);

                fSinTrkSinV0A[centrCode]->Fill(aodTrk1->Pt(), sinTRKsinV0A);
                fCosTrkCosV0A[centrCode]->Fill(aodTrk1->Pt(), cosTRKcosV0A);

                fSinTrkSinV0C[centrCode]->Fill(aodTrk1->Pt(), sinTRKsinV0C);
                fCosTrkCosV0C[centrCode]->Fill(aodTrk1->Pt(), cosTRKcosV0C);

            } else {

                fSinTrkCosV0AB1C[iCentV0]->Fill(aodTrk1->Pt(), sinTRKcosV0A);
                fCosTrkSinV0AB1C[iCentV0]->Fill(aodTrk1->Pt(), cosTRKsinV0A);

                fSinTrkCosV0CB1C[iCentV0]->Fill(aodTrk1->Pt(), sinTRKcosV0C);
                fCosTrkSinV0CB1C[iCentV0]->Fill(aodTrk1->Pt(), cosTRKsinV0C);

                fSinTrkSinV0AB1C[iCentV0]->Fill(aodTrk1->Pt(), sinTRKsinV0A);
                fCosTrkCosV0AB1C[iCentV0]->Fill(aodTrk1->Pt(), cosTRKcosV0A);

                fSinTrkSinV0CB1C[iCentV0]->Fill(aodTrk1->Pt(), sinTRKsinV0C);
                fCosTrkCosV0CB1C[iCentV0]->Fill(aodTrk1->Pt(), cosTRKcosV0C);

            }

        }

        if (fQA){

            if (!fBin1Cent){

                Double_t allQA[5] = {Double_t(centrCode),
                    aodTrk1->Pt(),
                    Double_t(aodTrk1->Charge()),
                    aodTrk1->Eta(),
                    aodTrk1->Phi()};
                fAllQA->Fill(allQA);

            } else {

                Double_t allQAB1C[5] = {Double_t(iCentV0),
                    aodTrk1->Pt(),
                    Double_t(aodTrk1->Charge()),
                    aodTrk1->Eta(),
                    aodTrk1->Phi()};
                fAllQAB1C->Fill(allQAB1C);
            }

        }



        if (aodTrk1->GetTPCsignalN() < fNoClusPid)
            continue;


        Double_t piDEDX = fPIDResponse->GetTPCResponse().GetExpectedSignal(aodTrk1, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, kTRUE);
        Double_t Dpi = aodTrk1->GetTPCsignal() - piDEDX;


        //pi high pT
        if ((Dpi > fMinPiCut) && (Dpi < fMaxPiCut) && (aodTrk1->Pt() >= 3.) && (aodTrk1->GetTPCmomentum() > 3.)){

            Double_t rapPiHPt = GetRapidity(0.139570, aodTrk1->Pt(), aodTrk1->Eta());

            if (TMath::Abs(rapPiHPt) < 0.5){

                if (!fBin1Cent){
                    fVnPihighPtA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPihighPtC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);

                } else {
                    fVnPihighPtAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPihighPtCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);
                }

            }

        }

        //p high pT
        if ((Dpi > fMinPCut) && (Dpi < fMaxPCut) && (aodTrk1->Pt() >= 3.) && (aodTrk1->GetTPCmomentum() > 3.)){

            Double_t rapPHPt = GetRapidity(0.938272, aodTrk1->Pt(), aodTrk1->Eta());

            if (TMath::Abs(rapPHPt) < 0.5){

                if (!fBin1Cent){
                    fVnPhighPtA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPhighPtC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);

                } else {
                    fVnPhighPtAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPhighPtCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);
                }

            }

        }



        Double_t nSigPiTPC = fPIDResponse->NumberOfSigmasTPC(aodTrk1, AliPID::kPion);
        Double_t nSigKTPC  = fPIDResponse->NumberOfSigmasTPC(aodTrk1, AliPID::kKaon);
        Double_t nSigPTPC  = fPIDResponse->NumberOfSigmasTPC(aodTrk1, AliPID::kProton);

        Double_t nSigPiTOF = fPIDResponse->NumberOfSigmasTOF(aodTrk1, AliPID::kPion);
        Double_t nSigKTOF  = fPIDResponse->NumberOfSigmasTOF(aodTrk1, AliPID::kKaon);
        Double_t nSigPTOF  = fPIDResponse->NumberOfSigmasTOF(aodTrk1, AliPID::kProton);

        Double_t nSigmaPi = 99999.;
        Double_t nSigmaK  = 99999.;
        Double_t nSigmaP  = 99999.;

        Float_t intL    = aodTrk1->GetIntegratedLength();
        Float_t timeTOF = aodTrk1->GetTOFsignal();
        Double_t betaPiK = -0.05;
        Double_t betaP = -0.05;

        //pi+K
        if ((aodTrk1->Pt() >= 0.4) && (aodTrk1->GetFlags()&AliESDtrack::kTOFout) && (aodTrk1->GetFlags()&AliESDtrack::kTIME) && (!(aodTrk1->GetFlags()&AliESDtrack::kTOFmismatch)) && (intL > 0) && (timeTOF > 0)){

            betaPiK = intL/2.99792458e-2/timeTOF;

            nSigmaPi = TMath::Sqrt(nSigPiTPC*nSigPiTPC + nSigPiTOF*nSigPiTOF);
            nSigmaK = TMath::Sqrt(nSigKTPC*nSigKTPC + nSigKTOF*nSigKTOF);

        }
        if (aodTrk1->Pt() < 0.4){
            if (aodTrk1->GetTPCsignal() <= 60)
                nSigmaPi = TMath::Abs(nSigPiTPC);
            if (aodTrk1->GetTPCsignal() >= 110)
                nSigmaK = TMath::Abs(nSigKTPC);
        }

        //p
        if ((aodTrk1->Pt() >= 0.5) && (aodTrk1->GetFlags()&AliESDtrack::kTOFout) && (aodTrk1->GetFlags()&AliESDtrack::kTIME) && (!(aodTrk1->GetFlags()&AliESDtrack::kTOFmismatch)) && (intL > 0) && (timeTOF > 0)){

            betaP = intL/2.99792458e-2/timeTOF;

            nSigmaP = TMath::Sqrt(nSigPTPC*nSigPTPC + nSigPTOF*nSigPTOF);

        }
        if (aodTrk1->Pt() < 0.5) {
            if (aodTrk1->GetTPCsignal() >= 110)
                nSigmaP = TMath::Abs(nSigPTPC);
        }


        //exclusive PID
        if (fExclPID){
            if ( (nSigmaPi < fNsigCut && nSigmaK < fNsigCut) || (nSigmaPi < fNsigCut && nSigmaP < fNsigCut) || (nSigmaK < fNsigCut && nSigmaP < fNsigCut) )
                continue;
        }


        Short_t minSigma = FindMinNSigma(nSigmaPi, nSigmaK, nSigmaP);

        if (minSigma == 0)
            continue;

        if ((nSigmaK == nSigmaPi) && (nSigmaK == nSigmaP))
            continue;


        Short_t pidFl = -1;

        //pi
        if (aodTrk1->Pt() < 4. && minSigma == 1 && !GetDoubleCountingPi(nSigmaPi, minSigma) && ((aodTrk1->Pt() >= 0.4 && betaPiK > 0.4) || (aodTrk1->Pt() < 0.4 && betaPiK < 0))){

            Double_t rapPi = GetRapidity(0.139570, aodTrk1->Pt(), aodTrk1->Eta());

            if (TMath::Abs(rapPi) < 0.5){

                if (!fBin1Cent){
                    fVnPiA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPiC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);
                } else {
                    fVnPiAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnPiCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);
                }

                pidFl = 0;
            }

        }

        //K
        if (aodTrk1->Pt() < 4. && minSigma == 2 && !GetDoubleCountingK(nSigmaK, minSigma) && ((aodTrk1->Pt() >= 0.4 && betaPiK > 0.4) || (aodTrk1->Pt() < 0.4 && betaPiK < 0))){

            Double_t rapK = GetRapidity(0.493667, aodTrk1->Pt(), aodTrk1->Eta());

            if (TMath::Abs(rapK) < 0.5){

                if (!fBin1Cent){
                    fVnKA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnKC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);
                } else {
                    fVnKAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
                    fVnKCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);
                }

                pidFl = 1;

            }

        }

        //p
        if (aodTrk1->Pt() < 4. && minSigma == 3 && !GetDoubleCountingP(nSigmaP, minSigma) && ((aodTrk1->Pt() >= 0.5 && betaP > 0.4) || (aodTrk1->Pt() < 0.5 && betaP < 0))){

            if ((aodTrk1->Charge() < 0 && aodTrk1->Pt() < 2.) || (aodTrk1->Pt() >= 2.)){

                Double_t rapP = GetRapidity(0.938272, aodTrk1->Pt(), aodTrk1->Eta());

                if (TMath::Abs(rapP) < 0.5){

                    if (!fBin1Cent){
                        fVnAntiPA[centrCode]->Fill(aodTrk1->Pt(), vnSPV0A);
                        fVnAntiPC[centrCode]->Fill(aodTrk1->Pt(), vnSPV0C);
                    } else {
                        fVnAntiPAB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0A);
                        fVnAntiPCB1C[iCentV0]->Fill(aodTrk1->Pt(), vnSPV0C);
                    }

                    pidFl = 2;

                }

            }

        }


        if (fQA && pidFl >= 0 && aodTrk1->Pt() < 4.){

            Float_t dedx = aodTrk1->GetTPCsignal();

            Double_t beta = 0;
            if (pidFl <  2)
                beta = betaPiK;
            else
                beta = betaP;

            if (!fBin1Cent){
                Double_t pidQA[6] = {Double_t(centrCode),
                    dedx,
                    beta,
                    Double_t(aodTrk1->Charge()),
                    Double_t(pidFl),
                    aodTrk1->P()};
                fPidQA->Fill(pidQA);

            } else {
                Double_t pidQAB1C[6] = {Double_t(iCentV0),
                    dedx,
                    beta,
                    Double_t(aodTrk1->Charge()),
                    Double_t(pidFl),
                    aodTrk1->P()};
                fPidQAB1C->Fill(pidQAB1C);

            }

        }

    }




    const Int_t nv0s = aod->GetNumberOfV0s();

    Double_t  lPrimaryVtxPosition[3];
    aod->GetPrimaryVertex()->GetXYZ(lPrimaryVtxPosition);

    //v0 vn
    for (Int_t jV0 = 0; jV0 < nv0s; jV0++) {

        AliAODv0* v0 = (AliAODv0*) aod->GetV0(jV0);

        if (!v0){
            delete v0;
            continue;
        }


        if (fEtaRange == 0){
            if (v0->Eta() >= 0)
                continue;
        }

        if (fEtaRange == 1){
            if (v0->Eta() < 0)
                continue;
        }


        // V0 cuts
        Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
        Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
        Double_t lV0CosineOfPointingAngle = v0->CosPointingAngle(aod->GetPrimaryVertex());
        Double_t lDcaV0Daughters = v0->DcaV0Daughters();
        Double_t lV0Radius = v0->RadiusV0();

        Double_t  lV0Position[3];
        v0->GetXYZ(lV0Position);
        Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
                TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
                TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2));

        if ( (lDcaPosToPrimVertex < fDCADghtPV) || (lDcaNegToPrimVertex < fDCADghtPV) || (lDcaV0Daughters > fMaxDCADght) || (lV0CosineOfPointingAngle < fCosPA) || (lV0Radius < fMinRad) || (lV0Radius > fMaxRad) || (v0->Pt() < fMinPt) || (v0->Pt() >= fMaxPt) || (TMath::Abs(v0->Eta()) >= fEtaCut) )
            continue;



        if (fArmPodCut){

            Double_t ptArm = v0->PtArmV0();
            Double_t angAlpha = v0->AlphaV0();

            if (ptArm <= 0.2*TMath::Abs(angAlpha))
                continue;

        }



        // Tracks quality cuts
        const AliAODTrack* negTrack = (AliAODTrack *)v0->GetDaughter(1);
        if (!negTrack){
            delete negTrack;
            continue;
        }

        const AliAODTrack* posTrack =(AliAODTrack *)v0->GetDaughter(0);
        if (!posTrack){
            delete posTrack;
            continue;
        }


        if (!posTrack->IsOn(AliAODTrack::kTPCrefit))
            continue;

        if (!negTrack->IsOn(AliAODTrack::kTPCrefit))
            continue;



        if (posTrack->Charge() == negTrack->Charge()){
            //cout<< "like sign, continue"<< endl;
            continue;
        }


        if ((posTrack->GetTPCNcls() < fNoClus) || (negTrack->GetTPCNcls() < fNoClus))
            continue;


        Float_t nCrossedRowsTPCPos = posTrack->GetTPCNCrossedRows();
        if (nCrossedRowsTPCPos < fNoClus)
            continue;
        Int_t findablePos = posTrack->GetTPCNclsF();
        if (findablePos <= 0)
            continue;
        if (nCrossedRowsTPCPos/findablePos < fNcrFind)
            continue;


        Float_t nCrossedRowsTPCNeg = negTrack->GetTPCNCrossedRows();
        if (nCrossedRowsTPCNeg < fNoClus)
            continue;
        Int_t findableNeg = negTrack->GetTPCNclsF();
        if (findableNeg <= 0)
            continue;
        if (nCrossedRowsTPCNeg/findableNeg < fNcrFind)
            continue;


        if (TMath::Abs(posTrack->Eta()) > fEtaCut || TMath::Abs(negTrack->Eta()) > fEtaCut)
            continue;


        if (fMinPtDght){
            if (posTrack->Pt() < fMinPt || negTrack->Pt() < fMinPt)
                continue;
        }



        if (negTrack->GetTPCsignalN() < fNoClusPid || posTrack->GetTPCsignalN() < fNoClusPid)
            continue;

        Double_t nSigmaPipos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion));
        Double_t nSigmaPineg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion));

        Double_t nSigmaPpos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton));
        Double_t nSigmaPneg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kProton));



        Int_t iPtV = GetPtBin(v0->Pt());
        if (iPtV < 0)
            continue;


        Double_t lInvMassK0s = v0->MassK0Short();
        Double_t lInvMassL = v0->MassLambda();

        Short_t flagV0 = -1;

        //K0
        if (nSigmaPipos < fNsigCut && nSigmaPineg < fNsigCut && lInvMassK0s > 0.4 && lInvMassK0s < 0.6){

            Double_t rapK0 = GetRapidity(0.497648, v0->Pt(), v0->Eta());

            if (TMath::Abs(rapK0) < 0.5){

                Double_t vnK0V0C = TMath::Cos(fNHarm*v0->Phi())*QxcnCor + TMath::Sin(fNHarm*v0->Phi())*QycnCor;
                Double_t vnK0V0A = TMath::Cos(fNHarm*v0->Phi())*QxanCor + TMath::Sin(fNHarm*v0->Phi())*QyanCor;

                fInvMassK0[iPtV][centrCode]->Fill(lInvMassK0s);
                fVnK0A[iPtV][centrCode]->Fill(lInvMassK0s, vnK0V0A);
                fVnK0C[iPtV][centrCode]->Fill(lInvMassK0s, vnK0V0C);

                if (iCentV0 == 0){

                    fInvMassK0B1C[iPtV][iCentV0]->Fill(lInvMassK0s);
                    fVnK0AB1C[iPtV][iCentV0]->Fill(lInvMassK0s, vnK0V0A);
                    fVnK0CB1C[iPtV][iCentV0]->Fill(lInvMassK0s, vnK0V0C);

                }

                flagV0 = 0;

            }

        }


        //Lambda
        if ( ((nSigmaPpos < fNsigCut && nSigmaPineg < fNsigCut) || (nSigmaPneg < fNsigCut || nSigmaPipos < fNsigCut)) && lInvMassL > 1.07 && lInvMassL < 1.17){

            Double_t rapL = GetRapidity(1.115683, v0->Pt(), v0->Eta());

            if (TMath::Abs(rapL) < 0.5){

                Double_t vnLV0C = TMath::Cos(fNHarm*v0->Phi())*QxcnCor + TMath::Sin(fNHarm*v0->Phi())*QycnCor;
                Double_t vnLV0A = TMath::Cos(fNHarm*v0->Phi())*QxanCor + TMath::Sin(fNHarm*v0->Phi())*QyanCor;

                fInvMassL[iPtV][centrCode]->Fill(lInvMassL);
                fVnLA[iPtV][centrCode]->Fill(lInvMassL, vnLV0A);
                fVnLC[iPtV][centrCode]->Fill(lInvMassL, vnLV0C);

                if (iCentV0 == 0){

                    fInvMassLB1C[iPtV][iCentV0]->Fill(lInvMassL);
                    fVnLAB1C[iPtV][iCentV0]->Fill(lInvMassL, vnLV0A);
                    fVnLCB1C[iPtV][iCentV0]->Fill(lInvMassL, vnLV0C);

                }

                flagV0 = 1;

            }

        }

        if (fQA && flagV0 >= 0){

            Double_t v0Data[8] = {Double_t(centrCode),
                v0->Pt(),
                lDcaPosToPrimVertex,
                lDcaNegToPrimVertex,
                lV0CosineOfPointingAngle,
                lDcaV0Daughters,
                lV0Radius,
                (Double_t)flagV0};
            fV0QA->Fill(v0Data);

        }


    }


    Double_t corV0AV0Cvn = QxanCor*QxcnCor + QyanCor*QycnCor;
    Double_t corV0ATPCvn = QxanCor*Qxtn + QyanCor*Qytn;
    Double_t corV0CTPCvn = QxcnCor*Qxtn + QycnCor*Qytn;

    if (!fBin1Cent){

        fV0AV0Cvn->Fill(v0Centr, corV0AV0Cvn);
        fV0ATPCvn->Fill(v0Centr, corV0ATPCvn);
        fV0CTPCvn->Fill(v0Centr, corV0CTPCvn);

        fV0AV0Cvnsq->Fill(v0Centr, corV0AV0Cvn*corV0AV0Cvn);
        fV0ATPCvnsq->Fill(v0Centr, corV0ATPCvn*corV0ATPCvn);
        fV0CTPCvnsq->Fill(v0Centr, corV0CTPCvn*corV0CTPCvn);

    } else {

        fV0AV0CvnB1C->Fill(v0Centr, corV0AV0Cvn);
        fV0ATPCvnB1C->Fill(v0Centr, corV0ATPCvn);
        fV0CTPCvnB1C->Fill(v0Centr, corV0CTPCvn);

        fV0AV0CvnB1Csq->Fill(v0Centr, corV0AV0Cvn*corV0AV0Cvn);
        fV0ATPCvnB1Csq->Fill(v0Centr, corV0ATPCvn*corV0ATPCvn);
        fV0CTPCvnB1Csq->Fill(v0Centr, corV0CTPCvn*corV0CTPCvn);
    }

}

//_____________________________________________________________________________
void AliAnalysisTaskPiKpK0Lamba::OpenInfoCalbration(Int_t run)
{

    if (!gGrid) {
        TGrid::Connect("alien://");
    }


    TFile* foadb = 0;
    if (!fRemChV0A)
        foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIR.root");
    else
        foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIRNoCh46V0A.root");


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



    AliOADBContainer* contQxnam = 0;
    if (fNHarm == 2.)
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
    else
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");

    if(!contQxnam){
        printf("OADB object fqxanm is not available in the file\n");
        return;
    }
    if(!(contQxnam->GetObject(run))){
        printf("OADB object fqxanm is not available for run %i\n", run);
        return;
    }
    fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));



    AliOADBContainer* contQynam = 0;
    if (fNHarm == 2.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
    else if (fNHarm == 3.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
    else if (fNHarm == 4.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya4m");

    if(!contQynam){
        printf("OADB object fqyanm is not available in the file\n");
        return;
    }
    if(!(contQynam->GetObject(run))){
        printf("OADB object fqyanm is not available for run %i\n", run);
        return;
    }
    fQynmV0A = ((TH1D*) contQynam->GetObject(run));



    AliOADBContainer* contQxnas = 0;
    if (fNHarm == 2.)
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
    else
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");

    if(!contQxnas){
        printf("OADB object fqxans is not available in the file\n");
        return;
    }
    if(!(contQxnas->GetObject(run))){
        printf("OADB object fqxans is not available for run %i\n", run);
        return;
    }
    fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));



    AliOADBContainer* contQynas = 0;
    if (fNHarm == 2.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
    else if (fNHarm == 3.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
    else if (fNHarm == 4.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya4s");

    if(!contQynas){
        printf("OADB object fqyans is not available in the file\n");
        return;
    }
    if(!(contQynas->GetObject(run))){
        printf("OADB object fqyans is not available for run %i\n", run);
        return;
    }
    fQynsV0A = ((TH1D*) contQynas->GetObject(run));



    AliOADBContainer* contQxncm = 0;
    if (fNHarm == 2.)
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
    else
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");

    if(!contQxncm){
        printf("OADB object fqxcnm is not available in the file\n");
        return;
    }
    if(!(contQxncm->GetObject(run))){
        printf("OADB object fqxcnm is not available for run %i\n", run);
        return;
    }
    fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));



    AliOADBContainer* contQyncm = 0;
    if (fNHarm == 2.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
    else if (fNHarm == 3.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
    else if (fNHarm == 4.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");

    if(!contQyncm){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
    if(!(contQyncm->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
    }
    fQynmV0C = ((TH1D*) contQyncm->GetObject(run));



    AliOADBContainer* contQxncs = 0;
    if (fNHarm == 2.)
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
    else
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");

    if(!contQxncs){
        printf("OADB object fqxc2s is not available in the file\n");
        return;
    }
    if(!(contQxncs->GetObject(run))){
        printf("OADB object fqxc2s is not available for run %i\n", run);
        return;
    }
    fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));



    AliOADBContainer* contQyncs = 0;
    if (fNHarm == 2.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
    else if (fNHarm == 3.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
    else if (fNHarm == 4.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");

    if(!contQyncs){
        printf("OADB object fqycnm is not available in the file\n");
        return;
    }
    if(!(contQyncs->GetObject(run))){
        printf("OADB object fqycns is not available for run %i\n", run);
        return;
    }
    fQynsV0C = ((TH1D*) contQyncs->GetObject(run));

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskPiKpK0Lamba::GetVertex(AliAODEvent* aod) const
{

    Float_t vtxz = -999;

    const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
    if (!trkVtx || trkVtx->GetNContributors()<=0) 
        return vtxz;

    vtxz = trkVtx->GetZ();

    return vtxz;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskPiKpK0Lamba::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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
Bool_t AliAnalysisTaskPiKpK0Lamba::plpMV(const AliVEvent *event)
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

//---------------------------------------------------
Short_t AliAnalysisTaskPiKpK0Lamba::FindMinNSigma(Double_t nSpi, Double_t nSk, Double_t nSp)
{

    Short_t kPID = 0;

    if((nSk == nSpi) && (nSk == nSp))
        return kPID;

    if((nSk < nSpi) && (nSk < nSp) && (nSk < fNsigCut))
        kPID = 2;

    if((nSpi < nSk) && (nSpi < nSp) && (nSpi < fNsigCut))
        kPID = 1;

    if((nSp < nSk) && (nSp < nSpi) && (nSp < fNsigCut))
        kPID = 3;

    return kPID;

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPiKpK0Lamba::GetDoubleCountingPi(Double_t nSpi, Short_t minNSigma)
{

    Bool_t piDC = kFALSE;

    if (nSpi < fNsigCut && minNSigma != 1)
        piDC = kTRUE;

    return piDC;

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPiKpK0Lamba::GetDoubleCountingK(Double_t nSk, Short_t minNSigma)
{

    Bool_t kDC = kFALSE;

    if (nSk < fNsigCut && minNSigma != 2)
        kDC = kTRUE;

    return kDC;

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPiKpK0Lamba::GetDoubleCountingP(Double_t nSp, Short_t minNSigma)
{

    Bool_t pDC = kFALSE;

    if (nSp < fNsigCut && minNSigma != 3)
        pDC = kTRUE;

    return pDC;

}

//__________________________________________________________
Double_t AliAnalysisTaskPiKpK0Lamba::GetRapidity(Double_t mass, Double_t Pt, Double_t Eta)
{

    Double_t rapid = TMath::Log( (TMath::Sqrt(mass*mass + Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta)) + Pt*TMath::SinH(Eta)) / TMath::Sqrt(mass*mass + Pt*Pt) );
    return rapid;

}

//____________________________________________________________________
Int_t AliAnalysisTaskPiKpK0Lamba::GetPtBin(Double_t valPt) const
{

    Int_t ptBin = -10;

    for (Int_t iPt = 0; iPt < fNPtBins; iPt++){

        if(valPt >= fPtBins[iPt] && valPt < fPtBins[iPt+1])
            ptBin = iPt;

    }

    return ptBin;

}

//_____________________________________________________________________________
void AliAnalysisTaskPiKpK0Lamba::Terminate(Option_t *)
{ 
    // Terminate loop
    Printf("Terminate()");
}
