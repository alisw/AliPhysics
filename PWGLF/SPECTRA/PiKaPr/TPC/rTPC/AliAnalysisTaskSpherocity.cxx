/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 *************************************************************************/

#include "AliAnalysisTaskSpherocity.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include "AliAnalysisTask.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
#include "AliAnalysisUtils.h"
//#include <AliSpherocityUtils.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include <AliCentrality.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliVParticle.h>
#include <AliPID.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>


// STL includes
#include <iostream>
using namespace std;


//
// Responsible:
// Antonio Ortiz (Lund)
// Peter Christiansen (Lund)
//


const Double_t AliAnalysisTaskSpherocity::fgkClight = 2.99792458e-2;
Float_t MAGF                    = 1;
const Int_t nHists              = 4;
const Int_t nCent               = 11;
const Double_t CentMin[nCent]   = {0.0,1.0,5.0 ,10.0,15.0,20.0,30.0,40.0,50.0,70.0,0.0};
const Double_t CentMax[nCent]   = {1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,100.0,100.0};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps
//const Int_t MultTrksMin[nCent]  = {51,41,36,31,26,21,16,11,6,0,0};
//const Int_t MultTrksMax[nCent]  = {100,50,40,35,30,25,20,15,10,5,100};

const Char_t *So[3]         = {"Jetty", "Isotropic", "Reference"};

const Double_t aPos[nCent]      = {49.9799,49.9659,0,0,0,0,0,0,0,0,0};
const Double_t bPos[nCent]      = {2.99619,2.91366,0,0,0,0,0,0,0,0,0};
const Double_t cPos[nCent]      = {-45.718,-45.5994,0,0,0,0,0,0,0,0,0};
const Double_t dPos[nCent]      = {290.013,290.042,0,0,0,0,0,0,0,0,0};
const Double_t ePos[nCent]      = {-1018.42,-1014.49,0,0,0,0,0,0,0,0,0};
const Double_t fPos[nCent]      = {1948.68,1931.84,0,0,0,0,0,0,0,0,0};
const Double_t gPos[nCent]      = {-1864.06,-1839.36,0,0,0,0,0,0,0,0,0};
const Double_t hPos[nCent]      = {692.752,680.421,0,0,0,0,0,0,0,0,0};

const Double_t aNeg[nCent]      = {50.078,50.046,0,0,0,0,0,0,0,0,0};
const Double_t bNeg[nCent]      = {6.67199,6.79992,0,0,0,0,0,0,0,0,0};
const Double_t cNeg[nCent]      = {103.662,109.86,0,0,0,0,0,0,0,0,0};
const Double_t dNeg[nCent]      = {611.034,668.241,0,0,0,0,0,0,0,0,0};
const Double_t eNeg[nCent]      = {1695.63,1916.44,0,0,0,0,0,0,0,0,0};
const Double_t fNeg[nCent]      = {2395.88,2815.04,0,0,0,0,0,0,0,0,0};
const Double_t gNeg[nCent]      = {1669.22,2057.21,0,0,0,0,0,0,0,0,0};
const Double_t hNeg[nCent]      = {455.362,595.391,0,0,0,0,0,0,0,0,0};

const Double_t aPosEl[nCent]    = {80.1263,79.9957,0,0,0,0,0,0,0,0,0};
const Double_t bPosEl[nCent]    = {5.28525,7.03079,0,0,0,0,0,0,0,0,0};
const Double_t cPosEl[nCent]    = {-32.7731,-42.9098,0,0,0,0,0,0,0,0,0};
const Double_t dPosEl[nCent]    = {68.4524,88.7057,0,0,0,0,0,0,0,0,0};
const Double_t ePosEl[nCent]    = {-44.1566,-56.6554,0,0,0,0,0,0,0,0,0};

const Double_t aNegEl[nCent]    = {79.8351,79.7387,0,0,0,0,0,0,0,0,0};
const Double_t bNegEl[nCent]    = {-8.46921,-8.60021,0,0,0,0,0,0,0,0,0};
const Double_t cNegEl[nCent]    = {-44.5947,-44.1718,0,0,0,0,0,0,0,0,0};
const Double_t dNegEl[nCent]    = {-86.2242,-84.4984,0,0,0,0,0,0,0,0,0};
const Double_t eNegEl[nCent]    = {-53.6285,-51.945,0,0,0,0,0,0,0,0,0};

ClassImp(AliAnalysisTaskSpherocity)
AliAnalysisTaskSpherocity::AliAnalysisTaskSpherocity():
AliAnalysisTaskSE(),
fESD(0x0),
fAOD(0x0),
fEventCuts(0x0),
fMC(0x0),
fMCStack(0x0),
fMCArray(0x0),
fPIDResponse(0x0),
fTrackFilter2015PbPb(0x0),
fTrackFilterGolden(0x0),
fTrackFilterTPC(0x0),
fTrackFilter(0x0),
utils(0x0),
//fSpheroUtils(0x0),
fAnalysisType("ESD"),
fAnalysisMC(kFALSE),
fAnalysisPbPb(kFALSE),
fisV0Mestimator(kFALSE),
fRandom(0x0),
fNcl(70),
fEtaCut(0.9),
fCentClass(3),
fSpherocity(2),
fMinCent(0.0),
fMaxCent(100.0),
fDeDxMIPMin(40),
fDeDxMIPMax(60),
fdEdxHigh(200),
fdEdxLow(40),
fJettyCutOff(0.503),
fIsotrCutOff(0.774),
fMinMult(10),
fNrec(0),
fSizeStep(0.1),
fMcProcessType(-999),
fTriggeredEventMB(-999),
fVtxStatus(-999),
fZvtx(-999),
fZvtxMC(-999),
fRun(-999),
fEventId(-999),
fListOfObjects(0),
fEvents(0x0),
hMult(0x0),
fTrcksVsTrklets(0x0),
fdEdxCalibrated(0x0),
fMakePid(0x0),
//fLowPt(0x0),
fLHC16l(0x0),
fcent(0x0),
fcentAfterPrimaries(0x0),
fcentAfterV0s(0x0),
hphiso(0x0),
hetaso(0x0),
hPtTruthVsPtRec(0x0),
hPtTruthVsPtRecJetty(0x0),
hPtTruthVsPtRecIsotr(0x0),
hTruthPhiSo(0x0),
hTruthEtaSo(0x0),
fEtaCalibrationNeg(0x0),
fEtaCalibration(0x0),
felededxfitPos(0x0),
felededxfitNeg(0x0),
fcutDCAxy(0x0),
fcutLow(0x0),
fcutHigh(0x0)

{
    
    for(Int_t j=0; j<nHists; ++j){
        hMIPVsV0M[j]=0;//TH2D, MIP vs V0M Mult. for different eta intervals
        pMIPVsV0M[j]=0;//TProfile, MIP vs V0M Mult. for different eta intervals
        hMIPVsNch[j]=0;//TH2D, MIP vs Nch for different eta intervals
        pMIPVsNch[j]=0;//TProfile, MIP vs Nch for different eta intervals
        
    }
    
    hSOrvsV0M  = 0;
    hSOrvsTrks = 0;
    hRefMultVsRefMultPer = 0;
    for(Int_t i = 0; i<nCent; ++i){
        
        // Histograms for PreCalibration
        
        for(Int_t so = 0; so < 3; ++so){
            hMIPVsEta[i][so]=0;
            pMIPVsEta[i][so]=0;
            hPlateauVsEta[i][so]=0;
            pPlateauVsEta[i][so]=0;
            hMIPVsEtaV0s[i][so]=0;
            pMIPVsEtaV0s[i][so]=0;
            
            hPtAll[i][so]=0;
            hPtpos_TPC[i][so]=0;
            hPtneg_TPC[i][so]=0;
            hPtpos_TOF[i][so]=0;
            hPtneg_TOF[i][so]=0;
            
        }
        
        hPhi[i]=0;
        
        // Histograms for PostCalibration
        
        for(Int_t j=0; j<nHists; ++j){
            for(Int_t so = 0; so < 3; ++so){
                hMIPVsPhi[i][j][so]=0;//TH2D, MIP vs phi for different eta intervals
                pMIPVsPhi[i][j][so]=0;//TProfile, MIP vs phi for different eta intervals
                hPlateauVsPhi[i][j][so]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
                pPlateauVsPhi[i][j][so]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
                
                hPtVsP[i][j][so]=0;//TH2D, Transverse momentum Vs momentum
                hDeDxVsP[i][j][so]=0;//TH2D, DeDx vs P
                
                histPiTof[i][j][so]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
                histPiV0[i][j][so]=0;//TH2D, dE/dx vs p, pi id by V0s
                histPV0[i][j][so]=0;// TH2D, dE/dx vs p, p id by V0s
                histEV0[i][j][so]=0;
                
                hnSigmaPiPos[i][j][so]=0;//TH2D, nSigmas vs Pt Pions
                hnSigmaKPos[i][j][so]=0;//TH2D, nSigmas vs Pt Kaons
                hnSigmaPPos[i][j][so]=0;//TH2D, nSigmas vs Pt Protons
                hnSigmaPiNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiPions
                hnSigmaKNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiKaons
                hnSigmaPNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiProtons
                
                hBetavsPpos[i][j][so]=0;
                hBetavsPneg[i][j][so]=0;
                
                hPtpos_TPC_Eta[i][j][so]=0;
                hPtneg_TPC_Eta[i][j][so]=0;
                hPtpos_TOF_Eta[i][j][so]=0;
                hPtneg_TOF_Eta[i][j][so]=0;
                hPpos_TOF_Eta[i][j][so]=0;
                hPneg_TOF_Eta[i][j][so]=0;
            }
        }
    }
    
    //	hTruthEtaSo = 0;
    //	hTruthPhiSo = 0;
    hSOtvsTrks  = 0;
    hSOtvsTrkst = 0;
    hSOtvsV0M   = 0;
    //default constructor
    for(Int_t cent=0;cent<nCent;++cent){
        hSOtVsSOm[cent]   = 0;
        for(Int_t pid=0;pid<7;++pid){
            for(Int_t so=0;so<3;++so){
                hMcIn[cent][pid][so]     = 0;
                hMcOut[cent][pid][so]    = 0;
                hMcInNeg[cent][pid][so]  = 0;
                hMcInPos[cent][pid][so]  = 0;
                hMcOutNeg[cent][pid][so] = 0;
                hMcOutPos[cent][pid][so] = 0;
            }
        }
    }
    
}

AliAnalysisTaskSpherocity::AliAnalysisTaskSpherocity(const char *name):
AliAnalysisTaskSE(name),
fESD(0x0),
fAOD(0x0),
fEventCuts(0x0),
fMC(0x0),
fMCStack(0x0),
fMCArray(0x0),
fPIDResponse(0x0),
fTrackFilter2015PbPb(0x0),
fTrackFilterGolden(0x0),
fTrackFilterTPC(0x0),
fTrackFilter(0x0),
utils(0x0),
//fSpheroUtils(0x0),
fAnalysisType("ESD"),
fAnalysisMC(kFALSE),
fAnalysisPbPb(kFALSE),
fisV0Mestimator(kFALSE),
fRandom(0x0),
fNcl(70),
fEtaCut(0.9),
fCentClass(3),
fSpherocity(2),
fMinCent(0.0),
fMaxCent(100.0),
fDeDxMIPMin(40),
fDeDxMIPMax(60),
fdEdxHigh(200),
fdEdxLow(40),
fJettyCutOff(0.503),
fIsotrCutOff(0.774),
fMinMult(10),
fNrec(0),
fSizeStep(0.1),
fMcProcessType(-999),
fTriggeredEventMB(-999),
fVtxStatus(-999),
fZvtx(-999),
fZvtxMC(-999),
fRun(-999),
fEventId(-999),
fListOfObjects(0), 
fEvents(0x0),
hMult(0x0),
fTrcksVsTrklets(0x0),
fdEdxCalibrated(0x0),
fMakePid(0x0),
//fLowPt(0x0),
fLHC16l(0x0),
fcent(0x0),
fcentAfterPrimaries(0x0),
fcentAfterV0s(0x0),
hphiso(0x0),
hetaso(0x0),
hPtTruthVsPtRec(0x0),
hPtTruthVsPtRecJetty(0x0),
hPtTruthVsPtRecIsotr(0x0),
hTruthPhiSo(0x0),
hTruthEtaSo(0x0),
fEtaCalibrationNeg(0x0),
fEtaCalibration(0x0),
felededxfitPos(0x0),
felededxfitNeg(0x0),
fcutDCAxy(0x0),
fcutLow(0x0),
fcutHigh(0x0)
{
    
    for(Int_t j=0; j<nHists; ++j){
        hMIPVsV0M[j]=0;//TH2D, MIP vs V0M Mult. for different eta intervals
        pMIPVsV0M[j]=0;//TProfile, MIP vs V0M Mult. for different eta intervals
        hMIPVsNch[j]=0;//TH2D, MIP vs Nch for different eta intervals
        pMIPVsNch[j]=0;//TProfile, MIP vs Nch for different eta intervals
    }
    
    hSOrvsV0M  = 0;
    hSOrvsTrks = 0;
    hRefMultVsRefMultPer = 0;
    
    for(Int_t i = 0; i<nCent; ++i){
        
        for(Int_t so = 0; so < 3; ++so){
            hMIPVsEta[i][so] = 0;
            pMIPVsEta[i][so]=0;
            hPlateauVsEta[i][so]=0;
            pPlateauVsEta[i][so]=0;
            hMIPVsEtaV0s[i][so]=0;
            pMIPVsEtaV0s[i][so]=0;
            
            hPtAll[i][so]=0;
            hPtpos_TPC[i][so]=0;
            hPtneg_TPC[i][so]=0;
            hPtpos_TOF[i][so]=0;
            hPtneg_TOF[i][so]=0;
        }
        
        hPhi[i]=0;
        
        for(Int_t j=0; j<nHists; ++j){
            
            for(Int_t so = 0; so < 3; ++so){
                hMIPVsPhi[i][j][so]=0;//TH2D, MIP vs phi for different eta intervals
                pMIPVsPhi[i][j][so]=0;//TProfile, MIP vs phi for different eta intervals
                hPlateauVsPhi[i][j][so]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
                pPlateauVsPhi[i][j][so]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
                
                hPtVsP[i][j][so]=0;//TH2D, Transverse momentum Vs momentum
                hDeDxVsP[i][j][so]=0;//TH2D, DeDx vs P
                
                histPiTof[i][j][so]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
                histPiV0[i][j][so]=0;//TH2D, dE/dx vs p, pi id by V0s
                histPV0[i][j][so]=0;// TH2D, dE/dx vs p, p id by V0s
                histEV0[i][j][so]=0;
                
                hnSigmaPiPos[i][j][so]=0;//TH2D, nSigmas vs Pt Pions
                hnSigmaKPos[i][j][so]=0;//TH2D, nSigmas vs Pt Kaons
                hnSigmaPPos[i][j][so]=0;//TH2D, nSigmas vs Pt Protons
                hnSigmaPiNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiPions
                hnSigmaKNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiKaons
                hnSigmaPNeg[i][j][so]=0;//TH2D, nSigmas vs Pt AntiProtons
                
                hBetavsPpos[i][j][so]=0;
                hBetavsPneg[i][j][so]=0;
                
                hPtpos_TPC_Eta[i][j][so]=0;
                hPtneg_TPC_Eta[i][j][so]=0;
                hPtpos_TOF_Eta[i][j][so]=0;
                hPtneg_TOF_Eta[i][j][so]=0;
                hPpos_TOF_Eta[i][j][so]=0;
                hPneg_TOF_Eta[i][j][so]=0;
            }
        }
        
    }
    
    //	hTruthEtaSo = 0;
    //	hTruthPhiSo = 0;
    hSOtvsTrks  = 0;
    hSOtvsTrkst = 0;
    hSOtvsV0M   = 0;
    // Default constructor (should not be used)
    for(Int_t cent=0; cent<nCent; ++cent){
        hSOtVsSOm[cent] = 0;
        for(Int_t pid=0; pid<7; ++pid){
            for(Int_t so=0; so<3; ++so){
                hMcIn[cent][pid][so]=0;
                hMcOut[cent][pid][so]=0;
                hMcInNeg[cent][pid][so]=0;
                hMcInPos[cent][pid][so]=0;
                hMcOutNeg[cent][pid][so]=0;
                hMcOutPos[cent][pid][so]=0;
            }
        }
    }
    
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskSpherocity::~AliAnalysisTaskSpherocity() {
    //
    // Destructor
    //
    
}
//______________________________________________________________________________
void AliAnalysisTaskSpherocity::UserCreateOutputObjects()
{
    // This method is called once per worker node
    // Here we define the output: histograms and debug tree if requested
    // We also create the random generator here so it might get different seeds...
    fRandom = new TRandom(0); // 0 means random seed
    
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man){
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
    }
    
    // Definition of trackcuts
    if(!fTrackFilter){
        fTrackFilter = new AliAnalysisFilter("trackFilter2015");
        SetTrackCutsSpherocity(fTrackFilter);
    }
    
    /*if(!fSpheroUtils){
     fSpheroUtils = new AliSpherocityUtils();
     fSpheroUtils->SetTrackFilter(fTrackFilter);
     fSpheroUtils->SetMinMult(10);
     fSpheroUtils->Init();
     
     printf("JettyCutOff  =  %f       IsotrCutOff  =  %f\n",fJettyCutOff,fIsotrCutOff);
     }*/
    
    //OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner(kTRUE);
    
    //
    // Histograms
    //
    
    fEvents = new TH2F( "fEvents", ";Evt. Sel.;Mult. Per.",12,0,12,13,0,13);
    fEvents->GetXaxis()->SetBinLabel(1, "Processed");
    fEvents->GetXaxis()->SetBinLabel(2, "Trigger");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(3, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(4, "DAQ");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(5, "BG");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(6, "INEL>0");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(7, "VtxRes&Proximity");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(9, "non-Selected So");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(10, "Selected So");//NotinVertexcut");
    fEvents->GetYaxis()->SetBinLabel(1,Form("%.2f-%.2f",CentMin[0],CentMax[0]));
    fEvents->GetYaxis()->SetBinLabel(2,Form("%.2f-%.2f",CentMin[1],CentMax[1]));
    fEvents->GetYaxis()->SetBinLabel(3,Form("%.2f-%.2f",CentMin[2],CentMax[2]));
    fEvents->GetYaxis()->SetBinLabel(4,Form("%.2f-%.2f",CentMin[3],CentMax[3]));
    fEvents->GetYaxis()->SetBinLabel(5,Form("%.2f-%.2f",CentMin[4],CentMax[4]));
    fEvents->GetYaxis()->SetBinLabel(6,Form("%.2f-%.2f",CentMin[5],CentMax[5]));
    fEvents->GetYaxis()->SetBinLabel(7,Form("%.2f-%.2f",CentMin[6],CentMax[6]));
    fEvents->GetYaxis()->SetBinLabel(8,Form("%.2f-%.2f",CentMin[7],CentMax[7]));
    fEvents->GetYaxis()->SetBinLabel(9,Form("%.2f-%.2f",CentMin[8],CentMax[8]));
    fEvents->GetYaxis()->SetBinLabel(10,Form("%.2f-%.2f",CentMin[9],CentMax[9]));
    fEvents->GetYaxis()->SetBinLabel(11,"0.0-100.0");
    fListOfObjects->Add(fEvents);
    
    hMult = new TH1F("hMult","Mult of events with SO measured",100,0,100);
    fListOfObjects->Add(hMult);
    
    fTrcksVsTrklets = new TH2F("fTrcksVsTrklets",";SPD_{Tracklets};Global Tracks",100,0,100,100,0,100);
    fListOfObjects->Add(fTrcksVsTrklets);
    
    fcent=new TH1F("fcent","fcent",13,0,13);
    fcentAfterPrimaries =new TH1F("fcentAfterPrimaries","fcentAfterPrimaries",13,0,13);
    fcentAfterV0s =new TH1F("fcentAfterV0s","fcentAfterV0s",13,0,13);
    fListOfObjects->Add(fcent);
    fListOfObjects->Add(fcentAfterPrimaries);
    fListOfObjects->Add(fcentAfterV0s);
    
    const Int_t nDeltaPiBins   = 80;
    const Double_t deltaPiLow  = 20;
    const Double_t deltaPiHigh = 100;
    
    const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
    
    const int nPtBins = 63;
    double ptBins[nPtBins+1] = {
        0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
        0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
        0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
        1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
        3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
        9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00,
        20.00,22.00,24.00,26.00,30.00};
        
    const int nPtBinsV0s = 25;
    double ptBinsV0s[nPtBinsV0s+1] = {
        0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
        1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
        9.0 , 12.0, 15.0, 20.0 };

//    const int nBinsPer = 10;
//    double perBins[nBinsPer+1] = {0.0,1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,100.0};
     printf("==============================\n");
	printf("Running the Fixed Code\n");
	printf("Jetty cut: %f\n",fJettyCutOff);
	printf("Isotropic cut: %f\n",fIsotrCutOff);
     printf("==============================\n");

    const int nBinsdEdx = fdEdxHigh-fdEdxLow;
    double binsdEdx[nBinsdEdx+1];

    for(int i = fdEdxLow; i <= fdEdxHigh; ++i){
        binsdEdx[i-fdEdxLow] = i;
//        printf("edges :: %f\n",binsdEdx[i-fdEdxLow]);
    }
    
    const Char_t* ending[nHists] = {"02", "24", "46", "68"};
    const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };
    
    fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
    fcutDCAxy->SetParameter(0,0.0105);
    fcutDCAxy->SetParameter(1,0.0350);
    fcutDCAxy->SetParameter(2,1.1);
    
    fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
    fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
    
    
    fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
    fEtaCalibration    = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
    
    felededxfitPos     = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
    felededxfitNeg     = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);
    
    hphiso = 0;
    hphiso = new TH1D("hphiso","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
    hphiso->Sumw2();
    fListOfObjects->Add(hphiso);
    
    hetaso = 0;
    hetaso = new TH1D("hetaso","spherocity; #eta; counts",40,-1.0,1.0);
    hetaso->Sumw2();
    fListOfObjects->Add(hetaso);
    
    Int_t nPhiBins = 36;
    
    hSOrvsV0M  = new TH2D("hSOrVsV0M","Measured SO vs V0M Per.;V0M Per.;#it{S}_{O} Reconstructed",100,0,100,1000,0,1);
    hSOrvsV0M->Sumw2();
    fListOfObjects->Add(hSOrvsV0M);
    
    hSOrvsTrks  = new TH2D("hSOrVsTrks","Measured SO vs Measured Ref. Mult. |#eta|<0.8;Reference mult. (|#eta|<0.8);#it{S}_{O} Reconstructed",100, 0, 100, 1000, 0, 1);
    hSOrvsTrks->Sumw2();
    fListOfObjects->Add(hSOrvsTrks);
    
    hRefMultVsRefMultPer = new TH2D("hRefMultVsRefMultPer","Ref Mult. vs Ref. Mult. Per. |#eta|<0.8;Ref. Mult Per.;Ref. Mult",100, 0, 100, 100, 0, 100);
    hRefMultVsRefMultPer->Sumw2();
    fListOfObjects->Add(hRefMultVsRefMultPer);
    
    for(Int_t i = 0; i<nCent; ++i){
        for(Int_t so = 0; so < 3; ++so){
            hMIPVsEta[i][so] = new TH2D(Form("hMIPVsEta%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
            pMIPVsEta[i][so] = new TProfile(Form("pMIPVsEta%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
            hPlateauVsEta[i][so] = new TH2D(Form("hPlateauVsEta%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
            pPlateauVsEta[i][so] = new TProfile(Form("pPlateauVsEta%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);
            hMIPVsEtaV0s[i][so] = new TH2D(Form("hMIPVsEtaV0s%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; dE/dx_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
            pMIPVsEtaV0s[i][so] = new TProfile(Form("pMIPVsEtaV0s%.2f-%.2f-%s",CentMin[i],CentMax[i],So[so]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMin, fDeDxMIPMax);
            hPtAll[i][so] = new TH1D(Form("hPt_%s_%.2f-%.2f",So[so],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtAll[i][so]->Sumw2();
            hPtpos_TPC[i][so] = new TH1D(Form("hPt_pos_TPC_%s_%.2f-%.2f",So[so],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtpos_TPC[i][so]->Sumw2();
            hPtneg_TPC[i][so] = new TH1D(Form("hPt_neg_TPC_%s_%.2f-%.2f",So[so],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtneg_TPC[i][so]->Sumw2();
            hPtpos_TOF[i][so] = new TH1D(Form("hPt_pos_TOF_%s_%.2f-%.2f",So[so],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtpos_TOF[i][so]->Sumw2();
            hPtneg_TOF[i][so] = new TH1D(Form("hPt_neg_TOF_%s_%.2f-%.2f",So[so],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtneg_TOF[i][so]->Sumw2();
        }
        
        hPhi[i] = new TH2D(Form("histPhi%.2f-%.2f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);
        
        for(Int_t j=0; j<nHists; j++) {            
            for(Int_t so = 0; so < 3; ++so){
                hDeDxVsP[i][j][so] = new TH2D(Form("hDeDxVsP%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, fdEdxHigh-fdEdxLow, fdEdxLow, fdEdxHigh);
                hDeDxVsP[i][j][so]->Sumw2();
                hMIPVsPhi[i][j][so] = new TH2D(Form("hMIPVsPhi%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
                hMIPVsPhi[i][j][so]->Sumw2();
                pMIPVsPhi[i][j][so] = new TProfile(Form("pMIPVsPhi%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMin, fDeDxMIPMax);
                pMIPVsPhi[i][j][so]->Sumw2();
                hPlateauVsPhi[i][j][so]  = new TH2D(Form("hPlateauVsPhi%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);
                hPlateauVsPhi[i][j][so]->Sumw2();
                pPlateauVsPhi[i][j][so] = new TProfile(Form("pPlateauVsPhi%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax, 95);
                pPlateauVsPhi[i][j][so]->Sumw2();
                hPtVsP[i][j][so] = new TH2D(Form("hPtVsP%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), ";#it{p} [GeV/c]; #it{p}_{T}", nPtBins, ptBins, nPtBins, ptBins);
                hPtVsP[i][j][so]->Sumw2();
                histPiV0[i][j][so]  = new TH2D(Form("hPiV0%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
                histPiV0[i][j][so]->Sumw2();
                histPV0[i][j][so]   = new TH2D(Form("hPV0%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
                histPV0[i][j][so]->Sumw2();
                histPiTof[i][j][so] = new TH2D(Form("hPiTOF%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
                histPiTof[i][j][so]->Sumw2();
                histEV0[i][j][so]   = new TH2D(Form("hEV0%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
                histEV0[i][j][so]->Sumw2();
                hnSigmaPiPos[i][j][so] = new TH2D(Form("hnSigma_Pion_pos_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p}_{T};nSigmaPiPos",nPtBins,ptBins,40,-10,10);
                hnSigmaPiPos[i][j][so]->Sumw2();
                hnSigmaKPos[i][j][so] = new TH2D(Form("hnSigma_Kaon_pos_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]), ";#it{p}_{T};nSigmaKPos",nPtBins,ptBins,40,-10, 10);
                hnSigmaKPos[i][j][so]->Sumw2();
                hnSigmaPPos[i][j][so] = new TH2D(Form("hnSigma_Proton_pos_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p}_{T};nSigmaPPos",nPtBins,ptBins,40,-10, 10);
                hnSigmaPPos[i][j][so]->Sumw2();
                hnSigmaPiNeg[i][j][so] = new TH2D(Form("hnSigma_Pion_neg_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p}_{T};nSigmaPiNeg",nPtBins,ptBins,40,-10,10);
                hnSigmaPiNeg[i][j][so]->Sumw2();
                hnSigmaKNeg[i][j][so] = new TH2D(Form("hnSigma_Kaon_neg_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p}_{T};nSigmaKNeg",nPtBins,ptBins,40,-10,10);
                hnSigmaKNeg[i][j][so]->Sumw2();
                hnSigmaPNeg[i][j][so] = new TH2D(Form("hnSigma_Proton_neg_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p}_{T};nSigmaPNeg",nPtBins,ptBins,40,-10,10);
                hnSigmaPNeg[i][j][so]->Sumw2();
                
                hBetavsPpos[i][j][so] = new TH2D(Form("hBetavsP_pos_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p} (GeV/#it{c});#beta",nPtBins,ptBins,300,0.3,1.1);
                hBetavsPpos[i][j][so]->Sumw2();
                hBetavsPneg[i][j][so] = new TH2D(Form("hBetavsP_neg_%.2f-%.2f-%s-%s",CentMin[i],CentMax[i],So[so],ending[j]),";#it{p} (GeV/#it{c});#beta",nPtBins,ptBins,300,0.3,1.1);
                hBetavsPneg[i][j][so]->Sumw2();
                
                hPtpos_TPC_Eta[i][j][so] = new TH1D(Form("hPt_pos_TPC_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
                hPtpos_TPC_Eta[i][j][so]->Sumw2();
                hPtneg_TPC_Eta[i][j][so] = new TH1D(Form("hPt_neg_TPC_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
                hPtneg_TPC_Eta[i][j][so]->Sumw2();
                
                hPtpos_TOF_Eta[i][j][so] = new TH1D(Form("hPt_pos_TOF_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
                hPtpos_TOF_Eta[i][j][so]->Sumw2();
                hPtneg_TOF_Eta[i][j][so] = new TH1D(Form("hPt_neg_TOF_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
                hPtneg_TOF_Eta[i][j][so]->Sumw2();

                hPpos_TOF_Eta[i][j][so] = new TH1D(Form("hP_pos_TOF_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p};Counts",nPtBins,ptBins);
                hPpos_TOF_Eta[i][j][so]->Sumw2();
                hPneg_TOF_Eta[i][j][so] = new TH1D(Form("hP_neg_TOF_%s_%.2f-%.2f-%s",So[so],CentMin[i],CentMax[i],ending[j]),";#it{p};Counts",nPtBins,ptBins);
                hPneg_TOF_Eta[i][j][so]->Sumw2();
            }
        }// eta loop
    } // centrality loop
    
    for(Int_t i = 0; i<nHists; ++i ){
        hMIPVsV0M[i] = new TH2D(Form("hMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
        hMIPVsV0M[i]->Sumw2();
        pMIPVsV0M[i] = new TProfile(Form("pMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMin, fDeDxMIPMax);
        pMIPVsV0M[i]->Sumw2();
        hMIPVsNch[i] = new TH2D(Form("hMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
        hMIPVsNch[i]->Sumw2();
        pMIPVsNch[i] = new TProfile(Form("pMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMin, fDeDxMIPMax);
        pMIPVsNch[i]->Sumw2();
    }
    
    if(!fAnalysisMC){
        
        /*for(Int_t j=0; j<nHists; ++j){
         fListOfObjects->Add(hMIPVsV0M[j]);
         fListOfObjects->Add(pMIPVsV0M[j]);
         fListOfObjects->Add(hMIPVsNch[j]);
         fListOfObjects->Add(pMIPVsNch[j]);
         }*/
        
        
        for(Int_t i=0; i<nCent; ++i ){
            
            if(i > 2)continue;
            /* for(Int_t so=0; so<3; ++so){
             fListOfObjects->Add(hMIPVsEta[i][so]);
             fListOfObjects->Add(pMIPVsEta[i][so]);
             fListOfObjects->Add(hPlateauVsEta[i][so]);
             fListOfObjects->Add(pPlateauVsEta[i][so]);
             fListOfObjects->Add(hMIPVsEtaV0s[i][so]);
             fListOfObjects->Add(pMIPVsEtaV0s[i][so]);
             
             for(Int_t j=0; j<nHists; ++j){
             fListOfObjects->Add(hMIPVsPhi[i][j][so]);
             fListOfObjects->Add(pMIPVsPhi[i][j][so]);
             fListOfObjects->Add(hPlateauVsPhi[i][j][so]);
             fListOfObjects->Add(pPlateauVsPhi[i][j][so]);
             }
             }*/
            fListOfObjects->Add(hPhi[i]);
            if(fMakePid){
                for(Int_t so=0; so<3; ++so){
                    
                    fListOfObjects->Add(hPtAll[i][so]);
                    fListOfObjects->Add(hPtneg_TPC[i][so]);
                    fListOfObjects->Add(hPtpos_TPC[i][so]);
                    fListOfObjects->Add(hPtneg_TOF[i][so]);
                    fListOfObjects->Add(hPtpos_TOF[i][so]);
                    
                    for(Int_t j=0; j<nHists; ++j){
                        fListOfObjects->Add(hPtVsP[i][j][so]);
                        //                        fListOfObjects->Add(histPiV0[i][j][so]);
                        //                        fListOfObjects->Add(histPiTof[i][j][so]);
                        //                        fListOfObjects->Add(histEV0[i][j][so]);
                        //                        fListOfObjects->Add(histPV0[i][j][so]);
                        fListOfObjects->Add(hDeDxVsP[i][j][so]);
                        fListOfObjects->Add(hnSigmaPiPos[i][j][so]);
                        fListOfObjects->Add(hnSigmaPiNeg[i][j][so]);
                        fListOfObjects->Add(hnSigmaKPos[i][j][so]);
                        fListOfObjects->Add(hnSigmaKNeg[i][j][so]);
                        fListOfObjects->Add(hnSigmaPPos[i][j][so]);
                        fListOfObjects->Add(hnSigmaPNeg[i][j][so]);
                        fListOfObjects->Add(hBetavsPneg[i][j][so]);
                        fListOfObjects->Add(hBetavsPpos[i][j][so]);
                        fListOfObjects->Add(hPtneg_TPC_Eta[i][j][so]);
                        fListOfObjects->Add(hPtpos_TPC_Eta[i][j][so]);
                        fListOfObjects->Add(hPtneg_TOF_Eta[i][j][so]);
                        fListOfObjects->Add(hPtpos_TOF_Eta[i][j][so]);
                        fListOfObjects->Add(hPneg_TOF_Eta[i][j][so]);
                        fListOfObjects->Add(hPpos_TOF_Eta[i][j][so]);
                    }
                }
            } //	if(MakePID)
        } //	Cent
    } //	!fAnalysisMC
    
    
    else{
        
        for(Int_t cent=0; cent<nCent; cent++) {
            for(Int_t pid=0; pid<7; pid++){
                for(Int_t so=0; so<3; so++){
                    hMcIn[cent][pid][so]=new TH1D(Form("hIn_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]), Form("MC in (pid %s)", Pid[pid]),nPtBins,ptBins);
                    hMcIn[cent][pid][so]->Sumw2();
                    hMcInNeg[cent][pid][so]=new TH1D(Form("hInNeg_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
                    hMcInNeg[cent][pid][so]->Sumw2();
                    hMcInPos[cent][pid][so]=new TH1D(Form("hInPos_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
                    hMcInPos[cent][pid][so]->Sumw2();
                    hMcOut[cent][pid][so]=new TH1D(Form("hMcOut_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
                    hMcOut[cent][pid][so]->Sumw2();
                    hMcOutNeg[cent][pid][so]=new TH1D(Form("hMcOutNeg_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
                    hMcOutNeg[cent][pid][so]->Sumw2();
                    hMcOutPos[cent][pid][so]=new TH1D(Form("hMcOutPos_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
                    hMcOutPos[cent][pid][so]->Sumw2();
                    
                    /*if(cent<3){
                     fListOfObjects->Add(hMcIn[cent][pid][so]);
                     fListOfObjects->Add(hMcInNeg[cent][pid][so]);
                     fListOfObjects->Add(hMcInPos[cent][pid][so]);
                     fListOfObjects->Add(hMcOut[cent][pid][so]);
                     fListOfObjects->Add(hMcOutNeg[cent][pid][so]);
                     fListOfObjects->Add(hMcOutPos[cent][pid][so]);
                     }*/
                }
            }	// pid Eff
            
            fListOfObjects->Add(hPhi[cent]);
            
            hSOtVsSOm[cent] = new TH2D(Form("hSOtVsSOm-%.2f-%.2f",CentMin[cent],CentMax[cent]),";#it{S}_{O} generated;#it{S}_{O} reconstructed",10, 0, 1, 10, 0, 1);
            hSOtVsSOm[cent]->Sumw2();
            fListOfObjects->Add(hSOtVsSOm[cent]);
            
        }	// cent Eff
        
        hPtTruthVsPtRec = new TH2D("hPtTruthVsPtRec","Pt Truth Vs Pt Rec;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
        hPtTruthVsPtRec->Sumw2();
        fListOfObjects->Add(hPtTruthVsPtRec);
        
        hPtTruthVsPtRecJetty = new TH2D("hPtTruthVsPtRecJetty","Pt Truth Vs Pt Rec Jetty Events;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
        hPtTruthVsPtRecJetty->Sumw2();
        fListOfObjects->Add(hPtTruthVsPtRecJetty);
        
        hPtTruthVsPtRecIsotr = new TH2D("hPtTruthVsPtRecIsotr","Pt Truth Vs Pt Rec Isotropic Events;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
        hPtTruthVsPtRecIsotr->Sumw2();
        fListOfObjects->Add(hPtTruthVsPtRecIsotr);
        
        hTruthEtaSo = new TH1D("hTruthEtaSo","spherocity; #eta; counts",40,-1.0,1.0);
        hTruthEtaSo->Sumw2();
        fListOfObjects->Add(hTruthEtaSo);
        
        hTruthPhiSo = new TH1D("hTruthPhiSo","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
        hTruthPhiSo->Sumw2();
        fListOfObjects->Add(hTruthPhiSo);
        
        hSOtvsTrks  = new TH2D("hSOtVsTrks","Truth SO vs Measured Ref. Mult.;Ref. mult. (|#eta|<0.8);#it{S}_{O} Truth ",100, 0, 100, 1000, 0, 1);
        hSOtvsTrks->Sumw2();
        fListOfObjects->Add(hSOtvsTrks);
        
        hSOtvsTrkst  = new TH2D("hSOtVsTrkst","Truth SO vs Truth. Mult.;Truth Mult. (|#eta|<0.8);#it{S}_{O} Truth ",100, 0, 100, 1000, 0, 1);
        hSOtvsTrkst->Sumw2();
        fListOfObjects->Add(hSOtvsTrkst);
        
        hSOtvsV0M  = new TH2D("hSOtVsV0M","Truth SO vs V0M Per.;V0M Per.;#it{S}_{O} Truth",100, 0, 100, 1000, 0, 1);;
        hSOtvsV0M->Sumw2();
        fListOfObjects->Add(hSOtvsV0M);
        
    }
    
    PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskSpherocity::UserExec(Option_t *)
{
    // Main loop
    
    //
    // First we make sure that we have valid input(s)!
    //
    
    AliVEvent *event = InputEvent();
    if (!event) {
        Error("UserExec", "Could not retrieve event");
        return;
    }
    
    if (fAnalysisType == "ESD"){
        fESD = dynamic_cast<AliESDEvent*>(event);
        if(!fESD){
            Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
            this->Dump();
            return;
        }
        
    } else{
        fAOD = dynamic_cast<AliAODEvent*>(event);
        if(!fAOD){
            Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
            this->Dump();
            return;
        }
    }
    
    if (fAnalysisMC){
        if (fAnalysisType == "ESD"){
            fMC = dynamic_cast<AliMCEvent*>(MCEvent());
            if(!fMC){
                Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
                this->Dump();
                return;
            }
            
            fMCStack = fMC->Stack();
            
            if(!fMCStack){
                Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
                this->Dump();
                return;
            }
        } else { // AOD
            
            fMC = dynamic_cast<AliMCEvent*>(MCEvent());
            if(fMC)
                fMC->Dump();
            
            fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
            if(!fMCArray){
                Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
                this->Dump();
                return;
            }
        }
    }
    
    utils = new AliAnalysisUtils();
    if ( !utils ){return;}
    
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
    if(!isINT7selected)
        return;
    
    fEvents->Fill(1.5,10);
    
    //--------------- Event Selection --------------------
    
    float V0MPercentile = -1;
    float RefPercentile = -1;
    int fnRefGlobal = -1;
    int IndxTrksMult = -1;
    int IndxV0MMult = -1;
    
    fnRefGlobal = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8 );
    
    AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
    if(!MultSelection){ return; }
    if( MultSelection-> IsEventSelected() ){
        V0MPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
        RefPercentile = MultSelection->GetMultiplicityPercentile("RefMult08",false);
        
        IndxV0MMult  = GetCentralityClass(V0MPercentile);
        IndxTrksMult = GetCentralityClass(RefPercentile);
//        printf("V0MPer === %f   RefPercentile == %f\n",V0MPercentile,RefPercentile);
    }
    else{
        return;}
    
    if(!fisV0Mestimator) {
      if(IndxTrksMult>=0)
	fCentClass = IndxTrksMult;
      else
	return;
    } else {
      if(IndxV0MMult>=0)
	fCentClass = IndxV0MMult;
      else
	return;
    }
    
    if( fESD->IsPileupFromSPDInMultBins() ){return;}
    fEvents->Fill(2.5,fCentClass);
    fEvents->Fill(2.5,10);
    
    if( fESD->IsIncompleteDAQ() ){return;}
    fEvents->Fill(3.5,fCentClass);
    fEvents->Fill(3.5,10);
    
    if( utils->IsSPDClusterVsTrackletBG(fESD) ){return;}
    fEvents->Fill(4.5,fCentClass);
    fEvents->Fill(4.5,10);
    
    if( !MultSelection->GetThisEventINELgtZERO() ){return;}
    fEvents->Fill(5.5,fCentClass);
    fEvents->Fill(5.5,10);
    
    if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) ){return;}
    fEvents->Fill(6.5,fCentClass);
    fEvents->Fill(6.5,10);
    
    if( !IsGoodZvertexPos(fESD) ){return;}
    fEvents->Fill(7.5,fCentClass);
    fEvents->Fill(7.5,10);
    
    float SOm = -1.0;
    SOm = GetSpherocity(hphiso, hetaso);
    if(fnRefGlobal < 10 && SOm>0)
        cout<<"-------------------------------   "<<fnRefGlobal<<	"      SOm  ==="<<SOm<<endl;
    
    double SOt = -1.0;
    //    if(fAnalysisMC)
    //        SOt = fSpheroUtils->GetEventShapeTrue(fMCStack,hTruthPhiSo,hTruthEtaSo);
    
    // Events with non-measured spherocity
    if( SOm < 0 ){
        fEvents->Fill(8.5,fCentClass);
        fEvents->Fill(8.5,10);
        return;}

    hSOrvsV0M->Fill(V0MPercentile,SOm);
    hSOrvsTrks->Fill(RefPercentile,SOm);
    hRefMultVsRefMultPer->Fill(RefPercentile,fnRefGlobal);
    
    if(fCentClass > 2)
        return;
    
    fEvents->Fill(9.5,fCentClass);
    fEvents->Fill(9.5,10);
    
    //	fTrcksVsTrklets->Fill(fnRefGlobal,nRec);
    
    ProduceArrayTrksESD(fESD,fCentClass,2);
    
    if(fAnalysisMC){
        
        int TruthMult = -1;
        TruthMult = GetMultiplicityParticles(0.8);
        
        if(SOt>0){
            hSOtvsTrks->Fill(fnRefGlobal,SOt);
            hSOtvsV0M->Fill(V0MPercentile,SOt);
            hSOtvsTrkst->Fill(TruthMult,SOt);
            hSOtVsSOm[fCentClass]->Fill(SOt,SOm);
            hSOtVsSOm[10]->Fill(SOt,SOm);
        }
    }
    
    if(fisV0Mestimator){
        if((0.0<SOm)&&(SOm<=fJettyCutOff)){
            ProduceArrayTrksESD(fESD,fCentClass,0);
//            ProduceArrayV0ESD(fESD,fCentClass,0);
            
            if(fAnalysisMC)
                PtRecVsPtTruth(fESD, kTRUE);
            
            fcent->Fill(fCentClass);
            fcent->Fill(10);
        }
        
        if((SOm>=fIsotrCutOff) && (SOm<1.0)){
            ProduceArrayTrksESD(fESD,fCentClass,1);
//            ProduceArrayV0ESD(fESD,fCentClass,1);
            
            if(fAnalysisMC)
                PtRecVsPtTruth(fESD, kFALSE);
            
            fcent->Fill(fCentClass);
            fcent->Fill(10);
        }
        
        if(fAnalysisMC){
            if((0<SOt) && (SOt<0.472))
                ProcessMCTruthESD(fCentClass,0);
            
            if((0.759<SOt) && (SOt<1.0))
                ProcessMCTruthESD(fCentClass,1);
            
        }
    }
    
    //Make PID using Mid-Rapidity estimator
    else{
        if(fCentClass > 2){return;}
        if((0<SOm) && (SOm<fJettyCutOff)){
            ProduceArrayTrksESD(fESD,fCentClass,0);
//            ProduceArrayV0ESD(fESD,fCentClass,0);
            fcent->Fill(fCentClass);
            fcent->Fill(10);}
        if((SOm>fIsotrCutOff) && (SOm<1.0)){
            ProduceArrayTrksESD(fESD,fCentClass,1);
//            ProduceArrayV0ESD(fESD,fCentClass,1);
            fcent->Fill(fCentClass);
            fcent->Fill(10);}
        
        if(fAnalysisMC){
            if((0<SOm) && (SOm<fJettyCutOff))
                ProcessMCTruthESD(fCentClass,0);
            if((fIsotrCutOff<SOm) && (SOm<1.0))
                ProcessMCTruthESD(fCentClass,1);}
    }
    
    PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSpherocity::GetCentralityClass(Float_t percentile)
{
    
    if((percentile<=0))return -1;
    
    Int_t Index = -1;
    if((percentile>70) && (percentile<=100))Index=9;
    else if((percentile>50) && (percentile<=70))Index=8;
    else if((percentile>40) && (percentile<=50))Index=7;
    else if((percentile>30) && (percentile<=40))Index=6;
    else if((percentile>20) && (percentile<=30))Index=5;
    else if((percentile>15) && (percentile<=20))Index=4;
    else if((percentile>10) && (percentile<=15))Index=3;
    else if((percentile>5) && (percentile<=10))Index=2;
    else if((percentile>1) && (percentile<=5))Index=1;
    else if((percentile>0) && (percentile<=1))Index=0;
    else Index=-1;
    
    return Index;
}
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::PtRecVsPtTruth( AliESDEvent *ESDevent, const Bool_t IsJetty )
{
    const Int_t nESDTracks = ESDevent->GetNumberOfTracks();
    for(Int_t iT = 0; iT < nESDTracks; iT++){
        
        AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
        if(!esdTrack){continue;}
        
        UInt_t selectDebug = 0;
        if(fTrackFilter){
            selectDebug = fTrackFilter->IsSelected(esdTrack);
            if (!selectDebug) {continue;}
        }
        
        if(TMath::Abs(esdTrack->Eta())>fEtaCut)
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      ==========  %f\n",esdTrack->Eta());
        
        
        Float_t dcaxy = 0.;
        Float_t dcaz = 0.;
        esdTrack->GetImpactParameters(dcaxy,dcaz);
        if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,esdTrack->Pt())){continue;}
        
        Short_t ncl = esdTrack->GetTPCNcls();
        if(ncl<50)
            printf("--------------------------------------------------         --------------- Ncl      ==========  %d\n",ncl);
        
        const Int_t label = TMath::Abs(esdTrack->GetLabel());
        TParticle* mcTrack = fMCStack->Particle(label);
        
        TParticlePDG* pdgPart = mcTrack->GetPDG();
        Double_t chargeMC = pdgPart->Charge();
        
        if (mcTrack){
            
            if((esdTrack->Charge()==0) || (TMath::Abs(chargeMC) < 0.1))
                continue;
            
            if((TMath::Abs(esdTrack->Eta())>fEtaCut) || (TMath::Abs(mcTrack->Eta())>fEtaCut))
                continue;
            
            if( fMCStack->IsPhysicalPrimary(label) ){
                hPtTruthVsPtRec->Fill(mcTrack->Pt(),esdTrack->Pt());
                
                if(IsJetty)
                    hPtTruthVsPtRecJetty->Fill(mcTrack->Pt(),esdTrack->Pt());
                if(!IsJetty)
                    hPtTruthVsPtRecIsotr->Fill(mcTrack->Pt(),esdTrack->Pt());
            }
        }
    }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSpherocity::GetMultiplicityParticles(Double_t etaCut)
{
    // Fill the special MC histogram with the MC truth info
    
    Int_t trackmult = 0;
    const Int_t nTracksMC = fMCStack->GetNtrack();
    
    for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
        
        TParticle* trackMC = fMCStack->Particle(iTracks);
        if(!trackMC)
            continue;
        
        if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
            continue;
        
        
        TParticlePDG* pdgPart = trackMC ->GetPDG();
        Double_t chargeMC = pdgPart->Charge();
        
        if( TMath::Abs(chargeMC) < 0.1 )
            continue;
        
        if ( TMath::Abs( trackMC->Eta() ) > etaCut )
            continue;
        
        trackmult++;
        
    }//MC track loop
    
    return trackmult;
    
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskSpherocity::GetPidCode(Int_t pdgCode) const
{
    // return our internal code for pions, kaons, and protons
    
    Short_t pidCode = 6;
    
    switch (TMath::Abs(pdgCode)) {
        case 211:
            pidCode = 1; // pion
            break;
        case 321:
            pidCode = 2; // kaon
            break;
        case 2212:
            pidCode = 3; // proton
            break;
        case 11:
            pidCode = 4; // electron
            break;
        case 13:
            pidCode = 5; // muon
            break;
        default:
            pidCode = 6;  // something else?
    };
    
    return pidCode;
}
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::ProcessMCTruthESD(const Int_t Cent, const Int_t Spherocity)
{
    // Fill the special MC histogram with the MC truth info
    
    cout<<"Cent Inside ProcessMCTruth ::: "<<Cent<<endl;
    const Int_t nTracksMC = fMCStack->GetNtrack();
    
    for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
        
        TParticle* trackMC = fMCStack->Particle(iTracks);
        if( !trackMC )
            continue;
        
        if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
            continue;
        
        TParticlePDG* pdgPart = trackMC ->GetPDG();
        Double_t chargeMC = pdgPart->Charge();
        
        if(chargeMC==0)
            continue;
        
        if ( TMath::Abs(trackMC->Eta()) > fEtaCut )
            continue;
        
        Int_t pdgCode = trackMC->GetPdgCode();
        Short_t pidCodeMC = 0;
        pidCodeMC = GetPidCode(pdgCode);
        
        hMcIn[Cent][0][Spherocity]->Fill(trackMC->Pt());
        hMcIn[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());
        
        if( chargeMC < 0 ){
            hMcInNeg[Cent][0][Spherocity]->Fill(trackMC->Pt());
            hMcInNeg[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());
        }
        else{
            hMcInPos[Cent][0][Spherocity]->Fill(trackMC->Pt());
            hMcInPos[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());
        }
        
    }//MC track loop
}
//____________________________________________________________________
TParticle* AliAnalysisTaskSpherocity::FindPrimaryMother(AliStack* stack, Int_t label)
{
    //
    // Finds the first mother among the primary particles of the particle identified by <label>,
    // i.e. the primary that "caused" this particle
    //
    // Taken from AliPWG0Helper class
    //
    
    Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
    if (motherLabel < 0)
        return 0;
    
    return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskSpherocity::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
    //
    // Finds the first mother among the primary particles of the particle identified by <label>,
    // i.e. the primary that "caused" this particle
    //
    // returns its label
    //
    // Taken from AliPWG0Helper class
    //
    const Int_t nPrim  = stack->GetNprimary();
    
    while (label >= nPrim) {
        
        //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));
        
        TParticle* particle = stack->Particle(label);
        if (!particle) {
            
            AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
            return -1;
        }
        
        // find mother
        if (particle->GetMother(0) < 0) {
            
            AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
            return -1;
        }
        
        label = particle->GetMother(0);
    }
    
    return label;
}

//____________________________________________________________________
TParticle* AliAnalysisTaskSpherocity::FindPrimaryMotherV0(AliStack* stack, Int_t label)
{
    //
    // Finds the first mother among the primary particles of the particle identified by <label>,
    // i.e. the primary that "caused" this particle
    //
    // Taken from AliPWG0Helper class
    //
    
    Int_t nSteps = 0;
    
    Int_t motherLabel = FindPrimaryMotherLabelV0(stack, label, nSteps);
    if (motherLabel < 0)
        return 0;
    
    return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskSpherocity::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
{
    //
    // Finds the first mother among the primary particles of the particle identified by <label>,
    // i.e. the primary that "caused" this particle
    //
    // returns its label
    //
    // Taken from AliPWG0Helper class
    //
    nSteps = 0;
    const Int_t nPrim  = stack->GetNprimary();
    
    while (label >= nPrim) {
        
        //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));
        
        nSteps++; // 1 level down
        
        TParticle* particle = stack->Particle(label);
        if (!particle) {
            
            AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
            return -1;
        }
        
        // find mother
        if (particle->GetMother(0) < 0) {
            
            AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
            return -1;
        }
        
        label = particle->GetMother(0);
    }
    
    return label;
}

//__________________________________________________________________
void AliAnalysisTaskSpherocity::ProduceArrayTrksESD( AliESDEvent *ESDevent, const int Cent, const int Spherocity ){
    
    const Int_t nESDTracks = ESDevent->GetNumberOfTracks();
    
    fcentAfterPrimaries->Fill(Cent);
    fcentAfterPrimaries->Fill(10);
    
    float V0MPer  = -1;
    
    AliMultSelection *MultSelection = (AliMultSelection*)ESDevent -> FindListObject("MultSelection");
    if(MultSelection-> IsEventSelected())
        V0MPer = MultSelection->GetMultiplicityPercentile("V0M",false);
    
    int multTPC = 0;
    for(Int_t iT = 0; iT < nESDTracks; iT++) {
        
        AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
        
        if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
            continue;
        
        UInt_t selectDebug = 0;
        if (fTrackFilterTPC) {
            selectDebug = fTrackFilterTPC->IsSelected(esdTrack);
            if (!selectDebug) {
                continue;
            }
        }
        
        multTPC++;
        
    }
    
    for(Int_t iT = 0; iT < nESDTracks; iT++) {
        
        AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
        if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
            continue;
        
        UInt_t selectDebug = 0;
        if(fTrackFilterGolden){
            selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
            if (!selectDebug) {
                continue;
            }
        }
        
        double eta      = esdTrack->Eta();
        double phi      = esdTrack->Phi();
        double momentum = esdTrack->P();
        double pt       = esdTrack->Pt();
        float  dedx     = esdTrack->GetTPCsignal();
        float  dedxUnc  = esdTrack->GetTPCsignal();        
        float dcaxy = 0.0;
        float dcaz = 0.0;
        esdTrack->GetImpactParameters(dcaxy,dcaz);
        
        int nh = -1;
        if(TMath::Abs(eta)<0.2)
            nh = 0;
        else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
            nh = 1;
        else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
            nh = 2;
        else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
            nh = 3;
        
        if(nh<0)
            continue;
        
        short ncl = esdTrack->GetTPCsignalN();
        if( ncl < fNcl )
            continue;
        
        if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
            continue;
        
        if( TOFPID(esdTrack) ){
            
            double trkLength = esdTrack->GetIntegratedLength();
            double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);
            
            if(esdTrack->Charge() < 0){
                hBetavsPneg[Cent][nh][Spherocity]->Fill(esdTrack->P(),beta);
                hPtneg_TOF_Eta[Cent][nh][Spherocity]->Fill(pt);
                hPneg_TOF_Eta[Cent][nh][Spherocity]->Fill(esdTrack->P());
                hPtneg_TOF[Cent][Spherocity]->Fill(pt);
            }else{
                hBetavsPpos[Cent][nh][Spherocity]->Fill(esdTrack->P(),beta);
                hPtpos_TOF_Eta[Cent][nh][Spherocity]->Fill(pt);
                hPpos_TOF_Eta[Cent][nh][Spherocity]->Fill(esdTrack->P());
                hPtpos_TOF[Cent][Spherocity]->Fill(pt);
            }
            
        }
        
        if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), MAGF, fcutLow, fcutHigh))
            continue;
        
        if(fdEdxCalibrated){
            if(eta < 0){
                if(fLHC16l == 1)
                    dedx *= 50/EtaCalibrationNeg(0,eta);
                else
                    dedx *= 50/EtaCalibrationNeg(1,eta);
            }
            else{
                if(fLHC16l == 1)
                    dedx *= 50/EtaCalibrationPos(0,eta);
                else
                    dedx *= 50/EtaCalibrationPos(1,eta);
            }
        }
        
        Short_t pidCode     = 0;
        if(fAnalysisMC) {
            
            const Int_t label = TMath::Abs(esdTrack->GetLabel());
            TParticle* mcTrack = 0;
            mcTrack = fMCStack->Particle(label);
            
            if (mcTrack){
                
                if( esdTrack->Charge()==0 )
                    continue;
                
                Int_t pdgCode = mcTrack->GetPdgCode();
                pidCode = GetPidCode(pdgCode);
                
                if( fMCStack->IsPhysicalPrimary(label) ){
                    if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
                        hMcOut[Cent][0][Spherocity]->Fill(esdTrack->Pt());
                        hMcOut[Cent][pidCode][Spherocity]->Fill(esdTrack->Pt());
                    }
                    
                    if( esdTrack->Charge() < 0.0 ){
                        if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {
                            hMcOutNeg[Cent][0][Spherocity]->Fill(esdTrack->Pt());
                            hMcOutNeg[Cent][pidCode][Spherocity]->Fill(esdTrack->Pt());
                        }
                    }
                    else{
                        if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
                            hMcOutPos[Cent][0][Spherocity]->Fill(esdTrack->Pt());
                            hMcOutPos[Cent][pidCode][Spherocity]->Fill(esdTrack->Pt());
                        }
                    }
                }	// Primary particles MC
            }	//mcTrack
        }	//fAnalysis MC
        
        
        //TOF
        Bool_t IsTOFout=kFALSE;
        if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
            IsTOFout=kTRUE;
        Float_t lengthtrack=esdTrack->GetIntegratedLength();
        Float_t timeTOF=esdTrack->GetTOFsignal();
        Double_t inttime[5]={0,0,0,0,0};
        esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
        Float_t beta = -99;
        if ( !IsTOFout ){
            if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
                beta = inttime[0] / timeTOF;
        }
        
        if(!fdEdxCalibrated){
            if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
                if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                    hMIPVsEta[Cent][Spherocity]->Fill(eta,dedxUnc);
                    pMIPVsEta[Cent][Spherocity]->Fill(eta,dedxUnc);
                    //hMIPVsEta[10][Spherocity]->Fill(eta,dedxUnc);
                    //pMIPVsEta[10][Spherocity]->Fill(eta,dedxUnc);
                }
                if( dedxUnc > 70 && dedxUnc < 90 ){
                    if(TMath::Abs(beta-1)<0.1){
                        hPlateauVsEta[Cent][Spherocity]->Fill(eta,dedxUnc);
                        pPlateauVsEta[Cent][Spherocity]->Fill(eta,dedxUnc);
                        //hPlateauVsEta[10][Spherocity]->Fill(eta,dedxUnc);
                        //pPlateauVsEta[10][Spherocity]->Fill(eta,dedxUnc);
                    }
                }
            }
        }
        else{
            if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
                if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                    hMIPVsEta[Cent][Spherocity]->Fill(eta,dedx);
                    pMIPVsEta[Cent][Spherocity]->Fill(eta,dedx);
                    //hMIPVsEta[10][Spherocity]->Fill(eta,dedx);
                    //pMIPVsEta[10][Spherocity]->Fill(eta,dedx);
                }
                if( dedxUnc > 70 && dedxUnc < 90 ){
                    if(TMath::Abs(beta-1)<0.1){
                        hPlateauVsEta[Cent][Spherocity]->Fill(eta,dedx);
                        pPlateauVsEta[Cent][Spherocity]->Fill(eta,dedx);
                        //hPlateauVsEta[10][Spherocity]->Fill(eta,dedx);
                        //pPlateauVsEta[10][Spherocity]->Fill(eta,dedx);
                    }
                }
            }
        }
        
        if(beta>1){
            histPiTof[Cent][nh][Spherocity]->Fill(momentum, dedx);
            //histPiTof[10][nh]->Fill(momentum, dedx);
        }
        
        if( momentum <= 0.6 && momentum >= 0.4  ){
            if( dedx < fDeDxMIPMax && dedx > fDeDxMIPMin ){
                hMIPVsPhi[Cent][nh][Spherocity]->Fill(phi,dedx);
                pMIPVsPhi[Cent][nh][Spherocity]->Fill(phi,dedx);
                hMIPVsNch[nh]->Fill(multTPC,dedx);
                pMIPVsNch[nh]->Fill(multTPC,dedx);
                hMIPVsV0M[nh]->Fill(V0MPer,dedx);  
                pMIPVsV0M[nh]->Fill(V0MPer,dedx);
            }
            if( dedx > 70 && dedx < 90 ){
                if(TMath::Abs(beta-1)<0.1){
                    hPlateauVsPhi[Cent][nh][Spherocity]->Fill(phi,dedx);
                    pPlateauVsPhi[Cent][nh][Spherocity]->Fill(phi,dedx);
                }
            }
        }
        
        hPtAll[Cent][Spherocity]->Fill(pt);
        hPtVsP[Cent][nh][Spherocity]->Fill(momentum,pt);
        hDeDxVsP[Cent][nh][Spherocity]->Fill(momentum,dedx);
        
        if(esdTrack->Charge() < 0.0){
            hnSigmaPiNeg[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
            hnSigmaKNeg[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
            hnSigmaPNeg[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
            hPtneg_TPC_Eta[Cent][nh][Spherocity]->Fill(pt);
            hPtneg_TPC[Cent][Spherocity]->Fill(pt);
            
        }else{
            hnSigmaPiPos[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
            hnSigmaKPos[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
            hnSigmaPPos[Cent][nh][Spherocity]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
            hPtpos_TPC_Eta[Cent][nh][Spherocity]->Fill(pt);
            hPtpos_TPC[Cent][Spherocity]->Fill(pt);
        }
        
    }//end of track loop
    
    
}

//----------------------------------------------------------------------------------

void AliAnalysisTaskSpherocity::ProduceArrayV0ESD( AliESDEvent *ESDevent, const Int_t Cent, const Int_t Spherocity ){
    
    Int_t nv0s = ESDevent->GetNumberOfV0s();
    
    fcentAfterV0s->Fill(Cent);
    fcentAfterV0s->Fill(10);
    
    const AliESDVertex *myBestPrimaryVertex = ESDevent->GetPrimaryVertex();
    
    if ( !myBestPrimaryVertex )
        return;
    if ( !(myBestPrimaryVertex->GetStatus()) )
        return;
    
    Double_t  lPrimaryVtxPosition[3];
    myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
    Double_t  lPrimaryVtxCov[6];
    myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
    Double_t  lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
    
    AliAODVertex* myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    
    //
    // LOOP OVER V0s, K0s, L, AL
    //
    
    for (Int_t iV0=0; iV0<nv0s; iV0++) {
        
        AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
        if ( !esdV0 ) continue;
        
        //check onfly status
        //		if( !esdV0->GetOnFlyStatus() )
        //			continue;
        
        if( esdV0->GetOnFlyStatus()!=0 )
            continue;
        
        // AliESDTrack (V0 Daughters)
        UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());
        
        AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);
        
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        
        // Remove like-sign
        if (pTrack->GetSign() == nTrack->GetSign())
            continue;
        
        // Eta cut on decay products
        if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
            continue;
        
        UInt_t selectDebug_p = 0;
        if ( fTrackFilterGolden ) {
            selectDebug_p = fTrackFilterGolden->IsSelected(pTrack);
            if (!selectDebug_p) {
                continue;
            }
        }
        
        UInt_t selectDebug_n = 0;
        if ( fTrackFilterGolden ) {
            selectDebug_n = fTrackFilterGolden->IsSelected(nTrack);
            if (!selectDebug_n) {
                continue;
            }
        }
        
        
        // Check if switch does anything!
        Bool_t isSwitched = kFALSE;
        if (pTrack->GetSign() < 0) { // switch
            isSwitched = kTRUE;
            AliESDtrack* helpTrack = nTrack;
            nTrack = pTrack;
            pTrack = helpTrack;
        }
        
        AliKFVertex primaryVtxKF( *myPrimaryVertex );
        AliKFParticle::SetField(ESDevent->GetMagneticField());
        
        // Also implement switch here!!!!!!
        AliKFParticle* negEKF  = 0; // e-
        AliKFParticle* posEKF  = 0; // e+
        AliKFParticle* negPiKF = 0; // pi -
        AliKFParticle* posPiKF = 0; // pi +
        AliKFParticle* posPKF  = 0; // p
        AliKFParticle* negAPKF = 0; // p-bar
        
        if(!isSwitched) {
            negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
            posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
            negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
            posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
            posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
            negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
        } else { // switch + and -
            negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
            posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
            negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
            posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
            posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
            negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
        }
        
        AliKFParticle v0GKF;  // Gamma e.g. from pi0
        v0GKF+=(*negEKF);
        v0GKF+=(*posEKF);
        v0GKF.SetProductionVertex(primaryVtxKF);
        
        AliKFParticle v0K0sKF; // K0 short
        v0K0sKF+=(*negPiKF);
        v0K0sKF+=(*posPiKF);
        v0K0sKF.SetProductionVertex(primaryVtxKF);
        
        AliKFParticle v0LambdaKF; // Lambda
        v0LambdaKF+=(*negPiKF);
        v0LambdaKF+=(*posPKF);
        v0LambdaKF.SetProductionVertex(primaryVtxKF);
        
        AliKFParticle v0AntiLambdaKF; // Lambda-bar
        v0AntiLambdaKF+=(*posPiKF);
        v0AntiLambdaKF+=(*negAPKF);
        v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);
        
        Double_t dmassG     = TMath::Abs(v0GKF.GetMass());
        Double_t dmassK     = TMath::Abs(v0K0sKF.GetMass()-0.498);
        Double_t dmassL     = TMath::Abs(v0LambdaKF.GetMass()-1.116);
        Double_t dmassAL    = TMath::Abs(v0AntiLambdaKF.GetMass()-1.116);
        
        if( dmassG  > 0.1 &&
           dmassK  > 0.1 &&
           dmassL  > 0.1 &&
           dmassAL > 0.1
           )
            continue;
        
        
        for( Int_t case_v0 = 0; case_v0 < 2; ++case_v0 ){
            switch(case_v0){
                case 0:{
                    
                    Bool_t fillPos = kFALSE;
                    Bool_t fillNeg = kFALSE;
                    
                    if(dmassG < 0.1)
                        continue;
                    
                    if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01){
                        continue;
                    }
                    
                    if(dmassL<0.01){
                        fillPos = kTRUE;
                    }
                    
                    if(dmassAL<0.01) {
                        if(fillPos)
                            continue;
                        fillNeg = kTRUE;
                    }
                    
                    if(dmassK<0.01) {
                        if(fillPos||fillNeg)
                            continue;
                        fillPos = kTRUE;
                        fillNeg = kTRUE;
                    }
                    
                    
                    for(Int_t j = 0; j < 2; j++) {
                        
                        AliESDtrack* track = 0;
                        
                        if(j==0) {
                            
                            if(fillNeg)
                                track = nTrack;
                            else
                                continue;
                        } else {
                            
                            if(fillPos)
                                track = pTrack;
                            else
                                continue;
                        }
                        
                        if(track->GetTPCsignalN()<fNcl)continue;
                        Double_t phi     = track->Phi();
                        
                        if(!PhiCut(track->Pt(), phi, track->Charge(), MAGF, fcutLow, fcutHigh))
                            continue;
                        
                        Double_t eta      = track->Eta();
                        Double_t momentum = track->P();
                        Double_t dedx     = track->GetTPCsignal();
                        Double_t dedxUnc  = track->GetTPCsignal();
                        
                        
                        if(fdEdxCalibrated){
                            if(eta < 0){
                                if(fLHC16l == 1)
                                    dedx   *= 50.0/EtaCalibrationNeg(0,eta);
                                else
                                    dedx   *= 50.0/EtaCalibrationNeg(1,eta);
                            }
                            else{
                                if(fLHC16l == 1)
                                    dedx   *= 50.0/EtaCalibrationPos(0,eta);
                                else
                                    dedx   *= 50.0/EtaCalibrationPos(1,eta);
                            }
                        }
                        
                        
                        if(fillPos&&fillNeg){
                            if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                                if(momentum<0.6&&momentum>0.4){
                                    hMIPVsEtaV0s[Cent][Spherocity]->Fill(eta,dedx);
                                    pMIPVsEtaV0s[Cent][Spherocity]->Fill(eta,dedx);
                                    //hMIPVsEtaV0s[10][Spherocity]->Fill(eta,dedx);
                                    //pMIPVsEtaV0s[10][Spherocity]->Fill(eta,dedx);
                                }
                            }
                        }
                        
                        Int_t nh = -1;
                        
                        if(TMath::Abs(eta)<0.2)
                            nh = 0;
                        else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
                            nh = 1;
                        else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
                            nh = 2;
                        else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
                            nh = 3;
                        
                        if(nh<0)
                            continue;
                        
                        if(fillPos&&fillNeg){
                            histPiV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
                            //histPiV0[10][nh]->Fill(momentum, dedx);
                        }
                        else{
                            histPV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
                            //histPV0[10][nh]->Fill(momentum, dedx);
                        }
                    }//end loop over two tracks
                };
                    break;
                    
                case 1:{//gammas
                    
                    Bool_t fillPos = kFALSE;
                    Bool_t fillNeg = kFALSE;
                    
                    
                    if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
                        if( dmassG<0.01 && dmassG>0.0001 ) {
                            if(nTrack->Eta() > 0){
                                if(fLHC16l == 1){
                                    if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationPosEl(0,nTrack->Eta())) < 5)
                                        fillPos = kTRUE;
                                }
                                else{
                                    if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationPosEl(1,nTrack->Eta())) < 5)
                                        fillPos = kTRUE;
                                    cout << "-------------------LHC16k"  << endl;
                                }
                                
                            }
                            if(nTrack->Eta() < 0){
                                if(fLHC16l == 1){
                                    if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationNegEl(0,nTrack->Eta())) < 5)
                                        fillPos = kTRUE;
                                }
                                else{
                                    if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationNegEl(1,nTrack->Eta())) < 5)
                                        fillPos = kTRUE;
                                    cout << "-------------------LHC16k"  << endl;
                                }
                            }
                            
                        } else {
                            continue;
                        }
                    }
                    
                    
                    
                    if(fillPos == kTRUE && fillNeg == kTRUE)
                        continue;
                    
                    AliESDtrack* track = 0;
                    if(fillNeg)
                        track = nTrack;
                    else if(fillPos)
                        track = pTrack;
                    else
                        continue;
                    
                    Double_t dedx     = track->GetTPCsignal();
                    Double_t eta      = track->Eta();
                    Double_t phi      = track->Phi();
                    Double_t momentum = track->P();
                    
                    if(fdEdxCalibrated){
                        if(eta < 0){
                            if(fLHC16l == 1)
                                dedx *= 50/EtaCalibrationNeg(0,eta);
                            else
                                dedx *= 50/EtaCalibrationNeg(1,eta);
                        }
                        else{
                            if(fLHC16l == 1)
                                dedx *= 50/EtaCalibrationPos(0,eta);
                            else
                                dedx *= 50/EtaCalibrationPos(1,eta);
                        }
                    }
                    
                    
                    //					       if(track->GetTPCsignalN()<=70)continue;
                    
                    if(!PhiCut(track->Pt(), phi, track->Charge(), MAGF, fcutLow, fcutHigh))
                        continue;
                    
                    Int_t nh = -1;
                    
                    if(TMath::Abs(eta)<0.2)
                        nh = 0;
                    else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
                        nh = 1;
                    else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
                        nh = 2;
                    else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
                        nh = 3;
                    
                    if(nh<0)
                        continue;
                    
                    
                    cout<<"Cent :: "<<Cent<<"Eta  :: "<<eta<<"dedx :: "<<dedx<<endl;
                    histEV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
                    //histEV0[10][nh]->Fill(momentum, dedx);
                };
                    break;
                    
                    
            }//end switch
            
        }//end loop over case V0
        
        
        // clean up loop over v0
        
        delete negPiKF;
        delete posPiKF;
        delete posPKF;
        delete negAPKF;
        
        
        
    }
    
    
    delete myPrimaryVertex;
    
    
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::selectVertex2015pp(AliESDEvent *esd,
                                                     Bool_t checkSPDres, //enable check on vtx resolution
                                                     Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
                                                     Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{
    
    if (!esd) return kFALSE;
    
    const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
    const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
    Bool_t hasSPD = spdVertex->GetStatus();
    Bool_t hasTrk = trkVertex->GetStatus();
    
    //Note that AliVertex::GetStatus checks that N_contributors is > 0
    //reject events if both are explicitly requested and none is available
    if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
    
    //reject events if none between the SPD or track verteces are available
    //if no trk vertex, try to fall back to SPD vertex;
    if (!hasTrk) {
        if (!hasSPD) return kFALSE;
        //on demand check the spd vertex resolution and reject if not satisfied
        if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
    } else {
        if (hasSPD) {
            //if enabled check the spd vertex resolution and reject if not satisfied
            //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
            if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
            if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
        }
    }
    
    //Cut on the vertex z position
    //const AliESDVertex * vertex = esd->GetPrimaryVertex();
    //if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{
    
    if( !spdVertex ) return kFALSE;
    if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::IsGoodZvertexPos(AliESDEvent *esd)
{
    
    if( !esd ) return kFALSE;
    //Cut on the vertex z position
    const AliESDVertex * vertex = esd->GetPrimaryVertex();
    if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSpherocity::SetTrackCutsSpherocity(AliAnalysisFilter* fTrackFilter){
    
    AliESDtrackCuts* esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    fTrackFilter->AddCuts(esdTrackCuts);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)
{
    if(pt < 2.0)
        return kTRUE;
    
    //Double_t phi = track->Phi();
    if(mag < 0)    // for negatve polarity field
        phi = TMath::TwoPi() - phi;
    if(q < 0) // for negatve charge
        phi = TMath::TwoPi()-phi;
    
    phi += TMath::Pi()/18.0; // to center gap in the middle
    phi = fmod(phi, TMath::Pi()/9.0);
    
    if(phi<phiCutHigh->Eval(pt)
       && phi>phiCutLow->Eval(pt))
        return kFALSE; // reject track
    
    hPhi[fCentClass]->Fill(pt, phi);
    
    return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskSpherocity::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){
    
    Double_t maxDCAxy = 10;
    maxDCAxy = fMaxDCAxy->Eval(ptI);
    return maxDCAxy;
    
}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::GetSpherocity( TH1D * hphi, TH1D *heta )
{
    vector<Float_t> pt;
    vector<Float_t> eta;
    vector<Float_t> phi;
    
    if(fESD)
        fNrec = ReadESDEvent( pt, eta, phi, hphi, heta );
    
    if( fNrec < fMinMult )
        return -0.5;
    
    float spherocity = AnalyseGetSpherocity( pt, eta, phi );
    
    return spherocity;
    
}
//_____________________________________________________________________
int AliAnalysisTaskSpherocity::ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta ){
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    
    int nTracks = fESD->GetNumberOfTracks();
    int nRec    = 0;
    
    for(int iT = 0; iT < nTracks; iT++) {
        AliESDtrack* Track = 0;
        Track = fESD->GetTrack(iT);
        if(!Track)
            continue;
        
        float eta  = Track->Eta();
        float pt   = Track->Pt();
        float phi  = Track->Phi();
        
        if( !(TMath::Abs(eta) < fEtaCut) )
            continue;
        
        //cuts in pt
        if( pt > 1E8 || pt < 0.15 )
            continue;
        
        //quality cuts
        if(!fTrackFilter->IsSelected(Track))
            continue;
        
        if(hphi)
            hphi->Fill(phi);
        if(heta)
            heta->Fill(eta);
        
        ptArray.push_back(pt);
        etaArray.push_back(eta);
        phiArray.push_back(phi);
        
        nRec++;
        
    } //close first loop on nTracks
    
    return nRec;
    
}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::AnalyseGetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi ){
    
    
    float spherocity = -10.0;
    float pFull = 0;
    float Spherocity = 2;
    
    //computing total pt
    float sumapt = 0;
    for(int i1 = 0; i1 < fNrec; ++i1){
    //    sumapt += pt[i1];       
        sumapt++;       
    }
    
    //Getting thrust
    for(int i = 0; i < 360/(fSizeStep); ++i){
        float numerador = 0;
        float phiparam  = 0;
        float nx = 0;
        float ny = 0;
        phiparam=( (TMath::Pi()) * i * fSizeStep ) / 180; // parametrization of the angle
        nx = TMath::Cos(phiparam);            // x component of an unitary vector n
        ny = TMath::Sin(phiparam);            // y component of an unitary vector n
        for(int i1 = 0; i1 < fNrec; ++i1){
            
//            float pxA = pt[i1] * TMath::Cos( phi[i1] );
//            float pyA = pt[i1] * TMath::Sin( phi[i1] );
            float pxA = 1.0 * TMath::Cos( phi[i1] );
            float pyA = 1.0 * TMath::Sin( phi[i1] );
            
            numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
        }
        pFull=TMath::Power( (numerador / sumapt),2 );
        if(pFull < Spherocity)//maximization of pFull
        {
            Spherocity = pFull;
        }
    }
    
    spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
    
    
    return spherocity;
    
}
//________________________________________________________________________
bool AliAnalysisTaskSpherocity::TOFPID(AliESDtrack * track)
{
    UInt_t status;
    status=track->GetStatus();
    
    if (!(status & AliESDtrack::kTOFout) || !(status & AliESDtrack::kTIME))
        return kFALSE;
    
    if (track->GetIntegratedLength() < 350.)
        return kFALSE;
    
    if (TMath::Abs(track->GetTOFsignalDx()) > 10.0 || TMath::Abs(track->GetTOFsignalDz()) > 10.0)
        return kFALSE;
    
    return kTRUE;
}
//________________________________________________________________________
Double_t AliAnalysisTaskSpherocity::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){
    
    
    for(Int_t i=0; i<8; ++i)
        fEtaCalibrationNeg->SetParameter(i,0);
    
    fEtaCalibrationNeg->SetParameter(0,aNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(1,bNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(2,cNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(3,dNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(4,eNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(5,fNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(6,gNeg[Cent]);
    fEtaCalibrationNeg->SetParameter(7,hNeg[Cent]);
    
    return fEtaCalibrationNeg->Eval(eta);
    
    
}
//________________________________________________________________________
Double_t AliAnalysisTaskSpherocity::EtaCalibrationPos( const Int_t Cent, const Double_t eta){
    
    
    for(Int_t i=0; i<8; ++i)
        fEtaCalibration->SetParameter(i,0);
    
    fEtaCalibration->SetParameter(0,aPos[Cent]);
    fEtaCalibration->SetParameter(1,bPos[Cent]);
    fEtaCalibration->SetParameter(2,cPos[Cent]);
    fEtaCalibration->SetParameter(3,dPos[Cent]);
    fEtaCalibration->SetParameter(4,ePos[Cent]);
    fEtaCalibration->SetParameter(5,fPos[Cent]);
    fEtaCalibration->SetParameter(6,gPos[Cent]);
    fEtaCalibration->SetParameter(7,hPos[Cent]);
    
    return fEtaCalibration->Eval(eta);
    
}
//________________________________________________________________________
Double_t AliAnalysisTaskSpherocity::EtaCalibrationNegEl(const Int_t Cent, const Double_t eta){
    
    
    for(Int_t i=0; i<5; ++i)
        felededxfitNeg->SetParameter(i,0);
    
    
    felededxfitNeg->SetParameter(0,aNegEl[Cent]);
    felededxfitNeg->SetParameter(1,bNegEl[Cent]);
    felededxfitNeg->SetParameter(2,cNegEl[Cent]);
    felededxfitNeg->SetParameter(3,dNegEl[Cent]);
    felededxfitNeg->SetParameter(4,eNegEl[Cent]);
    
    
    return felededxfitNeg->Eval(eta);
    
}
//________________________________________________________________________
Double_t AliAnalysisTaskSpherocity::EtaCalibrationPosEl(const Int_t Cent, const Double_t eta){
    
    
    for(Int_t i=0; i<5; ++i)
        felededxfitPos->SetParameter(i,0);
    
    felededxfitPos->SetParameter(0,aPosEl[Cent]);
    felededxfitPos->SetParameter(1,bPosEl[Cent]);
    felededxfitPos->SetParameter(2,cPosEl[Cent]);
    felededxfitPos->SetParameter(3,dPosEl[Cent]);
    felededxfitPos->SetParameter(4,ePosEl[Cent]);
    
    return felededxfitPos->Eval(eta);
    
}

