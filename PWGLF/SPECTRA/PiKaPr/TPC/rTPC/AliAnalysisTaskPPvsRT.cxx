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

#include "AliAnalysisTaskPPvsRT.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH3.h>
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

#include <iostream>
using namespace std;


//
// Responsible:
// Antonio Ortiz (Lund)
// Peter Christiansen (Lund)
//


static float Magf                = 1;
static const int nHists          = 4;
static const int nRt             = 5;
static const double CentMin[nRt] = {0.0,0.5,1.5,2.5,0.0};
static const double CentMax[nRt] = {0.5,1.5,2.5,5.0,5.0};

ClassImp(AliAnalysisTaskPPvsRT)
AliAnalysisTaskPPvsRT::AliAnalysisTaskPPvsRT():
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
fAnalysisType("ESD"),
fAnalysisMC(kFALSE),
fAnalysisPbPb(kFALSE),
fRandom(0x0),
fVtxCut(0x0),
fLeadingCut(5.0),
fNcl(70),
fEtaCut(0.9),
cent(3),
fMinCent(0.0),
fMaxCent(100.0),
fDeDxMIPMin(40),
fDeDxMIPMax(60),
fdEdxHigh(200),
fdEdxLow(40),
fPeriod("l"),
fMeanChT(7.11),
fMcProcessType(-999),
fTriggeredEventMB(-999),
fVtxStatus(-999),
fZvtx(-999),
fZvtxMC(-999),
fRun(-999),
fEventId(-999),
fListOfObjects(0),
fEvents(0x0),
fVtxMC(0x0),
fdEdxCalibrated(0x0),
fMakePid(0x0),
fLowPt(0x0),
fMultN(0x0),
fPtN(0x0),
fDphiN(0x0),
fMultA(0x0),
fPtA(0x0),
fDphiA(0x0),
fMultT(0x0),
fPtT(0x0),
fMultNMC(0x0),
fMultAMC(0x0),
fMultTMC(0x0),
fDphiT(0x0),
fDphi(0x0),
fRT(0x0),
fRTMC(0x0),
fPtLVsRt(0x0),
fPtLVsRtMC(0x0),
fEtaCalibration(0x0),
fEtaCalibrationEl(0x0),
fcutDCAxy(0x0),
fcutLow(0x0),
fcutHigh(0x0)

{
    
    for(Int_t i = 0; i<nRt; ++i){
        
        for(Int_t r= 0; r < 3; r++){
            hMIPVsEta[r][i]=0;
            pMIPVsEta[r][i]=0;
            hPtAll[r][i]=0;
        }
        
        hMIPVsEtaV0s[i]=0;
        pMIPVsEtaV0s[i]=0;
        hPlateauVsEta[i]=0;
        pPlateauVsEta[i]=0;
        hPhi[i]=0;
        hDCAxyVsPtPi[i]=0;
        hDCAxyVsPtp[i]=0;
        hDCAxyVsPtPiC[i]=0;
        hDCAxyVsPtpC[i]=0;
        
        for(Int_t j=0; j<nHists; ++j){
            
            for(Int_t r= 0; r < 3; r++)
                hDeDxVsP[r][i][j]=0;//TH2D, DeDx vs P
            
            hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
            pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
            hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
            pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
            hnSigmaPi[i][j]=0;//TH2D, nSigmas vs Pt Pions
            hnSigmak[i][j]=0;//TH2D, nSigmas vs Pt Kaons
            hnSigmap[i][j]=0;//TH2D, nSigmas vs Pt Protons
            histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
            histpPiV0[i][j]=0;//TH1D, pi id by V0s
            histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
            histpPV0[i][j]=0;// TH1D, p id by V0s
            histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
            histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
            histEV0[i][j]=0;
        }
    }
    
    for(Int_t j=0; j<nHists; ++j)
        hPtVsP[j]=0;
    
    //default constructor
    for(Int_t cent=0;cent<nRt;++cent){
        for(Int_t pid=0;pid<7;++pid){
            hMcIn[cent][pid]     = 0;
            hMcOut[cent][pid]    = 0;
            hDCApTPrim[cent][pid]  = 0;
            hDCApTWDec[cent][pid]  = 0;
            hDCApTMate[cent][pid]  = 0;
            hDCApTPrim2[cent][pid] = 0;
            hDCApTWDec2[cent][pid] = 0;
            hDCApTMate2[cent][pid] = 0;
        }
    }
    
}


AliAnalysisTaskPPvsRT::AliAnalysisTaskPPvsRT(const char *name):
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
fAnalysisType("ESD"),
fAnalysisMC(kFALSE),
fAnalysisPbPb(kFALSE),
fRandom(0x0),
fVtxCut(0x0),
fLeadingCut(5.0),
fNcl(70),
fEtaCut(0.9),
cent(3),
fMinCent(0.0),
fMaxCent(100.0),
fDeDxMIPMin(40),
fDeDxMIPMax(60),
fdEdxHigh(200),
fdEdxLow(40),
fPeriod("l"),
fMeanChT(7.11),
fMcProcessType(-999),
fTriggeredEventMB(-999),
fVtxStatus(-999),
fZvtx(-999),
fZvtxMC(-999),
fRun(-999),
fEventId(-999),
fListOfObjects(0),
fEvents(0x0),
fVtxMC(0x0),
fdEdxCalibrated(0x0),
fMakePid(0x0),
fLowPt(0x0),
fMultN(0x0),
fPtN(0x0),
fDphiN(0x0),
fMultA(0x0),
fPtA(0x0),
fDphiA(0x0),
fMultT(0x0),
fPtT(0x0),
fMultNMC(0x0),
fMultAMC(0x0),
fMultTMC(0x0),
fDphiT(0x0),
fDphi(0x0),
fRT(0x0),
fRTMC(0x0),
fPtLVsRt(0x0),
fPtLVsRtMC(0x0),
fEtaCalibration(0x0),
fEtaCalibrationEl(0x0),
fcutDCAxy(0x0),
fcutLow(0x0),
fcutHigh(0x0)

{
    
    for(Int_t i = 0; i<nRt; ++i){
        
        for(Int_t r= 0; r < 3; r++){
            hMIPVsEta[r][i]=0;
            pMIPVsEta[r][i]=0;
            hPtAll[r][i]=0;
        }
        
        hMIPVsEtaV0s[i]=0;
        pMIPVsEtaV0s[i]=0;
        hPlateauVsEta[i]=0;
        pPlateauVsEta[i]=0;
        hPhi[i]=0;
        hDCAxyVsPtPi[i]=0;
        hDCAxyVsPtp[i]=0;
        hDCAxyVsPtPiC[i]=0;
        hDCAxyVsPtpC[i]=0;
        
        for(Int_t j=0; j<nHists; ++j){
            
            for(Int_t r= 0; r < 3; r++)
                hDeDxVsP[r][i][j]=0;//TH2D, DeDx vs P
            
            hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
            pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
            hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
            pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
            hnSigmaPi[i][j]=0;//TH2D, nSigmas vs Pt Pions
            hnSigmak[i][j]=0;//TH2D, nSigmas vs Pt Kaons
            hnSigmap[i][j]=0;//TH2D, nSigmas vs Pt Protons
            histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
            histpPiV0[i][j]=0;//TH1D, pi id by V0s
            histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
            histpPV0[i][j]=0;// TH1D, p id by V0s
            histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
            histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
            histEV0[i][j]=0;
            
        }
        
    }

    for(Int_t j=0; j<nHists; ++j)
        hPtVsP[j]=0;
    
    for(Int_t cent=0; cent<nRt; ++cent){
        for(Int_t pid=0; pid<7; ++pid){
            hMcIn[cent][pid]=0;
            hMcOut[cent][pid]=0;
            hDCApTPrim[cent][pid]=0;
            hDCApTWDec[cent][pid]=0;
            hDCApTMate[cent][pid]=0;
            hDCApTPrim2[cent][pid]=0;
            hDCApTWDec2[cent][pid]=0;
            hDCApTMate2[cent][pid]=0;
        }
    }
    
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());//esto es nuevo
}

AliAnalysisTaskPPvsRT::~AliAnalysisTaskPPvsRT() {
    //
    // Destructor
    //
    
}
//______________________________________________________________________________
void AliAnalysisTaskPPvsRT::UserCreateOutputObjects()
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
        fTrackFilter = new AliAnalysisFilter("trackFilter");
        SetTrackCuts(fTrackFilter);
    }
    
    
    //OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner(kTRUE);
    
    //
    // Histograms
    //
    
    
    fEvents = new TH1F( "fEvents", "; Evt. Sel.",12,0,12);
    fEvents->GetXaxis()->SetBinLabel(1, "Processed");
    fEvents->GetXaxis()->SetBinLabel(2, "PhysSel+Trigger");
    fEvents->GetXaxis()->SetBinLabel(3, "INEL>0");
    fEvents->GetXaxis()->SetBinLabel(4, "BG");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(5, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(6, "Incom DAQ");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(7, "Res&Proximity");//NotinVertexcut");
    fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
    fListOfObjects->Add(fEvents);
    
    const Int_t nPtBins = 63;
    Double_t ptBins[nPtBins+1] = {
        0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
        0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
        0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
        1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
        3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
        9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00,
        20.00,22.00,24.00,26.00,30.00};
    
    const Int_t nPtBinsV0s = 25;
    Double_t ptBinsV0s[nPtBinsV0s+1] = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0,
        9.0 , 12.0, 15.0, 20.0 };
    
    const Int_t ndcaBins = 100;
    Double_t dcaBins[ndcaBins+1] = {
        -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1,
        -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1,
        -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1,
        -1.0, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6,
        -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15,
        -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
        2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4,
        3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
    
    fMultN = new TH1F("hMultN","hMult_Near; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultN->Sumw2();
    fListOfObjects->Add(fMultN);
    
    fMultA = new TH1F("hMultA","hMult_Away; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultA->Sumw2();
    fListOfObjects->Add(fMultA);
    
    fMultT = new TH1F("hMultT","hMult_Transverse; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultT->Sumw2();
    fListOfObjects->Add(fMultT);
    
    fPtN = new TH1F("hPtN","hPt_Near; #it{p}_{T}; Entries",nPtBins,ptBins);
    fPtN->Sumw2();
    fListOfObjects->Add(fPtN);
    
    fPtA = new TH1F("hPtA","hPt_Away; #it{p}_{T}; Entries",nPtBins,ptBins);
    fPtA->Sumw2();
    fListOfObjects->Add(fPtA);
    
    fPtT = new TH1F("hPtT","hPt_Transverse; #it{p}_{T}; Entries",nPtBins,ptBins);
    fPtT->Sumw2();
    fListOfObjects->Add(fPtT);
    
    fMultNMC = new TH1F("hMultNMC","hMult_Near; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultNMC->Sumw2();
    
    fMultAMC = new TH1F("hMultAMC","hMult_Away; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultAMC->Sumw2();
    
    fMultTMC = new TH1F("hMultTMC","hMult_Transverse; #it{N}_{ch}; Entries",131,-0.5,130.5);
    fMultTMC->Sumw2();
    
    fDphiN = new TH1F("hDphiN","hDphi_Near; #Delta#phi; Entries",2*64,-2*TMath::Pi(),2*TMath::Pi());
    fDphiN->Sumw2();
    fListOfObjects->Add(fDphiN);
    
    fDphiA = new TH1F("hDphiA","hDphi_Away; #Delta#phi; Entries",2*64,-2*TMath::Pi(),2*TMath::Pi());
    fDphiA->Sumw2();
    fListOfObjects->Add(fDphiA);
    
    fDphiT = new TH1F("hDphiT","hDphi_Transverse; #Delta#phi; Entries",2*64,-2*TMath::Pi(),2*TMath::Pi());
    fDphiT->Sumw2();
    fListOfObjects->Add(fDphiT);
    
    fDphi = new TH1F("hDphi","hDphi; #Delta#phi; Entries",2*64,-2*TMath::Pi(),2*TMath::Pi());
    fDphi->Sumw2();
    fListOfObjects->Add(fDphi);
    
   /*
   h: fMeanChT = 7.269
   i: fMeanChT = 7.257
   j: fMeanChT = 7.265
   k: fMeanChT = 7.261
   l: fMeanChT = 7.266
   o: fMeanChT = 7.211
   p: fMeanChT = 7.216
 
  */
    
    const int nBinsRT = 80;
    double binsRT[nBinsRT+1] = {0};
    
    for(int i = 0; i <= nBinsRT; ++i){
        binsRT[i] = (double)(i-1)/8.0;
    }
    
    fRT = new TH1D("hRT","hRT; R_{T}; Entries",nBinsRT,binsRT);
    fRT->Sumw2();
    fListOfObjects->Add(fRT);
    
    fPtLVsRt = new TH2D("hPtLVsRt", "; #it{p}^{L}_{T} (GeV/#it{c}); R_{T}",nPtBins,ptBins,nBinsRT,binsRT);
    fPtLVsRt->Sumw2();
    fListOfObjects->Add(fPtLVsRt);
    
    fRTMC = new TH1D("hRTMC","hRT; R_{T}; Entries",nBinsRT,binsRT);
    fRTMC->Sumw2();
    
    fPtLVsRtMC = new TH2D("hPtLVsRtMC", "; #it{p}^{L}_{T} (GeV/#it{c}); R_{T}",nPtBins,ptBins,nBinsRT,binsRT);
    fPtLVsRtMC->Sumw2();
    
    
    const Int_t nDeltaPiBins   = 80;
    const Double_t deltaPiLow  = 20;
    const Double_t deltaPiHigh = 100;
    
    const char *Region[3]      = {"Near","Away","Transverse"};
    const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
    
    
    const Char_t* ending[nHists] = {"02", "24", "46", "68"};
    const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };
    
    fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
    fcutDCAxy->SetParameter(0,0.0105);
    fcutDCAxy->SetParameter(1,0.0350);
    fcutDCAxy->SetParameter(2,1.1);
    
    fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
    fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);
    
    fEtaCalibration   = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
    fEtaCalibrationEl = new TF1("fDeDxVsEtaEl", "pol4", 0.0, 1.0);
    
    Int_t nPhiBins = 36;
    
    for(Int_t i = 0; i<nRt; ++i){
        
        for(Int_t r = 0; r<3; ++r){
            hMIPVsEta[r][i] = new TH2D(Form("hMIPVsEta_%s_%.1f_%.1f",Region[r],CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
            pMIPVsEta[r][i] = new TProfile(Form("pMIPVsEta_%s_%.1f_%.1f",Region[r],CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
            hPtAll[r][i] = new TH1D(Form("hPt_%s_%.1f_%.1f",Region[r],CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
            hPtAll[r][i]->Sumw2();
        }
        
        hMIPVsEtaV0s[i] = new TH2D(Form("hMIPVsEtaV0s_%.1f_%.1f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
        pMIPVsEtaV0s[i] = new TProfile(Form("pMIPVsEtaV0s_%.1f_%.1f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMin, fDeDxMIPMax);
        hPlateauVsEta[i] = new TH2D(Form("hPlateauVsEta_%.1f_%.1f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
        pPlateauVsEta[i] = new TProfile(Form("pPlateauVsEta_%.1f_%.1f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);
        hPhi[i] = new TH2D(Form("histPhi_%.1f_%.1f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);              
        hDCAxyVsPtPi[i] = new TH2D(Form("hDCAxyVsPtPi_%.1f_%.1f",CentMin[i],CentMax[i]), "hDCAxyVsPtPi; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
        hDCAxyVsPtPiC[i] = new TH2D(Form("hDCAxyVsPtPiC_%.1f_%.1f",CentMin[i],CentMax[i]), "hDCAxyVsPtPi Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
        hDCAxyVsPtPi[i]->Sumw2();
        hDCAxyVsPtPiC[i]->Sumw2();        
        hDCAxyVsPtp[i] = new TH2D(Form("hDCAxyVsPtp_%.1f_%.1f",CentMin[i],CentMax[i]), "hDCAxyVsPtp; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
        hDCAxyVsPtpC[i] = new TH2D(Form("hDCAxyVsPtpC_%.1f_%.1f",CentMin[i],CentMax[i]), "hDCAxyVsPtp Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
        hDCAxyVsPtp[i]->Sumw2();
        hDCAxyVsPtpC[i]->Sumw2();
        
        for(Int_t j=0; j<nHists; j++) {
                        
            for(Int_t r = 0; r<3; ++r){
                hDeDxVsP[r][i][j] = new TH2D( Form("hDeDxVsP_%s_%.1f_%.1f_%s",Region[r],CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, fdEdxHigh-fdEdxLow, fdEdxLow, fdEdxHigh);
                hDeDxVsP[r][i][j]->Sumw2();
            }
            
            hnSigmaPi[i][j] = new TH2D(Form("hnSigmaPi_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPi", nPtBins, ptBins, 40, -10, 10);
            hnSigmaPi[i][j]->Sumw2();
            
            hnSigmak[i][j] = new TH2D(Form("hnSigmak_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmak", nPtBins, ptBins, 40, -10, 10);
            hnSigmak[i][j]->Sumw2();
            
            hnSigmap[i][j] = new TH2D(Form("hnSigmap_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmap", nPtBins, ptBins, 40, -10, 10);
            hnSigmap[i][j]->Sumw2();
            
            hMIPVsPhi[i][j] = new TH2D(Form("hMIPVsPhi_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
            hMIPVsPhi[i][j]->Sumw2();
            
            pMIPVsPhi[i][j] = new TProfile(Form("pMIPVsPhi_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMin, fDeDxMIPMax);
            pMIPVsPhi[i][j]->Sumw2();
            
            hPlateauVsPhi[i][j] = new TH2D(Form("hPlateauVsPhi_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]),Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);
            hPlateauVsPhi[i][j]->Sumw2();
            
            pPlateauVsPhi[i][j] = new TProfile(Form("pPlateauVsPhi_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax, 95);
            pPlateauVsPhi[i][j]->Sumw2();
            
            histPiV0[i][j] = new TH2D(Form("hPiV0_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
            histPiV0[i][j]->Sumw2();
            
            histPV0[i][j] = new TH2D(Form("hPV0_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
            histPV0[i][j]->Sumw2();
            
            histPiTof[i][j] = new TH2D(Form("hPiTOF_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
            histPiTof[i][j]->Sumw2();
            
            histEV0[i][j] = new TH2D(Form("hEV0_%.1f_%.1f_%s",CentMin[i],CentMax[i],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
            histEV0[i][j]->Sumw2();
            
        }// eta loop
    } // centrality loop
    
    
    if(!fAnalysisMC){
        for(Int_t i=0; i<nRt; ++i ){
            
            for(Int_t r =0; r < 3; r++){
                fListOfObjects->Add(hMIPVsEta[r][i]);
                fListOfObjects->Add(pMIPVsEta[r][i]);
            }
            
            fListOfObjects->Add(hPlateauVsEta[i]);
            fListOfObjects->Add(pPlateauVsEta[i]);
            fListOfObjects->Add(hPhi[i]);
            
            for(Int_t j=0; j<nHists; ++j){
                fListOfObjects->Add(hMIPVsPhi[i][j]);
                fListOfObjects->Add(pMIPVsPhi[i][j]);
                fListOfObjects->Add(hPlateauVsPhi[i][j]);
                fListOfObjects->Add(pPlateauVsPhi[i][j]);
            }
            
            if(fMakePid){
                for(Int_t r =0; r < 3; r++)
                    fListOfObjects->Add(hPtAll[r][i]);
                if(fLowPt){
                    fListOfObjects->Add(hDCAxyVsPtPi[i]);
                    fListOfObjects->Add(hDCAxyVsPtp[i]);
                    
                    fListOfObjects->Add(hDCAxyVsPtPiC[i]);
                    fListOfObjects->Add(hDCAxyVsPtpC[i]);
                }
                
                for(Int_t j=0; j<nHists; ++j){
                    if(fLowPt){
                        fListOfObjects->Add(hnSigmaPi[i][j]);
                        fListOfObjects->Add(hnSigmak[i][j]);
                        fListOfObjects->Add(hnSigmap[i][j]);
                    }
                    
                    for(Int_t r =0; r < 3; r++)
                        fListOfObjects->Add(hDeDxVsP[r][i][j]);
                    
                    fListOfObjects->Add(histPiV0[i][j]);
                    fListOfObjects->Add(histPiTof[i][j]);
                    fListOfObjects->Add(histEV0[i][j]);
                    fListOfObjects->Add(histPV0[i][j]);
                }
            } //	if(MakePID)
        } //	Cent
        
        for(Int_t j=0; j<nHists; j++){
            hPtVsP[j] = new TH2D(Form("hPtVsP_%s",ending[j]),";#it{p} [GeV/c]; #it{p}_{T}",nPtBins,ptBins,nPtBins,ptBins);
            hPtVsP[j]->Sumw2();
            fListOfObjects->Add(hPtVsP[j]);
        }
    } //	!fAnalysisMC
    
    
    else{
        for(Int_t cent=0; cent<nRt; cent++) {
            for(Int_t pid=0; pid<7; pid++) {
                hMcIn[cent][pid] = new TH1D(Form("hIn_%s_%.1f_%.1f",Pid[pid],CentMin[cent],CentMax[cent]), Form("MC in (pid %s)", Pid[pid]),nPtBins,ptBins);
                hMcIn[cent][pid]->Sumw2();
                hMcOut[cent][pid] = new TH1D(Form("hMcOut_%s_%.1f_%.1f",Pid[pid],CentMin[cent],CentMax[cent]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
                hMcOut[cent][pid]->Sumw2();
                hDCApTPrim[cent][pid] = 0;
                hDCApTPrim[cent][pid] = new TH2D(Form("hDCA_%s_Prim_%d",Pid[pid],cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
                hDCApTPrim[cent][pid]->Sumw2();
                hDCApTWDec[cent][pid] = 0;
                hDCApTWDec[cent][pid] = new TH2D(Form("hDCA_%s_WDec_%d",Pid[pid],cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
                hDCApTWDec[cent][pid]->Sumw2();
                hDCApTMate[cent][pid] = 0;
                hDCApTMate[cent][pid] = new TH2D(Form("hDCA_%s_Mate_%d",Pid[pid],cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
                hDCApTMate[cent][pid]->Sumw2();
                
                hDCApTPrim2[cent][pid] = 0;
                hDCApTPrim2[cent][pid] = new TH2D(Form("hDCA2_%s_Prim_%d",Pid[pid],cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
                hDCApTPrim2[cent][pid]->Sumw2();
                hDCApTWDec2[cent][pid] = 0;
                hDCApTWDec2[cent][pid] = new TH2D(Form("hDCA2_%sWDec_%d",Pid[pid],cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
                hDCApTWDec2[cent][pid]->Sumw2();
                hDCApTMate2[cent][pid] = 0;
                hDCApTMate2[cent][pid] = new TH2D(Form("hDCA2_%sMate_%d",Pid[pid],cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
                hDCApTMate2[cent][pid]->Sumw2();
                
                fListOfObjects->Add(hMcIn[cent][pid]);
                fListOfObjects->Add(hMcOut[cent][pid]);
                fListOfObjects->Add(fMultNMC);
                fListOfObjects->Add(fMultAMC);
                fListOfObjects->Add(fMultTMC);
                //         fListOfObjects->Add(hDCApTPrim[cent][pid]);
                //         fListOfObjects->Add(hDCApTPrim2[cent][pid]);
                //         fListOfObjects->Add(hDCApTWDec[cent][pid]);
                //         fListOfObjects->Add(hDCApTWDec2[cent][pid]);
                //         fListOfObjects->Add(hDCApTMate[cent][pid]);
                //         fListOfObjects->Add(hDCApTMate2[cent][pid]);
                
            }	// pid MC
        }	//cent DCA MC
        
        fListOfObjects->Add(fRTMC);
        fListOfObjects->Add(fPtLVsRtMC);
        fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
        fListOfObjects->Add(fVtxMC);
        
    }
    
    
    fEventCuts.AddQAplotsToList(fListOfObjects);
    PostData(1, fListOfObjects);
    
}
//______________________________________________________________________________
void AliAnalysisTaskPPvsRT::UserExec(Option_t *)
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
    if (!utils)
    {
        cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
        return;
    }
    
    fEvents->Fill(0.5);
    
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
    if(!isINT7selected)
        return;
    fEvents->Fill(1.5);
    
    Int_t INEL = -1;
    INEL = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0);
    if( INEL < 1 )
        return;
    fEvents->Fill(2.5);
    
    if( utils->IsSPDClusterVsTrackletBG(fESD) )
        return;
    fEvents->Fill(3.5);
    
    if( fESD->IsPileupFromSPDInMultBins() )
        return;
    fEvents->Fill(4.5);
    
    if( fESD->IsIncompleteDAQ())
        return;
    fEvents->Fill(5.5);
    
    if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) )
        return;
    fEvents->Fill(6.5);
    
    if( !IsGoodZvertexPos(fESD) )
        return;
    fEvents->Fill(7.5);
    
    // Start Analysis
    
    if (fAnalysisType == "ESD"){
        
        
        AliESDtrack* LeadingTrack = GetLeadingTrack();
        if(!LeadingTrack)
            return;
        
        TObjArray* RegionsArray = new TObjArray(3);
        RegionsArray = SortRegions(LeadingTrack);

        TList* listToward = new TList();
        TList* listAway   = new TList();
        TList* listTrans  = new TList();

        listToward = (TList*)RegionsArray->At(0);
        listAway   = (TList*)RegionsArray->At(1);
        listTrans  = (TList*)RegionsArray->At(2);
        
        Double_t Rt = ((double)listTrans->GetEntries())/fMeanChT;
        
        fRT->Fill(Rt);
        fPtLVsRt->Fill(LeadingTrack->Pt(),Rt);
        
        int BinRT = GetBinRT(listTrans);
	if(BinRT < 0)
		return;
        ProduceArrayTrksESD(0,listToward,BinRT);
        ProduceArrayTrksESD(1,listAway,BinRT);
        ProduceArrayTrksESD(2,listTrans,BinRT);
        ProduceArrayV0ESD( fESD, 0);
        
        printf("********* Info about Event **********\n");
        printf(" %d Nch in the Near region\n",listToward->GetEntries());
        printf(" %d Nch in the Away region\n",listAway->GetEntries());
        printf(" %d Nch in Trans region\n",listTrans->GetEntries());
        printf("pT_Leading :: %f\n", LeadingTrack->Pt());
        cout << "Bin_RT    :: " << BinRT << endl;
        cout << "MeanChT   :: " << fMeanChT << endl;
        cout << "Period    :: " << fPeriod << endl;
        printf("*************************************\n");
        
        delete listToward;
        delete listAway;
        delete listTrans;
        delete RegionsArray;
        
        listToward   = 0x0;
        listAway     = 0x0;
        listTrans    = 0x0;
        RegionsArray = 0x0;
        
        if(fAnalysisMC)
            ProcessMCTruthESD();
        
    }
    
    PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
AliESDtrack* AliAnalysisTaskPPvsRT::GetLeadingTrack(){
    
    TObjArray* fTrks = new TObjArray();
    
    for(Int_t it = 0; it < fESD->GetNumberOfTracks(); ++it){
        
        AliESDtrack* esdTrack = fESD->GetTrack(it);
        if(!esdTrack)
            continue;
        
        if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
            continue;
        
        if(esdTrack->Pt() < 0.15)
            continue;
        
        UInt_t selectDebug = 0;
        if(fTrackFilterGolden){
            selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
            if (!selectDebug)
                continue;
        }
        
        Float_t dcaxy = 0.;
        Float_t dcaz = 0.;
        esdTrack->GetImpactParameters(dcaxy,dcaz);
        
        if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,esdTrack->Pt()) )
            continue;
        
        fTrks->Add(esdTrack);
        
    }
    
    Int_t nAccTrks = fTrks->GetEntriesFast();
    if(!nAccTrks)
        return 0x0;
    
    double *PtTrks = new double[nAccTrks];
    for(Int_t it=0; it<nAccTrks; it++){
        
        AliESDtrack* trk = (AliESDtrack*)(fTrks->At(it));
        if(!trk)
            continue;
        
        PtTrks[it] = trk->Pt();
    }
    
    Double_t LeadingPt = TMath::MaxElement(nAccTrks, PtTrks);
    AliESDtrack *LeadingTrk=0x0;
    
    for(Int_t it = 0; it < nAccTrks; ++it){
        
        AliESDtrack* trk = (AliESDtrack*)(fTrks->At(it));
        if(!trk)
            continue;
        
        if(trk->Pt() == LeadingPt){
            LeadingTrk = (AliESDtrack*)(fTrks->At(it));
            break;
        }
    }
     
    delete[] PtTrks;
    delete fTrks;
    PtTrks = 0x0;
    fTrks = 0x0;
    
    if(LeadingTrk->Pt() < fLeadingCut)
        return 0x0;
    
    else
        return LeadingTrk;
    
}
//_____________________________________________________________________________
TParticle* AliAnalysisTaskPPvsRT::GetLeadingTrackMC(TObjArray* TrksArray){
    
    if(!TrksArray)
        return 0x0;
    
    Int_t nAccTrks = TrksArray->GetEntriesFast();
    if(!nAccTrks)
        return 0x0;
    
    Double_t PtTrks[nAccTrks];
    for(Int_t it=0; it<nAccTrks; it++){
        
        TParticle* trk = (TParticle*)(TrksArray->At(it));
        if(!trk)
            continue;
        
        printf("pT = %f\n",trk->Pt());
        PtTrks[it] = trk->Pt();
    }
    
    Double_t LeadingPt = TMath::MaxElement(nAccTrks, PtTrks);
    printf("pT_Leading  = %f\n", LeadingPt);
    
    TParticle* LeadingTrk=0x0;
    for(Int_t it = 0; it < nAccTrks; ++it){
        TParticle* trk = (TParticle*)(TrksArray->At(it));
        if(!trk)
            continue;
        if(trk->Pt() == LeadingPt)
            LeadingTrk = (TParticle*)(TrksArray->At(it));
    }
    
    if(LeadingTrk->Pt() < fLeadingCut)
        return 0x0;
    
    else
        return LeadingTrk;
    
}
//_____________________________________________________________________________
TObjArray* AliAnalysisTaskPPvsRT::SortRegions(AliESDtrack* Ltrk)
{
    
    if(!Ltrk)
        return 0x0;
    
    TList *near       = new TList();
    TList *away       = new TList();
    TList *transverse = new TList();
    
    TObjArray* ArrayRegions = new TObjArray;
    ArrayRegions->AddLast(near);
    ArrayRegions->AddLast(away);
    ArrayRegions->AddLast(transverse);
    
    for(Int_t it = 0; it < fESD->GetNumberOfTracks(); ++it){
        
        AliESDtrack* esdTrack = fESD->GetTrack(it);
        if(!esdTrack)
            continue;
        
        if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
            continue;
        
        if(esdTrack->Pt() < 0.15)
            continue;
        
        if(Ltrk == esdTrack)
            continue;
        
        UInt_t selectDebug = 0;
        if(fTrackFilter){
            selectDebug = fTrackFilter->IsSelected(esdTrack);
            if (!selectDebug)
                continue;
        }
        
        Double_t dphi = DeltaPhi(esdTrack->Phi(),Ltrk->Phi());
        fDphi->Fill(dphi);
        
        if( TMath::Abs(dphi) < TMath::Pi()/3 ){
            near  ->Add(esdTrack);
            fPtN  ->Fill(esdTrack->Pt());
            fDphiN->Fill(dphi);
        }
        
        else if( TMath::Abs(dphi-TMath::Pi()) < TMath::Pi()/3 ){
            away  ->Add(esdTrack);
            fPtA  ->Fill(esdTrack->Pt());
            fDphiA->Fill(dphi);
        }
        else{
            transverse->Add(esdTrack);
            fPtT  ->Fill(esdTrack->Pt());
            fDphiT->Fill(dphi);
        }
        
        Int_t nh = -1;
        if(TMath::Abs(esdTrack->Eta())<0.2)
            nh = 0;
        else if(TMath::Abs(esdTrack->Eta())>=0.2 && TMath::Abs(esdTrack->Eta())<0.4)
            nh = 1;
        else if(TMath::Abs(esdTrack->Eta())>=0.4 && TMath::Abs(esdTrack->Eta())<0.6)
            nh = 2;
        else if(TMath::Abs(esdTrack->Eta())>=0.6 && TMath::Abs(esdTrack->Eta())<0.8)
            nh = 3;
        
        if(nh<0)
            continue;
        
        hPtVsP[nh]->Fill(esdTrack->P(),esdTrack->Pt());
    }
    
    fMultN->Fill(near->GetEntries());
    fMultA->Fill(away->GetEntries());
    fMultT->Fill(transverse->GetEntries());
    
    return ArrayRegions;
    
}
//_____________________________________________________________________________
TObjArray* AliAnalysisTaskPPvsRT::SortRegionsMC(TObjArray* TrksArray, TParticle* Ltrk){
    
    if(!TrksArray && !Ltrk)
        return 0x0;
    
    Int_t nTrks = TrksArray->GetEntriesFast();
    if(nTrks<=0)
        return 0x0;
    
    TList *near       = new TList();
    TList *away       = new TList();
    TList *transverse = new TList();
    
    TObjArray* ArrayRegions = new TObjArray;
    ArrayRegions->SetOwner(kTRUE);
    
    ArrayRegions->AddLast(near);
    ArrayRegions->AddLast(away);
    ArrayRegions->AddLast(transverse);
    
    for(Int_t it = 0; it < nTrks; ++it){
        
        TParticle* trk = (TParticle*)(TrksArray->At(it));
        if(!trk)
            continue;
        
        if(Ltrk == trk)
            continue;
        
        Double_t dphi = DeltaPhi(trk->Phi(),Ltrk->Phi());
        fDphi->Fill(dphi);
        
        if( TMath::Abs(dphi) < TMath::Pi()/3 ){
            near->Add(trk);
            fPtN  ->Fill(trk->Pt());
            fDphiN->Fill(dphi);
        }
        
        else if( TMath::Abs(dphi-TMath::Pi()) < TMath::Pi()/3 ){
            away->Add(trk);
            fPtA  ->Fill(trk->Pt());
            fDphiA->Fill(dphi);
        }
        else{
            transverse->Add(trk);
            fPtT  ->Fill(trk->Pt());
            fDphiT->Fill(dphi);
            //            printf("Dphi = %f\n",dphi);
        }
    }
    
    fMultNMC->Fill(near->GetEntries());
    fMultAMC->Fill(away->GetEntries());
    fMultTMC->Fill(transverse->GetEntries());
    
    return ArrayRegions;
    
    
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskPPvsRT::DeltaPhi(Double_t phi, Double_t Lphi,
                                           Double_t rangeMin, Double_t rangeMax)
{
    
    Double_t dphi = -999;
    Double_t pi = TMath::Pi();
    if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
    if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;
    
    if(phi < 0)          phi += 2*pi;
    else if(phi > 2*pi)  phi -= 2*pi;
    if(Lphi < 0)         Lphi += 2*pi;
    else if(Lphi > 2*pi) Lphi -= 2*pi;
    dphi = Lphi - phi;
    if (dphi < rangeMin)      dphi += 2*pi;
    else if (dphi > rangeMax) dphi -= 2*pi;
    
    return dphi;
}
//_____________________________________________________________________________
int AliAnalysisTaskPPvsRT::GetBinRT(TList* listT){
    
    double rt = -999.0;
    rt = (double)listT->GetEntries()/fMeanChT;
     int r = -10;
    
    if(rt < 0.5)
        r = 0;
    else if(rt >= 0.5 && rt < 1.5)
        r = 1;
    else if(rt >= 1.5 && rt < 2.5)
        r = 2;
    else if(rt >= 2.5 && rt < 5.0)
        r = 3;
    else
        r = -10;
    
    return r;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskPPvsRT::GetPidCode(Int_t pdgCode) const
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
void AliAnalysisTaskPPvsRT::ProcessMCTruthESD()
{
    
    TObjArray* fTrks        = new TObjArray();
    TObjArray* RegionsArray = new TObjArray();
    const Int_t nTracksMC = fMCStack->GetNtrack();
    
    for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
        
        TParticle* trackMC = fMCStack->Particle(iTracks);
        if(!trackMC)
            continue;
        
        if(!(fMCStack->IsPhysicalPrimary(iTracks)))
            continue;
        
        TParticlePDG* pdgPart = trackMC ->GetPDG();
        Double_t chargeMC = pdgPart->Charge();
        
        if(chargeMC==0)
            continue;
        
        if(TMath::Abs(trackMC->Eta()) > fEtaCut)
            continue;
        
        if(trackMC->Pt() < 0.15)
            continue;
        
        fTrks->Add(trackMC);
        
    }//MC track loop
    
    TParticle* Ltrk = GetLeadingTrackMC(fTrks);
    if(!Ltrk)
        return;
    
    RegionsArray = SortRegionsMC(fTrks,Ltrk);
    if(!RegionsArray)
        return;
    
    //    TList *listToward = (TList*)RegionsArray->At(0);
    //    TList *listAway   = (TList*)RegionsArray->At(1);
    TList *listTrans  = (TList*)RegionsArray->At(2);
    
    Double_t Rt = ((double)listTrans->GetEntries())/fMeanChT;
    
    fRTMC->Fill(Rt);
    fPtLVsRtMC->Fill(Ltrk->Pt(),Rt);
    
    Int_t BinRT = GetBinRT(listTrans);
    
    for (Int_t iTracks = 0; iTracks < listTrans->GetEntries(); iTracks++) {
        
        TParticle* trackMC = (TParticle*)listTrans->At(iTracks);
        if(!trackMC)
            continue;
        
        Int_t pdgCode = trackMC->GetPdgCode();
        Short_t pidCodeMC = 0;
        pidCodeMC = GetPidCode(pdgCode);
        
        hMcIn[BinRT][0]->Fill(trackMC->Pt());
        hMcIn[5][0]->Fill(trackMC->Pt());
        hMcIn[BinRT][pidCodeMC]->Fill(trackMC->Pt());
        hMcIn[5][pidCodeMC]->Fill(trackMC->Pt());
        
    }//MC track loop
    
    
    
}
//_____________________________________________________________________________
void AliAnalysisTaskPPvsRT::ProduceArrayTrksESD(const int& region, TList *lt, const int& Cent ){
    
    Int_t nESDTracks = lt->GetEntries();
    for(Int_t iT = 0; iT < nESDTracks; iT++) {
        
        AliESDtrack* esdTrack = (AliESDtrack*)lt->At(iT);
        
        if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
            continue;
        
        if(esdTrack->GetTPCsignalN() < fNcl)
            continue;
        
        if(esdTrack->Pt() < 0.15)
            continue;
        
        UInt_t selectDebug = 0;
        if(fTrackFilterGolden){
            selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
            if (!selectDebug)
                continue;
        }
        
        Double_t eta      = esdTrack->Eta();
        Double_t phi      = esdTrack->Phi();
        Double_t momentum = esdTrack->P();
        Double_t pt       = esdTrack->Pt();
        Float_t  dedx     = esdTrack->GetTPCsignal();
        Float_t  dedxUnc  = esdTrack->GetTPCsignal();
        
        Float_t dcaxy = 0.;
        Float_t dcaz = 0.;
        esdTrack->GetImpactParameters(dcaxy,dcaz);
        
        if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh))
            continue;
        
        if(fdEdxCalibrated){
            Int_t index = -1;
            index = GetIndex();
            dedx *= 50/EtaCalibration(index,eta);
        }
        
        Short_t pidCode     = 0;
        if(fAnalysisMC) {
            const Int_t label = TMath::Abs(esdTrack->GetLabel());
            TParticle* mcTrack = 0;
            mcTrack = fMCStack->Particle(label);
            
            if(mcTrack){
                if( esdTrack->Charge()==0 )
                    continue;
                
                Int_t pdgCode = mcTrack->GetPdgCode();
                pidCode = GetPidCode(pdgCode);
                
                if( fMCStack->IsPhysicalPrimary(label) ){
                    if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
                        hMcOut[Cent][0]->Fill(esdTrack->Pt());
                        hMcOut[4][0]->Fill(esdTrack->Pt());
                        hMcOut[Cent][pidCode]->Fill(esdTrack->Pt());
                        hMcOut[4][pidCode]->Fill(esdTrack->Pt());
                    }
                    
                    hDCApTPrim[Cent][0]->Fill(pt,dcaxy);
                    hDCApTPrim[4][0]->Fill(pt,dcaxy);
                    hDCApTPrim[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTPrim[4][pidCode]->Fill(pt,dcaxy);
                    
                    hDCApTPrim2[Cent][0]->Fill(pt,dcaxy);
                    hDCApTPrim2[4][0]->Fill(pt,dcaxy);
                    hDCApTPrim2[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTPrim2[4][pidCode]->Fill(pt,dcaxy);
                }	// Primary particles MC
                
                if( fMCStack->IsSecondaryFromWeakDecay(label) ){
                    hDCApTWDec[Cent][0]->Fill(pt,dcaxy);
                    hDCApTWDec[4][0]->Fill(pt,dcaxy);
                    hDCApTWDec[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTWDec[4][pidCode]->Fill(pt,dcaxy);
                    
                    hDCApTWDec2[Cent][0]->Fill(pt,dcaxy);
                    hDCApTWDec2[4][0]->Fill(pt,dcaxy);
                    hDCApTWDec2[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTWDec2[4][pidCode]->Fill(pt,dcaxy);
                }	// Weak Decay MC
                
                if( fMCStack->IsSecondaryFromMaterial(label) ){
                    hDCApTMate[Cent][0]->Fill(pt,dcaxy);
                    hDCApTMate[4][0]->Fill(pt,dcaxy);
                    hDCApTMate[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTMate[4][pidCode]->Fill(pt,dcaxy);
                    
                    hDCApTMate2[Cent][0]->Fill(pt,dcaxy);
                    hDCApTMate2[4][0]->Fill(pt,dcaxy);
                    hDCApTMate2[Cent][pidCode]->Fill(pt,dcaxy);
                    hDCApTMate2[4][pidCode]->Fill(pt,dcaxy);
                }	// Material Inte MC
            }	//mcTrack
        }	//fAnalysis MC
        
        // DCAxy cut
        if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
            continue;
        
        Bool_t IsTOFout=kFALSE;
        if((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
            IsTOFout=kTRUE;
        Float_t lengthtrack=esdTrack->GetIntegratedLength();
        Float_t timeTOF=esdTrack->GetTOFsignal();
        Double_t inttime[5]={0,0,0,0,0};
        esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
        Float_t beta = -99;
        if(!IsTOFout){
            if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
                beta = inttime[0] / timeTOF;
        }
        
        if(!fdEdxCalibrated){
            if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
                if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                    hMIPVsEta[region][Cent]->Fill(eta,dedxUnc);
                    hMIPVsEta[region][4]->Fill(eta,dedxUnc);
                    pMIPVsEta[region][Cent]->Fill(eta,dedxUnc);
                    pMIPVsEta[region][4]->Fill(eta,dedxUnc);
                }
                if( dedxUnc > 70 && dedxUnc < 90 ){
                    if(TMath::Abs(beta-1)<0.1){
                        hPlateauVsEta[Cent]->Fill(eta,dedxUnc);
                        hPlateauVsEta[4]->Fill(eta,dedxUnc);
                        pPlateauVsEta[Cent]->Fill(eta,dedxUnc);
                        pPlateauVsEta[4]->Fill(eta,dedxUnc);
                    }
                }
            }
        }
        else{
            if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
                if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                    hMIPVsEta[region][Cent]->Fill(eta,dedx);
                    hMIPVsEta[region][4]->Fill(eta,dedx);
                    pMIPVsEta[region][Cent]->Fill(eta,dedx);
                    pMIPVsEta[region][4]->Fill(eta,dedx);
                }
                if( dedxUnc > 70 && dedxUnc < 90 ){
                    if(TMath::Abs(beta-1)<0.1){
                        hPlateauVsEta[Cent]->Fill(eta,dedx);
                        hPlateauVsEta[4]->Fill(eta,dedx);
                        pPlateauVsEta[Cent]->Fill(eta,dedx);
                        pPlateauVsEta[4]->Fill(eta,dedx);
                    }
                }
            }
        }
        
        hPtAll[region][Cent]->Fill(pt);
        hPtAll[region][4]->Fill(pt);
        
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
        
        if(beta>1){
            histPiTof[Cent][nh]->Fill(momentum, dedx);
            histPiTof[4][nh]->Fill(momentum, dedx);
        }
        
        if( momentum <= 0.6 && momentum >= 0.4  ){
            if( dedx < fDeDxMIPMax && dedx > fDeDxMIPMin ){
                hMIPVsPhi[Cent][nh]->Fill(phi,dedx);
                hMIPVsPhi[4][nh]->Fill(phi,dedx);
                pMIPVsPhi[Cent][nh]->Fill(phi,dedx);
                pMIPVsPhi[4][nh]->Fill(phi,dedx);
            }
            if( dedx > 70 && dedx < 90 ){
                if(TMath::Abs(beta-1)<0.1){
                    hPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
                    hPlateauVsPhi[4][nh]->Fill(phi,dedx);
                    pPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
                    pPlateauVsPhi[4][nh]->Fill(phi,dedx);
                }
            }
        }
        
        hDeDxVsP[region][Cent][nh]->Fill(momentum,dedx);
        hDeDxVsP[region][4][nh]->Fill(momentum,dedx);
        
        hnSigmaPi[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
        hnSigmaPi[4][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
        hnSigmak[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
        hnSigmak[4][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
        hnSigmap[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
        hnSigmap[4][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
        
        
    }//end of track loop
    
    
}
//________________________________________________________________________
void AliAnalysisTaskPPvsRT::ProduceArrayV0ESD( AliESDEvent *ESDevent, const int& Cent ){
    
    Int_t nv0s = ESDevent->GetNumberOfV0s();
    
    //	fcentAfterV0s->Fill(Cent+1);
    //	fcentAfterV0s->Fill(11);
    
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
        if ( !esdV0 )
            continue;
        
        //check onfly status
        if( !esdV0->GetOnFlyStatus() )
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
                        
                        if(track->GetTPCsignalN()<70)continue;
                        Double_t phi     = track->Phi();
                        
                        if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh))
                            continue;
                        
                        Double_t eta      = track->Eta();
                        Double_t momentum = track->P();
                        Double_t dedx     = track->GetTPCsignal();
                        Double_t dedxUnc  = track->GetTPCsignal();
                        
                        
                        if(fdEdxCalibrated){
                            Int_t index = -1;
                            index = GetIndex();
                            dedx   *= 50.0/EtaCalibration(index,eta);
                        }
                        
                        if(fillPos&&fillNeg){
                            if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
                                if(momentum<0.6&&momentum>0.4){
                                    hMIPVsEtaV0s[Cent]->Fill(eta,dedx);
                                    hMIPVsEtaV0s[4]->Fill(eta,dedx);
                                    pMIPVsEtaV0s[Cent]->Fill(eta,dedx);
                                    pMIPVsEtaV0s[4]->Fill(eta,dedx);
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
                            histPiV0[Cent][nh]->Fill(momentum, dedx);
                            histPiV0[4][nh]->Fill(momentum, dedx);
                        }
                        else{
                            histPV0[Cent][nh]->Fill(momentum, dedx);
                            histPV0[4][nh]->Fill(momentum, dedx);
                        }
                    }//end loop over two tracks
                    
                };
                    break;
                    
                case 1:{//gammas
                    
                    Bool_t fillPos = kFALSE;
                    Bool_t fillNeg = kFALSE;
                    
                    
                    if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
                        if( dmassG<0.01 && dmassG>0.0001 ) {
                            Int_t index = -1;
                            index = GetIndex();
                            if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationEl(index,nTrack->Eta())) < 5)
                                fillPos = kTRUE;
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
                        Int_t index = -1;
                        index = GetIndex();
                        dedx *= 50/EtaCalibration(index,eta);
                    }
                    
                    if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh))
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
                    
                    histEV0[Cent][nh]->Fill(momentum, dedx);
                    histEV0[4][nh]->Fill(momentum, dedx);
                    
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
Bool_t AliAnalysisTaskPPvsRT::selectVertex2015pp(AliESDEvent *esd,
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
Bool_t AliAnalysisTaskPPvsRT::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{
    
    if( !spdVertex ) return kFALSE;
    if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPPvsRT::IsGoodZvertexPos(AliESDEvent *esd)
{
    
    if( !esd ) return kFALSE;
    //Cut on the vertex z position
    const AliESDVertex * vertex = esd->GetPrimaryVertex();
    if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPPvsRT::PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh)
{
    if(pt < 2.0)
        return kTRUE;
    
    if(mag < 0)    // for negatve polarity field
        phi = TMath::TwoPi() - phi;
    if(q < 0) // for negatve charge
        phi = TMath::TwoPi()-phi;
    
    phi += TMath::Pi()/18.0; // to center gap in the middle
    phi = fmod(phi, TMath::Pi()/9.0);
    
    if(phi<phiCutHigh->Eval(pt)
       && phi>phiCutLow->Eval(pt))
        return kFALSE; // reject track
    
    hPhi[4]->Fill(pt, phi);
    
    return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskPPvsRT::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){
    
    Double_t maxDCAxy = 10;
    maxDCAxy = fMaxDCAxy->Eval(ptI);
    return maxDCAxy;
    
}
//________________________________________________________________________
void AliAnalysisTaskPPvsRT::SetTrackCuts(AliAnalysisFilter* fTrackFilter){
    
    AliESDtrackCuts* esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    
    fTrackFilter->AddCuts(esdTrackCuts);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPPvsRT::EtaCalibration( const Int_t &indx, const Double_t &eta){
                                    //    h        i         j        l      k        o          p
    const Double_t aPos[nRt+2]      = {49.9044 ,50.0841  ,49.8419 ,49.9799 ,49.9659 ,50.0535 ,50.0649};
    const Double_t bPos[nRt+2]      = {4.05075 ,-0.743724,10.3952 ,2.99619 ,2.91366 ,-2.87404,-4.3589};
    const Double_t cPos[nRt+2]      = {-58.1027,-13.8508 ,-151.227,-45.718 ,-45.5994,31.159  ,66.268 };
    const Double_t dPos[nRt+2]      = {342.297 ,141.269  ,900.83  ,290.013 ,290.042 ,-151.257,-435.715};
    const Double_t ePos[nRt+2]      = {-1098.19,-567.054 ,-2833.25,-1018.42,-1014.49,282.703 ,1325.24};
    const Double_t fPos[nRt+2]      = {1944.32 ,1105.63  ,4884.59 ,1948.68 ,1931.84 ,-98.8756,-2030.25};
    const Double_t gPos[nRt+2]      = {-1749.14,-1022.19 ,-4318.99,-1864.06,-1839.36,-230.114,1542.09};
    const Double_t hPos[nRt+2]      = {617.929 ,355.158  ,1521.66 ,692.752 ,680.421 ,172.854 ,-468.577};
    
    const Double_t aNeg[nRt+2]      = {49.9261,50.0561,49.9583,50.078 ,50.046 ,49.9496,50.1258};
    const Double_t bNeg[nRt+2]      = {2.92422,4.68965,3.38038,6.67199,6.79992,2.45301,12.7977};
    const Double_t cNeg[nRt+2]      = {61.6661,65.891 ,53.1256,103.662,109.86 ,53.654 ,190.076};
    const Double_t dNeg[nRt+2]      = {421.545,394.542,314.489,611.034,668.241,363.689,1144.11};
    const Double_t eNeg[nRt+2]      = {1283.04,1100.56,825.296,1695.63,1916.44,1115.13,3411.98};
    const Double_t fNeg[nRt+2]      = {1944.85,1516.21,1021.01,2395.88,2815.04,1762.18,5402.99};
    const Double_t gNeg[nRt+2]      = {1442.98,989.24 ,548.44 ,1669.22,2057.21,1421.46,4379.16};
    const Double_t hNeg[nRt+2]      = {419.491,238.333,84.7945,455.362,595.391,469.45 ,1436.76};
    
    
    for(Int_t i=0; i<8; ++i)
        fEtaCalibration->SetParameter(i,0);
    
    if(eta<0){
        fEtaCalibration->SetParameter(0,aNeg[indx]);
        fEtaCalibration->SetParameter(1,bNeg[indx]);
        fEtaCalibration->SetParameter(2,cNeg[indx]);
        fEtaCalibration->SetParameter(3,dNeg[indx]);
        fEtaCalibration->SetParameter(4,eNeg[indx]);
        fEtaCalibration->SetParameter(5,fNeg[indx]);
        fEtaCalibration->SetParameter(6,gNeg[indx]);
        fEtaCalibration->SetParameter(7,hNeg[indx]);
    }
    else{
        fEtaCalibration->SetParameter(0,aPos[indx]);
        fEtaCalibration->SetParameter(1,bPos[indx]);
        fEtaCalibration->SetParameter(2,cPos[indx]);
        fEtaCalibration->SetParameter(3,dPos[indx]);
        fEtaCalibration->SetParameter(4,ePos[indx]);
        fEtaCalibration->SetParameter(5,fPos[indx]);
        fEtaCalibration->SetParameter(6,gPos[indx]);
        fEtaCalibration->SetParameter(7,hPos[indx]);
    }
    
    return fEtaCalibration->Eval(eta);
    
}
//________________________________________________________________________
Double_t AliAnalysisTaskPPvsRT::EtaCalibrationEl(const Int_t &indx, const Double_t &eta){
    
    const Double_t aPosEl[nRt+2]    = {79.8647 ,79.6737 ,80.3915 ,80.1263 ,79.9957 ,79.6537 ,80.6434 };
    const Double_t bPosEl[nRt+2]    = {6.50512 ,16.0745 ,9.53925 ,5.28525 ,7.03079 ,15.0221 ,0.40293 };
    const Double_t cPosEl[nRt+2]    = {-35.9277,-80.5639,-69.3773,-32.7731,-42.9098,-83.6391,-21.8162};
    const Double_t dPosEl[nRt+2]    = {73.1535 ,148.866 ,143.956 ,68.4524 ,88.7057 ,168.5   ,61.9147 };
    const Double_t ePosEl[nRt+2]    = {-47.1041,-90.3376,-89.5518,-44.1566,-56.6554,-107.999,-44.6593};
    
    const Double_t aNegEl[nRt+2]    = {79.6366 ,80.0767 ,79.6157 ,79.8351 ,79.7387 ,79.3638 ,79.9111 };
    const Double_t bNegEl[nRt+2]    = {-11.3437,-2.51009,-16.2468,-8.46921,-8.60021,-17.1977,-1.66066};
    const Double_t cNegEl[nRt+2]    = {-65.1353,-23.6188,-92.0783,-44.5947,-44.1718,-82.7998,-6.96109};
    const Double_t dNegEl[nRt+2]    = {-134.447,-65.5053,-180.753,-86.2242,-84.4984,-143.394,-16.0465};
    const Double_t eNegEl[nRt+2]    = {-87.7848,-51.1463,-112.997,-53.6285,-51.945 ,-81.3439,-10.3587};
    
    
    for(Int_t i=0; i<5; ++i)
        fEtaCalibrationEl->SetParameter(i,0);
    
    if(eta<0){
        fEtaCalibrationEl->SetParameter(0,aNegEl[indx]);
        fEtaCalibrationEl->SetParameter(1,bNegEl[indx]);
        fEtaCalibrationEl->SetParameter(2,cNegEl[indx]);
        fEtaCalibrationEl->SetParameter(3,dNegEl[indx]);
        fEtaCalibrationEl->SetParameter(4,eNegEl[indx]);
    }
    else{
        fEtaCalibrationEl->SetParameter(0,aPosEl[indx]);
        fEtaCalibrationEl->SetParameter(1,bPosEl[indx]);
        fEtaCalibrationEl->SetParameter(2,cPosEl[indx]);
        fEtaCalibrationEl->SetParameter(3,dPosEl[indx]);
        fEtaCalibrationEl->SetParameter(4,ePosEl[indx]);
    }
    
    return fEtaCalibrationEl->Eval(eta);
    
}
//________________________________________________________________________
Int_t AliAnalysisTaskPPvsRT::GetIndex()
{
    
    Int_t indx = -1;
    
    if(fPeriod=="h")
        indx = 0;
    else if(fPeriod=="i")
        indx = 1;
    else if(fPeriod=="j")
        indx = 2;
    else if(fPeriod=="l")
        indx = 3;
    else if(fPeriod=="k")
        indx = 4;
    else if(fPeriod=="o")
        indx = 5;
    else
        indx = 6;
    
    return indx;
    
}
