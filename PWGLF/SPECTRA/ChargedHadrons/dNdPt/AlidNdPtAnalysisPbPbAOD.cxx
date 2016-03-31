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
//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPbAOD class.
//
// Author: P. Luettig, 15.05.2013
// last modified: 10.06.2014
//------------------------------------------------------------------------------
/*
 * This task analysis measured data in PbPb collisions stored in AODs and extract
 * transverse momentum spectra for unidentified charged hadrons vs. centrality.
 * Based on MC the efficiency and secondary contamination are determined,
 * to correct the measured pT distribution.
 * Histograms for the pT resolution correction are also filled.
 *
 */


#include "AlidNdPtAnalysisPbPbAOD.h"

#include "AliAnalysisTaskSE.h"
#include <RVersion.h>
#include <iostream>

using namespace std;

ClassImp(AlidNdPtAnalysisPbPbAOD)

AlidNdPtAnalysisPbPbAOD::AlidNdPtAnalysisPbPbAOD(const char *name) : AliAnalysisTaskSE(name),
fOutputList(0),
// Histograms
fPt(0),
fMCPt(0),
fZvPtEtaCent(0),
fDeltaphiPtEtaPhiCent(0),
fDeltaphiPtEtaPhiZvCent(0),
fPtResptCent(0),
fPtResptptCent(0),
fPtEvent(0),
fMCRecPrimZvPtEtaCent(0),
fMCGenZvPtEtaCent(0),
fMCRecSecZvPtEtaCent(0),
fMCPtEtaPhiCent(0),
fMCRecPrimPtEtaPhiCent(0),
fMCGenPtEtaPhiCent(0),
fMCRecSecPtEtaPhiCent(0),
fMCPtEtaPhiZvCent(0),
fMCRecPrimPtEtaPhiZvCent(0),
fMCGenPtEtaPhiZvCent(0),
fMCRecSecPtEtaPhiZvCent(0),
fEventStatistics(0),
fEventStatisticsCentrality(0),
fMCEventStatisticsCentrality(0),
fAllEventStatisticsCentrality(0),
fEventStatisticsCentralityTrigger(0),
fZvMultCent(0),
fTriggerStatistics(0),
fCharge(0),
fMCCharge(0),
fDCAPtAll(0),
fDCAPtAccepted(0),
fMCDCAPtSecondary(0),
fMCDCAPtPrimary(0),
fCutPercClusters(0),
fCutPercCrossed(0),
fCrossCheckRowsLength(0),
fCrossCheckClusterLength(0),
fCrossCheckRowsLengthAcc(0),
fCrossCheckClusterLengthAcc(0),
fCrossCheckPtresLength(0),
fCrossCheckPtresRows(0),
fCutSettings(0),
fEventplaneDist(0),
fEventplaneRunDist(0),
fMCEventplaneDist(0),
fCorrelEventplaneMCDATA(0),
fCorrelEventplaneDefaultCorrected(0),
fEventplaneSubtractedPercentage(0),
fChargeOverPtRuns(0),
fVZEROMultCentrality(0),
fVEROMultRefMult(0),
// cross check for event plane resolution
fEPDistCent(0),
fPhiCent(0),
fPcosEPCent(0),
fPsinEPCent(0),
fPcosPhiCent(0),
fPsinPhiCent(0),
fEPContributionDifference(0),
// cross check for event plane determination
fDeltaPhiCent(0),
fDeltaPhiSymCent(0),
fMCRecTracksMult(0),
fMCGenTracksMult(0),
fCrossCheckFilterBitPhiCent(0),
fTriggerStringsFired(0),
fTriggerStringComplete(0),
//vertex histograms
fVertexZ(0),
fVertexZSPD(0),
fVertexZTPC(0),
fDeltaVertexZGlobalSPD(0),
fDeltaVertexZGlobalTPC(0),
fVertexContributors(0),
fVertexZAfterCuts(0),
fVertexZSPDAfterCuts(0),
fVertexZTPCAfterCuts(0),
fDeltaVertexZGlobalSPDAfterCuts(0),
fDeltaVertexZGlobalTPCAfterCuts(0),
fVertexContributorsAfterCuts(0),
fVertexContributorsAfterCutsCent(0),
fVertexContributorsAfterCutsSemi(0),
fVertexContributorsAfterCutsMB(0),
//global
fIsMonteCarlo(0),
fEventNumberForPtSpectra(0),
fFillEventPtSpectraHistogram(kFALSE),
fEPselector("Q"),
fCentEstimator("V0M"),
fDoMinBiasAnalysis(kFALSE),
fDisabledTriggerString(""),
// event cut variables
fCutMaxZVertex(10.),
fVertexMinContributors(1),
// track kinematic cut variables
fCutPtMin(0.15),
fCutPtMax(200.),
fCutEtaMin(-0.8),
fCutEtaMax(0.8),
// anchor point
fAnchorPointCorrectionFactor(1),
fDoAnchorPointSystStudy(kFALSE),
// track quality cut variables
fFilterBit(AliAODTrack::kTrkGlobal),
fHybridTracking(kFALSE),
fUseRelativeCuts(kFALSE),
fCutRequireTPCRefit(kTRUE),
fCutRequireITSRefit(kTRUE),
fCutMinNumberOfClusters(60),
fCutPercMinNumberOfClusters(0.2),
fCutMinNumberOfCrossedRows(120.),
fCutPercMinNumberOfCrossedRows(0.2),
fCutMinRatioCrossedRowsOverFindableClustersTPC(0.8),
fCutMaxChi2PerClusterTPC(4.),
fCutMaxFractionSharedTPCClusters(0.4),
fCutMaxDCAToVertexZ(3.0),
fCutMaxDCAToVertexXY(3.0),
fCutMaxChi2PerClusterITS(36.),
fCutDCAToVertex2D(kFALSE),
fCutRequireSigmaToVertex(kFALSE),
fCutMaxDCAToVertexXYPtDepPar0(0.0182),
fCutMaxDCAToVertexXYPtDepPar1(0.0350),
fCutMaxDCAToVertexXYPtDepPar2(1.01),
fCutAcceptKinkDaughters(kFALSE),
fCutMaxChi2TPCConstrainedGlobal(36.),
fCutLengthInTPCPtDependent(kFALSE),
fPrefactorLengthInTPCPtDependent(1),
fCrossCheckCorrelHisto(kFALSE),
// binning for THnSparse
fMultNbins(0),
fMultFineNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fPtCheckNbins(0),
fEtaNbins(0),
fEtaCheckNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fPhiNbins(0),
fPhiCorrNbins(0),
fDeltaphiNbins(0),
fRunNumberNbins(0),
fBinsMult(0),
fBinsMultFine(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsPtCheck(0),
fBinsEta(0),
fBinsEtaCheck(0),
fBinsZv(0),
fBinsCentrality(0),
fBinsPhi(0),
fBinsPhiCorr(0),
fBinsDeltaphi(0),
fBinsRunNumber(0)
{
  
  for(Int_t i = 0; i < cqMax; i++)
  {
    fCrossCheckAll[i] = 0;
    fCrossCheckAcc[i] = 0;
  }
  
  fMultNbins = 0;
  fMultFineNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fPtCheckNbins = 0;
  fEtaNbins = 0;
  fEtaCheckNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fPhiNbins = 0;
  fPhiCorrNbins = 0;
  fDeltaphiNbins = 0;
  fRunNumberNbins = 0;
  fBinsMult = 0;
  fBinsMultFine = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsPtCheck = 0;
  fBinsEta = 0;
  fBinsEtaCheck = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  fBinsPhi = 0;
  fBinsPhiCorr = 0;
  fBinsDeltaphi = 0;
  fBinsRunNumber = 0;
  
  fEventNumberForPtSpectra = 0;
  
  DefineOutput(1, TList::Class());
}

// destructor
AlidNdPtAnalysisPbPbAOD::~AlidNdPtAnalysisPbPbAOD()
{
  //
  //  because task is owner of the output list, all objects are deleted, when list->Clear() is called
  //
  if(fOutputList)
  {
    fOutputList->Clear();
    delete fOutputList;
  }
  fOutputList = 0;
}

void AlidNdPtAnalysisPbPbAOD::UserCreateOutputObjects()
{
  // create all output histograms here
  OpenFile(1, "RECREATE");
  
  fOutputList = new TList();
  fOutputList->SetOwner();
  
  //define default binning
  Double_t binsMultDefault[48] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,19.5, 20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5,200.5, 300.5, 400.5, 500.5, 600.5, 700.5, 800.5, 900.5, 1000.5, 2000.5, 3000.5, 4000.5, 5000.5, 6000.5, 7000.5, 8000.5, 9000.5, 10000.5 };
  Double_t binsMultDefaultFine[160] = {	-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,19.5, 20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5, 110.5, 120.5, 130.5, 140.5, 150.5, 160.5, 170.5, 180.5, 190.5, 200.5, 210.5, 220.5, 230.5, 240.5, 250.5, 260.5, 270.5, 280.5, 290.5, 300.5, 310.5, 320.5, 330.5, 340.5, 350.5, 360.5, 370.5, 380.5, 390.5, 400.5, 410.5, 420.5, 430.5, 440.5, 450.5, 460.5, 470.5, 480.5, 490.5, 500.5, 550.5, 600.5, 650.5, 700.5, 750.5, 800.5, 850.5, 900.5, 950.5, 1000.5, 1050.5, 1100.5, 1150.5, 1200.5, 1250.5, 1300.5, 1350.5, 1400.5, 1450.5, 1500.5, 1550.5, 1600.5, 1650.5, 1700.5, 1750.5, 1800.5, 1850.5, 1900.5, 1950.5, 2000.5, 2050.5, 2100.5, 2150.5, 2200.5, 2250.5, 2300.5, 2350.5, 2400.5, 2450.5, 2500.5, 2550.5, 2600.5, 2650.5, 2700.5, 2750.5, 2800.5, 2850.5, 2900.5, 2950.5, 3000.5, 3050.5, 3100.5, 3150.5, 3200.5, 3250.5, 3300.5, 3350.5, 3400.5, 3450.5, 3500.5, 3550.5, 3600.5, 3650.5, 3700.5, 3750.5, 3800.5, 3850.5, 3900.5, 3950.5, 4000.5, 4050.5, 4100.5, 4150.5, 4200.5, 4250.5, 4300.5, 4350.5, 4400.5, 4450.5, 4500.5, 4550.5, 4600.5, 4650.5, 4700.5, 4750.5, 4800.5, 4850.5, 4900.5, 4950.5, 5000.5};
  Double_t binsPtDefault[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  Double_t binsPtCorrDefault[51] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 200.0};
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[7] = {-30.,-10.,-5.,0.,5.,10.,30.};
  Double_t binsCentralityDefault[12] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
  
  Double_t binsPhiDefault[37] = { 0., 0.174533, 0.349066, 0.523599, 0.698132, 0.872665, 1.0472, 1.22173, 1.39626, 1.5708, 1.74533, 1.91986, 2.0944, 2.26893, 2.44346, 2.61799, 2.79253, 2.96706, 3.14159, 3.31613, 3.49066, 3.66519, 3.83972, 4.01426, 4.18879, 4.36332, 4.53786, 4.71239, 4.88692, 5.06145, 5.23599, 5.41052, 5.58505, 5.75959, 5.93412, 6.10865, 2.*TMath::Pi()};
  
  Double_t binsPhiCorrDefault[19] = { 0, 1./9.*TMath::Pi(), 2./9.*TMath::Pi(), 3./9.*TMath::Pi(), 4./9.*TMath::Pi(), 5./9.*TMath::Pi(), 6./9.*TMath::Pi(), 7./9.*TMath::Pi(), 8./9.*TMath::Pi(), 9./9.*TMath::Pi(), 10./9.*TMath::Pi(), 11./9.*TMath::Pi(), 12./9.*TMath::Pi(), 13./9.*TMath::Pi(), 14./9.*TMath::Pi(), 15./9.*TMath::Pi(), 16./9.*TMath::Pi(), 17./9.*TMath::Pi(), 18./9.*TMath::Pi() };
  
  Double_t binsDeltaphiDefault[9] = {  0, 1./16.*TMath::Pi(), 2./16.*TMath::Pi(), 3./16.*TMath::Pi(), 4./16.*TMath::Pi(), 5./16.*TMath::Pi(), 6./16.*TMath::Pi(), 7./16.*TMath::Pi(), 8./16.*TMath::Pi()};
  
  Double_t binsPtCheckDefault[8] = {0.,0.15,1.0,5.0, 10.0, 20.0, 50.0, 200.0};
  Double_t binsEtaCheckDefault[7] = {-1.0,-0.8,-0.4,0.,0.4,0.8,1.0};
  
  Double_t binsRunNumbers2011[186] = {
    167693, 167706, 167711, 167712, 167713, 167806, 167807, 167808, 167813, 167814, 167818, 167841, 167842, 167844, 167846, 167902, 167903, 167909, 167915, 167920, 167921, 167985, 167986, 167987, 167988, 168066, 168068, 168069, 168076, 168103, 168104, 168105, 168107, 168108, 168115, 168171, 168172, 168173, 168175, 168177, 168181, 168203, 168204, 168205, 168206, 168207, 168208, 168212, 168213, 168310, 168311, 168318, 168322, 168325, 168341, 168342, 168356, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168644, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169236, 169238, 169377, 169382, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169584, 169586, 169587, 169588, 169590, 169591, 169628, 169683, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169914, 169918, 169919, 169920, 169922, 169923, 169924, 169926, 169956, 169961, 169965, 169969, 169975, 169981, 170027, 170036, 170038, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170162, 170163, 170193, 170195, 170203, 170204, 170205, 170207, 170208, 170228, 170230, 170264, 170267, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170374, 170387, 170388, 170389, 170390, 170546, 170552, 170556, 170572, 170593, 170593+1
  };
  
  // if no binning is set, use the default
  if (!fBinsMult)		{ SetBinsMult(48,binsMultDefault); }
  if (!fBinsMultFine)	{ SetBinsMultFine(160,binsMultDefaultFine); }
  if (!fBinsPt)			{ SetBinsPt(82,binsPtDefault); }
  if (!fBinsPtCorr)		{ SetBinsPtCorr(51,binsPtCorrDefault); }
  if (!fBinsPtCheck)	{ SetBinsPtCheck(8,binsPtCheckDefault); }
  if (!fBinsEta)		{ SetBinsEta(31,binsEtaDefault); }
  if (!fBinsEtaCheck)	{ SetBinsEtaCheck(7,binsEtaCheckDefault); }
  if (!fBinsZv)			{ SetBinsZv(7,binsZvDefault); }
  if (!fBinsCentrality)	{ SetBinsCentrality(12,binsCentralityDefault); }
  if (!fBinsPhi)		{ SetBinsPhi(37,binsPhiDefault); }
  if (!fBinsPhiCorr)		{ SetBinsPhiCorr(19,binsPhiCorrDefault); }
  if (!fBinsDeltaphi)	{ SetBinsDeltaphi(9,binsDeltaphiDefault); }
  if (!fBinsRunNumber) 	{ SetBinsRunNumber(186, binsRunNumbers2011); }
  
  Int_t binsZvPtEtaCent[4]={fZvNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Double_t minbinsZvPtEtaCent[4]={-30.,0,-1.5,0};
  Double_t maxbinsZvPtEtaCent[4]={30  ,200,1.5,100};
  
  Int_t binsZvPtCorrEtaCent[4]={fZvNbins-1,fPtCorrNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Double_t minbinsZvPtCorrEtaCent[4]={-30.,0,-1.5,0};
  Double_t maxbinsZvPtCorrEtaCent[4]={30  ,200,1.5,100};
  
  Int_t binsPhiPtEtaCent[5]={fDeltaphiNbins-1,fPtNbins-1,fEtaNbins-1,fPhiNbins-1,fCentralityNbins-1};
  Double_t minbinsPhiPtEtaCent[5]={0,				0,-1.5,0 ,0};
  Double_t maxbinsPhiPtEtaCent[5]={TMath::Pi()/2.,200,1.5,2.*TMath::Pi(),100};
  
  Int_t binsDeltaphiPtEtaPhiZvCent[6] =       { fDeltaphiNbins-1, fPtNbins-1, fEtaNbins-1, fPhiNbins-1, fZvNbins-1, fCentralityNbins-1 };
  Double_t minbinsDeltaphiPtEtaPhiZvCent[6] = {  0              , 0     , -1.5        , 0            , -10.      ,   0};
  Double_t maxbinsDeltaphiPtEtaPhiZvCent[6] = {  TMath::Pi()/2.  , 200     , 1.5        , 2.*TMath::Pi()   , 10.      ,   100};
  
  Int_t binsPtCorrEtaPhiCorrZvCent[5] =       { fPtCorrNbins-1, fEtaNbins-1, fPhiCorrNbins-1, fZvNbins-1, fCentralityNbins-1 };
  Double_t minbinsPtCorrEtaPhiCorrZvCent[5] = { 0     , -1.5        , 0            , -10.      ,   0};
  Double_t maxbinsPtCorrEtaPhiCorrZvCent[5] = {  200     , 1.5        , 2.*TMath::Pi()   , 10.      ,   100};
  
  Int_t binsZvMultCent[3]={fZvNbins-1,fMultFineNbins-1,fCentralityNbins-1};
  
  Int_t binsPhiPtCorrEtaCent[4]={fPtCorrNbins-1,fEtaNbins-1,fPhiNbins-1,fCentralityNbins-1};
  Double_t minbinsPhiPtCorrEtaCent[4] = {0.0, -1.5, 0, 0};
  Double_t maxbinsPhiPtCorrEtaCent[4] = {200.0, 1.5, 2.*TMath::Pi(), 100};
  
  Int_t binsOneOverPtPtResCent[3]={400,300,11};
  Double_t minbinsOneOverPtPtResCent[3]={0,0,0};
  Double_t maxbinsOneOverPtPtResCent[3]={1,0.015,100};
  
  Int_t binsPtPtResCent[3]={fPtNbins-1,300,11};
  Double_t minbinsPtPtResCent[3]={0,0,0};
  Double_t maxbinsPtPtResCent[3]={200,0.2,100};
  
  // define Histograms
  fZvPtEtaCent = new THnSparseF("fZvPtEtaCent","Zv:Pt:Eta:Centrality",4,binsZvPtEtaCent);
  fZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fZvPtEtaCent->SetBinEdges(1,fBinsPt);
  fZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fZvPtEtaCent->GetAxis(0)->SetTitle("Zv (cm)");
  fZvPtEtaCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fZvPtEtaCent->GetAxis(2)->SetTitle("Eta");
  fZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fZvPtEtaCent->Sumw2();
  
  fDeltaphiPtEtaPhiCent = new THnSparseF("fDeltaphiPtEtaPhiCent","Deltaphi:Pt:Eta:Phi:Centrality",5,binsPhiPtEtaCent);
  fDeltaphiPtEtaPhiCent->SetBinEdges(0,fBinsDeltaphi);
  fDeltaphiPtEtaPhiCent->SetBinEdges(1,fBinsPt);
  fDeltaphiPtEtaPhiCent->SetBinEdges(2,fBinsEta);
  fDeltaphiPtEtaPhiCent->SetBinEdges(3,fBinsPhi);
  fDeltaphiPtEtaPhiCent->SetBinEdges(4,fBinsCentrality);
  fDeltaphiPtEtaPhiCent->GetAxis(0)->SetTitle("#Delta phi to ep");
  fDeltaphiPtEtaPhiCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fDeltaphiPtEtaPhiCent->GetAxis(2)->SetTitle("Eta");
  fDeltaphiPtEtaPhiCent->GetAxis(3)->SetTitle("Phi");
  fDeltaphiPtEtaPhiCent->GetAxis(4)->SetTitle("Centrality");
  fDeltaphiPtEtaPhiCent->Sumw2();
  
  fDeltaphiPtEtaPhiZvCent = new THnSparseF("fDeltaphiPtEtaPhiZvCent","Deltaphi:Pt:Eta:Phi:Zv:Centrality",6,binsDeltaphiPtEtaPhiZvCent);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(0,fBinsDeltaphi);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(1,fBinsPt);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(2,fBinsEta);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(3,fBinsPhi);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(4,fBinsZv);
  fDeltaphiPtEtaPhiZvCent->SetBinEdges(5,fBinsCentrality);
  fDeltaphiPtEtaPhiZvCent->GetAxis(0)->SetTitle("#Delta phi to ep");
  fDeltaphiPtEtaPhiZvCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fDeltaphiPtEtaPhiZvCent->GetAxis(2)->SetTitle("Eta");
  fDeltaphiPtEtaPhiZvCent->GetAxis(3)->SetTitle("Phi");
  fDeltaphiPtEtaPhiZvCent->GetAxis(4)->SetTitle("Zv");
  fDeltaphiPtEtaPhiZvCent->GetAxis(5)->SetTitle("Centrality");
  fDeltaphiPtEtaPhiZvCent->Sumw2();
  
  fPtResptCent = new THnSparseF("fPtResptCent","OneOverPt:PtRes:Centrality",3,binsOneOverPtPtResCent, minbinsOneOverPtPtResCent, maxbinsOneOverPtPtResCent);
  fPtResptCent->SetBinEdges(2, fBinsCentrality);
  fPtResptCent->GetAxis(0)->SetTitle("1/pT (GeV/c)^{-1}");
  fPtResptCent->GetAxis(1)->SetTitle("#sigma(1/pT)");
  fPtResptCent->GetAxis(2)->SetTitle("centrality");
  fPtResptCent->Sumw2();
  
  fPtResptptCent = new THnSparseF("fPtResptptCent","OneOverPt:PtRes*Pt:Centrality",3,binsPtPtResCent, minbinsPtPtResCent, maxbinsPtPtResCent);
  fPtResptptCent->SetBinEdges(0, fBinsPt);
  fPtResptptCent->SetBinEdges(2, fBinsCentrality);
  fPtResptptCent->GetAxis(0)->SetTitle("pT (GeV/c)");
  fPtResptptCent->GetAxis(1)->SetTitle("#sigma(1/pT)*pT");
  fPtResptptCent->GetAxis(2)->SetTitle("centrality");
  fPtResptptCent->Sumw2();
  
  
  fPtEvent = new TH2F("fPtEvent","fPtEvent;pT;eventnumber",fPtNbins-1, fBinsPt,200,0,200);
  fPtEvent->GetXaxis()->SetTitle("Pt (GeV/c)");
  fPtEvent->GetYaxis()->SetTitle("event number");
  fPtEvent->Sumw2();
  
  fMCRecPrimZvPtEtaCent = new THnSparseF("fMCRecPrimZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtCorrEtaCent, minbinsZvPtCorrEtaCent, maxbinsZvPtCorrEtaCent);
  fMCRecPrimZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCRecPrimZvPtEtaCent->SetBinEdges(1,fBinsPtCorr);
  fMCRecPrimZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecPrimZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecPrimZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  fMCRecPrimZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCRecPrimZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCRecPrimZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecPrimZvPtEtaCent->Sumw2();
  
  fMCGenZvPtEtaCent = new THnSparseF("fMCGenZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtCorrEtaCent, minbinsZvPtCorrEtaCent, maxbinsZvPtCorrEtaCent);
  fMCGenZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCGenZvPtEtaCent->SetBinEdges(1,fBinsPtCorr);
  fMCGenZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCGenZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCGenZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  fMCGenZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCGenZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCGenZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCGenZvPtEtaCent->Sumw2();
  
  fMCRecSecZvPtEtaCent = new THnSparseF("fMCRecSecZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtCorrEtaCent, minbinsZvPtCorrEtaCent, maxbinsZvPtCorrEtaCent);
  fMCRecSecZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCRecSecZvPtEtaCent->SetBinEdges(1,fBinsPtCorr);
  fMCRecSecZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecSecZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecSecZvPtEtaCent->GetAxis(0)->SetTitle("MC Sec Zv (cm)");
  fMCRecSecZvPtEtaCent->GetAxis(1)->SetTitle("MC Sec Pt (GeV/c)");
  fMCRecSecZvPtEtaCent->GetAxis(2)->SetTitle("MC Sec Eta");
  fMCRecSecZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecSecZvPtEtaCent->Sumw2();
  
  fMCPtEtaPhiCent = new THnF("fMCPtEtaPhiCent","Pt:Eta:Phi:Centrality",4,binsPhiPtCorrEtaCent, minbinsPhiPtCorrEtaCent, maxbinsPhiPtCorrEtaCent);
  fMCPtEtaPhiCent->SetBinEdges(0,fBinsPtCorr);
  fMCPtEtaPhiCent->SetBinEdges(1,fBinsEta);
  fMCPtEtaPhiCent->SetBinEdges(2,fBinsPhi);
  fMCPtEtaPhiCent->SetBinEdges(3,fBinsCentrality);
  fMCPtEtaPhiCent->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fMCPtEtaPhiCent->GetAxis(1)->SetTitle("Eta");
  fMCPtEtaPhiCent->GetAxis(2)->SetTitle("Phi");
  fMCPtEtaPhiCent->GetAxis(3)->SetTitle("Centrality");
  fMCPtEtaPhiCent->Sumw2();
  
  fMCRecPrimPtEtaPhiCent = new THnF("fMCRecPrimPtEtaPhiCent","mcPt:mcEta:mcPhi:Centrality",4,binsPhiPtCorrEtaCent, minbinsPhiPtCorrEtaCent, maxbinsPhiPtCorrEtaCent);
  fMCRecPrimPtEtaPhiCent->SetBinEdges(0,fBinsPtCorr);
  fMCRecPrimPtEtaPhiCent->SetBinEdges(1,fBinsEta);
  fMCRecPrimPtEtaPhiCent->SetBinEdges(2,fBinsPhi);
  fMCRecPrimPtEtaPhiCent->SetBinEdges(3,fBinsCentrality);
  fMCRecPrimPtEtaPhiCent->GetAxis(0)->SetTitle("MC Pt (GeV/c)");
  fMCRecPrimPtEtaPhiCent->GetAxis(1)->SetTitle("MC Eta");
  fMCRecPrimPtEtaPhiCent->GetAxis(2)->SetTitle("MC Phi");
  fMCRecPrimPtEtaPhiCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecPrimPtEtaPhiCent->Sumw2();
  
  fMCGenPtEtaPhiCent = new THnF("fMCGenPtEtaPhiCent","mcPt:mcEta:mcPhi:Centrality",4,binsPhiPtCorrEtaCent, minbinsPhiPtCorrEtaCent, maxbinsPhiPtCorrEtaCent);
  fMCGenPtEtaPhiCent->SetBinEdges(0,fBinsPtCorr);
  fMCGenPtEtaPhiCent->SetBinEdges(1,fBinsEta);
  fMCGenPtEtaPhiCent->SetBinEdges(2,fBinsPhi);
  fMCGenPtEtaPhiCent->SetBinEdges(3,fBinsCentrality);
  fMCGenPtEtaPhiCent->GetAxis(0)->SetTitle("MC Pt (GeV/c)");
  fMCGenPtEtaPhiCent->GetAxis(1)->SetTitle("MC Eta");
  fMCGenPtEtaPhiCent->GetAxis(2)->SetTitle("MC Phi");
  fMCGenPtEtaPhiCent->GetAxis(3)->SetTitle("Centrality");
  fMCGenPtEtaPhiCent->Sumw2();
  
  fMCRecSecPtEtaPhiCent = new THnF("fMCRecSecPtEtaPhiCent","mcPt:mcEta:mcPhi:Centrality",4,binsPhiPtCorrEtaCent, minbinsPhiPtCorrEtaCent, maxbinsPhiPtCorrEtaCent);
  fMCRecSecPtEtaPhiCent->SetBinEdges(0,fBinsPtCorr);
  fMCRecSecPtEtaPhiCent->SetBinEdges(1,fBinsEta);
  fMCRecSecPtEtaPhiCent->SetBinEdges(2,fBinsPhi);
  fMCRecSecPtEtaPhiCent->SetBinEdges(3,fBinsCentrality);
  fMCRecSecPtEtaPhiCent->GetAxis(0)->SetTitle("MC Sec Pt (GeV/c)");
  fMCRecSecPtEtaPhiCent->GetAxis(1)->SetTitle("MC Sec Eta");
  fMCRecSecPtEtaPhiCent->GetAxis(2)->SetTitle("MC Phi");
  fMCRecSecPtEtaPhiCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecSecPtEtaPhiCent->Sumw2();
  
  //
  // MC Histograms Pt:Eta:Phi:Zv:Centrality
  //
  
  fMCPtEtaPhiZvCent = new THnF("fMCPtEtaPhiZvCent","Pt:Eta:Phi:Zv:Centrality",5,binsPtCorrEtaPhiCorrZvCent, minbinsPtCorrEtaPhiCorrZvCent, maxbinsPtCorrEtaPhiCorrZvCent);
  fMCPtEtaPhiZvCent->SetBinEdges(0,fBinsPtCorr);
  fMCPtEtaPhiZvCent->SetBinEdges(1,fBinsEta);
  fMCPtEtaPhiZvCent->SetBinEdges(2,fBinsPhiCorr);
  fMCPtEtaPhiZvCent->SetBinEdges(3,fBinsZv);
  fMCPtEtaPhiZvCent->SetBinEdges(4,fBinsCentrality);
  fMCPtEtaPhiZvCent->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fMCPtEtaPhiZvCent->GetAxis(1)->SetTitle("Eta");
  fMCPtEtaPhiZvCent->GetAxis(2)->SetTitle("Phi");
  fMCPtEtaPhiZvCent->GetAxis(3)->SetTitle("Zv");
  fMCPtEtaPhiZvCent->GetAxis(4)->SetTitle("Centrality");
  fMCPtEtaPhiZvCent->Sumw2();
  
  fMCRecPrimPtEtaPhiZvCent = new THnF("fMCRecPrimPtEtaPhiZvCent","mcPt:mcEta:mcPhi:Zv:Centrality",5,binsPtCorrEtaPhiCorrZvCent, minbinsPtCorrEtaPhiCorrZvCent, maxbinsPtCorrEtaPhiCorrZvCent);
  fMCRecPrimPtEtaPhiZvCent->SetBinEdges(0,fBinsPtCorr);
  fMCRecPrimPtEtaPhiZvCent->SetBinEdges(1,fBinsEta);
  fMCRecPrimPtEtaPhiZvCent->SetBinEdges(2,fBinsPhiCorr);
  fMCRecPrimPtEtaPhiZvCent->SetBinEdges(3,fBinsZv);
  fMCRecPrimPtEtaPhiZvCent->SetBinEdges(4,fBinsCentrality);
  fMCRecPrimPtEtaPhiZvCent->GetAxis(0)->SetTitle("MC Pt (GeV/c)");
  fMCRecPrimPtEtaPhiZvCent->GetAxis(1)->SetTitle("MC Eta");
  fMCRecPrimPtEtaPhiZvCent->GetAxis(2)->SetTitle("MC Phi");
  fMCRecPrimPtEtaPhiZvCent->GetAxis(3)->SetTitle("Zv");
  fMCRecPrimPtEtaPhiZvCent->GetAxis(4)->SetTitle("Centrality");
  fMCRecPrimPtEtaPhiZvCent->Sumw2();
  
  fMCGenPtEtaPhiZvCent = new THnF("fMCGenPtEtaPhiZvCent","mcPt:mcEta:mcPhi:Zv:Centrality",5,binsPtCorrEtaPhiCorrZvCent, minbinsPtCorrEtaPhiCorrZvCent, maxbinsPtCorrEtaPhiCorrZvCent);
  fMCGenPtEtaPhiZvCent->SetBinEdges(0,fBinsPtCorr);
  fMCGenPtEtaPhiZvCent->SetBinEdges(1,fBinsEta);
  fMCGenPtEtaPhiZvCent->SetBinEdges(2,fBinsPhiCorr);
  fMCGenPtEtaPhiZvCent->SetBinEdges(3,fBinsZv);
  fMCGenPtEtaPhiZvCent->SetBinEdges(4,fBinsCentrality);
  fMCGenPtEtaPhiZvCent->GetAxis(0)->SetTitle("MC Pt (GeV/c)");
  fMCGenPtEtaPhiZvCent->GetAxis(1)->SetTitle("MC Eta");
  fMCGenPtEtaPhiZvCent->GetAxis(2)->SetTitle("MC Phi");
  fMCGenPtEtaPhiZvCent->GetAxis(3)->SetTitle("Zv");
  fMCGenPtEtaPhiZvCent->GetAxis(4)->SetTitle("Centrality");
  fMCGenPtEtaPhiZvCent->Sumw2();
  
  fMCRecSecPtEtaPhiZvCent = new THnF("fMCRecSecPtEtaPhiZvCent","mcPt:mcEta:mcPhi:Zv:Centrality",5,binsPtCorrEtaPhiCorrZvCent, minbinsPtCorrEtaPhiCorrZvCent, maxbinsPtCorrEtaPhiCorrZvCent);
  fMCRecSecPtEtaPhiZvCent->SetBinEdges(0,fBinsPtCorr);
  fMCRecSecPtEtaPhiZvCent->SetBinEdges(1,fBinsEta);
  fMCRecSecPtEtaPhiZvCent->SetBinEdges(2,fBinsPhiCorr);
  fMCRecSecPtEtaPhiZvCent->SetBinEdges(3,fBinsZv);
  fMCRecSecPtEtaPhiZvCent->SetBinEdges(4,fBinsCentrality);
  fMCRecSecPtEtaPhiZvCent->GetAxis(0)->SetTitle("MC Sec Pt (GeV/c)");
  fMCRecSecPtEtaPhiZvCent->GetAxis(1)->SetTitle("MC Sec Eta");
  fMCRecSecPtEtaPhiZvCent->GetAxis(2)->SetTitle("MC Phi");
  fMCRecSecPtEtaPhiZvCent->GetAxis(3)->SetTitle("Zv");
  fMCRecSecPtEtaPhiZvCent->GetAxis(4)->SetTitle("Centrality");
  fMCRecSecPtEtaPhiZvCent->Sumw2();
  
  fPt = new TH1F("fPt","fPt",2000,0,200);
  fPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fPt->GetYaxis()->SetTitle("dN/dp_{T}");
  fPt->Sumw2();
  
  fMCPt = new TH1F("fMCPt","fMCPt",2000,0,200);
  fMCPt->GetXaxis()->SetTitle("MC p_{T} (GeV/c)");
  fMCPt->GetYaxis()->SetTitle("dN/dp_{T}");
  fMCPt->Sumw2();
  
  fEventStatistics = new TH1F("fEventStatistics","fEventStatistics",10,0,10);
  fEventStatistics->GetYaxis()->SetTitle("number of events");
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fEventStatistics->SetBit(TH1::kCanRebin);
#endif
  
  fEventStatisticsCentrality = new TH1F("fEventStatisticsCentrality","fEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  fMCEventStatisticsCentrality = new TH1F("fMCEventStatisticsCentrality","fMCEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fMCEventStatisticsCentrality->GetYaxis()->SetTitle("number of MC events");
  
  fAllEventStatisticsCentrality = new TH1F("fAllEventStatisticsCentrality","fAllEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fAllEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  fEventStatisticsCentralityTrigger = new TH2F("fEventStatisticsCentralityTrigger","fEventStatisticsCentralityTrigger;centrality;trigger",100,0,100,3,0,3);
  fEventStatisticsCentralityTrigger->Sumw2();
  
  fZvMultCent = new THnSparseF("fZvMultCent","Zv:mult:Centrality",3,binsZvMultCent);
  fZvMultCent->SetBinEdges(0,fBinsZv);
  fZvMultCent->SetBinEdges(1,fBinsMultFine);
  fZvMultCent->SetBinEdges(2,fBinsCentrality);
  fZvMultCent->GetAxis(0)->SetTitle("Zv (cm)");
  fZvMultCent->GetAxis(1)->SetTitle("N_{acc}");
  fZvMultCent->GetAxis(2)->SetTitle("Centrality");
  fZvMultCent->Sumw2();
  
  fTriggerStatistics = new TH1F("fTriggerStatistics","fTriggerStatistics",10,0,10);
  fTriggerStatistics->GetYaxis()->SetTitle("number of events");
  
  fCharge = new TH1F("fCharge","fCharge",30, -5, 5);
  fCharge->GetXaxis()->SetTitle("Charge");
  fCharge->GetYaxis()->SetTitle("number of tracks");
  
  fMCCharge = new TH1F("fMCCharge","fMCCharge",30, -5, 5);
  fMCCharge->GetXaxis()->SetTitle("MC Charge");
  fMCCharge->GetYaxis()->SetTitle("number of tracks");
  
  Int_t binsDCAxyDCAzPtEtaPhi[6] =   { 10 , 10 , fPtCheckNbins-1, fEtaCheckNbins-1,             18, fCentralityNbins-1 };
  Double_t minDCAxyDCAzPtEtaPhi[6] = { -5 , -5 ,               0,             -1.5,             0.,                  0 };
  Double_t maxDCAxyDCAzPtEtaPhi[6] = {  5.,  5.,             100,              1.5, 2.*TMath::Pi(),                100 };
  
  fDCAPtAll = new THnSparseF("fDCAPtAll","fDCAPtAll",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fDCAPtAccepted = new THnSparseF("fDCAPtAccepted","fDCAPtAccepted",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fMCDCAPtSecondary = new THnSparseF("fMCDCAPtSecondary","fMCDCAPtSecondary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fMCDCAPtPrimary = new THnSparseF("fMCDCAPtPrimary","fMCDCAPtPrimary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  
  fDCAPtAll->SetBinEdges(2, fBinsPtCheck);
  fDCAPtAccepted->SetBinEdges(2, fBinsPtCheck);
  fMCDCAPtSecondary->SetBinEdges(2, fBinsPtCheck);
  fMCDCAPtPrimary->SetBinEdges(2, fBinsPtCheck);
  
  fDCAPtAll->SetBinEdges(3, fBinsEtaCheck);
  fDCAPtAccepted->SetBinEdges(3, fBinsEtaCheck);
  fMCDCAPtSecondary->SetBinEdges(3, fBinsEtaCheck);
  fMCDCAPtPrimary->SetBinEdges(3, fBinsEtaCheck);
  
  fDCAPtAll->SetBinEdges(5, fBinsCentrality);
  fDCAPtAccepted->SetBinEdges(5, fBinsCentrality);
  fMCDCAPtSecondary->SetBinEdges(5, fBinsCentrality);
  fMCDCAPtPrimary->SetBinEdges(5, fBinsCentrality);
  
  fDCAPtAll->Sumw2();
  fDCAPtAccepted->Sumw2();
  fMCDCAPtSecondary->Sumw2();
  fMCDCAPtPrimary->Sumw2();
  
  fDCAPtAll->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fDCAPtAll->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fDCAPtAll->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fDCAPtAll->GetAxis(3)->SetTitle("#eta");
  fDCAPtAll->GetAxis(4)->SetTitle("#phi");
  fDCAPtAll->GetAxis(5)->SetTitle("Centrality");
  
  fDCAPtAccepted->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fDCAPtAccepted->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fDCAPtAccepted->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fDCAPtAccepted->GetAxis(3)->SetTitle("#eta");
  fDCAPtAccepted->GetAxis(4)->SetTitle("#phi");
  fDCAPtAccepted->GetAxis(5)->SetTitle("Centrality");
  
  fMCDCAPtSecondary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fMCDCAPtSecondary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fMCDCAPtSecondary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fMCDCAPtSecondary->GetAxis(3)->SetTitle("#eta");
  fMCDCAPtSecondary->GetAxis(4)->SetTitle("#phi");
  fMCDCAPtSecondary->GetAxis(5)->SetTitle("Centrality");
  
  fMCDCAPtPrimary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fMCDCAPtPrimary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fMCDCAPtPrimary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fMCDCAPtPrimary->GetAxis(3)->SetTitle("#eta");
  fMCDCAPtPrimary->GetAxis(4)->SetTitle("#phi");
  fMCDCAPtPrimary->GetAxis(5)->SetTitle("Centrality");
  
  
  char cFullTempTitle[255];
  char cTempTitleAxis0All[255];
  char cTempTitleAxis0Acc[255];
  //   char cTempTitleAxis1[255];
  char cFullTempName[255];
  char cTempNameAxis0[255];
  //   char cTempNameAxis1[255];
  const Int_t iNbinRowsClusters = 21;
  //   Double_t dBinsRowsClusters[iNbinRowsClusters] = {0, 7.95, 15.9, 23.85, 31.8, 39.75, 47.7, 55.65, 63.6, 71.55, 79.5, 87.45, 95.4, 103.35, 111.3, 119.25, 127.2, 135.15, 143.1, 151.05, 159.};
  
  const Int_t iNbinChi = 51;
  const Int_t iNbinLength = 165;
  const Int_t iNbinRowsOverClusters = 60;
  //   Double_t dBinsChi[iNbinChi] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.6, 6.8, 7, 7.2, 7.4, 7.6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6, 9.8,10.};
  
  Int_t iNbin = 0;
  //   Double_t *dBins = 0x0;
  Double_t dBinMin = 0;
  Double_t dBinMax = 0;
  
  if(AreCrossCheckCorrelationHistosEnabled())
  {
    for(Int_t iCheckQuant = 0; iCheckQuant < cqMax; iCheckQuant++)
    {
      // iCheckQuant: 0 = CrossedRows, 1 = Nclusters, 2 = Chi^2/clusterTPC
      if(iCheckQuant == cqCrossedRows)
      {
        snprintf(cTempTitleAxis0All,255, "NcrossedRows before Cut");
        snprintf(cTempTitleAxis0Acc,255, "NcrossedRows after Cut");
        snprintf(cTempNameAxis0,255, "CrossedRows");
        iNbin = iNbinRowsClusters;
        dBinMin = 0;
        dBinMax = 159.;
      }
      else if(iCheckQuant == cqNcluster)
      {
        snprintf(cTempTitleAxis0All,255, "Nclusters before Cut");
        snprintf(cTempTitleAxis0Acc,255, "Nclusters after Cut");
        snprintf(cTempNameAxis0,255, "Clusters");
        iNbin = iNbinRowsClusters;
        dBinMin = 0;
        dBinMax = 159.;
      }
      else if(iCheckQuant == cqChi)
      {
        snprintf(cTempTitleAxis0All,255, "#Chi^{2}/cluster before Cut");
        snprintf(cTempTitleAxis0Acc,255, "#Chi^{2}/cluster after Cut");
        snprintf(cTempNameAxis0,255, "Chi");
        iNbin = iNbinChi;
        dBinMin = 0;
        dBinMax = 10.;
      }
      else if(iCheckQuant == cqLength)
      {
        snprintf(cTempTitleAxis0All,255, "Length in TPC before Cut (cm)");
        snprintf(cTempTitleAxis0Acc,255, "Length in TPC after Cut (cm)");
        snprintf(cTempNameAxis0,255, "Length");
        iNbin = iNbinLength;
        dBinMin = 0;
        dBinMax = 165.;
      }
      else if(iCheckQuant == cqRowsOverFindable)
      {
        snprintf(cTempTitleAxis0All,255, "Number of Crossed Rows / Number of Findable Clusters before Cut");
        snprintf(cTempTitleAxis0Acc,255, "Number of Crossed Rows / Number of Findable Clusters before Cut");
        snprintf(cTempNameAxis0,255, "RowsOverFindable");
        iNbin = iNbinRowsOverClusters;
        dBinMin = 0.6;
        dBinMax = 1.2;
      }
      
      
      Int_t binsCheckPtEtaPhi[5] = { iNbin, fPtCheckNbins-1, fEtaCheckNbins-1, 18, fCentralityNbins-1};
      //     Int_t binsCheckPtEtaPhi[5] = { iNbin, fPtNbins-1, fEtaCheckNbins-1, 18, fCentralityNbins-1};
      Double_t minCheckPtEtaPhi[5] = { dBinMin,  0, -1.5, 0., 0, };
      Double_t maxCheckPtEtaPhi[5] = { dBinMax, 100, 1.5, 2.*TMath::Pi(), 100};
      
      snprintf(cFullTempName, 255, "f%sPtEtaPhiAll",cTempNameAxis0);
      snprintf(cFullTempTitle, 255,"%s;%s;p_{T} (GeV/c);#eta;#phi;Centrality", cFullTempName, cTempTitleAxis0All);
      fCrossCheckAll[iCheckQuant] = new THnF(cFullTempName, cFullTempTitle, 5, binsCheckPtEtaPhi, minCheckPtEtaPhi, maxCheckPtEtaPhi);
      fCrossCheckAll[iCheckQuant]->SetBinEdges(1, fBinsPtCheck);
      fCrossCheckAll[iCheckQuant]->SetBinEdges(2, fBinsEtaCheck);
      fCrossCheckAll[iCheckQuant]->Sumw2();
      
      snprintf(cFullTempName, 255, "f%sPtEtaPhiAcc",cTempNameAxis0);
      snprintf(cFullTempTitle, 255,"%s;%s;p_{T} (GeV/c);#eta;#phi;Centrality", cFullTempName, cTempTitleAxis0Acc);
      fCrossCheckAcc[iCheckQuant] = new THnF(cFullTempName, cFullTempTitle, 5, binsCheckPtEtaPhi, minCheckPtEtaPhi, maxCheckPtEtaPhi);
      fCrossCheckAcc[iCheckQuant]->SetBinEdges(1, fBinsPtCheck);
      fCrossCheckAcc[iCheckQuant]->SetBinEdges(2, fBinsEtaCheck);
      fCrossCheckAcc[iCheckQuant]->Sumw2();
    } // end iCheckQuant
  } // AreCrossCheckCorrelationHistosEnabled
  
  fCutPercClusters = new TH1F("fCutPercClusters","fCutPercClusters;NclustersTPC;counts",160,0,160);
  fCutPercClusters->Sumw2();
  fCutPercCrossed = new TH1F("fCutPercCrossed","fCutPercCrossed;NcrossedRowsTPC;counts",160,0,160);
  fCutPercCrossed->Sumw2();
  
  fCrossCheckRowsLength = new TH3F("fCrossCheckRowsLength","fCrossCheckRowsLength;Length in TPC;NcrossedRows;Centrality",170,0,170,170,0,170, 100,0,100);
  fCrossCheckRowsLength->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckRowsLength->Sumw2();
  
  fCrossCheckClusterLength = new TH3F("fCrossCheckClusterLength","fCrossCheckClusterLength;Length in TPC;NClusters;Centrality",170,0,170,170,0,170, 100,0,100);
  fCrossCheckClusterLength->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckClusterLength->Sumw2();
  
  fCrossCheckRowsLengthAcc = new TH3F("fCrossCheckRowsLengthAcc","fCrossCheckRowsLengthAcc;Length in TPC;NcrossedRows;Centrality",170,0,170,170,0,170, 100,0,100);
  fCrossCheckRowsLengthAcc->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckRowsLengthAcc->Sumw2();
  
  fCrossCheckClusterLengthAcc = new TH3F("fCrossCheckClusterLengthAcc","fCrossCheckClusterLengthAcc;Length in TPC;NClusters;Centrality",170,0,170,170,0,170, 100,0,100);
  fCrossCheckClusterLengthAcc->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckClusterLengthAcc->Sumw2();
  
  fCrossCheckPtresLength = new TH3F("fCrossCheckPtresLength","fCrossCheckPtresLength;Length in TPC;#sigma(1/pT)*pT;Centrality",170,0,170,100,0,1, 100,0,100);
  fCrossCheckPtresLength->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckPtresLength->Sumw2();
  
  fCrossCheckPtresRows = new TH3F("fCrossCheckPtresRows","fCrossCheckPtresRows;NcrossedRows;#sigma(1/pT)*pT;Centrality",170,0,170,100,0,1, 100,0,100);
  fCrossCheckPtresRows->GetZaxis()->Set( fCentralityNbins-1, fBinsCentrality);
  fCrossCheckPtresRows->Sumw2();
  
  fCutSettings = new TH1F("fCutSettings","fCutSettings",100,0,10);
  fCutSettings->GetYaxis()->SetTitle("cut value");
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fCutSettings->SetBit(TH1::kCanRebin);
#endif
  
  fEventplaneDist = new TH1F("fEventplaneDist","fEventplaneDist",200, 0, 2.*TMath::Pi());
  fEventplaneDist->GetXaxis()->SetTitle("#phi (event plane)");
  fEventplaneDist->Sumw2();
  
  fEventplaneRunDist = new TH2F("fEventplaneRunDist","fEventplaneRunDist",200, 0, 2.*TMath::Pi(),fRunNumberNbins-1, fBinsRunNumber );
  fEventplaneRunDist->GetXaxis()->SetTitle("#phi (event plane)");
  fEventplaneRunDist->GetYaxis()->SetTitle("runnumber");
  fEventplaneRunDist->Sumw2();
  
  fMCEventplaneDist = new TH1F("fMCEventplaneDist","fMCEventplaneDist",20, -1.*TMath::Pi(), TMath::Pi());
  fMCEventplaneDist->GetXaxis()->SetTitle("#phi (MC event plane)");
  fMCEventplaneDist->Sumw2();
  
  fCorrelEventplaneMCDATA = new TH2F("fCorrelEventplaneMCDATA","fCorrelEventplaneMCDATA",40, -2.*TMath::Pi(), 2.*TMath::Pi(), 40, -2.*TMath::Pi(), 2.*TMath::Pi());
  fCorrelEventplaneMCDATA->GetXaxis()->SetTitle("#phi (event plane)");
  fCorrelEventplaneMCDATA->GetYaxis()->SetTitle("#phi (MC event plane)");
  fCorrelEventplaneMCDATA->Sumw2();
  
  Int_t binsCorrelPhiPhiCent[3] = { 40, 40, 10};
  Double_t minCorrelPhiPhiCent[3] = { -2.*TMath::Pi(), -2.*TMath::Pi(), 0};
  Double_t maxCorrelPhiPhiCent[3] = { 2.*TMath::Pi(), 2.*TMath::Pi(), 100};
  
  fCorrelEventplaneDefaultCorrected = new THnSparseF("fCorrelEventplaneDefaultCorrected","fCorrelEventplaneDefaultCorrected",3,binsCorrelPhiPhiCent, minCorrelPhiPhiCent, maxCorrelPhiPhiCent);
  fCorrelEventplaneDefaultCorrected->SetBinEdges(2, fBinsCentrality);
  fCorrelEventplaneDefaultCorrected->GetAxis(0)->SetTitle("#phi (event plane)");
  fCorrelEventplaneDefaultCorrected->GetAxis(1)->SetTitle("#phi (corrected event plane)");
  fCorrelEventplaneDefaultCorrected->GetAxis(2)->SetTitle("centrality");
  fCorrelEventplaneDefaultCorrected->Sumw2();
  
  fEventplaneSubtractedPercentage = new TH2F("fEventplaneSubtractedPercentage","fEventplaneSubtractedPercentage",100, 0,1, fCentralityNbins-1, fBinsCentrality);
  fEventplaneSubtractedPercentage->GetXaxis()->SetTitle("percentage of tracks, which have been subtracted during analysis");
  fEventplaneSubtractedPercentage->GetYaxis()->SetTitle("centrality");
  fEventplaneSubtractedPercentage->Sumw2();
  
  fChargeOverPtRuns = new TH2F("fChargeOverPtRuns","fChargeOverPtRuns",2000, -10, 10, fRunNumberNbins-1, fBinsRunNumber );
  fChargeOverPtRuns->GetXaxis()->SetTitle("q/#it{p}_{T}");
  fChargeOverPtRuns->GetYaxis()->SetTitle("runnumber");
  char cChargeOverPtTitle[255];
  snprintf(cChargeOverPtTitle, 255, "%.2f < pT < %.2f, %.2f < #eta < %.2f",GetCutPtMin(), GetCutPtMax(), GetCutEtaMin(), GetCutEtaMax());
  fChargeOverPtRuns->SetTitle(cChargeOverPtTitle);
  fChargeOverPtRuns->Sumw2();
  
  fVZEROMultCentrality = new TH2F("fVZEROMultCentrality","fVZEROMultCentrality",3000,0,30000, fCentralityNbins-1, fBinsCentrality);
  fVZEROMultCentrality->GetXaxis()->SetTitle("VZERO Multiplicity");
  fVZEROMultCentrality->GetYaxis()->SetTitle("centrality");
  fVZEROMultCentrality->Sumw2();
  
  fVEROMultRefMult = new TH2F("fVEROMultRefMult","fVEROMultRefMult",3000,0,30000,3000,0,30000);
  fVEROMultRefMult->GetXaxis()->SetTitle("VZERO Multiplicity");
  fVEROMultRefMult->GetYaxis()->SetTitle("reference multiplicity");
  fVEROMultRefMult->Sumw2();
  
  // cross check for event plane resolution
  fEPDistCent = new TH2F("fEPDistCent","fEPDistCent",20, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fEPDistCent->GetXaxis()->SetTitle("#phi (#Psi_{EP})");
  fEPDistCent->GetYaxis()->SetTitle("Centrality");
  fEPDistCent->Sumw2();
  
  fPhiCent = new TH2F("fPhiCent","fPhiCent",200, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fPhiCent->GetXaxis()->SetTitle("#phi");
  fPhiCent->GetYaxis()->SetTitle("Centrality");
  fPhiCent->Sumw2();
  
  fPcosEPCent = new TProfile("fPcosEPCent","fPcosEPCent", 100,0,100);
  fPcosEPCent->GetXaxis()->SetTitle("Centrality");
  fPcosEPCent->GetYaxis()->SetTitle("#LT cos 2 #Psi_{EP} #GT");
  fPcosEPCent->Sumw2();
  
  fPsinEPCent = new TProfile("fPsinEPCent","fPsinEPCent", 100,0,100);
  fPsinEPCent->GetXaxis()->SetTitle("Centrality");
  fPsinEPCent->GetYaxis()->SetTitle("#LT sin 2 #Psi_{EP} #GT");
  fPsinEPCent->Sumw2();
  
  fPcosPhiCent = new TProfile("fPcosPhiCent","fPcosPhiCent", 100,0,100);
  fPcosPhiCent->GetXaxis()->SetTitle("Centrality");
  fPcosPhiCent->GetYaxis()->SetTitle("#LT cos 2 #phi #GT");
  fPcosPhiCent->Sumw2();
  
  fPsinPhiCent = new TProfile("fPsinPhiCent","fPsinPhiCent", 100,0,100);
  fPsinPhiCent->GetXaxis()->SetTitle("Centrality");
  fPsinPhiCent->GetYaxis()->SetTitle("#LT sin 2 #phi #GT");
  fPsinPhiCent->Sumw2();
  
  fEPContributionDifference = new TH2F("fEPContributionDifference","fEPContributionDifference",200,-1,1,200,-1,1);
  fEPContributionDifference->GetXaxis()->SetTitle("difference Qx (own-lookup)/own");
  fEPContributionDifference->GetYaxis()->SetTitle("difference Qy (own-lookup)/own");
  fEPContributionDifference->Sumw2();
  
  fDeltaPhiCent = new TH2F("fDeltaPhiCent","fDeltaPhiCent",200, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fDeltaPhiCent->GetXaxis()->SetTitle("#Delta #phi");
  fDeltaPhiCent->GetYaxis()->SetTitle("Centrality");
  fDeltaPhiCent->Sumw2();
  
  fDeltaPhiSymCent = new TH2F("fDeltaPhiSymCent","fDeltaPhiSymCent",200, 0., 0.5*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fDeltaPhiSymCent->GetXaxis()->SetTitle("#Delta #phi");
  fDeltaPhiSymCent->GetYaxis()->SetTitle("Centrality");
  fDeltaPhiSymCent->Sumw2();
  
  fMCRecTracksMult = new TH1F("fMCRecTracksMult","fMCRecTracksMult",4000,-0.5,3999.5);
  fMCRecTracksMult->GetYaxis()->SetTitle("#recoTracks (MC)");
  fMCRecTracksMult->GetXaxis()->SetTitle("reference multiplicity");
  fMCRecTracksMult->Sumw2();
  
  fMCGenTracksMult = new TH1F("fMCGenTracksMult","fMCGenTracksMult",4000,-0.5,3999.5);
  fMCGenTracksMult->GetYaxis()->SetTitle("#genTracks (MC)");
  fMCGenTracksMult->GetXaxis()->SetTitle("reference multiplicity");
  fMCGenTracksMult->Sumw2();
  
  Int_t binsFilterBitPhiCent[3]={3,200,fCentralityNbins-1};
  Double_t minbinsFilterBitPhiCent[3]={0,-2.*TMath::Pi(),0};
  Double_t maxbinsFilterBitPhiCent[3]={3,2.*TMath::Pi(),100};
  
  fCrossCheckFilterBitPhiCent = new THnSparseF("fCrossCheckFilterBitPhiCent","fCrossCheckFilterBitPhiCent",3, binsFilterBitPhiCent, minbinsFilterBitPhiCent, maxbinsFilterBitPhiCent);
  fCrossCheckFilterBitPhiCent->SetBinEdges(2,fBinsCentrality);
  fCrossCheckFilterBitPhiCent->GetAxis(0)->SetTitle("FilterBit");
  fCrossCheckFilterBitPhiCent->GetAxis(1)->SetTitle("#phi");
  fCrossCheckFilterBitPhiCent->GetAxis(2)->SetTitle("Centrality");
  fCrossCheckFilterBitPhiCent->Sumw2();
  
  fTriggerStringsFired = new TH1F("fTriggerStringsFired","fTriggerStringsFired",15,0,15);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fTriggerStringsFired->SetBit(TH1::kCanRebin);
#endif
  fTriggerStringsFired->GetYaxis()->SetTitle("number of fired triggers");
  fTriggerStringsFired->Sumw2();
  
  fTriggerStringComplete = new TH1F("fTriggerStringComplete","fTriggerStringComplete",15,0,15);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fTriggerStringComplete->SetBit(TH1::kCanRebin);
#endif
  fTriggerStringComplete->GetYaxis()->SetTitle("number of events");
  fTriggerStringComplete->Sumw2();
  
  fVertexZ = new TH1F("fVertexZ","fVertexZ",600, -30, 30);
  fVertexZ->GetXaxis()->SetTitle("vtx_{z} (cm)");
  fVertexZ->GetYaxis()->SetTitle("number of events");
  fVertexZ->Sumw2();
  
  fVertexZSPD = new TH1F("fVertexZSPD","fVertexZSPD",600, -30, 30);
  fVertexZSPD->GetXaxis()->SetTitle("vtxSPD_{z} (cm)");
  fVertexZSPD->GetYaxis()->SetTitle("number of events");
  fVertexZSPD->Sumw2();
  
  fVertexZTPC = new TH1F("fVertexZTPC","fVertexZTPC",600, -30, 30);
  fVertexZTPC->GetXaxis()->SetTitle("vtxTPC_{z} (cm)");
  fVertexZTPC->GetYaxis()->SetTitle("number of events");
  fVertexZTPC->Sumw2();
  
  fDeltaVertexZGlobalSPD = new TH1F("fDeltaVertexZGlobalSPD","fDeltaVertexZGlobalSPD",6000, -30, 30);
  fDeltaVertexZGlobalSPD->GetXaxis()->SetTitle("vtx_{z} - vtxSPD_{z} (cm)");
  fDeltaVertexZGlobalSPD->GetYaxis()->SetTitle("number of events");
  fDeltaVertexZGlobalSPD->Sumw2();
  
  fDeltaVertexZGlobalTPC = new TH1F("fDeltaVertexZGlobalTPC","fDeltaVertexZGlobalTPC",6000, -30, 30);
  fDeltaVertexZGlobalTPC->GetXaxis()->SetTitle("vtx_{z} - vtxTPC_{z} (cm)");
  fDeltaVertexZGlobalTPC->GetYaxis()->SetTitle("number of events");
  fDeltaVertexZGlobalTPC->Sumw2();
  
  fVertexContributors = new TH1F("fVertexContributors","fVertexContributors",3000, 0 ,3000);
  fVertexContributors->GetXaxis()->SetTitle("contributors to primary vertex");
  fVertexContributors->GetYaxis()->SetTitle("number of events");
  fVertexContributors->Sumw2();
  
  fVertexZAfterCuts = new TH1F("fVertexZAfterCuts","fVertexZAfterCuts",600, -30, 30);
  fVertexZAfterCuts->GetXaxis()->SetTitle("vtx_{z} (cm) (with vertex cuts)");
  fVertexZAfterCuts->GetYaxis()->SetTitle("number of events");
  fVertexZAfterCuts->Sumw2();
  
  fVertexZSPDAfterCuts = new TH1F("fVertexZSPD","fVertexZSPD",600, -30, 30);
  fVertexZSPD->GetXaxis()->SetTitle("vtxSPD_{z} (cm) (with vertex cuts)");
  fVertexZSPD->GetYaxis()->SetTitle("number of events");
  fVertexZSPD->Sumw2();
  
  fVertexZTPCAfterCuts = new TH1F("fVertexZTPCAfterCuts","fVertexZTPCAfterCuts",600, -30, 30);
  fVertexZTPCAfterCuts->GetXaxis()->SetTitle("vtxTPC_{z} (cm) (with vertex cuts)");
  fVertexZTPCAfterCuts->GetYaxis()->SetTitle("number of events");
  fVertexZTPCAfterCuts->Sumw2();
  
  fDeltaVertexZGlobalSPDAfterCuts = new TH1F("fDeltaVertexZGlobalSPDAfterCuts","fDeltaVertexZGlobalSPDAfterCuts",6000, -30, 30);
  fDeltaVertexZGlobalSPDAfterCuts->GetXaxis()->SetTitle("vtx_{z} - vtxSPD_{z} (cm) (with vertex cuts)");
  fDeltaVertexZGlobalSPDAfterCuts->GetYaxis()->SetTitle("number of events");
  fDeltaVertexZGlobalSPDAfterCuts->Sumw2();
  
  fDeltaVertexZGlobalTPCAfterCuts = new TH1F("fDeltaVertexZGlobalTPCAfterCuts","fDeltaVertexZGlobalTPCAfterCuts",6000, -30, 30);
  fDeltaVertexZGlobalTPCAfterCuts->GetXaxis()->SetTitle("vtx_{z} - vtxTPC_{z} (cm) (with vertex cuts)");
  fDeltaVertexZGlobalTPCAfterCuts->GetYaxis()->SetTitle("number of events");
  fDeltaVertexZGlobalTPCAfterCuts->Sumw2();
  
  fVertexContributorsAfterCuts = new TH1F("fVertexContributorsAfterCuts","fVertexContributorsAfterCuts",3000, 0 ,3000);
  fVertexContributorsAfterCuts->GetXaxis()->SetTitle("contributors to primary vertex (with vertex cuts)");
  fVertexContributorsAfterCuts->GetYaxis()->SetTitle("number of events");
  fVertexContributorsAfterCuts->Sumw2();
  fVertexContributorsAfterCuts->SetMarkerStyle(21);
  fVertexContributorsAfterCuts->SetMarkerColor(kBlack);
  
  fVertexContributorsAfterCutsCent = new TH1F("fVertexContributorsAfterCutsCent","fVertexContributorsAfterCutsCent",3000, 0 ,3000);
  fVertexContributorsAfterCutsCent->GetXaxis()->SetTitle("contributors to primary vertex (with vertex cuts)");
  fVertexContributorsAfterCutsCent->GetYaxis()->SetTitle("number of events with central trigger");
  fVertexContributorsAfterCutsCent->Sumw2();
  fVertexContributorsAfterCutsCent->SetMarkerStyle(20);
  fVertexContributorsAfterCutsCent->SetMarkerColor(kBlue);
  
  fVertexContributorsAfterCutsSemi = new TH1F("fVertexContributorsAfterCutsSemi","fVertexContributorsAfterCutsSemi",3000, 0 ,3000);
  fVertexContributorsAfterCutsSemi->GetXaxis()->SetTitle("contributors to primary vertex (with vertex cuts)");
  fVertexContributorsAfterCutsSemi->GetYaxis()->SetTitle("number of events with semicentral trigger");
  fVertexContributorsAfterCutsSemi->Sumw2();
  fVertexContributorsAfterCutsSemi->SetMarkerStyle(20);
  fVertexContributorsAfterCutsSemi->SetMarkerColor(kRed);
  
  fVertexContributorsAfterCutsMB = new TH1F("fVertexContributorsAfterCutsMB","fVertexContributorsAfterCutsMB",3000, 0 ,3000);
  fVertexContributorsAfterCutsMB->GetXaxis()->SetTitle("contributors to primary vertex (with vertex cuts)");
  fVertexContributorsAfterCutsMB->GetYaxis()->SetTitle("number of events with MB trigger");
  fVertexContributorsAfterCutsMB->Sumw2();
  fVertexContributorsAfterCutsMB->SetMarkerStyle(20);
  fVertexContributorsAfterCutsMB->SetMarkerColor(kGreen+1);
  
  // Add Histos, Profiles etc to List
  fOutputList->Add(fZvPtEtaCent);
  fOutputList->Add(fDeltaphiPtEtaPhiCent);
  fOutputList->Add(fDeltaphiPtEtaPhiZvCent);
  fOutputList->Add(fPtResptCent);
  fOutputList->Add(fPtResptptCent);
  fOutputList->Add(fPtEvent);
  fOutputList->Add(fPt);
  fOutputList->Add(fMCRecPrimZvPtEtaCent);
  fOutputList->Add(fMCGenZvPtEtaCent);
  fOutputList->Add(fMCRecSecZvPtEtaCent);
  fOutputList->Add(fMCPtEtaPhiCent);
  fOutputList->Add(fMCRecPrimPtEtaPhiCent);
  fOutputList->Add(fMCGenPtEtaPhiCent);
  fOutputList->Add(fMCRecSecPtEtaPhiCent);
  
  fOutputList->Add(fMCPtEtaPhiZvCent);
  fOutputList->Add(fMCRecPrimPtEtaPhiZvCent);
  fOutputList->Add(fMCGenPtEtaPhiZvCent);
  fOutputList->Add(fMCRecSecPtEtaPhiZvCent);
  
  fOutputList->Add(fMCPt);
  fOutputList->Add(fEventStatistics);
  fOutputList->Add(fEventStatisticsCentrality);
  fOutputList->Add(fMCEventStatisticsCentrality);
  fOutputList->Add(fAllEventStatisticsCentrality);
  fOutputList->Add(fEventStatisticsCentralityTrigger);
  fOutputList->Add(fZvMultCent);
  fOutputList->Add(fTriggerStatistics);
  fOutputList->Add(fCharge);
  fOutputList->Add(fMCCharge);
  fOutputList->Add(fDCAPtAll);
  fOutputList->Add(fDCAPtAccepted);
  fOutputList->Add(fMCDCAPtSecondary);
  fOutputList->Add(fMCDCAPtPrimary);
  if(AreCrossCheckCorrelationHistosEnabled())
  {
    for(Int_t i = 0; i < cqMax; i++)
    {
      fOutputList->Add(fCrossCheckAll[i]);
      fOutputList->Add(fCrossCheckAcc[i]);
    }
  }
  fOutputList->Add(fCutPercClusters);
  fOutputList->Add(fCutPercCrossed);
  fOutputList->Add(fCrossCheckRowsLength);
  fOutputList->Add(fCrossCheckClusterLength);
  fOutputList->Add(fCrossCheckRowsLengthAcc);
  fOutputList->Add(fCrossCheckClusterLengthAcc);
  fOutputList->Add(fCrossCheckPtresLength);
  fOutputList->Add(fCrossCheckPtresRows);
  fOutputList->Add(fCutSettings);
  fOutputList->Add(fEventplaneDist);
  fOutputList->Add(fEventplaneRunDist);
  fOutputList->Add(fMCEventplaneDist);
  fOutputList->Add(fCorrelEventplaneMCDATA);
  fOutputList->Add(fCorrelEventplaneDefaultCorrected);
  fOutputList->Add(fEventplaneSubtractedPercentage);
  
  fOutputList->Add(fChargeOverPtRuns);
  fOutputList->Add(fVZEROMultCentrality);
  fOutputList->Add(fVEROMultRefMult);
  
  fOutputList->Add(fEPDistCent);
  fOutputList->Add(fPhiCent);
  fOutputList->Add(fPcosEPCent);
  fOutputList->Add(fPsinEPCent);
  fOutputList->Add(fPcosPhiCent);
  fOutputList->Add(fPsinPhiCent);
  fOutputList->Add(fEPContributionDifference);
  
  fOutputList->Add(fDeltaPhiCent);
  fOutputList->Add(fDeltaPhiSymCent);
  
  fOutputList->Add(fMCRecTracksMult);
  fOutputList->Add(fMCGenTracksMult);
  
  fOutputList->Add(fCrossCheckFilterBitPhiCent);
  
  fOutputList->Add(fTriggerStringsFired);
  fOutputList->Add(fTriggerStringComplete);
  
  fOutputList->Add(fVertexZ);
  fOutputList->Add(fVertexZSPD);
  fOutputList->Add(fVertexZTPC);
  fOutputList->Add(fDeltaVertexZGlobalSPD);
  fOutputList->Add(fDeltaVertexZGlobalTPC);
  fOutputList->Add(fVertexContributors);
  
  fOutputList->Add(fVertexZAfterCuts);
  fOutputList->Add(fVertexZSPDAfterCuts);
  fOutputList->Add(fVertexZTPCAfterCuts);
  fOutputList->Add(fDeltaVertexZGlobalSPDAfterCuts);
  fOutputList->Add(fDeltaVertexZGlobalTPCAfterCuts);
  fOutputList->Add(fVertexContributorsAfterCuts);
  
  fOutputList->Add(fVertexContributorsAfterCutsCent);
  fOutputList->Add(fVertexContributorsAfterCutsSemi);
  fOutputList->Add(fVertexContributorsAfterCutsMB);
  
  StoreCutSettingsToHistogram();
  
  PostData(1, fOutputList);
}

void AlidNdPtAnalysisPbPbAOD::UserExec(Option_t *option)
{
  //
  // Main Loop
  // called for each event
  //
  //cout << fBinsPhi[fPhiNbins-1] << endl;
  fEventStatistics->Fill("all events",1);
  
  // set ZERO pointers:
  AliInputEventHandler *inputHandler = NULL;
  AliAODTrack *track = NULL;
  AliAODMCParticle *mcPart = NULL;
  AliAODMCHeader *mcHdr = NULL;
  AliGenHijingEventHeader *genHijingHeader = NULL;
  //AliGenPythiaEventHeader *genPythiaHeader = NULL;
  AliEventplane *ep = NULL;
  
  TVector2 *epQvector = NULL;
  
  Bool_t bIsEventSelectedMB = kFALSE;
  Bool_t bIsEventSelectedSemi = kFALSE;
  Bool_t bIsEventSelectedCentral = kFALSE;
  Bool_t bIsEventSelected = kFALSE;
  Bool_t bIsPrimary = kFALSE;
  Bool_t bIsHijingParticle = kFALSE;
  Bool_t bMotherIsHijingParticle = kFALSE;
  //Bool_t bIsPythiaParticle = kFALSE;
  Bool_t bEventHasATrack = kFALSE;
  Bool_t bEventHasATrackInRange = kFALSE;
  Int_t nTriggerFired = 0;
  
  
  Double_t dMCTrackZvPtEtaCent[4] = {0};
  Double_t dTrackZvPtEtaCent[4] = {0};
  
  Double_t dMCTrackDeltaphiPtEtaPhiCent[5] = {0};
  Double_t dTrackDeltaphiPtEtaPhiCent[5] = {0};
  
  Double_t dMCTrackPtEtaPhiZvCent[5] = {0};
  Double_t dTrackDeltaphiPtEtaPhiZvCent[6] = {0};
  
  Double_t dMCTrackPtEtaPhiCent[4] = {0};
  Double_t dTrackPtEtaPhiCent[4] = {0};
  
  Double_t dDCA[2] = {0};
  
  Double_t dMCEventZv = -100;
  Double_t dEventZv = -100;
  Double_t dEventZvSPD = -100;
  Double_t dEventZvTPC = -100;
  Int_t iAcceptedMultiplicity = 0;
  Double_t dEventplaneAngle = -10;
  Double_t dEventplaneAngleCorrected = -10; // event plane angle, where tracks contributing to this angle have been subtracted
  Double_t dMCEventplaneAngle = -10;
  
  Double_t dReferenceMultiplicity = 0;
  Float_t dMCNRecTracks = 0;
  Float_t dMCNGenTracks = 0;
  
  Bool_t bHasEntriesInEventByEventPtSpectrum = kFALSE;
  
  fIsMonteCarlo = kFALSE;
  
  AliAODEvent *eventAOD = 0x0;
  eventAOD = dynamic_cast<AliAODEvent*>( InputEvent() );
  if (!eventAOD) {
    AliWarning("ERROR: eventAOD not available \n");
    return;
  }
  
  
  // check, which trigger has been fired
  inputHandler = (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  // only take tracks of events, which are triggered
  bIsEventSelected = ( inputHandler->IsEventSelected() & GetCollisionCandidates() );
  if(!bIsEventSelected) { return; }
  
  //   Bool_t isMB = (event->GetTriggerMask() & (ULong64_t(1)<<1));
  //   Bool_t isCentral = (event->GetTriggerMask() & (ULong64_t(1)<<4));
  //   Bool_t isSemiCentral = (event->GetTriggerMask() & (ULong64_t(1)<<7));
  
  //   bIsEventSelectedMB = ( inputHandler->IsEventSelected() & AliVEvent::kMB);
  //   bIsEventSelectedSemi = ( inputHandler->IsEventSelected() & AliVEvent::kSemiCentral);
  //   bIsEventSelectedCentral = ( inputHandler->IsEventSelected() & AliVEvent::kCentral);
  
  bIsEventSelectedMB = (eventAOD->GetTriggerMask() & (ULong64_t(1)<<1));
  bIsEventSelectedCentral = (eventAOD->GetTriggerMask() & (ULong64_t(1)<<4));
  bIsEventSelectedSemi = (eventAOD->GetTriggerMask() & (ULong64_t(1)<<7));
  
  if(bIsEventSelectedMB || bIsEventSelectedSemi || bIsEventSelectedCentral) fTriggerStatistics->Fill("all triggered events",1);
  if(bIsEventSelectedMB) { fTriggerStatistics->Fill("MB trigger",1); nTriggerFired++; }
  if(bIsEventSelectedSemi) { fTriggerStatistics->Fill("SemiCentral trigger",1); nTriggerFired++; }
  if(bIsEventSelectedCentral) { fTriggerStatistics->Fill("Central trigger",1); nTriggerFired++; }
  if(nTriggerFired == 0) { fTriggerStatistics->Fill("No trigger",1); }
  
  //   cout << "Fired Trigger Classes: " << eventAOD->GetFiredTriggerClasses().Data() << endl;
  
  TString sFiredTrigger = eventAOD->GetFiredTriggerClasses();
  
  // do not use some of the triggers
  TString sDisabledOnlineTrigger = GetDisabledOnlineTrigger();
  TObjArray *oaDisabledTrigger = sDisabledOnlineTrigger.Tokenize(" ");
  for(Int_t iString = 0; iString < oaDisabledTrigger->GetEntries(); iString++)
  {
    TObjString *os = (TObjString*)oaDisabledTrigger->At(iString);
    if(sFiredTrigger.Contains(os->GetString())) return;
  }
  
  // store trigger strings to histogram
  TObjArray *oaFiredTrigger = sFiredTrigger.Tokenize(" ");
  for(Int_t iString = 0; iString < oaFiredTrigger->GetEntries(); iString++)
  {
    TObjString *os = (TObjString*)oaFiredTrigger->At(iString);
    fTriggerStringsFired->Fill(os->GetString().Data(),1);
  }
  
  fTriggerStringComplete->Fill(sFiredTrigger.Data(), 1);
  
  //  if(eventAOD->IsPileupFromSPD()) cout << "Event is pileup " << endl;
  
  //   if(nTriggerFired == 0) { return; }
  //   if( !bIsEventSelected || nTriggerFired>1 ) return;
  //   fEventStatistics->Fill("events with only coll. cand.", 1);
  
  // check if there is a stack, if yes, then do MC loop
  TList *list = eventAOD->GetList();
  TClonesArray *stack = 0x0;
  stack = (TClonesArray*)list->FindObject(AliAODMCParticle::StdBranchName());
  
  if( stack )
  {
    fIsMonteCarlo = kTRUE;
    
    mcHdr = (AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());
    
    genHijingHeader = GetHijingEventHeader(mcHdr);
    //     genPythiaHeader = GetPythiaEventHeader(mcHdr);
    
    if(!genHijingHeader) { return; }
    
    //     if(!genPythiaHeader)  { return; }
    
    
    dMCEventZv = mcHdr->GetVtxZ();
    dMCTrackZvPtEtaCent[0] = dMCEventZv;
    dMCEventplaneAngle = genHijingHeader->ReactionPlaneAngle();//MoveEventplane(genHijingHeader->ReactionPlaneAngle());
    fEventStatistics->Fill("MC all events",1);
    fMCEventplaneDist->Fill(dMCEventplaneAngle);
  }
  
  AliCentrality* aCentrality = eventAOD->GetCentrality();
  //   Double_t dCentrality = aCentrality->GetCentralityPercentile("V0M");
  Double_t dCentrality = aCentrality->GetCentralityPercentile(GetCentralityEstimator().Data());
  
  if(GetDoMinBiasAnalysis()) { dCentrality = 1; }
  if( dCentrality < 0 ) return;
  
  // protection for bias on pt spectra if all triggers selected
  if( (bIsEventSelectedCentral)  && (dCentrality > 10) ) return;
  //   if( (bIsEventSelectedSemi) && ((dCentrality < 20) || (dCentrality > 50))) return;
  if( (bIsEventSelectedSemi) &&  (dCentrality > 50) ) return;
  
  if(DoAnchorPointSystStudy() && GetAnchorPointCorrectionFactor()>0)
  {
    dCentrality = dCentrality/GetAnchorPointCorrectionFactor(); // division of centrality value is same as moving bins upwards
  }
  
  fEventStatistics->Fill("after centrality selection",1);
  
  // get event plane Angle from AODHeader, default is Q
  ep = const_cast<AliAODEvent*>(eventAOD)->GetEventplane();
  if(ep) {
    dEventplaneAngle = ep->GetEventplane(GetEventplaneSelector().Data(),eventAOD);//MoveEventplane(ep->GetEventplane(GetEventplaneSelector().Data(),eventAOD));
    if(GetEventplaneSelector().CompareTo("Q") == 0)
    {
      epQvector = ep->GetQVector();
      if(epQvector) dEventplaneAngle = epQvector->Phi()/2.;//MoveEventplane(epQvector->Phi());
    }
  }
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //   AliEPSelectionTask *eptask = 0x0;
  //   eptask = dynamic_cast<AliEPSelectionTask *>(man->GetTask("EventplaneSelection"));
  
  if( (GetEventplaneSelector().CompareTo("Q") == 0) && !epQvector )
  {
    AliWarning("ERROR: epQvector not available \n");
    return;
  }
  
  //   cout << dEventplaneAngle << endl;
  fEventplaneDist->Fill(dEventplaneAngle);
  fEventplaneRunDist->Fill(dEventplaneAngle, (Double_t)eventAOD->GetRunNumber());
  
  // fill crosscheck histos
  fEPDistCent->Fill(dEventplaneAngle, dCentrality);
  fPcosEPCent->Fill(dCentrality, TMath::Cos(2.*dEventplaneAngle));
  fPsinEPCent->Fill(dCentrality, TMath::Sin(2.*dEventplaneAngle));
  
  // start with MC truth analysis
  if(fIsMonteCarlo)
  {
    
    //    TList *genHeaders = mcHdr->GetCocktailHeaders();
    //
    //    if(!genHeaders)     {
    //      AliError("ERROR: Could not retrieve genHeaders");
    //      return;
    //    }
    //
    //    Int_t nGenerators = genHeaders->GetEntries();
    //
    //    printf("N generators %d \n", nGenerators);
    //
    //    // igen==1 and igen==2 are pi0 and eta added signals, respectively, when looking at lhc13ef
    //    // check eventHeader2->GetName(); for your dataset
    //    for(Int_t igen = 0; igen < nGenerators; igen++) {
    //      AliGenEventHeader* eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
    //      TString name = eventHeader2->GetName();
    ////      cout << "Generator: " << igen << " " << name.Data() << endl;
    //    }
    
    if( dMCEventZv > GetCutMaxZVertex() )  { return; }
    
    dMCTrackZvPtEtaCent[0] = dMCEventZv;
    
    fEventStatistics->Fill("MC afterZv cut",1);
    
    for(Int_t iMCtrack = 0; iMCtrack < stack->GetEntriesFast(); iMCtrack++)
    {
      mcPart =(AliAODMCParticle*)stack->At(iMCtrack);
      // check for charge
      if( !(IsMCTrackAccepted(mcPart)) ) continue;
      
      if(!IsHijingParticle(mcPart, genHijingHeader)) { continue; }
      
      if(mcPart->IsPhysicalPrimary() )
      {
        // 	fMCHijingPrim->Fill("IsPhysicalPrimary",1);
      }
      else
      {
        // 	fMCHijingPrim->Fill("NOT a primary",1);
        continue;
      }
      
      
      //
      // ======================== fill histograms ========================
      dMCTrackZvPtEtaCent[1] = mcPart->Pt();
      dMCTrackZvPtEtaCent[2] = mcPart->Eta();
      dMCTrackZvPtEtaCent[3] = dCentrality;
      fMCGenZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
      
      dMCTrackDeltaphiPtEtaPhiCent[0] = RotatePhi(mcPart->Phi(), dEventplaneAngle, fBinsDeltaphi[fDeltaphiNbins-1]); // use eventplane and not reactionplan, similar to centrality vs impact paramter
      // 	  if( dMCTrackDeltaphiPtEtaPhiCent[0] < 0) dMCTrackDeltaphiPtEtaPhiCent[0] += 2.*TMath::Pi();
      // 	  else if( dMCTrackDeltaphiPtEtaPhiCent[0] > 2.*TMath::Pi()) dMCTrackDeltaphiPtEtaPhiCent[0] -= 2.*TMath::Pi();
      dMCTrackDeltaphiPtEtaPhiCent[1] = mcPart->Pt();
      dMCTrackDeltaphiPtEtaPhiCent[2] = mcPart->Eta();
      dMCTrackDeltaphiPtEtaPhiCent[3] = mcPart->Phi();
      dMCTrackDeltaphiPtEtaPhiCent[4] = dCentrality;
      
      dMCTrackPtEtaPhiCent[0] = mcPart->Pt();
      dMCTrackPtEtaPhiCent[1] = mcPart->Eta();
      dMCTrackPtEtaPhiCent[2] = mcPart->Phi();
      dMCTrackPtEtaPhiCent[3] = dCentrality;
      
      fMCGenPtEtaPhiCent->Fill(dMCTrackPtEtaPhiCent);
      
      
      dMCTrackPtEtaPhiZvCent[0] = mcPart->Pt();
      dMCTrackPtEtaPhiZvCent[1] = mcPart->Eta();
      dMCTrackPtEtaPhiZvCent[2] = mcPart->Phi();
      dMCTrackPtEtaPhiZvCent[3] = dMCEventZv;
      dMCTrackPtEtaPhiZvCent[4] = dCentrality;
      fMCGenPtEtaPhiZvCent->Fill(dMCTrackPtEtaPhiZvCent);
      
      bEventHasATrack = kTRUE;
      
      
      if( (dMCTrackZvPtEtaCent[1] > GetCutPtMin() ) &&
         (dMCTrackZvPtEtaCent[1] < GetCutPtMax() ) &&
         (dMCTrackZvPtEtaCent[2] > GetCutEtaMin() ) &&
         (dMCTrackZvPtEtaCent[2] < GetCutEtaMax() ) )
      {
        fMCPt->Fill(mcPart->Pt());
        fMCCharge->Fill(mcPart->Charge()/3.);
        bEventHasATrackInRange = kTRUE;
        dMCNGenTracks++;
      }
      
    }
  } // isMonteCarlo
  
  if(bEventHasATrack) { fEventStatistics->Fill("MC events with tracks",1); }
  if(bEventHasATrackInRange)
  {
    fEventStatistics->Fill("MC events with tracks in range",1);
    fMCEventStatisticsCentrality->Fill(dCentrality);
  }
  bEventHasATrack = kFALSE;
  bEventHasATrackInRange = kFALSE;
  
  
  //
  // Loop over recontructed tracks
  //
  
  dEventZv = eventAOD->GetPrimaryVertex()->GetZ();
  dEventZvSPD = eventAOD->GetPrimaryVertexSPD()->GetZ();
  dEventZvTPC = eventAOD->GetPrimaryVertexTPC()->GetZ();
  
  fVertexZ->Fill(dEventZv);
  fVertexZSPD->Fill(dEventZvSPD);
  fVertexZTPC->Fill(dEventZvTPC);
  fDeltaVertexZGlobalSPD->Fill(dEventZv - dEventZvSPD);
  fDeltaVertexZGlobalTPC->Fill(dEventZv - dEventZvTPC);
  fVertexContributors->Fill(eventAOD->GetPrimaryVertex()->GetNContributors());
  
  if( TMath::Abs(dEventZv) > GetCutMaxZVertex() ) return;
  if( eventAOD->GetPrimaryVertex()->GetNContributors() < GetNContributorsVertex() ) return;
  
  fVertexZAfterCuts->Fill(dEventZv);
  fVertexZSPDAfterCuts->Fill(dEventZvSPD);
  fVertexZTPCAfterCuts->Fill(dEventZvTPC);
  fDeltaVertexZGlobalSPDAfterCuts->Fill(dEventZv - dEventZvSPD);
  fDeltaVertexZGlobalTPCAfterCuts->Fill(dEventZv - dEventZvTPC);
  fVertexContributorsAfterCuts->Fill(eventAOD->GetPrimaryVertex()->GetNContributors());
  if(bIsEventSelectedCentral) { fVertexContributorsAfterCutsCent->Fill(eventAOD->GetPrimaryVertex()->GetNContributors()); }
  if(bIsEventSelectedSemi) { fVertexContributorsAfterCutsSemi->Fill(eventAOD->GetPrimaryVertex()->GetNContributors()); }
  if(bIsEventSelectedMB) { fVertexContributorsAfterCutsMB->Fill(eventAOD->GetPrimaryVertex()->GetNContributors()); }
  
  // count all events, which are within zv distribution
  fAllEventStatisticsCentrality->Fill(dCentrality/*, nTriggerFired*/);
  
  fEventStatistics->Fill("after Zv cut",1);
  
  Double_t dTotMultVZERO = -1.;
  for(Int_t iVZERObin = 0; iVZERObin < 64; iVZERObin++)
  {
    // 	dTotMultVZERO += eventAOD->GetVZEROEqMultiplicity(iVZERObin);
    dTotMultVZERO += eventAOD->GetVZEROData()->GetMultiplicity(iVZERObin);
  }
  fVZEROMultCentrality->Fill(dTotMultVZERO, dCentrality);
  
  
  dTrackZvPtEtaCent[0] = dEventZv;
  
  
  
  if(AreRelativeCutsEnabled())
  {
    if(!SetRelativeCuts(eventAOD)) return;
  }
  
  Int_t iSubtractedTracks = 0;
  
  for(Int_t itrack = 0; itrack < eventAOD->GetNumberOfTracks(); itrack++)
  {
    track = dynamic_cast<AliAODTrack*>(eventAOD->GetTrack(itrack));
    if(!track) AliFatal("Not a standard AOD");
    
    if(!track) continue;
    
    mcPart = NULL;
    dMCTrackZvPtEtaCent[1] = 0;
    dMCTrackZvPtEtaCent[2] = 0;
    dMCTrackZvPtEtaCent[3] = 0;
    
    dMCTrackDeltaphiPtEtaPhiCent[0] = 0;
    dMCTrackDeltaphiPtEtaPhiCent[1] = 0;
    dMCTrackDeltaphiPtEtaPhiCent[2] = 0;
    dMCTrackDeltaphiPtEtaPhiCent[3] = 0;
    dMCTrackDeltaphiPtEtaPhiCent[4] = 0;
    
    dMCTrackPtEtaPhiZvCent[0] = 0;
    dMCTrackPtEtaPhiZvCent[1] = 0;
    dMCTrackPtEtaPhiZvCent[2] = 0;
    dMCTrackPtEtaPhiZvCent[3] = 0;
    dMCTrackPtEtaPhiZvCent[4] = 0;
    
    dMCTrackPtEtaPhiCent[0] = 0;
    dMCTrackPtEtaPhiCent[1] = 0;
    dMCTrackPtEtaPhiCent[2] = 0;
    dMCTrackPtEtaPhiCent[3] = 0;
    
    bIsPrimary = kFALSE;
    
    GetDCA(track, eventAOD, dDCA);
    
    Double_t dDCAxyDCAzPt[5] = { dDCA[0], dDCA[1], track->Pt(), track->Eta(), track->Phi() };
    
    fDCAPtAll->Fill(dDCAxyDCAzPt);
    
    if( !(IsTrackAccepted(track, dCentrality, eventAOD->GetMagneticField())) ) continue;
    
    dTrackZvPtEtaCent[1] = track->Pt();
    dTrackZvPtEtaCent[2] = track->Eta();
    dTrackZvPtEtaCent[3] = dCentrality;
    
    
    
    if(GetEventplaneSelector().CompareTo("Q") == 0)
    {
      // subtract track contribution from eventplane
      Double_t dX = -1000;
      Double_t dY = -1000;
      
      dX = epQvector->X();
      dY = epQvector->Y();
      
      // 	  if(eptask)
      // 	  {
      // 		Double_t rms[2] ={1.,1.};
      // 		eptask->Recenter(1, rms);
      //
      // 		Double_t dContribX = (eptask->GetWeight(track) * TMath::Cos(2*track->Phi()) / rms[0]);
      // 		Double_t dContribY = (eptask->GetWeight(track) * TMath::Sin(2*track->Phi()) / rms[1]);
      // 		// 		dX -= (eptask->GetWeight(track) * TMath::Cos(2*track->Phi()) / rms[0]);
      // 		//         dY -= (eptask->GetWeight(track) * TMath::Sin(2*track->Phi()) / rms[1]);
      // 		dX -= dContribX;
      // 		dY -= dContribY;
      // 		iSubtractedTracks++;
      // 		if(track->GetID()>0)
      // 		{
      // 		  if( (dContribX != 0) && (dContribY != 0) ) {
      // 			fEPContributionDifference->Fill( (dContribX-ep->GetQContributionX(track))/dContribX, (dContribY-ep->GetQContributionY(track))/dContribY);
      // 		  }
      // 		}
      // 	  }
      if( (dX>-1000) && (dY>-1000) && (track->GetID()>0) ) // only subtract, if not default and track Id > 0 to avoid crash
      {
        dX -= ep->GetQContributionX(track);
        dY -= ep->GetQContributionY(track);
        iSubtractedTracks++;
      }
      
      TVector2 epCorrected(dX, dY);
      dEventplaneAngleCorrected = epCorrected.Phi()/2.; // see AlEPSelectionTask.cxx:354
    }
    else
    {
      dEventplaneAngleCorrected = dEventplaneAngle;
    }
    Double_t dFillEPCorrectionCheck[] = {dEventplaneAngle, dEventplaneAngleCorrected, dCentrality};
    fCorrelEventplaneDefaultCorrected->Fill(dFillEPCorrectionCheck);
    
    
    dTrackDeltaphiPtEtaPhiCent[0] = RotatePhi(track->Phi(), dEventplaneAngleCorrected, fBinsDeltaphi[fDeltaphiNbins-1]);
    
    dTrackDeltaphiPtEtaPhiCent[1] = track->Pt();
    dTrackDeltaphiPtEtaPhiCent[2] = track->Eta();
    dTrackDeltaphiPtEtaPhiCent[3] = track->Phi();
    dTrackDeltaphiPtEtaPhiCent[4] = dCentrality;
    
    dTrackDeltaphiPtEtaPhiZvCent[0] = dTrackDeltaphiPtEtaPhiCent[0];
    dTrackDeltaphiPtEtaPhiZvCent[1] = track->Pt();
    dTrackDeltaphiPtEtaPhiZvCent[2] = track->Eta();
    dTrackDeltaphiPtEtaPhiZvCent[3] = track->Phi();
    dTrackDeltaphiPtEtaPhiZvCent[4] = dEventZv;
    dTrackDeltaphiPtEtaPhiZvCent[5] = dCentrality;
    
    dTrackPtEtaPhiCent[0] = track->Pt();
    dTrackPtEtaPhiCent[1] = track->Eta();
    dTrackPtEtaPhiCent[2] = track->Phi();
    dTrackPtEtaPhiCent[3] = dCentrality;
    
    
    if( fIsMonteCarlo )
    {
      mcPart = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
      if( !mcPart ) { continue; }
      
      // check for charge
      // if( !(IsMCTrackAccepted(mcPart)) ) {  continue; }
      
      bIsHijingParticle = IsHijingParticle(mcPart, genHijingHeader);
      //       bIsPythiaParticle = IsPythiaParticle(mcPart, genPythiaHeader);
      
//      if(!bIsHijingParticle && mcPart->GetGeneratorIndex()==0)
//      cout << "IsHijing: " << bIsHijingParticle << " GeneratorIndex: " << mcPart->GetGeneratorIndex() <<  " label: " << mcPart->Label() << " NProduced-1: " <<(genHijingHeader->NProduced()-1) << endl;
      
      bIsPrimary = mcPart->IsPhysicalPrimary();
      
      dMCTrackZvPtEtaCent[1] = mcPart->Pt();
      dMCTrackZvPtEtaCent[2] = mcPart->Eta();
      dMCTrackZvPtEtaCent[3] = dCentrality;
      
      dMCTrackDeltaphiPtEtaPhiCent[0] = RotatePhi(mcPart->Phi(), dEventplaneAngle, fBinsDeltaphi[fDeltaphiNbins-1]); // use eventplane and not reactionplan, similar to centrality vs impact paramter
      
      dMCTrackDeltaphiPtEtaPhiCent[1] = mcPart->Pt();
      dMCTrackDeltaphiPtEtaPhiCent[2] = mcPart->Eta();
      dMCTrackDeltaphiPtEtaPhiCent[3] = mcPart->Phi();
      dMCTrackDeltaphiPtEtaPhiCent[4] = dCentrality;
      
      dMCTrackPtEtaPhiZvCent[0] = mcPart->Pt();
      dMCTrackPtEtaPhiZvCent[1] = mcPart->Eta();
      dMCTrackPtEtaPhiZvCent[2] = mcPart->Phi();
      dMCTrackPtEtaPhiZvCent[3] = dEventZv;
      dMCTrackPtEtaPhiZvCent[4] = dCentrality;
      
      dMCTrackPtEtaPhiCent[0] = mcPart->Pt();
      dMCTrackPtEtaPhiCent[1] = mcPart->Eta();
      dMCTrackPtEtaPhiCent[2] = mcPart->Phi();
      dMCTrackPtEtaPhiCent[3] = dCentrality;
      
      if(bIsPrimary && bIsHijingParticle)
      {
        fMCRecPrimZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
        fMCRecPrimPtEtaPhiCent->Fill(dMCTrackPtEtaPhiCent);
        fMCRecPrimPtEtaPhiZvCent->Fill(dMCTrackPtEtaPhiZvCent);
        fMCDCAPtPrimary->Fill(dDCAxyDCAzPt);
        if( (dMCTrackZvPtEtaCent[1] > GetCutPtMin()) &&
           (dMCTrackZvPtEtaCent[1] < GetCutPtMax()) &&
           (dMCTrackZvPtEtaCent[2] > GetCutEtaMin()) &&
           (dMCTrackZvPtEtaCent[2] < GetCutEtaMax()) )
        {
          dMCNRecTracks++;
        }
        
      }
      
      if(!bIsPrimary /*&& !bIsHijingParticle*/)
      {
        Int_t indexMoth = mcPart->GetMother();
        if(indexMoth >= 0)
        {
          AliAODMCParticle* moth = (AliAODMCParticle*)stack->At(indexMoth);
          bMotherIsHijingParticle = IsHijingParticle(moth, genHijingHeader);
          
          if(bMotherIsHijingParticle) // only store secondaries, which come from a not embedded signal!
          {
            fMCRecSecZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
            fMCRecSecPtEtaPhiCent->Fill(dMCTrackPtEtaPhiCent);
            fMCRecSecPtEtaPhiZvCent->Fill(dMCTrackPtEtaPhiZvCent);
            fMCDCAPtSecondary->Fill(dDCAxyDCAzPt);
            // 	  delete moth;
          }
        }
      }
    } // end isMonteCarlo
    
    // ======================== fill histograms ========================
    
    // only keep prim and sec from not embedded signal
    Bool_t bKeepMCTrack = kFALSE;
    if(fIsMonteCarlo)
    {
      if( (bIsHijingParticle && bIsPrimary) ^ (bMotherIsHijingParticle && !bIsPrimary) )
      {
        bKeepMCTrack = kTRUE;
      }
      else
      {
        continue;
      }
    }
    
    bEventHasATrack = kTRUE;
    
    fZvPtEtaCent->Fill(dTrackZvPtEtaCent);
    fDeltaphiPtEtaPhiCent->Fill(dTrackDeltaphiPtEtaPhiCent);
    fDeltaphiPtEtaPhiZvCent->Fill(dTrackDeltaphiPtEtaPhiZvCent);
    
    if(fIsMonteCarlo)
    {
      fMCPtEtaPhiCent->Fill(dTrackPtEtaPhiCent);
      fMCPtEtaPhiZvCent->Fill(dMCTrackPtEtaPhiZvCent);
    }
    
    fDCAPtAccepted->Fill(dDCAxyDCAzPt);
    
    if( (dTrackZvPtEtaCent[1] > GetCutPtMin()) &&
       (dTrackZvPtEtaCent[1] < GetCutPtMax()) &&
       (dTrackZvPtEtaCent[2] > GetCutEtaMin()) &&
       (dTrackZvPtEtaCent[2] < GetCutEtaMax()) )
    {
      iAcceptedMultiplicity++;
      bEventHasATrackInRange = kTRUE;
      fPt->Fill(track->Pt());
      fCharge->Fill(track->Charge());
      
      fPhiCent->Fill(track->Phi(), dCentrality);
      fPcosPhiCent->Fill(dCentrality, TMath::Cos(2.*track->Phi()));
      fPsinPhiCent->Fill(dCentrality, TMath::Sin(2.*track->Phi()));
      
      Double_t deltaphi = track->Phi() - dEventplaneAngleCorrected;
      // 	  if(deltaphi > TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      
      fDeltaPhiCent->Fill(deltaphi, dCentrality);
      fDeltaPhiSymCent->Fill(dTrackDeltaphiPtEtaPhiCent[0], dCentrality);
      
      fChargeOverPtRuns->Fill(track->Charge()/track->Pt(), (Double_t)eventAOD->GetRunNumber());
      
      if(GetFillEventPtSpectraHistogram() && (dCentrality < 5.)) {
        fPtEvent->Fill(track->Pt(), fEventNumberForPtSpectra);
        bHasEntriesInEventByEventPtSpectrum = kTRUE;
      }
    }
  } // end track loop
  Int_t iContributorsQVector = ep->GetQContributionXArray()->GetSize();
  if(iContributorsQVector) fEventplaneSubtractedPercentage->Fill((Double_t)iSubtractedTracks/(Double_t)iContributorsQVector, dCentrality);
  
  if(bEventHasATrack) { fEventStatistics->Fill("events with tracks",1); bEventHasATrack = kFALSE; }
  
  if(bEventHasATrackInRange)
  {
    fEventStatistics->Fill("events with tracks in range",1);
    fEventStatisticsCentrality->Fill(dCentrality);
    
    bEventHasATrackInRange = kFALSE;
  }
  
  if(bIsEventSelectedMB) fEventStatisticsCentralityTrigger->Fill(dCentrality, 0);
  if(bIsEventSelectedSemi) fEventStatisticsCentralityTrigger->Fill(dCentrality, 1);
  if(bIsEventSelectedCentral) fEventStatisticsCentralityTrigger->Fill(dCentrality, 2);
  
  Double_t dEventZvMultCent[3] = {dEventZv, static_cast<Double_t>(iAcceptedMultiplicity), dCentrality};
  fZvMultCent->Fill(dEventZvMultCent);
  
  AliAODHeader *aodHeader = (AliAODHeader*)eventAOD->GetHeader();
  dReferenceMultiplicity = aodHeader->GetRefMultiplicity();
  
  
  // store correlation between data and MC eventplane
  if(fIsMonteCarlo)
  {
    fCorrelEventplaneMCDATA->Fill(dEventplaneAngle, dMCEventplaneAngle);
    fMCRecTracksMult->Fill(dReferenceMultiplicity, dMCNRecTracks);
    fMCGenTracksMult->Fill(dReferenceMultiplicity, dMCNGenTracks);
  }
  
  PostData(1, fOutputList);
  
  fVEROMultRefMult->Fill(dTotMultVZERO, dReferenceMultiplicity);
  
  if(bHasEntriesInEventByEventPtSpectrum) fEventNumberForPtSpectra++;
  // delete pointers:
  //   delete [] iIndexAcceptedTracks;
}


Double_t AlidNdPtAnalysisPbPbAOD::RotatePhi(Double_t phiTrack, Double_t phiEP, Double_t dMaxDeltaPhi)
{
  Double_t dPhi = 0;
  dPhi = TMath::Abs(phiTrack - phiEP);
  
  //   if( dPhi <= TMath::Pi() )
  //   {
  // 	return dPhi;
  //   }
  if( dPhi > TMath::Pi() )
  {
    dPhi = 2.*TMath::Pi() - dPhi;
    // 	return dPhi;
  }
  
  if( dPhi > dMaxDeltaPhi)
  {
    dPhi = 2.*dMaxDeltaPhi - dPhi;
  }
  
  if(dPhi > dMaxDeltaPhi)
  {
    Printf("[E] dphi = %.4f , phiTrack = %.4f, phiEP = %.4f, maxDeltaPhi = %.4f", dPhi, phiTrack, phiEP, dMaxDeltaPhi);
  }
  
  //   return -9999.;
  
  return dPhi;
}

Bool_t AlidNdPtAnalysisPbPbAOD::SetRelativeCuts(AliAODEvent *event)
{
  //
  // this function determines the absolute cut event-by-event based on the
  // the percentage given from outside
  //  - cut set on Nclusters and NcrossedRows
  //
  
  if(!event) return kFALSE;
  
  AliAODTrack *tr = 0x0;
  TH1F *hCluster = new TH1F("hCluster","hCluster",160,0,160);
  TH1F *hCrossed = new TH1F("hCrossed","hCrossed",160,0,160);
  
  for(Int_t itrack = 0; itrack < event->GetNumberOfTracks(); itrack++)
  {
    tr = dynamic_cast<AliAODTrack*>(event->GetTrack(itrack));
    if(!tr) AliFatal("Not a standard AOD");
    if(!tr) continue;
    
    // do some selection already
    //if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { continue; }
    
    Double_t dNClustersTPC = tr->GetTPCNcls();
    Double_t dCrossedRowsTPC = tr->GetTPCClusterInfo(2,1);
    
    hCluster->Fill(dNClustersTPC);
    hCrossed->Fill(dCrossedRowsTPC);
  }
  
  // loop trough histogram to check, where percentage is reach
  Double_t dTotIntCluster = hCluster->Integral();
  Double_t dTotIntCrossed = hCrossed->Integral();
  Float_t dIntCluster = 0;
  Float_t dIntCrossed = 0;
  
  if(dTotIntCluster)
  {
    for(Int_t i = 0; i < hCluster->GetNbinsX(); i++)
    {
      if(hCluster->GetBinCenter(i) < 0) continue;
      dIntCluster += hCluster->GetBinContent(i);
      if(dIntCluster/dTotIntCluster > (1-GetCutPercMinNClustersTPC()))
      {
        SetCutMinNClustersTPC(hCluster->GetBinCenter(i));
        fCutPercClusters->Fill(hCluster->GetBinCenter(i));
        break;
      }
    }
  }
  
  if(dTotIntCrossed)
  {
    for(Int_t i = 0; i < hCrossed->GetNbinsX(); i++)
    {
      if(hCrossed->GetBinCenter(i) < 0) continue;
      dIntCrossed += hCrossed->GetBinContent(i);
      if(dIntCrossed/dTotIntCrossed > (1-GetCutPercMinNCrossedRowsTPC()))
      {
        SetCutMinNClustersTPC(hCrossed->GetBinCenter(i));
        fCutPercCrossed->Fill(hCrossed->GetBinCenter(i));
        break;
      }
    }
  }
  
  delete hCrossed;
  delete hCluster;
  return kTRUE;
  
}

Bool_t AlidNdPtAnalysisPbPbAOD::IsTrackAccepted(AliAODTrack *tr, Double_t dCentrality, Double_t bMagZ)
{
  //
  // this function checks the track parameters for quality
  // returns kTRUE if track is accepted
  //
  // - debug histograms (cuts vs pt,eta,phi) are filled in this function
  // - histogram for pt resolution correction are filled here as well
  //
  
  if(!tr) return kFALSE;
  
  if(tr->Charge()==0) { return kFALSE; }
  
  //
  // as done in AliAnalysisTaskFragmentationFunction
  //
  
  Short_t sign = tr->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[21];
  
  for(Int_t i = 0; i < 21; i++) cv[i] = 0;
  for(Int_t i = 0; i < 50; i++) xyz[i] = 0;
  for(Int_t i = 0; i < 50; i++) pxpypz[i] = 0;
  
  tr->GetXYZ(xyz);
  tr->GetPxPyPz(pxpypz);
  tr->GetCovarianceXYZPxPyPz(cv);
  
  // similar error occured as this one:
  // See https://savannah.cern.ch/bugs/?102721
  // which is one of the two 11h re-filtering follow-ups:
  // Andrea Dainese now first does the beam pipe
  // check and then copies from the vtrack (was the other
  // way around) to avoid the crash in the etp::Set()
  
  //   if(xyz[0]*xyz[0]+xyz[1]*xyz[1] > 3.*3.) { return kFALSE; }
  
  AliExternalTrackParam par(xyz, pxpypz, cv, sign);
  //   AliExternalTrackParam *par = new AliExternalTrackParam(xyz, pxpypz, cv, sign); // high mem consumption!!!!
  static AliESDtrack dummy;
  //   Double_t dLength = dummy.GetLengthInActiveZone(par,3,236, -5 ,0,0);
  //   Double_t dLengthInTPC = GetLengthInTPC(tr, 1.8, 220, bMagZ);
  
  Double_t dLengthInTPC = 0;
  if ( DoCutLengthInTPCPtDependent() ) { dLengthInTPC = dummy.GetLengthInActiveZone(&par,3,236, bMagZ ,0,0); }
  
  Double_t dNClustersTPC = tr->GetTPCNcls();
  Double_t dCrossedRowsTPC = tr->GetTPCNCrossedRows();//GetTPCClusterInfo(2,1);
  Double_t dFindableClustersTPC = tr->GetTPCNclsF();
  Double_t dChi2PerClusterTPC = (dNClustersTPC>0)?tr->Chi2perNDF()*(dNClustersTPC-5)/dNClustersTPC:-1.; // see AliDielectronVarManager.h
  Double_t dOneOverPt = tr->OneOverPt();
  Double_t dSigmaOneOverPt = TMath::Sqrt(par.GetSigma1Pt2());
  
  //   hAllCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  Double_t dCrossedRowsTPCOverFindableClustersTPC = 0;
  if(dFindableClustersTPC) dCrossedRowsTPCOverFindableClustersTPC = dCrossedRowsTPC/dFindableClustersTPC;
  Double_t dCheck[cqMax] = {dCrossedRowsTPC, dNClustersTPC, dChi2PerClusterTPC, dLengthInTPC, dCrossedRowsTPCOverFindableClustersTPC};// = new Double_t[cqMax];
  Double_t dKine[kqMax] = {tr->Pt(), tr->Eta(), tr->Phi()};// = new Double_t[kqMax];
  
  //   dKine[0] = tr->Pt();
  //   dKine[1] = tr->Eta();
  //   dKine[2] = tr->Phi();
  //
  //   dCheck[0] = dCrossedRowsTPC;
  //   dCheck[1] = dNClustersTPC;
  //   dCheck[2] = dChi2PerClusterTPC;
  
  
  FillDebugHisto(dCheck, dKine, dCentrality, kFALSE);
  
  fCrossCheckPtresLength->Fill(dLengthInTPC, dSigmaOneOverPt*tr->Pt(), dCentrality);
  fCrossCheckPtresRows->Fill(dCrossedRowsTPC, dSigmaOneOverPt*tr->Pt(), dCentrality);
  
  // filter bit 5
  //   if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { return kFALSE; }
  
  //   if(!(tr->TestFilterBit(GetFilterBit())) ) { return kFALSE; }
  
  if( RequireHybridTracking() && !(tr->IsHybridGlobalConstrainedGlobal()) ) { return kFALSE; }
  if( !(RequireHybridTracking()) && !(tr->TestFilterBit(GetFilterBit())) ) { return kFALSE; }
  
  
  
  // cut on length
  if( DoCutLengthInTPCPtDependent() && ( dLengthInTPC < GetPrefactorLengthInTPCPtDependent()*(130-5*TMath::Abs(1./tr->Pt())) )  ) { return kFALSE; }
  
  
  
  // filter bit 4
  //   if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) ) { return kFALSE; }
  
  //   hFilterCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  
  if(dFindableClustersTPC == 0) {return kFALSE; }
  if(dCrossedRowsTPC < GetCutMinNCrossedRowsTPC()) { return kFALSE; }
  if( (dCrossedRowsTPCOverFindableClustersTPC) < GetCutMinRatioCrossedRowsOverFindableClustersTPC() ) { return kFALSE; }
  if(dNClustersTPC < GetCutMinNClustersTPC()) { return kFALSE; }
  
  if (IsITSRefitRequired() && !(tr->GetStatus() & AliVTrack::kITSrefit)) { return kFALSE; } // no ITS refit
  
  // fill histogram for pT resolution correction
  Double_t dPtResolutionHisto[3] = { dOneOverPt, dSigmaOneOverPt, dCentrality };
  fPtResptCent->Fill(dPtResolutionHisto);
  
  Double_t dPtResolutionHisto2[3] = {tr->Pt(), dSigmaOneOverPt*tr->Pt(), dCentrality };
  fPtResptptCent->Fill(dPtResolutionHisto2);
  
  // fill debug histogram for all accepted tracks
  FillDebugHisto(dCheck, dKine, dCentrality, kTRUE);
  
  Double_t dFilterBitPhiCent[3] = {-10, -10, -10};
  if(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) dFilterBitPhiCent[0] = 0;
  else if(tr->TestFilterBit(AliAODTrack::kTrkGlobalSDD)) dFilterBitPhiCent[0] = 1;
  
  dFilterBitPhiCent[1] = tr->Phi();
  dFilterBitPhiCent[2] = dCentrality;
  fCrossCheckFilterBitPhiCent->Fill(dFilterBitPhiCent);
  
  // delete pointers
  
  return kTRUE;
}

Bool_t AlidNdPtAnalysisPbPbAOD::FillDebugHisto(Double_t *dCrossCheckVar, Double_t *dKineVar, Double_t dCentrality, Bool_t bIsAccepted){
  if(bIsAccepted)
  {
    if(AreCrossCheckCorrelationHistosEnabled())
    {
      for(Int_t iCrossCheck = 0; iCrossCheck < cqMax; iCrossCheck++)
      {
        Double_t dFillIt[5] = {dCrossCheckVar[iCrossCheck], dKineVar[0], dKineVar[1], dKineVar[2], dCentrality};
        fCrossCheckAcc[iCrossCheck]->Fill(dFillIt);
      }
    }
    
    fCrossCheckRowsLengthAcc->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqCrossedRows], dCentrality);
    fCrossCheckClusterLengthAcc->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqNcluster], dCentrality);
  }
  else
  {
    if(AreCrossCheckCorrelationHistosEnabled())
    {
      for(Int_t iCrossCheck = 0; iCrossCheck < cqMax; iCrossCheck++)
      {
        Double_t dFillIt[5] = {dCrossCheckVar[iCrossCheck], dKineVar[0], dKineVar[1], dKineVar[2], dCentrality};
        fCrossCheckAll[iCrossCheck]->Fill(dFillIt);
      }
    }
    
    fCrossCheckRowsLength->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqCrossedRows], dCentrality);
    fCrossCheckClusterLength->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqNcluster], dCentrality);
  }
  
  return kTRUE;
  
}

void AlidNdPtAnalysisPbPbAOD::StoreCutSettingsToHistogram()
{
  //
  // this function stores all cut settings to a histograms
  //
  
  fCutSettings->Fill("IsMonteCarlo",fIsMonteCarlo);
  
  fCutSettings->Fill("fCutMaxZVertex", fCutMaxZVertex);
  
  // kinematic cuts
  fCutSettings->Fill("fCutPtMin", fCutPtMin);
  fCutSettings->Fill("fCutPtMax", fCutPtMax);
  fCutSettings->Fill("fCutEtaMin", fCutEtaMin);
  fCutSettings->Fill("fCutEtaMax", fCutEtaMax);
  
  // track quality cut variables
  fCutSettings->Fill("fFilterBit", fFilterBit);
  if(fUseRelativeCuts) fCutSettings->Fill("fUseRelativeCuts", 1);
  if(fCutRequireTPCRefit) fCutSettings->Fill("fCutRequireTPCRefit", 1);
  if(fCutRequireITSRefit) fCutSettings->Fill("fCutRequireITSRefit", 1);
  if(fHybridTracking) fCutSettings->Fill("RequireHybridTracking", 1);
  
  fCutSettings->Fill("fCutMinNumberOfClusters", fCutMinNumberOfClusters);
  fCutSettings->Fill("fCutPercMinNumberOfClusters", fCutPercMinNumberOfClusters);
  fCutSettings->Fill("fCutMinNumberOfCrossedRows", fCutMinNumberOfCrossedRows);
  fCutSettings->Fill("fCutPercMinNumberOfCrossedRows", fCutPercMinNumberOfCrossedRows);
  
  fCutSettings->Fill("fCutMinRatioCrossedRowsOverFindableClustersTPC", fCutMinRatioCrossedRowsOverFindableClustersTPC);
  fCutSettings->Fill("fCutMaxFractionSharedTPCClusters", fCutMaxFractionSharedTPCClusters);
  fCutSettings->Fill("fCutMaxDCAToVertexXY", fCutMaxDCAToVertexXY);
  fCutSettings->Fill("fCutMaxChi2PerClusterITS", fCutMaxChi2PerClusterITS);
  
  if(fCutDCAToVertex2D) fCutSettings->Fill("fCutDCAToVertex2D", 1);
  if(fCutRequireSigmaToVertex) fCutSettings->Fill("fCutRequireSigmaToVertex",1);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar0", fCutMaxDCAToVertexXYPtDepPar0);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar1", fCutMaxDCAToVertexXYPtDepPar1);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar2", fCutMaxDCAToVertexXYPtDepPar2);
  
  if(fCutAcceptKinkDaughters) fCutSettings->Fill("fCutAcceptKinkDaughters", 1);
  fCutSettings->Fill("fCutMaxChi2TPCConstrainedGlobal", fCutMaxChi2TPCConstrainedGlobal);
  if(fCutLengthInTPCPtDependent) fCutSettings->Fill("fCutLengthInTPCPtDependent", 1);
  fCutSettings->Fill("fPrefactorLengthInTPCPtDependent", fPrefactorLengthInTPCPtDependent);
  fCutSettings->Fill(Form("EP selector %s", fEPselector.Data()), 1);
  
  
  
  fCutSettings->Fill("NContributorsVertex", fVertexMinContributors);
}

Bool_t AlidNdPtAnalysisPbPbAOD::GetDCA(const AliAODTrack *track, AliAODEvent *evt, Double_t d0z0[2])
{
  // function adapted from AliDielectronVarManager.h
  
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Bool_t ok=kFALSE;
  if(evt) {
    Double_t covd0z0[3];
    //AliAODTrack copy(*track);
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);
    
    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      d0z0[0]=-999.;
      d0z0[1]=-999.;
      //printf("This method can be used only for propagation inside the beam pipe \n");
      return kFALSE;
    }
    
    
    AliAODVertex *vtx =(AliAODVertex*)(evt->GetPrimaryVertex());
    Double_t fBzkG = evt->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
    //ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}


Bool_t AlidNdPtAnalysisPbPbAOD::IsMCTrackAccepted(AliAODMCParticle *part)
{
  if(!part) return kFALSE;
  
  Double_t charge = part->Charge()/3.;
  if (TMath::Abs(charge) < 0.001) return kFALSE;
  
  return kTRUE;
}

const char * AlidNdPtAnalysisPbPbAOD::GetParticleName(Int_t pdg)
{
  TParticlePDG * p1 = TDatabasePDG::Instance()->GetParticle(pdg);
  if(p1) return p1->GetName();
  return Form("%d", pdg);
}

AliGenHijingEventHeader* AlidNdPtAnalysisPbPbAOD::GetHijingEventHeader(AliAODMCHeader *header)
{
  //
  // inspired by PWGJE/AliPWG4HighPtSpectra.cxx
  //
  
  if(!header) return 0x0;
  AliGenHijingEventHeader* hijingGenHeader = NULL;
  
  TList* headerList = header->GetCocktailHeaders();
  
  for(Int_t i = 0; i < headerList->GetEntries(); i++)
  {
    hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
    if(hijingGenHeader) break;
  }
  
  if(!hijingGenHeader) return 0x0;
  
  return hijingGenHeader;
}

AliGenPythiaEventHeader* AlidNdPtAnalysisPbPbAOD::GetPythiaEventHeader(AliAODMCHeader *header)
{
  //
  // inspired by PWGJE/AliPWG4HighPtSpectra.cxx
  //
  
  if(!header) return 0x0;
  AliGenPythiaEventHeader* PythiaGenHeader = NULL;
  
  TList* headerList = header->GetCocktailHeaders();
  
  for(Int_t i = 0; i < headerList->GetEntries(); i++)
  {
    PythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
    if(PythiaGenHeader) break;
  }
  
  if(!PythiaGenHeader) return 0x0;
  
  return PythiaGenHeader;
}

//________________________________________________________________________
Bool_t AlidNdPtAnalysisPbPbAOD::IsHijingParticle(const AliAODMCParticle *part, AliGenHijingEventHeader* hijingGenHeader){
  
  // Check whether a particle is from Hijing or some injected
  // returns kFALSE if particle is injected
  
//  if(part->Label() > (hijingGenHeader->NProduced()-1)) return kFALSE;
//  return kTRUE;
  
  if(part->GetGeneratorIndex()>0) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AlidNdPtAnalysisPbPbAOD::IsPythiaParticle(const AliAODMCParticle *part, AliGenPythiaEventHeader* pythiaGenHeader){
  
  // Check whether a particle is from Pythia or some injected
  
  if(part->Label() > (pythiaGenHeader->NProduced()-1)) return kFALSE;
  return kTRUE;
}

Double_t* AlidNdPtAnalysisPbPbAOD::GetArrayClone(Int_t n, Double_t* source)
{
  if (!source || n==0) return 0;
  Double_t* dest = new Double_t[n];
  for (Int_t i=0; i<n ; i++) { dest[i] = source[i]; }
  return dest;
}

void AlidNdPtAnalysisPbPbAOD::Terminate(Option_t *)
{
  
}


