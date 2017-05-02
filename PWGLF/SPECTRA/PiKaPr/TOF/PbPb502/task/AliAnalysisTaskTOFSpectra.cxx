/**************************************************************************
*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                         *
*  Author: The ALICE Off-line Project.                                    *
*  Contributors are mentioned in the code where appropriate.              *
*                                                                         *
*  Permission to use, copy, modify and distribute this software and its   *
*  documentation strictly for non-commercial purposes is hereby granted   *
*  without fee, provided that the above copyright notice appears in all   *
*  copies and that both the copyright notice and this permission notice   *
*  appear in the supporting documentation. The authors make no claims     *
*  about the suitability of this software for any purpose. It is          *
*  provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified charged hadron spectra: TOF                   //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#define LOG_NO_INFO
#define LOG_NO_DEBUG
#include "AliLog.h"
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3I.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliPID.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisFilter.h"
#include "AliCFContainer.h"
#include "AliCFFrame.h"
#include "AliCFGridSparse.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliMultiplicity.h"
#include "AliAnalysisTaskTOFSpectra.h"
#include "AliESDUtils.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"
#include "AliPPVsMultUtils.h"
#include "AliTOFGeometry.h"
#include "AliCDBManager.h"
#ifdef USETREECLASS
#include "TClonesArray.h"
#include "AliAnTOFtrack.h"
#endif


ClassImp(AliAnalysisTaskTOFSpectra)

//________________________________________________________________________
AliAnalysisTaskTOFSpectra::AliAnalysisTaskTOFSpectra(const TString taskname, Bool_t hi, Bool_t mc, Bool_t tree, Bool_t chan, Bool_t cuts, Int_t simplecuts) :
AliAnalysisTaskSE(taskname),

fEvtMult(-999),
fEvtMultBin(-1),
fVertStatus(-1),
fNContrPrimVertex(-999),

fEvtMask(0),
fTrkMask(0),
fTPCPIDMask(0),
fTrkCutMask(0),
fMCTrkMask(0),
fDCAXYshift(0),
fDCAXY(-999),
fDCAZ(-999),
fRunNumber(0),
fEvtPhysSelected(kFALSE),
fEvtSelected(kFALSE),
fTOFout(kFALSE),
fTRDout(kFALSE),
fTime(kFALSE),
fMismatch(kFALSE),
fPassGoldenChi2(kFALSE),
fPassDCAxy(kFALSE),
fPassDCAz(kFALSE),
fITSTPCMatch(kFALSE),
fT0UsedMask(-999),
fBinPtIndex(-999),
fTPCClusters(0),
fTPCFindClusters(0),
fITSClusters(0),
fTPCCrossedRows(-999),
fTPCChi2(-999),
fTPCChi2PerNDF(-999),
fITSChi2(-999),
fITSChi2PerNDF(-999),
fGoldenChi2(-999),
fLengthActiveZone(-999),
fLength(-999),
fLengthRatio(-999),
fSign(kFALSE),
fTOFTime(-999),
fTOFMismatchTime(-999),
fTOFchan(-999),
fEta(-999),
fPhi(-999),
fPt(-999),
fTOFSignalTot(-999),
fP(-999),
fPTPC(-999),
fTOFImpactDZ(-999),
fTOFImpactDX(-999),
fT0TrkTime(-999),
fNTOFClusters(-999),
fTOFClusters(0x0),
fTPCSignal(-999),
fPhiout(-999),
fXout(-999),
fYout(-999),
fZout(-999),
//MC information
fNMCTracks(-999),
fMCPrimaries(-999),
fMCTOFMatch(-999),
fT0TrkSigma(-999),
fPtMC(-999),
fPMC(-999),
fPhiMC(-999),
fEtaMC(-999),
fProdInfo(-999),
fPdgcode(-999),
fPdgIndex(-999),
fRapidityMC(-999),
fSignMC(kFALSE),

fESD(0x0),
fMCEvt(0x0),
fMCStack(0x0),
fEventCut(0),
fESDtrackCuts(0x0),
fESDtrackCutsPrm(0x0),
fMultSel(0x0),
fTOFcls(0x0),
fTOFcalib(0x0),
fTOFT0maker(0x0),
fTimeResolution(60.),
fListHist(0x0),
fTreeTrack(0x0),
ArrayAnTrk(0x0),
fTreeTrackMC(0x0),

fHImode(hi),
fMCmode(mc),
fTreemode(tree),
fChannelmode(chan),
fCutmode(cuts),
fSimpleCutmode(simplecuts),
fBuilTPCTOF(kFALSE),
fBuilDCAchi2(kFALSE),
fUseTPCShift(kFALSE),
fPerformance(kFALSE),
fRecalibrateTOF(kFALSE),
fFineTOFReso(kFALSE),
fSelectBit(AliVEvent::kINT7),

tb(),//TBenchmark
hPerformanceTime(0x0),
hPerformanceCPUTime(0x0),

fPIDResponse(0x0),
fTOFPIDResponse(),
fESDpid(0x0),
hNEvt(0x0),
hEvtMult(0x0),
hEvtMultAftEvSel(0x0),
hEvtVtxXY(0x0),
hEvtVtxZ(0x0),
hNTrk(0x0),
hPadDist(0x0),
hTOFDist(0x0),
hBeta(0x0),
hBetaNoMismatch(0x0),
hBetaNoMismatchEtaCut(0x0),
hBetaNoMismatchEtaCutOut(0x0),
hBetaCentral(0x0),
hBetaNoMismatchCentral(0x0),
hBetaNoMismatchCentralEtaCut(0x0),
hBetaNoMismatchCentralEtaCutOut(0x0),
hTPCdEdx(0x0),
hCutVariation(0x0),
//TOF Geometrical information
hTOFResidualX(0x0),
hTOFResidualZ(0x0),
hTOFChannel(0x0),
//TOF and T0 distributions for resolution
hT0(0x0),
hT0Resolution(0x0),
hTimeOfFlightRes(0x0),
hTimeOfFlightTOFRes(0x0),
hTimeOfFlightGoodRes(0x0),
hTimeOfFlightResNoMismatch(0x0),
hTimeOfFlightResFine(0x0),
hTimeOfFlightResFinePerEvent(0x0),

fMultiplicityBin(kEvtMultBins+1),
fEtaRange(0.8),
fTOFmin(10000),
fTOFmax(80000),
fLengthmin(350),
fRapidityCut(0.5),
hChannelTime(0x0),
hTOFClusters(0x0),
hTOFClustersDCApass(0x0),
fChannelFirst(0),
fChannelLast(157248)
{
  
  //
  // standard constructur which should be used
  //
  AliInfo("**** CONSTRUCTOR CALLED ****");
  
  //Perform some tests
  //   TestDCAXYBinning();
  
  //Initialize all the variables
  Init();
  
  //Objects for cut variation
  for(Int_t cut = 0; cut < nCutVars; cut++) fCutVar[cut] = 0x0;
  
  //Objects for TOF calibration
  if(fRecalibrateTOF){
    fESDpid = new AliESDpid();
    fTOFcalib = new AliTOFcalib();
    fTOFT0maker = new AliTOFT0maker(fESDpid, fTOFcalib);
  }
  DefineInput(0, TChain::Class());
  DefineAllTheOutput();
  
  PrintStatus();
  AliDebug(2, "**** END OF CONSTRUCTOR ****");
}

//________________________________________________________________________
AliAnalysisTaskTOFSpectra::~AliAnalysisTaskTOFSpectra(){//Destructor
  AliInfo("**** DESTRUCTOR CALLED ****");
  
  if (fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts = 0;
  }
  
  if (fESDtrackCutsPrm) {
    delete fESDtrackCutsPrm;
    fESDtrackCutsPrm = 0;
  }
  
  for(Int_t cut = 0; cut < nCutVars; cut++){
    if(fCutVar[cut]){
      delete fCutVar[cut];
      fCutVar[cut] = 0;
    }
  }
  
  if (fTreeTrack) {
    delete fTreeTrack;
    fTreeTrack = 0;
  }
  
  if (ArrayAnTrk) {
    delete ArrayAnTrk;
    ArrayAnTrk = 0;
  }
  
  if (fTreeTrackMC) {
    delete fTreeTrackMC;
    fTreeTrackMC = 0;
  }
  
  if (fListHist) {
    delete fListHist;
    fListHist = 0;
  }
  
  if (hTimeOfFlightResFinePerEvent) {
    delete hTimeOfFlightResFinePerEvent;
    hTimeOfFlightResFinePerEvent = 0;
  }
  
  if (fTOFT0maker){
    delete fTOFT0maker;
    fTOFT0maker = 0;
  }
  
  if (fTOFcalib){
    delete fTOFcalib;
    fTOFcalib = 0;
  }
  
  if (fESDpid){
    delete fESDpid;
    fESDpid = 0;
  }
  
  AliDebug(2, "**** END OF DESTRUCTOR ****");
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::Init(){//Sets everything to default values
  AliInfo("Init()");
  
  for(Int_t i = 0; i < ntests; i++){
    tr[i] = 0;
    tc[i] = 0;
  }
  
  InitializeEventVar();
  
  //
  //Defining ranges and bins
  //
  AliInfo("Defining Multiplicity bins");
  for(Int_t i = 0; i <=  kEvtMultBins; i++){//Multiplicity
    // Multiplicity [0%, 100%] - not uniform
    //5% step from 0% to 10% -> 2 bins
    //10% step from 10% to 100% -> 9 bins
    //107% step from 100% to 207% -> 1 bins ---> NOT calibrated but put by default in bin 100.5  and also event selection
    //Total of 12 bins
    
    if(i == 0) fMultiplicityBin.AddAt(0., i);//Starting point
    else if(fMultiplicityBin.At(i-1)  < 10.) fMultiplicityBin.AddAt(fMultiplicityBin.At(i-1) + 5., i);
    else if(fMultiplicityBin.At(i-1)  < 100.) fMultiplicityBin.AddAt(fMultiplicityBin.At(i-1) + 10., i);
    else fMultiplicityBin.AddAt(fMultiplicityBin.At(i-1) + 107., i);
    if(i > 0) AliInfo(Form("Multiplicity Bin %i/%i [%.1f, %.1f]", i, kEvtMultBins, fMultiplicityBin.At(i-1), fMultiplicityBin.At(i)));
    
  }
  
  //
  //Shift to the TPC signal
  //
  AliInfo("Defining shifts to TPC signals");
  //Pass 1
  // const Double_t PionShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.495301, -2.849033, -3.038054, -3.106598, -3.105506, -3.075790, -2.991061, -2.913724, -2.852930, -2.803524, -2.757924, -2.696637, -2.632819, -2.568329, -2.496323, -2.390133, -2.274463, -2.166133, -2.063583, -1.965722, -1.870548, -1.775581, -1.687634, -1.606076, -1.533710, -1.469498, -1.421488, -1.389162, -1.367918, -1.350301, -1.337053, -1.326001, -1.310669, -1.291229, -1.269140, -1.226947, -1.172763, -1.103424, -1.039604, -0.978991, -0.925905, -0.879515, -0.837937, -0.803185, -0.765683};
  // const Double_t KaonShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 4.440578, 4.016308, 5.250707, 11.406953, 12.859423, 11.364870, 9.495074, 8.067584, 6.774985, 5.729376, 4.907116, 4.271800, 3.750499, 3.306340, 2.933824, 2.505940, 2.039071, 1.683853, 1.415126, 1.213763, 1.070441, 0.979520, 0.930015, 0.924567, 0.955567, 1.019278, 1.100308, 1.188906, 1.289394, 1.391640, 1.483409, 1.575449, 1.650027, 1.715081, 1.768910, 1.828771, 1.867502, 1.879254, 1.864644, 1.825089, 1.792100, 1.760907, 1.742615, 1.731851, 1.736609};
  // const Double_t ProtonShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 3.548278, 3.338830, 3.489038, 3.744695, 4.090173, 4.169596, 4.311995, 7.942016, 13.574720, 14.799053, 14.592191, 14.012937, 13.128907, 12.020871, 10.796905, 9.115128, 7.384271, 6.039855, 4.986625, 4.142134, 3.463985, 2.917085, 2.466358, 2.089473, 1.775816, 1.520427, 1.300421, 1.112384, 0.950658, 0.817063, 0.699761, 0.605128, 0.516923, 0.452642, 0.390945, 0.330683, 0.293475, 0.310124, 0.371988, 0.475335, 0.621825, 0.807812, 0.984995, 1.186980, 1.354701};
  
  //Pass 2
  const Double_t PionShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -4.106467, -4.416255, -4.686875, -4.813587, -4.881268, -4.895112, -4.882032, -4.849220, -4.776634, -4.703480, -4.647757, -4.600764, -4.552614, -4.507683, -4.472258, -4.416541, -4.333280, -4.251431, -4.173233, -4.098128, -4.026338, -3.958499, -3.896338, -3.839612, -3.792598, -3.757292, -3.734550, -3.724568, -3.724493, -3.727113, -3.729910, -3.733357, -3.734002, -3.726840, -3.718131, -3.691936, -3.653432, -3.603556, -3.564313, -3.526249, -3.498696, -3.474054, -3.453687, -3.430655, -3.406596};
  const Double_t KaonShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.534201, 1.272898, 3.331416, 8.894055, 8.481794, 6.733964, 5.201038, 4.018954, 3.010390, 2.192596, 1.534934, 1.020067, 0.594652, 0.240203, -0.058884, -0.425643, -0.804173, -1.088455, -1.303637, -1.462859, -1.576417, -1.647016, -1.683909, -1.685685, -1.656390, -1.609954, -1.545773, -1.474710, -1.402276, -1.334256, -1.273858, -1.221346, -1.180201, -1.146672, -1.122532, -1.100815, -1.096768, -1.114986, -1.155814, -1.205518, -1.247801, -1.278684, -1.293594, -1.296794, -1.295287};
  const Double_t ProtonShift[kPtBins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.130821, 0.914432, 0.923787, 1.031415, 1.133344, 1.489368, 6.903029, 12.047407, 13.346370, 13.112204, 11.835434, 10.205098, 8.706627, 7.433365, 6.318621, 4.920023, 3.516806, 2.428509, 1.575802, 0.895684, 0.350548, -0.095673, -0.459524, -0.764463, -1.022071, -1.241166, -1.424443, -1.586363, -1.725276, -1.842191, -1.951553, -2.040862, -2.122386, -2.185230, -2.234066, -2.292148, -2.327972, -2.315525, -2.264726, -2.172619, -2.061416, -1.917858, -1.767914, -1.623935, -1.490772};
  
  for(Int_t species = 0; species <  3; species++){//Particle loop
    for(Int_t ptbin = 0; ptbin <  kPtBins; ptbin++){//Pt loop
      switch (species) {
        case 0:
        fTPCShift[species][ptbin] = fUseTPCShift ? PionShift[ptbin] : 0;
        break;
        case 1:
        fTPCShift[species][ptbin] = fUseTPCShift ? KaonShift[ptbin] : 0;
        break;
        case 2:
        fTPCShift[species][ptbin] = fUseTPCShift ? ProtonShift[ptbin] : 0;
        break;
        default:
        break;
      }
      AliInfo(Form("TPC signal shift for species %i ptbin %i is %f", species, ptbin, fTPCShift[species][ptbin]));
      
    }
  }
  
  for(Int_t i = 0; i < kExpSpecies;i++) fTOFPIDProbability[i] = -999;
  
  for(Int_t i = 0; i < 3; i++) fRapidity[i] = -999;
  
  //
  //Histograms
  //
  for(Int_t i = 0; i < 2; i++){//After and before the cut loop
    hTrkTPCCls[i] = 0x0;
    hTrkTPCRows[i] = 0x0;
    hTrkTPCRatioRowsFindCls[i] = 0x0;
    hTrkTPCChi2NDF[i] = 0x0;
    hTrkITSChi2NDF[i] = 0x0;
    hTrkActiveLength[i] = 0x0;
    hTrkITSTPCMatch[i] = 0x0;
    hTrkDCAxy[i] = 0x0;
    hTrkDCAz[i] = 0x0;
    
    #ifdef CHECKTRACKCUTS
    for(Int_t j = 0; j < 3; j++){//Loop on the correlatio variable, i.e. pt [0] - eta [1] - phi [2]
      hTrkTPCClsCorr[i][j] = 0x0;
      hTrkTPCRowsCorr[i][j] = 0x0;
      hTrkTPCRatioRowsFindClsCorr[i][j] = 0x0;
      hTrkTPCChi2NDFCorr[i][j] = 0x0;
      hTrkITSChi2NDFCorr[i][j] = 0x0;
      hTrkActiveLengthCorr[i][j] = 0x0;
      hTrkITSTPCMatchCorr[i][j] = 0x0;
      hTrkDCAxyCorr[i][j] = 0x0;
      hTrkDCAzCorr[i][j] = 0x0;
      
    }
    #endif
    
  }
  
  
  for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
    for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
      for(Int_t species = 0; species < 3; species++){//Species loop
        
        hNumMatchMultTrk[charge][species][mult] = 0x0;
        hDenMatchMultTrk[charge][species][mult] = 0x0;
        
        hNumMatchMultTrkTRDOut[charge][species][mult] = 0x0;
        hDenMatchMultTrkTRDOut[charge][species][mult] = 0x0;
        
        hNumMatchMultTrkNoTRDOut[charge][species][mult] = 0x0;
        hDenMatchMultTrkNoTRDOut[charge][species][mult] = 0x0;
        
      }
      
      hNumMatchMultTrkInc[charge][mult] = 0x0;
      hDenMatchMultTrkInc[charge][mult] = 0x0;
      
    }
    
  }
  
  for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
    for(Int_t species = 0; species < 3; species++){//Species loop
      
      hDenTrkVertMultTrk[charge][species] = 0x0;
      hDenTrkTriggerMultTrk[charge][species] = 0x0;
      
      for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
        hDenPrimMCYCut[charge][species][mult] = 0x0;
        hDenPrimMCEtaCut[charge][species][mult] = 0x0;
        hDenPrimMCEtaYCut[charge][species][mult] = 0x0;
        hNumPrimMCTrueMatch[charge][species][mult] = 0x0;
        hNumPrimMCTrueMatchYCut[charge][species][mult] = 0x0;
        hNumPrimMCTrueMatchYCutTPC[charge][species][mult] = 0x0;
        hNumPrimMCConsistentMatchYCut[charge][species][mult] = 0x0;
      }
      
    }
  }
  
  for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
    hNumMatch[charge] = 0x0;
    hDenMatch[charge] = 0x0;
    
    hNumMatchPtEtaPhiout[charge] = 0x0;
    hDenMatchPtEtaPhiout[charge] = 0x0;
    
    hNumMatchEta[charge] = 0x0;
    hDenMatchEta[charge] = 0x0;
    
    hNumMatchphiOut[charge] = 0x0;
    hDenMatchphiOut[charge] = 0x0;
    
    hNumMatchEtaPtMa[charge] = 0x0;
    hDenMatchEtaPtMa[charge] = 0x0;
    
    hNumMatchphiOutPtMa[charge] = 0x0;
    hDenMatchphiOutPtMa[charge] = 0x0;
    
    hNumMatchTRDOut[charge] = 0x0;
    hDenMatchTRDOut[charge] = 0x0;
    
    hNumMatchNoTRDOut[charge] = 0x0;
    hDenMatchNoTRDOut[charge] = 0x0;
    
    for(Int_t species = 0; species < 3; species++){//Species loop
      hNumMatchTPC[charge][species] = 0x0;
      hDenMatchTPC[charge][species] = 0x0;
      
      hNumMatchMC[charge][species] = 0x0;
      hDenMatchMC[charge][species] = 0x0;
      
      hNumMatchPrimMC[charge][species] = 0x0;
      hDenMatchPrimMC[charge][species] = 0x0;
      for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
        
        hNumMatchPrimMCYCut[charge][species][mult] = 0x0;
        hDenMatchPrimMCYCut[charge][species][mult] = 0x0;
      }
    }
  }
  
  for(Int_t species = 0; species < kExpSpecies; species++){//Species loop
    fTOFExpSigma[species] = -999;
    fTOFExpTime[species] = -999;
    fTOFSigma[species] = -999;
    fTPCSigma[species] = -999;
    if(species < 3) fCombinedSigma[species] = -999;
    
    #ifdef CHECKCOMPUTEDVALUES
    hTOFExpectedComputed[species] = 0x0;
    hTOFExpectedComputedTPC[species] = 0x0;
    hTOFMomComputed[species] = 0x0;
    hTOFMomComputedTPC[species] = 0x0;
    #endif
    
    for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
      hTPCTOFSeparation[species][ptbin] = 0x0;
    }
  }
  
  #ifdef BUILDTOFDISTRIBUTIONS// Build TOF distributions only if requested
  for(Int_t species = 0; species < kSpecies; species++){
    for(Int_t charge = 0; charge < kCharges; charge++){
      for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
        hTOF[ptbin][charge][species] = 0x0;
        hTOFNoYCut[ptbin][charge][species] = 0x0;
        hTOFSigma[ptbin][charge][species] = 0x0;
        hTOFNoMismatch[ptbin][charge][species] = 0x0;
        hTOFSigmaNoMismatch[ptbin][charge][species] = 0x0;
      }
    }
  }
  #endif
  
  
  for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
    for(Int_t species = 0; species < 3; species++){//Species loop
      for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
        for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
          hDCAxy[charge][species][ptbin][mult] = 0x0;
          hDCAxyGoldenChi2[charge][species][ptbin][mult] = 0x0;
        }
        hDCAxyPrimMC[charge][species][ptbin] = 0x0;
        hDCAxySecStMC[charge][species][ptbin] = 0x0;
        hDCAxySecMatMC[charge][species][ptbin] = 0x0;
      }
    }
  }
  
  AliDebug(2, "Init()\t END");
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::UserCreateOutputObjects(){
  
  AliInfo("UserCreateOutputObjects");
  
  PrintStatus();
  SetTrackCuts();//Standard cuts
  if(fTreemode){
    SetCutVar();//Cut variation for tree
    PrintCutVariablesForTree();
  }
  SetSimpleCutVar();//Simple cut variation
  PrintCutVariables();
  
  //Define the DCA histogram binning
  Double_t DCAXYbin[fDCAXYbins+1];
  DCAXYbin[0] = -fDCAXYRange;
  for(Int_t i = 0; i < fDCAXYbins; i++){
    if(DCAXYbin[i] > -.1 && DCAXYbin[i] < .1) DCAXYbin[i+1] = DCAXYbin[i] + 0.0008;
    else DCAXYbin[i+1] = DCAXYbin[i] + 0.007733333333333333333333;
    // cout<<i<<"  ["<<DCAXYbin[i]<<","<<DCAXYbin[i+1]<<"] width = "<< DCAXYbin[i]-DCAXYbin[i+1]<<endl;
  }
  
  
  // Create histograms
  // Called once
  //List of the output histograms
  
  OpenFile(1);
  
  if (!fListHist) fListHist = new TList();
  fListHist->SetOwner();
  
  //Object for the event selection
  if(fHImode) fEventCut.SetupLHC15o();
  else{
    fEventCut.SetupRun2pp();
    fEventCut.fCentralityFramework = 0;
  }
  fEventCut.SetManualMode();
  
  //Histograms for event Selection
  fEventCut.AddQAplotsToList(fListHist);
  
  DefineTimePerformance();
  
  if(fChannelmode){
    hChannelTime = new TH2I("hChannelTime", "Histogram with Channel/Time correlation;Channel;TOF (ps)", kChannelBins, fChannelFirst -.5, fChannelLast -.5, kTimeBins, fTOFmin, fTOFmax);
    fListHist->AddLast(hChannelTime);
  }
  else{
    
    Int_t binstart = 1;
    
    hNEvt = new TH1F("hNEvt", "Number of processed events;Evt. Sel. Step;Counts", 15, -.5, 14.5);
    hNEvt->Sumw2();
    hNEvt->SetMinimum(0);
    hNEvt->GetXaxis()->SetBinLabel(binstart++, "Read from ESD");
    hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has AliESDtrackCuts");
    hNEvt->GetXaxis()->SetBinLabel(binstart++, "Pass Phys. Sel. + Trig");
    if(fHImode){
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has AliMultSelection");//Multiplicity estimator initialized
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has Calibrated Mult.");//kNoCalib: centrality not calibrated, this is the default value for the centrality code
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Pass Trigger");//kRejTrigger: do not pass the trigger
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has INEL>0");//kRejINELgtZERO: do not pass INEL>0 Cut
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Pass Vtx Cut");//pkRejVzCut: do not pass vertex Cut
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Pass Pile-up");//kRejPileupInMultBins: do not pass Pile-up Cut
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has consistent vertex");//kRejConsistencySPDandTrackVertices: do not pass consistency of vertex Cut
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "pass Trk.lets Vs Clusters");// kRejTrackletsVsClusters: do not pass Tracklets Vs Clusters Cut
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "Has Vertex Contributors");//kRejNonZeroNContribs: do not pass Contributors (to vertex) Cut
    }
    else{
      
      //AliPPVsMultUtils::IsMinimumBias(fESD))
      //AliPPVsMultUtils::IsAcceptedVertexPosition(fESD))
      //AliPPVsMultUtils::IsINELgtZERO(fESD))
      //AliPPVsMultUtils::IsNotPileupSPDInMultBins(fESD))
      //AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fESD))
      
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "IsMinimumBias");//AliPPVsMultUtils::IsMinimumBias(fESD))
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "IsAcceptedVertexPosition");//AliPPVsMultUtils::IsAcceptedVertexPosition(fESD))
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "IsINELgtZERO");//AliPPVsMultUtils::IsINELgtZERO(fESD))
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "IsNotPileupSPDInMultBins");//AliPPVsMultUtils::IsNotPileupSPDInMultBins(fESD))
      hNEvt->GetXaxis()->SetBinLabel(binstart++, "HasNoInconsistentSPDandTrackVertices");//AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fESD))
      
    }
    
    if(binstart > hNEvt->GetNbinsX() + 1) AliFatal(Form("binstart out of bounds!!"));
    fListHist->AddLast(hNEvt);
    
    hEvtMult = new TH1F("hEvtMult", "Event Multiplicity Before Event Selection;Multiplicity;Counts", 6002, fHImode ? -1. : -301, fHImode ? 6001. : 5701);
    if(hEvtMult->GetXaxis()->GetBinWidth(100) != 1.) AliWarning(Form("Bins have size %f which is different from one!", hEvtMult->GetXaxis()->GetBinWidth(100)));
    if(hEvtMult->GetXaxis()->GetBinCenter(hEvtMult->GetXaxis()->FindBin(-0.5)) != -0.5) AliWarning(Form("First bin has center in %f which is different from -0.5!", hEvtMult->GetXaxis()->GetBinCenter(1)));
    fListHist->AddLast(hEvtMult);
    
    hEvtMultAftEvSel = new TH1F("hEvtMultAftEvSel", "Event Multiplicity After Event Selection;Multiplicity;Counts", 6002, fHImode ? -1. : -301, fHImode ? 6001. : 5701);
    if(hEvtMultAftEvSel->GetXaxis()->GetBinWidth(100) != 1.) AliWarning(Form("Bins have size %f which is different from one!", hEvtMultAftEvSel->GetXaxis()->GetBinWidth(100)));
    if(hEvtMultAftEvSel->GetXaxis()->GetBinCenter(hEvtMultAftEvSel->GetXaxis()->FindBin(-0.5)) != -0.5) AliWarning(Form("First bin has center in %f which is different from -0.5!", hEvtMultAftEvSel->GetXaxis()->GetBinCenter(1)));
    fListHist->AddLast(hEvtMultAftEvSel);
    
    hEvtVtxXY = new TH1F("hEvtVtxXY", "XY primary vertex distance;(x^2+y^2)^(1/2) (cm);Counts", 100, -5., 5.);
    fListHist->AddLast(hEvtVtxXY);
    
    hEvtVtxZ = new TH1F("hEvtVtxZ", "Z primary vertex distance;Z (cm);Counts", 100, -50., 50.);
    fListHist->AddLast(hEvtVtxZ);
    
    binstart = 1;
    
    hNTrk = new TH1F("hNTrk", "Number of processed tracks;Trk. Sel. Step;Counts", 20, -0.5, 19.5);
    hNTrk->Sumw2();
    hNTrk->SetMinimum(0);
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Reconstructed tracks");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Pass Trk Cuts");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Tracks in |#eta| < %.2f", fEtaRange));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Non-zero charge");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "In pt range");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Pass DCAxy cut");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Has kTOFout");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Has kTIME");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Has kTRDout");
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fLength>%.0f", fLengthmin));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fTOFTime>%.0f", fTOFmin));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fTOFExpTime[AliPID::kPion]>%.0f", fTOFmin));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fTOFExpTime[AliPID::kKaon]>%.0f", fTOFmin));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fTOFExpTime[AliPID::kProton]>1%.0f", fTOFmin));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, Form("Has fTOFTime<%.0f", fTOFmax));
    hNTrk->GetXaxis()->SetBinLabel(binstart++, "Is selected");
    
    if(binstart > hNTrk->GetNbinsX() + 1) AliFatal(Form("binstart out of bounds!!"));
    fListHist->AddLast(hNTrk);
    
    
    //
    //Histograms for Cut values
    //
    const TString passtitle[2] = {"Not Checked", "Checked"};
    const TString passname[2] = {"_NoCheck", "_Check"};
    const Int_t nvar = 9;
    const Int_t varbins[nvar] =   { 200,  200, 200, 100, 100,  100,   2, 400, 400}; //TPCCls - TPCCrR - Ratio TPC crossed rows TPC findable clusters - TPC Chi2NDF - ITS Chi2NDF - Length in active region  - ITS/TPC match - DCAxy - DCAz
    const Double_t varmin[nvar] = {  0.,   0., 0.,  0.,  0.,   0., -.5,  -4., -4.};
    const Double_t varmax[nvar] = {200., 200., 2., 10., 50., 200., 1.5,   4.,  4.};
    for(Int_t i = 0; i < 2; i++){
      Int_t index = 0;
      hTrkTPCCls[i] = new TH1F(Form("hTrkTPCCls%s", passname[i].Data()) , Form("TPC Clusters %s (%i);TPC Clusters;Counts", passtitle[i].Data(), fESDtrackCuts->GetMinNClusterTPC()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkTPCCls[i]);
      
      index++;
      hTrkTPCRows[i] = new TH1F(Form("hTrkTPCRows%s", passname[i].Data()) , Form("TPC Rows %s (%f);TPC Crossed Rows;Counts", passtitle[i].Data(), fESDtrackCuts->GetMinNCrossedRowsTPC()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkTPCRows[i]);
      
      index++;
      hTrkTPCRatioRowsFindCls[i] = new TH1F(Form("hTrkTPCRatioRowsFindCls%s", passname[i].Data()) , Form("TPC Rows/Findable clusters %s (%f);TPC Crossed Rows/TPC Findable clusters;Counts", passtitle[i].Data(), fESDtrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkTPCRatioRowsFindCls[i]);
      
      index++;
      hTrkTPCChi2NDF[i] = new TH1F(Form("hTrkTPCChi2NDF%s", passname[i].Data()) , Form("TPC Chi2/NDF %s (%f);TPC Chi2/NDF;Counts", passtitle[i].Data(), fESDtrackCuts->GetMaxChi2PerClusterTPC()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkTPCChi2NDF[i]);
      
      index++;
      hTrkITSChi2NDF[i] = new TH1F(Form("hTrkITSChi2NDF%s", passname[i].Data()) , Form("ITS Chi2/NDF %s (%f);ITS Chi2/NDF;Counts", passtitle[i].Data(), fESDtrackCuts->GetMaxChi2PerClusterITS()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkITSChi2NDF[i]);
      
      index++;
      hTrkActiveLength[i] = new TH1F(Form("hTrkActiveLength%s", passname[i].Data()) , Form("Active Length %s (%f);Length in active region;Counts", passtitle[i].Data(), fESDtrackCuts->GetMinLengthActiveVolumeTPC()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkActiveLength[i]);
      
      index++;
      hTrkITSTPCMatch[i] = new TH1F(Form("hTrkITSTPCMatch%s", passname[i].Data()) , Form("ITS TPC match %s (%i);ITS/TPC match;Counts", passtitle[i].Data(), fESDtrackCuts->GetRequireITSRefit()), varbins[index], varmin[index], varmax[index]);
      hTrkITSTPCMatch[i]->GetXaxis()->SetBinLabel(1, "NoMatch");
      hTrkITSTPCMatch[i]->GetXaxis()->SetBinLabel(2, "Match");
      fListHist->AddLast(hTrkITSTPCMatch[i]);
      
      index++;
      hTrkDCAxy[i] = new TH1F(Form("hTrkDCAxy%s", passname[i].Data()) , Form("DCAxy %s (%s);DCAxy (cm);Counts", passtitle[i].Data(), fESDtrackCuts->GetMaxDCAToVertexXYPtDep()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkDCAxy[i]);
      
      index++;
      hTrkDCAz[i] = new TH1F(Form("hTrkDCAz%s", passname[i].Data()) , Form("DCAz %s (%f);DCAz (cm);Counts", passtitle[i].Data(), fESDtrackCuts->GetMaxDCAToVertexZ()), varbins[index], varmin[index], varmax[index]);
      fListHist->AddLast(hTrkDCAz[i]);
      
      if(index != nvar - 1) AliFatal("Index is different than designed");
      
      #ifdef CHECKTRACKCUTS
      for(Int_t j = 0; j < 3; j++){
        TString hname = "";
        Int_t bins;
        Double_t limits[2];
        switch (j) {
          case 0:
          hname = "Pt";
          bins = 3*kPtBins;
          limits[0] = 0;
          limits[1] = fBinPt[kPtBins];
          break;
          case 1:
          hname = "Eta";
          bins = 20;
          limits[0] = -1;
          limits[1] =  1;
          break;
          case 2:
          hname = "Phi";
          bins = 20;
          limits[0] = 0;
          limits[1] = TMath::TwoPi();
          break;
          default:
          AliFatal("index out of bound!");
          break;
        }
        
        index = 0;
        hTrkTPCClsCorr[i][j] = new TH2I(Form("hTrkTPCClsCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkTPCCls[i]->GetName(), hTrkTPCCls[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkTPCClsCorr[i][j]);
        
        index++;
        hTrkTPCRowsCorr[i][j] = new TH2I(Form("hTrkTPCRowsCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkTPCRows[i]->GetName(), hTrkTPCRows[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkTPCRowsCorr[i][j]);
        
        index++;
        hTrkTPCRatioRowsFindClsCorr[i][j] = new TH2I(Form("hTrkTPCRatioRowsFindClsCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkTPCRatioRowsFindCls[i]->GetName(), hTrkTPCRatioRowsFindCls[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkTPCRatioRowsFindClsCorr[i][j]);
        
        index++;
        hTrkTPCChi2NDFCorr[i][j] = new TH2I(Form("hTrkTPCChi2NDFCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkTPCChi2NDF[i]->GetName(), hTrkTPCChi2NDF[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkTPCChi2NDFCorr[i][j]);
        
        index++;
        hTrkITSChi2NDFCorr[i][j] = new TH2I(Form("hTrkITSChi2NDFCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkITSChi2NDF[i]->GetName(), hTrkITSChi2NDF[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkITSChi2NDFCorr[i][j]);
        
        index++;
        hTrkActiveLengthCorr[i][j] = new TH2I(Form("hTrkActiveLengthCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkActiveLength[i]->GetName(), hTrkActiveLength[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkActiveLengthCorr[i][j]);
        
        index++;
        hTrkITSTPCMatchCorr[i][j] = new TH2I(Form("hTrkITSTPCMatchCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkITSTPCMatch[i]->GetName(), hTrkITSTPCMatch[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        hTrkITSTPCMatchCorr[i][j]->GetXaxis()->SetBinLabel(1, "NoMatch");
        hTrkITSTPCMatchCorr[i][j]->GetXaxis()->SetBinLabel(2, "Match");
        fListHist->AddLast(hTrkITSTPCMatchCorr[i][j]);
        
        index++;
        hTrkDCAxyCorr[i][j] = new TH2I(Form("hTrkDCAxyCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkDCAxy[i]->GetName(), hTrkDCAxy[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkDCAxyCorr[i][j]);
        
        index++;
        hTrkDCAzCorr[i][j] = new TH2I(Form("hTrkDCAzCorr%s%s", hname.Data(), passname[i].Data()) , Form("%s;%s;%s", hTrkDCAz[i]->GetName(), hTrkDCAz[i]->GetXaxis()->GetTitle(), hname.Data()), varbins[index], varmin[index], varmax[index], bins, limits[0], limits[1]);
        fListHist->AddLast(hTrkDCAzCorr[i][j]);
        
        if(index != nvar - 1) AliFatal("Index is different than designed");
        
      }
      
      #endif
      
    }
    
    hPadDist = new TH2F("hPadDist", "Distribution in the pads;#DeltaX_{pad} (cm);#DeltaZ_{pad} (cm)", 400, -10, 10, 400, -10, 10);
    fListHist->AddLast(hPadDist);
    
    hTOFDist = new TH2F("hTOFDist", "Distribution in the TOF;Sector;Strip", 18, 0, 18, 91, 0, 91);
    fListHist->AddLast(hTOFDist);
    
    if(fPerformance){
      const Int_t Bnbins      = 4000;
      const Double_t Blim[2]  = {0., 1.5};
      const Double_t Bplim[2] = {0., 10.};
      
      hBeta = new TH2I("hBeta", Form("Distribution of the beta;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBeta);
      
      hBetaNoMismatch = new TH2I("hBetaNoMismatch", Form("Distribution of the beta w/o Mismatch;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatch);
      
      hBetaNoMismatchEtaCut = new TH2I("hBetaNoMismatchEtaCut", Form("Distribution of the beta w/o Mismatch and a |#eta| < 0.5;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatchEtaCut);
      
      hBetaNoMismatchEtaCutOut = new TH2I("hBetaNoMismatchEtaCutOut", Form("Distribution of the beta w/o Mismatch and a |#eta| > 0.2;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatchEtaCutOut);
      
      hBetaCentral = new TH2I("hBetaCentral", Form("Distribution of the beta Central Events;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaCentral);
      
      hBetaNoMismatchCentral = new TH2I("hBetaNoMismatchCentral", Form("Distribution of the beta w/o Mismatch Central Events;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatchCentral);
      
      hBetaNoMismatchCentralEtaCut = new TH2I("hBetaNoMismatchCentralEtaCut", Form("Distribution of the beta w/o Mismatch Central Events and a |#eta| < 0.5;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatchCentralEtaCut);
      
      hBetaNoMismatchCentralEtaCutOut = new TH2I("hBetaNoMismatchCentralEtaCutOut", Form("Distribution of the beta w/o Mismatch Central Events and a |#eta| > 0.2;%s;TOF #beta", pstring.Data()), Bnbins, Bplim[0], Bplim[1], Bnbins, Blim[0], Blim[1]);
      fListHist->AddLast(hBetaNoMismatchCentralEtaCutOut);
      
      hTPCdEdx = new TH2I("hTPCdEdx", Form("Distribution of the TPC dE/dx;%s;d#it{E}/d#it{x} in TPC (arb. units)", pstring.Data()), 1000, 0.25, 30., 1000, 0., 1000);
      fListHist->AddLast(hTPCdEdx);
    }
    
    hCutVariation = new TH1F("hCutVariation", "CutCounter;CutSet;Tracks", 20, 0. -.5, 20. -.5);
    fListHist->AddLast(hCutVariation);
    
    hTOFResidualX = new TH1F("hTOFResidualX", "X Residual;#DeltaX_{pad} (cm)", 150, -10, 10);
    fListHist->AddLast(hTOFResidualX);
    
    hTOFResidualZ = new TH1F("hTOFResidualZ", "Z Residual;#DeltaZ_{pad} (cm)", 150, -10, 10);
    fListHist->AddLast(hTOFResidualZ);
    
    hTOFChannel = new TH1F("hTOFChannel", "Channel;Channel;Counts", 500, 0., 170000);
    fListHist->AddLast(hTOFChannel);
    
    hT0 = new TH1F("hT0", "hT0;T0", 500, -250, 250);
    fListHist->AddLast(hT0);
    
    hT0Resolution = new TH1F("hT0Resolution", "T0Resolution;T0 #sigma (ps)", 250, 0, 250);
    fListHist->AddLast(hT0Resolution);
    
    hTimeOfFlightRes = new TH1F("hTimeOfFlightRes", "TOF Resolution in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500);
    fListHist->AddLast(hTimeOfFlightRes);
    
    hTimeOfFlightTOFRes = new TH1F("hTimeOfFlightTOFRes", "TOF Resolution with TOF T0 in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500);
    fListHist->AddLast(hTimeOfFlightTOFRes);
    
    hTimeOfFlightGoodRes = new TH1F("hTimeOfFlightGoodRes", "TOF Resolution Good in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500);
    hTimeOfFlightGoodRes->Sumw2();
    fListHist->AddLast(hTimeOfFlightGoodRes);
    
    hTimeOfFlightResNoMismatch = new TH1F("hTimeOfFlightResNoMismatch", "TOF Resolution Wo Mismatch in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 800, -4000, 4000);
    fListHist->AddLast(hTimeOfFlightResNoMismatch);
    
    if(fFineTOFReso){
      hTimeOfFlightResFine = new TH2F("hTimeOfFlightResFine", "TOF Resolution in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500, 100, 0, 100);
      fListHist->AddLast(hTimeOfFlightResFine);
      
      hTimeOfFlightResFinePerEvent = new TH1F("hTimeOfFlightResFinePerEvent", "", 100, -500, 500);//This histogram is not to be added as it is used only to store the information for each event, be sure to have it in the destructor!!!
      
    }
    
    hTOFClusters = new TH1F("hTOFClusters", "Number of clusters per track;TOF Clusters;Number of tracks", 50, -.5, 49.5);
    fListHist->AddLast(hTOFClusters);
    
    hTOFClustersDCApass = new TH1F("hTOFClustersDCApass", "Number of clusters per track that passed DCA cut;TOF Clusters;Number of tracks", 50, -.5, 49.5);
    fListHist->AddLast(hTOFClustersDCApass);
    
    for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
      const TString ptinterval =  Form("in pt bin %i [%.2f,%.2f] %s", ptbin, fBinPt[ptbin], fBinPt[ptbin+1], ptstringOnly.Data());
      if(fBuilTPCTOF){
        for(Int_t species = 0; species < kExpSpecies; species++){//Species Loop
          hTPCTOFSeparation[species][ptbin] = new TH2F(Form("hTPCTOFSeparation%s%i", pS_all[species].Data(), ptbin), Form("TPC TOF PID separation for %s %s;TOF separation (n#sigma);TPC separation (n#sigma)", pS_all[species].Data(), ptinterval.Data()), 200, -20.10, 19.90, 200, -20.10, 19.90);
          fListHist->AddLast(hTPCTOFSeparation[species][ptbin]);
        }
      }
      
      #ifdef BUILDTOFDISTRIBUTIONS// Build TOF distributions only if requested
      for(Int_t species = 0; species < kSpecies; species++){//Species Loop
        for(Int_t charge = 0; charge < kCharges; charge++){//Charge Loop
          const TString hname = Form("%s%s_%i", pC[charge].Data(), pS[species].Data(), ptbin);
          hTOF[ptbin][charge][species] = new TH1F(Form("hTOF%s", hname.Data()), Form("TOF %s%s %s;%s;Counts", pC[charge].Data(), pS[species].Data(), ptinterval.Data(), tofsignalstringSpecies[2+species].Data()), 1000, -10000, 10000);
          fListHist->AddLast(hTOF[ptbin][charge][species]);
          
          hTOFNoYCut[ptbin][charge][species] = new TH1F(Form("hTOFNoYCut%s", hname.Data()), Form("TOF w/o Y cut %s%s %s;%s;Counts", pC[charge].Data(), pS[species].Data(), ptinterval.Data(), tofsignalstringSpecies[2+species].Data()), 1000, -10000, 10000);
          fListHist->AddLast(hTOFNoYCut[ptbin][charge][species]);
          
          hTOFSigma[ptbin][charge][species] = new TH1F(Form("hTOFSigma%s", hname.Data()), Form("TOF Sigma %s%s %s;%s;Counts", pC[charge].Data(), pS[species].Data(), ptinterval.Data(), nsigmastringSpecies[2+species].Data()), 1000, -300, 300);
          fListHist->AddLast(hTOFSigma[ptbin][charge][species]);
          
          hTOFNoMismatch[ptbin][charge][species] = new TH1F(Form("hTOFNoMismatch%s", hname.Data()), Form("TOF No Mism. %s%s %s;%s;Counts", pC[charge].Data(), pS[species].Data(), ptinterval.Data(), tofsignalstringSpecies[2+species].Data()), 1000, -10000, 10000);
          fListHist->AddLast(hTOFNoMismatch[ptbin][charge][species]);
          
          hTOFSigmaNoMismatch[ptbin][charge][species] = new TH1F(Form("hTOFSigmaNoMismatch%s", hname.Data()), Form("TOF Sigma No Mism. %s%s %s;%s;Counts", pC[charge].Data(), pS[species].Data(), ptinterval.Data(), nsigmastringSpecies[2+species].Data()), 1000, -300, 300);
          fListHist->AddLast(hTOFSigmaNoMismatch[ptbin][charge][species]);
        }
      }
      #endif
      
    }
    
    #ifdef CHECKCOMPUTEDVALUES
    for(Int_t species = 0; species < kExpSpecies; species++){//Species Loop
      hTOFExpectedComputed[species] = new TH1F(Form("hTOFExpectedComputed%s", AliPID::ParticleShortName(species)), Form("T_{exp %s calc.}/T_{exp %s};T_{exp %s calc.}/T_{exp %s};Counts", AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species)), 2000, 0, 2);
      fListHist->AddLast(hTOFExpectedComputed[species]);
      
      hTOFExpectedComputedTPC[species] = new TH1F(Form("hTOFExpectedComputedTPC%s", AliPID::ParticleShortName(species)), Form("T_{exp %s calc. TPC}-T_{exp %s};T_{exp %s calc. TPC}-T_{exp %s};Counts", AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species)), 2000, 0, 2);
      fListHist->AddLast(hTOFExpectedComputedTPC[species]);
      
      hTOFMomComputed[species] = new TH1F(Form("hTOFMomComputed%s", AliPID::ParticleShortName(species)), Form("P_{exp %s calc.}/P;P_{exp %s calc.}/P;Counts", AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species)), 2000, 0, 2);
      fListHist->AddLast(hTOFMomComputed[species]);
      
      hTOFMomComputedTPC[species] = new TH1F(Form("hTOFMomComputedTPC%s", AliPID::ParticleShortName(species)), Form("P_{exp %s calc.}/P_{TPC};P_{exp %s calc.}/P_{TPC};Counts", AliPID::ParticleLatexName(species), AliPID::ParticleLatexName(species)), 2000, 0, 2);
      fListHist->AddLast(hTOFMomComputedTPC[species]);
      
    }
    #endif
    
    for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
      //*****
      //PT
      //*****
      
      //Numerator
      hNumMatch[charge] = new TH1F(Form("hNumMatch%s", pC[charge].Data()), Form("Matching Numerator %s;%s", pCharge[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
      hNumMatch[charge]->Sumw2();
      fListHist->AddLast(hNumMatch[charge]);
      
      //Denominator
      hDenMatch[charge] = new TH1F(Form("hDenMatch%s", pC[charge].Data()), Form("Matching Denominator %s;%s", pCharge[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
      hDenMatch[charge]->Sumw2();
      fListHist->AddLast(hDenMatch[charge]);
      
      //*****
      //ETA
      //*****
      
      //Numerator
      hNumMatchEta[charge] = new TH1F(Form("hNumMatch%sEta", pC[charge].Data()), Form("Matching Numerator in #eta %s;%s", pCharge[charge].Data(), etastring.Data()), kEtaBins, -fEtaRange, fEtaRange);
      hNumMatchEta[charge]->Sumw2();
      fListHist->AddLast(hNumMatchEta[charge]);
      
      //Denominator
      hDenMatchEta[charge] = new TH1F(Form("hDenMatch%sEta", pC[charge].Data()), Form("Matching Denominator in #eta %s;%s", pCharge[charge].Data(), etastring.Data()), kEtaBins, -fEtaRange, fEtaRange);
      hDenMatchEta[charge]->Sumw2();
      fListHist->AddLast(hDenMatchEta[charge]);
      
      //*****
      //PHI OUT
      //*****
      
      //Numerator
      hNumMatchphiOut[charge] = new TH1F(Form("hNumMatch%sphiOut", pC[charge].Data()), Form("Matching Numerator in #phi %s;%s", pCharge[charge].Data(), phistring.Data()), kPhiBins, 0, TMath::TwoPi());
      hNumMatchphiOut[charge]->Sumw2();
      fListHist->AddLast(hNumMatchphiOut[charge]);
      
      //Denominator
      hDenMatchphiOut[charge] = new TH1F(Form("hDenMatch%sphiOut", pC[charge].Data()), Form("Matching Denominator in #phi %s;%s", pCharge[charge].Data(), phistring.Data()), kPhiBins, 0, TMath::TwoPi());
      hDenMatchphiOut[charge]->Sumw2();
      fListHist->AddLast(hDenMatchphiOut[charge]);
      
      //*****
      //ETA (PT>0.5)
      //*****
      
      //Numerator
      hNumMatchEtaPtMa[charge] = new TH1F(Form("hNumMatch%sEtaPtMa", pC[charge].Data()), Form("Matching Numerator in #eta Max %s;%s", pCharge[charge].Data(), etastring.Data()), kEtaBins, -fEtaRange, fEtaRange);
      hNumMatchEtaPtMa[charge]->Sumw2();
      fListHist->AddLast(hNumMatchEtaPtMa[charge]);
      
      //Denominator
      hDenMatchEtaPtMa[charge] = new TH1F(Form("hDenMatch%sEtaPtMa", pC[charge].Data()), Form("Matching Denominator in #eta Max %s;%s", pCharge[charge].Data(), etastring.Data()), kEtaBins, -fEtaRange, fEtaRange);
      hDenMatchEtaPtMa[charge]->Sumw2();
      fListHist->AddLast(hDenMatchEtaPtMa[charge]);
      
      //*****
      //PHI OUT (PT>0.5)
      //*****
      
      //Numerator
      hNumMatchphiOutPtMa[charge] = new TH1F(Form("hNumMatch%sphiOutPtMa", pC[charge].Data()), Form("Matching Numerator in #phi Max %s;%s", pCharge[charge].Data(), phistring.Data()), kPhiBins, 0, TMath::TwoPi());
      hNumMatchphiOutPtMa[charge]->Sumw2();
      fListHist->AddLast(hNumMatchphiOutPtMa[charge]);
      
      //Denominator
      hDenMatchphiOutPtMa[charge] = new TH1F(Form("hDenMatch%sphiOutPtMa", pC[charge].Data()), Form("Matching Denominator in #phi Max %s;%s", pCharge[charge].Data(), phistring.Data()), kPhiBins, 0, TMath::TwoPi());
      hDenMatchphiOutPtMa[charge]->Sumw2();
      fListHist->AddLast(hDenMatchphiOutPtMa[charge]);
      
      //*****
      //TRD OUT
      //*****
      
      //Numerator
      hNumMatchTRDOut[charge] = new TH1F(Form("hNumMatch%sTRDOut", pC[charge].Data()), Form("Matching Numerator in TRD;%s", ptstring.Data()), kPtBins, fBinPt);
      hNumMatchTRDOut[charge]->Sumw2();
      fListHist->AddLast(hNumMatchTRDOut[charge]);
      
      //Denominator
      hDenMatchTRDOut[charge] = new TH1F(Form("hDenMatch%sTRDOut", pC[charge].Data()), Form("Matching Denominator in TRD;%s", ptstring.Data()), kPtBins, fBinPt);
      hDenMatchTRDOut[charge]->Sumw2();
      fListHist->AddLast(hDenMatchTRDOut[charge]);
      
      //*****
      //NO TRD OUT
      //*****
      
      //Numerator
      hNumMatchNoTRDOut[charge] = new TH1F(Form("hNumMatch%sNoTRDOut", pC[charge].Data()), Form("Matching Numerator with NO TRD;%s", ptstring.Data()), kPtBins, fBinPt);
      hNumMatchNoTRDOut[charge]->Sumw2();
      fListHist->AddLast(hNumMatchNoTRDOut[charge]);
      
      //Denominator
      hDenMatchNoTRDOut[charge] = new TH1F(Form("hDenMatch%sNoTRDOut", pC[charge].Data()), Form("Matching Denominator with NO TRD;%s", ptstring.Data()), kPtBins, fBinPt);
      hDenMatchNoTRDOut[charge]->Sumw2();
      fListHist->AddLast(hDenMatchNoTRDOut[charge]);
      
      //*****
      //PT, ETA and PHI
      //*****
      
      //Numerator
      hNumMatchPtEtaPhiout[charge] = new TH3I(Form("hNumMatchPtEtaPhiout%s", pC[charge].Data()), Form("Matching Numerator %s;%s;%s;%s", pCharge[charge].Data(), ptstring.Data(), etastring.Data(), phistring.Data()), kPtBins, 0, kPtBins, kEtaBins, -fEtaRange, fEtaRange, kPhiBins, 0, TMath::TwoPi());
      hNumMatchPtEtaPhiout[charge]->Sumw2();
      hNumMatchPtEtaPhiout[charge]->GetXaxis()->Set(kPtBins, fBinPt);
      fListHist->AddLast(hNumMatchPtEtaPhiout[charge]);
      
      //Denominator
      hDenMatchPtEtaPhiout[charge] = new TH3I(Form("hDenMatchPtEtaPhiout%s", pC[charge].Data()), Form("Matching Denominator %s;%s;%s;%s", pCharge[charge].Data(), ptstring.Data(), etastring.Data(), phistring.Data()), kPtBins, 0, kPtBins, kEtaBins, -fEtaRange, fEtaRange, kPhiBins, 0, TMath::TwoPi());
      hDenMatchPtEtaPhiout[charge]->Sumw2();
      hDenMatchPtEtaPhiout[charge]->GetXaxis()->Set(kPtBins, fBinPt);
      fListHist->AddLast(hDenMatchPtEtaPhiout[charge]);
      
      for(Int_t species = 0; charge < 2 && species < 3; species++){//Species loop Pi/Ka/Pr
        
        //*****
        //PID with TPC
        //*****
        
        //Numerator
        hNumMatchTPC[charge][species] = new TH1F(Form("hNumMatchTPC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Numerator with TPC PID for %s %s;%s", pCharge[charge].Data(), pSpecies[species].Data(), ptstring.Data()), kPtBins, fBinPt);
        hNumMatchTPC[charge][species]->Sumw2();
        fListHist->AddLast(hNumMatchTPC[charge][species]);
        
        //Denominator
        hDenMatchTPC[charge][species] = new TH1F(Form("hDenMatchTPC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Denominator with TPC PID for %s %s;%s", pCharge[charge].Data(), pSpecies[species].Data(), ptstring.Data()), kPtBins, fBinPt);
        hDenMatchTPC[charge][species]->Sumw2();
        fListHist->AddLast(hDenMatchTPC[charge][species]);
        
        //*****
        //Info from the DCAxy
        //*****
        for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
          for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
            hDCAxy[charge][species][ptbin][mult] = new TH1F(Form("hDCAxy%s%s_pt%i_mult%i", pC[charge].Data(), pS[species].Data(), ptbin, mult), Form("DCAxy Distribution of %s %s in pt %i [%.2f,%.2f] and mult %i;DCA_{xy} (cm);Counts", pCharge[charge].Data(), pSpecies[species].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin+1], mult), fDCAXYbins, -fDCAXYRange, fDCAXYRange);
            hDCAxy[charge][species][ptbin][mult]->Sumw2();
            fListHist->AddLast(hDCAxy[charge][species][ptbin][mult]);
            
            if(!fBuilDCAchi2) continue;
            hDCAxyGoldenChi2[charge][species][ptbin][mult] = new TH1F(Form("hDCAxyGoldenChi2%s%s_pt%i_mult%i", pC[charge].Data(), pS[species].Data(), ptbin, mult), Form("DCAxy w Golden Chi2 Distribution of %s %s in pt %i [%.2f,%.2f] and mult %i;DCA_{xy} (cm);Counts", pCharge[charge].Data(), pSpecies[species].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin+1], mult), fDCAXYbins, -fDCAXYRange, fDCAXYRange);
            hDCAxyGoldenChi2[charge][species][ptbin][mult]->Sumw2();
            fListHist->AddLast(hDCAxyGoldenChi2[charge][species][ptbin][mult]);
          }
        }
      }
    }
    
    
    if(fMCmode){
      for (Int_t charge = 0 ; charge < 2; charge++){//Charge loop
        for(Int_t species = 0; species < 3; species++){//Species loop
          hNumMatchMC[charge][species] = new TH1F(Form("hNumMatchMC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Numerator with MC PID %s %s;%s", pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
          hNumMatchMC[charge][species]->Sumw2();
          fListHist->AddLast(hNumMatchMC[charge][species]);
          
          hDenMatchMC[charge][species] = new TH1F(Form("hDenMatchMC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Denominator with MC PID %s %s;%s", pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
          hDenMatchMC[charge][species]->Sumw2();
          fListHist->AddLast(hDenMatchMC[charge][species]);
          
          hNumMatchPrimMC[charge][species] = new TH1F(Form("hNumMatchPrimMC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Numerator with MC PID and Primary %s %s;%s", pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
          hNumMatchPrimMC[charge][species]->Sumw2();
          fListHist->AddLast(hNumMatchPrimMC[charge][species]);
          
          hDenMatchPrimMC[charge][species] = new TH1F(Form("hDenMatchPrimMC%s%s", pC[charge].Data(), pS[species].Data()), Form("Matching Denominator with MC PID and Primary %s %s;%s", pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
          hDenMatchPrimMC[charge][species]->Sumw2();
          fListHist->AddLast(hDenMatchPrimMC[charge][species]);
          
          for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
            
            hNumMatchPrimMCYCut[charge][species][mult] = new TH1F(Form("hNumMatchPrimMCYCut%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), Form("Matching Numerator with MC PID and Primary |y| < %.1f %s %s;%s", fRapidityCut, pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
            hNumMatchPrimMCYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumMatchPrimMCYCut[charge][species][mult]);
            
            hDenMatchPrimMCYCut[charge][species][mult] = new TH1F(Form("hDenMatchPrimMCYCut%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), Form("Matching Denominator with MC PID and Primary |y| < %.1f %s %s;%s", fRapidityCut, pCharge[charge].Data(), pSpecies[charge].Data(), ptstring.Data()), kPtBins, fBinPt);
            hDenMatchPrimMCYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenMatchPrimMCYCut[charge][species][mult]);
            
          }
          
          //*****
          //Info from the DCAxy with MC information
          //*****
          
          for(Int_t ptbin = 0; ptbin < kPtBins; ptbin++){//Pt loop
            hDCAxyPrimMC[charge][species][ptbin] = new TH1F(Form("hDCAxyPrimMC%s%s_%i", pC[charge].Data(), pS[species].Data(), ptbin), Form("DCAxy Distribution of %s %s in pt %i [%.2f,%.2f];DCA_{xy} (cm);Counts", pCharge[charge].Data(), pSpecies[species].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin+1]), fDCAXYbins, -fDCAXYRange, fDCAXYRange);
            hDCAxyPrimMC[charge][species][ptbin]->Sumw2();
            fListHist->AddLast(hDCAxyPrimMC[charge][species][ptbin]);
            
            hDCAxySecStMC[charge][species][ptbin] = new TH1F(Form("hDCAxySecStMC%s%s_%i", pC[charge].Data(), pS[species].Data(), ptbin), Form("DCAxy Distribution of %s %s in pt %i [%.2f,%.2f];DCA_{xy} (cm);Counts", pCharge[charge].Data(), pSpecies[species].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin+1]), fDCAXYbins, -fDCAXYRange, fDCAXYRange);
            hDCAxySecStMC[charge][species][ptbin]->Sumw2();
            fListHist->AddLast(hDCAxySecStMC[charge][species][ptbin]);
            
            hDCAxySecMatMC[charge][species][ptbin] = new TH1F(Form("hDCAxySecMatMC%s%s_%i", pC[charge].Data(), pS[species].Data(), ptbin), Form("DCAxy Distribution of %s %s in pt %i [%.2f,%.2f];DCA_{xy} (cm);Counts", pCharge[charge].Data(), pSpecies[species].Data(), ptbin, fBinPt[ptbin], fBinPt[ptbin+1]), fDCAXYbins, -fDCAXYRange, fDCAXYRange);
            hDCAxySecMatMC[charge][species][ptbin]->Sumw2();
            fListHist->AddLast(hDCAxySecMatMC[charge][species][ptbin]);
          }
        }
      }
      
      //efficiency
      for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
        for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
          for(Int_t species = 0; species < 3; species++){//Species loop
            
            hNumMatchMultTrk[charge][species][mult] = new TH1F(Form("hNumMatch_MultTrk%i_%s%s", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hNumMatchMultTrk[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumMatchMultTrk[charge][species][mult]);
            
            hDenMatchMultTrk[charge][species][mult] = new TH1F(Form("hDenMatch_MultTrk%i_%s%s", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hDenMatchMultTrk[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenMatchMultTrk[charge][species][mult]);
            
            hNumMatchMultTrkTRDOut[charge][species][mult] = new TH1F(Form("hNumMatch_MultTrk%i_%s%s_TRDOut", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hNumMatchMultTrkTRDOut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumMatchMultTrkTRDOut[charge][species][mult]);
            
            hDenMatchMultTrkTRDOut[charge][species][mult] = new TH1F(Form("hDenMatch_MultTrk%i_%s%s_TRDOut", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hDenMatchMultTrkTRDOut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenMatchMultTrkTRDOut[charge][species][mult]);
            
            hNumMatchMultTrkNoTRDOut[charge][species][mult] = new TH1F(Form("hNumMatch_MultTrk%i_%s%s_NoTRDOut", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hNumMatchMultTrkNoTRDOut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumMatchMultTrkNoTRDOut[charge][species][mult]);
            
            hDenMatchMultTrkNoTRDOut[charge][species][mult] = new TH1F(Form("hDenMatch_MultTrk%i_%s%s_NoTRDOut", mult, pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
            hDenMatchMultTrkNoTRDOut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenMatchMultTrkNoTRDOut[charge][species][mult]);
            
          }
          
          hNumMatchMultTrkInc[charge][mult] = new TH1F(Form("hNumMatch_Inc_MultTrk%i_%s", mult, pC[charge].Data()), "", kPtBins, fBinPt);
          hNumMatchMultTrkInc[charge][mult]->Sumw2();
          fListHist->AddLast(hNumMatchMultTrkInc[charge][mult]);
          
          hDenMatchMultTrkInc[charge][mult] = new TH1F(Form("hDenMatch_Inc_MultTrk%i_%s", mult, pC[charge].Data()), "", kPtBins, fBinPt);
          hDenMatchMultTrkInc[charge][mult]->Sumw2();
          fListHist->AddLast(hDenMatchMultTrkInc[charge][mult]);
          
        }
        
      }
      
      for(Int_t charge = 0; charge < 2; charge++){//Charge loop Positive/Negative
        for(Int_t species = 0; species < 3; species++){//Species loop
          
          hDenTrkVertMultTrk[charge][species] = new TH1F(Form("hDenTrkVert_%s%s", pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
          hDenTrkVertMultTrk[charge][species]->Sumw2();
          fListHist->AddLast(hDenTrkVertMultTrk[charge][species]);
          
          hDenTrkTriggerMultTrk[charge][species] = new TH1F(Form("hDenTrkTrigger_%s%s", pC[charge].Data(), pS[species].Data()), "", kPtBins, fBinPt);
          hDenTrkTriggerMultTrk[charge][species]->Sumw2();
          fListHist->AddLast(hDenTrkTriggerMultTrk[charge][species]);
          
          for(Int_t mult = 0; mult < kEvtMultBins; mult++){//Multiplicity loop
            hDenPrimMCYCut[charge][species][mult] = new TH1F(Form("hDenPrimMCYCut_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "Primary particles", kPtBins, fBinPt);
            hDenPrimMCYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenPrimMCYCut[charge][species][mult]);
            
            hDenPrimMCEtaCut[charge][species][mult] = new TH1F(Form("hDenPrimMCEtaCut_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hDenPrimMCEtaCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenPrimMCEtaCut[charge][species][mult]);
            
            hDenPrimMCEtaYCut[charge][species][mult] = new TH1F(Form("hDenPrimMCEtaYCut_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hDenPrimMCEtaYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hDenPrimMCEtaYCut[charge][species][mult]);
            
            hNumPrimMCTrueMatch[charge][species][mult] = new TH1F(Form("hNumPrimMCTrueMatch_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hNumPrimMCTrueMatch[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumPrimMCTrueMatch[charge][species][mult]);
            
            hNumPrimMCTrueMatchYCut[charge][species][mult] = new TH1F(Form("hNumPrimMCTrueMatchYCut_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hNumPrimMCTrueMatchYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumPrimMCTrueMatchYCut[charge][species][mult]);
            
            hNumPrimMCTrueMatchYCutTPC[charge][species][mult] = new TH1F(Form("hNumPrimMCTrueMatchYCutTPC_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hNumPrimMCTrueMatchYCutTPC[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumPrimMCTrueMatchYCutTPC[charge][species][mult]);
            
            hNumPrimMCConsistentMatchYCut[charge][species][mult] = new TH1F(Form("hNumPrimMCConsistentMatchYCut_%s%s_%i", pC[charge].Data(), pS[species].Data(), mult), "", kPtBins, fBinPt);
            hNumPrimMCConsistentMatchYCut[charge][species][mult]->Sumw2();
            fListHist->AddLast(hNumPrimMCConsistentMatchYCut[charge][species][mult]);
          }
        }
      }
      
    }
    
    
    //TTree with the data information
    if(fTreemode){
      OpenFile(2);
      
      
      fTreeTrack = new TTree("fTreeTrack", "Track Properties");
      //
      //Event wide variables
      //fTreeTrack->Branch("fPrimVertex", &fPrimVertex, "fPrimVertex[3]/D");
      //fTreeTrack->Branch("fNContrPrimVertex", &fNContrPrimVertex, "fNContrPrimVertex/I");
      fTreeTrack->Branch("fEvtMultBin", &fEvtMultBin, "fEvtMultBin/S");
      //
      //Track variables
      #ifdef USETREECLASS
      ArrayAnTrk = new TClonesArray("AliAnTOFtrack");
      ArrayAnTrk->SetOwner();
      AliAnTOFtrack::Class()->IgnoreTObjectStreamer();
      fTreeTrack->Branch("AliAnTOFtrack", "TClonesArray", &ArrayAnTrk, 8000, 0);
      #else
      //       fTreeTrack->Branch("fDCAXY", &fDCAXY, "fDCAXY/F");
      //       fTreeTrack->Branch("fDCAZ", &fDCAZ, "fDCAZ/F");
      fTreeTrack->Branch("fTrkMask", &fTrkMask, "fTrkMask/b");//Track information in Mask
      fTreeTrack->Branch("fTPCPIDMask", &fTPCPIDMask, "fTPCPIDMask/b");//TPC PID information in Mask
      if(fCutmode) fTreeTrack->Branch("fTrkCutMask", &fTrkCutMask, "fTrkCutMask/s");//information for the cut variation
      //     fTreeTrack->Branch("fMismatch", &fMismatch, "fMismatch/O");
      //     fTreeTrack->Branch("fT0UsedMask", &fT0UsedMask, "fT0UsedMask/I");
      //fTreeTrack->Branch("TOFout", &fTOFout, "TOFout/O");
      //     fTreeTrack->Branch("fTRDout", &fTRDout, "fTRDout/O");
      //fTreeTrack->Branch("fTime", &fTime, "fTime/O");
      //     fTreeTrack->Branch("fP", &fP, "fP/D");
      //fTreeTrack->Branch("fTPCClusters", &fTPCClusters, "fTPCClusters/b");
      fTreeTrack->Branch("fLength", &fLength, "fLength/D");
      fTreeTrack->Branch("fLengthRatio", &fLengthRatio, "fLengthRatio/F");
      fTreeTrack->Branch("fSign", &fSign, "fSign/O");
      fTreeTrack->Branch("fTOFTime", &fTOFTime, "fTOFTime/D");
      fTreeTrack->Branch("fTOFMismatchTime", &fTOFMismatchTime, "fTOFMismatchTime/F");
      for(Int_t ipart = 0; ipart < kExpSpecies; ipart++) fTreeTrack->Branch(Form("fTOFExpTime%i",ipart), &fTOFExpTime[ipart], Form("fTOFExpTime%i/F",ipart));
      //fTreeTrack->Branch("fTOFImpactDZ", &fTOFImpactDZ, "fTOFImpactDZ/D");
      //fTreeTrack->Branch("fTOFImpactDX", &fTOFImpactDX, "fTOFImpactDX/D");
      fTreeTrack->Branch("fT0TrkTime", &fT0TrkTime, "fT0TrkTime/F");
      fTreeTrack->Branch("fTOFchan", &fTOFchan, "fTOFchan/I");
      fTreeTrack->Branch("fEta", &fEta, "fEta/D");
      fTreeTrack->Branch("fPhi", &fPhi, "fPhi/D");
      //fTreeTrack->Branch("fTOFSignalTot", &fTOFSignalTot, "fTOFSignalTot/F");
      fTreeTrack->Branch("fPt", &fPt, "fPt/D");
      for(Int_t ipart = 0; ipart < kExpSpecies; ipart++) fTreeTrack->Branch(Form("fTOFExpSigma%i",ipart), &fTOFExpSigma[ipart], Form("fTOFExpSigma%i/F",ipart));
      fTreeTrack->Branch("fNTOFClusters", &fNTOFClusters, "fNTOFClusters/I");
      
      fTreeTrack->Branch("fTPCSignal", &fTPCSignal, "fTPCSignal/F");
      //     fTreeTrack->Branch("fTPCSigma0", &fTPCSigma[ke], "fTPCSigma0/F");
      //     fTreeTrack->Branch("fTPCSigma1", &fTPCSigma[kmu], "fTPCSigma1/F");
      //     fTreeTrack->Branch("fTPCSigma2", &fTPCSigma[kpi], "fTPCSigma2/F");
      //     fTreeTrack->Branch("fTOFSigma0", &fTOFSigma[ke], "fTOFSigma0/F");
      //     fTreeTrack->Branch("fTOFSigma1", &fTOFSigma[kmu], "fTOFSigma1/F");
      //     fTreeTrack->Branch("fTOFSigma2", &fTOFSigma[kpi], "fTOFSigma2/F");
      //fTreeTrack->Branch("fTOFPIDProbability0", &fTOFPIDProbability[0], "fTOFPIDProbability0/D");
      //fTreeTrack->Branch("fTOFPIDProbability1", &fTOFPIDProbability[1], "fTOFPIDProbability1/D");
      //fTreeTrack->Branch("fTOFPIDProbability2", &fTOFPIDProbability[2], "fTOFPIDProbability2/D");
      //fTreeTrack->Branch("fTOFPIDProbability3", &fTOFPIDProbability[3], "fTOFPIDProbability3/D");
      //fTreeTrack->Branch("fTOFPIDProbability4", &fTOFPIDProbability[4], "fTOFPIDProbability4/D");
      
      //fTreeTrack->Branch("fPhiout", &fPhiout, "fPhiout/D");
      //fTreeTrack->Branch("fXout", &fXout, "fXout/D");
      //fTreeTrack->Branch("fYout", &fYout, "fYout/D");
      //fTreeTrack->Branch("fZout", &fZout, "fZout/D");
      
      if(fMCmode){//Ulterior branches for the MC
        fTreeTrack->Branch("fMCTrkMask", &fMCTrkMask, "fMCTrkMask/b");
      }
      #endif
      
      
      if(0 && fMCmode){
        fTreeTrackMC = new TTree("fTreeTrackMC", "Track MC Properties");
        fTreeTrackMC->Branch("fEvtMultBin", &fEvtMultBin, "fEvtMultBin/S");
        fTreeTrackMC->Branch("fDCAXY", &fDCAXY, "fDCAXY/F");
        fTreeTrackMC->Branch("fDCAZ", &fDCAZ, "fDCAZ/F");
        fTreeTrackMC->Branch("TOFout", &fTOFout, "TOFout/O");
        fTreeTrackMC->Branch("fTRDout", &fTRDout, "fTRDout/O");
        fTreeTrackMC->Branch("fTime", &fTime, "fTime/O");
        fTreeTrackMC->Branch("fP", &fP, "fP/D");
        fTreeTrackMC->Branch("fTPCClusters", &fTPCClusters, "fTPCClusters/b");
        fTreeTrackMC->Branch("fLength", &fLength, "fLength/D");
        fTreeTrackMC->Branch("fSign", &fSign, "fSign/O");
        fTreeTrackMC->Branch("fTOFTime", &fTOFTime, "fTOFTime/D");
        fTreeTrackMC->Branch("fTOFExpTime0", &fTOFExpTime[ke], "fTOFExpTime0/F");
        fTreeTrackMC->Branch("fTOFExpTime1", &fTOFExpTime[kmu], "fTOFExpTime1/F");
        fTreeTrackMC->Branch("fTOFExpTime2", &fTOFExpTime[kpi], "fTOFExpTime2/F");
        fTreeTrackMC->Branch("fTOFExpTime3", &fTOFExpTime[kK], "fTOFExpTime3/F");
        fTreeTrackMC->Branch("fTOFExpTime4", &fTOFExpTime[kp], "fTOFExpTime4/F");
        fTreeTrackMC->Branch("fTOFImpactDZ", &fTOFImpactDZ, "fTOFImpactDZ/D");
        fTreeTrackMC->Branch("fTOFImpactDX", &fTOFImpactDX, "fTOFImpactDX/D");
        fTreeTrackMC->Branch("fT0TrkTime", &fT0TrkTime, "fT0TrkTime/F");
        fTreeTrackMC->Branch("fTOFchan", &fTOFchan, "fTOFchan/I");
        fTreeTrackMC->Branch("fEta", &fEta, "fEta/D");
        fTreeTrackMC->Branch("fPhi", &fPhi, "fPhi/D");
        fTreeTrackMC->Branch("fTOFSignalTot", &fTOFSignalTot, "fTOFSignalTot/F");
        fTreeTrackMC->Branch("fPt", &fPt, "fPt/D");
        fTreeTrackMC->Branch("fTOFExpSigma0", &fTOFExpSigma[0], "fTOFExpSigma0/F");
        fTreeTrackMC->Branch("fTOFExpSigma1", &fTOFExpSigma[1], "fTOFExpSigma1/F");
        fTreeTrackMC->Branch("fTOFExpSigma2", &fTOFExpSigma[2], "fTOFExpSigma2/F");
        fTreeTrackMC->Branch("fTOFExpSigma3", &fTOFExpSigma[3], "fTOFExpSigma3/F");
        fTreeTrackMC->Branch("fTOFExpSigma4", &fTOFExpSigma[4], "fTOFExpSigma4/F");
        fTreeTrackMC->Branch("fTPCSignal", &fTPCSignal, "fTPCSignal/F");
        fTreeTrackMC->Branch("fTPCSigma0", &fTPCSigma[ke], "fTPCSigma0/F");
        fTreeTrackMC->Branch("fTPCSigma1", &fTPCSigma[kmu], "fTPCSigma1/F");
        fTreeTrackMC->Branch("fTPCSigma2", &fTPCSigma[kpi], "fTPCSigma2/F");
        //       fTreeTrackMC->Branch("fTOFSigma0", &fTOFSigma[ke], "fTOFSigma0/F");
        //       fTreeTrackMC->Branch("fTOFSigma1", &fTOFSigma[kmu], "fTOFSigma1/F");
        //       fTreeTrackMC->Branch("fTOFSigma2", &fTOFSigma[kpi], "fTOFSigma2/F");
        
        fTreeTrackMC->Branch("fTOFPIDProbability0", &fTOFPIDProbability[0], "fTOFPIDProbability0/D");
        fTreeTrackMC->Branch("fTOFPIDProbability1", &fTOFPIDProbability[1], "fTOFPIDProbability1/D");
        fTreeTrackMC->Branch("fTOFPIDProbability2", &fTOFPIDProbability[2], "fTOFPIDProbability2/D");
        fTreeTrackMC->Branch("fTOFPIDProbability3", &fTOFPIDProbability[3], "fTOFPIDProbability3/D");
        fTreeTrackMC->Branch("fTOFPIDProbability4", &fTOFPIDProbability[4], "fTOFPIDProbability4/D");
        
        fTreeTrackMC->Branch("fPhiout", &fPhiout, "fPhiout/D");
        fTreeTrackMC->Branch("fXout", &fXout, "fXout/D");
        fTreeTrackMC->Branch("fYout", &fYout, "fYout/D");
        fTreeTrackMC->Branch("fZout", &fZout, "fZout/D");
        
      }
      
      
    }
    
  }
  
  // Post output data.
  PostAllTheData();
  AliInfo("UserCreateOutputObjects\t END");
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::UserExec(Option_t *){
  AliInfo("UserExec");
  AliDebug(2, Form("Event: %.0f", hNEvt->GetBinContent(1)+1));
  
  
  StartTimePerformance(0);
  //Initialize all the variables to their default values
  InitializeEventVar();
  StopTimePerformance(0);
  
  //
  // main event loop
  //
  
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    AliError("fESD not available");
    // Post output data.
    PostAllTheData();
    return;
  }
  
  if(fChannelmode){
    RunTOFChannel();
    AliInfo("TOF Channel run terminated");
    // Post output data.
    PostAllTheData();
    return;
  }
  
  if(fRecalibrateTOF){
    if(!TOFCalibInitRun()){//Done once for all the run
      AliError("Required Run TOF re-calibration was not successful");
      // Post output data.
      PostAllTheData();
      return;
    }
    if(!TOFCalibInitEvent()){//Done for each event
      AliError("Required Event TOF re-calibration was not successful");
      // Post output data.
      PostAllTheData();
      return;
    }
  }
  
  Int_t EvtStart = 0;
  hNEvt->Fill(EvtStart++);//Number of events opened -->Read from ESD
  
  if (!fESDtrackCuts) {
    AliError("fESDtrackCuts not available");
    // Post output data.
    PostAllTheData();
    return;
  }
  
  hNEvt->Fill(EvtStart++);//Has AliESDtrackCuts
  
  // AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  // if (esdH) fESDpid = esdH->GetESDpid();
  
  //
  // create PID response
  //
  AliAnalysisManager *AnManager = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (AnManager->GetInputEventHandler());
  fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
  if(fRecalibrateTOF) fPIDResponse->SetTOFResponse(fESD, AliPIDResponse::kBest_T0);
  fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  if(fRecalibrateTOF) fTOFPIDResponse.SetTimeResolution(fTimeResolution);
  
  //
  //Physics Selection
  //
  //
  // check if event is selected by physics selection class
  //
  //   TString firedTriggerClass = fInputEvent->GetFiredTriggerClasses();
  //   AliInfo(Form("Fired Trigger Classes for Event %.0f\n%s", hNEvt->GetBinContent(1), firedTriggerClass.Data()));
  
  UInt_t PhysSelmask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  fEvtPhysSelected = (PhysSelmask & fSelectBit);//--> Physics Selection
  
  //
  //Physics Selection Cut
  //
  if (!fEvtPhysSelected) {
    AliDebug(2, Form("Event %.0f did not pass the physics selection", hNEvt->GetBinContent(1)));
    // Post output data.
    PostAllTheData();
    return;
  }
  
  //   AliInfo(Form("Event %.0f Passed the Phys. Sel. + Trig", hNEvt->GetBinContent(1)));
  hNEvt->Fill(EvtStart++); //-->Pass Phys. Sel. + Trig
  
  //
  //Event Selection
  //
  //
  //Check the Event Multiplicity, the event selection is embedded in the Multiplicity selection with the codes reported in the AliMultSelectionCuts
  //
  StartTimePerformance(1);
  if(fHImode){
    fMultSel = (AliMultSelection * ) fESD->FindListObject("MultSelection");
    if(!fMultSel) {
      //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
    }
    else{
      AliDebug(2, "Estimating centrality");
      fEvtMult = fMultSel->GetMultiplicityPercentile("V0M", kTRUE);//Event selection is embedded in the Multiplicity estimator so that the Multiplicity percentiles are well defined and refer to the same sample
    }
    
  }
  else{
    AliDebug(2, "Estimating Multiplicity at midrapidity");
    fEvtMult = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD, kTRUE);//fESDtrackCuts->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
    
  }
  StopTimePerformance(1);
  
  StartTimePerformance(2);
  fEvtSelected = SelectEvents(EvtStart);
  StopTimePerformance(2);
  
  StartTimePerformance(3);
  ComputeEvtMultiplicityBin();//Calculate in the handy binning the Multiplicity bin of the event
  StopTimePerformance(3);
  
  if(fMCmode) AnalyseMCParticles(); //First loop on stack Before the Physics Selection (and also after) Before the Event Selection (and also after)
  
  //
  //Filling the Multiplicity histogram before the event selection
  //
  hEvtMult->Fill(fEvtMult);
  
  //
  //Event Selection Cut
  //
  if(!fEvtSelected){//Event selection based on the Multiplicity selector
    AliDebug(2, "Did not pass the event selection");
    // Post the output data.
    PostAllTheData();
    return;
  }
  
  
  //
  //Now events are selected, there cannot be any inconsistency!!!
  //
  // monitor vertex position
  //
  const AliESDVertex* vertex = ObtainVertex();
  if (!vertex) {
    // Post output data.
    AliFatal(Form("Event selected for the analysis has vertex status %i (should be 3) and will be rejected!!!", fVertStatus));
    PostAllTheData();
    return;
  }
  
  AliDebug(2, Form("Event %.0f Has a good vertex", hNEvt->GetBinContent(1)));
  
  //Filling XY and Z distribution for vertex
  hEvtVtxXY->Fill(TMath::Sqrt(fPrimVertex[0]*fPrimVertex[0] + fPrimVertex[1]*fPrimVertex[1]));
  hEvtVtxZ->Fill(fPrimVertex[2]);
  
  //
  //Filling the Multiplicity histogram after the event selection
  //
  hEvtMultAftEvSel->Fill(fEvtMult);
  
  //
  // TOF geometrical variables
  //
  Int_t det[5];
  Float_t coord[3];
  
  //
  // track loop
  //
  StartTimePerformance(4);
  const Int_t nTrack = fESD->GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {
    if(iTrack == 0) AliDebug(2, "Starting loop on the tracks");
    if(iTrack == nTrack -1) AliDebug(2, "Ending loop on the tracks");
    
    if(iTrack == 0) SetEvtMaskBit(kIsNewEvent, 1);
    else if(iTrack == 1) SetEvtMaskBit(kIsNewEvent, 0);
    
    
    const AliESDtrack *track = fESD->GetTrack(iTrack);
    if (!track){continue;}
    hNTrk->Fill(0);// --> Tracks in Input
    
    //
    // start TOF analysis
    InitializeTrackVar();//Initialize the values needed at each track iteration
    
    //
    // Get the information on the track variables
    SetTrackValues(track, vertex);//Get the information on the track values
    
    //
    // Get the information on the track flags relative to the TOF but not only
    SetTrackFlags(track);
    
    //
    // Put the information into the histograms
    FillCutVariable(kFALSE);
    
    //
    //Track Cuts
    if (!fESDtrackCuts->AcceptTrack(track)) {//WARNING It is important that here the accepted tracks are not selected with the DCA cut, this will be done after filling the DCA histograms
      AliDebug(2, Form("Event %.0f Has track %i/%i  -> Not Accepted", hNEvt->GetBinContent(1), iTrack, nTrack));
      continue;
    }
    else AliDebug(2, Form("Event %.0f Has track %i/%i  -> Accepted", hNEvt->GetBinContent(1), iTrack, nTrack));
    hNTrk->Fill(1);// --> Tracks passes Track Selection
    
    if(TMath::Abs(fEta) > fEtaRange) continue;//Pseudorapidity cut
    hNTrk->Fill(2);// --> Eta cut
    
    if(track->GetSign() == 0) continue; //No neutral particles
    hNTrk->Fill(3);// --> Charged track
    
    FindPtBin();                                            // Find the corresponding Pt bin
    if(fBinPtIndex < 0 || fBinPtIndex >= kPtBins) continue; //Track has higher momentum than actual binning this cut makes sense because if you don't have the DCA templates you cannot correct for primaries!!
    hNTrk->Fill(4);// --> Track in pt range
    
    //
    //Compute the track Rapidity for the three hypothesis (pi/k/p)
    ComputeRapidity();
    
    //
    //Monte Carlo information of the track
    if(fMCmode && !GatherTrackMCInfo(track)){
      // Post output data.
      AliError("Track had a problem accessing the MC info");
      PostAllTheData();
      return;
    }
    
    //
    //TOF clusters
    //if(track->GetTOFHeader() && track->GetTOFHeader()->GetTriggerMask() && track->GetTOFHeader()->GetNumberOfTOFclusters() > -1) fNTOFClusters = track->GetTOFHeader()->GetNumberOfTOFclusters(); //all fired readout pads
    fNTOFClusters = track->GetNTOFclusters(); //All matchable clusters
    
    //
    //Fill the histogram with the distribution of TOF clusters
    //
    hTOFClusters->Fill(fNTOFClusters);
    
    //
    //Get the array of matchable clusters in TOF
    //
    fTOFClusters = track->GetTOFclusterArray();  //Index of the matchable cluster (there are fNTOFClusters of them)
    
    if(fNTOFClusters > 1){
      //Get the mismatch time from cluster which have poor quality, always taking the worst one! Doing so, it is possible to select when running localy the number of clusters to use at least! (e.g. 6/7/8 etc..)
      Int_t clsindex = -1;
      fTOFcls = (AliESDTOFCluster *) fESD->GetESDTOFClusters()->At(fTOFClusters[fNTOFClusters-1]);
      for(Int_t i = 0; i < fTOFcls->GetNMatchableTracks(); i++){
        if(fTOFcls->GetTrackIndex(i) == track->GetID()){
          clsindex = i;
          break;
        }
      }
      
      if(clsindex > -1){
        fTOFMismatchTime = fTOFcls->GetTime(clsindex);
        fLengthRatio = fLength > 0 ? fTOFcls->GetLength(clsindex)/fLength : -999;
      }
    }
    
    ////////////////
    //Track PID Info
    ////////////////
    
    //
    //PID with the TOF
    //
    fTOFTime = track->GetTOFsignal(); //Gets the TOF signal
    
    Double_t inttime[kExpSpecies];
    track->GetIntegratedTimes(inttime, kExpSpecies);// Returns the array with integrated times for each particle hypothesis
    fTOFExpTime[ke] = inttime[AliPID::kElectron];
    fTOFExpTime[kmu] = inttime[AliPID::kMuon];
    fTOFExpTime[kpi] = inttime[AliPID::kPion];
    fTOFExpTime[kK] = inttime[AliPID::kKaon];
    fTOFExpTime[kp] = inttime[AliPID::kProton];
    fTOFExpTime[kd] = inttime[AliPID::kDeuteron];
    
    #ifdef CHECKCOMPUTEDVALUES
    for(Int_t species = 0; species < kExpSpecies && fTOFout && fTime && (fP < 10.); species++){//Check on the expected time
      hTOFMomComputed[species]->Fill(ComputeExpectedMomentum(AliPID::ParticleMass(AliPID::kElectron+species), fTOFExpTime[species])/fP);
      hTOFMomComputedTPC[species]->Fill(ComputeExpectedMomentum(AliPID::ParticleMass(AliPID::kElectron+species), fTOFExpTime[species])/fPTPC);
      
      hTOFExpectedComputed[species]->Fill(ComputeExpectedTime(AliPID::ParticleMass(AliPID::kElectron+species), fP)/fTOFExpTime[species] );
      hTOFExpectedComputedTPC[species]->Fill(ComputeExpectedTime(AliPID::ParticleMass(AliPID::kElectron+species), fPTPC)/fTOFExpTime[species]);
    }
    #endif
    
    fTOFSignalTot = track->GetTOFsignalToT();// Time of threshold signal of the TOF
    fTOFImpactDZ = track->GetTOFsignalDz();  // local z  of track's impact on the TOF pad
    fTOFImpactDX = track->GetTOFsignalDx();  // local x  of track's impact on the TOF pad
    fTOFchan = track->GetTOFCalChannel();    // Channel Index of the TOF Signal
    
    fT0TrkTime = fTOFPIDResponse.GetStartTime(fP);      // T0best time
    fT0UsedMask = fTOFPIDResponse.GetStartTimeMask(fP); // T0best used time  ->  mask with the T0 used (0x1=T0-TOF,0x2=T0A,0x3=TOC) for p bins
    fT0TrkSigma = fTOFPIDResponse.GetStartTimeRes(fP);  // T0best resolution time
    
    //Fill histograms with start-time information
    hT0->Fill(fT0TrkTime);
    hT0Resolution->Fill(fT0TrkSigma);
    
    //Get TOF Expected Sigma
    for (Int_t ipart = 0; ipart < kExpSpecies; ipart++) fTOFExpSigma[ipart] = fTOFPIDResponse.GetExpectedSigma(fP, inttime[ipart], AliPID::ParticleMass(ipart));
    
    //Get TOF Separation in number of sigmas
    fTOFSigma[ke] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    fTOFSigma[kmu] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
    fTOFSigma[kpi] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    fTOFSigma[kK] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    fTOFSigma[kp] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    fTOFSigma[kd] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
    
    //Compute mismatch Probability
    if(fPIDResponse->ComputeTOFProbability(track, kExpSpecies, fTOFPIDProbability) ==  AliPIDResponse::kDetMismatch) fMismatch = kTRUE;
    else fMismatch = kFALSE;
    
    //
    //PID with the TPC
    //
    
    fTPCSignal = track->GetTPCsignal();
    
    //Get TPC Separation in number of sigmas
    fTPCSigma[ke] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fTPCSigma[kmu] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
    fTPCSigma[kpi] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) - fTPCShift[0][fBinPtIndex];
    fTPCSigma[kK] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon) - fTPCShift[1][fBinPtIndex];
    fTPCSigma[kp] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton) - fTPCShift[2][fBinPtIndex];
    fTPCSigma[kd] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    
    //Set the mask flags for TPC PID
    if(TMath::Abs(fTPCSigma[ke]) < 1.5) SetTPCPIDMaskBit(kIsTPCElectron, 1);//1.5 sigma cut for Electrons
    if(TMath::Abs(fTPCSigma[kmu]) < 1.5) SetTPCPIDMaskBit(kIsTPCMuon, 1);   //1.5 sigma cut for Muons
    if(TMath::Abs(fTPCSigma[kpi]) < 5) SetTPCPIDMaskBit(kIsTPCPion, 1);     //5.0 sigma cut for Pions
    if(TMath::Abs(fTPCSigma[kK]) < 5) SetTPCPIDMaskBit(kIsTPCKaon, 1);      //5.0 sigma cut for Kaons
    if(TMath::Abs(fTPCSigma[kp]) < 5) SetTPCPIDMaskBit(kIsTPCProton, 1);    //5.0 sigma cut for Protons
    if(TMath::Abs(fTPCSigma[kd]) < 1.5) SetTPCPIDMaskBit(kIsTPCDeuteron, 1);//1.5 sigma cut for Deuterons
    
    
    //Get Geometrical Parameters of the track
    AliExternalTrackParam* exttrack = (AliExternalTrackParam *)track->GetOuterParam();
    if(exttrack){
      fPhiout = exttrack->Phi();
      fXout = exttrack->GetX();
      fYout = exttrack->GetY();
      fZout = exttrack->GetZ();
    }
    
    //
    //PID with the Combined TPC/TOF, of course this part needs both good TPC and TOF signals
    //
    if(fTOFout && fTime){//Only if tracks are matched to a TOF signal, otherwise it shall not be in the histogram!
      //Test the TPC/TOF separation
      if(fBuilTPCTOF){
        for(Int_t species = 0; species < kExpSpecies; species++) hTPCTOFSeparation[species][fBinPtIndex]->Fill(fTOFSigma[species], fTPCSigma[species]);
      }
      
      for(Int_t species = 0; species < kSpecies; species++){//Species Loop
        //
        //Compute the combined separation with TPC and TOF
        //
        fCombinedSigma[species] = TMath::Sqrt(fTPCSigma[species + kpi]*fTPCSigma[species + kpi] + fTOFSigma[species + kpi]*fTOFSigma[species + kpi]);
        
      }
      
      //
      //DCA information
      //
      for(Int_t species = 0; species < 3; species++){//Combined information of the TOF and TPC signal to identify a pure sample
        if(fCombinedSigma[species] > 2.) continue;
        if(TMath::Abs(fRapidity[species]) >= fRapidityCut) continue;//Rapidity cut
        
        //Pure sample selected via TOF and TPC 2 sigma cut
        hDCAxy[fSign][species][fBinPtIndex][fEvtMultBin]->Fill(fDCAXY);
        
        //
        //DCA - with golden chi2 cut
        //
        if(!fBuilDCAchi2 || !fPassGoldenChi2) continue;
        hDCAxyGoldenChi2[fSign][species][fBinPtIndex][fEvtMultBin]->Fill(fDCAXY);
      }
      
      if(fPdgIndex != -999 && TMath::Abs(fRapidity[fPdgIndex]) < fRapidityCut){//MC info If the track is a Pion or a Kaon or a Proton
        if(fProdInfo == 0) hDCAxyPrimMC[fSign][fPdgIndex][fBinPtIndex]->Fill(fDCAXY);
        else if(fProdInfo == 1) hDCAxySecStMC[fSign][fPdgIndex][fBinPtIndex]->Fill(fDCAXY);
        else if(fProdInfo == 2) hDCAxySecMatMC[fSign][fPdgIndex][fBinPtIndex]->Fill(fDCAXY);
      }
      
    }
    //
    //DCAxy cut AND Ulterior cuts to match the requirements of the TOF
    if((fPassDCAxy = fESDtrackCutsPrm->AcceptTrack(track))){
      hNTrk->Fill(5);
      if(fTOFout == kTRUE){
        hNTrk->Fill(6);
        if(fTime == kTRUE){
          hNTrk->Fill(7);
          if(fTRDout == kTRUE){
            hNTrk->Fill(8);
            if(fLength > fLengthmin){
              hNTrk->Fill(9);
              if(fTOFTime > fTOFmin){
                hNTrk->Fill(10);
                if(fTOFExpTime[AliPID::kPion] > fTOFmin){
                  hNTrk->Fill(11);
                  if(fTOFExpTime[AliPID::kKaon] > fTOFmin){
                    hNTrk->Fill(12);
                    if(fTOFExpTime[AliPID::kProton] > fTOFmin){
                      hNTrk->Fill(13);
                      if(fTOFTime < fTOFmax) hNTrk->Fill(14);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    //
    //DCAxy cut
    //
    if(!fPassDCAxy){
      AliDebug(2, Form("Event %.0f Has track %i/%i  -> Not Accepted as primary", hNEvt->GetBinContent(1), iTrack, nTrack));
      continue;
    }
    AliDebug(2, Form("Event %.0f Has track %i/%i  -> Accepted as primary", hNEvt->GetBinContent(1), iTrack, nTrack));
    
    //
    //TOF clusters for particles which passed the DCA cut
    //
    hTOFClustersDCApass->Fill(fNTOFClusters);
    
    
    //
    //Check on the cut variables
    //
    FillCutVariable(kTRUE);
    
    //
    //Matching efficiency plots
    //Denominator
    for(Int_t i = 0; i < 3; i++){//Loop on pi/k/p
      Bool_t pass = kTRUE;
      for(Int_t j = 0; j < 3; j++){
        if(i == j){
          if(TMath::Abs(fTPCSigma[kpi + j]) > 3) pass = kFALSE;
        }
        else if(TMath::Abs(fTPCSigma[kpi + j]) < 3) pass = kFALSE;
      }
      if(pass) hDenMatchTPC[fSign][i]->Fill(fPt);
    }
    
    //
    //All Reconstructed tracks
    //
    hDenMatch[fSign]->Fill(fPt);                           //-> transverse momentum
    hDenMatchEta[fSign]->Fill(fEta);                       //-> eta
    hDenMatchphiOut[fSign]->Fill(fPhiout);                 //-> Phi
    hDenMatchPtEtaPhiout[fSign]->Fill(fPt, fEta, fPhiout); //-> transverse momentum, eta and phi
    
    if(fPt > 0.5){                                         //Lower cut on the transver momentum
      hDenMatchEtaPtMa[fSign]->Fill(fEta);                 //-> eta
      hDenMatchphiOutPtMa[fSign]->Fill(fPhiout);           //-> Phi
    }
    if(fTRDout) hDenMatchTRDOut[fSign]->Fill(fPt);         //Tracks with TRDout
    else hDenMatchNoTRDOut[fSign]->Fill(fPt);              //Tracks without TRDout
    
    
    if(fMCmode){//start to take MC data
      if(fPdgIndex != -999){//If the track is a Pion or a Kaon or a Proton
        hDenMatchMC[fSign][fPdgIndex]->Fill(fPt);
        if(fProdInfo == 0){
          hDenMatchPrimMC[fSign][fPdgIndex]->Fill(fPt);
          if(TMath::Abs(fRapidityMC) < fRapidityCut) hDenMatchPrimMCYCut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
        }
        
        //Setting the MC track information in the mask
        SetMCTrkMaskBit(kIsNegative, fSign);//Set the sign
        
        switch (fPdgIndex) {//Set the particle species
          case 0:
          SetMCTrkMaskBit(kIsPion, 1);
          break;
          case 1:
          SetMCTrkMaskBit(kIsKaon, 1);
          break;
          case 2:
          SetMCTrkMaskBit(kIsProton, 1);
          break;
          default:
          break;
        }
        
        switch (fProdInfo) {//Set the particle type
          case 0:
          SetMCTrkMaskBit(kIsPhysicalPrimary, 1);
          break;
          case 1:
          SetMCTrkMaskBit(kIsFromStrangeness, 1);
          break;
          case 2:
          SetMCTrkMaskBit(kIsFromMaterial, 1);
          break;
          default:
          break;
        }
        
      }
      //       AliDebug(2, "Gathered the track MC Info");
    }
    
    //Setting the track information in the Track mask
    SetTrkMaskBit(kNegTrk, fSign);//Set the sign of the track
    SetTrkMaskBit(kIsMismatch, fMismatch);//Set the mismatch
    
    // "FILL"        -> 0
    // "TOF"         -> 1
    // "T0A"         -> 2
    // "TOF.and.T0A" -> 3
    // "T0C"         -> 4
    // "TOF.and.T0C" -> 5
    // "T0AC;        -> 6
    // "TOF.and.T0AC"-> 7
    
    // cout<< "T0 used: ["<<fT0UsedMask<<"]"<<endl;
    if(fT0UsedMask == 1 || fT0UsedMask == 3 || fT0UsedMask == 5 || fT0UsedMask == 7){
      SetTrkMaskBit(kT0_0, 1);//Set the used T0
      // cout<<"Is one"<<endl;
    }
    if(fT0UsedMask == 2 || fT0UsedMask == 3 || fT0UsedMask == 6 || fT0UsedMask == 7){
      SetTrkMaskBit(kT0_1, 1);
      // cout<<"Is two"<<endl;
    }
    if(fT0UsedMask == 4 || fT0UsedMask == 5 || fT0UsedMask == 6 || fT0UsedMask == 7){
      SetTrkMaskBit(kT0_2, 1);
      // cout<<"Is three"<<endl;
    }
    
    
    SetTrkMaskBit(kIsTOFout, fTOFout);              //Set the flag for TOFout
    SetTrkMaskBit(kIsTOFTime, fTime);               //Set the flag for Time
    SetTrkMaskBit(kIsTRDout, fTRDout);              //Set the flag for TRDout
    SetTrkMaskBit(kPassGoldenChi2, fPassGoldenChi2);//Set the flag for GoldenChi2
    
    //Fill tree or storing the information in the TClonesArray, depending on the Flags
    if(fTOFout && fTime && fTreemode){
      if(fCutmode) AnalyseCutVariation(track);
      #ifdef USETREECLASS
      AliAnTOFtrack *AnTrk = (AliAnTOFtrack*)ArrayAnTrk->ConstructedAt(ArrayAnTrk->GetEntries());
      AnTrk->ComputeDCABin(fDCAXY, fDCAZ);//Convert Impact parameters to binning
      AnTrk->fTrkMask = fTrkMask;
      AnTrk->fTPCPIDMask = fTPCPIDMask;
      AnTrk->fTrkCutMask = fTrkCutMask;
      AnTrk->fLength = fLength;
      AnTrk->fLengthRatio = fLengthRatio;
      AnTrk->fTOFTime = fTOFTime;
      AnTrk->fTOFMismatchTime = fTOFMismatchTime;
      for(Int_t ipart = 0; ipart < kExpSpecies; ipart++) AnTrk->fTOFExpTime[ipart] = fTOFExpTime[ipart];
      for(Int_t ipart = 0; ipart < kExpSpecies; ipart++) AnTrk->fTOFExpSigma[ipart] = fTOFExpSigma[ipart];
      AnTrk->fT0TrkTime = fT0TrkTime;
      AnTrk->fTOFchan = fTOFchan;
      AnTrk->fEta = fEta;
      AnTrk->fPhi = fPhi;
      AnTrk->fPt = fPt;
      AnTrk->fPTPC = fPTPC;
      AnTrk->fNTOFClusters = fNTOFClusters;
      AnTrk->fTPCSignal = fTPCSignal;
      
      #else
      fTreeTrack->Fill();
      #endif
      
    }
    
    //Numerator (the ones that pass all the cuts)
    if((fTOFout == kTRUE) && (fTime == kTRUE) && (fLength > fLengthmin) && (fTOFTime > fTOFmin) && (fTOFExpTime[AliPID::kPion]>fTOFmin) && (fTOFExpTime[AliPID::kKaon] > fTOFmin) && (fTOFExpTime[AliPID::kProton] > fTOFmin) && (fTOFTime < fTOFmax)){
      hNTrk->Fill(15);
      
      //
      //If requested
      //Histograms with T-Texp-T0 with and without mismatch
      //
      #ifdef BUILDTOFDISTRIBUTIONS// Build TOF distributions only if requested
      for(Int_t species = 0; fTOFout && fTime && species < kSpecies; species++){
        const Double_t Tdiff = fTOFTime-fT0TrkTime-fTOFExpTime[species + kpi];
        const Double_t TdiffSigma = Tdiff/fTOFExpSigma[species + kpi];
        hTOFNoYCut[fBinPtIndex][fSign][species]->Fill(Tdiff);
        
        //Rapidity cut
        if(TMath::Abs(fRapidity[species]) >= fRapidityCut) continue;
        hTOF[fBinPtIndex][fSign][species]->Fill(Tdiff);
        hTOFSigma[fBinPtIndex][fSign][species]->Fill(TdiffSigma);
        
        //Skip if the TPC signal is not for pions/kaons or protons
        if(fTPCSigma[kpi] > 5 && fTPCSigma[kK] > 5 && fTPCSigma[kp] > 5) continue;
        
        hTOFNoMismatch[fBinPtIndex][fSign][species]->Fill(Tdiff);
        hTOFSigmaNoMismatch[fBinPtIndex][fSign][species]->Fill(TdiffSigma);
      }
      #endif
      
      //
      //Matching efficiency plots
      //Numerator
      for(Int_t i = 0; i < 3; i++){//Loop on pi/k/p
        Bool_t pass = kTRUE;
        for(Int_t j = 0; j < 3; j++){
          if(i == j){
            if(TMath::Abs(fTPCSigma[kpi + j]) > 3) pass = kFALSE;
          }
          else if(TMath::Abs(fTPCSigma[kpi + j]) < 3) pass = kFALSE;
        }
        if(pass) hNumMatchTPC[fSign][i]->Fill(fPt);
      }
      
      
      if(fPdgIndex != -999){//MC info
        hNumMatchMC[fSign][fPdgIndex]->Fill(fPt);
        if(fProdInfo == 0 ){//Primaries
          hNumMatchPrimMC[fSign][fPdgIndex]->Fill(fPt);
          if(TMath::Abs(fRapidityMC) < fRapidityCut) hNumMatchPrimMCYCut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
          if(fMCTOFMatch == 0){//True match in TOF
            hNumPrimMCTrueMatch[fSignMC][fPdgIndex][fEvtMultBin]->Fill(fPtMC);
            if(TMath::Abs(fRapidityMC) < fRapidityCut){
              hNumPrimMCTrueMatchYCut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
              if((TMath::Abs(fTPCSigma[kpi]) < 5) && (TMath::Abs(fTPCSigma[kK]) < 5) && (TMath::Abs(fTPCSigma[kp]) < 5)) hNumPrimMCTrueMatchYCutTPC[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
            }
            
          }
          else if(fMCTOFMatch > 0 && fTOFSigma[kpi + fPdgIndex] < 2.){//Consistent match in TOF
            if(TMath::Abs(fRapidityMC) < fRapidityCut){
              hNumPrimMCConsistentMatchYCut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
            }
          }
        }
      }
      
      //
      //All trakcs matched to TOF
      //
      hNumMatch[fSign]->Fill(fPt);                           //-> transverse momentum
      hNumMatchEta[fSign]->Fill(fEta);                       //-> eta
      hNumMatchphiOut[fSign]->Fill(fPhiout);                 //-> Phi
      hNumMatchPtEtaPhiout[fSign]->Fill(fPt, fEta, fPhiout); //-> transverse momentum, eta and phi
      
      if(fPt > 0.5){                                         //Lower cut on the transver momentum
        hNumMatchEtaPtMa[fSign]->Fill(fEta);                 //-> eta
        hNumMatchphiOutPtMa[fSign]->Fill(fPhiout);           //-> Phi
      }
      if(fTRDout) hNumMatchTRDOut[fSign]->Fill(fPt);         //Tracks with TRDout
      else hNumMatchNoTRDOut[fSign]->Fill(fPt);              //Tracks without TRDout
      
      
      //Performance plots
      if(fPerformance){
        //TOF
        const Double_t beta = fLength / ((fTOFTime - fT0TrkTime) * CSPEED);
        hBeta->Fill(fP, beta);
        if(fNTOFClusters < 2){
          hBetaNoMismatch->Fill(fP, beta);
          if(TMath::Abs(fEta) < 0.5) hBetaNoMismatchEtaCut->Fill(fP, beta);
          if(TMath::Abs(fEta) > 0.2) hBetaNoMismatchEtaCutOut->Fill(fP, beta);
          
        }
        
        if(fEvtMult <= 30.0 && fEvtMult >= 0.0){//Central collisions
          hBetaCentral->Fill(fP, beta);
          if(fNTOFClusters < 2){
            hBetaNoMismatchCentral->Fill(fP, beta);
            if(TMath::Abs(fEta) < 0.5) hBetaNoMismatchCentralEtaCut->Fill(fP, beta);
            if(TMath::Abs(fEta) > 0.2) hBetaNoMismatchCentralEtaCutOut->Fill(fP, beta);
          }
          
        }
        
        //TPC
        hTPCdEdx->Fill(fP, fTPCSignal);
        
      }
      
      //fine plot test matching efficiency
      hPadDist->Fill(fTOFImpactDX, fTOFImpactDZ);
      AliTOFGeometry::GetVolumeIndices(fTOFchan, det);
      coord[0] = AliTOFGeometry::GetX(det);
      coord[1] = AliTOFGeometry::GetY(det);
      coord[2] = AliTOFGeometry::GetZ(det);
      hTOFDist->Fill(AliTOFGeometry::GetSector(coord), AliTOFGeometry::GetStripNumberPerSM(AliTOFGeometry::GetPlate(coord), AliTOFGeometry::GetStrip(coord)));
      hTOFResidualX->Fill(fTOFImpactDX);
      hTOFResidualZ->Fill(fTOFImpactDZ);
      hTOFChannel->Fill(fTOFchan);
      if((fP>0.9) && (fP<1.1)){//P range selected for TOF resolution measurements
        Float_t deltat = fTOFTime-fT0TrkTime-fTOFExpTime[AliPID::kPion];
        hTimeOfFlightRes->Fill(deltat);
        if(fFineTOFReso) hTimeOfFlightResFinePerEvent->Fill(deltat);
        for(Int_t i = 0; i < 3; i++){//Loop on pi/k/p
          if(TMath::Abs(fTPCSigma[kpi+i]) < 5 && TMath::Abs(fTOFSigma[kpi+i]) < 5){
            hTimeOfFlightResNoMismatch->Fill(fTOFTime-fT0TrkTime-fTOFExpTime[AliPID::kPion]);
            break;
          }
        }
        if(fT0TrkTime != 0){//No T0 Fill
          hTimeOfFlightTOFRes->Fill(deltat);
          
        }
        if((TMath::Abs(fTOFImpactDX)<1.25) && (TMath::Abs(fTOFImpactDZ)<1.75)){
          hTimeOfFlightGoodRes->Fill(deltat);
        }
      }
      
    }
    
    if(fMCmode) AnalyseMCTracks();
    
    
  } // end of track loop
  StopTimePerformance(4);
  
  #ifdef USETREECLASS
  if(fTreemode){
    StartTimePerformance(5);
    fTreeTrack->Fill();
    StopTimePerformance(5);
    
    //
    //Prepare the array for new event
    StartTimePerformance(6);
    ArrayAnTrk->Clear();
    StopTimePerformance(6);
  }
  #endif
  
  if(fFineTOFReso){
    for (Int_t bin = 1; bin <= hTimeOfFlightResFinePerEvent->GetNbinsX(); bin++) {
      hTimeOfFlightResFine->Fill(hTimeOfFlightResFinePerEvent->GetXaxis()->GetBinCenter(bin), hTimeOfFlightResFinePerEvent->GetEntries(), hTimeOfFlightResFinePerEvent->GetBinContent(bin));
    }
    hTimeOfFlightResFinePerEvent->Reset();
  }
  
  FillTimePerformance();
  
  // Post output data.
  PostAllTheData();
  AliDebug(2, "UserExec\t END");
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::InitializeTrackVar(){
  //   AliInfo("InitializeTrackVar");
  
  fTOFcls = 0x0;
  fDCAXY = -999;
  fDCAZ = -999;
  fTOFout = kFALSE;
  fTRDout = kFALSE;
  fTime = kFALSE;
  fTPCClusters = 0;
  fTPCFindClusters = 0;
  fITSClusters = 0;
  fTPCCrossedRows = -999;
  fTPCChi2 = -999;
  fTPCChi2PerNDF = -999;
  fITSChi2 = -999;
  fITSChi2PerNDF = -999;
  fGoldenChi2 = -999;
  fLengthActiveZone = -999;
  fLength = -999;
  fLengthRatio = -999;
  fSign = kFALSE;
  fTOFTime = -999;
  fTOFMismatchTime = -999;
  fTOFchan = -999;
  fEta = -999;
  fPhi = -999;
  fPt = -999;
  fTOFSignalTot = -999;
  fP = -999;
  fPTPC = -999;
  fTOFImpactDZ = -999;
  fTOFImpactDX = -999;
  fT0TrkTime = -999;
  for (Int_t species = 0 ; species < kExpSpecies; species++) {
    fTOFExpSigma[species] = -999;
    fTOFExpTime[species] = -999;
    fTOFSigma[species] = -999;
    fTOFPIDProbability[species] = -999;
    fTPCSigma[species] = -999;
  }
  for (Int_t species = 0 ; species < 3; species++) {
    fCombinedSigma[species] = -999;
  }
  fNTOFClusters = -999;
  fTOFClusters = 0x0;
  fTPCSignal = -999;
  fXout = -999;
  fYout = -999;
  fZout = -999;
  fPhiout = -999;
  fMismatch = kFALSE;
  fPassGoldenChi2 = kFALSE;
  fPassDCAxy = kFALSE;
  fPassDCAz = kFALSE;
  fITSTPCMatch = kFALSE;
  fT0UsedMask = -999;
  fBinPtIndex = -999;
  fRapidity[0] = -999;
  fRapidity[1] = -999;
  fRapidity[2] = -999;
  
  ResetTrkMaskBit();//Set all bit to 0
  ResetTrkCutMaskBit();//Set all bit to 0
  ResetTPCPIDMaskBit();//Set all bit to 0
  
  if(fMCmode) InitializeMCTrackVar();
  //   AliDebug(2, "InitializeTrackVar\t END");
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::InitializeMCTrackVar(){
  //   AliInfo("InitializeMCTrackVar");
  
  fPMC = -999;
  fPtMC = -999;
  fPhiMC = -999;
  fEtaMC = -999;
  fProdInfo = -999;
  fRapidityMC = -999;
  fPdgcode = -999;
  fPdgIndex = -999;
  fSignMC = kFALSE;
  fMCTOFMatch = -999;
  
  ResetMCTrkMaskBit();//Set all bit to 0
  
  //   AliDebug(2, "InitializeMCTrackVar\t END");
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::InitializeEventVar(){
  //   AliInfo("InitializeEventVar");
  
  //Objects
  StartTimePerformance(7);
  fMCEvt = 0x0;
  fMCStack = 0x0;
  fMultSel = 0x0;
  StopTimePerformance(7);
  
  //Variables
  StartTimePerformance(8);
  fEvtMult = -999,
  fEvtMultBin = -1,
  fVertStatus = -1;
  for(Int_t i = 0; i < 3; i++) fPrimVertex[i] = -999;
  fNContrPrimVertex = -999;
  fNMCTracks = -999;
  fMCPrimaries = -999;
  fEvtPhysSelected = kFALSE;
  fEvtSelected = kFALSE;
  StopTimePerformance(8);
  
  StartTimePerformance(9);
  ResetEvtMaskBit();//Set all bit to 0
  
  InitializeTrackVar();
  StopTimePerformance(9);
  
  //   AliDebug(2, "InitializeEventVar\t END");
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::Terminate(Option_t *){
  AliInfo("Terminate");
  // Draw result to the screen if you want
  // Called once at the end of the query
  fListHist = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListHist) {
    printf("fListHist not available\n");
    PostAllTheData();
    return;
  }
  if(hNEvt){
    hNEvt = dynamic_cast<TH1F*>(fListHist->FindObject("hNEvt"));
    AliInfo("Printing Event Stats");
    for(Int_t bin = 1; bin  <=  hNEvt->GetNbinsX(); bin++) if(!((TString) hNEvt->GetXaxis()->GetBinLabel(bin)).EqualTo("")) printf("%s = %.0f\n", hNEvt->GetXaxis()->GetBinLabel(bin), hNEvt->GetBinContent(bin));
    
    //     printf("Number of Analyzed Events = %f\n", hNEvt->GetBinContent(1));
    //     printf("Number of Events that passed the trigger = %f\n", hNEvt->GetBinContent(3));
  }else{
    printf("hNEvt not available\n");
  }
  
  if(hNTrk){
    hNTrk = dynamic_cast<TH1F*>(fListHist->FindObject("hNTrk"));
    AliInfo("Printing Track Stats");
    for(Int_t bin = 1; bin  <=  hNTrk->GetNbinsX(); bin++) if(!((TString) hNTrk->GetXaxis()->GetBinLabel(bin)).EqualTo("")) printf("%s = %.0f\n", hNTrk->GetXaxis()->GetBinLabel(bin), hNTrk->GetBinContent(bin));
  }else{
    printf("hNTrk not available\n");
  }
  
  AliDebug(2, "Terminate\t END");
  return;
  
}

//*************************************************************
//********************Utility Methods**************************
//*************************************************************

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::ComputeEvtMultiplicityBin(){
  if(fEvtMultBin != -1){
    AliFatal(Form("Multiplicity bin already assigned to value %i!", fEvtMultBin));
    return;
  }
  
  if(fHImode){
    for(Int_t multbin = 0; multbin < kEvtMultBins; multbin++){        ///<  Computes the Multiplicity bin
      if(fEvtMult < fMultiplicityBin.At(multbin) || fEvtMult >= fMultiplicityBin.At(multbin+1) ) continue;
      AliDebug(2, Form("Found bin %i/%i : %f < fEvtMultBin %f < %f", multbin, kEvtMultBins, fMultiplicityBin.At(multbin), fEvtMult, fMultiplicityBin.At(multbin+1)));
      fEvtMultBin = multbin;
      break;
    }
  }
  else{
    
    if((fEvtMult>= 0) && (fEvtMult <= 5)){fEvtMultBin = 0;}
    if((fEvtMult>= 6) && (fEvtMult <= 9)){fEvtMultBin = 1;}
    if((fEvtMult>= 10) && (fEvtMult <= 14)){fEvtMultBin = 2;}
    if((fEvtMult>= 15) && (fEvtMult <= 22)){fEvtMultBin = 3;}
    if((fEvtMult>= 23) && (fEvtMult <= 32)){fEvtMultBin = 4;}
    if(fEvtMult>= 33){fEvtMultBin = 5;}
    if((fEvtMult>= -204) && (fEvtMult <= -200)){fEvtMultBin = 6;}
    
  }
  
  if(fEvtMultBin < 0 || fEvtMultBin >= kEvtMultBins) AliFatal(Form("Multiplicity bin wrongly assigned, for Multiplicity %f fEvtMultBin value: %hd!", fEvtMult, fEvtMultBin));
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::AnalyseMCParticles(){
  //   AliInfo("AnalyseMCParticles");
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    AliError("Could not retrieve MC event handler");
    return;
  }
  
  
  if (eventHandler) fMCEvt = eventHandler->MCEvent();
  if (!fMCEvt) {
    AliError("Could not retrieve MC event");
    return;
  }
  fPIDResponse->SetCurrentMCEvent(fMCEvt);
  
  fMCStack = fMCEvt->Stack();
  if (!fMCStack) return;
  
  fNMCTracks = fMCEvt->GetNumberOfTracks();
  fMCPrimaries = fMCEvt->GetNumberOfPrimaries();
  
  if(fEvtMultBin < 0 || fEvtMultBin > kEvtMultBins -1) AliFatal("The Multiplicity bin is not defined!!!");
  
  Bool_t passeta = kTRUE;
  Bool_t passy = kTRUE;
  TParticle * trackMC = nullptr;
  //loop on primary MC tracks Before Event Selection
  for(Int_t i = 0; i < fMCStack->GetNtrack(); i++) {
    if (!fMCStack->IsPhysicalPrimary(i)) continue;//Keep only primaries
    InitializeMCTrackVar();//Initialize all the variables to zero
    passeta = kTRUE;
    passy = kTRUE;
    
    trackMC = fMCStack->Particle(i);
    fPMC = trackMC->P();
    fPtMC = trackMC->Pt();
    fEtaMC = trackMC->Eta();
    //if(TMath::Abs(treeMCEtaBis)>= 0.9){continue;}
    fPhiMC = trackMC->Phi();
    fPdgcode = trackMC->GetPdgCode();
    if(TMath::Abs(trackMC->Y()) >= fRapidityCut) passy = kFALSE;//Rapidity cut
    if(TMath::Abs(trackMC->Eta()) >= fEtaRange) passeta = kFALSE;//Eta cut
    
    //Only pi/k/p are kept
    if((TMath::Abs(fPdgcode) == 211)) fPdgIndex = 0;      //Particle is a Pion
    else if((TMath::Abs(fPdgcode) == 321)) fPdgIndex = 1; //Particle is a Kaon
    else if((TMath::Abs(fPdgcode) == 2212)) fPdgIndex = 2;//Particle is a Proton
    else continue;
    
    if(fPdgcode > 0) fSignMC = kFALSE;
    else fSignMC = kTRUE;
    
    if(passy) hDenTrkTriggerMultTrk[fSignMC][fPdgIndex]->Fill(fPtMC);
    
    if(!fEvtPhysSelected) continue;//After Physics Selection
    //vertex efficiency correction+senza taglio in eta
    
    if(passy) hDenTrkVertMultTrk[fSignMC][fPdgIndex]->Fill(fPtMC);
    
    if(!fEvtSelected) continue;//After Event Selection
    
    if(passy) hDenPrimMCYCut[fSignMC][fPdgIndex][fEvtMultBin]->Fill(fPtMC);
    
    if(passeta) hDenPrimMCEtaCut[fSignMC][fPdgIndex][fEvtMultBin]->Fill(fPtMC);
    
    if(passeta && passy) hDenPrimMCEtaYCut[fSignMC][fPdgIndex][fEvtMultBin]->Fill(fPtMC);
    
    
  }
  
  InitializeMCTrackVar();
  
  //   AliDebug(2, "AnalyseMCParticles\t END");
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::AnalyseMCTracks(){//Returns kTRUE if it reaches the end
  //   AliInfo("AnalyseMCTracks");
  
  if(TMath::Abs(fRapidityMC) >= fRapidityCut){return kFALSE;}//Rapidity cut
  
  if((fTOFout == kTRUE) && (fTime == kTRUE)) hNumMatchMultTrkInc[fSign][fEvtMultBin]->Fill(fPt);
  hDenMatchMultTrkInc[fSign][fEvtMultBin]->Fill(fPt);
  
  if(fPdgIndex < 0) return kFALSE;//Pi/K/P only!
  hDenMatchMultTrk[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
  if(fTRDout == kTRUE) hDenMatchMultTrkTRDOut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
  else hDenMatchMultTrkNoTRDOut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
  if((fTOFout == kTRUE) && (fTime == kTRUE)){
    hNumMatchMultTrk[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
    if(fTRDout == kTRUE) hNumMatchMultTrkTRDOut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
    else hNumMatchMultTrkNoTRDOut[fSign][fPdgIndex][fEvtMultBin]->Fill(fPt);
  }
  
  
  //   fTreeTrackMC->Fill();
  
  /*
  *      kTOFout = TOF matching
  *      kTIME = good integrated time
  *      100000 > track->GetTOFsignal() > 12000 = TOF time reasanble range
  *      tracklength > 365 = should be greater than the TOF radius (370 cm)
  */
  
  
  //   AliDebug(2, "AnalyseMCTracks\t END");
  return kTRUE;
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::GatherTrackMCInfo(const AliESDtrack * trk){
  //       AliDebug(2, "Gathering the track MC Info");
  const Int_t TrkLabel = trk->GetLabel(); /*The Get*Label() getters return the label of the associated MC particle. The absolute value of this label is the index of the particle within the MC fMCStack. If the label is negative, this track was assigned a certain number of clusters that did not in fact belong to this track. */
  const Int_t AbsTrkLabel = TMath::Abs(TrkLabel);
  Int_t TOFTrkLabel[3] = {-1};//This can contain three particles wich occupy the same cluster
  Int_t mfl, uniqueID;
  
  trk->GetTOFLabel(TOFTrkLabel);//Gets the labels of the tracks matched to the TOF, this can be used to remove the mismatch and to compute the efficiency! The label to check is the first one, the others can come from different tracks
  
  if (TOFTrkLabel[0] == -1) {fMCTOFMatch = -1;}                // Track was not matched to any TOF hit.
  else if (AbsTrkLabel == TOFTrkLabel[0]) {fMCTOFMatch = 0;}   // Track was correctly matched to a TOF hit.
  else {fMCTOFMatch = 1;}                                      // Track was matched to a TOF hit but comes from mismatch!
  
  TParticle *part = (TParticle*) fMCStack->Particle(AbsTrkLabel);//Particle in the stack
  if(!part){
    AliError("Cannot find the TParticle!");
    return kFALSE;
  }
  
  //Define the MC truth on the track
  fPtMC = part->Pt();
  fPhiMC = part->Phi(); //angolo tra 0 e 2pi
  fEtaMC = part->Eta();
  fPdgcode = part->GetPdgCode();
  if(fSignMC != kFALSE) AliError("fSignMC already defined!!");
  if(fPdgcode < 0) fSignMC = kTRUE;//Negative
  
  if((TMath::Abs(fPdgcode) == 211)) fPdgIndex = 0;//Track is a Pion
  else if((TMath::Abs(fPdgcode) == 321)) fPdgIndex = 1;//Track is a Kaon
  else if((TMath::Abs(fPdgcode) == 2212)) fPdgIndex = 2;//Track is a Proton
  
  fRapidityMC = ComputeY(AliPID::ParticleMass(fPdgIndex + 2));
  
  //
  //Particle production
  if(fMCStack->IsPhysicalPrimary(AbsTrkLabel)) fProdInfo = 0;//Track is Physical Primary
  else{//If it not physical primary check the Origin
    Int_t indexMoth = part->GetFirstMother();
    if(indexMoth >= 0){
      TParticle* mother = fMCStack->Particle(indexMoth);
      Float_t PDGcodemother = TMath::Abs(mother->GetPdgCode());
      mfl = Int_t (PDGcodemother/ TMath::Power(10, Int_t(TMath::Log10(PDGcodemother))));
    }
    else mfl = 0;
    
    uniqueID = part->GetUniqueID();
    
    if(mfl == 3 && uniqueID == kPDecay) fProdInfo = 1;//Track comes from weak decay
    else fProdInfo = 2;//Track is from material
  }
  
  //       AliDebug(2, Form("Track is a %i %i", fSignMC, fPdgIndex));
  
  return kTRUE;
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::SelectEvents(Int_t &binstart){
  
  if (!fEventCut.AcceptEvent(fESD)) {
    return kFALSE;
  }
  else return kTRUE;
  
  
  
  if(fHImode){//Heavy Ion
    if(fEvtMult == -999){//Multiplicity estimator not initialized
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kNoCalib){//199: kNoCalib: centrality not calibrated, this is the default value for the centrality code
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejTrigger){//200: kRejTrigger: do not pass the trigger
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejINELgtZERO){//201: kRejINELgtZERO: do not pass INEL>0 Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejVzCut){//202: kRejVzCut: do not pass vertex Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejPileupInMultBins){//203: kRejPileupInMultBins: do not pass Pile-up Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejConsistencySPDandTrackVertices){//204: kRejConsistencySPDandTrackVertices: do not pass consistency of vertex Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejTrackletsVsClusters){//205: kRejTrackletsVsClusters: do not pass Tracklets Vs Clusters Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    if(fEvtMult == AliMultSelectionCuts::kRejNonZeroNContribs){//206: kRejNonZeroNContribs: do not pass Contributors (to vertex) Cut
      
      return kFALSE;
    }
    hNEvt->Fill(binstart++);
    
    
  }
  else{//pp
    
    //------------------------------------------------
    // Selection Investigation with AliPPVsMultUtils
    //------------------------------------------------
    Bool_t passall = kFALSE;
    //------------------------------------------------
    //Step 1: Check for Min-Bias Trigger
    //------------------------------------------------
    if(AliPPVsMultUtils::IsMinimumBias(fESD)){
      hNEvt->Fill(binstart++);
      //------------------------------------------------
      //Step 2: Check for INEL>0
      //------------------------------------------------
      
      if(AliPPVsMultUtils::IsAcceptedVertexPosition(fESD)){
        hNEvt->Fill(binstart++);
        //------------------------------------------------
        //Step 3: Check for Vertex-Z position
        //------------------------------------------------
        if(AliPPVsMultUtils::IsINELgtZERO(fESD)){
          hNEvt->Fill(binstart++);
          //------------------------------------------------
          //Step 4: Check for SPD Pileup
          //------------------------------------------------
          if(AliPPVsMultUtils::IsNotPileupSPDInMultBins(fESD)){
            hNEvt->Fill(binstart++);
            //------------------------------------------------
            //Step 5: Check for SPD / track vertex consistency
            //------------------------------------------------
            if(AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fESD)){
              hNEvt->Fill(binstart++);
              passall = kTRUE;
            }
          }
        }
      }
    }
    
    //------------------------------------------------
    //Check if event is selected
    //------------------------------------------------
    if(!AliPPVsMultUtils::IsEventSelected(fESD)){
      if(passall) AliFatal("The event selection for pp is inchoerent");
      return kFALSE;
    }
    else if(!passall) AliFatal("The event selection for pp is inchoerent");
    
    
  }
  
  return kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::RunTOFChannel(){
  
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = fESD->GetTrack(iTrack);
    if (!track){continue;}
    hChannelTime->Fill(track->GetTOFCalChannel(), track->GetTOFsignal());
    
  }
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::AnalyseCutVariation(const AliESDtrack *track){
  AliInfo("Checking cut variation for track");
  if(fMCmode || !fTreemode) AliFatal("The cut variation should be called only for Tree analysis in data!");
  //Check on all the demanded cuts and fill the mask accordingly
  //At this step all tracks have already passed the loose set of cuts and therefore it is not useful to check this again!
  
  if (!track){return kFALSE;}
  Int_t index = 0;
  hCutVariation->Fill(19);//All tracks
  for(UInt_t cut = 0; cut < nCuts; cut++){
    for(UInt_t i = 1; i < CutIndex[cut]; i++){//Loose cuts are already passed -> Skipping!
      if(fCutVar[index]->AcceptTrack(track)){
        hCutVariation->Fill(index);
        AliDebug(2, Form("Track passed the cut %s (%s) at %f this will switch the bit %i !", fCutVar[index]->GetName(), fCutVar[index]->GetTitle(), CutValues[cut][i], index));
        SetTrkCutMaskBit((fTrkCutMaskIndex) index, 1);
      }
      index++;
    }
  }
  
  if(index != nCutVars) AliFatal(Form("Wrong final index (%i should be %i)", index, nCutVars));
  
  return kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::SetCutVar(){
  AliInfo("SetCutVar");
  //Cut variations Looser cuts have to be in first position!
  if(fCutmode){
    if(fMCmode || !fTreemode) AliFatal("The cut variation should be called only for Tree analysis in data!");
    if(fSimpleCutmode >= 0) AliFatal("The cut variation should be called in alternative to the simple mode!");
    
    Int_t index = 0;
    for(UInt_t cut = 0; cut < nCuts; cut++){
      for(UInt_t i = 0; i < CutIndex[cut]; i++){//Starts from zero because the first loop is to set to loos cuts the default AliESDtrackCuts
        if(index < 0 || index >= nCutVars) AliFatal("Index exceding limits");
        if(i > 0 && !fCutVar[index]) fCutVar[index] = new AliESDtrackCuts(Form("fCutVar%i", index), Form("Cut variation for %s set %i", Cuts[cut].Data(), index));
        else if(i > 0) AliFatal("AliESDtrackCuts for cut variation already exists!");
        
        AliESDtrackCuts *cuts = 0x0;
        if(i == 0){
          if(cut == kDCAxy) cuts = fESDtrackCutsPrm;
          else cuts = fESDtrackCuts;
        }
        else cuts = fCutVar[index];
        
        TString c = "";
        switch (cut) {
          case kTPCrows:
          {
            cuts->SetMinNCrossedRowsTPC(CutValues[cut][i]);
            c += Form("%f", CutValues[cut][i]);
            break;
          }
          case kTrkChi2:
          {
            cuts->SetMaxChi2PerClusterTPC(CutValues[cut][i]);
            c += Form("%f", CutValues[cut][i]);
            break;
          }
          case kDCAz:
          {
            cuts->SetMaxDCAToVertexZ(CutValues[cut][i]);
            c += Form("%f", CutValues[cut][i]);
            break;
          }
          case kDCAxy:
          {
            const TString dep = Form("%f*(%s)", CutValues[cut][i], primfunct.Data());
            cuts->SetMaxDCAToVertexXYPtDep(dep);
            c += dep;
            break;
          }
          case kGeo:
          {
            cuts->SetCutGeoNcrNcl(CutValues[cut][5*i +0], CutValues[cut][5*i +1], CutValues[cut][5*i +2], CutValues[cut][5*i +3], CutValues[cut][5*i +4]);
            c += Form("(%f, %f, %f, %f, %f)", CutValues[cut][5*i +0], CutValues[cut][5*i +1], CutValues[cut][5*i +2], CutValues[cut][5*i +3], CutValues[cut][5*i +4]);
            break;
          }
          default:
          AliFatal("Requested set for track cuts not implemented!!");
          break;
        }
        
        if(i > 0){//New cut variable
          AliInfo(Form("Creating new Cut set %s (%s) #%i for %s value #%i %s%s", cuts->GetName(), cuts->GetTitle(), index, Cuts[cut].Data(), i, c.Data(), i == CutStdIndex[cut] ? " this is the standard one!": ""));
          index++;
        }
        else AliInfo(Form("Setting Cut %s (%s) to accept loose cuts for %s value: #%i %s", cuts->GetName(), cuts->GetTitle(), Cuts[cut].Data(), i, c.Data()));
        
        
      }
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::SetTrackCuts(){
  AliInfo("SetTrackCuts");
  
  //
  // create track cuts
  //
  //   AliESDtrackCuts* AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(Bool_t selPrimaries, Int_t clusterCut)
  //For now it is experimental the PbPB2015
  //   if(fHImode) fESDtrackCuts = (AliESDtrackCuts*) AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kTRUE, 1);
  //   else fESDtrackCuts = (AliESDtrackCuts*) AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
  
  if(!fESDtrackCuts){
    fESDtrackCuts = (AliESDtrackCuts*) AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE, 1);//WARNING KEEP THE DCA SO THAT YOU CAN USE THE SECONDARIES
    // fESDtrackCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
    fESDtrackCuts->SetCutGeoNcrNcl(GeoSetStd[0], GeoSetStd[1], GeoSetStd[2], GeoSetStd[3], GeoSetStd[4]);
    
    // TPC geometrical cut combined with cut on crossed rows and number of clusters
    //****************************************************************************
    //AliESDtrackCuts::SetCutGeoNcrNcl(Float_t deadZoneWidth,Float_t cutGeoNcrNclLength, Float_t cutGeoNcrNclGeom1Pt, Float_t cutGeoNcrNclFractionNcr,  Float_t cutGeoNcrNclFractionNcl)
    //****************************************************************************
    //Float_t fDeadZoneWidth;             // width of the TPC dead zone (missing pads + PRF +ExB)
    //Float_t fCutGeoNcrNclLength;        // cut on the geometical length  condition Ngeom(cm)>cutGeoNcrNclLength default=130
    //Float_t fCutGeoNcrNclGeom1Pt;       // 1/pt dependence slope  cutGeoNcrNclLength:=fCutGeoNcrNclLength-abs(1/pt)^fCutGeoNcrNclGeom1Pt
    //Float_t fCutGeoNcrNclFractionNcr;   // relative fraction cut Ncr  condition Ncr>cutGeoNcrNclFractionNcr*fCutGeoNcrNclLength
    //Float_t fCutGeoNcrNclFractionNcl;   // ralative fraction cut Ncr  condition Ncl>cutGeoNcrNclFractionNcl
    
    
    fESDtrackCuts->SetName("MainCuts");
  }
  else AliWarning("fESDtrackCuts already exists!");
  
  //THIS WILL BE USED TO CUT ON THE PRIMARIES
  //Taken directly from the GetStandardITSTPCTrackCuts2011
  if(!fESDtrackCutsPrm){
    fESDtrackCutsPrm = new AliESDtrackCuts("fESDtrackCutsPrm", "Cut For Primaries");
    fESDtrackCutsPrm->SetMaxDCAToVertexXYPtDep(primfunct);
    fESDtrackCutsPrm->SetName("PrimaryCuts");
  }
  else AliWarning("fESDtrackCutsPrm already exists!");
  
  //Check on track cuts
  
  const TString primstring = fESDtrackCuts->GetMaxDCAToVertexXYPtDep();
  if(!primstring.IsNull()) AliFatal(Form("The track cuts inserted contains a DCAxy cut (%s)! Please set it via SetESDtrackCutsPrm!!!", primstring.Data()));
  
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::SetSimpleCutVar(){
  AliInfo("SetSimpleCutVar");
  //Cut variation with simple parameters
  if(fSimpleCutmode >= 0){
    if(fTreemode) AliFatal("The cut variation should be called only for the analysis without Tree!");
    if(fCutmode) AliFatal("The cut variation should be called in alternative to the tree mode!");
    
    AliInfo(Form("Selecting cut option %i", fSimpleCutmode));
    
    Int_t index = 0;
    Bool_t set = kFALSE;
    for(UInt_t cut = 0; cut < nCuts; cut++){
      for(UInt_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the case of requiring the same cuts as the ones which are already applied from the standard case
        if(index == fSimpleCutmode){
          switch (cut) {
            case kTPCrows:
            {
              fESDtrackCuts->SetName(Form("MainCutsVersion%i", fSimpleCutmode));
              if(CutValues[cut][i] == fESDtrackCuts->GetMinNCrossedRowsTPC()) AliFatal("Value is already in the cut!!");
              fESDtrackCuts->SetMinNCrossedRowsTPC(CutValues[cut][i]);
              AliInfo(Form("Setting SetMinNCrossedRowsTPC(%f) : %f", CutValues[cut][i], fESDtrackCuts->GetMinNCrossedRowsTPC()));
              set = kTRUE;
              break;
            }
            case kTrkChi2:
            {
              fESDtrackCuts->SetName(Form("MainCutsVersion%i", fSimpleCutmode));
              if(CutValues[cut][i] == fESDtrackCuts->GetMaxChi2PerClusterTPC()) AliFatal("Value is already in the cut!!");
              fESDtrackCuts->SetMaxChi2PerClusterTPC(CutValues[cut][i]);
              AliInfo(Form("Setting SetMaxChi2PerClusterTPC(%f) : %f", CutValues[cut][i], fESDtrackCuts->GetMaxChi2PerClusterTPC()));
              set = kTRUE;
              break;
            }
            case kDCAz:
            {
              fESDtrackCuts->SetName(Form("MainCutsVersion%i", fSimpleCutmode));
              if(CutValues[cut][i] == fESDtrackCuts->GetMaxDCAToVertexZ()) AliFatal("Value is already in the cut!!");
              fESDtrackCuts->SetMaxDCAToVertexZ(CutValues[cut][i]);
              AliInfo(Form("Setting SetMaxDCAToVertexZ(%f) : %f", CutValues[cut][i], fESDtrackCuts->GetMaxDCAToVertexZ()));
              set = kTRUE;
              break;
            }
            case kDCAxy:
            {
              fESDtrackCutsPrm->SetName(Form("PrimaryCutsVersion%i", fSimpleCutmode));
              // 	      if(CutValues[cut][i] == fESDtrackCuts->GetMaxDCAToVertexZ()) AliFatal("Value is already in the cut!!");
              
              fESDtrackCutsPrm->SetMaxDCAToVertexXYPtDep(Form("%f*(%s)", CutValues[cut][i], primfunct.Data()));
              AliInfo(Form("Setting SetMaxDCAToVertexXYPtDep(%f) : %s", CutValues[cut][i], fESDtrackCutsPrm->GetMaxDCAToVertexXYPtDep()));
              set = kTRUE;
              break;
            }
            case kGeo:
            {
              fESDtrackCuts->SetName(Form("MainCutsVersion%i", fSimpleCutmode));
              fESDtrackCuts->SetCutGeoNcrNcl(CutValues[cut][5*i +0], CutValues[cut][5*i +1], CutValues[cut][5*i +2], CutValues[cut][5*i +3], CutValues[cut][5*i +4]);
              AliInfo(Form("Setting SetCutGeoNcrNcl(%f, %f, %f, %f, %f)", CutValues[cut][5*i +0], CutValues[cut][5*i +1], CutValues[cut][5*i +2], CutValues[cut][5*i +3], CutValues[cut][5*i +4]));
              set = kTRUE;
              break;
            }
            default:
            AliFatal(Form("Cutmode %i not yet implemented!!!", fSimpleCutmode));
            break;
          }
          break;
        }
        
        index++;
      }
      if(set) break;
      
    }
    if(!set) AliFatal(Form("No variation has been applied when instead variation %i has been asked!!!", fSimpleCutmode));
  }
  
}

//________________________________________________________________________
const AliESDVertex * AliAnalysisTaskTOFSpectra::ObtainVertex(){
  //vertex status (0) vertex does not exist, (1) not enough Contributors, (2) vertex exists but outside the fiducial volume, (3) Good vertex
  const AliESDVertex* vertex = fESD->GetPrimaryVertex/*GetPrimaryVertexTracks*/(); //! Primary vertex estimated using ESD tracks
  fVertStatus = 0;
  if (!vertex){
    AliError("Cannot find any vertex");
    return 0x0;
  }
  
  
  //Get the vertex Contributors
  fNContrPrimVertex = vertex->GetNContributors();
  fVertStatus++;
  if(fNContrPrimVertex < 1) { // # of tracklets/tracks used for the estimate
    AliError(Form("Vertex has %i Contributors %f bin %i", vertex->GetNContributors(), fEvtMult, fEvtMultBin));
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD(); //! Primary vertex estimated by the SPD
    fNContrPrimVertex = vertex->GetNContributors();
    if(fNContrPrimVertex < 1)
    {
      AliError(Form("SPD Vertex has %i Contributors %f bin %i", vertex->GetNContributors(), fEvtMult, fEvtMultBin));
      return 0x0;
    }
    
  }
  
  //Get the vertex X, Y and Z coordinate
  fPrimVertex[0] = vertex->GetX();
  fPrimVertex[1] = vertex->GetY();
  fPrimVertex[2] = vertex->GetZ();
  
  //Check the position of the vertex
  fVertStatus++;
  if(TMath::Abs(fPrimVertex[2]) > 10.){
    AliError(Form("Vertex is outside the confidence window along Z : %f", fPrimVertex[2]));
    return 0x0;
  }
  
  fVertStatus++;
  return vertex;
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::SetTrackFlags(const AliESDtrack *track){
  
  //
  //Track is both ITS and TPC refit
  if ((track->GetStatus()&AliESDtrack::kITSrefit) && (track->GetStatus()&AliESDtrack::kTPCrefit)){//Track has both the match in the ITS and TPC
    fITSTPCMatch = kTRUE;
  }
  else {fITSTPCMatch = kFALSE;}
  
  //
  //kTOFout flag
  if ((track->GetStatus()&AliESDtrack::kTOFout) !=  0){//Track has the kTOFout flag
    fTOFout = kTRUE;
  }
  else {fTOFout = kFALSE;}
  
  //
  //kTIME flag
  if ((track->GetStatus()&AliESDtrack::kTIME) !=  0){//Track has the kTIME flag
    fTime = kTRUE;
  }
  else {fTime = kFALSE;}
  
  //TRDout flag
  if ((track->GetStatus()&AliESDtrack::kTRDout) !=  0) {//Track has the kTRDout flag
    fTRDout = kTRUE;
  }
  else {fTRDout = kFALSE;}
  
  //
  //Check the Golden Chi2 and set the flag
  //
  if (fGoldenChi2 < primchi2) {//Track passed the Golden Chi2 cut
    fPassGoldenChi2 = kTRUE;
  }
  else {fPassGoldenChi2 = kFALSE;}
}

//________________________________________________________________________
void AliAnalysisTaskTOFSpectra::SetTrackValues(const AliESDtrack *track, const AliESDVertex *vertex){
  
  //
  //Track eta
  fEta = track->Eta();// This function return pseudorapidity
  
  //
  //Track phi in radians
  fPhi = track->Phi();
  
  //
  //Charge sign
  if(track->GetSign() > 0) fSign = kFALSE;//Positive case
  else if(track->GetSign() < 0) fSign = kTRUE;//Negative case
  
  //
  //Track Momentum
  fPt = track->Pt();                // This function returns the track transverse momentum
  fP = track->GetP();               // This function returns the track momentum
  fPTPC = track->GetTPCmomentum();  // This function returns the track momentum in the TPC
  
  //
  //Track impact parameters
  track->GetImpactParameters(fDCAXY, fDCAZ);//Impact parameters of the track to the vertex
  fDCAXY += fDCAXYshift;//If required shift the position of the DCAxy as to correct for the wrong values of the MC
  
  //
  //Clusters in TPC
  fTPCClusters = track->GetTPCNcls();
  
  //
  //Findable Clusters in TPC
  fTPCFindClusters = track->GetTPCNclsF();
  
  //
  //TPC crossed rows
  fTPCCrossedRows = track->GetTPCCrossedRows();
  
  //
  //Clusters in ITS
  fITSClusters = track->GetITSNcls();
  
  //
  //Track length
  fLength = track->GetIntegratedLength(); // Get the total track length
  
  //
  //Track length in active zone
  fLengthActiveZone = track->GetLengthInActiveZone(0, 3, 220., fESD->GetMagneticField());
  
  //
  //TPC chi2
  fTPCChi2 = track->GetTPCchi2();
  
  //
  //TPC chi2 per NDF
  fTPCChi2PerNDF = fTPCChi2/static_cast<Float_t>(fTPCClusters);
  
  //
  //ITS chi2
  fITSChi2 = track->GetITSchi2();
  
  //
  //ITS chi2 per NDF
  fITSChi2PerNDF = fITSChi2/static_cast<Float_t>(fITSClusters);
  
  //
  //Golden Chi2
  fGoldenChi2 = track->GetChi2TPCConstrainedVsGlobal(vertex);
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::TOFCalibInitRun() {
  if(!fRecalibrateTOF) AliFatal("Requiring TOF recalibration, without the fRecalibrateTOF flag on");
  
  // check run already initialized
  if (fRunNumber == fESD->GetRunNumber()) return kTRUE;//Skip if the run number of the event analysed is the same as the one alreay used for initialization
  else if(fRunNumber == 0) AliInfo(Form("First initialization of the run %i for TOF calibration", fESD->GetRunNumber()));//First time initialization
  else AliInfo(Form("Initialization from run %i to run %i for TOF calibration", fRunNumber, fESD->GetRunNumber()));//Already previous initialization
  fRunNumber = fESD->GetRunNumber();
  
  // init cdb
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(fRunNumber);
  
  // init TOF calib
  if (!fTOFcalib->Init(fRunNumber)) {
    AliError("cannot init TOF calib");
    return kFALSE;
  }
  
  AliInfo(Form("initialized for run %d", fRunNumber));
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskTOFSpectra::TOFCalibInitEvent() {
  if(!fRecalibrateTOF) AliFatal("Requiring TOF recalibration, without the fRecalibrateTOF flag on");
  
  // init TOF-T0 maker
  fTOFT0maker->SetTimeResolution(fTimeResolution);
  
  // calibrate ESD
  fTOFcalib->CalibrateESD(fESD);
  
  // compute T0-TOF and apply it
  fTOFT0maker->ComputeT0TOF(fESD);
  
  // Write information in ESD
  fTOFT0maker->WriteInESD(fESD);
  
  return kTRUE;
  
}
