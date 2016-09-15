#ifndef THNSPARSEDEFINITIONS_H
#define THNSPARSEDEFINITIONS_H

#ifndef __CINT__
#include "TCanvas.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TPaveText.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "AliPID.h"
#include "TH1.h"
#include "TAxis.h"
#endif

#define NEW_AXES

#ifdef NEW_AXES
  enum axesTHnSparseEta {
    kMCpid = 0,
    kSelectSpecies,
    kPtpcInner,
    kMultiplicity,
    kDeltaPrime,
    kEta
  };
#else
  enum axesTHnSparseEta {
    kMCpid = 0,
    kSelectSpecies,
    kPtpcInner,
    kPt,
    kDeDx,
    kMultiplicity,//kDelta,
    kDeltaPrime,
    kEta
  };
#endif

enum axesTHnSparsePID {
  kPidMCpid = 0,
  kPidSelectSpecies,
  kPidPt,
  //OLD kPidDelta,
  kPidDeltaPrime,
  kPidCentrality,
  kPidJetPt,
  kPidZ,
  kPidXi,
  kPidCharge,
  kPidTOFpidInfo,
  kPidDistance,
  kPidJt
};

/*OLD with TOF, p_TPC_Inner and p_vertex
enum axesTHnSparsePID {
  kPidMCpid = 0,
  kPidSelectSpecies,
  kPidPtpcInner,
  kPidPt,
  kPidPvertex,
  kPidDelta,
  kPidDeltaPrime,
  kPidDeltaTOF
};//*/

enum axesTHnSparsePIDgen {
  kPidGenMCpid = 0,
  kPidGenSelectSpecies,
  kPidGenPt,
  //OLD kPidGenDelta,
  kPidGenDeltaPrime,
  kPidGenCentrality,
  kPidGenJetPt,
  kPidGenZ,
  kPidGenXi,
  kPidGenCharge,
  kPidGenTOFpidInfo,
  kPidGenDistance,
  kPidGenJt
};

enum axesTHnSparsePIDgenYield {
  kPidGenYieldMCpid = 0,
  kPidGenYieldPt = 1,
  kPidGenYieldCentrality = 2,
  kPidGenYieldJetPt = 3,
  kPidGenYieldZ = 4,
  kPidGenYieldXi = 5,
  kPidGenYieldCharge = 6,
  kPidGenYieldDistance = 7,
  kPidGenYieldJt = 8,
  kPidGenYieldNumAxes = 9 
};

enum MCpid  {
  kEl = 1,
  kKa = 2,
  kMu = 3,
  kPi = 4,
  kPr = 5,
  kMuPlusPi = 10
};

enum PIDtype { 
  kMCid = 0, 
  kTPCid = 1, 
  kV0idNoTOF = 2, 
  kTPCandTOFid = 3,
  kV0idPlusTOFaccepted = 4,
  kV0idPlusTOFrejected = 5  
};

enum efficiencyAxes {
  kEffMCID = 0, 
  kEffTrackPt = 1, 
  kEffTrackEta = 2, 
  kEffTrackCharge = 3, 
  kEffCentrality = 4, 
  kEffJetPt = 5,
  kEffZ = 6,
  kEffXi = 7,
  kEffDistance = 8,
  kEffJt = 9,
  kEffNumAxes = 10
};

enum ptResolutionAxes { 
  kPtResJetPt = 0,
  kPtResGenPt = 1,
  kPtResRecPt = 2,
  kPtResCharge = 3,
  kPtResCentrality = 4, 
  kPtResNumAxes = 5
};

enum dEdxCheckAxes {
  kDeDxCheckPID = 0,
  kDeDxCheckP = 1,
  kDeDxCheckJetPt = 2,
  kDeDxCheckEtaAbs = 3,
  kDeDxCheckDeDx = 4,
  kDeDxCheckNumAxes = 5 
};
  
enum EffSteps {
  kStepGenWithGenCuts = 0, 
  kStepRecWithGenCuts = 1, 
  kStepRecWithGenCutsMeasuredObs = 2,
  kStepRecWithRecCutsMeasuredObs = 3, 
  kStepRecWithRecCutsMeasuredObsPrimaries = 4,
  kStepRecWithRecCutsMeasuredObsStrangenessScaled = 5,
  kStepRecWithRecCutsPrimaries = 6,
  kNumSteps = 7
};

enum chargeMode {
  kNegCharge = -1,
  kAllCharged = 0,
  kPosCharge = 1
};

enum TOFpidInfo {
  kNoTOFinfo = -2,
  kNoTOFpid = -1,
  kTOFpion = 0,
  kTOFkaon = 1,
  kTOFproton = 2,
  kNumTOFspecies = 3,
  kNumTOFpidInfoBins = 5
};

const TString partShortName[9] = { "El", "Ka", "Mu", "Pi", "Pr", "V0plusTOFel", "V0el", "V0pi", "V0pr" };

enum typePtBin {
  kPtBinTypeJets = 0,
  kPtBinTypePPMultOld = 1,
  kPtBinTypePPMult = 2
};


const Bool_t useNewMultiplicityRangeForMB = kTRUE;

// Coarser binning from ~5 GeV/c on to have reasonable statistics in high-mult bins
const Int_t nPtBinsPPmult = 53;
Double_t binsPtPPmult[nPtBinsPPmult + 1] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                           0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                            1.,  1.1, 1.2,  1.3, 1.4,  1.5, 1.6,  1.7, 1.8,  1.9,
                                            2.,  2.1, 2.2,  2.3, 2.4,  2.5, 2.6,  2.7, 2.8,  2.9, 
                                            3.,  3.2, 3.4,  3.6, 3.8,   4., 4.5,    5., 6.,   8., 
                                           10.,  15., 20.,  50. };

// Old binning for mult (used for thesis and results before 2015_04_21). 
// Much coarser binning at high pT and in general coarser binning to get reasonable weighting for regularisation
const Int_t nPtBinsPPmultOld = 40;
const Double_t binsPtPPmultOld[nPtBinsPPmultOld+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
             0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
             1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
             3.6, 3.8, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0,
             50.0 };

///*
//coarser binning at high pT and in general coarser binning to get reasonable weighting for regularisation
const Int_t nPtBinsJets = 53;
const Double_t binsPtJets[nPtBinsJets+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
             0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
             1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
             3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
             9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0,
             24.0, 30.0, 40.0, 50.0 };
//*/
             
/*
//coarser binning at high pT and, in addition, coarser binning around crossings for PbPb 2.76 ATeV FINAL version
const Int_t nPtBins = 50;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.95, 1.1,
             1.3, 1.4, 1.8, 2.2, 2.6, 3.0, 3.2, 3.4, 3.6, 3.8, 
             4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0, 
             50.0 };
*/



/*
//coarser binning at high pT and, in addition, coarser binning around crossings for pPb 5.023 ATeV FINAL version
const Int_t nPtBins = 52;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
             1.0, 1.1, 1.2, 1.3, 1.4, 1.8, 2.2, 2.6, 3.0, 3.2, 
             3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
             8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0,
             22.0, 30.0, 50.0 };
*/

/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 7 TeV FINAL version
const Int_t nPtBins = 47;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6,  0.65, 0.7, 0.75, 0.8, 1.2, 1.4, 1.6,
             1.8, 2.0, 2.6, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5,
             5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 
             13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 50.0 };
*/
/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 7 TeV VERY NEW version
const Int_t nPtBins = 46;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6,  0.65, 0.7, 0.75, 0.8, 1.2, 1.6,
             2.0, 2.5, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,
             5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
             14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0, 50.0 };
*/
/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 7 TeV NEW version
const Int_t nPtBins = 47;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6,  0.65, 0.7, 0.75, 0.9, 1.2, 1.4, 1.6,
             2.0, 2.5, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,
             5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
             14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0, 50.0 };
*/
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 7 TeV
/*
const Int_t nPtBins = 50;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5, 0.55, 0.6,  0.65, 0.7, 0.75, 0.9, 1.1, 1.2, 1.3,
             1.4, 1.6, 1.8, 2.0, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8,
             4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0,
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0,
             50.0 };
*/
/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 2.76 TeV FINAL version
const Int_t nPtBins = 45;
const Double_t binsPt[nPtBins+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
             0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.05,  
             1.2, 1.3, 1.4, 1.8, 2.9, 3.4, 3.6, 3.8, 4.0, 4.5, 
             5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 
             13.0, 14.0, 15.0, 16.0, 20.0, 50.0 };
*/
/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 2.76 TeV NEW version
const Int_t nPtBins = 49;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,
             1.05,  1.2,  1.3 , 1.4,  1.7,  2.4, 3.2 , 3.4 , 3.6,  3.8 ,
             4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0,
             50.0 };
*/
/*
// Coarser binning at high pT and, in addition, coarser binning around crossings for pp 2.76 TeV
const Int_t nPtBins = 50;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
             1.0,  1.1 , 1.2,  1.3 , 1.4, 1.8,  2.4, 3.2,  3.6,  3.8 ,
             4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0,
             50.0 };
*/

/*
//coarser binning at high pT
const Int_t nPtBins = 60;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
             0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
             1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
             2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
             4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 22.0, 26.0, 35.0,
             50.0 };
*/

/* OLD default as used in PID-Task to create THnSparses
const Int_t nPtBins = 68;
const Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
*/
/* "NEW" default as used in PID-task to create THnSparses (from 2015_04_21 on) - old binning is a true subset:
const Int_t nPtBins = 73;
Double_t binsPt[nPtBins + 1] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                   1.,  1.1, 1.2,  1.3, 1.4,  1.5, 1.6,  1.7, 1.8,  1.9,
                                   2.,  2.1, 2.2,  2.3, 2.4,  2.5, 2.6,  2.7, 2.8,  2.9, 
                                   3.,  3.2, 3.4,  3.6, 3.8,   4., 4.5,   5., 5.5,   6.,
                                  6.5,   7.,  8.,   9., 10.,  11., 12.,  13., 14.,  15.,
                                  16.,  18., 20.,  22., 24.,  26., 28.,  30., 32.,  34., 
                                  36.,  40., 45.,  50. };
*/

const Int_t nDeDxBins = 19;
const Double_t binsDeDx[nDeDxBins+1] = {50., 52., 54., 56., 58., 60., 65., 70., 75., 80., 
                                        85., 90., 100., 120., 160., 200., 250., 300., 400., 600. };	       
/*			       
const Int_t nDeDxBins = 35;
const Double_t binsDeDx[nDeDxBins+1] = {50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
             60., 62., 64., 66., 68., 70., 72., 74., 76., 78.,
             80., 85., 90., 95., 100., 120., 140., 160., 180., 200.,
             250., 300., 350., 400., 500, 650.};*/

//_____________________________________________________________________________________________________________________________
Bool_t CentralityHasDecimalsPlaces(Double_t cent)
{
  // Check, whether the centrality has decimals places or is an integer
  const Double_t temp1 = ((Int_t)cent)*1e6;
  const Double_t temp2 = cent*1e6;
  
  return TMath::Abs(temp1 - temp2) > 1;
}


//____________________________________________________________________________________________________________________
Int_t getLineColor(Int_t ID) {
  switch (ID) {
    case kEl:
      // El
      return kMagenta;
    case kKa:
      // Ka
      return kGreen;
    case kMu:
      // Mu
      return kOrange -3;
    case kPi:
      // Pi
      return kRed;
    case kPr:
      // Pr
      return kBlue;
    case kMuPlusPi:
      // Muons plus pions
      return kCyan;
    default:
      return 0;
  }
  
  return 0;
}


//____________________________________________________________________________________________________________________
Int_t getLineColorAliPID(Int_t ID) {
  switch (ID) {
    case AliPID::kElectron:
      // El
      return kMagenta;
    case AliPID::kKaon:
      // Ka
      return kGreen;
    case AliPID::kMuon:
      // Mu
      return kOrange - 3;
    case AliPID::kPion:
      // Pi
      return kRed;
    case AliPID::kProton:
      // Pr
      return kBlue;
    default:
      return 0;
  }
  
  return 0;
}


//____________________________________________________________________________________________________________________
void ClearTitleFromHistoInCanvas(TCanvas* c, Int_t padNum = -1)
{
  // Remove the title from a histogram plotted in the canvas without 
  // removing the title from the histogram itself.
  // If padNum is >= 0, this method will be applied to the corresponding
  // pad number
  
  c->Update();    // Update in order to have access to the title in the following

  TPaveText* paveTextTitle = (padNum >= 0) ? (TPaveText*)c->GetPad(padNum)->FindObject("title") : (TPaveText*)c->FindObject("title");
  if (paveTextTitle) 
    paveTextTitle->Clear();
  
  c->Modified();
}


//____________________________________________________________________________________________________________________
void ClearTitleFromHistoInPad(TPad* c)
{
  // Remove the title from a histogram plotted in the pad without 
  // removing the title from the histogram itself.
  
  c->Update();    // Update in order to have access to the title in the following

  TPaveText* paveTextTitle = (TPaveText*)c->FindObject("title");
  if (paveTextTitle) 
    paveTextTitle->Clear();
  
  c->Modified();
}


//____________________________________________________________________________________________________________________
Int_t GetAxisByTitle(const THnSparse* h, TString title)
{
  if (!h)
    return -1;
  
  for (Int_t iDim = 0; iDim < h->GetNdimensions(); iDim++) {
    if (!title.CompareTo(h->GetAxis(iDim)->GetTitle()))
      return iDim;
  }

  return -1;
}


//____________________________________________________________________________________________________________________
Double_t GetErrorRatioIndependent(Double_t num, Double_t numErr, Double_t den, Double_t denErr)
{
  if (num <= 0)
    return 0.;
  
  const Double_t den2 = den*den;
  return TMath::Sqrt((numErr*numErr * den2 + denErr*denErr * num*num) / (den2*den2));
}



//____________________________________________________________________________________________________________________
const Double_t* GetPtBins(Int_t type, Int_t& nBins)
{
  switch (type) {
    case kPtBinTypeJets:
      nBins = nPtBinsJets;
      return binsPtJets;
      break;
    case kPtBinTypePPMultOld:
      nBins = nPtBinsPPmultOld;
      return binsPtPPmultOld;
      break;
    case kPtBinTypePPMult:
      nBins = nPtBinsPPmult;
      return binsPtPPmult;
      break;
    default:
      break;
  }
  
  printf("ERROR: Requested unknown bin type!");
  nBins = -999;
  
  return 0x0;
}


//_____________________________________________________________________________________________________________________________
TH1D* GetRelErrorHisto(TH1D* h, TString name, TString yAxisTitle = "Rel. sys. error")
{
  if (!h)
    return 0x0;
  
  TH1D* hRelErr = new TH1D(*h);
  hRelErr->SetName(name.Data());
  hRelErr->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    const Double_t val = h->GetBinContent(i);
    const Double_t err = h->GetBinError(i);
    
    Double_t relErr = 0.;
    if (val > 0) 
      relErr = err / val;
    
    hRelErr->SetBinContent(i, relErr);
  }
  
  return hRelErr;
}


//_____________________________________________________________________________________________________________________________
TGraphAsymmErrors* GetRelErrorGraph(TGraphAsymmErrors* gr, TString name, TString yAxisTitle = "Rel. sys. error")
{
  if (!gr)
    return 0x0;
  
  TGraphAsymmErrors* grRelErr = new TGraphAsymmErrors(*gr);
  grRelErr->SetName(name.Data());
  grRelErr->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  for (Int_t i = 0; i < gr->GetN(); i++) {
    const Double_t val = gr->GetY()[i];
    const Double_t errLow = gr->GetEYlow()[i];
    const Double_t errHigh = gr->GetEYhigh()[i];
    
    Double_t relErrLow = 0.;
    Double_t relErrHigh = 0.;
    if (val > 0) {
      relErrLow = errLow / val;
      relErrHigh = errHigh / val;
    }
    
    grRelErr->SetPointEYlow(i, relErrLow);
    grRelErr->SetPointEYhigh(i, relErrHigh);
  }
  
  return grRelErr;
}


//____________________________________________________________________________________________________________________
const Char_t* GetLatexNameParticleAntiParticle(Int_t species)
{
  switch (species) {
    case AliPID::kElectron:
      return "e^{+}+e^{-}";
      break;
    case AliPID::kMuon:
      return "#mu^{+}+#mu^{-}";
      break;
    case AliPID::kPion:
      return "#pi^{+}+#pi^{-}";
      break;
    case AliPID::kKaon:
      return "K^{+}+K^{-}";
      break;
    case AliPID::kProton:
      return "p+#bar{p}";
      break;
  }
  
  return "UNKNOWN SPECIES";
}


//____________________________________________________________________________________________________________________
const Char_t* GetLatexNamePosParticle(Int_t species)
{
  switch (species) {
    case AliPID::kElectron:
      return "e^{+}";
      break;
    case AliPID::kMuon:
      return "#mu^{+}";
      break;
    case AliPID::kPion:
      return "#pi^{+}";
      break;
    case AliPID::kKaon:
      return "K^{+}";
      break;
    case AliPID::kProton:
      return "p";
      break;
  }
  
  return "UNKNOWN SPECIES";
}


//____________________________________________________________________________________________________________________
const Char_t* GetLatexNameNegParticle(Int_t species)
{
  switch (species) {
    case AliPID::kElectron:
      return "e^{-}";
      break;
    case AliPID::kMuon:
      return "#mu^{-}";
      break;
    case AliPID::kPion:
      return "#pi^{-}";
      break;
    case AliPID::kKaon:
      return "K^{-}";
      break;
    case AliPID::kProton:
      return "#bar{p}";
      break;
  }
  
  return "UNKNOWN SPECIES";
}


//____________________________________________________________________________________________________________________
const Char_t* GetLatexNameParticleAndOrAntiParticle(Int_t species, chargeMode charge, Bool_t addBrackets = kTRUE)
{
  if (charge == kPosCharge)
    return GetLatexNamePosParticle(species);
  else if (charge == kNegCharge)
    return GetLatexNameNegParticle(species);
  
  return Form("%s%s%s", addBrackets ? "(" : "", GetLatexNameParticleAntiParticle(species), addBrackets ? ")" : "");
}


//_____________________________________________________________________________________________________________________________
TH1D* GraphToHist(TGraphAsymmErrors* gr, TH1D* histRefBinning, TString histName, Bool_t checkConsistency)
{
  // Convert the graph "gr" to a histogram with the same bins as "histRefBinning".
  // If the consistency check is enabled and graph and histogram have different binning,
  // 0x0 is returned.
  
  
  if (!gr || !histRefBinning)
    return 0x0;
  
  // Check consistency
  if (checkConsistency) {
    for (Int_t i = 0; i < gr->GetN(); i++) {
      const Double_t x = gr->GetX()[i];
      const Double_t xLow = x - gr->GetEXlow()[i];
      const Double_t xHigh = x + gr->GetEXhigh()[i];
      
      const Int_t bin = histRefBinning->GetXaxis()->FindFixBin(x);
      const Double_t xBinCentre = histRefBinning->GetXaxis()->GetBinCenter(bin);
      const Double_t xBinLow = histRefBinning->GetXaxis()->GetBinLowEdge(bin);
      const Double_t xBinHigh = histRefBinning->GetXaxis()->GetBinUpEdge(bin);
      
      if (TMath::Abs(x - xBinCentre) > 1e-5 ||
          TMath::Abs(xLow - xBinLow) > 1e-5 ||
          TMath::Abs(xHigh - xBinHigh) > 1e-5) {
        printf("MultiplyHistAndGraph error: binning is not consistent; gr: %e (%e - %e) vs. hist: %e (%e - %e)!\n", x, 
               xLow, xHigh, xBinCentre, xBinLow, xBinHigh);
        return 0x0;
      }
    }
  }
  
  TH1D* hist = new TH1D(*histRefBinning);
  hist->Reset();
  hist->SetName(histName.Data());
  
  for (Int_t i = 0; i < gr->GetN(); i++) {
    const Double_t x = gr->GetX()[i];
    const Int_t bin = hist->GetXaxis()->FindFixBin(x);
    
    const Double_t grY = gr->GetY()[i];
    const Double_t grEYlow = gr->GetEYlow()[i];
    const Double_t grEYhigh = gr->GetEYhigh()[i];
    
    const Double_t grEY = (grEYlow + grEYhigh) / 2.;
    // Asymmetric errors not available for histogram - print warning!
    if (TMath::Abs(grEYlow - grEYhigh) > 1e-10)
      printf("Warning: asymmetric errors (%e vs. %e). Using mean...\n", grEYlow, grEYhigh);
    
    hist->SetBinContent(bin, grY);
    hist->SetBinError(bin, grEY);
  }
  
  return hist;
}


//____________________________________________________________________________________________________________________
TGraphAsymmErrors* HistToGraph(const TString grname, const TH1 *hh, const Double_t thres=0, const TH1 *herr=0x0, const Double_t xmin=-1e10, const Double_t xmax=1e10)
{
  if (!hh)
    return 0x0;
  
  const Int_t nbin = hh->GetNbinsX();
  Double_t xxs[nbin], yys[nbin], exs[nbin], eys[nbin];
  Int_t np=0;
  for(Int_t ii=1; ii<=nbin; ii++){
    const Double_t iyy = hh->GetBinContent(ii);
    if(iyy<=thres)
      continue;

    const Double_t iey = hh->GetBinError(ii);
    if(iey<1e-15){
      if(iyy>1e-15){
        printf("HistToGraph warning! should be fine if this is ratio %d %e %e\n", ii, iyy, iey); //exit(1);
      }
      //continue;
    }

    const Double_t ixx = hh->GetBinCenter(ii);
    if(ixx<xmin || ixx>xmax){
      //printf("test HistToGraph rejecting ixx %e xmin %e xmax %e\n", ixx, xmin, xmax);
      continue;
    }

    Double_t iex = 0;
    if(herr){
      iex = herr->GetBinContent(herr->GetXaxis()->FindBin(ixx));
    }
    else{
      iex = hh->GetBinWidth(ii)/2.;
    }

    xxs[np] = ixx;
    yys[np] = iyy;
    exs[np] = iex;
    eys[np] = iey;
    np++;
  }
  TGraphAsymmErrors * gr = new TGraphAsymmErrors(np, xxs, yys, exs, exs, eys, eys);
  gr->SetName(grname);
  gr->SetMaximum(hh->GetMaximum());
  gr->SetMinimum(hh->GetMinimum());
  gr->GetXaxis()->SetLimits(hh->GetXaxis()->GetXmin(), hh->GetXaxis()->GetXmax());
  
  gr->SetLineColor(hh->GetLineColor());
  gr->SetMarkerColor(hh->GetMarkerColor());
  gr->SetFillStyle(hh->GetFillStyle());
  gr->GetXaxis()->SetTitle(hh->GetXaxis()->GetTitle());
  gr->GetYaxis()->SetTitle(hh->GetYaxis()->GetTitle());
  return gr;
}


//____________________________________________________________________________________________________________________
void scaleGraph(TGraphAsymmErrors* g, Double_t scaleFactor, Double_t eps = 1e-5) {
  if (TMath::Abs(scaleFactor - 1.) < eps)
    return;
  
  if (!g)
    return;
  
  for (Int_t i = 0; i < g->GetN(); i++) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i] * scaleFactor);
    g->SetPointEYlow(i,  g->GetEYlow()[i] * scaleFactor);
    g->SetPointEYhigh(i, g->GetEYhigh()[i] * scaleFactor);
  }
}


//____________________________________________________________________________________________________________________
void setAxisTitlesItalic(TAxis* axis)
{
  if (!axis)
    return;
  
  TString temp = axis->GetTitle();
  temp = temp.ReplaceAll("p_", "#it{p}_");
  temp = temp.ReplaceAll("eV/c", "eV/#it{c}");
  temp = temp.ReplaceAll("z", "#it{z}^{ch}");
  temp = temp.ReplaceAll("#xi", "#it{#xi}^{ch}");
  temp = temp.ReplaceAll("N", "#it{N}");
  
  axis->SetTitle(temp.Data());
}


//____________________________________________________________________________________________________________________
void setAxisTitlesItalic(TH1* h)
{
  if (!h)
    return;
  
  TAxis* axis = 0;
  for (Int_t i = 0; i < 2; i++) {
    if (i == 0)
      axis = h->GetXaxis();
    else
      axis = h->GetYaxis();
    
    setAxisTitlesItalic(axis);
  }
}



//_____________________________________________________________________________________________________________________________
void SetErrorsGraphFromRelErrGraph(TGraphAsymmErrors* gr, TGraphAsymmErrors* grRelErrors)
{
  if (!gr || !grRelErrors) {
    printf("SetErrorsGraphFromRelErrGraph: Missing graph(s)!\n");
    return;
  }
  
  // Every point in gr needs a point in grRelErrors, but not vice-versa!
  if (gr->GetN() > grRelErrors->GetN()) {
    printf("SetErrorsGraphFromRelErrGraph: Graphs not consistent (different number of points)!\n");
    return;
  }
  
  for (Int_t i = 0; i < gr->GetN(); i++) {
    const Double_t val = gr->GetY()[i];
    
    if (val <= 0)
      continue;
    
    const Double_t x = gr->GetX()[i];
    const Double_t exLow = gr->GetEXlow()[i];
    const Double_t exHigh = gr->GetEXhigh()[i];
    
    Int_t iRelErrors = -1;
    
    for (Int_t j = 0; j < grRelErrors->GetN(); j++) {
      if (TMath::Abs(x - grRelErrors->GetX()[j]) < 1e-5 &&
          TMath::Abs(exLow - grRelErrors->GetEXlow()[j]) < 1e-5 &&
          TMath::Abs(exHigh - grRelErrors->GetEXhigh()[j]) < 1e-5) {
        iRelErrors = j;
        break;
      }
    }
    
    if (iRelErrors < 0) {
      printf("SetErrorsGraphFromRelErrGraph: Graphs not consistent (different different binning)!\n");
      return;
    }
    
    const Double_t relErrLow = grRelErrors->GetEYlow()[iRelErrors];
    const Double_t relErrHigh = grRelErrors->GetEYhigh()[iRelErrors];
    
    Double_t errLow = 0.;
    Double_t errHigh = 0.;
    if (val > 0) {
      errLow = relErrLow * val;
      errHigh = relErrHigh * val;
    }
    
    gr->SetPointEYlow(i, errLow);
    gr->SetPointEYhigh(i, errHigh);
  }
}


//_____________________________________________________________________________________________________________________________
void GetAxisRangeForMultiplicityAxisForMB(TAxis* axis, Int_t& low, Int_t& high, Bool_t noActionForOldRange = kTRUE)
{
  // Get range for the multiplicity axis for the MB range
  // Already checked: Works also for Integral(...)!
  
  // Old definition: MB = full range
  if (!useNewMultiplicityRangeForMB) {
    if (!noActionForOldRange) {
      low = 0;
      high = -1;
    }
  }
  else {
    // Get range to mult < 0 (including underflow)
    low = axis->FindFixBin(-9999);
    if (low == 0)
      low = -1; // To include underflow bin;
    high = axis->FindFixBin(-1e-5);
  }
}


//_____________________________________________________________________________________________________________________________
void PrintSettingsAxisRangeForMultiplicityAxisForMB()
{
  // Print which range is used for MB for the mult axis
  if (useNewMultiplicityRangeForMB) 
    printf("NEW range definition for the multiplicity in the MB case....\n");
  else
    printf("OLD range definition for the multiplicity in the MB case....\n");
}


//_____________________________________________________________________________________________________________________________
TH1D* rebinHistToRefHist(TH1D* h, TH1D* hRef, TString histName, Bool_t statError, Double_t xThreshold = -1, Bool_t normalisedTodPt = kTRUE,
                         Bool_t normalisedToPt = kTRUE, Double_t EPSILON = 1e-5)
{
  // If the binning of hRef is a real subset of that in h (comparison to the level of EPSILON) above some threshold xThreshold,
  // a pointer to a new histogram is returned which has the same binning above (below) xThreshold as hRef (h),
  // but with merged content from h. If xThreshold < 0, no threshold is applied.
  // A possible bin-dependent normalisation is undone before the "merging" and later re-applied.
  // Statistical errors are added in quadrature, systematic error are averaged.
  
  if (!h || !hRef)
    return 0x0;
  
  Int_t numBinsResult = 0;
  Double_t* resultBins = 0x0;
  // xThreshold must be a bin edge for both histograms
  if (xThreshold >= 0) {
    const Int_t binThrRef = hRef->GetXaxis()->FindFixBin(xThreshold);
    const Int_t binThr = h->GetXaxis()->FindFixBin(xThreshold);

    const Double_t binLowEdgeRef = hRef->GetXaxis()->GetBinLowEdge(binThrRef);
    const Double_t binLowEdge = h->GetXaxis()->GetBinLowEdge(binThr);
    
    if (TMath::Abs(xThreshold - binLowEdge) > EPSILON || TMath::Abs(xThreshold - binLowEdgeRef) > EPSILON) {
      printf("xThreshold is not a bin edge and cannot be used for rebinning!\n");
      return 0x0;
    }
    
    // Now construct the new axis
    numBinsResult = (binThr - 1) + (hRef->GetNbinsX() - (binThrRef - 1));
    resultBins = new Double_t[numBinsResult + 1];
    
    for (Int_t i = 1; i < binThr; i++)
      resultBins[i - 1] = h->GetXaxis()->GetBinLowEdge(i);
    for (Int_t i = binThr, j = binThrRef; i <= numBinsResult; i++, j++)
      resultBins[i - 1] = hRef->GetXaxis()->GetBinLowEdge(j);
    resultBins[numBinsResult] = hRef->GetXaxis()->GetBinUpEdge(hRef->GetNbinsX());
  }
  else {
    // Just take the binning from hRef
    numBinsResult = hRef->GetNbinsX();
    resultBins = new Double_t[numBinsResult + 1];
    
    for (Int_t i = 1; i <= numBinsResult; i++)
      resultBins[i - 1] = hRef->GetXaxis()->GetBinLowEdge(i);
    resultBins[numBinsResult] = hRef->GetXaxis()->GetBinUpEdge(numBinsResult);
  }
  
  
  // Create the new histo with the desired binning
  TH1D* hResult = new TH1D(*hRef);
  hResult->SetName(histName.Data());
  hResult->Reset();
  hResult->SetBins(numBinsResult, resultBins);
  
  delete resultBins;
  
  
  // Check, if binning is consistent (i.e. whether each bin edge of hResult has a corresponding one in h)
  for (Int_t binRef = 1; binRef <= hResult->GetNbinsX(); binRef++) {
    const Double_t binLowEdgeRef = hResult->GetXaxis()->GetBinLowEdge(binRef);
    
    const Int_t bin = h->GetXaxis()->FindFixBin(binLowEdgeRef);
    const Double_t binLowEdge = h->GetXaxis()->GetBinLowEdge(bin);
    if (TMath::Abs(binLowEdge - binLowEdgeRef) > EPSILON) {
      printf("Histograms not consistent. Cannot rebin!\n");
      delete hResult;
      return 0x0;
    }
  }
  
  // Also check the last upper bin edge.
  const Double_t binUpEdgeRef = hResult->GetXaxis()->GetBinUpEdge(hResult->GetNbinsX());
  
  // Note that FindFixBin returns the next higher bin if the value is exactly on the edge. If the last bin has the same
  // upper edge, it will return the overflow bin. Just subtract 1 in any case.
  const Int_t bin = h->GetXaxis()->FindFixBin(binUpEdgeRef) - 1;
  const Double_t binUpEdge = h->GetXaxis()->GetBinUpEdge(bin);
  if (TMath::Abs(binUpEdge - binUpEdgeRef) > EPSILON) {
    
    // It may happen that the ref histo has a larger range. This is OK, since the content can then just be set to zero.
    const Int_t binCheck = hResult->GetXaxis()->FindFixBin(binUpEdge);
      if (binCheck > hResult->GetNbinsX()) {
      printf("Histograms not consistent (last edge). Cannot rebin!\n");
      delete hResult;
      return 0x0;
    }
  }
  
  // Histograms are consistent! Do the merging
  for (Int_t bin = 1; bin <= hResult->GetNbinsX(); bin++) {
    Double_t cont = 0.;
    Double_t err = 0.;
    
    // Add up all bins from h that are within the range of the bin of hResult. Since the bin edges have already been checked against
    // consistency, it is safe to check only if the bin centre is within the boundaries.
    const Double_t xLowBound = hResult->GetXaxis()->GetBinLowEdge(bin);
    const Double_t xUpBound = hResult->GetXaxis()->GetBinUpEdge(bin);
    
    for (Int_t binToSum = 1; binToSum <= h->GetNbinsX(); binToSum++) {
      const Double_t xCentre = h->GetXaxis()->GetBinCenter(binToSum);
      
      if (xCentre < xLowBound)
        continue;
      
      if (xCentre > xUpBound)
        break;
      
      Double_t cont_tmp = h->GetBinContent(binToSum);
      Double_t err_tmp = h->GetBinError(binToSum);
      
      // Undo the bin-dependent normalisation, if neccessary
      Double_t undoNormFactor = 1.;
      
      if (normalisedToPt)
        undoNormFactor *= xCentre;
      
      if (normalisedTodPt)
        undoNormFactor *= h->GetXaxis()->GetBinWidth(binToSum);
      
      cont_tmp *= undoNormFactor;
      err_tmp *= undoNormFactor;
      
      cont += cont_tmp;
      
      // Stat. error: add errors quadratically; sys error: add errors
      if (statError)
        err += err_tmp * err_tmp;
      else
        err += err_tmp;
    }
    
    if (statError)
      err = TMath::Sqrt(err);
      
    // Re-apply normalisation, if neccessary
    Double_t normFactor = 1.;
    
    if (normalisedToPt)
      normFactor /= hResult->GetXaxis()->GetBinCenter(bin);
    
    if (normalisedTodPt)
      normFactor /= hResult->GetXaxis()->GetBinWidth(bin);
    
    cont *= normFactor;
    err *= normFactor;
    
    hResult->SetBinContent(bin, cont);
    hResult->SetBinError(bin, err);
  }
  
  return hResult;
}


//_____________________________________________________________________________________________________________________________
Double_t smartInterpolationPol2(Double_t xLow, Double_t x, Double_t xHigh, TH1D* hden, Double_t xHighThreshold = -999)
{
  // Instead of doin an interpolation, fit locally and integrate. xLow, x and xHigh are the lower bin edge, the centre and the
  // upper bin edge, respectively. hden is the histogram for the denominator. Return -999 in case of error.
  // If xHighThreshold is set, nothing will be interpolated above a certain value.
  
  const Double_t errValue = -999;
  // Reject invalid ranges
  if (xHigh <= xLow || x < xLow || !hden || (xHighThreshold > -990 && x >= xHighThreshold))
    return errValue;
  
  Int_t ibin = hden->GetXaxis()->FindFixBin(x);
  
  // If the binning is the same, no need for doing fitting and integration.
  if (TMath::Abs(hden->GetXaxis()->GetBinCenter(ibin) - x) < 1e-5 &&
      TMath::Abs(hden->GetXaxis()->GetBinLowEdge(ibin) - xLow) < 1e-5 &&
      TMath::Abs(hden->GetXaxis()->GetBinUpEdge(ibin) - xHigh) < 1e-5) {
    return hden->GetBinContent(ibin);
  }
  
  Int_t ibinLast = hden->GetNbinsX();
  if (xHighThreshold > -990) {
    Int_t iTemp = hden->GetXaxis()->FindFixBin(xHighThreshold);
    ibinLast = TMath::Min(iTemp, ibinLast);
  }
  
  // Do not extrapolate for central value
  if (x < hden->GetXaxis()->GetBinCenter(1) || x > hden->GetXaxis()->GetBinCenter(ibinLast))
    return errValue;
  
  
  // Check if a fit is possible
  
  // Integration normally from previous to next bin
  Int_t ibinIntLow = ibin - 1;
  Int_t ibinInt = ibin;
  Int_t ibinIntHigh = ibin + 1;
  
  // If central point in very first bin and right of the bin centre (checked above), then set integration window
  // from this bin to the next two
  if (ibinIntLow < 1) {
    ibinIntLow = 1;
    ibinInt = 2;
    ibinIntHigh = 3;
  }
  // If central point in very last bin and left of the bin centre (checked above), then set integration window
  // from this previous two bins to this bin
  else if (ibinIntHigh > ibinLast) {
    ibinIntHigh = ibinLast;
    ibinInt = ibinIntHigh - 1;
    ibinIntLow = ibinIntHigh - 2;
  }
  
  if (hden->GetBinContent(ibinIntLow)  < 0. ||
      hden->GetBinContent(ibinInt)     < 0. ||
      hden->GetBinContent(ibinIntHigh) < 0.) {
    return errValue;
  }
  
  TF1* f = (TF1*)gROOT->GetFunction("pol2");
  hden->Fit(f, "0QNWCI", "", hden->GetXaxis()->GetBinLowEdge(ibinIntLow), hden->GetXaxis()->GetBinUpEdge(ibinIntHigh));
  f->SetRange(TMath::Min(xLow, hden->GetXaxis()->GetBinLowEdge(ibinIntLow)), TMath::Max(xHigh, hden->GetXaxis()->GetBinUpEdge(ibinIntHigh)));
  Double_t ref = f->Integral(xLow, xHigh) / (xHigh - xLow);
  
  //printf("%e (%f), %e (%f), %e (%f) -> %e / %e = %e; ", hden->GetBinContent(ibin - 1), xLow, hden->GetBinContent(ibin), x, hden->GetBinContent(ibin + 1), xHigh, f->Integral(xLow, xHigh), (xHigh - xLow), ref); 
  //printf("%f, %f, %f\n", hden->GetXaxis()->GetBinLowEdge(ibin), hden->GetXaxis()->GetBinCenter(ibin), hden->GetXaxis()->GetBinUpEdge(ibin));
  
  return ref;
}


//_____________________________________________________________________________________________________________________________
Double_t smartInterpolationPowerLaw(Double_t xLow, Double_t x, Double_t xHigh, TH1D* hden, Double_t xHighThreshold = -999)
{
  // Instead of doin an interpolation, fit locally and integrate. xLow, x and xHigh are the lower bin edge, the centre and the
  // upper bin edge, respectively. hden is the histogram for the denominator. Return -999 in case of error.
  // If xHighThreshold is set, nothing will be interpolated above a certain value.
  
  const Double_t errValue = -999;
  // Reject invalid ranges
  if (xHigh <= xLow || x < xLow || !hden || (xHighThreshold > -990 && x >= xHighThreshold))
    return errValue;
  
  Int_t ibin = hden->GetXaxis()->FindFixBin(x);
  
  // If the binning is the same, no need for doing fitting and integration.
  if (TMath::Abs(hden->GetXaxis()->GetBinCenter(ibin) - x) < 1e-5 &&
      TMath::Abs(hden->GetXaxis()->GetBinLowEdge(ibin) - xLow) < 1e-5 &&
      TMath::Abs(hden->GetXaxis()->GetBinUpEdge(ibin) - xHigh) < 1e-5) {
    return hden->GetBinContent(ibin);
  }
  
  Int_t ibinLast = hden->GetNbinsX();
  if (xHighThreshold > -990) {
    Int_t iTemp = hden->GetXaxis()->FindFixBin(xHighThreshold);
    ibinLast = TMath::Min(iTemp, ibinLast);
  }
  
  // Do not extrapolate for central value
  if (x < hden->GetXaxis()->GetBinCenter(1) || x > hden->GetXaxis()->GetBinCenter(ibinLast))
    return errValue;
  
  
  // Check if a fit is possible
  
  // Integration normally from previous to next bin
  Int_t ibinIntLow = 0;
  Int_t ibinIntHigh = 0;
  
  if (x < hden->GetXaxis()->GetBinCenter(ibin)) {
    ibinIntLow = ibin - 1;
    ibinIntHigh = ibin;
  }
  else {
    ibinIntLow = ibin;
    ibinIntHigh = ibin + 1;
  }
  
  if (hden->GetBinContent(ibinIntLow)  < 0. ||
      hden->GetBinContent(ibinIntHigh) < 0.)
    return errValue;
  
  /*
  TF1* f = new TF1("f", "[0]*TMath::Power(x, [1])", TMath::Min(xLow, hden->GetXaxis()->GetBinLowEdge(ibinIntLow)),
                   TMath::Max(xHigh, hden->GetXaxis()->GetBinUpEdge(ibinIntHigh)));
  f->SetParLimits(1, -20, 20);
  f->SetParError(1, 0.5); // Step size!
  */
  TF1* f = new TF1("f", "expo", TMath::Min(xLow, hden->GetXaxis()->GetBinLowEdge(ibinIntLow)),
                   TMath::Max(xHigh, hden->GetXaxis()->GetBinUpEdge(ibinIntHigh)));
  hden->Fit(f, "0QBNCI", "", hden->GetXaxis()->GetBinLowEdge(ibinIntLow), hden->GetXaxis()->GetBinUpEdge(ibinIntHigh));
  Double_t ref = f->Integral(xLow, xHigh) / (xHigh - xLow);
  delete f;
  
  //printf("%e (%f), %e (%f), %e (%f) -> %e / %e = %e; ", hden->GetBinContent(ibin - 1), xLow, hden->GetBinContent(ibin), x, hden->GetBinContent(ibin + 1), xHigh, f->Integral(xLow, xHigh), (xHigh - xLow), ref); 
  //printf("%f, %f, %f\n", hden->GetXaxis()->GetBinLowEdge(ibin), hden->GetXaxis()->GetBinCenter(ibin), hden->GetXaxis()->GetBinUpEdge(ibin));
  
  return ref;
}

#endif
