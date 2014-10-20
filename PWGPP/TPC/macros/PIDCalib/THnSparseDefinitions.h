#ifndef THNSPARSEDEFINITIONS_H
#define THNSPARSEDEFINITIONS_H

#ifndef __CINT__
#include "TCanvas.h"
#include "THnSparse.h"
#include "TPaveText.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "AliPID.h"
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
  kPidXi
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
  kPidGenXi
};

enum axesTHnSparsePIDgenYield {
  kPidGenYieldMCpid = 0,
  kPidGenYieldPt = 1,
  kPidGenYieldCentrality = 2,
  kPidGenYieldJetPt = 3,
  kPidGenYieldZ = 4,
  kPidGenYieldXi = 5,
  kPidGenYieldNumAxes = 6 
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
  kEffNumAxes = 6 
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

///*
//coarser binning at high pT and in general coarser binning to get reasonable weighting for regularisation
const Int_t nPtBins = 53;
const Double_t binsPt[nPtBins+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
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
             
const Int_t nDeDxBins = 19;
const Double_t binsDeDx[nDeDxBins+1] = {50., 52., 54., 56., 58., 60., 65., 70., 75., 80., 
                                        85., 90., 100., 120., 160., 200., 250., 300., 400., 600. };	       
/*			       
const Int_t nDeDxBins = 35;
const Double_t binsDeDx[nDeDxBins+1] = {50., 51., 52., 53., 54., 55., 56., 57., 58., 59.,
             60., 62., 64., 66., 68., 70., 72., 74., 76., 78.,
             80., 85., 90., 95., 100., 120., 140., 160., 180., 200.,
             250., 300., 350., 400., 500, 650.};*/

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
  // Remove the title from a histogram plotted in the canvase without 
  // removing the title from the histogram itself.
  // If padNum is >= 0, this method will be applied to the corresponding
  // pad number
  
  c->Update();    // Update in order to have access to the title in the following

  TPaveText* paveTextTitle = (padNum >= 0) ? (TPaveText*)c->GetPad(padNum)->FindObject("title") : (TPaveText*)c->FindObject("title");
  if (paveTextTitle) 
    paveTextTitle->Clear();
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

#endif
