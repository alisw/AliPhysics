#include "TString.h"
#include "TMath.h"
#include "AliForwardFlowRun2Settings.h"
#include "TFile.h"

//________________________________________________________________________
AliForwardFlowRun2Settings::AliForwardFlowRun2Settings() :
  fDataType(-1),
  fPhiAcceptanceLowEdge(0),
  fPhiAcceptanceUpEdge(2*TMath::Pi()),
  fEtaLowEdge(-4.0),
  fEtaUpEdge(6.0),
  fNPhiBins(20),
  fZVtxAcceptanceLowEdge(-10),
  fZVtxAcceptanceUpEdge(10),
  fNZvtxBins(10),
  qctype("std"),
  fnoSamples(10),
  fNRefEtaBins(1),
  fNDiffEtaBins(50),
  fCentBins(10),
  nuacentral(),
  nuaforward(),
  doNUA(false),
  gap(0.0),
  mc(false)
{
}


Bool_t AliForwardFlowRun2Settings::ExtraEventCutFMD(TH2D forwarddNdedp, double cent, Bool_t mc){
  Bool_t useEvent = true;
  Int_t nBadBins = 0;
  Int_t phibins = forwarddNdedp.GetNbinsY();
  Double_t totalFMDpar = 0;

  for (Int_t etaBin = 1; etaBin <= forwarddNdedp.GetNbinsX(); etaBin++) {
  
    Double_t acceptance = 1.;
    Double_t eta = forwarddNdedp.GetXaxis()->GetBinCenter(etaBin);
    Double_t runAvg = 0;
    Double_t avgSqr = 0;
    Double_t max = 0;
    Int_t nInAvg = 0;


    for (Int_t phiBin = 0; phiBin <= phibins; phiBin++) {
      if (!mc){
        if ( fabs(eta) > 1.7) {
          if (phiBin == 0 && forwarddNdedp.GetBinContent(etaBin, 0) == 0) useEvent = false;
        }
      }
      Double_t weight = forwarddNdedp.GetBinContent(etaBin, phiBin);
      if (!weight){
        weight = 0;
      }
      totalFMDpar += weight;
      
      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) {
        max = weight;
      }
    } // End of phi loop

    // Outlier cut calculations
    double fSigmaCut = 4.0;
    if (nInAvg > 0) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = (nInAvg > 1 ? TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg) : 0);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent < 60) nBadBins++;
      else nBadBins = 0;
      // We still finish the loop, for fOutliers to make sense, 
      // but we do no keep the event for analysis 
      if (nBadBins > 3) useEvent = false;
     //if (nBadBins > 3) std::cout << "NUMBER OF BAD BINS > 3" << std::endl;
    }
  } // End of eta bin
  if (totalFMDpar < 10) useEvent = false;

  return useEvent;
}