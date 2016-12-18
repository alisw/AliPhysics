#include "THnSparse.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include "TObjArray.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TROOT.h"

#include <iostream>
#include <iomanip>

#include "AliPID.h"

#include "THnSparseDefinitions.h"
#include "histFitting/AliTPCPIDmathFit.h"

enum processMode { kPMpT = 0, kPMz = 1, kPMxi = 2, kPMdistance = 3, kPMjT = 4, kPMnum };
enum muonTreatment { kNoMuons = 0, kMuonFracEqualElFrac = 1, kMuonFracOverElFracTunedOnMCStandardTrackCuts = 2,
                     kMuonFracOverElFracTunedOnMCHybridTrackCuts = 3, kMuonFracOverElFracTunedOnMCHybridTrackCutsJets = 4,
                     kMuonFracOverElFracTunedOnMCStandardTrackCutsPPb = 5,
                     kNumHandlings = 6 };

const TString modeShortName[kPMnum] = { "Pt", "Z", "Xi", "R", "Jt" };
const TString modeLatexName[kPMnum] = { "p_{T}", "z", "#xi", "R", "j_{T}" };

const TString muonFractionHandlingShortName[kNumHandlings] =
  { "noMuons", "muonsEqualElectrons", "muonToElTunedOnMCStandardTrackCuts", "muonToElTunedOnMCHybridTrackCuts",
    "muonToElTunedOnMCHybridTrackCutsJets", "muonToElTunedOnMCStandardTrackCutsPPB" };

const Double_t epsilon = 1e-10;
const TString identifiedLabels[2] = { "Most Probable PID", "MC" };
Int_t isMC = 0;

const TString minimisationStrategy = "MIGRAD"; // "MINIMIZE"
Bool_t takeIntoAccountMuons = kTRUE;

// 0 = no muons, 1 = muonFrac=elFrac, 2(3) = muonFrac/elFrac tuned on MC for DefaultTrackCuts (hybridTrackCuts),
// 4 = muonFrac/elFrac tuned on MC for hybridTrackCuts for jet particles,
Int_t muonFractionHandling = 3; 

Int_t processingMode = kPMpT;

//TODO getErrorOf.... is COMPLETELY wrong now, since the parameter numbering has changed and the muons had come into play!!!!!!

// TODO CAREFUL: fitMethod == 1 adds errors of electrons to pions, but not to muons (should be added to electron error instead!)
const Bool_t muonContamination = kFALSE;//TODO CAREFUL: fitMethod == 1 takes into account the muon contamination in the error calculation!!!

const Bool_t normaliseResults = kTRUE; // Works only for fitMethod == 2

const Bool_t enableShift = kFALSE;
const Int_t dataAxis = kPidDeltaPrime;//kPidDelta; kPidDeltaPrime

const Int_t numSimultaneousFits = 4;

// Upper and lower axis bounds (y-axis) of (data - fit) / data QA histos
const Double_t fitQAaxisLowBound = -0.48;
const Double_t fitQAaxisUpBound = 0.48;

const Bool_t useDeltaPrime = (dataAxis == kPidDeltaPrime);

// Will be set later
Double_t muonFractionThresholdForFitting = -1.;
Double_t muonFractionThresholdBinForFitting = -1;
  
TF1 fMuonOverElFractionMC("fMuonOverElFractionMC", "[0]+[1]/TMath::Min(x, [4])+[2]*TMath::Min(x, [4])+[3]*TMath::Min(x, [4])*TMath::Min(x, [4])+[5]*TMath::Min(x, [4])*TMath::Min(x, [4])*TMath::Min(x, [4])+[6]*(x>[7])*TMath::Min(x-[7], [8]-[7])",
                          0.01, 50.);

TF1* fElectronFraction = 0x0;
Double_t lowFittingBoundElectronFraction = -1.; 

const Double_t gkLowFittingBoundElectronFraction = 3.0;
const Double_t gkElectronFractionThresholdForFitting = 10.0;

// Will be set later or might change in case of z and xi!
Double_t electronFractionThresholdForFitting = -1.;
Int_t electronFractionThresholdBinForFitting = -1;

TGraphErrors* gFractionElectronsData = 0x0;
Double_t lastPtForCallOfGetElectronFraction = -1;


// Most probable PID templates for low pT -> thresholds
// Values for z and xi will be calculated later according to jet pT
const Double_t mostProbPIDthresholdPt_pr = 0.45;
const Double_t mostProbPIDthresholdPt_ka = 0.3;

Double_t mostProbPIDthresholdZ_pr = -1.;
Double_t mostProbPIDthresholdZ_ka = -1.;

Double_t mostProbPIDthresholdXi_pr = 999.;
Double_t mostProbPIDthresholdXi_ka = 999.;

// There is no clear separation in any bin of R or jT, disable special templates
Double_t mostProbPIDthresholdR_pr = -1;
Double_t mostProbPIDthresholdR_ka = -1;

Double_t mostProbPIDthresholdJt_pr = -1;
Double_t mostProbPIDthresholdJt_ka = -1;

// Thresholds below which the regularisation is switched off for the corresponding species
// Values for z and xi will be calculated later according to jet pT
const Double_t regOffThresholdPt_pr = 0.45 - 1e-3;
const Double_t regOffThresholdPt_ka = 0.30 - 1e-3;

// Do not fit bins with less total entries than the following threshold
Double_t gYieldThresholdForFitting = 0;

// To be set at the beginning - pT binning
Int_t nPtBins = -1;
const Double_t* binsPt = 0x0;



//_____________________________________________________________________________________________________________________________
void ResetGlobals()
{
  // If 2 fits are run within same AliRoot, values must be reset. Otherwise, it might happen that things get mixed up!
  isMC = 0;
  takeIntoAccountMuons = kTRUE;
  muonFractionHandling = 3; 
  processingMode = kPMpT;
  muonFractionThresholdForFitting = -1.;
  muonFractionThresholdBinForFitting = -1;
  delete fElectronFraction;
  fElectronFraction = 0x0;
  lowFittingBoundElectronFraction = gkLowFittingBoundElectronFraction; 
  electronFractionThresholdForFitting = gkElectronFractionThresholdForFitting;
  electronFractionThresholdBinForFitting = -1;
  delete gFractionElectronsData;
  gFractionElectronsData = 0x0;
  lastPtForCallOfGetElectronFraction = -1;
  mostProbPIDthresholdZ_pr = -1.;
  mostProbPIDthresholdZ_ka = -1.;
  mostProbPIDthresholdXi_pr = 999.;
  mostProbPIDthresholdXi_ka = 999.;
  mostProbPIDthresholdR_pr = -1;
  mostProbPIDthresholdR_ka = -1;
  mostProbPIDthresholdJt_pr = -1;
  mostProbPIDthresholdJt_ka = -1;
  
  gYieldThresholdForFitting = 0.;
  
  nPtBins = -1;
  binsPt = 0x0;
}


//_____________________________________________________________________________________________________________________________
void RatioToRef(TH1* hToCompare, TH1* hRef)
{
  if (!hToCompare || !hRef)
    return;
  
  for (Int_t i = 1; i <= hToCompare->GetNbinsX(); i++) {
    const Double_t refValue = hRef->GetBinContent(i);
    const Double_t currValue = hToCompare->GetBinContent(i);
    
    if (refValue <= 0) {
      if (currValue > 0) {
        // Problem -> Not well defined, set error!
        hToCompare->SetBinContent(i, -999);
        hToCompare->SetBinError(i, 999);
      }
      else {
        hToCompare->SetBinContent(i, 0);
        hToCompare->SetBinError(i, 0);
      }
    }
    else {
      hToCompare->SetBinContent(i, currValue / refValue);
      hToCompare->SetBinError(i, hToCompare->GetBinError(i) / refValue);
    }
  }
}


//____________________________________________________________________________________________________________________
void PatchFractionWithTOF(Double_t& fractionTPC, Double_t& fractionErrorTPC, const Double_t yieldTPConlyTotal,
                          const Double_t yieldTOFspecies, const Double_t yieldTOFtotal)
{
  // Patch the fraction from TPC only with TOF in one single bin.
  // Only scale errors (assume TOF yields to be precise (and error of TPC yield is already "contained" in the fraction
  // by using the likelihood fit).

  if ((yieldTPConlyTotal + yieldTOFtotal) > 0.) {
    fractionTPC = (fractionTPC * yieldTPConlyTotal + yieldTOFspecies) / (yieldTPConlyTotal + yieldTOFtotal);
    fractionErrorTPC = fractionErrorTPC * yieldTPConlyTotal / (yieldTPConlyTotal + yieldTOFtotal);
  }
  else {
    if (fractionTPC > 1e-9)
      printf("Error TOF patching of fraction (denominator is zero): yieldTPConlyTotal %f, yieldTOFtotal %f, yieldTOFspecies %f, fractionTPC %f, fractionErrorTPC %f\n",
            yieldTPConlyTotal, yieldTOFtotal, yieldTOFspecies, fractionTPC, fractionErrorTPC);
  }
}


//____________________________________________________________________________________________________________________
inline Double_t PatchFractionWithTOFfast(const Double_t fractionTPC, const Double_t yieldTPConlyTotal, const Double_t yieldTOFspecies,
                                         const Double_t yieldTOFtotal)
{
  // Patch the fraction from TPC only with TOF in one single bin.

  const Double_t yieldTotalSummed = yieldTPConlyTotal + yieldTOFtotal;
  if (yieldTotalSummed > 0.) 
    return (fractionTPC * yieldTPConlyTotal + yieldTOFspecies) / yieldTotalSummed;
  
  return 0.; // Without yield, one cannot define a fraction at all
}


//____________________________________________________________________________________________________________________
inline Double_t UndoTOFpatchingForFraction(const Double_t patchedFractionTPC, const Double_t yieldTPConlyTotal, const Double_t yieldTOFspecies,
                                           const Double_t yieldTOFtotal)
{
  // Undo the TOF patching in a single bin.

  if (yieldTPConlyTotal > 0.)
    return (patchedFractionTPC * (yieldTPConlyTotal + yieldTOFtotal) -  yieldTOFspecies) / yieldTPConlyTotal;
    
  return patchedFractionTPC; // If no TPC yield, just return the TPC+TOF fraction as it is
}


//____________________________________________________________________________________________________________________
void PatchRatioWithTOF(Double_t& toPiRatioTPC, Double_t& toPiRatioErrorTPC, const Double_t yieldTPConlyPi,
                          const Double_t yieldTOFspecies, const Double_t yieldTOFpi)
{
  // Patch the to-pion-fraction from TPC only with TOF in one single bin.
  // This is simply using fraction_patched(x) / fraction_patched(pi) and deriving the error from the resulting formula.
  // Only scale errors (assume TOF yields to be precise (and error of TPC yield is already "contained" in the ratio
  // by using the likelihood fit).
  
  if (yieldTPConlyPi > 0.) {
    toPiRatioTPC = (toPiRatioTPC + yieldTOFspecies / yieldTPConlyPi) / (1. + yieldTOFpi / yieldTPConlyPi);
    toPiRatioErrorTPC = toPiRatioErrorTPC / (1. + yieldTOFpi / yieldTPConlyPi);
  }
  else {
    if (toPiRatioTPC > 1e-9)
      printf("Error TOF patching of to-pion-ratio (denominator is zero): yieldTPConlyPi %f, yieldTOFpi %f, yieldTOFspecies %f, toPiRatioTPC %f, toPiRatioErrorTPC %f\n",
            yieldTPConlyPi, yieldTOFpi, yieldTOFspecies, toPiRatioTPC, toPiRatioErrorTPC);
  }
}


//____________________________________________________________________________________________________________________
void ExtractTOFEfficiency(TH1* hYieldTOFKaons, TH1* hYieldTOFPions, TH1* hYieldTOFProtons,
                          TH1* hYieldKaons, TH1* hYieldPions, TH1* hYieldProtons,
                          TH1F** hTOFEfficiencyKaons, TH1F** hTOFEfficiencyPions, TH1F** hTOFEfficiencyProtons,
                          Bool_t fromMC)
{
  // Extract the efficiency of the TOF PID among the TPC tracks, i.e. Yield_TOF_k / (Yield_TOF_k + Yield_TPC_k)
  
  (*hTOFEfficiencyPions) = new TH1F(*((TH1F*)hYieldPions));
  (*hTOFEfficiencyPions)->SetName(Form("hTOFEfficiencyPions%s", fromMC ? "MC" : ""));
  (*hTOFEfficiencyPions)->GetYaxis()->SetTitle("TOF PID Efficiency");
  (*hTOFEfficiencyPions)->Add(hYieldTOFPions);
  // Binomial error, since efficiency like!
  (*hTOFEfficiencyPions)->Divide(hYieldTOFPions, (*hTOFEfficiencyPions), 1., 1., "B");
  (*hTOFEfficiencyPions)->SetDrawOption("histp0"); //Draw markers also for empty bins
  
  
  (*hTOFEfficiencyKaons) = new TH1F(*((TH1F*)hYieldKaons));
  (*hTOFEfficiencyKaons)->SetName(Form("hTOFEfficiencyKaons%s", fromMC ? "MC" : ""));
  (*hTOFEfficiencyKaons)->GetYaxis()->SetTitle("TOF PID Efficiency");
  (*hTOFEfficiencyKaons)->Add(hYieldTOFKaons);
  // Binomial error, since efficiency like!
  (*hTOFEfficiencyKaons)->Divide(hYieldTOFKaons, (*hTOFEfficiencyKaons), 1., 1., "B");
  (*hTOFEfficiencyKaons)->SetDrawOption("histp0"); //Draw markers also for empty bins
  
  (*hTOFEfficiencyProtons) = new TH1F(*((TH1F*)hYieldProtons));
  (*hTOFEfficiencyProtons)->SetName(Form("hTOFEfficiencyProtons%s", fromMC ? "MC" : ""));
  (*hTOFEfficiencyProtons)->GetYaxis()->SetTitle("TOF PID Efficiency");
  (*hTOFEfficiencyProtons)->Add(hYieldTOFProtons);
  // Binomial error, since efficiency like!
  (*hTOFEfficiencyProtons)->Divide(hYieldTOFProtons, (*hTOFEfficiencyProtons), 1., 1., "B");
  (*hTOFEfficiencyProtons)->SetDrawOption("histp0"); //Draw markers also for empty bins
  
}


//____________________________________________________________________________________________________________________
Double_t GetElectronFraction(const Double_t xCoordParameter, const Double_t *par, Int_t xBinParameter)
{
  // During the fit (both, simultaneous and non-simultaneous), the algorithm will always start off from
  // the low pT/z and go to higher pT/z (opposite direction for xi).
  // So, it is only necessary to do the fit once the first (last) fixed bin is reached (in case of xi).
  // Then the parameters for the electron fraction remain fixed until the next fit iteration.
  // Since only for the case of regularisation the electron fractions of all x bins are stored in mathFit,
  // the evaluation of this function is done here only in that case (only then the electron fraction will
  // be set to "-pT" (or "-z" or "-xi) and "-par" will be forwarded as xCoordParameter to this function.
  
  // NOTE 1: Electrons have index 3 per x bin
  // NOTE 2: In case of fitting vs. pT, xValue holds the LOG10 of pT!
 
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  // This function may only be called for enabled regularisation! Otherwise, statistical weights etc. are not available and it could crash
  // or give unpredictable results
  if (mathFit->GetNumXbinsRegularisation() <= 1) {
    printf("FATAL: GetElectronFraction called with NumXbinsRegularisation <= 1, i.e. regularisation is obviously off!\n");
    exit(-1);
  }
  if (xBinParameter < 0) {
    printf("FATAL: GetElectronFraction called with xBinParameter < 0, i.e. regularisation is obviously off!\n");
    exit(-1);
  }
  
  const Bool_t isPtMode = (processingMode == kPMpT);
  const Bool_t isXiMode = (processingMode == kPMxi);
  
  // lastPtForCallOfGetElectronFraction will be initialised with a value larger (smaller) than any pT/z (xi) during the fit.
  // So, if this function is called and the pT/z (xi) is smaller (larger) than lastPtForCallOfGetElectronFraction, the parameters
  // must have changed and the electron fit needs to be re-done (also see comment above)
  if ((isXiMode  && xCoordParameter > lastPtForCallOfGetElectronFraction) ||
      (!isXiMode && xCoordParameter < lastPtForCallOfGetElectronFraction)) {
    for (Int_t xBin = 0; xBin < mathFit->GetNumXbinsRegularisation(); xBin++) {
      
      // GetXvaluesForRegularisation holds (for pT mode) the Log10(pT)!
      const Double_t xCoord = isPtMode ? TMath::Power(10, mathFit->GetXvaluesForRegularisation()[xBin]) 
                                       : mathFit->GetXvaluesForRegularisation()[xBin];
      
      if ((isXiMode  && (xCoord > lowFittingBoundElectronFraction || xCoord < electronFractionThresholdForFitting)) ||
          (!isXiMode && (xCoord < lowFittingBoundElectronFraction || xCoord > electronFractionThresholdForFitting))) {
        gFractionElectronsData->SetPoint(xBin, -9999, 0);
        gFractionElectronsData->SetPointError(xBin, 0, 0);
        continue;
      }
      
      const Int_t parIndexWithFraction = 3 + xBin * mathFit->GetNumParametersPerXbin(); 
      
      // If TOF patching is not used, fractionTOFpatched is just the ordinary fraction
      Double_t fractionTOFpatched =  mathFit->GetApplyPatching()
                                        ? PatchFractionWithTOFfast(par[parIndexWithFraction],
                                                                   mathFit->GetXstatisticalWeight()[xBin],
                                                                   mathFit->GetParAdditional()[parIndexWithFraction],
                                                                   mathFit->GetXstatisticalWeight2()[xBin])
                                        : par[parIndexWithFraction];
      if (fractionTOFpatched > epsilon) { // Skip zero values (usually due to failed fits)
        gFractionElectronsData->SetPoint(xBin, xCoord, fractionTOFpatched);
        // Since the errors during the fitting are not reliable, use the following approximation on a statistical basis
        // (which indeed turns out to be rather good!)
        
        // Bin effective weight required for weighted data sets. In case of no weighting, the weight error is sqrt(weight),
        // i.e. effWeight is 1
        const Double_t effWeight = mathFit->GetXstatisticalWeightError()[xBin] * mathFit->GetXstatisticalWeightError()[xBin]
                                   / mathFit->GetXstatisticalWeight()[xBin];
        gFractionElectronsData->SetPointError(xBin, 0, effWeight * TMath::Sqrt(fractionTOFpatched 
                                                                               / mathFit->GetXstatisticalWeight()[xBin]));
      }
      else {
        gFractionElectronsData->SetPoint(xBin, -9999, 0);
        gFractionElectronsData->SetPointError(xBin, 0, 0);
      }
    }
    
    if (isXiMode)
      gFractionElectronsData->Fit(fElectronFraction, "Ex0NQ", "", electronFractionThresholdForFitting, lowFittingBoundElectronFraction);
    else
      gFractionElectronsData->Fit(fElectronFraction, "Ex0NQ", "", lowFittingBoundElectronFraction, electronFractionThresholdForFitting);
  }
  
  lastPtForCallOfGetElectronFraction = xCoordParameter;

  // If TOF patching is not used, the fraction from the fit is already without TOF patching.
  // If it is used, the fraction from the fit is WITH TOF patching and needs to be converted back to the TPC only fraction
  // to comply with the fit framework.
  Double_t outputFraction = fElectronFraction->Eval(xCoordParameter);

  if (mathFit->GetApplyPatching()) {
    const Int_t parIndexWithFraction = 3 +  xBinParameter * mathFit->GetNumParametersPerXbin(); 
    outputFraction = UndoTOFpatchingForFraction(outputFraction, mathFit->GetXstatisticalWeight()[xBinParameter],
                                                mathFit->GetParAdditional()[parIndexWithFraction],
                                                mathFit->GetXstatisticalWeight2()[xBinParameter]);
  }
  
  // Catch cases in which the fit function yields invalid fractions (i.e. < 0 or > 1)
  return TMath::Max(0.0, TMath::Min(1.0, outputFraction));
}


//____________________________________________________________________________________________________________________
Double_t GetElectronFractionError()
{
  // This function estimates the error of the electron fraction for the fixed values via using the parameter errors of
  // the electron fraction function. Note that the parameters (and errors) must be set before calling this function.
  
  // If there is no electron function set, then the value is most likely just fixed to zero. In this case, there is no valid
  // errors estimate possible, just return 0
  if (!fElectronFraction)
    return 0.;
  
  // Produce several values via setting the parameters to a random value, which is distributed with a gaussian with mean = parValue
  // and sigma = parError and then take the 2*RMS as the error
  const Int_t nGenValues = 1000;
  Double_t genValues[nGenValues];
  
  const Int_t nPars = fElectronFraction->GetNpar();
  Double_t par[nPars];
  
  TRandom3 rnd(0); // 0 means random seed
  
  const Double_t x = (processingMode == kPMxi)
                      ? electronFractionThresholdForFitting - 1.  // Some value below the threshold to obtain a fixed value
                      : electronFractionThresholdForFitting + 1.; // Some value above the threshold to obtain a fixed value
  for (Int_t i = 0 ; i < nGenValues; i++) {
    for (Int_t iPar = 0; iPar < nPars; iPar++)
      par[iPar] = rnd.Gaus(fElectronFraction->GetParameter(iPar), fElectronFraction->GetParError(iPar));
    
    genValues[i] = fElectronFraction->EvalPar(&x, &par[0]);
  }
  
  // NOTE: RMS is not really the root mean square, is it rather the sigma deviation, which is what is wanted here
  return 2. * TMath::RMS(nGenValues, &genValues[0]);
}


//____________________________________________________________________________________________________________________
Double_t GetMuonFractionFromElectronFractionAndPt(Double_t pT, Double_t elFrac)
{
  // NOTE: Since the electron fraction will be tuned on all (TOF + TPC), but the patching will be undone for the fit,
  // also the muons will get the correct (not TOF patched) fraction and will after TOF patching have the desired value.
  if (muonFractionHandling == kMuonFracOverElFracTunedOnMCStandardTrackCuts) {
//    return elFrac / (1. + 7.06909e+01 * TMath::Exp(-2.95078e+00 * TMath::Power(pT, 5.05016e-01)));
    return elFrac / (1. + 2.01840e+10 * TMath::Exp(-2.50480e+01 * TMath::Power(pT, 5.89044e-02)));
  }
  else if (muonFractionHandling == kMuonFracOverElFracTunedOnMCHybridTrackCuts) {
    fMuonOverElFractionMC.SetParameters(-6.87241e-01, 4.19528e-02, 4.52095e+00, -6.20026e+00, 5.16629e-01, 2.88604e+00, 3.68058e-02,
                                        2.21086e+00, 5.75003e+00);
    return elFrac * fMuonOverElFractionMC.Eval(pT);
  }
  else if (muonFractionHandling == kMuonFracOverElFracTunedOnMCHybridTrackCutsJets) {
    fMuonOverElFractionMC.SetParameters(-7.64548e-01, 2.47929e-02, 4.49057e+00, -2.06320e-01, 4.23339e-02, 1.19697e+02, 1.28832e-01,
                                        -1.71895e-01, 6.00000e+00);
    return elFrac * fMuonOverElFractionMC.Eval(pT);
  }
  else if (muonFractionHandling == kMuonFracOverElFracTunedOnMCStandardTrackCutsPPb) {
    // WITH PID cluster cut!
    fMuonOverElFractionMC.SetParameters(-6.62149e-01, 4.89591e-02, 4.58356e+00, -6.04319e+00, 6.25368e-01, 3.27191e+00, 1.69933e-01,
                                        1.00004e+00, 2.61438e+00);
    return elFrac * fMuonOverElFractionMC.Eval(pT);
  }
  else if (muonFractionHandling == kMuonFracEqualElFrac) {
    return elFrac;
  }
  
  return 0.;
}


//____________________________________________________________________________________________________________________
Double_t GetCorrelatedError(const Double_t x, const Double_t y, const Double_t cov00, const Double_t cov11, const Double_t cov01) 
{
  // Calculate the correlated error df of f:
  //                (cov00 cov01) (x)
  //df^2 = (x, y) * (cov01 cov11) (y) = x^2 * cov00 + y^2 * cov11 + 2 * x * y * cov01
  //
  // with  f = f(p1, p2) = p1 / p2
  // and (x, y) = (\partial f / \partial p1, \partial f / \partial p2)
  //            = (f / p1, -f / p2)

  const Double_t df2 = x * x * cov00 + y * y * cov11 + 2. * x * y * cov01;
  
  if (df2 < epsilon)
    return 0.;
  
  return TMath::Sqrt(df2);
}


//____________________________________________________________________________________________________________________
void GetRatioWithCorrelatedError(const Double_t fractionA, const Double_t fractionB,
                                 const Double_t fractionErrorA, const Double_t fractionErrorB,
                                 const Double_t covMatrixElementAB, Double_t& ratio, Double_t& ratioError)
{
  // Given fractions A and B with corresponding errors and the off-diagonal covariance matrix element of
  // these fractions, calculate the ratio A/B and the error taking into account the correlation.
  // The results are stored in ratio and ratioError.
  
  if (fractionB < epsilon) {
    ratio = -999.;
    ratioError = 999.;
    
    return;
  }
  
  /*
  if (fractionA < epsilon) {
    ratio = 0.;
    ratioError = 999.;
    
    return;
  }*/
  
  ratio = fractionA / fractionB;
  
  const Double_t x = 1. / fractionB;// = ratio / fractionA; -> New definition properly takes into account the case fractionA=0
  const Double_t y = -ratio / fractionB;
  
  // covMatrixElement(i, i) = error(i)^2
  ratioError = GetCorrelatedError(x, y, fractionErrorA * fractionErrorA, fractionErrorB * fractionErrorB, covMatrixElementAB); 
  
  //printf("frationA %e\nfractionB %e\nfractionErrorA %e\nfractionErrorB %e\ncovMatrixElementAB %e\nratio %e\nx %e\ny %e\nratioError %e\n\n",
  //       fractionA, fractionB, fractionErrorA, fractionErrorB, covMatrixElementAB, ratio, x, y, ratioError);
}


//____________________________________________________________________________________________________________________
TH1D* GetFitQAhisto(TH1D* hData, TF1* fFit, TString histName)
{
  // Calculate (data - fit) / data assuming fit to have no error
  if (!hData || !fFit)
    return 0x0;
 
  Double_t binContent = 0.;
  Double_t binError = 0.;
  Double_t eval = 0.;
  
  Double_t xx[3] = { 0., 0., 0. };
  Double_t *params = 0;
  fFit->InitArgs(xx, params);
 
  TH1D* hDeltaFitQA = (TH1D*)hData->Clone(histName.Data());
  
  for (Int_t binX = 1; binX <= hDeltaFitQA->GetNbinsX(); binX++) {
    xx[0] = hDeltaFitQA->GetXaxis()->GetBinCenter(binX);
    binContent = hDeltaFitQA->GetBinContent(binX);
    binError = hDeltaFitQA->GetBinError(binX);
    
    eval = fFit->EvalPar(xx);
    if (binContent) {
      hDeltaFitQA->SetBinContent(binX, (binContent - eval) / binContent);
      hDeltaFitQA->SetBinError(binX, TMath::Abs(eval / binContent * binError / binContent));
    }
    else {
      // As is done in root divide histo
      hDeltaFitQA->SetBinContent(binX, 0);
      hDeltaFitQA->SetBinError(binX, 0);
    }
  }
  
  hDeltaFitQA->GetYaxis()->SetTitle("(Data - Fit) / Data");
  hDeltaFitQA->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
  
  return hDeltaFitQA;
}
  

//____________________________________________________________________________________________________________________
TH1F* GetHistWithProperXbinning(THnSparse* hPIDdata, Int_t axisForMode, Int_t mode, const TString histName, const TString histTitle)
{
  // Create histogram with proper binning for the x-axis depending on the processing mode "mode"
  TH1F* h = 0x0;
  if (mode == kPMpT)
    h = new TH1F(histName.Data(), histTitle.Data(), nPtBins, binsPt);
  else {
    const TArrayD* histBins = hPIDdata->GetAxis(axisForMode)->GetXbins();
    if (histBins->fN == 0)
      h = new TH1F(histName.Data(), histTitle.Data(), hPIDdata->GetAxis(axisForMode)->GetNbins(),
                   hPIDdata->GetAxis(axisForMode)->GetXmin(), hPIDdata->GetAxis(axisForMode)->GetXmax());
    else
      h = new TH1F(histName.Data(), histTitle.Data(), hPIDdata->GetAxis(axisForMode)->GetNbins(), histBins->fArray);
  }
  
  return h;
}


//____________________________________________________________________________________________________________________
void SetStyleSpecialTemplate(TH1* h)
{
  if (!h)
    return;
  
  h->SetMarkerStyle(20);
}


//____________________________________________________________________________________________________________________
Double_t GetFractionErrorSpecialTemplate(Double_t totalYield, Double_t totalYieldError, TH1D* hTemplateUnnormalised)
{
  Double_t fracError = 0.;
  
  Double_t e1 = 0.;
  const Double_t b1 = hTemplateUnnormalised->IntegralAndError(hTemplateUnnormalised->GetXaxis()->GetFirst(),
                                                              hTemplateUnnormalised->GetXaxis()->GetLast(),
                                                              e1);
  const Double_t b2 = totalYield;
  if (b2 && b1 != b2) {
    const Double_t w = b1 / b2;
    const Double_t e2 = totalYieldError; 
    // In case of e1^2 = N1 (similar for e2), this is just the formula for the binomial error of the fraction f: delta_f = sqrt(f(1-f)/N_tot)
    fracError = TMath::Sqrt(TMath::Abs(((1. - 2. * w) * e1*e1 + w*w * e2*e2 ) / (b2*b2) ));
  }
  
  return fracError;
}


//____________________________________________________________________________________________________________________
Double_t GetFractionSpecialTemplate(Double_t totalYield, TH1D* hTemplateUnnormalised)
{
  // Calculate error for special templates. It is just the BINOMIAL error for yield_spec / yield_tot (true subsets in this case)
  // calcualted exactly as in root (TH1)
  if (!hTemplateUnnormalised)
    return 0;
  
  if (totalYield > 0)
    return (hTemplateUnnormalised->Integral() / totalYield);
  
  return 0.;
}


//____________________________________________________________________________________________________________________
void SetReasonableAxisRange(TAxis* axis, Int_t mode, Double_t pLow = -1, Double_t pHigh = -1)
{
  if (mode == kPMpT)
    axis->SetRangeUser(TMath::Max(0.15, pLow - 0.1), TMath::Min(50., pHigh + 0.1));
  else if (mode == kPMz)
    axis->SetRange(0, -1);
  else if (mode == kPMxi)
    axis->SetRange(0, -1);
  else if (mode == kPMdistance)
    axis->SetRange(0, -1);
  else if (mode == kPMjT)
    axis->SetRange(0, -1);
}

//____________________________________________________________________________________________________________________
void SetReasonableXaxisRange(TH1* h, Int_t& binLow, Int_t& binHigh)
{
  binLow = TMath::Max(1, h->FindFirstBinAbove(0));
  binHigh  = TMath::Min(h->GetNbinsX(), h->FindLastBinAbove(0));
  
  binLow = TMath::Min(binLow, h->GetXaxis()->FindFixBin(0.55));
  binHigh = TMath::Max(binHigh, h->GetXaxis()->FindFixBin(1.65));
  
  h->GetXaxis()->SetRange(binLow, binHigh);
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNoExponent(kTRUE);
}


//____________________________________________________________________________________________________________________
Int_t FindMomentumBin(const Double_t* pTbins, const Double_t value, const Int_t numPtBins = nPtBins)
{
  for (Int_t bin = 0; bin < numPtBins; bin++) {
    if (value >= pTbins[bin] && value < pTbins[bin + 1])
      return bin + 1; // Must be + 1 since bin of histogram starts with 1 and not with zero as the index of pTbins!
  }
  
  return -1;
}


//____________________________________________________________________________________________________________________
Double_t normaliseHist(TH1* h, Bool_t& hasNonVanishingIntegral, Double_t scaleFactor)
{
  // Scales the histogram with the scale factor. If the scale factor is < 0,
  // the histogram is normalised to it's integral.
  // In both cases, the normalisation factor is returned.
  // If scaleFactor < 0, hasNonVanishingIntegral is set to kTRUE in case of h->Integral() > 0.

  
  Double_t normFactor = 1.;
  hasNonVanishingIntegral = kFALSE;
  
  if (scaleFactor < 0) {
    Double_t integralTemp = h->Integral();
    if (integralTemp > 0) {
      hasNonVanishingIntegral = kTRUE;
      normFactor = 1.0 / integralTemp;
      h->Scale(normFactor);
    }
  }
  else {
    normFactor = scaleFactor;
    h->Scale(normFactor);
  }
  
  h->GetXaxis()->SetTitleOffset(1.0);
  
  return normFactor;
}


//____________________________________________________________________________________________________________________
Double_t normaliseHist(TH1* h, Double_t scaleFactor)
{
  // Function overload, which does not use "hasNonVanishingIntegral".
  
  Bool_t dummy = kFALSE;
  return normaliseHist(h, dummy, scaleFactor);
}


//____________________________________________________________________________________________________________________
void normaliseYieldHist(TH1* h, Double_t numEvents, Double_t deta, Bool_t normaliseToJets, Double_t numJets)
{
  // If not looking at jets (then it makes only sense to look at pT):
  // Yield histos are already normalised to dpT. Now normalise to 1/NeV 1/(2pi pT) 1/deta in addition.
  //
  // If looking at jets:
  // Yield histos are already normalised to bin width. Now just scale with 1/Njets.
  //
  // For data, the normalisation was done via inverseBinWidth in SetFractionsAndYields. For MC, this was done explicitely.
  
  if (normaliseToJets) {
    if (numJets > 0) {
      h->Scale(1. / numJets);
    }
    else
      printf("WARNING: Number of rec. jets not available. Histo \"%s\" will not be normalised to this number!\n", h->GetName());
  }
  else {
    if (numEvents <= 0) // Do not normalise
      numEvents = 1; 
    
    for (Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
      Double_t normFactor = 1. / (numEvents * 2 * TMath::Pi() * h->GetXaxis()->GetBinCenter(bin) * deta);
      h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
      h->SetBinError(bin, h->GetBinError(bin) * normFactor);
    }
  }
}


//____________________________________________________________________________________________________________________
void normaliseGenYieldMCtruthHist(TH1* h, Double_t numEvents, Double_t deta, Bool_t normaliseToJets, Double_t numJets)
{
  // If not looking at jets (then it makes only sense to look at pT):
  // Yield histos are NOT normalised to dpT yet. Now normalise to 1/NeV 1/(2pi pT) 1/deta 1/dpT.
  //
  // If looking at jets:
  // Yield histos are NOT normalised to bin width yet. Now normalise to 1/Njets 1/binWidth.
  
  if (normaliseToJets) {
    if (numJets <= 0) { // Do not normalise
      numJets = 1.;
      printf("WARNING: Number of gen. jets not available. Histo \"%s\" will not be normalised to this number!\n", h->GetName());
    }
    
    for (Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
      Double_t normFactor = 1. / (numJets * h->GetXaxis()->GetBinWidth(bin));
      h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
      h->SetBinError(bin, h->GetBinError(bin) * normFactor);
    }
  }
  else {
    if (numEvents <= 0) // Do not normalise
      numEvents = 1.; 
    
    for (Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
      Double_t normFactor = 1. / (numEvents * 2 * TMath::Pi() * h->GetXaxis()->GetBinCenter(bin) * h->GetXaxis()->GetBinWidth(bin) * deta);
      h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
      h->SetBinError(bin, h->GetBinError(bin) * normFactor);
    }
  }
}


//____________________________________________________________________________________________________________________
void setUpFitFunction(TF1* fitFunc, Int_t nBins, Bool_t noShift = kFALSE)
{
  fitFunc->SetLineColor(kGray + 1);
  fitFunc->SetLineWidth(2);
  fitFunc->SetLineStyle(1);
  fitFunc->SetNpx(nBins * 100);
  fitFunc->SetParName(0, "Pion fraction");
  fitFunc->SetParName(1, "Kaon fraction");
  fitFunc->SetParName(2, "Proton fraction");
  fitFunc->SetParName(3, "Electron fraction");
  fitFunc->SetParName(4, "Muon fraction");
  fitFunc->SetParName(5, "Total yield");
  if (noShift == kFALSE) {
    fitFunc->SetParName(6, "Shift of pion peak");
    fitFunc->SetParName(7, "Shift of kaon peak");
    fitFunc->SetParName(8, "Shift of proton peak");
    fitFunc->SetParName(9, "Shift of electron peak");
    fitFunc->SetParName(10, "Shift of muon peak");
  }
}


//____________________________________________________________________________________________________________________
inline Int_t findBinWithinRange(const TAxis* axis, Double_t value)
{
  Int_t bin = axis->FindFixBin(value);
  if (bin <= 0)
    bin = 1;
  if (bin > axis->GetNbins())
    bin = axis->GetNbins();
  
  return bin;
}


//____________________________________________________________________________________________________________________
Double_t linearInterpolation(const TH1* h, Double_t x, Double_t scaleFactor, Double_t shift, Double_t* error)
{  
  // Do linear interpolation between 2 bins to handle non-integer values of the shift parameters.
  // The shift also introduces some uncertainty, which is rather hard to estimate. Therefore, just take the maximum error of the involved bins.
  const Double_t xShifted = x - shift;

  // Just take value of bin, if beyond center of first/last bin
  if (xShifted <= h->GetBinCenter(1)) {
    if (error)
      *error = h->GetBinError(1) * scaleFactor;
    return h->GetBinContent(1) * scaleFactor;
  }
  else if(xShifted >= h->GetBinCenter(h->GetNbinsX())) {
    if (error)
      *error = h->GetBinError(h->GetNbinsX()) * scaleFactor;
    return h->GetBinContent(h->GetNbinsX()) * scaleFactor;
  }
  else {
    const Int_t xbin = h->FindFixBin(xShifted);
    Double_t x0, x1, y0, y1;
    
    if(xShifted <= h->GetBinCenter(xbin)) {
      y0 = h->GetBinContent(xbin - 1);
      x0 = h->GetBinCenter(xbin - 1);
      y1 = h->GetBinContent(xbin);
      x1 = h->GetBinCenter(xbin);
      
      if (error)
        *error = TMath::Max(h->GetBinError(xbin - 1), h->GetBinError(xbin)) * scaleFactor;
    } 
    else {
      y0 = h->GetBinContent(xbin);
      x0 = h->GetBinCenter(xbin);
      y1 = h->GetBinContent(xbin + 1);
      x1 = h->GetBinCenter(xbin + 1);
      
      if (error)
        *error = TMath::Max(h->GetBinError(xbin), h->GetBinError(xbin + 1)) * scaleFactor;
    }
    
    return scaleFactor * (y0 + (xShifted - x0) * ((y1 - y0) / (x1 - x0)));
  }
  
  return 0;
      
  /*Old version available for code bevor 03.05.2013*/
}


//____________________________________________________________________________________________________________________
void shiftHist(TH1D* h, Double_t shift, Bool_t useRegularisation = kFALSE)
{
  // Shift not available for regularisation. Just for convenience (can use the same code and only set one flag)
  // call this functions and then do nothing.
  // Actually, the shift is not availabe for simultaneous fitting also, but the parameter is just set to 0 there
  if (!h || useRegularisation)
    return;
  
  TString name = h->GetName();
  TH1D* hTemp = (TH1D*)h->Clone(Form("%s_clone", name.Data()));
  h->Reset();
  
  Double_t error = 0;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    // linearInterpolation with scaleFactor = 1.0, since histo is assumed to be properly scaled
    h->SetBinContent(i,  linearInterpolation(hTemp, h->GetXaxis()->GetBinCenter(i), 1.0, shift, &error));
    h->SetBinError(i, error);
  }
  
  delete hTemp;
}


//____________________________________________________________________________________________________________________
Double_t multiGaussFitForSimultaneousFitting(const Double_t *xx, const Double_t *par, const Int_t offset)
{
  // Offset for reference histos for delta_Species
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  // parXbinIndex (fixed) will be used my mathfit to hold the pT bin index (needed for regularisation)
  const Int_t xBinIndex = mathFit->GetXbinIndex();
  const Int_t numParsPerXbin = mathFit->GetNumParametersPerXbin();
  
  const Int_t numRefHistosPerFit = numSimultaneousFits + (takeIntoAccountMuons ? 1 : 0);
  const Int_t numRefHistosPerXbin = numRefHistosPerFit * numSimultaneousFits;
  
  const Int_t refHistOffset = offset * numRefHistosPerFit + xBinIndex * numRefHistosPerXbin;
  
  const TH1* hRefPi = mathFit->GetRefHisto(0 + refHistOffset);
  const TH1* hRefKa = mathFit->GetRefHisto(1 + refHistOffset);
  const TH1* hRefPr = mathFit->GetRefHisto(2 + refHistOffset);
  const TH1* hRefEl = mathFit->GetRefHisto(3 + refHistOffset);
  const TH1* hRefMu = takeIntoAccountMuons ? mathFit->GetRefHisto(4 + refHistOffset) : 0x0;
  
  if (!hRefEl || !hRefKa || !hRefPi || !hRefPr)
    return 0;
  
  if (takeIntoAccountMuons && !hRefMu)
    return 0;
  
  const Int_t parOffset = xBinIndex * numParsPerXbin;
  const Int_t parPi = 0 + parOffset;
  const Int_t parKa = 1 + parOffset;
  const Int_t parPr = 2 + parOffset;
  const Int_t parEl = 3 + parOffset;
  const Int_t parMu = 4 + parOffset;
  const Int_t parAll = 5 + parOffset;
  
  const Double_t scaleFactorPi = par[parAll] * (par[parPi] + (muonContamination ? par[parEl] : 0));
  const Double_t scaleFactorKa = par[parAll] * par[parKa];
  const Double_t scaleFactorPr = par[parAll] * par[parPr];
  //NOTE par[parEl] < 0 only in case of enabled regularisation. Otherwise, GetElectronFraction will NOT work!
  const Double_t parElFraction = (par[parEl] < 0) ? GetElectronFraction(-par[parEl], par, xBinIndex) : par[parEl];
  const Double_t scaleFactorEl = par[parAll] * parElFraction;
  // Fix muon fraction to electron fraction (or some modified electron fraction) if desired, i.e. corresponding par < 0
  const Double_t scaleFactorMu = (par[parMu] < 0)
                                    ? (par[parAll] * GetMuonFractionFromElectronFractionAndPt(-par[parMu], parElFraction)) 
                                    : (par[parAll] * par[parMu]);
  
  // Since one is looking at the same deltaSpecies for all considered species, the reference histograms have the same axes
  // => Only need to search for the bin once
  const Int_t binWithinRange = findBinWithinRange(hRefPi->GetXaxis(), xx[0]);
  const Double_t countPi = scaleFactorPi * hRefPi->GetBinContent(binWithinRange);
  const Double_t countKa = scaleFactorKa * hRefKa->GetBinContent(binWithinRange);
  const Double_t countPr = scaleFactorPr * hRefPr->GetBinContent(binWithinRange);
  const Double_t countEl = scaleFactorEl * hRefEl->GetBinContent(binWithinRange);
  const Double_t countMu = takeIntoAccountMuons ? scaleFactorMu * hRefMu->GetBinContent(binWithinRange) : 0;
  
  const Double_t res = countPi + countKa + countPr + countEl + countMu;
  
  
  return res;
}


//____________________________________________________________________________________________________________________
inline Double_t multiGaussFitDeltaPi(const Double_t *xx, const Double_t *par)
{
  return multiGaussFitForSimultaneousFitting(xx, par, 0);
}

//____________________________________________________________________________________________________________________
inline Double_t multiGaussFitDeltaKa(const Double_t *xx, const Double_t *par)
{
  return multiGaussFitForSimultaneousFitting(xx, par, 1);
}

//____________________________________________________________________________________________________________________
inline Double_t multiGaussFitDeltaPr(const Double_t *xx, const Double_t *par)
{
  return multiGaussFitForSimultaneousFitting(xx, par, 2);
}

//____________________________________________________________________________________________________________________
inline Double_t multiGaussFitDeltaEl(const Double_t *xx, const Double_t *par)
{
  return multiGaussFitForSimultaneousFitting(xx, par, 3);
}


//____________________________________________________________________________________________________________________
Double_t multiGaussFit(const Double_t *xx, const Double_t *par)
{
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  const TH1* hRefPi = mathFit->GetRefHisto(0);
  const TH1* hRefKa = mathFit->GetRefHisto(1);
  const TH1* hRefPr = mathFit->GetRefHisto(2);
  const TH1* hRefEl = mathFit->GetRefHisto(3);
  const TH1* hRefMu = takeIntoAccountMuons ? mathFit->GetRefHisto(4) : 0x0;
  
  if (!hRefEl || !hRefKa || !hRefPi || !hRefPr)
    return 0;
  
  if (takeIntoAccountMuons && !hRefMu)
    return 0;  
  
  // Do linear interpolation between 2 bins to handle non-integer values of the shift parameters
  const Double_t scaleFactorPi = par[5] * (par[0] + (muonContamination ? par[3] : 0));
  const Double_t scaleFactorKa = par[5] * par[1];
  const Double_t scaleFactorPr = par[5] * par[2];
  //NOTE par[3] < 0 will never happen because it is set only in case of regularisation which uses a different function call!
  const Double_t parElFraction = (par[3] < 0) ? GetElectronFraction(-par[3], par, -1) : par[3];
  const Double_t scaleFactorEl = par[5] * parElFraction;
  // Fix muon fraction to electron fraction (or some modified electron fraction) if desired, i.e. corresponding par < 0
  const Double_t scaleFactorMu = (par[4] < 0)
                                    ? (par[5] * GetMuonFractionFromElectronFractionAndPt(-par[4], parElFraction)) 
                                    : (par[5] * par[4]);

  const Double_t countPi = linearInterpolation(hRefPi, xx[0], scaleFactorPi, par[6], 0x0);
  const Double_t countKa = linearInterpolation(hRefKa, xx[0], scaleFactorKa, par[7], 0x0);
  const Double_t countPr = linearInterpolation(hRefPr, xx[0], scaleFactorPr, par[8], 0x0);
  const Double_t countEl = linearInterpolation(hRefEl, xx[0], scaleFactorEl, par[9], 0x0);
  const Double_t countMu = takeIntoAccountMuons ? linearInterpolation(hRefMu, xx[0], scaleFactorMu, par[10], 0x0) : 0;
  
  const Double_t res = countPi + countKa + countPr + countEl + countMu;
  
  /*
  const Double_t countPi = linearInterpolation(hRefPi, xx[0], par[6], 0x0);
  const Double_t countKa = linearInterpolation(hRefKa, xx[0], par[7], 0x0);
  const Double_t countPr = linearInterpolation(hRefPr, xx[0], par[8], 0x0);
  const Double_t countEl = linearInterpolation(hRefEl, xx[0], par[9], 0x0);
  const Double_t countMu = takeIntoAccountMuons ? linearInterpolation(hRefMu, xx[0], par[10], 0x0) : 0; 
  
  const Double_t res = par[5] * ((par[0] + (muonContamination ? par[3] : 0)) * countPi + par[1] * countKa 
                                 + par[2] * countPr + par[3] * countEl + par[4] * countMu);

  */
  
  return res;
}


//____________________________________________________________________________________________________________________
Double_t errorOfFitHistosForSimultaneousFitting(const Double_t *xx, const Double_t *par, const Int_t offset)
{
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  Double_t summedError = 0;
  
  // parXbinIndex (fixed) will be used my mathfit to hold the pT bin index (needed for regularisation)
  const Int_t xBinIndex = mathFit->GetXbinIndex();
  const Int_t numParsPerXbin = mathFit->GetNumParametersPerXbin();
  
  const Int_t numRefHistosPerFit = numSimultaneousFits + (takeIntoAccountMuons ? 1 : 0);
  const Int_t numRefHistosPerXbin = numRefHistosPerFit * numSimultaneousFits;
  
  const Int_t refHistOffset = offset * numRefHistosPerFit + xBinIndex * numRefHistosPerXbin;
  
  const TH1* hRefPi = mathFit->GetRefHisto(0 + refHistOffset);
  const TH1* hRefKa = mathFit->GetRefHisto(1 + refHistOffset);
  const TH1* hRefPr = mathFit->GetRefHisto(2 + refHistOffset);
  const TH1* hRefEl = mathFit->GetRefHisto(3 + refHistOffset);
  const TH1* hRefMu = takeIntoAccountMuons ? 
                        mathFit->GetRefHisto(4 + refHistOffset) 
                        : 0x0;
  
  if (!hRefEl || !hRefKa || !hRefPi || !hRefPr)
    return 0;
  
  if (takeIntoAccountMuons && !hRefMu)
    return 0;
  
  const Int_t parOffset = xBinIndex * numParsPerXbin;
  const Int_t parPi = 0 + parOffset;
  const Int_t parKa = 1 + parOffset;
  const Int_t parPr = 2 + parOffset;
  const Int_t parEl = 3 + parOffset;
  const Int_t parMu = 4 + parOffset;
  const Int_t parAll = 5 + parOffset;
  
  const Double_t scaleFactorPi = par[parAll] * (par[parPi] + (muonContamination ? par[parEl] : 0));
  const Double_t scaleFactorKa = par[parAll] * par[parKa];
  const Double_t scaleFactorPr = par[parAll] * par[parPr];
  //NOTE par[parEl] < 0 only in case of enabled regularisation. Otherwise, GetElectronFraction will NOT work!
  const Double_t scaleFactorEl = par[parAll] * ((par[parEl] < 0) ? GetElectronFraction(-par[parEl], par, xBinIndex) : par[parEl]);
  // Fix muon fraction to electron fraction (or some modified electron fraction) if desired, i.e. corresponding par < 0
  const Double_t scaleFactorMu = (par[parMu] < 0) ? 
                                    (scaleFactorEl * GetMuonFractionFromElectronFractionAndPt(-par[parMu], par[parEl]) / par[parEl]) 
                                    : (par[parAll] * par[parMu]);
  
  Double_t errorPi = 0, errorKa = 0, errorPr = 0, errorEl = 0, errorMu = 0;

  // Do linear interpolation between 2 bins to handle non-integer values of the shift parameters
  // Shift not implemented for simultaneous fit -> Just set all corresponding parameters to zero
  linearInterpolation(hRefPi, xx[0], scaleFactorPi, 0, &errorPi);
  linearInterpolation(hRefKa, xx[0], scaleFactorKa, 0, &errorKa);
  linearInterpolation(hRefPr, xx[0], scaleFactorPr, 0, &errorPr);
  linearInterpolation(hRefEl, xx[0], scaleFactorEl, 0, &errorEl);
  if (takeIntoAccountMuons)
    linearInterpolation(hRefMu, xx[0], scaleFactorMu, 0, &errorMu);
  
  summedError += TMath::Power(errorPi, 2);
  summedError += TMath::Power(errorKa, 2);
  summedError += TMath::Power(errorPr, 2);
  summedError += TMath::Power(errorEl, 2);
  if (takeIntoAccountMuons)
    summedError += TMath::Power(errorMu, 2);
  
  return summedError;
}


//____________________________________________________________________________________________________________________
inline Double_t errorOfFitHistosDeltaPi(const Double_t *xx, const Double_t *par)
{
  return errorOfFitHistosForSimultaneousFitting(xx, par, 0);
}


//____________________________________________________________________________________________________________________
inline Double_t errorOfFitHistosDeltaKa(const Double_t *xx, const Double_t *par)
{
  return errorOfFitHistosForSimultaneousFitting(xx, par, 1);
}


//____________________________________________________________________________________________________________________
inline Double_t errorOfFitHistosDeltaPr(const Double_t *xx, const Double_t *par)
{
  return errorOfFitHistosForSimultaneousFitting(xx, par, 2);
}


//____________________________________________________________________________________________________________________
inline Double_t errorOfFitHistosDeltaEl(const Double_t *xx, const Double_t *par)
{
  return errorOfFitHistosForSimultaneousFitting(xx, par, 3);
}


//____________________________________________________________________________________________________________________
Double_t errorOfFitHistos(const Double_t *xx, const Double_t *par)
{
  //TODO Error of shift is still not taken into account
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  Double_t summedError = 0;
  
  const TH1* hRefPi = mathFit->GetRefHisto(0);
  const TH1* hRefKa = mathFit->GetRefHisto(1);
  const TH1* hRefPr = mathFit->GetRefHisto(2);
  const TH1* hRefEl = mathFit->GetRefHisto(3);
  const TH1* hRefMu = takeIntoAccountMuons ? mathFit->GetRefHisto(4) : 0x0;

  if (!hRefEl || !hRefKa || !hRefPi || !hRefPr)
    return 0;

  if (takeIntoAccountMuons && !hRefMu)
    return 0;
  
  // Do linear interpolation between 2 bins to handle non-integer values of the shift parameters
  const Double_t scaleFactorPi = par[5] * (par[0] + (muonContamination ? par[3] : 0));
  const Double_t scaleFactorKa = par[5] * par[1];
  const Double_t scaleFactorPr = par[5] * par[2];
  //NOTE par[3] < 0 will never happen because it is set only in case of regularisation which uses a different function call!
  const Double_t scaleFactorEl = par[5] * ((par[3] < 0) ? GetElectronFraction(-par[3], par, -1) : par[3]);
  // Fix muon fraction to electron fraction (or some modified electron fraction) if desired, i.e. corresponding par < 0
  const Double_t scaleFactorMu = (par[4] < 0) ? (scaleFactorEl * GetMuonFractionFromElectronFractionAndPt(-par[4], par[3]) / par[3]) 
                                              : (par[5] * par[4]);
  
  Double_t errorPi = 0, errorKa = 0, errorPr = 0, errorEl = 0, errorMu = 0;

  linearInterpolation(hRefPi, xx[0], scaleFactorPi, par[6], &errorPi);
  linearInterpolation(hRefKa, xx[0], scaleFactorKa, par[7], &errorKa);
  linearInterpolation(hRefPr, xx[0], scaleFactorPr, par[8], &errorPr);
  linearInterpolation(hRefEl, xx[0], scaleFactorEl, par[9], &errorEl);
  if (takeIntoAccountMuons)
    linearInterpolation(hRefMu, xx[0], scaleFactorMu, par[10], &errorMu); // Assume same fraction as electron, i.e. same scale factor
  
  summedError += TMath::Power(errorPi, 2);
  summedError += TMath::Power(errorKa, 2);
  summedError += TMath::Power(errorPr, 2);
  summedError += TMath::Power(errorEl, 2);
  if (takeIntoAccountMuons)
    summedError += TMath::Power(errorMu, 2);

      
  /*
  for (Int_t index = 0; index < mathFit->GetNrefHistos(); index++)   {
    TH1* HREF = mathFit->GetRefHisto(index);
    Int_t bin = findBinWithinRange(HREF->GetXaxis(), xx[0]);
    summedError += TMath::Power(HREF->GetBinError(bin) * par[index] * par[mathFit->GetNrefHistos()], 2);
  }
  */
  return summedError;
}


//____________________________________________________________________________________________________________________
inline Double_t saveDivide(Double_t numerator, Double_t denominator) 
{
  return ((denominator != 0) ? numerator/denominator : 0 );
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfPionIntegral(TMatrixDSym covMat)
{
  return TMath::Sqrt(covMat(0, 0) + covMat(12, 12) + 2 * covMat(0, 12));
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfElectronFraction(Double_t* par, TMatrixDSym covMat) 
{
  Double_t g = saveDivide(par[3], (par[0] + par[1] + par[2] + 2 * par[3]));
  Double_t s1 = TMath::Power(g, 4) * (covMat(0, 0) + covMat(1, 1) + covMat(2, 2) + 4 * covMat(3, 3));
  Double_t s2 = TMath::Power(g, 2) * covMat(3, 3);
  Double_t s3 = (4 * TMath::Power(g, 4) - 2 * TMath::Power(g, 3)) * (covMat(3, 2) + covMat(3, 1) + covMat(3, 0));
  Double_t s4 = TMath::Power(g, 4) * 2 * (covMat(2, 1) + covMat(2, 0) +covMat(1, 0));
  
  return saveDivide(TMath::Sqrt(s1 + s2 + s3 + s4), par[3]);
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfKaonFraction(Double_t* par, TMatrixDSym covMat) 
{
  Double_t g = saveDivide(par[1], (par[0] + par[1] + par[2] + 2 * par[3]));
  Double_t s1 = TMath::Power(g, 4) * (covMat(0, 0) + covMat(1, 1) + covMat(2, 2) + 4 * covMat(3, 3));
  Double_t s2 = TMath::Power(g, 2) * covMat(1, 1);
  Double_t s3 = TMath::Power(g, 4) * (4 * covMat(3, 0) + 4 * covMat(3, 2) + 4 * covMat(3, 1) +
                                      2 * covMat(2, 1) + 2 * covMat(2, 0) + 2 * covMat(1, 0));
  Double_t s4 = TMath::Power(g, 3) * ((-4) * covMat(3, 1) - 2 * covMat(2, 1) - 2 * covMat(1, 0));
  
  return saveDivide(TMath::Sqrt(s1 + s2 + s3 + s4), par[1]);
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfPionFraction(Double_t* par, TMatrixDSym covMat) 
{
  Double_t g = saveDivide(par[0] + par[3], (par[0] + par[1] + par[2] + 2 * par[3]));
  Double_t s1 = TMath::Power(g, 4) * (covMat(0, 0) + covMat(1, 1) + covMat(2, 2) + 4 * covMat(3, 3));
  Double_t s2 = TMath::Power(g, 2) * (covMat(0, 0) + covMat(3, 3));
  Double_t s3 = TMath::Power(g, 4) * 2 * covMat(2, 1);
  Double_t s4 = (4 * TMath::Power(g, 4) - 2 * TMath::Power(g, 3)) * (covMat(3, 2) + covMat(3, 1));
  Double_t s5 = 2 * covMat(3, 0) * (2 * TMath::Power(g, 4) - 3 * TMath::Power(g, 3) + TMath::Power(g, 2));
  Double_t s6 = 2 * (covMat(2, 0) + covMat(1, 0)) * (TMath::Power(g, 4) - TMath::Power(g, 3));
  
  return saveDivide(TMath::Sqrt(s1 + s2 + s3 + s4 + s5 + s6), par[0] + par[3]);
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfProtonFraction(Double_t* par, TMatrixDSym covMat) 
{
  Double_t g = saveDivide(par[2], (par[0] + par[2] + par[1] + 2 * par[3]));
  Double_t s1 = TMath::Power(g, 4) * (covMat(0, 0) + covMat(1, 1) + covMat(2, 2) + 4 * covMat(3, 3));
  Double_t s2 = TMath::Power(g, 2) * covMat(2, 2);
  Double_t s3 = TMath::Power(g, 4) * (4 * covMat(3, 0) + 4 * covMat(3, 2) + 4 * covMat(3, 1) +
                                      2 * covMat(2, 1) + 2 * covMat(2, 0) + 2 * covMat(1, 0));
  Double_t s4 = TMath::Power(g, 3) * ((-4) * covMat(3, 2) - 2 * covMat(2, 1) - 2 * covMat(2, 0));
  
  return saveDivide(TMath::Sqrt(s1 + s2 + s3 + s4), par[2]);
}


//____________________________________________________________________________________________________________________
Double_t getErrorOfTotalIntegral(TMatrixDSym covMat) 
{
  Double_t s1 = covMat(0, 0) + covMat(1, 1) + covMat(2, 2) + 4 * covMat(3, 3);
  Double_t s2 = 4 * (covMat(3, 0) + covMat(3, 1) + covMat(3, 2));
  Double_t s3 = 2 * (covMat(2, 1) + covMat(2, 0) + covMat(1, 0));

  return TMath::Sqrt(s1 + s2 + s3);
}


//____________________________________________________________________________________________________________________
Double_t getMedianOfNonZeros(Double_t input[4])
{
  Double_t values[4] = {0,0,0,0};
  Int_t numNonZero = 0;
  if (input[0] > 0)  {
    values[numNonZero] = input[0];
    numNonZero++;
  }
  if (input[1] > 0)  {
    values[numNonZero] = input[1];
    numNonZero++;
  }
  if (input[2] > 0)  {
    values[numNonZero] = input[2];
    numNonZero++;
  }
  if (input[3] > 0)  {
    values[numNonZero] = input[3];
    numNonZero++;
  }
       
  return ((numNonZero > 0) ? TMath::Median(numNonZero, values) : 0);
}


//____________________________________________________________________________________________________________________
TCanvas* drawFinalFractions(Int_t mode, Double_t pLow, Double_t pHigh, Bool_t electronFixingUsed,
                            TH1* hFractionPions, TH1* hFractionKaons, TH1* hFractionProtons, TH1* hFractionElectrons, TH1* hFractionMuons,
                            TH1* hFractionSummed,
                            TH1* hFractionPionsMC, TH1* hFractionKaonsMC, TH1* hFractionProtonsMC, TH1* hFractionElectronsMC,
                            TH1* hFractionMuonsMC,
                            Bool_t plotIdentifiedSpectra,
                            Bool_t isMCtofPatchComparison = kFALSE)
{
  TCanvas* cFractions = new TCanvas(isMCtofPatchComparison ? "cFractionMCTOFpatched" : "cFractions",
                                    "Particle fractions", 100, 10, 1200, 800);
  cFractions->SetGridx(1);
  cFractions->SetGridy(1);
  cFractions->SetLogx(mode == kPMpT);
  hFractionPions->GetYaxis()->SetRangeUser(0.0, 1.0);
  SetReasonableAxisRange(hFractionPions->GetXaxis(), mode, pLow, pHigh);
  hFractionPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hFractionPions->GetXaxis()->SetNoExponent(kTRUE);
  hFractionPions->Draw("e p");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionPionsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionPionsMC->Draw("e p same");
  }
  
  SetReasonableAxisRange(hFractionKaons->GetXaxis(), mode, pLow, pHigh);
  hFractionKaons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionKaonsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionKaonsMC->Draw("e p same");
  }
  
  SetReasonableAxisRange(hFractionProtons->GetXaxis(), mode, pLow, pHigh);
  hFractionProtons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionProtonsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionProtonsMC->Draw("e p same");
  }
  
  SetReasonableAxisRange(hFractionElectrons->GetXaxis(), mode, pLow, pHigh);
  hFractionElectrons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionElectronsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionElectronsMC->Draw("e p same");
  }
  
  if (takeIntoAccountMuons) {
    SetReasonableAxisRange(hFractionMuons->GetXaxis(), mode, pLow, pHigh);
    if (muonFractionHandling != kNoMuons)
      hFractionMuons->Draw("e p same");
  }
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionMuonsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionMuonsMC->Draw("e p same");
  }
  
  hFractionSummed->Draw("e p same");
  
  if (fElectronFraction) {
    if (mode == kPMpT) 
      fElectronFraction->SetRange(lowFittingBoundElectronFraction, pHigh);
    else if (mode == kPMz)
      fElectronFraction->SetRange(lowFittingBoundElectronFraction, 1.0);
    else if (mode == kPMxi)
      fElectronFraction->SetRange(0., lowFittingBoundElectronFraction);
    else
      fElectronFraction->SetRange(0., 0.); // No special treatment of el for R and jT (seems to be the best option)
    
    if (electronFixingUsed)
      fElectronFraction->Draw("same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  if (plotIdentifiedSpectra)
    legend->SetNColumns(2);
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, isMCtofPatchComparison ? "TPC MC + TOF patch" : "Fit", "");
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, isMCtofPatchComparison ? "MC direct" : identifiedLabels[isMC].Data(), "");
  legend->AddEntry(hFractionPions, "#pi", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionPionsMC, "#pi", "p");
  legend->AddEntry(hFractionKaons, "K", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionKaonsMC, "K", "p");
  legend->AddEntry(hFractionProtons, "p", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionProtonsMC, "p", "p");
  legend->AddEntry(hFractionElectrons, "e", "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionElectronsMC, "e", "p");
  if (takeIntoAccountMuons && (muonFractionHandling != kNoMuons))
    legend->AddEntry(hFractionMuons, "#mu", "p");
  else
    legend->AddEntry((TObject*)0x0, "", "");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionMuonsMC, "#mu", "p");
  legend->AddEntry(hFractionSummed, "Total", "p");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(cFractions);
  
  return cFractions;
}


//____________________________________________________________________________________________________________________
TCanvas* drawFinalYields(Int_t mode, Double_t pLow, Double_t pHigh, 
                         TH1* hYieldPions, TH1* hYieldKaons, TH1* hYieldProtons, TH1* hYieldElectrons, TH1* hYieldMuons,
                         TH1* hYieldPionsMC, TH1* hYieldKaonsMC, TH1* hYieldProtonsMC, TH1* hYieldElectronsMC,
                         TH1* hYieldMuonsMC,
                         Bool_t plotIdentifiedSpectra)
{
  TCanvas* cYields = new TCanvas("cYields", "Particle yields",100,10,1200,800);
  cYields->SetGridx(1);
  cYields->SetGridy(1);
  cYields->SetLogx(mode == kPMpT);
  cYields->SetLogy(1);
  hYieldPions->GetYaxis()->SetRangeUser(hYieldElectrons->GetBinContent(hYieldElectrons->FindLastBinAbove(0.)) / 10.,
                                        hYieldPions->GetBinContent(hYieldPions->GetMaximumBin()) * 10.);
  SetReasonableAxisRange(hYieldPions->GetXaxis(), mode, pLow, pHigh);
  hYieldPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hYieldPions->GetXaxis()->SetNoExponent(kTRUE);
  hYieldPions->Draw("e p");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hYieldPionsMC->GetXaxis(), mode, pLow, pHigh);
    hYieldPionsMC->Draw("e p same");
  }
  
  SetReasonableAxisRange(hYieldKaons->GetXaxis(), mode, pLow, pHigh);
  hYieldKaons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hYieldKaonsMC->GetXaxis(), mode, pLow, pHigh);
    hYieldKaonsMC->Draw("e p same");
  }
  
  SetReasonableAxisRange(hYieldProtons->GetXaxis(), mode, pLow, pHigh);
  hYieldProtons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hYieldProtonsMC->GetXaxis(), mode, pLow, pHigh);
    hYieldProtonsMC->Draw("e p same");
  }
  
  if (takeIntoAccountMuons) {
    SetReasonableAxisRange(hYieldMuons->GetXaxis(), mode, pLow, pHigh);
    if (muonFractionHandling != kNoMuons)
      hYieldMuons->Draw("e p same");
    if (plotIdentifiedSpectra) {    
      SetReasonableAxisRange(hYieldMuonsMC->GetXaxis(), mode, pLow, pHigh);
      hYieldMuonsMC->Draw("e p same");
    }
  }
  
  SetReasonableAxisRange(hYieldElectrons->GetXaxis(), mode, pLow, pHigh);
  hYieldElectrons->Draw("e p same");
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hYieldElectronsMC->GetXaxis(), mode, pLow, pHigh);
    hYieldElectronsMC->Draw("e p same");
  }
  
  TLegend* legendYields = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legendYields->SetBorderSize(0);
  legendYields->SetFillColor(0);
  if (plotIdentifiedSpectra)
    legendYields->SetNColumns(2);
  if (plotIdentifiedSpectra)
    legendYields->AddEntry((TObject*)0x0, "Fit", "");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry((TObject*)0x0, identifiedLabels[isMC].Data(), "");
  legendYields->AddEntry(hYieldPions, "#pi", "p");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldPionsMC, "#pi", "p");
  legendYields->AddEntry(hYieldKaons, "K", "p");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldKaonsMC, "K", "p");
  legendYields->AddEntry(hYieldProtons, "p", "p");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldProtonsMC, "p", "p");
  legendYields->AddEntry(hYieldElectrons, "e", "p");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldElectronsMC, "e", "p");
  if (takeIntoAccountMuons && (muonFractionHandling != kNoMuons))
    legendYields->AddEntry(hYieldMuons, "#mu", "p");
  else
    legendYields->AddEntry((TObject*)0x0, "", "");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldMuonsMC, "#mu", "p");
  legendYields->Draw();
  
  ClearTitleFromHistoInCanvas(cYields);
  
  return cYields;
}
  
  
//____________________________________________________________________________________________________________________
TCanvas* drawFractionComparisonToMC(Int_t mode, Double_t pLow, Double_t pHigh, 
                                    TH1* hFractionPions, TH1* hFractionKaons, TH1* hFractionProtons, TH1* hFractionElectrons,
                                    TH1* hFractionMuons,
                                    TH1* hFractionPionsMC, TH1* hFractionKaonsMC, TH1* hFractionProtonsMC, TH1* hFractionElectronsMC,
                                    TH1* hFractionMuonsMC,
                                    TH1* hFractionComparisonPions, TH1* hFractionComparisonKaons, TH1* hFractionComparisonProtons,
                                    TH1* hFractionComparisonElectrons, TH1* hFractionComparisonMuons)
{
  // Compare data points with MC
  for (Int_t i = 1; i <= hFractionComparisonPions->GetNbinsX(); i++) {
    hFractionComparisonPions->SetBinContent(i, hFractionPions->GetBinContent(i));
    hFractionComparisonPions->SetBinError(i, hFractionPions->GetBinError(i));
    
    hFractionComparisonElectrons->SetBinContent(i, hFractionElectrons->GetBinContent(i));
    hFractionComparisonElectrons->SetBinError(i, hFractionElectrons->GetBinError(i));
    
    if (takeIntoAccountMuons) {
      hFractionComparisonMuons->SetBinContent(i, hFractionMuons->GetBinContent(i));
      hFractionComparisonMuons->SetBinError(i, hFractionMuons->GetBinError(i));
    }
    
    hFractionComparisonKaons->SetBinContent(i, hFractionKaons->GetBinContent(i));
    hFractionComparisonKaons->SetBinError(i, hFractionKaons->GetBinError(i));
    
    hFractionComparisonProtons->SetBinContent(i, hFractionProtons->GetBinContent(i));
    hFractionComparisonProtons->SetBinError(i, hFractionProtons->GetBinError(i));
  }
  
  RatioToRef(hFractionComparisonPions, hFractionPionsMC);
  RatioToRef(hFractionComparisonElectrons, hFractionElectronsMC);
  if (takeIntoAccountMuons)
    RatioToRef(hFractionComparisonMuons, hFractionMuonsMC);
  RatioToRef(hFractionComparisonKaons, hFractionKaonsMC);
  RatioToRef(hFractionComparisonProtons, hFractionProtonsMC);
  
  
  TCanvas* cFractionComparisons = new TCanvas("cFractionComparisons", "Particle fraction comparisons",100,10,1200,800);
  cFractionComparisons->SetGridx(1);
  cFractionComparisons->SetGridy(1);
  cFractionComparisons->SetLogx(mode == kPMpT);
  hFractionComparisonPions->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hFractionComparisonPions->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hFractionComparisonPions->GetXaxis()->SetNoExponent(kTRUE);
  hFractionComparisonPions->Draw("e p");
  
  hFractionComparisonElectrons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hFractionComparisonElectrons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonElectrons->Draw("e p same");
  
  if (takeIntoAccountMuons) {
    hFractionComparisonMuons->GetYaxis()->SetRangeUser(0.7, 1.3);
    SetReasonableAxisRange(hFractionComparisonMuons->GetXaxis(), mode, pLow, pHigh);
    if (muonFractionHandling != kNoMuons)
      hFractionComparisonMuons->Draw("e p same");
  }
  
  hFractionComparisonKaons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hFractionComparisonKaons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonKaons->Draw("e p same");
  
  hFractionComparisonProtons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hFractionComparisonProtons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonProtons->Draw("e p same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetNColumns(2);
  legend->AddEntry(hFractionComparisonPions, "#pi", "p");
  legend->AddEntry(hFractionComparisonKaons, "K", "p");
  legend->AddEntry(hFractionComparisonProtons, "p", "p");
  legend->AddEntry(hFractionComparisonElectrons, "e", "p");
  if (takeIntoAccountMuons && (muonFractionHandling != kNoMuons))
    legend->AddEntry(hFractionComparisonMuons, "#mu", "p");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(cFractionComparisons);
  
  return cFractionComparisons;
}


//____________________________________________________________________________________________________________________
TCanvas* drawTOFFractionComparisonToMC(Int_t mode, Double_t pLow, Double_t pHigh, TString fractionString,
                                       TH1* hYieldTOFTotal, TH1F** hYieldTOFTotalMC,
                                       TH1* hYieldTOFPions, TH1* hYieldTOFKaons, TH1* hYieldTOFProtons,
                                       TH1* hYieldTOFPionsMC, TH1* hYieldTOFKaonsMC, TH1* hYieldTOFProtonsMC, TH1* hYieldTOFElectronsMC,
                                       TH1* hYieldTOFMuonsMC,
                                       TH1F** hFractionTOFPions, TH1F** hFractionTOFKaons, TH1F** hFractionTOFProtons,
                                       TH1F** hFractionTOFPionsMC, TH1F** hFractionTOFKaonsMC, TH1F** hFractionTOFProtonsMC,
                                       TH1F** hFractionTOFElectronsMC, TH1F** hFractionTOFMuonsMC, TH1F** hFractionTOFPionsPlusMuonsMC,
                                       TH1F** hFractionTOFComparisonPions, TH1F** hFractionTOFComparisonKaons,
                                       TH1F** hFractionTOFComparisonProtons, TH1F** hFractionTOFComparisonPionsPlusMuons)
{
  TString fractionComparisonTitle = Form("Fraction fit / fraction %s", identifiedLabels[isMC].Data()); 
  
  (*hYieldTOFTotalMC) = new TH1F(*((TH1F*)hYieldTOFPionsMC));
  (*hYieldTOFTotalMC)->SetName("hYieldTOFTotalMC");
  (*hYieldTOFTotalMC)->GetYaxis()->SetTitle("Sum");
  (*hYieldTOFTotalMC)->SetLineColor(kBlack);
  (*hYieldTOFTotalMC)->SetMarkerColor(kBlack);
  (*hYieldTOFTotalMC)->Add(hYieldTOFKaonsMC);
  (*hYieldTOFTotalMC)->Add(hYieldTOFProtonsMC);
  (*hYieldTOFTotalMC)->Add(hYieldTOFMuonsMC);
  (*hYieldTOFTotalMC)->Add(hYieldTOFElectronsMC);
  
  // Binomial error for division, since subsets!
  (*hFractionTOFPions) = new TH1F(*((TH1F*)hYieldTOFPions));
  (*hFractionTOFPions)->SetName("hFractionTOFPions");
  (*hFractionTOFPions)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFPions)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFPions)->Divide(hYieldTOFPions, hYieldTOFTotal, 1., 1., "B");
  
  (*hFractionTOFKaons) = new TH1F(*((TH1F*)hYieldTOFKaons));
  (*hFractionTOFKaons)->SetName("hFractionTOFKaons");
  (*hFractionTOFKaons)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFKaons)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFKaons)->Divide(hYieldTOFKaons, hYieldTOFTotal, 1., 1., "B");
  
  (*hFractionTOFProtons) = new TH1F(*((TH1F*)hYieldTOFProtons));
  (*hFractionTOFProtons)->SetName("hFractionTOFProtons");
  (*hFractionTOFProtons)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFProtons)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFProtons)->Divide(hYieldTOFProtons, hYieldTOFTotal, 1., 1., "B");
  
  
  (*hFractionTOFPionsMC) = new TH1F(*((TH1F*)hYieldTOFPionsMC));
  (*hFractionTOFPionsMC)->SetName("hFractionTOFPionsMC");
  (*hFractionTOFPionsMC)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFPionsMC)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFPionsMC)->Divide(hYieldTOFPionsMC, (*hYieldTOFTotalMC), 1., 1., "B");
  
  (*hFractionTOFKaonsMC) = new TH1F(*((TH1F*)hYieldTOFKaonsMC));
  (*hFractionTOFKaonsMC)->SetName("hFractionTOFKaonsMC");
  (*hFractionTOFKaonsMC)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFKaonsMC)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFKaonsMC)->Divide(hYieldTOFKaonsMC, (*hYieldTOFTotalMC), 1., 1., "B");
  
  (*hFractionTOFProtonsMC) = new TH1F(*((TH1F*)hYieldTOFProtonsMC));
  (*hFractionTOFProtonsMC)->SetName("hFractionTOFProtonsMC");
  (*hFractionTOFProtonsMC)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFProtonsMC)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFProtonsMC)->Divide(hYieldTOFProtonsMC, (*hYieldTOFTotalMC), 1., 1., "B");
  
  (*hFractionTOFElectronsMC) = new TH1F(*((TH1F*)hYieldTOFElectronsMC));
  (*hFractionTOFElectronsMC)->SetName("hFractionTOFElectronsMC");
  (*hFractionTOFElectronsMC)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFElectronsMC)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFElectronsMC)->Divide(hYieldTOFElectronsMC, (*hYieldTOFTotalMC), 1., 1., "B");
  
  (*hFractionTOFMuonsMC) = new TH1F(*((TH1F*)hYieldTOFMuonsMC));
  (*hFractionTOFMuonsMC)->SetName("hFractionTOFMuonsMC");
  (*hFractionTOFMuonsMC)->GetYaxis()->SetTitle(fractionString.Data());
  (*hFractionTOFMuonsMC)->GetYaxis()->SetRangeUser(0.0, 1.0);
  (*hFractionTOFMuonsMC)->Divide(hYieldTOFMuonsMC, (*hYieldTOFTotalMC), 1., 1., "B");
  
  (*hFractionTOFPionsPlusMuonsMC) = new TH1F(*((TH1F*)(*hFractionTOFPionsMC)));
  (*hFractionTOFPionsPlusMuonsMC)->SetName("hFractionTOFPionsPlusMuonsMC");
  (*hFractionTOFPionsPlusMuonsMC)->SetTitle("#pi + #mu");
  (*hFractionTOFPionsPlusMuonsMC)->SetLineColor(kGray + 2);
  (*hFractionTOFPionsPlusMuonsMC)->SetMarkerColor(kGray + 2);
  (*hFractionTOFPionsPlusMuonsMC)->Add((*hFractionTOFMuonsMC));
  
  
  
  // Compare data points with MC
  (*hFractionTOFComparisonPions) = new TH1F(*((TH1F*)(*hFractionTOFPions)));
  (*hFractionTOFComparisonPions)->SetName("hFractionTOFComparisonPions");
  (*hFractionTOFComparisonPions)->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  RatioToRef((*hFractionTOFComparisonPions), (*hFractionTOFPionsMC));
  
  (*hFractionTOFComparisonKaons) = new TH1F(*((TH1F*)(*hFractionTOFKaons)));
  (*hFractionTOFComparisonKaons)->SetName("hFractionTOFComparisonKaons");
  (*hFractionTOFComparisonKaons)->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  RatioToRef((*hFractionTOFComparisonKaons), (*hFractionTOFKaonsMC));
  
  (*hFractionTOFComparisonProtons) = new TH1F(*((TH1F*)(*hFractionTOFProtons)));
  (*hFractionTOFComparisonProtons)->SetName("hFractionTOFComparisonProtons");
  (*hFractionTOFComparisonProtons)->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  RatioToRef((*hFractionTOFComparisonProtons), (*hFractionTOFProtonsMC));
  
  // NOTE: PionsPlusMuons means that ONLY on MC side these fractions are added (muons included in pions for data anyway)!
  (*hFractionTOFComparisonPionsPlusMuons) = new TH1F(*((TH1F*)(*hFractionTOFPions)));
  (*hFractionTOFComparisonPionsPlusMuons)->SetName("hFractionTOFComparisonPionsPlusMuons");
  (*hFractionTOFComparisonPionsPlusMuons)->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  RatioToRef((*hFractionTOFComparisonPionsPlusMuons), (*hFractionTOFPionsPlusMuonsMC));
  (*hFractionTOFComparisonPionsPlusMuons)->SetTitle((*hFractionTOFPionsPlusMuonsMC)->GetTitle());
  (*hFractionTOFComparisonPionsPlusMuons)->SetLineColor((*hFractionTOFPionsPlusMuonsMC)->GetLineColor());
  (*hFractionTOFComparisonPionsPlusMuons)->SetMarkerColor((*hFractionTOFPionsPlusMuonsMC)->GetMarkerColor());
  
  
  Double_t compRangeLow = 0.5;
  Double_t compRangeHigh = 1.5;
  TCanvas* cFractionTOFComparisons = new TCanvas("cFractionTOFComparisons", "Particle fraction (TOF) comparisons",100,10,1200,800);
  cFractionTOFComparisons->SetGridx(1);
  cFractionTOFComparisons->SetGridy(1);
  cFractionTOFComparisons->SetLogx(mode == kPMpT);
  (*hFractionTOFComparisonPions)->GetYaxis()->SetRangeUser(compRangeLow, compRangeHigh);
  SetReasonableAxisRange((*hFractionTOFComparisonPions)->GetXaxis(), mode, pLow, pHigh);
  (*hFractionTOFComparisonPions)->GetXaxis()->SetMoreLogLabels(kTRUE);
  (*hFractionTOFComparisonPions)->GetXaxis()->SetNoExponent(kTRUE);
  (*hFractionTOFComparisonPions)->Draw("e p");
  
  (*hFractionTOFComparisonKaons)->GetYaxis()->SetRangeUser(compRangeLow, compRangeHigh);
  SetReasonableAxisRange((*hFractionTOFComparisonKaons)->GetXaxis(), mode, pLow, pHigh);
  (*hFractionTOFComparisonKaons)->Draw("e p same");
  
  (*hFractionTOFComparisonProtons)->GetYaxis()->SetRangeUser(compRangeLow, compRangeHigh);
  SetReasonableAxisRange((*hFractionTOFComparisonProtons)->GetXaxis(), mode, pLow, pHigh);
  (*hFractionTOFComparisonProtons)->Draw("e p same");
  
  (*hFractionTOFComparisonPionsPlusMuons)->GetYaxis()->SetRangeUser(compRangeLow, compRangeHigh);
  SetReasonableAxisRange((*hFractionTOFComparisonPionsPlusMuons)->GetXaxis(), mode, pLow, pHigh);
  (*hFractionTOFComparisonPionsPlusMuons)->Draw("e p same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry((*hFractionTOFComparisonPions), "#pi", "p");
  legend->AddEntry((*hFractionTOFComparisonKaons), "K", "p");
  legend->AddEntry((*hFractionTOFComparisonProtons), "p", "p");
  legend->AddEntry((*hFractionTOFComparisonPionsPlusMuons), "#pi + #mu (for MC)", "p");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(cFractionTOFComparisons);
  
  return cFractionTOFComparisons;
}


//____________________________________________________________________________________________________________________
TCanvas* drawYieldComparisonToMC(Int_t mode, Double_t pLow, Double_t pHigh, 
                                 TH1* hYieldPions, TH1* hYieldKaons, TH1* hYieldProtons, TH1* hYieldElectrons,
                                 TH1* hYieldMuons,
                                 TH1* hYieldPionsMC, TH1* hYieldKaonsMC, TH1* hYieldProtonsMC, TH1* hYieldElectronsMC,
                                 TH1* hYieldMuonsMC,
                                 TH1* hYieldComparisonPions, TH1* hYieldComparisonKaons, TH1* hYieldComparisonProtons,
                                 TH1* hYieldComparisonElectrons, TH1* hYieldComparisonMuons) 
{
  // Compare data points with MC (yield)
  for (Int_t i = 1; i <= hYieldComparisonPions->GetNbinsX(); i++) {
    hYieldComparisonPions->SetBinContent(i, hYieldPions->GetBinContent(i));
    hYieldComparisonPions->SetBinError(i, hYieldPions->GetBinError(i));
    
    hYieldComparisonElectrons->SetBinContent(i, hYieldElectrons->GetBinContent(i));
    hYieldComparisonElectrons->SetBinError(i, hYieldElectrons->GetBinError(i));
    
    if (takeIntoAccountMuons) {
      hYieldComparisonMuons->SetBinContent(i, hYieldMuons->GetBinContent(i));
      hYieldComparisonMuons->SetBinError(i, hYieldMuons->GetBinError(i));
    }
    
    hYieldComparisonKaons->SetBinContent(i, hYieldKaons->GetBinContent(i));
    hYieldComparisonKaons->SetBinError(i, hYieldKaons->GetBinError(i));
    
    hYieldComparisonProtons->SetBinContent(i, hYieldProtons->GetBinContent(i));
    hYieldComparisonProtons->SetBinError(i, hYieldProtons->GetBinError(i));
  }
  
  RatioToRef(hYieldComparisonPions, hYieldPionsMC);
  RatioToRef(hYieldComparisonElectrons, hYieldElectronsMC);
  if (takeIntoAccountMuons)
    RatioToRef(hYieldComparisonMuons, hYieldMuonsMC);
  RatioToRef(hYieldComparisonKaons, hYieldKaonsMC);
  RatioToRef(hYieldComparisonProtons, hYieldProtonsMC);
  
  
  TCanvas* cYieldComparisons = new TCanvas("cYieldComparisons", "Particle yield comparisons",100,10,1200,800);
  cYieldComparisons->SetGridx(1);
  cYieldComparisons->SetGridy(1);
  cYieldComparisons->SetLogx(mode == kPMpT);
  hYieldComparisonPions->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hYieldComparisonPions->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hYieldComparisonPions->GetXaxis()->SetNoExponent(kTRUE);
  hYieldComparisonPions->Draw("e p");
  
  hYieldComparisonElectrons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hYieldComparisonElectrons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonElectrons->Draw("e p same");
  
  if (takeIntoAccountMuons) {
    hYieldComparisonMuons->GetYaxis()->SetRangeUser(0.7, 1.3);
    SetReasonableAxisRange(hYieldComparisonMuons->GetXaxis(), mode, pLow, pHigh);
    if (muonFractionHandling != kNoMuons)
      hYieldComparisonMuons->Draw("e p same");
  }
  
  hYieldComparisonKaons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hYieldComparisonKaons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonKaons->Draw("e p same");
  
  hYieldComparisonProtons->GetYaxis()->SetRangeUser(0.7, 1.3);
  SetReasonableAxisRange(hYieldComparisonProtons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonProtons->Draw("e p same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetNColumns(2);
  legend->AddEntry(hYieldComparisonPions, "#pi", "p");
  legend->AddEntry(hYieldComparisonKaons, "K", "p");
  legend->AddEntry(hYieldComparisonProtons, "p", "p");
  legend->AddEntry(hYieldComparisonElectrons, "e", "p");
  if (takeIntoAccountMuons && (muonFractionHandling != kNoMuons))
    legend->AddEntry(hYieldComparisonMuons, "#mu", "p");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(cYieldComparisons);
  
  return cYieldComparisons;
}


//____________________________________________________________________________________________________________________
TCanvas* drawFractionHistos(TString canvName, TString canvTitle, Int_t mode, Double_t pLow, Double_t pHigh, 
                            TH1* histDeltaPion, TH1* histDeltaElectron, TH1* histDeltaKaon, TH1* histDeltaProton, TH1* histMC,
                            Bool_t plotIdentifiedSpectra)
{
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  canv->SetLogx(mode == kPMpT);
  histDeltaPion->GetYaxis()->SetRangeUser(0.0, 1.0);
  histDeltaPion->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaPion->GetXaxis(), mode, pLow, pHigh);
  histDeltaPion->SetMarkerStyle(20);
  histDeltaPion->GetXaxis()->SetMoreLogLabels(kTRUE);
  histDeltaPion->GetXaxis()->SetNoExponent(kTRUE);
  histDeltaPion->Draw("e p");
  histDeltaElectron->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaElectron->GetXaxis(), mode, pLow, pHigh);
  histDeltaElectron->SetMarkerStyle(21);
  histDeltaElectron->Draw("e p same");
  histDeltaKaon->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaKaon->GetXaxis(), mode, pLow, pHigh);
  histDeltaKaon->SetMarkerStyle(22);
  histDeltaKaon->Draw("e p same");
  histDeltaProton->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaProton->GetXaxis(), mode, pLow, pHigh);
  histDeltaProton->SetMarkerStyle(29);
  histDeltaProton->Draw("e p same");
  if (plotIdentifiedSpectra) {
    histMC->GetYaxis()->SetTitle(canvTitle.Data());
    SetReasonableAxisRange(histMC->GetXaxis(), mode, pLow, pHigh);
    histMC->SetMarkerStyle(24);
    histMC->Draw("e p same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(histDeltaPion, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{#pi}}" : "_{#pi}"), "p");
  legend->AddEntry(histDeltaElectron, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{e}}" : "_{e}"), "p");
  legend->AddEntry(histDeltaKaon, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{K}}" : "_{K}"), "p");
  legend->AddEntry(histDeltaProton, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{p}}" : "_{p}"), "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(histMC, identifiedLabels[isMC].Data(), "p");
  legend->SetEntrySeparation(0.2);
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(canv);
  
  return canv;
}


//____________________________________________________________________________________________________________________
TCanvas* drawYieldHistos(TString canvName, TString canvTitle, Int_t mode, Double_t pLow, Double_t pHigh, 
                         TH1* histDeltaPion, TH1* histDeltaElectron, TH1* histDeltaKaon, TH1* histDeltaProton, TH1* histMC,
                         Bool_t plotIdentifiedSpectra)
{
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  canv->SetLogx(mode == kPMpT);
  canv->SetLogy(1);
  histDeltaPion->GetYaxis()->SetRangeUser(histDeltaPion->GetBinContent(histDeltaPion->FindLastBinAbove(0.)) / 10.,
                                          histDeltaPion->GetBinContent(histDeltaPion->GetMaximumBin()) * 10.);
  histDeltaPion->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaPion->GetXaxis(), mode, pLow, pHigh);
  histDeltaPion->SetMarkerStyle(20);
  histDeltaPion->GetXaxis()->SetMoreLogLabels(kTRUE);
  histDeltaPion->GetXaxis()->SetNoExponent(kTRUE);
  histDeltaPion->Draw("e p");
  histDeltaElectron->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaElectron->GetXaxis(), mode, pLow, pHigh);
  histDeltaElectron->SetMarkerStyle(21);
  histDeltaElectron->Draw("e p same");
  histDeltaKaon->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaKaon->GetXaxis(), mode, pLow, pHigh);
  histDeltaKaon->SetMarkerStyle(22);
  histDeltaKaon->Draw("e p same");
  histDeltaProton->GetYaxis()->SetTitle(canvTitle.Data());
  SetReasonableAxisRange(histDeltaProton->GetXaxis(), mode, pLow, pHigh);
  histDeltaProton->SetMarkerStyle(29);
  histDeltaProton->Draw("e p same");
  if (plotIdentifiedSpectra) {
    histMC->GetYaxis()->SetTitle(canvTitle.Data());
    SetReasonableAxisRange(histMC->GetXaxis(), mode, pLow, pHigh);
    histMC->SetMarkerStyle(24);
    histMC->Draw("e p same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  
  legend->AddEntry(histDeltaPion, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{#pi}}" : "_{#pi}"), "p");
  legend->AddEntry(histDeltaElectron, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{e}}" : "_{e}"), "p");
  legend->AddEntry(histDeltaKaon, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{K}}" : "_{K}"), "p");
  legend->AddEntry(histDeltaProton, Form("#Delta%s", useDeltaPrime ? "'_{#lower[-0.5]{p}}" : "_{p}"), "p");
  if (plotIdentifiedSpectra)
    legend->AddEntry(histMC, identifiedLabels[isMC].Data(), "p");
  legend->SetEntrySeparation(0.2);
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(canv);
  
  return canv;
}


//____________________________________________________________________________________________________________________
Int_t doSimultaneousFitRegularised(Int_t nPar, Double_t* gausParams, Double_t* parameterErrorsOut, Double_t* covMatrix,
                                   Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits, Double_t& reducedChiSquare,
                                   Double_t* specTemplateFractionErrorKa, Double_t* specTemplateFractionErrorPr)
{
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  Double_t chiSquare = -999;
  Int_t ndf = -1;
  
  AliTPCPIDmathFit::FitFunc_t* multiGaussFitArray = new AliTPCPIDmathFit::FitFunc_t[numSimultaneousFits];
  multiGaussFitArray[0] = multiGaussFitDeltaPi;  
  multiGaussFitArray[1] = multiGaussFitDeltaKa;  
  multiGaussFitArray[2] = multiGaussFitDeltaPr;  
  multiGaussFitArray[3] = multiGaussFitDeltaEl;  
  
  AliTPCPIDmathFit::FitFunc_t* errorOfFitHistosArray = new AliTPCPIDmathFit::FitFunc_t[numSimultaneousFits];
  errorOfFitHistosArray[0] = errorOfFitHistosDeltaPi;  
  errorOfFitHistosArray[1] = errorOfFitHistosDeltaKa;  
  errorOfFitHistosArray[2] = errorOfFitHistosDeltaPr;  
  errorOfFitHistosArray[3] = errorOfFitHistosDeltaEl;  
  
  //TODO errorFunction for bin errors of fit histos?  
  Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, 0x0, nPar, gausParams, parameterErrorsOut, covMatrix, 
                                     chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  //Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, errorOfFitHistosArray, nPar, gausParams, parameterErrorsOut, covMatrix, 
  //                                   chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  
  std::cout << std::endl;
  
  for (Int_t xBin = 0; xBin < mathFit->GetNumXbinsRegularisation(); xBin++) {
    std::cout << "x bin " << xBin << ":" << std::endl;
    
    Double_t sumFractions = 0;
    
    for (Int_t parIndex = xBin * mathFit->GetNumParametersPerXbin(); parIndex < (xBin + 1) * mathFit->GetNumParametersPerXbin();
         parIndex++) {
      Int_t parIndexModulo = parIndex % mathFit->GetNumParametersPerXbin();
    
      // NOTE: Covariance matrix is NOT set. But this doesn't matter since the parameter is fixed anyway, so
      // the error from the matrix would be zero.
      // parIndexModulo = 4 means muons, parIndexModulo = 3 means electrons, i.e. if parIndexModulo corresponds to muons,
      // then parIndexModulo - 1 corresponds to electrons.
      
      // Set electron fraction to value evaluated from a function above some threshold.
      // Fixed electron fraction < 0 does this job within the fitting functions
      if (parIndexModulo == 3 && gausParams[parIndex] < 0) {
        gausParams[parIndex]         = GetElectronFraction(-gausParams[parIndex], &gausParams[0], xBin);
        parameterErrorsOut[parIndex] = GetElectronFractionError();
      }
      // Set muon fraction equal to electron fraction (or some modified electron fraction) above some threshold,
      // which should be a reasonable approximation:
      // Fixed muon fraction < 0 does this job within the fitting functions
      else if (parIndexModulo == 4 && gausParams[parIndex] < 0) {
        gausParams[parIndex]         = GetMuonFractionFromElectronFractionAndPt(-gausParams[parIndex], gausParams[parIndex - 1]);
        parameterErrorsOut[parIndex] = parameterErrorsOut[parIndex - 1];
      }    
      
      // If special templates are used for K and/or p, set error accordingly
      if (parIndexModulo == 1 && TMath::Abs(lowParLimits[parIndex] - upParLimits[parIndex]) < epsilon) {
        // Kaons
        parameterErrorsOut[parIndex] = specTemplateFractionErrorKa[xBin];
      }
      if (parIndexModulo == 2 && TMath::Abs(lowParLimits[parIndex] - upParLimits[parIndex]) < epsilon) {
        // Protons
        parameterErrorsOut[parIndex] = specTemplateFractionErrorPr[xBin];
      }
      
      
      std::cout << "par[" << parIndex << "]: " << gausParams[parIndex] << " +- " << parameterErrorsOut[parIndex] << std::endl;
    
      if (parIndexModulo <= 3 || ((muonContamination || takeIntoAccountMuons) && parIndexModulo == 4))
        sumFractions += gausParams[parIndex];
    }
    
    std::cout << "Sum of fractions" << (muonContamination || takeIntoAccountMuons ? "(including muon contamination)" : "") << ": "
              << sumFractions; std::cout << std::endl;
    std::cout << std::endl << std::endl;
  }
  
  if (errFlag == 0) 
    std::cout << std::endl << "***Fit operation completed successfully***" << std::endl << std::endl;
  else
    std::cout << std::endl << "***Fit operation completed, but with errors***" << std::endl << std::endl;
  
  reducedChiSquare = (ndf > 0) ? chiSquare / ndf : -1;
  
  return errFlag;
}      


//____________________________________________________________________________________________________________________
Int_t doSimultaneousFit(TH1D** hDelta, Double_t xLow, Double_t xUp, Int_t nPar, Double_t* gausParams, Double_t* parameterErrorsOut, 
                        Double_t* covMatrix, Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits,
                        Double_t& reducedChiSquare, Double_t specTemplateFractionErrorKa, Double_t specTemplateFractionErrorPr)
{
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  Double_t chiSquare = -999;
  Int_t ndf = -1;
  
  //TODO:
  // Using no error on x (next TODO line) and no errorFunction (next after next TODO line) sometimes gives good results.
  // However, it can completely fail for low statistics for the fit histos.
  // Using either an error on x or the errorFunction both gives reasonable results, but might be slightly worse results in some cases
  // (shifted/distorted data). Maybe: Choose one method - the rest is for systematic errors of this fitting
  
  //TODO The next TODO marks are only relevant for chiSquare, but not for loglikelihood
  //TODO Use error in x also -> If reference histos have low statistics, this will be very important
  
  //TODO FOR INDIVIDUAL DeltaPrime fits: Replace input data and arrays in the following with the desired species and set numSimultaneousFits to 1
  for (Int_t i = 0; i < numSimultaneousFits; i++) {
    mathFit->InputData(hDelta[i], 0, i, xLow, xUp, -1., kFALSE); 
    //mathFit->InputData(hDelta[i], 0, i, xLow, xUp, -1., kTRUE); 
  }
  
  AliTPCPIDmathFit::FitFunc_t* multiGaussFitArray = new AliTPCPIDmathFit::FitFunc_t[numSimultaneousFits];
  multiGaussFitArray[0] = multiGaussFitDeltaPi;  
  multiGaussFitArray[1] = multiGaussFitDeltaKa;  
  multiGaussFitArray[2] = multiGaussFitDeltaPr;  
  multiGaussFitArray[3] = multiGaussFitDeltaEl;  
  
  AliTPCPIDmathFit::FitFunc_t* errorOfFitHistosArray = new AliTPCPIDmathFit::FitFunc_t[numSimultaneousFits];
  errorOfFitHistosArray[0] = errorOfFitHistosDeltaPi;  
  errorOfFitHistosArray[1] = errorOfFitHistosDeltaKa;  
  errorOfFitHistosArray[2] = errorOfFitHistosDeltaPr;  
  errorOfFitHistosArray[3] = errorOfFitHistosDeltaEl;  
  
  //TODO errorFunction for bin errors of fit histos?  
  Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, 0x0, nPar, gausParams, parameterErrorsOut, covMatrix, 
                                     chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  //Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, errorOfFitHistosArray, nPar, gausParams, parameterErrorsOut, covMatrix, 
  //                                   chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  
  std::cout << std::endl;
  
  // If the electron fraction is fixed, evaluate the error of the extrapolation of the fixed value
  if (TMath::Abs(lowParLimits[3] - upParLimits[3]) < epsilon) {
    // NOTE: Covariance matrix is NOT set. But this doesn't matter since the parameter is fixed anyway, so
    // the error from the matrix would be zero
    parameterErrorsOut[3] = GetElectronFractionError();
  }
  
  // Set muon fraction equal to electron fraction (or some modified electron fraction) above some threshold,
  // which should be a reasonable approximation:
  // Fixed muon fraction < 0 does this job within the fitting functions
  if (gausParams[4] < 0 ) {
    // NOTE: Covariance matrix is NOT set. But this doesn't matter since the parameter is fixed anyway, so
    // the error from the matrix would be zero
    gausParams[4] = GetMuonFractionFromElectronFractionAndPt(-gausParams[4], gausParams[3]);
    parameterErrorsOut[4] = parameterErrorsOut[3];
  }
  
  // If special templates are used for K and/or p, set error accordingly
  if (TMath::Abs(lowParLimits[1] - upParLimits[1]) < epsilon) {
    // Kaons
    parameterErrorsOut[1] = specTemplateFractionErrorKa;
  }
  if (TMath::Abs(lowParLimits[2] - upParLimits[2]) < epsilon) {
    // Protons
    parameterErrorsOut[2] = specTemplateFractionErrorPr;
  }
  
  Double_t sumFractions = 0;
  for (Int_t parIndex = 0; parIndex < nPar; parIndex++) {
    std::cout << "par[" << parIndex << "]: " << gausParams[parIndex] << " +- " << parameterErrorsOut[parIndex] << std::endl;
  }
  sumFractions = gausParams[0] + gausParams[1] + gausParams[2] + gausParams[3];
  // In case of muon contamination add muon fraction also
  if (muonContamination || takeIntoAccountMuons) {
    sumFractions += gausParams[4];
  }
  
  std::cout << "Sum of fractions" << (muonContamination || takeIntoAccountMuons ? "(including muon contamination)" : "") << ": " << sumFractions; std::cout << std::endl;
  
  if (errFlag == 0) 
    std::cout << std::endl << "***Fit operation completed successfully***" << std::endl << std::endl;
  else
    std::cout << std::endl << "***Fit operation completed, but with errors***" << std::endl << std::endl;
  
  reducedChiSquare = (ndf > 0) ? chiSquare / ndf : -1;
  
  return errFlag;
}        


//____________________________________________________________________________________________________________________
Int_t doFit(TH1D* hDelta, Double_t xLow, Double_t xUp, Int_t nPar, Double_t* gausParams, Double_t* parameterErrorsOut, Double_t* covMatrix,
            Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits, TF1* totalDeltaSpecies, Double_t& reducedChiSquare,
            Double_t specTemplateFractionErrorKa, Double_t specTemplateFractionErrorPr)
{
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  Double_t chiSquare = -999;
  Int_t ndf = -1;
  
  //TODO:
  // Using no error on x (next TODO line) and no errorFunction (next after next TODO line) sometimes gives good results.
  // However, it can completely fail for low statistics for the fit histos.
  // Using either an error on x or the errorFunction both gives reasonable results, but might be slightly worse results in some cases
  // (shifted/distorted data). Maybe: Choose one method - the rest is for systematic errors of this fitting
  
  //TODO The next TODO marks are only relevant for chiSquare, but not for loglikelihood
  //TODO Use error in x also -> If reference histos have low statistics, this will be very important
  mathFit->InputData(hDelta, 0, 0, xLow, xUp, -1., kFALSE); 
  //mathFit->InputData(hDelta, 0, 0, xLow, xUp, -1., kTRUE); 
  
  AliTPCPIDmathFit::FitFunc_t* multiGaussFitArray = new AliTPCPIDmathFit::FitFunc_t[1];
  multiGaussFitArray[0] = multiGaussFit; 
  
  AliTPCPIDmathFit::FitFunc_t* errorOfFitHistosArray = new AliTPCPIDmathFit::FitFunc_t[1];
  errorOfFitHistosArray[0] = errorOfFitHistos;
  
  //TODO errorFunction for bin errors of fit histos?
  Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, 0x0, nPar, gausParams, parameterErrorsOut, covMatrix, 
                                     chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  //Int_t errFlag = mathFit->MinuitFit(multiGaussFitArray, errorOfFitHistosArray, nPar, gausParams, parameterErrorsOut, covMatrix, 
  //                                   chiSquare, ndf, stepSize, lowParLimits, upParLimits);
  
  // If the electron fraction is fixed, evaluate the error of the extrapolation of the fixed value
  if (TMath::Abs(lowParLimits[3] - upParLimits[3]) < epsilon) {
    // NOTE: Covariance matrix is NOT set. But this doesn't matter since the parameter is fixed anyway, so
    // the error from the matrix would be zero
    parameterErrorsOut[3] = GetElectronFractionError();
  }
  
  // Set muon fraction equal to electron fraction (or some modified electron fraction) above some threshold, which should be a reasonable approximation:
  // Fixed muon fraction < 0 does this job within the fitting functions
  if (gausParams[4] < 0 ) {
    // NOTE: Covariance matrix is NOT set. But this doesn't matter since the parameter is fixed anyway, so
    // the error from the matrix would be zero
    gausParams[4] = GetMuonFractionFromElectronFractionAndPt(-gausParams[4], gausParams[3]);
    parameterErrorsOut[4] = parameterErrorsOut[3];
  }
  
  // If special templates are used for K and/or p, set error accordingly
  if (TMath::Abs(lowParLimits[1] - upParLimits[1]) < epsilon) {
    // Kaons
    parameterErrorsOut[1] = specTemplateFractionErrorKa;
  }
  if (TMath::Abs(lowParLimits[2] - upParLimits[2]) < epsilon) {
    // Protons
    parameterErrorsOut[2] = specTemplateFractionErrorPr;
  }
  
  Double_t sumFractions = 0;
  for (Int_t parIndex = 0; parIndex < nPar; parIndex++) {
    std::cout << totalDeltaSpecies->GetParName(parIndex) << ": " << gausParams[parIndex] << " +- " << parameterErrorsOut[parIndex] << std::endl;
  }
  sumFractions = gausParams[0] + gausParams[1] + gausParams[2] + gausParams[3];
  // In case of muon contamination add muon fraction also
  if (muonContamination || takeIntoAccountMuons) {
    sumFractions += gausParams[4];
  }
  
  std::cout << "Sum of fractions" << (muonContamination || takeIntoAccountMuons ? "(including muon contamination)" : "") << ": " << sumFractions;
  std::cout << std::endl;
  
  if (errFlag == 0) 
    std::cout << std::endl << "***Fit operation completed successfully***" << std::endl << std::endl;
  else
    std::cout << std::endl << "***Fit operation completed, but with errors***" << std::endl << std::endl;
  
  for (Int_t parIndex = 0; parIndex < nPar; parIndex++) {
    totalDeltaSpecies->SetParameter(parIndex, gausParams[parIndex]);
    totalDeltaSpecies->SetParError(parIndex, parameterErrorsOut[parIndex]);
  }
  
  reducedChiSquare = (ndf > 0) ? chiSquare / ndf : -1;
  
  return errFlag;
}        


//____________________________________________________________________________________________________________________
Double_t setFractionsAndYields(Int_t slice, Double_t inverseBinWidth, Int_t species, Double_t* parametersOut, 
                               Double_t* parameterErrorsOut, TH1* hFractionSpecies, TH1* hFractionPionsDeltaSpecies, 
                               TH1* hFractionElectronsDeltaSpecies, TH1* hFractionKaonsDeltaSpecies, TH1* hFractionProtonsDeltaSpecies,
                               TH1* hFractionMuonsDeltaSpecies, TH1* hYieldSpecies, TH1* hYieldPionsDeltaSpecies,
                               TH1* hYieldElectronsDeltaSpecies, TH1* hYieldKaonsDeltaSpecies, TH1* hYieldProtonsDeltaSpecies,
                               TH1* hYieldMuonsDeltaSpecies, 
                               Bool_t normaliseFractions = kFALSE)
{
  // Set fraction and yields in corresponding histograms. If normaliseFractions is kTRUE, the fractions will be normalised to unity
  // and the normalisation factor will be returned (i.e. 1./sumFraction)
  
  Double_t normalisationFactor = 1.0;
  
  // Since a log likelihood fit is anyway used, the normalisation should give a factor close to unity
  if (normaliseFractions) {
    Double_t sumFractions = parametersOut[0] + (muonContamination ? parametersOut[3] : 0) + parametersOut[1] + parametersOut[2] +
                            parametersOut[3] + (takeIntoAccountMuons ? parametersOut[4] : 0.);
    if (sumFractions > 0) {
      normalisationFactor = 1./sumFractions;
      for (Int_t i = 0; i < 5; i++) {
        parametersOut[i] *= normalisationFactor;
        
        // Do not introduce an error for the normalisation, i.e. just scale parameters and fractions with the same factor which is 
        // assumed to be exact.
        // Note that correlations should already be included in the parameterError        
        parameterErrorsOut[i] *= normalisationFactor;
      }
    }
  }
  
  Double_t sumOfParticles = inverseBinWidth * parametersOut[5];
  
  if (species == kPi) {
    hFractionSpecies->SetBinContent(slice + 1, (parametersOut[0]+(muonContamination ? parametersOut[3] : 0)));
    hFractionSpecies->SetBinError(slice + 1, parameterErrorsOut[0]);
  }
  else if (species == kEl) {
    hFractionSpecies->SetBinContent(slice + 1, parametersOut[3]);
    hFractionSpecies->SetBinError(slice + 1, parameterErrorsOut[3]);
  }
  else if (species == kKa) {    
    hFractionSpecies->SetBinContent(slice + 1, parametersOut[1]);
    hFractionSpecies->SetBinError(slice + 1, parameterErrorsOut[1]);
  }
  else if (species == kPr) {    
    hFractionSpecies->SetBinContent(slice + 1, parametersOut[2]);
    hFractionSpecies->SetBinError(slice + 1, parameterErrorsOut[2]);
  }
  else if (species == kMu) {
    if (takeIntoAccountMuons) {    
      hFractionSpecies->SetBinContent(slice + 1, parametersOut[4]);
      hFractionSpecies->SetBinError(slice + 1, parameterErrorsOut[4]);
      
      hYieldSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionSpecies->GetBinContent(slice + 1));
      hYieldSpecies->SetBinError(slice + 1, sumOfParticles * hFractionSpecies->GetBinError(slice + 1));
    }
    
    // Only set these histos for muons. The DeltaSpecies histos for muons will be set together with all other species
    return normalisationFactor;
  }
  
  hFractionPionsDeltaSpecies->SetBinContent(slice + 1, (parametersOut[0]+(muonContamination ? parametersOut[3] : 0)));
  hFractionPionsDeltaSpecies->SetBinError(slice + 1, parameterErrorsOut[0]);//TODO What about error of parOut[3]?
  hFractionElectronsDeltaSpecies->SetBinContent(slice + 1, parametersOut[3]);
  hFractionElectronsDeltaSpecies->SetBinError(slice + 1, parameterErrorsOut[3]);
  hFractionKaonsDeltaSpecies->SetBinContent(slice + 1, parametersOut[1]);
  hFractionKaonsDeltaSpecies->SetBinError(slice + 1, parameterErrorsOut[1]);
  hFractionProtonsDeltaSpecies->SetBinContent(slice + 1, parametersOut[2]);
  hFractionProtonsDeltaSpecies->SetBinError(slice + 1, parameterErrorsOut[2]);
  if (takeIntoAccountMuons) {
    hFractionMuonsDeltaSpecies->SetBinContent(slice + 1, parametersOut[4]);
    hFractionMuonsDeltaSpecies->SetBinError(slice + 1, parameterErrorsOut[4]);
  }
  
  hYieldSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionSpecies->GetBinContent(slice + 1));
  hYieldSpecies->SetBinError(slice + 1, sumOfParticles * hFractionSpecies->GetBinError(slice + 1));
  
  hYieldPionsDeltaSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionPionsDeltaSpecies->GetBinContent(slice + 1));
  hYieldPionsDeltaSpecies->SetBinError(slice + 1, sumOfParticles * hFractionPionsDeltaSpecies->GetBinError(slice + 1));
  hYieldElectronsDeltaSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionElectronsDeltaSpecies->GetBinContent(slice + 1));
  hYieldElectronsDeltaSpecies->SetBinError(slice + 1, sumOfParticles * hFractionElectronsDeltaSpecies->GetBinError(slice + 1));
  hYieldKaonsDeltaSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionKaonsDeltaSpecies->GetBinContent(slice + 1));
  hYieldKaonsDeltaSpecies->SetBinError(slice + 1, sumOfParticles * hFractionKaonsDeltaSpecies->GetBinError(slice + 1));
  hYieldProtonsDeltaSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionProtonsDeltaSpecies->GetBinContent(slice + 1));
  hYieldProtonsDeltaSpecies->SetBinError(slice + 1, sumOfParticles * hFractionProtonsDeltaSpecies->GetBinError(slice + 1));
  if (takeIntoAccountMuons) {
    hYieldMuonsDeltaSpecies->SetBinContent(slice + 1, sumOfParticles * hFractionMuonsDeltaSpecies->GetBinContent(slice + 1));
    hYieldMuonsDeltaSpecies->SetBinError(slice + 1, sumOfParticles * hFractionMuonsDeltaSpecies->GetBinError(slice + 1));
  }
  
  return normalisationFactor;
}


//____________________________________________________________________________________________________________________
void FitElContaminationForTOFpions(TAxis* pTaxis, TH2D* hGenDeltaUsed[][AliPID::kSPECIES + 1], TH1D* hMCmuonsAndPionsDummy, TH1D** hDeltaPi,
                                   TH1D** hDeltaEl, TH1D** hDeltaKa, TH1D** hDeltaPr, TH1D* hDeltaPiMC[][AliPID::kSPECIES], 
                                   TH1D* hDeltaElMC[][AliPID::kSPECIES], TH1D* hDeltaKaMC[][AliPID::kSPECIES], 
                                   TH1D* hDeltaPrMC[][AliPID::kSPECIES], Int_t mode, Bool_t useLogLikelihood, 
                                   Bool_t useWeightsForLogLikelihood, Bool_t plotIdentifiedSpectra, Int_t pSliceLow, Int_t pSliceHigh, 
                                   Int_t numSlices, Bool_t restrictJetPtAxis, Double_t actualUpperJetPt, Double_t xLow, Double_t xUp, 
                                   Int_t nBins, TString* speciesLabel, TFile* saveF, TH1D* hYieldTOFOrigBinningPions, 
                                   TH1D* hYieldTOFOrigBinningElectrons)
{
  // Take everything from the TOF pion bin: These are mainly pions, but also some electrons (also muons, but they need to be
  // treated separately anyway). The contamination from other particles should be negligable due to the choice of TOF cuts.
  // Fit the dEdx distributions of these particles with the corresponding el and pi templates (to use the same fit code, just
  // set the fraction of the other species fix to zero). Do NOT apply regularisation since the TOF efficiency might change. Also,
  // this is not required, since there is no el pi crossing in the region with finite TOF efficiency.
  //
  // Since the separation is good in the region where TOF is effective, the error of the fitting is negligable (or will be captured
  // when the parametrisations are changed).
  // This way, it is possible to benefit from the high TOF efficiency w/o el exclusion, but at the same time get rid of possible el
  // contamination
  //
  // NOTE: hYieldTOFOrigBinningPions has the same binning as hYieldPt
  
  std::cout << std::endl << std::endl << "Fitting electron contamination of TOF pions:" << std::endl
            << "Simultaneous fit" << std::endl << "No regularisation" << std::endl << "No TOF patching for this sample" << std::endl
            << "Fit method fixed to 2" << std::endl << "All other settings as for main fit" << std::endl << std::endl;
  
  // Fracs of each species + total yield in x bin
  const Int_t nParSimultaneousFit = AliPID::kSPECIES + 1; 
  
  Double_t parameterErrorsOut[nParSimultaneousFit] = { 0., };
  
  
  // Setup mathFit for simultaneous fit without regularisation and TOF patching; using fit method 2
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance(1, 4, 1810);
  mathFit->SetDebugLevel(0); 
  mathFit->SetEpsilon(5e-05);
  mathFit->SetMaxCalls(1e8);
  mathFit->SetMinimisationStrategy(minimisationStrategy);
  mathFit->SetUseLogLikelihood(useLogLikelihood);
  mathFit->SetUseWeightsForLogLikelihood(useWeightsForLogLikelihood);
  
  mathFit->SetScaleFactorError(2.);
  
  TH1D* hDeltaPiFitQA[numSlices];
  TH1D* hDeltaElFitQA[numSlices];
  TH1D* hDeltaKaFitQA[numSlices];
  TH1D* hDeltaPrFitQA[numSlices];
  
  for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < pTaxis->GetNbins(); slice++) {   
    if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
      continue; 
      
    // There won't (actually: shouldn't) be tracks with a pT larger than the jet pT
    if (mode == kPMpT && restrictJetPtAxis && binsPt[slice] >= actualUpperJetPt)
      continue;
      
    if (mode == kPMpT)
      std::cout << "Fitting range " << binsPt[slice] << " GeV/c < Pt < " << binsPt[slice + 1] << " GeV/c..." << std::endl;
    else {
      std::cout << "Fitting range " << pTaxis->GetBinLowEdge(slice + 1) << " < " << modeShortName[mode].Data() << " < ";
      std::cout << pTaxis->GetBinUpEdge(slice + 1) << "..." << std::endl;
    } 
    
    // Add/subtract some very small offset to be sure not to sit on the bin boundary, when looking for the integration/projection limits.
    const Int_t pBinLowProjLimit = (mode == kPMpT) ? hYieldTOFOrigBinningPions->GetXaxis()->FindFixBin(binsPt[slice] + 1e-5)    : slice + 1;
    const Int_t pBinUpProjLimit  = (mode == kPMpT) ? hYieldTOFOrigBinningPions->GetXaxis()->FindFixBin(binsPt[slice + 1]- 1e-5) : slice + 1;
      
    const Double_t totalYield = hYieldTOFOrigBinningPions->Integral(pBinLowProjLimit, pBinUpProjLimit);
    
    if (totalYield <= 0) {
      std::cout << "Skipped bin (yield is zero)!" << std::endl;
      continue;
    }
    
    if (totalYield < gYieldThresholdForFitting) {
      std::cout << "Skipped bin (yield (" << totalYield << ") is smaller than threshold (" << gYieldThresholdForFitting << "))!" << std::endl;
      continue;
    }
    
    TH1D *hGenDeltaElForElProj = 0x0, *hGenDeltaKaForElProj = 0x0, *hGenDeltaPiForElProj = 0x0, *hGenDeltaPrForElProj = 0x0;
    TH1D *hGenDeltaElForKaProj = 0x0, *hGenDeltaKaForKaProj = 0x0, *hGenDeltaPiForKaProj = 0x0, *hGenDeltaPrForKaProj = 0x0;
    TH1D *hGenDeltaElForPiProj = 0x0, *hGenDeltaKaForPiProj = 0x0, *hGenDeltaPiForPiProj = 0x0, *hGenDeltaPrForPiProj = 0x0;
    TH1D *hGenDeltaElForMuProj = 0x0, *hGenDeltaKaForMuProj = 0x0, *hGenDeltaPiForMuProj = 0x0, *hGenDeltaPrForMuProj = 0x0;
    TH1D *hGenDeltaElForPrProj = 0x0, *hGenDeltaKaForPrProj = 0x0, *hGenDeltaPiForPrProj = 0x0, *hGenDeltaPrForPrProj = 0x0;

    hGenDeltaElForElProj =(TH1D*)hGenDeltaUsed[kEl][kEl]->ProjectionY(Form("hGenDeltaElForElTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaKaForElProj =(TH1D*)hGenDeltaUsed[kKa][kEl]->ProjectionY(Form("hGenDeltaKaForElTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPiForElProj =(TH1D*)hGenDeltaUsed[kPi][kEl]->ProjectionY(Form("hGenDeltaPiForElTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPrForElProj =(TH1D*)hGenDeltaUsed[kPr][kEl]->ProjectionY(Form("hGenDeltaPrForElTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        
    hGenDeltaElForKaProj =(TH1D*)hGenDeltaUsed[kEl][kKa]->ProjectionY(Form("hGenDeltaElForKaTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaKaForKaProj =(TH1D*)hGenDeltaUsed[kKa][kKa]->ProjectionY(Form("hGenDeltaKaForKaTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPiForKaProj =(TH1D*)hGenDeltaUsed[kPi][kKa]->ProjectionY(Form("hGenDeltaPiForKaTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPrForKaProj =(TH1D*)hGenDeltaUsed[kPr][kKa]->ProjectionY(Form("hGenDeltaPrForKaTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
    hGenDeltaElForPiProj =(TH1D*)hGenDeltaUsed[kEl][kPi]->ProjectionY(Form("hGenDeltaElForPiTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaKaForPiProj =(TH1D*)hGenDeltaUsed[kKa][kPi]->ProjectionY(Form("hGenDeltaKaForPiTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPiForPiProj =(TH1D*)hGenDeltaUsed[kPi][kPi]->ProjectionY(Form("hGenDeltaPiForPiTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPrForPiProj =(TH1D*)hGenDeltaUsed[kPr][kPi]->ProjectionY(Form("hGenDeltaPrForPiTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        
    hGenDeltaElForMuProj =(TH1D*)hGenDeltaUsed[kEl][kMu]->ProjectionY(Form("hGenDeltaElForMuTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaKaForMuProj =(TH1D*)hGenDeltaUsed[kKa][kMu]->ProjectionY(Form("hGenDeltaKaForMuTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPiForMuProj =(TH1D*)hGenDeltaUsed[kPi][kMu]->ProjectionY(Form("hGenDeltaPiForMuTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPrForMuProj =(TH1D*)hGenDeltaUsed[kPr][kMu]->ProjectionY(Form("hGenDeltaPrForMuTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        
    hGenDeltaElForPrProj =(TH1D*)hGenDeltaUsed[kEl][kPr]->ProjectionY(Form("hGenDeltaElForPrTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaKaForPrProj =(TH1D*)hGenDeltaUsed[kKa][kPr]->ProjectionY(Form("hGenDeltaKaForPrTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPiForPrProj =(TH1D*)hGenDeltaUsed[kPi][kPr]->ProjectionY(Form("hGenDeltaPiForPrTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
    hGenDeltaPrForPrProj =(TH1D*)hGenDeltaUsed[kPr][kPr]->ProjectionY(Form("hGenDeltaPrForPrTOFpiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        
      
    // Normalise generated histos to TOTAL number of GENERATED particles for this species (i.e. including
    // entries that lie in the under or overflow bin), so that situations in which the generated spectra lie
    // at least partly outside the histo are treated properly. To find the total number of generated particle
    // species X, one can just take the integral of the generated histo for DeltaX (which should include all
    // generated entries) and apply the same normalisation factor to all other DeltaSpecies.
    // Also set some cosmetics

    
    // Generated electrons
    // If template is empty, it is empty for all delta_x!
    Bool_t nonEmptyTemplateEl = kFALSE;
    Double_t normEl = normaliseHist(hGenDeltaElForElProj, nonEmptyTemplateEl, -1);
    normaliseHist(hGenDeltaKaForElProj, normEl);
    normaliseHist(hGenDeltaPiForElProj, normEl);
    normaliseHist(hGenDeltaPrForElProj, normEl);
    
    
    // Generated kaons
    Double_t normKa = normaliseHist(hGenDeltaKaForKaProj, -1);
    normaliseHist(hGenDeltaElForKaProj, normKa);
    normaliseHist(hGenDeltaPiForKaProj, normKa);
    normaliseHist(hGenDeltaPrForKaProj, normKa);
    
    
    // Generated pions
    // If template is empty, it is empty for all delta_x!
    Bool_t nonEmptyTemplatePi = kFALSE;
    Double_t normPi = normaliseHist(hGenDeltaPiForPiProj, nonEmptyTemplatePi, -1);
    normaliseHist(hGenDeltaElForPiProj, normPi);
    normaliseHist(hGenDeltaKaForPiProj, normPi);
    normaliseHist(hGenDeltaPrForPiProj, normPi);
    

    Double_t normMu = 1;
    if (takeIntoAccountMuons) {
      // Generated pions
      // Since masses of muons and pions are so similar, the normalisation scheme should still work when looking at deltaPion instead
      normMu = normaliseHist(hGenDeltaPiForMuProj, -1);
      normaliseHist(hGenDeltaElForMuProj, normMu);
      normaliseHist(hGenDeltaKaForMuProj, normMu);
      normaliseHist(hGenDeltaPrForMuProj, normMu);
    }
    
    
    // Generated protons
    Double_t normPr = normaliseHist(hGenDeltaPrForPrProj, -1);
    normaliseHist(hGenDeltaElForPrProj, normPr);
    normaliseHist(hGenDeltaKaForPrProj, normPr);
    normaliseHist(hGenDeltaPiForPrProj, normPr);
      
    TF1* totalDeltaPion = 0x0;
    TF1* totalDeltaKaon = 0x0;
    TF1* totalDeltaProton = 0x0;
    TF1* totalDeltaElectron = 0x0;
    
    TLegend* legend = 0x0;

    totalDeltaPion = new TF1(Form("totalDeltaPion_TOFpi_slice%d", slice), multiGaussFitDeltaPi, xLow, xUp, nParSimultaneousFit);
    setUpFitFunction(totalDeltaPion, nBins);
    
    totalDeltaKaon = new TF1(Form("totalDeltaKaon_TOFpi_slice%d", slice), multiGaussFitDeltaKa, xLow, xUp, nParSimultaneousFit);
    setUpFitFunction(totalDeltaKaon, nBins);
    
    totalDeltaProton = new TF1(Form("totalDeltaProton_TOFpi_slice%d", slice), multiGaussFitDeltaPr, xLow, xUp, nParSimultaneousFit);
    setUpFitFunction(totalDeltaProton, nBins);
    
    totalDeltaElectron = new TF1(Form("totalDeltaElectron_TOFpi_slice%d", slice), multiGaussFitDeltaEl, xLow, xUp, nParSimultaneousFit);
    setUpFitFunction(totalDeltaElectron, nBins);
    
    // Legend is the same for all \Delta "species" plots
    legend = new TLegend(0.722126, 0.605932, 0.962069, 0.925932);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    if (plotIdentifiedSpectra)
      legend->SetNColumns(2);
    legend->AddEntry((TObject*)0x0, "Fit", "");
    if (plotIdentifiedSpectra)
      legend->AddEntry((TObject*)0x0, identifiedLabels[isMC].Data(), "");
    
    legend->AddEntry(hDeltaPi[slice], "Data", "Lp");
    if (plotIdentifiedSpectra)
      legend->AddEntry((TObject*)0x0, "", "");
    
    legend->AddEntry(totalDeltaPion, "Multi-template fit", "L");
    if (plotIdentifiedSpectra)
      legend->AddEntry(hMCmuonsAndPionsDummy, "#mu + #pi", "Lp");
    
    legend->AddEntry(hGenDeltaPiForPiProj, "#pi", "Lp");
    if (plotIdentifiedSpectra)
      legend->AddEntry(hDeltaPiMC[slice][kPi - 1], "#pi", "Lp");
    
    legend->AddEntry(hGenDeltaPiForElProj, "e", "Lp");
    if (plotIdentifiedSpectra) 
      legend->AddEntry(hDeltaPiMC[slice][kEl -1], "e", "Lp");
    
    // TOF selection rather clean, only small electron contamination (muons not taken into account);
    // Fix fraction of K,p,mu to zero (muons can anyway only be treated as an afterburner in case of TOF patching, fixing them
    // as function of pT,z,xi does not work because of unknown TOF efficiency
    Double_t fractionPions = 0.99;
    Double_t fractionErrorUpPions = 1.;
    Double_t fractionErrorLowPions = 0.;
    
    Double_t fractionKaons = 0.;
    Double_t fractionErrorUpKaons = 0.;
    Double_t fractionErrorLowKaons = 0.;
    
    Double_t fractionProtons = 0.;
    Double_t fractionErrorUpProtons = 0.;
    Double_t fractionErrorLowProtons = 0.;
    
    Double_t fractionElectrons = 0.01;
    Double_t fractionErrorUpElectrons = 1.;
    Double_t fractionErrorLowElectrons = 0.;
    
    Double_t fractionMuons = 0.;
    Double_t fractionErrorUpMuons = 0.;
    Double_t fractionErrorLowMuons = 0.;
    
    
    // If the templates are empty, they cannot be used for the fit (fractions do not change the chi^2 (except for regularisation)).
    // Fix fractions to zero in that case (usually only happens for protons (or kaons) at very low pT, where pre-PID works rather
    // perfectly) and set step size to zero.
    if (!nonEmptyTemplateEl) {
      fractionElectrons = 0.;
      fractionErrorUpElectrons = 0.;
      fractionErrorLowElectrons = 0.;
      
      printf("\nEl template empty in this bin - fixing fractions to zero!\n\n");
    }
    if (!nonEmptyTemplatePi) {
      fractionPions = 0.;
      fractionErrorUpPions = 0.;
      fractionErrorLowPions = 0.;
      
      printf("\nPi template empty in this bin - fixing fractions to zero!\n\n");
    }
    
    Double_t gausParamsSimultaneousFit[nParSimultaneousFit] = { 
      fractionPions,
      fractionKaons,
      fractionProtons,
      fractionElectrons,
      fractionMuons,
      totalYield
      // No shifts because they do not make too much sense (different eta + possible deviations from Bethe-Bloch in one x-Bin)
    };
    
    Double_t lowParLimitsSimultaneousFit[nParSimultaneousFit] = {
      fractionErrorLowPions,
      fractionErrorLowKaons,
      fractionErrorLowProtons,
      fractionErrorLowElectrons,
      fractionErrorLowMuons,
      totalYield
    };
    
    Double_t upParLimitsSimultaneousFit[nParSimultaneousFit] = {
      fractionErrorUpPions,
      fractionErrorUpKaons,
      fractionErrorUpProtons,
      fractionErrorUpElectrons,
      fractionErrorUpMuons,
      totalYield
    };
    
    Double_t stepSizeSimultaneousFit[nParSimultaneousFit] = {
      nonEmptyTemplatePi ? 0.1 : 0.,
      0.0,
      0.0,
      nonEmptyTemplateEl ? 0.1 : 0.,
      0.0,
      
      0.0
    };
    
        
    const TString binInfo = (mode == kPMpT) ? Form("%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                            : Form("%.2f_%s_%.2f", pTaxis->GetBinLowEdge(slice + 1), 
                                                    modeShortName[mode].Data(), pTaxis->GetBinUpEdge(slice + 1));
    
    const TString binInfoTitle = (mode == kPMpT) ? Form("%.2f < Pt <%.2f", binsPt[slice], binsPt[slice + 1])
                                                  : Form("%.2f < %s < %.2f", pTaxis->GetBinLowEdge(slice + 1), 
                                                        modeShortName[mode].Data(), 
                                                        pTaxis->GetBinUpEdge(slice + 1));
    
    const TString fitFuncSuffix = (mode == kPMpT) ? Form("%.3f_Pt_%.3f", binsPt[slice], binsPt[slice + 1])
                                                  : Form("%.3f_%s_%.3f", pTaxis->GetBinLowEdge(slice + 1), 
                                                          modeShortName[mode].Data(), 
                                                          pTaxis->GetBinUpEdge(slice + 1));
    
    TCanvas* cSingleFit[4] = { 0x0, };
    for (Int_t species = 0; species < 4; species++) {
      cSingleFit[species] = new TCanvas(Form("cSingleFitTOFpi_%s_%s", binInfo.Data(), speciesLabel[species].Data()), 
                                              Form("single fit (TOF #pi) for %s (%s)", binInfoTitle.Data(), speciesLabel[species].Data()),
                                              1366, 768);
      cSingleFit[species]->Divide(1, 2, 0.01, 0.);
      cSingleFit[species]->GetPad(1)->SetRightMargin(0.001);
      cSingleFit[species]->GetPad(2)->SetRightMargin(0.001);
      cSingleFit[species]->GetPad(1)->SetTopMargin(0.001);
      cSingleFit[species]->GetPad(2)->SetTopMargin(0.01);
      cSingleFit[species]->GetPad(1)->SetBottomMargin(0.01);
      
      cSingleFit[species]->GetPad(1)->SetGridx(kTRUE);
      cSingleFit[species]->GetPad(2)->SetGridx(kTRUE);
      cSingleFit[species]->GetPad(1)->SetGridy(kTRUE);
      cSingleFit[species]->GetPad(2)->SetGridy(kTRUE);
      
      cSingleFit[species]->GetPad(1)->SetLogy(kTRUE);
      cSingleFit[species]->GetPad(1)->SetLogx(kTRUE);
      cSingleFit[species]->GetPad(2)->SetLogx(kTRUE);
    }
      
    Int_t errFlag = 0;
    
    // Reset temp arrays for next slice
    for (Int_t ind = 0; ind < nParSimultaneousFit; ind++)
      parameterErrorsOut[ind] = 0;
    
    std::cout << "Fitting data simultaneously...." << std::endl << std::endl;
    
    // Add ref histos in initialisation step (w/ reg) or in the only loop (w/o reg)
    mathFit->ClearRefHistos();
    
    mathFit->AddRefHisto(hGenDeltaPiForPiProj);
    mathFit->AddRefHisto(hGenDeltaPiForKaProj);
    mathFit->AddRefHisto(hGenDeltaPiForPrProj);
    mathFit->AddRefHisto(hGenDeltaPiForElProj);
    if (takeIntoAccountMuons)
      mathFit->AddRefHisto(hGenDeltaPiForMuProj);
    
    mathFit->AddRefHisto(hGenDeltaKaForPiProj);
    mathFit->AddRefHisto(hGenDeltaKaForKaProj);
    mathFit->AddRefHisto(hGenDeltaKaForPrProj);
    mathFit->AddRefHisto(hGenDeltaKaForElProj);
    if (takeIntoAccountMuons)
      mathFit->AddRefHisto(hGenDeltaKaForMuProj);
    
    mathFit->AddRefHisto(hGenDeltaPrForPiProj);
    mathFit->AddRefHisto(hGenDeltaPrForKaProj);
    mathFit->AddRefHisto(hGenDeltaPrForPrProj);
    mathFit->AddRefHisto(hGenDeltaPrForElProj);
    if (takeIntoAccountMuons)
      mathFit->AddRefHisto(hGenDeltaPrForMuProj);
    
    mathFit->AddRefHisto(hGenDeltaElForPiProj);
    mathFit->AddRefHisto(hGenDeltaElForKaProj);
    mathFit->AddRefHisto(hGenDeltaElForPrProj);
    mathFit->AddRefHisto(hGenDeltaElForElProj);
    if (takeIntoAccountMuons)
      mathFit->AddRefHisto(hGenDeltaElForMuProj);

    TH1D* hDeltaSpecies[numSimultaneousFits] = { hDeltaPi[slice], hDeltaKa[slice], hDeltaPr[slice], hDeltaEl[slice] };
    Double_t reducedChiSquare = -1;
    
    errFlag = errFlag | 
              doSimultaneousFit(hDeltaSpecies, xLow, xUp, nParSimultaneousFit, gausParamsSimultaneousFit, parameterErrorsOut, 
                                0x0, stepSizeSimultaneousFit, lowParLimitsSimultaneousFit, upParLimitsSimultaneousFit, reducedChiSquare, 0, 0);
          
    // Forward parameters to single fits
    for (Int_t parIndex = 0; parIndex < nParSimultaneousFit; parIndex++) {
      // Fractions
      if (parIndex <= 4) {
        totalDeltaPion->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
        totalDeltaPion->SetParError(parIndex, parameterErrorsOut[parIndex]);
        
        totalDeltaKaon->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
        totalDeltaKaon->SetParError(parIndex, parameterErrorsOut[parIndex]);
        
        totalDeltaProton->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
        totalDeltaProton->SetParError(parIndex, parameterErrorsOut[parIndex]);
        
        totalDeltaElectron->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
        totalDeltaElectron->SetParError(parIndex, parameterErrorsOut[parIndex]);
      }
      // Total yield
      else if (parIndex == 5) {
        totalDeltaPion->SetParameter(parIndex, totalYield);
        totalDeltaPion->SetParError(parIndex, 0);
        
        totalDeltaKaon->SetParameter(parIndex, totalYield);
        totalDeltaKaon->SetParError(parIndex, 0);
        
        totalDeltaProton->SetParameter(parIndex, totalYield);
        totalDeltaProton->SetParError(parIndex, 0);
        
        totalDeltaElectron->SetParameter(parIndex, totalYield);
        totalDeltaElectron->SetParError(parIndex, 0);
      }
      // Hist shifts
      else {
        totalDeltaPion->SetParameter(parIndex, 0);
        totalDeltaPion->SetParError(parIndex, 0);
        
        totalDeltaKaon->SetParameter(parIndex, 0);
        totalDeltaKaon->SetParError(parIndex, 0);
        
        totalDeltaProton->SetParameter(parIndex, 0);
        totalDeltaProton->SetParError(parIndex, 0);
        
        totalDeltaElectron->SetParameter(parIndex, 0);
        totalDeltaElectron->SetParError(parIndex, 0);
      }
    }

    // Plot single fits
    
    Int_t binLow = -1;
    Int_t binHigh = -1;
    
    // DeltaPions
    cSingleFit[2]->cd(1);
    
    hDeltaPi[slice]->SetTitle("");
    hDeltaPi[slice]->GetYaxis()->SetTitle("Entries");
    SetReasonableXaxisRange(hDeltaPi[slice], binLow, binHigh);
    hDeltaPi[slice]->Draw("e");
    
    TF1* fitFuncTotalDeltaPion = (TF1*)totalDeltaPion->Clone(Form("Fit_Total_TOFpi_DeltaPion_%s", fitFuncSuffix.Data()));
    
    hDeltaPiFitQA[slice] = GetFitQAhisto(hDeltaPi[slice], fitFuncTotalDeltaPion, Form("hDeltaPiFitQA_TOFpi_%d", slice));
    
    hDeltaPi[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaPion);
    fitFuncTotalDeltaPion->Draw("same");   
    
    Double_t* parametersOut = &totalDeltaPion->GetParameters()[0];
    
    hGenDeltaPiForPiProj->Scale(parametersOut[5] * parametersOut[0]);
    hGenDeltaPiForPiProj->Draw("same");
    
    hGenDeltaPiForKaProj->Scale(parametersOut[5] * parametersOut[1]);
    hGenDeltaPiForKaProj->Draw("same");
    
    hGenDeltaPiForPrProj->Scale(parametersOut[5] * parametersOut[2]);
    hGenDeltaPiForPrProj->Draw("same");
    
    hGenDeltaPiForElProj->Scale(parametersOut[5] * parametersOut[3]);
    hGenDeltaPiForElProj->Draw("same");
    
    if (takeIntoAccountMuons) {
      hGenDeltaPiForMuProj->Scale(parametersOut[5] * parametersOut[4]);
      if (muonFractionHandling != kNoMuons)
        hGenDeltaPiForMuProj->Draw("same");
    }
    
    if (plotIdentifiedSpectra) {
      for (Int_t species = 0; species < 5; species++) 
        hDeltaPiMC[slice][species]->Draw("same");
      
      // Draw histo for sum of MC muons and pions
      TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPiMC[slice][kPi - 1]);
      hMCmuonsAndPions->Add(hDeltaPiMC[slice][kMu - 1]);
      hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPiMC[slice][kPi - 1]->GetName()));
      hMCmuonsAndPions->Draw("same");
    }
    
    hDeltaPi[slice]->Draw("esame");
    
    legend->Draw();
    
    cSingleFit[2]->cd(2);
    hDeltaPiFitQA[slice]->Draw("e");    
    
    // DeltaElectrons
    cSingleFit[0]->cd(1);
    
    hDeltaEl[slice]->SetTitle("");
    hDeltaEl[slice]->GetYaxis()->SetTitle("Entries");
    SetReasonableXaxisRange(hDeltaEl[slice], binLow, binHigh);
    hDeltaEl[slice]->Draw("e");
    
    TF1* fitFuncTotalDeltaElectron = (TF1*)totalDeltaElectron->Clone(Form("Fit_Total_TOFpi_DeltaElectron_%s", fitFuncSuffix.Data()));
    
    hDeltaElFitQA[slice] = GetFitQAhisto(hDeltaEl[slice], fitFuncTotalDeltaElectron, Form("hDeltaElFitQA_TOFpi_%d", slice));
    
    hDeltaEl[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaElectron);
    fitFuncTotalDeltaElectron->Draw("same");  
    
    parametersOut = &totalDeltaElectron->GetParameters()[0];
    
    hGenDeltaElForPiProj->Scale(parametersOut[5] * parametersOut[0]);
    hGenDeltaElForPiProj->Draw("same");
    
    hGenDeltaElForKaProj->Scale(parametersOut[5] * parametersOut[1]);
    hGenDeltaElForKaProj->Draw("same");
    
    hGenDeltaElForPrProj->Scale(parametersOut[5] * parametersOut[2]);
    hGenDeltaElForPrProj->Draw("same");
    
    hGenDeltaElForElProj->Scale(parametersOut[5] * parametersOut[3]);
    hGenDeltaElForElProj->Draw("same");
    
    if (takeIntoAccountMuons) {
      hGenDeltaElForMuProj->Scale(parametersOut[5] * parametersOut[4]);
      if (muonFractionHandling != kNoMuons)
        hGenDeltaElForMuProj->Draw("same");
    }
    
    if (plotIdentifiedSpectra) {
      for (Int_t species = 0; species < 5; species++) 
        hDeltaElMC[slice][species]->Draw("same");
      
      // Draw histo for sum of MC muons and pions
      TH1D* hMCmuonsAndPions = new TH1D(*hDeltaElMC[slice][kPi - 1]);
      hMCmuonsAndPions->Add(hDeltaElMC[slice][kMu - 1]);
      hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaElMC[slice][kPi - 1]->GetName()));
      hMCmuonsAndPions->Draw("same");
    }
    
    hDeltaEl[slice]->Draw("esame");
    
    legend->Draw();
    
    cSingleFit[0]->cd(2);
    hDeltaElFitQA[slice]->Draw("e");
    
    // DeltaKaons 
    cSingleFit[1]->cd(1);
    
    hDeltaKa[slice]->SetTitle("");
    hDeltaKa[slice]->GetYaxis()->SetTitle("Entries");
    SetReasonableXaxisRange(hDeltaKa[slice], binLow, binHigh);
    hDeltaKa[slice]->Draw("e");
    
    TF1* fitFuncTotalDeltaKaon = (TF1*)totalDeltaKaon->Clone(Form("Fit_Total_TOFpi_DeltaKaon_%s", fitFuncSuffix.Data()));
    
    hDeltaKaFitQA[slice] = GetFitQAhisto(hDeltaKa[slice], fitFuncTotalDeltaKaon, Form("hDeltaKaFitQA_TOFpi_%d", slice));
    
    hDeltaKa[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaKaon);
    fitFuncTotalDeltaKaon->Draw("same");  
    
    parametersOut = &totalDeltaKaon->GetParameters()[0];
    
    hGenDeltaKaForPiProj->Scale(parametersOut[5] * parametersOut[0]);
    hGenDeltaKaForPiProj->Draw("same");
    
    hGenDeltaKaForKaProj->Scale(parametersOut[5] * parametersOut[1]);
    hGenDeltaKaForKaProj->Draw("same");
    
    hGenDeltaKaForPrProj->Scale(parametersOut[5] * parametersOut[2]);
    hGenDeltaKaForPrProj->Draw("same");
    
    hGenDeltaKaForElProj->Scale(parametersOut[5] * parametersOut[3]);
    hGenDeltaKaForElProj->Draw("same");
    
    if (takeIntoAccountMuons) {
      hGenDeltaKaForMuProj->Scale(parametersOut[5] * parametersOut[4]);
      if (muonFractionHandling != kNoMuons)
        hGenDeltaKaForMuProj->Draw("same");
    }
    
    if (plotIdentifiedSpectra) {
      for (Int_t species = 0; species < 5; species++) 
        hDeltaKaMC[slice][species]->Draw("same");
      
      // Draw histo for sum of MC muons and pions
      TH1D* hMCmuonsAndPions = new TH1D(*hDeltaKaMC[slice][kPi - 1]);
      hMCmuonsAndPions->Add(hDeltaKaMC[slice][kMu - 1]);
      hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaKaMC[slice][kPi - 1]->GetName()));
      hMCmuonsAndPions->Draw("same");
    }
    
    hDeltaKa[slice]->Draw("esame");
    
    legend->Draw();
    
    cSingleFit[1]->cd(2);
    hDeltaKaFitQA[slice]->Draw("e");
    
    // DeltaProtons
    cSingleFit[3]->cd(1);
    
    hDeltaPr[slice]->SetTitle("");
    hDeltaPr[slice]->GetYaxis()->SetTitle("Entries");
    SetReasonableXaxisRange(hDeltaPr[slice], binLow, binHigh);
    hDeltaPr[slice]->Draw("e");
    
    TF1* fitFuncTotalDeltaProton = (TF1*)totalDeltaProton->Clone(Form("Fit_Total_TOFpi_DeltaProton_%s", fitFuncSuffix.Data()));
    
    hDeltaPrFitQA[slice] = GetFitQAhisto(hDeltaPr[slice], fitFuncTotalDeltaProton, Form("hDeltaPrFitQA_TOFpi_%d", slice));
    
    hDeltaPr[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaProton);
    
    fitFuncTotalDeltaProton->Draw("same");  
    
    parametersOut = &totalDeltaProton->GetParameters()[0];
    
    hGenDeltaPrForPiProj->Scale(parametersOut[5] * parametersOut[0]);
    hGenDeltaPrForPiProj->Draw("same");
    
    hGenDeltaPrForKaProj->Scale(parametersOut[5] * parametersOut[1]);
    hGenDeltaPrForKaProj->Draw("same");
    
    hGenDeltaPrForPrProj->Scale(parametersOut[5] * parametersOut[2]);
    hGenDeltaPrForPrProj->Draw("same");
    
    hGenDeltaPrForElProj->Scale(parametersOut[5] * parametersOut[3]);
    hGenDeltaPrForElProj->Draw("same");
    
    if (takeIntoAccountMuons) {
      hGenDeltaPrForMuProj->Scale(parametersOut[5] * parametersOut[4]);
      if (muonFractionHandling != kNoMuons)
        hGenDeltaPrForMuProj->Draw("same");
    }
    
    if (plotIdentifiedSpectra) {
      for (Int_t species = 0; species < 5; species++) 
        hDeltaPrMC[slice][species]->Draw("same");
      
      // Draw histo for sum of MC muons and pions
      TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPrMC[slice][kPi - 1]);
      hMCmuonsAndPions->Add(hDeltaPrMC[slice][kMu - 1]);
      hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
      hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPrMC[slice][kPi - 1]->GetName()));
      hMCmuonsAndPions->Draw("same");
    }
    
    hDeltaPr[slice]->Draw("esame");
    
    legend->Draw();
    
    cSingleFit[3]->cd(2);
    hDeltaPrFitQA[slice]->Draw("e");
    
    for (Int_t species = 0; species < 4; species++) {
      cSingleFit[species]->Modified();
      cSingleFit[species]->Update();
    }
    
    std::cout << std::endl << std::endl;
    
    // Update TOF yields. Nothing like "inverseBinWidth" used since this is also not taken into account in the original histograms.
    // It will be handled at another place, thus.
    
    // Parameters are the same for all anyway (simultaneous fit), just take those from pions.
    parametersOut = &totalDeltaPion->GetParameters()[0];
    
    // Force sum of fractions to be unity
    Double_t sumFractions = parametersOut[0] + parametersOut[3];
    // If something went wrong (sumFractions <= 0) (might happen for extremely low statistics?), just leave everything for the pions
    // and don't touch the histos.
    if (sumFractions > 0) {
      Double_t elFraction = parametersOut[3] / sumFractions;
      Double_t piFraction = parametersOut[0] / sumFractions;
      
      for (Int_t binX = pBinLowProjLimit; binX <= pBinUpProjLimit; binX++) {
        //NOTE: Since the pion yields are updated, but need to be used for the electron yields, they must be saved first!
        const Double_t piYield = hYieldTOFOrigBinningPions->GetBinContent(binX);
        const Double_t piYieldError = hYieldTOFOrigBinningPions->GetBinError(binX);
        
        // In each bin: Spit statistics between pi and el; assume that fit has no error, i.e. the error
        // is just scaled with sqrt(fraction), such that summing up both histos gives back the orriginal error
        hYieldTOFOrigBinningElectrons->SetBinContent(binX, piYield *  elFraction);
        hYieldTOFOrigBinningElectrons->SetBinError(binX, piYieldError * TMath::Sqrt(elFraction));
        
        hYieldTOFOrigBinningPions->SetBinContent(binX, piYield *  piFraction);
        hYieldTOFOrigBinningPions->SetBinError(binX, piYieldError * TMath::Sqrt(piFraction));
      }
    }
    
    // Save results
    
    TString saveDir = (mode == kPMpT) ? Form("TOFpiFits/SingleFit_%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                      : Form("TOFpiFits/SingleFit_%.2f_%s_%.2f", pTaxis->GetBinLowEdge(slice + 1), 
                                            modeShortName[mode].Data(), pTaxis->GetBinUpEdge(slice + 1));
    saveF->mkdir(saveDir.Data());
    saveF->cd(saveDir.Data());
    
    for (Int_t species = 0; species < 4; species++) {
      if (cSingleFit[species]) {
        cSingleFit[species]->Write();
        delete cSingleFit[species];
      }
    }
    
    if (hDeltaPi[slice])
      hDeltaPi[slice]->Write();
    
    if (hDeltaEl[slice])
      hDeltaEl[slice]->Write();
    
    if (hDeltaKa[slice])
      hDeltaKa[slice]->Write();
    
    if (hDeltaPr[slice])
      hDeltaPr[slice]->Write();
    
    
    if (hDeltaPiFitQA[slice])
      hDeltaPiFitQA[slice]->Write();
    delete hDeltaPiFitQA[slice];
    
    if (hDeltaElFitQA[slice])
      hDeltaElFitQA[slice]->Write();
    delete hDeltaElFitQA[slice];
    
    if (hDeltaKaFitQA[slice])
      hDeltaKaFitQA[slice]->Write();
    delete hDeltaKaFitQA[slice];
    
    if (hDeltaPrFitQA[slice])
      hDeltaPrFitQA[slice]->Write();
    delete hDeltaPrFitQA[slice];
    
    if (hGenDeltaElForElProj) 
      hGenDeltaElForElProj->Write();
    delete hGenDeltaElForElProj;
    
    if (hGenDeltaElForKaProj) 
      hGenDeltaElForKaProj->Write();
    delete hGenDeltaElForKaProj;
    
    if (hGenDeltaElForPiProj) 
      hGenDeltaElForPiProj->Write();
    delete hGenDeltaElForPiProj;
    
    if (hGenDeltaElForPrProj) 
      hGenDeltaElForPrProj->Write();
    delete hGenDeltaElForPrProj;
    
    if (hGenDeltaElForMuProj) 
      hGenDeltaElForMuProj->Write();
    delete hGenDeltaElForMuProj;
    
    delete fitFuncTotalDeltaElectron;
    
    if (hGenDeltaKaForElProj) 
      hGenDeltaKaForElProj->Write();
    delete hGenDeltaKaForElProj;
    
    if (hGenDeltaKaForKaProj) 
      hGenDeltaKaForKaProj->Write();
    delete hGenDeltaKaForKaProj;
    
    if (hGenDeltaKaForPiProj) 
      hGenDeltaKaForPiProj->Write();
    delete hGenDeltaKaForPiProj;
    
    if (hGenDeltaKaForPrProj) 
      hGenDeltaKaForPrProj->Write();
    delete hGenDeltaKaForPrProj;
    
    if (hGenDeltaKaForMuProj) 
      hGenDeltaKaForMuProj->Write();
    delete hGenDeltaKaForMuProj;
    
    delete fitFuncTotalDeltaKaon;
    

    if (hGenDeltaPiForElProj) 
      hGenDeltaPiForElProj->Write();
    delete hGenDeltaPiForElProj;
    
    if (hGenDeltaPiForKaProj) 
      hGenDeltaPiForKaProj->Write();
    delete hGenDeltaPiForKaProj;
    
    if (hGenDeltaPiForPiProj) 
      hGenDeltaPiForPiProj->Write();
    delete hGenDeltaPiForPiProj;
    
    if (hGenDeltaPiForPrProj) 
      hGenDeltaPiForPrProj->Write();
    delete hGenDeltaPiForPrProj;
    
    if (hGenDeltaPiForMuProj) 
      hGenDeltaPiForMuProj->Write();
    delete hGenDeltaPiForMuProj;
    
    delete fitFuncTotalDeltaPion;
    
    
    if (hGenDeltaPrForElProj) 
      hGenDeltaPrForElProj->Write();
    delete hGenDeltaPrForElProj;
    
    if (hGenDeltaPrForKaProj) 
      hGenDeltaPrForKaProj->Write();
    delete hGenDeltaPrForKaProj;
    
    if (hGenDeltaPrForPiProj) 
      hGenDeltaPrForPiProj->Write();
    delete hGenDeltaPrForPiProj;
    
    if (hGenDeltaPrForPrProj) 
      hGenDeltaPrForPrProj->Write();
    delete hGenDeltaPrForPrProj;
    
    if (hGenDeltaPrForMuProj) 
      hGenDeltaPrForMuProj->Write();
    delete hGenDeltaPrForMuProj;
    
    delete fitFuncTotalDeltaProton;
    
    delete totalDeltaElectron;
    delete totalDeltaKaon;
    delete totalDeltaPion;
    delete totalDeltaProton;
    
    delete legend;
    
    if (errFlag != 0)
      std::cout << "errFlag " << errFlag << std::endl << std::endl;
  }
  
  // Clean up - need new mathFit instance for main fit
  delete mathFit;
  
  std::cout << std::endl << std::endl << "Done!" << std::endl << std::endl;
}


//____________________________________________________________________________________________________________________
Int_t PID(TString fileName, Double_t deta, Double_t pLow, Double_t pHigh, Bool_t isMCdataSet, Int_t fitMethod, 
          Int_t muonFractionHandlingParameter, //0 = no muons, 1 = muonFrac=elFrac,
                                               //2(3) = muonFrac/elFrac tuned on MC for StandardTrackCuts(HybridTrackCuts)
          Bool_t useIdentifiedGeneratedSpectra, Bool_t plotIdentifiedSpectra, Int_t mode/*0=pT,1=z,2=xi,3=r,4=jT*/,
          Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
          Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
          Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/,
          Int_t rebin/* = 1 -> DON'T USE FOR PT (will not work since binsPt will and should be used!)*/,
          Int_t rebinDeltaPrime/* = 1*/,
          TString listName /* = "bhess_PID"*/,
          Bool_t useLogLikelihood /*= kTRUE*/, Bool_t useWeightsForLogLikelihood /*= kFALSE*/,
          Int_t regularisation /*= 0*/,
          Double_t regularisationFactor /*= 1*/,
          Bool_t applyTOFpatching,
          TString filePathNameFileWithInititalFractions /*= ""*/,
          TString* filePathNameResults /*= 0x0*/,
          Int_t binTypePt /*=kPtBinTypeJets*/,
          Double_t yieldThresholdForFitting /*= 0*/)
{
  // Do all the fitting
  
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  ResetGlobals();
  
  // Set the pT binning
  binsPt = GetPtBins(binTypePt, nPtBins);
  
  isMC = isMCdataSet;
  
  muonFractionHandling = muonFractionHandlingParameter;
  processingMode = mode;
  
  //NOTE/TODO: For mode != kPMpT, there is no well-defined pT per bin. So, one cannot use the usual function to extract
  // the muons from the electrons. This might be done, if one extracts the function vs. z,xi....
  // Therefore, just set the fraction to "noMuons" for the moment
  if (mode != kPMpT && muonFractionHandlingParameter != kNoMuons) {
    printf("WARNING: For mode != kPMpT the muon fraction handling is not yet properly implemented. Setting treatment to \"no muon\"!\n");
    muonFractionHandling = muonFractionHandlingParameter = kNoMuons;
  }
  
  Int_t genAxis = useDeltaPrime ? kPidGenDeltaPrime : 1000/*kPidGenDelta*/;
  if (!useDeltaPrime) {
    std::cout << "ERROR: delta plots no longer available!" << std::endl;
    return -1;
  }
  
  if (listName == "") {
    listName = fileName;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  
  if (rebin > 1 && mode == kPMpT) {
    std::cout << "ERROR: Requested re-binning of pT-axis! Since binsPt will be used, re-binning the data histo will lead to "
              << "unforeseen consequences!" << std::endl;
    return -1;
  }
    
  Int_t pSliceLow = -1;
  Int_t pSliceHigh = -1;
  
  Int_t axisForMode = kPidPt;
  Int_t axisGenForMode = kPidGenPt;
  Int_t axisGenYieldForMode = kPidGenYieldPt;
  
  std::cout << "Fitting \"" << fileName.Data() << "\" with settings:" << std::endl;
  
  std::cout << "Minimisation strategy: " << minimisationStrategy.Data() << std::endl;
  if (useLogLikelihood) 
    std::cout << "Binned loglikelihood fit" << (useWeightsForLogLikelihood ? " (weighted)" : "")  << std::endl;
  else
    std::cout << "ChiSquare fit" << std::endl;
  std::cout << "Processing mode: ";
  if (mode == kPMpT)
    std::cout << "pT" << ", binning scheme " << binTypePt << std::endl;
  else if (mode == kPMz) {
    std::cout << "z" << std::endl;
    axisForMode = kPidZ;
    axisGenForMode = kPidGenZ;
    axisGenYieldForMode = kPidGenYieldZ;
  }
  else if (mode == kPMxi) {
    std::cout << "xi" << std::endl;
    axisForMode = kPidXi;
    axisGenForMode = kPidGenXi;
    axisGenYieldForMode = kPidGenYieldXi;
  }
  else if (mode == kPMdistance) {
    std::cout << "r" << std::endl;
    axisForMode = kPidDistance;
    axisGenForMode = kPidGenDistance;
    axisGenYieldForMode = kPidGenYieldDistance;
  }
  else if (mode == kPMjT) {
    std::cout << "jT" << std::endl;
    axisForMode = kPidJt;
    axisGenForMode = kPidGenJt;
    axisGenYieldForMode = kPidGenYieldJt;
  }
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  std::cout << "Charge selection: ";
  if (chargeMode == kAllCharged)
    std::cout << "All charged particles" << std::endl;
  else if (chargeMode == kNegCharge)
    std::cout << "Negative particles only" << std::endl;
  else if (chargeMode == kPosCharge)
    std::cout << "Positive particles only" << std::endl;
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  if (regularisation > 0)
    std::cout << "Regularisation with +-" << regularisation << " bins and factor " << regularisationFactor << " for penalty term."
              << std::endl;
  else
    std::cout << "No regularisation" << std::endl;
  
  std::cout << "Assumption on muon fraction: ";
  if (muonFractionHandlingParameter >= 0 && muonFractionHandlingParameter < kNumHandlings)
    std::cout << muonFractionHandlingShortName[muonFractionHandlingParameter].Data() << std::endl;
  /*if (muonFractionHandlingParameter == kNoMuons)
    std::cout << "Identical zero" << std::endl;
  else if (muonFractionHandlingParameter == kMuonFracEqualElFrac)
    std::cout << "Equal electron fraction" << std::endl;
  else if (muonFractionHandlingParameter == kMuonFracOverElFracTunedOnMCStandardTrackCuts)
    std::cout << "Ratio to electron fraction tuned on MC for standard track cuts" << std::endl;
  else if (muonFractionHandlingParameter == kMuonFracOverElFracTunedOnMCHybridTrackCuts)
    std::cout << "Ratio to electron fraction tuned on MC for hybrid track cuts" << std::endl;
  else if (muonFractionHandlingParameter == kMuonFracOverElFracTunedOnMCHybridTrackCutsJets)
    std::cout << "Ratio to electron fraction tuned on MC for hybrid track cuts for jet particles" << std::endl;
  else if (muonFractionHandlingParameter == kMuonFracOverElFracTunedOnMCHybridTrackCutsJets)
    std::cout << "Ratio to electron fraction tuned on MC for hybrid track cuts for jet particles" << std::endl;*/
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  if (mode == kPMpT) {
    Int_t index = 0;
    while (pLow >= binsPt[index] && index < nPtBins)
      index++;
    pSliceLow = index - 1;
    
    index = 0;
    while (pHigh > binsPt[index] && index < nPtBins)
      index++;
    pSliceHigh = index - 1;
    
    Int_t numMomIntervals = pSliceHigh - pSliceLow + 1;
    
    if (numMomIntervals <= 0 || pSliceLow < 0 || pSliceHigh > nPtBins)  {
      std::cout << "Wrong choice of limits pLow/pHigh!" << std::endl;
      return -1;
    }
      
    pLow = binsPt[pSliceLow];
    pHigh =  binsPt[pSliceHigh + 1]; // need upper edge, but binsPt holds lower edge
    std::cout << "pLow/pHigh: ";
    std::cout << pLow << " / " << pHigh << std::endl;
  }
  
  std::cout << "TOF patching: " << applyTOFpatching << std::endl;
  
  gYieldThresholdForFitting = yieldThresholdForFitting;
  std::cout << "(Total) yield threshold for fitting: " << gYieldThresholdForFitting << std::endl;
  
  Bool_t initialiseWithFractionsFromFile = kFALSE;
  TFile* fInitialFractions = 0x0;
  TH1 *hInitFracEl = 0x0, *hInitFracKa = 0x0, *hInitFracPi = 0x0, *hInitFracMu = 0x0, *hInitFracPr = 0x0;
  
  if (filePathNameFileWithInititalFractions != "") {
    initialiseWithFractionsFromFile = kTRUE;
    
    std::cout << "Initialising fractions from file: " << filePathNameFileWithInititalFractions.Data() << std::endl;
  }
  else
    std::cout << "Not initialising fractions from file" << std::endl;

  if (initialiseWithFractionsFromFile) {
    fInitialFractions = TFile::Open(filePathNameFileWithInititalFractions.Data());
    if (!fInitialFractions) {
      std::cout << std::endl;
      std::cout << "Failed to open file with initial fractions \"" << filePathNameFileWithInititalFractions.Data() << "\"!"
                << std::endl;
      return -1;
    }
    
    hInitFracEl = (TH1*)fInitialFractions->Get("hFractionElectrons");
    hInitFracKa = (TH1*)fInitialFractions->Get("hFractionKaons");
    hInitFracPi = (TH1*)fInitialFractions->Get("hFractionPions");
    hInitFracMu = (TH1*)fInitialFractions->Get("hFractionMuons");
    hInitFracPr = (TH1*)fInitialFractions->Get("hFractionProtons");
    
    if (!hInitFracEl || ! hInitFracKa || ! hInitFracPi  || ! hInitFracMu  || ! hInitFracPr) {
      std::cout << std::endl;
      std::cout << "Failed to load initial fractions from file \"" << filePathNameFileWithInititalFractions.Data() << "\"!"
                << std::endl;
      
      fInitialFractions->Close();
      return -1;
    }
  }
  
  
  
  TObjArray* histList = 0x0;
  
  TFile* f = TFile::Open(fileName.Data());
  if (!f)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << fileName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  //TString listName = fileName;
  //listName = listName.ReplaceAll(".root", "");
  //listName = listName.Remove(1, listName.Last('/') + 1);
  histList = (TObjArray*)(f->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data histogram
  THnSparse* hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
  if (!hPIDdata) {
    std::cout << std::endl;
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  // If desired, rebin considered axis
  if (rebin > 1 || rebinDeltaPrime > 1) {
    const Int_t nDimensions = hPIDdata->GetNdimensions();
    Int_t rebinFactor[nDimensions];
    
    for (Int_t dim = 0; dim < nDimensions; dim++) {
      if (dim == axisForMode && rebin > 1)
        rebinFactor[dim] = rebin;
      else if (dim == kPidDeltaPrime && rebinDeltaPrime > 1)
        rebinFactor[dim] = rebinDeltaPrime;
      else
        rebinFactor[dim] = 1;
    }
    
    THnSparse* temp = hPIDdata->Rebin(&rebinFactor[0]);
    hPIDdata->Reset();
    hPIDdata = temp;
  }
  
  // Set proper errors, if not yet calculated
  if (!hPIDdata->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << hPIDdata->GetName() << "..." << std::endl;
    hPIDdata->Sumw2();
    Long64_t nBinsTHnSparse = hPIDdata->GetNbins();
    Double_t binContent = 0;
    
    for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
      binContent = hPIDdata->GetBinContent(bin);
      hPIDdata->SetBinError(bin, TMath::Sqrt(binContent));
    }
  }
  
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2;
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindFixBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindFixBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 0
        && upperCentralityBinLimit <= hPIDdata->GetAxis(kPidCentrality)->GetNbins() + 1) {
      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  if (!restrictCentralityAxis)
    GetAxisRangeForMultiplicityAxisForMB(hPIDdata->GetAxis(kPidCentrality), lowerCentralityBinLimit, upperCentralityBinLimit);
  
  hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  
  actualLowerCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinLowEdge(hPIDdata->GetAxis(kPidCentrality)->GetFirst());
  actualUpperCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinUpEdge(hPIDdata->GetAxis(kPidCentrality)->GetLast());
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis)
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
  else
    std::cout << "MB (" << actualLowerCentrality << " - " << actualUpperCentrality << ")" << std::endl;
  
  const Bool_t centralityHasDecimalsPlaces = CentralityHasDecimalsPlaces(actualLowerCentrality) ||
                                             CentralityHasDecimalsPlaces(actualUpperCentrality);
  
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -2;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindFixBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindFixBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 0
        && upperJetPtBinLimit <= hPIDdata->GetAxis(kPidJetPt)->GetNbins() + 1) {
      actualLowerJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "jet pT: ";
  if (restrictJetPtAxis) {
    std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  if (restrictJetPtAxis) {
    hPIDdata->GetAxis(kPidJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
  }
  
  
  // If desired, restrict charge axis
  const Int_t indexChargeAxisData = GetAxisByTitle(hPIDdata, "Charge (e_{0})");
  if (indexChargeAxisData < 0 && restrictCharge) {
    std::cout << "Error: Charge axis not found for data histogram!" << std::endl;
    return -1;
  }
  Int_t lowerChargeBinLimitData = -1;
  Int_t upperChargeBinLimitData = -2;
  Double_t actualLowerChargeData = -999;
  Double_t actualUpperChargeData = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindFixBin(-1. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindFixBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindFixBin(0. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindFixBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimitData <= upperChargeBinLimitData && lowerChargeBinLimitData >= 0
        && upperChargeBinLimitData <= hPIDdata->GetAxis(indexChargeAxisData)->GetNbins() + 1) {
      actualLowerChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinLowEdge(lowerChargeBinLimitData);
      actualUpperChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinUpEdge(upperChargeBinLimitData);
      
      std::cout << "Charge range data: " << actualLowerChargeData << " - " << actualUpperChargeData << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    hPIDdata->GetAxis(indexChargeAxisData)->SetRange(lowerChargeBinLimitData, upperChargeBinLimitData);
  }
  
  std::cout << std::endl;
 
  // Open file in which all the projections (= intermediate results) will be saved
  TString saveInterFName = fileName;
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  saveInterFName = Form("%s_Projections_%s_%d_%s%s%s%s%s%s.root", saveInterFName.ReplaceAll(".root", "").Data(), 
                        modeShortName[mode].Data(),
                        fitMethod, muonFractionHandlingShortName[muonFractionHandlingParameter].Data(),
                        useIdentifiedGeneratedSpectra ? "_idSpectra" : "",
                        restrictCentralityAxis ? Form(centralityHasDecimalsPlaces ? "_centrality%.0fem2_%.0fem2" : "_centrality%.0f_%.0f", 
                                                      centralityHasDecimalsPlaces ? 100.*actualLowerCentrality : actualLowerCentrality,
                                                      centralityHasDecimalsPlaces ? 100.*actualUpperCentrality : actualUpperCentrality) : "",
                        restrictJetPtAxis ? Form("_jetPt%.1f_%.1f", actualLowerJetPt, actualUpperJetPt) : "",
                        chargeString.Data(), applyTOFpatching ? "_TPConly" : "");
  TFile *saveInterF = TFile::Open(saveInterFName.Data(), "RECREATE");
  saveInterF->cd();
  
  
  const Int_t indexTOFinfoAxisData = GetAxisByTitle(hPIDdata, "TOF PID Info");
  if (indexTOFinfoAxisData < 0 && applyTOFpatching) {
    std::cout << "Error: TOF PID info axis not found for data histogram!" << std::endl;
    return -1;
  }
  
  
  TH1D* hYieldTOFOrigBinningPions = 0x0;
  TH1D* hYieldTOFOrigBinningKaons = 0x0;
  TH1D* hYieldTOFOrigBinningProtons = 0x0;
  
  
  TH1F* hYieldTOFPions = 0x0;
  TH1F* hYieldTOFKaons = 0x0;
  TH1F* hYieldTOFProtons = 0x0;
  TH1F* hYieldTOFElectrons = 0x0;
  TH1F* hYieldTOFTotal = 0x0;
  TH1F* hYieldTPConlyTotal = 0x0;
  TH2D* hTOFPurity = 0x0;
  
  TH2D* hMCdataTOF = 0x0;
  
  if (applyTOFpatching) {
    // Extract TOF yields (project to arbitrary selectSpecies to avoid multiple counting)
    hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1);
    Int_t tofBin = -1;
    
    tofBin = hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFpion);
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(tofBin, tofBin);
    hYieldTOFOrigBinningPions = hPIDdata->Projection(axisForMode, "e");
    hYieldTOFOrigBinningPions->SetName("hYieldTOFOrigBinningPions");
    hYieldTOFOrigBinningPions->SetTitle("#pi");
    hYieldTOFOrigBinningPions->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFOrigBinningPions->SetLineColor(getLineColor(kPi));
    hYieldTOFOrigBinningPions->SetMarkerColor(getLineColor(kPi));
    hYieldTOFOrigBinningPions->SetMarkerStyle(22);
    hYieldTOFOrigBinningPions->SetStats(kFALSE);
    
    tofBin = hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFkaon);
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(tofBin, tofBin);
    hYieldTOFOrigBinningKaons = hPIDdata->Projection(axisForMode, "e");
    hYieldTOFOrigBinningKaons->SetName("hYieldTOFOrigBinningKaons");
    hYieldTOFOrigBinningKaons->SetTitle("#K");
    hYieldTOFOrigBinningKaons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFOrigBinningKaons->SetLineColor(getLineColor(kKa));
    hYieldTOFOrigBinningKaons->SetMarkerColor(getLineColor(kKa));
    hYieldTOFOrigBinningKaons->SetMarkerStyle(22);
    hYieldTOFOrigBinningKaons->SetStats(kFALSE);
    
    tofBin = hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFproton);
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(tofBin, tofBin);
    hYieldTOFOrigBinningProtons = hPIDdata->Projection(axisForMode, "e");
    hYieldTOFOrigBinningProtons->SetName("hYieldTOFOrigBinningProtons");
    hYieldTOFOrigBinningProtons->SetTitle("p");
    hYieldTOFOrigBinningProtons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFOrigBinningProtons->SetLineColor(getLineColor(kPr));
    hYieldTOFOrigBinningProtons->SetMarkerColor(getLineColor(kPr));
    hYieldTOFOrigBinningProtons->SetMarkerStyle(22);
    hYieldTOFOrigBinningProtons->SetStats(kFALSE);
    
    // Reset TOF info axis range
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(0, -1);
    
    // Set up histos for the final binning
    hYieldTOFPions = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldTOFPions", "#pi");
    hYieldTOFPions->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFPions->SetLineColor(getLineColor(kPi));
    hYieldTOFPions->SetMarkerColor(getLineColor(kPi));
    hYieldTOFPions->SetMarkerStyle(22);
    hYieldTOFPions->SetStats(kFALSE);
    
    hYieldTOFKaons = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldTOFKaons", "K");
    hYieldTOFKaons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFKaons->SetLineColor(getLineColor(kKa));
    hYieldTOFKaons->SetMarkerColor(getLineColor(kKa));
    hYieldTOFKaons->SetMarkerStyle(22);
    hYieldTOFKaons->SetStats(kFALSE);
    
    hYieldTOFProtons = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldTOFProtons", "p");
    hYieldTOFProtons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFProtons->SetLineColor(getLineColor(kPr));
    hYieldTOFProtons->SetMarkerColor(getLineColor(kPr));
    hYieldTOFProtons->SetMarkerStyle(22);
    hYieldTOFProtons->SetStats(kFALSE);
    
    hYieldTOFElectrons = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldTOFElectrons", "e");
    hYieldTOFElectrons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTOFElectrons->SetLineColor(getLineColor(kEl));
    hYieldTOFElectrons->SetMarkerColor(getLineColor(kEl));
    hYieldTOFElectrons->SetMarkerStyle(22);
    hYieldTOFElectrons->SetStats(kFALSE);
    
    hYieldTPConlyTotal = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldTPConlyTotal", "TPC only total");
    hYieldTPConlyTotal->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    hYieldTPConlyTotal->SetLineColor(kBlack);
    hYieldTPConlyTotal->SetMarkerColor(kBlack);
    hYieldTPConlyTotal->SetStats(kFALSE);
    
    // Extract TOF purity in case of MC 
    if (isMC) {
      hTOFPurity = hPIDdata->Projection(indexTOFinfoAxisData, kPidMCpid, "e");
      hTOFPurity->SetName("hTOFPurity");
      hTOFPurity->SetTitle("TOF Purity");
      hTOFPurity->GetXaxis()->SetTitle("MC ID");
      hTOFPurity->GetYaxis()->SetTitle("TOF ID");
      
      for (Int_t binX = 1; binX <= hTOFPurity->GetNbinsX(); binX++)
        hTOFPurity->GetXaxis()->SetBinLabel(binX, hPIDdata->GetAxis(kPidMCpid)->GetBinLabel(binX));
      
      for (Int_t binY = 1; binY <= hTOFPurity->GetNbinsY(); binY++)
        hTOFPurity->GetYaxis()->SetBinLabel(binY, hPIDdata->GetAxis(indexTOFinfoAxisData)->GetBinLabel(binY));
      
      // Normalise rows to unity such that one can directly read off the purity
      for (Int_t binY = 1; binY <= hTOFPurity->GetNbinsY(); binY++) {
        Double_t sum = 0.;
        for (Int_t binX = 1; binX <= hTOFPurity->GetNbinsX(); binX++)
          sum += hTOFPurity->GetBinContent(binX, binY);
        if (sum > 0.) {
          for (Int_t binX = 1; binX <= hTOFPurity->GetNbinsX(); binX++) {
            hTOFPurity->SetBinContent(binX, binY, hTOFPurity->GetBinContent(binX, binY) / sum);
            hTOFPurity->SetBinError(binX, binY, hTOFPurity->GetBinError(binX, binY) / sum);
          }
        }
      }
    }
      
      
    // Count MC yield with TOF information separately (other will be TPC only), such that fit results can be compared to
    // MC TPC only truth and, finally, also TPC with TOF patching can be compared to full MC truth
    // Obtain MC information about particle yields.
    
    // First bin with TOF information is kTOFpion, last one is kTOFproton
    tofBin = hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFpion);
    Int_t tofBin2 = hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFproton);
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(tofBin, tofBin2);
    hMCdataTOF = (TH2D*)hPIDdata->Projection(kPidMCpid, axisForMode, "e");
    hMCdataTOF->SetName("hMCdataTOF");
    
    
    hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1);
    
    // In case of TOF patching: Restrict TOF info axis to no TOF region
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kNoTOFinfo),
                                                      hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kNoTOFpid));
  }

  // If jet axis is restricted, normalise to number of jets. If not, normalise to number of events
  Double_t numEvents = -1;
  Double_t numJetsRec = -1;
  Double_t numJetsGen = -1;
  TH1* hNumEvents = dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessed"));
  TH1* hNumEventsTriggerSel = dynamic_cast<TH1*>(histList->FindObject("fhEventsTriggerSel"));
  TH1* hNumEventsTriggerSelVtxCut = dynamic_cast<TH1*>(histList->FindObject("fhEventsTriggerSelVtxCut"));
  TH1* hNumEventsTriggerSelVtxCutNoPileUp = dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessedNoPileUpRejection"));
  
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;
  
  if (restrictJetPtAxis) {
    hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
    hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
    Bool_t createNewHistos = kFALSE;
    
    if (!hNjetsRec) {
      printf("Failed to load number of jets histo!\n");
      
      // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
      // period is considered; also: No multiplicity information)
      TString pathData = fileName;
      pathData.Replace(pathData.Last('/'), pathData.Length(), "");
      TString pathBackward = Form("%s/AnalysisResults.root", pathData.Data());
      TFile* fBackward = TFile::Open(pathBackward.Data());
    
     
      
      TString dirDataInFile = "";
      TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;

      TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

      TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCutsInc") : 0x0;
      TH1D* hFFJetPtGen = list ? (TH1D*)list->FindObject("fh1FFJetPtGenInc") : 0x0;
      
      if (hFFJetPtRec) {
        printf("***WARNING: For backward compatibility, using file \"%s\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n",
               pathBackward.Data());
        printf("ALSO: Using Njets for inclusive jets!!!!\n");
        
        createNewHistos = kTRUE;
        hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1,  hPIDdata->GetAxis(kPidJetPt)->GetNbins(),
                             hPIDdata->GetAxis(kPidJetPt)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
        }
      }
      
      if (!hNjetsRec)
        return -1;
      
      // Same for gen jets, if available
      if (hFFJetPtGen) {
        createNewHistos = kTRUE;
        hNjetsGen = new TH2D("fh2FFJetPtGen", "", 1, -1, 1,  hPIDdata->GetAxis(kPidJetPt)->GetNbins(),
                             hPIDdata->GetAxis(kPidJetPt)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsGen->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtGen->FindFixBin(hNjetsGen->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtGen->FindFixBin(hNjetsGen->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsGen->SetBinContent(1, iJet, hFFJetPtGen->Integral(lowerBin, upperBin));
        }
      }
    }
    
    numJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                 upperJetPtBinLimit) : -1.;
    numJetsRec = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                 upperJetPtBinLimit) : -1.;
    if (createNewHistos) {
      delete hNjetsGen;
      hNjetsGen = 0x0;
      
      delete hNjetsRec;
      hNjetsRec = 0x0;
    }
  }
  else {
    if (!hNumEvents) {
      std::cout << std::endl;
      std::cout << "Histo with number of processed events not found! Yields will NOT be normalised to this number!" << std::endl 
                << std::endl;
    }
    else {
      // Under- and overflow automatically taken into account for unrestricted axis by using -1 and -2 for limits
      numEvents = hNumEvents->Integral(lowerCentralityBinLimit, upperCentralityBinLimit);
      
      if (numEvents <= 0) {
        numEvents = -1;
        std::cout << std::endl;
        std::cout << "Number of processed events < 1 in selected range! Yields will NOT be normalised to this number!"
                  << std::endl << std::endl;
      }
    }
  }
  
  
  // TH1D hist with total yield per pT bin (-> project to arbitrary selectSpecies to avoid multiple counting)
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1);
  TH1D* hYieldPt = hPIDdata->Projection(axisForMode, "e");
  hYieldPt->SetName(Form("hYield%s", modeShortName[mode].Data()));
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1);
  
  
  // Fill \Delta\species histograms for each momentum slice
  Int_t nBins = hPIDdata->GetAxis(dataAxis)->GetNbins();
  Double_t xLow = hPIDdata->GetAxis(dataAxis)->GetXmin();
  Double_t xUp = hPIDdata->GetAxis(dataAxis)->GetXmax();
  
  const Int_t numSlices = (mode == kPMpT) ? nPtBins : hPIDdata->GetAxis(axisForMode)->GetNbins();
  
  TH1D* hDeltaPi[numSlices];
  TH1D* hDeltaEl[numSlices];
  TH1D* hDeltaKa[numSlices];
  TH1D* hDeltaPr[numSlices]; 
  
  TH1D* hDeltaPiFitQA[numSlices];
  TH1D* hDeltaElFitQA[numSlices];
  TH1D* hDeltaKaFitQA[numSlices];
  TH1D* hDeltaPrFitQA[numSlices];
    
  const Int_t nMCbins = 5;
  TH1D* hDeltaPiMC[numSlices][nMCbins];
  TH1D* hDeltaElMC[numSlices][nMCbins];
  TH1D* hDeltaKaMC[numSlices][nMCbins];
  TH1D* hDeltaPrMC[numSlices][nMCbins]; 
  
  
  TH2D* h2Delta[4];
  TH2D* h2DeltaMC[4][nMCbins];
  
  for (Int_t i = 0; i < 4; i++) {
    TString speciesLabel = hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(i + 1);
    
    hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(i + 1, i + 1);
    h2Delta[i] = hPIDdata->Projection(dataAxis, axisForMode, "e");
    h2Delta[i]->SetName(Form("h2Delta_%s", speciesLabel.Data()));
    h2Delta[i]->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
    h2Delta[i]->GetYaxis()->SetTitle(Form("#Delta%s = dE/dx %s <dE/dx>_{%s}", 
                                          useDeltaPrime ? Form("'_{#lower[-0.5]{%s}}", speciesLabel.Data()) 
                                                        : Form("_{%s}", speciesLabel.Data()),
                                          useDeltaPrime ? "/" : "-", speciesLabel.Data()));
    
    for (Int_t species = 0; species < nMCbins; species++) {
      hPIDdata->GetAxis(kPidMCpid)->SetRange(species + 1, species + 1); // Select MC species
      h2DeltaMC[i][species] = hPIDdata->Projection(dataAxis, axisGenForMode, "e");
      h2DeltaMC[i][species]->SetName(Form("h2Delta_MC_%s", speciesLabel.Data()));
      h2DeltaMC[i][species]->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisGenForMode)->GetTitle());
      h2DeltaMC[i][species]->GetYaxis()->SetTitle(h2Delta[i]->GetYaxis()->GetTitle());
    }
    hPIDdata->GetAxis(kPidMCpid)->SetRange(0, -1);
  }
  
  Int_t firstValidSlice = -1;
  for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < hPIDdata->GetAxis(axisForMode)->GetNbins(); slice++) {   
    if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
      continue; 
    
    if (firstValidSlice < 0)
      firstValidSlice = slice;
        
    // Add/subtract some very small offset to be sure not to sit on the bin boundary, when looking for the integration/projection limits.
    // For modes different from pT, just take 1 bin
    const Int_t pBinLowProjLimit = (mode == kPMpT) ? h2Delta[0]->GetXaxis()->FindFixBin(binsPt[slice] + 1e-5)    : slice + 1;
    const Int_t pBinUpProjLimit  = (mode == kPMpT) ? h2Delta[0]->GetXaxis()->FindFixBin(binsPt[slice + 1]- 1e-5) : slice + 1;
    
    const TString binInfo = (mode == kPMpT) ? Form("%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                            : Form("%.2f_%s_%.2f", hPIDdata->GetAxis(axisForMode)->GetBinLowEdge(pBinLowProjLimit), 
                                                                   modeShortName[mode].Data(),
                                                                   hPIDdata->GetAxis(axisForMode)->GetBinUpEdge(pBinUpProjLimit));
      
    hDeltaEl[slice] = h2Delta[0]->ProjectionY(Form("hDeltaEl_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
    hDeltaEl[slice]->GetXaxis()->SetTitle(h2Delta[0]->GetYaxis()->GetTitle());
    hDeltaEl[slice]->GetXaxis()->SetTitleOffset(1.0);
    hDeltaEl[slice]->SetStats(kFALSE);
    
    hDeltaKa[slice] = h2Delta[1]->ProjectionY(Form("hDeltaKa_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
    hDeltaKa[slice]->SetName(Form("hDeltaKa_%s", binInfo.Data()));
    hDeltaKa[slice]->GetXaxis()->SetTitle(h2Delta[1]->GetYaxis()->GetTitle());
    hDeltaKa[slice]->GetXaxis()->SetTitleOffset(1.0);
    hDeltaKa[slice]->SetStats(kFALSE);
    
    hDeltaPi[slice] = h2Delta[2]->ProjectionY(Form("hDeltaPi_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
    hDeltaPi[slice]->SetName(Form("hDeltaPi_%s", binInfo.Data()));
    hDeltaPi[slice]->GetXaxis()->SetTitle(h2Delta[2]->GetYaxis()->GetTitle());
    hDeltaPi[slice]->GetXaxis()->SetTitleOffset(1.0);
    hDeltaPi[slice]->SetStats(kFALSE);
    
    hDeltaPr[slice] = h2Delta[3]->ProjectionY(Form("hDeltaPr_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
    hDeltaPr[slice]->SetName(Form("hDeltaPr_%s", binInfo.Data()));
    hDeltaPr[slice]->GetXaxis()->SetTitle(h2Delta[3]->GetYaxis()->GetTitle());
    hDeltaPr[slice]->GetXaxis()->SetTitleOffset(1.0);
    hDeltaPr[slice]->SetStats(kFALSE);
    
    // For most prob PID low-pT templates this is also needed, therefore extract it always!
    //if (plotIdentifiedSpectra)  {
      // If identified spectra are available (mainly useful in the MC case) and shall be used,
      // create histos with signals from identified particles
  
      // DeltaEl
      for (Int_t species = 0; species < nMCbins; species++) {
        hDeltaElMC[slice][species] = h2DeltaMC[0][species]->ProjectionY(Form("hDeltaElMC_%s_species_%d", binInfo.Data(), species),
                                                                        pBinLowProjLimit, pBinUpProjLimit, "e");
        hDeltaElMC[slice][species]->SetLineColor(getLineColor(species + 1));
        hDeltaElMC[slice][species]->SetMarkerColor(getLineColor(species + 1));
        hDeltaElMC[slice][species]->SetMarkerStyle(24);
        hDeltaElMC[slice][species]->SetLineStyle(1);
        hDeltaElMC[slice][species]->GetXaxis()->SetTitle(h2DeltaMC[0][species]->GetYaxis()->GetTitle());
        hDeltaElMC[slice][species]->GetXaxis()->SetTitleOffset(1.0);
        hDeltaElMC[slice][species]->SetStats(kFALSE);
      }
      
      // DeltaKa
      for (Int_t species = 0; species < nMCbins; species++) {
        hDeltaKaMC[slice][species] = h2DeltaMC[1][species]->ProjectionY(Form("hDeltaKaMC_%s_species_%d", binInfo.Data(), species),
                                                                             pBinLowProjLimit, pBinUpProjLimit, "e");
        hDeltaKaMC[slice][species]->SetLineColor(getLineColor(species + 1));
        hDeltaKaMC[slice][species]->SetMarkerColor(getLineColor(species + 1));
        hDeltaKaMC[slice][species]->SetMarkerStyle(24);
        hDeltaKaMC[slice][species]->SetLineStyle(1);
        hDeltaKaMC[slice][species]->GetXaxis()->SetTitle(h2DeltaMC[1][species]->GetYaxis()->GetTitle());
        hDeltaKaMC[slice][species]->GetXaxis()->SetTitleOffset(1.0);
        hDeltaKaMC[slice][species]->SetStats(kFALSE);
      }
      
      // DeltaPi
      for (Int_t species = 0; species < nMCbins; species++) {
        hDeltaPiMC[slice][species] = h2DeltaMC[2][species]->ProjectionY(Form("hDeltaPiMC_%s_species_%d", binInfo.Data(), species),
                                                                             pBinLowProjLimit, pBinUpProjLimit, "e");
        hDeltaPiMC[slice][species]->SetLineColor(getLineColor(species + 1));
        hDeltaPiMC[slice][species]->SetMarkerColor(getLineColor(species + 1));
        hDeltaPiMC[slice][species]->SetMarkerStyle(24);
        hDeltaPiMC[slice][species]->SetLineStyle(1);
        hDeltaPiMC[slice][species]->GetXaxis()->SetTitle(h2DeltaMC[2][species]->GetYaxis()->GetTitle());
        hDeltaPiMC[slice][species]->GetXaxis()->SetTitleOffset(1.0);
        hDeltaPiMC[slice][species]->SetStats(kFALSE);
      }
      
      // DeltaPr
      for (Int_t species = 0; species < nMCbins; species++) {
        hDeltaPrMC[slice][species] = h2DeltaMC[3][species]->ProjectionY(Form("hDeltaPrMC_%s_species_%d", binInfo.Data(), species),
                                                                             pBinLowProjLimit, pBinUpProjLimit, "e");
        hDeltaPrMC[slice][species]->SetLineColor(getLineColor(species + 1));
        hDeltaPrMC[slice][species]->SetMarkerColor(getLineColor(species + 1));
        hDeltaPrMC[slice][species]->SetMarkerStyle(24);
        hDeltaPrMC[slice][species]->SetLineStyle(1);
        hDeltaPrMC[slice][species]->GetXaxis()->SetTitle(h2DeltaMC[3][species]->GetYaxis()->GetTitle());
        hDeltaPrMC[slice][species]->GetXaxis()->SetTitleOffset(1.0);
        hDeltaPrMC[slice][species]->SetStats(kFALSE);
      }
    //}
  }
  // In case of TOF patching, also get data points for the TOF pion bin
  TH1D* hDeltaPiTOFpi[numSlices];
  TH1D* hDeltaElTOFpi[numSlices];
  TH1D* hDeltaKaTOFpi[numSlices];
  TH1D* hDeltaPrTOFpi[numSlices]; 
  
  TH1D* hDeltaPiMCTOFpi[numSlices][nMCbins];
  TH1D* hDeltaElMCTOFpi[numSlices][nMCbins];
  TH1D* hDeltaKaMCTOFpi[numSlices][nMCbins];
  TH1D* hDeltaPrMCTOFpi[numSlices][nMCbins]; 
  
  for (Int_t i = 0; i < numSlices; i++) {
    hDeltaPiTOFpi[i] = 0x0;
    hDeltaElTOFpi[i] = 0x0;
    hDeltaKaTOFpi[i] = 0x0;
    hDeltaPrTOFpi[i] = 0x0;
    
    for (Int_t j = 0; j < nMCbins; j++) {
      hDeltaPiMCTOFpi[i][j] = 0x0;
      hDeltaElMCTOFpi[i][j] = 0x0;
      hDeltaKaMCTOFpi[i][j] = 0x0;
      hDeltaPrMCTOFpi[i][j] = 0x0;
    }
  }
  
  TH2D* h2DeltaTOFpi[4] = { 0x0, };
  TH2D* h2DeltaMCTOFpi[4][nMCbins];
  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < nMCbins; j++) {
      h2DeltaMCTOFpi[i][j] = 0x0;
    }
  }
  
  if (applyTOFpatching) {
    // Set to TOF pion range
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFpion),
                                                      hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kTOFpion));
    
    for (Int_t i = 0; i < 4; i++) {
      TString speciesLabel = hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(i + 1);
      
      hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(i + 1, i + 1);
      h2DeltaTOFpi[i] = hPIDdata->Projection(dataAxis, axisForMode, "e");
      h2DeltaTOFpi[i]->SetName(Form("h2DeltaTOFpi_%s", speciesLabel.Data()));
      h2DeltaTOFpi[i]->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
      h2DeltaTOFpi[i]->GetYaxis()->SetTitle(Form("#Delta%s = dE/dx %s <dE/dx>_{%s}", 
                                            useDeltaPrime ? Form("'_{#lower[-0.5]{%s}}", speciesLabel.Data()) 
                                                          : Form("_{%s}", speciesLabel.Data()),
                                            useDeltaPrime ? "/" : "-", speciesLabel.Data()));
      
      for (Int_t species = 0; species < nMCbins; species++) {
        hPIDdata->GetAxis(kPidMCpid)->SetRange(species + 1, species + 1); // Select MC species
        h2DeltaMCTOFpi[i][species] = hPIDdata->Projection(dataAxis, axisGenForMode, "e");
        h2DeltaMCTOFpi[i][species]->SetName(Form("h2DeltaTOFpi_MC_%s", speciesLabel.Data()));
        h2DeltaMCTOFpi[i][species]->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisGenForMode)->GetTitle());
        h2DeltaMCTOFpi[i][species]->GetYaxis()->SetTitle(h2Delta[i]->GetYaxis()->GetTitle());
      }
      hPIDdata->GetAxis(kPidMCpid)->SetRange(0, -1);
    }
    
    for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < hPIDdata->GetAxis(axisForMode)->GetNbins(); slice++) {   
      if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
        continue; 
      
      // Add/subtract some very small offset to be sure not to sit on the bin boundary, when looking for the integration/projection limits.
      // For modes different from pT, just take 1 bin
      const Int_t pBinLowProjLimit = (mode == kPMpT) ? h2DeltaTOFpi[0]->GetXaxis()->FindFixBin(binsPt[slice] + 1e-5)    : slice + 1;
      const Int_t pBinUpProjLimit  = (mode == kPMpT) ? h2DeltaTOFpi[0]->GetXaxis()->FindFixBin(binsPt[slice + 1]- 1e-5) : slice + 1;
      
      const TString binInfo = (mode == kPMpT) ? Form("%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                              : Form("%.2f_%s_%.2f", hPIDdata->GetAxis(axisForMode)->GetBinLowEdge(pBinLowProjLimit), 
                                                                    modeShortName[mode].Data(),
                                                                    hPIDdata->GetAxis(axisForMode)->GetBinUpEdge(pBinUpProjLimit));
        
      hDeltaElTOFpi[slice] = h2DeltaTOFpi[0]->ProjectionY(Form("hDeltaElTOFpi_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
      hDeltaElTOFpi[slice]->GetXaxis()->SetTitle(h2DeltaTOFpi[0]->GetYaxis()->GetTitle());
      hDeltaElTOFpi[slice]->GetXaxis()->SetTitleOffset(1.0);
      hDeltaElTOFpi[slice]->SetStats(kFALSE);
      
      hDeltaKaTOFpi[slice] = h2DeltaTOFpi[1]->ProjectionY(Form("hDeltaKaTOFpi_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
      hDeltaKaTOFpi[slice]->GetXaxis()->SetTitle(h2DeltaTOFpi[1]->GetYaxis()->GetTitle());
      hDeltaKaTOFpi[slice]->GetXaxis()->SetTitleOffset(1.0);
      hDeltaKaTOFpi[slice]->SetStats(kFALSE);
      
      hDeltaPiTOFpi[slice] = h2DeltaTOFpi[2]->ProjectionY(Form("hDeltaPiTOFpi_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
      hDeltaPiTOFpi[slice]->GetXaxis()->SetTitle(h2DeltaTOFpi[2]->GetYaxis()->GetTitle());
      hDeltaPiTOFpi[slice]->GetXaxis()->SetTitleOffset(1.0);
      hDeltaPiTOFpi[slice]->SetStats(kFALSE);
      
      hDeltaPrTOFpi[slice] = h2DeltaTOFpi[3]->ProjectionY(Form("hDeltaPrTOFpi_%s", binInfo.Data()), pBinLowProjLimit, pBinUpProjLimit, "e");
      hDeltaPrTOFpi[slice]->GetXaxis()->SetTitle(h2DeltaTOFpi[3]->GetYaxis()->GetTitle());
      hDeltaPrTOFpi[slice]->GetXaxis()->SetTitleOffset(1.0);
      hDeltaPrTOFpi[slice]->SetStats(kFALSE);
      
      // For most prob PID low-pT templates this is also needed, therefore extract it always!
      //if (plotIdentifiedSpectra)  {
        // If identified spectra are available (mainly useful in the MC case) and shall be used,
        // create histos with signals from identified particles
    
        // DeltaEl
        for (Int_t species = 0; species < nMCbins; species++) {
          hDeltaElMCTOFpi[slice][species] = h2DeltaMCTOFpi[0][species]->ProjectionY(Form("hDeltaElMCTOFpi_%s_species_%d", binInfo.Data(),
                                                                                         species), pBinLowProjLimit, pBinUpProjLimit, "e");
          hDeltaElMCTOFpi[slice][species]->SetLineColor(getLineColor(species + 1));
          hDeltaElMCTOFpi[slice][species]->SetMarkerColor(getLineColor(species + 1));
          hDeltaElMCTOFpi[slice][species]->SetMarkerStyle(24);
          hDeltaElMCTOFpi[slice][species]->SetLineStyle(1);
          hDeltaElMCTOFpi[slice][species]->GetXaxis()->SetTitle(h2DeltaMCTOFpi[0][species]->GetYaxis()->GetTitle());
          hDeltaElMCTOFpi[slice][species]->GetXaxis()->SetTitleOffset(1.0);
          hDeltaElMCTOFpi[slice][species]->SetStats(kFALSE);
        }
        
        // DeltaKa
        for (Int_t species = 0; species < nMCbins; species++) {
          hDeltaKaMCTOFpi[slice][species] = h2DeltaMCTOFpi[1][species]->ProjectionY(Form("hDeltaKaMCTOFpi_%s_species_%d", binInfo.Data(),
                                                                                         species), pBinLowProjLimit, pBinUpProjLimit, "e");
          hDeltaKaMCTOFpi[slice][species]->SetLineColor(getLineColor(species + 1));
          hDeltaKaMCTOFpi[slice][species]->SetMarkerColor(getLineColor(species + 1));
          hDeltaKaMCTOFpi[slice][species]->SetMarkerStyle(24);
          hDeltaKaMCTOFpi[slice][species]->SetLineStyle(1);
          hDeltaKaMCTOFpi[slice][species]->GetXaxis()->SetTitle(h2DeltaMCTOFpi[1][species]->GetYaxis()->GetTitle());
          hDeltaKaMCTOFpi[slice][species]->GetXaxis()->SetTitleOffset(1.0);
          hDeltaKaMCTOFpi[slice][species]->SetStats(kFALSE);
        }
        
        // DeltaPi
        for (Int_t species = 0; species < nMCbins; species++) {
          hDeltaPiMCTOFpi[slice][species] = h2DeltaMCTOFpi[2][species]->ProjectionY(Form("hDeltaPiMCTOFpi_%s_species_%d", binInfo.Data(), 
                                                                                         species), pBinLowProjLimit, pBinUpProjLimit, "e");
          hDeltaPiMCTOFpi[slice][species]->SetLineColor(getLineColor(species + 1));
          hDeltaPiMCTOFpi[slice][species]->SetMarkerColor(getLineColor(species + 1));
          hDeltaPiMCTOFpi[slice][species]->SetMarkerStyle(24);
          hDeltaPiMCTOFpi[slice][species]->SetLineStyle(1);
          hDeltaPiMCTOFpi[slice][species]->GetXaxis()->SetTitle(h2DeltaMCTOFpi[2][species]->GetYaxis()->GetTitle());
          hDeltaPiMCTOFpi[slice][species]->GetXaxis()->SetTitleOffset(1.0);
          hDeltaPiMCTOFpi[slice][species]->SetStats(kFALSE);
        }
        
        // DeltaPr
        for (Int_t species = 0; species < nMCbins; species++) {
          hDeltaPrMCTOFpi[slice][species] = h2DeltaMCTOFpi[3][species]->ProjectionY(Form("hDeltaPrMCTOFpi_%s_species_%d", binInfo.Data(),
                                                                                         species), pBinLowProjLimit, pBinUpProjLimit, "e");
          hDeltaPrMCTOFpi[slice][species]->SetLineColor(getLineColor(species + 1));
          hDeltaPrMCTOFpi[slice][species]->SetMarkerColor(getLineColor(species + 1));
          hDeltaPrMCTOFpi[slice][species]->SetMarkerStyle(24);
          hDeltaPrMCTOFpi[slice][species]->SetLineStyle(1);
          hDeltaPrMCTOFpi[slice][species]->GetXaxis()->SetTitle(h2DeltaMCTOFpi[3][species]->GetYaxis()->GetTitle());
          hDeltaPrMCTOFpi[slice][species]->GetXaxis()->SetTitleOffset(1.0);
          hDeltaPrMCTOFpi[slice][species]->SetStats(kFALSE);
        }
      //}
    }
    
    // Set back to original range (i.e. no TOF)
    hPIDdata->GetAxis(indexTOFinfoAxisData)->SetRange(hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kNoTOFinfo),
                                                      hPIDdata->GetAxis(indexTOFinfoAxisData)->FindFixBin(kNoTOFpid));
  }
  
  
  
  hPIDdata->GetAxis(kPidMCpid)->SetRange(0, -1);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1);
  
  hPIDdata->GetAxis(axisForMode)->SetRange(0, -1);
  
  // Start fitting of slices
  TCanvas* cSingleFit[numSlices][4];
  
  TF1*     fitFuncTotalDeltaPion[numSlices];
  TF1*     fitFuncTotalDeltaElectron[numSlices];
  TF1*     fitFuncTotalDeltaKaon[numSlices];
  TF1*     fitFuncTotalDeltaProton[numSlices];
  
  // Histos for particle fractions
  TH1F* hFractionElectrons = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hFractionElectrons", "e");
  hFractionElectrons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
  hFractionElectrons->GetYaxis()->SetTitle("Fraction");
  hFractionElectrons->SetLineColor(getLineColor(kEl));
  hFractionElectrons->SetMarkerColor(getLineColor(kEl));
  hFractionElectrons->SetMarkerStyle(20);
  hFractionElectrons->Sumw2();
  hFractionElectrons->SetStats(kFALSE);
  
  TH1F* hFractionElectronsDeltaPion = (TH1F*)hFractionElectrons->Clone("hFractionElectronsDeltaPion");
  TH1F* hFractionElectronsDeltaElectron = (TH1F*)hFractionElectrons->Clone("hFractionElectronsDeltaElectron");
  TH1F* hFractionElectronsDeltaKaon = (TH1F*)hFractionElectrons->Clone("hFractionElectronsDeltaKaon");
  TH1F* hFractionElectronsDeltaProton = (TH1F*)hFractionElectrons->Clone("hFractionElectronsDeltaProton");
  
  TH1F* hFractionKaons = (TH1F*)hFractionElectrons->Clone("hFractionKaons");
  hFractionKaons->SetTitle("K");
  hFractionKaons->SetLineColor(getLineColor(kKa));
  hFractionKaons->SetMarkerColor(getLineColor(kKa));
  
  TH1F* hFractionKaonsDeltaPion = (TH1F*)hFractionKaons->Clone("hFractionKaonsDeltaPion");
  TH1F* hFractionKaonsDeltaElectron = (TH1F*)hFractionKaons->Clone("hFractionKaonsDeltaElectron");
  TH1F* hFractionKaonsDeltaKaon = (TH1F*)hFractionKaons->Clone("hFractionKaonsDeltaKaon");
  TH1F* hFractionKaonsDeltaProton = (TH1F*)hFractionKaons->Clone("hFractionKaonsDeltaProton");
  
  TH1F* hFractionPions = (TH1F*)hFractionElectrons->Clone("hFractionPions");
  hFractionPions->SetTitle("#pi");
  hFractionPions->SetLineColor(getLineColor(kPi));
  hFractionPions->SetMarkerColor(getLineColor(kPi));
  
  TH1F* hFractionPionsDeltaPion = (TH1F*)hFractionPions->Clone("hFractionPionsDeltaPion");
  TH1F* hFractionPionsDeltaElectron = (TH1F*)hFractionPions->Clone("hFractionPionsDeltaElectron");
  TH1F* hFractionPionsDeltaKaon = (TH1F*)hFractionPions->Clone("hFractionPionsDeltaKaon");
  TH1F* hFractionPionsDeltaProton = (TH1F*)hFractionPions->Clone("hFractionPionsDeltaProton");
  
  TH1F* hFractionProtons = (TH1F*)hFractionElectrons->Clone("hFractionProtons");
  hFractionProtons->SetTitle("p");
  hFractionProtons->SetLineColor(getLineColor(kPr));
  hFractionProtons->SetMarkerColor(getLineColor(kPr));
 
  TH1F* hFractionProtonsDeltaPion = (TH1F*)hFractionProtons->Clone("hFractionProtonsDeltaPion");
  TH1F* hFractionProtonsDeltaElectron = (TH1F*)hFractionProtons->Clone("hFractionProtonsDeltaElectron");
  TH1F* hFractionProtonsDeltaKaon = (TH1F*)hFractionProtons->Clone("hFractionProtonsDeltaKaon");
  TH1F* hFractionProtonsDeltaProton = (TH1F*)hFractionProtons->Clone("hFractionProtonsDeltaProton");
  
  TH1F* hFractionMuons = (TH1F*)hFractionElectrons->Clone("hFractionMuons");
  hFractionMuons->SetTitle("#mu");
  hFractionMuons->SetLineColor(getLineColor(kMu));
  hFractionMuons->SetMarkerColor(getLineColor(kMu));
 
  TH1F* hFractionMuonsDeltaPion = (TH1F*)hFractionMuons->Clone("hFractionMuonsDeltaPion");
  TH1F* hFractionMuonsDeltaElectron = (TH1F*)hFractionMuons->Clone("hFractionMuonsDeltaElectron");
  TH1F* hFractionMuonsDeltaKaon = (TH1F*)hFractionMuons->Clone("hFractionMuonsDeltaKaon");
  TH1F* hFractionMuonsDeltaProton = (TH1F*)hFractionMuons->Clone("hFractionMuonsDeltaProton");

  TH1F* hFractionSummed = (TH1F*)hFractionProtons->Clone("hFractionSummed");
  hFractionSummed->SetTitle("Sum");
  hFractionSummed->SetLineColor(kBlack);
  hFractionSummed->SetMarkerColor(kBlack);
  
  
  // MC fractions
  TH1F* hFractionElectronsMC = (TH1F*)hFractionElectrons->Clone("hFractionElectronsMC");
  hFractionElectronsMC->SetMarkerStyle(24);
  TH1F* hFractionKaonsMC = (TH1F*)hFractionKaons->Clone("hFractionKaonsMC");
  hFractionKaonsMC->SetMarkerStyle(24);
  TH1F* hFractionPionsMC = (TH1F*)hFractionPions->Clone("hFractionPionsMC");
  hFractionPionsMC->SetMarkerStyle(24);
  TH1F* hFractionMuonsMC = (TH1F*)hFractionMuons->Clone("hFractionMuonsMC");
  hFractionMuonsMC->SetMarkerStyle(24);
  TH1F* hFractionProtonsMC = (TH1F*)hFractionProtons->Clone("hFractionProtonsMC");
  hFractionProtonsMC->SetMarkerStyle(24);
  
  
  // Comparison fit result<->MC
  TString fractionComparisonTitle = Form("Fraction fit / fraction %s", identifiedLabels[isMC].Data()); 
  TH1F* hFractionComparisonElectrons = (TH1F*)hFractionElectrons->Clone("hFractionComparisonElectrons");
  hFractionComparisonElectrons->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  TH1F* hFractionComparisonMuons = (TH1F*)hFractionMuons->Clone("hFractionComparisonMuons");
  hFractionComparisonMuons->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  TH1F* hFractionComparisonKaons = (TH1F*)hFractionKaons->Clone("hFractionComparisonKaons");
  hFractionComparisonKaons->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  TH1F* hFractionComparisonPions = (TH1F*)hFractionPions->Clone("hFractionComparisonPions");
  hFractionComparisonPions->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  TH1F* hFractionComparisonProtons = (TH1F*)hFractionProtons->Clone("hFractionComparisonProtons");
  hFractionComparisonProtons->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  TH1F* hFractionComparisonTotal = (TH1F*)hFractionSummed->Clone("hFractionComparisonTotal");
  hFractionComparisonTotal->GetYaxis()->SetTitle(fractionComparisonTitle.Data());
  
  
  
  // Histos for particle yields
  TH1F* hYieldElectrons = GetHistWithProperXbinning(hPIDdata, axisForMode, mode, "hYieldElectrons", "e");
  
  TString axisTitleYield = "";
  // Other title and normalisation when looking at jets
  if (restrictJetPtAxis) {
    axisTitleYield = Form("%sdN/d%s%s", numJetsRec > 0 ? "1/N_{Jets} " : "", modeLatexName[mode].Data(),
                          mode == kPMpT ? " (GeV/c)^{-1}" : "");
  }
  else {
    axisTitleYield = Form("%s1/(2#pi%s) d^{2}N/d#etad%s%s", numEvents > 0 ? "1/N_{ev} " : "",
                          modeLatexName[mode].Data(), modeLatexName[mode].Data(),
                          mode == kPMpT ? " (GeV/c)^{-2}" : "");
  }
  
  if (applyTOFpatching) {
    hYieldTOFOrigBinningPions->GetYaxis()->SetTitle(axisTitleYield.Data());
    hYieldTOFOrigBinningKaons->GetYaxis()->SetTitle(axisTitleYield.Data());
    hYieldTOFOrigBinningProtons->GetYaxis()->SetTitle(axisTitleYield.Data());
    hYieldTOFPions->GetYaxis()->SetTitle(axisTitleYield.Data());
    hYieldTOFKaons->GetYaxis()->SetTitle(axisTitleYield.Data());
    hYieldTOFProtons->GetYaxis()->SetTitle(axisTitleYield.Data());
  }
  
  hYieldElectrons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
  hYieldElectrons->GetYaxis()->SetTitle(axisTitleYield.Data());
  hYieldElectrons->SetLineColor(getLineColor(kEl));
  hYieldElectrons->SetMarkerColor(getLineColor(kEl));
  hYieldElectrons->SetMarkerStyle(20);
  hYieldElectrons->Sumw2();
  hYieldElectrons->SetStats(kFALSE);
  
  TH1F* hYieldElectronsDeltaPion = (TH1F*)hYieldElectrons->Clone("hYieldElectronsDeltaPion");
  TH1F* hYieldElectronsDeltaElectron = (TH1F*)hYieldElectrons->Clone("hYieldElectronsDeltaElectron");
  TH1F* hYieldElectronsDeltaKaon = (TH1F*)hYieldElectrons->Clone("hYieldElectronsDeltaKaon");
  TH1F* hYieldElectronsDeltaProton = (TH1F*)hYieldElectrons->Clone("hYieldElectronsDeltaProton");
  
  TH1F* hYieldKaons = (TH1F*)hYieldElectrons->Clone("hYieldKaons");
  hYieldKaons->SetTitle("K");
  hYieldKaons->SetLineColor(getLineColor(kKa));
  hYieldKaons->SetMarkerColor(getLineColor(kKa));
  
  TH1F* hYieldKaonsDeltaPion = (TH1F*)hYieldKaons->Clone("hYieldKaonsDeltaPion");
  TH1F* hYieldKaonsDeltaElectron = (TH1F*)hYieldKaons->Clone("hYieldKaonsDeltaElectron");
  TH1F* hYieldKaonsDeltaKaon = (TH1F*)hYieldKaons->Clone("hYieldKaonsDeltaKaon");
  TH1F* hYieldKaonsDeltaProton = (TH1F*)hYieldKaons->Clone("hYieldKaonsDeltaProton");
  
  TH1F* hYieldPions = (TH1F*)hYieldElectrons->Clone("hYieldPions");
  hYieldPions->SetTitle("#pi");
  hYieldPions->SetLineColor(getLineColor(kPi));
  hYieldPions->SetMarkerColor(getLineColor(kPi));
  
  TH1F* hYieldPionsDeltaPion = (TH1F*)hYieldPions->Clone("hYieldPionsDeltaPion");
  TH1F* hYieldPionsDeltaElectron = (TH1F*)hYieldPions->Clone("hYieldPionsDeltaElectron");
  TH1F* hYieldPionsDeltaKaon = (TH1F*)hYieldPions->Clone("hYieldPionsDeltaKaon");
  TH1F* hYieldPionsDeltaProton = (TH1F*)hYieldPions->Clone("hYieldPionsDeltaProton");
  
  TH1F* hYieldProtons = (TH1F*)hYieldElectrons->Clone("hYieldProtons");
  hYieldProtons->SetTitle("p");
  hYieldProtons->SetLineColor(getLineColor(kPr));
  hYieldProtons->SetMarkerColor(getLineColor(kPr));
  
  TH1F* hYieldProtonsDeltaPion = (TH1F*)hYieldProtons->Clone("hYieldProtonsDeltaPion");
  TH1F* hYieldProtonsDeltaElectron = (TH1F*)hYieldProtons->Clone("hYieldProtonsDeltaElectron");
  TH1F* hYieldProtonsDeltaKaon = (TH1F*)hYieldProtons->Clone("hYieldProtonsDeltaKaon");
  TH1F* hYieldProtonsDeltaProton = (TH1F*)hYieldProtons->Clone("hYieldProtonsDeltaProton");
  
  TH1F* hYieldMuons = (TH1F*)hYieldElectrons->Clone("hYieldMuons");
  hYieldMuons->SetTitle("#mu");
  hYieldMuons->SetLineColor(getLineColor(kMu));
  hYieldMuons->SetMarkerColor(getLineColor(kMu));
  
  TH1F* hYieldMuonsDeltaPion = (TH1F*)hYieldMuons->Clone("hYieldMuonsDeltaPion");
  TH1F* hYieldMuonsDeltaElectron = (TH1F*)hYieldMuons->Clone("hYieldMuonsDeltaElectron");
  TH1F* hYieldMuonsDeltaKaon = (TH1F*)hYieldMuons->Clone("hYieldMuonsDeltaKaon");
  TH1F* hYieldMuonsDeltaProton = (TH1F*)hYieldMuons->Clone("hYieldMuonsDeltaProton");
  
  // MC yields
  TH1F* hYieldElectronsMC = (TH1F*)hYieldElectrons->Clone("hYieldElectronsMC");
  hYieldElectronsMC->SetMarkerStyle(24);
  TH1F* hYieldMuonsMC = (TH1F*)hYieldElectrons->Clone("hYieldMuonsMC");
  hYieldMuonsMC->SetMarkerStyle(24);
  hYieldMuonsMC->SetLineColor(getLineColor(kMu));
  hYieldMuonsMC->SetMarkerColor(getLineColor(kMu));
  TH1F* hYieldKaonsMC = (TH1F*)hYieldKaons->Clone("hYieldKaonsMC");
  hYieldKaonsMC->SetMarkerStyle(24);
  TH1F* hYieldPionsMC = (TH1F*)hYieldPions->Clone("hYieldPionsMC");
  hYieldPionsMC->SetMarkerStyle(24);
  TH1F* hYieldProtonsMC = (TH1F*)hYieldProtons->Clone("hYieldProtonsMC");
  hYieldProtonsMC->SetMarkerStyle(24);
  TH1F* hYieldSummedMC = (TH1F*)hYieldProtonsMC->Clone("hYieldSummedMC");
  hYieldSummedMC->SetTitle("Sum");
  hYieldSummedMC->SetLineColor(kBlack);
  hYieldSummedMC->SetMarkerColor(kBlack);
  
  TH1F* hYieldTOFElectronsMC = 0x0;
  TH1F* hYieldTOFMuonsMC = 0x0;
  TH1F* hYieldTOFKaonsMC = 0x0;
  TH1F* hYieldTOFPionsMC = 0x0;
  TH1F* hYieldTOFProtonsMC = 0x0;
  
  if (applyTOFpatching) {
    hYieldTOFElectronsMC = (TH1F*)hYieldElectronsMC->Clone("hYieldTOFElectronsMC");
    hYieldTOFElectronsMC->SetMarkerStyle(26);
    hYieldTOFMuonsMC = (TH1F*)hYieldMuonsMC->Clone("hYieldTOFMuonsMC");
    hYieldTOFMuonsMC->SetMarkerStyle(26);
    hYieldTOFKaonsMC = (TH1F*)hYieldKaonsMC->Clone("hYieldTOFKaonsMC");
    hYieldTOFKaonsMC->SetMarkerStyle(26);
    hYieldTOFPionsMC = (TH1F*)hYieldPionsMC->Clone("hYieldTOFPionsMC");
    hYieldTOFPionsMC->SetMarkerStyle(26);
    hYieldTOFProtonsMC = (TH1F*)hYieldProtonsMC->Clone("hYieldTOFProtonsMC");
    hYieldTOFProtonsMC->SetMarkerStyle(26);
  }
  
  // Comparison fit result<->MC (yields)
  TString yieldComparisonTitle = Form("Yield fit / yield %s", identifiedLabels[isMC].Data()); 
  TH1F* hYieldComparisonElectrons = (TH1F*)hYieldElectrons->Clone("hYieldComparisonElectrons");
  hYieldComparisonElectrons->GetYaxis()->SetTitle(yieldComparisonTitle.Data());
  TH1F* hYieldComparisonMuons = (TH1F*)hYieldMuons->Clone("hYieldComparisonMuons");
  hYieldComparisonMuons->GetYaxis()->SetTitle(yieldComparisonTitle.Data());
  TH1F* hYieldComparisonKaons = (TH1F*)hYieldKaons->Clone("hYieldComparisonKaons");
  hYieldComparisonKaons->GetYaxis()->SetTitle(yieldComparisonTitle.Data());
  TH1F* hYieldComparisonPions = (TH1F*)hYieldPions->Clone("hYieldComparisonPions");
  hYieldComparisonPions->GetYaxis()->SetTitle(yieldComparisonTitle.Data());
  TH1F* hYieldComparisonProtons = (TH1F*)hYieldProtons->Clone("hYieldComparisonProtons");
  hYieldComparisonProtons->GetYaxis()->SetTitle(yieldComparisonTitle.Data());
  
  
  // To-pi ratios
  TString electronString[3] = { "e^{-}", "e^{+}+e^{-}", "e^{+}" };
  TString muonString[3] = { "#mu^{-}", "#mu^{+}+#mu^{-}", "#mu^{+}" };
  TString kaonString[3] = { "K^{-}", "K^{+}+K^{-}", "K^{+}" };
  TString pionString[3] = { "#pi^{-}", "#pi^{+}+#pi^{-}", "#pi^{+}" };
  TString protonString[3] = { "#bar{p}", "p+#bar{p}", "p" };
  
  TH1F* hRatioToPiElectrons = (TH1F*)hYieldElectrons->Clone("hRatioToPiElectrons");
  hRatioToPiElectrons->GetYaxis()->SetTitle(Form("d%sN_{%s}/%sd%s / d%sN_{%s}/%sd%s",
                                                 restrictJetPtAxis ? "" : "^{2}",
                                                 electronString[chargeMode+1].Data(),
                                                 restrictJetPtAxis ? "" : "d#eta",
                                                 modeLatexName[mode].Data(),
                                                 restrictJetPtAxis ? "" : "^{2}",
                                                 pionString[chargeMode+1].Data(),
                                                 restrictJetPtAxis ? "" : "d#eta",
                                                 modeLatexName[mode].Data()));
  hRatioToPiElectrons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", electronString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", electronString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  
  TH1F* hRatioToPiMuons = (TH1F*)hYieldMuons->Clone("hRatioToPiMuons");
  hRatioToPiMuons->GetYaxis()->SetTitle(Form("d%sN_{%s}/%sd%s / d%sN_{%s}/%sd%s",
                                             restrictJetPtAxis ? "" : "^{2}",
                                             muonString[chargeMode+1].Data(),
                                             restrictJetPtAxis ? "" : "d#eta",
                                             modeLatexName[mode].Data(),
                                             restrictJetPtAxis ? "" : "^{2}",
                                             pionString[chargeMode+1].Data(),
                                             restrictJetPtAxis ? "" : "d#eta",
                                             modeLatexName[mode].Data()));
  hRatioToPiMuons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", muonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", muonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  TH1F* hRatioToPiKaons = (TH1F*)hYieldKaons->Clone("hRatioToPiKaons");
  hRatioToPiKaons->GetYaxis()->SetTitle(Form("d%sN_{%s}/%sd%s / d%sN_{%s}/%sd%s",
                                             restrictJetPtAxis ? "" : "^{2}",
                                             kaonString[chargeMode+1].Data(),
                                             restrictJetPtAxis ? "" : "d#eta",
                                             modeLatexName[mode].Data(),
                                             restrictJetPtAxis ? "" : "^{2}",
                                             pionString[chargeMode+1].Data(),
                                             restrictJetPtAxis ? "" : "d#eta",
                                             modeLatexName[mode].Data()));
  hRatioToPiKaons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", kaonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", kaonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  TH1F* hRatioToPiProtons = (TH1F*)hYieldProtons->Clone("hRatioToPiProtons");
  hRatioToPiProtons->GetYaxis()->SetTitle(Form("d%sN_{%s}/%sd%s / d%sN_{%s}/%sd%s",
                                               restrictJetPtAxis ? "" : "^{2}",
                                               protonString[chargeMode+1].Data(),
                                               restrictJetPtAxis ? "" : "d#eta",
                                               modeLatexName[mode].Data(),
                                               restrictJetPtAxis ? "" : "^{2}",
                                               pionString[chargeMode+1].Data(),
                                               restrictJetPtAxis ? "" : "d#eta",
                                               modeLatexName[mode].Data()));
  hRatioToPiProtons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", protonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", protonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  // MC to-pi ratios
  TH1F* hRatioToPiElectronsMC = (TH1F*)hRatioToPiElectrons->Clone("hRatioToPiElectronsMC");
  hRatioToPiElectronsMC->SetMarkerStyle(24);
  TH1F* hRatioToPiMuonsMC = (TH1F*)hRatioToPiMuons->Clone("hRatioToPiMuonsMC");
  hRatioToPiMuonsMC->SetMarkerStyle(24);
  hRatioToPiMuonsMC->SetLineColor(getLineColor(kMu));
  hRatioToPiMuonsMC->SetMarkerColor(getLineColor(kMu));
  TH1F* hRatioToPiKaonsMC = (TH1F*)hRatioToPiKaons->Clone("hRatioToPiKaonsMC");
  hRatioToPiKaonsMC->SetMarkerStyle(24);
  TH1F* hRatioToPiProtonsMC = (TH1F*)hRatioToPiProtons->Clone("hRatioToPiProtonsMC");
  hRatioToPiProtonsMC->SetMarkerStyle(24);
  
  // Reduced Chi^2 of fits vs. pT for all Delta_Species
  TH2F* hReducedChiSquarePt = 0x0;
  if (mode == kPMpT)
    hReducedChiSquarePt = new TH2F("hReducedChiSquarePt", "e", nPtBins, binsPt, 4, 0, 4);
  else {
    const TArrayD* histBins = hPIDdata->GetAxis(axisForMode)->GetXbins();
    if (histBins->fN == 0)
      hReducedChiSquarePt = new TH2F(Form("hReducedChiSquare%s", modeShortName[mode].Data()), "e",
                                     hPIDdata->GetAxis(axisForMode)->GetNbins(), hPIDdata->GetAxis(axisForMode)->GetXmin(),
                                     hPIDdata->GetAxis(axisForMode)->GetXmax(), 4, 0, 4);
    else
      hReducedChiSquarePt = new TH2F(Form("hReducedChiSquare%s", modeShortName[mode].Data()), "e",
                                     hPIDdata->GetAxis(axisForMode)->GetNbins(), histBins->fArray, 4, 0, 4);
  }
  
  hReducedChiSquarePt->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
  hReducedChiSquarePt->GetYaxis()->SetTitle("Delta_{Species}");
  hReducedChiSquarePt->GetYaxis()->SetBinLabel(1, "e");
  hReducedChiSquarePt->GetYaxis()->SetBinLabel(2, "K");
  hReducedChiSquarePt->GetYaxis()->SetBinLabel(3, "#pi");
  hReducedChiSquarePt->GetYaxis()->SetBinLabel(4, "p");
  hReducedChiSquarePt->SetMarkerColor(kRed);
  hReducedChiSquarePt->SetMarkerStyle(20);
  hReducedChiSquarePt->SetStats(kFALSE);
  
  // Obtain MC information about particle yields
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1); // Do not count each particle more than once
  TH2D* hMCdata = (TH2D*)hPIDdata->Projection(kPidMCpid, axisForMode, "e");
  hMCdata->SetName("hMCdata");
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1); // Reset range
  
  // Free memory (strangely only works like this...)
  hPIDdata->Reset();

  
  
  // Extract the MC truth generated primary yields
  f->cd();
  THnSparse* hMCgeneratedYieldsPrimaries = isMCdataSet ? dynamic_cast<THnSparse*>(histList->FindObject("fhMCgeneratedYieldsPrimaries"))
                                                       : 0x0;
  saveInterF->cd();
  TH1D* hMCgenYieldsPrimSpecies[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    hMCgenYieldsPrimSpecies[i] = 0x0;
  
  if (hMCgeneratedYieldsPrimaries) {
    // If desired, rebin considered axis
    if (rebin > 1) {
      const Int_t nDimensions = hMCgeneratedYieldsPrimaries->GetNdimensions();
      Int_t rebinFactor[nDimensions];
      
      for (Int_t dim = 0; dim < nDimensions; dim++) {
        if (dim == axisGenYieldForMode && rebin > 1)
          rebinFactor[dim] = rebin;
        else
          rebinFactor[dim] = 1;
      }
      
      THnSparse* temp = hMCgeneratedYieldsPrimaries->Rebin(&rebinFactor[0]);
      hMCgeneratedYieldsPrimaries->Reset();
      hMCgeneratedYieldsPrimaries = temp;
    }
    
    // Set proper errors, if not yet calculated
    if (!hMCgeneratedYieldsPrimaries->GetCalculateErrors()) {
      std::cout << "Re-calculating errors of " << hMCgeneratedYieldsPrimaries->GetName() << "..." << std::endl;
      
      hMCgeneratedYieldsPrimaries->Sumw2();
      
      Long64_t nBinsTHnSparseGenYield = hMCgeneratedYieldsPrimaries->GetNbins();
      Double_t binContentGenYield = 0;
      for (Long64_t bin = 0; bin < nBinsTHnSparseGenYield; bin++) {
        binContentGenYield = hMCgeneratedYieldsPrimaries->GetBinContent(bin);
        hMCgeneratedYieldsPrimaries->SetBinError(bin, TMath::Sqrt(binContentGenYield));
      }
    }
      
    if (restrictJetPtAxis) 
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
    
    // Not needed anymore, since values have been adapted: if (restrictCentralityAxis)
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    
    
    if (restrictCharge) {
      const Int_t indexChargeAxisGenYield = GetAxisByTitle(hMCgeneratedYieldsPrimaries, "Charge (e_{0})");
      if (indexChargeAxisGenYield < 0) {
        std::cout << "Error: Charge axis not found for gen yield histogram!" << std::endl;
        return -1;
      }
  
      Int_t lowerChargeBinLimitGenYield = -1;
      Int_t upperChargeBinLimitGenYield = -2;
      Double_t actualLowerChargeGenYield = -999;
      Double_t actualUpperChargeGenYield = -999;
  
      // Add subtract a very small number to avoid problems with values right on the border between to bins
      if (chargeMode == kNegCharge) {
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindFixBin(-1. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindFixBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindFixBin(0. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindFixBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimitGenYield <= upperChargeBinLimitGenYield && lowerChargeBinLimitGenYield >= 0
          && upperChargeBinLimitGenYield <= hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetNbins() + 1) {
        actualLowerChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinLowEdge(lowerChargeBinLimitGenYield);
        actualUpperChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinUpEdge(upperChargeBinLimitGenYield);
        
        if (TMath::Abs(actualLowerChargeGenYield - actualLowerChargeData) > 1e-4 ||
            TMath::Abs(actualUpperChargeGenYield - actualUpperChargeData) > 1e-4) {
          std::cout << std::endl;
          std::cout << "Error: Charge range gen yield: " << actualLowerChargeGenYield << " - " << actualUpperChargeGenYield
                    << std::endl << "differs from that of data: " << actualLowerChargeData << " - " << actualUpperChargeData
                    << std::endl;
          return -1;
        }
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested charge range (gen yield) out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
      
      hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->SetRange(lowerChargeBinLimitGenYield, upperChargeBinLimitGenYield);
    }
    
    for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(MCid + 1, MCid + 1);
      
      hMCgenYieldsPrimSpecies[MCid] = hMCgeneratedYieldsPrimaries->Projection(axisGenYieldForMode, "e");
      hMCgenYieldsPrimSpecies[MCid]->SetName(Form("hMCgenYieldsPrimSpecies_%s", AliPID::ParticleShortName(MCid)));
      hMCgenYieldsPrimSpecies[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      
       
      // Choose the same binning as for the fitted yields, i.e. rebin the histogram (Rebin will create a clone!).
      // This is only needed for pT! z and xi are already treated by the re-binning of the THnSparse itself above!
      if (axisGenYieldForMode == kPidGenYieldPt) {
        TH1D* temp = (TH1D*)hMCgenYieldsPrimSpecies[MCid]->Rebin(nPtBins, hMCgenYieldsPrimSpecies[MCid]->GetName(), binsPt);
        // Delete the old binned histo and take the new binned one
        delete hMCgenYieldsPrimSpecies[MCid];
        hMCgenYieldsPrimSpecies[MCid] = temp;
      }
      
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(0, -1);
    }
    
    // Free memory (strangely only works like this...)
    hMCgeneratedYieldsPrimaries->Reset();
  }
  
  
  // Get expected shapes for pT bins
  TString Ytitle = "";
  
  // Array index 0 as unused dummy
  TH2D* hGenDelta[6][6]; // DeltaSpecies (first index) for species (second index)
  TH2D* hGenDeltaMCid[6][6]; // DeltaSpecies (first index) for species (second index)
  TH2D* hGenDeltaTOFpi[6][6]; // DeltaSpecies (first index) for species (second index) for the TOF pion bin
  
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 6; j++) {
      hGenDelta[i][j] = 0x0;
      hGenDeltaMCid[i][j] = 0x0;
      hGenDeltaTOFpi[i][j] = 0x0;
    }
  }
  
  THnSparse* current = 0x0;
  f->cd();
  THnSparse* hGenEl = dynamic_cast<THnSparse*>(histList->FindObject("hGenEl"));
  if (!hGenEl)  {
    std::cout << "Failed to load expected dEdx signal shape for: Electrons!" << std::endl;
    return -1;
  }
  
  THnSparse* hGenKa = dynamic_cast<THnSparse*>(histList->FindObject("hGenKa"));
  if (!hGenKa)  {
    std::cout << "Failed to load expected dEdx signal shape for: Kaons!" << std::endl;
    return -1;
  }
  
  THnSparse* hGenPi = dynamic_cast<THnSparse*>(histList->FindObject("hGenPi"));
  if (!hGenPi)  {
    std::cout << "Failed to load expected dEdx signal shape for: Pions!" << std::endl;
    return -1;
  }
  
  THnSparse* hGenMu = dynamic_cast<THnSparse*>(histList->FindObject("hGenMu"));
  if (!hGenMu)  {
    if (muonFractionHandling != kNoMuons)
      std::cout << "Failed to load expected dEdx signal shape for: Muons! Treated muons as pions in the following." << std::endl;
    takeIntoAccountMuons = kFALSE; 
  }
  
  THnSparse* hGenPr = dynamic_cast<THnSparse*>(histList->FindObject("hGenPr"));
  if (!hGenPr)  {
    std::cout << "Failed to load expected dEdx signal shape for: Protons!" << std::endl;
    return -1;
  }
  saveInterF->cd();
  for (Int_t MCid = kEl; MCid <= kPr; MCid++) {
    if (MCid == kEl)
      current = hGenEl;
    else if (MCid == kKa)
      current = hGenKa;
    else if (MCid == kMu) {
      if (takeIntoAccountMuons)
        current = hGenMu;
      else
       continue; // No histo for muons in this case
    }
    else if (MCid == kPi)
      current = hGenPi;
    else if (MCid == kPr)
      current = hGenPr;
    else
      break;
    
    // If desired, rebin considered axis
    if (rebin > 1 || rebinDeltaPrime > 1) {
      const Int_t nDimensions = current->GetNdimensions();
      Int_t rebinFactor[nDimensions];
      
      for (Int_t dim = 0; dim < nDimensions; dim++) {
        if (dim == axisGenForMode)
          rebinFactor[dim] = rebin;
        else if (dim == kPidGenDeltaPrime && rebinDeltaPrime > 1)
        rebinFactor[dim] = rebinDeltaPrime;
        else
          rebinFactor[dim] = 1;
      }
      
      THnSparse* temp = current->Rebin(&rebinFactor[0]);
      current->Reset();
      current = temp;
    }
    
    // Set proper errors, if not yet calculated
    if (!current->GetCalculateErrors()) {
      std::cout << "Re-calculating errors of " << current->GetName() << "..." << std::endl;
      
      current->Sumw2();
      
      Long64_t nBinsTHnSparseGen = current->GetNbins();
      Double_t binContentGen = 0;
      for (Long64_t bin = 0; bin < nBinsTHnSparseGen; bin++) {
        binContentGen = current->GetBinContent(bin);
        current->SetBinError(bin, TMath::Sqrt(binContentGen));
      }
    }
    
    // If desired, restrict centrality range
    // Not needed anymore, since values have been adapted: if (restrictCentralityAxis)
      current->GetAxis(kPidGenCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    
    // If desired, restrict jet pT range
    if (restrictJetPtAxis) {
      current->GetAxis(kPidGenJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
    }
    
    // If desired, restrict charge range
    if (restrictCharge) {
      const Int_t indexChargeAxisGen = GetAxisByTitle(current, "Charge (e_{0})");
      if (indexChargeAxisGen < 0) {
        std::cout << "Error: Charge axis not found for gen histogram!" << std::endl;
        return -1;
      }
  
      Int_t lowerChargeBinLimitGen = -1;
      Int_t upperChargeBinLimitGen = -2;
      Double_t actualLowerChargeGen = -999;
      Double_t actualUpperChargeGen = -999;
  
      // Add subtract a very small number to avoid problems with values right on the border between to bins
      if (chargeMode == kNegCharge) {
        lowerChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindFixBin(-1. + 0.001);
        upperChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindFixBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindFixBin(0. + 0.001);
        upperChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindFixBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimitGen <= upperChargeBinLimitGen && lowerChargeBinLimitGen >= 0
          && upperChargeBinLimitGen <= current->GetAxis(indexChargeAxisGen)->GetNbins() + 1) {
        actualLowerChargeGen = current->GetAxis(indexChargeAxisGen)->GetBinLowEdge(lowerChargeBinLimitGen);
        actualUpperChargeGen = current->GetAxis(indexChargeAxisGen)->GetBinUpEdge(upperChargeBinLimitGen);
        
        if (TMath::Abs(actualLowerChargeGen - actualLowerChargeData) > 1e-4 ||
            TMath::Abs(actualUpperChargeGen - actualUpperChargeData) > 1e-4) {
          std::cout << std::endl;
          std::cout << "Error: Charge range gen: " << actualLowerChargeGen << " - " << actualUpperChargeGen
                    << std::endl << "differs from that of data: " << actualLowerChargeData << " - " << actualUpperChargeData
                    << std::endl;
          return -1;
        }
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested charge range (gen) out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
      
      current->GetAxis(indexChargeAxisGen)->SetRange(lowerChargeBinLimitGen, upperChargeBinLimitGen);
    }
    
    const Int_t indexTOFinfoAxisGen = GetAxisByTitle(current, "TOF PID Info");
    if (indexTOFinfoAxisGen < 0 && applyTOFpatching) {
      std::cout << "Error: TOF PID info axis not found for data histogram!" << std::endl;
      return -1;
    }

    for (Int_t selectBin = 1; selectBin <= 4; selectBin++)  {
      Int_t selectMCid = (selectBin >= 3) ? (selectBin+1) : selectBin;

      current->GetAxis(kPidGenSelectSpecies)->SetRange(selectBin, selectBin);
      
      Ytitle = Form("#Delta%s = dE/dx %s <dE/dx>_{%s}",
                    useDeltaPrime ? Form("'_{#lower[-0.5]{%s}}", partShortName[selectMCid - 1].Data()) 
                                  : Form("_{%s}", partShortName[selectMCid - 1].Data()),
                    useDeltaPrime ? "/" : "-", partShortName[selectMCid - 1].Data());
      
      TH2* hGenCurrent = 0x0;
      
      if (applyTOFpatching) {
        // In case of TOF patching: Restrict TOF info axis to no TOF region.
        // This is necessary to get proper templates. Extreme example: Assuming that at low pT only protons
        // can reach the TOF (for whatever reason), but no pions, kaons etc. Then the templates would be wrong, if one
        // abandons the Pre-PID and produces templates for every species. Because there will be the contribution from the
        // TOF protons which are not in the data points
        current->GetAxis(indexTOFinfoAxisGen)->SetRange(current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kNoTOFinfo),
                                                        current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kNoTOFpid));
      }
      if (!useIdentifiedGeneratedSpectra)   {
        hGenDelta[selectMCid][MCid] = current->Projection(genAxis, axisGenForMode, "e");
        hGenDelta[selectMCid][MCid]->SetName(Form("hGenDelta%sFor%s", partShortName[selectMCid - 1].Data(),
                                                  partShortName[MCid - 1].Data()));
        hGenCurrent = hGenDelta[selectMCid][MCid];
      }
      else  {
        current->GetAxis(kPidGenMCpid)->SetRange(MCid, MCid);
        hGenDeltaMCid[selectMCid][MCid] = current->Projection(genAxis, axisGenForMode, "e");
        hGenDeltaMCid[selectMCid][MCid]->SetName(Form("hGenDelta%sForMCid%s", partShortName[selectMCid - 1].Data(), 
                                                      partShortName[MCid - 1].Data()));
        
        hGenCurrent = hGenDeltaMCid[selectMCid][MCid];
        current->GetAxis(kPidGenMCpid)->SetRange(0, -1);
      }
      
      hGenCurrent->GetYaxis()->SetTitle(Ytitle.Data());
      hGenCurrent->SetLineColor(getLineColor(MCid));
      hGenCurrent->SetMarkerColor(getLineColor(MCid));
      hGenCurrent->SetLineWidth(2);
      hGenCurrent->SetLineStyle(2);
      hGenCurrent->SetMarkerStyle(20);
      hGenCurrent->GetXaxis()->SetTitleOffset(1.0);
      
      if (applyTOFpatching) {
        // Set to TOF pion range
        current->GetAxis(indexTOFinfoAxisGen)->SetRange(current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kTOFpion),
                                                        current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kTOFpion));
        
        // In case of TOF patching: Also obtain generated signals for the TOF pion bin
         if (!useIdentifiedGeneratedSpectra)   {
          hGenDeltaTOFpi[selectMCid][MCid] = current->Projection(genAxis, axisGenForMode, "e");
          hGenDeltaTOFpi[selectMCid][MCid]->SetName(Form("hGenDelta%sFor%sTOFpi", partShortName[selectMCid - 1].Data(),
                                                         partShortName[MCid - 1].Data()));
        }
        else  {
          current->GetAxis(kPidGenMCpid)->SetRange(MCid, MCid);
          hGenDeltaTOFpi[selectMCid][MCid] = current->Projection(genAxis, axisGenForMode, "e");
          hGenDeltaTOFpi[selectMCid][MCid]->SetName(Form("hGenDelta%sForMCid%sTOFpi", partShortName[selectMCid - 1].Data(), 
                                                         partShortName[MCid - 1].Data()));
          
          current->GetAxis(kPidGenMCpid)->SetRange(0, -1);
        }
        
        hGenDeltaTOFpi[selectMCid][MCid]->GetYaxis()->SetTitle(Ytitle.Data());
        hGenDeltaTOFpi[selectMCid][MCid]->SetLineColor(getLineColor(MCid));
        hGenDeltaTOFpi[selectMCid][MCid]->SetMarkerColor(getLineColor(MCid));
        hGenDeltaTOFpi[selectMCid][MCid]->SetLineWidth(2);
        hGenDeltaTOFpi[selectMCid][MCid]->SetLineStyle(2);
        hGenDeltaTOFpi[selectMCid][MCid]->SetMarkerStyle(20);
        hGenDeltaTOFpi[selectMCid][MCid]->GetXaxis()->SetTitleOffset(1.0);
        
        // Set back to original range (i.e. no TOF)
        current->GetAxis(indexTOFinfoAxisGen)->SetRange(current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kNoTOFinfo),
                                                        current->GetAxis(indexTOFinfoAxisGen)->FindFixBin(kNoTOFpid));
      }
    }
    
    current->GetAxis(kPidGenSelectSpecies)->SetRange(0, -1);
    
    // Free memory (strangely only works like this...)
    current->Reset();
  }
  
  // Free a lot of memory for the following procedure. Histogram is not needed anymore (only its projections)
  // -> Seems not to help!
  //delete f;
  f->Close();
  
  
  // Save intermediate results
  saveInterF->cd();
  
  if (hMCdata)
    hMCdata->Write();
  
  if (hMCdataTOF)
    hMCdataTOF->Write();
  
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 6; j++) {
      if (hGenDelta[i][j])
        hGenDelta[i][j]->Write();
        
      if (hGenDeltaMCid[i][j])
        hGenDeltaMCid[i][j]->Write();
      
      if (hGenDeltaTOFpi[i][j])
        hGenDeltaTOFpi[i][j]->Write();
    }
  }
  
  for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < hFractionPions->GetXaxis()->GetNbins(); slice++) {   
    if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
      continue; 
    
    if(hDeltaPi[slice])  
      hDeltaPi[slice]->Write();
    if(hDeltaEl[slice])
      hDeltaEl[slice]->Write();
    if(hDeltaKa[slice])
      hDeltaKa[slice]->Write();
    if(hDeltaPr[slice])
      hDeltaPr[slice]->Write();
      
    //if (plotIdentifiedSpectra)  {
      for (Int_t species = 0; species < 5; species++) {
        if (hDeltaElMC[slice][species])
          hDeltaElMC[slice][species]->Write();
        if (hDeltaKaMC[slice][species])
          hDeltaKaMC[slice][species]->Write();
        if (hDeltaPiMC[slice][species])
          hDeltaPiMC[slice][species]->Write();
        if (hDeltaPrMC[slice][species])
          hDeltaPrMC[slice][species]->Write(); 
      }
    //}   
    
    if(hDeltaPiTOFpi[slice])  
      hDeltaPiTOFpi[slice]->Write();
    if(hDeltaElTOFpi[slice])
      hDeltaElTOFpi[slice]->Write();
    if(hDeltaKaTOFpi[slice])
      hDeltaKaTOFpi[slice]->Write();
    if(hDeltaPrTOFpi[slice])
      hDeltaPrTOFpi[slice]->Write();
      
    //if (plotIdentifiedSpectra)  {
      for (Int_t species = 0; species < 5; species++) {
        if (hDeltaElMCTOFpi[slice][species])
          hDeltaElMCTOFpi[slice][species]->Write();
        if (hDeltaKaMCTOFpi[slice][species])
         hDeltaKaMCTOFpi[slice][species]->Write();
        if (hDeltaPiMCTOFpi[slice][species])
         hDeltaPiMCTOFpi[slice][species]->Write();
        if (hDeltaPrMCTOFpi[slice][species])
         hDeltaPrMCTOFpi[slice][species]->Write(); 
      }
    //}   
  }
  
  if (hYieldTOFOrigBinningPions)
    hYieldTOFOrigBinningPions->Write();
  
  if (hYieldTOFOrigBinningKaons)
    hYieldTOFOrigBinningKaons->Write();
  
  if (hYieldTOFOrigBinningProtons)
    hYieldTOFOrigBinningProtons->Write();
  
  // File may not be closed because the projections are needed in the following!
  //saveInterF->Close();
    
  // Create histo for TOF electrons and reset everything to zero. If TOF patching is applied, it can be used to take into account
  // the electron contamination
  TH1D* hYieldTOFOrigBinningElectrons = 0x0;
  
  if (applyTOFpatching) {
    hYieldTOFOrigBinningElectrons = new TH1D(*hYieldTOFOrigBinningPions);
    hYieldTOFOrigBinningElectrons->Reset();
    
    hYieldTOFOrigBinningElectrons->SetName("hYieldTOFOrigBinningElectrons");
    hYieldTOFOrigBinningElectrons->SetTitle(AliPID::ParticleLatexName(AliPID::kElectron));
    hYieldTOFOrigBinningElectrons->SetLineColor(getLineColor(kEl));
    hYieldTOFOrigBinningElectrons->SetMarkerColor(getLineColor(kEl));
  }
  
  // Save some first results for the final output
  TString saveFName = fileName;
  saveFName = Form("%s_results_%s__%s_%d_reg%d_regFac%.2f_%s%s%s%s%s%s.root", saveFName.ReplaceAll(".root", "").Data(), 
                   useLogLikelihood ? (useWeightsForLogLikelihood ? "weightedLLFit" : "LLFit") : "ChiSquareFit",
                   modeShortName[mode].Data(), fitMethod, regularisation, regularisationFactor,
                   muonFractionHandlingShortName[muonFractionHandlingParameter].Data(),
                   useIdentifiedGeneratedSpectra ? "_idSpectra" : "",
                   restrictCentralityAxis ? Form(centralityHasDecimalsPlaces ? "_centrality%.0fem2_%.0fem2" : "_centrality%.0f_%.0f", 
                                                 centralityHasDecimalsPlaces ? 100.*actualLowerCentrality : actualLowerCentrality,
                                                 centralityHasDecimalsPlaces ? 100.*actualUpperCentrality : actualUpperCentrality) : "",
                   restrictJetPtAxis ? Form("_jetPt%.1f_%.1f", actualLowerJetPt, actualUpperJetPt) : "",
                   chargeString.Data(), applyTOFpatching ? "_TPConly" : "");
  TFile *saveF = TFile::Open(saveFName.Data(), "RECREATE");
  saveF->cd();
  
  if (hFractionElectrons)
      hFractionElectrons->Write(0, TObject::kWriteDelete);
      
  if (hFractionKaons)
    hFractionKaons->Write(0, TObject::kWriteDelete);
    
  if (hFractionPions)
    hFractionPions->Write(0, TObject::kWriteDelete);
    
  if (hFractionProtons)
    hFractionProtons->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuons)
    hFractionMuons->Write(0, TObject::kWriteDelete);
  
  if (hFractionSummed)
    hFractionSummed->Write(0, TObject::kWriteDelete);
  
  if (hFractionElectronsMC)
    hFractionElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionKaonsMC)
    hFractionKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionPionsMC)
    hFractionPionsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsMC)
    hFractionMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionProtonsMC)
    hFractionProtonsMC->Write(0, TObject::kWriteDelete);

  
  if (hYieldElectrons)
    hYieldElectrons->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaons)
    hYieldKaons->Write(0, TObject::kWriteDelete);
  
  if (hYieldPions)
    hYieldPions->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtons)
    hYieldProtons->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuons)
    hYieldMuons->Write(0, TObject::kWriteDelete);
  
  if (hYieldElectronsMC)
    hYieldElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsMC)
    hYieldMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsMC)
    hYieldKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsMC)
    hYieldPionsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsMC)
    hYieldProtonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldSummedMC)
    hYieldSummedMC->Write(0, TObject::kWriteDelete);
  
  
  if (hRatioToPiElectrons)
    hRatioToPiElectrons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiMuons)
    hRatioToPiMuons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiKaons)
    hRatioToPiKaons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiProtons)
    hRatioToPiProtons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiElectronsMC)
    hRatioToPiElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiMuonsMC)
    hRatioToPiMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiKaonsMC)
    hRatioToPiKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiProtonsMC)
    hRatioToPiProtonsMC->Write(0, TObject::kWriteDelete);
  
  // Dummy histo to create generic legend entries from
  TH1D* hMCmuonsAndPionsDummy = 0x0;
  if (plotIdentifiedSpectra && firstValidSlice >= 0) {
    hMCmuonsAndPionsDummy = new TH1D(*hDeltaPiMC[firstValidSlice][kPi]);
    hMCmuonsAndPionsDummy->SetLineColor(getLineColor(kMuPlusPi));
    hMCmuonsAndPionsDummy->SetMarkerColor(getLineColor(kMuPlusPi));
    hMCmuonsAndPionsDummy->SetName("hMCmuonsAndPionsDummy");
  }
  
  
  // No different muon treatment for z and xi (would require special functions - not worth the work!)
  muonFractionThresholdForFitting = 0.;//OLD 0.295;
  muonFractionThresholdBinForFitting = (mode == kPMpT) ? FindMomentumBin(binsPt, muonFractionThresholdForFitting) : -1;
  
  if (mode == kPMpT) {
    // For pT, move threshold slightly to the left to also fix the bin if the threshold is right on the lower edge
    electronFractionThresholdForFitting = electronFractionThresholdForFitting - 1e-4;
    // Move one bin to the right, such that the fitting and fixing is done from that bin on
    electronFractionThresholdBinForFitting = FindMomentumBin(binsPt, electronFractionThresholdForFitting) + 1;
    if (electronFractionThresholdBinForFitting < 1) {
      printf("Fatal error: electronFractionThresholdBinForFitting < 1!\n");
      exit(-1);
    }
    // To make the fit work during regularisation, set the threshold to the lower end of this bin
    electronFractionThresholdForFitting = binsPt[TMath::Min(nPtBins + 1, electronFractionThresholdBinForFitting - 1)];// -1 because lower threshold is wanted here!
    
    lastPtForCallOfGetElectronFraction = pHigh + 10.; // Make sure that this value is higher than in any call during the fit
    
    fElectronFraction = new TF1("fElectronFraction", Form("[0]+(x<%f)*[1]*(x-%f)", electronFractionThresholdForFitting, 
                                                          electronFractionThresholdForFitting), 
                                pLow, pHigh);
  }
  else if (mode == kPMz) {
    // Use threshold in pT to get threshold in z -> Turns out that it makes sense to use the middle of the jet pT range
    // as an estimate for z
    
    Double_t effectiveJetPt = (actualUpperJetPt + actualLowerJetPt) / 2.;
    // For z, move threshold slightly to the left to also fix the bin if the threshold is right on the lower edge
    Double_t effectiveThreshold = electronFractionThresholdForFitting / effectiveJetPt - 1e-4;
    // Move one bin to the right, such that the fitting and fixing is done from that bin on
    electronFractionThresholdBinForFitting = hYieldPt->GetXaxis()->FindFixBin(effectiveThreshold) + 1;
    // To make the fit work during regularisation, set the threshold to the lower end of this bin
    if (electronFractionThresholdBinForFitting > hYieldPt->GetNbinsX())
      electronFractionThresholdForFitting = effectiveThreshold = hYieldPt->GetXaxis()->GetBinUpEdge(hYieldPt->GetNbinsX());
    else
      electronFractionThresholdForFitting = effectiveThreshold = hYieldPt->GetXaxis()->GetBinLowEdge(electronFractionThresholdBinForFitting);
    
    // For the lower threshold of the fitting (in case of xi the upper actually), it turns out to yield reasonable results
    // if the lowest jet pT is chosen.
    effectiveJetPt = actualLowerJetPt;
    lowFittingBoundElectronFraction = lowFittingBoundElectronFraction / effectiveJetPt;
    
    lastPtForCallOfGetElectronFraction = 100.; // Make sure that this value is higher than in any call during the fit
    
    fElectronFraction = new TF1("fElectronFraction", Form("[0]+(x<%f)*[1]*(x-%f)", electronFractionThresholdForFitting,
                                                          electronFractionThresholdForFitting), 
                                0, 1.1);
  }
  else if (mode == kPMxi) {
    // Likewise as for z
    Double_t effectiveJetPt = (actualUpperJetPt + actualLowerJetPt) / 2.;
    // For xi, move threshold slightly to the right to also fix the bin if the threshold is right on the upper edge
    Double_t effectiveThreshold = TMath::Log(effectiveJetPt / electronFractionThresholdForFitting) + 1e-4;
    // Move one bin to the left, such that the fitting and fixing is done from that bin on
    electronFractionThresholdBinForFitting = hYieldPt->GetXaxis()->FindFixBin(effectiveThreshold) - 1;
    // To make the fit work during regularisation, set the threshold to the upper end of this bin
    if (electronFractionThresholdBinForFitting < 1)
      electronFractionThresholdForFitting = effectiveThreshold = hYieldPt->GetXaxis()->GetBinLowEdge(1);
    else
      electronFractionThresholdForFitting = effectiveThreshold = hYieldPt->GetXaxis()->GetBinUpEdge(electronFractionThresholdBinForFitting);
    
    effectiveJetPt = actualLowerJetPt;
    lowFittingBoundElectronFraction = TMath::Log(effectiveJetPt / lowFittingBoundElectronFraction);
    
    lastPtForCallOfGetElectronFraction = -999.; // Make sure that this value is smaller than in any call during the fit
    
    // NOTE: Fix fraction is BELOW the threshold for xi!
    fElectronFraction = new TF1("fElectronFraction", Form("[0]+(x>%f)*[1]*(x-%f)", electronFractionThresholdForFitting, 
                                                          electronFractionThresholdForFitting), 
                                0, 9.);
  }
  else {
    // No special treatment of el for R and jT (seems to be the best option), i.e., set arbitrary high threshold which is never reached
    fElectronFraction = 0x0;
    electronFractionThresholdBinForFitting = 99999;
  }
  
  if (fElectronFraction)
    fElectronFraction->SetParameters(0.01, 0.0);

  
  
  // For special templates at high dEdx, calculate the thresholds based on the thresholds for the track pT.
  // To be sure that clean PID is possible, take the most extreme jet pT in the considered bin to calculate the threshold.
  if (mode == kPMz) {
    mostProbPIDthresholdZ_ka = mostProbPIDthresholdPt_ka / actualUpperJetPt;
    mostProbPIDthresholdZ_pr = mostProbPIDthresholdPt_pr / actualUpperJetPt;
  }
  else if (mode == kPMxi) {
    mostProbPIDthresholdXi_ka = TMath::Log(actualUpperJetPt / mostProbPIDthresholdPt_ka);
    mostProbPIDthresholdXi_pr = TMath::Log(actualUpperJetPt / mostProbPIDthresholdPt_pr);
  }
  
  
  TString speciesLabel[4] = {"El", "Ka", "Pi", "Pr" };
  
  if (applyTOFpatching) {
    FitElContaminationForTOFpions(hFractionPions->GetXaxis(), hGenDeltaTOFpi, hMCmuonsAndPionsDummy, &hDeltaPiTOFpi[0], &hDeltaElTOFpi[0], &hDeltaKaTOFpi[0], &hDeltaPrTOFpi[0], hDeltaPiMCTOFpi, hDeltaElMCTOFpi, hDeltaKaMCTOFpi, hDeltaPrMCTOFpi, mode, useLogLikelihood, useWeightsForLogLikelihood, plotIdentifiedSpectra, pSliceLow, pSliceHigh, numSlices, restrictJetPtAxis, actualUpperJetPt, xLow, xUp, nBins,  speciesLabel, saveF, hYieldTOFOrigBinningPions, hYieldTOFOrigBinningElectrons);
    
    TH1D* hYieldTOFOrigBinningElectronsContaminationCleaned = new TH1D(*hYieldTOFOrigBinningElectrons);
    hYieldTOFOrigBinningElectronsContaminationCleaned->SetName("hYieldTOFOrigBinningElectronsContaminationCleaned");
    
    TH1D* hYieldTOFOrigBinningPionsContaminationCleaned = new TH1D(*hYieldTOFOrigBinningPions);
    hYieldTOFOrigBinningPionsContaminationCleaned->SetName("hYieldTOFOrigBinningPionsContaminationCleaned");
    
    saveInterF->cd();
    hYieldTOFOrigBinningElectronsContaminationCleaned->Write();
    hYieldTOFOrigBinningPionsContaminationCleaned->Write();
    saveF->cd();
  }
  
  
  // In case of regularisation, the actual number of x bins and the (for pT: logs of their) bin centres are required
  Int_t numXBins = 0;
  Double_t* xBinCentres = 0x0;
  Double_t* xBinStatisticalWeight = 0x0;
  Double_t* xBinStatisticalWeightError = 0x0;
  
  Double_t* xBinStatisticalWeightTOF = 0x0;
  Double_t* xBinStatisticalWeightErrorTOF = 0x0;
  Double_t* xTOFyield = 0x0;
  
  // Set the number of parameters per x bin:
  // Regularisation only implemented for simultaneous fit.
  const Int_t numParamsPerXbin = AliPID::kSPECIES + 1; // Fractions of each species + total yield in x bin
  
  // Construct the array of all the parameters that are to be regularised, i.e. only the FREE fractions
  // and NOT the total yields or the x bin
  Int_t nParToRegulariseSimultaneousFit = 0;
  Int_t* indexParametersToRegularise = 0x0;
  Int_t* lastNotFixedIndexOfParameters = 0x0;
  
  if (regularisation > 0) {
    Int_t xBinIndexTemp = 0;
    Int_t internalParIndexTemp = 0;
    
    // Switch off regularisation for very low pT (or low z and high xi) for kaons and protons.
    // Reason: Due to "kink" in efficiency for K and p (and efficiency close to zero) the assumption of
    // smooth raw fractions is not a good approximation any more. Anyway, statistics is sufficient in this
    // region and there is no crossing (e and pi have still quite different templates). This means there are
    // no outlier from this side expected, but only real jumps.
    // Also, there is a jump for the pion efficiency vs. z in the first bin (whereas it is smooth vs. pT!).
    
    // Find the bin with the threshold. This bin and the next (in case of xi, previous) one are switched off for regularisation.
    // Reason for also next bin: Don't want to bias this neighbour, since the other bin contributes to interpolation
    // (mainly important for z, where only one bin would be affected otherwise).
    
    Int_t binEffectiveRegThreshold_pi = -1;
    Int_t binEffectiveRegThreshold_pr = -1;
    Int_t binEffectiveRegThreshold_ka = -1;
    if (mode == kPMpT) {
      // No threshold for pi, since smooth vs pT
      binEffectiveRegThreshold_pr = hFractionPions->GetXaxis()->FindFixBin(regOffThresholdPt_pr) + 1;
      binEffectiveRegThreshold_ka = hFractionPions->GetXaxis()->FindFixBin(regOffThresholdPt_ka) + 1;
    }
    else if (mode == kPMz) {
      /*OLD// For z (and xi) it turns out that the threshold can be estimated from the pT threshold best if
      // the middle of the jet pT bin is used for the jet pT value
      Double_t effectiveJetPt = (actualUpperJetPt + actualLowerJetPt) / 2.;*/
      
      
      // For z (and xi) it turns out that the threshold can be estimated from the pT threshold best if
      // the lower edge of the jet pT bin is used for the jet pT value (ensures that all tracks have passed
      // the pT threshold)
      Double_t effectiveJetPt = actualLowerJetPt;
      Double_t regThresholdZ_pr = regOffThresholdPt_pr / effectiveJetPt;
      Double_t regThresholdZ_ka = regOffThresholdPt_ka / effectiveJetPt;
      
      // For pions: There is always a jump between first and second bin, independent of jetPt. Take out the first two bins.
      binEffectiveRegThreshold_pi = 1 + 1;
      binEffectiveRegThreshold_pr = hFractionPions->GetXaxis()->FindFixBin(regThresholdZ_pr) + 1;
      binEffectiveRegThreshold_ka = hFractionPions->GetXaxis()->FindFixBin(regThresholdZ_ka) + 1;
    }
    else if (mode == kPMxi) {
      //OLD Double_t effectiveJetPt = (actualUpperJetPt + actualLowerJetPt) / 2.;

      // Same comment as for z
      Double_t effectiveJetPt = actualLowerJetPt;
      // No threshold for pi, since smooth vs xi
      Double_t regThresholdXi_pr = TMath::Log(effectiveJetPt / regOffThresholdPt_pr);
      Double_t regThresholdXi_ka = TMath::Log(effectiveJetPt / regOffThresholdPt_ka);
      
      // Minus(!) 1, since low pT part sits at high xi
      binEffectiveRegThreshold_pr = hFractionPions->GetXaxis()->FindFixBin(regThresholdXi_pr) - 1;
      binEffectiveRegThreshold_ka = hFractionPions->GetXaxis()->FindFixBin(regThresholdXi_ka) - 1;
    }
    
    // Loop twice over data: Determine the number of bins in the first iteration, allocate the memory and fill in the 2. iteration
    for (Int_t i = 0; i < 2; i++) {
      if (i == 1) {
        if (numXBins == 0) {
          printf("No bins for fitting! Exiting...\n");
          
          return -1;
        }
        if (nParToRegulariseSimultaneousFit == 0) {
          printf("No parameters to regularise! Exiting...\n");
          
          return -1;
        }
        
        xBinCentres = new Double_t[numXBins];
        xBinStatisticalWeight = new Double_t[numXBins];
        xBinStatisticalWeightError = new Double_t[numXBins];
        
        if (applyTOFpatching) {
          xBinStatisticalWeightTOF = new Double_t[numXBins];
          xBinStatisticalWeightErrorTOF = new Double_t[numXBins];
          xTOFyield = new Double_t[numParamsPerXbin * numXBins];
        }
        
        indexParametersToRegularise = new Int_t[nParToRegulariseSimultaneousFit];
        
        lastNotFixedIndexOfParameters = new Int_t[numParamsPerXbin];
        
        // Set last not fixed index of parameter to numXBins, i.e. a index larger than any existing index.
        // This will not restrict the parameter regularisation range. In the following, the range for electrons
        // and muons will be restricted
        for (Int_t iPar = 0; iPar < numParamsPerXbin; iPar++) 
          lastNotFixedIndexOfParameters[iPar] = numXBins;
      }
      
      // For xi, start at highest xi and move to smaller xi, i.e. highest xi get the smallest x index. This way all the code
      // can be kept as it is!
      const Bool_t isXiMode = (processingMode == kPMxi);
      const Bool_t isPtMode = (processingMode == kPMpT);
      const Int_t sliceStart = isXiMode ? hFractionPions->GetXaxis()->GetNbins() : 0;
      const Int_t sliceEnd   = isXiMode ? 0. : (isPtMode ? nPtBins : hFractionPions->GetXaxis()->GetNbins());
      
      for (Int_t slice = sliceStart; isXiMode ? slice >= sliceEnd : slice < sliceEnd; isXiMode ? slice-- : slice++) { 
        if (isPtMode && (slice < pSliceLow || slice > pSliceHigh))
          continue; 
        
        // There won't (actually: shouldn't) be tracks with a pT larger than the jet pT
        if (isPtMode && restrictJetPtAxis && binsPt[slice] >= actualUpperJetPt)
          continue;
        
        const Int_t pBinLowProjLimit = isPtMode ? hYieldPt->GetXaxis()->FindFixBin(binsPt[slice] + 1e-5)    : slice + 1;
        const Int_t pBinUpProjLimit  = isPtMode ? hYieldPt->GetXaxis()->FindFixBin(binsPt[slice + 1] - 1e-5) : slice + 1;
        
        // NOTE: In case of regularisation, only the simultaneous fit values will be used, i.e. totalYield and not allDeltaSpecies!
        
        Double_t totalYieldError = 0;
        const Double_t totalYield = hYieldPt->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, totalYieldError);
        
        // NOTE if there is no TPC yield, there should also be no TOF yield (and at high pT, where this must not be true,
        // the TOF efficiency is zero, so this is also not an issue)
        if (totalYield <= 0) 
          continue;
        
        if (i == 1) {
          if (mode == kPMpT)
            // Take the logarithm (base 10!) in case of pT
            xBinCentres[xBinIndexTemp] = TMath::Log10((binsPt[slice + 1] + binsPt[slice]) / 2.);
          else
            xBinCentres[xBinIndexTemp] = hYieldPt->GetXaxis()->GetBinCenter(slice + 1);
          
          // To switch off the regularisation (compare comment above), just set the corresponding indexParametersToRegularise to
          // a negative value.
          Int_t factorForSwitchingOffReg_pi = 1;
          Int_t factorForSwitchingOffReg_pr = 1;
          Int_t factorForSwitchingOffReg_ka = 1;
          
          if (mode == kPMpT || mode == kPMz) {
            if ((slice + 1) <= binEffectiveRegThreshold_pi)
              factorForSwitchingOffReg_pi = -1;
            if ((slice + 1) <= binEffectiveRegThreshold_pr)
              factorForSwitchingOffReg_pr = -1;
            if ((slice + 1) <= binEffectiveRegThreshold_ka)
              factorForSwitchingOffReg_ka = -1;
          }
          else if (mode == kPMxi) {
            if ((slice + 1) >= binEffectiveRegThreshold_pi)
              factorForSwitchingOffReg_pi = -1;
            if ((slice + 1) >= binEffectiveRegThreshold_pr)
              factorForSwitchingOffReg_pr = -1;
            if ((slice + 1) >= binEffectiveRegThreshold_ka)
              factorForSwitchingOffReg_ka = -1;
          }
         
          xBinStatisticalWeight[xBinIndexTemp] = totalYield;
          
          // NOTE: The total yield is a fact - a number w/o error. However, one assigns this error here to use it
          // to calculate the effective weighting for the weighted likelihood fit (and the error is only used for this).
          // So, it is more like a weighting than an error...
          xBinStatisticalWeightError[xBinIndexTemp] = totalYieldError;
          
          // Mark the fractions for all species except for electrons and muons in this bin for regularisation
          // (negative value in special cases, when regularisation is switched off by hand)
          indexParametersToRegularise[internalParIndexTemp++] = factorForSwitchingOffReg_pi > 0 ? (numParamsPerXbin * xBinIndexTemp + 0)
                                                                                            : -1;// Pi
          indexParametersToRegularise[internalParIndexTemp++] = factorForSwitchingOffReg_ka > 0 ? (numParamsPerXbin * xBinIndexTemp + 1)
                                                                                            : -1;// Ka
          indexParametersToRegularise[internalParIndexTemp++] = factorForSwitchingOffReg_pr > 0 ? (numParamsPerXbin * xBinIndexTemp + 2)
                                                                                            : -1;// Pr
          
          printf("indexParametersToRegularise (< 0 if reg off) for bin %d (centre %f): pi %d, ka %d, pr %d", slice + 1,
                 isPtMode ? TMath::Power(10, xBinCentres[xBinIndexTemp]) : xBinCentres[xBinIndexTemp], 
                 indexParametersToRegularise[internalParIndexTemp - 3], indexParametersToRegularise[internalParIndexTemp - 2], 
                 indexParametersToRegularise[internalParIndexTemp - 1]);
          
          // Also mark electrons for regularisation in this bin, if not fixed
          if ( !(mode != kPMxi && ((slice + 1) >= electronFractionThresholdBinForFitting)) &&
               !(mode == kPMxi && ((slice + 1) <= electronFractionThresholdBinForFitting)) ) {
            // Same "switching off" as for pions, since fractions correlated
            indexParametersToRegularise[internalParIndexTemp++] = factorForSwitchingOffReg_pi > 0 ? (numParamsPerXbin * xBinIndexTemp + 3)
                                                                                                  : -1; 
            printf(", el %d", indexParametersToRegularise[internalParIndexTemp - 1]);
          }
          else {
            // Set the index of the last x bin in which the parameter is not fixed.
            // If the parameter is fixed in all x bins, this index will be -1.
            // NOTE: The xi case is also treated properly since the loop starts at the highest xi and moves to lower xi,
            // i.e. the indizes start at highest xi. The index is only used inside AliTPCPIDmathFit to stop the regularisation
            // for the corresponding species if the neighbouring bins are fixed. Reversing the direction of the loop renders
            // it possible to keep the code in AliTPCPIDmathFit as it is.
            if (xBinIndexTemp - 1 < lastNotFixedIndexOfParameters[3])
              lastNotFixedIndexOfParameters[3] = xBinIndexTemp - 1;
          }
          
          // Also mark muons for regularisation in this bin, if not fixed
          if( !(mode != kPMpT || (slice + 1) >= muonFractionThresholdBinForFitting) ) {
            // Same "switching off" as for pions, since fractions correlated
            indexParametersToRegularise[internalParIndexTemp++] = factorForSwitchingOffReg_pi > 0 ? (numParamsPerXbin * xBinIndexTemp + 4)
                                                                                                  : -1; 
            printf(", mu %d", indexParametersToRegularise[internalParIndexTemp - 1]);
          }
          else {
            // Set the index of the last x bin in which the parameter is not fixed.
            // If the parameter is fixed in all x bins, this index will be -1.
            if (xBinIndexTemp - 1 < lastNotFixedIndexOfParameters[4])
              lastNotFixedIndexOfParameters[4] = xBinIndexTemp - 1;
          }
          
          printf("\n");
          
          // TOF yield (binning is the same as for hYieldPt)
          if (applyTOFpatching) {
            // NOTE: TOF pions have already been split into pi and el
            Double_t tofYieldKaonsError = 0, tofYieldPionsError = 0, tofYieldProtonsError = 0, tofYielElectronsError = 0;
            const Double_t tofYieldPions   = hYieldTOFOrigBinningPions->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                         tofYieldPionsError);
            const Double_t tofYieldKaons   = hYieldTOFOrigBinningKaons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                         tofYieldKaonsError);
            const Double_t tofYieldProtons = hYieldTOFOrigBinningProtons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 
                                                                                           tofYieldProtonsError);
            const Double_t tofYieldElectrons = hYieldTOFOrigBinningElectrons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 
                                                                                               tofYielElectronsError);
            
            const Double_t tofYieldTotal = tofYieldPions + tofYieldKaons + tofYieldProtons + tofYieldElectrons;
            const Double_t tofYieldTotalError = TMath::Sqrt(tofYieldPionsError*tofYieldPionsError + tofYieldKaonsError*tofYieldKaonsError
                                                            + tofYieldProtonsError*tofYieldProtonsError
                                                            + tofYielElectronsError*tofYielElectronsError);
            
            // Apply same downscaling as for TPC and ALSO apply this downscaling for the identified yields! Reason:
            // If the electron fraction is calculated or the regularisation term, the TOF patching is applied and/or undone.
            // Since also the TPC fraction has this downscaling, scaling the TOF fraction makes the scale factor to drop
            // out for the (un-)patching as desired, but weakens the regularisation (where also the TOF stat weight contributes).
            xBinStatisticalWeightTOF[xBinIndexTemp] = tofYieldTotal;
            xBinStatisticalWeightErrorTOF[xBinIndexTemp] = tofYieldTotalError;
            
            const Int_t parOffset = xBinIndexTemp * numParamsPerXbin;
            for (Int_t par = 0; par < numParamsPerXbin; par++) {
              // Set all TOF yields to zero (just use the same structure as the parameters itself, although not all of them are used),
              // in particular, non-TOF species like muons and electrons must be set to zero!
              if (par == 0)
                xTOFyield[parOffset + par] = tofYieldPions;
              else if (par == 1)
                xTOFyield[parOffset + par] = tofYieldKaons;
              else if (par == 2)
                xTOFyield[parOffset + par] = tofYieldProtons;
              else if (par == 3)
                xTOFyield[parOffset + par] = tofYieldElectrons;
              else 
                xTOFyield[parOffset + par] = 0.;
            }
          }
          
          xBinIndexTemp++;
        }
        
        if (i == 0) {
          nParToRegulariseSimultaneousFit += AliPID::kSPECIES - 2; // Fracs for all species in this bin except for electrons and muons

          if ( !(mode != kPMxi && ((slice + 1) >= electronFractionThresholdBinForFitting)) &&
               !(mode == kPMxi && ((slice + 1) <= electronFractionThresholdBinForFitting)) )
            nParToRegulariseSimultaneousFit++; // Also regularise electrons in this bin (not fixed)
          
          if ( !(mode != kPMpT || (slice + 1) >= muonFractionThresholdBinForFitting) )
            nParToRegulariseSimultaneousFit++; // Also regularise muons in this bin (not fixed)

          numXBins++;
        }
      }
    }
  }
  AliTPCPIDmathFit* mathFit = 0x0;
  
  if (regularisation > 0) {
    mathFit = (fitMethod == 2) ? AliTPCPIDmathFit::Instance(numXBins, 4, 1810)
                               : AliTPCPIDmathFit::Instance(numXBins, 1, 1810);
  }
  else {
    //TODO FOR INDIVIDUAL DeltaPrime fits: Replace input data and arrays in the following with the desired species and set numSimultaneousFits to 1
    mathFit = (fitMethod == 2) ? AliTPCPIDmathFit::Instance(1, 4, 1810)
                               : AliTPCPIDmathFit::Instance(1, 1, 1810);
  }
  
  mathFit->SetDebugLevel(0); 
  mathFit->SetEpsilon(5e-05);
  mathFit->SetMaxCalls(1e8);
  
  mathFit->SetMinimisationStrategy(minimisationStrategy);
  
  mathFit->SetUseLogLikelihood(useLogLikelihood);
  mathFit->SetUseWeightsForLogLikelihood(useWeightsForLogLikelihood);
  
  if (fitMethod == 2) {
    // If the deltaPrime range is large enough, we artificially get a factor 4 in statistics by looking at the four
    // different deltaPrimeSpecies, which have (except for binning effects) the same information. 
    // Therefore, to get the "real" statistical error, we need to multiply the obtained error by sqrt(4) = 2
    
    //TODO FOR INDIVIDUAL DeltaPrime fits: set to 1.
    mathFit->SetScaleFactorError(2.);
  }
  
  mathFit->SetRegularisation(regularisation, regularisationFactor);
  
  // Number of parameters for fitting
  const Int_t nPar = 11;
  
  // Fracs of each species + total yield in x bin
  const Int_t nParSimultaneousFit = AliPID::kSPECIES + 1; 
  
  // Fracs of each species in x bin + tot yield in x bin
  const Int_t nParSimultaneousFitRegularised = numXBins * (AliPID::kSPECIES + 1); 
  
  // This call is in principle only relevant in case of regularisation. Otherwise, it has no effect.
  mathFit->SetApplyPatching(applyTOFpatching);
  
  if (regularisation > 0) {
    if (!mathFit->SetParametersToRegularise(nParToRegulariseSimultaneousFit, numParamsPerXbin, indexParametersToRegularise,
                                            lastNotFixedIndexOfParameters, xBinCentres, xBinStatisticalWeight, 
                                            xBinStatisticalWeightError, xBinStatisticalWeightTOF, xBinStatisticalWeightErrorTOF, xTOFyield))
      return -1;
  }
  
  delete [] xBinCentres;
  xBinCentres = 0x0;
  
  delete [] xBinStatisticalWeight;
  xBinStatisticalWeight = 0x0;
  
  delete [] xBinStatisticalWeightError;
  xBinStatisticalWeight = 0x0;
  
  delete [] xBinStatisticalWeightTOF;
  xBinStatisticalWeightTOF = 0x0;
  
  delete [] xBinStatisticalWeightErrorTOF;
  xBinStatisticalWeightTOF = 0x0;
  
  delete [] xTOFyield;
  xTOFyield = 0x0;
  
  delete [] indexParametersToRegularise;
  indexParametersToRegularise = 0x0;
  
  delete [] lastNotFixedIndexOfParameters;
  lastNotFixedIndexOfParameters = 0x0;
  
  
  
  gFractionElectronsData = new TGraphErrors(numXBins);
  
  // Fit each slice with sum of 4/5 shapes with means and sigmas fixed from last fitting step
  // For electrons: Fit up to certain pT bin and use constant value for higher momenta
  
  // Two iterations required for regularisation
  Bool_t regularisedFitDone = kFALSE;
  Double_t reducedChiSquareRegularisation = -1;
  
  Double_t gausParamsSimultaneousFitRegularised[nParSimultaneousFitRegularised];
  Double_t parameterErrorsOutRegularised[nParSimultaneousFitRegularised];
  Double_t lowParLimitsSimultaneousFitRegularised[nParSimultaneousFitRegularised];
  Double_t upParLimitsSimultaneousFitRegularised[nParSimultaneousFitRegularised];
  Double_t stepSizeSimultaneousFitRegularised[nParSimultaneousFitRegularised];
  
  for (Int_t i = 0; i < nParSimultaneousFitRegularised; i++) {
    gausParamsSimultaneousFitRegularised[i] = 0;
    parameterErrorsOutRegularised[i] = 0;
    lowParLimitsSimultaneousFitRegularised[i] = 0;
    upParLimitsSimultaneousFitRegularised[i] = 0;
    stepSizeSimultaneousFitRegularised[i] = 0;
  }
  
  Double_t fitRegularisedSpecTemplateFractionErrorKa[numXBins];
  Double_t fitRegularisedSpecTemplateFractionErrorPr[numXBins];
  
  for (Int_t i = 0; i < numXBins; i++) {
    fitRegularisedSpecTemplateFractionErrorKa[i] = 0;
    fitRegularisedSpecTemplateFractionErrorPr[i] = 0;
  }
  
  mathFit->ClearRefHistos();
  
  
  const Int_t nParUsed = (fitMethod == 2) ? ((regularisation <= 0) ? nParSimultaneousFit: nParSimultaneousFitRegularised) : nPar;
  Double_t parameterErrorsOut[nParUsed];
  Double_t covMatrix[nParUsed][nParUsed];
  
  const Bool_t isXiMode = (processingMode == kPMxi);
  const Bool_t isPtMode = (processingMode == kPMpT);
  const Int_t sliceStart = isXiMode ? hFractionPions->GetXaxis()->GetNbins() : 0;
  const Int_t sliceEnd   = isXiMode ? 0. : (isPtMode ? nPtBins : hFractionPions->GetXaxis()->GetNbins());
  
  Bool_t isFirstTimeToWrite = kTRUE;
  Bool_t electronFixingUsed = kFALSE;
  
  for (Int_t iter = 0; iter < 2; iter++) {
    if (regularisation <= 0 && iter == 0)
      continue; // Only one iteration w/o regularisation
  
    Int_t currXbin = 0;
    
    for (Int_t slice = sliceStart; isXiMode ? slice >= sliceEnd : slice < sliceEnd; isXiMode ? slice-- : slice++) { 
      if (isPtMode && (slice < pSliceLow || slice > pSliceHigh))
        continue; 
      
      // There won't (actually: shouldn't) be tracks with a pT larger than the jet pT
      if (isPtMode && restrictJetPtAxis && binsPt[slice] >= actualUpperJetPt)
        continue;
      
      if (regularisation <= 0) {
        if (isPtMode)
          std::cout << "Fitting range " << binsPt[slice] << " GeV/c < Pt < " << binsPt[slice + 1] << " GeV/c..." << std::endl;
        else {
          std::cout << "Fitting range " << hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) << " < " << modeShortName[mode].Data() << " < ";
          std::cout << hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) << "..." << std::endl;
        }
      }
      
      // inverseBinWidth = 1.0, if the raw yield for each bin is requested.
      // If divided by the bin size, the histograms give "yield per unit pT in the corresponding bin" or dN/dpT
      Double_t inverseBinWidth = isPtMode ? 1.0 / (binsPt[slice + 1] - binsPt[slice])
                                          : 1.0 / hYieldPt->GetBinWidth(slice + 1); 
      
      // Add/subtract some very small offset to be sure not to sit on the bin boundary, when looking for the integration/projection limits.
      const Int_t pBinLowProjLimit = isPtMode ? hYieldPt->GetXaxis()->FindFixBin(binsPt[slice] + 1e-5)    : slice + 1;
      const Int_t pBinUpProjLimit  = isPtMode ? hYieldPt->GetXaxis()->FindFixBin(binsPt[slice + 1]- 1e-5) : slice + 1;
      
      Double_t totalYieldError = 0;
      const Double_t totalYield = hYieldPt->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, totalYieldError);
      
      // TOF yield (binning is the same as for hYieldPt)
      if (applyTOFpatching) {
        Double_t tofYieldPionsError = 0;
        const Double_t tofYieldPions = hYieldTOFOrigBinningPions->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                  tofYieldPionsError);
        Double_t tofYieldKaonsError = 0;
        const Double_t tofYieldKaons = hYieldTOFOrigBinningKaons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                  tofYieldKaonsError);
        Double_t tofYieldProtonsError = 0;
        const Double_t tofYieldProtons = hYieldTOFOrigBinningProtons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                      tofYieldProtonsError);
        Double_t tofYieldElectronsError = 0;
        const Double_t tofYieldElectrons = hYieldTOFOrigBinningElectrons->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit,
                                                                                           tofYieldElectronsError);
        
        hYieldTOFPions->SetBinContent(slice + 1, tofYieldPions * inverseBinWidth);
        hYieldTOFPions->SetBinError(slice + 1, tofYieldPionsError * inverseBinWidth);
        
        hYieldTOFKaons->SetBinContent(slice + 1, tofYieldKaons * inverseBinWidth);
        hYieldTOFKaons->SetBinError(slice + 1, tofYieldKaonsError * inverseBinWidth);
        
        hYieldTOFProtons->SetBinContent(slice + 1, tofYieldProtons * inverseBinWidth);
        hYieldTOFProtons->SetBinError(slice + 1, tofYieldProtonsError * inverseBinWidth);
        
        hYieldTOFElectrons->SetBinContent(slice + 1, tofYieldElectrons * inverseBinWidth);
        hYieldTOFElectrons->SetBinError(slice + 1, tofYieldElectronsError * inverseBinWidth);
        
        hYieldTPConlyTotal->SetBinContent(slice + 1, totalYield * inverseBinWidth);
        hYieldTPConlyTotal->SetBinError(slice + 1, totalYieldError * inverseBinWidth);
        
        Double_t MCTOFtotal = -1, MCTOFelectrons = -1, MCTOFkaons = -1, MCTOFmuons = -1, MCTOFpions = -1, MCTOFprotons = -1;
        Double_t MCTOFelectronsErr = 0, MCTOFkaonsErr = 0, MCTOFmuonsErr = 0, MCTOFpionsErr = 0, MCTOFprotonsErr = 0;
        
        MCTOFelectrons = hMCdataTOF->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 1, 1, MCTOFelectronsErr) * inverseBinWidth;
        MCTOFkaons     = hMCdataTOF->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 2, 2, MCTOFkaonsErr)     * inverseBinWidth;
        MCTOFmuons     = hMCdataTOF->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 3, 3, MCTOFmuonsErr)     * inverseBinWidth;
        MCTOFpions     = hMCdataTOF->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 4, 4, MCTOFpionsErr)     * inverseBinWidth;
        MCTOFprotons   = hMCdataTOF->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 5, 5, MCTOFprotonsErr)   * inverseBinWidth;

        MCTOFelectronsErr *= inverseBinWidth;
        MCTOFkaonsErr     *= inverseBinWidth;
        MCTOFmuonsErr     *= inverseBinWidth;
        MCTOFpionsErr     *= inverseBinWidth;
        MCTOFprotonsErr   *= inverseBinWidth;
        
        MCTOFtotal = MCTOFelectrons + MCTOFkaons + MCTOFpions + MCTOFprotons + MCTOFmuons;
        
        if (MCTOFtotal > 0)  {
          hYieldTOFElectronsMC->SetBinContent(slice + 1, MCTOFelectrons);
          hYieldTOFElectronsMC->SetBinError(slice + 1, MCTOFelectronsErr);
          
          hYieldTOFMuonsMC->SetBinContent(slice + 1, MCTOFmuons);
          hYieldTOFMuonsMC->SetBinError(slice + 1, MCTOFmuonsErr);
          
          hYieldTOFKaonsMC->SetBinContent(slice + 1, MCTOFkaons);
          hYieldTOFKaonsMC->SetBinError(slice + 1, MCTOFkaonsErr);
          
          hYieldTOFPionsMC->SetBinContent(slice + 1, MCTOFpions);
          hYieldTOFPionsMC->SetBinError(slice + 1, MCTOFpionsErr);
          
          hYieldTOFProtonsMC->SetBinContent(slice + 1, MCTOFprotons);
          hYieldTOFProtonsMC->SetBinError(slice + 1, MCTOFprotonsErr);
        }
      }
      
      if (totalYield <= 0) {
        std::cout << "Skipped bin (yield is zero)!" << std::endl;
        continue;
      }
      
      if (totalYield < gYieldThresholdForFitting) {
        std::cout << "Skipped bin (yield (" << totalYield << ") is smaller than threshold (" << gYieldThresholdForFitting << "))!" << std::endl;
        continue;
      }
      
      Double_t allDeltaPionError = 0;
      Double_t allDeltaElectronError = 0;
      Double_t allDeltaKaonError = 0;
      Double_t allDeltaProtonError = 0;
      
      const Double_t allDeltaPion = hDeltaPi[slice]->IntegralAndError(hDeltaPi[slice]->GetXaxis()->GetFirst(),
                                                                      hDeltaPi[slice]->GetXaxis()->GetLast(),
                                                                      allDeltaPionError);
      const Double_t allDeltaElectron = hDeltaEl[slice]->IntegralAndError(hDeltaEl[slice]->GetXaxis()->GetFirst(),
                                                                          hDeltaEl[slice]->GetXaxis()->GetLast(),
                                                                          allDeltaElectronError);
      const Double_t allDeltaKaon = hDeltaKa[slice]->IntegralAndError(hDeltaKa[slice]->GetXaxis()->GetFirst(),
                                                                      hDeltaKa[slice]->GetXaxis()->GetLast(),
                                                                      allDeltaKaonError);
      const Double_t allDeltaProton = hDeltaPr[slice]->IntegralAndError(hDeltaPr[slice]->GetXaxis()->GetFirst(),
                                                                        hDeltaPr[slice]->GetXaxis()->GetLast(),
                                                                        allDeltaProtonError);
      
      TH1D *hGenDeltaElForElProj = 0x0, *hGenDeltaKaForElProj = 0x0, *hGenDeltaPiForElProj = 0x0, *hGenDeltaPrForElProj = 0x0;
      TH1D *hGenDeltaElForKaProj = 0x0, *hGenDeltaKaForKaProj = 0x0, *hGenDeltaPiForKaProj = 0x0, *hGenDeltaPrForKaProj = 0x0;
      TH1D *hGenDeltaElForPiProj = 0x0, *hGenDeltaKaForPiProj = 0x0, *hGenDeltaPiForPiProj = 0x0, *hGenDeltaPrForPiProj = 0x0;
      TH1D *hGenDeltaElForMuProj = 0x0, *hGenDeltaKaForMuProj = 0x0, *hGenDeltaPiForMuProj = 0x0, *hGenDeltaPrForMuProj = 0x0;
      TH1D *hGenDeltaElForPrProj = 0x0, *hGenDeltaKaForPrProj = 0x0, *hGenDeltaPiForPrProj = 0x0, *hGenDeltaPrForPrProj = 0x0;
      
      
      TH2D* hGenDeltaUsed[6][6];
      if (useIdentifiedGeneratedSpectra) { 
        for (Int_t i = 0; i < 6; i++) {
          for (Int_t j = 0; j < 6; j++) {
            hGenDeltaUsed[i][j] = hGenDeltaMCid[i][j];
          }
        }
      }
      else {
        for (Int_t i = 0; i < 6; i++) {
          for (Int_t j = 0; j < 6; j++) {
            hGenDeltaUsed[i][j] = hGenDelta[i][j];
          }
        }
      }
      
      Bool_t mostProbPIDLowPtTemplatesForKa = kFALSE;
      Bool_t mostProbPIDLowPtTemplatesForPr = kFALSE;
      
      if (!isMCdataSet) {
        // Do not use special templates in MC (instead of the most probable PID, they have the MC ID stored!), since the thresholds
        // are different to data due to different dEdx shapes
        if (mode == kPMpT) {
          mostProbPIDLowPtTemplatesForKa = (hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) <= mostProbPIDthresholdPt_ka);
          mostProbPIDLowPtTemplatesForPr = (hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) <= mostProbPIDthresholdPt_pr);
        }
        else if (mode == kPMz) {
          mostProbPIDLowPtTemplatesForKa = (hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) <= mostProbPIDthresholdZ_ka);
          mostProbPIDLowPtTemplatesForPr = (hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) <= mostProbPIDthresholdZ_pr);
        }
        else if (mode == kPMxi) {
          mostProbPIDLowPtTemplatesForKa = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) >= mostProbPIDthresholdXi_ka);
          mostProbPIDLowPtTemplatesForPr = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) >= mostProbPIDthresholdXi_pr);
        }
        else if (mode == kPMdistance) {
          mostProbPIDLowPtTemplatesForKa = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) <= mostProbPIDthresholdR_ka);
          mostProbPIDLowPtTemplatesForPr = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) <= mostProbPIDthresholdR_pr);
        }
        else if (mode == kPMjT) {
          mostProbPIDLowPtTemplatesForKa = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) <= mostProbPIDthresholdJt_ka);
          mostProbPIDLowPtTemplatesForPr = (hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) <= mostProbPIDthresholdJt_pr);
        }
      }
      
      
      
      hGenDeltaElForElProj =(TH1D*)hGenDeltaUsed[kEl][kEl]->ProjectionY(Form("hGenDeltaElForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaKaForElProj =(TH1D*)hGenDeltaUsed[kKa][kEl]->ProjectionY(Form("hGenDeltaKaForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPiForElProj =(TH1D*)hGenDeltaUsed[kPi][kEl]->ProjectionY(Form("hGenDeltaPiForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPrForElProj =(TH1D*)hGenDeltaUsed[kPr][kEl]->ProjectionY(Form("hGenDeltaPrForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
      // ALWAYS replace templates in desired region with those from most prob ID - one could take difference to fits as sys error,
      // but this is kind of artificial since there is clean ID possible!
      if (mostProbPIDLowPtTemplatesForKa) {
        hGenDeltaElForKaProj = (TH1D*)hDeltaElMC[slice][1]->Clone(Form("hGenDeltaElForKaProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaElForKaProj);
        hGenDeltaKaForKaProj = (TH1D*)hDeltaKaMC[slice][1]->Clone(Form("hGenDeltaKaForKaProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaKaForKaProj);
        hGenDeltaPiForKaProj = (TH1D*)hDeltaPiMC[slice][1]->Clone(Form("hGenDeltaPiForKaProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaPiForKaProj);
        hGenDeltaPrForKaProj = (TH1D*)hDeltaPrMC[slice][1]->Clone(Form("hGenDeltaPrForKaProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaPrForKaProj);
      }
      else {
        hGenDeltaElForKaProj =(TH1D*)hGenDeltaUsed[kEl][kKa]->ProjectionY(Form("hGenDeltaElForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaKaForKaProj =(TH1D*)hGenDeltaUsed[kKa][kKa]->ProjectionY(Form("hGenDeltaKaForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPiForKaProj =(TH1D*)hGenDeltaUsed[kPi][kKa]->ProjectionY(Form("hGenDeltaPiForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPrForKaProj =(TH1D*)hGenDeltaUsed[kPr][kKa]->ProjectionY(Form("hGenDeltaPrForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      }
      
      hGenDeltaElForPiProj =(TH1D*)hGenDeltaUsed[kEl][kPi]->ProjectionY(Form("hGenDeltaElForPiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaKaForPiProj =(TH1D*)hGenDeltaUsed[kKa][kPi]->ProjectionY(Form("hGenDeltaKaForPiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPiForPiProj =(TH1D*)hGenDeltaUsed[kPi][kPi]->ProjectionY(Form("hGenDeltaPiForPiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPrForPiProj =(TH1D*)hGenDeltaUsed[kPr][kPi]->ProjectionY(Form("hGenDeltaPrForPiProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
      if (takeIntoAccountMuons) {
        hGenDeltaElForMuProj =(TH1D*)hGenDeltaUsed[kEl][kMu]->ProjectionY(Form("hGenDeltaElForMuProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaKaForMuProj =(TH1D*)hGenDeltaUsed[kKa][kMu]->ProjectionY(Form("hGenDeltaKaForMuProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPiForMuProj =(TH1D*)hGenDeltaUsed[kPi][kMu]->ProjectionY(Form("hGenDeltaPiForMuProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPrForMuProj =(TH1D*)hGenDeltaUsed[kPr][kMu]->ProjectionY(Form("hGenDeltaPrForMuProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      }
      
      // ALWAYS replace templates in desired region with those from most prob ID - one could take difference to fits as sys error,
      // but this is kind of artificial since there is clean ID possible!
      if (mostProbPIDLowPtTemplatesForPr) {
        hGenDeltaElForPrProj = (TH1D*)hDeltaElMC[slice][4]->Clone(Form("hGenDeltaElForPrProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaElForPrProj);
        hGenDeltaKaForPrProj = (TH1D*)hDeltaKaMC[slice][4]->Clone(Form("hGenDeltaKaForPrProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaKaForPrProj);
        hGenDeltaPiForPrProj = (TH1D*)hDeltaPiMC[slice][4]->Clone(Form("hGenDeltaPiForPrProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaPiForPrProj);
        hGenDeltaPrForPrProj = (TH1D*)hDeltaPrMC[slice][4]->Clone(Form("hGenDeltaPrForPrProj%d", slice));
        SetStyleSpecialTemplate(hGenDeltaPrForPrProj);
      }
      else {
        hGenDeltaElForPrProj =(TH1D*)hGenDeltaUsed[kEl][kPr]->ProjectionY(Form("hGenDeltaElForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaKaForPrProj =(TH1D*)hGenDeltaUsed[kKa][kPr]->ProjectionY(Form("hGenDeltaKaForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPiForPrProj =(TH1D*)hGenDeltaUsed[kPi][kPr]->ProjectionY(Form("hGenDeltaPiForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
        hGenDeltaPrForPrProj =(TH1D*)hGenDeltaUsed[kPr][kPr]->ProjectionY(Form("hGenDeltaPrForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      }
      
      // Errors for kaon and proton fraction in case of special templates - separately for each fit method
      Double_t specTemplateFractionErrorKa = 0.;
      Double_t specTemplateFractionErrorKaDeltaEl = 0.;
      Double_t specTemplateFractionErrorKaDeltaPi = 0.;
      Double_t specTemplateFractionErrorKaDeltaKa = 0.;
      Double_t specTemplateFractionErrorKaDeltaPr = 0.;
      Double_t specTemplateFractionErrorPr = 0.;
      Double_t specTemplateFractionErrorPrDeltaEl = 0.;
      Double_t specTemplateFractionErrorPrDeltaPi = 0.;
      Double_t specTemplateFractionErrorPrDeltaKa = 0.;
      Double_t specTemplateFractionErrorPrDeltaPr = 0.;
      
      Bool_t nonEmptyTemplateEl = kFALSE;
      Bool_t nonEmptyTemplatePi = kFALSE;
      Bool_t nonEmptyTemplateKa = kFALSE;
      Bool_t nonEmptyTemplateMu = kFALSE;
      Bool_t nonEmptyTemplatePr = kFALSE;
      
      if (fitMethod == 2) {
        // Calculate error for special templates - here, no need to scale the errors, since no artificial grow in statistics!
        if (mostProbPIDLowPtTemplatesForKa) {
          specTemplateFractionErrorKa = GetFractionErrorSpecialTemplate(totalYield, totalYieldError, hDeltaKaMC[slice][1]);
        }
        if (mostProbPIDLowPtTemplatesForPr) {
          specTemplateFractionErrorPr = GetFractionErrorSpecialTemplate(totalYield, totalYieldError, hDeltaPrMC[slice][4]);
        }
        
        
        // Normalise generated histos to TOTAL number of GENERATED particles for this species (i.e. including
        // entries that lie in the under or overflow bin), so that situations in which the generated spectra lie
        // at least partly outside the histo are treated properly. To find the total number of generated particle
        // species X, one can just take the integral of the generated histo for DeltaX (which should include all
        // generated entries) and apply the same normalisation factor to all other DeltaSpecies.
        // Also set some cosmetics
        
        // NOTE: Just the same for low-pT templates! Here, the number of generated entries is replaced by the actual measured
        // number of tracks = templates
        
        // Generated electrons
        // If template is empty, it is empty for all delta_x!
        Double_t normEl = normaliseHist(hGenDeltaElForElProj, nonEmptyTemplateEl, -1);
        normaliseHist(hGenDeltaKaForElProj, normEl);
        normaliseHist(hGenDeltaPiForElProj, normEl);
        normaliseHist(hGenDeltaPrForElProj, normEl);
        
        
        // Generated kaons
        // If template is empty, it is empty for all delta_x!
        Double_t normKa = normaliseHist(hGenDeltaKaForKaProj, nonEmptyTemplateKa, -1);
        normaliseHist(hGenDeltaElForKaProj, normKa);
        normaliseHist(hGenDeltaPiForKaProj, normKa);
        normaliseHist(hGenDeltaPrForKaProj, normKa);
        
        
        // Generated pions
        // If template is empty, it is empty for all delta_x!
        Double_t normPi = normaliseHist(hGenDeltaPiForPiProj, nonEmptyTemplatePi, -1);
        normaliseHist(hGenDeltaElForPiProj, normPi);
        normaliseHist(hGenDeltaKaForPiProj, normPi);
        normaliseHist(hGenDeltaPrForPiProj, normPi);
        

        Double_t normMu = 1;
        if (takeIntoAccountMuons) {
          // Generated pions
          // If template is empty, it is empty for all delta_x!
          // Since masses of muons and pions are so similar, the normalisation scheme should still work when looking at deltaPion instead
          normMu = normaliseHist(hGenDeltaPiForMuProj, nonEmptyTemplateMu, -1);
          normaliseHist(hGenDeltaElForMuProj, normMu);
          normaliseHist(hGenDeltaKaForMuProj, normMu);
          normaliseHist(hGenDeltaPrForMuProj, normMu);
        }
        
        
        // Generated protons
        // If template is empty, it is empty for all delta_x!
        Double_t normPr = normaliseHist(hGenDeltaPrForPrProj, nonEmptyTemplatePr, -1);
        normaliseHist(hGenDeltaElForPrProj, normPr);
        normaliseHist(hGenDeltaKaForPrProj, normPr);
        normaliseHist(hGenDeltaPiForPrProj, normPr);
      }
      else {
        // Calculate error for special templates - here, no need to scale the errors, since no artificial grow in statistics!
        if (mostProbPIDLowPtTemplatesForKa) {
          specTemplateFractionErrorKaDeltaEl = GetFractionErrorSpecialTemplate(allDeltaElectron, allDeltaElectronError, hDeltaElMC[slice][1]);
          specTemplateFractionErrorKaDeltaKa = GetFractionErrorSpecialTemplate(allDeltaKaon, allDeltaKaonError, hDeltaKaMC[slice][1]);
          specTemplateFractionErrorKaDeltaPi = GetFractionErrorSpecialTemplate(allDeltaPion, allDeltaPionError, hDeltaPiMC[slice][1]);
          specTemplateFractionErrorKaDeltaPr = GetFractionErrorSpecialTemplate(allDeltaProton, allDeltaProtonError, hDeltaPrMC[slice][1]);
        }
        if (mostProbPIDLowPtTemplatesForPr) {
          specTemplateFractionErrorPrDeltaEl = GetFractionErrorSpecialTemplate(allDeltaElectron, allDeltaElectronError, hDeltaElMC[slice][4]);
          specTemplateFractionErrorPrDeltaKa = GetFractionErrorSpecialTemplate(allDeltaKaon, allDeltaKaonError, hDeltaKaMC[slice][4]);
          specTemplateFractionErrorPrDeltaPi = GetFractionErrorSpecialTemplate(allDeltaPion, allDeltaPionError, hDeltaPiMC[slice][4]);
          specTemplateFractionErrorPrDeltaPr = GetFractionErrorSpecialTemplate(allDeltaProton, allDeltaProtonError, hDeltaPrMC[slice][4]);
        }
        
        
        // Normalise generated histos to total number of particles for this delta
        // and also set some cosmetics
        
       // NOTE: Just the same normalisation for low-pT templates! 
        
        // DeltaEl
        normaliseHist(hGenDeltaElForElProj, -1);
        normaliseHist(hGenDeltaElForKaProj, -1);
        normaliseHist(hGenDeltaElForPiProj, -1);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaElForMuProj, -1);
        normaliseHist(hGenDeltaElForPrProj, -1);
        
        // DeltaKa
        normaliseHist(hGenDeltaKaForElProj, -1);
        normaliseHist(hGenDeltaKaForKaProj, -1);
        normaliseHist(hGenDeltaKaForPiProj, -1);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaKaForMuProj, -1);
        normaliseHist(hGenDeltaKaForPrProj, -1);
        
        // DeltaPi
        normaliseHist(hGenDeltaPiForElProj, -1);
        normaliseHist(hGenDeltaPiForKaProj, -1);
        normaliseHist(hGenDeltaPiForPiProj, -1);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaPiForMuProj, -1);
        normaliseHist(hGenDeltaPiForPrProj, -1);
        
        // DeltaPr
        normaliseHist(hGenDeltaPrForElProj, -1);
        normaliseHist(hGenDeltaPrForKaProj, -1);
        normaliseHist(hGenDeltaPrForPiProj, -1);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaPrForMuProj, -1);
        normaliseHist(hGenDeltaPrForPrProj, -1);
        
        // If template is empty, it is empty for all delta_x! Look at delta_species for species to ensure that the entries are not just
        // outside of the range! Take delta_pi for mu.
        nonEmptyTemplateEl = hGenDeltaElForElProj->Integral() > 0;
        nonEmptyTemplatePi = hGenDeltaPiForPiProj->Integral() > 0;
        nonEmptyTemplateKa = hGenDeltaKaForKaProj->Integral() > 0;
        nonEmptyTemplatePr = hGenDeltaPrForPrProj->Integral() > 0;
        if (takeIntoAccountMuons)
          nonEmptyTemplateMu = hGenDeltaPiForMuProj->Integral() > 0;
      }

      
      TF1* totalDeltaPion = 0x0;
      TF1* totalDeltaKaon = 0x0;
      TF1* totalDeltaProton = 0x0;
      TF1* totalDeltaElectron = 0x0;
      
      TLegend* legend = 0x0;
      
      if (iter == 1) { // Only needed for second iteration (= the only iteration w/o regularisation)
        // The number of parameters and their values will always be adjusted, such that using nPar parameters is fine
        totalDeltaPion = new TF1(Form("totalDeltaPion_slice%d", slice), (fitMethod == 2) ? multiGaussFitDeltaPi : multiGaussFit, 
                                 xLow, xUp, nPar);
        setUpFitFunction(totalDeltaPion, nBins);
        
        totalDeltaKaon = new TF1(Form("totalDeltaKaon_slice%d", slice), (fitMethod == 2) ? multiGaussFitDeltaKa : multiGaussFit,
                                 xLow, xUp, nPar);
        setUpFitFunction(totalDeltaKaon, nBins);
        
        totalDeltaProton = new TF1(Form("totalDeltaProton_slice%d", slice), (fitMethod == 2) ? multiGaussFitDeltaPr : multiGaussFit, 
                                   xLow, xUp, nPar);
        setUpFitFunction(totalDeltaProton, nBins);
        
        totalDeltaElectron = new TF1(Form("totalDeltaElectron_slice%d", slice), 
                                     (fitMethod == 2) ? multiGaussFitDeltaEl : multiGaussFit,
                                     xLow, xUp, nPar);
        setUpFitFunction(totalDeltaElectron, nBins);
        
        // Legend is the same for all \Delta "species" plots
        legend = new TLegend(0.722126, 0.605932, 0.962069, 0.925932);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        if (plotIdentifiedSpectra)
          legend->SetNColumns(2);
        legend->AddEntry((TObject*)0x0, "Fit", "");
        if (plotIdentifiedSpectra)
          legend->AddEntry((TObject*)0x0, identifiedLabels[isMC].Data(), "");
        
        legend->AddEntry(hDeltaPi[slice], "Data", "Lp");
        if (plotIdentifiedSpectra)
          legend->AddEntry((TObject*)0x0, "", "");
        
        legend->AddEntry(totalDeltaPion, "Multi-template fit", "L");
        if (plotIdentifiedSpectra)
          legend->AddEntry(hMCmuonsAndPionsDummy, "#mu + #pi", "Lp");
        
        if (takeIntoAccountMuons && (muonFractionHandling != kNoMuons))
          legend->AddEntry(hGenDeltaPiForMuProj, "#mu", "Lp");
        else if (plotIdentifiedSpectra)
          legend->AddEntry((TObject*)0x0, "", "");
        if (plotIdentifiedSpectra)
          legend->AddEntry(hDeltaPiMC[slice][kMu - 1], "#mu", "Lp");
        
        legend->AddEntry(hGenDeltaPiForPiProj, takeIntoAccountMuons ? "#pi" : "#pi + #mu", "Lp");
        if (plotIdentifiedSpectra)
          legend->AddEntry(hDeltaPiMC[slice][kPi - 1], "#pi", "Lp");
        
        legend->AddEntry(hGenDeltaPiForKaProj, "K", "Lp");
        if (plotIdentifiedSpectra)
          legend->AddEntry(hDeltaPiMC[slice][kKa - 1], "K", "Lp");
        
        legend->AddEntry(hGenDeltaPiForPrProj, "p", "Lp");
        if (plotIdentifiedSpectra)
          legend->AddEntry(hDeltaPiMC[slice][kPr - 1], "p", "Lp");
        
        legend->AddEntry(hGenDeltaPiForElProj, "e", "Lp");
        if (plotIdentifiedSpectra) 
          legend->AddEntry(hDeltaPiMC[slice][kEl -1], "e", "Lp");
      }
      
      
      // Allow tolerance of +-2% (for delta -> assume dEdx = 80 and take +-2%)
      //const Double_t peakTolerance = (useDeltaPrime ? 0.8 : 1.0) / hGenDeltaElForElProj->GetXaxis()->GetBinWidth(1);
      //const Double_t shiftStepSize = 0.01;
      const Double_t peakTolerance = (useDeltaPrime ? 0.02 : 1.6);
      const Double_t shiftStepSize = 0.01;
      
      // Assume fractions vs. pT to be smooth. Allow 1 sigma variations from bin to bin. For small pT, the error will be very small.
      // Therefore, allow at least a change of some percent.
      const Double_t nSigma = 1.;
      const Double_t minChange = 1.0; // This disables the sigma restriction
      
      Double_t fractionPions = (muonContamination ? 0.87 : 0.88);
      
      Double_t fractionErrorUpPions = 1.;
      Double_t fractionErrorLowPions = 0.;
      
      Int_t xBinInit = 0;
      if (initialiseWithFractionsFromFile) {
        Double_t xBinCentre = isPtMode ? (binsPt[slice + 1] + binsPt[slice]) / 2.
                                       : hYieldPt->GetXaxis()->GetBinCenter(slice + 1); 
        xBinInit = hInitFracPi->GetXaxis()->FindFixBin(xBinCentre);
        fractionPions = hInitFracPi->GetBinContent(xBinInit) + (muonContamination ? hInitFracMu->GetBinContent(xBinInit) : 0.);
      }
      else {
        // Set found fraction from last slice, if available. Note: Current bin for slice = slice + 1
        // => Bin for last slice = slice
        if (hFractionPions->GetBinContent(slice) > 0 && hFractionPions->GetBinError(slice) > 0) {
          fractionPions = hFractionPions->GetBinContent(slice);
          fractionErrorUpPions = TMath::Min(1.0, fractionPions + TMath::Max(minChange, nSigma * hFractionPions->GetBinError(slice)));
          fractionErrorLowPions = TMath::Max(0.0, fractionPions - TMath::Max(minChange, nSigma * hFractionPions->GetBinError(slice)));
        }
      }
      
      Double_t fractionKaons = 0.08;
      Double_t fractionErrorUpKaons = 1.;
      Double_t fractionErrorLowKaons = 0.;
      
      if (initialiseWithFractionsFromFile) {
        fractionKaons = hInitFracKa->GetBinContent(xBinInit);
      }
      else {
        if (hFractionKaons->GetBinContent(slice) > 0 && hFractionKaons->GetBinError(slice) > 0) {
          fractionKaons = hFractionKaons->GetBinContent(slice);
          fractionErrorUpKaons = TMath::Min(1.0, fractionKaons + TMath::Max(minChange, nSigma * hFractionKaons->GetBinError(slice)));
          fractionErrorLowKaons = TMath::Max(0.0, fractionKaons - TMath::Max(minChange, nSigma * hFractionKaons->GetBinError(slice)));
        }
      }
       
      Double_t fractionProtons = 0.02;
      Double_t fractionErrorUpProtons = 1.;
      Double_t fractionErrorLowProtons = 0.;
      
      if (initialiseWithFractionsFromFile) {
        fractionProtons = hInitFracPr->GetBinContent(xBinInit);
      }
      else {
        if (hFractionProtons->GetBinContent(slice) > 0 && hFractionProtons->GetBinError(slice) > 0) {
          fractionProtons = hFractionProtons->GetBinContent(slice);
          fractionErrorUpProtons = TMath::Min(1.0, fractionProtons +
                                                   TMath::Max(minChange, nSigma * hFractionProtons->GetBinError(slice)));
          fractionErrorLowProtons = TMath::Max(0.0, fractionProtons -
                                                    TMath::Max(minChange, nSigma * hFractionProtons->GetBinError(slice)));
        }
      }
        
      Double_t fractionElectrons = (takeIntoAccountMuons ? 0.01 : 0.02);
      Double_t fractionErrorUpElectrons = 1.;
      Double_t fractionErrorLowElectrons = 0.;
      
      if (initialiseWithFractionsFromFile) {
        fractionElectrons = hInitFracEl->GetBinContent(xBinInit);
      }
      else {
        if (hFractionElectrons->GetBinContent(slice) > 0 && hFractionElectrons->GetBinError(slice) > 0) {
          fractionElectrons = hFractionElectrons->GetBinContent(slice);
          fractionErrorUpElectrons = TMath::Min(1.0, fractionElectrons + 
                                                     TMath::Max(minChange, nSigma * hFractionElectrons->GetBinError(slice)));
          fractionErrorLowElectrons = TMath::Max(0.0, fractionElectrons -
                                                      TMath::Max(minChange, nSigma * hFractionElectrons->GetBinError(slice)));
        }
      }
      
      Double_t fractionMuons = (takeIntoAccountMuons ? 0.01 : 0.);
      Double_t fractionErrorUpMuons = 1.;
      Double_t fractionErrorLowMuons = 0.;
      if (!takeIntoAccountMuons) {
        fractionErrorUpMuons = 0.;
        fractionErrorLowMuons = 0.;
      }
      else {
        if (initialiseWithFractionsFromFile) {
          fractionMuons = hInitFracMu->GetBinContent(xBinInit);
        }
        else {
          if (hFractionMuons->GetBinContent(slice) > 0 && hFractionMuons->GetBinError(slice) > 0) {
            fractionMuons = hFractionMuons->GetBinContent(slice);
            fractionErrorUpMuons = TMath::Min(1.0, fractionMuons + TMath::Max(minChange, nSigma * hFractionMuons->GetBinError(slice)));
            fractionErrorLowMuons = TMath::Max(0.0, fractionMuons - TMath::Max(minChange, nSigma * hFractionMuons->GetBinError(slice)));
          }
        }
      }
      
      // If the templates are empty, they cannot be used for the fit (fractions do not change the chi^2 (except for regularisation)).
      // Fix fractions to zero in that case (usually only happens for protons (or kaons) at very low pT, where pre-PID works rather
      // perfectly); stepSize needs to be set to 0 as well!
      // Note that K and p will be overwritten anyway by the special templates. But empty templates should imply that also the special
      // templates are empty, such that it is equivalent.
      if (!nonEmptyTemplateEl) {
        fractionElectrons = 0.;
        fractionErrorUpElectrons = 0.;
        fractionErrorLowElectrons = 0.;
        
        printf("\nEl template empty in this bin - fixing fractions to zero!\n\n");
      }
      if (!nonEmptyTemplatePi) {
        fractionPions = 0.;
        fractionErrorUpPions = 0.;
        fractionErrorLowPions = 0.;
        
        printf("\nPi template empty in this bin - fixing fractions to zero!\n\n");
      }
      if (!nonEmptyTemplateKa) {
        fractionKaons = 0.;
        fractionErrorUpKaons = 0.;
        fractionErrorLowKaons = 0.;
        
        printf("\nKa template empty in this bin - fixing fractions to zero!\n\n");
      }
      if (!nonEmptyTemplatePr) {
        fractionProtons = 0.;
        fractionErrorUpProtons = 0.;
        fractionErrorLowProtons = 0.;
        
        printf("\nPr template empty in this bin - fixing fractions to zero!\n\n");
      }
      if (!nonEmptyTemplateMu && takeIntoAccountMuons) {
        fractionMuons = 0.;
        fractionErrorUpMuons = 0.;
        fractionErrorLowMuons = 0.;
        
        printf("\nMu template empty in this bin - fixing fractions to zero!\n\n");
      }
      
      
      // In case of special (low-pT) templates, fix the fraction to the correct value (i.e. set mean, error limits and step size
      // accordingly). Note that the unnormalised templates are required to extrat the proper fraction
      
      Double_t gausParamsPi[nPar] = { 
        fractionPions,
        mostProbPIDLowPtTemplatesForKa ? GetFractionSpecialTemplate(allDeltaPion, hDeltaPiMC[slice][1]) : fractionKaons,
        mostProbPIDLowPtTemplatesForPr ? GetFractionSpecialTemplate(allDeltaPion, hDeltaPiMC[slice][4]) : fractionProtons,
        fractionElectrons,
        fractionMuons,
        allDeltaPion,
        0,
        0,
        0,
        0,
        0
      };
      
      Double_t gausParamsEl[nPar] = { 
        fractionPions,
        mostProbPIDLowPtTemplatesForKa ? GetFractionSpecialTemplate(allDeltaElectron, hDeltaElMC[slice][1]) : fractionKaons,
        mostProbPIDLowPtTemplatesForPr ? GetFractionSpecialTemplate(allDeltaElectron, hDeltaElMC[slice][4]) : fractionProtons,
        fractionElectrons,
        fractionMuons,
        allDeltaElectron,
        0,
        0,
        0,
        0,
        0
      };
      
      Double_t gausParamsKa[nPar] = { 
        fractionPions,
        mostProbPIDLowPtTemplatesForKa ? GetFractionSpecialTemplate(allDeltaKaon, hDeltaKaMC[slice][1]) : fractionKaons,
        mostProbPIDLowPtTemplatesForPr ? GetFractionSpecialTemplate(allDeltaKaon, hDeltaKaMC[slice][4]) : fractionProtons,
        fractionElectrons,
        fractionMuons,
        allDeltaKaon,
        0,
        0,
        0,
        0,
        0
      };
      
      Double_t gausParamsPr[nPar] = { 
        fractionPions,
        mostProbPIDLowPtTemplatesForKa ? GetFractionSpecialTemplate(allDeltaProton, hDeltaPrMC[slice][1]) : fractionKaons,
        mostProbPIDLowPtTemplatesForPr ? GetFractionSpecialTemplate(allDeltaProton, hDeltaPrMC[slice][4]) : fractionProtons,
        fractionElectrons,
        fractionMuons,
        allDeltaProton,
        0,
        0,
        0,
        0,
        0
      };
      
      Double_t lowParLimitsPi[nPar] = {
        fractionErrorLowPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsPi[1] : fractionErrorLowKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsPi[2] : fractionErrorLowProtons,
        fractionErrorLowElectrons,
        fractionErrorLowMuons,
        allDeltaPion,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance
      };
      
      Double_t lowParLimitsEl[nPar] = {
        fractionErrorLowPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsEl[1] : fractionErrorLowKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsEl[2] : fractionErrorLowProtons,
        fractionErrorLowElectrons,
        fractionErrorLowMuons,
        allDeltaElectron,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance
      };
      
      Double_t lowParLimitsKa[nPar] = {
        fractionErrorLowPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsKa[1] : fractionErrorLowKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsKa[2] : fractionErrorLowProtons,
        fractionErrorLowElectrons,
        fractionErrorLowMuons,
        allDeltaKaon,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance
      };
      
      Double_t lowParLimitsPr[nPar] = {
        fractionErrorLowPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsPr[1] : fractionErrorLowKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsPr[2] : fractionErrorLowProtons,
        fractionErrorLowElectrons,
        fractionErrorLowMuons,
        allDeltaProton,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance,
        -peakTolerance
        -peakTolerance
      };
      
      
      Double_t upParLimitsPi[nPar] = {
        fractionErrorUpPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsPi[1] : fractionErrorUpKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsPi[2] : fractionErrorUpProtons,
        fractionErrorUpElectrons,
        fractionErrorUpMuons,
        allDeltaPion,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance
      };
      
      Double_t upParLimitsEl[nPar] = {
        fractionErrorUpPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsEl[1] : fractionErrorUpKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsEl[2] : fractionErrorUpProtons,
        fractionErrorUpElectrons,
        fractionErrorUpMuons,
        allDeltaElectron,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance
      };
      
      Double_t upParLimitsKa[nPar] = {
        fractionErrorUpPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsKa[1] : fractionErrorUpKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsKa[2] : fractionErrorUpProtons,
        fractionErrorUpElectrons,
        fractionErrorUpMuons,
        allDeltaKaon,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance
      };
      
      Double_t upParLimitsPr[nPar] = {
        fractionErrorUpPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsPr[1] : fractionErrorUpKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsPr[2] : fractionErrorUpProtons,
        fractionErrorUpElectrons,
        fractionErrorUpMuons,
        allDeltaProton,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance,
        peakTolerance,
      };
      
      Double_t stepSize[nPar] = {
        nonEmptyTemplatePi ? 0.1 : 0.,
        mostProbPIDLowPtTemplatesForKa ? 0. : (nonEmptyTemplateKa ? 0.1 : 0.),
        mostProbPIDLowPtTemplatesForPr ? 0. : (nonEmptyTemplatePr ? 0.1 : 0.),
        nonEmptyTemplateEl ? 0.1 : 0.,
        (takeIntoAccountMuons ? (nonEmptyTemplateMu ? 0.1 : 0.) : 0.),
        
        0.0,
        
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        (enableShift && takeIntoAccountMuons) ? shiftStepSize : 0.
      };
      
      
      Double_t gausParamsSimultaneousFit[nParSimultaneousFit] = { 
        fractionPions,
        mostProbPIDLowPtTemplatesForKa ? GetFractionSpecialTemplate(totalYield, hDeltaKaMC[slice][1]) : fractionKaons,
        mostProbPIDLowPtTemplatesForPr ? GetFractionSpecialTemplate(totalYield, hDeltaPrMC[slice][4]) : fractionProtons,
        fractionElectrons,
        fractionMuons,
        totalYield
        // No shifts because they do not make too much sense (different eta + possible deviations from Bethe-Bloch in one x-Bin)
      };
      
      Double_t lowParLimitsSimultaneousFit[nParSimultaneousFit] = {
        fractionErrorLowPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsSimultaneousFit[1] : fractionErrorLowKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsSimultaneousFit[2] : fractionErrorLowProtons,
        fractionErrorLowElectrons,
        fractionErrorLowMuons,
        totalYield
      };
      
      Double_t upParLimitsSimultaneousFit[nParSimultaneousFit] = {
        fractionErrorUpPions,
        mostProbPIDLowPtTemplatesForKa ? gausParamsSimultaneousFit[1] : fractionErrorUpKaons,
        mostProbPIDLowPtTemplatesForPr ? gausParamsSimultaneousFit[2] : fractionErrorUpProtons,
        fractionErrorUpElectrons,
        fractionErrorUpMuons,
        totalYield
      };
      
      Double_t stepSizeSimultaneousFit[nParSimultaneousFit] = {
        nonEmptyTemplatePi ? 0.1 : 0.,
        mostProbPIDLowPtTemplatesForKa ? 0. : (nonEmptyTemplateKa ? 0.1 : 0.),
        mostProbPIDLowPtTemplatesForPr ? 0. : (nonEmptyTemplatePr ? 0.1 : 0.),
        nonEmptyTemplateEl ? 0.1 : 0.,
        (takeIntoAccountMuons ? (nonEmptyTemplateMu ? 0.1 : 0.) : 0.),
        
        0.0
      };
      
      
      
      if (regularisation <= 0 && iter == 1) {
        // In case of no regularisation, do the fit of the electron fraction here (compare comment below)
        if ((slice + 1) == electronFractionThresholdBinForFitting) {
          // If TOF patching is to be applied, the TOF patched fraction should be described by the fitted curve.
          // Thus, patch the fraction, apply the fit and afterwards undo the TOF patching to have to proper fraction for the fit
          // (that will be converted to the final desired results after applying TOF patching again).
          
          if (applyTOFpatching) {
            // At this point, temporarily clone the electron histogram and patch its content.
            TH1F* hFractionElectronsTOFpatchedTemp = new TH1F(*hFractionElectrons);
            
            for (Int_t iBin = 1; iBin <= hFractionElectronsTOFpatchedTemp->GetNbinsX(); iBin++) {
              Double_t fracTPC = hFractionElectrons->GetBinContent(iBin);
              Double_t fracErrorTPC = hFractionElectrons->GetBinError(iBin);
              
              // NOTE: All yields have inverseBinWidth factor included -> It doesn't matter, if it is really for all yields because
              // it just drops out then.
              const Double_t yieldTOFpions   = hYieldTOFPions->GetBinContent(slice + 1);
              const Double_t yieldTOFkaons   = hYieldTOFKaons->GetBinContent(slice + 1);
              const Double_t yieldTOFprotons = hYieldTOFProtons->GetBinContent(slice + 1);
              
              // No TOF electrons used at the moment, but histogram filled from contamination of TOF pions
              const Double_t yieldTOFelectrons = hYieldTOFElectrons->GetBinContent(slice + 1); 
              
              const Double_t yieldTOFtotal = yieldTOFelectrons + yieldTOFpions + yieldTOFkaons + yieldTOFprotons;
              
              const Double_t yieldTPCtotal = hYieldTPConlyTotal->GetBinContent(slice + 1);
              
              PatchFractionWithTOF(fracTPC, fracErrorTPC, yieldTPCtotal, yieldTOFelectrons, yieldTOFtotal);
              
              hFractionElectronsTOFpatchedTemp->SetBinContent(iBin, fracTPC);
              hFractionElectronsTOFpatchedTemp->SetBinError(iBin, fracErrorTPC);
            }
            
            if (processingMode == kPMxi)
              hFractionElectronsTOFpatchedTemp->Fit(fElectronFraction, "N", "", electronFractionThresholdForFitting,
                                                    lowFittingBoundElectronFraction);
            else
              hFractionElectronsTOFpatchedTemp->Fit(fElectronFraction, "N", "", lowFittingBoundElectronFraction,
                                                    electronFractionThresholdForFitting);
            delete hFractionElectronsTOFpatchedTemp;
          }
          else {
            if (processingMode == kPMxi)
              hFractionElectrons->Fit(fElectronFraction, "N", "", electronFractionThresholdForFitting, lowFittingBoundElectronFraction);
            else if (processingMode == kPMdistance || processingMode == kPMjT) {
              // No special treatment of el for R and jT (seems to be the best option), i.e. do nothing here
            }
            else
              hFractionElectrons->Fit(fElectronFraction, "N", "", lowFittingBoundElectronFraction, electronFractionThresholdForFitting);
          }
          
        }
      }
      
      if ((regularisation > 0 && iter == 0) || (regularisation <= 0 && iter == 1)) {
        // Set the electron fraction to the negative pT -> A function will be used
        // to evaluate the electron fraction for each bin above the threshold
        if((mode == kPMxi && (slice + 1) <= electronFractionThresholdBinForFitting) || 
           (mode != kPMxi && (slice + 1) >= electronFractionThresholdBinForFitting)) {
          // In case of no regularisation, mathFit has no information about the fraction of other x bins.
          // Thus, the electron fraction is evaluated and set here. For the case w/ regularisation,
          // just "-pT" (or "-z" or "-xi") is set and the electron fraction will be evaluated during the fit.
          
          electronFixingUsed = kTRUE;
          
          Double_t xCoordinate = isPtMode ? (binsPt[slice + 1] + binsPt[slice]) / 2.
                                          :  hYieldPt->GetXaxis()->GetBinCenter(slice + 1);
          Double_t fixElectronFraction;
          if (regularisation <= 0) {
            // For the case w/o regularisation and with TOF patching, undo the TOF patching as discussed above
            fixElectronFraction = fElectronFraction->Eval(xCoordinate);
            
            if (applyTOFpatching) {
              // NOTE: All yields have inverseBinWidth factor included -> It doesn't matter, if it is really for all yields because
              // it just drops out then.
              const Double_t yieldTOFpions   = hYieldTOFPions->GetBinContent(slice + 1);
              const Double_t yieldTOFkaons   = hYieldTOFKaons->GetBinContent(slice + 1);
              const Double_t yieldTOFprotons = hYieldTOFProtons->GetBinContent(slice + 1);
              
              // No TOF electrons used at the moment, but histogram filled from contamination of TOF pions
              const Double_t yieldTOFelectrons = hYieldTOFElectrons->GetBinContent(slice + 1); 
              
              const Double_t yieldTOFtotal = yieldTOFelectrons + yieldTOFpions + yieldTOFkaons + yieldTOFprotons;
              
              const Double_t yieldTPCtotal = hYieldTPConlyTotal->GetBinContent(slice + 1);
              
              fixElectronFraction = UndoTOFpatchingForFraction(fixElectronFraction, yieldTPCtotal, yieldTOFelectrons, yieldTOFtotal);
            }
          }
          else
            fixElectronFraction = -xCoordinate;
          
          if (regularisation <= 0) {
            fixElectronFraction = TMath::Min(1.0, fixElectronFraction);
            fixElectronFraction = TMath::Max(0.0, fixElectronFraction);
          }
          
          gausParamsPi[3] = fixElectronFraction;
          lowParLimitsPi[3] = fixElectronFraction;
          upParLimitsPi[3] = fixElectronFraction;
          
          gausParamsEl[3] = fixElectronFraction;
          lowParLimitsEl[3] = fixElectronFraction;
          upParLimitsEl[3] = fixElectronFraction;
          
          gausParamsKa[3] = fixElectronFraction;
          lowParLimitsKa[3] = fixElectronFraction;
          upParLimitsKa[3] = fixElectronFraction;
          
          gausParamsPr[3] = fixElectronFraction;
          lowParLimitsPr[3] = fixElectronFraction;
          upParLimitsPr[3] = fixElectronFraction;
          
          stepSize[3] = 0.0;
          
          gausParamsSimultaneousFit[3] = fixElectronFraction;
          lowParLimitsSimultaneousFit[3] = fixElectronFraction;
          upParLimitsSimultaneousFit[3] = fixElectronFraction;
          
          stepSizeSimultaneousFit[3] = 0.0;
        }   
        
        
        // Set muon fraction equal to (some modified) electron fraction above some threshold, which should be a reasonable approximation:
        // Fixed muon fraction < 0 does this job within the fitting functions
        if(!isPtMode || (slice + 1) >= muonFractionThresholdBinForFitting) {
          // "Abuse" the muon fraction to forward the pT (z,xi), which can then be used to get some modified electron fraction
          const Double_t fixedValue = isPtMode ? -(binsPt[slice] + binsPt[slice + 1]) / 2.
                                               : -hYieldPt->GetXaxis()->GetBinCenter(slice + 1);
          
          gausParamsPi[4] = fixedValue;
          lowParLimitsPi[4] = fixedValue;
          upParLimitsPi[4] = fixedValue;
          
          gausParamsEl[4] = fixedValue;
          lowParLimitsEl[4] = fixedValue;
          upParLimitsEl[4] = fixedValue;
          
          gausParamsKa[4] = fixedValue;
          lowParLimitsKa[4] = fixedValue;
          upParLimitsKa[4] = fixedValue;
          
          gausParamsPr[4] = fixedValue;
          lowParLimitsPr[4] = fixedValue;
          upParLimitsPr[4] = fixedValue;
        
          stepSize[4] = 0.;
          
          gausParamsSimultaneousFit[4] = fixedValue;
          lowParLimitsSimultaneousFit[4] = fixedValue;
          upParLimitsSimultaneousFit[4] = fixedValue;
          
          stepSizeSimultaneousFit[4] = 0.0;
        }
      }
      
      // iter 0 used for initialisation
      if (regularisation > 0 && iter == 0) {
        const Int_t offset = currXbin * mathFit->GetNumParametersPerXbin();
        for (Int_t i = 0; i < mathFit->GetNumParametersPerXbin(); i++) {
          gausParamsSimultaneousFitRegularised[offset + i] = gausParamsSimultaneousFit[i];
          lowParLimitsSimultaneousFitRegularised[offset + i] = lowParLimitsSimultaneousFit[i];
          upParLimitsSimultaneousFitRegularised[offset + i] = upParLimitsSimultaneousFit[i];
          stepSizeSimultaneousFitRegularised[offset + i] = stepSizeSimultaneousFit[i];
        }
        
        // Store error of special templates in array
        fitRegularisedSpecTemplateFractionErrorKa[currXbin] = specTemplateFractionErrorKa;
        fitRegularisedSpecTemplateFractionErrorPr[currXbin] = specTemplateFractionErrorPr;
      }
      
      
      if (iter == 1) {
        // The parameters are only used for fitMethod < 2. Thus, they can be set for these methods,
        // although a different method is used
        totalDeltaPion->SetParameters(gausParamsPi);
        totalDeltaElectron->SetParameters(gausParamsEl);
        totalDeltaKaon->SetParameters(gausParamsKa);
        totalDeltaProton->SetParameters(gausParamsPr);
      }
      
      const TString binInfo = isPtMode ? Form("%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                       : Form("%.2f_%s_%.2f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                              modeShortName[mode].Data(), hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
      
      const TString binInfoTitle = isPtMode ? Form("%.2f < Pt <%.2f", binsPt[slice], binsPt[slice + 1])
                                            : Form("%.2f < %s < %.2f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                                   modeShortName[mode].Data(), 
                                                   hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
      
      const TString fitFuncSuffix = isPtMode ? Form("%.3f_Pt_%.3f", binsPt[slice], binsPt[slice + 1])
                                             : Form("%.3f_%s_%.3f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                                    modeShortName[mode].Data(), 
                                                    hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
      
      if (iter == 1) {
        for (Int_t species = 0; species < 4; species++) {
          cSingleFit[slice][species] = new TCanvas(Form("cSingleFit_%s_%s", binInfo.Data(), speciesLabel[species].Data()), 
                                                  Form("single fit for %s (%s)", binInfoTitle.Data(), speciesLabel[species].Data()),
                                                  1366, 768);
          cSingleFit[slice][species]->Divide(1, 2, 0.01, 0.);
          cSingleFit[slice][species]->GetPad(1)->SetRightMargin(0.001);
          cSingleFit[slice][species]->GetPad(2)->SetRightMargin(0.001);
          cSingleFit[slice][species]->GetPad(1)->SetTopMargin(0.001);
          cSingleFit[slice][species]->GetPad(2)->SetTopMargin(0.01);
          cSingleFit[slice][species]->GetPad(1)->SetBottomMargin(0.01);
          
          cSingleFit[slice][species]->GetPad(1)->SetGridx(kTRUE);
          cSingleFit[slice][species]->GetPad(2)->SetGridx(kTRUE);
          cSingleFit[slice][species]->GetPad(1)->SetGridy(kTRUE);
          cSingleFit[slice][species]->GetPad(2)->SetGridy(kTRUE);
          
          cSingleFit[slice][species]->GetPad(1)->SetLogy(kTRUE);
          cSingleFit[slice][species]->GetPad(1)->SetLogx(kTRUE);
          cSingleFit[slice][species]->GetPad(2)->SetLogx(kTRUE);
        }
      }
    
      // Problem: For p < 0.5 GeV/c, the fractions cannot be simply taken from the parameters because
      // not all entries of the histogram are inside the considered range.
      // Also: Small deviations of summed fractions from one if taking the fractions from different Delta species histos.
      // Therefore: Add up the integrals of the individual fits (\Delta species) and take the fraction of the sum
      Double_t integralTotal = 0;
      Double_t integralPions = 0, integralKaons = 0, integralProtons = 0, integralElectrons = 0, integralMuons = 0;
      
      Double_t integralPionsDeltaPion = 0, integralPionsDeltaElectron = 0, integralPionsDeltaKaon = 0, integralPionsDeltaProton = 0;
      Double_t integralElectronsDeltaPion = 0, integralElectronsDeltaElectron = 0, integralElectronsDeltaKaon = 0, 
              integralElectronsDeltaProton = 0;
      Double_t integralKaonsDeltaPion = 0, integralKaonsDeltaElectron = 0, integralKaonsDeltaKaon = 0, integralKaonsDeltaProton = 0;
      Double_t integralProtonsDeltaPion = 0, integralProtonsDeltaElectron = 0, integralProtonsDeltaKaon = 0, 
              integralProtonsDeltaProton = 0;
      Double_t integralMuonsDeltaPion = 0, integralMuonsDeltaElectron = 0, integralMuonsDeltaKaon = 0, integralMuonsDeltaProton = 0;
      
      /*
      Double_t integralErrorPions = 0, integralErrorKaons = 0, integralErrorProtons = 0, integralErrorElectrons = 0;
      
      Double_t integralErrorPionsDeltaPion = 0, integralErrorPionsDeltaElectron = 0, integralErrorPionsDeltaKaon = 0, 
              integralErrorPionsDeltaProton = 0;
      Double_t integralErrorElectronsDeltaPion = 0, integralErrorElectronsDeltaElectron = 0, integralErrorElectronsDeltaKaon = 0, 
              integralErrorElectronsDeltaProton = 0;
      Double_t integralErrorKaonsDeltaPion = 0, integralErrorKaonsDeltaElectron = 0, integralErrorKaonsDeltaKaon = 0, 
              integralErrorKaonsDeltaProton = 0;
      Double_t integralErrorProtonsDeltaPion = 0, integralErrorProtonsDeltaElectron = 0, integralErrorProtonsDeltaKaon = 0, 
              integralErrorProtonsDeltaProton = 0;
      
      Double_t integralErrorTotalDeltaPion = 0, integralErrorTotalDeltaElectron = 0, integralErrorTotalDeltaKaon = 0;   
      Double_t integralErrorTotalDeltaProton = 0;
      */
    
      Int_t errFlag = 0;
      
      // Reset temp arrays for next slice
      for (Int_t ind = 0; ind < nParUsed; ind++)
        parameterErrorsOut[ind] = 0;
      
      // Do not reset, if regularisation is on and the fit is done because the covariance matrix
      // will not be changed anymore in this case. On the other hand it will only be calculated once,
      // so resetting it would mean that is not available anymore.
      if (regularisation <= 0 || !regularisedFitDone) {
        for (Int_t i = 0; i < nParUsed; i++) {
          for (Int_t j = 0; j < nParUsed; j++) {
            covMatrix[i][j] = 0;
          }
        }
      }
      
      Double_t reducedChiSquare = -1;
      
      if (fitMethod == 2) {
        if (regularisation <= 0 && iter == 1)
          std::cout << "Fitting data simultaneously...." << std::endl << std::endl;
         
        // Add ref histos in initialisation step (w/ reg) or in the only loop (w/o reg)
        if ((regularisation > 0 && iter == 0) || (regularisation <= 0 && iter == 1)) {
          
          if (regularisation <= 0)
            mathFit->ClearRefHistos();
          
          mathFit->AddRefHisto(hGenDeltaPiForPiProj);
          mathFit->AddRefHisto(hGenDeltaPiForKaProj);
          mathFit->AddRefHisto(hGenDeltaPiForPrProj);
          mathFit->AddRefHisto(hGenDeltaPiForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaPiForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaKaForPiProj);
          mathFit->AddRefHisto(hGenDeltaKaForKaProj);
          mathFit->AddRefHisto(hGenDeltaKaForPrProj);
          mathFit->AddRefHisto(hGenDeltaKaForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaKaForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaPrForPiProj);
          mathFit->AddRefHisto(hGenDeltaPrForKaProj);
          mathFit->AddRefHisto(hGenDeltaPrForPrProj);
          mathFit->AddRefHisto(hGenDeltaPrForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaPrForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaElForPiProj);
          mathFit->AddRefHisto(hGenDeltaElForKaProj);
          mathFit->AddRefHisto(hGenDeltaElForPrProj);
          mathFit->AddRefHisto(hGenDeltaElForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaElForMuProj);
          
          // In reg case, fill in the data for this bin and continue with the nex bin
          if (regularisation > 0) {
            TH1D* hDeltaSpecies[numSimultaneousFits] = { hDeltaPi[slice], hDeltaKa[slice], hDeltaPr[slice], hDeltaEl[slice] };
            
            for (Int_t i = 0; i < numSimultaneousFits; i++) {
              mathFit->InputData(hDeltaSpecies[i], currXbin, i, xLow, xUp, -1., kFALSE); 
            }
            
            currXbin++;
            continue;
          }
        }
        
        if (regularisation > 0 && iter == 1 && !regularisedFitDone) {
          std::cout << "Fitting data simultaneously with regularisation...." << std::endl << std::endl;

          errFlag =  errFlag | doSimultaneousFitRegularised(nParSimultaneousFitRegularised, gausParamsSimultaneousFitRegularised, 
                                                            parameterErrorsOutRegularised, &covMatrix[0][0],
                                                            stepSizeSimultaneousFitRegularised, 
                                                            lowParLimitsSimultaneousFitRegularised, 
                                                            upParLimitsSimultaneousFitRegularised, reducedChiSquare,
                                                            fitRegularisedSpecTemplateFractionErrorKa,
                                                            fitRegularisedSpecTemplateFractionErrorPr);
          if (errFlag != 0)
            std::cout << "errFlag " << errFlag << std::endl << std::endl;
          
          reducedChiSquareRegularisation = reducedChiSquare;
          
          // Since everything is fitted in one go, only do this for the first x bin
          // (more convenient to put the fitting in the x bin loop in order to intialise
          // the parameters in the same way they are initialised for the fit w/o regularisation.
          regularisedFitDone = kTRUE;
        }
        
        if (regularisation > 0 && iter == 1) {
          // To allow for an identical processing, just forward the parameter results for the current xBin to the
          // array used by the standard simultaneous fit. The rest of the code is then the same for regularisation on and off
          
          for (Int_t i = 0; i < mathFit->GetNumParametersPerXbin(); i++) {
            const Int_t iRegularised = i + currXbin * mathFit->GetNumParametersPerXbin();
            
            gausParamsSimultaneousFit[i] = gausParamsSimultaneousFitRegularised[iRegularised];
            parameterErrorsOut[i]        = parameterErrorsOutRegularised[iRegularised];
          }
          
          // Same reducedChiSquare for all bins, since only one fit
          reducedChiSquare = reducedChiSquareRegularisation;
          
          
          // Also clear reference histos and load those for the current bin
          mathFit->ClearRefHistos();
          
          mathFit->AddRefHisto(hGenDeltaPiForPiProj);
          mathFit->AddRefHisto(hGenDeltaPiForKaProj);
          mathFit->AddRefHisto(hGenDeltaPiForPrProj);
          mathFit->AddRefHisto(hGenDeltaPiForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaPiForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaKaForPiProj);
          mathFit->AddRefHisto(hGenDeltaKaForKaProj);
          mathFit->AddRefHisto(hGenDeltaKaForPrProj);
          mathFit->AddRefHisto(hGenDeltaKaForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaKaForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaPrForPiProj);
          mathFit->AddRefHisto(hGenDeltaPrForKaProj);
          mathFit->AddRefHisto(hGenDeltaPrForPrProj);
          mathFit->AddRefHisto(hGenDeltaPrForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaPrForMuProj);
          
          mathFit->AddRefHisto(hGenDeltaElForPiProj);
          mathFit->AddRefHisto(hGenDeltaElForKaProj);
          mathFit->AddRefHisto(hGenDeltaElForPrProj);
          mathFit->AddRefHisto(hGenDeltaElForElProj);
          if (takeIntoAccountMuons)
            mathFit->AddRefHisto(hGenDeltaElForMuProj);
        }
      
      
        if (regularisation <= 0) {
          TH1D* hDeltaSpecies[numSimultaneousFits] = { hDeltaPi[slice], hDeltaKa[slice], hDeltaPr[slice], hDeltaEl[slice] };
          errFlag = errFlag | 
                    doSimultaneousFit(hDeltaSpecies, xLow, xUp, nParSimultaneousFit, gausParamsSimultaneousFit, parameterErrorsOut, 
                                      &covMatrix[0][0], stepSizeSimultaneousFit, lowParLimitsSimultaneousFit,
                                      upParLimitsSimultaneousFit, reducedChiSquare, specTemplateFractionErrorKa, specTemplateFractionErrorPr);
        }
        
        // Forward parameters to single fits
        for (Int_t parIndex = 0; parIndex < nPar; parIndex++) {
          // Fractions
          if (parIndex <= 4) {
            totalDeltaPion->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
            totalDeltaPion->SetParError(parIndex, parameterErrorsOut[parIndex]);
            
            totalDeltaKaon->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
            totalDeltaKaon->SetParError(parIndex, parameterErrorsOut[parIndex]);
            
            totalDeltaProton->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
            totalDeltaProton->SetParError(parIndex, parameterErrorsOut[parIndex]);
            
            totalDeltaElectron->SetParameter(parIndex, gausParamsSimultaneousFit[parIndex]);
            totalDeltaElectron->SetParError(parIndex, parameterErrorsOut[parIndex]);
          }
          // Total yield
          else if (parIndex == 5) {
            totalDeltaPion->SetParameter(parIndex, totalYield);
            totalDeltaPion->SetParError(parIndex, 0);
            
            totalDeltaKaon->SetParameter(parIndex, totalYield);
            totalDeltaKaon->SetParError(parIndex, 0);
            
            totalDeltaProton->SetParameter(parIndex, totalYield);
            totalDeltaProton->SetParError(parIndex, 0);
            
            totalDeltaElectron->SetParameter(parIndex, totalYield);
            totalDeltaElectron->SetParError(parIndex, 0);
          }
          // Hist shifts
          else {
            totalDeltaPion->SetParameter(parIndex, 0);
            totalDeltaPion->SetParError(parIndex, 0);
            
            totalDeltaKaon->SetParameter(parIndex, 0);
            totalDeltaKaon->SetParError(parIndex, 0);
            
            totalDeltaProton->SetParameter(parIndex, 0);
            totalDeltaProton->SetParError(parIndex, 0);
            
            totalDeltaElectron->SetParameter(parIndex, 0);
            totalDeltaElectron->SetParError(parIndex, 0);
          }
        }
        
        const Bool_t useRegularisation = regularisation > 0;
      
        // Plot single fits
        
        Int_t binLow = -1;
        Int_t binHigh = -1;
        
        // DeltaPions
        cSingleFit[slice][2]->cd(1);
        
        hDeltaPi[slice]->SetTitle("");
        hDeltaPi[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaPi[slice], binLow, binHigh);
        hDeltaPi[slice]->Draw("e");
        
        fitFuncTotalDeltaPion[slice] = (TF1*)totalDeltaPion->Clone(Form("Fit_Total_DeltaPion_%s", fitFuncSuffix.Data()));
        
        hDeltaPiFitQA[slice] = GetFitQAhisto(hDeltaPi[slice], fitFuncTotalDeltaPion[slice], Form("hDeltaPiFitQA_%d", slice));
        
        hDeltaPi[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaPion[slice]);
        fitFuncTotalDeltaPion[slice]->Draw("same");   
        
        Double_t* parametersOut = &totalDeltaPion->GetParameters()[0];
        
        hGenDeltaPiForPiProj->Scale(parametersOut[5] * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0)));
        shiftHist(hGenDeltaPiForPiProj, parametersOut[6], useRegularisation);
        hGenDeltaPiForPiProj->Draw("same");
        
        hGenDeltaPiForKaProj->Scale(parametersOut[5] * parametersOut[1]);
        shiftHist(hGenDeltaPiForKaProj, parametersOut[7], useRegularisation);
        hGenDeltaPiForKaProj->Draw("same");
        
        hGenDeltaPiForPrProj->Scale(parametersOut[5] * parametersOut[2]);
        shiftHist(hGenDeltaPiForPrProj, parametersOut[8], useRegularisation);
        hGenDeltaPiForPrProj->Draw("same");
        
        hGenDeltaPiForElProj->Scale(parametersOut[5] * parametersOut[3]);
        shiftHist(hGenDeltaPiForElProj, parametersOut[9], useRegularisation);
        hGenDeltaPiForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaPiForMuProj->Scale(parametersOut[5] * parametersOut[4]);
          shiftHist(hGenDeltaPiForMuProj, parametersOut[10], useRegularisation);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaPiForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaPiMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPiMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaPiMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPiMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaPi[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][2]->cd(2);
        hDeltaPiFitQA[slice]->Draw("e");    
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 3, reducedChiSquare);
        
       // TMatrixDSym covMatrixPi(nParUsed, &covMatrix[0][0]);   
       
        setFractionsAndYields(slice, inverseBinWidth, kPi, parametersOut, parameterErrorsOut, hFractionPions,
                              hFractionPionsDeltaPion, hFractionElectronsDeltaPion, hFractionKaonsDeltaPion,
                              hFractionProtonsDeltaPion, hFractionMuonsDeltaPion, hYieldPions, hYieldPionsDeltaPion, hYieldElectronsDeltaPion,
                              hYieldKaonsDeltaPion, hYieldProtonsDeltaPion, hYieldMuonsDeltaPion, normaliseResults);
        
        // Also set specific muon fractions and yields -> The deltaSpecies histos are not needed here: They will be set together with
        // the fraction and yields for all other species
        setFractionsAndYields(slice, inverseBinWidth, kMu, parametersOut, parameterErrorsOut, hFractionMuons,
                              0x0, 0x0, 0x0, 0x0, 0x0, hYieldMuons, 0x0, 0x0, 0x0, 0x0, 0x0, normaliseResults);
        
        
        // DeltaElectrons
        cSingleFit[slice][0]->cd(1);
        
        hDeltaEl[slice]->SetTitle("");
        hDeltaEl[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaEl[slice], binLow, binHigh);
        hDeltaEl[slice]->Draw("e");
        
        fitFuncTotalDeltaElectron[slice] = (TF1*)totalDeltaElectron->Clone(Form("Fit_Total_DeltaElectron_%s", fitFuncSuffix.Data()));
        
        hDeltaElFitQA[slice] = GetFitQAhisto(hDeltaEl[slice], fitFuncTotalDeltaElectron[slice], Form("hDeltaElFitQA_%d", slice));
        
        hDeltaEl[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaElectron[slice]);
        fitFuncTotalDeltaElectron[slice]->Draw("same");  
        
        parametersOut = &totalDeltaElectron->GetParameters()[0];
        
        hGenDeltaElForPiProj->Scale(parametersOut[5] * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0)));
        shiftHist(hGenDeltaElForPiProj, parametersOut[6], useRegularisation);
        hGenDeltaElForPiProj->Draw("same");
        
        hGenDeltaElForKaProj->Scale(parametersOut[5] * parametersOut[1]);
        shiftHist(hGenDeltaElForKaProj, parametersOut[7], useRegularisation);
        hGenDeltaElForKaProj->Draw("same");
        
        hGenDeltaElForPrProj->Scale(parametersOut[5] * parametersOut[2]);
        shiftHist(hGenDeltaElForPrProj, parametersOut[8], useRegularisation);
        hGenDeltaElForPrProj->Draw("same");
        
        hGenDeltaElForElProj->Scale(parametersOut[5] * parametersOut[3]);
        shiftHist(hGenDeltaElForElProj, parametersOut[9], useRegularisation);
        hGenDeltaElForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaElForMuProj->Scale(parametersOut[5] * parametersOut[4]);
          shiftHist(hGenDeltaElForMuProj, parametersOut[10], useRegularisation);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaElForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaElMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaElMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaElMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaElMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaEl[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][0]->cd(2);
        hDeltaElFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 1, reducedChiSquare);
        
        //TMatrixDSym covMatrixEl(nParUsed, &covMatrix[0][0]);
        
        setFractionsAndYields(slice, inverseBinWidth, kEl, parametersOut, parameterErrorsOut, hFractionElectrons,
                              hFractionPionsDeltaElectron, hFractionElectronsDeltaElectron, hFractionKaonsDeltaElectron,
                              hFractionProtonsDeltaElectron, hFractionMuonsDeltaElectron, hYieldElectrons, hYieldPionsDeltaElectron, 
                              hYieldElectronsDeltaElectron, hYieldKaonsDeltaElectron, hYieldProtonsDeltaElectron, hYieldMuonsDeltaElectron, 
                              normaliseResults);
        
        
        
        // DeltaKaons 
        cSingleFit[slice][1]->cd(1);
        
        hDeltaKa[slice]->SetTitle("");
        hDeltaKa[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaKa[slice], binLow, binHigh);
        hDeltaKa[slice]->Draw("e");
        
        fitFuncTotalDeltaKaon[slice] = (TF1*)totalDeltaKaon->Clone(Form("Fit_Total_DeltaKaon_%s", fitFuncSuffix.Data()));
        
        hDeltaKaFitQA[slice] = GetFitQAhisto(hDeltaKa[slice], fitFuncTotalDeltaKaon[slice], Form("hDeltaKaFitQA_%d", slice));
        
        hDeltaKa[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaKaon[slice]);
        fitFuncTotalDeltaKaon[slice]->Draw("same");  
        
        parametersOut = &totalDeltaKaon->GetParameters()[0];
        
        hGenDeltaKaForPiProj->Scale(parametersOut[5] * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0)));
        shiftHist(hGenDeltaKaForPiProj, parametersOut[6], useRegularisation);
        hGenDeltaKaForPiProj->Draw("same");
        
        hGenDeltaKaForKaProj->Scale(parametersOut[5] * parametersOut[1]);
        shiftHist(hGenDeltaKaForKaProj, parametersOut[7], useRegularisation);
        hGenDeltaKaForKaProj->Draw("same");
        
        hGenDeltaKaForPrProj->Scale(parametersOut[5] * parametersOut[2]);
        shiftHist(hGenDeltaKaForPrProj, parametersOut[8], useRegularisation);
        hGenDeltaKaForPrProj->Draw("same");
        
        hGenDeltaKaForElProj->Scale(parametersOut[5] * parametersOut[3]);
        shiftHist(hGenDeltaKaForElProj, parametersOut[9], useRegularisation);
        hGenDeltaKaForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaKaForMuProj->Scale(parametersOut[5] * parametersOut[4]);
          shiftHist(hGenDeltaKaForMuProj, parametersOut[10], useRegularisation);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaKaForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaKaMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaKaMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaKaMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaKaMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaKa[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][1]->cd(2);
        hDeltaKaFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 2, reducedChiSquare);
        
        //TMatrixDSym covMatrixKa(nParUsed, &covMatrix[0][0]);
        
        setFractionsAndYields(slice, inverseBinWidth, kKa, parametersOut, parameterErrorsOut, hFractionKaons,
                              hFractionPionsDeltaKaon, hFractionElectronsDeltaKaon, hFractionKaonsDeltaKaon, hFractionProtonsDeltaKaon,
                              hFractionMuonsDeltaKaon, hYieldKaons, hYieldPionsDeltaKaon, hYieldElectronsDeltaKaon, hYieldKaonsDeltaKaon,
                              hYieldProtonsDeltaKaon, hYieldMuonsDeltaKaon, normaliseResults);
        
        
        
        // DeltaProtons
        cSingleFit[slice][3]->cd(1);
        
        hDeltaPr[slice]->SetTitle("");
        hDeltaPr[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaPr[slice], binLow, binHigh);
        hDeltaPr[slice]->Draw("e");
        
        fitFuncTotalDeltaProton[slice] = (TF1*)totalDeltaProton->Clone(Form("Fit_Total_DeltaProton_%s", fitFuncSuffix.Data()));
        
        hDeltaPrFitQA[slice] = GetFitQAhisto(hDeltaPr[slice], fitFuncTotalDeltaProton[slice], Form("hDeltaPrFitQA_%d", slice));
        
        hDeltaPr[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaProton[slice]);
        
        fitFuncTotalDeltaProton[slice]->Draw("same");  
        
        parametersOut = &totalDeltaProton->GetParameters()[0];
        
        hGenDeltaPrForPiProj->Scale(parametersOut[5] * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0)));
        shiftHist(hGenDeltaPrForPiProj, parametersOut[6], useRegularisation);
        hGenDeltaPrForPiProj->Draw("same");
        
        hGenDeltaPrForKaProj->Scale(parametersOut[5] * parametersOut[1]);
        shiftHist(hGenDeltaPrForKaProj, parametersOut[7], useRegularisation);
        hGenDeltaPrForKaProj->Draw("same");
        
        hGenDeltaPrForPrProj->Scale(parametersOut[5] * parametersOut[2]);
        shiftHist(hGenDeltaPrForPrProj, parametersOut[8], useRegularisation);
        hGenDeltaPrForPrProj->Draw("same");
        
        hGenDeltaPrForElProj->Scale(parametersOut[5] * parametersOut[3]);
        shiftHist(hGenDeltaPrForElProj, parametersOut[9], useRegularisation);
        hGenDeltaPrForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaPrForMuProj->Scale(parametersOut[5] * parametersOut[4]);
          shiftHist(hGenDeltaPrForMuProj, parametersOut[10], useRegularisation);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaPrForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaPrMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPrMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaPrMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPrMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaPr[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][3]->cd(2);
        hDeltaPrFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 4, reducedChiSquare);
        
        //TMatrixDSym covMatrixPr(nParUsed, &covMatrix[0][0]);
        
        Double_t normalisationFactor = 1.0;
        normalisationFactor = setFractionsAndYields(slice, inverseBinWidth, kPr, parametersOut, parameterErrorsOut, 
                                                    hFractionProtons, hFractionPionsDeltaProton, hFractionElectronsDeltaProton, 
                                                    hFractionKaonsDeltaProton, hFractionProtonsDeltaProton, hFractionMuonsDeltaProton, 
                                                    hYieldProtons, hYieldPionsDeltaProton, hYieldElectronsDeltaProton,
                                                    hYieldKaonsDeltaProton, hYieldProtonsDeltaProton, hYieldMuonsDeltaProton, 
                                                    normaliseResults);
        
        // Fractions are the same for all plots -> just take deltaPion as default
        Double_t sumFractions = hFractionPionsDeltaPion->GetBinContent(slice + 1) + 
            hFractionElectronsDeltaPion->GetBinContent(slice + 1) + (takeIntoAccountMuons ? hFractionMuonsDeltaPion->GetBinContent(slice + 1) : 0.) +
            hFractionKaonsDeltaPion->GetBinContent(slice + 1) + hFractionProtonsDeltaPion->GetBinContent(slice + 1);
        
        hFractionSummed->SetBinContent(slice + 1, sumFractions);
        hFractionSummed->SetBinError(slice + 1, 
                                    TMath::Sqrt(TMath::Power(hFractionPionsDeltaPion->GetBinError(slice + 1), 2) + 
                                                TMath::Power(hFractionElectronsDeltaPion->GetBinError(slice + 1), 2) +
                                                (takeIntoAccountMuons ? TMath::Power(hFractionMuonsDeltaPion->GetBinError(slice + 1), 2) : 0.) +
                                                TMath::Power(hFractionKaonsDeltaPion->GetBinError(slice + 1), 2) +
                                                TMath::Power(hFractionProtonsDeltaPion->GetBinError(slice + 1), 2)));
          
        for (Int_t species = 0; species < 4; species++) {
          cSingleFit[slice][species]->Modified();
          cSingleFit[slice][species]->Update();
        }
        
        
        // Compute the to-pi ratios with proper error for the current slice
        // NOTE: error and covariance matrix are already scaled for the simultaneous fit
        // by mathFit (it was checked that all (i.e. also off-diagonal) matrix elements grow by fScaleError^2
        // NOTE 2: take the fractions and error from the histogram (takes correct muon and electrons fractions with errors set manually 
        // in case of fixed fraction; the parameters are fixed, so the elements of the off-diagonal elements of the covariance matrix 
        // remain zero!). The fractions are then also scaled to sum up to 1 (but correction factor usually close to unity).
        // The covariance matrix is NOT scaled like this. Therefore, scale the matrix elements accordingly.
        // If the normalisation is not done for the fractions, then this factor is unity by construction.
        // NOTE 3: Fixed parameters are already taken into account by mathFit, i.e. they APPEAR in the covariance matrix, but the
        // corresponding elements are zero as it should be
        
        
        
        Double_t covMatrixElementToPiForEl = 0.;
        Double_t covMatrixElementToPiForMu = 0.;
        Double_t covMatrixElementToPiForKa = 0.;
        Double_t covMatrixElementToPiForPr = 0.;
        
        // Get the correct covariance matrix elements and apply the normalisation factor
        Int_t parOffset = 0;
        
        // In case of regularisation, there is an offset with respect to the current slice
        if (useRegularisation)
          parOffset = currXbin * numParamsPerXbin;
        
        covMatrixElementToPiForEl = covMatrix[3 + parOffset][0 + parOffset] * normalisationFactor * normalisationFactor;
        covMatrixElementToPiForMu = covMatrix[4 + parOffset][0 + parOffset] * normalisationFactor * normalisationFactor;
        covMatrixElementToPiForKa = covMatrix[1 + parOffset][0 + parOffset] * normalisationFactor * normalisationFactor;
        covMatrixElementToPiForPr = covMatrix[2 + parOffset][0 + parOffset] * normalisationFactor * normalisationFactor;
        
        Double_t ratio = -999.;
        Double_t ratioError = 999.;
        Double_t currFractionSpecies = 0.;
        Double_t currFractionPions = 0.;
        Double_t currFractionErrorSpecies = 0.;
        Double_t currFractionErrorPions = 0.;
        Double_t covMatrixElementAB = 0.; // NOTE that there is only one covariance matrix (simultaneous fit!)
        
        currFractionPions = hFractionPions->GetBinContent(slice + 1);
        currFractionErrorPions = hFractionPions->GetBinError(slice + 1);
        
        // NOTE: Even in case of regularisation, when fractions of different bins become correlated, this does NOT change
        // the formula. Only the covariance matrix element for the considered fraction in the SAME slice needs to be taken
        // into account. Explanation: f = f(fracA_slice, fracB_slice), this means that \dell f / \dell fracA_slice+-1 = 0 (etc.).
        // So, the formula is the same, although the correlation between different slices is contained in the covariance matrix.
        
        // el-to-pi ratio
        currFractionSpecies = hFractionElectrons->GetBinContent(slice + 1);
        currFractionErrorSpecies = hFractionElectrons->GetBinError(slice + 1);
        covMatrixElementAB = covMatrixElementToPiForEl;
        
        GetRatioWithCorrelatedError(currFractionSpecies, currFractionPions, currFractionErrorSpecies, currFractionErrorPions,
                                    covMatrixElementAB, ratio, ratioError);
        
        hRatioToPiElectrons->SetBinContent(slice + 1, ratio);
        hRatioToPiElectrons->SetBinError(slice + 1, ratioError);
        
        // mu-to-pi ratio
        currFractionSpecies = hFractionMuons->GetBinContent(slice + 1);
        currFractionErrorSpecies = hFractionMuons->GetBinError(slice + 1);
        covMatrixElementAB = covMatrixElementToPiForMu;
        
        GetRatioWithCorrelatedError(currFractionSpecies, currFractionPions, currFractionErrorSpecies, currFractionErrorPions,
                                    covMatrixElementAB, ratio, ratioError);
        
        hRatioToPiMuons->SetBinContent(slice + 1, ratio);
        hRatioToPiMuons->SetBinError(slice + 1, ratioError);
        
        
        // K-to-pi ratio
        currFractionSpecies = hFractionKaons->GetBinContent(slice + 1);
        currFractionErrorSpecies = hFractionKaons->GetBinError(slice + 1);
        covMatrixElementAB = covMatrixElementToPiForKa;
        
        GetRatioWithCorrelatedError(currFractionSpecies, currFractionPions, currFractionErrorSpecies, currFractionErrorPions,
                                    covMatrixElementAB, ratio, ratioError);
        
        hRatioToPiKaons->SetBinContent(slice + 1, ratio);
        hRatioToPiKaons->SetBinError(slice + 1, ratioError);
        
        
        // p-to-pi ratio
        currFractionSpecies = hFractionProtons->GetBinContent(slice + 1);
        currFractionErrorSpecies = hFractionProtons->GetBinError(slice + 1);
        covMatrixElementAB = covMatrixElementToPiForPr;
        
        GetRatioWithCorrelatedError(currFractionSpecies, currFractionPions, currFractionErrorSpecies, currFractionErrorPions,
                                    covMatrixElementAB, ratio, ratioError);
        
        hRatioToPiProtons->SetBinContent(slice + 1, ratio);
        hRatioToPiProtons->SetBinError(slice + 1, ratioError);
        
        /*
        for (Int_t i = 0; i < nParUsed; i++) {
          for (Int_t j = 0; j < nParUsed; j++) {
            printf("\t%e", covMatrix[i][j]);
          }
          printf("\n");
        }
        */
        
        currXbin++;
      }
      //_____________________________________________________________________
      // Other methods without simultaneous fitting
      else {
        Int_t binLow = -1;
        Int_t binHigh = -1;
        
        // DeltaPions 
        
        std::cout << "Fitting deltaPion...." << std::endl << std::endl;
          
        cSingleFit[slice][2]->cd(1);
        
        mathFit->ClearRefHistos();
        mathFit->AddRefHisto(hGenDeltaPiForPiProj);
        mathFit->AddRefHisto(hGenDeltaPiForKaProj);
        mathFit->AddRefHisto(hGenDeltaPiForPrProj);
        mathFit->AddRefHisto(hGenDeltaPiForElProj);
        if (takeIntoAccountMuons)
          mathFit->AddRefHisto(hGenDeltaPiForMuProj);
        
        errFlag = errFlag |
                  doFit(hDeltaPi[slice], xLow, xUp, nPar, gausParamsPi, parameterErrorsOut, &covMatrix[0][0],
                        stepSize, lowParLimitsPi, upParLimitsPi, totalDeltaPion, reducedChiSquare,
                        specTemplateFractionErrorKaDeltaPi, specTemplateFractionErrorPrDeltaPi);
        
        hDeltaPi[slice]->SetTitle("");
        hDeltaPi[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaPi[slice], binLow, binHigh);
        hDeltaPi[slice]->Draw("e");
        
        fitFuncTotalDeltaPion[slice] = (TF1*)totalDeltaPion->Clone(Form("Fit_Total_DeltaPion_%s", fitFuncSuffix.Data()));
        
        hDeltaPiFitQA[slice] = GetFitQAhisto(hDeltaPi[slice], fitFuncTotalDeltaPion[slice], Form("hDeltaPiFitQA_%d", slice));
        
        hDeltaPi[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaPion[slice]);
        fitFuncTotalDeltaPion[slice]->Draw("same");   
        
        Double_t* parametersOut = &gausParamsPi[0];
        
        hGenDeltaPiForPiProj->Scale(gausParamsPi[5] * (gausParamsPi[0] + (muonContamination ? gausParamsPi[3] : 0)));
        shiftHist(hGenDeltaPiForPiProj, gausParamsPi[6]);
        hGenDeltaPiForPiProj->Draw("same");
        
        hGenDeltaPiForKaProj->Scale(gausParamsPi[5] * gausParamsPi[1]);
        shiftHist(hGenDeltaPiForKaProj, gausParamsPi[7]);
        hGenDeltaPiForKaProj->Draw("same");
        
        hGenDeltaPiForPrProj->Scale(gausParamsPi[5] * gausParamsPi[2]);
        shiftHist(hGenDeltaPiForPrProj, gausParamsPi[8]);
        hGenDeltaPiForPrProj->Draw("same");
        
        hGenDeltaPiForElProj->Scale(gausParamsPi[5] * gausParamsPi[3]);
        shiftHist(hGenDeltaPiForElProj, gausParamsPi[9]);
        hGenDeltaPiForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaPiForMuProj->Scale(gausParamsPi[5] * gausParamsPi[4]);
          shiftHist(hGenDeltaPiForMuProj, gausParamsPi[10]);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaPiForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaPiMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPiMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaPiMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPiMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaPi[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][2]->cd(2);
        hDeltaPiFitQA[slice]->Draw("e");    
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 3, reducedChiSquare);
        
        TMatrixDSym covMatrixPi(nParUsed, &covMatrix[0][0]);    
      
        if (fitMethod == 1)  {
          // Histos are normalised => expression equals integral
          integralPions = allDeltaPion * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0)); 
          integralTotal += integralPions;
          if (takeIntoAccountMuons) {
            integralMuons = allDeltaPion * parametersOut[4];
            integralTotal += integralMuons;
          }
          
          /*
          integralErrorTotalDeltaPion = getErrorOfTotalIntegral(covMatrixPi) * allDeltaPion;
          
          integralErrorPions = getErrorOfPionIntegral(covMatrixPi) * allDeltaPion;
          */
          
          integralPionsDeltaPion = integralPions;
          
          // Compare comment above
          integralElectronsDeltaPion = allDeltaPion * parametersOut[3];
          integralKaonsDeltaPion = allDeltaPion * parametersOut[1];
          integralProtonsDeltaPion = allDeltaPion * parametersOut[2];
          integralMuonsDeltaPion = allDeltaPion * parametersOut[4];
          
          /*
          integralErrorPionsDeltaPion = integralErrorPions;
          
          integralErrorElectronsDeltaPion = allDeltaPion * parameterErrorsOut[3];
          integralErrorKaonsDeltaPion = allDeltaPion * parameterErrorsOut[1];
          integralErrorProtonsDeltaPion = allDeltaPion * parameterErrorsOut[2];
          */
        }
        else  {
          setFractionsAndYields(slice, inverseBinWidth, kPi, parametersOut, parameterErrorsOut, hFractionPions,
                                hFractionPionsDeltaPion, hFractionElectronsDeltaPion, hFractionKaonsDeltaPion,
                                hFractionProtonsDeltaPion, hFractionMuonsDeltaPion, hYieldPions, hYieldPionsDeltaPion, hYieldElectronsDeltaPion,
                                hYieldKaonsDeltaPion, hYieldProtonsDeltaPion, hYieldMuonsDeltaPion);
          
          // Also set specific muon fractions and yields -> The deltaSpecies histos are not needed here: They will be set together with
          // the fraction and yields for all other species
          setFractionsAndYields(slice, inverseBinWidth, kMu, parametersOut, parameterErrorsOut, hFractionMuons,
                                0x0, 0x0, 0x0, 0x0, 0x0, hYieldMuons, 0x0, 0x0, 0x0, 0x0, 0x0);
        }
        
        
        std::cout << std::endl << std::endl;
        

        // DeltaElectrons
        
        std::cout << "Fitting deltaElectron...." << std::endl << std::endl;
        
        cSingleFit[slice][0]->cd(1);
        
        mathFit->ClearRefHistos();
        mathFit->AddRefHisto(hGenDeltaElForPiProj);
        mathFit->AddRefHisto(hGenDeltaElForKaProj);
        mathFit->AddRefHisto(hGenDeltaElForPrProj);
        mathFit->AddRefHisto(hGenDeltaElForElProj);
        if (takeIntoAccountMuons)
          mathFit->AddRefHisto(hGenDeltaElForMuProj);
        
        errFlag = errFlag |
                  doFit(hDeltaEl[slice], xLow, xUp, nPar, gausParamsEl, parameterErrorsOut, &covMatrix[0][0],
                        stepSize, lowParLimitsEl, upParLimitsEl, totalDeltaElectron, reducedChiSquare,
                        specTemplateFractionErrorKaDeltaEl, specTemplateFractionErrorPrDeltaEl);
        
        hDeltaEl[slice]->SetTitle("");
        hDeltaEl[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaEl[slice], binLow, binHigh);
        hDeltaEl[slice]->Draw("e");
        
        fitFuncTotalDeltaElectron[slice] = (TF1*)totalDeltaElectron->Clone(Form("Fit_Total_DeltaElectron_%s", fitFuncSuffix.Data()));
        
        hDeltaElFitQA[slice] = GetFitQAhisto(hDeltaEl[slice], fitFuncTotalDeltaElectron[slice], Form("hDeltaElFitQA_%d", slice));
        
        hDeltaEl[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaElectron[slice]);
        fitFuncTotalDeltaElectron[slice]->Draw("same");  
        
        parametersOut = &gausParamsEl[0];
        
        hGenDeltaElForPiProj->Scale(gausParamsEl[5] * (gausParamsEl[0] + (muonContamination ? gausParamsEl[3] : 0)));
        shiftHist(hGenDeltaElForPiProj, gausParamsEl[6]);
        hGenDeltaElForPiProj->Draw("same");
        
        hGenDeltaElForKaProj->Scale(gausParamsEl[5] * gausParamsEl[1]);
        shiftHist(hGenDeltaElForKaProj, gausParamsEl[7]);
        hGenDeltaElForKaProj->Draw("same");
        
        hGenDeltaElForPrProj->Scale(gausParamsEl[5] * gausParamsEl[2]);
        shiftHist(hGenDeltaElForPrProj, gausParamsEl[8]);
        hGenDeltaElForPrProj->Draw("same");
        
        hGenDeltaElForElProj->Scale(gausParamsEl[5] * gausParamsEl[3]);
        shiftHist(hGenDeltaElForElProj, gausParamsEl[9]);
        hGenDeltaElForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaElForMuProj->Scale(gausParamsEl[5] * gausParamsEl[4]);
          shiftHist(hGenDeltaElForMuProj, gausParamsEl[10]);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaElForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaElMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaElMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaElMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaElMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaEl[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][0]->cd(2);
        hDeltaElFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 1, reducedChiSquare);
        
        TMatrixDSym covMatrixEl(nParUsed, &covMatrix[0][0]);
        
        if (fitMethod == 1)  {
          integralElectrons = allDeltaElectron * parametersOut[3]; // Histos are normalised => expression equals integral
          integralTotal += integralElectrons;
          
          /*                                                                   
          integralErrorTotalDeltaElectron = getErrorOfTotalIntegral(covMatrixEl) * allDeltaElectron;
          
          integralErrorElectrons = allDeltaElectron * parameterErrorsOut[3];
          */
                          
          // Factor 2 in case of takeIntoAccountMuons will be applied below
          integralElectronsDeltaElectron = integralElectrons;

          // Compare comment above
          integralPionsDeltaElectron = allDeltaElectron * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0));
          integralKaonsDeltaElectron = allDeltaElectron * parametersOut[1];
          integralProtonsDeltaElectron = allDeltaElectron * parametersOut[2];
          integralMuonsDeltaElectron = allDeltaElectron * parametersOut[4];
          
          
          /*
          integralErrorElectronsDeltaElectron = integralErrorElectrons;
          
          integralErrorPionsDeltaElectron = getErrorOfPionIntegral(covMatrixEl) * allDeltaElectron;
          integralErrorKaonsDeltaElectron = allDeltaElectron * parameterErrorsOut[1];
          integralErrorProtonsDeltaElectron = allDeltaElectron * parameterErrorsOut[2];
          */
        }
        else  {
          setFractionsAndYields(slice, inverseBinWidth, kEl, parametersOut, parameterErrorsOut, hFractionElectrons,
                                hFractionPionsDeltaElectron, hFractionElectronsDeltaElectron, hFractionKaonsDeltaElectron,
                                hFractionProtonsDeltaElectron, hFractionMuonsDeltaElectron, hYieldElectrons, hYieldPionsDeltaElectron, 
                                hYieldElectronsDeltaElectron, hYieldKaonsDeltaElectron, hYieldProtonsDeltaElectron, hYieldMuonsDeltaElectron);
        }
        
        std::cout << std::endl << std::endl;
        
        // DeltaKaons 
        
        std::cout << "Fitting deltaKaon...." << std::endl << std::endl;
        
        cSingleFit[slice][1]->cd(1);
        
        mathFit->ClearRefHistos();
        mathFit->AddRefHisto(hGenDeltaKaForPiProj);
        mathFit->AddRefHisto(hGenDeltaKaForKaProj);
        mathFit->AddRefHisto(hGenDeltaKaForPrProj);
        mathFit->AddRefHisto(hGenDeltaKaForElProj);
        if (takeIntoAccountMuons)
          mathFit->AddRefHisto(hGenDeltaKaForMuProj);
        
        errFlag = errFlag |
                  doFit(hDeltaKa[slice], xLow, xUp, nPar, gausParamsKa, parameterErrorsOut, &covMatrix[0][0],
                        stepSize, lowParLimitsKa, upParLimitsKa, totalDeltaKaon, reducedChiSquare,
                        specTemplateFractionErrorKaDeltaKa, specTemplateFractionErrorPrDeltaKa);
        
        hDeltaKa[slice]->SetTitle("");
        hDeltaKa[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaKa[slice], binLow, binHigh);
        hDeltaKa[slice]->Draw("e");
        
        fitFuncTotalDeltaKaon[slice] = (TF1*)totalDeltaKaon->Clone(Form("Fit_Total_DeltaKaon_%s", fitFuncSuffix.Data()));
        
        hDeltaKaFitQA[slice] = GetFitQAhisto(hDeltaKa[slice], fitFuncTotalDeltaKaon[slice], Form("hDeltaKaFitQA_%d", slice));
        
        hDeltaKa[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaKaon[slice]);
        fitFuncTotalDeltaKaon[slice]->Draw("same");  
        
        parametersOut = &gausParamsKa[0];
        
        hGenDeltaKaForPiProj->Scale(gausParamsKa[5] * (gausParamsKa[0] + (muonContamination ? gausParamsKa[3] : 0)));
        shiftHist(hGenDeltaKaForPiProj, gausParamsKa[6]);
        hGenDeltaKaForPiProj->Draw("same");
        
        hGenDeltaKaForKaProj->Scale(gausParamsKa[5] * gausParamsKa[1]);
        shiftHist(hGenDeltaKaForKaProj, gausParamsKa[7]);
        hGenDeltaKaForKaProj->Draw("same");
        
        hGenDeltaKaForPrProj->Scale(gausParamsKa[5] * gausParamsKa[2]);
        shiftHist(hGenDeltaKaForPrProj, gausParamsKa[8]);
        hGenDeltaKaForPrProj->Draw("same");
        
        hGenDeltaKaForElProj->Scale(gausParamsKa[5] * gausParamsKa[3]);
        shiftHist(hGenDeltaKaForElProj, gausParamsKa[9]);
        hGenDeltaKaForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaKaForMuProj->Scale(gausParamsKa[5] * gausParamsKa[4]);
          shiftHist(hGenDeltaKaForMuProj, gausParamsKa[10]);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaKaForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaKaMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaKaMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaKaMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaKaMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaKa[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][1]->cd(2);
        hDeltaKaFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 2, reducedChiSquare);
        
        TMatrixDSym covMatrixKa(nParUsed, &covMatrix[0][0]);
        
        if (fitMethod == 1)  {
          integralKaons = allDeltaKaon * parametersOut[1]; // Histos are normalised => expression equals integral
          integralTotal += integralKaons;
          /*
          integralErrorTotalDeltaKaon = getErrorOfTotalIntegral(covMatrixKa) * allDeltaKaon;
                                                                      
          integralErrorKaons = allDeltaKaon * parameterErrorsOut[1];
          */
          
          
          integralKaonsDeltaKaon = integralKaons;
          
          // Compare comment above
          integralPionsDeltaKaon = allDeltaKaon * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0));
          integralElectronsDeltaKaon = allDeltaKaon * parametersOut[3];
          integralProtonsDeltaKaon = allDeltaKaon * parametersOut[2];
          integralMuonsDeltaKaon = allDeltaKaon * parametersOut[4];
          
          /*
          integralErrorKaonsDeltaKaon = integralErrorKaons;                                                            
                                  
          integralErrorPionsDeltaKaon = getErrorOfPionIntegral(covMatrixKa) * allDeltaKaon;
          integralErrorElectronsDeltaKaon = allDeltaKaon * parameterErrorsOut[3];
          integralErrorProtonsDeltaKaon = allDeltaKaon * parameterErrorsOut[2];
          */
        }
        else  {
          setFractionsAndYields(slice, inverseBinWidth, kKa, parametersOut, parameterErrorsOut, hFractionKaons,
                                hFractionPionsDeltaKaon, hFractionElectronsDeltaKaon, hFractionKaonsDeltaKaon, hFractionProtonsDeltaKaon,
                                hFractionMuonsDeltaKaon, hYieldKaons, hYieldPionsDeltaKaon, hYieldElectronsDeltaKaon, hYieldKaonsDeltaKaon, 
                                hYieldProtonsDeltaKaon, hYieldMuonsDeltaKaon);
        }
        
        std::cout << std::endl << std::endl;
        
        
        // DeltaProtons
        
        std::cout << "Fitting deltaProton...." << std::endl << std::endl;
        
        cSingleFit[slice][3]->cd(1);
        
        mathFit->ClearRefHistos();
        mathFit->AddRefHisto(hGenDeltaPrForPiProj);
        mathFit->AddRefHisto(hGenDeltaPrForKaProj);
        mathFit->AddRefHisto(hGenDeltaPrForPrProj);
        mathFit->AddRefHisto(hGenDeltaPrForElProj);
        if (takeIntoAccountMuons)
          mathFit->AddRefHisto(hGenDeltaPrForMuProj);
        
        errFlag = errFlag | 
                  doFit(hDeltaPr[slice], xLow, xUp, nPar, gausParamsPr, parameterErrorsOut, &covMatrix[0][0],
                        stepSize, lowParLimitsPr, upParLimitsPr, totalDeltaProton, reducedChiSquare,
                        specTemplateFractionErrorKaDeltaPr, specTemplateFractionErrorPrDeltaPr);
        
        hDeltaPr[slice]->SetTitle("");
        hDeltaPr[slice]->GetYaxis()->SetTitle("Entries");
        SetReasonableXaxisRange(hDeltaPr[slice], binLow, binHigh);
        hDeltaPr[slice]->Draw("e");
        
        fitFuncTotalDeltaProton[slice] = (TF1*)totalDeltaProton->Clone(Form("Fit_Total_DeltaProton_%s", fitFuncSuffix.Data()));
        
        hDeltaPrFitQA[slice] = GetFitQAhisto(hDeltaPr[slice], fitFuncTotalDeltaProton[slice], Form("hDeltaPrFitQA_%d", slice));
        
        hDeltaPr[slice]->GetListOfFunctions()->Add(fitFuncTotalDeltaProton[slice]);
        
        fitFuncTotalDeltaProton[slice]->Draw("same");  
        
        parametersOut = &gausParamsPr[0];
        
        hGenDeltaPrForPiProj->Scale(gausParamsPr[5] * (gausParamsPr[0] + (muonContamination ? gausParamsPr[3] : 0)));
        shiftHist(hGenDeltaPrForPiProj, gausParamsPr[6]);
        hGenDeltaPrForPiProj->Draw("same");
        
        hGenDeltaPrForKaProj->Scale(gausParamsPr[5] * gausParamsPr[1]);
        shiftHist(hGenDeltaPrForKaProj, gausParamsPr[7]);
        hGenDeltaPrForKaProj->Draw("same");
        
        hGenDeltaPrForPrProj->Scale(gausParamsPr[5] * gausParamsPr[2]);
        shiftHist(hGenDeltaPrForPrProj, gausParamsPr[8]);
        hGenDeltaPrForPrProj->Draw("same");
        
        hGenDeltaPrForElProj->Scale(gausParamsPr[5] * gausParamsPr[3]);
        shiftHist(hGenDeltaPrForElProj, gausParamsPr[9]);
        hGenDeltaPrForElProj->Draw("same");
        
        if (takeIntoAccountMuons) {
          hGenDeltaPrForMuProj->Scale(gausParamsPr[5] * gausParamsPr[4]);
          shiftHist(hGenDeltaPrForMuProj, gausParamsPr[10]);
          if (muonFractionHandling != kNoMuons)
            hGenDeltaPrForMuProj->Draw("same");
        }
        
        if (plotIdentifiedSpectra) {
          for (Int_t species = 0; species < 5; species++) 
            hDeltaPrMC[slice][species]->Draw("same");
          
          // Draw histo for sum of MC muons and pions
          TH1D* hMCmuonsAndPions = new TH1D(*hDeltaPrMC[slice][kPi - 1]);
          hMCmuonsAndPions->Add(hDeltaPrMC[slice][kMu - 1]);
          hMCmuonsAndPions->SetLineColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetMarkerColor(getLineColor(kMuPlusPi));
          hMCmuonsAndPions->SetName(Form("%s_muonsAdded", hDeltaPrMC[slice][kPi - 1]->GetName()));
          hMCmuonsAndPions->Draw("same");
        }
        
        hDeltaPr[slice]->Draw("esame");
        
        legend->Draw();
        
        cSingleFit[slice][3]->cd(2);
        hDeltaPrFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 4, reducedChiSquare);
        
        TMatrixDSym covMatrixPr(nParUsed, &covMatrix[0][0]);
        
        if (fitMethod == 1)  {
          integralProtons = allDeltaProton * parametersOut[2]; // Histos are normalised => expression equals integral
          integralTotal += integralProtons;
          /*
          integralErrorTotalDeltaProton = getErrorOfTotalIntegral(covMatrixPr) * allDeltaProton;
          
          integralErrorProtons = allDeltaProton * parameterErrorsOut[2];
          */                                       
          
          integralProtonsDeltaProton = integralProtons;
          
          // Compare comment above
          integralPionsDeltaProton = allDeltaProton * (parametersOut[0] + (muonContamination ? parametersOut[3] : 0));
          integralElectronsDeltaProton = allDeltaProton * parametersOut[3];
          integralKaonsDeltaProton = allDeltaProton * parametersOut[1];
          integralMuonsDeltaProton = allDeltaProton * parametersOut[4];
          
          
          /*
          integralErrorProtonsDeltaProton = integralErrorProtons;                                                            
                                  
          integralErrorPionsDeltaProton = getErrorOfPionIntegral(covMatrixPr) * allDeltaProton;
          integralErrorElectronsDeltaProton = allDeltaProton * parameterErrorsOut[3];
          integralErrorKaonsDeltaProton = allDeltaProton * parameterErrorsOut[1];
          */
        }
        else  {
          setFractionsAndYields(slice, inverseBinWidth, kPr, parametersOut, parameterErrorsOut, hFractionProtons,
                                hFractionPionsDeltaProton, hFractionElectronsDeltaProton, hFractionKaonsDeltaProton,
                                hFractionProtonsDeltaProton, hFractionMuonsDeltaProton, hYieldProtons, hYieldPionsDeltaProton,
                                hYieldElectronsDeltaProton, hYieldKaonsDeltaProton, hYieldProtonsDeltaProton, hYieldMuonsDeltaProton);
        }
        
        std::cout << std::endl << std::endl;
      
        
        if (fitMethod == 1)  { 
          // Calculate fractions and yields for method 1
          if (integralTotal > 0)  {   
            
            Double_t sumOfParticles = 0;
            
            // Check fraction and yield determination for systematics
            // DeltaPion
            Double_t integralTotalDeltaPion = integralPionsDeltaPion + integralElectronsDeltaPion +
                                              (takeIntoAccountMuons ? integralMuonsDeltaPion : 0.) + 
                                              integralKaonsDeltaPion + integralProtonsDeltaPion;
            totalDeltaPion->GetParameters(parametersOut);
            
            Double_t pionFractionDeltaPion = saveDivide(integralPionsDeltaPion, integralTotalDeltaPion);
            Double_t pionFractionErrorDeltaPion = getErrorOfPionFraction(parametersOut, covMatrixPi);
            hFractionPionsDeltaPion->SetBinContent(slice + 1, pionFractionDeltaPion);
            hFractionPionsDeltaPion->SetBinError(slice + 1, pionFractionErrorDeltaPion);
            
            Double_t electronFractionDeltaPion = saveDivide(integralElectronsDeltaPion, integralTotalDeltaPion);
            Double_t electronFractionErrorDeltaPion = getErrorOfElectronFraction(parametersOut, covMatrixPi);
            hFractionElectronsDeltaPion->SetBinContent(slice + 1, electronFractionDeltaPion);
            hFractionElectronsDeltaPion->SetBinError(slice + 1, electronFractionErrorDeltaPion);
            
            Double_t kaonFractionDeltaPion = saveDivide(integralKaonsDeltaPion, integralTotalDeltaPion);
            Double_t kaonFractionErrorDeltaPion = getErrorOfKaonFraction(parametersOut, covMatrixPi);
            hFractionKaonsDeltaPion->SetBinContent(slice + 1, kaonFractionDeltaPion);
            hFractionKaonsDeltaPion->SetBinError(slice + 1, kaonFractionErrorDeltaPion);
            
            Double_t protonFractionDeltaPion = saveDivide(integralProtonsDeltaPion, integralTotalDeltaPion);
            Double_t protonFractionErrorDeltaPion = getErrorOfProtonFraction(parametersOut, covMatrixPi);
            hFractionProtonsDeltaPion->SetBinContent(slice + 1, protonFractionDeltaPion);
            hFractionProtonsDeltaPion->SetBinError(slice + 1, protonFractionErrorDeltaPion);
            
            Double_t muonFractionDeltaPion = saveDivide(integralMuonsDeltaPion, integralTotalDeltaPion);
            // TODO Error is anyway not implemented correctly. Just take electron error as an approximation
            Double_t muonFractionErrorDeltaPion = getErrorOfElectronFraction(parametersOut, covMatrixPi);
            hFractionMuonsDeltaPion->SetBinContent(slice + 1, muonFractionDeltaPion);
            hFractionMuonsDeltaPion->SetBinError(slice + 1, muonFractionErrorDeltaPion);
            
            sumOfParticles = inverseBinWidth * gausParamsPi[5]; 
            
            hYieldPionsDeltaPion->SetBinContent(slice + 1, sumOfParticles * hFractionPionsDeltaPion->GetBinContent(slice + 1));
            hYieldPionsDeltaPion->SetBinError(slice + 1, sumOfParticles * hFractionPionsDeltaPion->GetBinError(slice + 1));
            hYieldElectronsDeltaPion->SetBinContent(slice + 1, sumOfParticles * hFractionElectronsDeltaPion->GetBinContent(slice + 1));
            hYieldElectronsDeltaPion->SetBinError(slice + 1, sumOfParticles * hFractionElectronsDeltaPion->GetBinError(slice + 1));
            hYieldKaonsDeltaPion->SetBinContent(slice + 1, sumOfParticles * hFractionKaonsDeltaPion->GetBinContent(slice + 1));
            hYieldKaonsDeltaPion->SetBinError(slice + 1, sumOfParticles * hFractionKaonsDeltaPion->GetBinError(slice + 1));
            hYieldProtonsDeltaPion->SetBinContent(slice + 1, sumOfParticles * hFractionProtonsDeltaPion->GetBinContent(slice + 1));
            hYieldProtonsDeltaPion->SetBinError(slice + 1, sumOfParticles * hFractionProtonsDeltaPion->GetBinError(slice + 1));
            hYieldMuonsDeltaPion->SetBinContent(slice + 1, sumOfParticles * hFractionMuonsDeltaPion->GetBinContent(slice + 1));
            hYieldMuonsDeltaPion->SetBinError(slice + 1, sumOfParticles * hFractionMuonsDeltaPion->GetBinError(slice + 1));
            
            
            // DeltaElectron
            Double_t integralTotalDeltaElectron = integralPionsDeltaElectron + integralElectronsDeltaElectron + 
                                                  (takeIntoAccountMuons ? integralMuonsDeltaElectron : 0.) + 
                                                  integralKaonsDeltaElectron + integralProtonsDeltaElectron;
            totalDeltaElectron->GetParameters(parametersOut);
            
            Double_t pionFractionDeltaElectron = saveDivide(integralPionsDeltaElectron, integralTotalDeltaElectron);
            Double_t pionFractionErrorDeltaElectron = getErrorOfPionFraction(parametersOut, covMatrixEl);
            hFractionPionsDeltaElectron->SetBinContent(slice + 1, pionFractionDeltaElectron);
            hFractionPionsDeltaElectron->SetBinError(slice + 1, pionFractionErrorDeltaElectron);
            
            Double_t electronFractionDeltaElectron = saveDivide(integralElectronsDeltaElectron, integralTotalDeltaElectron);
            Double_t electronFractionErrorDeltaElectron = getErrorOfElectronFraction(parametersOut, covMatrixEl);
            hFractionElectronsDeltaElectron->SetBinContent(slice + 1, electronFractionDeltaElectron);
            hFractionElectronsDeltaElectron->SetBinError(slice + 1, electronFractionErrorDeltaElectron);
            
            Double_t kaonFractionDeltaElectron = saveDivide(integralKaonsDeltaElectron, integralTotalDeltaElectron);
            Double_t kaonFractionErrorDeltaElectron = getErrorOfKaonFraction(parametersOut, covMatrixEl);
            hFractionKaonsDeltaElectron->SetBinContent(slice + 1, kaonFractionDeltaElectron);
            hFractionKaonsDeltaElectron->SetBinError(slice + 1, kaonFractionErrorDeltaElectron);
            
            Double_t protonFractionDeltaElectron = saveDivide(integralProtonsDeltaElectron, integralTotalDeltaElectron);
            Double_t protonFractionErrorDeltaElectron = getErrorOfProtonFraction(parametersOut, covMatrixEl);
            hFractionProtonsDeltaElectron->SetBinContent(slice + 1, protonFractionDeltaElectron);
            hFractionProtonsDeltaElectron->SetBinError(slice + 1, protonFractionErrorDeltaElectron);
            
            Double_t muonFractionDeltaElectron = saveDivide(integralMuonsDeltaElectron, integralTotalDeltaElectron);
            // TODO Error is anyway not implemented correctly. Just take electron error as an approximation
            Double_t muonFractionErrorDeltaElectron = getErrorOfElectronFraction(parametersOut, covMatrixEl);
            hFractionMuonsDeltaElectron->SetBinContent(slice + 1, muonFractionDeltaElectron);
            hFractionMuonsDeltaElectron->SetBinError(slice + 1, muonFractionErrorDeltaElectron);
            
            sumOfParticles = inverseBinWidth * gausParamsEl[5];
            
            hYieldPionsDeltaElectron->SetBinContent(slice + 1, sumOfParticles * hFractionPionsDeltaElectron->GetBinContent(slice + 1));
            hYieldPionsDeltaElectron->SetBinError(slice + 1, sumOfParticles * hFractionPionsDeltaElectron->GetBinError(slice + 1));
            hYieldElectronsDeltaElectron->SetBinContent(slice + 1, sumOfParticles * hFractionElectronsDeltaElectron->GetBinContent(slice + 1));
            hYieldElectronsDeltaElectron->SetBinError(slice + 1, sumOfParticles * hFractionElectronsDeltaElectron->GetBinError(slice + 1));
            hYieldKaonsDeltaElectron->SetBinContent(slice + 1, sumOfParticles * hFractionKaonsDeltaElectron->GetBinContent(slice + 1));
            hYieldKaonsDeltaElectron->SetBinError(slice + 1, sumOfParticles * hFractionKaonsDeltaElectron->GetBinError(slice + 1));
            hYieldProtonsDeltaElectron->SetBinContent(slice + 1, sumOfParticles * hFractionProtonsDeltaElectron->GetBinContent(slice + 1));
            hYieldProtonsDeltaElectron->SetBinError(slice + 1, sumOfParticles * hFractionProtonsDeltaElectron->GetBinError(slice + 1));
            hYieldMuonsDeltaElectron->SetBinContent(slice + 1, sumOfParticles * hFractionMuonsDeltaElectron->GetBinContent(slice + 1));
            hYieldMuonsDeltaElectron->SetBinError(slice + 1, sumOfParticles * hFractionMuonsDeltaElectron->GetBinError(slice + 1));
            
            
            // DeltaKaon
            Double_t integralTotalDeltaKaon = integralPionsDeltaKaon + integralElectronsDeltaKaon +
                                              (takeIntoAccountMuons ? integralMuonsDeltaKaon : 0.) + 
                                              integralKaonsDeltaKaon + integralProtonsDeltaKaon;
            totalDeltaKaon->GetParameters(parametersOut);
            
            Double_t pionFractionDeltaKaon = saveDivide(integralPionsDeltaKaon, integralTotalDeltaKaon);
            Double_t pionFractionErrorDeltaKaon = getErrorOfPionFraction(parametersOut, covMatrixKa);
            hFractionPionsDeltaKaon->SetBinContent(slice + 1, pionFractionDeltaKaon);
            hFractionPionsDeltaKaon->SetBinError(slice + 1, pionFractionErrorDeltaKaon);
            
            Double_t electronFractionDeltaKaon = saveDivide(integralElectronsDeltaKaon, integralTotalDeltaKaon);
            Double_t electronFractionErrorDeltaKaon = getErrorOfElectronFraction(parametersOut, covMatrixKa);
            hFractionElectronsDeltaKaon->SetBinContent(slice + 1, electronFractionDeltaKaon);
            hFractionElectronsDeltaKaon->SetBinError(slice + 1, electronFractionErrorDeltaKaon);
            
            Double_t kaonFractionDeltaKaon = saveDivide(integralKaonsDeltaKaon, integralTotalDeltaKaon);
            Double_t kaonFractionErrorDeltaKaon = getErrorOfKaonFraction(parametersOut, covMatrixKa);
            hFractionKaonsDeltaKaon->SetBinContent(slice + 1, kaonFractionDeltaKaon);
            hFractionKaonsDeltaKaon->SetBinError(slice + 1, kaonFractionErrorDeltaKaon);
            
            Double_t protonFractionDeltaKaon = saveDivide(integralProtonsDeltaKaon, integralTotalDeltaKaon);
            Double_t protonFractionErrorDeltaKaon = getErrorOfProtonFraction(parametersOut, covMatrixKa);
            hFractionProtonsDeltaKaon->SetBinContent(slice + 1, protonFractionDeltaKaon);
            hFractionProtonsDeltaKaon->SetBinError(slice + 1, protonFractionErrorDeltaKaon);
            
            Double_t muonFractionDeltaKaon = saveDivide(integralMuonsDeltaKaon, integralTotalDeltaKaon);
            // TODO Error is anyway not implemented correctly. Just take electron error as an approximation
            Double_t muonFractionErrorDeltaKaon = getErrorOfElectronFraction(parametersOut, covMatrixKa);
            hFractionMuonsDeltaKaon->SetBinContent(slice + 1, muonFractionDeltaKaon);
            hFractionMuonsDeltaKaon->SetBinError(slice + 1, muonFractionErrorDeltaKaon);
            
            sumOfParticles = inverseBinWidth * gausParamsKa[5];
            
            hYieldPionsDeltaKaon->SetBinContent(slice + 1, sumOfParticles * hFractionPionsDeltaKaon->GetBinContent(slice + 1));
            hYieldPionsDeltaKaon->SetBinError(slice + 1, sumOfParticles * hFractionPionsDeltaKaon->GetBinError(slice + 1));
            hYieldElectronsDeltaKaon->SetBinContent(slice + 1, sumOfParticles * hFractionElectronsDeltaKaon->GetBinContent(slice + 1));
            hYieldElectronsDeltaKaon->SetBinError(slice + 1, sumOfParticles * hFractionElectronsDeltaKaon->GetBinError(slice + 1));
            hYieldKaonsDeltaKaon->SetBinContent(slice + 1, sumOfParticles * hFractionKaonsDeltaKaon->GetBinContent(slice + 1));
            hYieldKaonsDeltaKaon->SetBinError(slice + 1, sumOfParticles * hFractionKaonsDeltaKaon->GetBinError(slice + 1));
            hYieldProtonsDeltaKaon->SetBinContent(slice + 1, sumOfParticles * hFractionProtonsDeltaKaon->GetBinContent(slice + 1));
            hYieldProtonsDeltaKaon->SetBinError(slice + 1, sumOfParticles * hFractionProtonsDeltaKaon->GetBinError(slice + 1));
            hYieldMuonsDeltaKaon->SetBinContent(slice + 1, sumOfParticles * hFractionMuonsDeltaKaon->GetBinContent(slice + 1));
            hYieldMuonsDeltaKaon->SetBinError(slice + 1, sumOfParticles * hFractionMuonsDeltaKaon->GetBinError(slice + 1));
            
            
            
            // DeltaProton
            Double_t integralTotalDeltaProton = integralPionsDeltaProton + integralElectronsDeltaProton +
                                                (takeIntoAccountMuons ? integralMuonsDeltaProton : 0.) + 
                                                integralKaonsDeltaProton + integralProtonsDeltaProton;
            totalDeltaProton->GetParameters(parametersOut);
            
            Double_t pionFractionDeltaProton = saveDivide(integralPionsDeltaProton, integralTotalDeltaProton);
            Double_t pionFractionErrorDeltaProton = getErrorOfPionFraction(parametersOut, covMatrixPr);
            hFractionPionsDeltaProton->SetBinContent(slice + 1, pionFractionDeltaProton);
            hFractionPionsDeltaProton->SetBinError(slice + 1, pionFractionErrorDeltaProton);
            
            Double_t electronFractionDeltaProton = saveDivide(integralElectronsDeltaProton, integralTotalDeltaProton);
            Double_t electronFractionErrorDeltaProton = getErrorOfElectronFraction(parametersOut, covMatrixPr);
            hFractionElectronsDeltaProton->SetBinContent(slice + 1, electronFractionDeltaProton);
            hFractionElectronsDeltaProton->SetBinError(slice + 1, electronFractionErrorDeltaProton);
            
            Double_t kaonFractionDeltaProton = saveDivide(integralKaonsDeltaProton, integralTotalDeltaProton);
            Double_t kaonFractionErrorDeltaProton = getErrorOfKaonFraction(parametersOut, covMatrixPr);
            hFractionKaonsDeltaProton->SetBinContent(slice + 1, kaonFractionDeltaProton);
            hFractionKaonsDeltaProton->SetBinError(slice + 1, kaonFractionErrorDeltaProton);
            
            Double_t protonFractionDeltaProton = saveDivide(integralProtonsDeltaProton, integralTotalDeltaProton);
            Double_t protonFractionErrorDeltaProton = getErrorOfProtonFraction(parametersOut, covMatrixPr);
            hFractionProtonsDeltaProton->SetBinContent(slice + 1, protonFractionDeltaProton);
            hFractionProtonsDeltaProton->SetBinError(slice + 1, protonFractionErrorDeltaProton);
            
            Double_t muonFractionDeltaProton = saveDivide(integralMuonsDeltaProton, integralTotalDeltaProton);
            // TODO Error is anyway not implemented correctly. Just take electron error as an approximation
            Double_t muonFractionErrorDeltaProton = getErrorOfElectronFraction(parametersOut, covMatrixPr);
            hFractionMuonsDeltaProton->SetBinContent(slice + 1, muonFractionDeltaProton);
            hFractionMuonsDeltaProton->SetBinError(slice + 1, muonFractionErrorDeltaProton);
            
            sumOfParticles = inverseBinWidth * gausParamsPr[5];
            
            hYieldPionsDeltaProton->SetBinContent(slice + 1, sumOfParticles * hFractionPionsDeltaProton->GetBinContent(slice + 1));
            hYieldPionsDeltaProton->SetBinError(slice + 1, sumOfParticles * hFractionPionsDeltaProton->GetBinError(slice + 1));
            hYieldElectronsDeltaProton->SetBinContent(slice + 1, sumOfParticles * hFractionElectronsDeltaProton->GetBinContent(slice + 1));
            hYieldElectronsDeltaProton->SetBinError(slice + 1, sumOfParticles * hFractionElectronsDeltaProton->GetBinError(slice + 1));
            hYieldKaonsDeltaProton->SetBinContent(slice + 1, sumOfParticles * hFractionKaonsDeltaProton->GetBinContent(slice + 1));
            hYieldKaonsDeltaProton->SetBinError(slice + 1, sumOfParticles * hFractionKaonsDeltaProton->GetBinError(slice + 1));
            hYieldProtonsDeltaProton->SetBinContent(slice + 1, sumOfParticles * hFractionProtonsDeltaProton->GetBinContent(slice + 1));
            hYieldProtonsDeltaProton->SetBinError(slice + 1, sumOfParticles * hFractionProtonsDeltaProton->GetBinError(slice + 1));
            hYieldMuonsDeltaProton->SetBinContent(slice + 1, sumOfParticles * hFractionMuonsDeltaProton->GetBinContent(slice + 1));
            hYieldMuonsDeltaProton->SetBinError(slice + 1, sumOfParticles * hFractionMuonsDeltaProton->GetBinError(slice + 1));
            
            
            
            // Take for XXXXfractionError the median of XXXXfractionErrorYYYY and do not take into account errors
            // with value zero, since the should correspond to a failed fit (but the other fits can still converge).
            // Same for the yields
            Double_t pionFraction = saveDivide(integralPions, integralTotal);
            Double_t errorsPions[4] = { pionFractionErrorDeltaPion, pionFractionErrorDeltaElectron, 
                pionFractionErrorDeltaKaon, pionFractionErrorDeltaProton };
            Double_t pionFractionError = getMedianOfNonZeros(errorsPions);
            
            Double_t electronFraction = saveDivide(integralElectrons, integralTotal);
            Double_t errorsElectrons[4] = { electronFractionErrorDeltaPion, electronFractionErrorDeltaElectron, 
                electronFractionErrorDeltaKaon, electronFractionErrorDeltaProton };
            Double_t electronFractionError = getMedianOfNonZeros(errorsElectrons);
            
            Double_t kaonFraction = saveDivide(integralKaons, integralTotal);
            Double_t errorsKaons[4] = { kaonFractionErrorDeltaPion, kaonFractionErrorDeltaElectron, 
                kaonFractionErrorDeltaKaon, kaonFractionErrorDeltaProton };
            Double_t kaonFractionError = getMedianOfNonZeros(errorsKaons);
            
            Double_t protonFraction = saveDivide(integralProtons, integralTotal);
            Double_t errorsProtons[4] = { protonFractionErrorDeltaPion, protonFractionErrorDeltaElectron, 
                protonFractionErrorDeltaKaon, protonFractionErrorDeltaProton };
            Double_t protonFractionError = getMedianOfNonZeros(errorsProtons);
            
            Double_t muonFraction = saveDivide(integralMuons, integralTotal);
            Double_t errorsMuons[4] = { muonFractionErrorDeltaPion, muonFractionErrorDeltaElectron, 
                muonFractionErrorDeltaKaon, muonFractionErrorDeltaProton };
            Double_t muonFractionError = getMedianOfNonZeros(errorsMuons);
            
            hFractionPions->SetBinContent(slice + 1, pionFraction);
            hFractionPions->SetBinError(slice + 1, pionFractionError);
            hFractionElectrons->SetBinContent(slice + 1, electronFraction);
            hFractionElectrons->SetBinError(slice + 1, electronFractionError);
            hFractionKaons->SetBinContent(slice + 1, kaonFraction);
            hFractionKaons->SetBinError(slice + 1, kaonFractionError);
            hFractionProtons->SetBinContent(slice + 1, protonFraction);
            hFractionProtons->SetBinError(slice + 1, protonFractionError);
            hFractionMuons->SetBinContent(slice + 1, muonFraction);
            hFractionMuons->SetBinError(slice + 1, muonFractionError);
            
            hFractionSummed->SetBinContent(slice + 1, pionFraction + electronFraction + (takeIntoAccountMuons ? muonFraction : 0.) +
                                                      kaonFraction + protonFraction);
            hFractionSummed->SetBinError(slice + 1, 
                                        TMath::Sqrt(TMath::Power(pionFractionError, 2) +
                                                    TMath::Power(electronFractionError, 2)  +
                                                    (takeIntoAccountMuons ? TMath::Power(muonFractionError, 2) : 0.) +
                                                    TMath::Power(kaonFractionError, 2) +
                                                    TMath::Power(protonFractionError, 2)));
            
            sumOfParticles = inverseBinWidth * integralTotal;
            
            hYieldPions->SetBinContent(slice + 1, sumOfParticles * hFractionPions->GetBinContent(slice + 1));
            hYieldPions->SetBinError(slice + 1, sumOfParticles * hFractionPions->GetBinError(slice + 1));
            hYieldElectrons->SetBinContent(slice + 1, sumOfParticles * hFractionElectrons->GetBinContent(slice + 1));
            hYieldElectrons->SetBinError(slice + 1, sumOfParticles * hFractionElectrons->GetBinError(slice + 1));
            hYieldKaons->SetBinContent(slice + 1, sumOfParticles * hFractionKaons->GetBinContent(slice + 1));
            hYieldKaons->SetBinError(slice + 1, sumOfParticles * hFractionKaons->GetBinError(slice + 1));
            hYieldProtons->SetBinContent(slice + 1, sumOfParticles * hFractionProtons->GetBinContent(slice + 1));
            hYieldProtons->SetBinError(slice + 1, sumOfParticles * hFractionProtons->GetBinError(slice + 1));
            hYieldMuons->SetBinContent(slice + 1, sumOfParticles * hFractionMuons->GetBinContent(slice + 1));
            hYieldMuons->SetBinError(slice + 1, sumOfParticles * hFractionMuons->GetBinError(slice + 1));
          }
        }
        else  {  
          Double_t SumFractionsDeltaElectron = hFractionPionsDeltaElectron->GetBinContent(slice + 1) + 
              hFractionElectronsDeltaElectron->GetBinContent(slice + 1) + 
              (takeIntoAccountMuons ? hFractionMuonsDeltaElectron->GetBinContent(slice + 1) : 0.) +
              hFractionKaonsDeltaElectron->GetBinContent(slice + 1) + hFractionProtonsDeltaElectron->GetBinContent(slice + 1);
          
          Double_t SumFractionsDeltaKaon = hFractionPionsDeltaKaon->GetBinContent(slice + 1) + 
              hFractionElectronsDeltaKaon->GetBinContent(slice + 1) +
              (takeIntoAccountMuons ? hFractionMuonsDeltaKaon->GetBinContent(slice + 1) : 0.) +
              hFractionKaonsDeltaKaon->GetBinContent(slice + 1) + hFractionProtonsDeltaKaon->GetBinContent(slice + 1);
          
          Double_t SumFractionsDeltaPion = hFractionPionsDeltaPion->GetBinContent(slice + 1) + 
              hFractionElectronsDeltaPion->GetBinContent(slice + 1) +
              (takeIntoAccountMuons ? hFractionMuonsDeltaPion->GetBinContent(slice + 1) : 0.) +
              hFractionKaonsDeltaPion->GetBinContent(slice + 1) + hFractionProtonsDeltaPion->GetBinContent(slice + 1);
          
          Double_t SumFractionsDeltaProton = hFractionPionsDeltaProton->GetBinContent(slice + 1) + 
              hFractionElectronsDeltaProton->GetBinContent(slice + 1) +
              (takeIntoAccountMuons ? hFractionMuonsDeltaProton->GetBinContent(slice + 1) : 0.) +
              hFractionKaonsDeltaProton->GetBinContent(slice + 1) + hFractionProtonsDeltaProton->GetBinContent(slice + 1);
          
          Double_t SumFractionsUsed = hFractionPionsDeltaPion->GetBinContent(slice + 1) + 
              hFractionElectronsDeltaElectron->GetBinContent(slice + 1) +
              (takeIntoAccountMuons ? hFractionMuonsDeltaPion->GetBinContent(slice + 1) : 0.) +
              hFractionKaonsDeltaKaon->GetBinContent(slice + 1) + hFractionProtonsDeltaProton->GetBinContent(slice + 1);
          
          hFractionSummed->SetBinContent(slice + 1, SumFractionsUsed);
          hFractionSummed->SetBinError(slice + 1, 
                                      TMath::Sqrt(TMath::Power(hFractionPionsDeltaPion->GetBinError(slice + 1), 2) + 
                                                  TMath::Power(hFractionElectronsDeltaElectron->GetBinError(slice + 1), 2) +
                                                  (takeIntoAccountMuons ? TMath::Power(hFractionMuonsDeltaPion->GetBinError(slice + 1),   
                                                                                        2) : 0.) +
                                                  TMath::Power(hFractionKaonsDeltaKaon->GetBinError(slice + 1), 2) +
                                                  TMath::Power(hFractionProtonsDeltaProton->GetBinError(slice + 1), 2)));
          
          
          std::cout << "Sum Fractions DeltaElectron: " << SumFractionsDeltaElectron;
          std::cout << (TMath::Abs(SumFractionsDeltaElectron - 1) >= 0.001 ? " WARNING: Deviation >= 0.001" : "") << std::endl;
          
          std::cout << "Sum Fractions DeltaKaon: " << SumFractionsDeltaKaon;
          std::cout << (TMath::Abs(SumFractionsDeltaKaon - 1) >= 0.001 ? " WARNING: Deviation >= 0.001" : "") << std::endl;
          
          std::cout << "Sum Fractions DeltaPion: " << SumFractionsDeltaPion;
          std::cout << (TMath::Abs(SumFractionsDeltaPion - 1) >= 0.001 ? " WARNING: Deviation >= 0.001" : "") << std::endl;
          
          std::cout << "Sum fractions DeltaProton: " << SumFractionsDeltaProton;
          std::cout << (TMath::Abs(SumFractionsDeltaProton - 1) >= 0.001 ? " WARNING: Deviation >= 0.001" : "") << std::endl;
          
          std::cout << "Sum fractions used: " << SumFractionsUsed;
          std::cout << (TMath::Abs(SumFractionsUsed - 1) >= 0.001 ? " WARNING: Deviation >= 0.001" : "") << std::endl;
        }
        
        for (Int_t species = 0; species < 4; species++) {
          cSingleFit[slice][species]->Modified();
          cSingleFit[slice][species]->Update();
        }
        
      
      }
      
      if (regularisation <= 0)
        std::cout << std::endl << std::endl;
      
      
      // MC results
      Double_t MCtotal = -1, MCelectrons = -1, MCkaons = -1, MCmuons = -1, MCpions = -1, MCprotons = -1;
      Double_t MCelectronsErr = 0, MCkaonsErr = 0, MCmuonsErr = 0, MCpionsErr = 0, MCprotonsErr = 0;
      
      MCelectrons = hMCdata->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 1, 1, MCelectronsErr) * inverseBinWidth;
      MCkaons     = hMCdata->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 2, 2, MCkaonsErr)     * inverseBinWidth;
      MCmuons     = hMCdata->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 3, 3, MCmuonsErr)     * inverseBinWidth;
      MCpions     = hMCdata->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 4, 4, MCpionsErr)     * inverseBinWidth;
      MCprotons   = hMCdata->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, 5, 5, MCprotonsErr)   * inverseBinWidth;

      MCelectronsErr *= inverseBinWidth;
      MCkaonsErr *= inverseBinWidth;
      MCmuonsErr *= inverseBinWidth;
      MCpionsErr *= inverseBinWidth;
      MCprotonsErr *= inverseBinWidth;
      
      MCtotal = MCelectrons + MCkaons + MCpions + MCprotons + MCmuons;
      
      if (MCtotal > 0)  {
        hYieldElectronsMC->SetBinContent(slice + 1, MCelectrons);
        hYieldElectronsMC->SetBinError(slice + 1, MCelectronsErr);
        
        hYieldMuonsMC->SetBinContent(slice + 1, MCmuons);
        hYieldMuonsMC->SetBinError(slice + 1, MCmuonsErr);
        
        hYieldKaonsMC->SetBinContent(slice + 1, MCkaons);
        hYieldKaonsMC->SetBinError(slice + 1, MCkaonsErr);
        
        hYieldPionsMC->SetBinContent(slice + 1, MCpions);
        hYieldPionsMC->SetBinError(slice + 1, MCpionsErr);
        
        hYieldProtonsMC->SetBinContent(slice + 1, MCprotons);
        hYieldProtonsMC->SetBinError(slice + 1, MCprotonsErr);
        
        hYieldSummedMC->SetBinContent(slice + 1, hYieldElectronsMC->GetBinContent(slice + 1) +
                                                hYieldKaonsMC->GetBinContent(slice + 1) +                           
                                                hYieldPionsMC->GetBinContent(slice + 1) +
                                                hYieldProtonsMC->GetBinContent(slice + 1) +
                                                hYieldMuonsMC->GetBinContent(slice + 1));
        hYieldSummedMC->SetBinError(slice + 1, TMath::Sqrt(TMath::Power(hYieldPionsMC->GetBinError(slice + 1), 2) + 
                                                          TMath::Power(hYieldElectronsMC->GetBinError(slice + 1), 2) +
                                                          TMath::Power(hYieldKaonsMC->GetBinError(slice + 1), 2) +
                                                          TMath::Power(hYieldProtonsMC->GetBinError(slice + 1), 2) +
                                                          TMath::Power(hYieldMuonsMC->GetBinError(slice + 1), 2)));
        
        // MCspecies and MCtotal are correlated. This can be taken into account via using the binomial error in the division
        hFractionElectronsMC->Divide(hYieldElectronsMC, hYieldSummedMC, 1., 1., "B");
        hFractionMuonsMC->Divide(hYieldMuonsMC, hYieldSummedMC, 1., 1., "B");
        hFractionKaonsMC->Divide(hYieldKaonsMC, hYieldSummedMC, 1., 1., "B");
        hFractionPionsMC->Divide(hYieldPionsMC, hYieldSummedMC, 1., 1., "B");
        hFractionProtonsMC->Divide(hYieldProtonsMC, hYieldSummedMC, 1., 1., "B");
      }
      
      // Save further results
      if (isFirstTimeToWrite) { // || (isXiMode && (sliceEnd - slice) % 18 == 0) || (!isXiMode && slice % 18 == 0)) {
        isFirstTimeToWrite = kFALSE;
        saveF->cd();  

        if (hFractionElectrons)
          hFractionElectrons->Write(0, TObject::kWriteDelete);
          
        if (hFractionKaons)
          hFractionKaons->Write(0, TObject::kWriteDelete);
        
        if (hFractionPions)
          hFractionPions->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtons)
          hFractionProtons->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuons)
          hFractionMuons->Write(0, TObject::kWriteDelete);
        
        if (hFractionSummed)
          hFractionSummed->Write(0, TObject::kWriteDelete);
          
          
        if (hFractionElectronsDeltaElectron)
          hFractionElectronsDeltaElectron->Write(0, TObject::kWriteDelete);
          
        if (hFractionKaonsDeltaElectron)
          hFractionKaonsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hFractionPionsDeltaElectron)
          hFractionPionsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtonsDeltaElectron)
          hFractionProtonsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuonsDeltaElectron)
          hFractionMuonsDeltaElectron->Write(0, TObject::kWriteDelete);
          
          
        if (hFractionElectronsDeltaPion)
          hFractionElectronsDeltaPion->Write(0, TObject::kWriteDelete);
          
        if (hFractionKaonsDeltaPion)
          hFractionKaonsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hFractionPionsDeltaPion)
          hFractionPionsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtonsDeltaPion)
          hFractionProtonsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuonsDeltaPion)
          hFractionMuonsDeltaPion->Write(0, TObject::kWriteDelete);
          
        
        if (hFractionElectronsDeltaKaon)
          hFractionElectronsDeltaKaon->Write(0, TObject::kWriteDelete);
          
        if (hFractionKaonsDeltaKaon)
          hFractionKaonsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hFractionPionsDeltaKaon)
          hFractionPionsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtonsDeltaKaon)
          hFractionProtonsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuonsDeltaKaon)
          hFractionMuonsDeltaKaon->Write(0, TObject::kWriteDelete);
          
        
        if (hFractionElectronsDeltaProton)
          hFractionElectronsDeltaProton->Write(0, TObject::kWriteDelete);
          
        if (hFractionKaonsDeltaProton)
          hFractionKaonsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hFractionPionsDeltaProton)
          hFractionPionsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtonsDeltaProton)
          hFractionProtonsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuonsDeltaProton)
          hFractionMuonsDeltaProton->Write(0, TObject::kWriteDelete);
          
          
        if (hFractionElectronsMC)
          hFractionElectronsMC->Write(0, TObject::kWriteDelete);
        
        if (hFractionKaonsMC)
          hFractionKaonsMC->Write(0, TObject::kWriteDelete);
        
        if (hFractionPionsMC)
          hFractionPionsMC->Write(0, TObject::kWriteDelete);
        
        if (hFractionMuonsMC)
          hFractionMuonsMC->Write(0, TObject::kWriteDelete);
        
        if (hFractionProtonsMC)
          hFractionProtonsMC->Write(0, TObject::kWriteDelete);
        
        
        
        
        if (hYieldElectrons)
          hYieldElectrons->Write(0, TObject::kWriteDelete);
          
        if (hYieldKaons)
          hYieldKaons->Write(0, TObject::kWriteDelete);
        
        if (hYieldPions)
          hYieldPions->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtons)
          hYieldProtons->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuons)
          hYieldMuons->Write(0, TObject::kWriteDelete);
          
          
        if (hYieldElectronsDeltaElectron)
          hYieldElectronsDeltaElectron->Write(0, TObject::kWriteDelete);
          
        if (hYieldKaonsDeltaElectron)
          hYieldKaonsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hYieldPionsDeltaElectron)
          hYieldPionsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtonsDeltaElectron)
          hYieldProtonsDeltaElectron->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuonsDeltaElectron)
          hYieldMuonsDeltaElectron->Write(0, TObject::kWriteDelete);
          
          
        if (hYieldElectronsDeltaPion)
          hYieldElectronsDeltaPion->Write(0, TObject::kWriteDelete);
          
        if (hYieldKaonsDeltaPion)
          hYieldKaonsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hYieldPionsDeltaPion)
          hYieldPionsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtonsDeltaPion)
          hYieldProtonsDeltaPion->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuonsDeltaPion)
          hYieldMuonsDeltaPion->Write(0, TObject::kWriteDelete);
          
        
        if (hYieldElectronsDeltaKaon)
          hYieldElectronsDeltaKaon->Write(0, TObject::kWriteDelete);
          
        if (hYieldKaonsDeltaKaon)
          hYieldKaonsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hYieldPionsDeltaKaon)
          hYieldPionsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtonsDeltaKaon)
          hYieldProtonsDeltaKaon->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuonsDeltaKaon)
          hYieldMuonsDeltaKaon->Write(0, TObject::kWriteDelete);
          
        
        if (hYieldElectronsDeltaProton)
          hYieldElectronsDeltaProton->Write(0, TObject::kWriteDelete);
          
        if (hYieldKaonsDeltaProton)
          hYieldKaonsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hYieldPionsDeltaProton)
          hYieldPionsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtonsDeltaProton)
          hYieldProtonsDeltaProton->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuonsDeltaProton)
          hYieldMuonsDeltaProton->Write(0, TObject::kWriteDelete);
          
          
        if (hYieldElectronsMC)
          hYieldElectronsMC->Write(0, TObject::kWriteDelete);
        
        if (hYieldKaonsMC)
          hYieldKaonsMC->Write(0, TObject::kWriteDelete);
        
        if (hYieldPionsMC)
          hYieldPionsMC->Write(0, TObject::kWriteDelete);
        
        if (hYieldMuonsMC)
          hYieldMuonsMC->Write(0, TObject::kWriteDelete);
        
        if (hYieldProtonsMC)
          hYieldProtonsMC->Write(0, TObject::kWriteDelete);
        
        if (hYieldSummedMC)
          hYieldSummedMC->Write(0, TObject::kWriteDelete);
      }
      
      TString saveDir = isPtMode ? Form("SingleFit_%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                 : Form("SingleFit_%.2f_%s_%.2f", hYieldPt->GetXaxis()->GetBinLowEdge(slice + 1), 
                                        modeShortName[mode].Data(), hYieldPt->GetXaxis()->GetBinUpEdge(slice + 1));
      saveF->mkdir(saveDir.Data());
      saveF->cd(saveDir.Data());
      
      for (Int_t species = 0; species < 4; species++) {
        if (cSingleFit[slice][species]) {
          cSingleFit[slice][species]->Write();
          delete cSingleFit[slice][species];
        }
      }
      
      if (hDeltaPi[slice])
        hDeltaPi[slice]->Write();
      
      if (hDeltaEl[slice])
        hDeltaEl[slice]->Write();
      
      if (hDeltaKa[slice])
        hDeltaKa[slice]->Write();
      
      if (hDeltaPr[slice])
        hDeltaPr[slice]->Write();
      
      
      if (hDeltaPiFitQA[slice])
        hDeltaPiFitQA[slice]->Write();
      delete hDeltaPiFitQA[slice];
      
      if (hDeltaElFitQA[slice])
        hDeltaElFitQA[slice]->Write();
      delete hDeltaElFitQA[slice];
      
      if (hDeltaKaFitQA[slice])
        hDeltaKaFitQA[slice]->Write();
      delete hDeltaKaFitQA[slice];
      
      if (hDeltaPrFitQA[slice])
        hDeltaPrFitQA[slice]->Write();
      delete hDeltaPrFitQA[slice];
      
      if (hGenDeltaElForElProj) 
        hGenDeltaElForElProj->Write();
      delete hGenDeltaElForElProj;
      
      if (hGenDeltaElForKaProj) 
        hGenDeltaElForKaProj->Write();
      delete hGenDeltaElForKaProj;
      
      if (hGenDeltaElForPiProj) 
        hGenDeltaElForPiProj->Write();
      delete hGenDeltaElForPiProj;
      
      if (hGenDeltaElForPrProj) 
        hGenDeltaElForPrProj->Write();
      delete hGenDeltaElForPrProj;
      
      if (hGenDeltaElForMuProj) 
        hGenDeltaElForMuProj->Write();
      delete hGenDeltaElForMuProj;
      
      //if (fitFuncTotalDeltaElectron[slice]) 
      //  fitFuncTotalDeltaElectron[slice]->Write();
      delete fitFuncTotalDeltaElectron[slice];
      
      if (hGenDeltaKaForElProj) 
        hGenDeltaKaForElProj->Write();
      delete hGenDeltaKaForElProj;
      
      if (hGenDeltaKaForKaProj) 
        hGenDeltaKaForKaProj->Write();
      delete hGenDeltaKaForKaProj;
      
      if (hGenDeltaKaForPiProj) 
        hGenDeltaKaForPiProj->Write();
      delete hGenDeltaKaForPiProj;
      
      if (hGenDeltaKaForPrProj) 
        hGenDeltaKaForPrProj->Write();
      delete hGenDeltaKaForPrProj;
      
      if (hGenDeltaKaForMuProj) 
        hGenDeltaKaForMuProj->Write();
      delete hGenDeltaKaForMuProj;
      
      //if (fitFuncTotalDeltaKaon[slice]) 
      //  fitFuncTotalDeltaKaon[slice]->Write();
      delete fitFuncTotalDeltaKaon[slice];
      
        
      if (hGenDeltaPiForElProj) 
        hGenDeltaPiForElProj->Write();
      delete hGenDeltaPiForElProj;
      
      if (hGenDeltaPiForKaProj) 
        hGenDeltaPiForKaProj->Write();
      delete hGenDeltaPiForKaProj;
      
      if (hGenDeltaPiForPiProj) 
        hGenDeltaPiForPiProj->Write();
      delete hGenDeltaPiForPiProj;
      
      if (hGenDeltaPiForPrProj) 
        hGenDeltaPiForPrProj->Write();
      delete hGenDeltaPiForPrProj;
      
      if (hGenDeltaPiForMuProj) 
        hGenDeltaPiForMuProj->Write();
      delete hGenDeltaPiForMuProj;
      
      //if (fitFuncTotalDeltaPion[slice]) 
      //  fitFuncTotalDeltaPion[slice]->Write();
      delete fitFuncTotalDeltaPion[slice];
      
      
      if (hGenDeltaPrForElProj) 
        hGenDeltaPrForElProj->Write();
      delete hGenDeltaPrForElProj;
      
      if (hGenDeltaPrForKaProj) 
        hGenDeltaPrForKaProj->Write();
      delete hGenDeltaPrForKaProj;
      
      if (hGenDeltaPrForPiProj) 
        hGenDeltaPrForPiProj->Write();
      delete hGenDeltaPrForPiProj;
      
      if (hGenDeltaPrForPrProj) 
        hGenDeltaPrForPrProj->Write();
      delete hGenDeltaPrForPrProj;
      
      if (hGenDeltaPrForMuProj) 
        hGenDeltaPrForMuProj->Write();
      delete hGenDeltaPrForMuProj;
      
      //if (fitFuncTotalDeltaProton[slice]) 
      //  fitFuncTotalDeltaProton[slice]->Write();
      delete fitFuncTotalDeltaProton[slice];
      
      delete totalDeltaElectron;
      delete totalDeltaKaon;
      delete totalDeltaPion;
      delete totalDeltaProton;
      
      delete legend;
      
      if (errFlag != 0)
        std::cout << "errFlag " << errFlag << std::endl << std::endl;
    }
  }
  
  if (applyTOFpatching) {
    // Total TOF yield is sum of pi, k, p
    hYieldTOFTotal = new TH1F(*hYieldTOFPions);
    hYieldTOFTotal->SetName("hYieldTOFTotal");
    hYieldTOFTotal->SetTitle("Sum");
    hYieldTOFTotal->SetLineColor(kBlack);
    hYieldTOFTotal->SetMarkerColor(kBlack);
    
    hYieldTOFTotal->Add(hYieldTOFKaons);
    hYieldTOFTotal->Add(hYieldTOFProtons);
    hYieldTOFTotal->Add(hYieldTOFElectrons);
  }
  
  
  // Calculate MC to-pi ratios -> In MC the yields are uncorrelated, so just divide the histos to get the correct result
  hRatioToPiElectronsMC->Divide(hYieldElectronsMC, hYieldPionsMC);
  hRatioToPiMuonsMC->Divide(hYieldMuonsMC, hYieldPionsMC);
  hRatioToPiKaonsMC->Divide(hYieldKaonsMC, hYieldPionsMC);
  hRatioToPiProtonsMC->Divide(hYieldProtonsMC, hYieldPionsMC);
  
  
  TCanvas* cFractions = drawFinalFractions(mode, pLow, pHigh, electronFixingUsed, hFractionPions, hFractionKaons, hFractionProtons, 
                                           hFractionElectrons, hFractionMuons, hFractionSummed, hFractionPionsMC, hFractionKaonsMC, 
                                           hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, plotIdentifiedSpectra);
  

  // Compare data points with MC
  TCanvas* cFractionComparisons = drawFractionComparisonToMC(mode, pLow, pHigh, hFractionPions, hFractionKaons, hFractionProtons,
                                                             hFractionElectrons, hFractionMuons, hFractionPionsMC, hFractionKaonsMC,
                                                             hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, 
                                                             hFractionComparisonPions, hFractionComparisonKaons, 
                                                             hFractionComparisonProtons, hFractionComparisonElectrons, 
                                                             hFractionComparisonMuons);
  
  // Normalise the yields
  normaliseYieldHist(hYieldPions, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldPionsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldPionsDeltaElectron, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldPionsDeltaPion, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldPionsDeltaKaon, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldPionsDeltaProton, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  normaliseYieldHist(hYieldElectrons, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldElectronsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldElectronsDeltaElectron, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldElectronsDeltaPion, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldElectronsDeltaKaon, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldElectronsDeltaProton, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  normaliseYieldHist(hYieldMuons, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldMuonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldMuonsDeltaElectron, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldMuonsDeltaPion, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldMuonsDeltaKaon, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldMuonsDeltaProton, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  normaliseYieldHist(hYieldKaons, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldKaonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldKaonsDeltaElectron, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldKaonsDeltaPion, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldKaonsDeltaKaon, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldKaonsDeltaProton, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  normaliseYieldHist(hYieldProtons, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldProtonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldProtonsDeltaElectron, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldProtonsDeltaPion, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldProtonsDeltaKaon, numEvents, deta, restrictJetPtAxis, numJetsRec);
  normaliseYieldHist(hYieldProtonsDeltaProton, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  normaliseYieldHist(hYieldSummedMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  
  if (applyTOFpatching) {
    normaliseYieldHist(hYieldTOFPions, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFKaons, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFProtons, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFElectrons, numEvents, deta, restrictJetPtAxis, numJetsRec);
    
    normaliseYieldHist(hYieldTOFTotal, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTPConlyTotal, numEvents, deta, restrictJetPtAxis, numJetsRec);
    
    normaliseYieldHist(hYieldTOFElectronsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFMuonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFKaonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFPionsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
    normaliseYieldHist(hYieldTOFProtonsMC, numEvents, deta, restrictJetPtAxis, numJetsRec);
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenYieldsPrimSpecies[i]) {
      Int_t color = kBlack;
      
      switch (i) {
        case AliPID::kElectron:
          color = getLineColor(kEl);
          break;
        case AliPID::kKaon:
          color = getLineColor(kKa);
          break;
        case AliPID::kMuon:
          color = getLineColor(kMu);
          break;
        case AliPID::kPion:
          color = getLineColor(kPi);
          break;
        case AliPID::kProton:
          color = getLineColor(kPr);
          break;
      }
      
      hMCgenYieldsPrimSpecies[i]->SetLineColor(color);
      hMCgenYieldsPrimSpecies[i]->SetMarkerColor(color);
      hMCgenYieldsPrimSpecies[i]->SetMarkerStyle(28);
      hMCgenYieldsPrimSpecies[i]->SetLineStyle(1);
      hMCgenYieldsPrimSpecies[i]->GetXaxis()->SetTitleOffset(1.0);
      hMCgenYieldsPrimSpecies[i]->SetStats(kFALSE);
      
      SetReasonableAxisRange(hMCgenYieldsPrimSpecies[i]->GetXaxis(), mode, pLow, pHigh);
      normaliseGenYieldMCtruthHist(hMCgenYieldsPrimSpecies[i], numEvents, deta, restrictJetPtAxis, numJetsGen);
    }
  }
  
  
  // Compare data points with MC (yield)
  TCanvas* cYieldComparisons = drawYieldComparisonToMC(mode, pLow, pHigh, hYieldPions, hYieldKaons, hYieldProtons,
                                                       hYieldElectrons, hYieldMuons, hYieldPionsMC, hYieldKaonsMC,
                                                       hYieldProtonsMC, hYieldElectronsMC, hYieldMuonsMC, 
                                                       hYieldComparisonPions, hYieldComparisonKaons, 
                                                       hYieldComparisonProtons, hYieldComparisonElectrons, 
                                                       hYieldComparisonMuons);
  
  
  TCanvas* cFractionsPions = drawFractionHistos("cFractionsPions", "Pion fractions", mode, pLow, pHigh, hFractionPionsDeltaPion, 
                                                hFractionPionsDeltaElectron, hFractionPionsDeltaKaon, hFractionPionsDeltaProton,
                                                hFractionPionsMC, plotIdentifiedSpectra);
  
  
  TCanvas* cFractionsElectrons = drawFractionHistos("cFractionsElectrons", "Electron fractions", mode, pLow, pHigh, 
                                                    hFractionElectronsDeltaPion, hFractionElectronsDeltaElectron,
                                                    hFractionElectronsDeltaKaon, hFractionElectronsDeltaProton, hFractionElectronsMC, 
                                                    plotIdentifiedSpectra);
  
  TCanvas* cFractionsKaons = drawFractionHistos("cFractionsKaons", "Kaon fractions", mode, pLow, pHigh, hFractionKaonsDeltaPion, 
                                                hFractionKaonsDeltaElectron, hFractionKaonsDeltaKaon, hFractionKaonsDeltaProton,
                                                hFractionKaonsMC, plotIdentifiedSpectra);
                                                
  TCanvas* cFractionsProtons = drawFractionHistos("cFractionsProtons", "Proton fractions", mode, pLow, pHigh, hFractionProtonsDeltaPion, 
                                                hFractionProtonsDeltaElectron, hFractionProtonsDeltaKaon, hFractionProtonsDeltaProton,
                                                hFractionProtonsMC, plotIdentifiedSpectra);
  
  TCanvas* cFractionsMuons = drawFractionHistos("cFractionsMuons", "Muon fractions", mode, pLow, pHigh, hFractionMuonsDeltaPion, 
                                                hFractionMuonsDeltaElectron, hFractionMuonsDeltaKaon, hFractionMuonsDeltaProton,
                                                hFractionMuonsMC, plotIdentifiedSpectra);
                                              
  
  
  TCanvas* cYields = drawFinalYields(mode, pLow, pHigh, hYieldPions, hYieldKaons, hYieldProtons, hYieldElectrons, hYieldMuons,
                                     hYieldPionsMC, hYieldKaonsMC, hYieldProtonsMC, hYieldElectronsMC, hYieldMuonsMC,
                                     plotIdentifiedSpectra);
  
  
  TCanvas* cYieldsPions = drawYieldHistos("cYieldsPions", "Pion yields", mode, pLow, pHigh, hYieldPionsDeltaPion, hYieldPionsDeltaElectron,
                                          hYieldPionsDeltaKaon, hYieldPionsDeltaProton, hYieldPionsMC, plotIdentifiedSpectra);
  
  
  TCanvas* cYieldsElectrons = drawYieldHistos("cYieldsElectrons", "Electron yields", mode, pLow, pHigh, hYieldElectronsDeltaPion,
                                              hYieldElectronsDeltaElectron, hYieldElectronsDeltaKaon, hYieldElectronsDeltaProton, hYieldElectronsMC,
                                              plotIdentifiedSpectra);
  
  TCanvas* cYieldsKaons = drawYieldHistos("cYieldsKaons", "Kaon yields", mode, pLow, pHigh, hYieldKaonsDeltaPion, hYieldKaonsDeltaElectron,
                                          hYieldKaonsDeltaKaon, hYieldKaonsDeltaProton, hYieldKaonsMC, plotIdentifiedSpectra);
  
  TCanvas* cYieldsProtons = drawYieldHistos("cYieldsProtons", "Proton yields", mode, pLow, pHigh, hYieldProtonsDeltaPion, hYieldProtonsDeltaElectron,
                                            hYieldProtonsDeltaKaon, hYieldProtonsDeltaProton, hYieldProtonsMC, plotIdentifiedSpectra);
  
  TCanvas* cYieldsMuons = drawYieldHistos("cYieldsMuons", "Muon yields", mode, pLow, pHigh, hYieldMuonsDeltaPion, hYieldMuonsDeltaElectron,
                                            hYieldMuonsDeltaKaon, hYieldMuonsDeltaProton, hYieldMuonsMC, plotIdentifiedSpectra);
  
  
  // Extract TOF efficiency -> NOTE: Must be done as long as the yields are still TPC only!
  TH1F* hTOFEfficiencyKaons = 0x0;
  TH1F* hTOFEfficiencyPions = 0x0;
  TH1F* hTOFEfficiencyProtons = 0x0;
  
  // Same for MC
  TH1F* hTOFEfficiencyKaonsMC = 0x0;
  TH1F* hTOFEfficiencyPionsMC = 0x0;
  TH1F* hTOFEfficiencyProtonsMC = 0x0;
    
  if (applyTOFpatching) {
    // Use the TOF pions AFTER cleanup for the efficiency (is anyway an effect < 1% for the pions)
    ExtractTOFEfficiency(hYieldTOFKaons, hYieldTOFPions, hYieldTOFProtons, hYieldKaons, hYieldPions, hYieldProtons, 
                         &hTOFEfficiencyKaons, &hTOFEfficiencyPions, &hTOFEfficiencyProtons, kFALSE);
    if (isMC) {
      ExtractTOFEfficiency(hYieldTOFKaonsMC, hYieldTOFPionsMC, hYieldTOFProtonsMC, hYieldKaonsMC, hYieldPionsMC, hYieldProtonsMC, 
                           &hTOFEfficiencyKaonsMC, &hTOFEfficiencyPionsMC, &hTOFEfficiencyProtonsMC, kTRUE);
      
    }
  }
  
  TH1F* hYieldTOFTotalMC = 0x0;
  TH1F* hFractionTOFPions = 0x0;
  TH1F* hFractionTOFKaons = 0x0;
  TH1F* hFractionTOFProtons = 0x0;
  TH1F* hFractionTOFPionsMC = 0x0;
  TH1F* hFractionTOFKaonsMC = 0x0;
  TH1F* hFractionTOFProtonsMC = 0x0;
  TH1F* hFractionTOFElectronsMC = 0x0;
  TH1F* hFractionTOFMuonsMC = 0x0;
  TH1F* hFractionTOFPionsPlusMuonsMC = 0x0;
  TH1F* hFractionTOFComparisonPions = 0x0;
  TH1F* hFractionTOFComparisonKaons = 0x0;
  TH1F* hFractionTOFComparisonProtons = 0x0;
  TH1F* hFractionTOFComparisonPionsPlusMuons = 0x0;
  
  TCanvas* cFractionTOFComparisons = 0x0;
  
  // Also for the comparison: Take the cleaned up TOF pions "without" electrons. Electrons are not interesting here.
  if (applyTOFpatching)
    cFractionTOFComparisons = drawTOFFractionComparisonToMC(mode, pLow, pHigh, hFractionPions->GetYaxis()->GetTitle(),
                                                            hYieldTOFTotal, &hYieldTOFTotalMC,
                                                            hYieldTOFPions, hYieldTOFKaons, hYieldTOFProtons,
                                                            hYieldTOFPionsMC, hYieldTOFKaonsMC, hYieldTOFProtonsMC,
                                                            hYieldTOFElectronsMC, hYieldTOFMuonsMC,
                                                            &hFractionTOFPions, &hFractionTOFKaons, &hFractionTOFProtons,
                                                            &hFractionTOFPionsMC, &hFractionTOFKaonsMC, &hFractionTOFProtonsMC,
                                                            &hFractionTOFElectronsMC, &hFractionTOFMuonsMC, 
                                                            &hFractionTOFPionsPlusMuonsMC,
                                                            &hFractionTOFComparisonPions, &hFractionTOFComparisonKaons,
                                                            &hFractionTOFComparisonProtons, &hFractionTOFComparisonPionsPlusMuons);
  
  // Save final results
  saveF->cd();
  
  if (fElectronFraction)
    fElectronFraction->Write();
  
  if (hFractionElectrons)
    hFractionElectrons->Write(0, TObject::kWriteDelete);
    
  if (hFractionKaons)
    hFractionKaons->Write(0, TObject::kWriteDelete);
  
  if (hFractionPions)
    hFractionPions->Write(0, TObject::kWriteDelete);
  
  if (hFractionProtons)
    hFractionProtons->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuons)
    hFractionMuons->Write(0, TObject::kWriteDelete);
  
  if (hFractionSummed)
    hFractionSummed->Write(0, TObject::kWriteDelete);
    
  
  if (hFractionElectronsDeltaElectron)
    hFractionElectronsDeltaElectron->Write(0, TObject::kWriteDelete);
      
  if (hFractionKaonsDeltaElectron)
    hFractionKaonsDeltaElectron->Write(0, TObject::kWriteDelete);
    
  if (hFractionPionsDeltaElectron)
    hFractionPionsDeltaElectron->Write(0, TObject::kWriteDelete);
    
  if (hFractionProtonsDeltaElectron)
    hFractionProtonsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsDeltaElectron)
    hFractionMuonsDeltaElectron->Write(0, TObject::kWriteDelete);
      
      
  if (hFractionElectronsDeltaPion)
    hFractionElectronsDeltaPion->Write(0, TObject::kWriteDelete);
      
  if (hFractionKaonsDeltaPion)
    hFractionKaonsDeltaPion->Write(0, TObject::kWriteDelete);
    
  if (hFractionPionsDeltaPion)
    hFractionPionsDeltaPion->Write(0, TObject::kWriteDelete);
    
  if (hFractionProtonsDeltaPion)
    hFractionProtonsDeltaPion->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsDeltaPion)
    hFractionMuonsDeltaPion->Write(0, TObject::kWriteDelete);
      
    
  if (hFractionElectronsDeltaKaon)
    hFractionElectronsDeltaKaon->Write(0, TObject::kWriteDelete);
      
  if (hFractionKaonsDeltaKaon)
    hFractionKaonsDeltaKaon->Write(0, TObject::kWriteDelete);
    
  if (hFractionPionsDeltaKaon)
    hFractionPionsDeltaKaon->Write(0, TObject::kWriteDelete);
    
  if (hFractionProtonsDeltaKaon)
    hFractionProtonsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsDeltaKaon)
    hFractionMuonsDeltaKaon->Write(0, TObject::kWriteDelete);
      
    
  if (hFractionElectronsDeltaProton)
    hFractionElectronsDeltaProton->Write(0, TObject::kWriteDelete);
      
  if (hFractionKaonsDeltaProton)
    hFractionKaonsDeltaProton->Write(0, TObject::kWriteDelete);
    
  if (hFractionPionsDeltaProton)
    hFractionPionsDeltaProton->Write(0, TObject::kWriteDelete);
    
  if (hFractionProtonsDeltaProton)
    hFractionProtonsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsDeltaProton)
    hFractionMuonsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hNumEvents)
    hNumEvents->Write();
  
  if (hNumEventsTriggerSel)
    hNumEventsTriggerSel->Write();
  
  if (hNumEventsTriggerSelVtxCut)
    hNumEventsTriggerSelVtxCut->Write();
  
  if (hNumEventsTriggerSelVtxCutNoPileUp)
    hNumEventsTriggerSelVtxCutNoPileUp->Write();
  
  if (hNjetsGen)
    hNjetsGen->Write();
  
  if (hNjetsRec)
    hNjetsRec->Write();
  
  if (cFractions)
    cFractions->Write();
  if (cFractionsPions)
    cFractionsPions->Write();
  if (cFractionsElectrons)
    cFractionsElectrons->Write();
  if (cFractionsKaons)
    cFractionsKaons->Write();
  if (cFractionsProtons)
    cFractionsProtons->Write();
  if (cFractionsMuons)
    cFractionsMuons->Write();
  
  
  if (hFractionElectronsMC)
    hFractionElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionKaonsMC)
    hFractionKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionPionsMC)
    hFractionPionsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionMuonsMC)
    hFractionMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hFractionProtonsMC)
    hFractionProtonsMC->Write(0, TObject::kWriteDelete);
  
  
  if (hFractionComparisonElectrons)
    hFractionComparisonElectrons->Write(0, TObject::kWriteDelete);
  
  if (hFractionComparisonMuons)
    hFractionComparisonMuons->Write(0, TObject::kWriteDelete);
  
  if (hFractionComparisonKaons)
    hFractionComparisonKaons->Write(0, TObject::kWriteDelete);
  
  if (hFractionComparisonPions)
    hFractionComparisonPions->Write(0, TObject::kWriteDelete);
  
  if (hFractionComparisonProtons)
    hFractionComparisonProtons->Write(0, TObject::kWriteDelete);
  
  if (hFractionComparisonTotal)
    hFractionComparisonTotal->Write(0, TObject::kWriteDelete);
  
  if (cFractionComparisons)
    cFractionComparisons->Write();
  
  
  if (hYieldComparisonElectrons)
    hYieldComparisonElectrons->Write(0, TObject::kWriteDelete);
  
  if (hYieldComparisonMuons)
    hYieldComparisonMuons->Write(0, TObject::kWriteDelete);
  
  if (hYieldComparisonKaons)
    hYieldComparisonKaons->Write(0, TObject::kWriteDelete);
  
  if (hYieldComparisonPions)
    hYieldComparisonPions->Write(0, TObject::kWriteDelete);
  
  if (hYieldComparisonProtons)
    hYieldComparisonProtons->Write(0, TObject::kWriteDelete);
  
  if (cYieldComparisons)
    cYieldComparisons->Write();
  
  
  if (hYieldElectrons)
    hYieldElectrons->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaons)
    hYieldKaons->Write(0, TObject::kWriteDelete);
  
  if (hYieldPions)
    hYieldPions->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtons)
    hYieldProtons->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuons)
    hYieldMuons->Write(0, TObject::kWriteDelete);
  
  
  if (hYieldElectronsDeltaElectron)
    hYieldElectronsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsDeltaElectron)
    hYieldKaonsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsDeltaElectron)
    hYieldPionsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsDeltaElectron)
    hYieldProtonsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsDeltaElectron)
    hYieldMuonsDeltaElectron->Write(0, TObject::kWriteDelete);
  
  
  if (hYieldElectronsDeltaPion)
    hYieldElectronsDeltaPion->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsDeltaPion)
    hYieldKaonsDeltaPion->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsDeltaPion)
    hYieldPionsDeltaPion->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsDeltaPion)
    hYieldProtonsDeltaPion->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsDeltaPion)
    hYieldMuonsDeltaPion->Write(0, TObject::kWriteDelete);
  
  
  if (hYieldElectronsDeltaKaon)
    hYieldElectronsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsDeltaKaon)
    hYieldKaonsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsDeltaKaon)
    hYieldPionsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsDeltaKaon)
    hYieldProtonsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsDeltaKaon)
    hYieldMuonsDeltaKaon->Write(0, TObject::kWriteDelete);
  
  
  if (hYieldElectronsDeltaProton)
    hYieldElectronsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsDeltaProton)
    hYieldKaonsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsDeltaProton)
    hYieldPionsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsDeltaProton)
    hYieldProtonsDeltaProton->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsDeltaProton)
    hYieldMuonsDeltaProton->Write(0, TObject::kWriteDelete);
  
  
  if (hYieldElectronsMC)
    hYieldElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldKaonsMC)
    hYieldKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldPionsMC)
    hYieldPionsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldMuonsMC)
    hYieldMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldProtonsMC)
    hYieldProtonsMC->Write(0, TObject::kWriteDelete);
  
  if (hYieldSummedMC)
    hYieldSummedMC->Write(0, TObject::kWriteDelete);
  
  
  if (hRatioToPiElectrons)
    hRatioToPiElectrons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiMuons)
    hRatioToPiMuons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiKaons)
    hRatioToPiKaons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiProtons)
    hRatioToPiProtons->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiElectronsMC)
    hRatioToPiElectronsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiMuonsMC)
    hRatioToPiMuonsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiKaonsMC)
    hRatioToPiKaonsMC->Write(0, TObject::kWriteDelete);
  
  if (hRatioToPiProtonsMC)
    hRatioToPiProtonsMC->Write(0, TObject::kWriteDelete);
  
  
  
  if (hReducedChiSquarePt)
    hReducedChiSquarePt->Write(0, TObject::kWriteDelete);
  
  if (cYields)
    cYields->Write();
  if (cYieldsPions)
    cYieldsPions->Write();
  if (cYieldsElectrons)
    cYieldsElectrons->Write();
  if (cYieldsKaons)
    cYieldsKaons->Write();
  if (cYieldsProtons)
    cYieldsProtons->Write();
  if (cYieldsMuons)
    cYieldsMuons->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenYieldsPrimSpecies[i])
      hMCgenYieldsPrimSpecies[i]->Write();
  }
  
  if (hTOFPurity)
    hTOFPurity->Write();
  
  if (hTOFEfficiencyKaons)
    hTOFEfficiencyKaons->Write();
  
  if (hTOFEfficiencyPions)
    hTOFEfficiencyPions->Write();
  
  if (hTOFEfficiencyProtons)
    hTOFEfficiencyProtons->Write();
  
  if (hTOFEfficiencyKaonsMC)
    hTOFEfficiencyKaonsMC->Write();
  
  if (hTOFEfficiencyPionsMC)
    hTOFEfficiencyPionsMC->Write();
  
  if (hTOFEfficiencyProtonsMC)
    hTOFEfficiencyProtonsMC->Write();
  
  if (hYieldTOFElectrons)
    hYieldTOFElectrons->Write();
  
  if (hYieldTOFPions)
    hYieldTOFPions->Write();
  
  if (hYieldTOFKaons)
    hYieldTOFKaons->Write();
  
  if (hYieldTOFProtons)
    hYieldTOFProtons->Write();
  
  if (hYieldTOFTotal)
    hYieldTOFTotal->Write();
  
  if (hYieldTPConlyTotal)
    hYieldTPConlyTotal->Write();
  
  if (hYieldTOFElectronsMC)
    hYieldTOFElectronsMC->Write();
  
  if (hYieldTOFMuonsMC)
    hYieldTOFMuonsMC->Write();
  
  if (hYieldTOFKaonsMC)
    hYieldTOFKaonsMC->Write();
  
  if (hYieldTOFPionsMC)
    hYieldTOFPionsMC->Write();
  
  if (hYieldTOFProtonsMC)
    hYieldTOFProtonsMC->Write();
  
  if (hYieldTOFTotalMC)
    hYieldTOFTotalMC->Write();
  
  if (hFractionTOFPions)
    hFractionTOFPions->Write();
  
  if (hFractionTOFKaons)
    hFractionTOFKaons->Write();
  
  if (hFractionTOFProtons)
    hFractionTOFProtons->Write();  
  
  if (hFractionTOFPionsMC)
    hFractionTOFPionsMC->Write();
  
  if (hFractionTOFElectronsMC)
    hFractionTOFElectronsMC->Write();
  
  if (hFractionTOFMuonsMC)
    hFractionTOFMuonsMC->Write();
  
  if (hFractionTOFKaonsMC)
    hFractionTOFKaonsMC->Write();
  
  if (hFractionTOFProtonsMC)
    hFractionTOFProtonsMC->Write();
  
  if (hFractionTOFPionsPlusMuonsMC)
    hFractionTOFPionsPlusMuonsMC->Write();
  
  if (hFractionTOFComparisonPions)
    hFractionTOFComparisonPions->Write();
  
  if (hFractionTOFComparisonKaons)
    hFractionTOFComparisonKaons->Write();
  
  if (hFractionTOFComparisonProtons)
    hFractionTOFComparisonProtons->Write();
  
  if (hFractionTOFComparisonPionsPlusMuons)
    hFractionTOFComparisonPionsPlusMuons->Write();
  
  if (cFractionTOFComparisons)
    cFractionTOFComparisons->Write();
  
  
  if (filePathNameResults)
    *filePathNameResults = saveFName;
  
  if (TMath::Abs(mathFit->GetScaleFactorError() - 1.) > 1e-6) {
    // If the deltaPrime range is large enough, we artificially get a factor 4 in statistics by looking at the four
    // different deltaPrimeSpecies, which have (except for binning effects) the same information. 
    // Therefore, to get the "real" statistical error, we need to multiply the obtained error by sqrt(4) = 2
    std::cout << "NOTE: Errors multiplied by " << mathFit->GetScaleFactorError() 
              << " to take into account artificially higher statistics (by factor of 4) due to same information "
              << "for all deltaPrimeSpecies (except for binning effects), if deltaPrimeRange sufficiently large!" << std::endl
              << std::endl;
  }
  
  if (fitMethod < 2) {
    std::cout << "WARNING: Errors might be wrong! Especially, for the to-pi ratios there are no correlations taken into account!"
              << std::endl;
  }
  
  if (isMCdataSet) {
    std::cout << "Special templates have been deactivated for MC!" << std::endl;
  }
  
  delete cFractions;
  delete cFractionComparisons;
  delete cYieldComparisons;
  delete cFractionsPions;
  delete cFractionsElectrons;
  delete cFractionsKaons;
  delete cFractionsProtons;
  delete cFractionsMuons;
  delete cYields;
  delete cYieldsPions;
  delete cYieldsKaons;
  delete cYieldsMuons;
  delete cYieldsProtons;
  delete cYieldsElectrons;
  
  if (applyTOFpatching) {
    TString saveFTOFName = saveFName;
    saveFTOFName = saveFTOFName.ReplaceAll("_TPConly.root", "_TOFpatched.root");
    TFile* saveFTOF = TFile::Open(saveFTOFName.Data(), "RECREATE");
    saveFTOF->cd();
    
    if (hTOFPurity)
      hTOFPurity->Write();
    
    if (hTOFEfficiencyKaons)
      hTOFEfficiencyKaons->Write();
    
    if (hTOFEfficiencyPions)
      hTOFEfficiencyPions->Write();
    
    if (hTOFEfficiencyProtons)
      hTOFEfficiencyProtons->Write();
    
    if (hTOFEfficiencyKaonsMC)
      hTOFEfficiencyKaonsMC->Write();
    
    if (hTOFEfficiencyPionsMC)
      hTOFEfficiencyPionsMC->Write();
    
    if (hTOFEfficiencyProtonsMC)
      hTOFEfficiencyProtonsMC->Write();
    
    if (hYieldTOFElectrons)
      hYieldTOFElectrons->Write();
    
    if (hYieldTOFPions)
      hYieldTOFPions->Write();
    
    if (hYieldTOFKaons)
      hYieldTOFKaons->Write();
    
    if (hYieldTOFProtons)
      hYieldTOFProtons->Write();
    
    if (hYieldTOFTotal)
      hYieldTOFTotal->Write();
    
    if (hYieldTPConlyTotal)
      hYieldTPConlyTotal->Write();
    
    if (hYieldTOFElectronsMC)
      hYieldTOFElectronsMC->Write();
    
    if (hYieldTOFMuonsMC)
      hYieldTOFMuonsMC->Write();
    
    if (hYieldTOFKaonsMC)
      hYieldTOFKaonsMC->Write();
    
    if (hYieldTOFPionsMC)
      hYieldTOFPionsMC->Write();
    
    if (hYieldTOFProtonsMC)
      hYieldTOFProtonsMC->Write();
    
    if (hYieldTOFTotalMC)
      hYieldTOFTotalMC->Write();
    
    if (hFractionTOFPions)
      hFractionTOFPions->Write();
    
    if (hFractionTOFKaons)
      hFractionTOFKaons->Write();
    
    if (hFractionTOFProtons)
      hFractionTOFProtons->Write();  
    
    if (hFractionTOFPionsMC)
      hFractionTOFPionsMC->Write();
    
    if (hFractionTOFElectronsMC)
      hFractionTOFElectronsMC->Write();
    
    if (hFractionTOFMuonsMC)
      hFractionTOFMuonsMC->Write();
    
    if (hFractionTOFKaonsMC)
      hFractionTOFKaonsMC->Write();
    
    if (hFractionTOFProtonsMC)
      hFractionTOFProtonsMC->Write();
    
    if (hFractionTOFPionsPlusMuonsMC)
      hFractionTOFPionsPlusMuonsMC->Write();
    
    if (hFractionTOFComparisonPions)
      hFractionTOFComparisonPions->Write();
    
    if (hFractionTOFComparisonKaons)
      hFractionTOFComparisonKaons->Write();
    
    if (hFractionTOFComparisonProtons)
      hFractionTOFComparisonProtons->Write();
    
    if (hFractionTOFComparisonPionsPlusMuons)
      hFractionTOFComparisonPionsPlusMuons->Write();
    
    if (cFractionTOFComparisons)
      cFractionTOFComparisons->Write();
    
    // Note: The normalisation of the yields is arbitrary for the patching of fractions and to-pion-ratios.
    // But to be consistent with the yield (and the patch for this), the same normalisation is used for the TOF and
    // TPC-only yields (inverseBinWidth, dEta, numEvents/numJets, etc).
    // Note also: Binning always the same
      
    Double_t fractionTPC = 0, fractionErrorTPC = 0, toPiRatioTPC = 0, toPiRatioErrorTPC = 0, yieldTPConlyPi = 0, 
             yieldTPConlyTotal = 0, yieldTOFspecies = 0, yieldTOFpi = 0, yieldTOFtotal = 0;
    
    for (Int_t binX = 1; binX <= hFractionPions->GetNbinsX(); binX++) { 
      // For the yields, just patch the fractions and then scale with the total yield (TPC only + TOF).
      // Since the TOF yields are assumed to be precise, this is the same as just adding up TPC only yield and
      // TOF yield and only taking into account the error of the TPC only yield.
      
      yieldTPConlyPi = hYieldPions->GetBinContent(binX);
      yieldTPConlyTotal = hYieldTPConlyTotal->GetBinContent(binX);
      yieldTOFtotal = hYieldTOFTotal->GetBinContent(binX);
      yieldTOFpi = hYieldTOFPions->GetBinContent(binX);
      
      // Electrons
      fractionTPC = hFractionElectrons->GetBinContent(binX);
      fractionErrorTPC = hFractionElectrons->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiElectrons->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiElectrons->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFElectrons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionElectrons->SetBinContent(binX, fractionTPC);
      hFractionElectrons->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiElectrons->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiElectrons->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldElectrons->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldElectrons->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Muons
      fractionTPC = hFractionMuons->GetBinContent(binX);
      fractionErrorTPC = hFractionMuons->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiMuons->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiMuons->GetBinError(binX);
      
      yieldTOFspecies = 0.;
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionMuons->SetBinContent(binX, fractionTPC);
      hFractionMuons->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiMuons->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiMuons->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldMuons->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldMuons->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Kaons
      fractionTPC = hFractionKaons->GetBinContent(binX);
      fractionErrorTPC = hFractionKaons->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiKaons->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiKaons->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFKaons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionKaons->SetBinContent(binX, fractionTPC);
      hFractionKaons->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiKaons->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiKaons->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldKaons->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldKaons->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Protons
      fractionTPC = hFractionProtons->GetBinContent(binX);
      fractionErrorTPC = hFractionProtons->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiProtons->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiProtons->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFProtons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionProtons->SetBinContent(binX, fractionTPC);
      hFractionProtons->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiProtons->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiProtons->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldProtons->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldProtons->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Pions
      fractionTPC = hFractionPions->GetBinContent(binX);
      fractionErrorTPC = hFractionPions->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFPions->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      
      hFractionPions->SetBinContent(binX, fractionTPC);
      hFractionPions->SetBinError(binX, fractionErrorTPC);
      
      hYieldPions->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldPions->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
    }
    
    // Do the same patching, but take the MC truth for TPC only instead of the fit result
    const Int_t markerStyleForMCTOFpatched = 29;
    TH1F* hYieldPionsMCTOFpatched = new TH1F(*hYieldPionsMC);
    hYieldPionsMCTOFpatched->SetName("hYieldPionsMCTOFpatched");
    hYieldPionsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hYieldKaonsMCTOFpatched = new TH1F(*hYieldKaonsMC);
    hYieldKaonsMCTOFpatched->SetName("hYieldKaonsMCTOFpatched");
    hYieldKaonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hYieldProtonsMCTOFpatched = new TH1F(*hYieldProtonsMC);
    hYieldProtonsMCTOFpatched->SetName("hYieldProtonsMCTOFpatched");
    hYieldProtonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hYieldMuonsMCTOFpatched = new TH1F(*hYieldMuonsMC);
    hYieldMuonsMCTOFpatched->SetName("hYieldMuonsMCTOFpatched");
    hYieldMuonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hYieldElectronsMCTOFpatched = new TH1F(*hYieldElectronsMC);
    hYieldElectronsMCTOFpatched->SetName("hYieldElectronsMCTOFpatched");
    hYieldElectronsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    
    TH1F* hFractionPionsMCTOFpatched = new TH1F(*hFractionPionsMC);
    hFractionPionsMCTOFpatched->SetName("hFractionPionsMCTOFpatched");
    hFractionPionsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hFractionKaonsMCTOFpatched = new TH1F(*hFractionKaonsMC);
    hFractionKaonsMCTOFpatched->SetName("hFractionKaonsMCTOFpatched");
    hFractionKaonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hFractionProtonsMCTOFpatched = new TH1F(*hFractionProtonsMC);
    hFractionProtonsMCTOFpatched->SetName("hFractionProtonsMCTOFpatched");
    hFractionProtonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hFractionMuonsMCTOFpatched = new TH1F(*hFractionMuonsMC);
    hFractionMuonsMCTOFpatched->SetName("hFractionMuonsMCTOFpatched");
    hFractionMuonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hFractionElectronsMCTOFpatched = new TH1F(*hFractionElectronsMC);
    hFractionElectronsMCTOFpatched->SetName("hFractionElectronsMCTOFpatched");
    hFractionElectronsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    
    TH1F* hRatioToPiKaonsMCTOFpatched = new TH1F(*hRatioToPiKaonsMC);
    hRatioToPiKaonsMCTOFpatched->SetName("hRatioToPiKaonsMCTOFpatched");
    hRatioToPiKaonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hRatioToPiProtonsMCTOFpatched = new TH1F(*hRatioToPiProtonsMC);
    hRatioToPiProtonsMCTOFpatched->SetName("hRatioToPiProtonsMCTOFpatched");
    hRatioToPiProtonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hRatioToPiMuonsMCTOFpatched = new TH1F(*hRatioToPiMuonsMC);
    hRatioToPiMuonsMCTOFpatched->SetName("hRatioToPiMuonsMCTOFpatched");
    hRatioToPiMuonsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    TH1F* hRatioToPiElectronsMCTOFpatched = new TH1F(*hRatioToPiElectronsMC);
    hRatioToPiElectronsMCTOFpatched->SetName("hRatioToPiElectronsMCTOFpatched");
    hRatioToPiElectronsMCTOFpatched->SetMarkerStyle(markerStyleForMCTOFpatched);
    
    
    for (Int_t binX = 1; binX <= hFractionPions->GetNbinsX(); binX++) { 
      yieldTPConlyPi = hYieldPionsMC->GetBinContent(binX);
      yieldTPConlyTotal = hYieldTPConlyTotal->GetBinContent(binX);
      yieldTOFtotal = hYieldTOFTotal->GetBinContent(binX);
      yieldTOFpi = hYieldTOFPions->GetBinContent(binX);
      
      // Electrons
      fractionTPC = hFractionElectronsMC->GetBinContent(binX);
      fractionErrorTPC = hFractionElectronsMC->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiElectronsMC->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiElectronsMC->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFElectrons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionElectronsMCTOFpatched->SetBinContent(binX, fractionTPC);
      hFractionElectronsMCTOFpatched->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiElectronsMCTOFpatched->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiElectronsMCTOFpatched->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldElectronsMCTOFpatched->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldElectronsMCTOFpatched->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Muons
      fractionTPC = hFractionMuonsMC->GetBinContent(binX);
      fractionErrorTPC = hFractionMuonsMC->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiMuonsMC->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiMuonsMC->GetBinError(binX);
      
      yieldTOFspecies = 0.;
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionMuonsMCTOFpatched->SetBinContent(binX, fractionTPC);
      hFractionMuonsMCTOFpatched->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiMuonsMCTOFpatched->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiMuonsMCTOFpatched->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldMuonsMCTOFpatched->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldMuonsMCTOFpatched->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Kaons
      fractionTPC = hFractionKaonsMC->GetBinContent(binX);
      fractionErrorTPC = hFractionKaonsMC->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiKaonsMC->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiKaonsMC->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFKaons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionKaonsMCTOFpatched->SetBinContent(binX, fractionTPC);
      hFractionKaonsMCTOFpatched->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiKaonsMCTOFpatched->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiKaonsMCTOFpatched->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldKaonsMCTOFpatched->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldKaonsMCTOFpatched->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Protons
      fractionTPC = hFractionProtonsMC->GetBinContent(binX);
      fractionErrorTPC = hFractionProtonsMC->GetBinError(binX);
      
      toPiRatioTPC = hRatioToPiProtonsMC->GetBinContent(binX);
      toPiRatioErrorTPC = hRatioToPiProtonsMC->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFProtons->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      PatchRatioWithTOF(toPiRatioTPC, toPiRatioErrorTPC, yieldTPConlyPi, yieldTOFspecies, yieldTOFpi);
      
      hFractionProtonsMCTOFpatched->SetBinContent(binX, fractionTPC);
      hFractionProtonsMCTOFpatched->SetBinError(binX, fractionErrorTPC);
      
      hRatioToPiProtonsMCTOFpatched->SetBinContent(binX, toPiRatioTPC);
      hRatioToPiProtonsMCTOFpatched->SetBinError(binX, toPiRatioErrorTPC);
      
      hYieldProtonsMCTOFpatched->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldProtonsMCTOFpatched->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
      
      // Pions
      fractionTPC = hFractionPionsMC->GetBinContent(binX);
      fractionErrorTPC = hFractionPionsMC->GetBinError(binX);
      
      yieldTOFspecies = hYieldTOFPions->GetBinContent(binX);
      
      PatchFractionWithTOF(fractionTPC, fractionErrorTPC, yieldTPConlyTotal, yieldTOFspecies, yieldTOFtotal);
      
      hFractionPionsMCTOFpatched->SetBinContent(binX, fractionTPC);
      hFractionPionsMCTOFpatched->SetBinError(binX, fractionErrorTPC);
      
      hYieldPionsMCTOFpatched->SetBinContent(binX, fractionTPC * (yieldTPConlyTotal + yieldTOFtotal));
      hYieldPionsMCTOFpatched->SetBinError(binX, fractionErrorTPC * (yieldTPConlyTotal + yieldTOFtotal));
    }
    
    // Check again total fraction
    hFractionSummed->Reset();
    hFractionSummed->Add(hFractionElectrons);
    hFractionSummed->Add(hFractionMuons);
    hFractionSummed->Add(hFractionKaons);
    hFractionSummed->Add(hFractionPions);
    hFractionSummed->Add(hFractionProtons);
    
    TH1F* hFractionSummedMCTOFpatched = new TH1F(*hFractionSummed);
    hFractionSummedMCTOFpatched->SetName("hFractionSummedMCTOFpatched");
    hFractionSummedMCTOFpatched->Reset();
    hFractionSummedMCTOFpatched->Add(hFractionElectronsMCTOFpatched);
    hFractionSummedMCTOFpatched->Add(hFractionMuonsMCTOFpatched);
    hFractionSummedMCTOFpatched->Add(hFractionPionsMCTOFpatched);
    hFractionSummedMCTOFpatched->Add(hFractionKaonsMCTOFpatched);
    hFractionSummedMCTOFpatched->Add(hFractionProtonsMCTOFpatched);
    
    // Sum up MC fractions for TPC only and TOF (independent!) and calculate new fractions and to-pion-ratios
    hYieldElectronsMC->Add(hYieldTOFElectronsMC);
    hYieldMuonsMC->Add(hYieldTOFMuonsMC);
    hYieldKaonsMC->Add(hYieldTOFKaonsMC);
    hYieldPionsMC->Add(hYieldTOFPionsMC);
    hYieldProtonsMC->Add(hYieldTOFProtonsMC);
    
    hYieldSummedMC->Reset();
    hYieldSummedMC->Add(hYieldElectronsMC);
    hYieldSummedMC->Add(hYieldMuonsMC);
    hYieldSummedMC->Add(hYieldKaonsMC);
    hYieldSummedMC->Add(hYieldPionsMC);
    hYieldSummedMC->Add(hYieldProtonsMC);
    
    // MCspecies and MCtotal are correlated. This can be taken into account via using the binomial error in the division
    hFractionElectronsMC->Divide(hYieldElectronsMC, hYieldSummedMC, 1., 1., "B");
    hFractionMuonsMC->Divide(hYieldMuonsMC, hYieldSummedMC, 1., 1., "B");
    hFractionKaonsMC->Divide(hYieldKaonsMC, hYieldSummedMC, 1., 1., "B");
    hFractionPionsMC->Divide(hYieldPionsMC, hYieldSummedMC, 1., 1., "B");
    hFractionProtonsMC->Divide(hYieldProtonsMC, hYieldSummedMC, 1., 1., "B");
    
    // Calculate MC to-pi ratios -> In MC the yields are uncorrelated, so just divide the histos to get the correct result
    hRatioToPiElectronsMC->Divide(hYieldElectronsMC, hYieldPionsMC);
    hRatioToPiMuonsMC->Divide(hYieldMuonsMC, hYieldPionsMC);
    hRatioToPiKaonsMC->Divide(hYieldKaonsMC, hYieldPionsMC);
    hRatioToPiProtonsMC->Divide(hYieldProtonsMC, hYieldPionsMC);
    
    TCanvas* cFractionMCTOFpatched = drawFinalFractions(mode, pLow, pHigh, electronFixingUsed,
                                                        hFractionPionsMCTOFpatched, hFractionKaonsMCTOFpatched,  
                                                        hFractionProtonsMCTOFpatched, hFractionElectronsMCTOFpatched,
                                                        hFractionMuonsMCTOFpatched, hFractionSummedMCTOFpatched, hFractionPionsMC, 
                                                        hFractionKaonsMC, hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, 
                                                        plotIdentifiedSpectra, kTRUE);

    TString compYtitle = "Fraction (MC w/ TOF) / fraction (MC direct)";
    TH1F* hFractionMCTOFpatchedComparisonPions = new TH1F(*hFractionComparisonPions);
    hFractionMCTOFpatchedComparisonPions->SetName("hFractionMCTOFpatchedComparisonPions");
    hFractionMCTOFpatchedComparisonPions->GetYaxis()->SetTitle(compYtitle.Data());
    hFractionMCTOFpatchedComparisonPions->GetYaxis()->SetTitleSize(0.4);
    hFractionMCTOFpatchedComparisonPions->GetYaxis()->SetTitleOffset(1.0);
    hFractionMCTOFpatchedComparisonPions->GetYaxis()->SetLabelSize(0.03);
    
    TH1F* hFractionMCTOFpatchedComparisonKaons = new TH1F(*hFractionComparisonKaons);
    hFractionMCTOFpatchedComparisonKaons->SetName("hFractionMCTOFpatchedComparisonKaons");
    hFractionMCTOFpatchedComparisonKaons->GetYaxis()->SetTitle(compYtitle.Data());
    hFractionMCTOFpatchedComparisonKaons->GetYaxis()->SetTitleSize(0.4);
    hFractionMCTOFpatchedComparisonKaons->GetYaxis()->SetTitleOffset(1.0);
    hFractionMCTOFpatchedComparisonKaons->GetYaxis()->SetLabelSize(0.03);

    TH1F* hFractionMCTOFpatchedComparisonProtons = new TH1F(*hFractionComparisonProtons);
    hFractionMCTOFpatchedComparisonProtons->SetName("hFractionMCTOFpatchedComparisonProtons");
    hFractionMCTOFpatchedComparisonProtons->GetYaxis()->SetTitle(compYtitle.Data());
    hFractionMCTOFpatchedComparisonProtons->GetYaxis()->SetTitleSize(0.4);
    hFractionMCTOFpatchedComparisonProtons->GetYaxis()->SetTitleOffset(1.0);
    hFractionMCTOFpatchedComparisonProtons->GetYaxis()->SetLabelSize(0.03);

    TH1F* hFractionMCTOFpatchedComparisonMuons = new TH1F(*hFractionComparisonMuons);
    hFractionMCTOFpatchedComparisonMuons->SetName("hFractionMCTOFpatchedComparisonMuons");
    hFractionMCTOFpatchedComparisonMuons->GetYaxis()->SetTitle(compYtitle.Data());
    hFractionMCTOFpatchedComparisonMuons->GetYaxis()->SetTitleSize(0.4);
    hFractionMCTOFpatchedComparisonMuons->GetYaxis()->SetTitleOffset(1.0);
    hFractionMCTOFpatchedComparisonMuons->GetYaxis()->SetLabelSize(0.03);

    TH1F* hFractionMCTOFpatchedComparisonElectrons = new TH1F(*hFractionComparisonElectrons);
    hFractionMCTOFpatchedComparisonElectrons->SetName("hFractionMCTOFpatchedComparisonElectrons");
    hFractionMCTOFpatchedComparisonElectrons->GetYaxis()->SetTitle(compYtitle.Data());
    hFractionMCTOFpatchedComparisonElectrons->GetYaxis()->SetTitleSize(0.4);
    hFractionMCTOFpatchedComparisonElectrons->GetYaxis()->SetTitleOffset(1.0);
    hFractionMCTOFpatchedComparisonElectrons->GetYaxis()->SetLabelSize(0.03);

    
    TCanvas* cFractionMCTOFpatchedComparisons = drawFractionComparisonToMC(mode, pLow, pHigh, hFractionPionsMCTOFpatched,
                                                                           hFractionKaonsMCTOFpatched, hFractionProtonsMCTOFpatched,
                                                                           hFractionElectronsMCTOFpatched, hFractionMuonsMCTOFpatched, hFractionPionsMC, hFractionKaonsMC,
                                                                           hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, 
                                                                           hFractionMCTOFpatchedComparisonPions,
                                                                           hFractionMCTOFpatchedComparisonKaons, 
                                                                           hFractionMCTOFpatchedComparisonProtons,
                                                                           hFractionMCTOFpatchedComparisonElectrons, 
                                                                           hFractionMCTOFpatchedComparisonMuons);
    cFractionMCTOFpatchedComparisons->SetName("cFractionMCTOFpatchedComparisons");
    
    cFractions = drawFinalFractions(mode, pLow, pHigh, electronFixingUsed, hFractionPions, hFractionKaons, hFractionProtons, 
                                    hFractionElectrons, hFractionMuons, hFractionSummed, hFractionPionsMC, hFractionKaonsMC, 
                                    hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, plotIdentifiedSpectra);
    
    cYields = drawFinalYields(mode, pLow, pHigh, hYieldPions, hYieldKaons, hYieldProtons, hYieldElectrons, hYieldMuons,
                              hYieldPionsMC, hYieldKaonsMC, hYieldProtonsMC, hYieldElectronsMC, hYieldMuonsMC,
                              plotIdentifiedSpectra);
    
    cFractionComparisons = drawFractionComparisonToMC(mode, pLow, pHigh, hFractionPions, hFractionKaons, hFractionProtons,
                                                      hFractionElectrons, hFractionMuons, hFractionPionsMC, hFractionKaonsMC,
                                                      hFractionProtonsMC, hFractionElectronsMC, hFractionMuonsMC, 
                                                      hFractionComparisonPions, hFractionComparisonKaons, 
                                                      hFractionComparisonProtons, hFractionComparisonElectrons, 
                                                      hFractionComparisonMuons);
    
    cYieldComparisons = drawYieldComparisonToMC(mode, pLow, pHigh, hYieldPions, hYieldKaons, hYieldProtons,
                                                hYieldElectrons, hYieldMuons, hYieldPionsMC, hYieldKaonsMC,
                                                hYieldProtonsMC, hYieldElectronsMC, hYieldMuonsMC, 
                                                hYieldComparisonPions, hYieldComparisonKaons, 
                                                hYieldComparisonProtons, hYieldComparisonElectrons, 
                                                hYieldComparisonMuons);
    
    if (hYieldPionsMCTOFpatched)
      hYieldPionsMCTOFpatched->Write();
    
    if (hYieldKaonsMCTOFpatched)
      hYieldKaonsMCTOFpatched->Write();
    
    if (hYieldProtonsMCTOFpatched)
      hYieldProtonsMCTOFpatched->Write();
    
    if (hYieldMuonsMCTOFpatched)
      hYieldMuonsMCTOFpatched->Write();
    
    if (hYieldElectronsMCTOFpatched)
      hYieldElectronsMCTOFpatched->Write();
    
    
    if (hFractionPionsMCTOFpatched)
      hFractionPionsMCTOFpatched->Write();
    
    if (hFractionKaonsMCTOFpatched)
      hFractionKaonsMCTOFpatched->Write();
    
    if (hFractionProtonsMCTOFpatched)
      hFractionProtonsMCTOFpatched->Write();
    
    if (hFractionMuonsMCTOFpatched)
      hFractionMuonsMCTOFpatched->Write();
    
    if (hFractionElectronsMCTOFpatched)
      hFractionElectronsMCTOFpatched->Write();
    
    
    if (hRatioToPiKaonsMCTOFpatched)
      hRatioToPiKaonsMCTOFpatched->Write();
    
    if (hRatioToPiProtonsMCTOFpatched)
      hRatioToPiProtonsMCTOFpatched->Write();
    
    if (hRatioToPiMuonsMCTOFpatched)
      hRatioToPiMuonsMCTOFpatched->Write();
    
    if (hRatioToPiElectronsMCTOFpatched)
      hRatioToPiElectronsMCTOFpatched->Write();
    
    if (hFractionSummedMCTOFpatched)
      hFractionSummedMCTOFpatched->Write();
    
    
    if (hFractionMCTOFpatchedComparisonElectrons)
      hFractionMCTOFpatchedComparisonElectrons->Write();
    
    if (hFractionMCTOFpatchedComparisonMuons)
      hFractionMCTOFpatchedComparisonMuons->Write();
    
    if (hFractionMCTOFpatchedComparisonKaons)
      hFractionMCTOFpatchedComparisonKaons->Write();
    
    if (hFractionMCTOFpatchedComparisonPions)
      hFractionMCTOFpatchedComparisonPions->Write();
    
    if (hFractionMCTOFpatchedComparisonProtons)
      hFractionMCTOFpatchedComparisonProtons->Write();
    
    
    if (cFractionMCTOFpatched)
      cFractionMCTOFpatched->Write();
    
    if (cFractionMCTOFpatchedComparisons)
      cFractionMCTOFpatchedComparisons->Write();
    
    if (cFractions)
      cFractions->Write();
    
    if (cYields)
      cYields->Write();
    
    if (cFractionComparisons)
      cFractionComparisons->Write();
    
    if (cYieldComparisons)
      cYieldComparisons->Write();
    
    if (hFractionComparisonElectrons)
      hFractionComparisonElectrons->Write();
    
    if (hFractionComparisonMuons)
      hFractionComparisonMuons->Write();
    
    if (hFractionComparisonKaons)
      hFractionComparisonKaons->Write();
    
    if (hFractionComparisonPions)
      hFractionComparisonPions->Write();
    
    if (hFractionComparisonProtons)
      hFractionComparisonProtons->Write();
    
    if (hYieldComparisonElectrons)
      hYieldComparisonElectrons->Write();
    
    if (hYieldComparisonMuons)
      hYieldComparisonMuons->Write();
    
    if (hYieldComparisonKaons)
      hYieldComparisonKaons->Write();
    
    if (hYieldComparisonPions)
      hYieldComparisonPions->Write();
    
    if (hYieldComparisonProtons)
      hYieldComparisonProtons->Write();
    
    
    if (hFractionSummed)
      hFractionSummed->Write();
    
    if (hFractionElectrons)
      hFractionElectrons->Write();
      
    if (hFractionKaons)
      hFractionKaons->Write();
    
    if (hFractionPions)
      hFractionPions->Write();
    
    if (hFractionProtons)
      hFractionProtons->Write();
    
    if (hFractionMuons)
      hFractionMuons->Write();
    
    if (hFractionElectronsMC)
      hFractionElectronsMC->Write();
    
    if (hFractionKaonsMC)
      hFractionKaonsMC->Write();
    
    if (hFractionPionsMC)
      hFractionPionsMC->Write();
    
    if (hFractionMuonsMC)
      hFractionMuonsMC->Write();
    
    if (hFractionProtonsMC)
      hFractionProtonsMC->Write();
    
    if (hYieldElectrons)
      hYieldElectrons->Write();
    
    if (hYieldKaons)
      hYieldKaons->Write();
    
    if (hYieldPions)
      hYieldPions->Write();
    
    if (hYieldProtons)
      hYieldProtons->Write();
    
    if (hYieldMuons)
      hYieldMuons->Write();
    
    if (hYieldElectronsMC)
      hYieldElectronsMC->Write();
    
    if (hYieldKaonsMC)
      hYieldKaonsMC->Write();
    
    if (hYieldPionsMC)
      hYieldPionsMC->Write();
    
    if (hYieldMuonsMC)
      hYieldMuonsMC->Write();
    
    if (hYieldProtonsMC)
      hYieldProtonsMC->Write();
    
    if (hYieldSummedMC)
      hYieldSummedMC->Write();
    
    if (hRatioToPiElectrons)
      hRatioToPiElectrons->Write();
    
    if (hRatioToPiMuons)
      hRatioToPiMuons->Write();
    
    if (hRatioToPiKaons)
      hRatioToPiKaons->Write();
    
    if (hRatioToPiProtons)
      hRatioToPiProtons->Write();
    
    if (hRatioToPiElectronsMC)
      hRatioToPiElectronsMC->Write();
    
    if (hRatioToPiMuonsMC)
      hRatioToPiMuonsMC->Write();
    
    if (hRatioToPiKaonsMC)
      hRatioToPiKaonsMC->Write();
    
    if (hRatioToPiProtonsMC)
      hRatioToPiProtonsMC->Write();
    
     for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hMCgenYieldsPrimSpecies[i])
        hMCgenYieldsPrimSpecies[i]->Write();
    }
    
    if (hNumEvents)
      hNumEvents->Write();
    
    if (hNumEventsTriggerSel)
      hNumEventsTriggerSel->Write();
    
    if (hNumEventsTriggerSelVtxCut)
      hNumEventsTriggerSelVtxCut->Write();
    
    if (hNumEventsTriggerSelVtxCutNoPileUp)
      hNumEventsTriggerSelVtxCutNoPileUp->Write();
    
    if (hNjetsGen)
      hNjetsGen->Write();
    
    if (hNjetsRec)
      hNjetsRec->Write();
    
    saveFTOF->Close();
    
    // Result path should be set to TOF patched results, if TOF patching is used
    if (filePathNameResults)
      *filePathNameResults = saveFTOFName;
    
    delete cFractions;
    delete cYields;
    delete cFractionComparisons;
    delete cYieldComparisons;
    delete cFractionTOFComparisons;
    delete cFractionMCTOFpatched;
    delete cFractionMCTOFpatchedComparisons;
  }
  
  // Clean up
  delete gFractionElectronsData;
  gFractionElectronsData = 0x0;
  
  delete fElectronFraction;
  fElectronFraction = 0x0;
  
  delete mathFit;
  mathFit = 0x0;
  
  saveF->Close();
  
  delete histList;
  
  TIter next(gDirectory->GetList());
  TObject* obj = 0x0;
  while ( (obj = (TObject*)next()) ) {
    if (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TF1::Class()) || obj->InheritsFrom(THnBase::Class()))
      delete obj;
  }
  
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  
  return 0; 
}
