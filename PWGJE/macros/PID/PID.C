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
#include "AliTPCPIDmathFit.h"

enum processMode { kPMpT = 0, kPMz = 1, kPMxi = 2 };
enum muonTreatment { kNoMuons = 0, kMuonFracEqualElFrac = 1, kMuonFracOverElFracTunedOnMCStandardTrackCuts = 2,
                     kMuonFracOverElFracTunedOnMCHybridTrackCuts = 3, kMuonFracOverElFracTunedOnMCHybridTrackCutsJets = 4,
                     kMuonFracOverElFracTunedOnMCStandardTrackCutsPPb = 5,
                     kNumHandlings = 6 };

const TString modeShortName[3] = { "Pt", "Z", "Xi" };
const TString modeLatexName[3] = { "P_{T}", "z", "#xi" };

const TString muonFractionHandlingShortName[kNumHandlings] =
  { "noMuons", "muonsEqualElectrons", "muonToElTunedOnMCStandardTrackCuts", "muonToElTunedOnMCHybridTrackCuts",
    "muonToElTunedOnMCHybridTrackCutsJets", "muonToElTunedOnMCStandardTrackCutsPPB" };

const Double_t epsilon = 1e-10;
const TString identifiedLabels[2] = { "Most Probable PID", "MC" };
Int_t isMC = 0;

TString minimisationStrategy = "MIGRAD"; // "MINIMIZE"
Bool_t takeIntoAccountMuons = kTRUE;

// 0 = no muons, 1 = muonFrac=elFrac, 2(3) = muonFrac/elFrac tuned on MC for DefaultTrackCuts (hybridTrackCuts),
// 4 = muonFrac/elFrac tuned on MC for hybridTrackCuts for jet particles,
Int_t muonFractionHandling = 3; 


//TODO getErrorOf.... is COMPLETELY wrong now, since the parameter numbering has changed and the muons had come into play!!!!!!

// TODO CAREFUL: fitMethod == 1 adds errors of electrons to pions, but not to muons (should be added to electron error instead!)
const Bool_t muonContamination = kFALSE;//TODO CAREFUL: fitMethod == 1 takes into account the muon contamination in the error calculation!!!

const Bool_t normaliseResults = kTRUE; // Works only for fitMethod == 2

const Bool_t enableShift = kFALSE;
const Int_t dataAxis = kPidDeltaPrime;//kPidDelta; kPidDeltaPrime

const Int_t numSimultaneousFits = 4;

// Upper and lower axis bounds (y-axis) of (data - fit) / data QA histos
const Double_t fitQAaxisLowBound = -0.5;
const Double_t fitQAaxisUpBound = 0.5;

Bool_t useDeltaPrime = (dataAxis == kPidDeltaPrime);

// Will be set later
Double_t muonFractionThresholdForFitting = -1.;
Double_t muonFractionThresholdBinForFitting = -1;
  
Double_t electronFractionThresholdForFitting = -1.;
Double_t electronFractionThresholdBinForFitting = -1;


TF1 fMuonOverElFractionMC("fMuonOverElFractionMC", "[0]+[1]/TMath::Min(x, [4])+[2]*TMath::Min(x, [4])+[3]*TMath::Min(x, [4])*TMath::Min(x, [4])+[5]*TMath::Min(x, [4])*TMath::Min(x, [4])*TMath::Min(x, [4])+[6]*(x>[7])*TMath::Min(x-[7], [8]-[7])",
                          0.01, 50.);

TF1* fElectronFraction = 0x0;
const Double_t lowFittingBoundElectronFraction = 3.0; 

TGraphErrors* gFractionElectronsData = 0x0;
Double_t lastPtForCallOfGetElectronFraction = -1;


//____________________________________________________________________________________________________________________
Double_t GetElectronFraction(const Double_t pT, const Double_t *par)
{
  // During the fit (both, simultaneous and non-simultaneous), the algorithm will always start off from
  // the low pT and go to higher pT. So, it is only necessary to do the fit once the first fixed bin is reached.
  // Then the parameters for the electron fraction remain fixed until the next fit iteration.
  // Since only for the case of regularisation the electron fractions of all x bins are stored in mathFit,
  // the evaluation of this function is done here only in that case (only then the electron fraction will
  // be set to "-pT".
  
  // NOTE 1: Electrons have index 3 per x bin
  // NOTE 2: This function is only called for fitting vs. pT. In that case, xValue holds the LOG of pT!
  
  AliTPCPIDmathFit* mathFit = AliTPCPIDmathFit::Instance();
  
  // lastPtForCallOfGetElectronFraction will be initialised with a value larger than any pT during the fit.
  // So, if this function is called and the pT is smaller than lastPtForCallOfGetElectronFraction, the parameters
  // must have changed and the electron fit needs to be re-done (also see comment above)
  if (pT < lastPtForCallOfGetElectronFraction) {
    for (Int_t xBin = 0; xBin < mathFit->GetNumXbinsRegularisation(); xBin++) {
      
      const Double_t xCoord = TMath::Exp(mathFit->GetXvaluesForRegularisation()[xBin]);
      const Int_t parIndexWithFraction = 3 + xBin * mathFit->GetNumParametersPerXbin(); 
      
      if (xCoord >= lowFittingBoundElectronFraction && xCoord <= electronFractionThresholdForFitting
          && par[parIndexWithFraction] > epsilon) { // Skip zero values (usually due to failed fits)
        gFractionElectronsData->SetPoint(xBin, TMath::Exp(mathFit->GetXvaluesForRegularisation()[xBin]), par[parIndexWithFraction]);
        // Since the errors during the fitting are not reliable, use the following approximation on a statistical basis
        // (which indeed turns out to be rather good!)
        
        // Bin effective weight required for weighted data sets. In case of no weighting, the weight error is sqrt(weight),
        // i.e. effWeight is 1
        const Double_t effWeight = mathFit->GetXstatisticalWeightError()[xBin] * mathFit->GetXstatisticalWeightError()[xBin]
                                   / mathFit->GetXstatisticalWeight()[xBin];
        gFractionElectronsData->SetPointError(xBin, 0, effWeight * TMath::Sqrt(par[parIndexWithFraction] 
                                                                               / mathFit->GetXstatisticalWeight()[xBin]));
      }
      else {
        gFractionElectronsData->SetPoint(xBin, -1, 0);
        gFractionElectronsData->SetPointError(xBin, 0, 0);
      }
    }
    
    gFractionElectronsData->Fit(fElectronFraction, "Ex0NQ", "", lowFittingBoundElectronFraction, electronFractionThresholdForFitting);
  }
  
  lastPtForCallOfGetElectronFraction = pT;
  
  // Catch cases in which the fit function yields invalid fractions (i.e. < 0 or > 1)
  return TMath::Max(0.0, TMath::Min(1.0, fElectronFraction->Eval(pT)));
}


//____________________________________________________________________________________________________________________
Double_t GetElectronFractionError()
{
  // This function estimates the error of the electron fraction for the fixed values via using the parameter errors of
  // the electron fraction function. Note that the parameters (and errors) must be set before calling this function.
  
  // Produce several values via setting the parameters to a random value, which is distributed with a gaussian with mean = parValue
  // and sigma = parError and then take the 2*RMS as the error
  const Int_t nGenValues = 1000;
  Double_t genValues[nGenValues];
  
  const Int_t nPars = fElectronFraction->GetNpar();
  Double_t par[nPars];
  
  TRandom3 rnd(0); // 0 means random seed
  
  const Double_t x = electronFractionThresholdForFitting + 1.; // Some value above the threshold to obtain a fixed value
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
  
  if (fractionA < epsilon) {
    ratio = 0.;
    ratioError = 999.;
    
    return;
  }
  
  ratio = fractionA / fractionB;
  
  const Double_t x = ratio / fractionA;
  const Double_t y = -ratio / fractionB;
  
  // covMatrixElement(i, i) = error(i)^2
  ratioError = GetCorrelatedError(x, y, fractionErrorA * fractionErrorA, fractionErrorB * fractionErrorB, covMatrixElementAB); 
  
  //printf("frationA %e\nfractionB %e\nfractionErrorA %e\nfractionErrorB %e\ncovMatrixElementAB %e\nratio %e\nx %e\ny %e\nratioError %e\n\n",
  //       fractionA, fractionB, fractionErrorA, fractionErrorB, covMatrixElementAB, ratio, x, y, ratioError);
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
}

//____________________________________________________________________________________________________________________
void SetReasonableXaxisRange(TH1* h, Int_t& binLow, Int_t& binHigh)
{
  binLow = TMath::Max(1, h->FindFirstBinAbove(0));
  binHigh  = TMath::Min(h->GetNbinsX(), h->FindLastBinAbove(0));
  
  h->GetXaxis()->SetRange(binLow, binHigh);
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNoExponent(kTRUE);
}


//____________________________________________________________________________________________________________________
Int_t FindMomentumBin(const Double_t* pTbins, const Double_t value, const Int_t numPtBins = nPtBins)
{
  for (Int_t bin = 0; bin < numPtBins; bin++) {
    if (value >= pTbins[bin] && value < pTbins[bin + 1])
      return bin;
  }
  
  return -1;
}


//____________________________________________________________________________________________________________________
Double_t normaliseHist(TH1* h, Double_t scaleFactor = -1)
{
  // Scales the histogram with the scale factor. If the scale factor is < 0,
  // the histogram is normalised to it's integral.
  // In both cases, the normalisation factor is returned.
  
  Double_t normFactor = 1.;
  
  if (scaleFactor < 0) {
    Double_t integralTemp = h->Integral();
    if (integralTemp > 0) {
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
void normaliseYieldHist(TH1* h, Double_t numEvents, Double_t deta)
{
  // Yield histos are already normalised to dpT. Now normalise to 1/NeV 1/(2pi pT) 1/deta in addition
  
  if (numEvents <= 0) // Do not normalise
    numEvents = 1; 
  
  for (Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
    Double_t normFactor = 1. / (numEvents * 2 * TMath::Pi() * h->GetXaxis()->GetBinCenter(bin) * deta);
    h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
    h->SetBinError(bin, h->GetBinError(bin) * normFactor);
  }
}


//____________________________________________________________________________________________________________________
void normaliseGenYieldMCtruthHist(TH1* h, Double_t numEvents, Double_t deta)
{
  // Yield histos are NOT normalised to dpT. Now normalise to 1/NeV 1/(2pi pT) 1/deta 1/dpT!
  
  if (numEvents <= 0) // Do not normalise
    numEvents = 1; 
  
  for (Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
    Double_t normFactor = 1. / (numEvents * 2 * TMath::Pi() * h->GetXaxis()->GetBinCenter(bin) * h->GetXaxis()->GetBinWidth(bin) * deta);
    h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
    h->SetBinError(bin, h->GetBinError(bin) * normFactor);
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
  const Double_t parElFraction = (par[parEl] < 0) ? GetElectronFraction(-par[parEl], par) : par[parEl];
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
  const Double_t parElFraction = (par[3] < 0) ? GetElectronFraction(-par[3], par) : par[3];
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
  const Double_t scaleFactorEl = par[parAll] * ((par[parEl] < 0) ? GetElectronFraction(-par[parEl], par) : par[parEl]);
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
  const Double_t scaleFactorEl = par[5] * ((par[3] < 0) ? GetElectronFraction(-par[3], par) : par[3]);
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
                                   Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits, Double_t& reducedChiSquare)
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
        gausParams[parIndex]         = GetElectronFraction(-gausParams[parIndex], &gausParams[0]);
        parameterErrorsOut[parIndex] = GetElectronFractionError();
      }
      // Set muon fraction equal to electron fraction (or some modified electron fraction) above some threshold,
      // which should be a reasonable approximation:
      // Fixed muon fraction < 0 does this job within the fitting functions
      else if (parIndexModulo == 4 && gausParams[parIndex] < 0) {
        gausParams[parIndex]         = GetMuonFractionFromElectronFractionAndPt(-gausParams[parIndex], gausParams[parIndex - 1]);
        parameterErrorsOut[parIndex] = parameterErrorsOut[parIndex - 1];
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
                        Double_t* covMatrix, Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits, Double_t& 
                        reducedChiSquare)
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
            Double_t* stepSize, Double_t* lowParLimits, Double_t* upParLimits, TF1* totalDeltaSpecies, Double_t& reducedChiSquare)
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
Double_t setFractionsAndYields(Int_t slice, Double_t inverseBinWidth, Double_t binWidthFitHisto, Int_t species, Double_t* parametersOut, 
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
  
  Double_t sumOfParticles = inverseBinWidth * parametersOut[5] / binWidthFitHisto; // Divide by binWidthFitHisto since parametersOut includes this width
  
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
Int_t PID(TString fileName, Double_t deta, Double_t pLow, Double_t pHigh, Bool_t isMCdataSet, Int_t fitMethod, 
          Int_t muonFractionHandlingParameter, //0 = no muons, 1 = muonFrac=elFrac,
                                               //2(3) = muonFrac/elFrac tuned on MC for StandardTrackCuts(HybridTrackCuts)
          Bool_t useIdentifiedGeneratedSpectra, Bool_t plotIdentifiedSpectra, Int_t mode/*0=pT,1=z,2=xi*/,
          Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
          Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
          Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/,
          Int_t rebin/* = 1 -> DON'T USE FOR PT (will not work since binsPt will and should be used!)*/,
          Int_t rebinDeltaPrime/* = 1*/,
          TString listName /* = "bhess_PID"*/,
          Bool_t useLogLikelihood /*= kTRUE*/, Bool_t useWeightsForLogLikelihood /*= kFALSE*/,
          Int_t regularisation /*= 0*/,
          Double_t regularisationFactor /*= 1*/,
          TString filePathNameFileWithInititalFractions /*= ""*/,
          TString* filePathNameResults /*= 0x0*/) 
{
  // Do all the fitting
  
  isMC = isMCdataSet;
  
  muonFractionHandling = muonFractionHandlingParameter;
  
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
  
  std::cout << "Fitting \"" << fileName.Data() << "\" with settings:" << std::endl;
  
  std::cout << "Minimisation strategy: " << minimisationStrategy.Data() << std::endl;
  if (useLogLikelihood) 
    std::cout << "Binned loglikelihood fit" << (useWeightsForLogLikelihood ? " (weighted)" : "")  << std::endl;
  else
    std::cout << "ChiSquare fit" << std::endl;
  std::cout << "Processing mode: ";
  if (mode == kPMpT)
    std::cout << "pT" << std::endl;
  else if (mode == kPMz) {
    std::cout << "z" << std::endl;
    axisForMode = kPidZ;
    axisGenForMode = kPidGenZ;
  }
  else if (mode == kPMxi) {
    std::cout << "xi" << std::endl;
    axisForMode = kPidXi;
    axisGenForMode = kPidGenXi;
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
  Int_t upperCentralityBinLimit = -1;
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= hPIDdata->GetAxis(kPidCentrality)->GetNbins()) {
      actualLowerCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
  }
  else {
    std::cout << "All" << std::endl;
  }
    
  if (restrictCentralityAxis) {
    hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  }
  
  
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -1;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 && upperJetPtBinLimit <= hPIDdata->GetAxis(kPidJetPt)->GetNbins()) {
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
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(-1. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimitData <= upperChargeBinLimitData && lowerChargeBinLimitData >= 1
        && upperChargeBinLimitData <= hPIDdata->GetAxis(indexChargeAxisData)->GetNbins()) {
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
  
  saveInterFName = Form("%s_Projections_%s_%d_%s%s%s%s%s.root", saveInterFName.ReplaceAll(".root", "").Data(), 
                        modeShortName[mode].Data(),
                        fitMethod, muonFractionHandlingShortName[muonFractionHandlingParameter].Data(),
                        useIdentifiedGeneratedSpectra ? "_idSpectra" : "",
                        restrictCentralityAxis ? Form("_centrality%.0f_%.0f", actualLowerCentrality, actualUpperCentrality) : "",
                        restrictJetPtAxis ? Form("_jetPt%.1f_%.1f", actualLowerJetPt, actualUpperJetPt) : "",
                        chargeString.Data());
  TFile *saveInterF = TFile::Open(saveInterFName.Data(), "RECREATE");
  saveInterF->cd();

  // TH1 hist with number of processed events
  Double_t numEvents = -1;
  TH1* hNumEvents = dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessed"));
  if (!hNumEvents) {
    std::cout << std::endl;
    std::cout << "Histo with number of processed events not found! Yields will NOT be normalised to this number!" << std::endl 
              << std::endl;
  }
  else {
    numEvents = restrictCentralityAxis ? hNumEvents->Integral(lowerCentralityBinLimit, upperCentralityBinLimit) : 
                                         hNumEvents->Integral();
    
    if (numEvents <= 0) {
      numEvents = -1;
      std::cout << std::endl;
      std::cout << "Number of processed events < 1 in selected range! Yields will NOT be normalised to this number!"
                << std::endl << std::endl;
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
    h2Delta[i]->GetYaxis()->SetTitle(Form("#Delta%s_{%s} = dE/dx %s <dE/dx>_{%s} (arb. units)", useDeltaPrime ? "'" : "", speciesLabel.Data(),
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
    const Int_t pBinLowProjLimit = (mode == kPMpT) ? h2Delta[0]->GetXaxis()->FindBin(binsPt[slice] + 1e-5)    : slice + 1;
    const Int_t pBinUpProjLimit  = (mode == kPMpT) ? h2Delta[0]->GetXaxis()->FindBin(binsPt[slice + 1]- 1e-5) : slice + 1;
    
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
    
    if (plotIdentifiedSpectra)  {
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
    }
  }
  hPIDdata->GetAxis(kPidMCpid)->SetRange(0, -1);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1);
  
  hPIDdata->GetAxis(axisForMode)->SetRange(0, -1);

  /*
  // TOF TODO
  TCanvas* cTOF = new TCanvas("cTOF", "TOF PID",100,10,1200,800);
  cTOF->Divide(4,1);
  cTOF->cd(1);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1);
  TH2D* h2TOFel = hPIDdata->Projection(kPidDeltaTOF, kPidPvertex);
  h2TOFel->SetName("h2TOFel");
  h2TOFel->GetYaxis()->SetTitle(Form("#DeltaTOF_{%s} (ps)", hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(1)));
  h2TOFel->Draw("colz");
  
  cTOF->cd(2);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(2, 2);
  TH2D* h2TOFka = hPIDdata->Projection(kPidDeltaTOF, kPidPvertex);
  h2TOFka->SetName("h2TOFka");
  h2TOFka->GetYaxis()->SetTitle(Form("#DeltaTOF_{%s} (ps)", hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(2)));
  h2TOFka->Draw("colz");
  
  cTOF->cd(3);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(3, 3);
  TH2D* h2TOFpi = hPIDdata->Projection(kPidDeltaTOF, kPidPvertex);
  h2TOFpi->SetName("h2TOFpi");
  h2TOFpi->GetYaxis()->SetTitle(Form("#DeltaTOF_{%s} (ps)", hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(3)));
  h2TOFpi->Draw("colz");
  
  cTOF->cd(4);
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(4, 4);
  TH2D* h2TOFpr = hPIDdata->Projection(kPidDeltaTOF, kPidPvertex);
  h2TOFpr->SetName("h2TOFpr");
  h2TOFpr->GetYaxis()->SetTitle(Form("#DeltaTOF_{%s} (ps)", hPIDdata->GetAxis(kPidSelectSpecies)->GetBinLabel(4)));
  h2TOFpr->Draw("colz");
  
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(0, -1);
  */
  
  // Start fitting of slices
  TCanvas* cSingleFit[numSlices][4];
  
  TF1*     fitFuncTotalDeltaPion[numSlices];
  TF1*     fitFuncTotalDeltaElectron[numSlices];
  TF1*     fitFuncTotalDeltaKaon[numSlices];
  TF1*     fitFuncTotalDeltaProton[numSlices];
  
  // Histos for particle fractions
  TH1F* hFractionElectrons = 0x0;
  if (mode == kPMpT)
    hFractionElectrons = new TH1F("hFractionElectrons", "e", nPtBins, binsPt);
  else {
    const TArrayD* histBins = hPIDdata->GetAxis(axisForMode)->GetXbins();
    if (histBins->fN == 0)
      hFractionElectrons = new TH1F("hFractionElectrons", "e", hPIDdata->GetAxis(axisForMode)->GetNbins(), hPIDdata->GetAxis(axisForMode)->GetXmin(),
                                    hPIDdata->GetAxis(axisForMode)->GetXmax());
    else
      hFractionElectrons = new TH1F("hFractionElectrons", "e", hPIDdata->GetAxis(axisForMode)->GetNbins(), histBins->fArray);
  }
  
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
  TH1F* hYieldElectrons = 0x0;
  if (mode == kPMpT)
    hYieldElectrons = new TH1F("hYieldElectrons", "e", nPtBins, binsPt);
  else {
    const TArrayD* histBins = hPIDdata->GetAxis(axisForMode)->GetXbins();
    if (histBins->fN == 0)
      hYieldElectrons = new TH1F("hYieldElectrons", "e", hPIDdata->GetAxis(axisForMode)->GetNbins(), hPIDdata->GetAxis(axisForMode)->GetXmin(),
                                    hPIDdata->GetAxis(axisForMode)->GetXmax());
    else
      hYieldElectrons = new TH1F("hYieldElectrons", "e", hPIDdata->GetAxis(axisForMode)->GetNbins(), histBins->fArray);
  }
  
  hYieldElectrons->GetXaxis()->SetTitle(hPIDdata->GetAxis(axisForMode)->GetTitle());
  hYieldElectrons->GetYaxis()->SetTitle(Form("%s1/(2#pi%s) d^{2}N/d#etad%s%s", numEvents > 0 ? "1/N_{ev} " : "",
                                             modeLatexName[mode].Data(), modeLatexName[mode].Data(),
                                             mode == kPMpT ? " (GeV/c)^{-2}" : 0));
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
  hRatioToPiElectrons->GetYaxis()->SetTitle(Form("d^{2}N_{%s}/d#etad%s / d^{2}N_{%s}/d#etad%s",
                                                 electronString[chargeMode+1].Data(),
                                                 modeLatexName[mode].Data(),
                                                 pionString[chargeMode+1].Data(),
                                                 modeLatexName[mode].Data()));
  hRatioToPiElectrons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", electronString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", electronString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  
  TH1F* hRatioToPiMuons = (TH1F*)hYieldMuons->Clone("hRatioToPiMuons");
  hRatioToPiMuons->GetYaxis()->SetTitle(Form("d^{2}N_{%s}/d#etad%s / d^{2}N_{%s}/d#etad%s",
                                             muonString[chargeMode+1].Data(),
                                             modeLatexName[mode].Data(),
                                             pionString[chargeMode+1].Data(),
                                             modeLatexName[mode].Data()));
  hRatioToPiMuons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", muonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", muonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  TH1F* hRatioToPiKaons = (TH1F*)hYieldKaons->Clone("hRatioToPiKaons");
  hRatioToPiKaons->GetYaxis()->SetTitle(Form("d^{2}N_{%s}/d#etad%s / d^{2}N_{%s}/d#etad%s",
                                             kaonString[chargeMode+1].Data(),
                                             modeLatexName[mode].Data(),
                                             pionString[chargeMode+1].Data(),
                                             modeLatexName[mode].Data()));
  hRatioToPiKaons->SetTitle(Form("%s", chargeMode == 0
                                             ? Form("(%s)/(%s)", kaonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())
                                             : Form("%s/%s", kaonString[chargeMode+1].Data(), pionString[chargeMode+1].Data())));
  
  TH1F* hRatioToPiProtons = (TH1F*)hYieldProtons->Clone("hRatioToPiProtons");
  hRatioToPiProtons->GetYaxis()->SetTitle(Form("d^{2}N_{%s}/d#etad%s / d^{2}N_{%s}/d#etad%s",
                                             protonString[chargeMode+1].Data(),
                                             modeLatexName[mode].Data(),
                                             pionString[chargeMode+1].Data(),
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

  
  
  // Extract the MC truth generated primary yields
  THnSparse* hMCgeneratedYieldsPrimaries = isMCdataSet ? dynamic_cast<THnSparse*>(histList->FindObject("fhMCgeneratedYieldsPrimaries"))
                                                       : 0x0;
  
  TH1D* hMCgenYieldsPrimSpecies[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    hMCgenYieldsPrimSpecies[i] = 0x0;
  
  if (hMCgeneratedYieldsPrimaries) {
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
    
    if (restrictCentralityAxis)
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
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(-1. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimitGenYield <= upperChargeBinLimitGenYield && lowerChargeBinLimitGenYield >= 1
          && upperChargeBinLimitGenYield <= hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetNbins()) {
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
      
      hMCgenYieldsPrimSpecies[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldPt, "e");
      hMCgenYieldsPrimSpecies[MCid]->SetName(Form("hMCgenYieldsPrimSpecies_%s", AliPID::ParticleShortName(MCid)));
      hMCgenYieldsPrimSpecies[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      
      // Choose the same binning as for the fitted yields, i.e. rebin the histogram (Rebin will create a clone!)
      TH1D* temp = (TH1D*)hMCgenYieldsPrimSpecies[MCid]->Rebin(nPtBins, hMCgenYieldsPrimSpecies[MCid]->GetName(), binsPt);
      // Delete the old binned histo and take the new binned one
      delete hMCgenYieldsPrimSpecies[MCid];
      hMCgenYieldsPrimSpecies[MCid] = temp;
      
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(0, -1);
    }
  }
  
  
  // Get expected shapes for pT bins
  TString Ytitle = "";
  
  // Array index 0 as unused dummy
  TH2D* hGenDelta[6][6]; // DeltaSpecies (first index) for species (second index)
  TH2D* hGenDeltaMCid[6][6]; // DeltaSpecies (first index) for species (second index)
  
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 6; j++) {
      hGenDelta[i][j] = 0x0;
      hGenDeltaMCid[i][j] = 0x0;
    }
  }
  
  THnSparse* current = 0x0;
  
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
    std::cout << "Failed to load expected dEdx signal shape for: Muons! Treated muons as pions in the following." << std::endl;
    takeIntoAccountMuons = kFALSE; 
  }
  
  THnSparse* hGenPr = dynamic_cast<THnSparse*>(histList->FindObject("hGenPr"));
  if (!hGenPr)  {
    std::cout << "Failed to load expected dEdx signal shape for: Protons!" << std::endl;
    return -1;
  }

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
    if (restrictCentralityAxis) {
      current->GetAxis(kPidGenCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    }
    
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
        lowerChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindBin(-1. + 0.001);
        upperChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindBin(0. + 0.001);
        upperChargeBinLimitGen = current->GetAxis(indexChargeAxisGen)->FindBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimitGen <= upperChargeBinLimitGen && lowerChargeBinLimitGen >= 1
          && upperChargeBinLimitGen <= current->GetAxis(indexChargeAxisGen)->GetNbins()) {
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
    
    
  
    for (Int_t selectBin = 1; selectBin <= 4; selectBin++)  {
      Int_t selectMCid = (selectBin >= 3) ? (selectBin+1) : selectBin;

      current->GetAxis(kPidGenSelectSpecies)->SetRange(selectBin, selectBin);
      
      Ytitle = Form("#Delta%s_{%s} = dE/dx %s <dE/dx>_{%s} (arb. units)", useDeltaPrime ? "'" : "",
                    partShortName[selectMCid - 1].Data(), useDeltaPrime ? "/" : "-",
                    partShortName[selectMCid - 1].Data());
      
      TH2* hGenCurrent = 0x0;
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
    }
    
    current->GetAxis(kPidGenSelectSpecies)->SetRange(0, -1);
  }    
  
  // Free a lot of memory for the following procedure. Histogram is not needed anymore (only its projections)
  delete f;
  
  // Save intermediate results
  //TODO save intermediate TOF results
  saveInterF->cd();
  
  if (hMCdata)
    hMCdata->Write();
  
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 6; j++) {
      if (hGenDelta[i][j])
        hGenDelta[i][j]->Write();
        
      if (hGenDeltaMCid[i][j])
        hGenDeltaMCid[i][j]->Write();
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
      
    if (plotIdentifiedSpectra)  {
      for (Int_t species = 0; species < 5; species++) {
        hDeltaElMC[slice][species]->Write();
        hDeltaKaMC[slice][species]->Write();
        hDeltaPiMC[slice][species]->Write();
        hDeltaPrMC[slice][species]->Write(); 
      }
    }   
  }
  
  // File may not be closed because the projections are needed in the following!
  //saveInterF->Close();
  
  // Save some first results for the final output
  TString saveFName = fileName;
  saveFName = Form("%s_results_%s__%s_%d_reg%d_regFac%.2f_%s%s%s%s%s.root", saveFName.ReplaceAll(".root", "").Data(), 
                   useLogLikelihood ? (useWeightsForLogLikelihood ? "weightedLLFit" : "LLFit") : "ChiSquareFit",
                   modeShortName[mode].Data(), fitMethod, regularisation, regularisationFactor,
                   muonFractionHandlingShortName[muonFractionHandlingParameter].Data(),
                   useIdentifiedGeneratedSpectra ? "_idSpectra" : "",
                   restrictCentralityAxis ? Form("_centrality%.0f_%.0f", actualLowerCentrality, actualUpperCentrality) : "",
                   restrictJetPtAxis ? Form("_jetPt%.1f_%.1f", actualLowerJetPt, actualUpperJetPt) : "",
                   chargeString.Data());
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
  
  
  muonFractionThresholdForFitting = 0.;//OLD 0.295;
  muonFractionThresholdBinForFitting = (mode == kPMpT) ? FindMomentumBin(binsPt, muonFractionThresholdForFitting) : -1;
  
  electronFractionThresholdForFitting = 9.;
  electronFractionThresholdBinForFitting = (mode == kPMpT) ? FindMomentumBin(binsPt, electronFractionThresholdForFitting) : -1;
  
  lastPtForCallOfGetElectronFraction = pHigh + 10.; // Make sure that this value is higher than in any call during the fit
  
  fElectronFraction = new TF1("fElectronFraction", Form("[0]+(x<%f)*[1]*(x-%f)", electronFractionThresholdForFitting, 
                                                        electronFractionThresholdForFitting), 
                              pLow, pHigh);
  fElectronFraction->SetParameters(0.01, 0.0);
  
  
  
  TString speciesLabel[4] = {"El", "Ka", "Pi", "Pr" };
  
  const Double_t binWidthFitHisto = 1.0; // Not used any longer
  
  // In case of regularisation, the actual number of x bins and the (for pT: logs of their) bin centres are required
  Int_t numXBins = 0;
  Double_t* xBinCentres = 0x0;
  Double_t* xBinStatisticalWeight = 0x0;
  Double_t* xBinStatisticalWeightError = 0x0;
  
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
        
        indexParametersToRegularise = new Int_t[nParToRegulariseSimultaneousFit];
        
        lastNotFixedIndexOfParameters = new Int_t[numParamsPerXbin];
        
        // Set last not fixed index of parameter to numXBins, i.e. a index larger than any existing index.
        // This will not restrict the parameter regularisation range. In the following, the range for electrons
        // and muons will be restricted
        for (Int_t iPar = 0; iPar < numParamsPerXbin; iPar++) 
          lastNotFixedIndexOfParameters[iPar] = numXBins;
      }
      
        
      for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < hFractionPions->GetXaxis()->GetNbins(); slice++) {   
        if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
          continue; 
        
        // There won't (actually: shouldn't) be tracks with a pT larger than the jet pT
        if (mode == kPMpT && restrictJetPtAxis && binsPt[slice] >= actualUpperJetPt)
          continue;
        
        const Int_t pBinLowProjLimit = (mode == kPMpT) ? hYieldPt->GetXaxis()->FindBin(binsPt[slice] + 1e-5)    : slice + 1;
        const Int_t pBinUpProjLimit  = (mode == kPMpT) ? hYieldPt->GetXaxis()->FindBin(binsPt[slice + 1]- 1e-5) : slice + 1;
        
        // NOTE: In case of regularisation, only the simultaneous fit values will be used, i.e. totalYield and not allDeltaSpecies!
        
        // Also take into account bin width in delta(Prime) plots -> Multiply by binWidthFitHisto
        Double_t totalYieldError = 0;
        const Double_t totalYield = binWidthFitHisto * hYieldPt->IntegralAndError(pBinLowProjLimit, pBinUpProjLimit, totalYieldError);
        totalYieldError *= binWidthFitHisto;
        
        if (totalYield <= 0) 
          continue;
        
        if (i == 1) {
          if (mode == kPMpT)
            // Take the logarithm in case of pT
            xBinCentres[xBinIndexTemp] = TMath::Log((binsPt[slice + 1] + binsPt[slice]) / 2.);
          else
            xBinCentres[xBinIndexTemp] = hFractionPions->GetXaxis()->GetBinCenter(slice + 1);
          
          xBinStatisticalWeight[xBinIndexTemp] = totalYield;
          
          // NOTE: The total yield is a fact - a number w/o error. However, one assigns this error here to use it
          // to calculate the effective weighting for the weighted likelihood fit (and the error is only used for this).
          // So, it is more like a weighting than an error...
          xBinStatisticalWeightError[xBinIndexTemp] = totalYieldError;
          
          
          // Mark the fractions for all species except for electrons and muons in this bin for regularisation
          for (Int_t speciesIndex = 0; speciesIndex < AliPID::kSPECIES - 2; speciesIndex++)
            indexParametersToRegularise[internalParIndexTemp++] = numParamsPerXbin * xBinIndexTemp + speciesIndex;
          
          // Also mark electrons for regularisation in this bin, if not fixed
          if( !(mode == kPMpT && slice >= electronFractionThresholdBinForFitting) ) {
            indexParametersToRegularise[internalParIndexTemp++] = numParamsPerXbin * xBinIndexTemp + 3; 
          }
          else {
            // Set the index of the last x bin in which the parameter is not fixed.
            // If the parameter is fixed in all x bins, this index will be -1.
            if (xBinIndexTemp - 1 < lastNotFixedIndexOfParameters[3])
              lastNotFixedIndexOfParameters[3] = xBinIndexTemp - 1;
          }
          
          // Also mark muons for regularisation in this bin, if not fixed
          if( !(mode != kPMpT || slice >= muonFractionThresholdBinForFitting) ) {
            indexParametersToRegularise[internalParIndexTemp++] = numParamsPerXbin * xBinIndexTemp + 4; 
          }
          else {
            // Set the index of the last x bin in which the parameter is not fixed.
            // If the parameter is fixed in all x bins, this index will be -1.
            if (xBinIndexTemp - 1 < lastNotFixedIndexOfParameters[4])
              lastNotFixedIndexOfParameters[4] = xBinIndexTemp - 1;
          }
          
          xBinIndexTemp++;
        }
        
        if (i == 0) {
          nParToRegulariseSimultaneousFit += AliPID::kSPECIES - 2; // Fracs for all species in this bin except for electrons and muons

          if( !(mode == kPMpT && slice >= electronFractionThresholdBinForFitting) )
            nParToRegulariseSimultaneousFit++; // Also regularise electrons in this bin (not fixed)
          
          if( !(mode != kPMpT || slice >= muonFractionThresholdBinForFitting) )
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
    mathFit->SetScaleFactorError(2.);
  }
  
  mathFit->SetRegularisation(regularisation, regularisationFactor);
  
  // Number of parameters for fitting
  const Int_t nPar = 11;
  
  // Fracs of each species + total yield in x bin
  const Int_t nParSimultaneousFit = AliPID::kSPECIES + 1; 
  
  // Fracs of each species in x bin + tot yield in x bin
  const Int_t nParSimultaneousFitRegularised = numXBins * (AliPID::kSPECIES + 1); 
  
  if (regularisation > 0) {
    if (!mathFit->SetParametersToRegularise(nParToRegulariseSimultaneousFit, numParamsPerXbin, indexParametersToRegularise,
                                            lastNotFixedIndexOfParameters, xBinCentres, xBinStatisticalWeight, 
                                            xBinStatisticalWeightError))
      return -1;
  }
  
  delete xBinCentres;
  xBinCentres = 0x0;
  
  delete xBinStatisticalWeight;
  xBinStatisticalWeight = 0x0;
  
  delete xBinStatisticalWeightError;
  xBinStatisticalWeight = 0x0;
  
  delete indexParametersToRegularise;
  indexParametersToRegularise = 0x0;
  
  delete lastNotFixedIndexOfParameters;
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
  
  mathFit->ClearRefHistos();
  
  
  const Int_t nParUsed = (fitMethod == 2) ? ((regularisation <= 0) ? nParSimultaneousFit: nParSimultaneousFitRegularised) : nPar;
  Double_t parameterErrorsOut[nParUsed];
  Double_t covMatrix[nParUsed][nParUsed];
  
  for (Int_t iter = 0; iter < 2; iter++) {
    if (regularisation <= 0 && iter == 0)
      continue; // Only one iteration w/o regularisation
  
    Int_t currXbin = 0;
    
    for (Int_t slice = 0; (mode == kPMpT) ? slice < nPtBins : slice < hFractionPions->GetXaxis()->GetNbins(); slice++) {   
      if (mode == kPMpT && (slice < pSliceLow || slice > pSliceHigh))
        continue; 
      
      // There won't (actually: shouldn't) be tracks with a pT larger than the jet pT
      if (mode == kPMpT && restrictJetPtAxis && binsPt[slice] >= actualUpperJetPt)
        continue;
      
      if (regularisation <= 0) {
        if (mode == kPMpT)
          std::cout << "Fitting range " << binsPt[slice] << " GeV/c < Pt < " << binsPt[slice + 1] << " GeV/c..." << std::endl;
        else {
          std::cout << "Fitting range " << hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1) << " < " << modeShortName[mode].Data() << " < ";
          std::cout << hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1) << "..." << std::endl;
        }
      }
      
      // Add/subtract some very small offset to be sure not to sit on the bin boundary, when looking for the integration/projection limits.
      const Int_t pBinLowProjLimit = (mode == kPMpT) ? hYieldPt->GetXaxis()->FindBin(binsPt[slice] + 1e-5)    : slice + 1;
      const Int_t pBinUpProjLimit  = (mode == kPMpT) ? hYieldPt->GetXaxis()->FindBin(binsPt[slice + 1]- 1e-5) : slice + 1;
      
      // Also take into account bin width in delta(Prime) plots -> Multiply by binWidthFitHisto
      const Double_t totalYield = binWidthFitHisto * hYieldPt->Integral(pBinLowProjLimit, pBinUpProjLimit);
      
      if (totalYield <= 0) {
        std::cout << "Skipped bin (yield is zero)!" << std::endl;
        continue;
      }
      
      const Double_t allDeltaPion = hDeltaPi[slice]->Integral();
      const Double_t allDeltaElectron = hDeltaEl[slice]->Integral();
      const Double_t allDeltaKaon = hDeltaKa[slice]->Integral();
      const Double_t allDeltaProton = hDeltaPr[slice]->Integral();
      
      // inverseBinWidth = 1.0, if the raw yield for each bin is requested.
      // If divided by the bin size, the histograms give "yield per unit pT in the corresponding bin" or dN/dpT
      Double_t inverseBinWidth = (mode == kPMpT) ? 1.0 / (binsPt[slice + 1] - binsPt[slice])
                                                : 1.0 / hYieldPt->GetBinWidth(slice + 1); 
      
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
      
      hGenDeltaElForElProj =(TH1D*)hGenDeltaUsed[kEl][kEl]->ProjectionY(Form("hGenDeltaElForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaKaForElProj =(TH1D*)hGenDeltaUsed[kKa][kEl]->ProjectionY(Form("hGenDeltaKaForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPiForElProj =(TH1D*)hGenDeltaUsed[kPi][kEl]->ProjectionY(Form("hGenDeltaPiForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPrForElProj =(TH1D*)hGenDeltaUsed[kPr][kEl]->ProjectionY(Form("hGenDeltaPrForElProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
      hGenDeltaElForKaProj =(TH1D*)hGenDeltaUsed[kEl][kKa]->ProjectionY(Form("hGenDeltaElForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaKaForKaProj =(TH1D*)hGenDeltaUsed[kKa][kKa]->ProjectionY(Form("hGenDeltaKaForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPiForKaProj =(TH1D*)hGenDeltaUsed[kPi][kKa]->ProjectionY(Form("hGenDeltaPiForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPrForKaProj =(TH1D*)hGenDeltaUsed[kPr][kKa]->ProjectionY(Form("hGenDeltaPrForKaProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
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
      
      hGenDeltaElForPrProj =(TH1D*)hGenDeltaUsed[kEl][kPr]->ProjectionY(Form("hGenDeltaElForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaKaForPrProj =(TH1D*)hGenDeltaUsed[kKa][kPr]->ProjectionY(Form("hGenDeltaKaForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPiForPrProj =(TH1D*)hGenDeltaUsed[kPi][kPr]->ProjectionY(Form("hGenDeltaPiForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      hGenDeltaPrForPrProj =(TH1D*)hGenDeltaUsed[kPr][kPr]->ProjectionY(Form("hGenDeltaPrForPrProj%d", slice), pBinLowProjLimit, pBinUpProjLimit, "e");
      
      
      
      if (fitMethod == 2) {
        // Normalise generated histos to TOTAL number of GENERATED particles for this species (i.e. including
        // entries that lie in the under or overflow bin), so that situations in which the generated spectra lie
        // at least partly outside the histo are treated properly. To find the total number of generated particle
        // species X, one can just take the integral of the generated histo for DeltaX (which should include all
        // generated entries) and apply the same normalisation factor to all other DeltaSpecies.
        // Also set some cosmetics
        
        // Generated electrons
        Double_t normEl = normaliseHist(hGenDeltaElForElProj, -1);
        normaliseHist(hGenDeltaKaForElProj, normEl);
        normaliseHist(hGenDeltaPiForElProj, normEl);
        normaliseHist(hGenDeltaPrForElProj, normEl);
        
        
        // Generated kaons
        Double_t normKa = normaliseHist(hGenDeltaKaForKaProj, -1);
        normaliseHist(hGenDeltaElForKaProj, normKa);
        normaliseHist(hGenDeltaPiForKaProj, normKa);
        normaliseHist(hGenDeltaPrForKaProj, normKa);
        
        
        // Generated pions
        Double_t normPi = normaliseHist(hGenDeltaPiForPiProj, -1);
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
      }
      else {
        // Normalise generated histos to total number of particles for this delta
        // and also set some cosmetics
        
        // DeltaEl
        normaliseHist(hGenDeltaElForElProj);
        normaliseHist(hGenDeltaElForKaProj);
        normaliseHist(hGenDeltaElForPiProj);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaElForMuProj);
        normaliseHist(hGenDeltaElForPrProj);    
        
        // DeltaKa
        normaliseHist(hGenDeltaKaForElProj);
        normaliseHist(hGenDeltaKaForKaProj);
        normaliseHist(hGenDeltaKaForPiProj);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaKaForMuProj);
        normaliseHist(hGenDeltaKaForPrProj);
        
        // DeltaPi
        normaliseHist(hGenDeltaPiForElProj);
        normaliseHist(hGenDeltaPiForKaProj);
        normaliseHist(hGenDeltaPiForPiProj);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaPiForMuProj);
        normaliseHist(hGenDeltaPiForPrProj);
        
        // DeltaPr
        normaliseHist(hGenDeltaPrForElProj);
        normaliseHist(hGenDeltaPrForKaProj);
        normaliseHist(hGenDeltaPrForPiProj);
        if (takeIntoAccountMuons)
          normaliseHist(hGenDeltaPrForMuProj);
        normaliseHist(hGenDeltaPrForPrProj);
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
        
        if (takeIntoAccountMuons)
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
        Double_t xBinCentre = (mode == kPMpT) ? (binsPt[slice + 1] + binsPt[slice]) / 2.
                                              : hYieldPt->GetXaxis()->GetBinCenter(slice + 1); 
        xBinInit = hInitFracPi->GetXaxis()->FindBin(xBinCentre);
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
      
      Double_t gausParamsPi[nPar] = { 
        fractionPions,
        fractionKaons,
        fractionProtons,
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
        fractionKaons,
        fractionProtons,
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
        fractionKaons,
        fractionProtons,
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
        fractionKaons,
        fractionProtons,
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
        fractionErrorLowKaons,
        fractionErrorLowProtons,
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
        fractionErrorLowKaons,
        fractionErrorLowProtons,
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
        fractionErrorLowKaons,
        fractionErrorLowProtons,
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
        fractionErrorLowKaons,
        fractionErrorLowProtons,
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
        fractionErrorUpKaons,
        fractionErrorUpProtons,
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
        fractionErrorUpKaons,
        fractionErrorUpProtons,
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
        fractionErrorUpKaons,
        fractionErrorUpProtons,
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
        fractionErrorUpKaons,
        fractionErrorUpProtons,
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
        0.1,
        0.1,
        0.1,
        0.1,
        (takeIntoAccountMuons ? 0.1 : 0.),
        
        0.0,
        
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        enableShift ? shiftStepSize : 0.,
        (enableShift && takeIntoAccountMuons) ? shiftStepSize : 0.
      };
      
      
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
        0.1,
        0.1,
        0.1,
        0.1,
        (takeIntoAccountMuons ? 0.1 : 0.),
        
        0.0
      };
      
      if (regularisation <= 0 && iter == 1) {
        // In case of no regularisation, do the fit of the electron fraction here (compare comment below)
        if (mode == kPMpT && slice == electronFractionThresholdBinForFitting) 
          hFractionElectrons->Fit(fElectronFraction, "N", "", lowFittingBoundElectronFraction, electronFractionThresholdForFitting);
      }
      
      if ((regularisation > 0 && iter == 0) || (regularisation <= 0 && iter == 1)) {
        // Set the electron fraction to the negative pT -> A function will be used
        // to evaluate the electron fraction for each bin above the threshold
        if(mode == kPMpT && slice >= electronFractionThresholdBinForFitting) {
          // In case of no regularisation, mathFit has no information about the fraction of other x bins.
          // Thus, the electron fraction is evaluated and set here. For the case w/ regularisation,
          // just "-pT" is set and the electron fraction will be evaluated during the fit.
          Double_t fixElectronFraction = (regularisation <= 0) ? fElectronFraction->Eval((binsPt[slice + 1] + binsPt[slice]) / 2.)
                                                               : -(binsPt[slice + 1] + binsPt[slice]) / 2.;
          
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
        if(mode != kPMpT || slice >= muonFractionThresholdBinForFitting) {
          // "Abuse" the muon fraction to forward the pT, which can then be used to get some modified electron fraction
          const Double_t fixedValue = -(binsPt[slice] + binsPt[slice + 1]) / 2.;
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
      }
      
      
      if (iter == 1) {
        // The parameters are only used for fitMethod < 2. Thus, they can be set for these methods,
        // although a different method is used
        totalDeltaPion->SetParameters(gausParamsPi);
        totalDeltaElectron->SetParameters(gausParamsEl);
        totalDeltaKaon->SetParameters(gausParamsKa);
        totalDeltaProton->SetParameters(gausParamsPr);
      }
      
      const TString binInfo = (mode == kPMpT) ? Form("%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                              : Form("%.2f_%s_%.2f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                                     modeShortName[mode].Data(), hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
      
      const TString binInfoTitle = (mode == kPMpT) ? Form("%.2f < Pt <%.2f", binsPt[slice], binsPt[slice + 1])
                                                   : Form("%.2f < %s < %.2f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                                          modeShortName[mode].Data(), 
                                                          hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
      
      const TString fitFuncSuffix = (mode == kPMpT) ? Form("%.3f_Pt_%.3f", binsPt[slice], binsPt[slice + 1])
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
        
          // TODO At the moment, the covariance matrix is NOT forwarded from TMinuit (has completely different dimensions)
          // -> Since it is not used at the moment, this is not necessary. If it is to be used in future,
          // this has to be implemented! But one has to be careful, since parameters from different bins then
          // depend on each other! Furthermore, there will be no errors for fixed parameters like muon fraction or electron fraction
          // above the corresponding threshold in the covariance matrix, but an estimated error will be set manually.
          errFlag =  errFlag | doSimultaneousFitRegularised(nParSimultaneousFitRegularised, gausParamsSimultaneousFitRegularised, 
                                                            parameterErrorsOutRegularised, &covMatrix[0][0],
                                                            stepSizeSimultaneousFitRegularised, 
                                                            lowParLimitsSimultaneousFitRegularised, 
                                                            upParLimitsSimultaneousFitRegularised, reducedChiSquare);
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
                                      upParLimitsSimultaneousFit, reducedChiSquare);
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
        SetReasonableXaxisRange(hDeltaPi[slice], binLow, binHigh);
        hDeltaPi[slice]->Draw("e");
        
        fitFuncTotalDeltaPion[slice] = (TF1*)totalDeltaPion->Clone(Form("Fit_Total_DeltaPion_%s", fitFuncSuffix.Data()));
        
        hDeltaPiFitQA[slice] = (TH1D*)hDeltaPi[slice]->Clone(Form("hDeltaPiFitQA_%d", slice));
        hDeltaPiFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaPiFitQA[slice]->Add(fitFuncTotalDeltaPion[slice], -1);
        hDeltaPiFitQA[slice]->Divide(hDeltaPi[slice]);
        
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
        hDeltaPiFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
        hDeltaPiFitQA[slice]->Draw("e");    
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 3, reducedChiSquare);
        
       // TMatrixDSym covMatrixPi(nParUsed, &covMatrix[0][0]);   
       
        setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kPi, parametersOut, parameterErrorsOut, hFractionPions,
                              hFractionPionsDeltaPion, hFractionElectronsDeltaPion, hFractionKaonsDeltaPion,
                              hFractionProtonsDeltaPion, hFractionMuonsDeltaPion, hYieldPions, hYieldPionsDeltaPion, hYieldElectronsDeltaPion,
                              hYieldKaonsDeltaPion, hYieldProtonsDeltaPion, hYieldMuonsDeltaPion, normaliseResults);
        
        // Also set specific muon fractions and yields -> The deltaSpecies histos are not needed here: They will be set together with
        // the fraction and yields for all other species
        setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kMu, parametersOut, parameterErrorsOut, hFractionMuons,
                              0x0, 0x0, 0x0, 0x0, 0x0, hYieldMuons, 0x0, 0x0, 0x0, 0x0, 0x0, normaliseResults);
        
        
        // DeltaElectrons
        cSingleFit[slice][0]->cd(1);
        
        hDeltaEl[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaEl[slice], binLow, binHigh);
        hDeltaEl[slice]->Draw("e");
        
        fitFuncTotalDeltaElectron[slice] = (TF1*)totalDeltaElectron->Clone(Form("Fit_Total_DeltaElectron_%s", fitFuncSuffix.Data()));
        
        hDeltaElFitQA[slice] = (TH1D*)hDeltaEl[slice]->Clone(Form("hDeltaElFitQA_%d", slice));
        hDeltaElFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaElFitQA[slice]->Add(fitFuncTotalDeltaElectron[slice], -1);
        hDeltaElFitQA[slice]->Divide(hDeltaEl[slice]);
        
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
        hDeltaElFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
        hDeltaElFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 1, reducedChiSquare);
        
        //TMatrixDSym covMatrixEl(nParUsed, &covMatrix[0][0]);
        
        setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kEl, parametersOut, parameterErrorsOut, hFractionElectrons,
                              hFractionPionsDeltaElectron, hFractionElectronsDeltaElectron, hFractionKaonsDeltaElectron,
                              hFractionProtonsDeltaElectron, hFractionMuonsDeltaElectron, hYieldElectrons, hYieldPionsDeltaElectron, 
                              hYieldElectronsDeltaElectron, hYieldKaonsDeltaElectron, hYieldProtonsDeltaElectron, hYieldMuonsDeltaElectron, 
                              normaliseResults);
        
        
        
        // DeltaKaons 
        cSingleFit[slice][1]->cd(1);
        
        hDeltaKa[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaKa[slice], binLow, binHigh);
        hDeltaKa[slice]->Draw("e");
        
        fitFuncTotalDeltaKaon[slice] = (TF1*)totalDeltaKaon->Clone(Form("Fit_Total_DeltaKaon_%s", fitFuncSuffix.Data()));
        
        hDeltaKaFitQA[slice] = (TH1D*)hDeltaKa[slice]->Clone(Form("hDeltaKaFitQA_%d", slice));
        hDeltaKaFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaKaFitQA[slice]->Add(fitFuncTotalDeltaKaon[slice], -1);
        hDeltaKaFitQA[slice]->Divide(hDeltaKa[slice]);
        
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
        hDeltaKaFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
        hDeltaKaFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 2, reducedChiSquare);
        
        //TMatrixDSym covMatrixKa(nParUsed, &covMatrix[0][0]);
        
        setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kKa, parametersOut, parameterErrorsOut, hFractionKaons,
                              hFractionPionsDeltaKaon, hFractionElectronsDeltaKaon, hFractionKaonsDeltaKaon, hFractionProtonsDeltaKaon,
                              hFractionMuonsDeltaKaon, hYieldKaons, hYieldPionsDeltaKaon, hYieldElectronsDeltaKaon, hYieldKaonsDeltaKaon,
                              hYieldProtonsDeltaKaon, hYieldMuonsDeltaKaon, normaliseResults);
        
        
        
        // DeltaProtons
        cSingleFit[slice][3]->cd(1);
        
        hDeltaPr[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaPr[slice], binLow, binHigh);
        hDeltaPr[slice]->Draw("e");
        
        fitFuncTotalDeltaProton[slice] = (TF1*)totalDeltaProton->Clone(Form("Fit_Total_DeltaProton_%s", fitFuncSuffix.Data()));
        
        hDeltaPrFitQA[slice] = (TH1D*)hDeltaPr[slice]->Clone(Form("hDeltaPrFitQA_%d", slice));
        hDeltaPrFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaPrFitQA[slice]->Add(fitFuncTotalDeltaProton[slice], -1);
        hDeltaPrFitQA[slice]->Divide(hDeltaPr[slice]);
        
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
        hDeltaPrFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
        hDeltaPrFitQA[slice]->Draw("e");
        
        hReducedChiSquarePt->SetBinContent(slice + 1, 4, reducedChiSquare);
        
        //TMatrixDSym covMatrixPr(nParUsed, &covMatrix[0][0]);
        
        Double_t normalisationFactor = 1.0;
        normalisationFactor = setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kPr, parametersOut, parameterErrorsOut, 
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
                        stepSize, lowParLimitsPi, upParLimitsPi, totalDeltaPion, reducedChiSquare);
        
        hDeltaPi[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaPi[slice], binLow, binHigh);
        hDeltaPi[slice]->Draw("e");
        
        fitFuncTotalDeltaPion[slice] = (TF1*)totalDeltaPion->Clone(Form("Fit_Total_DeltaPion_%s", fitFuncSuffix.Data()));
        
        hDeltaPiFitQA[slice] = (TH1D*)hDeltaPi[slice]->Clone(Form("hDeltaPiFitQA_%d", slice));
        hDeltaPiFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaPiFitQA[slice]->Add(fitFuncTotalDeltaPion[slice], -1);
        hDeltaPiFitQA[slice]->Divide(hDeltaPi[slice]);
        
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
        hDeltaPiFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
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
          setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kPi, parametersOut, parameterErrorsOut, hFractionPions,
                                hFractionPionsDeltaPion, hFractionElectronsDeltaPion, hFractionKaonsDeltaPion,
                                hFractionProtonsDeltaPion, hFractionMuonsDeltaPion, hYieldPions, hYieldPionsDeltaPion, hYieldElectronsDeltaPion,
                                hYieldKaonsDeltaPion, hYieldProtonsDeltaPion, hYieldMuonsDeltaPion);
          
          // Also set specific muon fractions and yields -> The deltaSpecies histos are not needed here: They will be set together with
          // the fraction and yields for all other species
          setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kMu, parametersOut, parameterErrorsOut, hFractionMuons,
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
                        stepSize, lowParLimitsEl, upParLimitsEl, totalDeltaElectron, reducedChiSquare);
        
        hDeltaEl[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaEl[slice], binLow, binHigh);
        hDeltaEl[slice]->Draw("e");
        
        fitFuncTotalDeltaElectron[slice] = (TF1*)totalDeltaElectron->Clone(Form("Fit_Total_DeltaElectron_%s", fitFuncSuffix.Data()));
        
        hDeltaElFitQA[slice] = (TH1D*)hDeltaEl[slice]->Clone(Form("hDeltaElFitQA_%d", slice));
        hDeltaElFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaElFitQA[slice]->Add(fitFuncTotalDeltaElectron[slice], -1);
        hDeltaElFitQA[slice]->Divide(hDeltaEl[slice]);
        
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
        hDeltaElFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
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
          setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kEl, parametersOut, parameterErrorsOut, hFractionElectrons,
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
                        stepSize, lowParLimitsKa, upParLimitsKa, totalDeltaKaon, reducedChiSquare);
        
        hDeltaKa[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaKa[slice], binLow, binHigh);
        hDeltaKa[slice]->Draw("e");
        
        fitFuncTotalDeltaKaon[slice] = (TF1*)totalDeltaKaon->Clone(Form("Fit_Total_DeltaKaon_%s", fitFuncSuffix.Data()));
        
        hDeltaKaFitQA[slice] = (TH1D*)hDeltaKa[slice]->Clone(Form("hDeltaKaFitQA_%d", slice));
        hDeltaKaFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaKaFitQA[slice]->Add(fitFuncTotalDeltaKaon[slice], -1);
        hDeltaKaFitQA[slice]->Divide(hDeltaKa[slice]);
        
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
        hDeltaKaFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
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
          setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kKa, parametersOut, parameterErrorsOut, hFractionKaons,
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
                        stepSize, lowParLimitsPr, upParLimitsPr, totalDeltaProton, reducedChiSquare);
        
        hDeltaPr[slice]->SetTitle("");
        SetReasonableXaxisRange(hDeltaPr[slice], binLow, binHigh);
        hDeltaPr[slice]->Draw("e");
        
        fitFuncTotalDeltaProton[slice] = (TF1*)totalDeltaProton->Clone(Form("Fit_Total_DeltaProton_%s", fitFuncSuffix.Data()));
        
        hDeltaPrFitQA[slice] = (TH1D*)hDeltaPr[slice]->Clone(Form("hDeltaPrFitQA_%d", slice));
        hDeltaPrFitQA[slice]->GetYaxis()->SetTitle("(Data - Fit) / Data");
        hDeltaPrFitQA[slice]->Add(fitFuncTotalDeltaProton[slice], -1);
        hDeltaPrFitQA[slice]->Divide(hDeltaPr[slice]);
        
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
        hDeltaPrFitQA[slice]->GetYaxis()->SetRangeUser(fitQAaxisLowBound, fitQAaxisUpBound);
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
          setFractionsAndYields(slice, inverseBinWidth, binWidthFitHisto, kPr, parametersOut, parameterErrorsOut, hFractionProtons,
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
            
            sumOfParticles = inverseBinWidth * gausParamsPi[5] / binWidthFitHisto; // Divide by binWidthFitHisto, since gausParamsXX includes this width
            
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
            
            sumOfParticles = inverseBinWidth * gausParamsEl[5] / binWidthFitHisto; // Divide by binWidthFitHisto, since gausParamsXX includes this width
            
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
            
            sumOfParticles = inverseBinWidth * gausParamsKa[5] / binWidthFitHisto; // Divide by binWidthFitHisto, since gausParamsXX includes this width
            
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
            
            sumOfParticles = inverseBinWidth * gausParamsPr[5] / binWidthFitHisto; // Divide by binWidthFitHisto, since gausParamsXX includes this width
            
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
            
            sumOfParticles = inverseBinWidth * integralTotal / binWidthFitHisto; // Divide by binWidthFitHisto, since integralTotal includes this width
            
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
      if (slice % 18 == 0 || slice == pSliceLow) {
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
      
      TString saveDir = (mode == kPMpT) ? Form("SingleFit_%.2f_Pt_%.2f", binsPt[slice], binsPt[slice + 1])
                                        : Form("SingleFit_%.2f_%s_%.2f", hFractionPions->GetXaxis()->GetBinLowEdge(slice + 1), 
                                              modeShortName[mode].Data(), hFractionPions->GetXaxis()->GetBinUpEdge(slice + 1));
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
  
  // Calculate MC to-pi ratios -> In MC the yields are uncorrelated, so just divide the histos to get the correct result
  hRatioToPiElectronsMC->Divide(hYieldElectronsMC, hYieldPionsMC);
  hRatioToPiMuonsMC->Divide(hYieldMuonsMC, hYieldPionsMC);
  hRatioToPiKaonsMC->Divide(hYieldKaonsMC, hYieldPionsMC);
  hRatioToPiProtonsMC->Divide(hYieldProtonsMC, hYieldPionsMC);
  
  
  TCanvas* cFractions = new TCanvas("cFractions", "Particle fractions",100,10,1200,800);
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
    hFractionMuons->Draw("e p same");
  }
  if (plotIdentifiedSpectra) {
    SetReasonableAxisRange(hFractionMuonsMC->GetXaxis(), mode, pLow, pHigh);
    hFractionMuonsMC->Draw("e p same");
  }
  
  hFractionSummed->Draw("e p same");
  
  if (mode == kPMpT) {
    fElectronFraction->SetRange(lowFittingBoundElectronFraction, pHigh);
    fElectronFraction->Draw("same");
  }
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  if (plotIdentifiedSpectra)
    legend->SetNColumns(2);
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, "Fit", "");
  if (plotIdentifiedSpectra)
    legend->AddEntry((TObject*)0x0, identifiedLabels[isMC].Data(), "");
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
  if (takeIntoAccountMuons)
    legend->AddEntry(hFractionMuons, "#mu", "p");
  else
    legend->AddEntry((TObject*)0x0, "", "");
  if (plotIdentifiedSpectra)
    legend->AddEntry(hFractionMuonsMC, "#mu", "p");
  legend->AddEntry(hFractionSummed, "Total", "p");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(cFractions);
  

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
    
    hFractionComparisonTotal->SetBinContent(i, hFractionSummed->GetBinContent(i));
    hFractionComparisonTotal->SetBinError(i, hFractionSummed->GetBinError(i));
  }
  
  hFractionComparisonPions->Divide(hFractionPionsMC);
  hFractionComparisonElectrons->Divide(hFractionElectronsMC);
  if (takeIntoAccountMuons)
    hFractionComparisonMuons->Divide(hFractionMuonsMC);
  hFractionComparisonKaons->Divide(hFractionKaonsMC);
  hFractionComparisonProtons->Divide(hFractionProtonsMC);
  
  
  TCanvas* cFractionComparisons = new TCanvas("cFractionComparisons", "Particle fraction comparisons",100,10,1200,800);
  cFractionComparisons->SetGridx(1);
  cFractionComparisons->SetGridy(1);
  cFractionComparisons->SetLogx(mode == kPMpT);
  hFractionComparisonPions->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hFractionComparisonPions->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hFractionComparisonPions->GetXaxis()->SetNoExponent(kTRUE);
  hFractionComparisonPions->Draw("e p");
  
  hFractionComparisonElectrons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hFractionComparisonElectrons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonElectrons->Draw("e p same");
  
  if (takeIntoAccountMuons) {
    hFractionComparisonMuons->GetYaxis()->SetRangeUser(0.0, 10.0);
    SetReasonableAxisRange(hFractionComparisonMuons->GetXaxis(), mode, pLow, pHigh);
    hFractionComparisonMuons->Draw("e p same");
  }
  
  hFractionComparisonKaons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hFractionComparisonKaons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonKaons->Draw("e p same");
  
  hFractionComparisonProtons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hFractionComparisonProtons->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonProtons->Draw("e p same");
  
  hFractionComparisonTotal->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hFractionComparisonTotal->GetXaxis(), mode, pLow, pHigh);
  hFractionComparisonTotal->Draw("e p same");
  
  TLegend* legend2 = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetNColumns(2);
  legend2->AddEntry(hFractionComparisonPions, "#pi", "p");
  legend2->AddEntry(hFractionComparisonKaons, "K", "p");
  legend2->AddEntry(hFractionComparisonProtons, "p", "p");
  legend2->AddEntry(hFractionComparisonElectrons, "e", "p");
  if (takeIntoAccountMuons)
    legend2->AddEntry(hFractionComparisonMuons, "#mu", "p");
  legend2->AddEntry(hFractionComparisonTotal, "Total", "p");
  legend2->Draw();
  
  ClearTitleFromHistoInCanvas(cFractionComparisons);
  
  // Normalise the yields
  normaliseYieldHist(hYieldPions, numEvents, deta);
  normaliseYieldHist(hYieldPionsMC, numEvents, deta);
  normaliseYieldHist(hYieldPionsDeltaElectron, numEvents, deta);
  normaliseYieldHist(hYieldPionsDeltaPion, numEvents, deta);
  normaliseYieldHist(hYieldPionsDeltaKaon, numEvents, deta);
  normaliseYieldHist(hYieldPionsDeltaProton, numEvents, deta);
  
  normaliseYieldHist(hYieldElectrons, numEvents, deta);
  normaliseYieldHist(hYieldElectronsMC, numEvents, deta);
  normaliseYieldHist(hYieldElectronsDeltaElectron, numEvents, deta);
  normaliseYieldHist(hYieldElectronsDeltaPion, numEvents, deta);
  normaliseYieldHist(hYieldElectronsDeltaKaon, numEvents, deta);
  normaliseYieldHist(hYieldElectronsDeltaProton, numEvents, deta);
  
  normaliseYieldHist(hYieldMuons, numEvents, deta);
  normaliseYieldHist(hYieldMuonsMC, numEvents, deta);
  normaliseYieldHist(hYieldMuonsDeltaElectron, numEvents, deta);
  normaliseYieldHist(hYieldMuonsDeltaPion, numEvents, deta);
  normaliseYieldHist(hYieldMuonsDeltaKaon, numEvents, deta);
  normaliseYieldHist(hYieldMuonsDeltaProton, numEvents, deta);
  
  normaliseYieldHist(hYieldKaons, numEvents, deta);
  normaliseYieldHist(hYieldKaonsMC, numEvents, deta);
  normaliseYieldHist(hYieldKaonsDeltaElectron, numEvents, deta);
  normaliseYieldHist(hYieldKaonsDeltaPion, numEvents, deta);
  normaliseYieldHist(hYieldKaonsDeltaKaon, numEvents, deta);
  normaliseYieldHist(hYieldKaonsDeltaProton, numEvents, deta);
  
  normaliseYieldHist(hYieldProtons, numEvents, deta);
  normaliseYieldHist(hYieldProtonsMC, numEvents, deta);
  normaliseYieldHist(hYieldProtonsDeltaElectron, numEvents, deta);
  normaliseYieldHist(hYieldProtonsDeltaPion, numEvents, deta);
  normaliseYieldHist(hYieldProtonsDeltaKaon, numEvents, deta);
  normaliseYieldHist(hYieldProtonsDeltaProton, numEvents, deta);
  
  normaliseYieldHist(hYieldSummedMC, numEvents, deta);
  
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
      
      SetReasonableAxisRange(hMCgenYieldsPrimSpecies[i]->GetXaxis(), kPMpT, pLow, pHigh);
      normaliseGenYieldMCtruthHist(hMCgenYieldsPrimSpecies[i], numEvents, deta);
    }
  }
  
  
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
  
  hYieldComparisonPions->Divide(hYieldPionsMC);
  hYieldComparisonElectrons->Divide(hYieldElectronsMC);
  if (takeIntoAccountMuons)
    hYieldComparisonMuons->Divide(hYieldMuonsMC);
  hYieldComparisonKaons->Divide(hYieldKaonsMC);
  hYieldComparisonProtons->Divide(hYieldProtonsMC);
  
  
  TCanvas* cYieldComparisons = new TCanvas("cYieldComparisons", "Particle yield comparisons",100,10,1200,800);
  cYieldComparisons->SetGridx(1);
  cYieldComparisons->SetGridy(1);
  cYieldComparisons->SetLogx(mode == kPMpT);
  hYieldComparisonPions->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hYieldComparisonPions->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonPions->GetXaxis()->SetMoreLogLabels(kTRUE);
  hYieldComparisonPions->GetXaxis()->SetNoExponent(kTRUE);
  hYieldComparisonPions->Draw("e p");
  
  hYieldComparisonElectrons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hYieldComparisonElectrons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonElectrons->Draw("e p same");
  
  if (takeIntoAccountMuons) {
    hYieldComparisonMuons->GetYaxis()->SetRangeUser(0.0, 10.0);
    SetReasonableAxisRange(hYieldComparisonMuons->GetXaxis(), mode, pLow, pHigh);
    hYieldComparisonMuons->Draw("e p same");
  }
  
  hYieldComparisonKaons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hYieldComparisonKaons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonKaons->Draw("e p same");
  
  hYieldComparisonProtons->GetYaxis()->SetRangeUser(0.0, 10.0);
  SetReasonableAxisRange(hYieldComparisonProtons->GetXaxis(), mode, pLow, pHigh);
  hYieldComparisonProtons->Draw("e p same");
  
  TLegend* legend3 = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetNColumns(2);
  legend3->AddEntry(hYieldComparisonPions, "#pi", "p");
  legend3->AddEntry(hYieldComparisonKaons, "K", "p");
  legend3->AddEntry(hYieldComparisonProtons, "p", "p");
  legend3->AddEntry(hYieldComparisonElectrons, "e", "p");
  if (takeIntoAccountMuons)
    legend3->AddEntry(hYieldComparisonMuons, "#mu", "p");
  legend3->Draw();
  
  ClearTitleFromHistoInCanvas(cYieldComparisons);
  
  
  
  
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
  if (takeIntoAccountMuons)
    legendYields->AddEntry(hYieldMuons, "#mu", "p");
  else
    legendYields->AddEntry((TObject*)0x0, "", "");
  if (plotIdentifiedSpectra)
    legendYields->AddEntry(hYieldMuonsMC, "#mu", "p");
  legendYields->Draw();
  
  ClearTitleFromHistoInCanvas(cYields);
  
  
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
  
  delete gFractionElectronsData;
  delete fElectronFraction;
  
  delete mathFit;
  
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
  
  saveF->Close();

  return 0; 
}
