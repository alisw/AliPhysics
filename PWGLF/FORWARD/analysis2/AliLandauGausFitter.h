#ifndef ALILANDAUGAUSFITTER_H
#define ALILANDAUGAUSFITTER_H
/**
 * @file   AliLandauGausFitter.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Mar 11 08:56:03 2014
 * 
 * @brief Declaration and implementation of fitter of Landau-Gauss
 * distributions to energy loss spectra.
 * 
 * @ingroup pwglf_forward 
 */
#include <AliLandauGaus.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>
#include <TArray.h>
#include <TFitResult.h>
#include <TError.h>

/**
 * Fit Landau-Gauss distributions to energy loss distributions 
 * 
 * @see AliLandauGaus 
 *
 * @ingroup pwglf_forward 
 */
class AliLandauGausFitter 
{
public:
  /** 
   * Enumeration of parameters 
   */
  enum { 
    /** Index of pre-constant @f$ C@f$ */
    kC		= AliLandauGaus::kC,
    /** Index of most probable value @f$ \Delta_p@f$ */
    kDelta	= AliLandauGaus::kDelta, 
    /** Index of Landau width @f$ \xi@f$ */
    kXi		= AliLandauGaus::kXi, 
    /** Index of Gaussian width @f$ \sigma@f$ */
    kSigma	= AliLandauGaus::kSigma, 
    /** Index of Gaussian additional width @f$ \sigma_n@f$ */
    kSigmaN	= AliLandauGaus::kSigmaN,
    /** Index of Number of particles @f$ N@f$ */
    kN		= AliLandauGaus::kN, 
    /** Base index of particle strengths @f$ a_i@f$ for 
	@f$i=2,\ldots,N@f$ */
    kA		= AliLandauGaus::kA
  };
  /** 
   * Get the fit options to use 
   * 
   * @return Fit options used through-out
   */
  static const char* GetFitOptions() { return "RNS"; }
  /** 
   * Constructor 
   * 
   * @param lowCut     Lower cut of spectrum - data below this cuts is ignored
   * @param maxRange   Maximum range to fit to 
   * @param minusBins  The number of bins below maximum to use 
   */
  AliLandauGausFitter(Double_t lowCut, Double_t maxRange, UShort_t minusBins)
    : fLowCut(lowCut), fMaxRange(maxRange), fMinusBins(minusBins), 
      fFitResults(0), fFunctions(0), fDebug(false) 
  {
    fFitResults.SetOwner();
    fFunctions.SetOwner();
  }
  /** 
   * Destructor
   * 
   */
  virtual ~AliLandauGausFitter()
  {
    fFitResults.Delete();
    fFunctions.Delete();
  }
  /** 
   * Enable/disable debugging output. 
   * 
   * @param debug If true, enable debugging output
   */
  void SetDebug(Bool_t debug=true) { fDebug = debug; }
  /** 
   * Clear internal arrays
   * 
   */
  void Clear() 
  {
    fFitResults.Clear();
    fFunctions.Clear();
  }
  /** 
   * Fit a 1-particle signal to the passed energy loss distribution 
   * 
   * Note that this function clears the internal arrays first 
   *
   * @param dist    Data to fit the function to 
   * @param sigman If larger than zero, the initial guess of the
   *               detector induced noise. If zero or less, then this 
   *               parameter is ignored in the fit (fixed at 0)
   * 
   * @return The function fitted to the data 
   */
  TF1* Fit1Particle(TH1* dist, Double_t sigman=-1);
  /** 
   * Fit a N-particle signal to the passed energy loss distribution 
   *
   * If there's no 1-particle fit present, it does that first 
   *
   * @param dist   Data to fit the function to 
   * @param n      Number of particle signals to fit 
   * @param sigman If larger than zero, the initial guess of the
   *               detector induced noise. If zero or less, then this 
   *               parameter is ignored in the fit (fixed at 0)
   * 
   * @return The function fitted to the data 
   */
  TF1* FitNParticle(TH1* dist, UShort_t n, Double_t sigman=-1);
  /** 
   * Fit a composite distribution of energy loss from both primaries
   * and secondaries
   * 
   * @param dist   Distribution 
   * @param sigman If larger than zero, the initial guess of the
   *                detector included noise.  If zero or less this
   *                parameter is fixed to 0.
   * 
   * @return Function fitted to the data 
   */
  TF1* FitComposite(TH1* dist, Double_t sigman);
  /**
   * Get Lower cut on data 
   *
   * @return Lower cut on data 
   */
  Double_t GetLowCut() const { return fLowCut; }
  /**
   * Get Maximum range to fit 
   *
   * @return Maximum range to fit 
   */
  Double_t GetMaxRange() const { return fMaxRange; }
  /**
   * Get Number of bins from maximum to fit 1st peak
   *
   * @return Number of bins from maximum to fit 1st peak
   */
  UShort_t GetMinusBins() const { return fMinusBins; }
  /**
   * Get Array of fit results 
   *
   * @return Array of fit results 
   */
  const TObjArray& GetFitResults() const { return fFitResults; }
  /** 
   * Get Array of fit results  
   *
   * @return Array of fit results 
   */
  TObjArray& GetFitResults() { return fFitResults; }
  /**
   * Get Array of functions 
   *
   * @return Array of functions 
   */
  const TObjArray& GetFunctions() const { return fFunctions; }
  /** 
   * Get Array of functions  
   *
   * @return Array of functions 
   */
  TObjArray& GetFunctions() { return fFunctions; }  
protected:
  /** 
   * Set parameter limits on the function @a f.  The limits are only
   * set if @f$ low \le test \le high@f$
   * 
   * @param f     Function object 
   * @param iPar  Parameter number 
   * @param test  Initial guess
   * @param low   Low limit 
   * @param high  high limit 
   */
  void SetParLimits(TF1* f, Int_t iPar, Double_t test, 
		    Double_t low, Double_t high) 
  {
    if (test < low || test > high) {
      ::Warning("","Initial value of %s=%f not in [%f,%f]", 
		f->GetParName(iPar), test, low, high);
      return;
    }
    if (fDebug) 
      ::Info(/*"SetParLimits"*/"", "Set par limits on %-12s=%f: [%f,%f]",
	     f->GetParName(iPar), test, low, high);
    f->SetParLimits(iPar, low, high);
  }
  const Double_t fLowCut;     // Lower cut on data 
  const Double_t fMaxRange;   // Maximum range to fit 
  const UShort_t fMinusBins;  // Number of bins from maximum to fit 1st peak
  TObjArray fFitResults;      // Array of fit results 
  TObjArray fFunctions;       // Array of functions 
  Bool_t    fDebug;           // Debug flag
};


//____________________________________________________________________
inline TF1*
AliLandauGausFitter::Fit1Particle(TH1* dist, Double_t sigman)
{
  // Clear the cache 
  Clear();

  // Find the fit range 
  Int_t    cutBin  = TMath::Max(dist->GetXaxis()->FindBin(fLowCut),3);
  Int_t    maxBin  = TMath::Min(dist->GetXaxis()->FindBin(fMaxRange),
				dist->GetNbinsX());
  dist->GetXaxis()->SetRange(cutBin, maxBin);
  // dist->GetXaxis()->SetRangeUser(fLowCut, fMaxRange);
  
  // Get the bin with maximum 
  Int_t    peakBin = dist->GetMaximumBin();
  Double_t peakE   = dist->GetBinLowEdge(peakBin);
  Double_t rmsE    = dist->GetRMS();
  
  // Get the low edge 
  // dist->GetXaxis()->SetRangeUser(fLowCut, peakE);
  Int_t    minBin = peakBin - fMinusBins; // dist->GetMinimumBin();
  Double_t minE   = TMath::Max(dist->GetBinCenter(minBin),fLowCut);
  Double_t maxE   = dist->GetBinCenter(peakBin+2*fMinusBins);

  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxE);
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("Fit1Particle", 
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxE, minEb, maxEb, intg);
    return 0;
  }
    
  // Restore the range 
  dist->GetXaxis()->SetRange(1, maxBin);
  
  // Define the function to fit 
  TF1* f = AliLandauGaus::MakeF1(intg,peakE,peakE/10,peakE/5,sigman,minE,maxE);
  SetParLimits(f, kDelta, peakE,   minE, fMaxRange);
  SetParLimits(f, kXi,    peakE,   0,    2*rmsE); // 0.1
  SetParLimits(f, kSigma, peakE/5, 1e-5, rmsE); // 0.1
  if (sigman <= 0)  
    f->FixParameter(kSigmaN, 0);
  else 
    SetParLimits(f, kSigmaN, peakE, 0, rmsE);
  

  TString opts(Form("%s%s", GetFitOptions(), fDebug ? "" : "Q"));
  // Do the fit, getting the result object 
  if (fDebug) 
    ::Info(/*"Fit1Particle"*/"", "Fitting in the range %f,%f", minE, maxE);
  TFitResultPtr r = dist->Fit(f, opts, "", minE, maxE);
  if (!r.Get()) { 
    ::Warning("Fit1Particle", 
	      "No fit returned when processing %s in the range [%f,%f] "
	      "options %s", dist->GetName(), minE, maxE, opts.Data());
    return 0;
  }
  // f->SetRange(minE, fMaxRange);
  fFitResults.AddAtAndExpand(new TFitResult(*r), 0);
  fFunctions.AddAtAndExpand(f, 0);

  return f;
}
//____________________________________________________________________
inline TF1*
AliLandauGausFitter::FitNParticle(TH1* dist, UShort_t n, Double_t sigman)
{
  // Get the seed fit result 
  TFitResult* r = static_cast<TFitResult*>(fFitResults.At(0));
  TF1*        f = static_cast<TF1*>(fFunctions.At(0));
  if (!r || !f) { 
    f = Fit1Particle(dist, sigman);
    r = static_cast<TFitResult*>(fFitResults.At(0));
    if (!r || !f) { 
      ::Warning("FitNLandau", "No first shot at landau fit");
      return 0;
    }
  }

  // Get some parameters from seed fit 
  Double_t delta1  = r->Parameter(kDelta);
  Double_t xi1     = r->Parameter(kXi);
  Double_t maxEi   = n * (delta1 + xi1 * TMath::Log(n)) + 2 * n * xi1;
  Double_t minE    = f->GetXmin();

  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxEi);
  Double_t rmsE  = dist->GetRMS();
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("FitNParticle",
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxEi, minEb, maxEb, intg);
    return 0;
  }

  // Array of weights 
  TArrayD a(n-1);
  for (UShort_t i = 2; i <= n; i++) 
    a.fArray[i-2] = (n == 2 ? 0.05 : 0.000001);
  // Make the fit function 
  f = AliLandauGaus::MakeFn(r->Parameter(kC),
			    r->Parameter(kDelta),
			    r->Parameter(kXi),
			    r->Parameter(kSigma),
			    r->Parameter(kSigmaN),
			    n, a.fArray, minE, maxEi);
  SetParLimits(f,kDelta, f->GetParameter(kDelta), minE, fMaxRange);
  SetParLimits(f,kXi,    f->GetParameter(kXi),    0,    2*rmsE); //0.1
  SetParLimits(f,kSigma, f->GetParameter(kSigma), 1e-5, 2*rmsE); //0.1
  if (sigman <= 0)  f->FixParameter(kSigmaN, 0);
  else 
    SetParLimits(f,kSigmaN, f->GetParameter(kSigmaN), 0, rmsE);

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    SetParLimits(f,kA+i-2, a[i-2], 0, 1);
  }

  // Do the fit 
  TString opts(Form("%s%s", GetFitOptions(), fDebug ? "" : "Q"));
  if (fDebug) 
    ::Info(/*"FitNParticle"*/"", 
	   "Fitting in the range %f,%f (%d)", minE, maxEi, n);
  TFitResultPtr tr = dist->Fit(f, opts, "", minE, maxEi);
  
  // f->SetRange(minE, fMaxRange);
  fFitResults.AddAtAndExpand(new TFitResult(*tr), n-1);
  fFunctions.AddAtAndExpand(f, n-1);
  
  return f;
}  
//____________________________________________________________________
inline TF1*
AliLandauGausFitter::FitComposite(TH1* dist, Double_t sigman)
{
  // Find the fit range 
  Int_t    cutBin  = TMath::Max(dist->GetXaxis()->FindBin(fLowCut),3);
  Int_t    maxBin  = TMath::Min(dist->GetXaxis()->FindBin(fMaxRange),
				dist->GetNbinsX());
  dist->GetXaxis()->SetRange(cutBin, maxBin);
  
  // Get the bin with maximum 
  Int_t    peakBin = dist->GetMaximumBin();
  Double_t peakE   = dist->GetBinLowEdge(peakBin);
  
  // Get the low edge 
  // dist->GetXaxis()->SetRangeUser(fLowCut, peakE);
  Int_t    minBin = peakBin - fMinusBins; // dist->GetMinimumBin();
  Double_t minE   = TMath::Max(dist->GetBinCenter(minBin),fLowCut);
  Double_t maxE   = dist->GetBinCenter(peakBin+2*fMinusBins);

  // Get the range in bins and the integral of that range 
  Int_t    minEb = dist->GetXaxis()->FindBin(minE);
  Int_t    maxEb = dist->GetXaxis()->FindBin(maxE);
  Double_t intg  = dist->Integral(minEb, maxEb);
  if (intg <= 0) {
    ::Warning("Fit1Particle", 
	      "Integral of %s between [%f,%f] [%03d,%03d] = %f < 0", 
	      dist->GetName(), minE, maxE, minEb, maxEb, intg);
    return 0;
  }
    
  // Restore the range 
  dist->GetXaxis()->SetRange(1, maxBin);
  
  // Define the function to fit 
  TF1* seed = AliLandauGaus::MakeF1(1,peakE,peakE/10,peakE/5,sigman,minE,maxE);

  // Set initial guesses, parameter names, and limits  
  seed->SetParLimits(kDelta, minE, fMaxRange);
  seed->SetParLimits(kXi,    0.00, 0.1);
  seed->SetParLimits(kSigma, 1e-5, 0.1);
  if (sigman <= 0)  seed->FixParameter(kSigmaN, 0);
  else              seed->SetParLimits(kSigmaN, 0, fMaxRange);

  // Do the fit, getting the result object 
  if (fDebug) 
    ::Info(/*"FitComposite"*/"", "Fitting seed in the range %f,%f", minE, maxE);
  /* TFitResultPtr r = */ dist->Fit(seed, GetFitOptions(), "", minE, maxE);

  maxE = dist->GetXaxis()->GetXmax();
  TF1* comp = 
    AliLandauGaus::MakeComposite(0.8 * seed->GetParameter(kC),
				 seed->GetParameter(kDelta),
				 seed->GetParameter(kDelta)/10,
				 seed->GetParameter(kDelta)/5, 
				 1.20 * seed->GetParameter(kC),
				 seed->GetParameter(kXi),
				 minE, maxE);
  // comp->SetParLimits(kC,       minE, fMaxRange); // C
  comp->SetParLimits(kDelta,      minE, fMaxRange); // Delta
  comp->SetParLimits(kXi,         0.00, fMaxRange); // Xi 
  comp->SetParLimits(kSigma,      1e-5, fMaxRange); // Sigma
  // comp->SetParLimits(kSigma+1, minE, fMaxRange); // C
  comp->SetParLimits(kSigma+2,    0.00, fMaxRange); // Xi'
  comp->SetLineColor(kRed+1);
  comp->SetLineWidth(3);
  
  // Do the fit, getting the result object 
  TString opts(Form("%s%s", GetFitOptions(), fDebug ? "" : "Q"));
  if (fDebug) 
    ::Info("FitComposite", "Fitting composite in the range %f,%f", minE, maxE);
  /* TFitResultPtr r = */ dist->Fit(comp, opts, "", minE, maxE);

#if 0
  // This is to store the two components with the output
  TF1* part1 = static_cast<TF1*>(seed->Clone("part1"));
  part1->SetLineColor(kGreen+1);
  part1->SetLineWidth(4);
  part1->SetRange(minE, maxE);
  part1->SetParameters(comp->GetParameter(0), // C 
		       comp->GetParameter(1), // Delta
		       comp->GetParameter(2), // Xi
		       comp->GetParameter(3), // sigma
		       0);
  part1->Save(minE,maxE,0,0,0,0);
  dist->GetListOfFunctions()->Add(part1);

  TF1* part2 = static_cast<TF1*>(seed->Clone("part2"));
  part2->SetLineColor(kBlue+1);
  part2->SetLineWidth(4);
  part2->SetRange(minE, maxE);
  part2->SetParameters(comp->GetParameter(4), // C 
		       comp->GetParameter(1), // Delta
		       comp->GetParameter(5), // Xi
		       comp->GetParameter(3), // sigma
		       0);
  part2->Save(minE,maxE,0,0,0,0);
  dist->GetListOfFunctions()->Add(part2);
#endif
  return comp;
}

#endif
// Local Variables: 
//  mode: C++
// End:
