/* $Id$ */

#ifndef ALIUNFOLDING_H
#define ALIUNFOLDING_H

//
// class that implements several unfolding methods
// I.e. chi2 minimization and bayesian unfolding
// The whole class is static and not thread-safe (due to the fact that MINUIT unfolding is not thread-safe)
//

// TMatrixD, TVectorD defined here, because it does not seem possible to predeclare these (or i do not know how)
// -->
// $ROOTSYS/include/TVectorDfwd.h:21: conflicting types for `typedef struct TVectorT<Double_t> TVectorD'
// PWG0/AliUnfolding.h:21: previous declaration as `struct TVectorD'

#include "TObject.h"
#include <TMatrixD.h>
#include <TVectorD.h>

class TH1;
class TH2;
class TF1;
class TCanvas;
class TVirtualPad;
class TAxis;

class AliUnfolding : public TObject
{
  public:
  enum RegularizationType { kNone = 0, kPol0, kPol1, kLog, kEntropy, kCurvature, kRatio, kPowerLaw, kLogLog };
    enum MethodType { kInvalid = -1, kChi2Minimization = 0, kBayesian = 1, kFunction = 2};

    virtual ~AliUnfolding() {};

    static void SetUnfoldingMethod(MethodType methodType);
    static void SetCreateOverflowBin(Float_t overflowBinLimit);
    static void SetSkipBinsBegin(Int_t nBins);
    static void SetNbins(Int_t nMeasured, Int_t nUnfolded);
    static void SetChi2Regularization(RegularizationType type, Float_t weight);
    static void SetMinuitStepSize(Float_t stepSize) { fgMinuitStepSize = stepSize; }
    static void SetMinuitPrecision(Float_t pres) {fgMinuitPrecision = pres;}
    static void SetMinuitMaxIterations(Int_t iter) {fgMinuitMaxIterations = iter;}
    static void SetMinimumInitialValue(Bool_t flag, Float_t value = -1) { fgMinimumInitialValue = flag; fgMinimumInitialValueFix = value; }
    static void SetNormalizeInput(Bool_t flag) { fgNormalizeInput = flag; }
    static void SetNotFoundEvents(Float_t notFoundEvents) { fgNotFoundEvents = notFoundEvents; }
    static void SetSkip0BinInChi2(Bool_t flag) { fgSkipBin0InChi2 = flag; }
    static void SetBayesianParameters(Float_t smoothing, Int_t nIterations);
    static void SetFunction(TF1* function);
    static void SetDebug(Bool_t flag) { fgDebug = flag; }
    static void SetPowern(Int_t n) {fgPowern = -1*n;}

    static Int_t Unfold(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check = kFALSE);
    static Int_t UnfoldGetBias(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result);

    static TH1* GetPenaltyPlot(Double_t* params);
    static TH1* GetPenaltyPlot(TH1* corrected);

    static TH1* GetResidualsPlot(Double_t* params);
    static TH1* GetResidualsPlot(TH1* corrected);

    static Double_t GetChi2FromFit() {return fChi2FromFit;}
    static Double_t GetPenaltyVal()  {return fPenaltyVal;}
    static Double_t GetAvgResidual() {return fAvgResidual;}

    static void DrawResults(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TCanvas *canvas = 0, Int_t reuseHists = kFALSE,TH1 *unfolded=0);
    static void InteractiveUnfold(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions);
    static void RedrawInteractive();  

  protected:
    AliUnfolding() {};

    static Int_t UnfoldWithMinuit(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check);
    static Int_t UnfoldWithBayesian(TH2* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* aResult);
    static Int_t UnfoldWithFunction(TH2* correlation, TH1* efficiency, TH1* measured, TH1* /* initialConditions */, TH1* aResult);

    static void CreateOverflowBin(TH2* correlation, TH1* measured); 
    static void SetStaticVariables(TH2* correlation, TH1* measured, TH1* efficiency);

    static void MakePads();
    static void DrawGuess(Double_t *params, TVirtualPad *pfolded=0, TVirtualPad *pres=0, TVirtualPad *ppen=0, Int_t reuseHists = kFALSE, TH1* unfolded=0);

    static Double_t RegularizationPol0(TVectorD& params);
    static Double_t RegularizationPol1(TVectorD& params);
    static Double_t RegularizationTotalCurvature(TVectorD& params);
    static Double_t RegularizationEntropy(TVectorD& params);
    static Double_t RegularizationLog(TVectorD& params);
    static Double_t RegularizationRatio(TVectorD& params);
    static Double_t RegularizationPowerLaw(TVectorD& params);
    static Double_t RegularizationLogLog(TVectorD& params);

    static void Chi2Function(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t);
    static void TF1Function(Int_t& unused1, Double_t* unused2, Double_t& chi2, Double_t *params, Int_t unused3);

    // static variable to be accessed by MINUIT
    static TMatrixD* fgCorrelationMatrix;            // contains fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationMatrixSquared;     // contains squared fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationCovarianceMatrix;  // contains the errors of fCurrentESD
    static TVectorD* fgCurrentESDVector;             // contains fCurrentESD
    static TVectorD* fgEntropyAPriori;               // a-priori distribution for entropy regularization
    static TVectorD* fgEfficiency;                   // efficiency
    /*
    static TVectorD* fgBinWidths;                    // bin widths to be taken into account in regularization
    static TVectorD* fgBinPos;                       // bin positions of unfolded
    */
    static TAxis *fgUnfoldedAxis;                    // bin widths and positions for unfolded
    static TAxis *fgMeasuredAxis;                    // bin widths and positions for measured

    static TF1* fgFitFunction;                       // fit function

    // --- configuration params follow ---
    static MethodType fgMethodType;                  // unfolding method to be used
    static Int_t fgMaxParams;                        // bins in unfolded histogram = number of fit params
    static Int_t fgMaxInput;                         // bins in measured histogram
    static Float_t fgOverflowBinLimit;               // to fix fluctuations at high multiplicities, all entries above the limit are summarized in one bin

    static RegularizationType fgRegularizationType;  // regularization that is used during Chi2 method
    static Float_t fgRegularizationWeight;           // factor for regularization term
    static Int_t fgSkipBinsBegin;                    // (optional) skip the given number of bins in the regularization
    static Float_t fgMinuitStepSize;                 // (usually not needed to be changed) step size in minimization
    static Float_t fgMinuitPrecision;                // precision used by minuit. default = 1e-6
    static Int_t   fgMinuitMaxIterations;            // maximum number of iterations used by minuit. default = 5000
    static Bool_t fgMinimumInitialValue;             // set all initial values at least to the smallest value among the initial values
    static Float_t fgMinimumInitialValueFix;         // use this as the minimum initial value instead of determining it automatically
    static Bool_t fgNormalizeInput;                  // normalize input spectrum
    static Float_t fgNotFoundEvents;                 // constraint on the total number of not found events sum(guess * (1/eff -1))
    static Bool_t fgSkipBin0InChi2;                  // skip bin 0 (= 0 measured) in chi2 function

    static Float_t fgBayesianSmoothing;              // smoothing parameter (0 = no smoothing)
    static Int_t   fgBayesianIterations;             // number of iterations in Bayesian method

    static Bool_t fgDebug;                           // debug flag
    // --- end of configuration ---
    
    static Int_t fgCallCount;                        // call count to chi2 function

    static Int_t fgPowern;                           // power of power law for regularization with power law

    static Double_t fChi2FromFit;                    // Chi2 from fit at current iteration
    static Double_t fPenaltyVal;                     // Penalty value at current iteration (\beta * PU)
    static Double_t fAvgResidual;                    // Sum residuals / nbins

    static Int_t fgPrintChi2Details;                 // debug for chi2 calc

    // Pointers for interactive unfolder
    static TCanvas *fgCanvas;                        // Canvas for interactive unfolder
    static TH1 *fghUnfolded;                         // Unfolding result for interactive unfolder
    static TH2 *fghCorrelation;                      // Response matrix for interactive unfolder
    static TH1 *fghEfficiency;                       // Efficiency histo for interactive unfolder
    static TH1 *fghMeasured;                         // Measured distribution for interactive unfolder

private:
    AliUnfolding(const AliUnfolding&);
    AliUnfolding& operator=(const AliUnfolding&);

  ClassDef(AliUnfolding, 0);
};

#endif

