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

class AliUnfolding : public TObject
{
  public:
    enum RegularizationType { kNone = 0, kPol0, kPol1, kLog, kEntropy, kCurvature, kRatio };
    enum MethodType { kInvalid = -1, kChi2Minimization = 0, kBayesian = 1, kFunction = 2};

    virtual ~AliUnfolding() {};

    static void SetUnfoldingMethod(MethodType methodType);
    static void SetCreateOverflowBin(Float_t overflowBinLimit);
    static void SetSkipBinsBegin(Int_t nBins);
    static void SetNbins(Int_t nMeasured, Int_t nUnfolded);
    static void SetChi2Regularization(RegularizationType type, Float_t weight);
    static void SetMinuitStepSize(Float_t stepSize) { fgMinuitStepSize = stepSize; }
    static void SetMinimumInitialValue(Bool_t flag, Float_t value = -1) { fgMinimumInitialValue = flag; fgMinimumInitialValueFix = value; }
    static void SetNormalizeInput(Bool_t flag) { fgNormalizeInput = flag; }
    static void SetNotFoundEvents(Float_t notFoundEvents) { fgNotFoundEvents = notFoundEvents; }
    static void SetSkip0BinInChi2(Bool_t flag) { fgSkipBin0InChi2 = flag; }
    static void SetBayesianParameters(Float_t smoothing, Int_t nIterations);
    static void SetFunction(TF1* function);
    static void SetDebug(Bool_t flag) { fgDebug = flag; }

    static Int_t Unfold(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check = kFALSE);
    static Int_t UnfoldGetBias(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result);

    static TH1* GetPenaltyPlot(Double_t* params);
    static TH1* GetPenaltyPlot(TH1* corrected);
  
  protected:
    AliUnfolding() {};

    static Int_t UnfoldWithMinuit(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check);
    static Int_t UnfoldWithBayesian(TH2* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* aResult);
    static Int_t UnfoldWithFunction(TH2* correlation, TH1* efficiency, TH1* measured, TH1* /* initialConditions */, TH1* aResult);

    static void DrawGuess(Double_t *params);
    static void CreateOverflowBin(TH2* correlation, TH1* measured); 
    static void SetStaticVariables(TH2* correlation, TH1* measured, TH1* efficiency);

    static Double_t RegularizationPol0(TVectorD& params);
    static Double_t RegularizationPol1(TVectorD& params);
    static Double_t RegularizationTotalCurvature(TVectorD& params);
    static Double_t RegularizationEntropy(TVectorD& params);
    static Double_t RegularizationLog(TVectorD& params);
    static Double_t RegularizationRatio(TVectorD& params);

    static void Chi2Function(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t);
    static void TF1Function(Int_t& unused1, Double_t* unused2, Double_t& chi2, Double_t *params, Int_t unused3);

    // static variable to be accessed by MINUIT
    static TMatrixD* fgCorrelationMatrix;            // contains fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationMatrixSquared;     // contains squared fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationCovarianceMatrix;  // contains the errors of fCurrentESD
    static TVectorD* fgCurrentESDVector;             // contains fCurrentESD
    static TVectorD* fgEntropyAPriori;               // a-priori distribution for entropy regularization
    static TVectorD* fgEfficiency;                   // efficiency
    static TVectorD* fgBinWidths;                    // bin widths to be taken into account in regularization

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

private:
    AliUnfolding(const AliUnfolding&);
    AliUnfolding& operator=(const AliUnfolding&);

  ClassDef(AliUnfolding, 0);
};

#endif

