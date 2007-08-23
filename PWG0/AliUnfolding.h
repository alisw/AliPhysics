/* $Id$ */

#ifndef ALIUNFOLDING_H
#define ALIUNFOLDING_H

//
// class that implements several unfolding methods
// E.g. chi2 minimization and bayesian unfolding
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
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TCollection;

class AliUnfolding : public TObject
{
  public:
    enum RegularizationType { kNone = 0, kPol0, kPol1, kLog, kEntropy, kCurvature };
    enum MethodType { kChi2Minimization = 0, kBayesian = 1 };

    AliUnfolding();
    virtual ~AliUnfolding();

    void SetInput(TH2* correlationMatrix, TH1* efficiency, TH1* measured) { fCurrentCorrelation = correlationMatrix; fCurrentEfficiency = efficiency; fCurrentESD = measured; }
    void SetInitialConditions(TH1* initialConditions) { fInitialConditions = initialConditions; }
    const TH1* GetResult() const { return fResult; }

    static void SetParameters(Int_t measuredBins, Int_t unfoldedBins, Bool_t bigbin) { fMaxInput = measuredBins; fMaxParams = unfoldedBins; fgCreateBigBin = bigbin; }
    static void SetChi2MinimizationParameters(RegularizationType type, Float_t weight) { fgRegularizationType = type; fgRegularizationWeight = weight; }
    static void SetRegularizationRange(Int_t start, Int_t end) { fgRegularizationRangeStart = start; fgRegularizationRangeEnd = end; }
    static void SetBayesianParameters(Float_t smoothing, Int_t nIterations) { fgBayesianSmoothing = smoothing; fgBayesianIterations = nIterations; }

    Int_t ApplyMinuitFit(Bool_t check = kFALSE);
    Int_t ApplyBayesianMethod(Bool_t determineError = kTRUE);
    Int_t ApplyNBDFit();
    Int_t ApplyLaszloMethod();

    TH1* StatisticalUncertainty(MethodType methodType, Bool_t randomizeMeasured, Bool_t randomizeResponse, TH1* compareTo = 0);

  protected:
    static Double_t RegularizationPol0(TVectorD& params);
    static Double_t RegularizationPol1(TVectorD& params);
    static Double_t RegularizationTotalCurvature(TVectorD& params);
    static Double_t RegularizationEntropy(TVectorD& params);
    static Double_t RegularizationLog(TVectorD& params);

    static void MinuitFitFunction(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t);
    static void MinuitNBD(Int_t& unused1, Double_t* unused2, Double_t& chi2, Double_t *params, Int_t unused3);

    void SetupCurrentHists();

    Int_t UnfoldWithBayesian(* aEfficiency, TH1* measured, TH1* initialConditions, TH1* aResult, Float_t regPar, Int_t nIterations);
    Int_t UnfoldWithMinuit(TH1* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check);

    Float_t BayesCovarianceDerivate(Float_t matrixM[251][251], TH2* hResponse, Int_t k, Int_t i, Int_t r, Int_t u);

    TH1* fCurrentESD;         //! measured spectrum
    TH2* fCurrentCorrelation; //! correlation matrix
    TH1* fCurrentEfficiency;  //! efficiency
    TH1* fInitialConditions;  //! initial conditions
    TH1* fResult;             //! unfolding result

    // static variable to be accessed by MINUIT
    static TMatrixD* fgCorrelationMatrix;            //! contains fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationCovarianceMatrix;  //! contains the errors of fCurrentESD
    static TVectorD* fgCurrentESDVector;             //! contains fCurrentESD
    static TVectorD* fgEntropyAPriori;               //! a-priori distribution for entropy regularization

    static TF1* fgNBD;   //! negative binomial distribution

    static Int_t fgMaxParams;  //! bins in unfolded histogram = number of fit params
    static Int_t fgMaxInput;   //! bins in measured histogram

    // configuration params follow
    static RegularizationType fgRegularizationType; //! regularization that is used during Chi2 method
    static Float_t fgRegularizationWeight;          //! factor for regularization term
    static Int_t fgRegularizationRangeStart;        //! first bin where regularization is applied
    static Int_t fgRegularizationRangeEnd;          //! last bin + 1 where regularization is applied
    static Bool_t  fgCreateBigBin;                  //! to fix fluctuations at high multiplicities, all entries above a certain limit are summarized in one bin

    static Float_t fgBayesianSmoothing;             //! smoothing parameter (0 = no smoothing)
    static Int_t   fgBayesianIterations;            //! number of iterations in Bayesian method
    // end of configuration

 private:
    AliUnfolding(const AliUnfolding&);
    AliUnfolding& operator=(const AliUnfolding&);

  ClassDef(AliUnfolding, 0);
};

#endif

