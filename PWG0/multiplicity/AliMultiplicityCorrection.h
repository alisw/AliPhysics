/* $Id$ */

#ifndef ALIMULTIPLICITYCORRECTION_H
#define ALIMULTIPLICITYCORRECTION_H

#include "TNamed.h"

//
// class that contains the correction matrix and the functions for
// correction the multiplicity spectrum
// implements a several unfolding methods: e.g. chi2 minimization and bayesian unfolding
//

class TH1;
class TH2;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TCollection;

// defined here, because it does not seem possible to predeclare these (or i do not know how)
// -->
// $ROOTSYS/include/TVectorDfwd.h:21: conflicting types for `typedef struct TVectorT<Double_t> TVectorD'
// PWG0/dNdEta/AliMultiplicityCorrection.h:21: previous declaration as `struct TVectorD'

#include <TMatrixD.h>
#include <TVectorD.h>
#include <AliPWG0Helper.h>

class AliMultiplicityCorrection : public TNamed {
  public:
    enum EventType { kTrVtx = 0, kMB, kINEL, kNSD };
    enum RegularizationType { kNone = 0, kPol0, kPol1, kLog, kEntropy, kCurvature };
    enum MethodType { kChi2Minimization = 0, kBayesian = 1 };
    enum { kESDHists = 4, kMCHists = 5, kCorrHists = 8, kQualityRegions = 3 };

    AliMultiplicityCorrection();
    AliMultiplicityCorrection(const Char_t* name, const Char_t* title);
    virtual ~AliMultiplicityCorrection();
    
    static AliMultiplicityCorrection* Open(const char* fileName, const char* folderName = "Multiplicity");

    virtual Long64_t Merge(TCollection* list);

    void FillMeasured(Float_t vtx, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20);
    void FillGenerated(Float_t vtx, Bool_t triggered, Bool_t vertex, AliPWG0Helper::MCProcessType processType, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll);

    void FillCorrection(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20);

    Bool_t LoadHistograms(const Char_t* dir = 0);
    void SaveHistograms(const char* dir = 0);
    void DrawHistograms();
    void DrawComparison(const char* name, Int_t inputRange, Bool_t fullPhaseSpace, Bool_t normalizeESD, TH1* mcHist, Bool_t simple = kFALSE);

    Int_t ApplyMinuitFit(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Bool_t check = kFALSE, TH1* initialConditions = 0);

    static void SetRegularizationParameters(RegularizationType type, Float_t weight, Int_t minuitParams = -1);
    static void SetBayesianParameters(Float_t smoothing, Int_t nIterations);
    static void SetCreateBigBin(Bool_t flag) { fgCreateBigBin = flag; }

    void ApplyNBDFit(Int_t inputRange, Bool_t fullPhaseSpace);

    void ApplyBayesianMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Float_t regPar = 1, Int_t nIterations = 100, TH1* initialConditions = 0, Bool_t determineError = kTRUE);

    static TH1* CalculateStdDev(TH1** results, Int_t max);
    TH1* StatisticalUncertainty(MethodType methodType, Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Bool_t randomizeMeasured, Bool_t randomizeResponse, TH1* compareTo = 0);

    void ApplyGaussianMethod(Int_t inputRange, Bool_t fullPhaseSpace);

    void ApplyLaszloMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType);

    TH2F* GetMultiplicityESD(Int_t i) { return fMultiplicityESD[i]; }
    TH2F* GetMultiplicityVtx(Int_t i) { return fMultiplicityVtx[i]; }
    TH2F* GetMultiplicityMB(Int_t i) { return fMultiplicityMB[i]; }
    TH2F* GetMultiplicityINEL(Int_t i) { return fMultiplicityINEL[i]; }
    TH2F* GetMultiplicityNSD(Int_t i) { return fMultiplicityNSD[i]; }
    TH2F* GetMultiplicityMC(Int_t i, EventType eventType);
    TH3F* GetCorrelation(Int_t i) { return fCorrelation[i]; }
    TH1F* GetMultiplicityESDCorrected(Int_t i) { return fMultiplicityESDCorrected[i]; }

    void SetMultiplicityESD(Int_t i, TH2F* hist)  { fMultiplicityESD[i]  = hist; }
    void SetMultiplicityVtx(Int_t i, TH2F* hist)  { fMultiplicityVtx[i]  = hist; }
    void SetMultiplicityMB(Int_t i, TH2F* hist)   { fMultiplicityMB[i]   = hist; }
    void SetMultiplicityINEL(Int_t i, TH2F* hist) { fMultiplicityINEL[i] = hist; }
    void SetMultiplicityNSD(Int_t i, TH2F* hist) { fMultiplicityNSD[i] = hist; }
    void SetMultiplicityMC(Int_t i, EventType eventType, TH2F* hist);
    void SetCorrelation(Int_t i, TH3F* hist) { fCorrelation[i] = hist; }
    void SetMultiplicityESDCorrected(Int_t i, TH1F* hist) { fMultiplicityESDCorrected[i] = hist; }

    void SetGenMeasFromFunc(TF1* inputMC, Int_t id);
    TH2F* CalculateMultiplicityESD(TH1* inputMC, Int_t correlationMap);

    void GetComparisonResults(Float_t* mc = 0, Int_t* mcLimit = 0, Float_t* residuals = 0, Float_t* ratioAverage = 0) const;

    TH1* GetEfficiency(Int_t inputRange, EventType eventType);
    TH1* GetTriggerEfficiency(Int_t inputRange);

    static void SetQualityRegions(Bool_t SPDStudy);
    Float_t GetQuality(Int_t region) const { return fQuality[region]; }

    void FFT(Int_t dir, Int_t m, Double_t *x, Double_t *y);

  protected:
    static const Int_t fgkMaxParams;  //! bins in unfolded histogram = number of fit params
    static const Int_t fgkMaxInput;   //! bins in measured histogram

    static Double_t RegularizationPol0(TVectorD& params);
    static Double_t RegularizationPol1(TVectorD& params);
    static Double_t RegularizationTotalCurvature(TVectorD& params);
    static Double_t RegularizationEntropy(TVectorD& params);
    static Double_t RegularizationLog(TVectorD& params);

    static void MinuitFitFunction(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t);
    static void MinuitNBD(Int_t& unused1, Double_t* unused2, Double_t& chi2, Double_t *params, Int_t unused3);

    static void DrawGuess(Double_t *params);

    void SetupCurrentHists(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Bool_t createBigBin);

    Float_t BayesCovarianceDerivate(Float_t matrixM[251][251], TH2* hResponse, Int_t k, Int_t i, Int_t r, Int_t u);
    static Int_t UnfoldWithBayesian(TH1* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* aResult, Float_t regPar, Int_t nIterations);
    static Int_t UnfoldWithMinuit(TH1* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check);

    TH1* fCurrentESD;         //! current input esd
    TH1* fCurrentCorrelation; //! current correlation
    TH1* fCurrentEfficiency;  //! current efficiency

    // static variable to be accessed by MINUIT
    static TMatrixD* fgCorrelationMatrix;            //! contains fCurrentCorrelation in matrix form
    static TMatrixD* fgCorrelationCovarianceMatrix;  //! contains the errors of fCurrentESD
    static TVectorD* fgCurrentESDVector;             //! contains fCurrentESD
    static TVectorD* fgEntropyAPriori;               //! a-priori distribution for entropy regularization

    static TF1* fgNBD;   //! negative binomial distribution

    // configuration params follow
    static RegularizationType fgRegularizationType; //! regularization that is used during Chi2 method
    static Float_t fgRegularizationWeight;          //! factor for regularization term
    static Bool_t  fgCreateBigBin;                  //! to fix fluctuations at high multiplicities, all entries above a certain limit are summarized in one bin
    static Int_t   fgNParamsMinuit;                 //! number of parameters minuit uses for unfolding (todo: to be merged w/ fgkMaxParams that has to be const. for the moment)

    static Float_t fgBayesianSmoothing;             //! smoothing parameter (0 = no smoothing)
    static Int_t   fgBayesianIterations;            //! number of iterations in Bayesian method
    // end of configuration

    TH2F* fMultiplicityESD[kESDHists]; // multiplicity histogram: vtx vs multiplicity; array: |eta| < 0.5, 1.0, 1.5, 2 (0..3)

    TH2F* fMultiplicityVtx[kMCHists];  // multiplicity histogram of events that have a reconstructed vertex : vtx vs multiplicity; array: |eta| < 0.5, 1.0, 1.5, 2, inf (0..4)
    TH2F* fMultiplicityMB[kMCHists];   // multiplicity histogram of triggered events                        : vtx vs multiplicity; array: |eta| < 0.5, 1.0, 1.5, 2, inf (0..4)
    TH2F* fMultiplicityINEL[kMCHists]; // multiplicity histogram of all (inelastic) events                  : vtx vs multiplicity; array: |eta| < 0.5, 1.0, 1.5, 2, inf (0..4)
    TH2F* fMultiplicityNSD[kMCHists]; // multiplicity histogram of NSD events                  : vtx vs multiplicity; array: |eta| < 0.5, 1.0, 1.5, 2, inf (0..4)

    TH3F* fCorrelation[kCorrHists];              // vtx vs. (gene multiplicity (trig+vtx)) vs. (meas multiplicity); array: |eta| < 0.5, 1, 1.5, 2 (0..3 and 4..7), the first corrects to the eta range itself, the second to full phase space

    TH1F* fMultiplicityESDCorrected[kCorrHists]; // corrected histograms

    Int_t fLastBinLimit;        //! last bin limit, determined in SetupCurrentHists()
    Float_t fLastChi2MC;        //! last Chi2 between MC and unfolded ESD (calculated in DrawComparison)
    Int_t   fLastChi2MCLimit;   //! bin where the last chi2 breached a certain threshold, used to evaluate the multiplicity reach (calc. in DrawComparison)
    Float_t fLastChi2Residuals; //! last Chi2 of the ESD and the folded unfolded ESD (calculated in DrawComparison)
    Float_t fRatioAverage;      //! last average of |ratio-1| where ratio = unfolded / mc (bin 2..150)

    static Int_t   fgQualityRegionsB[kQualityRegions]; //! begin, given in multiplicity units
    static Int_t   fgQualityRegionsE[kQualityRegions]; //! end
    Float_t fQuality[kQualityRegions];                 //! stores the quality of the last comparison (calculated in DrawComparison). Contains 3 values that are averages of (MC - unfolded) / e(MC) in 3 regions, these are defined in fQualityRegionB,E

 private:
    AliMultiplicityCorrection(const AliMultiplicityCorrection&);
    AliMultiplicityCorrection& operator=(const AliMultiplicityCorrection&);

  ClassDef(AliMultiplicityCorrection, 4);
};

#endif

