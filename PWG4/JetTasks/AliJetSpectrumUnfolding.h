#ifndef ALIJETSPECTRUMUNFOLDING_H
#define ALIJETSPECTRUMUNFOLDING_H

#include "TNamed.h"
//
// class that contains the correction matrix and the functions for
// correction the jet spectrum
// implements 2-dim bayesian method
//

class TH1;
class TH2;
class TH3;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TF2;
class TCollection;

#include <TMatrixD.h>
#include <TVectorD.h>
#include <THnSparse.h>

class AliJetSpectrumUnfolding : public TNamed {
  public:
    AliJetSpectrumUnfolding();
    AliJetSpectrumUnfolding(const Char_t* name, const Char_t* title);
    virtual ~AliJetSpectrumUnfolding();

    virtual Long64_t Merge(TCollection* list);

    Bool_t LoadHistograms(const Char_t* dir = 0);
    void SaveHistograms();
    void DrawHistograms();
    void DrawComparison(const char* name, TH2* genSpec);

    void SetBayesianParameters(Float_t smoothing, Int_t nIterations);
    void ApplyBayesianMethod(Float_t regPar = 1, Int_t nIterations = 100, TH2* initialConditions = 0, Bool_t determineError = kTRUE);

    TH2F* GetRecSpectrum()      { return fRecSpectrum; }
    TH2F* GetGenSpectrum()      { return fGenSpectrum; }
    TH2F* GetUnfSpectrum()      { return fUnfSpectrum; }
    THnSparseF* GetCorrelation(){ return fCorrelation; }

    void SetRecSpectrum(TH2F* hist)      { fRecSpectrum  = hist; }
    void SetGenSpectrum(TH2F* hist)      { fGenSpectrum  = hist; }
    void SetUnfSpectrum(TH2F* hist)      { fUnfSpectrum  = hist; }
    void SetCorrelation(THnSparseF* hist){ fCorrelation  = hist; }

    void SetGenRecFromFunc(TF2* inputGen);
    TH2F* CalculateRecSpectrum(TH2* inputGen);

    static void NormalizeToBinWidth(TH2* hist);

    void GetComparisonResults(Float_t* gen = 0, Int_t* genLimit = 0, Float_t* residuals = 0, Float_t* ratioAverage = 0) const;

  protected:
    void SetupCurrentHists(Bool_t createBigBin);

    static Double_t BayesCov(THnSparseF* M, THnSparseF* correlation, Int_t* binTM, Int_t* binTM1);
    static Double_t BayesUncertaintyTerms(THnSparseF *M, THnSparseF *C, Int_t* binTM, Int_t* binTM1, Double_t Nt);
    static Int_t UnfoldWithBayesian(THnSparseF* correlation, TH2* measured, TH2* initialConditions, TH2* aResult, Float_t regPar, Int_t nIterations, Bool_t calculateErrors = kFALSE);

    static Float_t fgBayesianSmoothing;             //! smoothing parameter (0 = no smoothing)
    static Int_t   fgBayesianIterations;            //! number of iterations in Bayesian method

    TH2F* fCurrentRec;
    THnSparseF* fCurrentCorrelation;

    TH2F* fRecSpectrum;
    TH2F* fGenSpectrum;
    TH2F* fUnfSpectrum;
    THnSparseF* fCorrelation;

    Float_t fLastChi2MC;        //! last Chi2 between MC and unfolded ESD (calculated in DrawComparison)
    Int_t   fLastChi2MCLimit;   //! bin where the last chi2 breached a certain threshold
    Float_t fLastChi2Residuals; //! last Chi2 of the ESD and the folded unfolded ESD (calculated in DrawComparison)
    Float_t fRatioAverage;      //! last average of |ratio-1| where ratio = unfolded / mc (bin 2..150)

 private:
    AliJetSpectrumUnfolding(const AliJetSpectrumUnfolding&);
    AliJetSpectrumUnfolding& operator=(const AliJetSpectrumUnfolding&);

  ClassDef(AliJetSpectrumUnfolding, 2);
};

#endif

