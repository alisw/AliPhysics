#ifndef ALIJETSPECTRUMUNFOLDING_H
#define ALIJETSPECTRUMUNFOLDING_H
//
// class that contains the correction matrix and the functions for
// correction the jet spectrum
// implements 2-dim bayesian method
//




class TH1;
class TH2;
class TH3;
class TH1F;
class TH3F;
class TF1;
class TF2;
class TCollection;

#include "TNamed.h"
#include <TH2F.h>
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
    void DrawComparison(const char* name, TH2* const genSpec);

    void SetBayesianParameters(Float_t smoothing, Int_t nIterations);
    void ApplyBayesianMethod(Float_t regPar = 1, Int_t nIterations = 100, TH2* initialConditions = 0, Bool_t determineError = kTRUE);

    TH2F* GetRecSpectrum() const      { return fRecSpectrum; }
    TH2F* GetGenSpectrum() const     { return fGenSpectrum; }
    TH2F* GetUnfSpectrum() const     { return fUnfSpectrum; }
    THnSparseF* GetCorrelation() const { return fCorrelation; }

    void SetRecSpectrum(TH2F* const hist)   { if(fRecSpectrum) delete fRecSpectrum;
      fRecSpectrum  = hist; }
    void SetGenSpectrum(TH2F* const hist)      { if(fGenSpectrum)delete fGenSpectrum;
      fGenSpectrum  = hist; }
    void SetUnfSpectrum(TH2F*  const hist)      { if(fUnfSpectrum)delete fUnfSpectrum;
      fUnfSpectrum  = hist; }
    void SetCorrelation(THnSparseF* const hist){ if(fCorrelation)delete fCorrelation;
      fCorrelation  = hist; }

    void SetGenRecFromFunc(TF2* inputGen);
    TH2F* CalculateRecSpectrum(TH2* inputGen);

    static void NormalizeToBinWidth(TH2* const hist);

    void GetComparisonResults(Float_t* const gen = 0, Int_t* const genLimit = 0, Float_t* const residuals = 0, Float_t* const ratioAverage = 0) const;

  protected:
    void SetupCurrentHists(Bool_t createBigBin);

    static Double_t BayesCov(THnSparseF* const M, THnSparseF* const correlation, Int_t* const binTM, Int_t* const binTM1);
    static Double_t BayesUncertaintyTerms(THnSparseF* const M, THnSparseF *const C, Int_t* const binTM, Int_t* const binTM1, Double_t nt);
    static Int_t UnfoldWithBayesian(THnSparseF* const correlation, TH2* const measured, TH2* const initialConditions, TH2* const aResult, Float_t regPar, Int_t nIterations, Bool_t calculateErrors = kFALSE);

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

    static const Int_t fgkNBINSE;
    static const Int_t fgkNBINSZ;
    static const Int_t fgkNEVENTS;
    static const Double_t fgkaxisLowerLimitE;
    static const Double_t fgkaxisLowerLimitZ;
    static const Double_t fgkaxisUpperLimitE;
    static const Double_t fgkaxisUpperLimitZ;



    AliJetSpectrumUnfolding(const AliJetSpectrumUnfolding&);
    AliJetSpectrumUnfolding& operator=(const AliJetSpectrumUnfolding&);

  ClassDef(AliJetSpectrumUnfolding, 2);
};

#endif

