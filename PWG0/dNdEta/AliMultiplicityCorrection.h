/* $Id$ */

#ifndef ALIMULTIPLICITYCORRECTION_H
#define ALIMULTIPLICITYCORRECTION_H

#include "TNamed.h"

//
// class that contains the correction matrix and the functions for
// correction the multiplicity spectrum
//

class TH1;
class TH2;
class TH1F;
class TH2F;
class TH3F;

class AliMultiplicityCorrection : public TNamed {
  public:
    AliMultiplicityCorrection();
    AliMultiplicityCorrection(const Char_t* name, const Char_t* title);
    virtual ~AliMultiplicityCorrection();

    virtual Long64_t Merge(TCollection* list);

    void FillMeasured(Float_t vtx, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20);
    void FillGenerated(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll);

    void FillCorrection(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20);

    Bool_t LoadHistograms(const Char_t* dir);
    void SaveHistograms();
    void DrawHistograms();
    void DrawComparison(const char* name, Int_t mcID, Int_t esdCorrId, Bool_t normalizeESD = kTRUE);

    void ApplyMinuitFit(Int_t inputRange, Bool_t fullPhaseSpace, Bool_t check = kFALSE);
    void ApplyMinuitFitAll();

    void ApplyBayesianMethod(Int_t inputRange, Bool_t fullPhaseSpace);

    void ApplyGaussianMethod(Int_t inputRange, Bool_t fullPhaseSpace);

    TH2F* GetMultiplicityESD(Int_t i) { return fMultiplicityESD[i]; }
    TH2F* GetMultiplicityMC(Int_t i) { return fMultiplicityMC[i]; }
    TH3F* GetCorrelation(Int_t i) { return fCorrelation[i]; }

    void SetMultiplicityESD(Int_t i, TH2F* hist) { fMultiplicityESD[i] = hist; }
    void SetMultiplicityMC(Int_t i, TH2F* hist) { fMultiplicityMC[i] = hist; }
    void SetCorrelation(Int_t i, TH3F* hist) { fCorrelation[i] = hist; }

    TH1F* GetMultiplicityESDCorrected(Int_t i) { return fMultiplicityESDCorrected[i]; }

    static void NormalizeToBinWidth(TH1* hist);
    static void NormalizeToBinWidth(TH2* hist);

  protected:
    enum { kESDHists = 4, kMCHists = 5, kCorrHists = 8 };

    static const Int_t fgMaxParams; // number of fit params

    static Double_t RegularizationPol0(Double_t *params);
    static Double_t RegularizationPol1(Double_t *params);
    static Double_t RegularizationTotalCurvature(Double_t *params);

    static void MinuitFitFunction(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t);

    void SetupCurrentHists(Int_t inputRange, Bool_t fullPhaseSpace);

    static TH1* fCurrentESD;         // static variable to be accessed by MINUIT
    static TH1* fCurrentCorrelation; // static variable to be accessed by MINUIT

    TH2F* fMultiplicityESD[kESDHists]; // multiplicity histogram: vtx vs multiplicity; array: |eta| < 0.5, 1, 1.5, 2 (0..3)
    TH2F* fMultiplicityMC[kMCHists];   // multiplicity histogram: vtx vs multiplicity; array: |eta| < 0.5, 1, 1.5, 2, inf (0..4)

    TH3F* fCorrelation[kCorrHists];    // vtx vs. (gene multiplicity) vs. (meas multiplicity); array: |eta| < 0.5, 1, 1.5, 2 (0..3 and 4..7), the first corrects to the eta range itself, the second to full phase space
    TH1F* fMultiplicityESDCorrected[kCorrHists]; // corrected histograms

 private:
    AliMultiplicityCorrection(const AliMultiplicityCorrection&);
    AliMultiplicityCorrection& operator=(const AliMultiplicityCorrection&);

  ClassDef(AliMultiplicityCorrection, 1);
};

#endif
