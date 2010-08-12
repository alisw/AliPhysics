//
// Class AliRsnFitResult
//
// It contains all results from a fit to an invariant mass spectrum
//

#ifndef ALIRSNFITRESULT
#define ALIRSNFITRESULT

#include <TMath.h>

class AliRsnFitResult : public TNamed
{
  public:
  
    enum
    {
      kSgAmp,
      kSgCenter,
      kSgWidth,
      kSgIntegral,
      kBgIntegral,
      kSg2Bg,
      kSignificance,
      kTotalIntegral,
      kSumIntegral,
      kSumIntegralError,
      kChi2,
      kNDF,
      kValues
    };
    
    enum ESignalType
    {
      kBreitWigner,
      kGaus
    };
    
    enum EBackgroundType
    {
      kSquareRoot,
      kLinear,
      kPoly2
    };
  
    AliRsnFitResult(const char *name = "none", ESignalType sType = kGaus, EBackgroundType bType = kPoly2);
    //AliRsnFitResult(const AliRsnFitResult& copy);
    //AliRsnFitResult& operator=(cont AliRsnFitResult& copy);
    
    Double_t NormSquareRoot(Double_t *x, Double_t *par);
    Double_t NormPoly2(Double_t *x, Double_t *par);
    Double_t NormLinear(Double_t *x, Double_t *par);
    Double_t NormBreitWigner(Double_t *x, Double_t *par);
    Double_t NormGaus(Double_t *x, Double_t *par);
    
    Double_t Signal(Double_t *x, Double_t *par);
    Double_t Background(Double_t *x, Double_t *par);
    Double_t Sum(Double_t *x, Double_t *par);
    
    Int_t    GetNParBackground();
    Bool_t   InitFunctions();
    Bool_t   SingleFit(const char *fitOpts, Double_t mass, Double_t width);
    
    void     SetHistogram(TH1F *histogram);
    void     SetPeakRange(Double_t min, Double_t max) {fPeakRange[0] = TMath::Min(min, max); fPeakRange[1] = TMath::Max(min, max);}
    void     SetFullRange(Double_t min, Double_t max) {fFullRange[0] = TMath::Min(min, max); fFullRange[1] = TMath::Max(min, max);}
  
  private:
  
    ESignalType      fSignalType;          // signal function
    EBackgroundType  fBackgroundType;      // background function
    
    Double_t         fFullRange[2];        // full fit range
    Double_t         fPeakRange[2];        // peak fit range
    Double_t         fFitResult[kValues];  // array of fit results
    
    TH1D            *fHistogram;           // histogram to be fitted
    TF1             *fSignalFcn;           // function for signal
    TF1             *fBackgroundFcn;       // function for background
    TF1             *fSumFcn;              // total function
  
    ClassDef(AliRsnFitResult, 1)
};

#endif
