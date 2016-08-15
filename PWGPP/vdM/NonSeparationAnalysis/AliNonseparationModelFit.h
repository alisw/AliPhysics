// -*- C++ -*-

#ifndef _ALI_NONSEPARATION_MODEL_FIT_H_
#define _ALI_NONSEPARATION_MODEL_FIT_H_

#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>

class TObjArray;
class TGraph;
class TTree;
class TString;

class AliNonseparationModelFit : public TObject {
public:
  AliNonseparationModelFit();
  virtual ~AliNonseparationModelFit();

  void Add(Int_t idx, TTree* filename, TGraph *gRate);

  static const char* GetParName(Int_t idx);

  Double_t MinuitFunction(const Double_t *par);
  
  void SetFitToRates(Bool_t b) { fFitToRates = b; }
  Bool_t FitToRates() const { return fFitToRates; }

  void SetScaleRateError(Double_t s) { fScaleRateError = s; }

  void DoFit(TVectorD& startPar, const TString &saveFileName);

  Double_t GetNDF()  const { return GetNDFMoments()  + GetNDFRates();  }
  Double_t GetChi2() const { return GetChi2Moments() + GetChi2Rates(); }

  Double_t GetNDFMoments()  const;
  Double_t GetNDFRates()    const;
  Double_t GetChi2Moments() const;
  Double_t GetChi2Rates()   const;

  static Bool_t IsGoodLumiRegionFit(Int_t k, Int_t scanType, const TVectorD& sep,
				    const TVectorD& mom, const TMatrixDSym& cov);

protected:
private:
  AliNonseparationModelFit(const AliNonseparationModelFit& );
  AliNonseparationModelFit& operator=(const AliNonseparationModelFit& );

  TObjArray *fMoments;       //!
  TObjArray *fRates;         //!

  TVectorD  fPar;            //!

  Int_t     fNRates;         //!
  Bool_t    fFitToRates;     //!
  Double_t  fScaleRateError; //!

  Bool_t    fFitToOffsetScans; //!

  TVectorD  fMuOffsetsX;     //!
  TVectorD  fMuOffsetsY;     //!

  Double_t  fChi2Moments[6]; //!
  Double_t  fChi2Rates[6];   //!

  Double_t  fNDFMoments[6]; //!
  Double_t  fNDFRates[6];   //!

  ClassDef(AliNonseparationModelFit, 1);
} ;

#endif // _ALI_NONSEPARATION_MODEL_FIT_H_
