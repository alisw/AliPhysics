// -*- C++ -*-

#ifndef _ALI_NONSEPARATION_MODEL_FIT_H_
#define _ALI_NONSEPARATION_MODEL_FIT_H_

#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMinuitMinimizer.h>
#include <Math/Functor.h>
#include <TMinuit.h>
#include <TFitter.h>
#include <Math/MinimizerOptions.h>

class TObjArray;
class TCut;
class TGraph;
class TTree;
class TString;

class AliNonseparationModelFit : public TObject {
public:
  AliNonseparationModelFit();
  virtual ~AliNonseparationModelFit();

  void Add(Int_t idx, TTree* filename, const TCut& cutMoments, TGraph *gRate);

  static const char* GetParName(Int_t idx);

  Double_t MinuitFunction(const Double_t *par);

  void SetFitToRates(Bool_t b) { fFitToRates = b; }
  Bool_t FitToRates() const { return fFitToRates; }

  void SetScaleRateError(Double_t s) { fScaleRateError = s; }

  void SetVar(Int_t idx, Double_t val, Double_t step, Double_t min, Double_t max);

  void DoFit( const TString &saveFileName);

  Double_t GetNDF()  const { return GetNDFMoments()  + GetNDFRates();  }
  Double_t GetChi2() const { return GetChi2Moments() + GetChi2Rates(); }

  Double_t GetNDFMoments()  const;
  Double_t GetNDFRates()    const;
  Double_t GetChi2Moments() const;
  Double_t GetChi2Rates()   const;

  TMinuitMinimizer& GetMinimizer() { return fMinimizer; }

protected:
private:
  AliNonseparationModelFit(const AliNonseparationModelFit& );
  AliNonseparationModelFit& operator=(const AliNonseparationModelFit& );

  TObjArray fMoments;           //!
  TObjArray fRates;             //!

  TVectorD  fPar;               //!

  Int_t     fNRates;            //!
  Bool_t    fFitToRates;        //!
  Bool_t    fFixCrossingAngles; //!
  Double_t  fScaleRateError;    //!

  Bool_t    fFitToOffsetScans;  //!

  TVectorD  fMuOffsetsX;        //!
  TVectorD  fMuOffsetsY;        //!

  TVectorD  fChi2Moments;       //!
  TVectorD  fChi2Rates;         //!

  TVectorD  fNDFMoments;        //!
  TVectorD  fNDFRates;          //!

  TMinuitMinimizer    fMinimizer; //!
  ROOT::Math::Functor fFcn;       //!

  ClassDef(AliNonseparationModelFit, 3);
} ;

#endif // _ALI_NONSEPARATION_MODEL_FIT_H_
