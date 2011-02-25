#ifndef ALISIGNIFICANCECALCULATOR_H
#define ALISIGNIFICANCECALCULATOR_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to calculate the statistical significance from          //
// AliMultiVeector objects for signal and background             //
// Origin: Francesco Prino (prino@to.infn.it)                    //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliMultiDimVector.h"
#include "TObject.h"

class AliSignificanceCalculator : public TObject{
 public:
  AliSignificanceCalculator();
  AliSignificanceCalculator(AliMultiDimVector* sig, AliMultiDimVector* bkg, Float_t normsig=1., Float_t normbkg=1.);
  AliSignificanceCalculator(AliMultiDimVector* sig, AliMultiDimVector* bkg, AliMultiDimVector* err2sig, AliMultiDimVector* err2bkg, Float_t normsig=1., Float_t normbkg=1.);

  ~AliSignificanceCalculator();

  void SetSignal(AliMultiDimVector* sig, Float_t norm=1.){
    if(fSignal) delete fSignal;
    fSignal=sig;
    fNormSig=norm;
    if(fSignal && fBackground) CalculateSignificance();
  }
  void SetBackground(AliMultiDimVector* bac, Float_t norm=1.){
    if(fBackground) delete fBackground;
    fBackground=bac;
    fNormBkg=norm;
    if(fSignal && fBackground) CalculateSignificance();
  }
  void SetErrSquareSignal(AliMultiDimVector* err2sig, Float_t norm=1.){
    if(fErrSquareSignal) delete fErrSquareSignal;
    fErrSquareSignal=err2sig;
    fNormSig=norm;
    if(fSignal && fBackground) CalculateSignificance();
  }
  void SetErrSquareBackground(AliMultiDimVector* err2bkg, Float_t norm=1.){
    if(fErrSquareBackground) delete fErrSquareBackground;
    fErrSquareBackground=err2bkg;
    fNormBkg=norm;
    if(fSignal && fBackground) CalculateSignificance();
  }
  
  void SetNormalizations(Float_t normSig, Float_t normBkg){
    fNormSig=normSig;
    fNormBkg=normBkg;
    if(fSignal && fBackground) CalculateSignificance();
  }

  AliMultiDimVector* GetSignal() const {return fSignal;}
  AliMultiDimVector* GetBackground() const {return fBackground;}
  AliMultiDimVector* GetSignificance() const {return fSignificance;}
  AliMultiDimVector* GetSignificanceError() const {return fErrSignificance;}

  void CalculateSignificance();
  Float_t GetMaxSignificance(Int_t* cutIndices, Int_t ptbin) const{
    Float_t sigMax=0;
    if(fSignificance) fSignificance->FindMaximum(sigMax,cutIndices,ptbin);
    return sigMax;
  }
  AliMultiDimVector* CalculatePurity() const;
  AliMultiDimVector* CalculatePurityError() const;
  AliMultiDimVector* CalculateSOverB() const;
  AliMultiDimVector* CalculateSOverBError() const;

 private:
  Bool_t Check() const;
  AliSignificanceCalculator(const AliSignificanceCalculator& c);
  AliSignificanceCalculator& operator=(const AliSignificanceCalculator& c);

  AliMultiDimVector* fSignal;              // signal matrix
  AliMultiDimVector* fErrSquareSignal;     // matrix with err^2 for signal
  AliMultiDimVector* fBackground;          // background matrix
  AliMultiDimVector* fErrSquareBackground; // matrix with err^2 for background
  AliMultiDimVector* fSignificance;        // significance matrix
  AliMultiDimVector* fErrSignificance;     // matrix with error on significance
  Float_t fNormSig;                        // signal normalization
  Float_t fNormBkg;                        // background normalization

  ClassDef(AliSignificanceCalculator,1); // class to compute and maximise significance

};

#endif
