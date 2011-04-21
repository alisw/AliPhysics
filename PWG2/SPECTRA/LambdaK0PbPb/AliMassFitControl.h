#ifndef ALIMASSFITCONTROL_H
#define ALIMASSFITCONTROL_H
#include "TMath.h"

//-------------------------------------------------------------------------
// AliMassFitControl
// This class controls the parameters (fit range, polynomial order, rebinning fact)
// for fitting an invariant mass histogram.
// A typical use is for fitting a number of pt bins projected from a 2D pt-M histo.
//-------------------------------------------------------------------------
class AliMassFitControl : public TObject {
 public:
  AliMassFitControl(){
    fPtUpper=0.0; fPtLower=0.0; fPolyOrder = 2; fRebinFactor = 1;
    fBinLower=0;fBinUpper=0;
  };
  AliMassFitControl(Float_t ptLo, Float_t ptUp, Int_t polyO, Int_t RBF){
    fPtUpper=ptUp; fPtLower=ptLo; fPolyOrder = polyO; fRebinFactor = RBF;
    fBinLower=0;fBinUpper=0;
  };
  AliMassFitControl(Double_t ptLo, Double_t ptUp, Int_t polyO, Int_t RBF){
    fPtUpper=ptUp; fPtLower=ptLo; fPolyOrder = polyO; fRebinFactor = RBF;
    fMinMass=1.085; fMaxMass=1.17;
    fBinLower=0;fBinUpper=0;
  };
  ///// the minimum value of parameter have to be zero with this constructor
  AliMassFitControl(Double_t ptLo, Double_t ptUp, Int_t polyO, Int_t RBF, Double_t m1, Double_t m2){
    fPtUpper=ptUp; fPtLower=ptLo; fPolyOrder = polyO; fRebinFactor = RBF;
    fMinMass=m1; fMaxMass=m2; fhistMin = 0;
    fBinLower=0;fBinUpper=0;
  };
  ///// the minimum value of (QA) parameter don't have to be zero. can by arbitrary value "histMin"
  AliMassFitControl(Double_t histMin, Double_t ptLo, Double_t ptUp, Int_t polyO, Int_t RBF, Double_t m1, Double_t m2){
    fhistMin = histMin; fPtUpper=ptUp; fPtLower=ptLo; fPolyOrder = polyO; fRebinFactor = RBF;
    fMinMass=m1; fMaxMass=m2;
    fBinLower=0;fBinUpper=0;
  };
  AliMassFitControl(Int_t BinLo, Int_t BinUp, Int_t polyO, Int_t RBF){
    fPtUpper=0.0; fPtLower=0.0; fPolyOrder = polyO; fRebinFactor = RBF;
    fBinLower=BinLo;fBinUpper=BinUp;
  };
  ~AliMassFitControl(){;};

  // Functions for sorting
  Bool_t IsEqual(const TObject *obj) const {return fPtLower == ((AliMassFitControl*)obj)->fPtLower;}; //Not sure whether this one reqd for sorting
  Bool_t IsSortable() const { return kTRUE; };
  Int_t Compare(const TObject *obj) const
    {
    if ( fPtLower < ((AliMassFitControl*)obj)->fPtLower) return -1;
    if (fPtLower > ((AliMassFitControl*)obj)->fPtLower) return 1;
    return 0;
  };

  Int_t RebinFactor() const {return fRebinFactor;};
  Int_t BinLower() const {return fBinLower;};
  Int_t BinUpper() const {return fBinUpper;};
  Double_t PtUpper() const {return fPtUpper;};
  Double_t PtLower() const {return fPtLower;};
  Double_t MinMass() const {return fMinMass;};
  Double_t MaxMass() const {return fMaxMass;};


  Bool_t FixedQuad() const {if(fPolyOrder < 2) {return kTRUE;}
  else {return kFALSE;};
  };
  Bool_t FixedLin() const {if(fPolyOrder < 1) {return kTRUE;}
  else {return kFALSE;};
  };
  Double_t DPt() const {return fPtUpper-fPtLower;};
  void CalcBinLimits(Int_t BinsPerGeV){
//    fBinLower = TMath::Nint(1.+fPtLower*BinsPerGeV); //this assumes that histogram allways start from zero. not true for QA
    fBinLower = TMath::Nint(1.+(fPtLower-fhistMin)*BinsPerGeV); //this is valid if hist. doesn't start from zero but from number: histMin
    fBinUpper = TMath::Nint((fPtUpper-fhistMin)*BinsPerGeV); //this is valid if hist. doesn't start from zero but from number: histMin    
  };

protected:
  Double_t fhistMin; //value of the first bin for QA parameter
  Double_t fPtUpper; // Upper pt limit
  Double_t fPtLower; // Lower pt limit. This was previously public due to use in sorting functions but it seems OK here.
  Int_t fBinUpper; // Upper bin limit
  Int_t fBinLower; // Lower bin limit
  Int_t fPolyOrder; // Polynomial order - 0=constant, 1=linear, 2=quadratic, others not supported
  Int_t fRebinFactor; // Rebinning factor
  Double_t fMinMass; // Minimum to use as fitting range
  Double_t fMaxMass; // Maximum to use as fitting range
    
  
  ClassDef(AliMassFitControl,1);
};
#endif
