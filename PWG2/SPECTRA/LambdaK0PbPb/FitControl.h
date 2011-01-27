#include "TMath.h"

class FitControl : public TObject {
 protected:
  Double_t mPtUpper; // Upper and lower pt limits
  Int_t mBinUpper; // Upper and lower bin limits
  Int_t mBinLower;
  Int_t mPolyOrder; // Polynomial order - 0=constant, 1=linear, 2=quadratic, others not supported
  Int_t mRebinFactor; // Rebinning factor
  Double_t mMinMass; // Min and max to use as fitting range
  Double_t mMaxMass;

 public:
  // This data member moved here due to use in sorting functions
  Double_t mPtLower;

  FitControl(){
    mPtUpper=0.0; mPtLower=0.0; mPolyOrder = 2; mRebinFactor = 1;
    mBinLower=0;mBinUpper=0;
  };
  FitControl(Float_t ptLo, Float_t ptUp, Int_t polyO, Int_t RBF){
    mPtUpper=ptUp; mPtLower=ptLo; mPolyOrder = polyO; mRebinFactor = RBF;
    mBinLower=0;mBinUpper=0;
  };
  FitControl(Double_t ptLo, Double_t ptUp, Int_t polyO, Int_t RBF){
    mPtUpper=ptUp; mPtLower=ptLo; mPolyOrder = polyO; mRebinFactor = RBF;
    mMinMass=1.085; mMaxMass=1.17;
    mBinLower=0;mBinUpper=0;
  };
  FitControl(Double_t ptLo, Double_t ptUp, Int_t polyO, Int_t RBF, Double_t m1, Double_t m2){
    mPtUpper=ptUp; mPtLower=ptLo; mPolyOrder = polyO; mRebinFactor = RBF;
    mMinMass=m1; mMaxMass=m2;
    mBinLower=0;mBinUpper=0;
  };
  FitControl(Int_t BinLo, Int_t BinUp, Int_t polyO, Int_t RBF){
    mPtUpper=0.0; mPtLower=0.0; mPolyOrder = polyO; mRebinFactor = RBF;
    mBinLower=BinLo;mBinUpper=BinUp;
  };
  ~FitControl(){;};

  // Functions for sorting
  Bool_t IsEqual(const TObject *obj) const {return mPtLower == ((FitControl*)obj)->mPtLower;}; //Not sure whether this one reqd for sorting
  Bool_t IsSortable() const { return kTRUE; };
  Int_t Compare(const TObject *obj) const
    {
    if ( mPtLower < ((FitControl*)obj)->mPtLower) return -1;
    if (mPtLower > ((FitControl*)obj)->mPtLower) return 1;
    return 0;
  };

  Int_t RebinFactor(){return mRebinFactor;};
  Int_t BinLower(){return mBinLower;};
  Int_t BinUpper(){return mBinUpper;};
  Double_t PtUpper(){return mPtUpper;};
  Double_t PtLower(){return mPtLower;};
  Double_t MinMass(){return mMinMass;};
  Double_t MaxMass(){return mMaxMass;};


  Bool_t FixedQuad(){if(mPolyOrder < 2) {return kTRUE;}
  else {return kFALSE;};
  };
  Bool_t FixedLin(){if(mPolyOrder < 1) {return kTRUE;}
  else {return kFALSE;};
  };
  Double_t DPt(){return mPtUpper-mPtLower;};
  void CalcBinLimits(Int_t BinsPerGeV){
    mBinLower = TMath::Nint(1.+mPtLower*BinsPerGeV);
    mBinUpper = TMath::Nint(mPtUpper*BinsPerGeV);
  };

  ClassDef(FitControl,1);
};
