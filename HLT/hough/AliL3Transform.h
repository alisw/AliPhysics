#ifndef ALIL3TRANSFORM_H
#define ALIL3TRANSFORM_H



//#include "AliTPCParam.h"
#include "AliL3RootTypes.h"
class TFile;

class AliL3Transform {
 private:
//  AliTPCParam *fParam;
  Int_t fNTimeBins;
  Int_t fNRowLow;
  Int_t fNRowUp;
  Int_t fNSectorLow;
  Int_t fNSectorUp;
  Double_t fPadPitchWidthLow;
  Double_t fPadPitchWidthUp;
  Double_t fZWidth;
  Double_t fZSigma;
  Int_t fNSector;
  Int_t fNSlice;
  Int_t fNRow;
  Double_t fPi;
  Double_t fCos[36];
  Double_t fSin[36];
  Double_t fX[174];
  Int_t fNPads[174];
 public:
  AliL3Transform();
  virtual ~AliL3Transform();
  void Init();

  Double_t GetPadPitchWidthLow() {return fPadPitchWidthLow;}
  Double_t GetPadPitchWidthUp() {return fPadPitchWidthUp;}

  Bool_t Slice2Sector(Int_t slice, Int_t slicerow, Int_t & sector, Int_t &row) const;

  Bool_t Sector2Slice(Int_t & slice, Int_t  sector) const;
  Bool_t Sector2Slice(Int_t & slice, Int_t & slicerow,Int_t  sector, Int_t row) const;
  
  Double_t Row2X(Int_t slicerow);
  Int_t GetNPads(Int_t row){return (row<174)?fNPads[row]:0;}
  Int_t GetNTimeBins(){return fNTimeBins;}

  Double_t GetEta(Float_t *xyz);
  Double_t GetPhi(Float_t *xyz);
  Double_t GetMaxY(Int_t slicerow);
  void Local2Global(Float_t *xyz,Int_t slice);
  void Local2GlobalAngle(Float_t *angle,Int_t slice);
  void Global2LocalAngle(Float_t *angle,Int_t slice);

  void Raw2Local(Float_t *xyz,Int_t sector,Int_t row,Float_t pad,Float_t time);
  void Local2Global(Float_t *xyz,Int_t sector,Int_t row);
  void Global2Local(Float_t *xyz,Int_t sector);
  void Raw2Global(Float_t *xyz,Int_t sector,Int_t row,Float_t pad,Float_t time);
  void Local2Raw(Float_t *xyz,Int_t sector,Int_t row);
  void Global2Raw(Float_t *xyz,Int_t sector,Int_t row);
  
  ClassDef(AliL3Transform,1)
};


#endif
