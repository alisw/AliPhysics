#ifndef ALIL3TRANSFORM_H
#define ALIL3TRANSFORM_H

#include "AliL3RootTypes.h"
class TFile;

class AliL3Transform {
 private:
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
  Double_t fNRotShift;
  Double_t fZLength;
  Double_t fZOffset;
  Double_t fCos[36]; //fill this following Init
  Double_t fSin[36]; //fill this following Init
  Double_t fX[176];  //fill this following Init
  Int_t fNPads[176]; //fill this following Init
  Int_t fVersion; //flags which version one is using
  void Init(); //old init used by Anders for AliRoot <= 3.06
 public:
  AliL3Transform();
  AliL3Transform(const char *pathname);
  virtual ~AliL3Transform();
  Int_t GetVersion(){return fVersion;}
  void Init(const Char_t* path); //new init for all AliRoot versions

  Double_t GetPadPitchWidthLow() {return fPadPitchWidthLow;}
  Double_t GetPadPitchWidthUp() {return fPadPitchWidthUp;}
  Double_t GetPadPitchWidth(Int_t patch) {if(patch<=2) return fPadPitchWidthLow; else return fPadPitchWidthUp;}
  Double_t GetZWidth() {return fZWidth;}
  
  Bool_t Slice2Sector(Int_t slice, Int_t slicerow, Int_t &sector, Int_t &row) const;

  Bool_t Sector2Slice(Int_t &slice, Int_t sector) const;
  Bool_t Sector2Slice(Int_t &slice, Int_t &slicerow, Int_t sector, Int_t row) const;
  
  Double_t Row2X(Int_t slicerow);
  Int_t GetNPads(Int_t row){return (row<fNRow)?fNPads[row]:0;}
  Int_t GetNTimeBins(){return fNTimeBins;}

  Double_t GetEta(Float_t *xyz);
  Double_t GetEta(Int_t row, Int_t pad, Int_t time);
  Double_t GetPhi(Float_t *xyz);
  Double_t GetMaxY(Int_t slicerow);
  void XYZtoRPhiEta(Float_t *rpe, Float_t *xyz);
  void Local2Global(Float_t *xyz, Int_t slice);
  void Local2GlobalAngle(Float_t *angle, Int_t slice);
  void Global2LocalAngle(Float_t *angle, Int_t slice);

  void Raw2Local(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  void Local2Global(Float_t *xyz, Int_t sector, Int_t row);
  void Global2Local(Float_t *xyz, Int_t sector, Bool_t isSlice=kFALSE);
  void Raw2Global(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  void Local2Raw(Float_t *xyz, Int_t sector, Int_t row);
  void Global2Raw(Float_t *xyz, Int_t sector, Int_t row);
  
  ClassDef(AliL3Transform,1) //Transformation class for ALICE TPC
};


#endif





