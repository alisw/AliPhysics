#ifndef ALIL3TRANSFORM_H
#define ALIL3TRANSFORM_H

#include "AliL3RootTypes.h"

class AliL3Transform {

 private:
  const static Double_t fBFACT;
  
  static Double_t fBField;
  static Int_t fBFieldFactor;
  static Int_t fNTimeBins;
  static Int_t fNRowLow;
  static Int_t fNRowUp;
  static Int_t fNRowUp1;
  static Int_t fNRowUp2;
  static Int_t fNSectorLow;
  static Int_t fNSectorUp;
  static Double_t fPadPitchWidthLow;
  static Double_t fPadPitchWidthUp;
  static Double_t fZWidth;
  static Double_t fZSigma;
  static Int_t fNSector;
  static Int_t fNSlice;
  static Int_t fNRow;
  static Double_t fPi;
  static Double_t fNRotShift;
  static Double_t fZLength;
  static Double_t fZOffset;
  static Int_t fNPads[159]; //fill this following Init
  static Double_t fX[159];  //fill this following Init
  static Int_t fVersion;    //flags the version
  static Int_t fRows[6][2];
  static Int_t fNRows[6];
  static Int_t fNPatches;
  static Double_t fDiffT; //Transversal diffusion constant
  static Double_t fDiffL; //Longitudinal diffusion constant
  static Double_t fAnodeWireSpacing; 
  static Double_t fInnerPadLength;
  static Double_t fOuter1PadLength;
  static Double_t fOuter2PadLength;
  static Double_t fInnerPRFSigma;
  static Double_t fOuter1PRFSigma;
  static Double_t fOuter2PRFSigma;
  static Double_t fTimeSigma; //Minimal longitudinal width
  static Int_t fADCSat;
  
 public:

  static Bool_t Init(Char_t* path,Bool_t UseAliTPCParam=kFALSE); //init for different AliRoot versions
  static Bool_t MakeInitFile(Char_t *filename,Char_t *path);     //create the init file in path
  static Bool_t ReadInit(Char_t *path);                          //read it possibly from root file 

  static void SetBField(Double_t f) {fBField = f;}

  static const Char_t* GetParamName() {return "75x40_100x60_150x60";}

  static Int_t GetFirstRow(Int_t patch);
  static Int_t GetLastRow(Int_t patch);
  static Int_t GetNRows(Int_t patch);
  static Int_t GetNRows() {return fNRow;}
  static Int_t GetNPatches() {return fNPatches;}
  static Double_t GetBField() {return fBField;}
  static Double_t GetBFact() {return fBFACT;}
  static Double_t GetBFieldValue() {return (fBField*fBFACT);}
  static Float_t Deg2Rad(Float_t angle) {return angle*fPi/180;}
  static Double_t Pi() {return fPi;}
  static Int_t GetVersion(){return fVersion;}
  static Double_t GetPadPitchWidthLow() {return fPadPitchWidthLow;}
  static Double_t GetPadPitchWidthUp() {return fPadPitchWidthUp;}
  static Double_t GetPadPitchWidth(Int_t patch) {return patch < 2 ? fPadPitchWidthLow : fPadPitchWidthUp;}  
  static Double_t GetZWidth() {return fZWidth;}
  static Double_t GetZLength() {return fZLength;}
  static Double_t GetZOffset() {return fZOffset;}
  static Double_t GetDiffT() {return fDiffT;}
  static Double_t GetDiffL() {return fDiffL;}
  static Double_t GetAnodeWireSpacing() {return fAnodeWireSpacing;}
  static Double_t GetPadLength(Int_t padrow);
  static Double_t GetPRFSigma(Int_t padrow);
  static Double_t GetTimeSigma() {return fTimeSigma;}
  static Int_t GetADCSat() {return fADCSat;}
  static Int_t GetNSectorLow() {return fNSectorLow;}
  static Int_t GetNSectorUp() {return fNSectorUp;}
  
  static Bool_t Slice2Sector(Int_t slice, Int_t slicerow, Int_t &sector, Int_t &row);
  static Bool_t Sector2Slice(Int_t &slice, Int_t sector);
  static Bool_t Sector2Slice(Int_t &slice, Int_t &slicerow, Int_t sector, Int_t row);

  static Int_t GetNPads(Int_t row){return (row<fNRow)?fNPads[row]:0;}
  static Int_t GetNTimeBins(){return fNTimeBins;}
  static Double_t Row2X(Int_t slicerow);
  static Double_t GetMaxY(Int_t slicerow);
  static Double_t GetEta(Float_t *xyz);
  static Double_t GetEta(Int_t slice,Int_t padrow, Int_t pad, Int_t time);
  static Double_t GetPhi(Float_t *xyz);

  static void XYZtoRPhiEta(Float_t *rpe, Float_t *xyz);
  static void Local2Global(Float_t *xyz, Int_t slice);
  static void Local2GlobalAngle(Float_t *angle, Int_t slice);
  static void Global2LocalAngle(Float_t *angle, Int_t slice);

  static void Raw2Local(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  static void Local2Global(Float_t *xyz, Int_t sector, Int_t row);
  static void Global2Local(Float_t *xyz, Int_t sector, Bool_t isSlice=kFALSE);
  static void Raw2Global(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  static void Local2Raw(Float_t *xyz, Int_t sector, Int_t row);
  static void Global2Raw(Float_t *xyz, Int_t sector, Int_t row);
  
  ClassDef(AliL3Transform,1)
};
#endif
