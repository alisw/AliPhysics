// @(#) $Id$

#ifndef ALIL3TRANSFORM_H
#define ALIL3TRANSFORM_H

#ifdef use_aliroot
  class AliRunLoader;
#endif

#include "AliL3RootTypes.h"

class AliL3Transform {

 public:
  enum VersionType { fV_default=0, fV_deprecated=1, fV_aliroot=10, fV_cosmics=100};

 private:
  static const Double_t fBFACT;
  static const Double_t fPi;
  static const Double_t fPi2;
  static const Double_t f2Pi;
  static const Double_t fAnodeWireSpacing; 
  static const Double_t fToDeg;

  static Int_t fNPatches; //6 (dont change this) 
  static Int_t fRows[6][2];
  static Int_t fNRows[6];

  static Double_t fBField;
  static Double_t fBFieldFactor;
  static Double_t fSolenoidBField;
  static Int_t fNTimeBins;
  static Int_t fNRowLow;
  static Int_t fNRowUp;
  static Int_t fNRowUp1;
  static Int_t fNRowUp2;
  static Int_t fNSectorLow;
  static Int_t fNSectorUp;
  static Int_t fSlice2Sector[36][2];
  static Int_t fSector2Slice[72];
  static Int_t fSectorLow[72];
  static Double_t fPadPitchWidthLow;
  static Double_t fPadPitchWidthUp;
  static Double_t fZWidth;
  static Double_t fZSigma;
  static Double_t fZLength;
  static Double_t fZOffset;
  static Int_t fNSector; //72  (dont change this)
  static Int_t fNSlice;  //36  (dont change this)
  static Int_t fNRow;    //159 (dont change this)
  static Double_t fNRotShift; //Rotation shift (eg. 0.5 for 10 degrees)
  static Int_t fNPads[159]; //fill this following Init and fVersion
  static Double_t fX[159];  //X position in local coordinates
  static Int_t fVersion;  //flags the version
  static Double_t fDiffT; //Transversal diffusion constant
  static Double_t fDiffL; //Longitudinal diffusion constant
  static Double_t fOmegaTau; //ExB effects
  static Double_t fInnerPadLength;
  static Double_t fOuter1PadLength;
  static Double_t fOuter2PadLength;
  static Double_t fInnerPRFSigma;
  static Double_t fOuter1PRFSigma;
  static Double_t fOuter2PRFSigma;
  static Double_t fTimeSigma; //Minimal longitudinal width
  static Int_t fADCSat; //ADC Saturation (1024 = 10 bit)
  static Int_t fZeroSup; //Zero suppression threshold
  static Double_t fCos[36]; //stores the cos value for local to global rotations  
  static Double_t fSin[36]; //stores the sin value for local to global rotations  

 public:
#ifdef use_aliroot
  static Bool_t Init(AliRunLoader *runLoader); //init transformer params using a run loader
#endif
  static Bool_t Init(Char_t* path,Bool_t UseAliTPCParam=kFALSE); //init transformer settings (versions)
  static Bool_t MakeInitFile(Char_t *rootfilename,Char_t *filename); //create the init file from rootfile
  static Bool_t ReadInit(Char_t *path);         //read init (possibly from root file)
  static Bool_t ReadInitFile(Char_t *path);     //read init from text file 
  static Bool_t SaveInitFile(Char_t *filename); //save parameters in init file

  //setters
  static void SetNPatches(Int_t i){fNPatches = i;}
  static void SetNRows(Int_t s[6]){for(Int_t i=0;i<fNPatches;i++) fNRows[i] = s[i];}
  static void SetRows(Int_t s[6][2]){
    for(Int_t i=0;i<fNPatches;i++){
      fRows[i][0] = s[i][0];
      fRows[i][1] = s[i][1];
    }
  }
  static void SetBField(Double_t f) {fBField = f;} //careful, these 3 are not independent!
  static void SetBFieldFactor(Double_t f) {fBFieldFactor = f;fBField=fBFieldFactor*fSolenoidBField*0.1;}
  static void SetSolenoidBField(Double_t f){fSolenoidBField = f;fBField=fBFieldFactor*fSolenoidBField*0.1;}
  static void SetNTimeBins(Int_t i){fNTimeBins = i;}
  static void SetNRowLow(Int_t i){fNRowLow = i;}
  static void SetNRowUp(Int_t i){fNRowUp = i;}
  static void SetNRowUp1(Int_t i){fNRowUp1 = i;}
  static void SetNRowUp2(Int_t i){fNRowUp2 = i;}
  static void SetSlice2Sector(Int_t s[36][2]){
    for(Int_t i=0;i<fNSlice;i++){
      fSlice2Sector[i][0] = s[i][0];
      fSlice2Sector[i][1] = s[i][1];
    }
  }
  static void SetSector2Slice(Int_t s[72]){for(Int_t i=0;i<fNSector;i++) fSector2Slice[i] = s[i];}
  static void SetSectorLow(Int_t s[72]){for(Int_t i=0;i<fNSector;i++) fSectorLow[i] = s[i];}
  static void SetNSectorLow(Int_t i){fNSectorLow = i;}
  static void SetNSectorUp(Int_t i){fNSectorUp = i;}
  static void SetPadPitchWidthLow(Double_t f){fPadPitchWidthLow = f;}
  static void SetPadPitchWidthUp(Double_t f){fPadPitchWidthUp = f;}
  static void SetZWidth(Double_t f){fZWidth = f;}
  static void SetZSigma(Double_t f){fZSigma = f;}
  static void SetZLength(Double_t f){fZLength = f;}
  static void SetZOffset(Double_t f){fZOffset = f;}
  static void SetNSector(Int_t i){fNSector = i;}
  static void SetNSlice(Int_t i){fNSlice = i;}
  static void SetNRow(Int_t i){fNRow = i;}
  static void SetNRotShift(Double_t f){fNRotShift = f;}
  static void SetNPads(Int_t pads[159]){for(Int_t i=0;i<fNRow;i++) fNPads[i] = pads[i];}
  static void SetX(Double_t xs[159]){for(Int_t i=0;i<fNRow;i++) fX[i] = xs[i];}
  static void SetVersion(Int_t i){fVersion = i;}
  static void SetDiffT(Double_t f){fDiffT = f;}
  static void SetDiffL(Double_t f){fDiffL = f;}
  static void SetOmegaTau(Double_t f){fOmegaTau = f;}
  static void SetInnerPadLength(Double_t f){fInnerPadLength = f;}
  static void SetOuter1PadLength(Double_t f){fOuter1PadLength = f;}
  static void SetOuter2PadLength(Double_t f){fOuter2PadLength = f;}
  static void SetInnerPRFSigma(Double_t f){fInnerPRFSigma = f;}
  static void SetOuter1PRFSigma(Double_t f){fOuter1PRFSigma = f;}
  static void SetOuter2PRFSigma(Double_t f){fOuter2PRFSigma = f;}
  static void SetTimeSigma(Double_t f){fTimeSigma = f;}
  static void SetADCSat(Int_t i) {fADCSat = i;}
  static void SetZeroSup(Int_t i) {fZeroSup = i;}

  //getters
  static const Char_t* GetParamName() {return "75x40_100x60_150x60";}
  static const Double_t Pi()     {return fPi;}
  static const Double_t PiHalf() {return fPi2;}
  static const Double_t TwoPi()  {return f2Pi;}
  static const Double_t GetAnodeWireSpacing() {return fAnodeWireSpacing;}
  static const Double_t GetBFact() {return fBFACT;}
  static const Double_t ToRad() {return 1./fToDeg;}
  static const Double_t ToDeg() {return fToDeg;}

  static Int_t GetFirstRow(Int_t patch);
  static Int_t GetLastRow(Int_t patch);
  static Int_t GetNRows(Int_t patch);
  static Int_t GetPatch(Int_t padrow);
  static Int_t GetNRows() {return fNRow;}
  static Int_t GetNRowLow() {return fNRowLow;}
  static Int_t GetNRowUp1() {return fNRowUp1;}
  static Int_t GetNRowUp2() {return fNRowUp2;}
  static Int_t GetPadRow(Float_t x);
  static Int_t GetNPatches() {return fNPatches;}
  static Int_t GetNPads(Int_t row);
  static Int_t GetNTimeBins(){return fNTimeBins;}
  static Double_t GetBField() {return fBField;}
  static Double_t GetSolenoidField() {return fSolenoidBField;}
  static Double_t GetBFactFactor() {return fBFieldFactor;}
  static Double_t GetBFieldValue() {return (fBField*fBFACT);}
  static Float_t Deg2Rad(Float_t angle) {return angle/fToDeg;}
  static Float_t Rad2Deg(Float_t angle) {return angle*fToDeg;}
  static Int_t GetVersion(){return fVersion;}
  static Double_t GetPadPitchWidthLow() {return fPadPitchWidthLow;}
  static Double_t GetPadPitchWidthUp() {return fPadPitchWidthUp;}
  static Double_t GetPadPitchWidth(Int_t patch);
  static Double_t GetZWidth() {return fZWidth;}
  static Double_t GetZLength() {return fZLength;}
  static Double_t GetZOffset() {return fZOffset;}
  static Double_t GetDiffT() {return fDiffT;}
  static Double_t GetDiffL() {return fDiffL;}
  static Double_t GetParSigmaY2(Int_t padrow,Float_t z,Float_t angle);
  static Double_t GetParSigmaZ2(Int_t padrow,Float_t z,Float_t tgl);
  static Double_t GetOmegaTau() {return fOmegaTau;}
  static Double_t GetPadLength(Int_t padrow);
  static Double_t GetPRFSigma(Int_t padrow);
  static Double_t GetTimeSigma() {return fTimeSigma;}
  static Double_t GetZSigma() {return fZSigma;}
  static Int_t GetADCSat() {return fADCSat;}
  static Int_t GetZeroSup() {return fZeroSup;}
  static Int_t GetNSlice() {return fNSlice;}
  static Int_t GetNSector() {return fNSector;}
  static Int_t GetNSectorLow() {return fNSectorLow;}
  static Int_t GetNSectorUp() {return fNSectorUp;}
  
  static Bool_t Slice2Sector(Int_t slice, Int_t slicerow, Int_t &sector, Int_t &row);
  static Bool_t Sector2Slice(Int_t &slice, Int_t sector);
  static Bool_t Sector2Slice(Int_t &slice, Int_t &slicerow, Int_t sector, Int_t row);

  static Double_t Row2X(Int_t slicerow);
  static Double_t GetMaxY(Int_t slicerow);
  static Double_t GetEta(Float_t *xyz);
  static Double_t GetEta(Int_t slice,Int_t padrow, Int_t pad, Int_t time);
  static Double_t GetPhi(Float_t *xyz);
  static Double_t GetZFast(Int_t slice, Int_t time, Float_t vertex=0.);

  static void XYZtoRPhiEta(Float_t *rpe, Float_t *xyz);
  static void Local2Global(Float_t *xyz, Int_t slice);
  static void Local2GlobalAngle(Float_t *angle, Int_t slice);
  static void Global2LocalAngle(Float_t *angle, Int_t slice);

  //we have 3 different system: Raw   : row, pad, time
  //                            Local : x,y and global z
  //                            Global: global x,y and global z
  //the methods with HLT in the name differ from the other
  //as you specify slice and slicerow, instead of sector
  //and sector row. In that way we safe "a few ifs"
  static void Raw2Local(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  static void RawHLT2Local(Float_t *xyz,Int_t slice,Int_t slicerow,Float_t pad,Float_t time);
  static void Raw2Local(Float_t *xyz, Int_t sector, Int_t row, Int_t pad, Int_t time);
  static void RawHLT2Local(Float_t *xyz,Int_t slice,Int_t slicerow,Int_t pad,Int_t time);
  static void Local2Global(Float_t *xyz, Int_t sector, Int_t row);
  static void LocHLT2Global(Float_t *xyz, Int_t slice, Int_t slicerow);
  static void Global2Local(Float_t *xyz, Int_t sector);
  static void Global2LocHLT(Float_t *xyz, Int_t slice);
  static void Raw2Global(Float_t *xyz, Int_t sector, Int_t row, Float_t pad, Float_t time);
  static void RawHLT2Global(Float_t *xyz, Int_t slice, 
                            Int_t slicerow, Float_t pad, Float_t time);
  static void Raw2Global(Float_t *xyz, Int_t sector, Int_t row, Int_t pad, Int_t time);
  static void RawHLT2Global(Float_t *xyz, Int_t slice, 
                            Int_t slicerow, Int_t pad, Int_t time);
  static void Local2Raw(Float_t *xyz, Int_t sector, Int_t row);
  static void LocHLT2Raw(Float_t *xyz, Int_t slice, Int_t slicerow);
  static void Global2Raw(Float_t *xyz, Int_t sector, Int_t row);
  static void Global2HLT(Float_t *xyz, Int_t slice, Int_t slicerow);

  static void PrintCompileOptions();
  
  ClassDef(AliL3Transform,1)
};
#endif
