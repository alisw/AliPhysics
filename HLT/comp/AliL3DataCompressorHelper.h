// @(#) $Id$

#ifndef AliL3_DataCompressorHelper
#define AliL3_DataCompressorHelper

#include "AliL3RootTypes.h"

class AliL3DataCompressorHelper {
  
 private:
  static Int_t fNumPadBits;
  static Int_t fNumTimeBits;
  static Int_t fNumChargeBits;
  static Int_t fNumShapeBits;
  
  static Float_t fXYResidualStep1;
  static Float_t fXYResidualStep2;
  static Float_t fXYResidualStep3;
  static Float_t fZResidualStep1;
  static Float_t fZResidualStep2;
  static Float_t fZResidualStep3;
  static Float_t fXYWidthStep;
  static Float_t fZWidthStep;
  static Int_t fClusterCharge;
  
 protected:
  
 public:
  static void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape);
  static void SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  static void SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  static const Int_t GetNPadBits() {return fNumPadBits;}
  static const Int_t GetNTimeBits() {return fNumTimeBits;}
  static const Int_t GetNChargeBits() {return fNumChargeBits;}
  static const Int_t GetNShapeBits() {return fNumShapeBits;}
  static const Float_t GetXYWidthStep() {return fXYWidthStep;}
  static const Float_t GetZWidthStep() {return fZWidthStep;}
  static const Int_t GetClusterCharge() {return fClusterCharge;}
  static const Float_t GetXYResidualStep(Int_t row);
  static const Float_t GetZResidualStep(Int_t row);

  ClassDef(AliL3DataCompressorHelper,1) 

};

#endif
