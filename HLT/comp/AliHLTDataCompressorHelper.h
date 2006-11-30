// @(#) $Id$

#ifndef AliHLT_DataCompressorHelper
#define AliHLT_DataCompressorHelper

#include "AliHLTRootTypes.h"

class AliHLTDataCompressorHelper {
  
 public:
  virtual ~AliHLTDataCompressorHelper() {}

  static void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape);
  static void SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  static void SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  static void SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape);
  static Int_t GetNPadBits() {return fgNumPadBits;}
  static Int_t GetNTimeBits() {return fgNumTimeBits;}
  static Int_t GetNChargeBits() {return fgNumChargeBits;}
  static Int_t GetNShapeBits() {return fgNumShapeBits;}
  static Float_t GetXYWidthStep() {return fgXYWidthStep;}
  static Float_t GetZWidthStep() {return fgZWidthStep;}
  static Int_t GetClusterCharge() {return fgClusterCharge;}
  static Float_t GetXYResidualStep(Int_t row);
  static Float_t GetZResidualStep(Int_t row);
  static Int_t GetNPadBitsRemaining() {return fgNumPadBitsRemaining;}
  static Int_t GetNTimeBitsRemaining() {return fgNumTimeBitsRemaining;}
  static Int_t GetNShapeBitsRemaining() {return fgNumShapeBitsRemaining;}
  static Float_t GetPadPrecisionFactor();
  static Float_t GetTimePrecisionFactor();

  //taken from TMath
  static Int_t Nint(Double_t x); 
  static Int_t Abs(Int_t d) { return (d > 0) ? d : -d; }
  static Double_t Abs(Double_t d) { return (d > 0) ? d : -d; }

 private:
  static Int_t fgNumPadBits; // Number of pad bits
  static Int_t fgNumTimeBits; // Number of time bits
  static Int_t fgNumChargeBits; // Number of charge bits
  static Int_t fgNumShapeBits; // Number of shape bits
  static Int_t fgNumPadBitsRemaining; // Number of remaining pad bits
  static Int_t fgNumTimeBitsRemaining; // Number of remaining time bits
  static Int_t fgNumShapeBitsRemaining; // Number of remaining shape bits

  static Float_t fgXYResidualStep1; // XY resbual at step 1
  static Float_t fgXYResidualStep2; // XY residual at step 2
  static Float_t fgXYResidualStep3; // XY resudual at step 3
  static Float_t fgZResidualStep1; // Z residual at step 1
  static Float_t fgZResidualStep2; // Z resudual at step 2
  static Float_t fgZResidualStep3; // Z resudual at step 3
  static Float_t fgXYWidthStep; // Width of XY step
  static Float_t fgZWidthStep;  // Width of Z step
  static Int_t fgClusterCharge; // Cluster charge


  ClassDef(AliHLTDataCompressorHelper,1) 

};

typedef AliHLTDataCompressorHelper AliL3DataCompressorHelper; // for backward compatibility

#endif
