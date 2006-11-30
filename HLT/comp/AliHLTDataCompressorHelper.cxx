// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTTransform.h"

#include "AliHLTDataCompressorHelper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTDataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliHLTDataCompressorHelper)


Int_t AliHLTDataCompressorHelper::fgNumTimeBits = 12;
Int_t AliHLTDataCompressorHelper::fgNumPadBits = 12;
Int_t AliHLTDataCompressorHelper::fgNumChargeBits = 14;
Int_t AliHLTDataCompressorHelper::fgNumShapeBits = 14;
Float_t AliHLTDataCompressorHelper::fgXYResidualStep1 = 0.03;
Float_t AliHLTDataCompressorHelper::fgXYResidualStep2 = 0.03;
Float_t AliHLTDataCompressorHelper::fgXYResidualStep3 = 0.03;
Float_t AliHLTDataCompressorHelper::fgZResidualStep1 = 0.05;
Float_t AliHLTDataCompressorHelper::fgZResidualStep2 = 0.05;
Float_t AliHLTDataCompressorHelper::fgZResidualStep3 = 0.05;
Float_t AliHLTDataCompressorHelper::fgXYWidthStep = 0.005;
Float_t AliHLTDataCompressorHelper::fgZWidthStep = 0.005;
Int_t AliHLTDataCompressorHelper::fgClusterCharge = 100;
Int_t AliHLTDataCompressorHelper::fgNumPadBitsRemaining = 18;
Int_t AliHLTDataCompressorHelper::fgNumTimeBitsRemaining = 19;
Int_t AliHLTDataCompressorHelper::fgNumShapeBitsRemaining = 11;

void AliHLTDataCompressorHelper::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  // sets the numbers of bits
  fgNumPadBits = pad;
  fgNumTimeBits = time;
  fgNumChargeBits = charge;
  fgNumShapeBits = shape;
}

void AliHLTDataCompressorHelper::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the transverse resolution
  fgXYResidualStep1 = res1;
  fgXYResidualStep2 = res2;
  fgXYResidualStep3 = res3;
  fgXYWidthStep = width;
}

void AliHLTDataCompressorHelper::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the longitudinal resolution
  fgZResidualStep1 = res1;
  fgZResidualStep2 = res2;
  fgZResidualStep3 = res3;
  fgZWidthStep = width;
}

void AliHLTDataCompressorHelper::SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape)
{
  // sets the numbers of remaining bits
  fgNumPadBitsRemaining = pad;
  fgNumTimeBitsRemaining = time;
  fgNumShapeBitsRemaining = shape;
}

Float_t AliHLTDataCompressorHelper::GetXYResidualStep(Int_t row) 
{
  // gets the XY residual step
  if(row < AliHLTTransform::GetNRowLow())
    return fgXYResidualStep1;
  else if(row < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1())
    return fgXYResidualStep2;
  else if(row < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1() + AliHLTTransform::GetNRowUp2())
    return fgXYResidualStep3;
  else
    {
      cerr<<"AliHLTDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTDataCompressorHelper::GetZResidualStep(Int_t row) 
{
  // gets the Z residual step
  if(row < AliHLTTransform::GetNRowLow())
    return fgZResidualStep1;
  else if(row < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1())
    return fgZResidualStep2;
  else if(row < AliHLTTransform::GetNRowLow() + AliHLTTransform::GetNRowUp1() + AliHLTTransform::GetNRowUp2())
    return fgZResidualStep3;
  else
    {
      cerr<<"AliHLTDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTDataCompressorHelper::GetPadPrecisionFactor()
{
  // gets pad precision factor
  Int_t nbits = fgNumPadBitsRemaining;
  if(nbits >=21)
    return 10000;
  if(nbits >= 18)
    return 1000;
  if(nbits >= 14) 
    return 100;
  if(nbits >= 11)
    return 10;
  if(nbits >= 8)
    return 1;
  else 
    {
      cerr<<"AliHLTDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}

Float_t AliHLTDataCompressorHelper::GetTimePrecisionFactor()
{
  // gest time precision factor
  Int_t nbits = fgNumTimeBitsRemaining;
  if(nbits >=23)
    return 10000;
  if(nbits >= 19)
    return 1000;
  if(nbits >= 16) 
    return 100;
  if(nbits >= 13)
    return 10;
  if(nbits >= 9)
    return 1;
  else 
    {
      cerr<<"AliHLTDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}


Int_t AliHLTDataCompressorHelper::Nint(Double_t x)
{
   // Round to nearest integer. Rounds half integers 
   // to the nearest even integer.

   Int_t i=0;
   if (x >= 0) {
      i = Int_t(x + 0.5);
      if (x + 0.5 == Double_t(i) && i & 1) i--;
   } else {
      i = Int_t(x - 0.5);
      if (x - 0.5 == Double_t(i) && i & 1) i++;

   }
   return i;
}
