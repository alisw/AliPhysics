// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3RootTypes.h"
#include "AliL3Transform.h"

#include "AliL3DataCompressorHelper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliL3DataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliL3DataCompressorHelper)


Int_t AliL3DataCompressorHelper::fgNumTimeBits = 12;
Int_t AliL3DataCompressorHelper::fgNumPadBits = 12;
Int_t AliL3DataCompressorHelper::fgNumChargeBits = 14;
Int_t AliL3DataCompressorHelper::fgNumShapeBits = 14;
Float_t AliL3DataCompressorHelper::fgXYResidualStep1 = 0.03;
Float_t AliL3DataCompressorHelper::fgXYResidualStep2 = 0.03;
Float_t AliL3DataCompressorHelper::fgXYResidualStep3 = 0.03;
Float_t AliL3DataCompressorHelper::fgZResidualStep1 = 0.05;
Float_t AliL3DataCompressorHelper::fgZResidualStep2 = 0.05;
Float_t AliL3DataCompressorHelper::fgZResidualStep3 = 0.05;
Float_t AliL3DataCompressorHelper::fgXYWidthStep = 0.005;
Float_t AliL3DataCompressorHelper::fgZWidthStep = 0.005;
Int_t AliL3DataCompressorHelper::fgClusterCharge = 100;
Int_t AliL3DataCompressorHelper::fgNumPadBitsRemaining = 18;
Int_t AliL3DataCompressorHelper::fgNumTimeBitsRemaining = 19;
Int_t AliL3DataCompressorHelper::fgNumShapeBitsRemaining = 11;

void AliL3DataCompressorHelper::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  // sets the numbers of bits
  fgNumPadBits = pad;
  fgNumTimeBits = time;
  fgNumChargeBits = charge;
  fgNumShapeBits = shape;
}

void AliL3DataCompressorHelper::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the transverse resolution
  fgXYResidualStep1 = res1;
  fgXYResidualStep2 = res2;
  fgXYResidualStep3 = res3;
  fgXYWidthStep = width;
}

void AliL3DataCompressorHelper::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the longitudinal resolution
  fgZResidualStep1 = res1;
  fgZResidualStep2 = res2;
  fgZResidualStep3 = res3;
  fgZWidthStep = width;
}

void AliL3DataCompressorHelper::SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape)
{
  // sets the numbers of remaining bits
  fgNumPadBitsRemaining = pad;
  fgNumTimeBitsRemaining = time;
  fgNumShapeBitsRemaining = shape;
}

Float_t AliL3DataCompressorHelper::GetXYResidualStep(Int_t row) 
{
  // gets the XY residual step
  if(row < AliL3Transform::GetNRowLow())
    return fgXYResidualStep1;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1())
    return fgXYResidualStep2;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1() + AliL3Transform::GetNRowUp2())
    return fgXYResidualStep3;
  else
    {
      cerr<<"AliL3DataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliL3DataCompressorHelper::GetZResidualStep(Int_t row) 
{
  // gets the Z residual step
  if(row < AliL3Transform::GetNRowLow())
    return fgZResidualStep1;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1())
    return fgZResidualStep2;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1() + AliL3Transform::GetNRowUp2())
    return fgZResidualStep3;
  else
    {
      cerr<<"AliL3DataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliL3DataCompressorHelper::GetPadPrecisionFactor()
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
      cerr<<"AliL3DataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}

Float_t AliL3DataCompressorHelper::GetTimePrecisionFactor()
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
      cerr<<"AliL3DataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}


Int_t AliL3DataCompressorHelper::Nint(Double_t x)
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
