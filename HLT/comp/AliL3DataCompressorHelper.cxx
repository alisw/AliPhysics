// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3RootTypes.h"
#include "AliL3Transform.h"

#include "AliL3DataCompressorHelper.h"

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliL3DataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliL3DataCompressorHelper)


Int_t AliL3DataCompressorHelper::fNumTimeBits = 12;
Int_t AliL3DataCompressorHelper::fNumPadBits = 12;
Int_t AliL3DataCompressorHelper::fNumChargeBits = 14;
Int_t AliL3DataCompressorHelper::fNumShapeBits = 14;
Float_t AliL3DataCompressorHelper::fXYResidualStep1 = 0.03;
Float_t AliL3DataCompressorHelper::fXYResidualStep2 = 0.03;
Float_t AliL3DataCompressorHelper::fXYResidualStep3 = 0.03;
Float_t AliL3DataCompressorHelper::fZResidualStep1 = 0.05;
Float_t AliL3DataCompressorHelper::fZResidualStep2 = 0.05;
Float_t AliL3DataCompressorHelper::fZResidualStep3 = 0.05;
Float_t AliL3DataCompressorHelper::fXYWidthStep = 0.005;
Float_t AliL3DataCompressorHelper::fZWidthStep = 0.005;
Int_t AliL3DataCompressorHelper::fClusterCharge = 100;
Int_t AliL3DataCompressorHelper::fNumPadBitsRemaining = 18;
Int_t AliL3DataCompressorHelper::fNumTimeBitsRemaining = 19;
Int_t AliL3DataCompressorHelper::fNumShapeBitsRemaining = 11;

void AliL3DataCompressorHelper::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  fNumPadBits = pad;
  fNumTimeBits = time;
  fNumChargeBits = charge;
  fNumShapeBits = shape;
}

void AliL3DataCompressorHelper::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  fXYResidualStep1 = res1;
  fXYResidualStep2 = res2;
  fXYResidualStep3 = res3;
  fXYWidthStep = width;
}

void AliL3DataCompressorHelper::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  fZResidualStep1 = res1;
  fZResidualStep2 = res2;
  fZResidualStep3 = res3;
  fZWidthStep = width;
}

void AliL3DataCompressorHelper::SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape)
{
  fNumPadBitsRemaining = pad;
  fNumTimeBitsRemaining = time;
  fNumShapeBitsRemaining = shape;
}

const Float_t AliL3DataCompressorHelper::GetXYResidualStep(Int_t row) 
{
  if(row < AliL3Transform::GetNRowLow())
    return fXYResidualStep1;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1())
    return fXYResidualStep2;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1() + AliL3Transform::GetNRowUp2())
    return fXYResidualStep3;
  else
    {
      cerr<<"AliL3DataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

const Float_t AliL3DataCompressorHelper::GetZResidualStep(Int_t row) 
{
  if(row < AliL3Transform::GetNRowLow())
    return fZResidualStep1;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1())
    return fZResidualStep2;
  else if(row < AliL3Transform::GetNRowLow() + AliL3Transform::GetNRowUp1() + AliL3Transform::GetNRowUp2())
    return fZResidualStep3;
  else
    {
      cerr<<"AliL3DataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

const Float_t AliL3DataCompressorHelper::GetPadPrecisionFactor()
{
  Int_t nbits = fNumPadBitsRemaining;
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

const Float_t AliL3DataCompressorHelper::GetTimePrecisionFactor()
{
  Int_t nbits = fNumTimeBitsRemaining;
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
