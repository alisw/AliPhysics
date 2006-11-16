// @(#) $Id$
// Original: AliL3DataCompressorHelper.cxx,v 1.5 2004/06/15 10:26:57 hristov Exp $

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStdIncludes.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCTransform.h"

#include "AliHLTTPCDataCompressorHelper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTTPCDataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliHLTTPCDataCompressorHelper)

AliHLTTPCDataCompressorHelper::AliHLTTPCDataCompressorHelper() 
{
}

AliHLTTPCDataCompressorHelper::~AliHLTTPCDataCompressorHelper() 
{
}

Int_t AliHLTTPCDataCompressorHelper::fgNumTimeBits = 12;
Int_t AliHLTTPCDataCompressorHelper::fgNumPadBits = 12;
Int_t AliHLTTPCDataCompressorHelper::fgNumChargeBits = 14;
Int_t AliHLTTPCDataCompressorHelper::fgNumShapeBits = 14;
Float_t AliHLTTPCDataCompressorHelper::fgXYResidualStep1 = 0.03;
Float_t AliHLTTPCDataCompressorHelper::fgXYResidualStep2 = 0.03;
Float_t AliHLTTPCDataCompressorHelper::fgXYResidualStep3 = 0.03;
Float_t AliHLTTPCDataCompressorHelper::fgZResidualStep1 = 0.05;
Float_t AliHLTTPCDataCompressorHelper::fgZResidualStep2 = 0.05;
Float_t AliHLTTPCDataCompressorHelper::fgZResidualStep3 = 0.05;
Float_t AliHLTTPCDataCompressorHelper::fgXYWidthStep = 0.005;
Float_t AliHLTTPCDataCompressorHelper::fgZWidthStep = 0.005;
Int_t AliHLTTPCDataCompressorHelper::fgClusterCharge = 100;
Int_t AliHLTTPCDataCompressorHelper::fgNumPadBitsRemaining = 18;
Int_t AliHLTTPCDataCompressorHelper::fgNumTimeBitsRemaining = 19;
Int_t AliHLTTPCDataCompressorHelper::fgNumShapeBitsRemaining = 11;

void AliHLTTPCDataCompressorHelper::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  // sets the numbers of bits
  fgNumPadBits = pad;
  fgNumTimeBits = time;
  fgNumChargeBits = charge;
  fgNumShapeBits = shape;
}

void AliHLTTPCDataCompressorHelper::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the transverse resolution
  fgXYResidualStep1 = res1;
  fgXYResidualStep2 = res2;
  fgXYResidualStep3 = res3;
  fgXYWidthStep = width;
}

void AliHLTTPCDataCompressorHelper::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // sets the longitudinal resolution
  fgZResidualStep1 = res1;
  fgZResidualStep2 = res2;
  fgZResidualStep3 = res3;
  fgZWidthStep = width;
}

void AliHLTTPCDataCompressorHelper::SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape)
{
  // sets the numbers of remaining bits
  fgNumPadBitsRemaining = pad;
  fgNumTimeBitsRemaining = time;
  fgNumShapeBitsRemaining = shape;
}

Float_t AliHLTTPCDataCompressorHelper::GetXYResidualStep(Int_t row) 
{
  // gets the XY residual step
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fgXYResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fgXYResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fgXYResidualStep3;
  else
    {
      cerr<<"AliHLTTPCDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTTPCDataCompressorHelper::GetZResidualStep(Int_t row) 
{
  // gets the Z residual step
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fgZResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fgZResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fgZResidualStep3;
  else
    {
      cerr<<"AliHLTTPCDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTTPCDataCompressorHelper::GetPadPrecisionFactor()
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
      cerr<<"AliHLTTPCDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}

Float_t AliHLTTPCDataCompressorHelper::GetTimePrecisionFactor()
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
      cerr<<"AliHLTTPCDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}


Int_t AliHLTTPCDataCompressorHelper::Nint(Double_t x)
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
