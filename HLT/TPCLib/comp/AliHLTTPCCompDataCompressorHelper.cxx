// $Id$

/**************************************************************************
 * TPCCompModelAnalysisright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Anders Vestbo <mailto:vestbo@fi.uib.no>                       *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCCompDataCompressorHelper.cxx
    @author Anders Vestbo
    @date   30-11-2006
    @brief  The Data Compressor helper for the Vestbo-compression for the TPC
*/

#include "AliHLTStdIncludes.h"

#include "AliHLTTPCTransform.h"

#include "AliHLTTPCCompDataCompressorHelper.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTTPCCompDataCompression
//
// Interface class; binary <-> AliROOT handling of TPC data compression classes.
//


ClassImp(AliHLTTPCCompDataCompressorHelper)

AliHLTTPCCompDataCompressorHelper::AliHLTTPCCompDataCompressorHelper()
{
}

// see header file for more documentation
// used variables for compression
Int_t AliHLTTPCCompDataCompressorHelper::fgNumTimeBits = 12;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumPadBits = 12;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumChargeBits = 14;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumShapeBits = 14;
Float_t AliHLTTPCCompDataCompressorHelper::fgXYResidualStep1 = 0.03;
Float_t AliHLTTPCCompDataCompressorHelper::fgXYResidualStep2 = 0.03;
Float_t AliHLTTPCCompDataCompressorHelper::fgXYResidualStep3 = 0.03;
Float_t AliHLTTPCCompDataCompressorHelper::fgZResidualStep1 = 0.05;
Float_t AliHLTTPCCompDataCompressorHelper::fgZResidualStep2 = 0.05;
Float_t AliHLTTPCCompDataCompressorHelper::fgZResidualStep3 = 0.05;
Float_t AliHLTTPCCompDataCompressorHelper::fgXYWidthStep = 0.005;
Float_t AliHLTTPCCompDataCompressorHelper::fgZWidthStep = 0.005;
Int_t AliHLTTPCCompDataCompressorHelper::fgClusterCharge = 100;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumPadBitsRemaining = 18;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumTimeBitsRemaining = 19;
Int_t AliHLTTPCCompDataCompressorHelper::fgNumShapeBitsRemaining = 11;

void AliHLTTPCCompDataCompressorHelper::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  // see header file for class documentation
  // sets the numbers of bits
  fgNumPadBits = pad;
  fgNumTimeBits = time;
  fgNumChargeBits = charge;
  fgNumShapeBits = shape;
}

void AliHLTTPCCompDataCompressorHelper::SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // see header file for class documentation
  // sets the transverse resolution
  fgXYResidualStep1 = res1;
  fgXYResidualStep2 = res2;
  fgXYResidualStep3 = res3;
  fgXYWidthStep = width;
}

void AliHLTTPCCompDataCompressorHelper::SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width)
{
  // see header file for class documentation
  // sets the longitudinal resolution
  fgZResidualStep1 = res1;
  fgZResidualStep2 = res2;
  fgZResidualStep3 = res3;
  fgZWidthStep = width;
}

void AliHLTTPCCompDataCompressorHelper::SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape)
{
  // see header file for class documentation
  // sets the numbers of remaining bits
  fgNumPadBitsRemaining = pad;
  fgNumTimeBitsRemaining = time;
  fgNumShapeBitsRemaining = shape;
}

Float_t AliHLTTPCCompDataCompressorHelper::GetXYResidualStep(Int_t row) 
{
  // see header file for class documentation
  // gets the XY residual step
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fgXYResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fgXYResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fgXYResidualStep3;
  else
    {
      cerr<<"AliHLTTPCCompDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTTPCCompDataCompressorHelper::GetZResidualStep(Int_t row) 
{
  // see header file for class documentation 
  // gets the Z residual step
  if(row < AliHLTTPCTransform::GetNRowLow())
    return fgZResidualStep1;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1())
    return fgZResidualStep2;
  else if(row < AliHLTTPCTransform::GetNRowLow() + AliHLTTPCTransform::GetNRowUp1() + AliHLTTPCTransform::GetNRowUp2())
    return fgZResidualStep3;
  else
    {
      cerr<<"AliHLTTPCCompDataCompressorHelper::GetXYResidualStep : Wrong row number "<<row<<endl;
      return -1;
    }
}

Float_t AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor()
{
  // see header file for class documentation
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
      cerr<<"AliHLTTPCCompDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}

Float_t AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor()
{
  // see header file for class documentation
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
      cerr<<"AliHLTTPCCompDataCompressorHelper::GetRemainingPadFactor : Too few bits for the pad direction: "<<nbits<<endl;
      return 1;
    }
}


Int_t AliHLTTPCCompDataCompressorHelper::Nint(Double_t x)
{
  // see header file for class documentation 
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
