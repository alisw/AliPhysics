/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONPadInfo
///
/// Class to summarize ESD data at pad
///
/// \author Philippe Pillot, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONPadInfo.h"

#include "AliLog.h"

#include <Riostream.h>

using std::cout;
using std::endl;
/// \cond CLASSIMP
ClassImp(AliMUONPadInfo)
/// \endcond

//_____________________________________________________________________________
AliMUONPadInfo::AliMUONPadInfo()
: TObject(),
  fPadId(0),
  fPadPlaneType(0),
  fPadX(0.),
  fPadY(0.),
  fPadDimX(0.),
  fPadDimY(0.),
  fPadCharge(0.),
  fPadADC(0),
  fPadSaturated(0),
  fPadCalibrated(0),
  fPedMean(0.),
  fPedSigma(0.)
{
  /// default constructor
}

//_____________________________________________________________________________
AliMUONPadInfo::AliMUONPadInfo (const AliMUONPadInfo& padInfo)
: TObject(padInfo),
  fPadId(padInfo.fPadId),
  fPadPlaneType(padInfo.fPadPlaneType),
  fPadX(padInfo.fPadX),
  fPadY(padInfo.fPadY),
  fPadDimX(padInfo.fPadDimX),
  fPadDimY(padInfo.fPadDimY),
  fPadCharge(padInfo.fPadCharge),
  fPadADC(padInfo.fPadADC),
  fPadSaturated(padInfo.fPadSaturated),
  fPadCalibrated(padInfo.fPadCalibrated),
  fPedMean(padInfo.fPedMean),
  fPedSigma(padInfo.fPedSigma)
{
  /// Copy constructor
}

//_____________________________________________________________________________
AliMUONPadInfo& AliMUONPadInfo::operator=(const AliMUONPadInfo& padInfo)
{
  /// Equal operator
  if (this == &padInfo) return *this;
  
  TObject::operator=(padInfo); // don't forget to invoke the base class' assignment operator
  
  fPadId = padInfo.fPadId;
  fPadPlaneType = padInfo.fPadPlaneType;
  fPadX = padInfo.fPadX;
  fPadY = padInfo.fPadY;
  fPadDimX = padInfo.fPadDimX;
  fPadDimY = padInfo.fPadDimY;
  fPadCharge = padInfo.fPadCharge;
  fPadADC = padInfo.fPadADC;
  fPadSaturated = padInfo.fPadSaturated;
  fPadCalibrated = padInfo.fPadCalibrated;
  fPedMean = padInfo.fPedMean;
  fPedSigma = padInfo.fPedSigma;
  
  return *this;
}

//__________________________________________________________________________
AliMUONPadInfo::~AliMUONPadInfo()
{
  /// Destructor
}

//_____________________________________________________________________________
void AliMUONPadInfo::Print(Option_t* option) const
{
  /// print pad info content
  /// also print calibration parameters if option=FULL
    
  cout<<Form("- padID=%u (det=%d, manuI=%d, manuC=%d, cath=%d)",
	     GetPadId(), GetDetElemId(), GetManuId(), GetManuChannel(), GetCathode())<<endl;
  
  cout<<Form("    position=(%5.2f, %5.2f), dimension=(%5.2f, %5.2f)",
	     GetPadX(), GetPadY(), GetPadDimX(), GetPadDimY())<<endl;
  
  cout<<Form("    charge=%5.2f, ADC=%d", GetPadCharge(), GetPadADC())<<endl;
    
  if (strstr(option,"FULL")) {
    cout<<Form("    pedestal (mean=%5.2f, sigma=%5.2f)", GetPedMean(), GetPedSigma())<<endl;
  }
  
}

