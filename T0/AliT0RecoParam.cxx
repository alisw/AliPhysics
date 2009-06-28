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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with T0 reconstruction parameters                                  //
//    fMeanAmplitude -
//           for low flux time-amplitude correction equalize time to amplitude 1 MIP; 
//           for high flux - to 15MIP   
//    To have nice time spectra after reconstruction we need to know 
//    reference point to write t(i) - RefPoint. 
//    It can be apparatus RefPoint or one of PMT                       //  
//    fRefPoint - number of channel with RF
//
//       Alla.Maevskaya@cern.ch
/////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
 
#include "AliT0RecoParam.h"
#include "Riostream.h"

ClassImp(AliT0RecoParam)




//_____________________________________________________________________________
AliT0RecoParam::AliT0RecoParam():
  AliDetectorRecoParam(),
   fRefAmp(7),
    fRefPoint(0)
{
  //
  // constructor
  //
  SetName("T0");
  SetTitle("T0");
}

//_____________________________________________________________________________
AliT0RecoParam::~AliT0RecoParam() 
{
  //
  // destructor
  //  
}

//_____________________________________________________________________________

AliT0RecoParam::AliT0RecoParam(const AliT0RecoParam &p):
  AliDetectorRecoParam(p),       
   fRefAmp(p.fRefAmp),
   fRefPoint(p.fRefPoint)
{
 
 //copy constructor

}
//_____________________________________________________________________________

AliT0RecoParam& AliT0RecoParam:: operator=(const AliT0RecoParam &p)
{
  //
  // assign. operator
  //

  if (this == &p)
    return *this;
  
  AliDetectorRecoParam::operator=(p);
  fRefAmp = p.fRefAmp;
  fRefPoint = p.fRefPoint;
  return *this;

}
//_____________________________________________________________________________
 
AliT0RecoParam *AliT0RecoParam::GetLowFluxParam()
{
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliT0RecoParam *param = new AliT0RecoParam();
  param->fRefAmp = 1;
  param->fRefPoint = 0;
  param->SetName("Low Flux");
  param->SetTitle("Low Flux");
  return param;
}

//_____________________________________________________________________________

AliT0RecoParam *AliT0RecoParam::GetHighFluxParam()
{
  //
  // make reco parameters for high flux env.
  //

  AliT0RecoParam *param = new AliT0RecoParam();
  param->fRefAmp = 5;
  param->fRefPoint = 0;
  //
  param->SetName("High Flux");
  param->SetTitle("High Flux");
  return param;
}


//_____________________________________________________________________________

AliT0RecoParam *AliT0RecoParam::GetLaserTestParam()
{
  //
  // special setting for laser
  //
  AliT0RecoParam *param = new AliT0RecoParam();
  param->fRefAmp = 1;
  param->fRefPoint = 1;
  //
  param->SetName("Laser Flux");
  param->SetTitle("Laser Flux");
  return param;
}
//_____________________________________________________________________________
void AliT0RecoParam::PrintParameters() const
{
  //
  // Printing of the used T0 reconstruction parameters
  //
  AliInfo(Form(" Reference amplitude for walk corerection : %f", fRefAmp));
  AliInfo(Form(" Reference point in channel  : %i", fRefPoint));
 
}
