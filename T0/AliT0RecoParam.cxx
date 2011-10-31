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
   fRefPoint(0),
   fLatencyL1(0),  
   fLatencyL1A(0),  
   fLatencyL1C(0),  
   fLatencyHPTDC(0),      
  fVertexShift(0),
  fEqualised(0)
{
  //
  // constructor
  SetName("T0");
  SetTitle("T0");

  fSatelliteThresholds[0] =  -15;
  fSatelliteThresholds[1] =  -1.5;

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
  fRefPoint(p.fRefPoint),
  fLatencyL1(p.fLatencyL1),  
  fLatencyL1A(p.fLatencyL1A),  
  fLatencyL1C(p.fLatencyL1C),  
  fLatencyHPTDC(p.fLatencyHPTDC),      
  fVertexShift(p.fVertexShift), 
  fEqualised( p.fEqualised)
{ 
 //copy constructor
  fSatelliteThresholds[0] = (p.fSatelliteThresholds[0]);
  fSatelliteThresholds[1] = (p.fSatelliteThresholds[1]);

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
  fLatencyL1 = p.fLatencyL1;
  fLatencyL1A = p.fLatencyL1A;
  fLatencyL1C = p.fLatencyL1C;

  fLatencyHPTDC = p.fLatencyHPTDC;
  fVertexShift = p.fVertexShift;

  fSatelliteThresholds[0] = (p.fSatelliteThresholds[0]);
  fSatelliteThresholds[1] = (p.fSatelliteThresholds[1]);
  fEqualised = p.fEqualised;

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
  param->fLatencyL1 = 7782.05;
  param->fLatencyL1A = 7781.90;
  param->fLatencyL1C =  7782.19;
  param->fLatencyHPTDC = 22000;
  param->fVertexShift = 0;
  for (Int_t i=0; i<500; i++)
    {
     param-> fLow[i]=0.;
     param-> fHigh[i]=10000.;
    }
  param->SetName("Low Flux");
  param->SetTitle("Low Flux");
  param->SetSatelliteThresholds(-15, -1.5);
  param->SetEq(0);
  return param;
}

//_____________________________________________________________________________

AliT0RecoParam *AliT0RecoParam::GetHighFluxParam()
{
  //
  // make reco parameters for high flux env.
  //

  AliT0RecoParam *param = new AliT0RecoParam();
  param->fRefAmp = 10;
  param->fRefPoint = 0;
  param->fLatencyL1 = 7782.05;
  param->fLatencyL1A = 7781.90;
  param->fLatencyL1C =  7782.19;
  param->fVertexShift = 0;
  param->fLatencyHPTDC = 22000;
  for (Int_t i=0; i<500; i++)
    {
      param-> fLow[i]=0.;
      param-> fHigh[i]=20000.;
    }
  //
  param->SetSatelliteThresholds(-15, -1.5);
  param->SetEq(0);

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
  param->fLatencyL1 = 7782.05;
  param->fLatencyL1A = 7781.90;
  param->fLatencyL1C =  7782.19;
  param->fLatencyHPTDC = 22000;
  param->fVertexShift = 0;
  param->SetSatelliteThresholds(-15, -1.5);
  param->SetEq(0);
  
  for (Int_t i=0; i<500; i++)
    {
     param-> fLow[i]=0.;
     param-> fHigh[i]=12000.;
    }
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
  AliInfo(Form(" Bad channel in channel  : %i", fRefPoint));
  cout<<" AliT0RecoParam::PrintParameters() "<<endl;
  for (Int_t i=0; i<105; i++) cout<<i<<" "<<fLow[i]<<" "<<fHigh[i]<<endl; 
}
