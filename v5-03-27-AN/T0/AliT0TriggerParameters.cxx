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

/* $Id: AliT0TriggerParameters.cxx 28275 2008-08-28 06:09:21Z alla $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0TriggerParameters.h"
#include "AliLog.h"

#include "Riostream.h"

ClassImp(AliT0TriggerParameters)

//________________________________________________________________
  AliT0TriggerParameters::AliT0TriggerParameters():TObject(),
  fSwtPmt(0),        
  fAmpCentr(5), 
  fAmpSemiCentr(1), 
  fTimeWindowLow(0),
  fTimeWindowHigh(1)

{
  //
  // for (Int_t i=0; i<sizeof(fThreshold)/sizeof(*fThreshold); i++) fThreshold[i] = 0;
 for (Int_t i=0; i<24; i++) fThreshold[i] = 0;
}
//_____________________________________________________________________________

AliT0TriggerParameters::AliT0TriggerParameters(const AliT0TriggerParameters &r):
  TObject(),       
  fSwtPmt(0),        
  fAmpCentr(5), 
  fAmpSemiCentr(1), 
  fTimeWindowLow(0),
  fTimeWindowHigh(1)
{
 
 //copy constructor
  ((AliT0TriggerParameters &) r).Copy(*this);
 
}

//_____________________________________________________________________________

AliT0TriggerParameters& AliT0TriggerParameters:: operator=(const AliT0TriggerParameters &p)
{
  //
  // assign. operator
  //

  if (this == &p)
    return *this;
  AliT0TriggerParameters:: operator=(p);
  fSwtPmt = p.fSwtPmt;
  fAmpCentr = p.fAmpCentr;
  fAmpSemiCentr = p.fAmpSemiCentr;
  fTimeWindowLow = p.fTimeWindowLow;
  fTimeWindowHigh = p.fTimeWindowHigh;
  return *this;
 
}

//________________________________________________________________
AliT0TriggerParameters::~AliT0TriggerParameters()
{
  //
  // destrictor
}
//________________________________________________________________
void AliT0TriggerParameters::Reset()
{
  //reset values

  memset(fThreshold,0,24*sizeof(Int_t));
 
}


//________________________________________________________________
void  AliT0TriggerParameters::Print(Option_t*) const
{
  // print time values

  printf("\n	----	Threshold	----\n\n");
  printf(" Switched on/off\n");
  for (Int_t i=0; i<24; i++) {
    printf(" Threshold  %i status %i ", fThreshold[i], GetPMTstatus(i));
  }
  AliInfo(Form(" Time window around vertex : %f %f",fTimeWindowLow, fTimeWindowHigh ));
  AliInfo(Form(" Amplitude threshold: central  %i semi-central %i", fAmpCentr,fAmpSemiCentr));

} 


//________________________________________________________________
void AliT0TriggerParameters::SetPMTstatus(Int_t i, Int_t val)
{
  if(val)fSwtPmt |= 1<<i;
  else fSwtPmt &= ~(1<<i);

}
//________________________________________________________________
Int_t AliT0TriggerParameters::GetPMTstatus(Int_t i) const
{

  return (1<<i)&fSwtPmt;

}


