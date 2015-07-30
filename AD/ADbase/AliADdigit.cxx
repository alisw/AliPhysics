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

#include "AliADdigit.h"
#include "AliADConst.h"

ClassImp(AliADdigit)

//__________________________________________________________________________
AliADdigit::AliADdigit()
   :AliDigit(),
    fPMNumber(0),
    fTime(0.),
    fWidth(0.),
    fIntegrator(0),
    fBBflag(0),
    fBGflag(0)

{
  // Standard default
  // constructor 
  for(Int_t iClock = 0; iClock < kADNClocks; ++iClock) fChargeADC[iClock] = 0;
}

//__________________________________________________________________________
AliADdigit::AliADdigit(Int_t   PMnumber, Float_t time, 
                             Float_t width,
			     Bool_t integrator,
			     Short_t *chargeADC,
			     Bool_t  BBflag,
		  	     Bool_t  BGflag,
			     Int_t *labels)
:AliDigit(),
fPMNumber(PMnumber),
fTime(time),
fWidth(width),
fIntegrator(integrator),
fBBflag(BBflag),
fBGflag(BGflag)
{  
  // Constructor
  // Used in the digitizer
  if (chargeADC) {
    for(Int_t iClock = 0; iClock < kADNClocks; ++iClock)
      fChargeADC[iClock] = chargeADC[iClock];
  }
  else {
    for(Int_t iClock = 0; iClock < kADNClocks; ++iClock)
      fChargeADC[iClock] = 0;
  }

  if (labels)
    for(Int_t iTrack = 0; iTrack < 3; ++iTrack) fTracks[iTrack] = labels[iTrack];
}
//__________________________________________________________________________
Bool_t AliADdigit::GetIntegratorFlag(Int_t clock)
{
if (clock >= 0 && clock < kADNClocks){
	if(clock%2 == 0) return fIntegrator;
	else return !fIntegrator;
	}
	
else return kFALSE;
}
//__________________________________________________________________________
void AliADdigit::Print(const Option_t*) const
{
    // Dumps digit object
    Dump();
}
