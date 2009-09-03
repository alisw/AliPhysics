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

#include "AliVZEROdigit.h"

ClassImp(AliVZEROdigit)

//__________________________________________________________________________
AliVZEROdigit::AliVZEROdigit()
   :AliDigit(),
    fTrack(0),
    fEvent(0),
    fPMNumber(0),
    fADC(0.),
    fTime(0.),
    fWidth(0.),
    fBBFlag(0),
    fBGFlag(0),
    fIntegrator(0)

{
  // Standard default constructor 
}

//__________________________________________________________________________
AliVZEROdigit::AliVZEROdigit(Int_t* tracks, Int_t *digits)
   :AliDigit(tracks),
   fTrack(0),
   fEvent(0),
   fPMNumber(0),
   fADC(0.),
   fTime(0.),
   fWidth(0.),
   fBBFlag(0),
   fBGFlag(0),
   fIntegrator(0)
  
{
  // Creates VZERO digits  
  // Fills the AliVZEROdigit data members from the array digits. 
    
  fTrack      =  tracks[0];
  fEvent      =  digits[0];  
}

//__________________________________________________________________________
AliVZEROdigit::AliVZEROdigit(Int_t PMnumber, Float_t adc, Float_t time)
   :AliDigit(),
   fTrack(0),
   fEvent(0),
   fPMNumber(0),
   fADC(0.),
   fTime(0.),
   fWidth(0.),
   fBBFlag(0),
   fBGFlag(0),
   fIntegrator(0)
{  
   fPMNumber   = PMnumber;
   fADC        = adc;
   fTime       = time;
}

//__________________________________________________________________________
AliVZEROdigit::AliVZEROdigit(Int_t  PMnumber, Float_t adc, Float_t time, 
                             Float_t width, Bool_t BeamBeamFlag, Bool_t BeamGasFlag)
   :AliDigit(),
   fTrack(0),
   fEvent(0),
   fPMNumber(PMnumber),
   fADC(adc),
   fTime(time),
   fWidth(width),
   fBBFlag(BeamBeamFlag),
   fBGFlag(BeamGasFlag),
   fIntegrator(0)
{  

}
//__________________________________________________________________________
AliVZEROdigit::AliVZEROdigit(Int_t   PMnumber, Float_t adc, Float_t time, 
                             Float_t width, Bool_t BeamBeamFlag, Bool_t BeamGasFlag, Bool_t integrator)
:AliDigit(),
fTrack(0),
fEvent(0),
fPMNumber(PMnumber),
fADC(adc),
fTime(time),
fWidth(width),
fBBFlag(BeamBeamFlag),
fBGFlag(BeamGasFlag),
fIntegrator(integrator)
{  
	
}

//__________________________________________________________________________
void AliVZEROdigit::Print(const Option_t*) const
{
    // Dumps digit object
    
    Dump();
}


