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

/* */

///////////////////////////////////////////////////
//
// Reconstructed Trigger track
// in
// ALICE
// dimuon
// spectrometer
// note: equivalent to AliMUONTriggerTrack for tracking,
// no need for a AliMUONTriggerTrackParam
///////////////////////////////////////////////////

#include <Riostream.h> // for cout
#include <stdlib.h> // for exit()

#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TObjArray.h>

#include "AliMUONEventReconstructor.h" 
#include "AliMUONTriggerTrack.h"


//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack()
{
    fEventReconstructor = 0;
    fx11 = 0.;
    fy11 = 0.;
    fthetax = 0.;
    fthetay = 0.;
    fGTPattern = 0;
}
//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t thetax, Float_t thetay, Long_t theGTPattern,  AliMUONEventReconstructor* EventReconstructor)
{
    fEventReconstructor = EventReconstructor; // link back to EventReconstructor
    fx11 = x11;
    fy11 = y11;
    fthetax = thetax;
    fthetay = thetay;
    fGTPattern = theGTPattern;

}

//__________________________________________________________________________
AliMUONTriggerTrack::~AliMUONTriggerTrack()
{
  // Destructor
    ;
    
}

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack (const AliMUONTriggerTrack& MUONTriggerTrack):TObject(MUONTriggerTrack)
{
  fEventReconstructor = new AliMUONEventReconstructor(*MUONTriggerTrack.fEventReconstructor);
  fx11 = MUONTriggerTrack.fx11;
  fy11 = MUONTriggerTrack.fy11;
  fthetax = MUONTriggerTrack.fthetax;
  fthetay = MUONTriggerTrack.fthetay;
  fGTPattern = MUONTriggerTrack.fGTPattern;
}
      
//__________________________________________________________________________
AliMUONTriggerTrack & AliMUONTriggerTrack::operator=(const AliMUONTriggerTrack&
MUONTriggerTrack)
{
    if (this == &MUONTriggerTrack)
	return *this;
    
    fEventReconstructor = new AliMUONEventReconstructor(*MUONTriggerTrack.fEventReconstructor);
    fx11 = MUONTriggerTrack.fx11;
    fy11 = MUONTriggerTrack.fy11;
    fthetax = MUONTriggerTrack.fthetax;
    fthetay = MUONTriggerTrack.fthetay;
    fGTPattern = MUONTriggerTrack.fGTPattern;

    return *this;
}
