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

/* $Id$*/

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

#include "AliMUONTriggerTrack.h"
#include "AliMUONTrackReconstructor.h" 


//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack()
  : TObject(),
    fx11(0),
    fy11(0),
    fthetax(0),
    fthetay(0),
    fGTPattern(0)

{
  /// default ctr
}
//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t thetax, Float_t thetay, Long_t theGTPattern)
    : TObject(),
      fx11(x11),
      fy11(y11),
      fthetax(thetax),
      fthetay(thetay),
      fGTPattern(theGTPattern)
{
/// ctor from local trigger output
}

//__________________________________________________________________________
AliMUONTriggerTrack::~AliMUONTriggerTrack()
{
  /// Destructor
    ;
    
}

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack (const AliMUONTriggerTrack& theMUONTriggerTrack)
    : TObject(theMUONTriggerTrack),
      fx11(theMUONTriggerTrack.fx11),
      fy11(theMUONTriggerTrack.fy11),
      fthetax(theMUONTriggerTrack.fthetax),
      fthetay(theMUONTriggerTrack.fthetay),
      fGTPattern(theMUONTriggerTrack.fGTPattern)    
{
///
/// copy ctor
///
}
      
//__________________________________________________________________________
AliMUONTriggerTrack & AliMUONTriggerTrack::operator=(const AliMUONTriggerTrack&
theMUONTriggerTrack)
{
/// Assignment operator

    // check assignement to self
    if (this == &theMUONTriggerTrack)
	return *this;
    
    /// base class assignement
    TObject::operator=(theMUONTriggerTrack);

    fx11 = theMUONTriggerTrack.fx11;
    fy11 = theMUONTriggerTrack.fy11;
    fthetax = theMUONTriggerTrack.fthetax;
    fthetay = theMUONTriggerTrack.fthetay;
    fGTPattern = theMUONTriggerTrack.fGTPattern;

    return *this;
}
