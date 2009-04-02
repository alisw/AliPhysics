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
// Class AliMUONTriggerTrack
//---------------------------
// Reconstructed Trigger track in ALICE dimuon spectrometer
// Note: equivalent to AliMUONTriggerTrack for tracking,
// No need for a AliMUONTriggerTrackParam
// Author: Philippe Crochet
//-----------------------------------------------------------------------------

#include "AliMUONTriggerTrack.h"
#include "AliMUONTrackReconstructor.h" 
#include <Riostream.h>
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerTrack)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack()
  : TObject(),
    fx11(0),
    fy11(0),
    fthetax(0),
    fthetay(0),
    floTrgNum(0),
    fGTPattern(0),
    fHitsPatternInTrigCh(0)

{
  /// default ctr
      AliDebug(1,Form("this=%p",this));
}
//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t thetax, Float_t thetay, Int_t loTrgNum, Long_t theGTPattern, UShort_t hitsPatternInTrigCh)
    : TObject(),
      fx11(x11),
      fy11(y11),
      fthetax(thetax),
      fthetay(thetay),
      floTrgNum(loTrgNum),
      fGTPattern(theGTPattern),
      fHitsPatternInTrigCh(fHitsPatternInTrigCh)
{
/// ctor from local trigger output
        AliDebug(1,Form("this=%p x11=%f y11=%f thetax=%f thetay=%f loTrgNum=%d GTPattern=%ld HitsPatternInTrigCh %i",
                        this,x11,y11,thetax,thetay,loTrgNum,theGTPattern,fHitsPatternInTrigCh));

}

//__________________________________________________________________________
AliMUONTriggerTrack::~AliMUONTriggerTrack()
{
  /// Destructor
  AliDebug(1,Form("this=%p",this));
}

//__________________________________________________________________________
AliMUONTriggerTrack::AliMUONTriggerTrack (const AliMUONTriggerTrack& theMUONTriggerTrack)
    : TObject(theMUONTriggerTrack),
      fx11(theMUONTriggerTrack.fx11),
      fy11(theMUONTriggerTrack.fy11),
      fthetax(theMUONTriggerTrack.fthetax),
      fthetay(theMUONTriggerTrack.fthetay),
      floTrgNum(theMUONTriggerTrack.floTrgNum),
      fGTPattern(theMUONTriggerTrack.fGTPattern),
      fHitsPatternInTrigCh(theMUONTriggerTrack.fHitsPatternInTrigCh)  
{
///
/// copy ctor
///
        AliDebug(1,Form("this=%p copy ctor",this));

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
    floTrgNum = theMUONTriggerTrack.floTrgNum;
    fGTPattern = theMUONTriggerTrack.fGTPattern;
    fHitsPatternInTrigCh = theMUONTriggerTrack.fHitsPatternInTrigCh;

    return *this;
}

//__________________________________________________________________________
void
AliMUONTriggerTrack::Print(Option_t*) const
{
/// Printing

  cout << Form("(X,Y)11=(%7.2f,%7.2f) Theta(X,Y)=(%7.2f,%7.2f) LocalBoard #%3d GlobalTriggerPattern %x HitsPatternInTrigCh %x",
               fx11,fy11,fthetax,fthetay,floTrgNum,fGTPattern,fHitsPatternInTrigCh) << endl;
}
