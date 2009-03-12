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

/* $Id: AliTRDtrackletMCM.cxx 28397 2008-09-02 09:33:00Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  MCM tracklet                                                          //
//                                                                        //
//  Author: J. Klein (Jochen.Klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrackletMCM.h"
#include "AliLog.h"

ClassImp(AliTRDtrackletMCM)

AliTRDtrackletMCM::AliTRDtrackletMCM(UInt_t trackletWord) :
  AliTRDtrackletBase(),
  fGeo(0x0),
  fHCId(-1),
  fTrackletWord(trackletWord), 
  fMCM(-1), 
  fROB(-1), 
  fLabel(-1)
{ 
    fGeo = new AliTRDgeometry();
}

AliTRDtrackletMCM::AliTRDtrackletMCM(UInt_t trackletWord, Int_t hcid) :
  AliTRDtrackletBase(),
  fGeo(0x0),
  fHCId(hcid),
  fTrackletWord(trackletWord), 
  fMCM(-1),
  fROB(-1),
  fLabel(-1)
{ 
    fGeo = new AliTRDgeometry();
}

AliTRDtrackletMCM::AliTRDtrackletMCM(UInt_t trackletWord, Int_t hcid, Int_t rob, Int_t mcm) :
  AliTRDtrackletBase(),
  fGeo(0x0),
  fHCId(hcid),
  fTrackletWord(trackletWord), 
  fMCM(mcm),
  fROB(rob),
  fLabel(-1)
{ 
    fGeo = new AliTRDgeometry();
}

AliTRDtrackletMCM::AliTRDtrackletMCM(const AliTRDtrackletMCM &rhs) :
  AliTRDtrackletBase(rhs),
  fGeo(0x0),
  fHCId(rhs.fHCId),
  fTrackletWord(rhs.fTrackletWord),
  fMCM(rhs.fMCM),
  fROB(rhs.fROB),
  fLabel(rhs.fLabel)
{
    fGeo = new AliTRDgeometry();
}

AliTRDtrackletMCM::~AliTRDtrackletMCM() 
{
    delete fGeo;
}

Int_t AliTRDtrackletMCM::GetYbin() const {
  // returns (signed) value of Y
  if (fTrackletWord & 0x1000) {
    return -((~(fTrackletWord-1)) & 0x1fff);
  }
  else {
    return (fTrackletWord & 0x1fff);
  }
}

Int_t AliTRDtrackletMCM::GetdY() const 
{
  // returns (signed) value of the deflection length
  if (fTrackletWord & (1 << 19)) {
    return -((~((fTrackletWord >> 13) - 1)) & 0x7f);
  }
  else {
    return ((fTrackletWord >> 13) & 0x7f);
  }
}
