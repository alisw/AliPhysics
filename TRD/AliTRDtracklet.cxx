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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  A TRD tracklet                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDtracklet.h"

ClassImp(AliTRDtracklet)

//_____________________________________________________________________________
AliTRDtracklet::AliTRDtracklet()
  :TObject()
  ,fY(0)
  ,fZ(0)
  ,fX(0)
  ,fAlpha(0)
  ,fSigma2(0)
  ,fP0(0)
  ,fP1(0)
  ,fNFound(0)
  ,fNCross(0)
  ,fPlane(0)
  ,fExpectedSigma2(0)
  ,fChi2(0)
  ,fTilt(0)
  ,fMaxPos(0)
  ,fMaxPos4(0)
  ,fMaxPos5(0) 
{ 
  //
  // Default contructor
  //

}

//_____________________________________________________________________________
AliTRDtracklet::~AliTRDtracklet()
{
  //
  // Destructor
  //

}
