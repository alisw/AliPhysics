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

/////////////////////////////////////////////////////////////////////////
//                                                                     //
//       Prototype ESD class                                           //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "AliESDvertex.h"

ClassImp(AliESDvertex)

//_______________________________________________________________________
AliESDvertex::AliESDvertex():
  fNPrimary(0),
  fCoordinates(3),
  fErrorMatrix(6),
  fPrimaryTracks(0),
  fEffectiveMass(0),
  fEffectiveMassError(0)
{
  cout << "ESDvertex def ctor" << endl;
}

//_______________________________________________________________________
AliESDvertex::AliESDvertex(const AliESDvertex &esdv):
  TObject(esdv),
  fNPrimary(0),
  fCoordinates(0),
  fErrorMatrix(0),
  fPrimaryTracks(0),
  fEffectiveMass(0),
  fEffectiveMassError(0)
{
}


