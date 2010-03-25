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

// *
// *
// *
// * this class defines the DeltaBCOffset object to be stored
// * in OCDB in order to apply DeltaBC correction during 
// * reconstruction. 
// *
// *
// *

#include "AliTOFDeltaBCOffset.h"

ClassImp(AliTOFDeltaBCOffset)

//_________________________________________________________

AliTOFDeltaBCOffset::AliTOFDeltaBCOffset() :
  TObject(),
  fDeltaBCOffset(0)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFDeltaBCOffset::~AliTOFDeltaBCOffset()
{
  /*
   * default destructor
   */
}

//_________________________________________________________

AliTOFDeltaBCOffset::AliTOFDeltaBCOffset(const AliTOFDeltaBCOffset &source) :
  TObject(source),
  fDeltaBCOffset(source.fDeltaBCOffset)
{
  /*
   * copy constructor
   */
}

//_________________________________________________________

AliTOFDeltaBCOffset &
AliTOFDeltaBCOffset::operator=(const AliTOFDeltaBCOffset &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  fDeltaBCOffset = source.fDeltaBCOffset;
  return *this;
}

