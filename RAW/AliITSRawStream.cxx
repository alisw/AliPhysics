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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a base class for providing access to ITS digits in raw data.
///
/// Derived class should implement the Next method.
///
/// It loops over all ITS digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliITSRawStream.h"

ClassImp(AliITSRawStream)


AliITSRawStream::AliITSRawStream(AliRawReader* rawReader)
{
// create an object to read ITS raw digits

  fRawReader = rawReader;
  fModuleID = fPrevModuleID = fCoord1 = fCoord2 = fSignal = -1;
}

AliITSRawStream::AliITSRawStream(const AliITSRawStream& stream) :
  TObject(stream)
{
  Fatal("AliITSRawStream", "copy constructor not implemented");
}

AliITSRawStream& AliITSRawStream::operator = (const AliITSRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

