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

//____________________________________________________________________
//                                                                          
// Buffer to read RAW ALTRO FMD format from a AliRawReader 
// 
// This class derives from AliAltroBuffer, but overloads the memer
// function Next to do some extra processing.  In particular, it tries
// to autodetect the sample rate.  If zero-suppression was used when
// writing the raw data, then the automatic discovery will not work,
// and the sample rate should be set explicitly. 
//
#include "AliFMDRawStream.h"		// ALIFMDRAWSTREAM_H
#include <AliRawReader.h>		// ALIRAWREADER_H

//____________________________________________________________________
ClassImp(AliFMDRawStream)

//____________________________________________________________________
AliFMDRawStream::AliFMDRawStream(AliRawReader* reader, UShort_t sampleRate) 
  : AliAltroRawStream(reader), 
    fSampleRate(sampleRate),
    fPrevTime(-1), 
    fExplicitSampleRate(kFALSE)
{
  if (fSampleRate > 0) fExplicitSampleRate = kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left
  fPrevTime = fTime;
  if (AliAltroRawStream::Next()) {
    if (!fExplicitSampleRate && fPrevPad != fPad) 
      fSampleRate = fTimeBunch / 128;
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
//
// EOF
//
