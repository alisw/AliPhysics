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
/// This class provides access to EMCAL digits in raw data.
///
/// It loops over all EMCAL digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
/// usage: 
/// root > AliRawReaderFile rawReader(#event) ; 
/// root > AliPHOSRawStream input(&rawReader) ; 
/// root > while (input.Next()) ..... 
///
///////////////////////////////////////////////////////////////////////////////

#include "AliEMCALRawStream.h"
#include "AliRawReader.h"

ClassImp(AliEMCALRawStream)


//_____________________________________________________________________________
AliEMCALRawStream::AliEMCALRawStream(AliRawReader* rawReader) :
  AliAltroRawStream(rawReader),
  fId(-1),
  fPrevId(-1),
  fModule(-1),
  fPrevModule(-1)
{
// create an object to read EMCAL raw digits

  SelectRawData("EMCAL");

  fNoAltroMapping = kTRUE;
}

//_____________________________________________________________________________
AliEMCALRawStream::AliEMCALRawStream(const AliEMCALRawStream& stream) :
  AliAltroRawStream(stream),
  fId(-1),
  fPrevId(-1),
  fModule(-1),
  fPrevModule(-1)
{
  Fatal("AliEMCALRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliEMCALRawStream& AliEMCALRawStream::operator = (const AliEMCALRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliEMCALRawStream::~AliEMCALRawStream()
{
// destructor
}

//_____________________________________________________________________________
void AliEMCALRawStream::Reset()
{
  // reset emcal raw stream params
  AliAltroRawStream::Reset();
  fId = fPrevId = fModule = fPrevModule = -1;
}

//_____________________________________________________________________________
Bool_t AliEMCALRawStream::Next()
{
  // Read next EMCAL signal
  // Apply the EMCAL altro mapping to get
  // the module and id indeces
  fPrevModule = fModule;
  fPrevId = fId;
  if (AliAltroRawStream::Next()) {
    //    if (IsNewHWAddress())
    ApplyAltroMapping();
    return kTRUE;
  }
  else
    return kFALSE;
}

//_____________________________________________________________________________
void AliEMCALRawStream::ApplyAltroMapping()
{
  // Take the DDL index, load
  // the corresponding altro mapping
  // object and fill the module and id indeces
  fModule = fSegmentation[0];
  fId = fSegmentation[2];
}
