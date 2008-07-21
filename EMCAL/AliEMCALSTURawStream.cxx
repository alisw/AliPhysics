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

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to EMCALSTU DDL raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliEMCALSTURawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliEMCALSTURawStream)

//_____________________________________________________________________________
AliEMCALSTURawStream::AliEMCALSTURawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fNum2x2Words(0)
{
  // reset arrays
  memset(fJetPatchWords, 0, sizeof(fJetPatchWords));

  // Reset RawStream
  //
  fRawReader->Reset();

  // select the raw data corresponding to
  // the EMCALSTU detector id
  fRawReader->Select("EMCAL", kEMCALSTUDDL);
}

//_____________________________________________________________________________
AliEMCALSTURawStream::AliEMCALSTURawStream(const AliEMCALSTURawStream& stream) :
  TObject(stream),
  fRawReader(NULL),
  fNum2x2Words(0)
{
  // Copy constructor
  AliFatal("Copy constructor not implemented");
}

//_____________________________________________________________________________
AliEMCALSTURawStream& AliEMCALSTURawStream::operator = (const AliEMCALSTURawStream& source
							/* stream */)
{
  // assignment operator; use copy ctor
  if (&source == this) return *this;

  new (this) AliEMCALSTURawStream(source);
  return *this;
}

//_____________________________________________________________________________
AliEMCALSTURawStream::~AliEMCALSTURawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliEMCALSTURawStream::Reset()
{
  // reset raw stream params
  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliEMCALSTURawStream::Next()
{
  // read the whole EMCALSTU raw data stream
  // return kFALSE in case of error

  UChar_t *data = NULL;

  // EMCALSTU raw data should contain CDH so we don't need any mods a la
  //  fRawReader->RequireHeader(kTRUE);

  // let's pick up the data block, if we can
  if (!fRawReader->ReadNextData(data)) {
    return kFALSE;
  }

  // check that the amount of data we picked up is what we expected
  // 20080721 (DS): the rest of this method needs to be verified as the exact data format 
  // is defined
  int nwordsExtra = fRawReader->GetDataSize() - (kNumJetPatchWords + kNumGammaJetPatchWords);
  if ( nwordsExtra < 0 ) {
    AliError(Form("Wrong EMCALSTU raw data size: %d", fRawReader->GetDataSize()));
    return kFALSE;
  }

  // unpack the data into our local arrays
  // jet-patch
  int iword = 0;
  int ioffset = 0;
  for (iword = 0; iword<kNumJetPatchWords; iword++) {
    fJetPatchWords[iword] = (data[iword] & 0xffffffff);
  }
  // gamma-jet patch
  ioffset += kNumJetPatchWords;
  for (iword = 0; iword<kNumGammaJetPatchWords; iword++) {
    fGammaJetPatchWords[iword] = (data[iword+ioffset] & 0xffffffff);
  }

  // also 2x2sum's, if there appears to be enough/more data..
  fNum2x2Words = 0;
  if ( nwordsExtra > kNum2x2Words ) {
    ioffset += kNumGammaJetPatchWords;
    for (iword = 0; iword<kNum2x2Words; iword++) {
      f2x2Words[iword] = (data[iword+ioffset] & 0xffffffff);
    }
    fNum2x2Words = kNum2x2Words;
  }

  return kTRUE;
}

