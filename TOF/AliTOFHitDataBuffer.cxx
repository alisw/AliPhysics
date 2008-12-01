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

/*
$Log$
Revision 1.1  2007/04/27 11:03:09  arcelli
container for TOF raw data

 authors: Roberto Preghenella, preghenella@bo.infn.it
          with contribution from Chiara Zampolli, zampolli@bo.infn.it 
*/


////////////////////////////////////////////////////////////////////////
//                                                                    //
//     This class provides access to TOF raw data in DDL files.       //
//                                                                    //
//      It loops over all TOF raw data given by the AliRawReader.     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

//#include "AliTOFHitData.h"
#include "AliTOFHitDataBuffer.h"

ClassImp(AliTOFHitDataBuffer)

AliTOFHitDataBuffer::AliTOFHitDataBuffer(Int_t bufferSize) :
  TObject(),
  fBufferSize(bufferSize),
  fBuffer(new AliTOFHitData[bufferSize]),
  fEntries(0)
{
}

//-----------------------------------------------------------------------------
AliTOFHitDataBuffer::AliTOFHitDataBuffer(const AliTOFHitDataBuffer &source):
  TObject(source),
  fBufferSize(source.fBufferSize),
  fBuffer(new AliTOFHitData[fBufferSize]),
  fEntries(source.fEntries)
{
  // copy ctr
  for (Int_t i = 0; i < fEntries; ++i)
    fBuffer[i] = source.fBuffer[i];
}

//-----------------------------------------------------------------------------
AliTOFHitDataBuffer& AliTOFHitDataBuffer::operator=(const AliTOFHitDataBuffer & source) 
{ 
  // ass operator
  if(this!=&source) {
    TObject::operator=(source);
    fEntries = source.fEntries < fBufferSize ? source.fEntries : fBufferSize;
    for (Int_t i = 0; i < fEntries; ++i) fBuffer[i]=source.fBuffer[i];
  }
  return *this;
}

//-----------------------------------------------------------------------------
AliTOFHitDataBuffer::~AliTOFHitDataBuffer()
{
  delete [] fBuffer;
}

//-----------------------------------------------------------------------------
Bool_t AliTOFHitDataBuffer::Add(AliTOFHitData &HitData) {
  // adding a new entry 
  if (fEntries >= fBufferSize){
    AliError("The buffer is completely full. ");
    return kTRUE;
  }
  fBuffer[fEntries++] = HitData;
  return kFALSE;
}


