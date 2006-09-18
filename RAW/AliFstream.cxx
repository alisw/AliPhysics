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

//-----------------------------------------------------------------
// This is the class which is to be used during the writing of
// simulated raw data (DDL files format).
// It is using the root functionality in order to deal correctly
// with little/big endian issue. By convention the detector raw
// data payload is stored always with little endian (this corresponds
// to the real life situation when the detector data is coming from
// the hardware). The implementation of this class is based on Root
// tobuf() method defined in Bytes.h
//-------------------------------------------------------------------------

#include <unistd.h>
#include <Riostream.h>

#include "AliFstream.h"
#include "AliLog.h"

ClassImp(AliFstream)

//______________________________________________________________________________
AliFstream::AliFstream():
  fFile(0x0),
  fBuffer(0x0),
  fBufferSize(0),
  fSwap(kFALSE)
{
  // Default constructor

  for (Int_t i = 0; i < 8; i++) fBuffer[i] = 0;
}

//______________________________________________________________________________
AliFstream::AliFstream(const char *fileName):
  fFile(0x0),
  fBuffer(0x0),
  fBufferSize(0),
  fSwap(kFALSE)
{
  // Constructor
  // Takes the input filename and
  // opens the output stream

#ifndef __DECCXX
  fFile = new fstream(fileName, ios::binary|ios::out);
#else
  fFile = new fstream(fileName, ios::out);
#endif

  // Check endianess
  UInt_t temp = 1;
  UChar_t *ptemp = (UChar_t *)&temp;
  if (!ptemp[0]) fSwap = kTRUE;
}

//______________________________________________________________________________
AliFstream::~AliFstream()
{
  // Destructor
  //
  if (fFile) {
    fFile->close();
    delete fFile;
  }
  if (fBuffer) delete [] fBuffer;
}

//______________________________________________________________________________
void AliFstream::Seekp(UInt_t position)
{
  // Go to a given position
  // inside the output stream
  if (fFile) fFile->seekp(position);
}

//______________________________________________________________________________
UInt_t AliFstream::Tellp()
{
  // Return the current
  // position inside the
  // output stream
  if (fFile) return fFile->tellp();
  else return 0;
}

//______________________________________________________________________________
void AliFstream::WriteBuffer(const char *buffer, UInt_t size, Bool_t force)
{
  // Write the buffer to a file
  // In case the user gives a 'force'
  // flag then the buffer is written
  // as it is. Otherwise, we check the
  // endianess and swap the buffer data
  // so that it is always stored using
  // little endian format.

  // The raw data payload size is always
  // 4 bytes aligned
  if ((size % 4) != 0)
    AliFatal(Form("Size of the buffer is not multiple of 4 (size = %d) !",size));

  if (force) {
    fFile->write(buffer,size);
  }
  else {
    if (!fSwap) {
      // Little endian - do nothing
      fFile->write(buffer,size);
    }
    else {
      // Big endian - swap the buffer contents
      if (size > fBufferSize) {
	if (fBuffer) delete [] fBuffer;
	fBuffer = new UChar_t[size];
	fBufferSize = size;
      }
      swab(buffer,fBuffer,size);
      fFile->write((const char *)fBuffer,size);
    }
  }
}
