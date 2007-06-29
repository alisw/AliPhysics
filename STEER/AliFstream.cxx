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
#include <stdio.h>
#include "AliFstream.h"
#include "AliLog.h"


ClassImp(AliFstream)

//______________________________________________________________________________
AliFstream::AliFstream():
  fFile(0x0),
  fBuffer(0x0),
  fBufferSize(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliFstream::AliFstream(const char *fileName):
  fFile(0x0),
  fBuffer(0x0),
  fBufferSize(0)
{
  // Constructor
  // Takes the input filename and
  // opens the output stream

#ifndef __DECCXX
  fFile = new fstream(fileName, ios::binary|ios::out);
#else
  fFile = new fstream(fileName, ios::out);
#endif
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
UInt_t AliFstream::Swap(UInt_t x)
{
   // Swap the endianess of the integer value 'x'

   return (((x & 0x000000ffU) << 24) | ((x & 0x0000ff00U) <<  8) |
           ((x & 0x00ff0000U) >>  8) | ((x & 0xff000000U) >> 24));
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
  
  if ((size % sizeof(UInt_t)) != 0)
    AliFatal(Form("Size of the buffer is not multiple of 4 (size = %d) !",size));
  
  if (force) {
    fFile->write(buffer,size);
  }
  else {
#ifdef R__BYTESWAP
    fFile->write(buffer,size);
#else
    size /= sizeof(UInt_t);

    if (size > fBufferSize) {
      if (fBuffer) delete [] fBuffer;
      fBuffer = new UInt_t[size];
      fBufferSize = size;
    }

    UInt_t *buf = (UInt_t *)buffer;
    for (UInt_t i = 0; i < size; i++, buf++) {
      UInt_t value = Swap(*buf);
      memcpy(fBuffer+i, &value, sizeof(UInt_t));
    }

    fFile->write((const char *)fBuffer,size*sizeof(UInt_t));
#endif
  }
}
