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
//
// This is a class for reading a raw data file and providing
// information about digits
//
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReaderFile.h"


ClassImp(AliRawReaderFile)


AliRawReaderFile::AliRawReaderFile(const char* fileName, Bool_t addNumber)
{
// create an object to read digits from the given input file(s)
// if addNumber is true, a number starting at 1 is appended to the file name

  fFileName = fileName;
  if (!addNumber) {
    fFileNumber = -1;
#ifndef __DECCXX
    fStream = new fstream(fileName, ios::binary|ios::in);
#else
    fStream = new fstream(fileName, ios::in);
#endif
  } else {
    fFileNumber = 0;
    fStream = NULL;
    OpenNextFile();
  }
  fMiniHeader = new AliMiniHeader;
  fBuffer = NULL;
  fBufferSize = 0;
}

AliRawReaderFile::~AliRawReaderFile()
{
// close the input file

  if (fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
  }
  delete fMiniHeader;
  if (fBuffer) delete[] fBuffer;
}


Bool_t AliRawReaderFile::OpenNextFile()
{
  if (fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
    fStream = NULL;
  }
  if (fFileNumber < 0) return kFALSE;

  fFileNumber++;
  char fileName[256];
  sprintf(fileName, "%s%d", fFileName.Data(), fFileNumber);
#ifndef __DECCXX 
  fStream = new fstream(fileName, ios::binary|ios::in);
#else
  fStream = new fstream(fileName, ios::in);
#endif
#if defined(__HP_aCC) || defined(__DECCXX)
  return (fStream->rdbuf()->is_open());
#else
  return (fStream->is_open());
#endif
}


Bool_t AliRawReaderFile::ReadMiniHeader()
{
// read a mini header at the current stream position
// returns kFALSE if the mini header could not be read

  if (!fStream) return kFALSE;
  do {
    if (fCount > 0) fStream->seekg(Int_t(fStream->tellg()) + fCount);
    while (!fStream->read((char*) fMiniHeader, sizeof(AliMiniHeader))) {
      if (!OpenNextFile()) return kFALSE;
    }
    CheckMiniHeader();
    fCount = fMiniHeader->fSize;
  } while (!IsSelected());
  return kTRUE;
}

Bool_t AliRawReaderFile::ReadNextData(UChar_t*& data)
{
// reads the next payload at the current stream position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadMiniHeader()) return kFALSE;
  }
  if (fBufferSize < fCount) {
    if (fBuffer) delete[] fBuffer;
    fBufferSize = Int_t(fCount*1.2);
    fBuffer = new UChar_t[fBufferSize];
  }
  if (!fStream->read((char*) fBuffer, fCount)) {
    Error("ReadNext", "could not read data!");
    return kFALSE;
  }
  fCount = 0;

  data = fBuffer;
  return kTRUE;
}

Bool_t AliRawReaderFile::ReadNext(UChar_t* data, Int_t size)
{
// reads the next block of data at the current stream position
// returns kFALSE if the data could not be read

  if (!fStream->read((char*) data, size)) {
    Error("ReadNext", "could not read data!");
    return kFALSE;
  }
  fCount -= size;
  return kTRUE;
}


Bool_t AliRawReaderFile::Reset()
{
// reset the current stream position to the beginning of the file

  if ((fFileNumber > 0) && fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
    fStream = NULL;
    fFileNumber = 0;
  }

  if (!fStream) {
    if (fFileNumber < 0) {
#ifndef __DECCXX
      fStream = new fstream(fFileName, ios::binary|ios::in);
#else
      fStream = new fstream(fFileName, ios::in);
#endif
    } else {
      if (!OpenNextFile()) return kFALSE;
    }
  }

  if (!fStream || !fStream->rdbuf()->is_open()) return kFALSE;
  fStream->seekg(0);
  fCount = 0;
  return kTRUE;
}

