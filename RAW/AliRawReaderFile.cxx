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
/// This is a class for reading raw data files.
///
/// The files of one event are expected to be in one directory. The name 
/// of the directory is "raw" + the event number. Each file contains
/// the raw data (with data header) of one DDL. The convention for the
/// file names is "DET_#DDL.ddl". "DET" is the name of the detector and
/// "#DDL" is the unique equipment ID.
///
/// The constructor of AliRawReaderFile takes the event number or the
/// directory name as argument.
/// 
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReaderFile.h"
#include <TSystem.h>


ClassImp(AliRawReaderFile)


AliRawReaderFile::AliRawReaderFile(Int_t eventNumber) :
  fDirName("raw"),
  fDirectory(NULL),
  fStream(NULL),
  fEquipmentId(-1),
  fBuffer(NULL),
  fBufferSize(0)
{
// create an object to read digits from the given event

  fDirName += eventNumber;
  fDirectory = gSystem->OpenDirectory(fDirName);
  if (!fDirectory) {
    Error("AliRawReaderFile", "could not open directory %s", fDirName.Data());
  }
  OpenNextFile();
  fHeader = new AliRawDataHeader;
}

AliRawReaderFile::AliRawReaderFile(const char* dirName) :
  fDirName(dirName),
  fDirectory(NULL),
  fStream(NULL),
  fEquipmentId(-1),
  fBuffer(NULL),
  fBufferSize(0)
{
// create an object to read digits from the given directory

  fDirectory = gSystem->OpenDirectory(fDirName);
  if (!fDirectory) {
    Error("AliRawReaderFile", "could not open directory %s", fDirName.Data());
  }
  OpenNextFile();
  fHeader = new AliRawDataHeader;
}

AliRawReaderFile::AliRawReaderFile(const AliRawReaderFile& rawReader) :
  AliRawReader(rawReader)
{
  Fatal("AliRawReaderFile", "copy constructor not implemented");
}

AliRawReaderFile& AliRawReaderFile::operator = (const AliRawReaderFile& 
					      /* rawReader */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

AliRawReaderFile::~AliRawReaderFile()
{
// close the input file

  if (fDirectory) gSystem->FreeDirectory(fDirectory);
  if (fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
  }
  delete fHeader;
  if (fBuffer) delete[] fBuffer;
}


Bool_t AliRawReaderFile::OpenNextFile()
{
// open the next file
// returns kFALSE if the current file is the last one

  if (fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
    fStream = NULL;
    fEquipmentId = -1;
  }

  if (!fDirectory) return kFALSE;
  TString entry;
  while (entry = gSystem->GetDirEntry(fDirectory)) {
    if (entry.IsNull()) return kFALSE;
    if (!entry.EndsWith(".ddl")) continue;
    char* fileName = gSystem->ConcatFileName(fDirName, entry);
#ifndef __DECCXX 
    fStream = new fstream(fileName, ios::binary|ios::in);
#else
    fStream = new fstream(fileName, ios::in);
#endif
    break;
  }

  if (!fStream) return kFALSE;
  entry.Remove(0, entry.Last('_')+1);
  entry.Remove(entry.Length()-4);
  fEquipmentId = atoi(entry.Data());
#if defined(__HP_aCC) || defined(__DECCXX)
  return (fStream->rdbuf()->is_open());
#else
  return (fStream->is_open());
#endif
}


Bool_t AliRawReaderFile::ReadHeader()
{
// read a data header at the current stream position
// returns kFALSE if the mini header could not be read

  if (!fStream) return kFALSE;
  do {
    if (fCount > 0) fStream->seekg(Int_t(fStream->tellg()) + fCount);
    while (!fStream->read((char*) fHeader, sizeof(AliRawDataHeader))) {
      if (!OpenNextFile()) return kFALSE;
    }
    if (fHeader->fSize != 0xFFFFFFFF) {
      fCount = fHeader->fSize - sizeof(AliRawDataHeader);
    } else {
      UInt_t currentPos = fStream->tellg();
      fStream->seekg(0, ios::end);
      fCount = UInt_t(fStream->tellg()) - currentPos;
      fStream->seekg(currentPos);
    }
  } while (!IsSelected());
  return kTRUE;
}

Bool_t AliRawReaderFile::ReadNextData(UChar_t*& data)
{
// reads the next payload at the current stream position
// returns kFALSE if the data could not be read

  while (fCount == 0) {
    if (!ReadHeader()) return kFALSE;
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

  void* directory = gSystem->OpenDirectory(fDirName);
  if (!directory) {
    Error("Reset", "could not open directory %s", fDirName.Data());
    return kFALSE;
  }

  if (fStream) {
#if defined(__HP_aCC) || defined(__DECCXX)
    if (fStream->rdbuf()->is_open()) fStream->close();
#else
    if (fStream->is_open()) fStream->close();
#endif
    delete fStream;
    fStream = NULL;
  }

  if (fDirectory) gSystem->FreeDirectory(fDirectory);
  fDirectory = directory;

  OpenNextFile();
  fCount = 0;
  return kTRUE;
}

