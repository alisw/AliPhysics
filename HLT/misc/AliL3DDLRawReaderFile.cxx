// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3RootTypes.h"
#include "AliL3StandardIncludes.h"
#include "AliL3Logging.h"

#include "AliL3DDLRawReaderFile.h"

#if __GNUC__ == 3
using namespace std;
#endif

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

/** \class AliL3DDLRawReaderFile
<pre>
//_____________________________________________________________
// AliL3DDLRawReaderFile (taken from the offline AliROOT code,
// original authors: D.Favretto and A.K.Mohanty)
//
// This is the base class for reading ddl raw data 
// and providing information about digits
</pre>
*/

ClassImp(AliL3DDLRawReaderFile)

AliL3DDLRawReaderFile::AliL3DDLRawReaderFile(const Char_t* name, Bool_t addnum)
{
  // create an object to read digits from the given input file(s)
  // if addNumber is true, a number starting at 1 is appended to the file name

  fFileName = new Char_t[1024];
  strcpy(fFileName,name);
  if (!addnum) {
    fFileNumber = -1;
#ifndef __DECCXX
    fStream = new fstream(fFileName, ios::binary|ios::in);
#else
    fStream = new fstream(fFileName, ios::in);
#endif
  } else {
    fFileNumber = 0;
    fStream = NULL;
    OpenNextFile();
  }
  fMiniHeader = new AliL3DDLMiniHeader;
  fBuffer = NULL;
  fBufferSize = 0;
}

AliL3DDLRawReaderFile::~AliL3DDLRawReaderFile()
{
  // close the input file
  if(fFileName) delete fFileName;

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

Bool_t AliL3DDLRawReaderFile::OpenNextFile()
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
  if (fFileNumber < 0) {
    LOG(AliL3Log::kError,"AliL3DDLRawReaderFile::OpenNextFile","File")
      <<"Could not open file, file number is negative."<<ENDLOG;
    return kFALSE;
  }

  fFileNumber++;
  Char_t fileName[1024];
  sprintf(fileName, "%s%d", fFileName, fFileNumber);

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

Bool_t AliL3DDLRawReaderFile::ReadMiniHeader()
{
  // read a mini header at the current stream position
  // returns kFALSE if the mini header could not be read

  if (!fStream) return kFALSE;
  do {
    if (fCount > 0) fStream->seekg(Int_t(fStream->tellg()) + fCount);
    while (!fStream->read((Char_t*) fMiniHeader, sizeof(AliL3DDLMiniHeader))) {
      if (!OpenNextFile()) return kFALSE;
    }
    //cout << fMiniHeader->fSize << " " << fMiniHeader->fDetectorID << " " << fMiniHeader->fVersion << " " << fMiniHeader->fCompressionFlag << " " << fMiniHeader->fDDLID << endl;
    //cout << "loop here " << (Int_t)fMiniHeader->fDDLID<< endl;
    CheckMiniHeader();
    fCount = fMiniHeader->fSize;
  } while (!IsSelected());
  return kTRUE;
}

Bool_t AliL3DDLRawReaderFile::ReadNextData(UChar_t*& data)
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
  if (!fStream->read((Char_t*) fBuffer, fCount)) {
    LOG(AliL3Log::kError,"AliL3DDLRawReaderFile::ReadNextData","Data")
      <<"Could not read next data!"<<ENDLOG;
    return kFALSE;
  }
  fCount = 0;

  data = fBuffer;
  return kTRUE;
}

Bool_t AliL3DDLRawReaderFile::ReadNext(UChar_t* data, Int_t size)
{
  // reads the next block of data at the current stream position
  // returns kFALSE if the data could not be read

  if (!fStream->read((Char_t*) data, size)) {
    LOG(AliL3Log::kError,"AliL3DDLRawReaderFile::ReadNext","Data")
      <<"Could not read next data!"<<ENDLOG;
    return kFALSE;
  }
  fCount -= size;
  return kTRUE;
}

Bool_t AliL3DDLRawReaderFile::Reset()
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
      if (!OpenNextFile()){
	LOG(AliL3Log::kError,"AliL3DDLRawReaderFile::Reset","Data")
	  <<"Could not reset data stream!"<<ENDLOG;
	return kFALSE;
      }
    }
  }

  if (!fStream || !fStream->rdbuf()->is_open()) return kFALSE;
  fStream->seekg(0);
  fCount = 0;
  return kTRUE;
}
