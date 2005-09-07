// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include <iostream>

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCStandardIncludes.h"
#include "AliHLTTPCLogging.h"

#include "AliHLTTPCDDLRawReaderFile.h"

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

/** \class AliHLTTPCDDLRawReaderFile
<pre>
//_____________________________________________________________
// AliHLTTPCDDLRawReaderFile (taken from the offline AliROOT code,
// original authors: D.Favretto and A.K.Mohanty)
//
// This is the base class for reading ddl raw data 
// and providing information about digits
</pre>
*/


ClassImp(AliHLTTPCDDLRawReaderFile)


AliHLTTPCDDLRawReaderFile::AliHLTTPCDDLRawReaderFile(const Char_t* name, Bool_t addnum)
{
  // create an object to read digits from the given input file(s)
  // if addNumber is true, a number starting at 1 is appended to the file name

  fFileName = new Char_t[1024];
  strcpy(fFileName,name);
  if (!addnum) {
    fFileNumber = -1;
    fStream = new fstream(fFileName, ios::binary|ios::in);
  } else {
    fFileNumber = 0;
    fStream = NULL;
    OpenNextFile();
  }
  fMiniHeader = new AliHLTTPCDDLMiniHeader;
  fBuffer = NULL;
  fBufferSize = 0;
}

AliHLTTPCDDLRawReaderFile::~AliHLTTPCDDLRawReaderFile()
{
  // close the input file
  if(fFileName) delete fFileName;

  if (fStream) {
    if (fStream->is_open()) fStream->close();
    delete fStream;
  }
  delete fMiniHeader;
  if (fBuffer) delete[] fBuffer;
}

Bool_t AliHLTTPCDDLRawReaderFile::OpenNextFile()
{
  if (fStream) {
    if (fStream->is_open()) fStream->close();
    delete fStream;
    fStream = NULL;
  }
  if (fFileNumber < 0) {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLRawReaderFile::OpenNextFile","File")
      <<"Could not open file, file number is negative."<<ENDLOG;
    return kFALSE;
  }

  fFileNumber++;
  Char_t fileName[1024];
  sprintf(fileName, "%s%d", fFileName, fFileNumber);

  fStream = new fstream(fileName, ios::binary|ios::in);
  return (fStream->is_open());
}

Bool_t AliHLTTPCDDLRawReaderFile::ReadMiniHeader()
{
  // read a mini header at the current stream position
  // returns kFALSE if the mini header could not be read

  if (!fStream) return kFALSE;
  do {
    if (fCount > 0) fStream->seekg(Int_t(fStream->tellg()) + fCount);
    while (!fStream->read((Char_t*) fMiniHeader, sizeof(AliHLTTPCDDLMiniHeader))) {
      if (!OpenNextFile()) return kFALSE;
    }
    //cout << fMiniHeader->fSize << " " << fMiniHeader->fDetectorID << " " << fMiniHeader->fVersion << " " << fMiniHeader->fCompressionFlag << " " << fMiniHeader->fDDLID << endl;
    //cout << "loop here " << (Int_t)fMiniHeader->fDDLID<< endl;
    CheckMiniHeader();
    fCount = fMiniHeader->fSize;
  } while (!IsSelected());
  return kTRUE;
}

Bool_t AliHLTTPCDDLRawReaderFile::ReadNextData(UChar_t*& data)
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
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLRawReaderFile::ReadNextData","Data")
      <<"Could not read next data!"<<ENDLOG;
    return kFALSE;
  }
  fCount = 0;

  data = fBuffer;
  return kTRUE;
}

Bool_t AliHLTTPCDDLRawReaderFile::ReadNext(UChar_t* data, Int_t size)
{
  // reads the next block of data at the current stream position
  // returns kFALSE if the data could not be read

  if (!fStream->read((Char_t*) data, size)) {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLRawReaderFile::ReadNext","Data")
      <<"Could not read next data!"<<ENDLOG;
    return kFALSE;
  }
  fCount -= size;
  return kTRUE;
}

Bool_t AliHLTTPCDDLRawReaderFile::Reset()
{
  // reset the current stream position to the beginning of the file

  if ((fFileNumber > 0) && fStream) {
    if (fStream->is_open()) fStream->close();
    delete fStream;
    fStream = NULL;
    fFileNumber = 0;
  }

  if (!fStream) {
    if (fFileNumber < 0) {
      fStream = new fstream(fFileName, ios::binary|ios::in);
    } else {
      if (!OpenNextFile()){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDDLRawReaderFile::Reset","Data")
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
