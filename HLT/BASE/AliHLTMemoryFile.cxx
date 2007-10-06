// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTMemoryFile.cxx
    @author Matthias Richter
    @date   
    @brief  ROOT file in memory.                                          */

#include "AliHLTMemoryFile.h"
#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMemoryFile);

AliHLTMemoryFile::AliHLTMemoryFile()
  :
  fpBuffer(NULL),
  fBufferSize(0),
  fPosition(0),
  fSize(0),
  fErrno(0),
  fbClosed(0),
  fHeaderSize(0),
  fTrailerSize(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTMemoryFile::AliHLTMemoryFile(void* pBuffer, int iSize)
  :
  TFile("/dev/null", "CREATE"),
  AliHLTLogging(),
  fpBuffer((char*)pBuffer),
  fBufferSize(iSize),
  fPosition(0),
  fSize(0),
  fErrno(0),
  fbClosed(0),
  fHeaderSize(0),
  fTrailerSize(0)
{
  // see header file for class documentation
  //HLTDebug("created memory file %p, capacity %d, ROOT version %d", this, fBufferSize, fVersion);
}

AliHLTMemoryFile::~AliHLTMemoryFile()
{
  // see header file for function documentation
  //HLTDebug("deleting file %p size %d", this, fSize);
  if (!fbClosed) {
    HLTWarning("memory file not closed, possible data loss");
  }
}

void AliHLTMemoryFile::Close(const Option_t*)
{
  CloseMemoryFile();
}

int AliHLTMemoryFile::CloseMemoryFile(int bFlush)
{
  fErrno=0;
  if (fbClosed) return 0;
  if (bFlush) {
    TFile::Close();
  }
  fpBuffer=NULL;
  fBufferSize=fPosition=0;
  fbClosed=1;
  if (fErrno==ENOSPC) {
    HLTError("error flushing memory file, buffer too small");
  } else if (fErrno>0) {
    HLTError("error flushing memory file");
  }
  return -fErrno;
}

Int_t    AliHLTMemoryFile::SysOpen(const char* /*pathname*/, Int_t /*flags*/, UInt_t /*mode*/)
{
  // see header file for function documentation
  if (fpBuffer==NULL || fSize==0) return 1;
  //HLTDebug("opening file %p capacity %d", this, fSize);
  fErrno=0;
  errno=fErrno=ENOSPC;
  return -1;
}

Int_t    AliHLTMemoryFile::SysClose(Int_t /*fd*/)
{
  // see header file for function documentation
  //HLTDebug("closing file %p size %d", this, fSize);
  return 0;
}

Int_t    AliHLTMemoryFile::SysRead(Int_t /*fd*/, void *buf, Int_t len)
{
  // see header file for function documentation
  if (buf==NULL) return 0;
  fErrno=0;
  //HLTDebug("reading buffer of size %d at position %d", len, fPosition);
  if (fpBuffer==NULL || fBufferSize==0) return 0;
  int read=len<fSize-fPosition?len:fSize-fPosition;
  memcpy(buf, fpBuffer+fPosition, read);
  fPosition+=read;
  if (fPosition>=fSize) fSize=fPosition+1;
  return read;
}

Int_t    AliHLTMemoryFile::SysWrite(Int_t /*fd*/, const void *buf, Int_t len)
{
  // see header file for function documentation
  if (buf==NULL) return 0;
  fErrno=0;
  //HLTDebug("writing buffer of size %d at position %d", len, fPosition);
  if (len<fBufferSize-fPosition) {
    memcpy(fpBuffer+fPosition, buf, len);
    fPosition+=len;
    if (fPosition>=fSize) fSize=fPosition+1;
    return len;
  }
  errno=fErrno=ENOSPC;
  return -1;
}

Long64_t AliHLTMemoryFile::SysSeek(Int_t /*fd*/, Long64_t offset, Int_t whence)
{
  // see header file for function documentation
  //HLTDebug("seek %d from %d", offset, whence);
  fErrno=0;
  int position=(int)offset;
  switch (whence) {
  case SEEK_SET:
    // nothing to do
    break;
  case SEEK_CUR:
    position+=fPosition;
    break;
  case SEEK_END:
    position+=fSize;
  default:
    position=-1;
    errno=EINVAL;
  }
  if (position>=0) {
    if (position<fBufferSize) {
      fPosition=position;
    } else {
      position=-1;
      errno=fErrno=ENOSPC;
    }
  }
  return position;
}

Int_t    AliHLTMemoryFile::SysStat(Int_t /*fd*/, Long_t */*id*/, Long64_t *size, Long_t */*flags*/, Long_t */*modtime*/)
{
  // see header file for function documentation
  if (size) *size=fSize;
  return 0;
}

Int_t    AliHLTMemoryFile::SysSync(Int_t /*fd*/)
{
  // see header file for function documentation
  return 0;
}

int AliHLTMemoryFile::WriteHeaderBuffer(const char* pHeader, int size)
{
  // see header file for function documentation
  fErrno=0;
  if (fHeaderSize==0) {
    if (fSize+size<fBufferSize) {
      if (fSize>0) {
	// move exiting data
	memcpy(fpBuffer+size, fpBuffer, fSize);
      }
      memcpy(fpBuffer, pHeader, size);
      fpBuffer+=size;
      fPosition+=size;
      fBufferSize-=size;
      fHeaderSize=size;
    } else {
      HLTError("no space left in memory file");
      fErrno=ENOSPC;
    }
  } else {
    HLTError("header exists");
    fErrno=EEXIST;
  }
  return -fErrno;
}

// int AliHLTMemoryFile::WriteTrailerBuffer(const char* pTrailer, int size)
// {
//   // see header file for function documentation
//   fErrno=0;
//   if (fD>0) {
//     HLTError("file must be closed to write trailer");
//     return EPERM;
//   }
//   if (fSize+size<fBufferSize) {
//     memcpy(fpBuffer+fSize, pTrailer, size);
//   } else {
//     HLTError("no space left in memory file");
//     fErrno=ENOSPC;
//   }
//   return fErrno;
// }
