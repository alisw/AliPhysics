// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCHWCFData.cxx
/// @author Matthias Richter
/// @date   2011-08-04
/// @brief  Decoder methods for the HWCF format
///
#include "AliHLTTPCHWCFData.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCHWCFEmulator.h"
#include "AliHLTTPCTransform.h"
#include "AliRawDataHeader.h"
#include "TFile.h"
#include <memory>
#include <ostream>

ClassImp(AliHLTTPCHWCFData) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCHWCFData::AliHLTTPCHWCFData(int forceVersion)
 : fpBuffer(NULL)
 , fBufferSize(0)
 , fVersion(-1)
 , fForcedVersion(forceVersion)
 , fRCUTrailerSize(0)
 , fpFileBuffer(NULL)
 , fIterator()
 , fIteratorEnd()
{
  // constructor
}

const unsigned AliHLTTPCHWCFData::fgkAliHLTTPCHWClusterSize=sizeof(AliHLTTPCHWCFData::AliHLTTPCHWClusterV1)/sizeof(AliHLTUInt32_t);

AliHLTTPCHWCFData::~AliHLTTPCHWCFData()
{
  // destructor
  if (fpFileBuffer) delete fpFileBuffer;
}

int AliHLTTPCHWCFData::Init(const AliHLTUInt8_t* pBuffer, int bufferSize)
{
  // init the internal pointer
  Reset();
  if (!pBuffer || (unsigned)bufferSize<sizeof(AliHLTUInt32_t)) return -EINVAL;
  if (bufferSize%sizeof(AliHLTUInt32_t)) {
    HLTError("invalid buffer size %d, expecting multiple of 4", bufferSize);
    return -EINVAL;
  }
  const AliHLTUInt32_t* rcuTrailer=reinterpret_cast<const AliHLTUInt32_t*>(pBuffer);
  rcuTrailer+=bufferSize/sizeof(AliHLTUInt32_t)-1;

  if ((*rcuTrailer >> 30) != 3) {
    // The HWCFEmulator does not write the RCU trailer at the moment
    // skip this error
    // HLTError("can not read last rcu trailer word");
    // return -ENODATA;
  } else {
    // number of 32 bit RCU trailer words can be found in the last word
    int nofRCUTrailerWords = (*rcuTrailer & 0x7F);
    if (nofRCUTrailerWords < 2) {
      HLTError("Invalid trailer size found (%d bytes)", nofRCUTrailerWords*4);
      return -ENODATA;
    }

    if (nofRCUTrailerWords*4>bufferSize) {
      HLTError("inconsistent RCU trailer size, exceeds buffer size");
      return -ENODATA;
    }

    // check if the first RCU trailer word starts with pattern '10'
    rcuTrailer-=nofRCUTrailerWords-1;
    if ((*rcuTrailer >> 30) != 2) {
      HLTError("inconsistent first RCU trailer word: can not find indicator pattern '10' in bit 31 and 30 (0x%08x): trailer size %d word(s), buffer size %d byte", *rcuTrailer, nofRCUTrailerWords, bufferSize);
      return -ENODATA;
    }

    fRCUTrailerSize = nofRCUTrailerWords*4;
  }

  fpBuffer=pBuffer;
  fBufferSize=bufferSize;

  return 0;
}

int AliHLTTPCHWCFData::Reset()
{
  // reset
  fpBuffer=NULL;
  fBufferSize=0;
  fRCUTrailerSize=0;
  return 0;
}

int AliHLTTPCHWCFData::CheckVersion()
{
  // check the data buffer for format version
  if (fVersion>=0) return fVersion;
  if (CheckAssumption(kHWCFDataV1, fpBuffer, fBufferSize-fRCUTrailerSize)) {
    fVersion=kHWCFDataV1;
  } else if (CheckAssumption(kHWCFDataV0, fpBuffer, fBufferSize-fRCUTrailerSize)) {
    fVersion=kHWCFDataV0;
  }
  return fVersion;
}

bool AliHLTTPCHWCFData::CheckAssumption(int format, const AliHLTUInt8_t* pData, int size) const
{
  // check the format assumption for data buffer
  int elementsize=GetElementSize(format);
  // size has to be known
  if (elementsize<0) return false;
  // buffer must be multiple of element size
  if (size<elementsize || (size%elementsize)!=0) return false;
  for (int trial=0; trial<10 && (trial+1)*elementsize<=size; trial++) {
    AliHLTUInt32_t header=AliHLTTPCHWCFEmulator::ReadBigEndian(*reinterpret_cast<const AliHLTUInt32_t*>(pData+trial*elementsize));
    // cluster header starts with 11 in bit 30 and 31
    if ((header&0xc0000000)!=0xc0000000) return false;
    // check that the padrow is within bounds
    if (((header >> 24) & 0x3f)>(unsigned)AliHLTTPCTransform::GetNRows(-1)) return false;
  }
  return true;
}

Int_t AliHLTTPCHWCFData::GetNumberOfClusters() const
{
  // get number of clusters
  if (fVersion<0) return 0;
  int elementsize=GetElementSize(fVersion);
  if (elementsize<0) return 0;
  if (!fpBuffer || fBufferSize==0 || fBufferSize<fRCUTrailerSize) return 0;
  return (fBufferSize-fRCUTrailerSize)/elementsize;
}

Int_t    AliHLTTPCHWCFData::GetPadRow(int i)  const
{
  // get raw coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetPadRow();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetPadRow();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -1;
}

Float_t  AliHLTTPCHWCFData::GetPad(int i)     const
{
  // get pad coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetPad();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetPad();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -10000.;
}

Float_t  AliHLTTPCHWCFData::GetTime(int i)    const
{
  // get time coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetTime();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetTime();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -10000.;
}

Float_t  AliHLTTPCHWCFData::GetSigmaY2(int i) const
{
  // get sigmaY2 coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetSigmaY2();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetSigmaY2();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -10000.;
}

Float_t  AliHLTTPCHWCFData::GetSigmaZ2(int i) const
{
  // get sigmaZ2 coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetSigmaZ2();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetSigmaZ2();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -10000.;
}

Int_t    AliHLTTPCHWCFData::GetCharge(int i)  const
{
  // get charge coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetCharge();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetCharge();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -1;
}

Int_t    AliHLTTPCHWCFData::GetQMax(int i)    const
{
  // get qmax coordinate
  if (fVersion>=0 && CheckBounds(i)) {
    switch (fVersion) {
    case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(Get(i))->GetQMax();
    case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(Get(i))->GetQMax();
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
  }
  return -1;
}

void AliHLTTPCHWCFData::Print(const char* option)
{
  // print info
  cout << "HWCF format version " << fVersion << endl;
  if (fVersion<0) return;

  cout << "    " << GetNumberOfClusters() << " cluster(s)" << endl;
  if (GetRCUTrailerSize()>0) {
    cout << "    RCU trailer: " << GetRCUTrailerSize() << " word(s)" << endl;
  }

  if (strcmp(option, "all")==0) {
    for (unsigned i=0; (int)i<GetNumberOfClusters(); i++) {
      cout << /* setw(5) <<*/ i << ":";
      cout << /* setw(8) <<*/ " " << GetPadRow(i);
      cout << /* setw(8) <<*/ " " << GetPad(i);
      cout << /* setw(8) <<*/ " " << GetTime(i);
      cout << /* setw(8) <<*/ " " << GetSigmaY2(i);
      cout << /* setw(8) <<*/ " " << GetSigmaZ2(i);
      cout << /* setw(8) <<*/ " " << GetCharge(i);
      cout << /* setw(8) <<*/ " " << GetQMax(i);
      cout << endl;
    }
  }
}

int AliHLTTPCHWCFData::Open(const char* filename)
{
  // open block from file and add to collection
  if (!filename) return -EINVAL;
  
  TString input=filename;
  input+="?filetype=raw";
  std::auto_ptr<TFile> pFile(new TFile(input));
  if (!pFile.get()) return -ENOMEM;
  if (pFile->IsZombie()) return -ENOENT;

  int iResult=0;
  if (pFile->GetSize()<(int)sizeof(AliRawDataHeader)) {
    HLTError("file %s to small", filename);
    return -ENODATA;
  }

  pFile->Seek(0);
  std::auto_ptr<TArrayC> buffer(new TArrayC);
  if (!buffer.get()) return -ENOMEM;

  buffer->Set(pFile->GetSize());
  if (pFile->ReadBuffer(buffer->GetArray(), buffer->GetSize())==0) {
  } else {
    HLTError("failed reading %d byte(s) from file %s", pFile->GetSize(), filename);
    iResult=-ENODATA;
  }

  AliHLTUInt8_t* pBuffer=reinterpret_cast<AliHLTUInt8_t*>(buffer->GetArray()+sizeof(AliRawDataHeader));
  unsigned bufferSize=buffer->GetSize()-sizeof(AliRawDataHeader);
  if ((iResult=Init(pBuffer, bufferSize))<0 ||
      (iResult=CheckVersion())<0) {
    Reset();
    return iResult;
  }

  fpFileBuffer=buffer.release();
  return GetNumberOfClusters();
}

Int_t    AliHLTTPCHWCFData::AliHLTTPCHWClusterV0::GetPadRow()  const
{
  // bit 24 to 29
  AliHLTUInt32_t header=AliHLTTPCHWCFEmulator::ReadBigEndian(fHeader);
  if ((header>>30) != 3) return -EBADMSG;
  return (header >> 24) & 0x3f;
}

Int_t    AliHLTTPCHWCFData::AliHLTTPCHWClusterV0::GetCharge()  const
{
  // 24 bit fixed point number with 6 bits after the point
  AliHLTUInt32_t header=AliHLTTPCHWCFEmulator::ReadBigEndian(fHeader);
  return (header & 0xFFFFFF )>>6;
}

Int_t    AliHLTTPCHWCFData::AliHLTTPCHWClusterV1::GetPadRow()  const
{
  // bit 24 to 29
  AliHLTUInt32_t header=AliHLTTPCHWCFEmulator::ReadBigEndian(fHeader);
  if ((header>>30) != 3) return -EBADMSG;
  return (header >> 24) & 0x3f;
}

Int_t    AliHLTTPCHWCFData::AliHLTTPCHWClusterV1::GetCharge()  const
{
  // 32 bit fixed point number with 12 bits after the point
  return AliHLTTPCHWCFEmulator::ReadBigEndian(fCharge)>>12;
}

Int_t    AliHLTTPCHWCFData::AliHLTTPCHWClusterV1::GetQMax()    const
{
  // 24 bit fixed point number with 12 bits after the point
  AliHLTUInt32_t header=AliHLTTPCHWCFEmulator::ReadBigEndian(fHeader);
  return (header & 0xFFFFFF )>>12;
}
