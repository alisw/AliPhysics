//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCHWCFDATA_H
#define ALIHLTTPCHWCFDATA_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCHWCFData.h
/// @author Matthias Richter
/// @date   2011-08-04
/// @brief  Decoder methods for the HWCF format
///

#include "AliHLTTPCRootTypes.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTErrorGuard.h"

class TArrayC;

/**
 * @class AliHLTTPCHWCFData
 * The class provides decoding functionality for the output format of the
 * TPC HW ClusterFinder
 *
 * Two formats have been defined in the past and can be detected:
 * version 0: 5 32bit words 20 Byte
 * <pre>
 *   word 0: header (big endian 32bit unsigned)
 *           bit 31-30: 0x11 indicates cluster
 *           bit 29-24: row number in partition
 *           bit 23-0: Qtot, fixed point number with 6 bits after the point
 *   word 1: pad (float)
 *   word 2: time (float)
 *   word 3: pad variance (float)
 *   word 4: time variance (float)
 * </pre>
 *
 * version 1: 6 32bit words 24 Byte
 * <pre>
 *   word 0: header (big endian 32bit unsigned)
 *           bit 31-30: 0x11 indicates cluster
 *           bit 29-24: row number in partition
 *           bit 23-0: Qmax, fixed point number with 6 bits after the point
 *   word 1: total charge 32bit  big endian, fixed point number with 12 bits after the point
 *   word 2: pad (float)
 *   word 3: time (float)
 *   word 4: pad variance (float)
 *   word 5: time variance (float)
 * </pre>
 */
class AliHLTTPCHWCFData : public AliHLTLogging {
 public:
  AliHLTTPCHWCFData(int forceVersion=-1);
  virtual ~AliHLTTPCHWCFData();

  int Init(const AliHLTUInt8_t* pBuffer, int bufferSize);
  int Reset();

  Int_t    GetNumberOfClusters() const;

  Int_t    GetPadRow(int i)  const;
  Float_t  GetPad(int i)     const;
  Float_t  GetTime(int i)    const;
  Float_t  GetSigmaY2(int i) const;
  Float_t  GetSigmaZ2(int i) const;
  Int_t    GetCharge(int i)  const;
  Int_t    GetQMax(int i)    const;

  int CheckVersion();
  bool CheckAssumption(int format, const AliHLTUInt8_t* pData, int size) const;

  // check if index is within bounds
  bool CheckBounds(int i) const {
    if (fVersion<0) {
      ALIHLTERRORGUARD(1, "");
      return false;
    }
    int elementsize=GetElementSize(fVersion);
    if (elementsize<0) return false;
    return ((i+1)*elementsize+fRCUTrailerSize<=fBufferSize);
  }

  // get the size of one element
  int GetElementSize(int version) const {
    switch (version) {
    case 0: return sizeof(AliHLTTPCHWClusterV0);
    case 1: return sizeof(AliHLTTPCHWClusterV1);
    default:
      ALIHLTERRORGUARD(1, "invalid format version %d", fVersion);
    }
    return -1;
  }

  // pointer to RCU trailer
  const AliHLTUInt8_t*  GetRCUTrailer() const
  {
    if (fRCUTrailerSize<=0 || fpBuffer==NULL || fBufferSize<fRCUTrailerSize) return NULL;
    return fpBuffer+(fBufferSize-fRCUTrailerSize);
  }

  // size of RCU trailer
  int GetRCUTrailerSize() const { return fRCUTrailerSize; }

  // print info
  void Print(const char* option);

  // open a file and init
  int Open(const char* filename);

  enum {
    kHWCFDataV0 = 0,
    kHWCFDataV1 = 1,
  };

  struct AliHLTTPCHWClusterV0 {
    AliHLTUInt32_t fHeader;
    Float_t        fPad;
    Float_t        fTime;
    Float_t        fSigmaY2;
    Float_t        fSigmaZ2;

    Int_t    GetPadRow()  const;
    Float_t  GetPad()     const {return fPad+0.5;}
    Float_t  GetTime()    const {return fTime;}
    Float_t  GetSigmaY2() const {
      Float_t sy2 = fSigmaY2 - fPad*fPad;
      return (sy2>0) ?sy2 :0.;
    }
    Float_t  GetSigmaZ2() const {
      Float_t sz2 = fSigmaZ2 - fTime*fTime;
      return (sz2>0) ?sz2 :0.;
    }
    Int_t    GetCharge()  const;
    Int_t    GetQMax()    const {return -1;}
  };

  struct AliHLTTPCHWClusterV1 {
    AliHLTUInt32_t fHeader;
    AliHLTUInt32_t fCharge;
    Float_t        fPad;
    Float_t        fTime;
    Float_t        fSigmaY2;
    Float_t        fSigmaZ2;

    Int_t    GetPadRow()  const;
    Float_t  GetPad()     const {return fPad+0.5;}
    Float_t  GetTime()    const {return fTime;}
    Float_t  GetSigmaY2() const {
      Float_t sy2 = fSigmaY2 - fPad*fPad;
      return (sy2>0) ?sy2 :0.;
    }
    Float_t  GetSigmaZ2() const {
      Float_t sz2 = fSigmaZ2 - fTime*fTime;
      return (sz2>0) ?sz2 :0.;
    }
    Int_t    GetCharge()  const;
    Int_t    GetQMax()    const;
  };

  template<typename T>
  class AliHLTTPCHWClusterDecoder {
  public:
    AliHLTTPCHWClusterDecoder(const T* pClusterArray, int entries);
    ~AliHLTTPCHWClusterDecoder();

    // i'th element, no bounds check for performance reasons
    const T& operator[](unsigned i) {
      return fpClusterArray[i];
    }

    Int_t    GetPadRow(int i)  const {return fpClusterArray[i]->GetPadRow();}
    Float_t  GetPad(int i)     const {return fpClusterArray[i]->GetPad();}
    Float_t  GetTime(int i)    const {return fpClusterArray[i]->GetTime();}
    Float_t  GetSigmaY2(int i) const {return fpClusterArray[i]->GetSigmaY2();}
    Float_t  GetSigmaZ2(int i) const {return fpClusterArray[i]->GetSigmaZ2();}
    Int_t    GetCharge(int i)  const {return fpClusterArray[i]->GetCharge();}
    Int_t    GetQMax(int i)    const {return fpClusterArray[i]->GetQMax();}

  private:
    const T* fpClusterArray; //! array of clusters
    int fEntries;            //! number of entries
  };

  class iterator {
  public:
    iterator()
      : fData(NULL), fVersion(-1), fElementSize(0) {}
    iterator(const AliHLTUInt8_t* pData, int version, int elementSize)
      : fData(pData), fVersion(version), fElementSize(elementSize) {}
    iterator(const iterator& i)
      : fData(i.fData), fVersion(i.fVersion), fElementSize(i.fElementSize) {}
    iterator& operator=(const iterator& i) {
      if (this==&i) return *this;
      fData=i.fData; fVersion=i.fVersion; fElementSize=i.fElementSize; return *this;
    }
    ~iterator() {fData=NULL;}

    bool operator==(const iterator& i) const  {return (fData!=NULL) && (fData==i.fData);}
    bool operator!=(const iterator& i) const  {return (fData!=NULL) && (fData!=i.fData);}
    // prefix operators
    iterator& operator++() {fData+=fElementSize; return *this;}
    iterator& operator--() {fData-=fElementSize; return *this;}
    // postfix operators
    iterator operator++(int) {iterator i(*this); fData+=fElementSize; return i;}
    iterator operator--(int) {iterator i(*this); fData-=fElementSize; return i;}

    iterator& operator+=(int step) {fData+=step*fElementSize; return *this;}

    Int_t    GetPadRow()  const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetPadRow();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetPadRow();
      } return -1;
    }
    Float_t  GetPad()     const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetPad();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetPad();
      } return -10000.;
    }
    Float_t  GetTime()    const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetTime();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetTime();
      } return -10000.;
    }
    Float_t  GetSigmaY2() const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetSigmaY2();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetSigmaY2();
      } return -10000.;
    }
    Float_t  GetSigmaZ2() const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetSigmaZ2();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetSigmaZ2();
      } return -10000.;
    }
    Int_t    GetCharge()  const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetCharge();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetCharge();
      } return -1;
    }
    Int_t    GetQMax()    const {
      switch (fVersion) {
      case 0: return reinterpret_cast<const AliHLTTPCHWClusterV0*>(fData)->GetQMax();
      case 1: return reinterpret_cast<const AliHLTTPCHWClusterV1*>(fData)->GetQMax();
      } return -1;
    }

  protected:
  private:
    const AliHLTUInt8_t* fData; //! data
    int fVersion; //! format version
    int fElementSize; //! element size
  };

  // prepare iterator and end marker
  iterator& begin() {
    fIterator.~iterator();
    new (&fIterator) iterator(fpBuffer, fVersion, GetElementSize(fVersion));
    fIteratorEnd=fIterator;
    fIteratorEnd+=GetNumberOfClusters();
    return fIterator;
  }

  // get loop end marker
  iterator& end() {
    return fIteratorEnd;
  }

  // find one single element
  iterator& find(int i) {
    fIterator.~iterator();
    fIteratorEnd.~iterator();
    if (i>=GetNumberOfClusters()) return fIterator;
    new (&fIterator) iterator(fpBuffer, fVersion, GetElementSize(fVersion));
    fIterator+=i;
    fIteratorEnd=fIterator;
    fIteratorEnd+=1;
    return fIterator;
  }

  static const unsigned  fgkAliHLTTPCHWClusterSize;
 protected:

 private:
  AliHLTTPCHWCFData(const AliHLTTPCHWCFData&);
  AliHLTTPCHWCFData& operator=(const AliHLTTPCHWCFData&);

  // get pointer to i'th element
  const AliHLTUInt8_t* Get(int i) const
  {
    if (!fpBuffer) return NULL;
    int elementsize=GetElementSize(fVersion);
    if (elementsize<0) return NULL;
    return fpBuffer+(i*elementsize);
  }


  const AliHLTUInt8_t* fpBuffer; //! pointer to data buffer
  int fBufferSize; //! size of data buffer

  int fVersion; //! format version
  int fForcedVersion; //! forced format version
  int fRCUTrailerSize; //! size of the RCU trailer in Byte

  TArrayC* fpFileBuffer; //! internal buffer for file content

  iterator fIterator; //! iterator
  iterator fIteratorEnd; //! end iterator

  ClassDef(AliHLTTPCHWCFData, 0)
};
#endif
