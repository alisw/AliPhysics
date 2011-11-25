// $Id$
#ifndef ALIHLTTPCRAWCLUSTER_H
#define ALIHLTTPCRAWCLUSTER_H

#include "AliHLTTPCRootTypes.h"

/**
 * @struct AliHLTTPCRawCluster
 * Primitive data of a TPC cluster in raw coordinates. The plan is to store the
 * data in a compressed format by limiting the resolution of the float values.
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTTPCRawCluster {
  AliHLTTPCRawCluster()
    : fPadRow(0)
    , fPad(0.)
    , fTime(0.)
    , fSigmaY2(0.)
    , fSigmaZ2(0.)
    , fCharge(0)
    , fQMax(0)
  {}

  AliHLTTPCRawCluster(short PadRow,
		      float Pad,
		      float Time,
		      float SigmaY2,
		      float SigmaZ2,
		      unsigned short Charge,
		      unsigned short QMax
		      )
    : fPadRow(PadRow)
    , fPad(Pad)
    , fTime(Time)
    , fSigmaY2(SigmaY2)
    , fSigmaZ2(SigmaZ2)
    , fCharge(Charge)
    , fQMax(QMax)
  {}

  AliHLTTPCRawCluster(const AliHLTTPCRawCluster& other)
    : fPadRow(other.fPadRow)
    , fPad(other.fPad)
    , fTime(other.fTime)
    , fSigmaY2(other.fSigmaY2)
    , fSigmaZ2(other.fSigmaZ2)
    , fCharge(other.fCharge)
    , fQMax(other.fQMax)
  {}

  AliHLTTPCRawCluster& operator=(const AliHLTTPCRawCluster& other) {
    if (this==&other) return *this;
    this->~AliHLTTPCRawCluster();
    new (this) AliHLTTPCRawCluster(other);
    return *this;
  }

  void Clear() {
    this->~AliHLTTPCRawCluster();
    new (this) AliHLTTPCRawCluster;
  }

  short fPadRow;
  float fPad;
  float fTime;
  float fSigmaY2;
  float fSigmaZ2;
  unsigned short fCharge;
  unsigned short fQMax;

  Int_t   GetPadRow()  const {return fPadRow;}
  Float_t GetPad()     const {return fPad;}
  Float_t GetTime()    const {return fTime;}
  Float_t GetSigmaY2() const {return fSigmaY2;}
  Float_t GetSigmaZ2() const {return fSigmaZ2;}
  Int_t   GetCharge()  const {return fCharge;}
  Int_t   GetQMax()    const {return fQMax;}

  void SetPadRow(Short_t padrow)  {fPadRow=padrow;}
  void SetPad(Float_t pad)     {fPad=pad;}
  void SetTime(Float_t time)    {fTime=time;}
  void SetSigmaY2(Float_t sigmaY2) {fSigmaY2=sigmaY2;}
  void SetSigmaZ2(Float_t sigmaZ2) {fSigmaZ2=sigmaZ2;}
  void SetCharge(UShort_t charge)  {fCharge=charge;}
  void SetQMax(UShort_t qmax)    {fQMax=qmax;}
};
typedef struct AliHLTTPCRawCluster AliHLTTPCRawCluster;

struct AliHLTTPCRawClusterData
{
  UInt_t fVersion; // version number
  UInt_t fCount;   // number of clusters
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCRawCluster  fClusters[1]; // array of clusters  
#else
  AliHLTTPCRawCluster  fClusters[0]; // array of clusters 
#endif
};
typedef struct AliHLTTPCRawClusterData AliHLTTPCRawClusterData;

#endif
