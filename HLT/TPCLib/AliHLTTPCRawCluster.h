// $Id$
#ifndef ALIHLTTPCRAWCLUSTER_H
#define ALIHLTTPCRAWCLUSTER_H

#include "Rtypes.h"

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
    , fSigmaPad2(0.)
    , fSigmaTime2(0.)
    , fCharge(0)
    , fQMax(0)
    , fFlags(0)
  {}

  AliHLTTPCRawCluster(short PadRow,
		      float Pad,
		      float Time,
		      float SigmaPad2,
		      float SigmaTime2,
		      unsigned short Charge,
		      unsigned short QMax
		      )
    : fPadRow(PadRow)
    , fPad(Pad)
    , fTime(Time)
    , fSigmaPad2(SigmaPad2)
    , fSigmaTime2(SigmaTime2)
    , fCharge(Charge)
    , fQMax(QMax)
    , fFlags(0)
  {}

  AliHLTTPCRawCluster(const AliHLTTPCRawCluster& other)
    : fPadRow(other.fPadRow)
    , fPad(other.fPad)
    , fTime(other.fTime)
    , fSigmaPad2(other.fSigmaPad2)
    , fSigmaTime2(other.fSigmaTime2)
    , fCharge(other.fCharge)
    , fQMax(other.fQMax)
    , fFlags(other.fFlags)
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
  unsigned short fFlags; //Flags: (1 << 0): Split in pad direction
                         //       (1 << 1): Split in time direction
  float fPad;
  float fTime;
  float fSigmaPad2;
  float fSigmaTime2;
  unsigned short fCharge;
  unsigned short fQMax;

  Int_t   GetPadRow()  const {return fPadRow;}
  Float_t GetPad()     const {return fPad;}
  Float_t GetTime()    const {return fTime;}
  Float_t GetSigmaPad2() const {return fSigmaPad2;}
  Float_t GetSigmaTime2() const {return fSigmaTime2;}
  Int_t   GetCharge()  const {return fCharge;}
  Int_t   GetQMax()    const {return fQMax;}
  Bool_t  GetFlagSplitPad() const {return (fFlags & (1 << 0));}
  Bool_t  GetFlagSplitTime() const {return (fFlags & (1 << 1));}
  Bool_t  GetFlagSplitAny() const {return (fFlags & 3);}

  void SetPadRow(Short_t padrow)  {fPadRow=padrow;}
  void SetPad(Float_t pad)     {fPad=pad;}
  void SetTime(Float_t time)    {fTime=time;}
  void SetSigmaPad2(Float_t sigmaPad2) {fSigmaPad2=sigmaPad2;}
  void SetSigmaTime2(Float_t sigmaTime2) {fSigmaTime2=sigmaTime2;}
  void SetCharge(UShort_t charge)  {fCharge=charge;}
  void SetQMax(UShort_t qmax)    {fQMax=qmax;}

  void ClearFlags() {fFlags = 0;}
  void SetFlagSplitPad() {fFlags |= (1 << 0);}
  void SetFlagSplitTime() {fFlags |= (1 << 1);}
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
