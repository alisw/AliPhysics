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
  /* UInt_t fP0;       // First  4 bytes of the packed data */
  /* UInt_t fP1;       // Second 4 bytes of the packed data */
  /* UInt_t fP2;       // Third  4 bytes of the packed data */
  /* UInt_t fP3;       // Last   4 bytes of the packed data */
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

  /* Int_t GetPadRow()    const {return        (fP0>>24) &     0xff;} */
  /* Float_t GetPad()     const {return (float( fP0      & 0xffffff)-8388608.)*1.e-4;} */
  /* Float_t GetTime()    const {return (float( fP1      & 0xffffff)-8388608.)*1.e-4;} */
  /* Float_t GetSigmaY2() const {return (float( fP2      &   0xffff)-32768.)*1.e-4;} */
  /* Float_t GetSigmaZ2() const {return (float((fP2>>16) &   0xffff)-32768.)*1.e-4;} */
  /* Int_t   GetCharge()  const {return         fP3      &  0xfffff;} */
  /* Int_t   GetQMax()    const {return        (fP3>>20) &    0x3ff;} */

  /* void SetPadRow(Int_t padrow)       { fP0&=0x00ffffff; fP0|=(padrow & 0xff)<<24;} */
  /* void SetPad(Float_t pad)         { fP0&=0xff000000; float v=pad*1.e4+8388608.; */
  /*                                       if( v<0 ) v=0; else if( v>0x00FFFFFF ) v=0x00FFFFFF; */
  /* 					UInt_t iv=(UInt_t) v; fP0|=iv; */
  /* } */
  /* void SetTime(Float_t time)       { fP1&=0xff000000; float v=time*1.e4+8388608.; */
  /*                                       if( v<0 ) v=0; else if( v>0x00FFFFFF ) v=0x00FFFFFF; */
  /* 					UInt_t iv=(UInt_t) v; fP1|=iv; */
  /* } */
  /* void SetSigmaY2(Float_t sigmaY2) { fP2&=0xffff0000; float v=sigmaY2*1.e4+32768.; */
  /*                                       if( v<0 ) v=0; else if( v>0x0000FFFF ) v=0x0000FFFF; */
  /* 					UInt_t iv=(UInt_t) v; fP2|=iv; */
  /* } */
  /* void SetSigmaZ2(Float_t sigmaZ2) { fP2&=0x0000ffff; float v=sigmaZ2*1.e4+32768.; */
  /*                                       if( v<0 ) v=0; else if( v>0x0000FFFF ) v=0x0000FFFF; */
  /* 					UInt_t iv=(UInt_t) v; fP2|=iv<<16; */
  /* } */
  /* void SetCharge(Int_t charge)     { fP3&=0xfff00000; fP3|=charge&0xfffff;} */
  /* void SetQMax(Int_t qmax)         { fP3&=0x000fffff; fP3|=(qmax&0x3ff)<<20;} */
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
