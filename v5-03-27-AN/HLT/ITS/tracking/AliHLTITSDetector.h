// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTITSDETECTOR_H
#define ALIHLTITSDETECTOR_H

#include "AliITSRecPoint.h"
#include "TMath.h"

class AliITSDetTypeRec;

/**
 * @class AliHLTITSDetector
 *
 * The AliHLTTPCCAHit class is the internal representation
 * of the TPC clusters for the AliHLTTPCCATracker algorithm.
 *
 */
class AliHLTITSDetector 
{ 
 public:
  AliHLTITSDetector():fR(0),fRmisal(0),fPhi(0),fSinPhi(0),fCosPhi(0),fYmin(0),fYmax(0),fZmin(0),fZmax(0){}
  AliHLTITSDetector(Double_t r,Double_t phi):fR(r),fRmisal(r),fPhi(phi),fSinPhi(TMath::Sin(phi)),fCosPhi(TMath::Cos(phi)),fYmin(10000),fYmax(-1000),fZmin(10000),fZmax(-1000){}
  ~AliHLTITSDetector() {}
  inline void GetGlobalXYZ( const AliITSRecPoint *cl, Double_t xyz[3]) const;
  Double_t GetR()   const {return fR;}
  Double_t GetRmisal()   const {return fRmisal;}
  Double_t GetPhi() const {return fPhi;}
  Double_t GetYmin() const {return fYmin;}
  Double_t GetYmax() const {return fYmax;}
  Double_t GetZmin() const {return fZmin;}
  Double_t GetZmax() const {return fZmax;}  
  void SetRmisal(Double_t rmisal) {fRmisal = rmisal;}
  void SetYmin(Double_t min) {fYmin = min;}
  void SetYmax(Double_t max) {fYmax = max;}
  void SetZmin(Double_t min) {fZmin = min;}
  void SetZmax(Double_t max) {fZmax = max;}

 private:

  AliHLTITSDetector(const AliHLTITSDetector& det);
  AliHLTITSDetector & operator=(const AliHLTITSDetector& det){
    this->~AliHLTITSDetector();new(this) AliHLTITSDetector(det);
    return *this;}
  Double_t fR;    // polar coordinates: r 
  Double_t fRmisal;    // polar coordinates: r, with misalignment 
  Double_t fPhi;  // polar coordinates: phi
  Double_t fSinPhi; // sin of phi;
  Double_t fCosPhi; // cos of phi 
  Double_t fYmin;   //  local y minimal
  Double_t fYmax;   //  local max y
  Double_t fZmin;   //  local z min
  Double_t fZmax;   //  local z max  
};

inline void  AliHLTITSDetector::GetGlobalXYZ(const AliITSRecPoint *cl, Double_t xyz[3]) const
{
  //
  // get cluster coordinates in global cooordinate 
  //
  xyz[2] = cl->GetZ();
  xyz[0] = fR*fCosPhi - cl->GetY()*fSinPhi;
  xyz[1] = fR*fSinPhi + cl->GetY()*fCosPhi;
}

#endif
