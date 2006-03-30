#ifndef ALIITSTRACKV2_H
#define ALIITSTRACKV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                       ITS Track Class
//
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------


/*****************************************************************************
 *                          December 18, 2000                                *
 *  Internal view of the ITS track parametrisation as well as the order of   *
 *           track parameters are subject for possible changes !             *
 *  Use GetExternalParameters() and GetExternalCovariance() to access ITS    *
 *      track information regardless of its internal representation.         *
 * This formation is now fixed in the following way:                         *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *****************************************************************************/

#include <AliKalmanTrack.h>

#include "AliITSrecoV2.h"

class AliESDtrack;
class AliITSStrLine;

//_____________________________________________________________________________
class AliITStrackV2 : public AliKalmanTrack {
public:
  AliITStrackV2();
  AliITStrackV2(AliESDtrack& t,Bool_t c=kFALSE) throw (const Char_t *);
  AliITStrackV2(const AliITStrackV2& t);
  Int_t PropagateToVertex(Double_t d=0., Double_t x0=0.);
  Int_t Propagate(Double_t alpha, Double_t xr);
  virtual Int_t CorrectForMaterial(Double_t d, Double_t x0=21.82);
  Int_t PropagateTo(Double_t xr, Double_t d, Double_t x0=21.82);
  Double_t PropagateToDCA(AliKalmanTrack *p, Double_t d=0., Double_t x0=0.); 
  Int_t Update(const AliCluster* cl,Double_t chi2,UInt_t i);
  Int_t Improve(Double_t x0,Double_t xyz[3],Double_t ers[3]);
  void SetdEdx(Double_t dedx) {fdEdx=dedx;}
  void SetSampledEdx(Float_t q, Int_t i);
  void CookdEdx(Double_t low=0., Double_t up=0.51);
  void SetDetectorIndex(Int_t i) {SetLabel(i);}
  void ResetCovariance();
  void ResetClusters() { SetChi2(0.); SetNumberOfClusters(0); }
  void UpdateESDtrack(ULong_t flags) const;
  
  AliESDtrack *GetESDtrack() const {return fESDtrack;}

  Int_t GetDetectorIndex() const {return GetLabel();}
  Double_t GetX()    const {return fX;}
  Double_t GetAlpha()const {return fAlpha;}
  Double_t GetdEdx() const {return fdEdx;}
  Double_t GetPIDsignal() const {return GetdEdx();}
  Double_t GetY()    const {return fP0;}
  Double_t GetZ()    const {return fP1;}
  Double_t GetSnp()  const {return fP2;}
  Double_t GetTgl()  const {return fP3;}
  Double_t GetC()    const {return fP4;}
  Double_t Get1Pt() const {
      return (TMath::Sign(1e-9,fP4) + fP4)*GetLocalConvConst();
  }
  Double_t GetD(Double_t x=0, Double_t y=0) const;
  Double_t GetZat(Double_t x=0) const;

  Double_t GetSigmaY2() const {return fC00;}
  Double_t GetSigmaZ2() const {return fC11;}
  Int_t Compare(const TObject *o) const;
  void GetExternalParameters(Double_t& xr, Double_t x[5]) const ;
  void GetExternalCovariance(Double_t cov[15]) const ;
  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}
  Int_t GetGlobalXYZat(Double_t r,Double_t &x,Double_t &y,Double_t &z) const;
  void ApproximateHelixWithLine(Double_t xk, AliITSStrLine *line);
  Double_t GetPredictedChi2(const AliCluster *cluster) const;
   Int_t Invariant() const;
 
protected:
  void GetXYZ(Float_t r[3]) const;

  Double_t fX;              // X-coordinate of this track (reference plane)
  Double_t fAlpha;          // rotation angle

  Double_t fdEdx;           // dE/dx

  Double_t fP0;             // Y-coordinate of a track 
  Double_t fP1;             // Z-coordinate of a track
  Double_t fP2;             // sine of the track momentum azimuthal angle
  Double_t fP3;             // tangent of the track momentum dip angle
  Double_t fP4;             // track curvature

  Double_t fC00;                         // covariance
  Double_t fC10, fC11;                   // matrix
  Double_t fC20, fC21, fC22;             // of the
  Double_t fC30, fC31, fC32, fC33;       // track
  Double_t fC40, fC41, fC42, fC43, fC44; // parameters 

  UInt_t fIndex[kMaxLayer]; // indices of associated clusters 
  Float_t fdEdxSample[4];   // array of dE/dx samples b.b.

  AliESDtrack *fESDtrack;   //! pointer to the connected ESD track
  ClassDef(AliITStrackV2,4)   //ITS reconstructed track
};

inline 
void AliITStrackV2::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  //---------------------------------------------------------------------
  // This function return external ITS track representation
  //---------------------------------------------------------------------
     xr=fX;          
     x[0]=GetY(); x[1]=GetZ(); x[2]=GetSnp(); x[3]=GetTgl(); x[4]=Get1Pt();
}

inline
void AliITStrackV2::SetSampledEdx(Float_t q, Int_t i) {
  //----------------------------------------------------------------------
  // This function stores dEdx sample corrected for the track segment length 
  // Origin: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
  //----------------------------------------------------------------------
  if (i<0) return;
  if (i>3) return;
  Double_t s=GetSnp(), t=GetTgl();
  q *= TMath::Sqrt((1-s*s)/(1+t*t));
  fdEdxSample[i]=q;
}

inline void AliITStrackV2::GetXYZ(Float_t r[3]) const {
  //---------------------------------------------------------------------
  // Returns the position of the track in the global coord. system 
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  r[0]=fX*cs - fP0*sn; r[1]=fX*sn + fP0*cs; r[2]=fP1;
}

#endif


