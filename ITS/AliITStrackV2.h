#ifndef ALIITSTRACKV2_H
#define ALIITSTRACKV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                       ITS Track Class
//
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
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

class AliTPCtrack;

//_____________________________________________________________________________
class AliITStrackV2  : public AliKalmanTrack {
public:
  AliITStrackV2():AliKalmanTrack(){}
  AliITStrackV2(const AliTPCtrack& t) throw (const Char_t *);
  AliITStrackV2(const AliITStrackV2& t);
  Int_t 
  PropagateToVertex(Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
  void SetdEdx(Float_t dedx) {fdEdx=dedx;}
  void SetDetectorIndex(Int_t i) {SetLabel(i);}
  
  void *operator new(size_t s,void *p) { return p; }
  void *operator new(size_t s) { return ::operator new(s); }
  
  Int_t GetDetectorIndex() const {return GetLabel();}
  Double_t GetX()    const {return fX;}
  Double_t GetAlpha()const {return fAlpha;}
  Float_t  GetdEdx() const {return fdEdx;}

  Double_t GetY()    const {return fP0;}
  Double_t GetZ()    const {return fP1;}
  Double_t GetSnp()  const {return fP2;}
  Double_t GetTgl()  const {return fP3;}
  Double_t Get1Pt()  const {return fP4*kConvConst;}


  Double_t GetD() const;


  Double_t GetSigmaY2() const {return fC00;}
  Double_t GetSigmaZ2() const {return fC11;}

  Int_t Compare(const TObject *o) const;

  void GetExternalParameters(Double_t& xr, Double_t x[5]) const ;
  void GetExternalCovariance(Double_t cov[15]) const ;

  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}
  Int_t GetGlobalXYZat(Double_t r,Double_t &x,Double_t &y,Double_t &z) const;

  Int_t Propagate(Double_t alpha,
                  Double_t xr,Double_t x0,Double_t rho,Double_t pm=0.139);

  Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Int_t Update(const AliCluster* cl,Double_t chi2,UInt_t i);

  Double_t GetPredictedChi2(const AliCluster *cluster, Double_t *m,
                  Double_t x0, Double_t pm=0.139) const;
  Int_t Update(const Double_t *m, Double_t chi2, UInt_t i);
  Int_t Improve(Double_t x0,Double_t yv,Double_t zv);

  Int_t Invariant() const;

  //protected:
Int_t 
PropagateTo(Double_t xr,Double_t x0=21.82,Double_t rho=2.33,Double_t pm=0.139);
 
private:
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

  ClassDef(AliITStrackV2,1)   //ITS reconstructed track
};

inline 
void AliITStrackV2::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  //---------------------------------------------------------------------
  // This function return external TPC track representation
  //---------------------------------------------------------------------
     xr=fX;          
     x[0]=GetY(); x[1]=GetZ(); x[2]=GetSnp(); x[3]=GetTgl(); x[4]=Get1Pt();
}

#endif


