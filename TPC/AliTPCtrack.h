#ifndef ALITPCTRACK_H
#define ALITPCTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

/*****************************************************************************
 *                          December 18, 2000                                *
 *  Internal view of the TPC track parametrisation as well as the order of   *
 *           track parameters are subject for possible changes !             *
 *  Use GetExternalParameters() and GetExternalCovariance() to access TPC    *
 *      track information regardless of its internal representation.         *
 * This formation is now fixed in the following way:                         *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuth angle    *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *****************************************************************************/

#include <AliKalmanTrack.h>
#include <TMath.h>

const Double_t kConversionConstant=100/0.299792458/0.2; 

class AliTPCClustersArray;
class AliTPCcluster;

//_____________________________________________________________________________
class AliTPCtrack : public AliKalmanTrack {
public:
  AliTPCtrack():AliKalmanTrack(){}
  AliTPCtrack(UInt_t index, const Double_t xx[5], 
              const Double_t cc[15], Double_t xr, Double_t alpha); 
  AliTPCtrack(const AliTPCtrack& t);
  Int_t PropagateToVertex(
                    Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
  Int_t Rotate(Double_t angle);
  void CookLabel(AliTPCClustersArray *carray);
  void SetdEdx(Float_t dedx) {fdEdx=dedx;}

  Double_t GetX()     const {return fX;}
  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetdEdx()  const {return fdEdx;}

  Double_t GetY()   const {return fP0;}
  Double_t GetZ()   const {return fP1;}
  Double_t GetSnp() const {return fX*fP3 - fP2;}             
  Double_t Get1Pt() const {return fP3*kConversionConstant;}             
  Double_t GetTgl() const {return fP4;}

  Double_t GetSigmaY2() const {return fC00;}
  Double_t GetSigmaZ2() const {return fC11;}

  Int_t Compare(const TObject *o) const;

  void GetExternalParameters(Double_t& xr, Double_t x[5]) const ;
  void GetExternalCovariance(Double_t cov[15]) const ;

//******** To be removed next release !!! **************
  Double_t GetEta() const {return fP2;}
  Double_t GetC()   const {return fP3;}
  void GetCovariance(Double_t cc[15]) const {
    cc[0 ]=fC00;
    cc[1 ]=fC10;  cc[2 ]=fC11;
    cc[3 ]=fC20;  cc[4 ]=fC21;  cc[5 ]=fC22;
    cc[6 ]=fC30;  cc[7 ]=fC31;  cc[8 ]=fC32;  cc[9 ]=fC33;
    cc[10]=fC40;  cc[11]=fC41;  cc[12]=fC42;  cc[13]=fC43;  cc[14]=fC44;
  }  
//****************************************************** 

  void GetCluster(Int_t i, Int_t &sec, Int_t &row, Int_t &ncl) const;

  virtual Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Int_t PropagateTo(Double_t xr,
                    Double_t x0=28.94,Double_t rho=0.9e-3,Double_t pm=0.139);
  Int_t Update(const AliCluster* c, Double_t chi2, UInt_t i);

private: 
  Double_t fX;              // X-coordinate of this track (reference plane)
  Double_t fAlpha;          // Rotation angle the local (TPC sector)
                            // coordinate system and the global ALICE one.

  Double_t fdEdx;           // dE/dx 

  Double_t fP0;             // Y-coordinate of a track
  Double_t fP1;             // Z-coordinate of a track
  Double_t fP2;             // C*x0
  Double_t fP3;             // track curvature
  Double_t fP4;             // tangent of the track momentum dip angle

  Double_t fC00;                         // covariance
  Double_t fC10, fC11;                   // matrix
  Double_t fC20, fC21, fC22;             // of the
  Double_t fC30, fC31, fC32, fC33;       // track
  Double_t fC40, fC41, fC42, fC43, fC44; // parameters
 
  UInt_t fIndex[200];       // indices of associated clusters 

  ClassDef(AliTPCtrack,1)   // Time Projection Chamber reconstructed tracks
};


inline 
void AliTPCtrack::GetCluster(Int_t i,Int_t &sec,Int_t &row,Int_t &ncl) const {
  //return sector, pad row and the index of the i-th cluster of this track 
  Int_t index=fIndex[i];
  sec=(index&0xff000000)>>24; 
  row=(index&0x00ff0000)>>16; 
  ncl=(index&0x0000ffff)>>00;
}

inline 
void AliTPCtrack::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  //---------------------------------------------------------------------
  // This function return external TPC track representation
  //---------------------------------------------------------------------
     xr=fX;          
     x[0]=GetY(); x[1]=GetZ(); x[2]=GetSnp(); x[3]=GetTgl(); x[4]=Get1Pt();
}

#endif


