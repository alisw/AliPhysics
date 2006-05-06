#ifndef ALITPCTRACK_H
#define ALITPCTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

//-------------------------------------------------------
//                    TPC Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//
//                        December 18, 2000                                
//  Internal view of the TPC track parametrisation as well as the order of   
//           track parameters are subject for possible changes !             
//  Use GetExternalParameters() and GetExternalCovariance() to access TPC    
//      track information regardless of its internal representation.         
// This formation is now fixed in the following way:                         
//      external param0:   local Y-coordinate of a track (cm)                
//      external param1:   local Z-coordinate of a track (cm)                
//      external param2:   local sine of the track momentum azimuth angle    
//      external param3:   tangent of the track momentum dip angle           
//      external param4:   1/pt (1/(GeV/c))                                  
//
#include <AliKalmanTrack.h>
#include <TMath.h>

#include "AliTPCreco.h"
#include "AliExternalTrackParam.h"
class AliESDtrack;

//_____________________________________________________________________________
class AliTPCtrack : public AliKalmanTrack {
public:
  AliTPCtrack();
  AliTPCtrack(UInt_t index, const Double_t xx[5], 
              const Double_t cc[15], Double_t xr, Double_t alpha); 
  AliTPCtrack(const AliESDtrack& t);
  AliTPCtrack(const AliTPCtrack& t);
  virtual ~AliTPCtrack() {}
  Int_t PropagateToVertex(Double_t x0=36.66,Double_t rho=1.2e-3);
  Int_t Rotate(Double_t angle);
  void SetdEdx(Double_t dedx) {fdEdx=dedx;}

  Double_t GetYat(Double_t x) const ;
  Double_t GetX()     const {return fX;}
  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetdEdx()  const {return fdEdx;}
  Double_t GetPIDsignal()  const {return GetdEdx();}

  Double_t GetY()   const {return fP0;}
  Double_t GetZ()   const {return fP1;}
  Double_t GetSnp() const {return fX*fP4 - fP2;}
  Double_t Get1Pt() const;
  Double_t GetTgl() const {return fP3;}

  Double_t GetSigmaY2() const {return fC00;}
  Double_t GetSigmaZ2() const {return fC11;}

// Some methods useful for users. Implementation is not
// optimized for speed but for minimal maintanance effort
  Double_t Phi() const;
  Double_t Theta() const {return TMath::Pi()/2.-TMath::ATan(GetTgl());}
  Double_t Px() const {return TMath::Cos(Phi())/TMath::Abs(Get1Pt());}
  Double_t Py() const {return TMath::Sin(Phi())/TMath::Abs(Get1Pt());}
  Double_t Pz() const {return GetTgl()/TMath::Abs(Get1Pt());}
  Double_t Pt() const {return 1./TMath::Abs(Get1Pt());}
  Double_t P() const {return TMath::Sqrt(Pt()*Pt()+Pz()*Pz());}

  Int_t Compare(const TObject *o) const;

  void GetExternalParameters(Double_t& xr, Double_t x[5]) const ;
  void GetExternalCovariance(Double_t cov[15]) const ;

  // [SR, 01.04.2003]
  Int_t GetNumber() const {return fNumber;}
  void  SetNumber(Int_t n) {fNumber = n;} 
  //

  Int_t GetClusterIndex(Int_t i) const {return fIndex[i];}

//******** To be removed next release !!! **************
  Double_t GetEta() const {return fP2;}
  Double_t GetC()   const {return fP4;}
  void GetCovariance(Double_t cc[15]) const {
    cc[0 ]=fC00;
    cc[1 ]=fC10;  cc[2 ]=fC11;
    cc[3 ]=fC20;  cc[4 ]=fC21;  cc[5 ]=fC22;
    cc[6 ]=fC40;  cc[7 ]=fC41;  cc[8 ]=fC42;  cc[9 ]=fC44;
    cc[10]=fC30;  cc[11]=fC31;  cc[12]=fC32;  cc[13]=fC43;  cc[14]=fC33;
  }
//****************************************************** 

  virtual Double_t GetPredictedChi2(const AliCluster *cluster) const;
  Int_t PropagateTo(Double_t xr,Double_t x0=28.94,Double_t rho=0.9e-3);
  Int_t Update(const AliCluster* c, Double_t chi2, UInt_t i);
  void ResetCovariance();
  void ResetClusters() {SetNumberOfClusters(0); SetChi2(0.);}
  void UpdatePoints();//update points 
  Float_t* GetPoints() {return fPoints;}
  //
  //  TClonesArray * fHelixIn;    //array of helixes  
  //TClonesArray * fHelixOut;   //array of helixes 
  //
  Float_t Density(Int_t row0, Int_t row1); //calculate cluster density
  Float_t Density2(Int_t row0, Int_t row1); //calculate cluster density
  Double_t GetZat0() const;
  Double_t GetD(Double_t x=0, Double_t y=0) const;
  AliExternalTrackParam & GetReference(){ return fReference;}
  void UpdateReference(){ new (&fReference) AliExternalTrackParam(*this);}
  Int_t   GetKinkIndex(Int_t i) const{ return fKinkIndexes[i];}
  Int_t*  GetKinkIndexes() { return fKinkIndexes;}
  Int_t   GetV0Index(Int_t i) const{ return fV0Indexes[i];}
  Int_t*  GetV0Indexes() { return fV0Indexes;}
protected: 
  void GetXYZ(Float_t r[3]) const;

  Double_t fX;              // X-coordinate of this track (reference plane)
  Double_t fAlpha;          // Rotation angle the local (TPC sector)
                            // coordinate system and the global ALICE one.

  Double_t fdEdx;           // dE/dx

  Double_t fP0;             // Y-coordinate of a track
  Double_t fP1;             // Z-coordinate of a track
  Double_t fP2;             // C*x0
  Double_t fP3;             // tangent of the track momentum dip angle
  Double_t fP4;             // track curvature

  Double_t fC00;                         // covariance
  Double_t fC10, fC11;                   // matrix
  Double_t fC20, fC21, fC22;             // of the
  Double_t fC30, fC31, fC32, fC33;       // track
  Double_t fC40, fC41, fC42, fC43, fC44; // parameters

  Int_t fIndex[kMaxRow];       // indices of associated clusters 
  Float_t fPoints[4];            //first, max dens row  end points of the track and max density
  //[SR, 01.04.2003]
  Int_t fNumber;         // magic number used for number of clusters
  // MI addition
  Float_t fSdEdx;           // sigma of dedx 
  //
  Int_t   fNFoundable;      //number of foundable clusters - dead zone taken to the account
  Bool_t  fBConstrain;   // indicate seeding with vertex constrain
  Int_t   fLastPoint;     // last  cluster position     
  Int_t   fFirstPoint;    // first cluster position
  Int_t fRemoval;         // removal factor
  Int_t fTrackType;       // track type - 0 - normal - 1 - kink -  2 -V0  3- double found
  Int_t fLab2;            // index of corresponding track (kink, V0, double)
  Int_t   fNShared;       // number of shared points 
  AliExternalTrackParam   fReference; // track parameters at the middle of the chamber
  Float_t  fKinkPoint[12];      //radius, of kink,  dfi and dtheta
  Int_t    fKinkIndexes[3];     // kink indexes - minus = mother + daughter
  Int_t    fV0Indexes[3];     // kink indexes - minus = mother + daughter

  ClassDef(AliTPCtrack,2)   // Time Projection Chamber reconstructed tracks
};

inline 
void AliTPCtrack::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  // This function return external TPC track representation
     xr=fX;
     x[0]=GetY(); x[1]=GetZ(); x[2]=GetSnp(); x[3]=GetTgl();
     x[4]=(TMath::Sign(1e-9,fP4) + fP4)*GetLocalConvConst();
}

inline void AliTPCtrack::GetXYZ(Float_t r[3]) const {
  //---------------------------------------------------------------------
  // Returns the position of the track in the global coord. system 
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  r[0]=fX*cs - fP0*sn; r[1]=fX*sn + fP0*cs; r[2]=fP1;
}


#endif

