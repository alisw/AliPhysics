#ifndef ALITRDTRACK_H
#define ALITRDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliKalmanTrack.h>
#include <TMath.h>

#include "AliTRDgeometry.h"
#include "AliESDtrack.h"
#include "TVector2.h"

class AliTRDcluster;
class AliTPCtrack;
class AliESDtrack;
class AliTrackReference;

const unsigned kMAX_CLUSTERS_PER_TRACK=210; 

class AliTRDtracklet :public TObject{
  friend class AliTRDtrack;
 public:
  AliTRDtracklet();
  void Set(Float_t x, Float_t y, Float_t z, Float_t alpha, Float_t error2){fX=x; fY=y; fZ=z; fAlpha=alpha; fSigma2= error2;}
  void SetP0(Float_t p0){fP0=p0;}
  void SetP1(Float_t p1){fP1=p1;}
  void SetN(Int_t n){fNFound=n;}
  void SetNCross(Int_t nc){fNCross=nc;}
  void SetPlane(Int_t plane){fPlane=plane;}
  void SetSigma2(Float_t sigma2){fExpectedSigma2=sigma2;}
  void SetChi2(Float_t chi2){fChi2=chi2;}
  void SetTilt(Float_t tilt){fTilt=tilt;}
  void SetMaxPos(Short_t pos, Short_t pos4, Short_t pos5){fMaxPos = pos; fMaxPos4 = pos4; fMaxPos5 = pos5;}
  Float_t GetX() const { return fX;}
  Float_t GetY() const { return fY;}
  Float_t GetZ() const {return fZ;}
  Float_t GetAlpha() const { return fAlpha;}
  Float_t GetTrackletSigma2() const { return fSigma2;}
  //
  Float_t GetP0() const {return fP0;}
  Float_t GetP1() const {return fP1;}
  Int_t GetN() const {return fNFound;}
  Int_t GetNCross() const {return fNCross;}  
  Int_t GetPlane() const {return fPlane;}
  Float_t GetClusterSigma2() const {return fExpectedSigma2;}
  Float_t GetChi2() const {return fChi2;}
  Float_t GetTilt() const {return fTilt;}
 protected:
  Float_t fY;                 // y position
  Float_t fZ;                 // z position
  Float_t fX;                 // x position
  Float_t fAlpha;             // rotation angle
  Float_t fSigma2;            // expected error of tracklet position
  Float_t fP0;                // offset in y
  Float_t fP1;                // offset in tangent
  Int_t   fNFound;            // number of found clusters
  Int_t   fNCross;            // number of crosses
  Int_t   fPlane;             // plane number
  Float_t fExpectedSigma2;    // expected sigma of residual distribution of clusters
  Float_t fChi2;              // chi2 of the tracklet
  Float_t fTilt;              // tilt factor 
  Short_t fMaxPos;            // time bin with max charge
  Short_t fMaxPos4;            // time bin with max charge
  Short_t fMaxPos5;            // time bin with max charge
  ClassDef(AliTRDtracklet,2)
};


class AliTRDtrack : public AliKalmanTrack {

// Represents reconstructed TRD track
  friend class AliTRDtracker;
public:

   AliTRDtrack():AliKalmanTrack(){fBackupTrack=0;}
   AliTRDtrack(const AliTRDcluster *c, UInt_t index, const Double_t xx[5],
               const Double_t cc[15], Double_t xr, Double_t alpha);  
   AliTRDtrack(const AliTRDtrack& t);    
   AliTRDtrack(const AliKalmanTrack& t, Double_t alpha); 
   AliTRDtrack(const AliESDtrack& t);    
   static AliTRDtrack * MakeTrack(const AliTrackReference *ref, Double_t mass);
   ~AliTRDtrack();
   Int_t    Compare(const TObject *o) const;
   void     CookdEdx(Double_t low=0.05, Double_t up=0.55);   
   Float_t    StatusForTOF();
   Double_t GetAlpha() const {return fAlpha;}
   Int_t    GetSector() const {
     //if (fabs(fAlpha) < AliTRDgeometry::GetAlpha()/2) return 0;
     return Int_t(TVector2::Phi_0_2pi(fAlpha)/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;}

   Double_t GetC()     const {return fC;}
   Int_t    GetClusterIndex(Int_t i) const {return fIndex[i];}    
   Float_t  GetClusterdQdl(Int_t i) const {return fdQdl[i];}    

   void     GetCovariance(Double_t cov[15]) const;  
   Double_t GetdEdx()  const {return fdEdx;}
   Double_t GetPIDsignal()  const {return GetdEdx();}
   Float_t GetPIDsignals(Int_t i) const {return fdEdxPlane[i];}
   Int_t  GetPIDTimBin(Int_t i) const {return fTimBinPlane[i];}
   Double_t GetEta()   const {return fE;}

   void     GetExternalCovariance(Double_t cov[15]) const ;   
   void     GetExternalParameters(Double_t& xr, Double_t x[5]) const ;

   Double_t GetLikelihoodElectron() const { return fLhElectron; };

   Double_t Get1Pt() const {
      return (TMath::Sign(1e-9,fC) + fC)*GetLocalConvConst();
   }
   Double_t GetP()     const {  
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetPredictedChi2(const AliTRDcluster*, Double_t h01) const ;
   Double_t GetPt()    const {return 1./Get1Pt();}   
   void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const ;
   void     GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const ;
   Int_t    GetSeedLabel() const { return fSeedLab; }
   Double_t GetSigmaC2()   const {return fCcc;}
   Double_t GetSigmaTgl2() const {return fCtt;}
   Double_t GetSigmaY2()   const {return TMath::Abs(fCyy);}
   Double_t GetSigmaZ2()   const {return TMath::Abs(fCzz);}
   Double_t GetSnp()  const {return fX*fC - fE;}
   Double_t GetTgl()  const {return fT;}
   Double_t GetX()    const {return fX;}
   Double_t GetY()    const {return fY;}
   Double_t GetZ()    const {return fZ;}
   UInt_t * GetBackupIndexes()  {return fIndexBackup;}
   UInt_t * GetIndexes()  {return fIndex;}
   Double_t GetYat(Double_t xk) const {     
//-----------------------------------------------------------------
// This function calculates the Y-coordinate of a track at the plane x=xk.
// Needed for matching with the TOF (I.Belikov)
//-----------------------------------------------------------------
      Double_t c1=fC*fX - fE, r1=TMath::Sqrt(1.- c1*c1);
      Double_t c2=fC*xk - fE, r2=TMath::Sqrt(1.- c2*c2);
      return fY + (xk-fX)*(c1+c2)/(r1+r2);
   }
   Int_t GetProlongation(Double_t xk, Double_t &y, Double_t &z);

   void SetStop(Bool_t stop) {fStopped=stop;}
   Bool_t GetStop() const {return fStopped;}

   Int_t    PropagateTo(Double_t xr, Double_t x0=8.72, Double_t rho=5.86e-3);
   Int_t    PropagateToX(Double_t xr, Double_t step);
   Int_t    PropagateToR(Double_t xr, Double_t step);
   void     ResetCovariance();   
   void     ResetCovariance(Float_t mult);   
   void ResetClusters() { SetChi2(0.); SetNumberOfClusters(0); }
   Int_t    Rotate(Double_t angle, Bool_t absolute=kFALSE);

   void     SetdEdx(Float_t dedx) {fdEdx=dedx;}  
   void SetPIDsignals(Float_t dedx, Int_t i) {fdEdxPlane[i]=dedx;}
   void  SetPIDTimBin(Int_t timbin, Int_t i) {fTimBinPlane[i]=timbin;}
   void     SetLikelihoodElectron(Float_t l) { fLhElectron = l; };  

   void     SetSampledEdx(Float_t q, Int_t i) {
               Double_t s=GetSnp(), t=GetTgl();
               q*= TMath::Sqrt((1-s*s)/(1+t*t));
               fdQdl[i]=q;
            }     
   void     SetSampledEdx(Float_t q) {
              Double_t s=GetSnp(), t=GetTgl();
              q*= TMath::Sqrt((1-s*s)/(1+t*t));
              fdQdl[fNdedx]=q;
	      fNdedx++;
            }     

   void     SetSeedLabel(Int_t lab) { fSeedLab=lab; }

   Int_t    Update(const AliTRDcluster* c, Double_t chi2, UInt_t i, 
                   Double_t h01);
   Int_t    UpdateMI(const AliTRDcluster* c, Double_t chi2, UInt_t i, 
                   Double_t h01, Int_t plane);
   Int_t    UpdateMI(const AliTRDtracklet & tracklet);

  //
  void AddNWrong() {fNWrong++;}
  
  Int_t GetNWrong() const {return fNWrong;}
  Int_t GetNRotate() const {return fNRotate;}
  Int_t GetNCross() const {return fNCross;}
  void  IncCross() {fNCross++; if (fBackupTrack) fBackupTrack->IncCross();}
  AliTRDtrack *  GetBackupTrack(){return fBackupTrack;}
  void    MakeBackupTrack();
  //


protected:
   void GetXYZ(Float_t r[3]) const;

   Double_t GetPredictedChi2(const AliCluster*/*c*/) const {return 0.;}
   Int_t Update(const AliCluster*/*c*/, Double_t /*chi2*/, UInt_t /*i*/) {
     return 0;
   }

   Int_t    fSeedLab;     // track label taken from seeding  
   Float_t  fdEdx;        // dE/dx 
   Float_t  fdEdxT;        // dE/dx  - truncated mean
   Float_t  fDE;          // integrated delta energy
   Float_t  fdEdxPlane[kNPlane];  // dE/dx from all 6 planes
   Int_t  fTimBinPlane[kNPlane];  // time bin of Max cluster from all 6 planes

   Double_t fAlpha;       // rotation angle
   Double_t fX;           // running local X-coordinate of the track (time bin)
   Bool_t   fStopped;     // track stop indication

   Double_t fY;             // Y-coordinate of the track
   Double_t fZ;             // Z-coordinate of the track
   Double_t fE;             // C*x0
   Double_t fT;             // tangent of the track momentum dip angle
   Double_t fC;             // track curvature

   Double_t fCyy;                         // covariance
   Double_t fCzy, fCzz;                   // matrix
   Double_t fCey, fCez, fCee;             // of the
   Double_t fCty, fCtz, fCte, fCtt;       // track
   Double_t fCcy, fCcz, fCce, fCct, fCcc; // parameters   
   
   UInt_t  fIndex[kMAX_CLUSTERS_PER_TRACK];  // global indexes of clusters  
   UInt_t  fIndexBackup[kMAX_CLUSTERS_PER_TRACK]; //backup indexes of clusters - used in iterations
   Float_t fdQdl[kMAX_CLUSTERS_PER_TRACK];   // cluster amplitudes corrected 
                                             // for track angles    
                           
   Float_t fLhElectron;    // Likelihood to be an electron    
   Int_t fNWrong;    // number of wrong clusters
   Int_t fNRotate;   // number of rotation
   Int_t fNCross;     // number of the cross materials
   Int_t fNExpected;  //expected number of cluster
   Int_t fNLast;      //number of clusters in last 2 layers
   Int_t fNExpectedLast; //number of expected clusters on last 2 layers
   Int_t      fNdedx;      //number of clusters for dEdx measurment
   Float_t fChi2Last;      //chi2 in the  last 2 layers
   AliTRDtracklet fTracklets[6]; //tracklets
   Float_t     fBudget[3];       // integrated material budget
   AliTRDtrack * fBackupTrack; //! backup track
   ClassDef(AliTRDtrack,3) // TRD reconstructed tracks
};                     

inline void AliTRDtrack::GetXYZ(Float_t r[3]) const {
  //---------------------------------------------------------------------
  // Returns the position of the track in the global coord. system 
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  r[0]=fX*cs - fY*sn; r[1]=fX*sn + fY*cs; r[2]=fZ;
}

#endif   
