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

const unsigned kMAXCLUSTERSPERTRACK=210; 

class AliTRDtrack : public AliKalmanTrack {

  //////////////////////////////////////////////////////////////////
  // Represents a reconstructed TRD track                         //
  //////////////////////////////////////////////////////////////////

  friend class AliTRDtracker;

public:

   AliTRDtrack():AliKalmanTrack(){fBackupTrack=0;}
   AliTRDtrack(const AliTRDcluster *c, UInt_t index, const Double_t xx[5],
               const Double_t cc[15], Double_t xr, Double_t alpha);  
   AliTRDtrack(const AliTRDtrack& t);    
   AliTRDtrack(const AliKalmanTrack& t, Double_t alpha); 
   AliTRDtrack(const AliESDtrack& t);    
   ~AliTRDtrack();
   Int_t    Compare(const TObject *o) const;
   void     CookdEdx(Double_t low=0.05, Double_t up=0.55);   
   Float_t  StatusForTOF();
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
   Float_t  GetPIDsignals(Int_t i) const {return fdEdxPlane[i];}
   Int_t    GetPIDTimBin(Int_t i) const {return fTimBinPlane[i];}
   Double_t GetEta()   const {return fE;}

   void     GetExternalCovariance(Double_t cov[15]) const ;   
   void     GetExternalParameters(Double_t& xr, Double_t x[5]) const ;

   Double_t GetLikelihoodElectron() const { return fLhElectron; };

   Double_t Get1Pt()   const {return (1e-9*TMath::Abs(fC)/fC + fC)*GetConvConst(); } 
   Double_t GetP()     const {  
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetPredictedChi2(const AliTRDcluster *c, Double_t h01) const ;
   Double_t GetPt()    const {return 1./Get1Pt();}   
   void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const ;
   void     GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const ;
   Int_t    GetSeedLabel() const { return fSeedLab; }
   Double_t GetSigmaC2()   const {return fCcc;}
   Double_t GetSigmaTgl2() const {return fCtt;}
   Double_t GetSigmaY2()   const {return fCyy;}
   Double_t GetSigmaZ2()   const {return fCzz;}
   Double_t GetSnp()  const {return fX*fC - fE;}
   Double_t GetTgl()  const {return fT;}
   Double_t GetX()    const {return fX;}
   Double_t GetY()    const {return fY;}
   Double_t GetZ()    const {return fZ;}
   UInt_t * GetBackupIndexes()  {return fIndexBackup;}
   UInt_t * GetIndexes()  {return fIndex;}
   //-----------------------------------------------------------------
   // This function calculates the Y-coordinate of a track at the plane x=xk.
   // Needed for matching with the TOF (I.Belikov)
   Double_t GetYat(Double_t xk) const {     
      Double_t c1=fC*fX - fE, r1=TMath::Sqrt(1.- c1*c1);
      Double_t c2=fC*xk - fE, r2=TMath::Sqrt(1.- c2*c2);
      return fY + (xk-fX)*(c1+c2)/(r1+r2);
   }
   //-----------------------------------------------------------------
   Int_t    GetProlongation(Double_t xk, Double_t &y, Double_t &z);

   void     SetStop(Bool_t stop) {fStopped=stop;}
   Bool_t   GetStop() const {return fStopped;}

   Int_t    PropagateTo(Double_t xr, Double_t x0=8.72, Double_t rho=5.86e-3);
   void     ResetCovariance();   
   void     ResetCovariance(Float_t mult);   
   void     ResetClusters() { SetChi2(0.); SetNumberOfClusters(0); }
   Int_t    Rotate(Double_t angle);

   void     SetdEdx(Float_t dedx) {fdEdx=dedx;}  
   void     SetPIDsignals(Float_t dedx, Int_t i) {fdEdxPlane[i]=dedx;}
   void     SetPIDTimBin(Int_t timbin, Int_t i) {fTimBinPlane[i]=timbin;}
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

  void      AddNWrong() {fNWrong++;}
  
  Int_t     GetNWrong() const {return fNWrong;}
  Int_t     GetNRotate() const {return fNRotate;}
  Int_t     GetNCross() const {return fNCross;}
  void      IncCross() {fNCross++; if (fBackupTrack) fBackupTrack->IncCross();}
  AliTRDtrack *  GetBackupTrack(){return fBackupTrack;}
  void      MakeBackupTrack();

protected:

   Int_t    fSeedLab;               // track label taken from seeding  
   Float_t  fdEdx;                  // dE/dx 
   Float_t  fdEdxPlane[kNPlane];    // dE/dx from all 6 planes
   Int_t    fTimBinPlane[kNPlane];  // time bin of Max cluster from all 6 planes

   Double_t fAlpha;         // rotation angle
   Double_t fX;             // running local X-coordinate of the track (time bin)
   Bool_t   fStopped;       // track stop indication

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
   
   UInt_t  fIndex[kMAXCLUSTERSPERTRACK];       // global indexes of clusters  
   UInt_t  fIndexBackup[kMAXCLUSTERSPERTRACK]; // backup indexes of clusters - used in iterations
   Float_t fdQdl[kMAXCLUSTERSPERTRACK];        // cluster amplitudes corrected 
                                               // for track angles    
                           
   Float_t fLhElectron;    // Likelihood to be an electron    
   Int_t   fNWrong;        // number of wrong clusters
   Int_t   fNRotate;       // number of rotation
   Int_t   fNCross;        // number of the cross materials
   Int_t   fNExpected;     // expected number of cluster
   Int_t   fNLast;         // number of clusters in last 2 layers
   Int_t   fNExpectedLast; // number of expected clusters on last 2 layers
   Int_t   fNdedx;         // number of clusters for dEdx measurment
   Float_t fChi2Last;      // chi2 in the  last 2 layers

   AliTRDtrack * fBackupTrack; //! backup track

   ClassDef(AliTRDtrack,2) // TRD reconstructed tracks

};                     


#endif   
