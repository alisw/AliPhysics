#ifndef ALITOFTRACK_H
#define ALITOFTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliKalmanTrack.h>
#include <TMath.h>
#include "AliTOFGeometry.h"
#include "TVector2.h"

class AliESDtrack;
class AliTOFtrack : public AliKalmanTrack {

// Represents reconstructed TOF track

public:

   AliTOFtrack():AliKalmanTrack(){}
   AliTOFtrack(const AliTOFtrack& t);    
   AliTOFtrack(const AliESDtrack& t);    


   Double_t GetAlpha() const {return fAlpha;}
   Int_t    GetSector() const {
     return Int_t(TVector2::Phi_0_2pi(fAlpha)/AliTOFGeometry::GetAlpha())%AliTOFGeometry::NSectors();}

   Double_t GetC()     const {return fC;}
   void     GetCovariance(Double_t cov[15]) const;  
   Double_t GetEta()   const {return fE;}

   void     GetExternalCovariance(Double_t cov[15]) const ;   
   void     GetExternalParameters(Double_t& xr, Double_t x[5]) const ;
   Double_t GetSigmaY2() const {return fCyy;}
   Double_t GetSigmaZ2() const {return fCzz;}

   Double_t Get1Pt()   const {return (1e-9*TMath::Abs(fC)/fC + fC)*GetConvConst();}
   Double_t GetP()     const {  
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetPt()    const {return 1./Get1Pt();}   
   void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const ;
   void     GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const ;
   Int_t    GetSeedLabel() const { return fSeedLab; }
   Int_t    GetSeedIndex() const { return fSeedInd; }
   void     SetSeedLabel(Int_t lab) { fSeedLab=lab; }
   void     SetSeedIndex(Int_t ind) { fSeedInd=ind; }
   Double_t GetSnp()  const {return fX*fC - fE;}
   Double_t GetTgl()  const {return fT;}
   Double_t GetX()    const {return fX;}
   Double_t GetY()    const {return fY;}
   Double_t GetZ()    const {return fZ;}
  
   Int_t Compare(const TObject *o) const;


   Double_t GetYat(Double_t xk, Bool_t skip) const {     
//-----------------------------------------------------------------
// This function calculates the Y-coordinate of a track at the plane x=xk.
// Needed for matching with the TOF (I.Belikov)
//-----------------------------------------------------------------
     skip=kFALSE;
     Double_t c1=fC*fX - fE, r1=TMath::Sqrt(TMath::Abs(1.- c1*c1));
     Double_t c2=fC*xk - fE, r2=TMath::Sqrt(TMath::Abs(1.- c2*c2));
      if( ((1.- c2*c2)<0) || ((1.- c1*c1)<0) ) skip=kTRUE;
      return fY + (xk-fX)*(c1+c2)/(r1+r2);
   }


   Int_t    PropagateTo(Double_t xr, Double_t x0=8.72, Double_t rho=5.86e-3);
   Int_t    PropagateToInnerTOF(Bool_t holes);
   void     ResetCovariance();   
   void     ResetCovariance(Float_t mult);   
   Int_t    Rotate(Double_t angle);



protected:

  
   Int_t    fSeedInd;     // ESD seed track index  
   Int_t    fSeedLab;     // track label taken from seeding  
   Double_t fAlpha;       // rotation angle
   Double_t fX;           // running local X-coordinate of the track (time bin)
  
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

 private:

   void GetPropagationParameters(Bool_t holes, Double_t *param);
   
   ClassDef(AliTOFtrack,0) // TOF reconstructed tracks

};                     

#endif   
