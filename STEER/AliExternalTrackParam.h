#ifndef ALIEXTERNALTRACKPARAM_H
#define ALIEXTERNALTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TObject.h"
#include "AliMagF.h"
#include "TVector3.h"

class AliKalmanTrack;

class AliExternalTrackParam: public TObject {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  AliExternalTrackParam(const AliKalmanTrack& track);

  virtual void SetMass(Double_t mass) {fMass=mass;}
  virtual Double_t GetMass() const {return fMass;}

  virtual const Double_t* GetParameter() const;
  virtual const Double_t* GetCovariance() const;
  virtual Double_t     X() const {return fX;};
  virtual Double_t     Alpha() const {return fAlpha;};
  virtual AliExternalTrackParam* CreateExternalParam() const;
  virtual void         ResetCovariance(Double_t factor = 10.,
				       Bool_t clearOffDiagonal = kTRUE);
  virtual Double_t     Y() const {return fParam[0];};
  virtual Double_t     Z() const {return fParam[1];};
  virtual void         GetXYZ(Float_t r[3]) const;
  virtual void         GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z ) const;
  virtual Bool_t       PropagateTo(Double_t x,  Double_t x0, Double_t rho);
  virtual Bool_t       PropagateToDCA(Double_t x, Double_t y,  Double_t x0, Double_t rho);
  virtual Bool_t       RotateTo(Double_t alpha);
  virtual Bool_t       CorrectForMaterial(Double_t d, Double_t x0, Double_t rho);
  virtual Bool_t       GetProlongationAt(Double_t x, Double_t& y, 
					 Double_t& z) const;
  virtual Double_t     GetXAtVertex(Double_t x = 0, Double_t y = 0) const;

  //  virtual Double_t     GetPredictedChi2(const AliCluster* cluster);
  //  virtual Bool_t       Update(const AliCluster* cluster);

  virtual Double_t     SigmaPhi() const;
  virtual Double_t     SigmaTheta() const;
  virtual Double_t     SigmaPt() const;
  virtual TVector3     Momentum() const;
  virtual TVector3     Position() const;

  virtual void         Print(Option_t* option = "") const;
  // local magnetic field manipulation 
  void     SaveLocalConvConst();
  Double_t GetLocalConvConst() const;
  
  static void SetFieldMap(const AliMagF *map) { fgkFieldMap=map; }
  static const AliMagF *GetFieldMap() { return fgkFieldMap; }

  static void SetUniformFieldTracking() {
     if (fgkFieldMap==0) {
        printf("AliKalmanTrack: Field map has not been set !\n"); 
        exit(1);
     } 
     fgConvConst=1000/0.299792458/(fgkFieldMap->SolenoidField()+1e-13);
  }
  static void SetNonuniformFieldTracking() { fgConvConst=0.; }

 private:
  Double_t             fMass;       // mass associated to the particle
  Double_t             fX;          // x coordinate for the parametrisation
  Double_t             fAlpha;      // azimuthal angle for the parametrisation
  Double_t             fParam[5];   // track parameter (y, z, sin(azimuthal angel), tan(dip angle), 1/pt)
  Double_t             fCovar[15];  // track parameter covariance
  //
  static const AliMagF *fgkFieldMap;//pointer to the magnetic field map
  static Double_t fgConvConst;      //conversion "curvature(1/cm) -> pt(GeV/c)"
  Double_t fLocalConvConst;         //local conversion "curvature(1/cm) -> pt(GeV/c)"


  ClassDef(AliExternalTrackParam, 3)
};


inline void AliExternalTrackParam::SaveLocalConvConst() {
  //---------------------------------------------------------------------
  // Saves local conversion constant "curvature (1/cm) -> pt (GeV/c)" 
  //---------------------------------------------------------------------
     if (fgConvConst > 0 || fgConvConst < 0) return; //uniform field tracking
     Float_t r[3]={0.,0.,0.}; GetXYZ(r);
     Float_t b[3]; fgkFieldMap->Field(r,b);
     fLocalConvConst=1000/0.299792458/(1e-13 - b[2]);
} 

inline Double_t AliExternalTrackParam::GetLocalConvConst() const {
  //---------------------------------------------------------------------
  // Returns conversion constant "curvature (1/cm) -> pt (GeV/c)" 
  //---------------------------------------------------------------------
     if (fgConvConst > 0 || fgConvConst < 0) return fgConvConst; //uniform field tracking
     return fLocalConvConst;
} 


inline void AliExternalTrackParam::GetXYZ(Float_t r[3]) const {
  //---------------------------------------------------------------------
  // Returns the position of the track in the global coord. system 
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  r[0]=fX*cs - fParam[0]*sn; r[1]=fX*sn + fParam[0]*cs; r[2]=fParam[1];
}

inline void AliExternalTrackParam::GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const {
  //---------------------------------------------------------------------
  // Returns the position of the track in the global coord. system 
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  x=fX*cs - fParam[0]*sn; y=fX*sn + fParam[0]*cs; z=fParam[1];
}




#endif
