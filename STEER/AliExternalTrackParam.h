#ifndef ALIEXTERNALTRACKPARAM_H
#define ALIEXTERNALTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/*****************************************************************************
 *              "External" track parametrisation class                       *
 *                                                                           *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *                                                                           *
 * The parameters are estimated at an exact position x in a local coord.     *
 * system rotated by angle alpha with respect to the global coord.system.    *
 *        Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                     *
 *****************************************************************************/
#include "TMath.h"

#include "AliVTrack.h"

const Double_t kVeryBig=1./kAlmost0;
const Double_t kMostProbablePt=0.35;

Double_t ApproximateBetheBloch(Double_t);

class AliVVertex;
class TPolyMarker3D; 

class AliExternalTrackParam: public AliVTrack {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(const AliExternalTrackParam &);
  AliExternalTrackParam& operator=(const AliExternalTrackParam & trkPar);
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  AliExternalTrackParam(const AliVTrack *vTrack);
  AliExternalTrackParam(Double_t xyz[3],Double_t pxpypz[3],
			Double_t cv[21],Short_t sign);
  virtual ~AliExternalTrackParam(){}

  void Set(Double_t x,Double_t alpha,
			const Double_t param[5], const Double_t covar[15]);
  void Set(Double_t xyz[3],Double_t pxpypz[3],Double_t cv[21],Short_t sign);

  static void SetMostProbablePt(Double_t pt) { fgMostProbablePt=pt; }
  static Double_t GetMostProbablePt() { return fgMostProbablePt; }

  void Reset();
  void ResetCovariance(Double_t s2);
  void AddCovariance(const Double_t cov[15]);

  const Double_t *GetParameter() const {return fP;}
  const Double_t *GetCovariance() const {return fC;}

  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetX() const {return fX;}
  Double_t GetY()    const {return fP[0];}
  Double_t GetZ()    const {return fP[1];}
  Double_t GetSnp()  const {return fP[2];}
  Double_t GetTgl()  const {return fP[3];}
  Double_t GetSigned1Pt()  const {return fP[4];}

  Double_t GetSigmaY2() const {return fC[0];}
  Double_t GetSigmaZY() const {return fC[1];}
  Double_t GetSigmaZ2() const {return fC[2];}
  Double_t GetSigmaSnpY() const {return fC[3];}
  Double_t GetSigmaSnpZ() const {return fC[4];}
  Double_t GetSigmaSnp2() const {return fC[5];}
  Double_t GetSigmaTglY() const {return fC[6];}
  Double_t GetSigmaTglZ() const {return fC[7];}
  Double_t GetSigmaTglSnp() const {return fC[8];}
  Double_t GetSigmaTgl2() const {return fC[9];}
  Double_t GetSigma1PtY() const {return fC[10];}
  Double_t GetSigma1PtZ() const {return fC[11];}
  Double_t GetSigma1PtSnp() const {return fC[12];}
  Double_t GetSigma1PtTgl() const {return fC[13];}
  Double_t GetSigma1Pt2() const {return fC[14];}

  // additional functions for AliVParticle
  Double_t Px() const;
  Double_t Py() const;
  Double_t Pz() const;
  Double_t Pt() const { return TMath::Abs(GetSignedPt()); }
  Double_t P() const { return GetP(); }
  Bool_t   PxPyPz(Double_t p[3]) const { return GetPxPyPz(p); }
  
  Double_t Xv() const;
  Double_t Yv() const;
  Double_t Zv() const;
  Bool_t   XvYvZv(Double_t x[3]) const { return GetXYZ(x); }

  Double_t OneOverPt() const { return 1./Pt(); }
  Double_t Phi() const;
  Double_t Theta() const;
  virtual Double_t E() const;
  virtual Double_t M() const;
  Double_t Eta() const;
  virtual Double_t Y() const;
  Short_t  Charge() const { return (Short_t)GetSign(); }
  virtual const Double_t *PID() const { return 0x0; }

  // additional functions from AliVTrack
  virtual Int_t    GetID() const { return -999; }
  virtual UChar_t  GetITSClusterMap() const {return 0; }
  virtual ULong_t  GetStatus() const { return 0; }

  Double_t GetSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t GetP() const;
  Double_t GetSignedPt() const {
    return (TMath::Abs(fP[4])>kAlmost0) ? 1./fP[4]:TMath::Sign(kVeryBig,fP[4]);
  }
  Double_t Get1P() const;
  Double_t GetC(Double_t b) const {return fP[4]*b*kB2C;}
  void GetDZ(Double_t x,Double_t y,Double_t z,Double_t b,Float_t dz[2]) const; 
  Double_t GetD(Double_t xv, Double_t yv, Double_t b) const; 
  Double_t GetLinearD(Double_t xv, Double_t yv) const; 
  Bool_t CorrectForMeanMaterial(Double_t xOverX0, Double_t xTimesRho, 
        Double_t mass,  Bool_t anglecorr=kFALSE,
	Double_t (*f)(Double_t)=ApproximateBetheBloch);
  Double_t GetPredictedChi2(Double_t p[2],Double_t cov[3]) const;

  Double_t 
    GetPredictedChi2(Double_t p[3],Double_t covyz[3],Double_t covxyz[3]) const;
  Bool_t 
    PropagateTo(Double_t p[3],Double_t covyz[3],Double_t covxyz[3],Double_t b);

  Double_t *GetResiduals(Double_t *p,Double_t *cov,Bool_t updated=kTRUE) const;
  Bool_t Update(Double_t p[2],Double_t cov[3]);
  Bool_t Rotate(Double_t alpha);
  Bool_t PropagateTo(Double_t x, Double_t b);
  Bool_t Propagate(Double_t alpha, Double_t x, Double_t b) {
    if (Rotate(alpha))
      if (PropagateTo(x,b)) return kTRUE;
    return kFALSE;
  }
  void   Propagate(Double_t len,Double_t x[3],Double_t p[3],Double_t bz) const;
  Bool_t Intersect(Double_t pnt[3], Double_t norm[3], Double_t bz) const;

  void GetHelixParameters(Double_t h[6], Double_t b) const;
  Double_t GetDCA(const AliExternalTrackParam *p, Double_t b,
    Double_t &xthis,Double_t &xp) const;
  Double_t PropagateToDCA(AliExternalTrackParam *p, Double_t b);
  Bool_t PropagateToDCA(const AliVVertex *vtx, Double_t b, Double_t maxd,
                        Double_t dz[2]=0, Double_t cov[3]=0);
  
  void GetDirection(Double_t d[3]) const;
  Bool_t GetPxPyPz(Double_t *p) const;
  Bool_t GetXYZ(Double_t *p) const;
  Bool_t GetCovarianceXYZPxPyPz(Double_t cv[21]) const;
  Bool_t GetPxPyPzAt(Double_t x, Double_t b, Double_t p[3]) const;
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t r[3]) const;
  Bool_t GetYAt(Double_t x,  Double_t b,  Double_t &y) const;
  Bool_t GetZAt(Double_t x,  Double_t b,  Double_t &z) const;
  void Print(Option_t* option = "") const;
  Double_t GetSnpAt(Double_t x,Double_t b) const;

  //Deprecated
  Bool_t CorrectForMaterial(Double_t d, Double_t x0, Double_t mass,
			    Double_t (*f)(Double_t)=ApproximateBetheBloch);

  Bool_t GetDistance(AliExternalTrackParam *param2, Double_t x, Double_t dist[3], Double_t b);
  Int_t GetIndex(Int_t i, Int_t j) const;
  Int_t GetLabel() const {return -1;} 
  //
  // visualization (M. Ivanov)
  //
  virtual void FillPolymarker(TPolyMarker3D *pol, Float_t magf, Float_t minR, Float_t maxR, Float_t stepR);
  virtual void DrawTrack(Float_t magF, Float_t minR, Float_t maxR, Float_t stepR);
 protected:
  Double_t &Par(Int_t i) {return fP[i];}
  Double_t &Cov(Int_t i) {return fC[i];}
  
 private:
  Double32_t           fX;     // X coordinate for the point of parametrisation
  Double32_t           fAlpha; // Local <-->global coor.system rotation angle
  Double32_t           fP[5];  // The track parameters
  Double32_t           fC[15]; // The track parameter covariance matrix

  static Double32_t    fgMostProbablePt; // "Most probable" pt
                                         // (to be used if Bz=0)
  ClassDef(AliExternalTrackParam, 8)
};

inline void AliExternalTrackParam::ResetCovariance(Double_t s2) {
  //
  // Reset the covarince matrix to "something big"
  //
    fC[0]*= s2;
    fC[1] = 0.;  fC[2]*= s2;
    fC[3] = 0.;  fC[4] = 0.;  fC[5]*= s2;
    fC[6] = 0.;  fC[7] = 0.;  fC[8] = 0.;  fC[9]*= s2;
    fC[10]= 0.;  fC[11]= 0.;  fC[12]= 0.;  fC[13]= 0.;  fC[14]*=s2;
}




#endif
