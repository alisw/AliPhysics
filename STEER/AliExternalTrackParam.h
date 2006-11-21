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
#include "TObject.h"
#include "TMath.h"

const Double_t kAlmost1=0.999;
const Double_t kAlmost0=1e-33;
const Double_t kVeryBig=1./kAlmost0;

const Double_t kB2C=0.299792458e-3;
const Double_t kAlmost0Field=1.e-13;
const Double_t kMostProbablePt=0.35;

class AliESDVertex;

Double_t ApproximateBetheBloch(Double_t);

class AliExternalTrackParam: public TObject {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(const AliExternalTrackParam &);
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  virtual ~AliExternalTrackParam(){}

  void Set(Double_t x,Double_t alpha,
			const Double_t param[5], const Double_t covar[15]);
  void Reset();
  void ResetCovariance(Double_t s2) {
    fC[0]*= s2;
    fC[1] = 0.;  fC[2]*= s2;
    fC[3] = 0.;  fC[4] = 0.;  fC[5]*= s2;
    fC[6] = 0.;  fC[7] = 0.;  fC[8] = 0.;  fC[9]*= s2;
    fC[10]= 0.;  fC[11]= 0.;  fC[12]= 0.;  fC[13]= 0.;  fC[14]*=s2;
  }

  const Double_t *GetParameter() const {return fP;}
  const Double_t *GetCovariance() const {return fC;}

  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetX() const {return fX;}
  Double_t GetY()    const {return fP[0];}
  Double_t GetZ()    const {return fP[1];}
  Double_t GetSnp()  const {return fP[2];}
  Double_t GetTgl()  const {return fP[3];}
  Double_t Get1Pt()  const {return fP[4];}

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

  Double_t GetSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t GetP() const;
  Double_t GetPt() const {
    return (TMath::Abs(fP[4])>kAlmost0) ? 1./fP[4]:TMath::Sign(kVeryBig,fP[4]);
  }
  Double_t Get1P() const;
  Double_t GetC(Double_t b) const {return fP[4]*b*kB2C;}
  void GetDZ(Double_t x,Double_t y,Double_t z,Double_t b,Float_t dz[2]) const; 
  Double_t GetD(Double_t xv, Double_t yv, Double_t b) const; 
  Double_t GetLinearD(Double_t xv, Double_t yv) const; 
  Bool_t CorrectForMaterial(Double_t d, Double_t x0, Double_t mass,
			    Double_t (*f)(Double_t)=ApproximateBetheBloch);
  Double_t GetPredictedChi2(Double_t p[2],Double_t cov[3]) const;
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
  Bool_t PropagateToDCA(const AliESDVertex *vtx, Double_t b, Double_t maxd);

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

protected:
  Double_t &Par(Int_t i) {return fP[i];}
  Double_t &Cov(Int_t i) {return fC[i];}

private:
  Double32_t           fX;     // X coordinate for the point of parametrisation
  Double32_t           fAlpha; // Local <-->global coor.system rotation angle
  Double32_t           fP[5];  // The track parameters
  Double32_t           fC[15]; // The track parameter covariance matrix

  ClassDef(AliExternalTrackParam, 5)
};

#endif
