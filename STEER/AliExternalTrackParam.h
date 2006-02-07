#ifndef ALIEXTERNALTRACKPARAM_H
#define ALIEXTERNALTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TObject.h"

class AliKalmanTrack;

class AliExternalTrackParam: public TObject {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  AliExternalTrackParam(const AliKalmanTrack& track);

  void Reset();
  void Set(const AliKalmanTrack& track);

  const Double_t* GetParameter() const {return fP;}
  const Double_t* GetCovariance() const {return fC;}
  virtual Double_t GetX() const {return fX;}
  virtual Double_t GetAlpha() const {return fAlpha;}
  Double_t GetSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t GetP() const;
  Double_t GetD(Double_t b, Double_t x=0, Double_t y=0) const; 
  Bool_t GetPxPyPz(Double_t *p) const;
  Bool_t GetXYZ(Double_t *p) const;
  Bool_t GetCovarianceXYZPxPyPz(Double_t cv[21]) const;
  Bool_t GetPxPyPzAt(Double_t x, Double_t b, Double_t p[3]) const;
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t r[3]) const;

  void Print(Option_t* option = "") const;

private:
  Double_t             fX;      // x coordinate for the parametrisation
  Double_t             fAlpha;  // azimuthal angle for the parametrisation
  Double_t             fP[5];   // track parameter (y, z, sin(azimuthal angel), tan(dip angle), 1/pt)
  Double_t             fC[15];  // track parameter covariance

  ClassDef(AliExternalTrackParam, 4)
};

#endif
