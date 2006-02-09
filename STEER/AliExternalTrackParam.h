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

const Double_t kB2C=0.299792458e-3;
const Double_t kAlmost1=0.999;
const Double_t kAlmost0=1e-33;
const Double_t kVeryBig=1./kAlmost0;

class AliKalmanTrack;

class AliExternalTrackParam: public TObject {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  AliExternalTrackParam(const AliKalmanTrack& track);

  void Reset();
  void Set(const AliKalmanTrack& track);

  const Double_t *GetParameter() const {return fP;}
  const Double_t *GetCovariance() const {return fC;}
  virtual Double_t GetX() const {return fX;}
  virtual Double_t GetAlpha() const {return fAlpha;}
  Double_t GetSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t GetP() const;
  Double_t GetD(Double_t b, Double_t xv, Double_t yv) const; 
  Double_t GetLinearD(Double_t xv, Double_t yv) const; 
  Double_t GetPredictedChi2(Double_t p[2],Double_t cov[3]) const;
  Bool_t Update(Double_t p[2],Double_t cov[3]);
  Bool_t Rotate(Double_t alpha);
  Bool_t PropagateTo(Double_t x, Double_t b);
  Bool_t Propagate(Double_t alpha, Double_t x, Double_t b) {
    if (Rotate(alpha))
      if (PropagateTo(x,b)) return kTRUE;
    return kFALSE;
  }

  Bool_t GetPxPyPz(Double_t *p) const;
  Bool_t GetXYZ(Double_t *p) const;
  Bool_t GetCovarianceXYZPxPyPz(Double_t cv[21]) const;
  Bool_t GetPxPyPzAt(Double_t x, Double_t b, Double_t p[3]) const;
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t r[3]) const;

  void Print(Option_t* option = "") const;

private:
  Double_t             fX;     // X coordinate for the point of parametrisation
  Double_t             fAlpha; // Local <-->global coor.system rotation angle
  Double_t             fP[5];  // The track parameters
  Double_t             fC[15]; // The track parameter covariance matrix

  ClassDef(AliExternalTrackParam, 4)
};

#endif
