#ifndef ALIEXTERNALTRACKPARAM_H
#define ALIEXTERNALTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliTrackParam.h"
class AliKalmanTrack;

class AliExternalTrackParam: public AliTrackParam {
 public:
  AliExternalTrackParam();
  AliExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  AliExternalTrackParam(const AliKalmanTrack& track);

  virtual const Double_t* GetParameter() const;
  virtual const Double_t* GetCovariance() const;
  virtual Double_t     X() const {return fX;};
  virtual Double_t     Alpha() const {return fAlpha;};
  virtual AliExternalTrackParam* CreateExternalParam() const;
  virtual void         ResetCovariance(Double_t factor = 10.,
				       Bool_t clearOffDiagonal = kTRUE);
  virtual Double_t     Y() const {return fParam[0];};
  virtual Double_t     Z() const {return fParam[1];};

  virtual Bool_t       PropagateTo(Double_t x, Double_t* length = NULL);
  virtual Bool_t       RotateTo(Double_t alpha);
  virtual Bool_t       CorrectForMaterial(Double_t dAngle2, Double_t dPrel);
  virtual Bool_t       GetProlongationAt(Double_t x, Double_t& y, 
					 Double_t& z) const;
  virtual Double_t     GetXAtVertex(Double_t x = 0, Double_t y = 0) const;

  virtual Double_t     GetPredictedChi2(const AliCluster* cluster);
  virtual Bool_t       Update(const AliCluster* cluster);

  virtual Double_t     SigmaPhi() const;
  virtual Double_t     SigmaTheta() const;
  virtual Double_t     SigmaPt() const;
  virtual TVector3     Momentum() const;
  virtual TVector3     Position() const;

  virtual void         Print(Option_t* option = "") const;

 private:
  Double_t             fX;          // x coordinate for the parametrisation
  Double_t             fAlpha;      // azimuthal angle for the parametrisation
  Double_t             fParam[5];   // track parameter (y, z, sin(azimuthal angel), tan(dip angle), 1/pt)
  Double_t             fCovar[15];  // track parameter covariance

  ClassDef(AliExternalTrackParam, 1)
};

#endif
