#ifndef ALITRACKPARAM_H
#define ALITRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TObject.h>
#include <TVector3.h>

class AliCluster;
class AliExternalTrackParam;


class FastMath {
public:
  FastMath();  
  static Double_t FastAsin(Double_t x);   
 private: 
  static Double_t fgFastAsin[20000];
};


class AliTrackParam: public TObject {
 public:
  virtual const Double_t* GetParameter() const = 0;
  virtual const Double_t* GetCovariance() const = 0;
  virtual Double_t     X() const = 0;
  virtual Double_t     Alpha() const = 0;
  virtual AliExternalTrackParam* CreateExternalParam() const = 0;
  virtual void         ResetCovariance(Double_t factor = 10.,
				       Bool_t clearOffDiagonal = kTRUE) = 0;

  virtual Double_t     Y() const = 0;
  virtual Double_t     Z() const = 0;

  virtual Bool_t       PropagateTo(Double_t x, Double_t* length = NULL) = 0;
  virtual Bool_t       RotateTo(Double_t alpha) = 0;
  Bool_t               RotateBy(Double_t dAlpha)
    {return RotateTo(Alpha() + dAlpha);};
  virtual Bool_t       RotateAndPropagateTo(Double_t alpha, Double_t x, 
					    Double_t* length);
  virtual Bool_t       CorrectForMaterial(Double_t dAngle2, 
					  Double_t dPrel) = 0;
  virtual Bool_t       GetProlongationAt(Double_t x, Double_t& y, 
					 Double_t& z) const = 0;
  virtual Double_t     GetXAtVertex(Double_t x = 0, Double_t y = 0) const = 0;
  virtual Double_t     GetDsdx() const;

  virtual Double_t     GetPredictedChi2(const AliCluster* cluster) = 0;
  virtual Bool_t       Update(const AliCluster* cluster) = 0;

  virtual Double_t     Phi() const;
  virtual Double_t     SigmaPhi() const;
  virtual Double_t     Theta() const;
  virtual Double_t     SigmaTheta() const;
  virtual Double_t     Eta() const;
  virtual Double_t     Px() const;
  virtual Double_t     Py() const;
  virtual Double_t     Pz() const;
  virtual Double_t     Pt() const;
  virtual Double_t     SigmaPt() const;
  virtual Double_t     P() const;
  virtual TVector3     Momentum() const;
  virtual TVector3     Position() const;
  virtual TVector3     PositionAt(Double_t x) const;

  ClassDef(AliTrackParam, 1)
};

#endif
