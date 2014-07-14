#ifndef ALIVVTRACKPOINT_H
#define ALIVVTRACKPOINT_H

#include "Rtypes.h"
class TGeoRotation;
class AliTrackPoint;


class AliVVtrackPoint {
  public:
  AliVVtrackPoint() {}
  virtual ~AliVVtrackPoint() {}


  virtual void     SetXYZ(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/, const Float_t *cov = 0) {if (cov) return;}
  virtual void     SetXYZ(const Float_t* /*xyz*/, const Float_t *cov = 0) {if (cov) return;}
  virtual void     SetCov(const Float_t* /*cov*/) {}
  virtual void     SetVolumeID(UShort_t /*volid*/) {}
  virtual void     SetCharge(Float_t /*charge*/) {}
  virtual void     SetDriftTime(Float_t /*time*/) {}
  virtual void     SetChargeRatio(Float_t /*ratio*/) {}
  virtual void     SetClusterType(Int_t /*clutype*/) {}
  virtual void     SetExtra(Bool_t flag=kTRUE) {if (flag) return;}

  virtual Float_t  GetX() const { return 0.; }
  virtual Float_t  GetY() const { return 0.; }
  virtual Float_t  GetZ() const { return 0.; }
  virtual void     GetXYZ(Float_t* /*xyz*/, Float_t* cov = 0) const {if (cov) return;}
  virtual const Float_t* GetCov() const { return NULL; }
  virtual UShort_t GetVolumeID() const { return 0; }
  virtual Float_t  GetCharge() const { return 0.; }
  virtual Float_t  GetDriftTime() const { return 0.;}
  virtual Float_t  GetChargeRatio() const { return 0.;}
  virtual Int_t    GetClusterType() const { return 0;}
  virtual Bool_t   IsExtra() const { return kFALSE;}

  virtual Float_t  GetResidual(const AliVVtrackPoint& /*p*/, Bool_t weighted = kFALSE) const {if (weighted) return 0.; else return 0.;}
  virtual Bool_t   GetPCA(const AliVVtrackPoint& /*p*/, AliVVtrackPoint& /*out*/) const {return kFALSE;}

  virtual Float_t  GetAngle() const {return 0.;}
  virtual Bool_t   GetRotMatrix(TGeoRotation& /*rot*/) const {return kFALSE;}
  //virtual void SetAlignCovMatrix(const TMatrixDSym& /*alignparmtrx*/) {}

  //virtual AliTrackPoint& Rotate(Float_t /*alpha*/) const {return *this;}
  //virtual AliTrackPoint& MasterToLocal() const {return *this;}

  virtual void     Print(Option_t*) const {}

  ClassDef(AliVVtrackPoint, 1);
};

#endif
