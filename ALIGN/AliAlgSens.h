#ifndef ALIALGSENS_H
#define ALIALGSENS_H

#include "AliAlgVol.h"
#include <TMath.h>

class AliAlgDet;
class AliAlgPoint;
class TObjArray;
class AliExternalTrackParam;


/*--------------------------------------------------------
  End-chain alignment volume in detector branch, where the actual measurement is done.
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgSens : public AliAlgVol
{
 public:
  //
  AliAlgSens(const char* name=0, Int_t vid=0, Int_t iid=0);
  virtual ~AliAlgSens();
  //
  virtual void AddChild(AliAlgVol*);
  //
  void       SetDetector(const AliAlgDet* det)         {fDet = det;}
  const AliAlgDet* GetDetector()                 const {return fDet;}
  //
  void  IncrementStat()                                {fNProcPoints++;}
  //
  // derivatives calculation
  virtual void DPosTraDParGeom(const double *tra, double* deriv,const AliAlgVol* parent=0) const;
  //
  virtual void DPosTraDParGeomLOC(const double *tra, double* deriv) const;
  virtual void DPosTraDParGeomTRA(const double *tra, double* deriv) const;
  virtual void DPosTraDParGeomLOC(const double *tra, double* deriv, const AliAlgVol* parent) const;
  virtual void DPosTraDParGeomTRA(const double *tra, double* deriv, const AliAlgVol* parent) const;
  //
  void GetModifiedMatrixT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const;
  void GetModifiedMatrixT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const;
  //
  void            SetAddError(double y, double z)            {fAddError[0]=y;fAddError[1]=z;}
  const Double_t* GetAddError()                        const {return fAddError;} 
  //
  virtual void   PrepareMatrixT2L();
  //
  virtual void   SetTrackingFrame();
  virtual Bool_t IsSensor()                       const {return kTRUE;}
  virtual void   Print(const Option_t *opt="")    const;
  //
  virtual void UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const;
  //
 protected:
  //
  virtual Bool_t  IsSortable()                         const {return kTRUE;}
  virtual Int_t   Compare(const TObject* a)            const;
  //
  // --------- dummies -----------
  AliAlgSens(const AliAlgSens&);
  AliAlgSens& operator=(const AliAlgSens&);
  //
 protected:
  //
  Double_t fAddError[2];              // additional error increment for measurement
  const AliAlgDet* fDet;              // pointer on detector
  //
  ClassDef(AliAlgSens,1)
};


#endif
