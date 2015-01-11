#ifndef ALIPMDMEANSM_H
#define ALIPMDMEANSM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;
class AliPMD;

class AliPMDMeanSm: public TNamed
{
 public:
  AliPMDMeanSm();
  AliPMDMeanSm(const char* name);
  AliPMDMeanSm(const AliPMDMeanSm &meanda);
  AliPMDMeanSm& operator= (const AliPMDMeanSm &meanda);
  virtual ~AliPMDMeanSm();
  void    Reset();
  void    SetMeanSm(Int_t det, Int_t smn,Float_t meansm);
  Float_t GetMeanSm(Int_t det, Int_t smn) const;
  virtual void Print(Option_t *) const;
  
 protected:

  enum
      {
	kDet = 2,        // Number of plane
	kModule = 24     // Modules per plane
      };
  Float_t fMeanSm[kDet][kModule];

  ClassDef(AliPMDMeanSm,0) // calibration class for gainfactors
};
#endif
