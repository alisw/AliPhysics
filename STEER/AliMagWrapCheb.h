#ifndef ALIMAGWRAPCHEB_H
#define ALIMAGWRAPCHEB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Wrapper for AliMagFCheb: set of magnetic field data + Tosca parameterization by Chebyshev polynomials
// 
// Author: ruben.shahoyan@cern.ch
//

#include "AliMagFC.h"
#include "AliMagFCheb.h"


class AliMagWrapCheb : public AliMagFC
{
public:
  enum constants {k2kG, k4kG, k5kG};
  AliMagWrapCheb();
  AliMagWrapCheb(const char *name, const char *title, Int_t integ,
		 Float_t factor=1, Float_t fmax=15, Int_t map = k2kG,
		 Bool_t dipoleON = kTRUE,
		 const char* path="$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  AliMagWrapCheb(const AliMagWrapCheb& maps);             
  AliMagWrapCheb& operator=(const AliMagWrapCheb& maps);
  virtual ~AliMagWrapCheb();
  //
  virtual void Field(float *x, float *b)                        const;
  virtual void Field(double *x, double *b)                      const;
  virtual void GetTPCInt(Float_t *xyz, Float_t *b)              const;
  virtual void GetTPCIntCyl(Float_t *rphiz, Float_t *b)         const;
  //
  AliMagFCheb* GetMeasuredMap()                                 const {return fMeasuredMap;}
  void SetMeasuredMap(AliMagFCheb* parm)                        {if (fMeasuredMap) delete fMeasuredMap; fMeasuredMap = parm;}
  virtual Float_t SolenoidField()                               const {return -Factor()*fSolenoid;}
  //
 protected:
  AliMagFCheb* fMeasuredMap;     // Measured part of the field map
  Float_t      fSolenoid;        // Solenoid field setting
  //   
  ClassDef(AliMagWrapCheb, 2)    // Class for all Alice MagField wrapper for measured data + Tosca parameterization
};


#endif
