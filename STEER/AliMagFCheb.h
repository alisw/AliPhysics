#ifndef ALIMAGFCHEB_H
#define ALIMAGFCHEB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Interface between the AliMagWrapCheb and AliMagF: set of magnetic field data + Tosca 
// parameterization by Chebyshev polynomials
// 
// Author: ruben.shahoyan@cern.ch
//

#include "AliMagFC.h"
class AliMagWrapCheb;


class AliMagFCheb : public AliMagFC
{
public:
  enum constants {k2kG, k4kG, k5kG};
  AliMagFCheb();
  AliMagFCheb(const char *name, const char *title, Int_t integ,
	      Float_t factor=1, Float_t fmax=15, Int_t map = k2kG,
	      Bool_t dipoleON = kTRUE,
	      const char* path="$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  AliMagFCheb(const AliMagFCheb& maps);             
  AliMagFCheb& operator=(const AliMagFCheb& maps);
  virtual ~AliMagFCheb();
  //
  virtual void Field(const Float_t *x, Float_t *b)              const;
  virtual void Field(const Double_t *x, Double_t *b)            const;
  virtual void GetTPCInt(const Float_t *xyz, Float_t *b)        const;
  virtual void GetTPCIntCyl(const Float_t *rphiz, Float_t *b)   const;
  //
  AliMagWrapCheb* GetMeasuredMap()                              const {return fMeasuredMap;}
  void SetMeasuredMap(AliMagWrapCheb* parm);
  virtual Float_t SolenoidField()                               const {return -Factor()*fSolenoid;}
  //
 protected:
  AliMagWrapCheb* fMeasuredMap;     // Measured part of the field map
  Float_t         fSolenoid;        // Solenoid field setting
  //   
  ClassDef(AliMagFCheb, 2)       // Class for all Alice MagField wrapper for measured data + Tosca parameterization
};


#endif
