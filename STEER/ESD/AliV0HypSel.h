#ifndef ALIV0HYPSEL_H
#define ALIV0HYPSEL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    V0 hypothesis selection params
// Used to check if given V0 is candidate for given hypothesis
//------------------------------------------------------------------

#include "TNamed.h"


//_____________________________________________________________________________
class AliV0HypSel : public TNamed {
public:

  AliV0HypSel();
  AliV0HypSel(const AliV0HypSel& src);
  AliV0HypSel(const char *name, float m0,float m1, float mass, float sigma, float nsig, float margin);
  void Validate();

  float GetM0()     const {return fM0;}
  float GetM1()     const {return fM1;}
  float GetMass()   const {return fMass;}
  float GetSigmaM() const {return fSigmaM;}
  float GetNSigma() const {return fNSigma;}
  float GetMarginAdd() const {return fMarginAdd;}
  float GetMassMargin(float pT) const {return fNSigma*fSigmaM*(1.+pT)+fMarginAdd;}

  virtual void Print(const Option_t *) const;
  
private:
  Float_t fM0;         // mass of the 1st prong
  Float_t fM1;         // mass of the 2nd prong
  Float_t fMass ;      // expected V0 mass
  Float_t fSigmaM;     // rough sigma estimate for sigmaMass = fSigmaM*(1+Pt) parameterization
  Float_t fNSigma;     // number fSigmaM to apply
  Float_t fMarginAdd;  // additional additive safety margin
  
  ClassDef(AliV0HypSel,1)  // V0 Hypothesis selection
};


#endif


