#ifndef ALITRDSIMPLEGEN_H
#define ALITRDSIMPLEGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Particle generator for the simplescopic TRD simulator                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TObject.h>

class AliTRDsimpleGen : public TObject {
 
 public:     

  AliTRDsimpleGen();
  AliTRDsimpleGen(const AliTRDsimpleGen &g); 
                                                                                
  virtual ~AliTRDsimpleGen();
  AliTRDsimpleGen &operator=(const AliTRDsimpleGen &g);    

  virtual void         Copy(TObject &g) const;
  virtual void         NewParticle(Int_t ievent);

  virtual void         SetMomentum(Double_t min, Double_t max) { fMomMin = min;
                                                                 fMomMax = max; };
  virtual void         SetPdg(Int_t pdg)                       { fPdg    = pdg; };

 protected:

  Int_t           fPdg;              //  Particle PDG code
  Double_t        fMomMin;           //  Particle minimum momentum
  Double_t        fMomMax;           //  Particle maximum momentum

  ClassDef(AliTRDsimpleGen,1)        //  Particle generator for the simplified TRD simulator
 
};
#endif                                                                          
