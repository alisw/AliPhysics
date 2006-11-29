#ifndef ALIACORDEHIT_H
#define ALIACORDEHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//                                            //
//  Hit class for ACORDE                         //
//  Interface                                 //
//  Getters, Setters and member variables     //
//  declared here                             //
//                                            //
////////////////////////////////////////////////

#include "AliHit.h"

class AliACORDEhit : public AliHit {  
public:
  AliACORDEhit();
  AliACORDEhit(Int_t shunt, Int_t track, Int_t* vol, Float_t *hits);
  AliACORDEhit(const AliACORDEhit& hit);
  virtual ~AliACORDEhit();

  AliACORDEhit& operator=(const AliACORDEhit& hit);
  Bool_t operator==(const AliACORDEhit& hit);
  Bool_t operator<(const AliACORDEhit& hit);

  Float_t ParticleId()   const {return fId;}
  Float_t Px()           const {return fPx;}
  Float_t Py()           const {return fPy;}
  Float_t Pz()           const {return fPz;}
  Float_t Eloss()        const {return fEloss;}
  Float_t Medium()       const {return fMedium;}
  Float_t Energy()       const;
  Float_t PolarAngle()   const;
  Float_t AzimuthAngle() const;

protected:
  Float_t fId;     //
  Float_t fPx;     //
  Float_t fPy;     //
  Float_t fPz;     //
  Float_t fEloss;  //
  Float_t fMedium; //

private:
  ClassDef(AliACORDEhit,1)  // Hit for ACORDE (ACORDE)
};
#endif /* ALIACORDEHIT_H */
