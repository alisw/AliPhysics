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
  virtual ~AliACORDEhit();

  Int_t GetModule() const {return fModule;}
  Int_t GetPlastic() const {return fPlastic;}

  Float_t TrackId()   const {return fTrackId;}
  Float_t GetTime() const {return fTime;}
  Float_t Px()           const {return fPx;}
  Float_t Py()           const {return fPy;}
  Float_t Pz()           const {return fPz;}
  Float_t Energy() const {return fEnergy;}
  Float_t Eloss()        const {return fEloss;}
  Float_t PolarAngle()   const;
  Float_t AzimuthAngle() const;
  Float_t TrkLength() const {return fTrkLength;}


protected:
  Int_t fModule;
  Int_t fPlastic;
  Float_t fTrackId;
  Float_t fTime;
  Float_t fPx;  
  Float_t fPy;  
  Float_t fPz;  
  Float_t fEloss;
  Float_t fEnergy;
  Float_t fTrkLength;

private:
  ClassDef(AliACORDEhit,1)  // Hit for ACORDE (ACORDE)
};

typedef AliACORDEhit AliCRThit; // for backward compatibility

#endif /* ALIACORDEHIT_H */
