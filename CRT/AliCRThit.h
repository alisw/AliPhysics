#ifndef ALICRTHIT_H
#define ALICRTHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//                                            //
//  Hit class for CRT                         //
//  Interface                                 //
//  Getters, Setters and member variables     //
//  declared here                             //
//                                            //
////////////////////////////////////////////////

#include "AliHit.h"

class AliCRThit : public AliHit {
  
public:
  AliCRThit();
  AliCRThit(Int_t shunt, Int_t track, Int_t* vol, Float_t *hits);
  AliCRThit(const AliCRThit & hit);
  AliCRThit& operator= (const AliCRThit& hit);
  virtual ~AliCRThit() {}

  // getters for AliCRThit object
  Float_t GetId()   const {return fId;}
  Float_t GetX()     const {return fX;}
  Float_t GetY()     const {return fY;}
  Float_t GetZ()     const {return fZ;}
  Float_t GetPx()     const {return fPx;}
  Float_t GetPy()     const {return fPy;}
  Float_t GetPz()     const {return fPz;}
  Float_t GetMedium()     const {return fMedium;}
  Float_t GetELoss()   const {return fELoss;}
  Float_t GetCRT() const {return fCRTh;}
  Float_t GetCRTMod() const {return fCRTMod;}
  Float_t GetCRTMag() const {return fCRTMag;}
  Float_t GetCRTRICH() const {return fCRTRICH;}
  Float_t GetCRTTPC() const {return fCRTTPC;}

  const Int_t* GetVolume() const {return fVolume;}

protected:
  Int_t   fVolume[5];
  Int_t   fCopy;
  Float_t fId;     
  Float_t fPx;          //
  Float_t fPy;          //
  Float_t fPz;          //
  Float_t fMedium;      //
  Float_t fELoss;       //
  Float_t fCRTh;
  Float_t fCRTMod;
  Float_t fCRTMag;
  Float_t fCRTRICH;
  Float_t fCRTTPC;

private:
  ClassDef(AliCRThit,1)  // Hit for CRT (ACORDE)

};

#endif /* ALICRTHIT_H */
