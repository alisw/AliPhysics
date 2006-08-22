#ifndef ALITOFHITT0_H
#define ALITOFHITT0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            //
//  Hit class for TOF                         //
//  Interface                                 //
//  Getters, Setters and member variables     //
//  declared here                             //
//                                            //
////////////////////////////////////////////////
 
/* $Id$ */

#include "AliHit.h"

class AliTOFhitT0 : public AliHit {  
public:
  AliTOFhitT0();
  AliTOFhitT0(Int_t shunt, Int_t track, Int_t* vol, 
            Float_t *hits);
  AliTOFhitT0(const AliTOFhitT0 & hit) ;
  virtual ~AliTOFhitT0() {}
  // getters for AliTOFhitT0 object
  Int_t   GetSector() const {return fSector;}
  Int_t   GetPlate()  const {return fPlate;}
  Int_t   GetPadx()   const {return fPadx;}
  Int_t   GetPadz()   const {return fPadz;}
  Int_t   GetStrip()  const {return fStrip;}
  Float_t GetTof()    const {return fTof;}
  Float_t GetLen()    const {return fLenTof;}
  Float_t GetMom()    const {return fPmom;}
  Float_t GetPx()     const {return fPx;}
  Float_t GetPy()     const {return fPy;}
  Float_t GetPz()     const {return fPz;}
  Float_t GetDx()     const {return fDx;}
  Float_t GetDy()     const {return fDy;}
  Float_t GetDz()     const {return fDz;}
  Float_t GetIncA()   const {return fIncA;}
  Float_t GetEdep()   const {return fEdep;}
  
 protected:
  Int_t      fSector;  // number of sector 
  Int_t      fPlate;   // number of plate
  Int_t      fStrip;   // number of strip
  Int_t      fPadx;    // number of pad along x
  Int_t      fPadz;    // number of pad along z
  // X, Y and Z coordinates of the hit are defined on mother class
  // AliHit
  Float_t    fPx;      // px in TOF
  Float_t    fPy;      // py in TOF
  Float_t    fPz;      // pz in TOF
  Float_t    fPmom;    // P in TOF
  Float_t    fTof;     // Time of Flight
  Float_t    fLenTof;  // track length at the TOF
  Float_t    fDx;      // x of impact point in pad r.s.
  Float_t    fDy;      // y of impact point in pad r.s.
  Float_t    fDz;      // z of impact point in pad r.s.
  Float_t    fIncA;    // Incidence angle
  Float_t    fEdep;    // Energy lost in TOF sensitive layer
  
  ClassDef(AliTOFhitT0,1)  // Hit for Time Of Flight
    };

#endif /* ALITOFHITT0_H */
