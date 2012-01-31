#ifndef ALIGLAUBERNUCLEON_H
#define ALIGLAUBERNUCLEON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//
//  AliGlauberNucleon
//  support class for Glauber MC
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

class TNamed;

class AliGlauberNucleon : public TNamed {
  
private:
   Double32_t fX;            //Position of nucleon
   Double32_t fY;            //Position of nucleon
   Double32_t fZ;            //Position of nucleon
   Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
   Int_t      fNColl;        //Number of binary collisions

public:
   AliGlauberNucleon();
   virtual   ~AliGlauberNucleon() {}

   void       Collide()            {fNColl++;}
   Int_t      GetNColl()     const {return fNColl;}
   Double_t   GetX()         const {return fX;}
   Double_t   GetY()         const {return fY;}
   Double_t   GetZ()         const {return fZ;}
   Bool_t     IsInNucleusA() const {return fInNucleusA;}
   Bool_t     IsInNucleusB() const {return !fInNucleusA;}
   Bool_t     IsSpectator()  const {return !fNColl;}
   Bool_t     IsWounded()    const {return fNColl;}
   void       Reset()              {fNColl=0;}
   void       SetInNucleusA()      {fInNucleusA=1;}
   void       SetInNucleusB()      {fInNucleusA=0;}
   void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}

   ClassDef(AliGlauberNucleon,1)
};

#endif
