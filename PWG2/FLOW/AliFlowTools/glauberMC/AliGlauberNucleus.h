#ifndef ALIGLAUBERNUCLEUS_H
#define ALIGLAUBERNUCLEUS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//
//  AliGlauberNucleus
//  support classfor Glauber MC
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

class TNamed;
class TObjArray;
class TF1;

class AliGlauberNucleus : public TNamed {
private:
   Int_t      fN;          //Number of nucleons
   Double_t   fR;          //Parameters of function
   Double_t   fA;          //Parameters of function
   Double_t   fW;          //Parameters of function
   Double_t   fMinDist;    //Minimum separation distance
   Int_t      fF;          //Type of radial distribution
   Int_t      fTrials;     //Store trials needed to complete nucleus
   TF1*       fFunction;   //Probability density function rho(r)
   TObjArray* fNucleons;   //Array of nucleons

   void       Lookup(Option_t* name);

public:
   AliGlauberNucleus(Option_t* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
   virtual ~AliGlauberNucleus();

   using      TObject::Draw;
   void       Draw(Double_t xs, Int_t col);
   Int_t      GetN()             const {return fN;}
   Double_t   GetR()             const {return fR;}
   Double_t   GetA()             const {return fA;}
   Double_t   GetW()             const {return fW;}
   TObjArray *GetNucleons()      const {return fNucleons;}
   Int_t      GetTrials()        const {return fTrials;}
   void       SetN(Int_t in)           {fN=in;}
   void       SetR(Double_t ir);
   void       SetA(Double_t ia);
   void       SetW(Double_t iw);
   void       SetMinDist(Double_t min) {fMinDist=min;}
   void       ThrowNucleons(Double_t xshift=0.);

   ClassDef(AliGlauberNucleus,1)
};

#endif
