#ifndef ALIRICHHIT_H
#define ALIRICHHIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliHit.h"

class AliRICHHit : public AliHit {
 public:
    Int_t     fChamber;       // Chamber number
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas   
    Float_t   fPHfirst;       // first padhit
    Float_t   fPHlast;        // last padhit
    Float_t   fLoss;          // did it hit the freon?
    Float_t   fMomX;          // x Momentum at photochatode entry point
    Float_t   fMomY;          // y Momentum at photochatode entry point
    Float_t   fMomZ;          // z Momentum at photochatode entry point
    Float_t   fNPads;         // Pads hit
    Float_t   fCerenkovAngle; // Dummy cerenkov angle
    Float_t   fMomFreoX;      // x Momentum at freon entry point
    Float_t   fMomFreoY;      // y Momentum at freon entry point
    Float_t   fMomFreoZ;      // z Momentum at freon entry point		   

 public:
    AliRICHHit() {}
    AliRICHHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliRICHHit() {}
    
    Int_t   Chamber()             {return fChamber;}
    Float_t Particle()            {return fParticle;}    
    Float_t Theta()               {return fTheta;}
    Float_t Phi()                 {return fPhi;}
    Float_t Tlength()             {return fTlength;}
    Float_t Eloss()               {return fEloss;}
    Float_t   PHfirst()           {return fPHfirst;}
    Float_t   PHlast()            {return fPHlast;}
    Float_t MomX()                {return fMomX;}
    Float_t MomY()                {return fMomY;}
    Float_t MomZ()                {return fMomZ;}
    Float_t CerenkovAngle()       {return fMomX;}
    Float_t MomFreoX()            {return fMomX;}
    Float_t MomFreoY()            {return fMomY;}
    Float_t MomFreoZ()            {return fMomZ;}
    
    
    ClassDef(AliRICHHit,1)  //Hits object for set:RICH
};
#endif
