#ifndef ALIVZEROHIT_H
#define ALIVZEROHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO  //
////////////////////////////////////////////////
 
#include "AliHit.h"
#include "TObjArray.h"
#include "TArrayF.h"
 
class AliVZEROhit : public AliHit {
 
public:
  AliVZEROhit() {}
  AliVZEROhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliVZEROhit() {};
  
  Int_t Volume()  {return fVolume;};
  Int_t Copy()    {return fCopy;};
  Float_t TrackPiD() {return fTrackPiD;};
  Float_t Tof()   {return fTof;};
  Float_t Charge() {return fCharge;};
  Float_t RingNumber() {return fRingNumber;};
  Float_t Pt()    {return fPt;};
  Float_t Pmom()  {return fPmom;};
  Float_t Px()    {return fPx;};
  Float_t Py()    {return fPy;};
  Float_t Pz()    {return fPz;};
  Float_t Vx()    {return fVx;};
  Float_t Vy()    {return fVy;};
  Float_t Vz()    {return fVz;};
  Float_t Eloss() {return fEloss;}
  Float_t Tleng() {return fTleng;}
 
private:
  Int_t   fVolume;                // Current volume ID
  Int_t   fCopy;                  // Copy number
  Float_t fTrackPiD;              // Root particle ID 
  Float_t fTof;                   // Time of flight wrt vertex
  Float_t fCharge;                // Charge of particle
  Float_t fTheta; 
  Float_t fPhi;
  Float_t fRingNumber;            // RingNumber
  
  Float_t fPt;
  Float_t fPmom;
  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Float_t fVx;            //  Vertex x coordinate  
  Float_t fVy;            //  Vertex y coordinate  
  Float_t fVz;            //  Vertex z coordinate    
  Float_t fEloss;         //  energy loss  in VZERO detector
  Float_t fTleng;         //  track length in VZERO detector
  
    
  ClassDef(AliVZEROhit,2)  //Hits for detector VZERO
};
#endif
