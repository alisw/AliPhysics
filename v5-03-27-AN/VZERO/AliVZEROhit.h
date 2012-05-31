#ifndef ALIVZEROHIT_H
#define ALIVZEROHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//                                            //
//  Manager and hits classes for set : VZERO  //
//                                            //
////////////////////////////////////////////////
 
#include "AliHit.h"
#include "TObjArray.h"
#include "TArrayF.h"
 
class AliVZEROhit : public AliHit {
 
public:
  AliVZEROhit();
  AliVZEROhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliVZEROhit() {};
  
  Int_t   Volume()  const {return fVolume;};
  Int_t   CopyNumber()    const {return fCopy;};
  Float_t TrackPiD() const {return fTrackPiD;};
  Float_t Tof()   const {return fTof;};
  Float_t Charge() const {return fCharge;};
  Float_t RingNumber() const {return fRingNumber;};
  Float_t Pt()    const {return fPt;};
  Float_t Pmom()  const {return fPmom;};
  Float_t Px()    const {return fPx;};
  Float_t Py()    const {return fPy;};
  Float_t Pz()    const {return fPz;};
  Float_t Vx()    const {return fVx;};
  Float_t Vy()    const {return fVy;};
  Float_t Vz()    const {return fVz;};
  Float_t Eloss() const {return fEloss;}
  Float_t Tleng() const {return fTleng;}
  Int_t   Nphot() const {return fNphot;}
  Int_t   Cell()  const {return fCell;}
 
private:
  Int_t   fVolume;        // Current volume ID
  Int_t   fCopy;          // Current copy number
  Float_t fTrackPiD;      // Track PiD 
  Float_t fTof;           // Particle time of flight wrt vertex
  Float_t fCharge;        // Particle charge
  Float_t fTheta;         // Incident theta angle in degrees 
  Float_t fPhi;           // Incident phi angle in degrees
  Float_t fRingNumber;    // RingNumber
  
  Float_t fPt;            // Local transverse momentum of the particle
  Float_t fPmom;          // Local P momentum of the particle
  Float_t fPx;            // Local Px momentum of the particle
  Float_t fPy;            // Local Py momentum of the particle
  Float_t fPz;            // Local Pz momentum of the particle
  Float_t fVx;            // Vertex x coordinate  
  Float_t fVy;            // Vertex y coordinate  
  Float_t fVz;            // Vertex z coordinate    
  Float_t fEloss;         // Energy loss  in VZERO detector
  Float_t fTleng;         // Track length in VZERO detector
  Int_t   fNphot;         // Number of photons created by current hit 
  Int_t   fCell;          // Scintillator cell number from 0 to 71 

  ClassDef(AliVZEROhit,2) //  Hits for detector VZERO
};
#endif
