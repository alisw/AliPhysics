#ifndef VZEROHIT_H
#define VZEROHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO  //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"
#include "TObjArray.h"
#include "TArrayF.h"
#include "TMath.h"
 
class AliVZEROhit : public AliHit {
 
public:
  AliVZEROhit() {}
  AliVZEROhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliVZEROhit() {}
  inline Int_t GetVolume() {return fVolume;};
  inline Int_t GetCopy() {return fCopy;};
  inline Float_t GetX() {return fX;};
  inline Float_t GetY() {return fY;};
  inline Float_t GetZ() {return fZ;};
  inline Float_t GetXloc() {return fXloc;};
  inline Float_t GetYloc() {return fYloc;};
  inline Float_t GetZloc() {return fZloc;};
  inline Float_t GetEdep() {return fEdep;};
  inline Float_t GetEtot() {return fEtot;};
  inline Float_t GetTrackPiD() {return fTrackPiD;};
  inline Float_t GetParticle() {return fParticle;};
  inline Float_t GetTof() {return fTof;};
  inline Float_t IsTrackEntering() {return fIsTrackEntering;};
  inline Float_t IsTrackExiting() {return fIsTrackExiting;};
  inline Float_t GetCharge() {return fCharge;};
  inline Float_t IsCerenkov() {return fIsCerenkov;};
  inline Float_t GetMultiplicity() {return fMulti;};
  inline Float_t GetTheta() {return fTheta;};
  inline Float_t GetPhi() {return fPhi;};
  inline Float_t GetNGCerenkovs() {return fNGCerenkovs;};
  
public:
  Int_t   fVolume;                // Current volume ID
  Int_t   fCopy;                  // Copy number
  Float_t fXloc;                  // x coordinate in STRT coord
  Float_t fYloc;                  // y coordinate in STRT coord 
  Float_t fZloc;                  // z coordinate in STRT coord 
  Float_t fEdep;                  // Energy loss
  Float_t fEtot;                  // Total energy of particle 
  Float_t fTrackPiD;              // Root particle ID 
  Float_t fParticle;              // Geant particle ID 
  Float_t fTof;                   // Time of flight wrt vertex
  Float_t fIsTrackEntering;       // Entrance flag
  Float_t fIsTrackExiting;        // Exit flag
  Float_t fCharge;                // Charge of particle
  Float_t fIsCerenkov;            // Particle is a cerenkov photon
  Float_t fMulti;                 // Multiplicity of entering charged particles
  Float_t fTheta; 
  Float_t fPhi;
  Float_t fNGCerenkovs;
    
  ClassDef(AliVZEROhit,1)  //Hits for detector VZERO
};
#endif
