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
 
class AliVZEROhit : public AliHit {
 
public:
  AliVZEROhit() {}
  AliVZEROhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliVZEROhit() {}
  virtual Int_t Volume() {return fVolume;};
  virtual Int_t Copy() {return fCopy;};
  virtual Float_t TrackPiD() {return fTrackPiD;};
  virtual Float_t Tof() {return fTof;};
  virtual Float_t Charge() {return fCharge;};
  virtual Float_t RingNumber() {return fRingNumber;};
  
  virtual Float_t Pt()   {return fPt;};
  virtual Float_t Pmom() {return fPmom;};
  virtual Float_t Px()   {return fPx;};
  virtual Float_t Py()   {return fPy;};
  virtual Float_t Pz()   {return fPz;};
  
  Float_t Eloss()        {return fEloss;}
  Float_t Tleng()        {return fTleng;}
  
private:
  Int_t   fVolume;                // Current volume ID
  Int_t   fCopy;                  // Copy number
  Float_t fTrackPiD;              // Root particle ID 
  Float_t fTof;                   // Time of flight wrt vertex
  Float_t fCharge;                // Charge of particle
  Float_t fTheta; 
  Float_t fPhi;
  Float_t fRingNumber;
  
  Float_t fPt;
  Float_t fPmom;
  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  
  Float_t fEloss;         //  energy loss  in VZERO detector
  Float_t fTleng;         //  track length in VZERO detector
  
    
  ClassDef(AliVZEROhit,1)  //Hits for detector VZERO
};
#endif
