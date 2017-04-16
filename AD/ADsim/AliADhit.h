// -*- C++ -*-
#ifndef ALIADHIT_H
#define ALIADHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliADhit.h  $ */

////////////////////////////////////////////////
//                                            //
//  Manager and hits classes for set : AD     //
//                                            //
////////////////////////////////////////////////

#include "AliHit.h"
#include "TObjArray.h"


class AliADhit : public AliHit {
public:
  AliADhit();
  AliADhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual          ~AliADhit();
  Int_t     GetModule()	   const { return fModule; }
  Float_t   GetEnergy()    const { return fEk; }
  Float_t   GetPt()        const { return fPt; }
  Float_t   GetPx()        const { return fPx; }
  Float_t   GetPy()        const { return fPy; }
  Float_t   GetPz()        const { return fPz; }
  Float_t   GetTof()       const { return fTof; };
  Float_t   GetTrackleng() const { return fTleng; }
  Float_t   GetEloss()     const { return fEloss; }
  Int_t     GetNphot()     const { return fNphot; }
  Int_t     GetCell()      const { return fCell; }
  Int_t     GetTrackStatus()    const { return fPrimary; }
  Int_t     GetTrackPDG()       const { return fPDG; }
  Int_t     GetTrackPDGMother() const { return fPDGMother; }

protected:
  Int_t fModule;
  Float_t   fEk;            // kinetic energy of the entering particle
  Float_t   fPt;            // Local transverse momentum of the particle
  Float_t   fPx;            // Local Px momentum of the particle
  Float_t   fPy;            // Local Py momentum of the particle
  Float_t   fPz;            // Local Pz momentum of the particle
  Float_t   fTof;           // Particle time of flight wrt vertex
  Float_t   fTleng;         // Track length in ADA or ADD detector
  Float_t   fEloss;         // Energy loss in AD detector
  Int_t     fNphot;         // Number of photons created by current hit
  Int_t     fCell;          // Scintillator cell (Sector Number) (ADA1 = 10-14, ADA2 = 20-24, ADC1 = 30-34, ADC2 = 40-44)

  // The following are for fast analysis, but nor really needed, can be retrive from AliStack with track ID on fTrack
  Int_t   fPrimary;       // flag (track primary or secondary)
  Int_t   fPDG;           // PDG code
  Int_t   fPDGMother;     // PDG code of the mother

private:
  ClassDef(AliADhit,1); //  Hits for detector AD
};

#endif
