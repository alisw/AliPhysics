#ifndef ALIGENHIJINGEVENTHEADER_H
#define ALIGENHIJINGEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenHijingEventHeader : public AliGenEventHeader
{
 public:

  AliGenHijingEventHeader(const char* name){;}
  virtual ~AliGenHijingEventHeader() {}
  // Getters
  Float_t TotalEnergy()  {return fTotalEnergy;} 
  Int_t   HardScatters() {return fNHardScatters;}
  Int_t   ProjectileParticipants()  {return fNProjectileParticipants;}
  Int_t   TargetParticipants()      {return fNTargetParticipants;}	  
  Int_t   NN()    {return fNNColl;}
  Int_t   NNw()   {return fNNwColl;}
  Int_t   NwN()   {return fNwNColl;}
  Int_t   NwNw()  {return fNwNwColl;}
  // Setters
  void SetTotalEnergy(Float_t energy)  {fTotalEnergy=energy;}
  void SetHardScatters(Int_t n)  {fNHardScatters=n;}
  void SetParticipants(Int_t np, Int_t nt)
      {fNProjectileParticipants=np, fNTargetParticipants=nt;}
  void SetCollisions(Int_t nn, Int_t nnw, Int_t nwn, Int_t nwnw)
      {fNNColl=nn, fNNwColl=nnw, fNwNColl=nwn,  fNwNwColl=nwnw;}
  
protected:
  Float_t fTotalEnergy;              // Total energy of produced particles
  Int_t   fNHardScatters;            // Number of hard scatterings
  Int_t   fNProjectileParticipants;  // Number of projectiles participants
  Int_t   fNTargetParticipants;      // Number of target participants
  Int_t   fNNColl;                   // Number of N-N collisions
  Int_t   fNNwColl;                  // Number of N-Nwounded collisions
  Int_t   fNwNColl;                  // Number of Nwounded-N collisons
  Int_t   fNwNwColl;                 // Number of Nwounded-Nwounded collisions
  
  
  ClassDef(AliGenHijingEventHeader,1) // Event header for hijing event
};

#endif
