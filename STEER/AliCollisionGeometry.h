#ifndef ALICOLLISIONGEOMETRY_H
#define ALICOLLISIONGEOMETRY_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>

class AliCollisionGeometry
{
public:
    AliCollisionGeometry();
    virtual ~AliCollisionGeometry(){;}
    // Getters
    Float_t ImpactParameter()   {return fImpactParameter;}
    Int_t   HardScatters() {return fNHardScatters;}
    Int_t   ProjectileParticipants()  {return fNProjectileParticipants;}
    Int_t   TargetParticipants()      {return fNTargetParticipants;}
    Int_t   Spectatorsn()	{return fSpecn;}
    Int_t   Spectatorsp()	{return fSpecp;}
    Int_t   NN()    {return fNNColl;}
    Int_t   NNw()   {return fNNwColl;}
    Int_t   NwN()   {return fNwNColl;}
    Int_t   NwNw()  {return fNwNwColl;}
    // Setters
    void SetImpactParameter(Float_t b)     {fImpactParameter=b;}
    void SetHardScatters(Int_t n)  {fNHardScatters=n;}
    void SetParticipants(Int_t np, Int_t nt)
	{fNProjectileParticipants=np, fNTargetParticipants=nt;}
    void SetCollisions(Int_t nn, Int_t nnw, Int_t nwn, Int_t nwnw)
	{fNNColl=nn, fNNwColl=nnw, fNwNColl=nwn,  fNwNwColl=nwnw;}
    void SetSpectators(Int_t nspecn, Int_t nspecp)
	{fSpecn=nspecn, fSpecp=nspecp;}
 protected:
    Int_t   fNHardScatters;            // Number of hard scatterings
    Int_t   fNProjectileParticipants;  // Number of projectiles participants
    Int_t   fNTargetParticipants;      // Number of target participants
    Int_t   fNNColl;                   // Number of N-N collisions
    Int_t   fNNwColl;                  // Number of N-Nwounded collisions
    Int_t   fNwNColl;                  // Number of Nwounded-N collisons
    Int_t   fNwNwColl;                 // Number of Nwounded-Nwounded collisions
    Int_t   fSpecn;                    // Number of spectators neutrons
    Int_t   fSpecp;                    // Number of spectators protons
    Float_t fImpactParameter;          // Impact Parameter

  ClassDef(AliCollisionGeometry,1)     // Collision Geometry
};
#endif






