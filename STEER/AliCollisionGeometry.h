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
    Float_t ReactionPlaneAngle() {return fReactionPlaneAngle;}
    Int_t   HardScatters() {return fNHardScatters;}
    Int_t   ProjectileParticipants()  {return fNProjectileParticipants;}
    Int_t   TargetParticipants()      {return fNTargetParticipants;}
    Int_t   ProjSpectatorsn()	{return fProjectileSpecn;}
    Int_t   ProjSpectatorsp()	{return fProjectileSpecp;}
    Int_t   TargSpectatorsn()	{return fTargetSpecn;	 }
    Int_t   TargSpectatorsp()	{return fTargetSpecp;	 }
    Int_t   NN()    {return fNNColl;}
    Int_t   NNw()   {return fNNwColl;}
    Int_t   NwN()   {return fNwNColl;}
    Int_t   NwNw()  {return fNwNwColl;}
    // Setters
    void SetImpactParameter(Float_t b)     {fImpactParameter=b;}
    void SetReactionPlaneAngle(Float_t phi)     {fReactionPlaneAngle = phi;}
    void SetHardScatters(Int_t n)  {fNHardScatters=n;}
    void SetParticipants(Int_t np, Int_t nt)
	{fNProjectileParticipants=np, fNTargetParticipants=nt;}
    void SetCollisions(Int_t nn, Int_t nnw, Int_t nwn, Int_t nwnw)
	{fNNColl=nn, fNNwColl=nnw, fNwNColl=nwn,  fNwNwColl=nwnw;}
    void SetSpectators(Int_t nprojspecn, Int_t nprojspecp, Int_t ntargspecn, Int_t ntargspecp)
	{fProjectileSpecn=nprojspecn, fProjectileSpecp=nprojspecp, 
	 fTargetSpecn=ntargspecn, fTargetSpecp=ntargspecp;}
 protected:
    Int_t   fNHardScatters;            // Number of hard scatterings
    Int_t   fNProjectileParticipants;  // Number of projectiles participants
    Int_t   fNTargetParticipants;      // Number of target participants
    Int_t   fNNColl;                   // Number of N-N collisions
    Int_t   fNNwColl;                  // Number of N-Nwounded collisions
    Int_t   fNwNColl;                  // Number of Nwounded-N collisons
    Int_t   fNwNwColl;                 // Number of Nwounded-Nwounded collisions
    Int_t   fProjectileSpecn;	       // Num. of spectator neutrons from projectile nucleus
    Int_t   fProjectileSpecp;	       // Num. of spectator protons from projectile nucleus
    Int_t   fTargetSpecn;    	       // Num. of spectator neutrons from target nucleus
    Int_t   fTargetSpecp;    	       // Num. of spectator protons from target nucleus
    Float_t fImpactParameter;          // Impact Parameter
    Float_t fReactionPlaneAngle;       // Reaction plane angle
    
  ClassDef(AliCollisionGeometry,3)     // Collision Geometry
};
#endif






