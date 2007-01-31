#ifndef ALIGENMC_H
#define ALIGENMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Base class for generators using external MC generators.
// For example AliGenPythia using Pythia.
// Provides basic functionality: setting of kinematic cuts on 
// decay products and particle selection.
// andreas.morsch@cern.ch

class TClonesArray;
class TParticle;
#include <TArrayI.h>   
#include <TString.h>

class AliGeometry;
#include "AliDecayer.h"
#include "AliGenerator.h"

class AliGenMC : public AliGenerator
{
 public:
    AliGenMC();
    AliGenMC(Int_t npart);
    virtual ~AliGenMC();
    virtual void Init();
    virtual void SetForceDecay(Decay_t decay = kAll) {fForceDecay = decay;}
    virtual void SetCutOnChild(Int_t flag = 0) {fCutOnChild = flag;}
    virtual void SetChildMomentumRange(Float_t pmin = 0, Float_t pmax = 1.e10)
	{fChildPMin = pmin; fChildPMax = pmax;}
    virtual void SetChildPtRange(Float_t ptmin = 0, Float_t ptmax = 20.)
	{fChildPtMin = ptmin; fChildPtMax = ptmax;}
    virtual void SetChildPhiRange(Float_t phimin = 0., Float_t phimax = 360.)
	{fChildPhiMin = TMath::Pi()*phimin/180;
	fChildPhiMax  = TMath::Pi()*phimax/180;}
    virtual void SetChildThetaRange(Float_t thetamin = 0, Float_t thetamax = 180)
	{fChildThetaMin = TMath::Pi()*thetamin/180;
	fChildThetaMax  = TMath::Pi()*thetamax/180;}
    virtual void SetChildYRange(Float_t ymin = -12, Float_t ymax = 12)
	{fChildYMin = ymin;
	fChildYMax  = ymax;}
    virtual void SetMaximumLifetime(Float_t time = 1.e-15) {fMaxLifeTime = time;}
   
    virtual void SetGeometryAcceptance(AliGeometry * GeometryAcceptance=0) {fGeometryAcceptance = GeometryAcceptance;}

    virtual void SetPdgCodeParticleforAcceptanceCut(Int_t PdgCodeParticleforAcceptanceCut=0) {fPdgCodeParticleforAcceptanceCut = PdgCodeParticleforAcceptanceCut;}

    virtual void SetNumberOfAcceptedParticles(Int_t NumberOfAcceptedParticles=2) {fNumberOfAcceptedParticles = NumberOfAcceptedParticles;}
    
    virtual Bool_t CheckAcceptanceGeometry(Int_t np, TClonesArray* particles);
    virtual void   SetProjectile(TString proj="P", Int_t a = 1, Int_t z = 1)
	{fProjectile = proj; fAProjectile = a; fZProjectile = z;}    
    virtual void   SetTarget(TString tar="P", Int_t a = 1, Int_t z = 1)
	{fTarget = tar; fATarget = a; fZTarget = z;}
    virtual void   SetCrossingAngle(Float_t phiX, Float_t phiY) {fXingAngleX = phiX; fXingAngleY = phiY;}
    virtual void Boost();

 protected:
    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip) const;
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip) const;
    // all kinematic selection cuts go here 
    Bool_t KinematicSelection(TParticle *particle, Int_t flag) const;
    Int_t  CheckPDGCode(Int_t pdgcode) const;

 protected:
    TClonesArray* fParticles;   //!Particle  List
    TArrayI     fParentSelect;  //!Parent particles to be selected 
    TArrayI     fChildSelect;   //!Decay products to be selected
    Int_t       fCutOnChild;    // Cuts on decay products (children)  are enabled/disabled
    Float_t     fChildPtMin;    // Children minimum pT
    Float_t     fChildPtMax;    // Children maximum pT
    Float_t     fChildPMin;     // Children minimum p
    Float_t     fChildPMax;     // Children maximum p
    Float_t     fChildPhiMin;   // Children minimum phi
    Float_t     fChildPhiMax;   // Children maximum phi
    Float_t     fChildThetaMin; // Children minimum theta
    Float_t     fChildThetaMax; // Children maximum theta
    Float_t     fChildYMin;     // Children minimum y
    Float_t     fChildYMax;     // Children maximum y
    Float_t     fXingAngleX;    // Crossing angle X
    Float_t     fXingAngleY;    // Crossing angle Y    
    Decay_t     fForceDecay;    // Decay channel forced
    Float_t     fMaxLifeTime;   // Maximum lifetime for unstable particles
    Int_t       fAProjectile;   // Projectile A
    Int_t       fZProjectile;   // Projectile Z
    Int_t       fATarget;       // Target A
    Int_t       fZTarget;       // Target Z
    TString     fProjectile;    // Projectile
    TString     fTarget;        // Target
    Double_t    fDyBoost;       // dy for boost into lab frame
    AliGeometry * fGeometryAcceptance; // Geometry to which particles must be simulated
    Int_t       fPdgCodeParticleforAcceptanceCut;  // Abs(PDG Code) of the particle to which the GeometryAcceptance must be applied
    Int_t       fNumberOfAcceptedParticles;  // Number of accepted particles in GeometryAcceptance with the right Abs(PdgCode) 

 private:
    AliGenMC(const AliGenMC &MC);
    AliGenMC & operator=(const AliGenMC & rhs);
    
    ClassDef(AliGenMC,5)       // AliGenerator implementation for generators using MC methods
};
#endif





