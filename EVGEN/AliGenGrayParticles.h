#ifndef ALIGENGRAYPARTICLES_H
#define ALIGENGRAYPARTICLES_H
/* Copyright(c) 198-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"
class AliGrayParticleModel;

class AliGenGrayParticles : public AliGenerator
{
public:
    AliGenGrayParticles();
    AliGenGrayParticles(Int_t npart);
    virtual ~AliGenGrayParticles();
    virtual void Init();
    virtual void Generate();
    virtual void SetPmax(Float_t pmax = 10.) {fPmax = pmax;}
    virtual void SetNominalCmsEnergy(Float_t energy = 14000.) {fCMS = energy;}
    virtual void SetTarget(Float_t a=208, Float_t z=82) {fATarget = a; fZTarget = z;}
    virtual void SetCharge(Int_t c = 1) {fCharge = c;}
    virtual void SetTemperature(Double_t t = 0.05) {fTemperature = t;}
    virtual void SetBetaSource(Double_t b = 0.05) {fBetaSource = b;}
    //
    virtual void SetGrayParticleModel(AliGrayParticleModel* model) 
	{fGrayParticleModel = model;}
    virtual Bool_t NeedsCollisionGeometry() {return kTRUE;}
    virtual void   SetCollisionGeometry(AliCollisionGeometry* geom)
	{fCollisionGeometry = geom;}
	    
 protected:
    void     GenerateSlow(Int_t charge, Double_t T, Double_t beta, Float_t* q);
    Double_t Maxwell(Double_t m, Double_t p, Double_t t);
    void     Lorentz(Double_t m, Double_t beta, Float_t* q);
 protected:
    Float_t  fCMS;         // Center of mass energy
    Float_t  fMomentum;    // Target nucleus momentum
    Float_t  fBeta;        // Target nucleus beta
    Float_t  fPmax;        // Maximum slow nucleon momentum
    Float_t  fATarget;     // Target nucleus mass number
    Float_t  fZTarget;     // Target nucleus charge number
    Int_t    fCharge;      // Slow nucleon charge
    Float_t  fTemperature; // Source Temperature
    Float_t  fBetaSource;  // Source beta
    //
    AliGrayParticleModel* fGrayParticleModel; // The gray particle model
  ClassDef(AliGenGrayParticles,1) // Gray Particle Generator
};
#endif






