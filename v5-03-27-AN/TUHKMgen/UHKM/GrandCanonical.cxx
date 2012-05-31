//////////////////////////////////////////////////////////////////////////////////       
//                                                                              //
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna     //
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru  //
//                           November. 2, 2005                                  //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <TError.h>
#include <TMath.h>
#include "GrandCanonical.h"
#include "HankelFunction.h"
#include "UKUtility.h"
#include "ParticlePDG.h"
#include "DatabasePDG.h"


//_______________________________________________________________________________
GrandCanonical::GrandCanonical():
  fTemperature(-1111),
  fBaryonPotential(-1111),
  fStrangePotential(-1111),
  fElectroPotential(-1111),
  fNMax(-1111),
  fInitialized(kFALSE)
{
  //
  // default constructor
  //
}

//_______________________________________________________________________________
GrandCanonical::GrandCanonical(Int_t nmax, Double_t temperature, Double_t baryonPotential, Double_t strangePotential, Double_t electroPotential):
  fTemperature(temperature),
  fBaryonPotential(baryonPotential),
  fStrangePotential(strangePotential),
  fElectroPotential(electroPotential),
  fNMax(nmax),
  fInitialized(kTRUE)
{
  //
  // constructor
  //
}

//_______________________________________________________________________________
GrandCanonical::~GrandCanonical() {
//
// destructor
//
}

//_______________________________________________________________________________
void GrandCanonical::Temperature(Double_t value) {
  //
  // set temperature
  //
  fTemperature = value;
  if(fNMax!=-1111 && TMath::Abs(fBaryonPotential+1111)>1.e-10 && 
    TMath::Abs(fStrangePotential+1111)>1.e-10 && TMath::Abs(fElectroPotential+1111)>1.e-10)
    fInitialized = kTRUE;
}

//_______________________________________________________________________________
void GrandCanonical::BaryonPotential(Double_t value) {
  //
  // set baryo chemical potential
  //
  fBaryonPotential = value;
  if(fNMax!=-1111 && TMath::Abs(fTemperature+1111)>1.e-10 && 
    TMath::Abs(fStrangePotential+1111)>1.e-10 && TMath::Abs(fElectroPotential+1111)>1.e-10)
    fInitialized = kTRUE;
}

//_______________________________________________________________________________
void GrandCanonical::StrangePotential(Double_t value) {
  //
  // set strange potential
  //
  fStrangePotential = value;
  if(fNMax!=-1111 && TMath::Abs(fTemperature+1111)>1.e-10 &&
     TMath::Abs(fBaryonPotential+1111)>1.e-10 && TMath::Abs(fElectroPotential+1111)>1.e-10)
    fInitialized = kTRUE;
}

//_______________________________________________________________________________
void GrandCanonical::ElectroPotential(Double_t value) {
  //
  // set electro chemical potential
  //
  fElectroPotential = value;
  if(fNMax!=-1111 && TMath::Abs(fTemperature+1111)>1.e-10 &&
     TMath::Abs(fBaryonPotential+1111)>1.e-10 && TMath::Abs(fStrangePotential+1111)>1.e-10)
    fInitialized = kTRUE;
}

//_______________________________________________________________________________
void GrandCanonical::NMax(Int_t value) {
  //
  // set the number of iterations
  //
  fNMax = value;
  if(fTemperature!=-1111 && fBaryonPotential!=-1111 && fStrangePotential!=-1111 && fElectroPotential!=-1111)
    if(fNMax!=-1111 && TMath::Abs(fBaryonPotential+1111)>1.e-10 && 
    TMath::Abs(fStrangePotential+1111)>1.e-10 && TMath::Abs(fElectroPotential+1111)>1.e-10)
    fInitialized = kTRUE;
}

//_______________________________________________________________________________
Double_t GrandCanonical::ParticleEnergyDensity(ParticlePDG *const particle) {
  //
  // compute the energy density for a given particle
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::ParticleEnergyDensity", "GrandCanonical object not fully initialized!!");
  
  // Compute the particle energy density
  Double_t degFactor = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
  Double_t mass = particle->GetMass();                                                // PDG table mass
  Double_t d = Int_t(2.*particle->GetSpin()) & 1 ? 1. : -1;                                   // Bose-Einstein/Fermi-Dirac factor
  Double_t preFactor = (degFactor*mass*mass*fTemperature*fTemperature/kHbarc/kHbarc/kHbarc)/(2.*TMath::Pi()*TMath::Pi()); 

  Double_t postFactor = 0.;
  //compute chemical potential
  Double_t potential = fBaryonPotential * particle->GetBaryonNumber() + 
                       fStrangePotential * particle->GetStrangeness() +
		       fElectroPotential * particle->GetElectricCharge();

  for(Int_t n = 1; n <= fNMax; ++n) {
    postFactor += TMath::Power(-d, n+1)/(n*n) *
                  TMath::Exp(n*potential/fTemperature) *
                  (3.*HankelKn(2, n*mass/fTemperature) + n*mass/fTemperature*HankelK1(n*mass/fTemperature));
  }
  return preFactor * postFactor;
}

//_______________________________________________________________________________
Double_t GrandCanonical::ParticleNumberDensity(ParticlePDG *const particle) {
  //
  // compute the particle number density
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::ParticleNumberDensity", "GrandCanonical object not fully initialized!!");

  Double_t degFactor = 2.*particle->GetSpin() + 1.;
  Double_t mass = particle->GetMass();     
  Double_t d = Int_t(2*particle->GetSpin()) & 1 ? 1. : -1.;
  Double_t preFactor = (degFactor*mass*mass*fTemperature/kHbarc/kHbarc/kHbarc)/(2.*TMath::Pi()*TMath::Pi());

  Double_t postFactor = 0.;
  Double_t potential = fBaryonPotential * particle->GetBaryonNumber() + 
                       fStrangePotential * particle->GetStrangeness() +
                       fElectroPotential * particle->GetElectricCharge();
                           
  for(Int_t n = 1; n <= fNMax; ++n) {
    postFactor += TMath::Power(-d, n+1)/n * 
                  TMath::Exp(n*potential/fTemperature) *
                  HankelKn(2, n*mass/fTemperature);        
  }
  return preFactor * postFactor;
}

//_______________________________________________________________________________
Double_t GrandCanonical::EnergyDensity(DatabasePDG *const database) {
  //
  // compute the total energy density
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::EnergyDensity", "GrandCanonical object not fully initialized!!");

  Double_t meanEnergyDensity = 0.;

  for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(currParticle);
    meanEnergyDensity += ParticleEnergyDensity(particle);
  }

  return meanEnergyDensity;
}

//_______________________________________________________________________________
Double_t GrandCanonical::BaryonDensity(DatabasePDG *const database) {
  //
  // compute the baryon density
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::BaryonDensity", "GrandCanonical object not fully initialized!!");

  Double_t meanBaryonDensity = 0.;

  for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(currParticle);
    meanBaryonDensity += ParticleNumberDensity(particle)*particle->GetBaryonNumber();
  }
  return meanBaryonDensity;
}

//_______________________________________________________________________________
Double_t GrandCanonical::StrangeDensity(DatabasePDG *const database) {
  //
  // compute the strangeness density
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::StrangeDensity", "GrandCanonical object not fully initialized!!");

  Double_t meanStrangeDensity = 0.;

  for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(currParticle);
    meanStrangeDensity += ParticleNumberDensity(particle)*particle->GetStrangeness();
  }

  return meanStrangeDensity;
}

//_______________________________________________________________________________
Double_t GrandCanonical::ElectroDensity(DatabasePDG *const database) {
  //
  // compute the electro number density
  //
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::ElectroDensity", "GrandCanonical object not fully initialized!!");

  Double_t meanElectroDensity = 0.;
  
  //hadrons
  for(Int_t currParticle = 0; currParticle<database->GetNParticles(); currParticle++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(currParticle);
    meanElectroDensity += ParticleNumberDensity(particle)*particle->GetElectricCharge();
  }

  return meanElectroDensity;
}
