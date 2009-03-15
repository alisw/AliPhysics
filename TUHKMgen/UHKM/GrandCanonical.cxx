/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/
#include <iostream>

#include <TError.h>
#include <TMath.h>
#ifndef GRANDCANONICAL_INCLUDED
#include "GrandCanonical.h"
#endif
#ifndef HANKELFUNCTION_INCLUDED
#include "HankelFunction.h"
#endif
#ifndef UKUTILITY_INCLUDED
#include "UKUtility.h"
#endif

GrandCanonical::GrandCanonical() {
  fInitialized = kFALSE;
  fNMax = -1111;
  fTemperature = -1111;
  fBaryonPotential = -1111;
  fStrangePotential = -1111;
  fElectroPotential = -1111;
}

GrandCanonical::GrandCanonical(Int_t nmax, Double_t temperature, Double_t baryonPotential, Double_t strangePotential, Double_t electroPotential) {
  fNMax = nmax;
  fTemperature = temperature;
  fBaryonPotential = baryonPotential;
  fStrangePotential = strangePotential;
  fElectroPotential = electroPotential;
  fInitialized = kTRUE;
}

GrandCanonical::~GrandCanonical() {}


void GrandCanonical::Temperature(Double_t value) {
  fTemperature = value;
  if(fNMax!=-1111 && fBaryonPotential!=-1111 && fStrangePotential!=-1111 && fElectroPotential!=-1111)
    fInitialized = kTRUE;
}

void GrandCanonical::BaryonPotential(Double_t value) {
  fBaryonPotential = value;
  if(fNMax!=-1111 && fTemperature!=-1111 && fStrangePotential!=-1111 && fElectroPotential!=-1111)
    fInitialized = kTRUE;
}

void GrandCanonical::StrangePotential(Double_t value) {
  fStrangePotential = value;
  if(fNMax!=-1111 && fTemperature!=-1111 && fBaryonPotential!=-1111 && fElectroPotential!=-1111)
    fInitialized = kTRUE;
}

void GrandCanonical::ElectroPotential(Double_t value) {
  fElectroPotential = value;
  if(fNMax!=-1111 && fTemperature!=-1111 && fBaryonPotential!=-1111 && fStrangePotential!=-1111)
    fInitialized = kTRUE;
}

void GrandCanonical::NMax(Int_t value) {
  fNMax = value;
  if(fTemperature!=-1111 && fBaryonPotential!=-1111 && fStrangePotential!=-1111 && fElectroPotential!=-1111)
    fInitialized = kTRUE;
}

Double_t GrandCanonical::ParticleEnergyDensity(ParticlePDG* particle) {
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::ParticleEnergyDensity", "GrandCanonical object not fully initialized!!");
  
  // Compute the particle energy density
  Double_t degFactor = 2.*particle->GetSpin() + 1.;                                    // degeneracy factor
  Double_t mass = particle->GetMass();                                                // PDG table mass
  Double_t d = Int_t(2.*particle->GetSpin()) & 1 ? 1. : -1;                                   // Bose-Einstein/Fermi-Dirac factor
  Double_t preFactor = (degFactor*mass*mass*fTemperature*fTemperature/hbarc/hbarc/hbarc)/(2.*TMath::Pi()*TMath::Pi()); 

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

Double_t GrandCanonical::ParticleNumberDensity(ParticlePDG* particle) {
  // Check if all the thermodinamic parameters are set
  if(!fInitialized)
    Fatal("GrandCanonical::ParticleNumberDensity", "GrandCanonical object not fully initialized!!");

  Double_t degFactor = 2.*particle->GetSpin() + 1.;
  Double_t mass = particle->GetMass();     
  Double_t d = Int_t(2*particle->GetSpin()) & 1 ? 1. : -1.;
  Double_t preFactor = (degFactor*mass*mass*fTemperature/hbarc/hbarc/hbarc)/(2.*TMath::Pi()*TMath::Pi());

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


Double_t GrandCanonical::EnergyDensity(DatabasePDG* database) {
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

Double_t GrandCanonical::BaryonDensity(DatabasePDG* database) {
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

Double_t GrandCanonical::StrangeDensity(DatabasePDG* database) {
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

Double_t GrandCanonical::ElectroDensity(DatabasePDG* database) {
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
