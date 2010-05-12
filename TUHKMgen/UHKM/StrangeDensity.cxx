//
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
//                           November. 2, 2005
//
//
//This class is used to obtain grand canonical description  of strange density
//by means of the temperature and chemical potentials (input). As for output

#include <TMath.h>
#include "StrangeDensity.h"
#include "DatabasePDG.h"
#include "ParticlePDG.h"
#include "UKUtility.h"

//__________________________________________________________
StrangeDensity::StrangeDensity():
  fTemperature(0.*kGeV),
  fBaryonPotential(0.*kGeV),
  fStrangePotential(0.*kGeV),
  fNMax(5)
{
  //
  // constructor
  //
}

//__________________________________________________________
Double_t StrangeDensity::StrangenessDensity(const DatabasePDG* database) {
  //
  // compute hadron system strangeness density
  //
  Double_t meanStrangenessDensity = 0.;
  for(Int_t particleIndex = 0; particleIndex < database->GetNParticles(); particleIndex++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(particleIndex);
    Double_t particleDensity = ParticleNumberDensity(particle);
    meanStrangenessDensity += particleDensity*particle->GetStrangeness();
  }
  return meanStrangenessDensity;
}

//__________________________________________________________
Double_t StrangeDensity::ParticleNumberDensity(ParticlePDG* pDef) {
  //
  // compute hadron number density
  //
  Double_t particleMass = pDef->GetMass();
  Int_t particleStrangeness = Int_t(pDef->GetStrangeness());
  Double_t particleBaryon = pDef->GetBaryonNumber();
  //compute chemical potential
  Double_t particleChemPotential = fBaryonPotential*particleBaryon + 
	                           fStrangePotential*particleStrangeness;
  //compute degeneracy factor
  Double_t particleDegFactor = 2*pDef->GetSpin() + 1.;     // IA: In ParticlePDG() GetSpin() returns spin not 2*spin !!
  Double_t d = 1.;//for fermions
  if(Int_t(2*pDef->GetSpin())%2 == 0)//it gives 0 for Spin = 0,2,4,.. and it gives 1 for Spin = 1,3,7,
    d = -1; //for bosons

  Double_t prefactor;
  Double_t postfactor;
  prefactor = (particleDegFactor*particleMass*particleMass*
	       fTemperature/kHbarc/kHbarc/kHbarc)/(2.*TMath::Pi()*TMath::Pi());  
  postfactor = 0.;
 
  for(Int_t n = 1; n <= fNMax; n++) {
    postfactor += pow(-d,n+1)/(n)*exp(n*particleChemPotential/fTemperature)*
      HankelKn(2,n*particleMass/fTemperature);
  }
  return prefactor*postfactor;
}

