#ifndef NAStrangeDensity_h
#include "StrangeDensity.h"
#endif

NAStrangeDensity::NAStrangeDensity() {
  fTemperature = 0.*GeV;
  fBaryonPotential = 0.*GeV;
  fStrangePotential = 0.*GeV;
  fNMax = 5;
}
// compute hadron system strangeness density
Double_t NAStrangeDensity::StrangenessDensity(DatabasePDG* database) {
  Double_t meanStrangenessDensity = 0.;
  for(Int_t particleIndex = 0; particleIndex < database->GetNParticles(); particleIndex++) {
    ParticlePDG *particle = database->GetPDGParticleByIndex(particleIndex);
    Double_t particleDensity = ParticleNumberDensity(particle);
    meanStrangenessDensity += particleDensity*particle->GetStrangeness();
  }
  return meanStrangenessDensity;
}

// compute hadron number density
Double_t NAStrangeDensity::ParticleNumberDensity(ParticlePDG* pDef) {
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
	       fTemperature/hbarc/hbarc/hbarc)/(2.*N_PI*N_PI);  
  postfactor = 0.;
 
  for(Int_t n = 1; n <= fNMax; n++) {
    postfactor += pow(-d,n+1)/(n)*exp(n*particleChemPotential/fTemperature)*
      HankelKn(2,n*particleMass/fTemperature);
  }
  return prefactor*postfactor;
}

