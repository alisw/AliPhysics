/*                                                                        
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#ifndef GRANDCANONICAL_INCLUDED
#define GRANDCANONICAL_INCLUDED

#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif
#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif

class GrandCanonical {

 private:

  Double_t    fTemperature;     
  Double_t    fBaryonPotential;	
  Double_t    fStrangePotential;
  Double_t    fElectroPotential;

  //  Number of terms for summation, if fNMax = 1 then 
  //  Maxwell-Boltzmann distribution will be recovered
  Int_t       fNMax;
  Bool_t fInitialized;

 public:
  GrandCanonical();
  GrandCanonical(Int_t nmax, Double_t temperature, Double_t baryonPotential, Double_t strangePotential, Double_t electroPotential);
  ~GrandCanonical();

  void     Temperature(Double_t value); 
  Double_t Temperature() { return fTemperature; }
  void     BaryonPotential(Double_t value);
  Double_t BaryonPotential() { return fBaryonPotential; }
  void     StrangePotential(Double_t value);
  Double_t StrangePotential() { return fStrangePotential; }
  void     ElectroPotential(Double_t value);
  Double_t ElectroPotential() { return fElectroPotential; }
  void     NMax(Int_t value); 
  Int_t    NMax() { return fNMax; }

  // compute of system baryon number, system strangeness, system charge and 
  // system energy
  // calculate system energy density
  Double_t EnergyDensity(DatabasePDG* database);
  // calculate system baryon density
  Double_t BaryonDensity(DatabasePDG* database);
  // calculate system strangeness density
  Double_t StrangeDensity(DatabasePDG* database);
  // calculate system electro density
  Double_t ElectroDensity(DatabasePDG* database);
  // compute of particle number density 
  Double_t ParticleNumberDensity(ParticlePDG* particle);
  // compute the particle energy density 
  Double_t ParticleEnergyDensity(ParticlePDG* particle); 
};

#endif
