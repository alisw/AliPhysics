//////////////////////////////////////////////////////////////////////////////////       
//                                                                              //
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna     //
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru  //
//                           November. 2, 2005                                  //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#ifndef GRANDCANONICAL_H
#define GRANDCANONICAL_H

class DatabasePDG;
class ParticlePDG;

class GrandCanonical {
 public:
  GrandCanonical();
  GrandCanonical(Int_t nmax, Double_t temperature, Double_t baryonPotential, Double_t strangePotential, Double_t electroPotential);
  ~GrandCanonical();

  void     Temperature(Double_t value); 
  Double_t Temperature() const { return fTemperature; }
  void     BaryonPotential(Double_t value);
  Double_t BaryonPotential() const { return fBaryonPotential; }
  void     StrangePotential(Double_t value);
  Double_t StrangePotential() const { return fStrangePotential; }
  void     ElectroPotential(Double_t value);
  Double_t ElectroPotential() const { return fElectroPotential; }
  void     NMax(Int_t value); 
  Int_t    NMax() const { return fNMax; }

  // compute of system baryon number, system strangeness, system charge and 
  // system energy
  // calculate system energy density
  Double_t EnergyDensity(DatabasePDG *const database);
  // calculate system baryon density
  Double_t BaryonDensity(DatabasePDG *const database);
  // calculate system strangeness density
  Double_t StrangeDensity(DatabasePDG *const database);
  // calculate system electro density
  Double_t ElectroDensity(DatabasePDG *const database);
  // compute of particle number density 
  Double_t ParticleNumberDensity(ParticlePDG *const particle);
  // compute the particle energy density 
  Double_t ParticleEnergyDensity(ParticlePDG *const particle); 

 private:

  Double_t    fTemperature;        // temperature
  Double_t    fBaryonPotential;	   // baryon chemical potential
  Double_t    fStrangePotential;   // strangeness chemical potential
  Double_t    fElectroPotential;   // electro chemical potential
  
  Int_t       fNMax;               //  Number of terms for summation, if fNMax = 1 then
  Bool_t fInitialized;             // flag
};

#endif
