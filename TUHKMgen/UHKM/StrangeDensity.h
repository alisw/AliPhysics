/*

        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
                           November. 2, 2005

*/

//This class is used to obtain grand canonical description  of strange density
//by means of the temperature and chemical potentials (input). As for output
//we get  strange density.

#ifndef NAStrangeDensity_h
#define NAStrangeDensity_h 1

#ifndef NAMathUtil_h
#include "MathUtil.h"
#endif
#ifndef HANKELFUNCTION_INCLUDED
#include "HankelFunction.h"
#endif
#ifndef PARTICLE_INCLUDED
#include "Particle.h"
#endif
#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif

class NAStrangeDensity {

 private:
  //input
  Double_t fTemperature;
  Double_t fBaryonPotential;	
  Double_t fStrangePotential;
  Int_t fNMax;   //number of terms for summation, if nMax = 1 then
                //Maxwell-Boltzmann distribution will be recovered	

  Double_t ParticleNumberDensity(ParticlePDG* particle);

 public:
  NAStrangeDensity();
  ~NAStrangeDensity(){};

  //for input
  void SetTemperature(Double_t value) {fTemperature = value;}
  void SetBaryonPotential(Double_t value) {fBaryonPotential = value;}
  void SetStrangePotential(Double_t value) {fStrangePotential = value;}
  void SetNMax(Int_t value) {
    fNMax = value; 
    if(fNMax < 1) fNMax = 1;
  }
  // compute hadron system strangeness density
  Double_t StrangenessDensity(DatabasePDG* database);
};

#endif
