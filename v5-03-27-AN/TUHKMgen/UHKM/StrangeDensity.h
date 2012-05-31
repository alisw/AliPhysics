//
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
//                           November. 2, 2005
//
//
//This class is used to obtain grand canonical description  of strange density
//by means of the temperature and chemical potentials (input). As for output
//we get  strange density.

#ifndef STRANGEDENSITY_H
#define STRANGEDENSITY_H

#include "HankelFunction.h"

class DatabasePDG;
class ParticlePDG;

class StrangeDensity {
 public:
  StrangeDensity();
  ~StrangeDensity(){};

  //for input
  void SetTemperature(Double_t value) {fTemperature = value;}
  void SetBaryonPotential(Double_t value) {fBaryonPotential = value;}
  void SetStrangePotential(Double_t value) {fStrangePotential = value;}
  void SetNMax(Int_t value) {
    fNMax = value; 
    if(fNMax < 1) fNMax = 1;
  }
  // compute hadron system strangeness density
  Double_t StrangenessDensity(const DatabasePDG* database);

 private:
  //input
  Double_t fTemperature;             // temperature
  Double_t fBaryonPotential;	     // baryon potential
  Double_t fStrangePotential;        // strange potential
  Int_t fNMax;   //number of terms for summation, if nMax = 1 then
                //Maxwell-Boltzmann distribution will be recovered	

  Double_t ParticleNumberDensity(ParticlePDG* particle);
};

#endif
