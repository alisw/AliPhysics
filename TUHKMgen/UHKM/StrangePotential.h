//
//
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
//                           November. 2, 2005
//
//
//This class is used to calculate strange potential from
//the known initial strange density = 0 at given temperature and baryon potential.

#ifndef STRANGEPOTENTIAL_H
#define STRANGEPOTENTIAL_H

#include "StrangeDensity.h"
#include "EquationSolver.h"
#include "DatabasePDG.h"
#include "UKUtility.h"
                                         
class StrangePotential {
 public:
  StrangePotential(const Double_t initialStrangeDensity=0, DatabasePDG* database=0x0) :
    fTemperature(0),
    fBaryonPotential(0),
    fStrangeDensity(initialStrangeDensity),
    fMinStrangePotential(0.0001*kGeV),
    fMaxStrangePotential(0.9*kGeV),
    fNIteration(100),
    fNSolverIteration(100),
    fTolerance(1.e-8),
    fDatabase(database),
    fGc()
    {};

  ~StrangePotential() {};
   
  Double_t operator()(const Double_t strangePotential) { 
    return (fStrangeDensity - this->CalculateStrangeDensity(strangePotential))/fStrangeDensity; 
  }	

  void SetTemperature(Double_t value) {fTemperature = value;}
  void SetBaryonPotential(Double_t value) {fBaryonPotential = value;}
  void SetMinStrangePotential(Double_t value) {fMinStrangePotential = value;}
  void SetMaxStrangePotential(Double_t value) {fMaxStrangePotential = value;}
  Double_t CalculateStrangePotential();

 private:
  StrangePotential();
  StrangePotential(const StrangePotential&);
  StrangePotential& operator=(const StrangePotential&);

  Double_t fTemperature;         // temperature
  Double_t fBaryonPotential;     // baryo-chemical potential
  Double_t fStrangeDensity;      // strangeness density
  Double_t fMinStrangePotential;//initial min value of strange potential 
  Double_t fMaxStrangePotential;//initial max value of strange potential
  Int_t fNIteration; //to find proper [minStrangePotential, maxStrangePotential] interval
  Int_t fNSolverIteration; //to find root in [minStrangePotential,maxStrangePotential] interval
  Double_t fTolerance;//to find root 
  DatabasePDG* fDatabase;        // PDG database
  StrangeDensity fGc;            // strangeness density object
  //compute hadron  system strange density through strange potential
  Double_t CalculateStrangeDensity(const Double_t strangePotential);
};

#endif
