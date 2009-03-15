#ifndef NAStrangePotential_h
#define NAStrangePotential_h
#include "StrangeDensity.h"
#include "EquationSolver.h"
/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

//This class is used to calculate strange potential from 
//the known initial strange density = 0 at given temperature and baryon potential.

class NAStrangePotential {
 public:
  NAStrangePotential(const Double_t initialStrangeDensity) :
    fStrangeDensity(initialStrangeDensity),
    fNIteration(100),
    fTolerance(1.e-8),
    fNSolverIteration(100),
    fMinStrangePotential(0.0001*GeV),
    fMaxStrangePotential(0.9*GeV)
    {};

  ~NAStrangePotential() {};
  
  Double_t operator()(const Double_t strangePotential) { 
    return (fStrangeDensity - this->CalculateStrangeDensity(strangePotential))/fStrangeDensity; 
  }	

 private:
  //  default constructor is not accesible
  NAStrangePotential(){};

 public:
  void SetTemperature(Double_t value) {fTemperature = value;}
  void SetBaryonPotential(Double_t value) {fBaryonPotential = value;}
  void SetMinStrangePotential(Double_t value) {fMinStrangePotential = value;}
  void SetMaxStrangePotential(Double_t value) {fMaxStrangePotential = value;}
  Double_t CalculateStrangePotential();

 private:
  //compute hadron  system strange density through strange potential
  Double_t CalculateStrangeDensity(const Double_t strangePotential);
  Double_t fTemperature;
  Double_t fBaryonPotential;
  Double_t fStrangeDensity;
  Double_t fMinStrangePotential;//initial min value of strange potential 
  Double_t fMaxStrangePotential;//initial max value of strange potential
  Int_t fNIteration; //to find proper [minStrangePotential, maxStrangePotential] interval
  Int_t fNSolverIteration; //to find root in [minStrangePotential,maxStrangePotential] interval
  Double_t fTolerance;//to find root 
  NAStrangeDensity fGc;
};

#endif
