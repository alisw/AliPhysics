/*

        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
                           November. 2, 2005

*/

//This class is used to calculate strange potential from
//the known initial strange density = 0 at given temperature and baryon potential.

#ifndef NAStrangePotential_h
#define NAStrangePotential_h 1

#ifndef NAStrangeDensity_h
#include "StrangeDensity.h"
#endif
#ifndef NAEquationSolver_h
#include "EquationSolver.h"
#endif
#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif
                                         
class NAStrangePotential {
 public:
  NAStrangePotential(const Double_t initialStrangeDensity=0, DatabasePDG* database=0x0) :
    fTemperature(0),
    fBaryonPotential(0),
    fStrangeDensity(initialStrangeDensity),
    fMinStrangePotential(0.0001*GeV),
    fMaxStrangePotential(0.9*GeV),
    fNIteration(100),
    fNSolverIteration(100),
    fTolerance(1.e-8),
    fDatabase(database),
    fGc()
    {};

  ~NAStrangePotential() {};
   
  Double_t operator()(const Double_t strangePotential) { 
    return (fStrangeDensity - this->CalculateStrangeDensity(strangePotential))/fStrangeDensity; 
  }	

  void SetTemperature(Double_t value) {fTemperature = value;}
  void SetBaryonPotential(Double_t value) {fBaryonPotential = value;}
  void SetMinStrangePotential(Double_t value) {fMinStrangePotential = value;}
  void SetMaxStrangePotential(Double_t value) {fMaxStrangePotential = value;}
  Double_t CalculateStrangePotential();

 private:
  NAStrangePotential(const NAStrangePotential&);
  NAStrangePotential& operator=(const NAStrangePotential&);

  Double_t fTemperature;
  Double_t fBaryonPotential;
  Double_t fStrangeDensity;
  Double_t fMinStrangePotential;//initial min value of strange potential 
  Double_t fMaxStrangePotential;//initial max value of strange potential
  Int_t fNIteration; //to find proper [minStrangePotential, maxStrangePotential] interval
  Int_t fNSolverIteration; //to find root in [minStrangePotential,maxStrangePotential] interval
  Double_t fTolerance;//to find root 
  DatabasePDG* fDatabase;
  NAStrangeDensity fGc;
  //compute hadron  system strange density through strange potential
  Double_t CalculateStrangeDensity(const Double_t strangePotential);
  //default constructor is not accesible
  NAStrangePotential();

};

#endif
