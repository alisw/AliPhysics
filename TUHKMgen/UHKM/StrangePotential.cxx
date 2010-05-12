//                                                                          
//                                                                            
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
//                           November. 2, 2006                                
//
//

#include <TMath.h>
#include "StrangePotential.h"

//_____________________________________________________________________
Double_t StrangePotential::CalculateStrangePotential() {
  //
  // compute strange potential
  //
  Double_t minFunction = this->operator()(fMinStrangePotential);
  Double_t maxFunction = this->operator()(fMaxStrangePotential); 
  
  Int_t iter = 0;  
  while(minFunction < 0.0 && iter++ < fNIteration) {
    fMinStrangePotential -= 0.5*fMinStrangePotential;
    minFunction = this->operator()(fMinStrangePotential);
  }
   
  iter = 0;  
  while(minFunction*maxFunction > 0.0 && iter++ < fNIteration) {
    fMaxStrangePotential += 1.5*TMath::Abs(fMaxStrangePotential-fMinStrangePotential);
    maxFunction = this->operator()(fMaxStrangePotential);
  }
	
  if(minFunction*maxFunction > 0.0) {
    Error("StrangePotential::CalculateStrangePotential", "minFunction*maxFunction is positive!\n");
    return 0.;
  }

  EquationSolver<StrangePotential> * theSolver = 
    new EquationSolver<StrangePotential>(fNSolverIteration, fTolerance);

  theSolver->SetIntervalLimits(fMinStrangePotential, fMaxStrangePotential);
  
  if (!theSolver->Brent(*this))
    Error("StrangePotential::CalculateStrangePotential", "the root is not found!\n");
  
  Double_t strangePotential = theSolver->GetRoot();
  delete theSolver;
  return strangePotential;
}

//_____________________________________________________________________
Double_t StrangePotential::CalculateStrangeDensity(const Double_t strangePotential)
{
  //
  //calculate hadron system strange density
  //
  fGc.SetStrangePotential(strangePotential);
  fGc.SetTemperature(fTemperature);
  fGc.SetBaryonPotential(fBaryonPotential);
  return fGc.StrangenessDensity(fDatabase);
}
