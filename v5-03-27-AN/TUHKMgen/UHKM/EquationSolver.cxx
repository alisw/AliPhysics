//////////////////////////////////////////////////////////////////////////////////       
//                                                                              //
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna     //
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru  //
//                           November. 2, 2005                                  //
//                                                                              //
//           This equation solver class is taken from GEANT4 and modified!!     //
//////////////////////////////////////////////////////////////////////////////////

#include "EquationSolver.h"

template <class Function> 
EquationSolver<Function>::EquationSolver(const EquationSolver & right) {
  fMaxIter = right.GetMaxIterations();
  fTolerance = right.GetTolerance();
  fA = right.GetIntervalLowerLimit();
  fB = right.GetIntervalUpperLimit();
  fRoot = right.GetRoot();
}

// operators
template <class Function> 
EquationSolver<Function> & EquationSolver<Function>::operator=(const EquationSolver & right) {
  fMaxIter = right.GetMaxIterations();
  fTolerance = right.GetTolerance();
  fA = right.GetIntervalLowerLimit();
  fB = right.GetIntervalUpperLimit();
  fRoot = right.GetRoot();
  return *this;
}

template <class Function> 
Bool_t EquationSolver<Function>::operator==(const EquationSolver & right) const {
  if (this == &right) return true;
  else return false;
}

template <class Function> 
Bool_t EquationSolver<Function>::operator!=(const EquationSolver & right) const {
  return !operator==(right);
}

