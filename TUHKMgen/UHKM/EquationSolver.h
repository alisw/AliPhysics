//////////////////////////////////////////////////////////////////////////////////       
//                                                                              //
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna     //
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru  //
//                           November. 2, 2005                                  //
//                                                                              //
//           This equation solver class is taken from GEANT4 and modified!!     //
//////////////////////////////////////////////////////////////////////////////////
#ifndef EQUATIONSOLVER_H
#define EQUATIONSOLVER_H
#include <Rtypes.h>

#define DefaultTolerance 5.0e-14

template <class Function> 
class EquationSolver {

 public:
  enum {kDefaultMaxIter = 100};  
  // default constructor
  EquationSolver() : fMaxIter(kDefaultMaxIter), fTolerance(DefaultTolerance),
    fA(0.0), fB(0.0), fRoot(0.0) {};
	
  EquationSolver(const Int_t iterations, const Double_t tol) :
    fMaxIter(iterations), fTolerance(tol),
    fA(0.0), fB(0.0), fRoot(0.0) {};

  // copy constructor	
  EquationSolver(const EquationSolver & right);

  // destructor
  ~EquationSolver() {};
	
  // operators
  EquationSolver & operator=(const EquationSolver & right);
  Bool_t operator==(const EquationSolver & right) const;
  Bool_t operator!=(const EquationSolver & right) const;
		
  Int_t GetMaxIterations(void) const {return fMaxIter;}
  void SetMaxIterations(const Int_t iterations) {fMaxIter=iterations;}
	
  Double_t GetTolerance(void) const {return fTolerance;}
  void SetTolerance(const Double_t epsilon) {fTolerance = epsilon;}
  
  Double_t GetIntervalLowerLimit(void) const {return fA;}
  Double_t GetIntervalUpperLimit(void) const {return fB;}
	
  void SetIntervalLimits(const Double_t Limit1, const Double_t Limit2);

  Double_t GetRoot(void) const {return fRoot;}	
	
  // Calculates the root by the Brent's method
  Bool_t Brent(Function& theFunction);

 private:
  Int_t fMaxIter;         // Maximum number of iterations
  Double_t fTolerance;    // tolerance (precision)
  Double_t fA;            // interval limits [a,b] which should bracket the root
  Double_t fB;            // interval limits [a,b] which should bracket the root
  Double_t fRoot;         // the equation's root
};

#include "EquationSolver.icc"

#endif
