
#ifndef NAEquationSolver_h
#define NAEquationSolver_h 1
#include <Rtypes.h>
#include "MathUtil.h"
/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/
//This equation solver class is taken from GEANT4 and modified!!


#define DefaultTolerance 5.0e-14

template <class Function> 
class NAEquationSolver {

 public:
  enum {DefaultMaxIter = 100};
	
 private:
  // Maximum number of iterations
  Int_t fMaxIter;
  Double_t fTolerance;
  // interval limits [a,b] which should bracket the root
  Double_t fA;
  Double_t fB;
  // root
  Double_t fRoot;

 public:    
  // default constructor
  NAEquationSolver() : fMaxIter(DefaultMaxIter), fTolerance(DefaultTolerance),
    fA(0.0), fB(0.0), fRoot(0.0) {};
	
  NAEquationSolver(const Int_t iterations, const Double_t tol) :
    fMaxIter(iterations), fTolerance(tol),
    fA(0.0), fB(0.0), fRoot(0.0) {};

  // copy constructor	
  NAEquationSolver(const NAEquationSolver & right);

  // destructor
  ~NAEquationSolver() {};
	
  // operators
  NAEquationSolver & operator=(const NAEquationSolver & right);
  Bool_t operator==(const NAEquationSolver & right) const;
  Bool_t operator!=(const NAEquationSolver & right) const;
		
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
};

#include "EquationSolver.icc"

#endif
