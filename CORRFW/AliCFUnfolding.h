
#ifndef ALICFUNFOLDING_H
#define ALICFUNFOLDING_H

//--------------------------------------------------------------------//
//                                                                    //
// AliCFUnfolding Class                                               //
// Class to handle general unfolding procedure using bayesian method  //
//                                                                    //
// Author : renaud.vernet@cern.ch                                     //
//--------------------------------------------------------------------//

#include "TNamed.h"
#include "THnSparse.h"
#include "AliLog.h"

class TF1;
class TRandom3;

class AliCFUnfolding : public TNamed {

 public :

  AliCFUnfolding();
  AliCFUnfolding(const Char_t* name, const Char_t* title, const Int_t nVar, 
		 const THnSparse* response, const THnSparse* efficiency, const THnSparse* measured, const THnSparse* prior=0x0, 
		 Double_t maxConvergencePerDOF = 1.e-06, UInt_t randomSeed = 0,
		 Int_t maxNumIterations = 10);
  ~AliCFUnfolding();
  void UnsetCorrelatedErrors()  {AliError("===================> DEPRECATED <=====================");}
  void SetUseCorrelatedErrors() {AliError("===================> DEPRECATED <=====================");}
  void SetMaxNumberOfIterations(Int_t n = 10)  {
    AliError("===================> DEPRECATED : should be set in constructor <=====================");
    AliError("===================> DEPRECATED : will be removed soon         <=====================");
    fMaxNumIterations = n;
  }

  void SetNRandomIterations(Int_t n = 100) {fNRandomIterations = n;};

  void UseSmoothing(TF1* fcn=0x0, Option_t* opt="iremn") { // if fcn=0x0 then smooth using neighbouring bins 
    fUseSmoothing=kTRUE;                                   // this function must NOT be used if fNVariables > 3
    fSmoothFunction=fcn;                                   // the option "opt" is used if "fcn" is specified
    fSmoothOption=opt;
  } 
                                                                                                
  void Unfold();

  const THnSparse* GetResponse()             const {return fResponseOrig;}
  const THnSparse* GetEfficiency()           const {return fEfficiencyOrig;}
  const THnSparse* GetMeasured()             const {return fMeasuredOrig;}
  const THnSparse* GetOriginalPrior()        const {return fPriorOrig;}
        THnSparse* GetInverseResponse()      const {return fInverseResponse;}
        THnSparse* GetPrior()                const {return fPrior;}
	THnSparse* GetUnfolded()             const {return fUnfoldedFinal;}
        THnSparse* GetEstMeasured()          const {return fMeasuredEstimate;}
	THnSparse* GetConditional()          const {return fConditional;}
	TF1*       GetSmoothFunction()       const {return fSmoothFunction;}
	THnSparse* GetDeltaUnfoldedProfile() const {return fDeltaUnfoldedP;}
	Int_t      GetDOF();                 // Returns number of degrees of freedom

  static Short_t  SmoothUsingNeighbours(THnSparse*); // smoothes the unfolded spectrum using the neighbouring cells

 private :
  AliCFUnfolding(const AliCFUnfolding& c);
  AliCFUnfolding& operator= (const AliCFUnfolding& c);
  
  //
  // user-related settings
  //
  const THnSparse     *fResponseOrig;     // Response matrix : dimensions must be 2N = 2 x (number of variables)
                                          // dimensions 0 ->  N-1 must be filled with reconstructed values
                                          // dimensions N -> 2N-1 must be filled with generated values
  const THnSparse     *fPriorOrig;        // This is the assumed generated distribution : dimensions must be N = number of variables
                                          // it will be used at the first step 
                                          // then will be updated automatically at each iteration
  const THnSparse     *fEfficiencyOrig;   // Efficiency map : dimensions must be N = number of variables (modified)
                                          // this map must be filled only with "true" values of the variables (should not do "bin smearing")
  const THnSparse     *fMeasuredOrig;     // Measured spectrum to be unfolded : dimensions must be N = number of variables (modified)

        Int_t          fMaxNumIterations; // Maximum  number of iterations to be performed
	Int_t          fNVariables;       // Number of variables used in analysis spectra (pt, y, ...)
        Bool_t         fUseSmoothing;     // Smooth the unfolded sectrum at each iteration; default is kFALSE
	TF1           *fSmoothFunction;   // Function used to smooth the unfolded spectrum
	Option_t      *fSmoothOption;     // Option to use during the fit (with fSmoothFunction) ; default is "iremn"

  //
  // internal settings
  //
	Double_t       fMaxConvergence;    // Convergence criterion in case of correlated error calculation
        Int_t          fNRandomIterations; // Number of random distributed measured spectra
	THnSparse     *fResponse;          // Copy of the original response matrix    (modified)
	THnSparse     *fPrior;             // Copy of the original prior spectrum     (modified)
	THnSparse     *fEfficiency;        // Copy of original efficiency             (modified)
	THnSparse     *fMeasured;          // Copy of the original measureed spectrum (modified)
        THnSparse     *fInverseResponse;   // Inverse response matrix
        THnSparse     *fMeasuredEstimate;  // Estimation of the measured (M) spectrum given the a priori (T) distribution
	THnSparse     *fConditional;       // Matrix holding the conditional probabilities P(M|T)
        THnSparse     *fUnfolded;          // Unfolded spectrum (modified before and during error calculation)
	THnSparse     *fUnfoldedFinal;     // Final unfolded spectrum
        Int_t         *fCoordinates2N;     // Coordinates in 2N (measured,true) space
	Int_t         *fCoordinatesN_M;    // Coordinates in measured space
	Int_t         *fCoordinatesN_T;    // Coordinates in true space


  /* correlated error calculation */
  THnSparse     *fRandomResponse;    // Randomized distribution for each bin of the response matrix     to calculate correlated errors
  THnSparse     *fRandomEfficiency;  // Randomized distribution for each bin of the efficiency spectrum to calculate correlated errors
  THnSparse     *fRandomMeasured;    // Randomized distribution for each bin of the measured   spectrum to calculate correlated errors
  TRandom3      *fRandom3;           // Object to get random number following Poisson distribution
  THnSparse     *fDeltaUnfoldedP;    // Profile of the delta-unfolded distribution
  THnSparse     *fDeltaUnfoldedN;    // Entries of the delta-unfolded distribution (count for each bin)
  Short_t        fNCalcCorrErrors;   // Book-keeping to prevend infinite loop
  UInt_t         fRandomSeed;        // Random seed


  // functions
  void     Init();                  // initialisation of the internal settings
  void     GetCoordinates();        // gets a cell coordinates in Measured and True space
  void     CreateConditional();     // creates the conditional matrix from the response matrix
  void     CreateEstMeasured();     // creates the measured spectrum estimation from the conditional matrix and the prior distribution
  void     CreateInvResponse();     // creates the inverse response function (Bayes Theorem) from the conditional matrix and the prior distribution
  void     CreateUnfolded();        // creates the unfolded spectrum from the inverse response matrix and the measured distribution
  void     CreateFlatPrior();       // creates a flat a priori distribution in case the one given in the constructor is null
  Double_t GetChi2();               // returns the chi2 between unfolded and prior spectra
  Short_t  Smooth();                // function calling smoothing methods
  Short_t  SmoothUsingFunction();   // smoothes the unfolded spectrum using a fit function

  /* correlated error calculation */
  Double_t GetConvergence();            // Returns convergence criterion
  void     CalculateCorrelatedErrors(); // Calculates correlated errors for the final unfolded spectrum
  void     CreateRandomizedDist();      // Create randomized dist from measured distribution
  void     FillDeltaUnfoldedProfile();  // Fills the fDeltaUnfoldedP profile
  void     SetMaxConvergencePerDOF (Double_t val);

  ClassDef(AliCFUnfolding,1);
};

#endif
