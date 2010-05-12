//                                                                            
//                                                                            
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
//                           November. 2, 2005                                
//
//
//      This class is taken from the GEANT4 tool kit  and changed!!!!!

#include <TError.h>
#include "RandArrayFunction.h"


RandArrayFunction::RandArrayFunction(const Double_t *aProbFunc, Int_t theProbSize, Int_t intType):
  fIntegralPdf(),
  fNBins(theProbSize),
  fOneOverNbins(0),
  fInterpolationType(intType)
{
  PrepareTable(aProbFunc);
}

RandArrayFunction::RandArrayFunction(Int_t theProbSize, Int_t intType):
  fIntegralPdf(),
  fNBins(theProbSize), 
  fOneOverNbins(0),
  fInterpolationType(intType)
{}

void RandArrayFunction::PrepareTable(const Double_t* aProbFunc) {
  //Prepares fIntegralPdf.
  if(fNBins < 1) {
    Error("RandArrayFunction::PrepareTable",
	  "RandArrayFunction constructed with no bins - will use flat distribution.");
    UseFlatDistribution();
    return;
  }

  fIntegralPdf.resize(fNBins + 1);
  fIntegralPdf[0] = 0;
  Int_t ptn;
  for (ptn = 0; ptn < fNBins; ++ptn ) {
    Double_t weight = aProbFunc[ptn];
    if (weight < 0.) {
      // We can't stomach negative bin contents, they invalidate the 
      // search algorithm when the distribution is fired.
      Warning("RandArrayFunction::PrepareTable",
	      "RandArrayFunction constructed with negative-weight bin %d == %f -- will substitute 0 weight",
	      ptn, weight);
      weight = 0.;
    }
    fIntegralPdf[ptn + 1] = fIntegralPdf[ptn] + weight;
  }

  if (fIntegralPdf[fNBins] <= 0.) {
    Warning("RandArrayFunction::PrepareTable",
	    "RandArrayFunction constructed with nothing in bins - will use flat distribution");
    UseFlatDistribution();
    return;
  }

  for (ptn = 0; ptn < fNBins + 1; ++ptn)
    fIntegralPdf[ptn] /= fIntegralPdf[fNBins];
  
  // And another useful variable is ...
  fOneOverNbins = 1.0 / fNBins;
  // One last chore:
  if (fInterpolationType && fInterpolationType != 1) {
    Info("RandArrayFunction::PrepareTable",
	 "RandArrayFunction does not recognize fInterpolationType %d \n"
	 "Will use type 0 (continuous linear interpolation)", fInterpolationType);
    fInterpolationType = 0;
  }
} 

void RandArrayFunction::UseFlatDistribution() {
  //Called only by PrepareTable in case of user error. 
  fNBins = 1;
  fIntegralPdf.resize(2);
  fIntegralPdf[0] = 0;
  fIntegralPdf[1] = 1;
  fOneOverNbins = 1.0;
} 

Double_t RandArrayFunction::MapRandom(Double_t rand) const {
  // Private method to take the random (however it is created) and map it
  // according to the distribution.

  Int_t nBelow = 0;	  // largest k such that I[k] is known to be <= rand
  Int_t nAbove = fNBins;  // largest k such that I[k] is known to be >  rand
  Int_t middle;
     
  while (nAbove > nBelow+1) {
    middle = (nAbove + nBelow+1)>>1;
    rand >= fIntegralPdf[middle] ? nBelow = middle : nAbove = middle;
  }// after this loop, nAbove is always nBelow+1 and they straddle rad:
       
  /*assert ( nAbove = nBelow+1 );
    assert ( fIntegralPdf[nBelow] <= rand );
    assert ( fIntegralPdf[nAbove] >= rand );*/  
  // If a defective engine produces rand=1, that will 
  // still give sensible results so we relax the > rand assertion

  if (fInterpolationType == 1) {
    return nBelow * fOneOverNbins;
  } 
  else {
    Double_t binMeasure = fIntegralPdf[nAbove] - fIntegralPdf[nBelow];
    // binMeasure is always aProbFunc[nBelow], 
    // but we don't have aProbFunc any more so we subtract.
    
    if (!binMeasure) { 
      // rand lies right in a bin of measure 0.  Simply return the center
      // of the range of that bin.  (Any value between k/N and (k+1)/N is 
      // equally good, in this rare case.)
      return (nBelow + .5) * fOneOverNbins;
    }

    Double_t binFraction = (rand - fIntegralPdf[nBelow]) / binMeasure;
    
    return (nBelow + binFraction) * fOneOverNbins;
  }
} 

void RandArrayFunction::FireArray(Int_t size, Double_t *vect) const {
  for (Int_t i = 0; i < size; ++i)
    vect[i] = Fire();
}
