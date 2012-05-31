/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/
//This class is taken from the GEANT4 tool kit and changed!!!!!

//========================================================================================
//RandArrayFunction defines several methods for shooting generally distributed random values, 
//given a user-defined probability distribution function.

//The probability distribution function Pdf must be provided by the user as an array of 
//positive real numbers. The array size must also be provided. Pdf doesn't need to be 
//normalized to 1.

  // if IntType = 0 ( default value ) a uniform random number is
  // generated using the StandardRand() engine. The uniform number is then transformed
  // to the user's distribution using the cumulative probability
  // distribution constructed from his histogram. The cumulative
  // distribution is inverted using a binary search for the nearest
  // bin boundary and a linear interpolation within the
  // bin. RandArrayFunction therefore generates a constant density within
  // each bin.
  // if IntType = 1 no interpolation is performed and the result is a
  // discrete distribution.

  //A speculate set of Shoot()/ShootArray() and Fire()/FireArray() methods is provided 
  //to Shoot random numbers via an instantiated RandArrayFunction object. These methods 
  //act directly on the flat distribution provided by a StandardRand() engine. 
  //An Operator () is also provided. 

//  example.
//      ...
//      Double_t* Pdf;
//      Int_t fNBins;
//      ...
//      RandArrayFunction FunctDist(Pdf,fNBins);
//      ... 
//      Double_t num = FunctDist.Shoot();//Shoot() provides the same functionality as Fire()

//  example.
//      ...
//      Double_t* Pdf;
//      Int_t fNBins;
//      ...
//      RandArrayFunction FunctDist(Pdf,fNBins);
//      ... 
//      Double_t num = FunctDist(); 

//  example.
//      ...
//      Double_t* Pdf;
//      Int_t fNBins;
//      ...
//      RandArrayFunction FunctDist(Pdf,fNBins);
//      ...
//	    Int_t size = 50;
//	    Double_t* vect = new Double_t[size];
//      FunctDist.FireArray (size, vect);

//========================================================================================

#ifndef RANDARRAYFUNCTION_H
#define RANDARRAYFUNCTION_H

#include <vector>
#include <TRandom.h>

class RandArrayFunction {
 public:
  RandArrayFunction(const Double_t *aProbFunc, Int_t theProbSize, Int_t interpolationType = 0);
  RandArrayFunction(Int_t probSize, Int_t interpolationType = 0);

  Double_t Shoot()const;
  Double_t Fire()const;
  Double_t operator()()const;
  void     ShootArray(Int_t size, Double_t *array)const;
  void     FireArray(Int_t size, Double_t *array)const;

  void     PrepareTable(const Double_t *aProbFunc);

 private:
  void     UseFlatDistribution();
  Double_t MapRandom(Double_t rand)const;
  Double_t StandardRand()const;

  std::vector<Double_t> fIntegralPdf;         //
  Int_t                 fNBins;               //
  Double_t              fOneOverNbins;        //
  Int_t                 fInterpolationType;   //

};

inline Double_t RandArrayFunction::StandardRand() const {
  return gRandom->Rndm();
}

inline Double_t RandArrayFunction::Fire() const {
  return MapRandom(StandardRand());
}

inline Double_t RandArrayFunction::Shoot() const {
  return Fire();
}

inline Double_t RandArrayFunction::operator()() const {
  return Fire();
}

inline void RandArrayFunction::ShootArray(Int_t size, Double_t *array) const {
  FireArray(size, array);
}

#endif
