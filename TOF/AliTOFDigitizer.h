#ifndef ALITOFDIGITIZER_H
#define ALITOFDIGITIZER_H

/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_______________________________________________________________________________//
//                                                                               //
//  Task Class for making Digits in TOF                                          //
// Class performs digitization of Summable digits (in the TOF case this is just  //
// sum of contributions of all signals into a given pad).                        //
// In addition it performs mixing of summable digits from different events.      //
//                                                                               //
// -- Author: Fabrizio Pierella (Bologna University)
//                                                                               //
//_______________________________________________________________________________//

/* $Id$ */

#include "AliDigitizer.h"

class AliDigitizationInput;
class AliTOFHitMap;
class AliTOFSDigit;
class AliTOFcalib;

class AliTOFDigitizer : public AliDigitizer {
 public:
  
  AliTOFDigitizer();
  AliTOFDigitizer(AliDigitizationInput * digInput);
  virtual ~AliTOFDigitizer();
  AliTOFDigitizer(const AliTOFDigitizer &source); // copy constructor
  AliTOFDigitizer& operator=(const AliTOFDigitizer &source); // ass. op.
  
  // Do the main work
  void Digitize(Option_t* option=0) ;
  TClonesArray* SDigits() const {return fSDigitsArray;}
  void ReadSDigit(Int_t inputFile);
  void CreateDigits();
  void InitDecalibration() const;
  void DecalibrateTOFSignal();
  
 private:
  void CollectSDigit(const AliTOFSDigit * const sdigit) ;
  Int_t PutNoise(Int_t /*charge*/)const {return 0;}; // not yet
						     // implemented
						     // due to the low
						     // noise expected
						     // level
  TClonesArray *fDigits;       //! array with digits
  TClonesArray *fSDigitsArray; //! List of summable digits; used as a
			       //container for all sdigits to be
			       //merged
  AliTOFHitMap *fhitMap ;      //! hit map used to perform the merging
  AliTOFcalib * fCalib;        //! calibration object
  
  ClassDef(AliTOFDigitizer,2)  // TOF/Merging/Digitization
};    
#endif
