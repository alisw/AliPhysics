#ifndef ALITOFDIGITIZER_H
#define ALITOFDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  Task Class for making Digits in TOF      
// Class performs digitization of Summable digits (in the TOF case this is just
// sum of contributions of all signals into a given pad). 
// In addition it performs mixing of summable digits from different events.
//                  
//*-- Author: Fabrizio Pierella (Bologna University)

#include "TH1F.h"
#include "TList.h"
#include "AliDigitizer.h"
#include "AliTOFcalib.h"
#include "AliTOFCal.h"
#include "AliTOFGeometry.h"

class AliRunDigitizer;
class AliTOFHitMap;
class AliTOFSDigit;

class AliTOFDigitizer : public AliDigitizer {
 public:
  
  AliTOFDigitizer();
  AliTOFDigitizer(AliRunDigitizer * manager);
  virtual ~AliTOFDigitizer();
  AliTOFDigitizer(const AliTOFDigitizer &source); // copy constructor
  AliTOFDigitizer& operator=(const AliTOFDigitizer &source); // ass. op.
  
  // Do the main work
  void Exec(Option_t* option=0) ;
  TClonesArray* SDigits() const {return fSDigitsArray;}
  void ReadSDigit(Int_t inputFile);
  void CreateDigits();
  void InitDecalibration(AliTOFcalib *calib) const;
  void DecalibrateTOFSignal(AliTOFcalib *calib);
  
 private:
  void CollectSDigit(AliTOFSDigit * sdigit) ;
  Int_t PutNoise(Int_t /*charge*/)const {return 0;}; // not yet
						     // implemented
						     // due to the low
						     // noise expected
						     // level
  AliTOFGeometry *fGeom;       // Pointer to the TOF geometry
  TClonesArray *fDigits;       //! array with digits
  TClonesArray *fSDigitsArray; //! List of summable digits; used as a
			       //container for all sdigits to be
			       //merged
  AliTOFHitMap *fhitMap ;      //! hit map used to perform the merging
  
  ClassDef(AliTOFDigitizer,0)  // TOF/Merging/Digitization
};    
#endif
