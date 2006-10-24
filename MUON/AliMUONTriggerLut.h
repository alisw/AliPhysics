#ifndef ALIMUONTRIGGERLUT_H
#define ALIMUONTRIGGERLUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONTriggerLut
/// \brief MUON trigger look up table class 
///
/// \author Philippe Crochet

#include <TNamed.h>

class TH3;

//----------------------------------------------
class AliMUONTriggerLut : public TNamed 
{
 public: 
  AliMUONTriggerLut();    // constructor
  virtual ~AliMUONTriggerLut();   // destructor

  void ReadFromFile(const char* filename);
  
  void GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev, Int_t ystrip, 
		    Int_t lutLpt[2], Int_t lutHpt[2]);

 protected:
  // assignment operator
  AliMUONTriggerLut& operator=(const AliMUONTriggerLut& AliMUONTriggerLut); 
		
 private:
  // copy constructor
  AliMUONTriggerLut (const AliMUONTriggerLut& AliMUONTriggerLut);
  Int_t GetMask(Int_t ystrip);

private:
  TH3 *fLptPlus; ///< 3-d histogram with 234x32x31 bins Low pt Plus  
  TH3 *fLptMinu; ///< 3-d histogram with 234x32x31 bins Low pt Minus
  TH3 *fLptUnde; ///< 3-d histogram with 234x32x31 bins Low pt Undefined
  TH3 *fHptPlus; ///< 3-d histogram with 234x32x31 bins High pt Plus
  TH3 *fHptMinu; ///< 3-d histogram with 234x32x31 bins High pt Minus 
  TH3 *fHptUnde; ///< 3-d histogram with 234x32x31 bins High pt Undefined 
  TH3 *fAptPlus; ///< 3-d histogram with 234x32x31 bins All pt Plus  
  TH3 *fAptMinu; ///< 3-d histogram with 234x32x31 bins All pt Minus  
  TH3 *fAptUnde; ///< 3-d histogram with 234x32x31 bins All pt Undefined    

  ClassDef(AliMUONTriggerLut,1) // Trigger Look up Table class

};
#endif






