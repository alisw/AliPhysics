#ifndef ALIMUONTRIGGERLUT_H
#define ALIMUONTRIGGERLUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONTriggerLut
/// \brief MUON trigger look up table class 

#include <TNamed.h>

class TH3S;

//----------------------------------------------
class AliMUONTriggerLut : public TNamed 
{
 public: 
  AliMUONTriggerLut();    // constructor
  ~AliMUONTriggerLut();   // destructor

  void LoadLut();
  
  void GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev, Int_t ystrip, 
		    Int_t lutLpt[2], Int_t lutHpt[2], Int_t lutApt[2]);

 protected:
  // copy constructor
  AliMUONTriggerLut (const AliMUONTriggerLut& AliMUONTriggerLut); 
  // assignment operator
  AliMUONTriggerLut& operator=(const AliMUONTriggerLut& AliMUONTriggerLut); 
		
 private:
  Int_t GetMask(Int_t ystrip);

  ClassDef(AliMUONTriggerLut,1) // Trigger Look up Table class

    private:
  TH3S *fLptPlus; //3-d histogram with 234x32x31 bins Low pt Plus  
  TH3S *fLptMinu; //3-d histogram with 234x32x31 bins Low pt Minus
  TH3S *fLptUnde; //3-d histogram with 234x32x31 bins Low pt Undefined
  TH3S *fHptPlus; //3-d histogram with 234x32x31 bins High pt Plus
  TH3S *fHptMinu; //3-d histogram with 234x32x31 bins High pt Minus 
  TH3S *fHptUnde; //3-d histogram with 234x32x31 bins High pt Undefined 
  TH3S *fAptPlus; //3-d histogram with 234x32x31 bins All pt Plus  
  TH3S *fAptMinu; //3-d histogram with 234x32x31 bins All pt Minus  
  TH3S *fAptUnde; //3-d histogram with 234x32x31 bins All pt Undefined    

};
#endif






