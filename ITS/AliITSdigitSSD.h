#ifndef AliITSDIGITSSD_H
#define AliITSDIGITSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliITS.h"

class AliITSdigitSSD: public AliITSdigit {

protected:       
  Bool_t  fSide;              //true if P side
  Int_t   fStripNumber;       //Number of Strip
  Int_t   fSignal;            //Signal on ADC 
    
public:
    
  AliITSdigitSSD() {};
  ~AliITSdigitSSD() {};
  AliITSdigitSSD(Int_t *tracks, Int_t *digits, Int_t strNo, Int_t s, Bool_t p);

  // Methods for accesing signal on strip				
  void   SetSignal(Int_t s) {fSignal = s;}           
  Int_t  GetSignal() const {return fSignal;}               		         
  
  // Methods for accesing strip number  
  Int_t  GetStripNumber() const {return fStripNumber;};
  void   SetStripNumber(Int_t s) {fStripNumber = s;};
    
  // Methods for accesing side of the strip P/N
  Bool_t IsSideP() const {return fSide;};              //returns true when side P
  void   SetSideP(Bool_t b) {fSide = b;};        //set side of digit
    	
  ClassDef(AliITSdigitSSD, 1)
};

#endif
