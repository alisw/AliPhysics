#ifndef ALIFITDIGIT_H
#define ALIFITDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/***********************************************************************
class for FIT Digits object
 Alla.Maevskaya@cern.ch 
***********************************************************************/
#include <AliDigit.h>
class AliFITDigit: public AliDigit {

 public:
  AliFITDigit();
  AliFITDigit(Int_t npmt, 
	      Int_t timeCFD, Int_t timeLED, Int_t timeQT0, Int_t timeQT1 , 
	      Int_t *labels=0);
  virtual ~AliFITDigit();
  Int_t TimeCFD() {return fTimeCFD;}
  Int_t TimeLED() {return fTimeLED;}
  Int_t TimeQT0 () {return fTimeQT0;}
  Int_t TimeQT1 () {return fTimeQT1;}
  Int_t NPMT ()  {return fNPMT;}

 private: 
  
  Int_t fTimeCFD;    // array's CFD
  Int_t fTimeLED;    // array's LED
  Int_t fTimeQT0;    // array's QT0 start QTC amplitude
  Int_t fTimeQT1;    // array's QT1 stop QTC amplitude
  Int_t fNPMT;       // number of channel [0,79] C side; [80,159] Aside
 
  AliFITDigit( const AliFITDigit& );
  AliFITDigit& operator=(const AliFITDigit&); 
  
 
  ClassDef(AliFITDigit,1)  //Digit (Header) object for set FIT
};


#endif




