#ifndef ALIDIGITSARRAY_H
#define ALIDIGITSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliDigitsArray        //
////////////////////////////////////////////////

#include "AliSegmentArray.h"
class AliDetectorParam;

class AliDigitsArray : public AliSegmentArray {
public:
  AliDigitsArray();
  AliDigitsArray(const AliDigitsArray &param); // copy constructor
 AliDigitsArray &operator = (const AliDigitsArray & param);
  virtual ~AliDigitsArray();  
  virtual   Bool_t Setup(AliDetectorParam *param);  //setup array according parameters
  const AliDetectorParam *  GetParam() {return fParam;} 
  virtual Bool_t SetParam(AliDetectorParam * param);
protected:  
  AliDetectorParam * fParam;      //pointer to detector parameters 
  ClassDef(AliDigitsArray,1) // Digits manager  
};
  
#endif //ALIDIGITSARRAY_H
