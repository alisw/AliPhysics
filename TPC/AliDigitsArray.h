#ifndef ALIDIGITSARRAY_H
#define ALIDIGITSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliDigitsArray        //
////////////////////////////////////////////////
class AliDetectorParam;

class AliDigitsArray : public AliSegmentArray {
public:
  AliDigitsArray();
  ~AliDigitsArray();  
  virtual   Bool_t Setup(AliDetectorParam *param);  //setup array according parameters
  const AliDetectorParam *  GetParam() {return fParam;} 
  virtual Bool_t SetParam(AliDetectorParam * param);
protected:  
  AliDetectorParam * fParam;      //pointer to detector parameters 
  ClassDef(AliDigitsArray,1) 
};
  
#endif //ALIDIGITSARRAY_H
