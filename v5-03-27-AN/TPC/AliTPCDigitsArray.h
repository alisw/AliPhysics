#ifndef ALITPCDIGITSARRAY_H
#define ALITPCDIGITSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC digits                   //
////////////////////////////////////////////////


#include "AliDigits.h" 
#include "AliDigitsArray.h"
#include "AliTPCParam.h" 

class AliDigits;
class AliDetectorParam;

class AliTPCDigitsArray : public AliDigitsArray {
public:
  AliTPCDigitsArray(Bool_t sim=kTRUE);
  virtual   ~AliTPCDigitsArray();
  AliDigits *  GetRow(Int_t sector,Int_t row); //return pointer to row from array
  AliDigits *  CreateRow(Int_t sector, Int_t row); //
  AliDigits *  LoadRow(Int_t sector,Int_t row);
  Bool_t StoreRow(Int_t sector,Int_t row);
  Bool_t ClearRow(Int_t sector,Int_t row);
  Bool_t Setup(AliDetectorParam *param);  
  
  Bool_t IsSimulated(){return fBSim;}
  Bool_t  Update(); 
private:  
  Bool_t fBSim;             //signalize if we have digits with track ID
  Int_t  fCompression;      //default compression for AliDigits - used in storing
  Int_t  fTrackLevel;        //default level for track ID storing
  ClassDef(AliTPCDigitsArray,1) // TPC digits manager
};
  
#endif //ALITPCCLUSTERSARRAY_H
