#ifndef ALITPCDIGITSARRAY_H
#define ALITPCDIGITSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigits.h" 
#include "AliSegmentArray.h"
#include "AliDigitsArray.h"
#include "AliTPCParam.h" 

#include <TMatrix.h>
#include <TTree.h>
#include <TClonesArray.h>


class TClonesArray;
class TObjArray;
class AliTPCPRF2D;
class AliTPCRF1D;

class AliTPCDigitsArray : public AliDigitsArray {
public:
  AliTPCDigitsArray(Bool_t sim=kTRUE);
  ~AliTPCDigitsArray();
  AliDigits *  GetRow(Int_t sector,Int_t row); //return pointer to row from array
  AliDigits *  CreateRow(Int_t sector, Int_t row); //
  AliDigits *  LoadRow(Int_t sector,Int_t row);
  Bool_t StoreRow(Int_t sector,Int_t row);
  Bool_t ClearRow(Int_t sector,Int_t row);
  Bool_t Setup(AliDetectorParam *param);  
  
  Bool_t IsSimulated(){return fBSim;}
  Bool_t  Update(); //
private:  
  //AliTPCPRF2D * fPRF;           //x and y pad response function object
  //AliTPCRF1D * fRF;             //z (time) response function object
  Bool_t fBSim;             //signalize if we have digits with track ID
  Int_t  fCompression;      //default compression for AliDigits - used in storing
  Int_t  fTrackLevel;        //default level for track ID storing
  ClassDef(AliTPCDigitsArray,1) 
};
  
#endif //ALITPCCLUSTERSARRAY_H
