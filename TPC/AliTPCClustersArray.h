#ifndef ALITPCCLUSTERSARRAY_H
#define ALITPCCLUSTERSARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////


#include "AliClustersArray.h"
 




class AliTPCClustersRow;

class AliTPCClustersArray : public AliClustersArray {
public:
  AliTPCClustersArray();
  virtual ~AliTPCClustersArray();
  AliTPCClustersRow * GetRow(Int_t sector,Int_t row);  
  AliTPCClustersRow * CreateRow(Int_t sector, Int_t row); //
  AliTPCClustersRow * LoadRow(Int_t sector,Int_t row);
  Bool_t StoreRow(Int_t sector,Int_t row);
  Bool_t ClearRow(Int_t sector,Int_t row);
  Bool_t Setup(const AliDetectorParam *param);     
  //construct array  according parameters in fParam   
  Bool_t  Update(); //blabla 
  AliSegmentID * NewSegment(); //create new segment - AliTPCClustersRow
protected:
  //void MakeTree(); 
 
private:   
  ClassDef(AliTPCClustersArray,1) // Cluster manager 
};
  
#endif //ALITPCCLUSTERSARRAY_H
