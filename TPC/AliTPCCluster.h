#ifndef TPCClusters_H
#define TPCClusters_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h" 
#include "AliSegmentArray.h"
#include "AliTPCParam.h" 

#include <TMatrix.h>
#include <TTree.h>
#include <TClonesArray.h>


class TClonesArray;
class TObjArray;


class AliTPCClustersRow : public AliSegmentID{
public:
  AliTPCClustersRow();
  AliTPCClustersRow(Int_t size);
  void InsertCluster(const AliTPCcluster* c ); //insert copy of cluster  
  const AliTPCcluster* operator[](Int_t i); 
  Int_t  Find(Double_t y) const;   //find nearest cluster in y direction
  void Sort();
public:
  TClonesArray * fClusters;  
  Int_t  fNclusters;  
  ClassDef(AliTPCClustersRow,1) // Cluster manager 
};


class AliTPCClustersArray : public AliSegmentArray {
public:
  AliTPCClustersArray();
  ~AliTPCClustersArray();
  const AliTPCClustersRow *  GetRow(Int_t sector,Int_t row); 
  Bool_t LoadRow(Int_t sector,Int_t row);
  Bool_t StoreRow(Int_t sector,Int_t row);
  Bool_t Setup(AliTPCParam *param);  
  //construct array  according parameters in fParam   
  const AliTPCParam & GetParam() {return *fParam;} 
private:  
  AliSegmentID * NewSegment(){ return (AliSegmentID*)new AliTPCClustersRow;}
  AliTPCParam * fParam;      //pointer to TPC parameters
  //AliTPCClustersRow ** fRow;  //pointer to array of pointers to cluster row
  ClassDef(AliTPCClustersArray,1) // Cluster manager
};
  
#endif

