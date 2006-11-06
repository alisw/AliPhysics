#ifndef ALICLUSTERS_H
#define ALICLUSTERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC   clusters          //
////////////////////////////////////////////////


#include "AliSegmentID.h"
#include <TClonesArray.h>

class TObjArray;


class AliClusters : public AliSegmentID{
public:
  AliClusters();
  AliClusters(const AliClusters &param); // copy constructor
  AliClusters &operator = (const AliClusters & param); 
  ~AliClusters();
  virtual TObject* InsertCluster(const TObject* c ); //insert copy of cluster  
  const TObject* operator[](Int_t i); 
  virtual Int_t  Find(Double_t y) const;   //find nearest cluster in y direction
  void Sort();
  TClonesArray * GetArray(){return fClusters;}
  void SetArray(Int_t length); //construct clonnes array of objects of type fClass
  void DrawClusters(Float_t shiftx, Float_t shifty, Int_t color, Int_t size, Int_t style);
  Bool_t SetClass(const Text_t *classname);
protected:
  TClonesArray * fClusters;  
  Int_t  fNclusters;  
  TClass * fClass; //!type of cluster class 
  ClassDef(AliClusters,1) // Cluster manager
};


#endif //ALICLUSTERS_H
