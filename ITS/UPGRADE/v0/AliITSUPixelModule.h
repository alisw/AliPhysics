#ifndef ALIITSUPIXELMODULE_H
#define ALIITSUPIXELMODULE_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////
//  ITS Cluster Finder Module  Container      //
//  Used in the ClusterFinder to keep         //
//  clusters in one module before storing     //
//  them in treeR                             //
////////////////////////////////////////////////

#include <TObject.h>

class AliITSUPixelModule :public TObject {

 public :

  enum {kMaxLab=12}; // maximum number of MC labels associated to the cluster
  AliITSUPixelModule();
  AliITSUPixelModule( UShort_t module, UInt_t col, UInt_t row, UInt_t charge, Int_t lab[kMaxLab]);

  virtual ~AliITSUPixelModule(){;}

  void SetPixelCoord( UInt_t col, UInt_t row) {fCol = col ; fRow = row; }
  
  void SetLabels(Int_t lab[kMaxLab]);
  void SetCharge(UInt_t charge) {fCharge = charge;}

  UShort_t GetModule() const {return fModule; }
  UInt_t GetCol() const {return fCol; }
  UInt_t GetRow() const {return fRow; }
  UInt_t GetCharge() const {return fCharge;}
  Int_t GetLabel(Int_t i) const {return fLabels[i];}
  void  PrintInfo();

 protected:
  UInt_t fCharge;
  UShort_t fModule;
  UInt_t fCol;
  UInt_t fRow;
  Int_t fLabels[kMaxLab];


  AliITSUPixelModule(const AliITSUPixelModule &p); // copy constructor
  AliITSUPixelModule& operator=(const AliITSUPixelModule &p);  // assignment operator


  ClassDef(AliITSUPixelModule,0)

    };
#endif 

