#ifndef ALITPCALIGN_H
#define ALITPCALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//     This class creates the alignment object from the surveyer data  //      
//     for the entire TPC in the magnet (ALICE R.F.)                   //
/////////////////////////////////////////////////////////////////////////
#include "AliAlignObjParams.h"
#include <TMatrixDfwd.h>
#include <TMatrixT.h>

class AliTPCAlign : public TObject {

 public:
  AliTPCAlign();
  AliTPCAlign(Int_t reportloc, Int_t reportglob);
  AliTPCAlign(const AliTPCAlign &align); // copy constructor
  AliTPCAlign &operator = (const AliTPCAlign &align); //assignment operator
  Bool_t LoadSurveyData();
  Double_t ComputeTransform();
  void CreateAlignObj();
  void Run();
  void SetDebug(Int_t debug){fDebug=debug;}
  void StoreAlignObj();
  virtual   ~AliTPCAlign();
  //
 private:

  Char_t *fFileLoc;                   // file with ideal points
  Char_t *fFileGlob;                  // file with surveyed points
  AliAlignObjParams *fTPCAlignObj;             // TPC alignment object
  TMatrixD fX;                      // transformation parameters
  TMatrixD fA;                     // coefficients
  TMatrixD fY;                      // "measurements"
  Int_t fDebug;                     // debug flag
  

  ClassDef(AliTPCAlign,0);
};
#endif
