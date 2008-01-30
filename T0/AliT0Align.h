#ifndef ALIT0ALIGN_H
#define ALIT0ALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//     This class creates the alignment object from the surveyer data  //      
//     for the T0 				                       //
/////////////////////////////////////////////////////////////////////////
#include "AliAlignObjParams.h"
#include <TMatrixDfwd.h>
#include <TMatrixT.h>

// Class creating the T0 aligmnent objects
// from the surveys done by surveyers at Point2.
// Position of T0 alignment objects is computed.

class AliT0Align : public TObject {

 public:
  AliT0Align();
  AliT0Align(Int_t reportloc, Int_t reportglob);
  AliT0Align(const AliT0Align &align); // copy constructor
  AliT0Align &operator = (const AliT0Align &align); //assignment operator
  Bool_t LoadSurveyData();
  Double_t ComputePosition();
  void CreateAlignObj();
  void Run();
  void SetDebug(Int_t debug){debug=fDebug;}
  void StoreAlignObj();
  virtual   ~AliT0Align();
  //
 private:

  Char_t *fFileGlob;                  // file with surveyed points
  AliAlignObjParams *fT0AAlignObj;  // T0-A alignment object
  AliAlignObjParams *fT0CAlignObj;  // T0-C alignment object
  Int_t fDebug;                     // debug flag
  Float_t fXPos;		    // "x" coordinate of T0-C with respect to Global Reference System  
  Float_t fYPos;  		    // "y" coordinate of T0-C with respect to Global Reference System
  Int_t fRepLoc;		    // Location of the report: 0 - DCDB (Grid), 1,2 ... - file on local disc    

  ClassDef(AliT0Align,0);
};
#endif
