#ifndef ALIFITALIGN_H
#define ALIFITALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//     This class creates the alignment object from the surveyer data  //      
//     for the FIT 				                       //
/////////////////////////////////////////////////////////////////////////
#include "AliAlignObjParams.h"
#include <TMatrixDfwd.h>
#include <TMatrixT.h>

// Class creating the FIT aligmnent objects
// from the surveys done by surveyers at Point2.
// Position of FIT alignment objects is computed.

class AliFITAlign : public TObject {

 public:
  AliFITAlign();
  AliFITAlign(Int_t reportloc, Int_t side, Int_t reportglob); 
  // reportloc - location of the report: 0 - DCDB (Grid), 1,2 ... - file on local disc
  // side - side in ALICE: 0 - A-side or 1 - C-side
  // reportglob - EDMS number of the survey report
  AliFITAlign(const AliFITAlign &align); // copy constructor
  AliFITAlign &operator = (const AliFITAlign &align); //assignment operator
  virtual   ~AliFITAlign();
  // void SetDebug(Int_t debug){debug=fDebug;}
  Bool_t LoadSurveyData();
  Double_t ComputePosition();
  void CreateAlignObj();
  void Run();
  void StoreAlignObj();
  //
 private:

  Char_t *fFileGlob;                  // file with surveyed points
  AliAlignObjParams *fFITAAlignObj;  // FIT-A alignment object
  AliAlignObjParams *fFITCAlignObj;  // FIT-C alignment object
  Int_t fDebug;                     // debug flag
  Float_t fXPosC;		    // "x" coordinate of FIT-C with respect to Global Reference System  
  Float_t fYPosC;  		    // "y" coordinate of FIT-C with respect to Global Reference System
  Float_t fXPosA;                    // "x" coordinate of FIT-A with respect to Global Reference System
  Float_t fYPosA;                    // "y" coordinate of FIT-A with respect to Global Reference System

  Int_t fRepLoc;		    // Location of the report: 0 - DCDB (Grid), 1,2 ... - file on local disc    
  Int_t fRepGlob;  		    // Number of the survey report	
  Int_t fSide;			    // Side in ALICE: 0: A-side or 1: C-side	
  Char_t *fUser;                    // AliEn user

  ClassDef(AliFITAlign,0);
};
#endif
