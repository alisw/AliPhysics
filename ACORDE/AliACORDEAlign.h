#ifndef ALIACORDEALIGN_H
#define ALIACORDEALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//     This class creates the alignment object from the surveyer data  //      
//     for the ACORDE 				                       //
/////////////////////////////////////////////////////////////////////////

#include "AliAlignObjParams.h"
#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include "AliSurveyToAlignObjs.h"

// Class creating the ACORDE aligmnent objects
// from the surveys done by surveyers at Point2.

class AliACORDEAlign : public  TObject{

 public:
  AliACORDEAlign();
  AliACORDEAlign(Int_t reportloc,Int_t reportglob);
  AliACORDEAlign(const AliACORDEAlign &align); // copy constructor
  AliACORDEAlign &operator = (const AliACORDEAlign &align); //assignment operator
  void  ComputePosition();
  void Run();
  void LoadSurveyData();
//  void CreateACORDEAlignObjs();
  void StoreAlignObj();
  void SetDebug(Int_t debug){fDebug=debug;}
  virtual ~AliACORDEAlign();
  //
 private:

  Char_t *fFileGlob;
  Int_t fRepLoc;
  Int_t fRepGlob;
  Char_t *fUser;
  TMatrixD fX;
//  AliAlignObjParams *fACORDEAlignObj;
  TObjArray  *fAlignACORDEObjArray;
  Int_t fDebug;                     // debug flag
                  //measurements
  ClassDef(AliACORDEAlign,0);
  
};
#endif
