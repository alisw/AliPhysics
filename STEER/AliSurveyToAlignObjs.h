#ifndef ALISURVEYTOALIGNOBJS_H
#define ALISURVEYTOALIGNOBJS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//     Base class for creating the alignment objects from survey       //
//     for a given subdetector                                         //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
#include <TSystem.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliAlignObjParams.h"
#include "AliSurveyObj.h"


class AliSurveyToAlignObjs : public TObject {

 public:
  AliSurveyToAlignObjs();
  AliSurveyToAlignObjs(const AliSurveyToAlignObjs &align); // copy constructor
  AliSurveyToAlignObjs &operator = (const AliSurveyToAlignObjs &align); //assignment operator

  Bool_t LoadSurveyFromLocalFile(const char* filename);
  Bool_t LoadSurveyFromAlienFile(const char* det, Int_t repNum, Int_t repVersion);

  virtual Bool_t CreateAlignObjs() = 0;
  virtual void Run() = 0;
  Bool_t StoreAlignObjToFile(const char* filename, const char* det);
  Bool_t StoreAlignObjToCDB(const char* cdbFolder, const char* det);
  virtual   ~AliSurveyToAlignObjs();
  //
 protected:
  AliSurveyObj      *fSurveyObj;         // survey object
  TObjArray         *fSurveyPoints;     // array of survey points
  TClonesArray      *fAlignObjArray;    // TClonesArray of alignment objects
  AliAlignObjParams *fAlignObj;         // alignment object

  ClassDef(AliSurveyToAlignObjs,0);
};
#endif
