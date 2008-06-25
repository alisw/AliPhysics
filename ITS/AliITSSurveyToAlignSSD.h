#ifndef ALIITSSURVEYTOALIGNSSD_H
#define ALIITSSURVEYTOALIGNSSD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//     This class creates the alignment object from the survey data    //      
//     for the SSD layers.                                             //
//                                                                     //
//     Author: Panos.Christakoglou@cern.ch - NIKHEF/UU                 //
/////////////////////////////////////////////////////////////////////////
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliAlignObjParams.h"


class AliITSSurveyToAlignSSD : public TObject {

 public:
  AliITSSurveyToAlignSSD();
  AliITSSurveyToAlignSSD(Int_t run, Int_t reportloc, Int_t reportglob);
  AliITSSurveyToAlignSSD(const AliITSSurveyToAlignSSD &align); // copy constructor
  AliITSSurveyToAlignSSD &operator = (const AliITSSurveyToAlignSSD &align); //assignment operator
  Bool_t LoadSurveyData();
  void CreateAlignObj();
  void Run();
  void SetDebug(Int_t debug){fDebug=debug;}
  void StoreAlignObj();
  virtual   ~AliITSSurveyToAlignSSD();
  //
 private:
  Int_t   fRun;                         // the run number for the OCDB
  Char_t *fFileLoc;                     // file with ideal points
  Char_t *fFileGlob;                    // file with surveyed points
  TObjArray *fSurveyPoints;             // array of survey points
  TClonesArray *fSSDAlignObj;           // TClonesArray of AliAlignObj
  AliAlignObjParams *fSSDAlignObjParam; // SSD alignment object param
  Int_t fDebug;                         // debug flag
  

  ClassDef(AliITSSurveyToAlignSSD,0);
};
#endif
