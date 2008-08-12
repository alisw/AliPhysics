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
#include <TSystem.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliAlignObjParams.h"


class AliITSSurveyToAlignSSD : public TObject {

 public:
  AliITSSurveyToAlignSSD();
  AliITSSurveyToAlignSSD(Int_t reportId, Int_t runNumber);
  AliITSSurveyToAlignSSD(const AliITSSurveyToAlignSSD &align); // copy constructor
  AliITSSurveyToAlignSSD &operator = (const AliITSSurveyToAlignSSD &align); //assignment operator

  void SetRun(Int_t runNumber) {fRun = runNumber;}
  void SetReportId(Int_t reportId) {
      fFileLocal = new Char_t[80];
      Char_t path[50];
      sprintf(path,gSystem->Getenv("ALICE_ROOT")); 
    sprintf(fFileLocal,"%s/ITS/Survey_SSD_%d.txt",path,reportId);  
  }
  Bool_t LoadSurveyData();
  void CreateAlignObj();
  void Run();
  void SetDebug(Int_t debug){fDebug=debug;}
  void StoreAlignObj();
  virtual   ~AliITSSurveyToAlignSSD();
  //
 private:
  Int_t fRun;                           // run number 
  Char_t *fFileLocal;                   // local file with surveyed points
  Char_t *fFileGrid;                    // GRID file with surveyed points
  TObjArray *fSurveyPoints;             // array of survey points
  TClonesArray *fSSDAlignObj;           // TClonesArray of AliAlignObj
  AliAlignObjParams *fSSDAlignObjParam; // SSD alignment object param
  Int_t fDebug;                         // debug flag
  

  ClassDef(AliITSSurveyToAlignSSD,0);
};
#endif
