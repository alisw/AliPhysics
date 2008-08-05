/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------

#ifndef ALIPROTONCORRECTIONTASK_H
#define ALIPROTONCORRECTIONTASK_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliESDtrack;
class AliVParticle;

class AliProtonCorrectionTask : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepReconstructible = 1,
    kStepReconstructed   = 2,
    kStepSelected        = 3
  };

  AliProtonCorrectionTask();
  AliProtonCorrectionTask(const Char_t* name);
  AliProtonCorrectionTask& operator= (const AliProtonCorrectionTask& c);
  AliProtonCorrectionTask(const AliProtonCorrectionTask& c);
  virtual ~AliProtonCorrectionTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManagerProtons(AliCFManager* io) {fCFManagerProtons = io;}   // global correction manager
  AliCFManager * GetCFManagerProtons() const {return fCFManagerProtons;}           // get corr manager
  void           SetCFManagerAntiProtons(AliCFManager* io) {fCFManagerAntiProtons = io;} // global correction manager
  AliCFManager * GetCFManagerAntiProtons() const {return fCFManagerAntiProtons;}         // get corr manager
  void           SetQAList(TList* list) {fQAHistList = list;}

  // Data types
  Bool_t IsReadTPCTracks() const {return fReadTPCTracks;}
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadTPCTracks (Bool_t flag=kTRUE) {fReadTPCTracks=flag;}
  void   SetReadAODData   (Bool_t flag=kTRUE) {fReadAODData=flag;}

 protected:
  Double_t Rapidity(Double_t px, Double_t py, Double_t pz);

  Bool_t          fReadTPCTracks;         // flag to loop on TPC tracks only
  Bool_t          fReadAODData;           // flag for AOD/ESD input files
  AliCFManager   *fCFManagerProtons;      // pointer to the CF manager
  AliCFManager   *fCFManagerAntiProtons;  // pointer to the CF manager
  TList          *fQAHistList;            // list of QA histograms

  // Histograms
  //Number of events
  TH1I  *fHistEventsProcessed; // simple histo for monitoring the number of events processed
  
  ClassDef(AliProtonCorrectionTask,1);
};

#endif
