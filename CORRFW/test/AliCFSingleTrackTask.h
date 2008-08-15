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

#ifndef ALICFSINGLETRACKTASK_H
#define ALICFSINGLETRACKTASK_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliESDtrack;
class AliVParticle;

class AliCFSingleTrackTask : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepReconstructible = 1,
    kStepReconstructed   = 2,
    kStepSelected        = 3
  };

  AliCFSingleTrackTask();
  AliCFSingleTrackTask(const Char_t* name);
  AliCFSingleTrackTask& operator= (const AliCFSingleTrackTask& c);
  AliCFSingleTrackTask(const AliCFSingleTrackTask& c);
  virtual ~AliCFSingleTrackTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}           // get corr manager
  void           SetQAList(TList* list) {fQAHistList = list;}

  // Data types
  Bool_t IsReadTPCTracks() const {return fReadTPCTracks;}
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadTPCTracks (Bool_t flag=kTRUE) {fReadTPCTracks=flag;}
  void   SetReadAODData   (Bool_t flag=kTRUE) {fReadAODData=flag;}

 protected:
  
  Bool_t          fReadTPCTracks ; // flag to loop on TPC tracks only (ESD mode only)
  Bool_t          fReadAODData ;   // flag for AOD/ESD input files
  AliCFManager   *fCFManager    ;  // pointer to the CF manager
  TList          *fQAHistList   ;  // list of QA histograms

  // Histograms
  //Number of events
  TH1I  *fHistEventsProcessed; // simple histo for monitoring the number of events processed
  
  ClassDef(AliCFSingleTrackTask,1);
};

#endif
