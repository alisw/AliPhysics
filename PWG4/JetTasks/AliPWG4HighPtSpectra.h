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
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTSPECTRA_H
#define ALIPWG4HIGHPTSPECTRA_H

#include "AliAnalysisTask.h"
#include "AliCFManager.h"

class TH1I;
class TH1F;
class TH1D;
class TFile ;
//class AliCFManager;
class AliESDtrackCuts;
class AliESDEvent;

class AliPWG4HighPtSpectra : public AliAnalysisTask {
 public:

  enum {
    kStepReconstructed        = 0,
    kStepReconstructedTPCOnly = 1,
    kStepSecondaries          = 2,
    kStepMCtrackable          = 3,
    kStepReconstructedMC      = 4,
    kStepMCAcceptance         = 5
  };

  AliPWG4HighPtSpectra();
  AliPWG4HighPtSpectra(const Char_t* name);
  // AliPWG4HighPtSpectra& operator= (const AliPWG4HighPtSpectra& c);
  //  AliPWG4HighPtSpectra(const AliPWG4HighPtSpectra& c);
  ~AliPWG4HighPtSpectra() {;};

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void     SetCFManagerPos(const AliCFManager* io1) {fCFManagerPos = io1;}   // global correction manager 
  const AliCFManager * GetCFManagerPos() const {return fCFManagerPos;}           // get corr manager 
  void     SetCFManagerNeg(const AliCFManager* io2) {fCFManagerNeg = io2;}   // global correction manager 
  const AliCFManager * GetCFManagerNeg() const {return fCFManagerNeg;}            // get corr manager
  
  //AliESDtrackCuts setters
  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  //Select the trigger
  void SelectTrigger(Int_t trig) { fTrigger = trig; } 

  // Data types
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadAODData(Bool_t flag=kTRUE) {fReadAODData=flag;}
  
 protected:
  Bool_t              fReadAODData ;       // flag for AOD/ESD input files
  const AliCFManager  *fCFManagerPos    ;  // pointer to the CF manager for positive charged particles
  const AliCFManager  *fCFManagerNeg    ;  // pointer to the CF manager for negative charged particles
 
  AliESDEvent *fESD;              //! ESD object
  //AliESDtrackCuts options. Must be setted in AddTaskPWG4HighPtQAMC.C. They correspond with different steps in container.
  AliESDtrackCuts *fTrackCuts;    // trackCuts applied
  Int_t fTrigger;                 //Trigger flag as defined in AliAnalysisHelperJetTasks.h 

 private:
 AliPWG4HighPtSpectra(const AliPWG4HighPtSpectra&);
 AliPWG4HighPtSpectra& operator=(const AliPWG4HighPtSpectra&);


  // Histograms
  //Number of events
  TList *fHistList;            //! List of output histograms
  TH1F *fNEventAll;            //! Event counter
  TH1F *fNEventSel;            //! Event counter: Selected events for analysis

  ClassDef(AliPWG4HighPtSpectra,1);
};

#endif
