//-*- Mode: C++ -*-

#ifndef ALIANALYSISTASKTRIGGER_H
#define ALIANALYSISTASKTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Study trigger efficiencies for high-pt trigger
// Author: Jochen Thaeder <jochen@thaeder.de> 

#include "AliAnalysisTaskSE.h"

class TH1F;
class AliESDEvent;
class AliMCEvent;
class AliKineTrackCuts;

#include "AliStack.h"
#include "TParticle.h"
#include "TRandom3.h"

class AliAnalysisTasktrigger : public AliAnalysisTaskSE {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisTasktrigger(const char *name = "AliAnalysisTasktrigger");
  virtual ~AliAnalysisTasktrigger();
  
  /*
   * ---------------------------------------------------------------------------------
   *                                    Methods
   * ---------------------------------------------------------------------------------
   */

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisTasktrigger(const AliAnalysisTasktrigger&); // not implemented
  AliAnalysisTasktrigger& operator=(const AliAnalysisTasktrigger&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                            Setup Methods - private
   * ---------------------------------------------------------------------------------
   */

  Bool_t SetupEvent();

  void   SetupESDTrackCuts();  

  void   SetupTrigHistograms();
  void   SetupPtHistograms();
  void   SetupMultHistograms();

  /*
   * ---------------------------------------------------------------------------------
   *                            Helper Methods - private
   * ---------------------------------------------------------------------------------
   */
  
  TParticle* GetChargedPhysicalPrimary( AliStack* stack, Int_t idx );
  Bool_t IsFindableMC(Int_t idx, Float_t length);

  void AddTriggerHist( TH1F* hist );
  void AddPtHist( TH1F* hist );
  void AddMultHist( TH1F* hist );
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Trigger Methods - private
   * ---------------------------------------------------------------------------------
   */
  
  void EvaluateTrigger();

  /*
   * ---------------------------------------------------------------------------------
   *                              Fill Methods - private
   * ---------------------------------------------------------------------------------
   */

  void FillCutStudies( Int_t mode );
  void FillCounters( Int_t mode );
  void FillTriggerHistograms();

  void FillTriggerStudies();
  void FillTriggerStudiesMC();

  /*
   * ---------------------------------------------------------------------------------
   *                         Static Members - private
   * ---------------------------------------------------------------------------------
   */

  static const Int_t     fgkNSettings;       // N Settings
  static const Int_t     fgkNTrigger;        // N Trigger 

  static const Double_t  fgkTriggerPt[];     // Array of Pt Settings
  static const Char_t   *fgkTrigger[];       // Array of Trigger Names

  static const Int_t     fgkNSelectionCuts;  // N Selection and Cut Types
  static const Char_t   *fgkSelectionCuts[]; // Selection and Cut Types

  enum mode {kOFF, kHLT, kMC, kNModes};

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  TRandom3        *fRandom;                  //! Random Generator

  AliMCEvent      *fMC;                      //! MC object
  AliESDEvent     *fESD;                     //! ESD object
  AliESDEvent     *fESDHLT;                  //! ESD - HLT object

  AliESDtrackCuts *fESDTrackCuts;            //! ESD cuts  
  AliESDtrackCuts *fESDHLTTrackCuts;         //! HLT adopted track cuts 
  AliKineTrackCuts*fMCTrackCuts;             //! MC track cuts

  Bool_t           fIsSelected;              //! Event Selected by physics selection + Primary Vertex
  Bool_t           fIsSelectedHLT;           //! HLT Event Selected by physics selection + Primary Vertex
  Bool_t           fIsSelectedMC;            //! MC Event Selected by physics selection + Primary Vertex

  Bool_t           fIsSelectedTask;          //! Event Selected by physics selection

  TObjArray       *fOutputContainer;         //! output data container

  // --------------------------------------------------------------------

  Int_t           *fPtCount;                 //! Pt count per setting [kNModes*fgkTrigger]
  Int_t           *fMultiplicity;            //! Multiplicty counters [fgkNSelectionCuts]

  // --------------------------------------------------------------------

  Bool_t          *fTrigger;                 // Trigger [kNModes*fgkNTrigger]

  ClassDef(AliAnalysisTasktrigger, 3);
};

#endif
