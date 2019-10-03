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

#ifndef ALICFV0TASK_H
#define ALICFV0TASK_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class AliMCEventHandler;
class AliCFManager;
class AliESDv0;
class AliESDEvent;
class AliAODv0;

class AliCFV0Task : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepReconstructible = 1,
    kStepReconstructed   = 2,
    kStepSelected        = 3
  };

  AliCFV0Task();
  AliCFV0Task(const Char_t* name);
  AliCFV0Task& operator= (const AliCFV0Task& c);
  AliCFV0Task(const AliCFV0Task& c);
  virtual ~AliCFV0Task();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager()                 {return fCFManager;} // get corr manager

  void     SetRebuildV0s(Bool_t flag)       {fRebuildV0s = flag;}       // setter for V0 on-the-fly reconstruction
  void     SetV0PDG(Int_t code)             {fV0PDG = code; }           // defines the PDG code of searched V0's
  Int_t    IsMcV0(AliESDv0*) const ;           // checks if the AliESDv0 can be associated, returns mother label
  Int_t    IsMcV0(AliAODv0*) const ;           // checks if the AliAODv0 can be associated, returns mother label
  Int_t    GetV0Label(UInt_t, UInt_t) const ;                // returns label of V0 given 2 daughter labels
  static   Double_t GetRapidity(Double_t, Double_t) ;  // returns the rapidity of the V0 (assuming PDG code)


 protected:
  void          RebuildV0s(AliESDEvent*) ;   // reconstructs V0's on-fly

  Bool_t          fRebuildV0s;   //  flag for on-the-fly V0 reconstruction
  Int_t           fV0PDG;        //  PDG code of searched V0's
  AliCFManager   *fCFManager  ;  //  pointer to the CF manager

  // Histograms
  //Number of events
  TH1I *fHistEventsProcessed; //! simple histo for monitoring the number of events processed
  
  ClassDef(AliCFV0Task,1);
};

#endif
