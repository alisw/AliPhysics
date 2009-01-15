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
// Author : C. Zampolli, on an example from R. Vernet
//-----------------------------------------------------------------------

#ifndef ALICFHEAVYFLAVOURTASK_H
#define ALICFHEAVYFLAVOURTASK_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class AliCFManager;
class AliAODRecoDecayHF2Prong;

class AliCFHeavyFlavourTask : public AliAnalysisTaskSE {
  public:

  enum {
    kStepGenerated       = 0,
    kStepReconstructed   = 1,
  };

  AliCFHeavyFlavourTask();
  AliCFHeavyFlavourTask(const Char_t* name);
  AliCFHeavyFlavourTask& operator= (const AliCFHeavyFlavourTask& c);
  AliCFHeavyFlavourTask(const AliCFHeavyFlavourTask& c);
  virtual ~AliCFHeavyFlavourTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void SetCFManager(AliCFManager* const io) {fCFManager = io;}   // global correction manager
  AliCFManager* GetCFManager() const {return fCFManager;} // get corr manager

  void SetPDG(Int_t code) {fPDG = code; } // defines the PDG code of searched HF
  Int_t IsMcVtx(AliAODRecoDecayHF2Prong* const vtx) const ; // checks if the AliAODRecoDecayHF2Prong can be associated to an MC particle, returns mother label
  Int_t GetVtxLabel(Int_t* labels) const ; // returns label of vertex given the daughter labels

  
 protected:
  Int_t fPDG;         //  PDG code of searched V0's
  AliCFManager* fCFManager  ; //  pointer to the CF manager
  TH1I* fHistEventsProcessed;   //! simple histo for monitoring the number of events processed
  Int_t fCountMC;               // MC particle found
  Int_t fEvents;                // n. of events
  
  ClassDef(AliCFHeavyFlavourTask,0);
};

#endif
