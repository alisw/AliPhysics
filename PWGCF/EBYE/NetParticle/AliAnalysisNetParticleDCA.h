//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEDCA_H
#define ALIANALYSISNETPARTICLEDCA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Efficiency and Contaminations for NetParticle Distributions
// Authors: Jochen Thaeder <jochen@thaeder.de>

#include "THnSparse.h"

#include "AliAnalysisNetParticleHelper.h"

class AliESDEvent;
class AliESDInputHandler;
class AliMCEvent;
class AliStack;

class AliAnalysisNetParticleDCA : public TNamed {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleDCA();
  virtual ~AliAnalysisNetParticleDCA();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize */
  void Initialize(AliESDtrackCuts *cuts, AliESDtrackCuts *cutsBkg, AliAnalysisNetParticleHelper* helper);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliMCEvent *mcEvent);

  /** Reset Event */
  void ResetEvent();

  /** Process Event */
  void Process();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to DCA THnSparse */
  THnSparseF* GetHnDCA() {return fHnDCA;}

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleDCA(const AliAnalysisNetParticleDCA&); // not implemented
  AliAnalysisNetParticleDCA& operator=(const AliAnalysisNetParticleDCA&); // not implemented


#if 0
  /*
   * ---------------------------------------------------------------------------------
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Create the efficiency / contamination THnSparseF */
  void CreateHistograms();

  /** Fill DCA ThnSparse */
  void FillDCA(); 

  /** Check if particle is contamination */
  void CheckDCATrack(Int_t label, Float_t sign, Int_t idxTrack);
      
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
#endif
  AliAnalysisNetParticleHelper *fHelper;      //! Ptr to helper class

  // --- ESD only ----------------------------------------------------------

  AliESDEvent        *fESD;                   //! ESD object
  AliESDtrackCuts    *fESDTrackCuts;          //! ESD cuts  
  AliESDtrackCuts    *fESDTrackCutsBkg;       //! ESD cuts  
  
  // --- MC only -----------------------------------------------------------

  AliStack           *fStack;                 //! Ptr to stack
  AliMCEvent         *fMCEvent;               //! Ptr to MC event

  // -----------------------------------------------------------------------

  THnSparseF         *fHnDCA;                 //  THnSparseF contamination DCA

  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleDCA, 1);
};

#endif
