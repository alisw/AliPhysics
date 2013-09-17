//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEQA_H
#define ALIANALYSISNETPARTICLEQA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/**
 * Class for for NetParticle QA
 * -- Create input for QA
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "THnSparse.h"
#include "TH1F.h"
#include "TF1.h"

#include "AliAnalysisNetParticleHelper.h"

class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDInputHandler;
class AliESDtrackCuts;
class AliAODInputHandler;

class AliAnalysisNetParticleQA : public TNamed {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleQA();
  virtual ~AliAnalysisNetParticleQA();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize */
  Int_t Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts, Bool_t isMC, Int_t trackCutBit);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent);

  /** Reset Event */
  void ResetEvent();

  /** Process NetParticle QAs */ 
  Int_t Process();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to efficiency THnSparse */
  THnSparseF* GetHnQA()  {return fHnQA;}


  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleQA(const AliAnalysisNetParticleQA&); // not implemented
  AliAnalysisNetParticleQA& operator=(const AliAnalysisNetParticleQA&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                           Process - Private
   * ---------------------------------------------------------------------------------
   */

  /** Create the QA THnSparseF */
  void CreateHistograms();

  /** Process ESD/AOD tracks and fill histograms */
  Int_t ProcessTracks();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleHelper *fHelper;        //! Ptr to helper class
  // =======================================================================
  Int_t                 fPdgCode;               // PDG code of particle to be found 
  // =======================================================================
  AliESDEvent          *fESD;                   //! Ptr to ESD event
  AliESDtrackCuts      *fESDTrackCuts;          //! ESD cuts  
  AliPIDResponse       *fPIDResponse;           //! Ptr to PID response Object
  // -----------------------------------------------------------------------
  AliAODEvent          *fAOD;                   //! Ptr to AOD event
  TClonesArray         *fArrayMC;               //! array of MC particles
  Int_t                 fAODtrackCutBit;        //  Track filter bit for AOD tracks
  // -----------------------------------------------------------------------
  Bool_t                fIsMC;                  //  Is MC event
  AliMCEvent           *fMCEvent;               //! Ptr to MC event
  AliStack             *fStack;                 //! Ptr to stack
  // =======================================================================
  Float_t               fCentralityBin;         //  Centrality of current event  
  Int_t                 fNTracks;               //  N Tracks in the current event
  // =======================================================================
  THnSparseF           *fHnQA;                  //! THnSparseF : tracks for QA
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleQA, 1);
};

#endif
