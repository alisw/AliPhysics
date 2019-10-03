//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEBASE_H
#define ALIANALYSISNETPARTICLEBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * Class for NetParticle Distributions
 * -- Base Class
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "THnSparse.h"

#include "AliAnalysisNetParticleHelper.h"

class AliESDEvent;
class AliESDInputHandler;
class AliMCEvent;
class AliStack;
class AliAODInputHandler;

class AliAnalysisNetParticleBase : public TNamed {

 public:

  AliAnalysisNetParticleBase();
  AliAnalysisNetParticleBase(const Char_t* name, const Char_t* title);
  virtual ~AliAnalysisNetParticleBase();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize */
  void Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts = NULL);

  /** Setup Event */
  Int_t SetupEvent();

  /** Reset Event */
  void ResetEvent();

  /** Process Event - To be implemented by every class */
  virtual void Process() = 0;

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleBase(const AliAnalysisNetParticleBase&); // not implemented
  AliAnalysisNetParticleBase& operator=(const AliAnalysisNetParticleBase&); // not implemented

 protected:
  /*
   * ---------------------------------------------------------------------------------
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Event-wise Initialization - Can be implemented by every class */
  virtual void Init() {};

  /** Create histograms - Can be implemented by every class */
  virtual void CreateHistograms() {};

  /** Event-wise Reset - Can be implemented by every class */
  virtual void Reset() {};

  /** Event-wise Setup - Can be implemented by every class */
  virtual Int_t Setup() { return 0;};

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleHelper *fHelper;      //! Ptr to helper class
  // =======================================================================
  Int_t               fPdgCode;               // PDG code of particle to be found 
  // =======================================================================
  AliESDEvent        *fESD;                   //! ESD object
  AliESDtrackCuts    *fESDTrackCuts;          //! ESD cuts  
  // -----------------------------------------------------------------------
  AliAODEvent        *fAOD;                   //! AOD object
  TClonesArray       *fArrayMC;               //! array of MC particles
  Int_t               fAODtrackCutBit;        //  Track filter bit for AOD tracks
  // -----------------------------------------------------------------------
  Bool_t              fIsMC;                  //  Is MC event
  AliMCEvent         *fMCEvent;               //! Ptr to MC event
  AliStack           *fStack;                 //! Ptr to stack
  // =======================================================================
  Float_t             fCentralityBin;         //  Centrality of current event  
  Int_t               fNTracks;               //  N Tracks in the current event
  // =======================================================================

  ClassDef(AliAnalysisNetParticleBase, 1);
};

#endif
