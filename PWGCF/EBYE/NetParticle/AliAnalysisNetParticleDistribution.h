//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEDISTRIBUTION_H
#define ALIANALYSISNETPARTICLEDISTRIBUTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Helper Class for for NetParticle Distributions
// Authors: Jochen Thaeder <jochen@thaeder.de>

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

class AliAnalysisNetParticleDistribution : public TNamed {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleDistribution();
  virtual ~AliAnalysisNetParticleDistribution();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize */
  Int_t Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts, Bool_t isMC, Float_t *ptRange, Float_t etaMax, Int_t trackCutBit);

  /** Add histograms to outlist */
  void CreateHistograms(TList *outList);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent);

  /** Resre Event */
  void ResetEvent();

  /** Process NetParticle Distributions */ 
  Int_t Process();

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleDistribution(const AliAnalysisNetParticleDistribution&); // not implemented
  AliAnalysisNetParticleDistribution& operator=(const AliAnalysisNetParticleDistribution&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                           Process - Private
   * ---------------------------------------------------------------------------------
   */
   /** Process ESD tracks and fill histograms */
  Int_t ProcessESDTracks();

   /** Process AOD tracks and fill histograms */
  Int_t ProcessAODTracks();

  /** Process primary particles from the stack and fill histograms */
  Int_t ProcessStackParticles();

  /** Process control particles from the stack and fill histograms */
  Int_t ProcessStackControlParticles();
  
 /*
   * ---------------------------------------------------------------------------------
   *                            Helper Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Add set of histograms */
  void AddHistSet(const Char_t *name, const Char_t *title);

  /** Fill set of histograms */
  void FillHistSet(const Char_t *name, Float_t *np, Int_t controlIdx = -1);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleHelper *fHelper;        //! Ptr to helper class

  TList                *fOutList;               //! Output data container
  // -----------------------------------------------------------------------
  AliESDInputHandler   *fESDHandler;            //! Ptr to ESD handler 
  AliPIDResponse       *fPIDResponse;           //! Ptr to PID response Object
  AliESDEvent          *fESD;                   //! Ptr to ESD event

  AliAODInputHandler   *fAODHandler;            //! Ptr to AOD handler 
  AliAODEvent          *fAOD;                   //! Ptr to AOD event

  Bool_t                fIsMC;                  //  Is MC event

  AliMCEvent           *fMCEvent;               //! Ptr to MC event
  AliStack             *fStack;                 //! Ptr to stack

  AliESDtrackCuts      *fESDTrackCuts;          //! ESD cuts  
  // -----------------------------------------------------------------------
  Float_t               fEtaMax;                //  Max, absolut eta
  Float_t              *fPtRange;               //  Array of pt [min,max]

  Int_t                 fAODtrackCutBit;        //  Track filter bit for AOD tracks
  // -----------------------------------------------------------------------
  Float_t              *fNp;                    //  Array of particle/anti-particle counts

  Int_t                 fNCorrNp;               //  N sets of arrays of corrected particle/anti-particle counts
  Float_t             **fCorrNp;                //  Array of corrected particle/anti-particle counts

  Int_t                 fNMCNp;                 //  N sets of arrays of MC particle/anti-particle counts
  Float_t             **fMCNp;                  //  Array of MC particle/anti-particle counts

  Int_t                 fNControlMCNp;          //  N sets of arrays of control MC particle/anti-particle counts
  Float_t             **fControlMCNp;           //  Array of control MC particle/anti-particle counts
  // -----------------------------------------------------------------------
  THnSparseF           *fHnTrackUnCorr;         //  THnSparseF : uncorrected probe particles
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleDistribution, 1);
};

#endif
