//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEEFFCONT_H
#define ALIANALYSISNETPARTICLEEFFCONT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/**
 * Class for for NetParticle Distributions
 * -- Efficiency and contaminations for netParticle distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "THnSparse.h"

#include "AliAnalysisNetParticleHelper.h"

class AliESDEvent;
class AliESDInputHandler;
class AliAODInputHandler;
class AliMCEvent;

class AliAnalysisNetParticleEffCont : public TNamed {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleEffCont();
  virtual ~AliAnalysisNetParticleEffCont();
  
  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Initialize */
  void Initialize(AliESDtrackCuts *cuts, AliAnalysisNetParticleHelper* helper, Int_t trackCutBit);

  /** Setup Event */
  Int_t SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent); 

  /** Reset Event */
  void ResetEvent();

  /** Process Event */
  void Process();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to efficiency THnSparse */
  THnSparseF* GetHnEff()  {return fHnEff;}

  /** Get Ptr to contaminiation THnSparse */
  THnSparseF* GetHnCont() {return fHnCont;}

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleEffCont(const AliAnalysisNetParticleEffCont&); // not implemented
  AliAnalysisNetParticleEffCont& operator=(const AliAnalysisNetParticleEffCont&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Create the efficiency / contamination THnSparseF */
  void CreateHistograms();

  /** Fill MC labels */
  void FillMCLabels(); 

  /** Fill efficiency THnSparse */
  void FillMCEffHist();
  void FillMCEffHistAOD();

  /** Check if particle is contamination */
  void CheckContTrack(Int_t label, Float_t sign, Int_t idxTrack);
  void CheckContTrackAOD(Int_t label, Float_t sign, Int_t idxTrack);
      
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  AliAnalysisNetParticleHelper *fHelper;      //! Ptr to helper class

  // -----------------------------------------------------------------------

  Int_t               fPdgCode;               // PDG code of particle to be found 

  // --- ESD only ----------------------------------------------------------

  AliESDEvent        *fESD;                   //! ESD object
  AliESDtrackCuts    *fESDTrackCuts;          //! ESD cuts  

  // --- AOD only ----------------------------------------------------------

  AliAODEvent        *fAOD;                   //! AOD object
  TClonesArray       *fArrayMC;               //! array of MC particles

  // -----------------------------------------------------------------------

  Float_t             fCentralityBin;         //  Centrality of current event  
  Int_t               fNTracks;               //  N Tracks in the current event
  
  Int_t               fAODtrackCutBit;        //  Track filter bit for AOD tracks

  // --- MC only -----------------------------------------------------------

  AliStack           *fStack;                 //! Ptr to stack
  AliMCEvent         *fMCEvent;               //! Ptr to MC event

  Int_t             **fLabelsRec;             //! 2x nTracks large array with labels for MC particles

  // -----------------------------------------------------------------------

  THnSparseF         *fHnEff;                 //  THnSparseF efficiency 
  THnSparseF         *fHnCont;                //  THnSparseF contamination

  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleEffCont, 1);
};

#endif
