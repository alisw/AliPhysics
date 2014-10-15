//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEDCA_H
#define ALIANALYSISNETPARTICLEDCA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * Class for NetParticle Distributions
 * -- DCA distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "THnSparse.h"

#include "AliAnalysisNetParticleBase.h"

class AliAnalysisNetParticleDCA : public AliAnalysisNetParticleBase {

 public:

  AliAnalysisNetParticleDCA();
  virtual ~AliAnalysisNetParticleDCA();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Process Event - implements purely virtual method */
  virtual void Process();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Setter / Getter
   * ---------------------------------------------------------------------------------
   */

  /** Set Background ESD Track Cuts */
  void SetESDTrackCutsBkg(AliESDtrackCuts *p) {fESDTrackCutsBkg = p;}

  /** Get Ptr to DCA THnSparse */
  THnSparseD* GetHnDCA() {return fHnDCA;}

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleDCA(const AliAnalysisNetParticleDCA&); // not implemented
  AliAnalysisNetParticleDCA& operator=(const AliAnalysisNetParticleDCA&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Create the efficiency / contamination THnSparseD  
   *    implements purely virtual method 
   */
  virtual void CreateHistograms();

  /** Get contamination index of track */
  Int_t GetContIdxTrack(Int_t label, Int_t sign);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // =======================================================================
  AliESDtrackCuts    *fESDTrackCutsBkg;       //! ESD cuts  
  // =======================================================================
  THnSparseD         *fHnDCA;                 //  THnSparseD contamination DCA
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleDCA, 1);
};

#endif
