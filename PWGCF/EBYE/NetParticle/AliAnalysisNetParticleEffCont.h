//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEEFFCONT_H
#define ALIANALYSISNETPARTICLEEFFCONT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/**
 * Class for NetParticle Distributions
 * -- Efficiency and contaminations for netParticle distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

class AliVTrack;

#include "THnSparse.h"

#include "AliAnalysisNetParticleBase.h"

class AliAnalysisNetParticleEffCont: public AliAnalysisNetParticleBase {

 public:

  AliAnalysisNetParticleEffCont();
  virtual ~AliAnalysisNetParticleEffCont();
  
  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Process Event - implements purely virtual method */
  virtual void Process();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Setter/Getter
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

  /** Event-wise Initialization - Can be implemented by every class */
  virtual void Init();

  /** Create the efficiency / contamination THnSparse */
  virtual void CreateHistograms();

  /** Event-wise Reset - Can be implemented by every class */
  virtual void Reset();

  /** Event-wise Setup - Can be implemented by every class */
  virtual Int_t Setup();

  // -----------------------------------------------------------------------

  /** Fill MC labels */
  void FillMCLabels(); 

  /** Fill efficiency THnSparse */
  void FillMCEffHist();

  /** Check if particle is contamination */
  void CheckContTrack(AliVTrack* track);
      
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // =======================================================================
  Int_t             **fLabelsRec;             //! 2x nTracks large array with labels for MC particles
  // =======================================================================
  THnSparseF         *fHnEff;                 //  THnSparseF efficiency 
  THnSparseF         *fHnCont;                //  THnSparseF contamination
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleEffCont, 1);
};

#endif
