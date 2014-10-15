//-*- Mode: C++ -*-

#ifndef ALIANALYSISNETPARTICLEDISTRIBUTION_H
#define ALIANALYSISNETPARTICLEDISTRIBUTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/**
 * Class for NetParticle Distributions
 * -- Create input for distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

#include "THnSparse.h"
#include "TList.h"

#include "AliAnalysisNetParticleBase.h"

class AliAnalysisNetParticleDistribution : public AliAnalysisNetParticleBase {

 public:

  AliAnalysisNetParticleDistribution();
  virtual ~AliAnalysisNetParticleDistribution();

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

  void SetOutList(TList* l) {fOutList = l;}

  ///////////////////////////////////////////////////////////////////////////////////

 private:

  AliAnalysisNetParticleDistribution(const AliAnalysisNetParticleDistribution&); // not implemented
  AliAnalysisNetParticleDistribution& operator=(const AliAnalysisNetParticleDistribution&); // not implemented

  /*
   * ---------------------------------------------------------------------------------
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Event-wise Initialization - implements virtual method */
  virtual void Init();

  /** Event-wise Reset - implements virtual method */
  virtual void Reset();

  /** HistSet-wise Reset */
 void ResetHistSet();

  /** Create the efficiency / contamination THnSparseD  - implements virtual method */
  virtual void CreateHistograms();

  // -----------------------------------------------------------------------

   /** Process ESD/AOD tracks and fill histograms */
  Int_t ProcessTracks();

  /** Process primary particles from the stack and fill histograms */
  Int_t ProcessParticles();
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Helper Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Add set of histograms */
  void AddHistSetCent(const Char_t *name, const Char_t *title);
  void AddHistSetCentPt(const Char_t *name, const Char_t *title);

  /** Fill set of histograms */
  void FillHistSetCent(const Char_t *name, Int_t idx, Bool_t isMC);
  void FillHistSetCentPt(const Char_t *name, Int_t idx, Bool_t isMC);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // =======================================================================
  TList                *fOutList;               //! Output data container
  // =======================================================================
  Int_t                 fOrder;                 //  Max order of higher order distributions
  // -----------------------------------------------------------------------
  Int_t                 fNNp;                   //  N sets of arrays of particle/anti-particle counts
  Int_t               **fNp;                    //  Array of particle/anti-particle counts
  Int_t              ***fNpPt;                  //  Array of particle/anti-particle per ptBin counts

  Int_t                 fNMCNp;                 //  N sets of arrays of MC particle/anti-particle counts
  Int_t               **fMCNp;                  //  Array of MC particle/anti-particle counts
  Int_t              ***fMCNpPt;                //  Array of MC particle/anti-particle per ptBin counts
  // -----------------------------------------------------------------------
  Double_t            **fRedFactp;              //  Array of particle/anti-particle reduced factorial
  // =======================================================================
  THnSparseD           *fHnTrackUnCorr;         //  THnSparseD : uncorrected probe particles
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleDistribution, 1);
};

#endif
