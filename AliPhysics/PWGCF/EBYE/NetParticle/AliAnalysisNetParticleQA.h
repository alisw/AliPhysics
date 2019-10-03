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

#include "AliAnalysisNetParticleBase.h"

class AliAnalysisNetParticleQA : public AliAnalysisNetParticleBase {

 public:

  AliAnalysisNetParticleQA();
  virtual ~AliAnalysisNetParticleQA();

  /*
   * ---------------------------------------------------------------------------------
   *                                 Public Methods
   * ---------------------------------------------------------------------------------
   */

  /** Process Event - implements purely virtual method */
  virtual void Process();

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
   *                                Methods - private
   * ---------------------------------------------------------------------------------
   */

  /** Create the efficiency / contamination THnSparseD  - implements virtual method */
  virtual void CreateHistograms();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // =======================================================================
  THnSparseF           *fHnQA;                  //! THnSparseF : tracks for QA
  // -----------------------------------------------------------------------

  ClassDef(AliAnalysisNetParticleQA, 1);
};

#endif
