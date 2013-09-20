/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/****************************************
 * analysis task for flow analysis with *
 *     multi-particle correlations      * 
 *                                      * 
 * author: Ante Bilandzic               *
 *         (abilandzic@gmail.com)       * 
 ***************************************/

#ifndef ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H
#define ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowAnalysisWithMultiparticleCorrelations.h"
#include "AliFlowEventSimple.h"

//================================================================================================================

class AliAnalysisTaskMultiparticleCorrelations : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskMultiparticleCorrelations();
  AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskMultiparticleCorrelations(){}; 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // ...

 private:
  AliAnalysisTaskMultiparticleCorrelations(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  AliAnalysisTaskMultiparticleCorrelations& operator=(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  
  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithMultiparticleCorrelations *fMPC; // "multi-particle correlations" object
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  
  ClassDef(AliAnalysisTaskMultiparticleCorrelations,0); 

};

//================================================================================================================

#endif











