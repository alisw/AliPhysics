#ifndef ALIANALYSISTASKJETCORREL_H
#define ALIANALYSISTASKJETCORREL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// Main class for two-particle correlations.
// Calls AliJetCorrelSelector and AliJetCorrelMaker for setup, then
// AliJetCorrelReader for ESD/AOD input reading into CorrelList_t lists, then
// AliJetCorrelMixer for event mixing and AliJetCorrelWriter for output histos
//-- Author: Paul Constantin

#include "AliJetCorrelReader.h"
#include "AliJetCorrelMixer.h"

class AliAnalysisTaskJetCorrel : public AliAnalysisTaskSE {
  
 public:
  AliAnalysisTaskJetCorrel();
  AliAnalysisTaskJetCorrel(AliJetCorrelSelector* s);
  virtual ~AliAnalysisTaskJetCorrel();
  
  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
 private:
  AliESDEvent *fjcESD;                               //! ESD event
  TList *fOutputContainer;                           // Histogram container
  AliJetCorrelSelector *fSelector;                   //! User selection object
  UInt_t fNumCorrel, fNumTrigg, fNumAssoc, fNumEvts; //! counters
  AliJetCorrelMaker *fMaker;                         //! Correlation maker object
  AliJetCorrelWriter *fWriter;                       //! Output writer object
  AliJetCorrelReader *fReader;                       //! Input reader object
  AliJetCorrelMixer *fMixer;                         //! Event mixing object
  CorrelList_t *fTriggList, *fAssocList;             //! Trigger&Associated particle lists
  
  void CrossCorrelate(CorrelList_t * const TriggList, CorrelList_t * const AssocList,
		      UInt_t cBin, UInt_t vBin, UInt_t iCor);
  
  // disable (make private) copy constructor and assignment operator:
  AliAnalysisTaskJetCorrel(const AliAnalysisTaskJetCorrel&);
  AliAnalysisTaskJetCorrel& operator=(const AliAnalysisTaskJetCorrel&);
  
  ClassDef(AliAnalysisTaskJetCorrel, 1);
};
 
#endif
