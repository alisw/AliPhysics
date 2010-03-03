#ifndef __ALIJETCORRELMIXER_H__
#define __ALIJETCORRELMIXER_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// Event mixing class. A 6-dimensinal pool fPool is maintained:
// type(fixed), type_idx(fixed), vertex(fixed), centrality(fixed), 
// and event(dynamic), particle(dynamic)
// fixed dimensions are fixed at task initialization time (via AliJetCorrelSelector)
// event&particle are allow to float during runtime (via linked lists TList=event & CorrelList_t=particle)
//-- Author: Paul Constantin

#include  "AliJetCorrelSelector.h"
#include  "AliJetCorrelWriter.h"

class AliJetCorrelMixer : public TObject {
  
 public:
  AliJetCorrelMixer();
  ~AliJetCorrelMixer();
  void Init(AliJetCorrelSelector * const s, AliJetCorrelMaker * const m, AliJetCorrelWriter * const w);
  
  // pool manipulation:
  void FillPool(CorrelList_t *partList, UInt_t pIdx, UInt_t vBin, UInt_t cBin);
  UInt_t PoolSize(PoolType_t pType, UInt_t vBin, UInt_t cBin);
  void CleanPool(PoolType_t pType);
  // mixing methods:  
  void Mix(UInt_t vBin, UInt_t cBin, UInt_t it, UInt_t ia, UInt_t ic);
  // print methods:
  void ShowPool(PoolType_t pType, UInt_t vBin, UInt_t cBin);
  void ShowSummary(PoolType_t pType, UInt_t pIdx, UInt_t vBin, UInt_t cBin);
  
 private:
  AliJetCorrelSelector* fSelector;       // user selection object
  AliJetCorrelMaker* fMaker;             // correlation maker object
  AliJetCorrelWriter* fWriter;           // output writer object
  CorrelList_t *TriggEvnt, *AssocEvnt;   // particle lists
  CorrelListIter_t AssocIter, TriggIter; // particle list iterators
  UInt_t fNumCentBins, fNumVertBins, fPoolDepth, fNumCorrel, fNumTriggs, fNumAssocs; // counters
  TList* fPool[2][kMAXNUMCORREL][kMAXVERTBIN][kMAXCENTBIN]; // the particle pools used for mixing
  
  // disable (make private) copy constructor and assignment operator:
  AliJetCorrelMixer(const AliJetCorrelMixer&);
  AliJetCorrelMixer* operator=(const AliJetCorrelMixer&);
  
  ClassDef(AliJetCorrelMixer, 1);
};

#endif

