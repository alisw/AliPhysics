#ifndef __ALIJETCORRELREADER_H__
#define __ALIJETCORRELREADER_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//______________________________________________________________________________________
// Class for input (ESD or AOD) reading and filling of Trigger&Associated particle lists
//-- Author: Paul Constantin
 
#include "AliJetCorrelWriter.h"
#include "AliJetCorrelSelector.h"

class AliJetCorrelReader : public TObject {
  
 public:
  AliJetCorrelReader(); 
  ~AliJetCorrelReader();
  
  void Init(AliJetCorrelSelector * const s, AliJetCorrelWriter * const w);    
  void SetEvent(AliESDEvent * const e) {jcESD=e;}
  
  Float_t GetMultiplicity(); 
  Float_t GetVertex();
  Bool_t VtxOutPipe();
  void FillLists(CorrelList_t* list1, CorrelList_t* list2);
  
 private:    
  AliESDEvent *jcESD;              // input event (ESD/AOD)
  AliJetCorrelSelector *fSelector; // user selection object
  AliJetCorrelWriter *fWriter;     // output writer object
  
  void FillList(CorrelList_t* list);
  void FillESDTrackLists(CorrelList_t* list1,CorrelList_t* list2);
  void FillESDTrackList(CorrelList_t* list);
  void FillESDPhotonList(CorrelList_t* list);
  void FillESDDiphotonList(CorrelList_t* list);
  void FillESDDielectronList(CorrelList_t* list);
  void FillParentList(CorrelList_t* list1, CorrelList_t* list2);
  
  // disable (make private) copy constructor and assignment operator:
  AliJetCorrelReader(const AliJetCorrelReader&);
  AliJetCorrelReader& operator=(const AliJetCorrelReader&);
  
  ClassDef(AliJetCorrelReader, 1);
};

#endif
