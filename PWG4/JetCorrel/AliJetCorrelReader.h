#ifndef ALIJETCORRELREADER_H
#define ALIJETCORRELREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________________________________________
// Class for input (ESD or AOD) reading.
// At the moment only ESD input is really implemented, AOD to be added later.
// Its products are the Trigger&Associated particle lists
//-- Author: Paul Constantin
 
#include "AliJetCorrelWriter.h"

class AliJetCorrelReader : public TObject {
  
 public:
  AliJetCorrelReader(); 
  ~AliJetCorrelReader();
  
  void Init(AliJetCorrelSelector * const s, AliJetCorrelWriter * const w);    
  void SetEvent(AliESDEvent * const e) {fjcESD=e;}
  
  Float_t GetMultiplicity() const; 
  Float_t GetVertex() const;
  Bool_t VtxOutPipe() const;
  void FillLists(CorrelList_t* list1, CorrelList_t* list2);
  
 private:    
  AliESDEvent *fjcESD;              //! input event (ESD/AOD)
  AliJetCorrelSelector *fSelector;  //! user selection object
  AliJetCorrelWriter *fWriter;      //! output writer object
  
  void FillList(CorrelList_t* list, Bool_t isTrigg);
  void FillESDTrackLists(CorrelList_t* list1,CorrelList_t* list2);
  void FillESDTrackList(CorrelList_t* list, Bool_t isTrigg);
  void FillESDPhotonList(CorrelList_t* list, Bool_t isTrigg);
  void FillESDDiphotonList(CorrelList_t* list, Bool_t isTrigg);
  void FillESDDielectronList(CorrelList_t* list, Bool_t isTrigg);
  void FillParentList(CorrelList_t* list1, CorrelList_t* list2, Bool_t isTrigg);
  
  // disable (make private) copy constructor and assignment operator:
  AliJetCorrelReader(const AliJetCorrelReader&);
  AliJetCorrelReader& operator=(const AliJetCorrelReader&);
  
  ClassDef(AliJetCorrelReader, 1);
};

#endif
