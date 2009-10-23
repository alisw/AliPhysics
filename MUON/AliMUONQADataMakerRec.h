#ifndef ALIMUONQADATAMAKERREC_H
#define ALIMUONQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONQADataMakerRec
/// \brief MUON Quality assurance data maker
///

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliMUONVQADataMakerRec;

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec(Bool_t tracker=kTRUE, Bool_t trigger=kTRUE);         
  virtual ~AliMUONQADataMakerRec();

  AliMUONVQADataMakerRec* Tracker() const { return fTracker; }
  AliMUONVQADataMakerRec* Trigger() const { return fTrigger; }
  
  virtual void InitDigits(); 
  virtual void InitESDs(); 
  virtual void InitRaws(); 
  virtual void InitRecPoints(); 

  virtual void StartOfDetectorCycle(); 
  
  void MakeDigits();
  
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list);

private:
  
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);

  AliMUONVQADataMakerRec* fTracker; ///< tracker sub-qadatamaker
  AliMUONVQADataMakerRec* fTrigger; ///< trigger sub-qadatamaker
  
  ClassDef(AliMUONQADataMakerRec,10)  // MUON Quality assurance data maker

};
#endif
