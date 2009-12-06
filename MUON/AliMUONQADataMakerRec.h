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

  /// Return tracker sub-qadatamaker
  AliMUONVQADataMakerRec* Tracker() const { return fTracker; }
  /// Return trigger sub-qadatamaker
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

  virtual void ResetDetector(AliQAv1::TASKINDEX_t task);
  
private:
  /// Not implemented
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  /// Not implemented
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);

  AliMUONVQADataMakerRec* fTracker; ///< tracker sub-qadatamaker
  AliMUONVQADataMakerRec* fTrigger; ///< trigger sub-qadatamaker
  
  ClassDef(AliMUONQADataMakerRec,10)  // MUON Quality assurance data maker

};
#endif
