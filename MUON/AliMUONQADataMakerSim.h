#ifndef ALIMUONQADATAMAKERSIM_H
#define ALIMUONQADATAMAKERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONQADataMakerSim
/// \brief MUON Quality assurance data maker
///
//  Author Christian Finck

// --- ROOT system ---
class TObjArray; 

// --- AliRoot header files ---
class AliMUONVHitStore;
class AliMUONVDigitStore;

#include "AliQADataMakerSim.h"

class AliMUONQADataMakerSim: public AliQADataMakerSim {

public:
  AliMUONQADataMakerSim();         
  AliMUONQADataMakerSim(const AliMUONQADataMakerSim& qadm);   
  AliMUONQADataMakerSim& operator=(const AliMUONQADataMakerSim& qadm);
  virtual ~AliMUONQADataMakerSim();
  
private:
  virtual void   StartOfDetectorCycle(); 
      /// init hits QA from Array 
  virtual void   InitHits(); 
     /// init SDigits QA from Array 
  virtual void   InitSDigits();  
     /// init SDigits QA from Array
  virtual void   InitDigits();
    
    /// make hits QA from Array (not implemented)
  virtual void   MakeHits(TClonesArray* /*hits*/) {return;}
    /// make hits QA from tree
  virtual void   MakeHits(TTree* hitsTree);
    /// make SDigits QA from Array (not implemented)
  virtual void   MakeSDigits(TClonesArray* /*sigits*/) {return;} 
    /// make SDigits QA from Tree
  virtual void   MakeSDigits(TTree* sigitsTree);
     /// make Digits QA from Array (not implemented)
  virtual void   MakeDigits(TClonesArray* /*digits*/)  {return;}
      /// make SDigits QA from Tree
  virtual void   MakeDigits(TTree* digitsTree);
  
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list);

  AliMUONVHitStore*   fHitStore;   //!< pointer to hit store
  AliMUONVDigitStore* fDigitStore; //!< pointer to digit store
                                    
  ClassDef(AliMUONQADataMakerSim,2)  // MUON Quality assurance data maker

};
#endif
