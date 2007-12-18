#ifndef AliMUONQADataMakerSim_H
#define AliMUONQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup sim
/// \class AliMUONQADataMakerSim
/// \brief MUON Quality assurance data maker
///
//  Author Christian Finck

// dummy function for simulation part
// to avoid circular dependencie

// --- ROOT system ---
class TObjArray; 

// --- AliRoot header files ---
class AliMUONVClusterStore;
class AliMUONVTrackStore;

#include "AliQADataMakerSim.h"

class AliMUONQADataMakerSim: public AliQADataMakerSim {

public:
  AliMUONQADataMakerSim();         
  AliMUONQADataMakerSim(const AliMUONQADataMakerSim& qadm);   
  AliMUONQADataMakerSim& operator=(const AliMUONQADataMakerSim& qadm);
  virtual ~AliMUONQADataMakerSim();
  
private:
  virtual void   StartOfDetectorCycle(); 
      /// init hits QA from Array (not implemented)
  virtual void   InitHits() {return;} 
     /// init SDigits QA from Array (not implemented)
  virtual void   InitSDigits() {return;}  
     /// init SDigits QA from Array (not implemented)
  virtual void   InitDigits() {return;} 
    
    /// make hits QA from Array (not implemented)
  virtual void   MakeHits(TClonesArray* /*hits*/) {return;}
    /// make hits QA from tree (not implemented)
  virtual void   MakeHits(TTree* /*hits*/)        {return;}
    /// make SDigits QA from Array (not implemented)
  virtual void   MakeSDigits(TClonesArray* /*sigits*/) {return;} 
    /// make SDigits QA from Tree (not implemented)
  virtual void   MakeSDigits(TTree* /*sigits*/)        {return;} 
     /// make Digits QA from Array (not implemented)
  virtual void   MakeDigits(TClonesArray* /*digits*/)  {return;}
      /// make SDigits QA from Tree (not implemented)
  virtual void   MakeDigits(TTree* /*digits*/)         {return;}
  
    virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray* list);

  ClassDef(AliMUONQADataMakerSim,1)  // MUON Quality assurance data maker

};
#endif
