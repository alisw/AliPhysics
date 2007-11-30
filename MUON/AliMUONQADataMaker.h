#ifndef ALIMUONQADATAMAKER_H
#define ALIMUONQADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup rec
/// \class AliMUONQADataMaker
/// \brief MUON Quality assurance data maker
///
//  Author Christian Finck

// dummy function for simulation part
// to avoid circular dependencie

// --- ROOT system ---
class TList; 

// --- AliRoot header files ---
class AliMUONVClusterStore;
class AliMUONVTrackStore;

#include "AliQADataMaker.h"

class AliMUONQADataMaker: public AliQADataMaker {

public:
  AliMUONQADataMaker();         
  AliMUONQADataMaker(const AliMUONQADataMaker& qadm);   
  AliMUONQADataMaker& operator=(const AliMUONQADataMaker& qadm);
  virtual ~AliMUONQADataMaker();
  
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
  
  virtual void   InitRaws(); 
  virtual void   InitRecPoints(); 
  virtual void   InitESDs(); 
  
  virtual void   MakeRaws(AliRawReader* rawReader); 
  virtual void   MakeRecPoints(TTree* recpo); 
  virtual void   MakeESDs(AliESDEvent* esd) ;
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TList* list);

  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store

  ClassDef(AliMUONQADataMaker,1)  // MUON Quality assurance data maker

};
#endif
