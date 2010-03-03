#ifndef __CORRELLIST_H__
#define __CORRELLIST_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// The particle list class CorrelList_t: implements a single-linked list
// Supporting classes:
// list node defined by class CorrelListNode_t and 
// list iterator defined by class CorrelListIter_t
//-- Author: Paul Constantin

#include "CorrelParticle.h"
#include "CorrelTrack.h"
#include "CorrelRecoParent.h"

class CorrelList_t;
class CorrelListIter_t;

class CorrelListNode_t {
 public:
  CorrelListNode_t();
  CorrelListNode_t(CorrelParticle_t* p, CorrelListNode_t* n);
  CorrelListNode_t(const CorrelListNode_t& rhs);
  ~CorrelListNode_t();
  const CorrelListNode_t& operator=(const CorrelListNode_t& rhs);
  CorrelParticle_t* GetData() const {return fPart;}
  CorrelListNode_t* GetNext() const {return fNext;}
    
 private:
    CorrelParticle_t* fPart;
    CorrelListNode_t* fNext;      
};

class CorrelListIter_t {
 public:
  CorrelListIter_t();
  CorrelListIter_t(CorrelListNode_t* theNode);
  CorrelListIter_t(const CorrelListIter_t& rhs);
  ~CorrelListIter_t();
  const CorrelListIter_t& operator=(const CorrelListIter_t& rhs);
    
  Bool_t HasEnded() const {return (fCurr==NULL);}
  void Check();
  void Move() {Check(); fCurr=fCurr->GetNext();}
  CorrelListNode_t* Node() {Check(); return fCurr;}
  CorrelParticle_t* Data() {Check(); return fCurr->GetData();}
  
 private:
  CorrelListNode_t* fCurr; // iterator "current node"
};

class CorrelList_t : public TObject {
 public:
  CorrelList_t();
  virtual ~CorrelList_t() {Reset();}
  void Reset();
  const CorrelList_t& operator=(const CorrelList_t& rhs); // makes shallow copy
  CorrelList_t* DeepCopy();                     // use this method to get a deep copy
  
  void Push(CorrelParticle_t* p);
  UInt_t Size() const {return fSize;}
  UInt_t EvtID() const {return fEvtID;}
  PartType_t PartID() const {return fPartID;}
  PoolType_t PoolID() const {return fPoolID;}
  Bool_t Filled() const {return fFilled;}
  void SetFilled(Bool_t f) {fFilled=f;}
  void Label(PartType_t p, PoolType_t l, UInt_t e) {fPartID=p; fPoolID=l; fEvtID=e;}
  CorrelListIter_t Head() const {return CorrelListIter_t(fHead);}
  void ShowHead();
  void Show();
  
 private:
  UInt_t fSize;             // list size
  UInt_t fEvtID;            // event ID
  Bool_t fFilled;           // is filled
  PartType_t fPartID;       // particle ID
  PoolType_t fPoolID;       // pool type
  CorrelListNode_t* fHead;  // list head
  
  CorrelList_t(const CorrelList_t& rhs);        // forbid copy constructor
};

#endif
