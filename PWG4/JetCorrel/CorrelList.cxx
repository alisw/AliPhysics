/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

//__________________________________________
// The particle list class CorrelList_t: implements a single-linked list
// Supporting classes:
// list node defined by class CorrelListNode_t and 
// list iterator defined by class CorrelListIter_t
//-- Author: Paul Constantin

#include "CorrelList.h"

using namespace std;

CorrelListNode_t::CorrelListNode_t() : fPart(NULL), fNext(NULL) {
}

CorrelListNode_t::CorrelListNode_t(CorrelParticle_t* p, CorrelListNode_t* n) : fPart(p), fNext(n) {
}

CorrelListNode_t::CorrelListNode_t(const CorrelListNode_t& rhs) : fPart(rhs.fPart), fNext(rhs.fNext) {
}

CorrelListNode_t::~CorrelListNode_t() {
  if(fPart) delete fPart;
  fNext = NULL;
}

const CorrelListNode_t& CorrelListNode_t::operator=(const CorrelListNode_t& rhs){
  fPart = rhs.fPart;
  fNext = rhs.fNext;
  return *this;
}

CorrelListIter_t::CorrelListIter_t() : fCurr(NULL) {
}

CorrelListIter_t::CorrelListIter_t(CorrelListNode_t* theNode) : fCurr(theNode) {
}

CorrelListIter_t::CorrelListIter_t(const CorrelListIter_t& rhs) : fCurr(rhs.fCurr) {
}

CorrelListIter_t::~CorrelListIter_t(){
}

void CorrelListIter_t::Check() const {
  // performs bounds check
  if(HasEnded()){
    std::cerr<<"CorrelListIter_t::Check() - ERROR: attempt to access null iterator!"<<std::endl; 
    exit(-1);
  }
}  

const CorrelListIter_t& CorrelListIter_t::operator=(const CorrelListIter_t& rhs){
  fCurr = rhs.fCurr;
  return *this;
}

CorrelList_t::CorrelList_t() : 
  fSize(0), fEvtID(0), fFilled(kFALSE), fPartID(t_unknown), fHead(NULL) {
  // constructor
}

const CorrelList_t& CorrelList_t::operator=(const CorrelList_t& rhs){
  fSize   = rhs.Size();
  fEvtID  = rhs.EvtID();
  fFilled = rhs.Filled();
  fPartID = rhs.PartID();
  fHead   = rhs.Head().Node();
  return *this;
}

void CorrelList_t::Push(CorrelParticle_t* p){
  // Push method creates new node, fills it with data object, hangs list to current fHead
  // (NULL on first insertion), then moves the fHead to this new node.

  if(!p){std::cerr<<"CorrelList_t::Push() - ERROR: cannot push undefined object!"<<std::endl; exit(-1);}
  CorrelListNode_t* newNode = new CorrelListNode_t(p,fHead);
  fHead = newNode;
  fSize++;
}

CorrelList_t* CorrelList_t::DeepCopy(){
  // returns deep copy of caller list
  // use it to store lists in memory for mixing pools 
  
  CorrelList_t *copy = new CorrelList_t;
  copy->Label(this->PartID(), this->EvtID());
  copy->SetFilled(this->Filled());
  // fHead and fSize are set by Push() below

  CorrelListIter_t iter = this->Head();
  while(!iter.HasEnded()){
    CorrelParticle_t*   gener = iter.Data();
    CorrelTrack_t*      track = dynamic_cast<CorrelTrack_t*>(gener);
    CorrelRecoParent_t* recon = dynamic_cast<CorrelRecoParent_t*>(gener);
    if(track)      copy->Push(track->Copy());
    else if(recon) copy->Push(recon->Copy());
    else           copy->Push(gener->Copy());
    iter.Move();
  } // iterator loop over caller list

  return copy;
}

void CorrelList_t::Reset(){
  // deep delete
  if(fSize>0){
    CorrelListIter_t iter = Head();
    while(!iter.HasEnded()){
      CorrelListIter_t current = iter; // get current node
      iter.Move();                     // move iterator to next node
      if(current.Node())               // redundant deletes happen in same part corr (dihadron)
	delete current.Node();         // with overlapping pT bins
      fSize--;
    }
  }
  fEvtID = 0;
  fFilled = kFALSE;
  fPartID = t_unknown;
  fHead = NULL;
}

void CorrelList_t::ShowHead() const {
  // top printout
  std::cout<<" CorrelList_t("<<this<<") head="<<fHead<<" size="<<Size()<<" filled="<<Filled()<<" evt="<<EvtID()<<" part="<<PartID()<<std::endl;
}

void CorrelList_t::Show() const {
  // full printout
  ShowHead();
  CorrelListIter_t iter = Head();
  while(!iter.HasEnded()){
    std::cout<<"At node("<<iter.Node()<<")";
    CorrelParticle_t* part = iter.Data();
    std::cout<<" is particle("<<part<<"):";
    if(part) part->Show(); else std::cout<<std::endl;
    iter.Move();
  }
}
