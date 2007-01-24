//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// This class implements the use of a map of integers.              //
// For simplicity, this version is just a sorted linked list,       //
// but it may be rewritten later for better performance if required.//
//////////////////////////////////////////////////////////////////////  

#include "AliITSIntMap.h"
#include "AliITSIntMapNode.h"

AliITSIntMap::AliITSIntMap():
  fNrEntries(0),
  fFirst(0)
{}

AliITSIntMap::AliITSIntMap(const AliITSIntMap& imap):
  fNrEntries(0),
  fFirst(0)
{
  // copy constructor
  for (Int_t index=0; index<imap.fNrEntries; index++) {
    this->Insert(imap.GetKey(index),imap.GetVal(index));
  }
}

AliITSIntMap::~AliITSIntMap() 
{}

AliITSIntMap& AliITSIntMap::operator=(const AliITSIntMap& imap) {
  // assignment operator
  if (this!=&imap) {
    this->Clear();
    for (Int_t index=0; index<imap.fNrEntries; index++) {
      this->Insert(imap.GetKey(index),imap.GetVal(index));
    }
  }
  return *this;
}

void AliITSIntMap::Clear() {
  // clear the whole map
  while (fNrEntries>0) {
    Remove(fFirst->GetKey());
  }
}

Bool_t AliITSIntMap::Insert(Int_t key, Int_t val) {
  // insert a new node into the map (returns true if the node was not present before)
  if (fNrEntries==0) {
    fFirst = new AliITSIntMapNode(key, val, NULL);
    fNrEntries++;
    return kTRUE;
  }
  else {
    AliITSIntMapNode* it1 = fFirst;
    AliITSIntMapNode* it2 = fFirst->GetNext();
    while (it2!=NULL && key > it2->GetKey()) {
      it1 = it2;
      it2 = it1->GetNext();
    }
    if (key < it1->GetKey()) {
      fFirst = new AliITSIntMapNode(key, val, it1);
      fNrEntries++;
      return kTRUE;
    }
    else if (key > it1->GetKey() && 
	     ( (it2!=NULL && key < it2->GetKey()) || (it2==NULL) )
	     ) {
      it1->SetNext(new AliITSIntMapNode(key, val, it2));
      fNrEntries++;
      return kTRUE;
    }
  }
  return kFALSE;
}

Bool_t AliITSIntMap::Remove(Int_t key) {
  // remove a node from the map (returns true if the node was found)
  if (fNrEntries>0) {
    AliITSIntMapNode* it1 = fFirst;
    AliITSIntMapNode* it2 = fFirst->GetNext();
    while (it2!=NULL && key > it2->GetKey()) {
      it1 = it2;
      it2 = it1->GetNext();
    }
    if (key == it1->GetKey()) {
      fFirst = it1->GetNext();
      delete it1;
      fNrEntries--;
      return kTRUE;
    }
    else if (it2 != NULL && key == it2->GetKey()) {
      it1->SetNext(it2->GetNext());
      delete it2;
      fNrEntries--;
      return kTRUE;
    }
  }
  return kFALSE; 
}

AliITSIntMapNode* AliITSIntMap::Find(Int_t key) const {
  // finds a node and returns it (returns NULL if not found)
  if (fNrEntries>0) {
    AliITSIntMapNode* it1 = fFirst;
    AliITSIntMapNode* it2 = fFirst->GetNext();
    while (it2!=NULL && key > it2->GetKey()) {
      it1 = it2;
      it2 = it1->GetNext();
    }
    if (key == it1->GetKey()) {
      return it1;
    }
    else if (it2 != NULL && key == it2->GetKey()) {
      return it2;
    }
  }
  return NULL; 
}

Int_t AliITSIntMap::GetKey(UInt_t index) const {
  // returns the key of the node at position 'index' in the map
  // returns -1 if out of bounds
  if (index < fNrEntries) {
    AliITSIntMapNode* it1 = fFirst;
    for (UInt_t i=0; i<index; i++) {
      it1 = it1->GetNext();
    }
    return it1->GetKey();
  }
  else {
    return -1;
  }
}

Int_t AliITSIntMap::GetVal(UInt_t index) const {
  // returns the value of the node at position 'index' in the map
  // returns -1 if out of bounds
  if (index < fNrEntries) {
    AliITSIntMapNode* it1 = fFirst;
    for (UInt_t i=0; i<index; i++) {
      it1 = it1->GetNext();
    }
    return it1->GetVal();
  }
  else {
    return -1;
  }
}

void AliITSIntMap::PrintEntries() const {
  // prints all the entries (key,value pairs) of the map
  printf("*** Map Entries: (key , value)***\n");
  AliITSIntMapNode* it1 = fFirst;  
  for (UInt_t i=0; i<fNrEntries; i++) {
    printf("%d , %d\n",it1->GetKey(),it1->GetVal());
    it1 = it1->GetNext();
  }
  printf("*********************************\n");
}
