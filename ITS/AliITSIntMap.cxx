//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// This class implements the use of a map of integers.              //
// The values are kept in a binary tree, which is automatically     //
// reordered to be more balanced when the tree height gets too large//
//////////////////////////////////////////////////////////////////////  

#include "AliITSIntMap.h"
#include "AliITSIntMapNode.h"
#include <TMath.h>
#include <TError.h>

AliITSIntMap::AliITSIntMap():
  fNrEntries(0),
  fRoot(NULL),
  fFastAccess(kFALSE),
  fFastAccessSerialize(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{}

AliITSIntMap::AliITSIntMap(AliITSIntMapNode* root, UInt_t nrEntries):
  fNrEntries(nrEntries),
  fRoot(root),
  fFastAccess(kFALSE),
  fFastAccessSerialize(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{}

AliITSIntMap::AliITSIntMap(const AliITSIntMap& imap):
  fNrEntries(0),
  fRoot(NULL),
  fFastAccess(kFALSE),
  fFastAccessSerialize(kFALSE),
  fFastAccessArray(NULL),
  fDummyIndex(0)
{
  // copy constructor
  *this = imap;
}

AliITSIntMap::~AliITSIntMap() {
  Clear();
}

AliITSIntMap& AliITSIntMap::operator=(const AliITSIntMap& imap) {
  // assignment operator
  if (this!=&imap) {
    this->Clear();
    fRoot = CloneNode(imap.fRoot);
    fFastAccess=kFALSE;
    fFastAccessSerialize=kFALSE;
    fFastAccessArray=NULL;
    fDummyIndex=0;
  }
  return *this;
}

void AliITSIntMap::Clear() {
  // clear the whole map
  ClearFastAccess();
  ClearNode(fRoot);
}

void AliITSIntMap::ClearNode(AliITSIntMapNode* &node) {
  // clear this node and all children nodes
  if (node==NULL) return;
  ClearNode(node->Left());
  ClearNode(node->Right());
  delete node;
  fNrEntries--;
  node = NULL;
  fFastAccess=kFALSE;
  fFastAccessSerialize=kFALSE;
}

AliITSIntMap* AliITSIntMap::Clone() const {
  // returns a clone of the map
  AliITSIntMapNode* newRoot;
  newRoot = CloneNode(fRoot);
  AliITSIntMap* newMap = new AliITSIntMap(newRoot,fNrEntries);
  return newMap;
}

AliITSIntMapNode* AliITSIntMap::CloneNode(AliITSIntMapNode* node) const {
  if (node==NULL) return NULL;
  else return new AliITSIntMapNode(node->Key(),node->Val(),CloneNode(node->Left()),CloneNode(node->Right()));
}

Bool_t AliITSIntMap::Insert(Int_t key, Int_t val) {
  // insert a new node into the map (returns true if the node was not present before)
  UInt_t entriesBefore = fNrEntries;
  InsertNode(key,val,fRoot,0);
  if (fNrEntries>entriesBefore) return kTRUE;
  else return kFALSE;
}

void AliITSIntMap::InsertNode(Int_t key, Int_t val, AliITSIntMapNode* &node, UInt_t height) {
  // method to insert a node in the tree (used recursively)
  height++;
  if (node==NULL) {
    node = new AliITSIntMapNode(key,val,NULL,NULL);
    fNrEntries++;
    fFastAccess=kFALSE;
    fFastAccessSerialize=kFALSE;
    UInt_t balanceHeight = (UInt_t) (log(fNrEntries+1)/log(2)+1);
    if ( (height-balanceHeight)*(height-balanceHeight) > fNrEntries ) {
      Balance();
    }
  }
  else if (key < node->Key()) {
    InsertNode(key,val,node->Left(),height);
  }
  else if (key > node->Key()) {
    InsertNode(key,val,node->Right(),height);
  }
  else { // (key==node->Key()): do nothing (avoid duplicates)
    //    Warning("AliITSIntMap::InsertNode","Node with key %d already in map. Not inserted.",key);
  }
}

Bool_t AliITSIntMap::Remove(Int_t key) {
  // remove a node from the map (returns true if the node was found)
  UInt_t entriesBefore = fNrEntries;
  RemoveNode(key,fRoot);
  if (fNrEntries<entriesBefore) return kTRUE;
  else return kFALSE;
}

void AliITSIntMap::RemoveNode(Int_t key, AliITSIntMapNode* &node) {
  // method to remove a node in the tree (used recursively)
  if (node == NULL) return;   // node not found; do nothing
  if (key < node->Key()) {
    RemoveNode(key,node->Left());
  }
  else if (key > node->Key()) {
    RemoveNode(key,node->Right());
  }
  else if (node->Left()!=NULL && node->Right()!=NULL) { // Two children
    if (fNrEntries%2==0) { // for better balance, remove from left or right sub tree
      AliITSIntMapNode* moveNode = FindMinNode(node->Right());
      node->SetKey(moveNode->Key());
      node->SetVal(moveNode->Val());
      RemoveNode(moveNode->Key(),node->Right());
    }
    else {
      AliITSIntMapNode* moveNode = FindMaxNode(node->Left());
      node->SetKey(moveNode->Key());
      node->SetVal(moveNode->Val());
      RemoveNode(moveNode->Key(),node->Left());
    }
  }
  else {
    AliITSIntMapNode* oldNode = node;
    node = (node->Left()!=NULL) ? node->Left() : node->Right();
    fNrEntries--;
    delete oldNode;
    fFastAccess=kFALSE;
    fFastAccessSerialize=kFALSE;
  }
}

Bool_t AliITSIntMap::Pop(Int_t& key, Int_t& val) {
  // removes one entry (root) from tree, giving its key,val pair
  if (fRoot!=NULL) {
    key = fRoot->Key();
    val = fRoot->Val();
    return Remove(key);
  }
  else return kFALSE;
}

AliITSIntMapNode* AliITSIntMap::FindMinNode(AliITSIntMapNode* node) const {
  // returns the node with smallest key in the sub tree starting from node
  if (node==NULL) return NULL;
  else if (node->Left()==NULL) return node;
  else return FindMinNode(node->Left());
}

AliITSIntMapNode* AliITSIntMap::FindMaxNode(AliITSIntMapNode* node) const {
  // returns the node with largest key in the sub tree starting from node
  if (node==NULL) return NULL;
  else if (node->Right()==NULL) return node;
  else return FindMaxNode(node->Right());
}

AliITSIntMapNode* AliITSIntMap::Find(Int_t key) const {
  // finds a node and returns it (returns NULL if not found)
  return FindNode(key,fRoot,0);
}

AliITSIntMapNode*  AliITSIntMap::FindNode(Int_t key, AliITSIntMapNode* node, UInt_t height) const {
  // method to find a node in the tree (used recursively)
  if (node==NULL) return NULL;
  height++;
  if (key<node->Key()) return FindNode(key,node->Left(),height);
  else if (key>node->Key()) return FindNode(key,node->Right(),height);
  else { // Match
//    //*** balance if height too high. const above have to be removed if this is needed ***
//    UInt_t balanceHeight = (UInt_t) (log(fNrEntries+1)/log(2)+1);
//    if ( (height-balanceHeight)*(height-balanceHeight) > fNrEntries ) {
//      Balance();
//    }
    return node;
  }
}

Int_t AliITSIntMap::GetVal(Int_t key) const {
  // returns the value for the node with key
  AliITSIntMapNode* node = Find(key);
  if (node!=NULL) {
    return node->Val();
  }
  else {
    Warning("AliITSIntMap::GetVal","Node with key %d not found in map. Returning -999.",key);
    return -999;
  }
}

void AliITSIntMap::Balance() {
  // method to balance the tree
  //  printf("balance H=%d --> ",GetTreeHeight());
  if (fNrEntries==0) return;
  if (!fFastAccess) InitFastAccess();
  fRoot = BalanceNode(0,fNrEntries-1);
  //  printf("H=%d\n",GetTreeHeight());
}

AliITSIntMapNode* AliITSIntMap::BalanceNode(Int_t lowInd, Int_t highInd) {
  // balances the tree by selecting the center of an index range 
  // (used recursively)
  if (lowInd>highInd) return NULL;
  Int_t thisInd = lowInd+(highInd-lowInd)/2;
  fFastAccessArray[thisInd]->Left() = BalanceNode(lowInd,thisInd-1);
  fFastAccessArray[thisInd]->Right() = BalanceNode(thisInd+1,highInd);
  return fFastAccessArray[thisInd];
}

void AliITSIntMap::ClearFastAccess(){
  // clears the fast access array of pointers
  if (fFastAccessArray!=NULL) {
    delete [] fFastAccessArray;
    fFastAccessArray=NULL;
  }
  fFastAccess=kFALSE;
  fFastAccessSerialize=kFALSE;
}

void AliITSIntMap::InitFastAccess(){
  // initializes the fast access array
  if (fFastAccess) return;
  ClearFastAccess();
  if (fNrEntries>0) {
    fFastAccessArray = new AliITSIntMapNode*[fNrEntries];
    fDummyIndex=0;
    InitFastAccessNode(fRoot);
    fFastAccess=kTRUE;
  }
}

void AliITSIntMap::InitFastAccessNode(AliITSIntMapNode* node) {
  // initializes the fast access array starting from node (used recursively)
  if (node==NULL) return;
  InitFastAccessNode(node->Left());
  fFastAccessArray[fDummyIndex++] = node;
  InitFastAccessNode(node->Right());
}

void AliITSIntMap::InitFastAccessSerialize(){
  // initializes the fast access array
  if (fFastAccessSerialize) return;
  ClearFastAccess();
  if (fNrEntries>0) {
    fFastAccessArray = new AliITSIntMapNode*[fNrEntries];
    fDummyIndex=0;
    InitFastAccessSerializeNode(fRoot);
    fFastAccessSerialize=kTRUE;
  }
}

void AliITSIntMap::InitFastAccessSerializeNode(AliITSIntMapNode* node) {
  // initializes the fast access array for tree ordering starting from node (used recursively)
  if (node==NULL) return;
  fFastAccessArray[fDummyIndex++] = node;
  InitFastAccessSerializeNode(node->Left());
  InitFastAccessSerializeNode(node->Right());
}

Int_t AliITSIntMap::GetKeyIndex(UInt_t index) {
  // returns the key of the node at position 'index' in the map
  // returns -1 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess && !fFastAccessSerialize) InitFastAccess();
    return fFastAccessArray[index]->Key();
  }
  return -1;
}

Int_t AliITSIntMap::GetValIndex(UInt_t index) {
  // returns the value of the node at position 'index' in the map
  // returns -1 if out of bounds
  if (index<fNrEntries) {
    if (!fFastAccess && !fFastAccessSerialize) InitFastAccess();
    return fFastAccessArray[index]->Val();
  }
  return -1;
}

AliITSIntMapNode* AliITSIntMap::FindNodeIndex(UInt_t index, AliITSIntMapNode* node) const {
  // method to find the index:th node in the tree (used recursively)
  // this method should not be needed anymore, since GetKeyIndex/GetValIndex is faster
  static UInt_t fTmpInd;
  if (node==fRoot) fTmpInd=0;
  if (node->Left()!=NULL) {
    AliITSIntMapNode* tmpResult = FindNodeIndex(index,node->Left());
    if (tmpResult != NULL) {
      return tmpResult;
    }
  }
  if (fTmpInd==index) return node;
  fTmpInd++;
  if (node->Right()!=NULL) {
    AliITSIntMapNode* tmpResult = FindNodeIndex(index,node->Right());
    if (tmpResult != NULL) {
      return tmpResult;
    }
  }
  return NULL;
}

void AliITSIntMap::PrintEntries() const {
  // prints all the entries (key,value pairs) of the map
  printf("*** Map Entries: (key , value)***\n");
  PrintNode(fRoot);
  printf("*********************************\n");
}

void AliITSIntMap::PrintNode(AliITSIntMapNode* node) const {
  // method to print node entry (key,value) (used recursively)
  if (node==NULL) return;
  if (node->Left()!=NULL) PrintNode(node->Left());
  printf("%d , %d\n",node->Key(),node->Val());
  if (node->Right()!=NULL) PrintNode(node->Right());
}

UInt_t AliITSIntMap::GetTreeHeight() const {
  // returns the height of the tree
  return GetTreeHeightNode(fRoot);
}

UInt_t AliITSIntMap::GetTreeHeightNode(AliITSIntMapNode* node) const {
  // returns tree height for the sub tree starting from node (used recursively)
  if (node==NULL) return 0;
  UInt_t leftH = GetTreeHeightNode(node->Left());
  UInt_t rightH = GetTreeHeightNode(node->Right());
  if (leftH>=rightH) return leftH+1;
  else               return rightH+1;
}
