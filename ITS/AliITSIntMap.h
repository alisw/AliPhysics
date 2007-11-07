#ifndef ALI_ITS_INTMAP_H
#define ALI_ITS_INTMAP_H

//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// This class implements the use of a map of integers.              //
//                                                                  //
//////////////////////////////////////////////////////////////////////  

#include <Rtypes.h>

class AliITSIntMapNode;

class AliITSIntMap {

 public:
  AliITSIntMap();
  AliITSIntMap(AliITSIntMapNode* root, UInt_t nrEntries);
  AliITSIntMap(const AliITSIntMap& imap);
  virtual ~AliITSIntMap();
  AliITSIntMap& operator=(const AliITSIntMap& imap);

  void               Clear();
  AliITSIntMap*      Clone() const;
  Bool_t             Insert(Int_t key, Int_t val);
  Bool_t             Remove(Int_t key);
  Bool_t             Pop(Int_t& key, Int_t& val);
  AliITSIntMapNode*  Find(Int_t key) const;
  Int_t              GetVal(Int_t key) const;

  UInt_t             GetNrEntries() const {return fNrEntries;}
  Int_t              GetKeyIndex(UInt_t index);
  Int_t              GetValIndex(UInt_t index);

  void               PrintEntries() const;
  void               Balance();
  void               PrepareSerialize() {InitFastAccessSerialize();}
  void               PrepareSerializeOrdered() {InitFastAccess();}
  UInt_t             GetTreeHeight() const;

 private:
  UInt_t             fNrEntries;          // nr of entries in map
  AliITSIntMapNode*  fRoot;               // link to first node of tree
  Bool_t             fFastAccess;         // is fast access array initialized (key ordered)?
  Bool_t             fFastAccessSerialize;// is fast access array initialized (tree ordered)?
  AliITSIntMapNode** fFastAccessArray;    // array of pointers to nodes
  UInt_t             fDummyIndex;         // dummy index used when traversing tree

  void               ClearNode(AliITSIntMapNode* &node); // delete this node and all below
  AliITSIntMapNode*  CloneNode(AliITSIntMapNode* node) const;
  void               InsertNode(Int_t key, Int_t val, AliITSIntMapNode* &node, UInt_t height);
  void               RemoveNode(Int_t key, AliITSIntMapNode* &node);
  AliITSIntMapNode*  FindNode(Int_t key, AliITSIntMapNode* node, UInt_t height) const;
  AliITSIntMapNode*  FindNodeIndex(UInt_t index, AliITSIntMapNode* node) const;
  AliITSIntMapNode*  FindMinNode(AliITSIntMapNode* node) const;
  AliITSIntMapNode*  FindMaxNode(AliITSIntMapNode* node) const;
  void               PrintNode(AliITSIntMapNode* node) const;
  UInt_t             GetTreeHeightNode(AliITSIntMapNode* node) const;

  AliITSIntMapNode*  BalanceNode(Int_t lowInd, Int_t highInd);
  void               ClearFastAccess();
  void               InitFastAccess();
  void               InitFastAccessNode(AliITSIntMapNode* node);
  void               InitFastAccessSerialize();
  void               InitFastAccessSerializeNode(AliITSIntMapNode* node);

};

#endif
