#ifndef ALIITSINTMAPNODE_H
#define ALIITSINTMAPNODE_H

//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// Class for the nodes to put in the integer map (AliITSIntMap)     //
//                                                                  //
//////////////////////////////////////////////////////////////////////  

#include <Rtypes.h>

class AliITSIntMapNode {

 public:
  AliITSIntMapNode();
  AliITSIntMapNode(Int_t key, Int_t val, AliITSIntMapNode* left, AliITSIntMapNode* right);
  AliITSIntMapNode(const AliITSIntMapNode& obj);
  virtual ~AliITSIntMapNode();
  AliITSIntMapNode& operator=(const AliITSIntMapNode& obj);

  Int_t             Key() const {return fKey;}
  Int_t             Val() const {return fVal;}
  AliITSIntMapNode*& Left() {return fLeft;}
  AliITSIntMapNode*& Right() {return fRight;}

  void              SetKey(Int_t key) {fKey=key;}
  void              SetVal(Int_t val) {fVal=val;}
  void              SetLeft(AliITSIntMapNode* obj) {fLeft = obj;}
  void              SetRight(AliITSIntMapNode* obj) {fRight = obj;}

 private:
  Int_t fKey;               // key (for sorting)
  Int_t fVal;               // value
  AliITSIntMapNode* fLeft;  // pointer to left object in bin tree
  AliITSIntMapNode* fRight; // pointer to right object in bin tree

};

#endif
