#ifndef ALI_ITS_INTMAPNODE_H
#define ALI_ITS_INTMAPNODE_H

//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// Class for the nodes to put in the integer map (AliITSIntMap)     //
//                                                                  //
//////////////////////////////////////////////////////////////////////  

#include <Rtypes.h>

class AliITSIntMapNode {

 public:
  AliITSIntMapNode();
  AliITSIntMapNode(Int_t key, Int_t val, AliITSIntMapNode* next);
  AliITSIntMapNode(const AliITSIntMapNode& obj);
  virtual ~AliITSIntMapNode();
  AliITSIntMapNode& operator=(const AliITSIntMapNode& obj);

  Int_t             GetKey() const {return fKey;}
  Int_t             GetVal() const {return fVal;}
  AliITSIntMapNode* GetNext() {return fNext;}
  void              SetKey(Int_t key) {fKey=key;}
  void              SetVal(Int_t val) {fVal=val;}
  void              SetNext(AliITSIntMapNode* obj) {fNext = obj;}

 private:
  Int_t fKey;               // key (for sorting)
  Int_t fVal;               // value
  AliITSIntMapNode* fNext;  // pointer to next object in list

};

#endif
