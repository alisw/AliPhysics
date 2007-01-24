#ifndef ALI_ITS_INTMAP_H
#define ALI_ITS_INTMAP_H

//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// This class implements the use of a map of integers.              //
// For simplicity, this version is just a sorted linked list,       //
// but it may be rewritten later for better performance if required.//
//////////////////////////////////////////////////////////////////////  

#include <Rtypes.h>

class AliITSIntMapNode;

class AliITSIntMap {

 public:
  AliITSIntMap();
  AliITSIntMap(const AliITSIntMap& imap);
  virtual ~AliITSIntMap();
  AliITSIntMap& operator=(const AliITSIntMap& imap);

  void               Clear();
  Bool_t             Insert(Int_t key, Int_t val);
  Bool_t             Remove(Int_t key);
  AliITSIntMapNode*  Find(Int_t key) const;
  Int_t              GetKey(UInt_t index) const;
  Int_t              GetVal(UInt_t index) const;
  UInt_t             GetNrEntries() const {return fNrEntries;}
  void               PrintEntries() const;

 private:
  UInt_t             fNrEntries;   // nr of entries in map
  AliITSIntMapNode*  fFirst;       // link to first node of map

};

#endif
