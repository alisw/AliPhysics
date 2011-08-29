#ifndef ALIEMCALFIXEDWINDOWCLUSTERINFO_H
#define ALIEMCALFIXEDWINDOWCLUSTERINFO_H

// $Id$

#include <TNamed.h>

class TArrayI;

class AliEMCALFixedWindowClusterInfo : public TNamed {
  
public:
  AliEMCALFixedWindowClusterInfo();
  AliEMCALFixedWindowClusterInfo(const char* name, Int_t size = 1000);
  AliEMCALFixedWindowClusterInfo(const TString& name, Int_t size = 1000);
  virtual ~AliEMCALFixedWindowClusterInfo();
  
  Bool_t GetInfoFromId(Int_t idclus, Int_t &index, Int_t &eta, Int_t &phi);
  Bool_t GetInfoFromIndex(Int_t index, Int_t &idclus, Int_t &eta, Int_t &phi);
  
  void Add(Int_t idclus, Int_t index, Int_t eta, Int_t phi);
  Bool_t SetIndexFromId(Int_t idclus, Int_t index);
  Bool_t RemoveId(Int_t idclus);
  Bool_t RemoveIndex(Int_t index);
  void Expand(Int_t size);
  Int_t GetSize();
  Bool_t ContainsId(Int_t idclus);
  Bool_t ContainsIndex(Int_t index);
  Int_t GetLastElementId();
  
  virtual void Clear(Option_t* option = "");
  
protected:
  Int_t GetPositionFromId(Int_t idclus);
  Int_t GetPositionFromIndex(Int_t index);
  
private:
  Int_t     fSize;
  Int_t     lastElement;
  TArrayI  *fIds;
  TArrayI  *fIndexes;
  TArrayI  *fPhi;
  TArrayI  *fEta;
  
  ClassDef(AliEMCALFixedWindowClusterInfo, 1);
};
#endif //ALIEMCALFIXEDWINDOWCLUSTERINFO_H
