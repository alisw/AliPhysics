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
  
public:
  Bool_t  GetInfoFromId(Int_t idclus, Int_t &index, Int_t &eta, Int_t &phi)     const;
  Bool_t  GetInfoFromIndex(Int_t index, Int_t &idclus, Int_t &eta, Int_t &phi)  const;
  Int_t   GetSize()                                                             const;
  Bool_t  ContainsId(Int_t idclus)                                              const;
  Bool_t  ContainsIndex(Int_t index)                                            const;
  Int_t   GetLastElementPosition()                                              const;
  void    Add(Int_t idclus, Int_t index, Int_t eta, Int_t phi);
  Bool_t  SetIndexFromId(Int_t idclus, Int_t index);
  Bool_t  RemoveId(Int_t idclus);
  Bool_t  RemoveIndex(Int_t index);
  void    Expand(Int_t size);
  
  virtual void Clear(Option_t* option = "");
  
protected:
  Int_t GetPositionFromId(Int_t idclus) const;
  Int_t GetPositionFromIndex(Int_t index) const;
  
  Int_t     fLastPos;         // Last non-empty position in the arrays
  TArrayI  *fIds;             // Array containing unique IDs of clusters
  TArrayI  *fIndexes;         // Array containing corresponding indexes in the cluster collection (usually a TClonesObject)
  TArrayI  *fPhi;             // Array containing corresponding Phi index
  TArrayI  *fEta;             // Array containing corresponding Eta index
  
private:
  AliEMCALFixedWindowClusterInfo(const AliEMCALFixedWindowClusterInfo&);            // not implemented
  AliEMCALFixedWindowClusterInfo &operator=(const AliEMCALFixedWindowClusterInfo&); // not implemented
  
  ClassDef(AliEMCALFixedWindowClusterInfo, 1);
};
#endif //ALIEMCALFIXEDWINDOWCLUSTERINFO_H
