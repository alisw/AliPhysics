// $Id$

#ifndef ALIJFCLUSTERDIFFERENCEH
#define ALIJFCLUSTERDIFFERENCEH

#include <Riostream.h>
#include "AliJFCluster.h"

class AliJFClusterDifference
{
 public:
  AliJFClusterDifference();
  AliJFClusterDifference(const AliJFClusterDifference &copy);
  AliJFClusterDifference(AliJFCluster *i, AliJFCluster *j);
  AliJFClusterDifference(Float_t dij,AliJFCluster *i, AliJFCluster *j);
  virtual ~AliJFClusterDifference() {}

  AliJFClusterDifference& operator=(const AliJFClusterDifference &copy);
  friend ostream& operator<<(ostream &o, const AliJFClusterDifference &j);
  friend bool operator< (const AliJFClusterDifference &a, const AliJFClusterDifference &b);
  friend bool operator==(const AliJFClusterDifference &a, const AliJFClusterDifference &b);

  inline Int_t GetNLastMerge() const {return fNLastMerge;}
  inline AliJFCluster* GetI()     const {return fI;}
  inline AliJFCluster* GetJ()     const {return fJ;}
  inline Float_t GetDij()      const {return fDij;}
  
  inline Bool_t IsDiagonal()      const {return (fI==fJ ? kTRUE:kFALSE);}
  inline Bool_t IsValid()         const {return (fI->IsValid()&&fJ->IsValid());}
  inline Bool_t IsValidDiagonal() const {return ((fI->IsValid())&&(fI->GetNMerge()==fNLastMerge));}
  inline Bool_t IsValidEntry()    const {
    if(IsDiagonal()) return IsValidDiagonal(); 
    else return IsValid();
  }
  inline Bool_t IsValidPointer()  const {return (fI!=NULL && fJ!=NULL ? kTRUE:kFALSE);}

  Float_t SetValues(AliJFCluster *i,AliJFCluster *j);

 protected:
  Int_t fNLastMerge;
  AliJFCluster *fI;
  AliJFCluster *fJ;
  Float_t fDij;

  ClassDef(AliJFClusterDifference,0) //AliJFClusterDifference Class
};

inline bool operator<(const AliJFClusterDifference &a, const AliJFClusterDifference &b){
  return (a.GetDij()<b.GetDij());
}

inline bool operator==(const AliJFClusterDifference &a, const AliJFClusterDifference &b){
  if((a.GetI()==b.GetI())&&(a.GetJ()==b.GetJ())) return true;
  return false;
}

#endif /*ALIJFCLUSTERDIFFERENCEH*/
