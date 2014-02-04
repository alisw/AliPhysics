#ifndef AliClusterContainer_H
#define AliClusterContainer_H

// $Id$

class TLorentzVector;

class AliVEvent;
class AliVCluster;

#include "AliEmcalContainer.h"

class AliClusterContainer : public AliEmcalContainer {
 public:
  AliClusterContainer();
  AliClusterContainer(const char *name); 
  virtual ~AliClusterContainer(){;}

  Bool_t                      AcceptCluster(AliVCluster         *vp) const;
  AliVCluster                *GetAcceptCluster(Int_t i)              const;
  AliVCluster                *GetAcceptClusterWithLabel(Int_t lab)   const;
  AliVCluster                *GetCluster(Int_t i)                    const;
  AliVCluster                *GetClusterWithLabel(Int_t lab)         const;
  AliVCluster                *GetLeadingCluster(const char* opt="")       ;
  void                        GetMomentum(TLorentzVector &mom, Int_t i) const;
  AliVCluster                *GetNextAcceptCluster(Int_t i=-1)            ;
  AliVCluster                *GetNextCluster(Int_t i=-1)                  ;
  Int_t                       GetNClusters()                         const { return GetNEntries();   }
  void                        SetClassName(const char *clname);
  void                        SetClusECut(Double_t cut)                    { fClusECut      = cut ; }
  void                        SetClusPtCut(Double_t cut)                   { fClusPtCut      = cut ; }
  void                        SetClusTimeCut(Double_t min, Double_t max)   { fClusTimeCutLow = min ; fClusTimeCutUp = max ; }
  void                        SetClusterBitMap(UInt_t m)                   { fClusterBitMap     = m ; }
  void                        SetMCClusterBitMap(UInt_t m)                 { fMCClusterBitMap   = m ; }
  void                        SetMinMCLabel(Int_t s)                       { fMinMCLabel        = s ; }

 protected:
  Double_t                    fClusPtCut;                  // cut on cluster pt
  Double_t                    fClusECut;                   // cut on cluster E
  Double_t                    fClusTimeCutLow;             // low time cut for clusters
  Double_t                    fClusTimeCutUp;              // up time cut for clusters
  UInt_t                      fClusterBitMap;              // bit map of accepted clusters (non MC)
  UInt_t                      fMCClusterBitMap;            // bit map of accepted MC clusters
  Int_t                       fMinMCLabel;                 // minimum MC label value for the tracks/clusters being considered MC particles

 private:
  AliClusterContainer(const AliClusterContainer& obj); // copy constructor
  AliClusterContainer& operator=(const AliClusterContainer& other); // assignment

  ClassDef(AliClusterContainer,2);
};

#endif

