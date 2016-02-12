#ifndef AliClusterContainer_H
#define AliClusterContainer_H

class TLorentzVector;

class AliVEvent;

#include <AliVCluster.h>

#include "AliEmcalContainer.h"

class AliClusterContainer : public AliEmcalContainer {
 public:
  typedef enum AliVCluster::VCluUserDefEnergy_t VCluUserDefEnergy_t;
  
  AliClusterContainer();
  AliClusterContainer(const char *name); 
  virtual ~AliClusterContainer(){;}

  Bool_t                      AcceptCluster(AliVCluster* vp)              ;
  AliVCluster                *GetAcceptCluster(Int_t i)                   ;
  AliVCluster                *GetAcceptClusterWithLabel(Int_t lab)        ;
  AliVCluster                *GetCluster(Int_t i)                    const;
  AliVCluster                *GetClusterWithLabel(Int_t lab)         const;
  AliVCluster                *GetLeadingCluster(const char* opt="")       ;
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetNextMomentum(TLorentzVector &mom, Int_t i=-1);
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom, Int_t i=-1);
  AliVCluster                *GetNextAcceptCluster(Int_t i=-1)            ;
  AliVCluster                *GetNextCluster(Int_t i=-1)                  ;
  Int_t                       GetNClusters()                         const { return GetNEntries();   }
  Int_t                       GetNAcceptedClusters()                      ;
  void                        SetClassName(const char* clname);
  void                        SetClusECut(Double_t cut)                    { fClusECut        = cut ; }
  void                        SetClusPtCut(Double_t cut)                   { fClusPtCut       = cut ; }
  void                        SetClusTimeCut(Double_t min, Double_t max)   { fClusTimeCutLow  = min ; fClusTimeCutUp = max ; }
  void                        SetMinMCLabel(Int_t s)                       { fMinMCLabel      = s   ; }
  void                        SetMaxMCLabel(Int_t s)                       { fMaxMCLabel      = s   ; }
  void                        SetMCLabelRange(Int_t min, Int_t max)        { SetMinMCLabel(min)     ; SetMaxMCLabel(max)    ; }
  void                        SetExoticCut(Bool_t e)                       { fExoticCut       = e   ; }

  void                        SetClusUserDefEnergyCut(Int_t t, Double_t cut);
  Double_t                    GetClusUserDefEnergyCut(Int_t t) const;

  void                        SetClusNonLinCorrEnergyCut(Double_t cut)                     { SetClusUserDefEnergyCut(AliVCluster::kNonLinCorr, cut); }
  void                        SetClusHadCorrEnergyCut(Double_t cut)                        { SetClusUserDefEnergyCut(AliVCluster::kHadCorr, cut)   ; }
  void                        SetDefaultClusterEnergy(Int_t d)                             { fDefaultClusterEnergy = d                             ; }

  Int_t                       GetDefaultClusterEnergy() const                              { return fDefaultClusterEnergy                          ; }
  Double_t                    GetClusPtCut() const                                         { return fClusPtCut                                     ; }

  const char*                 GetTitle() const;

 protected:
  
  Double_t         fClusPtCut;                  // cut on cluster pt
  Double_t         fClusECut;                   // cut on cluster E
  
  Double_t         fClusTimeCutLow;             // low time cut for clusters
  Double_t         fClusTimeCutUp;              // up time cut for clusters
  Int_t            fMinMCLabel;                 // minimum MC label
  Int_t            fMaxMCLabel;                 // maximum MC label

  Bool_t           fExoticCut;                  // reject clusters marked as "exotic"
  Double_t         fUserDefEnergyCut[AliVCluster::kLastUserDefEnergy+1]; // cut on the energy of the cluster after higher level corrections (see AliVCluster.h)
  Int_t            fDefaultClusterEnergy;       // default cluster energy: -1 for clus->E(); otherwise clus->GetUserDefEnergy(fDefaultClusterEnergy)

 private:
  AliClusterContainer(const AliClusterContainer& obj); // copy constructor
  AliClusterContainer& operator=(const AliClusterContainer& other); // assignment

  ClassDef(AliClusterContainer,4);
};

#endif

