#ifndef AliClusterContainer_H
#define AliClusterContainer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TLorentzVector;

class AliVEvent;

#include <TArrayI.h>
#include <AliVCluster.h>

#include "AliEmcalContainer.h"

typedef AliEmcalIterableContainerT<AliVCluster> AliClusterIterableContainer;

/**
 * @class AliClusterContainer
 * @brief Container structure for EMCAL clusters
 * @ingroup EMCALCOREFW
 * @author Martha Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * Container with name, TClonesArray and cuts for calo clusters
 */
class AliClusterContainer : public AliEmcalContainer {
 public:
  typedef enum AliVCluster::VCluUserDefEnergy_t VCluUserDefEnergy_t;

  AliClusterContainer();
  AliClusterContainer(const char *name); 
  virtual ~AliClusterContainer(){;}

  virtual TObject *operator[] (int index) const { return GetCluster(index); }

  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const              { return AcceptCluster(i, rejectionReason);}
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const   { return AcceptCluster(dynamic_cast<const AliVCluster*>(obj), rejectionReason);}
  virtual Bool_t              AcceptCluster(Int_t i, UInt_t &rejectionReason)                 const;
  virtual Bool_t              AcceptCluster(const AliVCluster* vp, UInt_t &rejectionReason)   const;
  virtual Bool_t              ApplyClusterCuts(const AliVCluster* clus, UInt_t &rejectionReason) const;
  AliVCluster                *GetAcceptCluster(Int_t i)              const;
  AliVCluster                *GetAcceptClusterWithLabel(Int_t lab)   const;
  void                        SetClusECut(Double_t cut)                    { SetMinE(cut)     ; }
  void                        SetClusPtCut(Double_t cut)                   { SetMinPt(cut)    ; }
  Double_t                    GetClusPtCut()                         const { return GetMinPt(); }
  AliVCluster                *GetCluster(Int_t i)                    const;
  AliVCluster                *GetClusterWithLabel(Int_t lab)         const;
  AliVCluster                *GetLeadingCluster(const char* opt="")       ;
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVCluster* vc, Double_t mass) const;
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVCluster* clus) const;
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      GetNextMomentum(TLorentzVector &mom);
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom);
  AliVCluster                *GetNextAcceptCluster();
  AliVCluster                *GetNextCluster();
  Int_t                       GetNClusters()                         const { return GetNEntries();   }
  Int_t                       GetNAcceptedClusters()                 const;
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

  const char*                 GetTitle() const;

  const AliClusterIterableContainer      all() const;
  const AliClusterIterableContainer      accepted() const;

  AliClusterIterableContainer::iterator  accept_begin()  const { return accepted().begin()   ; }
  AliClusterIterableContainer::iterator  accept_end()    const { return accepted().end()     ; }
  AliClusterIterableContainer::iterator  accept_rbegin() const { return accepted().rbegin()  ; }
  AliClusterIterableContainer::iterator  accept_rend()   const { return accepted().rend()    ; }

  AliClusterIterableContainer::iterator  begin()         const { return all().begin()        ; }
  AliClusterIterableContainer::iterator  end()           const { return all().end()          ; }
  AliClusterIterableContainer::iterator  rbegin()        const { return all().rbegin()       ; }
  AliClusterIterableContainer::iterator  rend()          const { return all().rend()         ; }

 protected:
  
  Double_t         fClusTimeCutLow;             ///< low time cut for clusters
  Double_t         fClusTimeCutUp;              ///< up time cut for clusters
  Bool_t           fExoticCut;                  ///< reject clusters marked as "exotic"
  Double_t         fUserDefEnergyCut[AliVCluster::kLastUserDefEnergy+1]; ///< cut on the energy of the cluster after higher level corrections (see AliVCluster.h)
  Int_t            fDefaultClusterEnergy;       ///< default cluster energy: -1 for clus->E(); otherwise clus->GetUserDefEnergy(fDefaultClusterEnergy)

 private:
  AliClusterContainer(const AliClusterContainer& obj); // copy constructor
  AliClusterContainer& operator=(const AliClusterContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliClusterContainer,5);
  /// \endcond
};

/**
 * Unit test for the iterators. Comparing iterators against for-loop of clusters.
 * All clusters selected in the for-loop must be found in order to pass the test.
 * @ingroup EMCALCOREFW
 * @param cont Cluster container used for the test.
 * @param iteratorType type of the iterator (0 = accept_iterator, 1 = all_iterator)
 * @param verbose Switch on verbosity in case of true
 * @return Result of the unit test (0 - passed, 1 - clusters missing, 2 - excess clusters)
 */
int TestClusterContainerIterator(const AliClusterContainer *const cont, int iteratorType = 0, bool verbose = false);

#endif

