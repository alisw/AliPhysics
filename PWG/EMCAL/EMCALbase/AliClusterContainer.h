#ifndef AliClusterContainer_H
#define AliClusterContainer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TLorentzVector;

class AliVEvent;

#include <TArrayI.h>
#include <AliVCluster.h>

#include "AliEmcalContainer.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVCluster, EMCALIterableContainer::operator_star_object<AliVCluster> > AliClusterIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVCluster, EMCALIterableContainer::operator_star_pair<AliVCluster> > AliClusterIterableMomentumContainer;
#endif

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
  void                        SetIncludePHOS(Bool_t b)                     { fIncludePHOS = b       ; }
  void                        SetIncludePHOSonly(Bool_t b)                 { fIncludePHOSonly = b   ; }
  void                        SetPhosMinNcells(Int_t n)                    { fPhosMinNcells = n; }
  void                        SetPhosMinM02(Double_t m)                    { fPhosMinM02 = m; }
  void                        SetArray(const AliVEvent * event);
  void                        SetClusUserDefEnergyCut(Int_t t, Double_t cut);
  Double_t                    GetClusUserDefEnergyCut(Int_t t) const;

  void                        SetClusNonLinCorrEnergyCut(Double_t cut)                     { SetClusUserDefEnergyCut(AliVCluster::kNonLinCorr, cut); }
  void                        SetClusHadCorrEnergyCut(Double_t cut)                        { SetClusUserDefEnergyCut(AliVCluster::kHadCorr, cut)   ; }
  void                        SetDefaultClusterEnergy(Int_t d)                             { fDefaultClusterEnergy = d                             ; }

  Int_t                       GetDefaultClusterEnergy() const                              { return fDefaultClusterEnergy                          ; }

  const char*                 GetTitle() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Get the EMCal container utils associated with particle containers
  static const AliEmcalContainerIndexMap <TClonesArray, AliVCluster>& GetEmcalContainerIndexMap() { return fgEmcalContainerIndexMap; }

  const AliClusterIterableContainer      all() const;
  const AliClusterIterableContainer      accepted() const;

  const AliClusterIterableMomentumContainer      all_momentum() const;
  const AliClusterIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  /**
   * Create default array name for the cluster container. The
   * default array name will be
   * - *caloClusters* in case of AOD event
   * - *CaloClusters* in case of ESD event
   * @param[in] ev Input event, used for data type selection
   * @return Appropriate default array name
   */
  virtual TString             GetDefaultArrayName(const AliVEvent * const ev) const;

  
#if !(defined(__CINT__) || defined(__MAKECINT__))
  static AliEmcalContainerIndexMap <TClonesArray, AliVCluster> fgEmcalContainerIndexMap; //!<! Mapping from containers to indices
#endif

  Double_t         fClusTimeCutLow;             ///< low time cut for clusters
  Double_t         fClusTimeCutUp;              ///< up time cut for clusters
  Bool_t           fExoticCut;                  ///< reject clusters marked as "exotic"
  Double_t         fUserDefEnergyCut[AliVCluster::kLastUserDefEnergy+1]; ///< cut on the energy of the cluster after higher level corrections (see AliVCluster.h)
  Int_t            fDefaultClusterEnergy;       ///< default cluster energy: -1 for clus->E(); otherwise clus->GetUserDefEnergy(fDefaultClusterEnergy)
  Bool_t           fIncludePHOS;                ///< flag to accept PHOS clusters in addition to EMCal clusters
  Bool_t           fIncludePHOSonly;            ///< flag to accept only PHOS clusters (and reject EMCal clusters)
  Int_t            fPhosMinNcells;              ///< min number of phos cells per cluster
  Double_t         fPhosMinM02;                 ///< min value of M02 for phos clusters

 private:
  AliClusterContainer(const AliClusterContainer& obj); // copy constructor
  AliClusterContainer& operator=(const AliClusterContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliClusterContainer,9);
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

