/************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef AliClusterContainer_H
#define AliClusterContainer_H

class TLorentzVector;

class AliVCaloCells;
class AliVEvent;

#include <map>
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

  /// Relates string to the cluster energy enumeration for %YAML configuration
  static const std::map <std::string, VCluUserDefEnergy_t> fgkClusterEnergyTypeMap; //!<!

  /**
   * @brief Default constructor.
   */
  AliClusterContainer();

  /**
   * @brief Standard constructor.
   * @param name Name of the array connected to this container
   */
  AliClusterContainer(const char *name); 

  /**
   * @brief Destructor
   */
  virtual ~AliClusterContainer(){;}

  /**
   * @brief Access to cluster at a given index
   * @param index Index in the container
   * @return Cluster at the index (nullptr if the index is larger than the container size)
   */
  virtual TObject *operator[] (int index) const { return GetCluster(index); }

  /**
   * @brief Check whether the cluster at a given index is accepted
   * @param i Index of the cluster
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the cluster is accepted, false otherwise
   * 
   * The function is used by AliEmcalContainer in order to select the cluster at a given position
   */
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const              { return AcceptCluster(i, rejectionReason);}

  /**
   * @brief Check whether a object is accepted by the cluster cuts defined for this container
   * @param obj Object to be checked (must inherit from AliVCluster)
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the cluster is accepted, false otherwise 
   * 
   * The function is used by AliEmcalContainer in order to select the cluster
   */
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const   { return AcceptCluster(dynamic_cast<const AliVCluster*>(obj), rejectionReason);}

  /**
   * @brief Check whether the cluster at a given index is accepted
   * @param i Index of the cluster
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the cluster is accepted, false otherwise
   */
  virtual Bool_t              AcceptCluster(Int_t i, UInt_t &rejectionReason)                 const;

  /**
   * @brief Check whether a cluster is accepted by the cluster cuts defined for this container
   * @param vp Cluster to be checked
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the cluster is accepted, false otherwise 
   */
  virtual Bool_t              AcceptCluster(const AliVCluster* vp, UInt_t &rejectionReason)   const;

  /**
   * @brief Apply cluster selection cuts
   * @param clus The cluster to which the cuts will be applied
   * @param[out] rejectionReason Contains the bit specifying the rejection reason
   * @return True if the cluster is accepted, false otherwise
   * 
   * Cluster selection cuts contain time, energy, exoticity, ... . 
   * Kinematic cuts are handled by the AliEmcalContainer directly
   */
  virtual Bool_t              ApplyClusterCuts(const AliVCluster* clus, UInt_t &rejectionReason) const;

  /**
   * @brief Check if the particle passed the special cut on the number of cells
   * @param clus Cluster to be checked
   * @return True if the cluster is accepted, false otherwise
   * 
   * Check if cluster has only one cell but should be accepted for the analysis
   * Only below a certain energy threshold (recommended = 4 GeV)
   * Based on the output of AliEmcalCorrectionClusterLowEnergyEfficiency.cxx
   * Chi2() value is used to flag clusters that should pass the cut
   * Clusters will only be flagged in MC, data does not change
   */
  virtual Bool_t              GetPassedSpecialNCell(const AliVCluster* clus) const;

  /**
   * @brief Access to cluster at position i if the cluster is accepted
   * @param i Index of the custer in the array
   * @return Cluster at index i (nullptr if index is larger than array size or cluster is not accepted)
   */
  AliVCluster                *GetAcceptCluster(Int_t i)              const;

  /**
   * @brief Get cluster based on cluster label in array if the cluster is accepted
   * @param lab Cluster label
   * @return Cluster associated with label (nullptr if cluster is not fould or cluster is not accepted)
   */
  AliVCluster                *GetAcceptClusterWithLabel(Int_t lab)   const;

  /**
   * @brief Set the default cluster energy cut
   * @param cut Cluster energy cut
   * @deprecated Use SetClusUserDefEnergyCut in order to specify which energy type is used for cluster selection
   * 
   * The cut is applied on the default cluster energy
   */
  void                        SetClusECut(Double_t cut)                    { SetMinE(cut)     ; }

  /**
   * @brief Set the cluster \f$p_{t}\f$-cut
   * @param \f$p_{t}\f$-cut of the cluster
   */
  void                        SetClusPtCut(Double_t cut)                   { SetMinPt(cut)    ; }

  /**
   * @brief Get the cluster \f$p_{t}\f$-cut
   * @return \f$p_{t}\f$-cut of the cluster
   */
  Double_t                    GetClusPtCut()                         const { return GetMinPt(); }

  /**
   * Get \f$i^{th|}\f$ cluster in array
   * @param i Position in array
   * @return Cluster at positon (nullptr if index > array size)
   */
  AliVCluster                *GetCluster(Int_t i)                    const;

  /**
   * Get cluster based on cluster label in array
   * @param lab Cluster label
   * @return Cluster associated with label (nullptr if cluster is not fould)
   */
  AliVCluster                *GetClusterWithLabel(Int_t lab)         const;

  /**
  * @brief Get the leading cluster 
  * @param opt Option: e - use cluster energy (default: \f$E_{t}\f$)
  * @return Leading cluster 
  */
  AliVCluster                *GetLeadingCluster(const char* opt="")       ;

  /**
   * @brief Fill momnetum vector for a given cluster under a user-defined mass hypothesis
   * @param[out] mom Cluster momentum vector
   * @param vc Cluster for which to calculate the momentum vector
   * @param mass Mass hypothesis (default: 0 if mass hypothesis is negative)
   * @return True if the cluster position is valid, false if the cluster position is invalid
   */
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVCluster* vc, Double_t mass) const;

  /**
   * @brief Fill momnetum vector for a given cluster
   * @param[out] mom Cluster momentum vector
   * @param clus Cluster for which to calculate the momentum vector
   * @return True if the cluster position is valid, false if the cluster position is invalid
   * 
   * As mass hypothesis the mass hypothesis specified in AliEmcalContainer 
   * is used (default: -1). In case of a negative mass hypothesis masses
   * particles (photons) are used instead.
   */
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVCluster* clus) const;

  /**
   * @brief Fill momentum vector of the \f$i^{th}\f$ cluster in the array
   * @param[out] mom Cluster momentum vector
   * @param i Index of the cluster in the arrays
   * @return True if access was successful, false if the index is beyond the container size
   */
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Fill momentum vector of the i^th cluster in the array in case the cluster was accepted
   * @param[out] mom Cluster momentum vector
   * @param i Index of the cluster in the arrays
   * @return True if access was successful, false if the index is beyond the container size or the cluster was not accepted
   */
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Fill momentum vector of the next cluster in the array
   * @param[out] mom Cluster momentum vector
   * @return True if access was successful, false if there is no more cluster in the array
   * @deprecated Use all_momentum() in order to iterate over all cluster momentum vectors in the container
   * 
   * Internal iterator over cluster momenta. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  Bool_t                      GetNextMomentum(TLorentzVector &mom);

  /**
   * Fill momentum vector of the next accepted cluster in the array
   * @param[out] mom Cluster momentum vector
   * @return True if access was successful, false if there is no more accepted cluster in the array
   * @deprecated Use accepted_momentum() in order to iterate over accepted cluster momentum vectors in the container
   * 
   * Internal iterator over accepted cluster momenta. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom);

  /**
   * @brief Get next accepted cluster in the container
   * @return Next accepted cluster (nullptr if no more accepted clusters are available)
   * @deprecated Use accepted() to interate over the accepted clusters in the container
   * 
   * Internal iterator over accepted clusters. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  AliVCluster                *GetNextAcceptCluster();

  /**
   * @brief Get next cluster in the container
   * @return Next custer (nullptr if no more clusters are available)
   * @deprecated Use all() to interate over the clusters in the container
   * 
   * Internal iterator over accepted clusters. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  AliVCluster                *GetNextCluster();

  /**
   * @brief Get the number of clusters in the array
   * @return Number of clusters in the array 
   */
  Int_t                       GetNClusters()                         const { return GetNEntries();   }

  /**
   * @brief Get number of accepted clusters in the container
   * @return Number of accepted clusters
   */
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
  void 						            SetEmcalM02Range(Double_t min, Double_t max) { fEmcalMinM02 = min; fEmcalMaxM02 = max; }
  void                        SetEmcalMaxM02Energy(Double_t max)           { fEmcalMaxM02CutEnergy = max; }
  void                        SetMaxFractionEnergyLeadingCell(Double_t max)  { fMaxFracEnergyLeadingCell = max; }
  void                        SetEmcalMaxNCellEfficiencyEnergy(Double_t max) { fEmcalMaxNCellEffCutEnergy = max; }

  /**
   * @brief Connect the container to the array with content stored inside the virtual event.
   * @param event Input event containing the array with content.
   * 
   * The object name in the event must match the name given in the constructor.
   * Additionally register the array into the index map and connect EMCAL cells
   */
  void                        SetArray(const AliVEvent * event);

  /**
   * @brief Set the energy cut of the applied on cluster energy of type t
   * @param t Cluster energy type (base energy, non-linearity corrected energy, hadronically corrected energy)
   * @param cut Cluster energy cut
   */
  void                        SetClusUserDefEnergyCut(Int_t t, Double_t cut);

  /**
   * @brief Get the energy cut of the applied on cluster energy of type t
   * @param t Cluster energy type (base energy, non-linearity corrected energy, hadronically corrected energy)
   * @return Cluster energy cut
   */
  Double_t                    GetClusUserDefEnergyCut(Int_t t) const;

  void                        SetClusNonLinCorrEnergyCut(Double_t cut)                     { SetClusUserDefEnergyCut(AliVCluster::kNonLinCorr, cut); }
  void                        SetClusHadCorrEnergyCut(Double_t cut)                        { SetClusUserDefEnergyCut(AliVCluster::kHadCorr, cut)   ; }
  void                        SetDefaultClusterEnergy(Int_t d)                             { fDefaultClusterEnergy = d                             ; }

  Int_t                       GetDefaultClusterEnergy() const                              { return fDefaultClusterEnergy                          ; }

  const char*                 GetTitle() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Get the EMCal container utils associated with particle containers
  static const AliEmcalContainerIndexMap <TClonesArray, AliVCluster>& GetEmcalContainerIndexMap() { return fgEmcalContainerIndexMap; }

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliClusterIterableContainer      all() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliClusterIterableContainer      accepted() const;

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliClusterIterableMomentumContainer      all_momentum() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliClusterIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  /**
   * @brief Create default array name for the cluster container. 
   * @param[in] ev Input event, used for data type selection
   * @return Appropriate default array name
   * 
   * The default array name will be
   * - *caloClusters* in case of AOD event
   * - *CaloClusters* in case of ESD event
   */
  virtual TString             GetDefaultArrayName(const AliVEvent * const ev) const;


#if !(defined(__CINT__) || defined(__MAKECINT__))
  static AliEmcalContainerIndexMap <TClonesArray, AliVCluster> fgEmcalContainerIndexMap; //!<! Mapping from containers to indices
#endif

  AliVCaloCells   *fEMCALCells;                 ///< pointer to EMCAL cells object
  Double_t         fClusTimeCutLow;             ///< low time cut for clusters
  Double_t         fClusTimeCutUp;              ///< up time cut for clusters
  Bool_t           fExoticCut;                  ///< reject clusters marked as "exotic"
  Double_t         fUserDefEnergyCut[AliVCluster::kLastUserDefEnergy+1]; ///< cut on the energy of the cluster after higher level corrections (see AliVCluster.h)
  Int_t            fDefaultClusterEnergy;       ///< default cluster energy: -1 for clus->E(); otherwise clus->GetUserDefEnergy(fDefaultClusterEnergy)
  Bool_t           fIncludePHOS;                ///< flag to accept PHOS clusters in addition to EMCal clusters
  Bool_t           fIncludePHOSonly;            ///< flag to accept only PHOS clusters (and reject EMCal clusters)
  Int_t            fPhosMinNcells;              ///< min number of phos cells per cluster
  Double_t         fPhosMinM02;                 ///< min value of M02 for phos clusters
  Double_t		     fEmcalMinM02;				   ///< min value of M02 for EMCAL clusters
  Double_t 		     fEmcalMaxM02;				   ///< max value of M02 for EMCAL clusters
  Double_t         fEmcalMaxM02CutEnergy;       ///< max EMCal cluster energy for which to apply M02 cut
  Double_t         fMaxFracEnergyLeadingCell;   ///< max fraction of energy in the leading cell
  Double_t         fEmcalMaxNCellEffCutEnergy;  ///< maximum energy that a 1 cell cluster can have to be considered for NCell efficiency

 private:
  AliClusterContainer(const AliClusterContainer& obj); // copy constructor
  AliClusterContainer& operator=(const AliClusterContainer& other); // assignment

  ClassDef(AliClusterContainer,13);
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
