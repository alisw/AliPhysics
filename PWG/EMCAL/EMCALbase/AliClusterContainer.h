#ifndef AliClusterContainer_H
#define AliClusterContainer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TLorentzVector;

class AliVEvent;

#include <iterator>
#include <TArrayI.h>
#include <AliVCluster.h>

#include "AliEmcalContainer.h"

/**
 * @class AliClusterContainer
 * @brief Container structure for EMCAL clusters
 * @ingroup EMCALCOREFW
 * @author Martha Verweij
 * @author Salvatore Aiola
 *
 * Container with name, TClonesArray and cuts for calo clusters
 */
class AliClusterContainer : public AliEmcalContainer {
 public:
  typedef enum AliVCluster::VCluUserDefEnergy_t VCluUserDefEnergy_t;
  
  /**
   * @class accept_iterator
   * @brief stl iterator over accepted clusters
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 16, 2016
   *
   * Bi-directional iterator over all accepted clusters in the array. The
   * iterator is used in the standard way a c++ iterator is used. The
   * container is responsible creating it. Afterwards one uses the prefix
   * or postfix increment or decrement operator, depending on the direction
   * of the iteration.
   *
   * ~~~{.cxx}
   * AliClusterContainer cont;
   * for(AliClusterContainer::iterator it = cont.accept_begin(); it != cont.accept_end(); ++cont){
   *   std::cout << "Cluster energy: " << (*it)->E() std::endl;      // do something with the cluster inside
   * }
   * ~~~
   *
   * As the postfix operator needs to perform a copy of the iterator,
   * the prefix operator will have a better performance and should
   * be preferred.
   *
   * Note: This class is using a list of accepted clusters in the
   * backend. This list is copied in the copy constructor and
   * the assignment operator. Can be done in a more performant way
   * using c++11 shared_ptr.
   */
  class accept_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                               AliVCluster,std::ptrdiff_t,
                                               AliVCluster **, AliVCluster *&>
  {
  public:
    accept_iterator(const AliClusterContainer *cont, int startpos, bool forward  = true);
    accept_iterator(const accept_iterator &other);
    virtual ~accept_iterator() {}
    accept_iterator &operator=(const accept_iterator &other);
    bool operator!=(const accept_iterator &other) const;

    accept_iterator &operator++();
    accept_iterator operator++(int);
    accept_iterator &operator--();
    accept_iterator operator--(int);

    AliVCluster *operator*() const;

  private:
    accept_iterator();

    const AliClusterContainer       *fkContainer;       ///< Container iterated over
    TArrayI                          fAcceptIndices;    ///< Indices of accepted clusters
    int                              fCurrentPos;       ///< Current position inside the container
    bool                             fForward;          ///< Direction, expressed in forward direction
  };

  /**
   * @class all_iterator
   * @brief stl iterator over all clusters
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 16, 2016
   *
   * Bi-directional iterator over all clusters in the array. The
   * iterator is used in the standard way a c++ iterator is used.
   * The container is responsible creating it. Afterwards one uses
   * the prefix or postfix increment or decrement operator, depending
   * on the direction of the iteration.
   *
   * ~~~{.cxx}
   * AliClusterContainer cont;
   * for(AliClusterContainer::iterator it = cont.begin(); it != cont.end(); ++cont){
   *   std::cout << "Cluster energy: " << (*it)->E() std::endl;      // do something with the cluster inside
   * }
   * ~~~
   *
   * If you are using c++11 this reduces to
   *
   * ~~~{.cxx}
   * for(auto it : cont){
   *   std::cout << "Cluster energy: " << it->E() std::endl;      // do something with the cluster inside
   * }
   * ~~~
   *
   * using c++11 range-based iterations.
   *
   * As the postfix operator needs to perform a copy of the iterator,
   * the prefix operator will have a better performance and should
   * be preferred.
   */
  class all_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                            AliVCluster,std::ptrdiff_t,
                                            AliVCluster **, AliVCluster *&>
  {
  public:
    all_iterator(const AliClusterContainer *cont, int startpos, bool forward  = true);
    all_iterator(const all_iterator &other);
    virtual ~all_iterator() {}
    all_iterator &operator=(const all_iterator &other);
    bool operator!=(const all_iterator &other) const;

    all_iterator &operator++();
    all_iterator operator++(int);
    all_iterator &operator--();
    all_iterator operator--(int);

    AliVCluster *operator*() const;

  private:
    all_iterator();

    const AliClusterContainer       *fkContainer;       ///< Container iterated over
    int                              fCurrentPos;       ///< Current position inside the container
    bool                             fForward;          ///< Direction, expressed in forward direction
  };

  AliClusterContainer();
  AliClusterContainer(const char *name); 
  virtual ~AliClusterContainer(){;}

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

  accept_iterator             accept_begin() const;
  accept_iterator             accept_end() const;
  accept_iterator             accept_rbegin() const;
  accept_iterator             accept_rend() const;

  all_iterator                begin() const;
  all_iterator                end() const;
  all_iterator                rbegin() const;
  all_iterator                rend() const;

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

