/************************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                 *
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
/// \cond CLASSIMP
#include "AliEmcalIterableContainer.h"
/// \endcond

#ifndef ALIEMCALCONTAINER_H
#define ALIEMCALCONTAINER_H

class TLorentzVector;
class AliTLorentzVector;
class AliVEvent;
class AliNamedArrayI;
class AliVParticle;

#include <TNamed.h>
#include <TClonesArray.h>

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<TObject, EMCALIterableContainer::operator_star_object<TObject> > AliEmcalIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<TObject, EMCALIterableContainer::operator_star_pair<TObject> > AliEmcalIterableMomentumContainer;
#endif

/**
 * @class AliEmcalContainer
 * @brief Base class for container structures within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 *
 * This class is the base class for container object used in the EMCAL framework.
 * The main purpose is to connect this to objects stored as list objects in the
 * input event, among them particles, EMCAL clusters, or jets. The core of the container
 * is a pointer to a TClonesArray representing the the content in the input event.
 *
 * Containers can be connected. For example, jet containers need access to the track
 * container and the cluster container in case constituent information is required.
 *
 * In addition, kinematical cuts can be applied, accessing only content which is selected
 * using the selection cuts to be applied.
 *
 * Iterator functionality is provided by the AliEmcalIterableContainer. Iterators are
 * available both for all entries and accepted entries in the container
 * - The function all creates an interable container over all entries in the EMCAL container
 * - The function accept creates an iterable container over accepted entries in the EMCal container
 * The following example demonstrates the usage of the iterator
 *
 * ~~~{.cxx}
 * AliEmcalContainer *cont;
 * AliEmcalIterableContainer alliter = cont->all(),
 *                           acceptiter = cont->accepted();
 * for(AliEmcalIterableContainer::iterator it = alliter.begin(); it != alliter.end(); ++it){
 *   // Iterating over all objects in the EMCAL container
 *   // Do something with the object
 * }
 * for(AliEmcalIterableContainer::iterator it = acceptiter.begin(); it != acceptiter.end(); ++it){
 *   // Iterating over all objects in the EMCAL container
 *   // Do something with the object
 * }
 * ~~~
 *
 * Using c++11 and range-based iteration this code can be simplified to
 *
 * ~~~{.cxx}
 * AliEmcalContainer *cont;
 * for(auto en : cont->all()) { // Iterating over all objects in the container
 *   // Do something with the object
 * }
 * for(auto en : cont->accepted()) { // Iterating over accepted objects in the container
 *   // Do something with the object
 * }
 * ~~~
 *
 * The usage of EMCAL containers is described under \subpage EMCALcontainers
 */
class AliEmcalContainer : public TObject {
 public:
  /**
   * @enum RejectionReason
   * @brief Bit definition for the reason a particle was rejected
   */
  enum RejectionReason {
    // General
    kNullObject = 1<<0,                  ///< Object is NULL
    kPtCut = 1<<1,                       ///< \f$ p_{t} \f$ cut
    kAcceptanceCut = 1<<2,               ///< particle not in acceptance in \f$ \eta \f$ and/or \f$ \phi \f$
    kMCLabelCut = 1<<3,                  ///< Invalid MC label
    kBitMapCut = 1<<4,                   ///< kBitMapCut
    kHFCut = 1<<5,                       ///< kHFCut
    // leave bit 6 free for future implementations
    
    // AliParticleContainer
    kNotHybridTrack = 1<<7,              ///< Track did not pass the hybrid track cuts
    kMCFlag = 1<<8,                      ///< Cut on the MC flag
    kMCGeneratorCut = 1<<9,              ///< Generator flag mismatch
    kChargeCut = 1<<10,                  ///< Particle charge did not match
    kMinDistanceTPCSectorEdgeCut = 1<<11,///< Track too close to the TPC sector boundary
    // leave bit 12 free for future implementations

    // AliClusterContainer
    kIsEMCalCut = 1<<13,                 ///< Cluster not in the EMCAL
    kTimeCut = 1<<14,                    ///< Cell time cut not passed
    kEnergyCut = 1<<15,                  ///< Energy below threshold
    kExoticCut = 1<<16,                  ///< Cluster is exotic cluster
    // leave bit 17 free for future implementations

    // AliJetContainer
    kAreaCut = 1<<18,                    ///< Cut on the jet area
    kAreaEmcCut = 1<<19,                 ///< Cut on the jet area in the EMCAL
    kZLeadingChCut = 1<<20,              ///< Cut on the z of the leading charged constituent
    kZLeadingEmcCut = 1<<21,             ///< Cut on the z of the leading particle in the EMCAL
    kNEFCut = 1<<22,                     ///< Cut on the neutral energy fraction
    kMinLeadPtCut = 1<<23,               ///< Cut on the minimum \f$ p_{t} \f$ of the leading particle
    kMaxTrackPtCut = 1<<24,              ///< Cut on the maximum track \f$ p_{t} \f$
    kMaxClusterPtCut = 1<<25,            ///< Cut on the maximum cluster \f$ p_{t} \f$
    kFlavourCut = 1<<26,                 ///< Cut on flavour content in the jet
    kTagStatus = 1<<27,                  ///< Cut on jet tag status
    kMinNConstituents = 1<<28,            ///< Cut on the minimum number of constituents
    kOverlapTpcHole = 1<<29             ///<Cut  on the regions of acceptance with bad sectors 
  };

  /**
   * @brief Default constructor. 
   * 
   * This constructor is only for ROOT I/O and not to be used by users. 
   * The container will not connect to an array in the input event.
   */
  AliEmcalContainer();

  /**
   * @brief Standard (named) constructor. 
   * 
   * The name provided must match the name of the array inside the list 
   * objects in the input event the EMCAL container connects to. The 
   * EMCAL container can get a different name, to be specified in the 
   * function SetEvent.
   * @param name Name of the container in the input event.
   */
  AliEmcalContainer(const char *name); 

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalContainer(){;}

  /**
   * @brief Index operator.
   * 
   * Accessing object in the container at a given index. Operates on all 
   * objects inside the container. 
   * @param index Index of the object to access
   * @return Object at the given index (NULL if out of range)
   */
  virtual TObject *operator[](int index) const = 0;

  /**
   * @brief Apply kinematical selection to the momentum vector provided. 
   * 
   * Selection is done in
   * - \f$ p_{t} \f$ (E)
   * - \f$ \eta \f$
   * - \f$ \phi \f$
   * @param[in] mom Momentum vector to select
   * @param[out] rejectionReason Bitmap for reason why object is rejected
   * @return True if the momentum vector is selected, false otherwise
   */
  virtual Bool_t              ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const;
  TClonesArray               *GetArray()                      const { return fClArray                   ; }
  const TString&              GetArrayName()                  const { return fClArrayName               ; }
  const TString&              GetClassName()                  const { return fClassName                 ; }
  TClass                     *GetClass()                      const { return fClArray == 0 ? 0 : fClArray->GetClass(); }
  Double_t                    GetMinE()                       const { return fMinE  ; }
  Double_t                    GetMaxE()                       const { return fMaxE  ; }
  Double_t                    GetMinPt()                      const { return fMinPt  ; }
  Double_t                    GetMaxPt()                      const { return fMaxPt  ; }
  Double_t                    GetMinEta()                     const { return fMinEta ; }
  Double_t                    GetMaxEta()                     const { return fMaxEta ; }
  Double_t                    GetMinPhi()                     const { return fMinPhi ; }
  Double_t                    GetMaxPhi()                     const { return fMaxPhi ; }
  Double_t                    GetEtaSwing()                   const { return fMaxEta - fMinEta; }
  Double_t                    GetPhiSwing()                   const { return fMaxPhi - fMinPhi; }
  Double_t                    GetAcceptance()                 const { return GetEtaSwing() * GetPhiSwing(); }
  Int_t                       GetCurrentID()                  const { return fCurrentID                 ; }
  Bool_t                      GetIsParticleLevel()            const { return fIsParticleLevel           ; }

  /**
   * @brief Get the index in the container from a given label
   * @param lab Label to check
   * @return Index (-1 if not found)
   */
  Int_t                       GetIndexFromLabel(Int_t lab)    const;

  Int_t                       GetNEntries()                   const { return fClArray ? fClArray->GetEntriesFast() : 0 ; }
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) const = 0;
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) const = 0;
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom) = 0;
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom) = 0;
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const = 0;
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const = 0;

  /**
   * @brief Count accepted entries in the container
   * @return Number of accepted events in the container
   */
  Int_t                       GetNAcceptEntries() const;

  /**
   * @brief Reset the iterator to a given index
   * 
   * Resetting the internal iterator to a new starting position. This can be
   * within the range of indices of the EMCAL container. To start from the beginning
   * one needs to select -1 as iterator index.
   * 
   * @param[in] Index to be stated from (must be 1 lower than the desired starting point)
   */
  void                        ResetCurrentID(Int_t i=-1)            { fCurrentID = i                    ; }

  /**
   * Connect the container to the array with content stored inside the virtual event.
   * The object name in the event must match the name given in the constructor
   * @param event Input event containing the array with content.
   */
  virtual void                SetArray(const AliVEvent *event);
  void                        SetArrayName(const char *n)           { fClArrayName = n                  ; }
  void                        SetVertex(Double_t *vtx)              { memcpy(fVertex, vtx, sizeof(Double_t) * 3); }
  void                        SetBitMap(UInt_t m)                   { fBitMap = m                       ; }
  void                        SetIsParticleLevel(Bool_t b)          { fIsParticleLevel = b              ; }
  void                        SortArray()                           { fClArray->Sort()                  ; }

  TClass*                     GetLoadedClass()                      { return fLoadedClass               ; }

  /**
   * @brief Preparation for the next event.
   * @param[in] event The event to be processed.
   */
  virtual void                NextEvent(const AliVEvent *event);
  void                        SetMinMCLabel(Int_t s)                            { fMinMCLabel      = s   ; }
  void                        SetMaxMCLabel(Int_t s)                            { fMaxMCLabel      = s   ; }
  void                        SetMCLabelRange(Int_t min, Int_t max)             { SetMinMCLabel(min)     ; SetMaxMCLabel(max)    ; }
  void                        SetELimits(Double_t min, Double_t max)    { fMinE   = min ; fMaxE   = max ; }
  void                        SetMinE(Double_t min)                     { fMinE   = min ; }
  void                        SetMaxE(Double_t max)                     { fMaxE   = max ; }
  void                        SetPtLimits(Double_t min, Double_t max)   { fMinPt  = min ; fMaxPt  = max ; }
  void                        SetMinPt(Double_t min)                    { fMinPt  = min ; }
  void                        SetMaxPt(Double_t max)                    { fMaxPt  = max ; }
  void                        SetEtaLimits(Double_t min, Double_t max)  { fMaxEta = max ; fMinEta = min ; }
  void                        SetPhiLimits(Double_t min, Double_t max)  { fMaxPhi = max ; fMinPhi = min ; }
  void                        SetMassHypothesis(Double_t m)             { fMassHypothesis         = m   ; }
  void                        SetClassName(const char *clname);

  /**
   * @brief Set embedding status
   * 
   * Embedding means that the container consists only of tracks from the embedded event.
   * @param[in] b If true the container handles the embedded event
   */
  void                        SetIsEmbedding(Bool_t b)                  { fIsEmbedding = b ; }

  /**
   * @brief Get embedding status
   * 
   * Embedding means that the container consists only of tracks/clusters from the embedded event.
   * @return True if the container handles an embeded event, false otherwise
   */
  Bool_t                      GetIsEmbedding() const                    { return fIsEmbedding; }

  const char*                 GetName()                       const { return fName.Data()               ; }

  /**
   * @brief Set the name of the class of the objets inside the underlying array.
   * @param[in] clname Name of the class of the object inside the underlying array.
   */
  void                        SetName(const char* n)                { fName = n                         ; }

  /**
   * @brief Calculates the relative phi between two angle values and returns it in [-Pi, +Pi] range.
   * @param mphi First angle value
   * @param vphi Second angle value
   * @return Difference between mphi and vphi
   */
  static Double_t             RelativePhi(Double_t ang1, Double_t ang2);

  /**
   * @brief Helper function to calculate the distance between two jets or a jet and a particle
   * @param part1 First particle in the check
   * @param part2 Second particle to compare to
   * @param dist Maximum distance under which partices are considered as "same" in \f$ p_{t} \f$, \f$ \eta \f$ and \f$ \phi \f$
   * @return True if the particles are considered as the same, false otherwise
   */
  static Bool_t               SamePart(const AliVParticle* part1, const AliVParticle* part2, Double_t dist = 1.e-4);
  
  /**
   * @brief Returns the highest bit in the rejection map as reason why the object
   * was rejected.
   * @return Highest bit of rejection reason
   */
  static UShort_t             GetRejectionReasonBitPosition(UInt_t rejectionReason);

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliEmcalIterableContainer   all() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliEmcalIterableContainer   accepted() const;

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliEmcalIterableMomentumContainer   all_momentum() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliEmcalIterableMomentumContainer   accepted_momentum() const;
#endif

 protected:
  /**
   * @brief Handling default Array names. 
   * 
   * Default names might differ
   * based on the input event type. Therefore it is determined
   * for the first time and event is handled.
   *
   * @param ev Input event used to read the input array
   * @return Default array name
   */
  virtual TString             GetDefaultArrayName(const AliVEvent * const ev) const { return ""; }

  /**
   * @brief Retrieve the vertex from the given event. 
   * 
   * It sets fVertex to the vertex of the current event.
   * @param[in] event Input event containing the vertex.
   */
  void                        GetVertexFromEvent(const AliVEvent * event);

  TString                     fName;                    ///< object name
  TString                     fClArrayName;             ///< name of branch
  TString                     fBaseClassName;           ///< name of the base class that this container can handle
  Bool_t                      fIsParticleLevel;         ///< whether or not it is a particle level object collection
  UInt_t                      fBitMap;                  ///< bitmap mask
  Double_t                    fMinPt;                   ///< Min. cut on particle \f$ p_{t} \f$
  Double_t                    fMaxPt;                   ///< Max. cut on particle \f$ p_{t} \f$
  Double_t                    fMaxE;                    ///< Min. cut on particle energy
  Double_t                    fMinE;                    ///< Max. cut on particle energy
  Double_t                    fMinEta;                  ///< Min. cut on particle \f$ \eta \f$
  Double_t                    fMaxEta;                  ///< Max. cut on particle \f$ \eta \f$
  Double_t                    fMinPhi;                  ///< Min. cut on particle \f$ \phi \f$
  Double_t                    fMaxPhi;                  ///< Max. cut on particle \f$ \phi \f$
  Int_t                       fMinMCLabel;              ///< minimum MC label
  Int_t                       fMaxMCLabel;              ///< maximum MC label
  Double_t                    fMassHypothesis;          ///< if < 0 it will use a PID mass when available
  Bool_t                      fIsEmbedding;             ///< if true, this container will connect to an external event
  TClonesArray               *fClArray;                 //!<! Pointer to array in input event
  Int_t                       fCurrentID;               //!<! current ID for automatic loops
  AliNamedArrayI             *fLabelMap;                //!<! Label-Index map
  Double_t                    fVertex[3];               //!<! event vertex array
  TClass                     *fLoadedClass;             //!<! Class of the objects contained in the TClonesArray

 private:
  TString                     fClassName;               ///< name of the class in the TClonesArray

  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,9);
};
#endif
