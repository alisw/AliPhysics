/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALITRACKCONTAINER_H
#define ALITRACKCONTAINER_H

class AliVEvent;
class AliVParticle;
class AliVCuts;
class AliTLorentzVector;

#include <map>

#include <TArrayC.h>

#include "AliVTrack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelResultHybrid.h"
#include "AliParticleContainer.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVTrack, EMCALIterableContainer::operator_star_object<AliVTrack> > AliTrackIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVTrack, EMCALIterableContainer::operator_star_pair<AliVTrack> > AliTrackIterableMomentumContainer;
#endif

/**
 * @class AliTrackContainer
 * @brief Container with name, TClonesArray and cuts for particles
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 */
class AliTrackContainer : public AliParticleContainer {
 public:
 
  /**
   * @class TrackOwnerHandler
   * @brief Unique_ptr implementation for ROOT5 compatibility
   * @ingroup EMCALCOREFW
   * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
   * @since Dec. 20, 2017
   * 
   * This class treats ownership over a shared object in the way a unique_ptr would 
   * do this: Ownership can only be transfered.
   */
  class TrackOwnerHandler : public TObject {
  public:

    /**
     * @brief Default constructor
     */
    TrackOwnerHandler();

    /**
     * @brief Constructor
     * @param managedobject Object to be handled by this manager
     * @param ownership Ownership status
     */
    TrackOwnerHandler(TObjArray *managedobject, Bool_t ownership);

    /**
     * @brief Copy constructor
     * @param other Other ownership handler to copy from
     */
    TrackOwnerHandler(const TrackOwnerHandler &other);

    /**
     * @brief Assignment operator
     * @param other Other ownership handler to assign from
     * @return This ownership handler 
     */
    TrackOwnerHandler &operator=(const TrackOwnerHandler &other);

    /**
     * @brief Destructor
     */
    virtual ~TrackOwnerHandler();

    /**
     * @brief Set object to manage by the ownership handler
     * @param obj Object to be managed
     */
    void SetObject(TObjArray *obj);

    /**
     * @brief Set ownership of the object to be handled
     * @param owner 
     */
    void SetOwner(Bool_t owner);

    /**
     * @brief Access to data the ownership handler manages
     * @return TObjArray* 
     */
    TObjArray *GetData() const { return fManagedObject; }

    /**
     * @brief Check if this handler is the owner of the object it handles
     * @return If true the object is owned by this handler
     */
    Bool_t IsOwner() const { return fOwnership; }

    /**
     * @brief Transfer owership status over this object to other handler
     * @param target Handler that will receive the ownership status
     */
    void TransferOwnershipTo(TrackOwnerHandler &target);

    /**
     * @brief Receive ownership status from other handler
     * @param source Handler from which to receive the ownership status
     */
    void ReceiveOwnershipFrom(TrackOwnerHandler &source);

  private:
    TObjArray          *fManagedObject;            ///< Object managed by the handler
    Bool_t              fOwnership;                 ///< Ownership implementation
    
    /// @cond CLASSIMP
    ClassDef(TrackOwnerHandler, 1);
    /// @endcond
  };

  typedef AliEmcalTrackSelection::ETrackFilterType_t ETrackFilterType_t;

  /// Relates string to the track filter enumeration for %YAML configuration
  static const std::map <std::string, AliEmcalTrackSelection::ETrackFilterType_t> fgkTrackFilterTypeMap; //!<!

  /**
   * @enum ETrackType_t
   * @brief Status of a track after track selection
   */
  enum ETrackType_t {
    kRejected = -1,                  ///< Track rejected
    kUndefined = 0,                  ///< Track status undefined
    kHybridGlobal = 0,               ///< Track selected under the global hybrid track cuts
    kHybridConstrainedTrue = 1,      ///< Track selected under the constrained hybrid track cuts (true constrained)
    kHybridConstrainedFake = 2,      ///< Track selected under the constrained hybrid track cuts (fake constrained)
    kHybridConstrainedNoITSrefit = 3,///< Track selected under the constrained hybrid track cuts without ITS refit
  };

  /**
   * @brief Default constructor.
   */
  AliTrackContainer();

  /**
   * @brief Standard constructor.
   * @param[in] name Name of the container (= name of the array operated on)
   * @param[in] period Name of the period, needed for track cut selection
   */
  AliTrackContainer(const char *name, const char *period = "");

  /**
   * @brief Destructor
   */
  virtual ~AliTrackContainer(){;}

  /**
   * @brief Select track based on track cuts in container
   * @param[in] vp Track to be checked
   * @param[out] rejectionReason Bitmap encoding the reason why the
   * track was rejected. Note: The variable is not set to NULL
   * inside this function before changing its value.
   * @return True if the particle was accepted, false otherwise
   * 
   * Perform track quality selection of the track vp using
   * the track cuts assigned to this container.
   */
  virtual Bool_t              ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const;

  /**
   * @brief Check whether the track at a given index is accepted
   * @param i Index of the track
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the track is accepted, false otherwise
   * 
   * The function is used by AliEmcalContainer in order to select the cluster at a given position
   */
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const                        { return AcceptTrack(i, rejectionReason)        ; }

  /**
   * @brief Check whether the track at a given index is accepted
   * @param obj Object to be checked (must inherit from AliVTrack)
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the track is accepted, false otherwise
   * 
   * The function is used by AliEmcalContainer in order to select the cluster at a given position
   */
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const             { return AcceptTrack(dynamic_cast<const AliVTrack*>(obj), rejectionReason); }

  /**
   * @brief Check whether the track at a given index is accepted
   * @param i Index of the track
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the track is accepted, false otherwise
   */
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const                      { return AcceptTrack(i, rejectionReason); }

  /**
   * @brief Check whether the track at a given index is accepted
   * @param obj Object to be checked (must inherit from AliVTrack)
   * @param[out] rejectionReason Bitmap containing selections criteria which were not passed
   * @return True if the track is accepted, false otherwise
   */
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const       { return AcceptTrack(dynamic_cast<const AliVTrack*>(vp), rejectionReason); }

  /**
   * @brief Get the particle at a given index in the container
   * @param i Index of the particle in the container
   * @return Particle at the index (nullptr if the index exceeds the number of particles in the container)
   */
  virtual AliVParticle       *GetParticle(Int_t i=-1)                const { return GetTrack(i)           ; }

  /**
   * @brief Get the particle at a given index in the container if accepted
   * @param i Index of the particle in the container
   * @return Particle at the index (nullptr if particle is not accepted or the index exceeds the number of particles in the container)
   */
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)          const { return GetAcceptTrack(i)     ; }

  /**
   * @brief Get the next accepted particle in the container
   * @return Next accepted particle in the container (null if there are no more accepted particles)
   * @deprecated Use accepted() in order to iterate over particles 
   * 
   * Internal iterator over accepted tracks. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  virtual AliVParticle       *GetNextAcceptParticle()                      { return GetNextAcceptTrack()  ; }
  
  /**
   * @brief Get the next particle in the container
   * @return Next particle in the container (null if there are no more particles)
   * @deprecated Use all() in order to iterate over particles 
   * 
   * Internal iterator over all tracks. Container must be reset
   * at the beginning of the iteration via Reset().
   */
  virtual AliVParticle       *GetNextParticle()                            { return GetNextTrack()        ; }

  /**
   * @brief Select track based on track cuts in container
   * @param[in] vp Track to be checked
   * @param[in] rejectionReason Bitmap encoding the reason why the
   * track was rejected. Note: The variable is not set to NULL
   * inside this function before changing its value.
   * @return True if the track is accepted, false otherwise
   * 
   * Perform full track selection for the particle vp, consisting
   * of kinematical track selection and track quality
   * cut provided by the cuts assigned to this container
   */
  virtual Bool_t              AcceptTrack(const AliVTrack* vp, UInt_t &rejectionReason)  const;

  /**
   * @brief Select track at a given position based on track cuts in container
   * @param[in] i Index of the track to check
   * @param[in] rejectionReason Bitmap encoding the reason why the
   * track was rejected. Note: The variable is not set to NULL
   * inside this function before changing its value.
   * @return True if the track is accepted, false otherwise
   * 
   * Perform full track selection for the particle in
   * the container stored at position i, consisting
   * of kinematical track selection and track quality
   * cut provided by the cuts assigned to this container
   */
  virtual Bool_t              AcceptTrack(Int_t i, UInt_t &rejectionReason) const;

  /**
   * @brief Get the leading track (highest energetic track) in the container 
   * @param opt Option for the leading track selection 
   * @return Leading track in the container
   */
  virtual AliVTrack          *GetLeadingTrack(const char* opt="")          { return static_cast<AliVTrack*>(GetLeadingParticle(opt)); }

  /**
   * @brief Get track at index in the container
   * @param[in] i Index of the particle in the container
   * @return pointer to particle if particle is accepted, NULL otherwise
   */
  virtual AliVTrack          *GetTrack(Int_t i=-1)                   const;

  /**
   * @brief Get track at index in the container if accepted by the track selection provided
   * @param[in] i Index of the particle in the container
   * @return pointer to particle if particle is accepted, NULL otherwise
   */
  virtual AliVTrack          *GetAcceptTrack(Int_t i=-1)             const;

  /**
   * @brief Get next accepted particle in the container selected using the track cuts provided.
   * @return Next accepted particle (NULL if the end of the array is reached)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::accept_iterator instead
   */
  virtual AliVTrack          *GetNextAcceptTrack()                        ;

  /**
   * @brief Get next particle in the container
   * @return Next track in the container (NULL if end of the container is reached)
   * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::all_iterator instead
   */
   virtual AliVTrack          *GetNextTrack()                              ;

  /**
   * @brief Retrieve the track type using the IndexOf method of the TClonesArray
   * to retrieve the index of the provided track.
   * @param track Track for which the type is requested.
   * @return The track type from fTrackTypes
   */
  Char_t                      GetTrackType(const AliVTrack* track)   const;

  /**
   * @brief Calculate momentum vector for the input track under a user-defined mass hypothesis
   * @param[out] mom Momentum vector to be filled
   * @param[in] track Track from which the momentum information is obtained.
   * @param[in] mass (Optional) Mass hypothesis
   * @return True in case the track is valid, false otherwise
   * 
   * Retrieve momentum information of a track and fill a TLorentzVector
   * with it. In case the optional parameter mass is provided, it is used as mass
   * hypothesis, otherwise the mass hypothesis from the particle itself is used.
   */
  virtual Bool_t              GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track, Double_t mass) const;

  /**
   * @brief Fills a TLorentzVector with the momentum information of the track provided
   * under a global mass hypothesis.
   * @param[out] mom Momentum vector of the particle provided
   * @param[in] track Track from which to obtain the momentum information
   * @return Always true
   */
  virtual Bool_t              GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track) const;

  /**
   * @brief Get momentum vector of a track a given position in the container
   * @param[out] mom Momentum vector of the \f$ i^{th} \f$ particle in the array
   * @param[in] i Index of the particle to check
   * @return True if the request was successful, false otherwise
   * 
   * Fills a TLorentzVector with the momentum information of the
   * \f$ i^{th} \f$ particle in the container, using a global
   * mass hypothesis. In case the provided index is out of
   * range, false is returned as return value.
   */
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Get momentum vector of the track at given index in case the track was accepted
   * @param[out] mom Momentum vector of the accepted particle
   * @param[in] i Index to check
   * @return True if the request was successfull, false otherwise
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * \f$ i^{th} \f$ accepted particle in the container, using a
   * global mass hypothesis. In case the provided index is out of
   * range, or the particle under the index is not accepted, false
   * is returned as return value.
   */
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;

  /**
   * @brief Get momentum vector of the next track in the container
   * @param[out] mom Momentum vector of the next particle
   * @return True if the request was successful, false otherwise
   * @deprecated Old style iterator - use all_iterator instead
   * 
   * Fills a TLorentzVector with the momentum information of the
   * next particle in the container, using a global mass hypothesis.
   * In case the iterator reached the end of the array, false
   * is returned as return value.
   */
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);

  /**
   * @brief Get momentum vector of the next accepted track
   * @param[out] mom Momentum vector of the next particle in the array
   * @return True if the request was successfull, false (no more entries) otherwise
   * @deprecated Old style iterator - use accept_iterator instead
   * 
   * Fills a TLorentzVector with the monentum infomation of the
   * next accepted particle in the container, using a global
   * mass hypothesis. In case the iteration reached the end of
   * the array, false is returned as return value.
   */
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);

  /**
   * @brief Get the number of tracks handled by the container
   * @return Number of tracks in the container
   */
  Int_t                       GetNTracks()                              const   { return GetNParticles()         ; }

  /**
   * @brief Get the number of accepted tracks in the container
   * @return Number of tracks in the container 
   */
  Int_t                       GetNAcceptedTracks()                              { return GetNAcceptedParticles() ; }

  ETrackFilterType_t          GetTrackFilterType()                      const   { return fTrackFilterType; }
  Char_t                      GetTrackType(Int_t i)                     const   { return i >= 0 && i < fTrackTypes.GetSize() ? fTrackTypes[i] : (Char_t)kUndefined ; }

  /**
   * @brief Get array from event 
   * @param[in] event Event from which the data is read
   * 
   * Also creating the virtual track selection
   * for the period provided in the constructor.
   */
  void                        SetArray(const AliVEvent *event);

  void                        SetTrackFilterType(ETrackFilterType_t f)          { fTrackFilterType = f; }
  void                        SetFilterHybridTracks(Bool_t f)                   { if (f) fTrackFilterType = AliEmcalTrackSelection::kHybridTracks; else fTrackFilterType = AliEmcalTrackSelection::kNoTrackFilter; }   // legacy method
  void                        SetITSHybridTrackDistinction(Bool_t doUse)        { fITSHybridTrackDistinction = doUse; }

  void                        SetTrackCutsPeriod(const char* period)            { fTrackCutsPeriod = period; }

  /**
   * @brief Add new track cuts to the container.
   * @param[in] cuts Cuts to be  added
   */
  void                        AddTrackCuts(AliVCuts *cuts);

  /**
   * @brief Get number of track cut objects assigned to this container.
   * @return Number of track cut objects
   */
  Int_t                       GetNumberOfCutObjects() const;

  /**
   * @brief Get the cut object at index (icut) assigned to this container.
   * @param[in] icut Index of the cut in the container
   * @return Cut object at the index if existing, NULL otherwise
   */
  AliVCuts                   *GetTrackCuts(Int_t icut);
  void                        SetAODFilterBits(UInt_t bits)                     { fAODFilterBits   = bits  ; }
  void                        AddAODFilterBit(UInt_t bit)                       { fAODFilterBits  |= bit   ; }
  UInt_t                      GetAODFilterBits()                          const { return fAODFilterBits    ; }
  Bool_t                      IsHybridTrackSelection() const;

  void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }
  void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

  /**
   * @brief Preparation for the next event
   * @param event Next event to be processed
   * 
   * Run the track selection of all bit and store the pointers to
   * selected tracks in a separate array.
   */
  void                        NextEvent(const AliVEvent* event);

  static void                 SetDefTrackCutsPeriod(const char* period)       { fgDefTrackCutsPeriod = period; }
  static TString              GetDefTrackCutsPeriod()                         { return fgDefTrackCutsPeriod  ; }

  /**
   * @brief Get the title of the track container
   * @return Title of the container
   * 
   * Build title of the container consisting of the container name
   * and a string encoding the minimum \f$ p_{t} \f$ cut applied
   * in the kinematic track selection.
   */
  const char*                 GetTitle() const;

  /**
   * @brief Test function checking whether the entries in the track array are the same as in the input array
   * @return If true entries are consistent 
   */
  Bool_t CheckArrayConsistency() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliTrackIterableContainer      all() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliTrackIterableContainer      accepted() const;

  /**
   * @brief Create an iterable container interface over all objects in the
   * EMCAL container.
   * @return iterable container over all objects in the EMCAL container
   */
  const AliTrackIterableMomentumContainer      all_momentum() const;

  /**
   * @brief Create an iterable container interface over accepted objects in the
   * EMCAL container.
   * @return iterable container over accepted objects in the EMCAL container
   */
  const AliTrackIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  /**
   * @brief Create default array name for the track container. 
   * @param[in] ev Input event, used for data type selection
   * @return Appropriate default array name
   * 
   * The default array name will be
   * - *tracks* in case of AOD event
   * - *Tracks* in case of ESD event
   */
  virtual TString             GetDefaultArrayName(const AliVEvent * const ev) const;

  PWG::EMCAL::AliEmcalTrackSelResultHybrid::HybridType_t  GetHybridDefinition(const PWG::EMCAL::AliEmcalTrackSelResultPtr &selectionResult) const;

  static TString              fgDefTrackCutsPeriod;           //!<! default period string used to generate track cuts

  ETrackFilterType_t          fTrackFilterType;               ///< track filter type
  TObjArray                  *fListOfCuts;                    ///< list of track cut objects
  Bool_t                      fSelectionModeAny;              ///< accept track if any of the cuts is fulfilled
  Bool_t                      fITSHybridTrackDistinction;     ///< Distinct hybrid tracks via SPD information
  UInt_t                      fAODFilterBits;                 ///< track filter bits
  TString                     fTrackCutsPeriod;               ///< period string used to generate track cuts
  AliEmcalTrackSelection     *fEmcalTrackSelection;  //!<! track selection object
  TrackOwnerHandler           fFilteredTracks;                //!<! tracks filtered using fEmcalTrackSelection
  TArrayC                     fTrackTypes;                    //!<! track types

 private:
  AliTrackContainer(const AliTrackContainer& obj); // copy constructor
  AliTrackContainer& operator=(const AliTrackContainer& other); // assignment

  ClassDef(AliTrackContainer,1);
};

#endif

