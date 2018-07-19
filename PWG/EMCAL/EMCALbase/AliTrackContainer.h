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
    TrackOwnerHandler();
    TrackOwnerHandler(TObjArray *managedobject, Bool_t ownership);
    TrackOwnerHandler(const TrackOwnerHandler &other);
    TrackOwnerHandler &operator=(const TrackOwnerHandler &other);
    virtual ~TrackOwnerHandler();

    void SetObject(TObjArray *obj);
    void SetOwner(Bool_t owner);

    TObjArray *GetData() const { return fManagedObject; }
    Bool_t IsOwner() const { return fOwnership; }

    void TransferOwnershipTo(TrackOwnerHandler &target);
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

  AliTrackContainer();
  AliTrackContainer(const char *name, const char *period = "");
  virtual ~AliTrackContainer(){;}

  virtual Bool_t              ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const;
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const                        { return AcceptTrack(i, rejectionReason)        ; }
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const             { return AcceptTrack(dynamic_cast<const AliVTrack*>(obj), rejectionReason); }
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const                      { return AcceptTrack(i, rejectionReason); }
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const       { return AcceptTrack(dynamic_cast<const AliVTrack*>(vp), rejectionReason); }
  virtual AliVParticle       *GetParticle(Int_t i=-1)                const { return GetTrack(i)           ; }
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)          const { return GetAcceptTrack(i)     ; }
  virtual AliVParticle       *GetNextAcceptParticle()                      { return GetNextAcceptTrack()  ; }
  virtual AliVParticle       *GetNextParticle()                            { return GetNextTrack()        ; }
  virtual Bool_t              AcceptTrack(const AliVTrack* vp, UInt_t &rejectionReason)  const;
  virtual Bool_t              AcceptTrack(Int_t i, UInt_t &rejectionReason) const;
  virtual AliVTrack          *GetLeadingTrack(const char* opt="")          { return static_cast<AliVTrack*>(GetLeadingParticle(opt)); }
  virtual AliVTrack          *GetTrack(Int_t i=-1)                   const;
  virtual AliVTrack          *GetAcceptTrack(Int_t i=-1)             const;
  virtual AliVTrack          *GetNextAcceptTrack()                        ;
  virtual AliVTrack          *GetNextTrack()                              ;
  Char_t                      GetTrackType(const AliVTrack* track)   const;
  virtual Bool_t              GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track, Double_t mass) const;
  virtual Bool_t              GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track) const;
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);
  Int_t                       GetNTracks()                              const   { return GetNParticles()         ; }
  Int_t                       GetNAcceptedTracks()                              { return GetNAcceptedParticles() ; }
  ETrackFilterType_t          GetTrackFilterType()                      const   { return fTrackFilterType; }
  Char_t                      GetTrackType(Int_t i)                     const   { return i >= 0 && i < fTrackTypes.GetSize() ? fTrackTypes[i] : (Char_t)kUndefined ; }

  void                        SetArray(const AliVEvent *event);

  void                        SetTrackFilterType(ETrackFilterType_t f)          { fTrackFilterType = f; }
  void                        SetFilterHybridTracks(Bool_t f)                   { if (f) fTrackFilterType = AliEmcalTrackSelection::kHybridTracks; else fTrackFilterType = AliEmcalTrackSelection::kNoTrackFilter; }   // legacy method
  void                        SetITSHybridTrackDistinction(Bool_t doUse)        { fITSHybridTrackDistinction = doUse; }

  void                        SetTrackCutsPeriod(const char* period)            { fTrackCutsPeriod = period; }
  void                        AddTrackCuts(AliVCuts *cuts);
  Int_t                       GetNumberOfCutObjects() const;
  AliVCuts                   *GetTrackCuts(Int_t icut);
  void                        SetAODFilterBits(UInt_t bits)                     { fAODFilterBits   = bits  ; }
  void                        AddAODFilterBit(UInt_t bit)                       { fAODFilterBits  |= bit   ; }
  UInt_t                      GetAODFilterBits()                          const { return fAODFilterBits    ; }
  Bool_t                      IsHybridTrackSelection() const;

  void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }
  void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

  void                        NextEvent(const AliVEvent* event);

  static void                 SetDefTrackCutsPeriod(const char* period)       { fgDefTrackCutsPeriod = period; }
  static TString              GetDefTrackCutsPeriod()                         { return fgDefTrackCutsPeriod  ; }

  const char*                 GetTitle() const;

  /**
   * @brief Test function checking whether the entries in the track array are the same as in the input array
   * 
   * @return Bool_t CheckArrayConsistency 
   */
  Bool_t CheckArrayConsistency() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))
  const AliTrackIterableContainer      all() const;
  const AliTrackIterableContainer      accepted() const;

  const AliTrackIterableMomentumContainer      all_momentum() const;
  const AliTrackIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  /**
   * Create default array name for the track container. The
   * default array name will be
   * - *tracks* in case of AOD event
   * - *Tracks* in case of ESD event
   * @param[in] ev Input event, used for data type selection
   * @return Appropriate default array name
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

  /// \cond CLASSIMP
  ClassDef(AliTrackContainer,1);
  /// \endcond
};

#endif

