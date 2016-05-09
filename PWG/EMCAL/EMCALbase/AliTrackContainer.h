#ifndef ALITRACKCONTAINER_H
#define ALITRACKCONTAINER_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliVEvent;
class AliVParticle;
class AliVCuts;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliVTrack.h"
#include "AliEmcalTrackSelection.h"
#include "AliParticleContainer.h"

typedef AliEmcalIterableContainerT<AliVTrack> AliTrackIterableContainer;

/**
 * @class AliTrackContainer
 * @brief Container with name, TClonesArray and cuts for particles
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 */
class AliTrackContainer : public AliParticleContainer {
 public:

  typedef AliEmcalTrackSelection::ETrackFilterType_t ETrackFilterType_t;

  /**
   * @enum ETrackType_t
   * @brief Status of a track after track selection
   */
  enum ETrackType_t {
    kRejected = -1,                  ///< Track rejected
    kUndefined = 0,                  ///< Track status undefined
    kHybridGlobal = 0,               ///< Track selected under the global hybrid track cuts
    kHybridConstrained = 1,          ///< Track selected under the constrained hybrid track cuts
    kHybridConstrainedNoITSrefit = 2,///< Track selected under the constrained hybrid track cuts without ITS refit
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
  Char_t                      GetTrackType(Int_t i)                     const   { return i >= 0 && i < fTrackTypes.GetSize() ? fTrackTypes[i] : kUndefined ; }

  void                        SetArray(const AliVEvent *event);

  void                        SetTrackFilterType(ETrackFilterType_t f)          { fTrackFilterType = f; }
  void                        SetFilterHybridTracks(Bool_t f)                   { if (f) fTrackFilterType = AliEmcalTrackSelection::kHybridTracks; else fTrackFilterType = AliEmcalTrackSelection::kNoTrackFilter; }   // legacy method

  void                        SetTrackCutsPeriod(const char* period)            { fTrackCutsPeriod = period; }
  void                        AddTrackCuts(AliVCuts *cuts);
  Int_t                       GetNumberOfCutObjects() const;
  AliVCuts                   *GetTrackCuts(Int_t icut);
  void                        SetAODFilterBits(UInt_t bits)                     { fAODFilterBits   = bits  ; }
  void                        AddAODFilterBit(UInt_t bit)                       { fAODFilterBits  |= bit   ; }
  UInt_t                      GetAODFilterBits()                          const { return fAODFilterBits    ; }

  void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }
  void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

  void                        NextEvent();

  static void                 SetDefTrackCutsPeriod(const char* period)       { fgDefTrackCutsPeriod = period; }
  static TString              GetDefTrackCutsPeriod()                         { return fgDefTrackCutsPeriod  ; }

  const char*                 GetTitle() const;

  const AliTrackIterableContainer      all() const;
  const AliTrackIterableContainer      accepted() const;

  AliTrackIterableContainer::iterator  accept_begin()  const { return accepted().begin()   ; }
  AliTrackIterableContainer::iterator  accept_end()    const { return accepted().end()     ; }
  AliTrackIterableContainer::iterator  accept_rbegin() const { return accepted().rbegin()  ; }
  AliTrackIterableContainer::iterator  accept_rend()   const { return accepted().rend()    ; }

  AliTrackIterableContainer::iterator  begin()         const { return all().begin()        ; }
  AliTrackIterableContainer::iterator  end()           const { return all().end()          ; }
  AliTrackIterableContainer::iterator  rbegin()        const { return all().rbegin()       ; }
  AliTrackIterableContainer::iterator  rend()          const { return all().rend()         ; }

 protected:
  static TString              fgDefTrackCutsPeriod;           //!<! default period string used to generate track cuts

  ETrackFilterType_t          fTrackFilterType;               ///< track filter type
  TObjArray                  *fListOfCuts;                    ///< list of track cut objects
  Bool_t                      fSelectionModeAny;              ///< accept track if any of the cuts is fulfilled
  UInt_t                      fAODFilterBits;                 ///< track filter bits
  TString                     fTrackCutsPeriod;               ///< period string used to generate track cuts
  AliEmcalTrackSelection     *fEmcalTrackSelection;           //!<! track selection object
  TObjArray                  *fFilteredTracks;                //!<! tracks filtered using fEmcalTrackSelection
  TArrayC                     fTrackTypes;                    //!<! track types

 private:
  AliTrackContainer(const AliTrackContainer& obj); // copy constructor
  AliTrackContainer& operator=(const AliTrackContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliTrackContainer,1);
  /// \endcond
};

#endif

