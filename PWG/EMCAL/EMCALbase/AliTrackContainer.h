#ifndef ALITRACKCONTAINER_H
#define ALITRACKCONTAINER_H

class AliVEvent;
class AliVParticle;
class AliVCuts;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliVTrack.h"
#include "AliEmcalTrackSelection.h"
#include "AliParticleContainer.h"

class AliTrackContainer : public AliParticleContainer {
 public:

  typedef AliEmcalTrackSelection::ETrackFilterType_t ETrackFilterType_t;

  enum ETrackType_t {
    kRejected = -1,
    kUndefined = 0,
    kHybridGlobal = 0,
    kHybridConstrained = 1,
    kHybridConstrainedNoITSrefit = 2,
  };

  AliTrackContainer();
  AliTrackContainer(const char *name, const char *period = "");
  virtual ~AliTrackContainer(){;}

  virtual Bool_t              ApplyTrackCuts(const AliVTrack* vp);
  virtual Bool_t              AcceptObject(Int_t i)                        { return AcceptTrack(i)        ; }
  virtual Bool_t              AcceptObject(const TObject* obj)             { return AcceptTrack(dynamic_cast<const AliVTrack*>(obj)) ; }
  virtual Bool_t              AcceptParticle(Int_t i)                      { return AcceptTrack(i)        ; }
  virtual Bool_t              AcceptParticle(const AliVParticle* vp)       { return AcceptTrack(dynamic_cast<const AliVTrack*>(vp))  ; }
  virtual AliVParticle       *GetParticle(Int_t i=-1)                const { return GetTrack(i)           ; }
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)                { return GetAcceptTrack(i)     ; }
  virtual AliVParticle       *GetNextAcceptParticle()                      { return GetNextAcceptTrack()  ; }
  virtual AliVParticle       *GetNextParticle()                            { return GetNextTrack()        ; }
  virtual Bool_t              AcceptTrack(const AliVTrack* vp)            ;
  virtual Bool_t              AcceptTrack(Int_t i)                        ;
  virtual AliVTrack          *GetLeadingTrack(const char* opt="")          { return static_cast<AliVTrack*>(GetLeadingParticle(opt)); }
  virtual AliVTrack          *GetTrack(Int_t i=-1)                   const;
  virtual AliVTrack          *GetAcceptTrack(Int_t i=-1)                  ;
  virtual AliVTrack          *GetNextAcceptTrack()                        ;
  virtual AliVTrack          *GetNextTrack()                              ;
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliVTrack* part, Double_t mass);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, const AliVTrack* part);
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);
  Int_t                       GetNTracks()                              const   { return GetNParticles()         ; }
  Int_t                       GetNAcceptedTracks()                              { return GetNAcceptedParticles() ; }
  ETrackFilterType_t          GetTrackFilterType()                      const   { return fTrackFilterType; }
  Char_t                      GetTrackType(Int_t i)                     const   { return i >= 0 && i < fTrackTypes.GetSize() ? fTrackTypes[i] : kUndefined ; }

  void                        SetArray(AliVEvent *event);

  void                        SetClassName(const char *clname);

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

 protected:
  static TString              fgDefTrackCutsPeriod;           //!default period string used to generate track cuts

  ETrackFilterType_t          fTrackFilterType;               // track filter type
  TObjArray                  *fListOfCuts;                    // list of track cut objects
  Bool_t                      fSelectionModeAny;              // accept track if any of the cuts is fulfilled
  UInt_t                      fAODFilterBits;                 // track filter bits
  TString                     fTrackCutsPeriod;               // period string used to generate track cuts
  AliEmcalTrackSelection     *fEmcalTrackSelection;           //!track selection object
  TObjArray                  *fFilteredTracks;                //!tracks filtered using fEmcalTrackSelection
  TArrayC                     fTrackTypes;                    //!track types

 private:
  AliTrackContainer(const AliTrackContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliTrackContainer& other); // assignment

  ClassDef(AliTrackContainer,1);
};

#endif

