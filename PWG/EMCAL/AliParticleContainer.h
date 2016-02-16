#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H

class AliVEvent;
class AliVParticle;
class AliVCuts;
class AliTLorentzVector;

#include <TArrayC.h>

#include "AliAODMCParticle.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalContainer.h"

class AliParticleContainer : public AliEmcalContainer {
 public:

  typedef AliEmcalTrackSelection::ETrackFilterType_t ETrackFilterType_t;

  enum ETrackType_t {
    kRejected = -1,
    kUndefined = 0,
    kHybridGlobal = 0,
    kHybridConstrained = 1,
    kHybridConstrainedNoITSrefit = 2,
  };

  AliParticleContainer();
  AliParticleContainer(const char *name, const char *period = "");
  virtual ~AliParticleContainer(){;}

  virtual Bool_t              ApplyParticleCuts(const AliVParticle* vp);
  virtual Bool_t              ApplyKinematicCuts(const AliTLorentzVector& mom);
  virtual Bool_t              AcceptObject(Int_t i)              { return AcceptParticle(i);}
  virtual Bool_t              AcceptObject(const TObject* obj)   { return AcceptParticle(dynamic_cast<const AliVParticle*>(obj));}
  virtual Bool_t              AcceptParticle(const AliVParticle* vp)         ;
  virtual Bool_t              AcceptParticle(Int_t i)                        ;
  Double_t                    GetParticlePtCut()                        const   { return GetMinPt()     ; }
  Double_t                    GetParticleEtaMin()                       const   { return GetMinEta()    ; }
  Double_t                    GetParticleEtaMax()                       const   { return GetMaxEta()    ; }
  Double_t                    GetParticlePhiMin()                       const   { return GetMinPhi()    ; }
  Double_t                    GetParticlePhiMax()                       const   { return GetMaxPhi()    ; }
  void                        SetParticlePtCut(Double_t cut)                    { SetMinPt(cut)         ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { SetEtaLimits(min, max); }
  void                        SetParticlePhiLimits(Double_t min, Double_t max)  { SetPhiLimits(min, max); }
  AliVParticle               *GetLeadingParticle(const char* opt="")         ;
  AliVParticle               *GetParticle(Int_t i=-1)                   const;
  AliVParticle               *GetAcceptParticle(Int_t i=-1)                  ;
  AliVParticle               *GetParticleWithLabel(Int_t lab)           const;
  AliVParticle               *GetAcceptParticleWithLabel(Int_t lab)          ;
  AliVParticle               *GetNextAcceptParticle(Int_t i=-1)              ;
  AliVParticle               *GetNextParticle(Int_t i=-1)                    ;
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVParticle* part, Double_t mass);
  Bool_t                      GetMomentum(TLorentzVector &mom, const AliVParticle* part);
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetNextMomentum(TLorentzVector &mom, Int_t i=-1);
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom, Int_t i=-1);
  Int_t                       GetNParticles()                           const   {return GetNEntries();}
  Int_t                       GetNAcceptedParticles()                   ;
  ETrackFilterType_t          GetTrackFilterType()                      const   { return fTrackFilterType; }
  Char_t                      GetTrackType(Int_t i)                     const   { return i >= 0 && i < fTrackTypes.GetSize() ? fTrackTypes[i] : kUndefined ; }

  void                        SetArray(AliVEvent *event);

  void                        SetClassName(const char *clname);
  void                        SetMinDistanceTPCSectorEdge(Double_t min)         { fMinDistanceTPCSectorEdge = min; }
  void                        SetMCFlag(UInt_t m)                               { fMCFlag          = m ; }
  void                        SelectHIJING(Bool_t s)                            { if (s) fGeneratorIndex = 0; else fGeneratorIndex = -1; }
  void                        SetGeneratorIndex(Short_t i)                      { fGeneratorIndex = i  ; }
  void                        SelectPhysicalPrimaries(Bool_t s)                 { if (s) fMCFlag |=  AliAODMCParticle::kPhysicalPrim ; 
                                                                                  else   fMCFlag &= ~AliAODMCParticle::kPhysicalPrim ; }
  void                        SetCharge(Short_t c)                              { fCharge = c         ; }
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

  Double_t                    fMinDistanceTPCSectorEdge;      // require minimum distance to edge of TPC sector edge
  UInt_t                      fMCFlag;                        // select MC particles with flags
  Short_t                     fGeneratorIndex;                // select MC particles with generator index (default = -1 = switch off selection)
  Short_t                     fCharge;                        // select particles with charge=fCharge
  ETrackFilterType_t          fTrackFilterType;               // track filter type
  TObjArray                  *fListOfCuts;                    // list of track cut objects
  Bool_t                      fSelectionModeAny;              // accept track if any of the cuts is fulfilled
  UInt_t                      fAODFilterBits;                 // track filter bits
  TString                     fTrackCutsPeriod;               // period string used to generate track cuts
  AliEmcalTrackSelection     *fEmcalTrackSelection;           //!track selection object
  TObjArray                  *fFilteredTracks;                //!tracks filtered using fEmcalTrackSelection
  TArrayC                     fTrackTypes;                    //!track types

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  ClassDef(AliParticleContainer,8);

};

#endif

