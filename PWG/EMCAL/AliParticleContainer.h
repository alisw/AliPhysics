#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H

class AliVEvent;
class AliVParticle;
class AliVCuts;

#include "AliAODMCParticle.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalContainer.h"

class AliParticleContainer : public AliEmcalContainer {
 public:

  typedef AliEmcalTrackSelection::ETrackFilterType_t ETrackFilterType_t;

  enum ETrackType_t {
    kUndefined = 0,
    kHybridGlobal = 0,
    kHybridConstrained = 1,
    kHybridConstrainedNoITSrefit = 2,
  };

  AliParticleContainer();
  AliParticleContainer(const char *name, const char *period = "");
  virtual ~AliParticleContainer(){;}

  Bool_t                      AcceptParticle(AliVParticle         *vp)       ;
  Double_t                    GetParticlePtCut()                        const   { return fParticlePtCut  ; }
  Double_t                    GetParticleEtaMin()                       const   { return fParticleMinEta ; }
  Double_t                    GetParticleEtaMax()                       const   { return fParticleMaxEta ; }
  Double_t                    GetParticlePhiMin()                       const   { return fParticleMinPhi ; }
  Double_t                    GetParticlePhiMax()                       const   { return fParticleMaxPhi ; }
  AliVParticle               *GetLeadingParticle(const char* opt="")         ;
  AliVParticle               *GetParticle(Int_t i)                      const;
  AliVParticle               *GetAcceptParticle(Int_t i)                     ;
  AliVParticle               *GetParticleWithLabel(Int_t lab)           const;
  AliVParticle               *GetAcceptParticleWithLabel(Int_t lab)          ;
  AliVParticle               *GetNextAcceptParticle(Int_t i=-1)              ;
  AliVParticle               *GetNextParticle(Int_t i=-1)                    ;
  void                        GetMomentum(TLorentzVector &mom, Int_t i) const;
  Int_t                       GetNParticles()                           const   {return GetNEntries();}
  Int_t                       GetNAcceptedParticles()                   ;
  ETrackFilterType_t          GetTrackFilterType()                      const   { return fTrackFilterType; }
  ETrackType_t                GetTrackType()                            const   { return fTrackType      ; }

  void                        SetArray(AliVEvent *event);

  void                        SetClassName(const char *clname);
  void                        SetMCTrackBitMap(UInt_t m)                        { fMCTrackBitMap   = m ; }
  void                        SetMinMCLabel(Int_t s)                            { fMinMCLabel      = s ; }
  void                        SetMinMCLabelAccept(Int_t s)                      { fMinMCLabelAccept= s ; }
  void                        SetParticlePtCut(Double_t cut)                    { fParticlePtCut = cut ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { fParticleMaxEta = max ; fParticleMinEta = min ; }
  void                        SetParticlePhiLimits(Double_t min, Double_t max, Double_t offset=0.)  { fParticleMaxPhi = max ; fParticleMinPhi = min ; fPhiOffset = offset;}
  void                        SetMinDistanceTPCSectorEdge(Double_t min)         { fMinDistanceTPCSectorEdge = min; }
  void                        SetTrackBitMap(UInt_t m)                          { fTrackBitMap     = m ; }
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

  void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }
  void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

  static void                 SetDefTrackCutsPeriod(const char* period)       { fgDefTrackCutsPeriod = period; }
  static TString              GetDefTrackCutsPeriod()                         { return fgDefTrackCutsPeriod  ; }

 protected:
  static TString              fgDefTrackCutsPeriod;           //!default period string used to generate track cuts

  Double_t                    fParticlePtCut;                 // cut on particle pt
  Double_t                    fParticleMinEta;                // cut on particle eta
  Double_t                    fParticleMaxEta;                // cut on particle eta
  Double_t                    fParticleMinPhi;                // cut on particle phi
  Double_t                    fParticleMaxPhi;                // cut on particle phi
  Double_t                    fPhiOffset;                     // phi offset
  Double_t                    fMinDistanceTPCSectorEdge;      // require minimum distance to edge of TPC sector edge
  UInt_t                      fTrackBitMap;                   // bit map of accepted tracks (non MC)
  UInt_t                      fMCTrackBitMap;                 // bit map of accepted MC tracks
  Int_t                       fMinMCLabel;                    // minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                       fMinMCLabelAccept;              // minimum MC label value to accept particle
  UInt_t                      fMCFlag;                        // select MC particles with flags
  Short_t                     fGeneratorIndex;                // select MC particles with generator index (default = -1 = switch off selection)
  Short_t                     fCharge;                        // select particles with charge=fCharge
  ETrackFilterType_t          fTrackFilterType;               // track filter type
  TObjArray                  *fListOfCuts;                    // list of track cut objects
  Bool_t                      fSelectionModeAny;              // accept track if any of the cuts is fulfilled
  UInt_t                      fAODFilterBits;                 // track filter bits
  TString                     fTrackCutsPeriod;               // period string used to generate track cuts
  AliEmcalTrackSelection     *fEmcalTrackSelection;           //!track selection object
  ETrackType_t                fTrackType;                     //!last accepted track type

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  ClassDef(AliParticleContainer,8);

};

#endif

