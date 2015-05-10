#ifndef ALIEMCALMCTRAKCSELECTOR_H
#define ALIEMCALMCTRAKCSELECTOR_H

class TClonesArray;
class TString;
class AliVEvent;
class AliMCEvent;
class AliNamedArrayI;
class AliAODMCParticle;

#include "AliAnalysisTaskSE.h"

class AliEmcalMCTrackSelector : public AliAnalysisTaskSE {
 public:
  AliEmcalMCTrackSelector();
  AliEmcalMCTrackSelector(const char *name);
  virtual ~AliEmcalMCTrackSelector();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  void SetOnlyPhysPrim(Bool_t s)                        { fOnlyPhysPrim     = s    ; }  
  void SetChargedMC(Bool_t c = kTRUE)                   { fChargedMC        = c    ; }
  void SetEtaMax(Double_t e)                            { fEtaMax           = e    ; }
  void SetRejectNK(Bool_t r = kTRUE)                    { fRejectNK         = r    ; }
  void SetOnlyHIJING(Bool_t s)                          { fOnlyHIJING       = s    ; }
  void SetParticlesOutName(const char *name)            { fParticlesOutName = name ; }
  void SetSpecialPDG(Int_t pdg)                         { fSpecialPDG       = pdg  ; }
  
 protected:
  void                      ConvertMCParticles();    // for ESD analysis
  void                      CopyMCParticles();       // for AOD analysis

  TString                   fParticlesOutName;     // name of output particle array
  Bool_t                    fOnlyPhysPrim;         // true = only physical primary particles
  Bool_t                    fRejectNK;             // true = reject K_0^L and neutrons
  Bool_t                    fChargedMC;            // true = only charged particles
  Bool_t                    fOnlyHIJING;           // true = only HIJING particles
  Double_t                  fEtaMax;               // maximum eta to accept particles
  Int_t                     fSpecialPDG;           // include particles with this PDG code even if they are not primary particles (and exclude their daughters)
  TString                   fParticlesMapName;     //!name of the particle map
  Bool_t                    fInit;                 //!true = task initialized
  TClonesArray             *fParticlesIn;          //!particle array in (AOD)
  TClonesArray             *fParticlesOut;         //!particle array out
  AliNamedArrayI           *fParticlesMap;         //!particle index/label
  AliVEvent                *fEvent;                //!event
  AliMCEvent               *fMC;                   //!MC event (ESD)
  Bool_t                    fIsESD;                //!ESD or AOD analysis
  Bool_t                    fDisabled;             //!Disable task if a problem occurs at initialization

 private:
  AliEmcalMCTrackSelector(const AliEmcalMCTrackSelector&);            // not implemented
  AliEmcalMCTrackSelector &operator=(const AliEmcalMCTrackSelector&); // not implemented

  ClassDef(AliEmcalMCTrackSelector, 4); // Task to select particle in MC events
};
#endif
