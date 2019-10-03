#ifndef ALIDIELECTRONMC_H
#define ALIDIELECTRONMC_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   #
//#              Class AliDielectronMC                #
//#       Cut Class for Jpsi->e+e- analysis           #
//#                                                   #
//#   by WooJin J. Park, GSI / W.J.Park@gsi.de        #
//#                                                   #
//#####################################################

#ifndef ROOT_TObject
#include <TObject.h>
#include <TMCProcess.h>
#endif

class AliESDEvent;
class AliHFEpid;
class AliMCEvent;
class AliESDtrack;
class AliAODTrack;
class TParticle;
class AliMCParticle;
class AliAODMCParticle;
class AliAODMCHeader;
class AliGenHijingEventHeader;

#include "AliDielectronSignalMC.h"
#include "AliDielectronPair.h"

#include <iostream>

class AliDielectronMC : public TObject{

public:
  enum AnalysisType {kUNSET=0, kESD, kAOD};

  AliDielectronMC(AnalysisType type=kUNSET);
  virtual ~AliDielectronMC();

  void SetHasMC(Bool_t hasMC) { fHasMC=hasMC; }
  Bool_t HasMC() const { return fHasMC; }

  void SetCheckHF(Bool_t checkHF) { fCheckHF=checkHF; }
  Bool_t CheckHF() const { return fCheckHF; }

  static AliDielectronMC* Instance();

  void Initialize();                              // initialization
  Int_t GetNMCTracks();                                     // return number of generated tracks
  Int_t GetNMCTracksFromStack();                            // return number of generated tracks from MC event
  Int_t GetNPrimary() const;                                      // return number of primary tracks
  Int_t GetNPrimaryFromStack();                                   // return number of primary tracks from MC event
  Int_t GetMCPID(const AliESDtrack* _track);                      // return MC PID
  Int_t GetMCPID(const AliAODTrack* _track);                      // return MC PID for AODtrack
  Int_t GetMCPIDFromStack(const AliESDtrack* _track);             // return MC PID from MC event
  Int_t GetMotherPDG(const AliESDtrack* _track);                  // return mother PID from the MC event
  Int_t GetMotherPDG(const AliAODTrack* _track);                  // return mother PID from the MC event
  Int_t GetMotherPDG(const AliMCParticle* _track);                  // return mother PID from the MC event
  Int_t GetMotherPDG(const AliAODMCParticle* _track);                  // return mother PID from the MC event
  //Int_t GetMotherPDGFromStack(const AliESDtrack* _track);         // return mother PID from the MC event
  Int_t GetMCProcess(const AliESDtrack* _track);                  // return process number
  //Int_t GetMCProcessFromStack(const AliESDtrack* _track);         // return process number
  Int_t GetMCProcessMother(const AliESDtrack* _track);            // return process number of the mother track
  //Int_t GetMCProcessMotherFromStack(const AliESDtrack* _track);   // return process number of the mother track

  Bool_t ConnectMCEvent();

  Bool_t IsMotherPdg(const AliDielectronPair* pair, Int_t pdgMother);
  Bool_t IsMotherPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Bool_t IsMCMotherToEE(const AliVParticle *particle, Int_t pdgMother);
  Bool_t IsMCTruth(const AliDielectronPair* pair, const AliDielectronSignalMC* signalMC) const;
  Bool_t IsMCTruth(const AliDielectronPair* pair1, const AliDielectronPair* pair2, const AliDielectronSignalMC* signalMC1, const AliDielectronSignalMC* signalMC2) const;
  Bool_t IsMCTruth(AliVParticle* mcD1, AliVParticle* mcD2, const AliDielectronSignalMC* signalMC) const;
  Bool_t IsMCTruth(AliVParticle* mcD1, AliVParticle* mcD2, AliVParticle* mcD3, AliVParticle* mcD4, const AliDielectronSignalMC* signalMC1, const AliDielectronSignalMC* signalMC2) const;
  Bool_t IsMCTruth(Int_t label, AliDielectronSignalMC* signalMC, Int_t branch) const;
  Int_t GetMothersLabel(Int_t daughterLabel) const;
  Int_t GetFirstMothersLabelInChain(Int_t daughterLabel) const;
  Int_t GetMinLabelParticlePrimaryAround(Int_t label) const;
  Int_t GetMaxLabelParticlePrimaryAround(Int_t label) const;
  Int_t GetPdgFromLabel(Int_t label) const;
  Int_t GetHFProcess(Int_t label);

  Bool_t IsPrimary(Int_t label) const;
  Bool_t IsPhysicalPrimary(Int_t label) const;  // checks if a particle is physical primary
  Bool_t IsSecondary(Int_t label) const;
  Bool_t CheckGEANTProcess(Int_t label, TMCProcess process) const;
  Bool_t IsSecondaryFromWeakDecay(Int_t label) const;
  Bool_t IsSecondaryFromMaterial(Int_t label) const;
  //Bool_t IsEleFromInjectedSignal(Int_t label) const;
  Bool_t IsFromBGEvent(Int_t label) const;
  Bool_t CheckHijingHeader() const;

  Bool_t HaveSameMother(const AliDielectronPair *pair) const;

  Int_t GetLabelMotherWithPdg(const AliDielectronPair* pair, Int_t pdgMother);
  Int_t GetLabelMotherWithPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);

//   AliVParticle* GetMCTrackFromMCEvent(const AliVParticle *track);   // return MC track directly from MC event
  AliVParticle* GetMCTrackFromMCEvent(Int_t label) const;           // return MC track directly from MC event
  TParticle* GetMCTrackFromStack(const AliESDtrack* _track);        // return MC track from MC event
  AliVParticle* GetMCTrack(const AliVParticle* _track); // return MC track from MC event
  AliMCParticle* GetMCTrack(const AliESDtrack* _track);             // return MC track associated with reco track
  AliAODMCParticle* GetMCTrack( const AliAODTrack* _track);          // return MC track associated with reco AOD track

  //TParticle* GetMCTrackMotherFromStack(const AliESDtrack* _track);  // return MC mother track from MC event
  AliMCParticle* GetMCTrackMother(const AliESDtrack* _track);       // return MC mother track from MC event
  AliAODMCParticle* GetMCTrackMother(const AliAODTrack* _track);       // return MC mother fot track AODTrack
  AliMCParticle* GetMCTrackMother(const AliMCParticle* _particle);       // return MC mother track from stack
  AliAODMCParticle* GetMCTrackMother(const AliAODMCParticle* _particle);       // return MC mother track from stack

  Int_t NumberOfDaughters(const AliESDtrack* track);                 // return number of daughters
  Int_t NumberOfDaughters(const AliAODTrack* track);                 // return number of daughters for AOD track
  Int_t NumberOfDaughters(const AliMCParticle* particle);                 // return number of daughters
  Int_t NumberOfDaughters(const AliAODMCParticle* particle);                 // return number of daughters

  void GetDaughters(const TObject *mother, AliVParticle* &d1, AliVParticle* &d2);
  Int_t IsJpsiPrimary(const AliDielectronPair * pair);
  Int_t IsJpsiPrimary(const AliVParticle * pair);
  Bool_t CheckParticleSource(Int_t label, AliDielectronSignalMC::ESource source) const;
  Bool_t CheckParticleSource(const AliAODMCParticle *mcPart, AliDielectronSignalMC::ESource source) const;

  Bool_t GetPrimaryVertex(Double_t &primVtxX, Double_t &primVtxY, Double_t &primVtxZ);

  AliMCEvent* GetMCEvent() { return fMCEvent; }         // return the AliMCEvent

private:
  AliMCEvent    *fMCEvent;  // MC event object
  AliAODMCHeader *fAODMCHeader; // MC AOD Header

  mutable AliGenHijingEventHeader *fGenCocktailHeader; //! Generated cocktail header

  AnalysisType fAnaType;    // Analysis type
  Bool_t fHasMC;            // Do we have an MC handler?
  Bool_t fCheckHF;          // Do we look for HF correlated pairs?

  std::map<Int_t,Int_t> fhfproc; // quark label and HF process

  mutable Int_t  fHasHijingHeader;  //! //mutable needed to change it in a const function.

  static AliDielectronMC* fgInstance; //! singleton pointer
  TClonesArray* fMcArray; //mcArray for AOD MC particles


  AliDielectronMC(const AliDielectronMC &c);
  AliDielectronMC &operator=(const AliDielectronMC &c);

  Bool_t IsMCMotherToEEesd(const AliMCParticle *particle, Int_t pdgMother);
  Bool_t IsMCMotherToEEaod(const AliAODMCParticle *particle, Int_t pdgMother);



  Int_t GetLabelMotherWithPdgESD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Int_t GetLabelMotherWithPdgAOD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);

  Bool_t ComparePDG(Int_t particlePDG, Int_t requiredPDG, Bool_t pdgExclusion, Bool_t checkBothCharges) const;
  Bool_t CheckIsRadiative(Int_t label) const;
  Bool_t CheckRadiativeDecision(Int_t mLabel, const AliDielectronSignalMC * const signalMC) const;
  Bool_t MotherIsGrandmother(int labelM1, int labelM2, int labelG1, int labelG2, bool motherIsGrandmother) const;
  Bool_t CheckStackParticle(Int_t labelPart, Int_t requiredPDG) const;
  Bool_t CompareDaughterPDG(Int_t labelM, Int_t reqPDG, Bool_t PDGexclusion, Bool_t CheckBothChargesDaughter) const;

  //MC - Heavy Flavour related methods
  Bool_t LoadHFPairs();
  Int_t  IsaBhadron(Int_t pdg) const;

  ClassDef(AliDielectronMC, 2)
};

//
// inline functions
//
inline Bool_t AliDielectronMC::IsMotherPdg(const AliDielectronPair* pair, Int_t pdgMother)
{
  return IsMotherPdg(pair->GetFirstDaughterP(),pair->GetSecondDaughterP(),pdgMother);
}
//___________________________________________________________
inline Bool_t AliDielectronMC::IsMotherPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother){
  return GetLabelMotherWithPdg(particle1,particle2,pdgMother)>=0;
}
//___________________________________________________________
inline Int_t AliDielectronMC::GetLabelMotherWithPdg(const AliDielectronPair* pair, Int_t pdgMother){
  return GetLabelMotherWithPdg(pair->GetFirstDaughterP(),pair->GetSecondDaughterP(),pdgMother);
}

#endif
