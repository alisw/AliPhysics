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
#endif
class AliESDEvent;
class AliHFEpid;
class AliStack;
class AliMCEvent;
class AliESDtrack;
class TParticle;
class AliMCParticle;
class AliAODMCParticle;

#include "AliDielectronPair.h"

class AliDielectronMC : public TObject{
  
public:
  enum AnalysisType {kUNSET=0, kESD, kAOD};
  
  AliDielectronMC(AnalysisType type=kUNSET);
  virtual ~AliDielectronMC();

  void SetHasMC(Bool_t hasMC) { fHasMC=hasMC; }
  Bool_t HasMC() const { return fHasMC; }
  
  static AliDielectronMC* Instance();
  
  void Initialize();                              // initialization
  Int_t GetNMCTracks();                                     // return number of generated tracks
  Int_t GetNMCTracksFromStack();                            // return number of generated tracks from stack
  Int_t GetMCPID(const AliESDtrack* _track);                      // return MC PID
  Int_t GetMCPIDFromStack(const AliESDtrack* _track);             // return MC PID
  Int_t GetMotherPDG(const AliESDtrack* _track);                  // return mother PID from the MC stack
  Int_t GetMotherPDGFromStack(const AliESDtrack* _track);         // return mother PID from the MC stack
  Int_t GetMCProcess(const AliESDtrack* _track);                  // return process number
  Int_t GetMCProcessFromStack(const AliESDtrack* _track);         // return process number
  Int_t GetMCProcessMother(const AliESDtrack* _track);            // return process number of the mother track
  Int_t GetMCProcessMotherFromStack(const AliESDtrack* _track);   // return process number of the mother track
  
  Bool_t ConnectMCEvent();
  Bool_t UpdateStack();

  Bool_t IsMotherPdg(const AliDielectronPair* pair, Int_t pdgMother);
  Bool_t IsMotherPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Bool_t IsMCMotherToEE(const AliVParticle *particle, Int_t pdgMother);

  Bool_t HaveSameMother(const AliDielectronPair *pair);
  
  Int_t GetLabelMotherWithPdg(const AliDielectronPair* pair, Int_t pdgMother);
  Int_t GetLabelMotherWithPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  
  AliVParticle* GetMCTrackFromMCEvent(AliVParticle *track);   // return MC track directly from MC event
  AliVParticle* GetMCTrackFromMCEvent(Int_t itrk);           // return MC track directly from MC event
  TParticle* GetMCTrackFromStack(const AliESDtrack* _track);        // return MC track from stack
  AliMCParticle* GetMCTrack(const AliESDtrack* _track);             // return MC track associated with reco track
  
  TParticle* GetMCTrackMotherFromStack(const AliESDtrack* _track);  // return MC mother track from stack
  AliMCParticle* GetMCTrackMother(const AliESDtrack* _track);       // return MC mother track from stack
  AliMCParticle* GetMCTrackMother(const AliMCParticle* _particle);       // return MC mother track from stack
  AliAODMCParticle* GetMCTrackMother(const AliAODMCParticle* _particle);       // return MC mother track from stack

  Int_t NumberOfDaughters(const AliESDtrack* track);                 // return number of daughters
  Int_t NumberOfDaughters(const AliMCParticle* particle);                 // return number of daughters
  Int_t NumberOfDaughters(const AliAODMCParticle* particle);                 // return number of daughters

  void GetDaughters(const TObject *mother, AliVParticle* &d1, AliVParticle* &d2);
  
private:
  AliMCEvent    *fMCEvent;  // MC event object
  AliStack      *fStack;    // MC stack

  AnalysisType fAnaType;    // Analysis type
  Bool_t fHasMC;            // Do we have an MC handler?
  
  static AliDielectronMC* fgInstance; //! singleton pointer
  
  AliDielectronMC(const AliDielectronMC &c);
  AliDielectronMC &operator=(const AliDielectronMC &c);

  Bool_t IsMCMotherToEEesd(const AliMCParticle *particle, Int_t pdgMother);
  Bool_t IsMCMotherToEEaod(const AliAODMCParticle *particle, Int_t pdgMother);

  Int_t GetLabelMotherWithPdgESD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Int_t GetLabelMotherWithPdgAOD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  
  ClassDef(AliDielectronMC, 0)
};

//
// inline functions
//
inline Bool_t AliDielectronMC::IsMotherPdg(const AliDielectronPair* pair, Int_t pdgMother)
{
  return IsMotherPdg(pair->GetFirstDaughter(),pair->GetSecondDaughter(),pdgMother);
}
//___________________________________________________________
inline Bool_t AliDielectronMC::IsMotherPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother){
  return GetLabelMotherWithPdg(particle1,particle2,pdgMother)>=0;
}
//___________________________________________________________
inline Int_t AliDielectronMC::GetLabelMotherWithPdg(const AliDielectronPair* pair, Int_t pdgMother){
  return GetLabelMotherWithPdg(pair->GetFirstDaughter(),pair->GetSecondDaughter(),pdgMother);
}

#endif 
