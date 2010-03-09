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

  static AliDielectronMC* Instance();
  
  void Initialize();                              // initialization
  Int_t GetNMCTracks();                                     // return number of generated tracks
  Int_t GetNMCTracksFromStack();                            // return number of generated tracks from stack
  Int_t GetMCPID(AliESDtrack* _track);                      // return MC PID
  Int_t GetMCPIDFromStack(AliESDtrack* _track);             // return MC PID
  Int_t GetMotherPDG(AliESDtrack* _track);                  // return mother PID from the MC stack
  Int_t GetMotherPDGFromStack(AliESDtrack* _track);         // return mother PID from the MC stack
  Int_t GetMCProcess(AliESDtrack* _track);                  // return process number
  Int_t GetMCProcessFromStack(AliESDtrack* _track);         // return process number
  Int_t GetMCProcessMother(AliESDtrack* _track);            // return process number of the mother track
  Int_t GetMCProcessMotherFromStack(AliESDtrack* _track);   // return process number of the mother track
  
  Bool_t ConnectMCEvent();
  Bool_t UpdateStack();

  Bool_t IsMotherPdg(const TObject* particle, Int_t pdgMother);
  Bool_t IsMotherPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Bool_t IsMCMotherToEE(const AliVParticle *particle, Int_t pdgMother);
  
  AliVParticle* GetMCTrackFromMCEvent(AliVParticle *track);   // return MC track directly from MC event
  AliVParticle* GetMCTrackFromMCEvent(Int_t _itrk);           // return MC track directly from MC event
  TParticle* GetMCTrackFromStack(AliESDtrack* _track);        // return MC track from stack
  AliMCParticle* GetMCTrack(AliESDtrack* _track);             // return MC track associated with reco track
  TParticle* GetMCTrackMotherFromStack(AliESDtrack* _track);  // return MC mother track from stack
  AliMCParticle* GetMCTrackMother(AliESDtrack* _track);       // return MC mother track from stack

  void GetDaughters(const TObject *mother, AliVParticle* &d1, AliVParticle* &d2);
  
private:
  AliMCEvent    *fMCEvent;  // MC event object
  AliStack      *fStack;    // MC stack

  AnalysisType fAnaType;    // Analysis type
  
  static AliDielectronMC* fgInstance; //! singleton pointer
  
  AliDielectronMC(const AliDielectronMC &c);
  AliDielectronMC &operator=(const AliDielectronMC &c);

  Bool_t IsMCMotherToEEesd(const AliMCParticle *particle, Int_t pdgMother);
  Bool_t IsMCMotherToEEaod(const AliAODMCParticle *particle, Int_t pdgMother);

  Bool_t IsMotherPdgESD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  Bool_t IsMotherPdgAOD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother);
  
  ClassDef(AliDielectronMC, 0)
};

inline Bool_t AliDielectronMC::IsMotherPdg(const TObject* particle, Int_t pdgMother)
{
  //
  //
  //
  if (particle->IsA()!=AliDielectronPair::Class()) return kFALSE;
  const AliDielectronPair *pair=static_cast<const AliDielectronPair*>(particle);
  return IsMotherPdg(pair->GetFirstDaughter(),pair->GetSecondDaughter(),pdgMother);
}

#endif 
