/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//#####################################################
//#                                                   #
//#              Class AliDielectronMC                #
//#       Cut Class for Jpsi->e+e- analysis           #
//#                                                   #
//#   by WooJin J. Park, GSI / W.J.Park@gsi.de        #
//#                                                   #
//#####################################################

#include <AliAnalysisManager.h>
#include <AliAODHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include <AliAODMCHeader.h>
#include <AliStack.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliLog.h>

#include <AliGenCocktailEventHeader.h>
#include <AliGenHijingEventHeader.h>

#include <TClonesArray.h>
#include <TParticle.h>
#include <TMCProcess.h>

#include "AliDielectronSignalMC.h"
#include "AliDielectronMC.h"


ClassImp(AliDielectronMC)

AliDielectronMC* AliDielectronMC::fgInstance=0x0;

//____________________________________________________________
AliDielectronMC* AliDielectronMC::Instance()
{
  //
  // return pointer to singleton implementation
  //
  if (fgInstance) return fgInstance;

  AnalysisType type=kUNSET;
  Bool_t hasMC=kFALSE;
  if (AliAnalysisManager::GetAnalysisManager()){
    if (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()) type=kESD;
    else if (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()) type=kAOD;

    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(type == kESD) hasMC=mcHandler!=0x0;
    }

  fgInstance=new AliDielectronMC(type);

  fgInstance->SetHasMC(hasMC);

  return fgInstance;
}

//____________________________________________________________
AliDielectronMC::AliDielectronMC(AnalysisType type):
  fMCEvent(0x0),
  fStack(0x0),
  fAnaType(type),
  fHasMC(kTRUE),
  fHasHijingHeader(-1),
  fMcArray(0x0)
{
  //
  // default constructor
  //
}


//____________________________________________________________
AliDielectronMC::~AliDielectronMC()
{
  //
  // default destructor
  //

}

//____________________________________________________________
void AliDielectronMC::Initialize()
{
  //
  // initialize MC class
  //
  if (!ConnectMCEvent()) AliError("Initialization of MC object failed!");
}

//____________________________________________________________
Int_t AliDielectronMC::GetNMCTracks()
{
  //
  //  return the number of generated tracks from MC event
  //
  if(fAnaType == kESD){
    if (!fMCEvent){ AliError("No fMCEvent"); return 0; }
    return fMCEvent->GetNumberOfTracks();}
  else if(fAnaType == kAOD){
    if(!fMcArray) { AliError("No fMcArray"); return 0; }
    return fMcArray->GetEntriesFast();
  }
  return 0;
}

//____________________________________________________________
Int_t AliDielectronMC::GetNMCTracksFromStack()
{
  //
  //  return the number of generated tracks from stack
  //
  if (!fStack){ AliError("No fStack"); return -999; }
  return fStack->GetNtrack();
}

//____________________________________________________________
Int_t AliDielectronMC::GetNPrimary() const
{
  //
  //  return the number of primary track from MC event
  //
  if (!fMCEvent){ AliError("No fMCEvent"); return 0; }
  return fMCEvent->GetNumberOfPrimaries();
}

//____________________________________________________________
Int_t AliDielectronMC::GetNPrimaryFromStack()
{
  //
  //  return the number of primary track from stack
  //
  if (!fStack){ AliError("No fStack"); return -999; }
  return fStack->GetNprimary();
}

//____________________________________________________________
AliVParticle* AliDielectronMC::GetMCTrackFromMCEvent(Int_t label) const
{
  //
  // return MC track directly from MC event
  // used not only for tracks but for mothers as well, therefore do not use abs(label)
  //
  if (label<0) return NULL;
  AliVParticle * track=0x0;
  if(fAnaType == kESD){
    if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
    track = fMCEvent->GetTrack(label); //  tracks from MC event (ESD)
  } else if(fAnaType == kAOD) {
    if (!fMcArray){ AliError("No fMcArray"); return NULL;}
    if (label>fMcArray->GetEntriesFast()) { AliWarning(Form("track %d out of array size %d",label,fMcArray->GetEntriesFast())); return NULL;}
    track = (AliVParticle*)fMCEvent->GetTrack(label);
  }
  return track;
}

//____________________________________________________________
Bool_t AliDielectronMC::ConnectMCEvent()
{
  //
  // connect stack object from the mc handler
  //
  fMcArray = 0x0;
  fMCEvent = 0x0;
  fHasHijingHeader=-1;

  if(fAnaType == kESD){
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcHandler){ /*AliError("Could not retrive MC event handler!");*/ return kFALSE; }
    if (!mcHandler->InitOk() ){ return kFALSE;}
    if (!mcHandler->TreeK() ) { return kFALSE;}
//    if (!mcHandler->TreeTR() ){ Printf("ERROR: !mcHandler->TreeTR()"); return kFALSE;}

    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if (!mcEvent){ AliError("Could not retrieve MC event!"); return kFALSE; }
    fMCEvent = mcEvent;

    if (!UpdateStack()) return kFALSE;
  }
  else if(fAnaType == kAOD)
  {
    AliAODInputHandler* aodHandler=(AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodHandler) return kFALSE;
    AliAODEvent *aod=aodHandler->GetEvent();
    if (!aod) return kFALSE;

    fMCEvent = aodHandler->MCEvent();
    if (!fMCEvent) AliError("No MCEvent available");

    fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMcArray){ /*AliError("Could not retrieve MC array!");*/ return kFALSE; }
    else fHasMC=kTRUE;
  }
  return kTRUE;
}

//____________________________________________________________
Bool_t AliDielectronMC::UpdateStack()
{
  //
  // update stack with new event
  //
  if (!fMCEvent){ AliError("No fMCEvent"); return kFALSE;}
  AliStack* stack = fMCEvent->Stack();
  if (!stack){ AliError("Could not retrive stack!"); return kFALSE; }
  fStack = stack;
  return kTRUE;
}

//____________________________________________________________
AliMCParticle* AliDielectronMC::GetMCTrack( const AliESDtrack* _track)
{
  //
  // return MC track
  //
  if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
  Int_t nStack = fMCEvent->GetNumberOfTracks();
  Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
  if(label>nStack)return NULL;
  AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
  return mctrack;
}


//____________________________________________________________
AliAODMCParticle* AliDielectronMC::GetMCTrack( const AliAODTrack* _track)
{
  //
  // return MC track
  //
 if(!fMcArray) { AliError("No fMcArray"); return NULL;}
 Int_t nStack = fMcArray->GetEntriesFast();
 Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
 if(label > nStack) return NULL;
 AliAODMCParticle *mctrack = (AliAODMCParticle*)fMcArray->At(label);
 return mctrack;
}

//____________________________________________________________
TParticle* AliDielectronMC::GetMCTrackFromStack(const AliESDtrack* _track)
{
  //
  // return MC track from stack
  //
  Int_t label = TMath::Abs(_track->GetLabel());
  if (!fStack) AliWarning("fStack is not available. Update stack first.");
  TParticle* mcpart = fStack->Particle(label);
  if (!mcpart) return NULL;
  return mcpart;
}

//____________________________________________________________
AliMCParticle* AliDielectronMC::GetMCTrackMother(const AliESDtrack* _track)
{
  //
  // return MC track mother
  //
  AliMCParticle* mcpart = GetMCTrack(_track);
  if (!mcpart) return NULL;
  if(mcpart->GetMother()<0) return NULL;
  AliMCParticle* mcmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(mcpart->GetMother()));
  if (!mcmother) return NULL;
  return mcmother;
}

//______________________________________________________________
AliAODMCParticle* AliDielectronMC::GetMCTrackMother(const AliAODTrack* _track)
{
 //
 // return MC track mother
 //
 AliAODMCParticle* mcpart = GetMCTrack(_track);
 if (!mcpart) return NULL;
 if(mcpart->GetMother() < 0) return NULL;
 AliAODMCParticle* mcmother = dynamic_cast<AliAODMCParticle *>(fMcArray->At(mcpart->GetMother()));
 if (!mcmother) return NULL;
 return mcmother;
}
//____________________________________________________________
AliMCParticle* AliDielectronMC::GetMCTrackMother(const AliMCParticle* _particle){
  //
  // return MC track mother
  //
  if(_particle->GetMother() < 0) return NULL;
  AliMCParticle* mcmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(_particle->GetMother()));
  return mcmother;
}

//____________________________________________________________
AliAODMCParticle* AliDielectronMC::GetMCTrackMother(const AliAODMCParticle* _particle){
  //
  // return MC track mother
  //
  if( _particle->GetMother() < 0) return NULL;
  AliAODMCParticle* mcmother = dynamic_cast<AliAODMCParticle *>(fMcArray->At(_particle->GetMother()));
  return mcmother;
}

//____________________________________________________________
TParticle* AliDielectronMC::GetMCTrackMotherFromStack(const AliESDtrack* _track)
{
  //
  // return MC track mother from stack
  //
  TParticle* mcpart = GetMCTrackFromStack(_track);
  if ( !mcpart || mcpart->GetFirstMother()<=0 ) return NULL;
  TParticle* mcmother = fStack->Particle(mcpart->GetFirstMother());
  if (!mcmother) return NULL;
  return mcmother;
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCPID(const AliESDtrack* _track)
{
  //
  // return PDG code of the track from the MC truth info
  //
  AliMCParticle* mcpart = GetMCTrack(_track);
  if (!mcpart) return -999;
  return mcpart->PdgCode();
}

//__________________________________________________________
Int_t AliDielectronMC::GetMCPID(const AliAODTrack* _track)
{
 //
 // return PDG code of the track from the MC truth info
 //
 AliAODMCParticle* mcpart = GetMCTrack(_track);
 if (!mcpart) return -999;
 return mcpart->PdgCode();
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCPIDFromStack(const AliESDtrack* _track)
{
  //
  // return MC PDG code from stack
  //
  TParticle* mcpart = GetMCTrackFromStack(_track);
  if (!mcpart) return -999;
  return mcpart->GetPdgCode();
}

//____________________________________________________________
Int_t AliDielectronMC::GetMotherPDG( const AliESDtrack* _track)
{
  //
  // return PDG code of the mother track from the MC truth info
  //
  AliMCParticle* mcmother = GetMCTrackMother(_track);
  if (!mcmother) return -999;
  return mcmother->PdgCode();
}

//________________________________________________________
Int_t AliDielectronMC::GetMotherPDG( const AliAODTrack* _track)
{
  //
  // return PDG code of the mother track from the MC truth info
  //
  AliAODMCParticle* mcmother = GetMCTrackMother(_track);
  if (!mcmother) return -999;
  return mcmother->PdgCode();
}

//________________________________________________________
Int_t AliDielectronMC::GetMotherPDG( const AliMCParticle* _track)
{
  //
  // return PDG code of the mother track from the MC truth info
  //
  AliMCParticle* mcmother = GetMCTrackMother(_track);
  if (!mcmother) return -999;
  return mcmother->PdgCode();
}

//________________________________________________________
Int_t AliDielectronMC::GetMotherPDG( const AliAODMCParticle* _track)
{
  //
  // return PDG code of the mother track from the MC truth info
  //
  AliAODMCParticle* mcmother = GetMCTrackMother(_track);
  if (!mcmother) return -999;
  return mcmother->PdgCode();
}

//____________________________________________________________
Int_t AliDielectronMC::GetMotherPDGFromStack(const AliESDtrack* _track)
{
  //
  // return PDG code of the mother track from stack
  //
  TParticle* mcmother = GetMCTrackMotherFromStack(_track);
  if (!mcmother) return -999;
  return mcmother->GetPdgCode();
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCProcess(const AliESDtrack* _track)
{
  //
  // return process number of the track
  //
  AliMCParticle* mcpart = GetMCTrack(_track);
  if (!mcpart) return -999;
  return 0;
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCProcessFromStack(const AliESDtrack* _track)
{
  //
  // return process number of the track
  //
  TParticle* mcpart = GetMCTrackFromStack(_track);
  if (!mcpart) return -999;
  return mcpart->GetUniqueID();
}

//____________________________________________________________
Int_t AliDielectronMC::NumberOfDaughters(const AliESDtrack* track)
{
  //
  // returns the number of daughters
  //
  AliMCParticle *mcmother=GetMCTrackMother(track);
  if(!mcmother||!mcmother->Particle()) return -999;
  //   return mcmother->GetFirstDaughter()>0?mcmother->GetLastDaughter()-mcmother->GetFirstDaughter()+1:0;
  return mcmother->Particle()->GetNDaughters();
}

//_________________________________________________________
Int_t AliDielectronMC::NumberOfDaughters(const AliAODTrack* track)
{
  //
  // returns the number of daughters
  //
  AliAODMCParticle *mcmother=GetMCTrackMother(track);
  if(!mcmother) return -999;
  return NumberOfDaughters(mcmother);

}

//____________________________________________________________
Int_t AliDielectronMC::NumberOfDaughters(const AliMCParticle* particle)
{
  //
  // returns the number of daughters
  //
  AliMCParticle *mcmother=GetMCTrackMother(particle);
  if(!mcmother||!mcmother->Particle()) return -999;
  //return mcmother->GetFirstDaughter()>0?mcmother->GetLastDaughter()-mcmother->GetFirstDaughter()+1:0;
  return mcmother->Particle()->GetNDaughters();
}

//____________________________________________________________
Int_t AliDielectronMC::NumberOfDaughters(const AliAODMCParticle* particle)
{
  //
  // returns the number of daughters
  //
  AliAODMCParticle *mcmother=GetMCTrackMother(particle);
  if(!mcmother) return -999;
  return mcmother->GetNDaughters();
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCProcessMother(const AliESDtrack* _track)
{
  //
  // return process number of the mother of the track
  //
  AliMCParticle* mcmother = GetMCTrackMother(_track);
  if (!mcmother) return -999;
  return 0;
}

//____________________________________________________________
Int_t AliDielectronMC::GetMCProcessMotherFromStack(const AliESDtrack* _track)
{
  //
  // return process number of the mother of the track
  //
  TParticle* mcmother = GetMCTrackMotherFromStack(_track);
  if (!mcmother) return -999;
  return mcmother->GetUniqueID();
}

//____________________________________________________________
Bool_t AliDielectronMC::IsMCMotherToEE(const AliVParticle *particle, Int_t pdgMother)
{
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  //
  if (fAnaType==kESD && !fMCEvent) return kFALSE;
  if (fAnaType==kAOD && !fMcArray) return kFALSE;
  if (!particle) return kFALSE;

  if (particle->IsA()==AliMCParticle::Class()){
    return IsMCMotherToEEesd(static_cast<const AliMCParticle*>(particle),pdgMother);
  } else if (particle->IsA()==AliAODMCParticle::Class()){
   return IsMCMotherToEEaod(static_cast<const AliAODMCParticle*>(particle),pdgMother);
  } else {
    AliError("Unknown particle type");
  }
  return kFALSE;
}

//____________________________________________________________
Bool_t AliDielectronMC::IsMCMotherToEEesd(const AliMCParticle *particle, Int_t pdgMother)
{
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  // ESD case
  //

  //check pdg code
  if (particle->PdgCode()!=pdgMother) return kFALSE;
  Int_t ifirst = particle->GetFirstDaughter();
  Int_t ilast  = particle->GetLastDaughter();

  //check number of daughters
  if ((ilast-ifirst)!=1) return kFALSE;
  AliMCParticle *firstD=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(ifirst));
  AliMCParticle *secondD=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(ilast));

  //TODO: check how you can get rid of the hardcoded numbers. One should make use of the PdgCodes set in AliDielectron!!!
  if (firstD->Charge()>0){
    if (firstD->PdgCode()!=-11) return kFALSE;
    if (secondD->PdgCode()!=11) return kFALSE;
  }else{
    if (firstD->PdgCode()!=11) return kFALSE;
    if (secondD->PdgCode()!=-11) return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________
Bool_t AliDielectronMC::IsMCMotherToEEaod(const AliAODMCParticle *particle, Int_t pdgMother)
{
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  // AOD case
  //

  if (particle->GetPdgCode()!=pdgMother) return kFALSE;
  if (particle->GetNDaughters()!=2) return kFALSE;

  Int_t ifirst = particle->GetDaughter(0);
  Int_t ilast  = particle->GetDaughter(1);

  //check number of daughters
  if ((ilast-ifirst)!=1) return kFALSE;

  AliAODMCParticle *firstD=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(ifirst));
  AliAODMCParticle *secondD=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(ilast));

  //TODO: check how you can get rid of the hardcoded numbers. One should make use of the PdgCodes set in AliDielectron!!!

  if (firstD->Charge()>0){
    if (firstD->GetPdgCode()!=-11) return kFALSE;
    if (secondD->GetPdgCode()!=11) return kFALSE;
  }else{
    if (firstD->GetPdgCode()!=11) return kFALSE;
    if (secondD->GetPdgCode()!=-11) return kFALSE;
  }
  return kTRUE;
}

//____________________________________________________________
Int_t AliDielectronMC::GetLabelMotherWithPdg(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother)
{
  //
  // test if mother of particle 1 and 2 has pdgCode pdgMother and is the same;
  //
  if (fAnaType==kESD){
  if (!fMCEvent) return -1;
  return GetLabelMotherWithPdgESD(particle1, particle2, pdgMother);
  }
  else if (fAnaType==kAOD)
  {
  if (!fMcArray) return -1;
  return GetLabelMotherWithPdgAOD(particle1, particle2, pdgMother);
  }

  return -1;
}

//____________________________________________________________
Int_t AliDielectronMC::GetLabelMotherWithPdgESD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother)
{
  //
  // test if mother of particle 1 and 2 has pdgCode +-11 (electron),
  //    have the same mother and the mother had pdg code pdgMother
  // ESD case
  //TODO: check how you can get rid of the hardcoded numbers. One should make use of the PdgCodes set in AliDielectron!!!
  //
  // negative label indicate poor matching quality
  Int_t lblPart1 = TMath::Abs(particle1->GetLabel());
  Int_t lblPart2 = TMath::Abs(particle2->GetLabel());
  AliMCParticle *mcPart1=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(lblPart1));
  AliMCParticle *mcPart2=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(lblPart2));

  if (!mcPart1||!mcPart2) return -1;

  Int_t lblMother1=mcPart1->GetMother();
  Int_t lblMother2=mcPart2->GetMother();
  AliMCParticle *mcMother1=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(lblMother1));

  if (!mcMother1) return -1;
  if (lblMother1!=lblMother2) return -1;
  if (TMath::Abs(mcPart1->PdgCode())!=11) return -1;
  if (mcPart1->PdgCode()!=-mcPart2->PdgCode()) return -1;
  if (mcMother1->PdgCode()!=pdgMother) return -1;

  return lblMother1;
}

//____________________________________________________________
Int_t AliDielectronMC::GetLabelMotherWithPdgAOD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother)
{
  //
  // test if mother of particle 1 and 2 has pdgCode +-11 (electron),
  //    have the same mother and the mother had pdg code pdgMother
  // AOD case
  //TODO: check how you can get rid of the hardcoded numbers. One should make use of the PdgCodes set in AliDielectron!!!
  //
  // negative label indicate poor matching quality
  Int_t lblPart1 = TMath::Abs(particle1->GetLabel());
  Int_t lblPart2 = TMath::Abs(particle2->GetLabel());
  AliAODMCParticle *mcPart1=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(lblPart1));
  AliAODMCParticle *mcPart2=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(lblPart2));

  if (!mcPart1||!mcPart2) return -1;

  Int_t lblMother1=mcPart1->GetMother();
  Int_t lblMother2=mcPart2->GetMother();
  AliAODMCParticle *mcMother1=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(lblMother1));

  if (!mcMother1) return -1;
  if (lblMother1!=lblMother2) return -1;
  if (TMath::Abs(mcPart1->GetPdgCode())!=11) return -1;
  if (mcPart1->GetPdgCode()!=-mcPart2->GetPdgCode()) return -1;
  if (mcMother1->GetPdgCode()!=pdgMother) return -1;

  return lblMother1;
}

//____________________________________________________________
void AliDielectronMC::GetDaughters(const TObject *mother, AliVParticle* &d1, AliVParticle* &d2)
{
  //
  // Get First two daughters of the mother
  //
  Int_t lblD1=-1;
  Int_t lblD2=-1;
  d1=0;
  d2=0;
  if (fAnaType==kAOD){
    if(!fMcArray) return;
    const AliAODMCParticle *aodMother=static_cast<const AliAODMCParticle*>(mother);
    lblD1=aodMother->GetDaughter(0);
    lblD2=aodMother->GetDaughter(1);
    d1 = (AliVParticle*)fMcArray->At(lblD1);
    d2 = (AliVParticle*)fMcArray->At(lblD2);
   } else if (fAnaType==kESD){
    if (!fMCEvent) return;
    const AliMCParticle *aodMother=static_cast<const AliMCParticle*>(mother);
    lblD1=aodMother->GetFirstDaughter();
    lblD2=aodMother->GetLastDaughter();
    d1=fMCEvent->GetTrack(lblD1);
    d2=fMCEvent->GetTrack(lblD2);
   }
}


//________________________________________________________________________________
Int_t AliDielectronMC::GetMothersLabel(Int_t daughterLabel) const {
  //
  //  Get the label of the mother for particle with label daughterLabel
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(daughterLabel<0) return -1;
  if (fAnaType==kAOD) {
    if(!fMcArray) return -1;
    if(GetMCTrackFromMCEvent(daughterLabel))
      return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(daughterLabel)))->GetMother();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return -1;
    if(GetMCTrackFromMCEvent(daughterLabel))
      return (static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(daughterLabel)))->GetMother();
  }
  return -1;
}


//________________________________________________________________________________
Int_t AliDielectronMC::GetPdgFromLabel(Int_t label) const {
  //
  //  Get particle code using the label from stack
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(label<0) return 0;
  if(fAnaType==kAOD) {
    if(!fMcArray) return 0;
    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->PdgCode();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return 0;
    return (static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(label)))->PdgCode();
  }
  return 0;
}


//________________________________________________________________________________
Bool_t AliDielectronMC::ComparePDG(Int_t particlePDG, Int_t requiredPDG, Bool_t pdgExclusion, Bool_t checkBothCharges) const {
  //
  //  Test the PDG codes of particles with the required ones
  //
  Bool_t result = kTRUE;
  Int_t absRequiredPDG = TMath::Abs(requiredPDG);

  switch(absRequiredPDG) {
  case 0:
    result = kTRUE;    // PDG not required (any code will do fine)
    break;
  case 100:     // light flavoured mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=100 && TMath::Abs(particlePDG)<=199;
    else {
      if(requiredPDG>0) result = particlePDG>=100 && particlePDG<=199;
      if(requiredPDG<0) result = particlePDG>=-199 && particlePDG<=-100;
    }
    break;
  case 1000:     // light flavoured baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=1000 && TMath::Abs(particlePDG)<=1999;
    else {
      if(requiredPDG>0) result = particlePDG>=1000 && particlePDG<=1999;
      if(requiredPDG<0) result = particlePDG>=-1999 && particlePDG<=-1000;
    }
    break;
  case 200:     // light flavoured mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=200 && TMath::Abs(particlePDG)<=299;
    else {
      if(requiredPDG>0)result = particlePDG>=200 && particlePDG<=299;
      if(requiredPDG<0)result = particlePDG>=-299 && particlePDG<=-200;
    }
    break;
  case 2000:     // light flavoured baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=2000 && TMath::Abs(particlePDG)<=2999;
    else {
      if(requiredPDG>0) result = particlePDG>=2000 && particlePDG<=2999;
      if(requiredPDG<0) result = particlePDG>=-2999 && particlePDG<=-2000;
    }
    break;
  case 300:     // all strange mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=300 && TMath::Abs(particlePDG)<=399;
    else {
      if(requiredPDG>0) result = particlePDG>=300 && particlePDG<=399;
      if(requiredPDG<0) result = particlePDG>=-399 && particlePDG<=-300;
    }
    break;
  case 3000:     // all strange baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=3000 && TMath::Abs(particlePDG)<=3999;
    else {
      if(requiredPDG>0) result = particlePDG>=3000 && particlePDG<=3999;
      if(requiredPDG<0) result = particlePDG>=-3999 && particlePDG<=-3000;
    }
    break;
  case 400:     // all charmed mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=400 && TMath::Abs(particlePDG)<=499;
    else {
      if(requiredPDG>0) result = particlePDG>=400 && particlePDG<=499;
      if(requiredPDG<0) result = particlePDG>=-499 && particlePDG<=-400;
    }
    break;
  case 401:     // open charm mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=400 && TMath::Abs(particlePDG)<=439;
    else {
      if(requiredPDG>0) result = particlePDG>=400 && particlePDG<=439;
      if(requiredPDG<0) result = particlePDG>=-439 && particlePDG<=-400;
    }
    break;
  case 402:     // open charm mesons and baryons together
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=400 && TMath::Abs(particlePDG)<=439) ||
	       (TMath::Abs(particlePDG)>=4000 && TMath::Abs(particlePDG)<=4399);
    else {
      if(requiredPDG>0) result = (particlePDG>=400 && particlePDG<=439) ||
			         (particlePDG>=4000 && particlePDG<=4399);
      if(requiredPDG<0) result = (particlePDG>=-439 && particlePDG<=-400) ||
 			         (particlePDG>=-4399 && particlePDG<=-4000);
    }
    break;
  case 403:     // all charm hadrons
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=400 && TMath::Abs(particlePDG)<=499) ||
	       (TMath::Abs(particlePDG)>=4000 && TMath::Abs(particlePDG)<=4999);
    else {
      if(requiredPDG>0) result = (particlePDG>=400 && particlePDG<=499) ||
			         (particlePDG>=4000 && particlePDG<=4999);
      if(requiredPDG<0) result = (particlePDG>=-499 && particlePDG<=-400) ||
 			         (particlePDG>=-4999 && particlePDG<=-4000);
    }
    break;
  case 404:     // charged open charmed mesons NO s quark
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=410 && TMath::Abs(particlePDG)<=419);
    else {
      if(requiredPDG>0) result = (particlePDG>=410 && particlePDG<=419);
      if(requiredPDG<0) result = (particlePDG>=-419 && particlePDG<=-410);
    }
    break;
  case 405:     // neutral open charmed mesons
    if(checkBothCharges)
      result =TMath::Abs(particlePDG)>=420 && TMath::Abs(particlePDG)<=429;
    else {
      if(requiredPDG>0) result = particlePDG>=420 && particlePDG<=429 ;
      if(requiredPDG<0) result = particlePDG>=-429 && particlePDG<=-420;
    }
    break;

    case 406:     // charged open charmed mesons with s quark
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=430 && TMath::Abs(particlePDG)<=439);
    else {
      if(requiredPDG>0) result = (particlePDG>=430 && particlePDG<=439);
      if(requiredPDG<0) result = (particlePDG>=-439 && particlePDG<=-430);
    }
    break;

  case 4000:     // all charmed baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=4000 && TMath::Abs(particlePDG)<=4999;
    else {
      if(requiredPDG>0) result = particlePDG>=4000 && particlePDG<=4999;
      if(requiredPDG<0) result = particlePDG>=-4999 && particlePDG<=-4000;
    }
    break;
  case 4001:     // open charm baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=4000 && TMath::Abs(particlePDG)<=4399;
    else {
      if(requiredPDG>0) result = particlePDG>=4000 && particlePDG<=4399;
      if(requiredPDG<0) result = particlePDG>=-4399 && particlePDG<=-4000;
    }
    break;
  case 500:      // all beauty mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=500 && TMath::Abs(particlePDG)<=599;
    else {
      if(requiredPDG>0) result = particlePDG>=500 && particlePDG<=599;
      if(requiredPDG<0) result = particlePDG>=-599 && particlePDG<=-500;
    }
    break;
  case 501:      // open beauty mesons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=500 && TMath::Abs(particlePDG)<=549;
    else {
      if(requiredPDG>0) result = particlePDG>=500 && particlePDG<=549;
      if(requiredPDG<0) result = particlePDG>=-549 && particlePDG<=-500;
    }
    break;
  case 502:      // open beauty mesons and baryons
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=500 && TMath::Abs(particlePDG)<=549) ||
	       (TMath::Abs(particlePDG)>=5000 && TMath::Abs(particlePDG)<=5499);
    else {
      if(requiredPDG>0) result = (particlePDG>=500 && particlePDG<=549) ||
			         (particlePDG>=5000 && particlePDG<=5499);
      if(requiredPDG<0) result = (particlePDG>=-549 && particlePDG<=-500) ||
                                 (particlePDG>=-5499 && particlePDG<=-5000);
    }
    break;
  case 503:      // all beauty hadrons
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=500 && TMath::Abs(particlePDG)<=599) ||
	       (TMath::Abs(particlePDG)>=5000 && TMath::Abs(particlePDG)<=5999);
    else {
      if(requiredPDG>0) result = (particlePDG>=500 && particlePDG<=599) ||
			         (particlePDG>=5000 && particlePDG<=5999);
      if(requiredPDG<0) result = (particlePDG>=-599 && particlePDG<=-500) ||
                                 (particlePDG>=-5999 && particlePDG<=-5000);
    }
    break;

  case 504:     // neutral open beauty mesons NO s quark
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=510 && TMath::Abs(particlePDG)<=519);
    else {
      if(requiredPDG>0) result = (particlePDG>=510 && particlePDG<=519);
      if(requiredPDG<0) result = (particlePDG>=-519 && particlePDG<=-510);
    }
    break;
  case 505:     // charged open beauty mesons
    if(checkBothCharges)
      result =TMath::Abs(particlePDG)>=520 && TMath::Abs(particlePDG)<=529;
    else {
      if(requiredPDG>0) result = particlePDG>=520 && particlePDG<=529 ;
      if(requiredPDG<0) result = particlePDG>=-529 && particlePDG<=-520;
    }
    break;

  case 506:     // charged open beauty mesons with s quark
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=530 && TMath::Abs(particlePDG)<=539);
    else {
      if(requiredPDG>0) result = (particlePDG>=530 && particlePDG<=539);
      if(requiredPDG<0) result = (particlePDG>=-539 && particlePDG<=-530);
    }
    break;

  case 5000:      // all beauty baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=5000 && TMath::Abs(particlePDG)<=5999;
    else {
      if(requiredPDG>0) result = particlePDG>=5000 && particlePDG<=5999;
      if(requiredPDG<0) result = particlePDG>=-5999 && particlePDG<=-5000;
    }
    break;
  case 5001:      // open beauty baryons
    if(checkBothCharges)
      result = TMath::Abs(particlePDG)>=5000 && TMath::Abs(particlePDG)<=5499;
    else {
      if(requiredPDG>0) result = particlePDG>=5000 && particlePDG<=5499;
      if(requiredPDG<0) result = particlePDG>=-5499 && particlePDG<=-5000;
    }
    break;
  case 902:      // // open charm,beauty  mesons and baryons together
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=400 && TMath::Abs(particlePDG)<=439) ||
      (TMath::Abs(particlePDG)>=4000 && TMath::Abs(particlePDG)<=4399) ||
      (TMath::Abs(particlePDG)>=500 && TMath::Abs(particlePDG)<=549) ||
      (TMath::Abs(particlePDG)>=5000 && TMath::Abs(particlePDG)<=5499);
    else {
      if(requiredPDG>0) result = (particlePDG>=400 && particlePDG<=439) ||
			  (particlePDG>=4000 && particlePDG<=4399)      ||
			  (particlePDG>=500 && particlePDG<=549)        ||
			  (particlePDG>=5000 && particlePDG<=5499);
      if(requiredPDG<0) result = (particlePDG>=-439 && particlePDG<=-400) ||
			  (particlePDG>=-4399 && particlePDG<=-4000)      ||
			  (particlePDG>=-549 && particlePDG<=-500)        ||
			  (particlePDG>=-5499 && particlePDG<=-5000);
    }
    break;
  case 903:      // // all hadrons in the code range 100-599, 1000-5999
    if(checkBothCharges)
      result = (TMath::Abs(particlePDG)>=100 && TMath::Abs(particlePDG)<=599) ||
      (TMath::Abs(particlePDG)>=1000 && TMath::Abs(particlePDG)<=5999);
    else {
      if(requiredPDG>0) result = (particlePDG>=100 && particlePDG<=599) ||
        (particlePDG>=1000 && particlePDG<=5999);
      if(requiredPDG<0) result = (particlePDG>=-599 && particlePDG<=-100) ||
        (particlePDG>=-5999 && particlePDG<=-1000);
    }
    break;
  default:          // all specific cases
    if(checkBothCharges)
      result = (absRequiredPDG==TMath::Abs(particlePDG));
    else
      result = (requiredPDG==particlePDG);
  }

  if(absRequiredPDG!=0 && pdgExclusion) result = !result;
  return result;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsPhysicalPrimary(Int_t label) const {
  //
  // Check if the particle with label "label" is a physical primary according to the
  // definition in AliStack::IsPhysicalPrimary(Int_t label)
  // Convention for being physical primary:
  // 1.) particles produced in the collision
  // 2.) stable particles with respect to strong and electromagnetic interactions
  // 3.) excludes initial state particles
  // 4.) includes products of directly produced Sigma0 hyperon decay
  // 5.) includes products of directly produced pi0 decays
  // 6.) includes products of directly produced beauty hadron decays
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsPhysicalPrimary();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    return fStack->IsPhysicalPrimary(label);
  }
  return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsPrimary(Int_t label) const {
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;

    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsPrimary();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    return (label>=0 && label<=GetNPrimary());
  }
  return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondary(Int_t label) const {
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    AliAODMCParticle* mctrack = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label));
    Bool_t isSecondary = mctrack->IsSecondaryFromMaterial() || mctrack->IsSecondaryFromWeakDecay();
    return isSecondary;
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    return (label>=GetNPrimary() && !IsPhysicalPrimary(label));
  }
  return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondaryFromWeakDecay(Int_t label) const {
  //
  // Check if the particle with label "label" is a physical secondary from weak decay according to the
  // definition in AliStack::IsSecondaryFromWeakDecay(Int_t label)
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsSecondaryFromWeakDecay();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    return fStack->IsSecondaryFromWeakDecay(label);
  }
  return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondaryFromMaterial(Int_t label) const {
  //
  // Check if the particle with label "label" is a physical secondary from weak decay according to the
  // definition in AliStack::IsSecondaryFromMaterial(Int_t label)
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsSecondaryFromMaterial();
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    return fStack->IsSecondaryFromMaterial(label);
  }
  return kFALSE;
}


//________________________________________________________________________________
Bool_t AliDielectronMC::CheckGEANTProcess(Int_t label, TMCProcess process) const {
  //
  //  Check the GEANT process for the particle
  //  NOTE: for tracks the absolute label should be passed
  //
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    UInt_t processID = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label))->GetMCProcessCode();
    //    printf("process: id %d --> %s \n",processID,TMCProcessName[processID]);
    return (process==processID);
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    AliError(Form("return of GEANT process not implemented for ESD "));
    return kFALSE;
  }
  return kFALSE;

}

//________________________________________________________________________________
Bool_t AliDielectronMC::CheckParticleSource(Int_t label, AliDielectronSignalMC::ESource source) const {
  //
  //  Check the source for the particle
  //  NOTE: for tracks the absolute label should be passed
  //

  switch (source) {
    case AliDielectronSignalMC::kDontCare :
      return kTRUE;
    break;
    case AliDielectronSignalMC::kPrimary :
      // true if label is in the list of particles from physics generator
      // NOTE: This includes all physics event history (initial state particles,
      //       exchange bosons, quarks, di-quarks, strings, un-stable particles, final state particles)
      //       Only the final state particles make it to the detector!!
      return IsPrimary(label);
    break;
    case AliDielectronSignalMC::kFinalState :
      // primary particles created in the collision which reach the detectors
      // These would be:
      // 1.) particles produced in the collision
      // 2.) stable particles with respect to strong and electromagnetic interactions
      // 3.) excludes initial state particles
      // 4.) includes products of directly produced Sigma0 hyperon decay
      // 5.) includes products of directly produced pi0 decays
      // 6.) includes products of directly produced beauty hadron decays
      return IsPhysicalPrimary(label);
    break;
    case AliDielectronSignalMC::kDirect :
      // Primary particles which do not have any mother
      // This is the case for:
      // 1.) Initial state particles (the 2 protons in Pythia pp collisions)
      // 2.) In some codes, with sudden freeze-out, all particles generated from the fireball are direct.
      //     There is no history for these particles.
      // 3.) Certain particles added via MC generator cocktails (e.g. J/psi added to pythia MB events)
      return (label>=0 && GetMothersLabel(label)<0);
      break;
    case AliDielectronSignalMC::kNoCocktail :
      // Particles from the HIJING event and NOT from the AliGenCocktail
      return (label>=0 && GetMothersLabel(label)>=0);
      break;
    case AliDielectronSignalMC::kSecondary :
      // particles which are created by the interaction of final state primaries with the detector
      // or particles from strange weakly decaying particles (e.g. lambda, kaons, etc.)
      return IsSecondary(label);
      // return (label>=GetNPrimary() && !IsPhysicalPrimary(label)); // old definition
    break;
    case AliDielectronSignalMC::kSecondaryFromWeakDecay :
      // secondary particle from weak decay
      // or particles from strange weakly decaying particles (e.g. lambda, kaons, etc.)
      return (IsSecondaryFromWeakDecay(label));
    break;
    case AliDielectronSignalMC::kSecondaryFromMaterial :
      // secondary particle from material
      return (IsSecondaryFromMaterial(label));
    break;
    case AliDielectronSignalMC::kFromBGEvent :
      // NOT implemented for AODs
      // used to select electrons which are not from injected signals.
      return (IsFromBGEvent(label));
      break;
    case AliDielectronSignalMC::kFinalStateFromBGEvent :
      // NOT implemented for AODs
      // used to select electrons which are not from injected signals.
      return (IsPhysicalPrimary(label) && IsFromBGEvent(label));
      break;
    default :
      return kFALSE;
  }
  return kFALSE;

}

/*
// (please keep this for reference...)
//________________________________________________________________________________
Bool_t AliDielectronMC::IsEleFromInjectedSignal(Int_t label) const {
  ///
  /// Function to check if the particle with label "label" originates from an injected signal.
  /// used criteria:
  /// isPhysPrim:   in a chain of decaying particles, denotes the first which is stable under strong and EM force.
  /// - injected:   the physical primary will NOT be fromBGEvent! (<= that is the relevant criterion!)
  /// - Geant:      the physical primary will be fromBGEvent!
  /// fromBGEvent:  particle comes from the generated MC event (?). true for most pi, K, p, some e.
  /// !fromBGEvent: particle could be injected or produced in Geant (from weak decay / photon conv).
  ///               weak decays and photon conv are handled by Geant -> will NOT be fromBGEvent!
  /// examples:     for injected J/psi -> ee: the e+- will be phys prim (because of EM decay), but not fromBGEvent!
  ///               for photon conversions: the photon will be phys prim. (not the e+-!), and will be fromBGEvent!
  ///               the same applies for weak decays.
  ///
  /// The following procedure is the most correct, and works fine when looping over the MC-ESD event,
  /// but for stack electrons coming from scattering etc (which are not physical primary) this leads to
  /// very long loops and may crash the task.
  /// One would have to apply at least some kinematic cuts before, to reject all/most of these cases.
  ///
  /// The short procedure, with the difference that electrons from conversions and weak decays are also rejected, is just:
  /// return !(IsFromBGEvent(label));
  ///
  if(label<0) return kFALSE;
  //if(label<0) label *= -1;  // not sure what is more correct

  Bool_t isPhysPrim  = IsPhysicalPrimary(label);
  Bool_t fromBGEvent = IsFromBGEvent(label);
  // if the particle isn't already the physical primary, then iteratively go back to it:
  Int_t labelMother = GetMCTrackFromMCEvent(label)->GetMother(); // works for AOD and ESD, no explicit cast needed.
  Int_t counter=0;
  while (!isPhysPrim) {
    Printf(Form("IsInjectedSignal(): label=%d pdg=%d labelMother=%d pdgmother=%d", label, GetPdgFromLabel(label), labelMother, GetPdgFromLabel(labelMother)));
    if (counter>10) {
      AliWarning(Form("probably infinite loop! label=%d pdg=%d labelMother=%d pdgmother=%d", label, GetPdgFromLabel(label), labelMother, GetPdgFromLabel(labelMother)));
      return kFALSE; // electrons coming from a long cascade of scattering particles are probably not injected...
    } counter++;
    if (labelMother<0) return kTRUE; // to avoid infinite loop
    isPhysPrim  = IsPhysicalPrimary(labelMother); // save property of mother
    fromBGEvent = IsFromBGEvent(labelMother);     // save property of mother
     // prepare to get properties of grandmother
    labelMother = GetMCTrackFromMCEvent(labelMother)->GetMother(); // works for AOD and ESD, no explicit cast needed.
  }
  // now check if the physical primary is from background event. if not, it was injected.
  if (!fromBGEvent) return kTRUE;

  return kFALSE;
}
*/


//________________________________________________________________________________
Bool_t AliDielectronMC::IsFromBGEvent(Int_t label) const {
  ///
  /// Check if the particle with label "label" is from the background MC event,
  /// which means that it is not injected.
  ///
  if(label<0) return kFALSE;
  if(fAnaType==kAOD) {
    AliWarning("IsFromBGEvent() not implemented for AOD!");
    return kFALSE;
//    if(!fMcArray) return kFALSE;
//    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsFromBGEvent(); // IsFromBGEvent() does not exist for AliAODMCParticle.
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    if (CheckHijingHeader()) return fMCEvent->IsFromBGEvent(label); // Works for HIJING inside Cocktail
    //else if (CheckSomeOtherHeader()) return ...;
    else {
      AliWarning("No headers to make decision! Assuming no injected signals are present.");
      return kTRUE;
    }
  }
  return kFALSE;
}


//________________________________________________________________________________
Bool_t AliDielectronMC::CheckHijingHeader() const {

//  if (fHasHijingHeader > -1) return Bool_t(fHasHijingHeader); // avoid many calls of the code below.

  if(fAnaType==kAOD) {
    AliWarning("CheckHijingHeader() not implemented for AOD!");
    return (fHasHijingHeader=0); //return kFALSE;
    //    AliAODInputHandler* aodHandler=(AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    //    if (!aodHandler) return kFALSE;
    //    AliAODEvent *aod=aodHandler->GetEvent();
    //    if (!aod) return kFALSE;
    //    AliAODHeader* header = aod->GetHeader(); // not sure if correct
  }
  else if(fAnaType==kESD) {
    // taken from AliMCEvent::IsFromBGEvent()
    if (!fMCEvent) return (fHasHijingHeader=0); //return kFALSE;
    AliGenCocktailEventHeader* coHeader = dynamic_cast<AliGenCocktailEventHeader*> (fMCEvent->GenEventHeader());
    if (!coHeader) return (fHasHijingHeader=0); //return kFALSE;
    TList* list = coHeader->GetHeaders();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(list->FindObject("Hijing"));
    if (hijingH) return (fHasHijingHeader=1); //return kTRUE;
  }
  return (fHasHijingHeader=0);
}


//________________________________________________________________________________
Bool_t AliDielectronMC::CheckIsRadiative(Int_t label) const
{
  //
  // Check if the particle has a three body decay, one being a photon
  //
  if(label<0) return kFALSE;


  if(fAnaType==kAOD) {
    if(!fMcArray) return kFALSE;
    AliAODMCParticle *mother=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label));
    if (!mother) return kFALSE;
    const Int_t nd=mother->GetNDaughters();
    if (nd==2) return kFALSE;
    for (Int_t i=2; i<nd; ++i)
      if (GetMCTrackFromMCEvent(mother->GetDaughter(0)+i)->PdgCode()!=22) return kFALSE; //last daughter is photon
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    AliMCParticle *mother=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(label));
    if (!mother) return kFALSE;
    const Int_t nd=(mother->GetLastDaughter()-mother->GetFirstDaughter()+1);
    if (nd==2) return kFALSE;
    for (Int_t i=2; i<nd; ++i)
      if (GetMCTrackFromMCEvent(mother->GetFirstDaughter()+i)->PdgCode()!=22) return kFALSE; //last daughters are photons
  }
  return kTRUE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::CheckRadiativeDecision(Int_t mLabel, const AliDielectronSignalMC * const signalMC) const
{
  //
  // Check for the decision of the radiative type request
  //

  if (!signalMC) return kFALSE;

  if (signalMC->GetJpsiRadiative()==AliDielectronSignalMC::kAll) return kTRUE;

  Bool_t isRadiative=CheckIsRadiative(mLabel);
  if ((signalMC->GetJpsiRadiative()==AliDielectronSignalMC::kIsRadiative) && !isRadiative) return kFALSE;
  if ((signalMC->GetJpsiRadiative()==AliDielectronSignalMC::kIsNotRadiative) && isRadiative) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsMCTruth(Int_t label, AliDielectronSignalMC* signalMC, Int_t branch) const {
  //
  // Check if the particle corresponds to the MC truth in signalMC in the branch specified
  //

  // NOTE:  Some particles have the sign of the label flipped. It is related to the quality of matching
  //        between the ESD and the MC track. The negative labels indicate a poor matching quality
  //if(label<0) return kFALSE;

  if(label<0) label *= -1;

  AliVParticle* part = GetMCTrackFromMCEvent(label);
  if (!part) {
    AliError(Form("Could not find MC particle with label %d",label));
    return kFALSE;
  }

  // check geant process if set
  if(signalMC->GetCheckGEANTProcess() && !CheckGEANTProcess(label,signalMC->GetGEANTProcess())) return kFALSE;

  // check the leg
  if(!ComparePDG(part->PdgCode(),signalMC->GetLegPDG(branch),signalMC->GetLegPDGexclude(branch),signalMC->GetCheckBothChargesLegs(branch))) return kFALSE;
  if(!CheckParticleSource(label, signalMC->GetLegSource(branch))) return kFALSE;

  // check the mother
  AliVParticle* mcMother=0x0;
  Int_t mLabel = -1;
  if(signalMC->GetMotherPDG(branch)!=0 || signalMC->GetMotherSource(branch)!=AliDielectronSignalMC::kDontCare) {
    if(part) {
      mLabel = GetMothersLabel(label);
      mcMother = GetMCTrackFromMCEvent(mLabel);
    }
    if(!mcMother && !signalMC->GetMotherPDGexclude(branch)) return kFALSE;

    if(!ComparePDG((mcMother ? mcMother->PdgCode() : 0),signalMC->GetMotherPDG(branch),signalMC->GetMotherPDGexclude(branch),signalMC->GetCheckBothChargesMothers(branch))) return kFALSE;
    if(!CheckParticleSource(mLabel, signalMC->GetMotherSource(branch))) return kFALSE;

    //check for radiative deday
    if (!CheckRadiativeDecision(mLabel, signalMC)) return kFALSE;
  }

  // check the grandmother
  AliVParticle* mcGrandMother=0x0;
  if(signalMC->GetGrandMotherPDG(branch)!=0 || signalMC->GetGrandMotherSource(branch)!=AliDielectronSignalMC::kDontCare) {
    Int_t gmLabel = -1;
    if(mcMother) {
      gmLabel = GetMothersLabel(mLabel);
      mcGrandMother = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(gmLabel));
    }
    if(!mcGrandMother && !signalMC->GetGrandMotherPDGexclude(branch)) return kFALSE;

    if(!ComparePDG((mcGrandMother ? mcGrandMother->PdgCode() : 0),signalMC->GetGrandMotherPDG(branch),signalMC->GetGrandMotherPDGexclude(branch),signalMC->GetCheckBothChargesGrandMothers(branch))) return kFALSE;
    if(!CheckParticleSource(gmLabel, signalMC->GetGrandMotherSource(branch))) return kFALSE;
  }

  return kTRUE;
}



//________________________________________________________________________________
Bool_t AliDielectronMC::IsMCTruth(const AliDielectronPair* pair, const AliDielectronSignalMC* signalMC) const {
  //
  // Check if the pair corresponds to the MC truth in signalMC
  //

  // legs (daughters)
  const AliVParticle * mcD1 = pair->GetFirstDaughterP();
  const AliVParticle * mcD2 = pair->GetSecondDaughterP();
  Int_t labelD1 = (mcD1 ? TMath::Abs(mcD1->GetLabel()) : -1);
  Int_t labelD2 = (mcD2 ? TMath::Abs(mcD2->GetLabel()) : -1);
  Int_t d1Pdg = GetPdgFromLabel(labelD1);
  Int_t d2Pdg = GetPdgFromLabel(labelD2);

  // mothers
  AliVParticle* mcM1=0x0;
  AliVParticle* mcM2=0x0;

  // grand-mothers
  AliVParticle* mcG1 = 0x0;
  AliVParticle* mcG2 = 0x0;

  // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
  Bool_t directTerm = kTRUE;
  // daughters
  directTerm = directTerm && mcD1 && ComparePDG(d1Pdg,signalMC->GetLegPDG(1),signalMC->GetLegPDGexclude(1),signalMC->GetCheckBothChargesLegs(1))
               && CheckParticleSource(labelD1, signalMC->GetLegSource(1));

  directTerm = directTerm && mcD2 && ComparePDG(d2Pdg,signalMC->GetLegPDG(2),signalMC->GetLegPDGexclude(2),signalMC->GetCheckBothChargesLegs(2))
               && CheckParticleSource(labelD2, signalMC->GetLegSource(2));

  // mothers
  Int_t labelM1 = -1;
  if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    labelM1 = GetMothersLabel(labelD1);
    if(labelD1>-1 && labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
    directTerm = directTerm && (mcM1 || signalMC->GetMotherPDGexclude(1))
                 && ComparePDG((mcM1 ? mcM1->PdgCode() : 0),signalMC->GetMotherPDG(1),signalMC->GetMotherPDGexclude(1),signalMC->GetCheckBothChargesMothers(1))
                 && CheckParticleSource(labelM1, signalMC->GetMotherSource(1))
                 && CheckRadiativeDecision(labelM1,signalMC);
  }

  Int_t labelM2 = -1;
  if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    labelM2 = GetMothersLabel(labelD2);
    if(labelD2>-1 && labelM2>-1) mcM2 = GetMCTrackFromMCEvent(labelM2);
    directTerm = directTerm && (mcM2 || signalMC->GetMotherPDGexclude(2))
                 && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC->GetMotherPDG(2),signalMC->GetMotherPDGexclude(2),signalMC->GetCheckBothChargesMothers(2))
                 && CheckParticleSource(labelM2, signalMC->GetMotherSource(2))
                 && CheckRadiativeDecision(labelM2,signalMC);
  }

  // grand-mothers
  Int_t labelG1 = -1;
  if(signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    labelG1 = GetMothersLabel(labelM1);
    if(mcM1 && labelG1>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
    directTerm = directTerm && (mcG1 || signalMC->GetGrandMotherPDGexclude(1))
                 && ComparePDG((mcG1 ? mcG1->PdgCode() : 0),signalMC->GetGrandMotherPDG(1),signalMC->GetGrandMotherPDGexclude(1),signalMC->GetCheckBothChargesGrandMothers(1))
                 && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(1));
  }

  Int_t labelG2 = -1;
  if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    labelG2 = GetMothersLabel(labelM2);
    if(mcM2 && labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
    directTerm = directTerm && (mcG2 || signalMC->GetGrandMotherPDGexclude(2))
                 && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC->GetGrandMotherPDG(2),signalMC->GetGrandMotherPDGexclude(2),signalMC->GetCheckBothChargesGrandMothers(2))
                 && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(2));
  }

  // Cross term
  Bool_t crossTerm = kTRUE;
  // daughters
  crossTerm = crossTerm && mcD2
              && ComparePDG(d2Pdg,signalMC->GetLegPDG(1),signalMC->GetLegPDGexclude(1),signalMC->GetCheckBothChargesLegs(1))
              && CheckParticleSource(labelD2, signalMC->GetLegSource(1));

  crossTerm = crossTerm && mcD1
              && ComparePDG(d1Pdg,signalMC->GetLegPDG(2),signalMC->GetLegPDGexclude(2),signalMC->GetCheckBothChargesLegs(2))
              && CheckParticleSource(labelD1, signalMC->GetLegSource(2));

  // mothers
  if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    if(!mcM2 && labelD2>-1) {
      labelM2 = GetMothersLabel(labelD2);
      if(labelM2>-1) mcM2 = GetMCTrackFromMCEvent(labelM2);
    }
    crossTerm = crossTerm && (mcM2 || signalMC->GetMotherPDGexclude(1))
                && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC->GetMotherPDG(1),signalMC->GetMotherPDGexclude(1),signalMC->GetCheckBothChargesMothers(1))
                && CheckParticleSource(labelM2, signalMC->GetMotherSource(1))
                && CheckRadiativeDecision(labelM2,signalMC);
  }

  if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    if(!mcM1 && labelD1>-1) {
      labelM1 = GetMothersLabel(labelD1);
      if(labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
    }
    crossTerm = crossTerm && (mcM1 || signalMC->GetMotherPDGexclude(2))
                && ComparePDG((mcM1 ? mcM1->PdgCode() : 0),signalMC->GetMotherPDG(2),signalMC->GetMotherPDGexclude(2),signalMC->GetCheckBothChargesMothers(2))
                && CheckParticleSource(labelM1, signalMC->GetMotherSource(2))
                && CheckRadiativeDecision(labelM1,signalMC);
  }

  // grand-mothers
  if(signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    if(!mcG2 && mcM2) {
      labelG2 = GetMothersLabel(labelM2);
      if(labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
    }
    crossTerm = crossTerm && (mcG2 || signalMC->GetGrandMotherPDGexclude(1))
                && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC->GetGrandMotherPDG(1),signalMC->GetGrandMotherPDGexclude(1),signalMC->GetCheckBothChargesGrandMothers(1))
                && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(1));
  }

  if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    if(!mcG1 && mcM1) {
      labelG1 = GetMothersLabel(labelM1);
      if(labelG1>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
    }
    crossTerm = crossTerm && (mcG1 || signalMC->GetGrandMotherPDGexclude(2))
                && ComparePDG((mcG1 ? mcG1->PdgCode() : 0),signalMC->GetGrandMotherPDG(2),signalMC->GetGrandMotherPDGexclude(2),signalMC->GetCheckBothChargesGrandMothers(2))
                && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(2));
  }

  Bool_t motherRelation = kTRUE;
  if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kSame) {
    motherRelation = motherRelation && HaveSameMother(pair);
  }
  if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
    motherRelation = motherRelation && !HaveSameMother(pair);
  }
  // check geant process if set
  Bool_t processGEANT = kTRUE;
  if(signalMC->GetCheckGEANTProcess()) {
    if(!CheckGEANTProcess(labelD1,signalMC->GetGEANTProcess()) &&
       !CheckGEANTProcess(labelD2,signalMC->GetGEANTProcess())   ) processGEANT= kFALSE;
  }

  // check particle stack for pdg code
  Bool_t pdgInStack = kTRUE;
  if (signalMC->GetCheckStackForPDG()) {
    pdgInStack = CheckStackParticle(labelM1,signalMC->GetStackPDG()) && CheckStackParticle(labelM2,signalMC->GetStackPDG());
  }
  // check if a mother is also a grandmother
  Bool_t motherIsGrandmother = kTRUE;
  if (signalMC->GetCheckMotherGrandmotherRelation()) {
    motherIsGrandmother = kFALSE;
    motherIsGrandmother = MotherIsGrandmother(labelM1,labelM2,labelG1,labelG2, signalMC->GetMotherIsGrandmother());
  }

  return ((directTerm || crossTerm) && motherRelation && processGEANT && motherIsGrandmother && pdgInStack);

}

//___________________________________________________________
Bool_t AliDielectronMC::MotherIsGrandmother(int labelM1, int labelM2, int labelG1, int labelG2, bool motherIsGrandmother) const
{
  //
  // Check if the mother of one particle is the grandmother of the other
  //
  Bool_t result = ((labelM1 == labelG2) || (labelM2 == labelG1));
  if (motherIsGrandmother) return result;
  else return !result;

}

//____________________________________________________________
Bool_t AliDielectronMC::CheckStackParticle(Int_t labelPart, Int_t requiredPDG) const
{
  //
  // Check the stack of a particle and exclude if there is a certain pdg code found
  //
  Bool_t result = kTRUE;
  Int_t labelMother = GetMothersLabel(labelPart);
  Int_t motherPDG = GetPdgFromLabel(labelMother);
  Int_t i = 0;
  while(TMath::Abs(motherPDG) > 10 && labelMother > -1 && result && i<10){
    if (motherPDG == 0) return kFALSE;
    result = ComparePDG(motherPDG, requiredPDG, kTRUE, kTRUE);
    labelMother = GetMothersLabel(labelMother);
    motherPDG = GetPdgFromLabel(labelMother);
    i++;
  }
  return result;
}

//___________________________________________________________
Bool_t AliDielectronMC::CompareDaughterPDG(Int_t labelM, Int_t reqPDG, Bool_t PDGexclusion, Bool_t CheckBothChargesDaughter) const
{
  //
  // Check if one of the daughters has reqPDG
  //
  Bool_t result = kFALSE;
  //Get Mother from label
  AliMCParticle *mother=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM));
  //Get number auf daughters, so you can go through them
  const Int_t nd=(mother->GetLastDaughter()-mother->GetFirstDaughter()+1);
  //Loop over daughters and get stop loop if PDG is found.
  for (Int_t i=0; i<nd; ++i) {
    //Go through the daughters, set result true and break if PDG code fits
    if (ComparePDG((GetMCTrackFromMCEvent(mother->GetFirstDaughter()+i)->PdgCode()),reqPDG,kFALSE,CheckBothChargesDaughter)) {
      result = kTRUE;
      break;
    }
  }
  //if exclusion is needed do so
  if (PDGexclusion) return !result;
  return result;
}

//____________________________________________________________
Bool_t AliDielectronMC::HaveSameMother(const AliDielectronPair * pair) const
{
  //
  // Check whether two particles have the same mother
  //

  const AliVParticle * daughter1 = pair->GetFirstDaughterP();
  const AliVParticle * daughter2 = pair->GetSecondDaughterP();
  if (!daughter1 || !daughter2) return 0;

  AliVParticle *mcDaughter1=GetMCTrackFromMCEvent(daughter1->GetLabel());
  AliVParticle *mcDaughter2=GetMCTrackFromMCEvent(daughter2->GetLabel());
  if (!mcDaughter1 || !mcDaughter2) return 0;

  Int_t labelMother1=-1;
  Int_t labelMother2=-1;

  if (mcDaughter1->IsA()==AliMCParticle::Class()){
    labelMother1=(static_cast<AliMCParticle*>(mcDaughter1))->GetMother();
    labelMother2=(static_cast<AliMCParticle*>(mcDaughter2))->GetMother();
  } else if (mcDaughter1->IsA()==AliAODMCParticle::Class()) {
    labelMother1=(static_cast<AliAODMCParticle*>(mcDaughter1))->GetMother();
    labelMother2=(static_cast<AliAODMCParticle*>(mcDaughter2))->GetMother();
  }

  Bool_t sameMother=(labelMother1>-1)&&(labelMother2>-1)&&(labelMother1==labelMother2);

  return sameMother;
}

//________________________________________________________________
Int_t AliDielectronMC::IsJpsiPrimary(const AliDielectronPair * pair)
{
 // return: "0" for primary jpsi
 //         "1" for secondary jpsi (from beauty)
 //         "2" for background
 if(!HaveSameMother(pair)) return 2;
 AliVParticle *mcDaughter1=GetMCTrackFromMCEvent((pair->GetFirstDaughterP())->GetLabel());
 Int_t labelMother=-1;

  if (mcDaughter1->IsA()==AliMCParticle::Class()){
     labelMother=(static_cast<AliMCParticle*>(mcDaughter1))->GetMother();
     } else if (mcDaughter1->IsA()==AliAODMCParticle::Class()) {
     labelMother=(static_cast<AliAODMCParticle*>(mcDaughter1))->GetMother();
     }

 AliVParticle* mcMother=GetMCTrackFromMCEvent(labelMother);
 if(!IsMCMotherToEE(mcMother,443)) return 2;
 return IsJpsiPrimary(mcMother);
}

//______________________________________________________________
Int_t AliDielectronMC::IsJpsiPrimary(const AliVParticle * particle)
{
  // return: "0" for primary jpsi
  //         "1" for secondary jpsi (come from B decay)
 Int_t labelMoth=-1;
 Int_t pdgCode;

 if (particle->IsA()==AliMCParticle::Class()){
     labelMoth = (static_cast<const AliMCParticle*>(particle))->GetMother();
     while(labelMoth>0){
       particle = GetMCTrackFromMCEvent(labelMoth);
       pdgCode = TMath::Abs((static_cast<const AliMCParticle*>(particle))->PdgCode());
       if((pdgCode>500 && pdgCode<600) || (pdgCode>5000 && pdgCode<6000)) return 1;
       labelMoth = (static_cast<const AliMCParticle*>(particle))->GetMother();
       }
    }
 else if (particle->IsA()==AliAODMCParticle::Class()){
     labelMoth = (static_cast<const AliAODMCParticle*>(particle))->GetMother();
     while(labelMoth>0){
     particle = GetMCTrackFromMCEvent(labelMoth);
     pdgCode = TMath::Abs((static_cast<const AliAODMCParticle*>(particle))->PdgCode());
     if((pdgCode>500 && pdgCode<600) || (pdgCode>5000 && pdgCode<6000)) return 1;
     labelMoth = (static_cast<const AliAODMCParticle*>(particle))->GetMother();
     }
  }
  return 0;
}

//______________________________________________________________
Bool_t AliDielectronMC::GetPrimaryVertex(Double_t &primVtxX, Double_t &primVtxY, Double_t &primVtxZ)
{
  if(fAnaType == kESD){
    const AliVVertex* mcVtx =  fMCEvent->GetPrimaryVertex();
    if(!mcVtx) return kFALSE;
    primVtxX = mcVtx->GetX();
    primVtxY = mcVtx->GetY();
    primVtxZ = mcVtx->GetZ();
  }
  else if(fAnaType == kAOD){
    AliAODEvent *aod=((AliAODInputHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetEvent();
    if(!aod) return kFALSE;
    AliAODMCHeader *mcHead = dynamic_cast<AliAODMCHeader*>(aod->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!mcHead) return kFALSE;
    primVtxX = mcHead->GetVtxX();
    primVtxY = mcHead->GetVtxY();
    primVtxZ = mcHead->GetVtxZ();
  }
  return kTRUE;
}
