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

#include <TClonesArray.h>
#include <TParticle.h>

#include "AliDielectronSignalMC.h"
#include "AliDielectronMC.h"

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
AliVParticle* AliDielectronMC::GetMCTrackFromMCEvent(Int_t itrk) const
{
  //
  // return MC track directly from MC event
  //
  if (itrk<0) return NULL;
  AliVParticle * track=0x0;
  if(fAnaType == kESD){
    if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
    track = fMCEvent->GetTrack(itrk); //  tracks from MC event (ESD)
  } else if(fAnaType == kAOD) {
    if (!fMcArray){ AliError("No fMcArray"); return NULL;}
    if (itrk>fMcArray->GetEntriesFast()) { AliWarning(Form("track %d out of array size %d",itrk,fMcArray->GetEntriesFast())); return NULL;}
    track = (AliVParticle*)fMcArray->At(itrk); //  tracks from MC event (AOD)
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
  
  if(fAnaType == kESD){
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcHandler){ /*AliError("Could not retrive MC event handler!");*/ return kFALSE; }
    if (!mcHandler->InitOk() ) return kFALSE;
    if (!mcHandler->TreeK() )  return kFALSE;
    if (!mcHandler->TreeTR() ) return kFALSE;
    
    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if (!mcEvent){ /*AliError("Could not retrieve MC event!");*/ return kFALSE; }
    fMCEvent = mcEvent;
    
    if (!UpdateStack()) return kFALSE;
  }
  else if(fAnaType == kAOD)
  {
    AliAODInputHandler* aodHandler=(AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodHandler) return kFALSE;
    AliAODEvent *aod=aodHandler->GetEvent();
    if (!aod) return kFALSE;

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
  Int_t label = TMath::Abs(_track->GetLabel());
  if(label>nStack)return NULL;
  if(label<0) return NULL;
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
 Int_t label = _track->GetLabel();
 if(label < 0 || label > fMcArray->GetEntriesFast()) return NULL;
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
  
  AliMCParticle *mcPart1=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(particle1->GetLabel()));
  AliMCParticle *mcPart2=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(particle2->GetLabel()));
  
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
  AliAODMCParticle *mcPart1=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(particle1->GetLabel()));
  AliAODMCParticle *mcPart2=static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(particle2->GetLabel()));
  
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
Bool_t AliDielectronMC::CheckParticleSource(Int_t label, AliDielectronSignalMC::ESource source) const {
  //
  //  Check the source for the particle 
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
      return (label>=0 && label<=GetNPrimary());      
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
      return (label>=GetNPrimary() && !IsPhysicalPrimary(label));
    break;
    default :
      return kFALSE;
  }
  return kFALSE;
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
  const AliVParticle * mcD1 = pair->GetFirstDaughter();
  const AliVParticle * mcD2 = pair->GetSecondDaughter();
  Int_t labelD1 = (mcD1 ? mcD1->GetLabel() : -1);
  Int_t labelD2 = (mcD2 ? mcD2->GetLabel() : -1);
  if(labelD1<0) labelD1 *= -1;
  if(labelD2<0) labelD2 *= -1;
  Int_t d1Pdg = 0;
  Int_t d2Pdg = 0;
  d1Pdg=GetPdgFromLabel(labelD1);
  d2Pdg=GetPdgFromLabel(labelD2);
  
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
      if(labelG2>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
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
 
  return ((directTerm || crossTerm) && motherRelation);
}



//____________________________________________________________
Bool_t AliDielectronMC::HaveSameMother(const AliDielectronPair * pair) const
{
  //
  // Check whether two particles have the same mother
  //

  const AliVParticle * daughter1 = pair->GetFirstDaughter();
  const AliVParticle * daughter2 = pair->GetSecondDaughter();
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
 AliVParticle *mcDaughter1=GetMCTrackFromMCEvent((pair->GetFirstDaughter())->GetLabel());
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


Bool_t AliDielectronMC::GetPrimaryVertex(Double_t &primVtxX, Double_t &primVtxY, Double_t &primVtxZ){

     if(fAnaType == kESD){
     const AliVVertex* mcVtx =  fMCEvent->GetPrimaryVertex();
     if(!mcVtx) return kFALSE;
     primVtxX = mcVtx->GetX();
     primVtxY = mcVtx->GetY();
     primVtxZ = mcVtx->GetZ();
     }else if(fAnaType == kAOD){
     AliAODEvent *aod=((AliAODInputHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetEvent();
     AliAODMCHeader *mcHead = dynamic_cast<AliAODMCHeader*>(aod->FindListObject(AliAODMCHeader::StdBranchName()));
     primVtxX = mcHead->GetVtxX();
     primVtxY = mcHead->GetVtxY();
     primVtxZ = mcHead->GetVtxZ();
     }

return kTRUE;
}
