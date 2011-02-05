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
#include <AliStack.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliLog.h>

#include <TClonesArray.h>
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
Int_t AliDielectronMC::GetNPrimary()
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
AliVParticle* AliDielectronMC::GetMCTrackFromMCEvent(Int_t itrk)
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
    if (!fMcArray){ AliError("No fMCEvent"); return NULL;}
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
  if(fAnaType == kESD){
    fMCEvent=0x0;
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcHandler){ AliError("Could not retrive MC event handler!"); return kFALSE; }
    if (!mcHandler->InitOk() ) return kFALSE;
    if (!mcHandler->TreeK() )  return kFALSE;
    if (!mcHandler->TreeTR() ) return kFALSE;
    
    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if (!mcEvent){ AliError("Could not retrieve MC event!"); return kFALSE; }
    fMCEvent = mcEvent;
    
    if (!UpdateStack()) return kFALSE;
  }
  else if(fAnaType == kAOD)
  {
    fMcArray = 0x0;
    AliAODEvent *aod=((AliAODInputHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetEvent();
    fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMcArray) return kFALSE;
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

  AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
  return mctrack;
}


//____________________________________________________________
AliAODMCParticle* AliDielectronMC::GetMCTrack( const AliAODTrack* _track)
{
  //
  // return MC track
  //
 if(!fMcArray) { AliError("No fMCArray"); return NULL;}
 Int_t label = _track->GetLabel();
 if(label < 0) return NULL;
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

//____________________________________________________________
Bool_t AliDielectronMC::HaveSameMother(const AliDielectronPair * pair)
{
  //
  // Check whether two particles have the same mother
  //

  const AliVParticle * daughter1 = pair->GetFirstDaughter();
  const AliVParticle * daughter2 = pair->GetSecondDaughter();

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
