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
  Bool_t checkHF=kFALSE;
  if (AliAnalysisManager::GetAnalysisManager()){
    if (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()) type=kESD;
    else if (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()) type=kAOD;

    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(type == kESD) hasMC=mcHandler!=0x0;
    else if (type == kAOD) hasMC=mcHandler!=0x0;
    }

  fgInstance=new AliDielectronMC(type);

  fgInstance->SetHasMC(hasMC);

  fgInstance->SetCheckHF(checkHF);

  return fgInstance;
}

//____________________________________________________________
AliDielectronMC::AliDielectronMC(AnalysisType type):
  fMCEvent(0x0),
  fAODMCHeader(0x0),
  fGenCocktailHeader(0x0),
  fAnaType(type),
  fHasMC(kTRUE),
  fCheckHF(kFALSE),
  fhfproc(),
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
  // if(fAnaType == kESD){
  //   if (!fMCEvent){ AliError("No fMCEvent"); return 0; }
  //   return fMCEvent->GetNumberOfTracks();}
  // else if(fAnaType == kAOD){
  //   if(!fMcArray) { AliError("No fMcArray"); return 0; }
  //   return fMcArray->GetEntriesFast();
  // }
  // return 0;
  // Splitting between ESDs and AODs not needed anymore should work fine this way. PD 2017-10-05
  if(!fMCEvent){
    AliError("No fMCEvent");
    return -999;
  }
  else return fMCEvent->GetNumberOfTracks();
}

//____________________________________________________________
Int_t AliDielectronMC::GetNMCTracksFromStack()
{
  //
  //  return the number of generated tracks from stack
  //
  if (!fMCEvent){ AliError("No fMCEvent"); return -999; }
  return fMCEvent->GetNumberOfTracks();
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
  if (!fMCEvent){ AliError("No fMCEvent"); return -999; }
  return fMCEvent->GetNumberOfPrimaries();
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
  if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
  track = fMCEvent->GetTrack(label); //  tracks from MC event
  //Splitting between ESD and AOD should be obsolete 2017-10-05 PD
  // if(fAnaType == kESD){
  //   if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
  //   track = fMCEvent->GetTrack(label); //  tracks from MC event (ESD)
  // } else if(fAnaType == kAOD) {
  //   if (!fMcArray){ AliError("No fMcArray"); return NULL;}
  //   if (label>fMcArray->GetEntriesFast()) { AliWarning(Form("track %d out of array size %d",label,fMcArray->GetEntriesFast())); return NULL;}
  //   track = (AliVParticle*)fMCEvent->GetTrack(label);
  // }
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

    if (fCheckHF){
      fhfproc.clear();
      fCheckHF=LoadHFPairs(); // So far only compatible with ESD
    }
  }
  else if(fAnaType == kAOD)
  {
    AliAODInputHandler* aodHandler=(AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!aodHandler) return kFALSE;
    AliAODEvent *aod=aodHandler->GetEvent();
    if (!aod) return kFALSE;

    fMCEvent = aodHandler->MCEvent();
    fAODMCHeader = (AliAODMCHeader*) aod->FindListObject(AliAODMCHeader::StdBranchName());
    if (!fMCEvent) AliError("No MCEvent available");

    fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMcArray){ /*AliError("Could not retrieve MC array!");*/ return kFALSE; }
    else fHasMC=kTRUE;
  }
  return kTRUE;
}

//____________________________________________________________
AliVParticle* AliDielectronMC::GetMCTrack(const AliVParticle* _track)
{
  //
  // return MC track without in/output casting
  //
  if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
  Int_t nStack = fMCEvent->GetNumberOfTracks();
  Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
  if(label > nStack)  return NULL;
  AliVParticle *mctrack = (fMCEvent->GetTrack(label));
  return mctrack;
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
  if(!fMCEvent) { AliError("No fMCEvent"); return NULL;}
  Int_t nStack = fMCEvent->GetNumberOfTracks();
  Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
  if(label > nStack) return NULL;
  AliAODMCParticle *mctrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
  return mctrack;
}

//____________________________________________________________
TParticle* AliDielectronMC::GetMCTrackFromStack(const AliESDtrack* _track)
{
  //
  // return MC track from MC event
  //
   if (!fMCEvent){ AliError("No fMCEvent"); return NULL;}
   Int_t nStack = fMCEvent->GetNumberOfTracks();
   Int_t label  = TMath::Abs(_track->GetLabel());
   if(label > nStack) return NULL;
   TParticle* mcpart = fMCEvent->Particle(label);
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

  return AliDielectronMC::GetMCTrackMother(mcpart);
  // if(mcpart->GetMother() < 0) return NULL;
  // AliAODMCParticle* mcmother = dynamic_cast<AliAODMCParticle *>(fMcArray->At(mcpart->GetMother()));
  // if (!mcmother) return NULL;
  // return mcmother;
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
  AliAODMCParticle* mcmother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(_particle->GetMother()));
  return mcmother;
}

//____________________________________________________________
//TParticle* AliDielectronMC::GetMCTrackMotherFromStack(const AliESDtrack* _track)
//{
  //
  // return MC track mother from MC event
  //
//  TParticle* mcpart = GetMCTrackFromStack(_track);
//  if ( !mcpart || mcpart->GetFirstMother()<=0 ) return NULL;
//  TParticle* mcmother = fMCEvent->Particle(mcpart->GetFirstMother());
//  if (!mcmother) return NULL;
//  return mcmother;
//}

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
  // return MC PDG code from MC event
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
//Int_t AliDielectronMC::GetMotherPDGFromStack(const AliESDtrack* _track)
//{
  //
  // return PDG code of the mother track from MC event
  //
//  TParticle* mcmother = GetMCTrackMotherFromStack(_track);
//  if (!mcmother) return -999;
//  return mcmother->GetPdgCode();
//}

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
//Int_t AliDielectronMC::GetMCProcessFromStack(const AliESDtrack* _track)
//{
  //
  // return process number of the track
  //
  //TParticle* mcpart = GetMCTrackFromStack(_track);
  //if (!mcpart) return -999;
  //return mcpart->GetUniqueID();
//}

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
//Int_t AliDielectronMC::GetMCProcessMotherFromStack(const AliESDtrack* _track)
//{
  //
  // return process number of the mother of the track
  //
//  TParticle* mcmother = GetMCTrackMotherFromStack(_track);
//  if (!mcmother) return -999;
//  return mcmother->GetUniqueID();
//}

//____________________________________________________________
Bool_t AliDielectronMC::IsMCMotherToEE(const AliVParticle *particle, Int_t pdgMother)
{
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  //
  if(!fMCEvent) return kFALSE;
  // if (fAnaType==kESD && !fMCEvent) return kFALSE;
  // if (fAnaType==kAOD && !fMcArray) return kFALSE;
  if (!particle) return kFALSE;

  if (!(particle->IsA()==AliMCParticle::Class() || particle->IsA()==AliAODMCParticle::Class())){
    AliError("Unknown particle type");
    return kFALSE;
  }
  //check pdg code
  if (particle->PdgCode()!=pdgMother) return kFALSE;
  Int_t ifirst = particle->GetDaughterFirst();
  Int_t ilast  = particle->GetDaughterLast();

  //check number of daughters
  if ((ilast-ifirst)!=1) return kFALSE;
  AliVParticle *firstD = (AliVParticle*) (GetMCTrackFromMCEvent(ifirst));
  AliVParticle *secondD = (AliVParticle*) (GetMCTrackFromMCEvent(ilast));

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
Bool_t AliDielectronMC::IsMCMotherToEEesd(const AliMCParticle *particle, Int_t pdgMother)
{
  // Obsolete function case splitting not needed 2017-10-05 PD
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  // ESD case
  //

  //check pdg code
  if (particle->PdgCode()!=pdgMother) return kFALSE;
  Int_t ifirst = particle->GetDaughterFirst();
  Int_t ilast  = particle->GetDaughterLast();

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
  // Obsolete function case splitting not needed 2017-10-05 PD
  //
  // Check if the Mother 'particle' is of type pdgMother and decays to e+e-
  // AOD case
  //

  if (particle->GetPdgCode()!=pdgMother) return kFALSE;
  if (particle->GetNDaughters()!=2) return kFALSE;

  Int_t ifirst = particle->GetDaughterFirst();
  Int_t ilast  = particle->GetDaughterLast();

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
  //TODO: check how you can get rid of the hardcoded numbers. One should make use of the PdgCodes set in AliDielectron!!!
  //
  if (!fMCEvent) return -1;

  Int_t lblMother1=particle1->GetMother();
  Int_t lblMother2=particle2->GetMother();

  AliVParticle *mcMother1 = (AliVParticle*) GetMCTrackFromMCEvent(lblMother1);

  if (!mcMother1) return -1;
  if (lblMother1!=lblMother2) return -1;
  if (TMath::Abs(particle1->PdgCode())!=11) return -1;
  if (particle1->PdgCode()!=-particle2->PdgCode()) return -1;
  if (mcMother1->PdgCode()!=pdgMother) return -1;

  return lblMother1;

  // AOD ESD case splitting is obsolete 2017-10-05 PD
  // if (fAnaType==kESD) return GetLabelMotherWithPdgESD(particle1, particle2, pdgMother);
  // if (fAnaType==kAOD) return GetLabelMotherWithPdgAOD(particle1, particle2, pdgMother);

  // return -1;
}

//____________________________________________________________
Int_t AliDielectronMC::GetLabelMotherWithPdgESD(const AliVParticle *particle1, const AliVParticle *particle2, Int_t pdgMother)
{
  //Obsolete function 2017-10-05 PD
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
  //Obsolete function 2017-10-05 PD
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

  if (!fMCEvent){
    AliError("No fMCEvent - break!");
    return;
  }
  AliVParticle *mcMother = (AliVParticle*) mother;
  lblD1 = mcMother->GetDaughterFirst();
  lblD2 = mcMother->GetDaughterLast();
  d1 = fMCEvent->GetTrack(lblD1);
  d2 = fMCEvent->GetTrack(lblD2);

}


//________________________________________________________________________________
Int_t AliDielectronMC::GetMothersLabel(Int_t daughterLabel) const {
  //
  //  Get the label of the mother for particle with label daughterLabel
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(daughterLabel<0) return -1;
  if (!fMCEvent) return -1;
  Int_t momsLabel = GetMCTrackFromMCEvent(daughterLabel)->GetMother();
  return momsLabel;
}

//________________________________________________________________________________
Int_t AliDielectronMC::GetFirstMothersLabelInChain(Int_t daughterLabel) const {
  //
  //  Get the label of the very first mother in the decay chain for particle with label daughterLabel
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(daughterLabel<0) return -1;
  if (!fMCEvent) return -1;

  Int_t labelfirstmother = -1;
  Int_t firstMotherIndex = GetMCTrackFromMCEvent(daughterLabel)->GetMother(); // get label of the mother

  while(firstMotherIndex>0){
    labelfirstmother = firstMotherIndex;
    firstMotherIndex = GetMCTrackFromMCEvent(labelfirstmother)->GetMother(); // get label of mother of mother....
  }

  if(firstMotherIndex==0) labelfirstmother = firstMotherIndex; // In case of PYTHIA first mother is 0 or 1 but should not happen in Hijing
  // labelfirstmother is -1 if particle has no mother (primary) or equaled to the very first primary particle in the chain
  return labelfirstmother;

}
//________________________________________________________________________________
Int_t AliDielectronMC::GetMinLabelParticlePrimaryAround(Int_t label) const {
  //
  //  Get the label of the last primary particle around the primary particle at the label
  //  NOTE: for tracks, the absolute label should be passed
  //

  if(label<0) return -1;
  if (!fMCEvent) return -1;

  Int_t labelminfirstMother = label; // start at the particle label
  // get label of the mother: should be -1 since we assume that this is the first mother in the chain
  Int_t firstMotherIndex = GetMCTrackFromMCEvent(label)->GetMother();
  if(firstMotherIndex>-1) return -1;

  while(firstMotherIndex<0){
    labelminfirstMother--;
    if(labelminfirstMother<0){
      firstMotherIndex = 0;
    }
    else{
      firstMotherIndex = GetMCTrackFromMCEvent(labelminfirstMother)->GetMother(); // get the mother label of the neighborhood particle
    }
  }

  labelminfirstMother ++; // set back by one to the last primary particle in the - neighborhood

  return labelminfirstMother;

}
//________________________________________________________________________________
Int_t AliDielectronMC::GetMaxLabelParticlePrimaryAround(Int_t label) const {
  //
  //  Get the label of the last primary particle around the primary particle at the label
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(label<0) return -1;
  if (!fMCEvent) return -1;

  Int_t labelmaxfirstMother = label; // start at the particle label
  // get label of the mother: should be -1 since we assume that this is the first mother in the chain
  Int_t firstMotherIndex = GetMCTrackFromMCEvent(label)->GetMother();
  if(firstMotherIndex>-1) return -1;

  while(firstMotherIndex<0){
    labelmaxfirstMother++;
    if(labelmaxfirstMother > fMCEvent->GetNumberOfTracks()){
      firstMotherIndex = 0;
    }
    else{
      firstMotherIndex = GetMCTrackFromMCEvent(labelmaxfirstMother)->GetMother(); // get the mother label of the neighborhood particle
    }
  }

  labelmaxfirstMother --; // set back by one to the last primary particle in the + neighborhood

  return labelmaxfirstMother;

}
//________________________________________________________________________________
Int_t AliDielectronMC::GetPdgFromLabel(Int_t label) const {
  //
  //  Get particle code using the label from MC event
  //  NOTE: for tracks, the absolute label should be passed
  //
  if(label<0) return 0;
  if (!fMCEvent) return 0;
  return GetMCTrackFromMCEvent(label)->PdgCode();
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
  case 600:     // all dielectron from same mother sources
    result = /*TMath::Abs(particlePDG)==22  ||*/ // photon
             TMath::Abs(particlePDG)==111 || // pion
             TMath::Abs(particlePDG)==113 || // rho
             TMath::Abs(particlePDG)==221 || // eta
             TMath::Abs(particlePDG)==331 || // eta'
             TMath::Abs(particlePDG)==223 || // omega
             TMath::Abs(particlePDG)==333 || // phi
             TMath::Abs(particlePDG)==443 || // jpsi
             TMath::Abs(particlePDG)==100443 // psi 2S
             ;
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

  if (!fMCEvent) return kFALSE;
  if(label<0) return kFALSE;
  return GetMCTrackFromMCEvent(label)->IsPhysicalPrimary();
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsPrimary(Int_t label) const {
  //
  if(label<0) return kFALSE;
  if(!fMCEvent) return kFALSE;
  if(fAnaType==kAOD) {
    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsPrimary();
  } else if(fAnaType==kESD) {
    return (label>=0 && label<=GetNPrimary());
  }
  return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondary(Int_t label) const {
  //
  if(label<0) return kFALSE;
  if (!fMCEvent) return kFALSE;
  if(IsPrimary(label)) return kFALSE;
  if(IsPhysicalPrimary(label)) return kFALSE;

  return kTRUE;

  // Old implementation different approach for AODs and ESDs not needed anymore -- 2017-10-05 PD
  // if(fAnaType==kAOD) {
  //   if(!fMcArray) return kFALSE;
  //   AliAODMCParticle* mctrack = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label));
  //   Bool_t isSecondary = mctrack->IsSecondaryFromMaterial() || mctrack->IsSecondaryFromWeakDecay();
  //   return isSecondary;
  // } else if(fAnaType==kESD) {
  //   if (!fMCEvent) return kFALSE;
  //   return (label>=GetNPrimary() && !IsPhysicalPrimary(label));
  // }
  // return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondaryFromWeakDecay(Int_t label) const {
  //
  // Check if the particle with label "label" is a physical secondary from weak decay according to the
  // definition in AliStack::IsSecondaryFromWeakDecay(Int_t label)
  //
  if(label<0) return kFALSE;
  if(!fMCEvent) return kFALSE;

  return ((AliVParticle*) GetMCTrackFromMCEvent(label))->IsSecondaryFromWeakDecay();

  // Old implementation different approach for AODs and ESDs not needed anymore -- 2017-10-05 PD

  // if(fAnaType==kAOD) {
  //   if(!fMcArray) return kFALSE;
  //   return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsSecondaryFromWeakDecay();
  // } else if(fAnaType==kESD) {
  //   if(!fMCEvent) return kFALSE;
  //   return fMCEvent->IsSecondaryFromWeakDecay(label);
  // }
  // return kFALSE;
}

//________________________________________________________________________________
Bool_t AliDielectronMC::IsSecondaryFromMaterial(Int_t label) const {
  //
  // Check if the particle with label "label" is a physical secondary from weak decay according to the
  // definition in AliStack::IsSecondaryFromMaterial(Int_t label)
  //
  if(label<0) return kFALSE;
  if(!fMCEvent) return kFALSE;

  return ((AliVParticle*) GetMCTrackFromMCEvent(label))->IsSecondaryFromMaterial();

  // Old implementation different approach for AODs and ESDs not needed anymore -- 2017-10-05 PD
  // if(fAnaType==kAOD) {
  //   if(!fMcArray) return kFALSE;
  //   return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsSecondaryFromMaterial();
  // } else if(fAnaType==kESD) {
  //   if (!fMCEvent) return kFALSE;
  //   return fMCEvent->IsSecondaryFromMaterial(label);
  // }
  // return kFALSE;
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
      // 3.) Certain particles added via MC generator cocktails (e.g. J/psi added to pythia MB events
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

//___________________________________________________________________
Bool_t AliDielectronMC::CheckParticleSource(const AliAODMCParticle *mcPart, AliDielectronSignalMC::ESource source) const {
  //
  //  Check the source for the particle mcPart
  //  NOTE: preferable function for AODs, since one anyhow needs the particle to check it's source and getting the correct label is not always do able.
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
      return mcPart->IsPrimary();
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
      return mcPart->IsPhysicalPrimary();
    break;
    case AliDielectronSignalMC::kDirect :
      // Primary particles which do not have any mother
      // This is the case for:
      // 1.) Initial state particles (the 2 protons in Pythia pp collisions)
      // 2.) In some codes, with sudden freeze-out, all particles generated from the fireball are direct.
      //     There is no history for these particles.
      // 3.) Certain particles added via MC generator cocktails (e.g. J/psi added to pythia MB events)
      return (mcPart->GetMother() <= 0);
      break;
    case AliDielectronSignalMC::kNoCocktail :
      // Particles from the HIJING event and NOT from the AliGenCocktail
      return (mcPart->GetMother() > 0);
      break;
    case AliDielectronSignalMC::kSecondary :
      // particles which are created by the interaction of final state primaries with the detector
      // or particles from strange weakly decaying particles (e.g. lambda, kaons, etc.)
      return (!mcPart->IsPrimary() && !mcPart->IsPhysicalPrimary());
      // return (mcPart->IsSecondaryFromMaterial() || mcPart->IsSecondaryFromWeakDecay()); // potential alternative definition of "secondary". Should be the same but is not tested
      // return (label>=GetNPrimary() && !IsPhysicalPrimary(label)); // old definition
    break;
    case AliDielectronSignalMC::kSecondaryFromWeakDecay :
      // secondary particle from weak decay
      // or particles from strange weakly decaying particles (e.g. lambda, kaons, etc.)
      return (mcPart->IsSecondaryFromWeakDecay());
    break;
    case AliDielectronSignalMC::kSecondaryFromMaterial :
      // secondary particle from material
      return (mcPart->IsSecondaryFromMaterial());
    break;
    case AliDielectronSignalMC::kFromBGEvent :
      // Use original esd track label to compare with nProduced particles
      // used to select electrons which are not from injected signals.
      return (IsFromBGEvent(mcPart->GetLabel()));
      break;
    case AliDielectronSignalMC::kFinalStateFromBGEvent :
      // Use original esd track label to compare with nProduced particles
      // used to select electrons which are not from injected signals.
      return (mcPart->IsPhysicalPrimary() && IsFromBGEvent(mcPart->GetLabel()));
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

  if (!fMCEvent) return kFALSE;
  if(!CheckHijingHeader()){
    AliError("No Hijing generator header - selection of IsFromBGEvent not possible!");
    return kFALSE;
  }
  Int_t nBgrdParticles = fGenCocktailHeader->NProduced();
  return (label < nBgrdParticles);


// Splitting between AOD-ESD analysis is obsolete. 2017-10-05 PD

//   if(fAnaType==kAOD) {
//     AliWarning("IsFromBGEvent() not implemented for AOD!");
//     return kFALSE;
// //    if(!fMcArray) return kFALSE;
// //    return (static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(label)))->IsFromBGEvent(); // IsFromBGEvent() does not exist for AliAODMCParticle.
//   } else if(fAnaType==kESD) {
//     if (!fMCEvent) return kFALSE;
//     if (CheckHijingHeader()) return fMCEvent->IsFromBGEvent(label); // Works for HIJING inside Cocktail
//     //else if (CheckSomeOtherHeader()) return ...;
//     else {
//       AliWarning("No headers to make decision! Assuming no injected signals are present.");
//       return kTRUE;
//     }
//   }
//   return kFALSE;
}


//________________________________________________________________________________
Bool_t AliDielectronMC::CheckHijingHeader() const {

//  if (fHasHijingHeader > -1) return Bool_t(fHasHijingHeader); // avoid many calls of the code below.

  if(fAnaType==kAOD) {
    if(!fAODMCHeader) AliError("No AODMCHeader - Break!");
    TList* headerList = fAODMCHeader->GetCocktailHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      if ((headerList->At(i))->IsA() == AliGenHijingEventHeader::Class()){
        fGenCocktailHeader = (AliGenHijingEventHeader*) headerList->At(i);
        return (fHasHijingHeader = kTRUE);
      }
    }
    return (fHasHijingHeader = kFALSE);
  }

  else if(fAnaType==kESD) {
    // taken from AliMCEvent::IsFromBGEvent()
    if (!fMCEvent) return (fHasHijingHeader=0); //return kFALSE;
    AliGenCocktailEventHeader* coHeader = dynamic_cast<AliGenCocktailEventHeader*> (fMCEvent->GenEventHeader());
    if (!coHeader) return (fHasHijingHeader=0); //return kFALSE;
    TList* list = coHeader->GetHeaders();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(list->FindObject("Hijing"));
    fGenCocktailHeader = hijingH;
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
      if (GetMCTrackFromMCEvent(mother->GetDaughterLabel(0)+i)->PdgCode()!=22) return kFALSE; //last daughter is photon
  } else if(fAnaType==kESD) {
    if (!fMCEvent) return kFALSE;
    AliMCParticle *mother=static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(label));
    if (!mother) return kFALSE;
    const Int_t nd=(mother->GetDaughterLast()-mother->GetDaughterFirst()+1);
    if (nd==2) return kFALSE;
    for (Int_t i=2; i<nd; ++i)
      if (GetMCTrackFromMCEvent(mother->GetDaughterFirst()+i)->PdgCode()!=22) return kFALSE; //last daughters are photons
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
  if(!ComparePDG(part->PdgCode(), signalMC->GetLegPDG(branch), signalMC->GetLegPDGexclude(branch), signalMC->GetCheckBothChargesLegs(branch))) return kFALSE;
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

    if(!ComparePDG((mcMother ? mcMother->PdgCode() : 0), signalMC->GetMotherPDG(branch), signalMC->GetMotherPDGexclude(branch), signalMC->GetCheckBothChargesMothers(branch))) return kFALSE;
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
  directTerm = directTerm && mcD1 && ComparePDG(d1Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD1, signalMC->GetLegSource(1));

  directTerm = directTerm && mcD2 && ComparePDG(d2Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2))
               && CheckParticleSource(labelD2, signalMC->GetLegSource(2));

  // mothers
  Int_t labelM1 = -1;
  if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    labelM1 = GetMothersLabel(labelD1);
    if(labelD1>-1 && labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
    directTerm = directTerm && (mcM1 || signalMC->GetMotherPDGexclude(1))
                 && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1),signalMC->GetCheckBothChargesMothers(1))
                 && CheckParticleSource(labelM1, signalMC->GetMotherSource(1))
                 && CheckRadiativeDecision(labelM1, signalMC);
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
                 && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1))
                 && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(1));
  }

  Int_t labelG2 = -1;
  if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    labelG2 = GetMothersLabel(labelM2);
    if(mcM2 && labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
    directTerm = directTerm && (mcG2 || signalMC->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC->GetGrandMotherPDG(2),signalMC->GetGrandMotherPDGexclude(2),signalMC->GetCheckBothChargesGrandMothers(2))
                 && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(2));
  }

  // Cross term
  Bool_t crossTerm = kTRUE;
  // daughters
  crossTerm = crossTerm && mcD2 && ComparePDG(d2Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD2, signalMC->GetLegSource(1));
  crossTerm = crossTerm && mcD1 && ComparePDG(d1Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2))
              && CheckParticleSource(labelD1, signalMC->GetLegSource(2));

  // mothers
  if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    if(!mcM2 && labelD2>-1) {
      labelM2 = GetMothersLabel(labelD2);
      if(labelM2>-1) mcM2 = GetMCTrackFromMCEvent(labelM2);
    }
    crossTerm = crossTerm && (mcM2 || signalMC->GetMotherPDGexclude(1)) && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1), signalMC->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM2, signalMC->GetMotherSource(1)) && CheckRadiativeDecision(labelM2,signalMC);
  }

  if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    if(!mcM1 && labelD1>-1) {
      labelM1 = GetMothersLabel(labelD1);
      if(labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
    }
    crossTerm = crossTerm && (mcM1 || signalMC->GetMotherPDGexclude(2)) && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(2), signalMC->GetMotherPDGexclude(2), signalMC->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM1, signalMC->GetMotherSource(2)) && CheckRadiativeDecision(labelM1,signalMC);
  }

  // grand-mothers
  if(signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
    if(!mcG2 && mcM2) {
      labelG2 = GetMothersLabel(labelM2);
      if(labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
    }
    crossTerm = crossTerm && (mcG2 || signalMC->GetGrandMotherPDGexclude(1)) && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(1));
  }

  if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
    if(!mcG1 && mcM1) {
      labelG1 = GetMothersLabel(labelM1);
      if(labelG1>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
    }
    crossTerm = crossTerm && (mcG1 || signalMC->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(2), signalMC->GetGrandMotherPDGexclude(2), signalMC->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(2));
  }

  Bool_t motherRelation = kTRUE;
  if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kSame) {
    motherRelation = motherRelation && HaveSameMother(pair);
  }
  if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
    motherRelation = motherRelation && !HaveSameMother(pair);
  }

  Bool_t grandMotherRelation = kTRUE;
  if(signalMC->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
    grandMotherRelation = grandMotherRelation && HaveSameGrandMother(pair);
  }
  if(signalMC->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
    grandMotherRelation = grandMotherRelation && !HaveSameGrandMother(pair);
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
  // check is first mothers are closed (correlated charm/beauty)
  Bool_t correlated = kTRUE;
  if (signalMC->GetCheckCorrelatedHF()){
    Int_t labelfirstmother1 = GetFirstMothersLabelInChain(labelD1);
    Int_t labelminfirstmother1 = GetMinLabelParticlePrimaryAround(labelfirstmother1);
    Int_t labelmaxfirstmother1 = GetMaxLabelParticlePrimaryAround(labelfirstmother1);
    Int_t labelfirstmother2 = GetFirstMothersLabelInChain(labelD2);
    Int_t labelminfirstmother2 = GetMinLabelParticlePrimaryAround(labelfirstmother2);
    Int_t labelmaxfirstmother2 = GetMaxLabelParticlePrimaryAround(labelfirstmother2);
    if(labelfirstmother1==-1 || labelminfirstmother1==-1 || labelmaxfirstmother1==-1) correlated = kFALSE; // security
    if(labelfirstmother2==-1 || labelminfirstmother2==-1 || labelmaxfirstmother2==-1) correlated = kFALSE; // security
    if(labelfirstmother1<labelminfirstmother2 || labelfirstmother1>labelmaxfirstmother2) correlated = kFALSE;
    if(labelfirstmother2<labelminfirstmother1 || labelfirstmother2>labelmaxfirstmother1) correlated = kFALSE;
  }

  // check if a mother is also a grandmother
  Bool_t motherIsGrandmother = kTRUE;
  if (signalMC->GetCheckMotherGrandmotherRelation()) {
    motherIsGrandmother = kFALSE;
    motherIsGrandmother = MotherIsGrandmother(labelM1,labelM2,labelG1,labelG2, signalMC->GetMotherIsGrandmother());
  }

  return ((directTerm || crossTerm) && motherRelation && grandMotherRelation && processGEANT && motherIsGrandmother && pdgInStack && correlated);

}



  //________________________________________________________________________________
  // Define IsMCTruth for 2 pairs
Bool_t AliDielectronMC::IsMCTruth(const AliDielectronPair* pair1, const AliDielectronPair* pair2, const AliDielectronSignalMC* signalMC1, const AliDielectronSignalMC* signalMC2) const {
    //
    // Check if the pair corresponds to the MC truth in signalMC
    //
    // legs (daughters)
    const AliVParticle * mcD1 = pair1->GetFirstDaughterP();
    const AliVParticle * mcD2 = pair1->GetSecondDaughterP();
    const AliVParticle * mcD3 = pair2->GetFirstDaughterP();
    const AliVParticle * mcD4 = pair2->GetSecondDaughterP();
    Int_t labelD1 = (mcD1 ? TMath::Abs(mcD1->GetLabel()) : -1);
    Int_t labelD2 = (mcD2 ? TMath::Abs(mcD2->GetLabel()) : -1);
    Int_t labelD3 = (mcD3 ? TMath::Abs(mcD3->GetLabel()) : -1);
    Int_t labelD4 = (mcD4 ? TMath::Abs(mcD4->GetLabel()) : -1);
    Int_t d1Pdg = GetPdgFromLabel(labelD1);
    Int_t d2Pdg = GetPdgFromLabel(labelD2);
    Int_t d3Pdg = GetPdgFromLabel(labelD3);
    Int_t d4Pdg = GetPdgFromLabel(labelD4);

    // mothers
    AliVParticle* mcM1=0x0;
    AliVParticle* mcM2=0x0;
    AliVParticle* mcM3=0x0;
    AliVParticle* mcM4=0x0;

    // grand-mothers
    AliVParticle* mcG1 = 0x0;
    AliVParticle* mcG2 = 0x0;
    AliVParticle* mcG3 = 0x0;
    AliVParticle* mcG4 = 0x0;

    // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
    Bool_t directTerm = kTRUE;
    // daughters
    directTerm = directTerm && mcD1 && ComparePDG(d1Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD1, signalMC1->GetLegSource(1));
    directTerm = directTerm && mcD2 && ComparePDG(d2Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD2, signalMC1->GetLegSource(2));
    directTerm = directTerm && mcD3 && ComparePDG(d3Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD3, signalMC2->GetLegSource(1));
    directTerm = directTerm && mcD4 && ComparePDG(d4Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD4, signalMC2->GetLegSource(2));


    // mothers
    Int_t labelM1 = -1;
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM1 = GetMothersLabel(labelD1);
      if(labelD1>-1 && labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
      directTerm = directTerm && (mcM1 || signalMC1->GetMotherPDGexclude(1))
                   && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1),signalMC1->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(labelM1, signalMC1->GetMotherSource(1))
                   && CheckRadiativeDecision(labelM1, signalMC1);
    }

    Int_t labelM2 = -1;
    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM2 = GetMothersLabel(labelD2);
      if(labelD2>-1 && labelM2>-1) mcM2 = GetMCTrackFromMCEvent(labelM2);
      directTerm = directTerm && (mcM2 || signalMC1->GetMotherPDGexclude(2))
                   && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC1->GetMotherPDG(2),signalMC1->GetMotherPDGexclude(2),signalMC1->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(labelM2, signalMC1->GetMotherSource(2))
                   && CheckRadiativeDecision(labelM2,signalMC1);
    }

    Int_t labelM3 = -1;
    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM3 = GetMothersLabel(labelD3);
      if(labelD3>-1 && labelM3>-1) mcM3 = GetMCTrackFromMCEvent(labelM3);
      directTerm = directTerm && (mcM3 || signalMC2->GetMotherPDGexclude(1))
                   && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1),signalMC2->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(labelM3, signalMC2->GetMotherSource(1))
                   && CheckRadiativeDecision(labelM3, signalMC2);
    }

    Int_t labelM4 = -1;
    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM4 = GetMothersLabel(labelD4);
      if(labelD4>-1 && labelM4>-1) mcM4 = GetMCTrackFromMCEvent(labelM4);
      directTerm = directTerm && (mcM4 || signalMC2->GetMotherPDGexclude(2))
                   && ComparePDG((mcM4 ? mcM4->PdgCode() : 0),signalMC2->GetMotherPDG(2),signalMC2->GetMotherPDGexclude(2),signalMC2->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(labelM4, signalMC2->GetMotherSource(2))
                   && CheckRadiativeDecision(labelM4,signalMC2);
    }

    // grand-mothers
    Int_t labelG1 = -1;
    if(signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelG1 = GetMothersLabel(labelM1);
      if(mcM1 && labelG1>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
      directTerm = directTerm && (mcG1 || signalMC1->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(labelG1, signalMC1->GetGrandMotherSource(1));
    }

    Int_t labelG2 = -1;
    if(signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelG2 = GetMothersLabel(labelM2);
      if(mcM2 && labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
      directTerm = directTerm && (mcG2 || signalMC1->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC1->GetGrandMotherPDG(2),signalMC1->GetGrandMotherPDGexclude(2),signalMC1->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(labelG2, signalMC1->GetGrandMotherSource(2));
    }

    Int_t labelG3 = -1;
    if(signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelG3 = GetMothersLabel(labelM3);
      if(mcM3 && labelG3>-1) mcG3 = GetMCTrackFromMCEvent(labelG3);
      directTerm = directTerm && (mcG3 || signalMC2->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(labelG3, signalMC2->GetGrandMotherSource(1));
    }

    Int_t labelG4 = -1;
    if(signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelG4 = GetMothersLabel(labelM4);
      if(mcM4 && labelG4>-1) mcG4 = GetMCTrackFromMCEvent(labelG4);
      directTerm = directTerm && (mcG4 || signalMC2->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG4 ? mcG4->PdgCode() : 0),signalMC2->GetGrandMotherPDG(2),signalMC2->GetGrandMotherPDGexclude(2),signalMC2->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(labelG4, signalMC2->GetGrandMotherSource(2));
    }

    // Cross term
    Bool_t crossTerm = kTRUE;
    // daughters
    crossTerm = crossTerm && mcD2 && ComparePDG(d2Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD2, signalMC1->GetLegSource(1));
    crossTerm = crossTerm && mcD1 && ComparePDG(d1Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD1, signalMC1->GetLegSource(2));
    crossTerm = crossTerm && mcD4 && ComparePDG(d4Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD4, signalMC2->GetLegSource(1));
    crossTerm = crossTerm && mcD3 && ComparePDG(d3Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD3, signalMC2->GetLegSource(2));

    // mothers
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM2 && labelD2>-1) {
        labelM2 = GetMothersLabel(labelD2);
        if(labelM2>-1) mcM2 = GetMCTrackFromMCEvent(labelM2);
      }
      crossTerm = crossTerm && (mcM2 || signalMC1->GetMotherPDGexclude(1)) && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1), signalMC1->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM2, signalMC1->GetMotherSource(1)) && CheckRadiativeDecision(labelM2,signalMC1);
    }

    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM1 && labelD1>-1) {
        labelM1 = GetMothersLabel(labelD1);
        if(labelM1>-1) mcM1 = GetMCTrackFromMCEvent(labelM1);
      }
      crossTerm = crossTerm && (mcM1 || signalMC1->GetMotherPDGexclude(2)) && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(2), signalMC1->GetMotherPDGexclude(2), signalMC1->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM1, signalMC1->GetMotherSource(2)) && CheckRadiativeDecision(labelM1,signalMC1);
    }

    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM4 && labelD4>-1) {
        labelM4 = GetMothersLabel(labelD4);
        if(labelM4>-1) mcM4 = GetMCTrackFromMCEvent(labelM4);
      }
      crossTerm = crossTerm && (mcM4 || signalMC2->GetMotherPDGexclude(1)) && ComparePDG((mcM4 ? mcM4->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1), signalMC2->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM4, signalMC2->GetMotherSource(1)) && CheckRadiativeDecision(labelM4,signalMC2);
    }

    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM3 && labelD3>-1) {
        labelM3 = GetMothersLabel(labelD3);
        if(labelM3>-1) mcM1 = GetMCTrackFromMCEvent(labelM3);
      }
      crossTerm = crossTerm && (mcM3 || signalMC2->GetMotherPDGexclude(2)) && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(2), signalMC2->GetMotherPDGexclude(2), signalMC2->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM3, signalMC2->GetMotherSource(2)) && CheckRadiativeDecision(labelM3,signalMC2);
    }

    // grand-mothers
    if(signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG2 && mcM2) {
        labelG2 = GetMothersLabel(labelM2);
        if(labelG2>-1) mcG2 = GetMCTrackFromMCEvent(labelG2);
      }
      crossTerm = crossTerm && (mcG2 || signalMC1->GetGrandMotherPDGexclude(1)) && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG2, signalMC1->GetGrandMotherSource(1));
    }

    if(signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG1 && mcM1) {
        labelG1 = GetMothersLabel(labelM1);
        if(labelG1>-1) mcG1 = GetMCTrackFromMCEvent(labelG1);
      }
      crossTerm = crossTerm && (mcG1 || signalMC1->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(2), signalMC1->GetGrandMotherPDGexclude(2), signalMC1->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG1, signalMC1->GetGrandMotherSource(2));
    }

    if(signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG4 && mcM4) {
        labelG4 = GetMothersLabel(labelM2);
        if(labelG4>-1) mcG4 = GetMCTrackFromMCEvent(labelG4);
      }
      crossTerm = crossTerm && (mcG4 || signalMC2->GetGrandMotherPDGexclude(1)) && ComparePDG((mcG4 ? mcG4->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG4, signalMC2->GetGrandMotherSource(1));
    }

    if(signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG3 && mcM3) {
        labelG3 = GetMothersLabel(labelM3);
        if(labelG3>-1) mcG3 = GetMCTrackFromMCEvent(labelG3);
      }
      crossTerm = crossTerm && (mcG3 || signalMC2->GetGrandMotherPDGexclude(2)) && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(2), signalMC2->GetGrandMotherPDGexclude(2), signalMC2->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG3, signalMC2->GetGrandMotherSource(2));
    }

    Bool_t motherRelation1 = kTRUE;
    Bool_t motherRelation2 = kTRUE;
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      motherRelation1 = motherRelation1 && HaveSameMother(pair1);
    }
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      motherRelation1 = motherRelation1 && !HaveSameMother(pair1);
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      motherRelation2 = motherRelation2 && HaveSameMother(pair2);
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      motherRelation2 = motherRelation2 && !HaveSameMother(pair2);
    }

    Bool_t grandMotherRelation1 = kTRUE;
    Bool_t grandMotherRelation2 = kTRUE;
    if(signalMC1->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      grandMotherRelation1 = grandMotherRelation1 && HaveSameGrandMother(pair1);
    }
    if(signalMC1->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      grandMotherRelation1 = grandMotherRelation1 && !HaveSameGrandMother(pair1);
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      grandMotherRelation2 = grandMotherRelation2 && HaveSameGrandMother(pair2);
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      grandMotherRelation2 = grandMotherRelation2 && !HaveSameGrandMother(pair2);
    }

    // check geant process if set
    Bool_t processGEANT1 = kTRUE;
    Bool_t processGEANT2 = kTRUE;
    if(signalMC1->GetCheckGEANTProcess()) {
      if(!CheckGEANTProcess(labelD1,signalMC1->GetGEANTProcess()) &&
         !CheckGEANTProcess(labelD2,signalMC1->GetGEANTProcess())   ) processGEANT1= kFALSE;
    }
    if(signalMC2->GetCheckGEANTProcess()) {
      if(!CheckGEANTProcess(labelD3,signalMC2->GetGEANTProcess()) &&
         !CheckGEANTProcess(labelD4,signalMC2->GetGEANTProcess())   ) processGEANT2= kFALSE;
    }

    // check particle stack for pdg code
    Bool_t pdgInStack1 = kTRUE;
    Bool_t pdgInStack2 = kTRUE;
    if (signalMC1->GetCheckStackForPDG()) {
      pdgInStack1 = CheckStackParticle(labelM1,signalMC1->GetStackPDG()) && CheckStackParticle(labelM2,signalMC1->GetStackPDG());
    }
    if (signalMC2->GetCheckStackForPDG()) {
      pdgInStack2 = CheckStackParticle(labelM3,signalMC2->GetStackPDG()) && CheckStackParticle(labelM4,signalMC2->GetStackPDG());
    }
    // check is first mothers are closed (correlated charm/beauty)
    Bool_t correlated1 = kTRUE;
    Bool_t correlated2 = kTRUE;
    if (signalMC1->GetCheckCorrelatedHF()){
      Int_t labelfirstmother1 = GetFirstMothersLabelInChain(labelD1);
      Int_t labelminfirstmother1 = GetMinLabelParticlePrimaryAround(labelfirstmother1);
      Int_t labelmaxfirstmother1 = GetMaxLabelParticlePrimaryAround(labelfirstmother1);
      Int_t labelfirstmother2 = GetFirstMothersLabelInChain(labelD2);
      Int_t labelminfirstmother2 = GetMinLabelParticlePrimaryAround(labelfirstmother2);
      Int_t labelmaxfirstmother2 = GetMaxLabelParticlePrimaryAround(labelfirstmother2);
      if(labelfirstmother1==-1 || labelminfirstmother1==-1 || labelmaxfirstmother1==-1) correlated1 = kFALSE; // security
      if(labelfirstmother2==-1 || labelminfirstmother2==-1 || labelmaxfirstmother2==-1) correlated1 = kFALSE; // security
      if(labelfirstmother1<labelminfirstmother2 || labelfirstmother1>labelmaxfirstmother2) correlated1 = kFALSE;
      if(labelfirstmother2<labelminfirstmother1 || labelfirstmother2>labelmaxfirstmother1) correlated1 = kFALSE;
    }
    if (signalMC2->GetCheckCorrelatedHF()){
      Int_t labelfirstmother1 = GetFirstMothersLabelInChain(labelD3);
      Int_t labelminfirstmother1 = GetMinLabelParticlePrimaryAround(labelfirstmother1);
      Int_t labelmaxfirstmother1 = GetMaxLabelParticlePrimaryAround(labelfirstmother1);
      Int_t labelfirstmother2 = GetFirstMothersLabelInChain(labelD4);
      Int_t labelminfirstmother2 = GetMinLabelParticlePrimaryAround(labelfirstmother2);
      Int_t labelmaxfirstmother2 = GetMaxLabelParticlePrimaryAround(labelfirstmother2);
      if(labelfirstmother1==-1 || labelminfirstmother1==-1 || labelmaxfirstmother1==-1) correlated2 = kFALSE; // security
      if(labelfirstmother2==-1 || labelminfirstmother2==-1 || labelmaxfirstmother2==-1) correlated2 = kFALSE; // security
      if(labelfirstmother1<labelminfirstmother2 || labelfirstmother1>labelmaxfirstmother2) correlated2 = kFALSE;
      if(labelfirstmother2<labelminfirstmother1 || labelfirstmother2>labelmaxfirstmother1) correlated2 = kFALSE;
    }

    // check if a mother is also a grandmother
    Bool_t motherIsGrandmother1 = kTRUE;
    Bool_t motherIsGrandmother2 = kTRUE;
    if (signalMC1->GetCheckMotherGrandmotherRelation()) {
      motherIsGrandmother1 = kFALSE;
      motherIsGrandmother1 = MotherIsGrandmother(labelM1,labelM2,labelG1,labelG2, signalMC1->GetMotherIsGrandmother());
    }
    if (signalMC2->GetCheckMotherGrandmotherRelation()) {
      motherIsGrandmother2 = kFALSE;
      motherIsGrandmother2 = MotherIsGrandmother(labelM3,labelM4,labelG3,labelG4, signalMC2->GetMotherIsGrandmother());
    }

    // check if the mother from one pair is the grandmother of the other pair
    Bool_t motherIsGrandmotherfromDiffPair1 = kTRUE;
    Bool_t motherIsGrandmotherfromDiffPair2 = kTRUE;
    if (signalMC1->GetCheckMotherGrandmotherDiffPairRelation()) {
      labelM1 = GetMothersLabel(labelD1);
      labelM3 = GetMothersLabel(labelD3);
      labelG1 = GetMothersLabel(labelM1);
      labelG3 = GetMothersLabel(labelM3);
      motherIsGrandmotherfromDiffPair1 = kFALSE;
      motherIsGrandmotherfromDiffPair1 = MotherIsGrandmother(labelM1,labelM3,labelG1,labelG3, signalMC1->GetMotherIsGrandmotherDiffPair());
    }
    else if (signalMC2->GetCheckMotherGrandmotherDiffPairRelation()) {
      labelM1 = GetMothersLabel(labelD1);
      labelM3 = GetMothersLabel(labelD3);
      labelG1 = GetMothersLabel(labelM1);
      labelG3 = GetMothersLabel(labelM3);
      motherIsGrandmotherfromDiffPair2 = kFALSE;
      motherIsGrandmotherfromDiffPair2 = MotherIsGrandmother(labelM1,labelM3,labelG1,labelG3, signalMC2->GetMotherIsGrandmotherDiffPair());
    }

    return ((directTerm || crossTerm) && motherRelation1 && motherRelation2 && grandMotherRelation1 && grandMotherRelation2 && processGEANT1 && processGEANT2 && motherIsGrandmother1 && motherIsGrandmother2 && motherIsGrandmotherfromDiffPair1 && motherIsGrandmotherfromDiffPair2 && pdgInStack1 && pdgInStack2 && correlated1 && correlated2);

}



//________________________________________________________________________________
// Define IsMCTruth for 4 particles
// This IsMCTruth function is checking the correlation between particle 1&2 with the first MCSignal and the correlation between particle 3&4 with the second MCSignal
Bool_t AliDielectronMC::IsMCTruth(AliVParticle* mcD1, AliVParticle* mcD2, AliVParticle* mcD3, AliVParticle* mcD4, const AliDielectronSignalMC* signalMC1, const AliDielectronSignalMC* signalMC2) const {
  //
  // Check if the pair corresponds to the MC truth in signalMC
  //

  if ((mcD1->IsA()) != (mcD2->IsA()) && (mcD3->IsA()) != (mcD4->IsA()) ) {
    AliError("AliDielectron::IsMCTruth(): Not same particle types");
  }

  if (mcD1->IsA() == AliMCParticle::Class())     {
    mcD1 = static_cast<AliMCParticle*>(mcD1);
    mcD2 = static_cast<AliMCParticle*>(mcD2);
    mcD3 = static_cast<AliMCParticle*>(mcD3);
    mcD4 = static_cast<AliMCParticle*>(mcD4);

    // legs (daughters)
    Int_t labelD1 = (mcD1 ? TMath::Abs(mcD1->GetLabel()) : -1);
    Int_t labelD2 = (mcD2 ? TMath::Abs(mcD2->GetLabel()) : -1);
    Int_t labelD3 = (mcD3 ? TMath::Abs(mcD3->GetLabel()) : -1);
    Int_t labelD4 = (mcD4 ? TMath::Abs(mcD4->GetLabel()) : -1);
    Int_t d1Pdg = mcD1->PdgCode();
    Int_t d2Pdg = mcD2->PdgCode();
    Int_t d3Pdg = mcD3->PdgCode();
    Int_t d4Pdg = mcD4->PdgCode();


    // mothers
    AliMCParticle* mcM1 = 0x0;
    AliMCParticle* mcM2 = 0x0;
    AliMCParticle* mcM3 = 0x0;
    AliMCParticle* mcM4 = 0x0;

    // grand-mothers
    AliMCParticle* mcG1 = 0x0;
    AliMCParticle* mcG2 = 0x0;
    AliMCParticle* mcG3 = 0x0;
    AliMCParticle* mcG4 = 0x0;

    // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
    Bool_t directTerm = kTRUE;
    // daughters
    directTerm = directTerm && mcD1 && ComparePDG(d1Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD1, signalMC1->GetLegSource(1));
    directTerm = directTerm && mcD2 && ComparePDG(d2Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD2, signalMC1->GetLegSource(2));
    directTerm = directTerm && mcD3 && ComparePDG(d3Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD3, signalMC2->GetLegSource(1));
    directTerm = directTerm && mcD4 && ComparePDG(d4Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD4, signalMC2->GetLegSource(2));

    // mothers
    Int_t labelM1 = -1;
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM1 = mcD1->GetMother();
      if(labelD1>-1 && labelM1>-1) mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      directTerm = directTerm
                   && (mcM1 || signalMC1->GetMotherPDGexclude(1))
                   && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1),signalMC1->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(labelM1, signalMC1->GetMotherSource(1));
    }

    Int_t labelM2 = -1;
    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM2 = mcD2->GetMother();
      if(labelD2>-1 && labelM2>-1) mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      directTerm = directTerm
                   && (mcM2 || signalMC1->GetMotherPDGexclude(2))
                   && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC1->GetMotherPDG(2),signalMC1->GetMotherPDGexclude(2),signalMC1->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(labelM2, signalMC1->GetMotherSource(2));
    }

    Int_t labelM3 = -1;
    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM3 = mcD3->GetMother();
      if(labelD1>-1 && labelM1>-1) mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      directTerm = directTerm
                   && (mcM3 || signalMC2->GetMotherPDGexclude(1))
                   && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1),signalMC2->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(labelM3, signalMC2->GetMotherSource(1));
    }

    Int_t labelM4 = -1;
    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM4 = mcD4->GetMother();
      if(labelD4>-1 && labelM4>-1) mcM4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      directTerm = directTerm
                   && (mcM4 || signalMC2->GetMotherPDGexclude(2))
                   && ComparePDG((mcM4 ? mcM4->PdgCode() : 0),signalMC2->GetMotherPDG(2),signalMC2->GetMotherPDGexclude(2),signalMC2->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(labelM4, signalMC2->GetMotherSource(2));
    }

    // grand-mothers
    Int_t labelG1 = -1;
    if((signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM1 != 0x0) {
      labelG1 = mcM1->GetMother();
      if(mcM1 && labelG1>-1) mcG1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      directTerm = directTerm
                   && (mcG1 || signalMC1->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(labelG1, signalMC1->GetGrandMotherSource(1));
    }

    Int_t labelG2 = -1;
    if((signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM2 != 0x0) {
      labelG2 = mcM2->GetMother();
      if(mcM2 && labelG2>-1) mcG2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      directTerm = directTerm
                   && (mcG2 || signalMC1->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC1->GetGrandMotherPDG(2),signalMC1->GetGrandMotherPDGexclude(2),signalMC1->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(labelG2, signalMC1->GetGrandMotherSource(2));
    }

    Int_t labelG3 = -1;
    if((signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM3 != 0x0) {
      labelG3 = mcM3->GetMother();
      if(mcM3 && labelG3>-1) mcG3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG3));
      directTerm = directTerm
                   && (mcG3 || signalMC2->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(labelG3, signalMC2->GetGrandMotherSource(1));
    }

    Int_t labelG4 = -1;
    if((signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM4 != 0x0) {
      labelG4 = mcM4->GetMother();
      if(mcM4 && labelG4>-1) mcG4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG4));
      directTerm = directTerm
                   && (mcG4 || signalMC2->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG4 ? mcG4->PdgCode() : 0),signalMC2->GetGrandMotherPDG(2),signalMC2->GetGrandMotherPDGexclude(2),signalMC2->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(labelG4, signalMC2->GetGrandMotherSource(2));
    }

    // Cross term
    Bool_t crossTerm = kTRUE;
    // daughters
    crossTerm = crossTerm && mcD2 && ComparePDG(d2Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD2, signalMC1->GetLegSource(1));
    crossTerm = crossTerm && mcD1 && ComparePDG(d1Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD1, signalMC1->GetLegSource(2));
    crossTerm = crossTerm && mcD4 && ComparePDG(d4Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD4, signalMC2->GetLegSource(1));
    crossTerm = crossTerm && mcD3 && ComparePDG(d3Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD3, signalMC2->GetLegSource(2));

    // mothers
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM2 && labelD2>-1) {
        labelM2 = mcD2->GetMother();
        if(labelM2>-1) mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      }
      crossTerm = crossTerm
                  && (mcM2 || signalMC1->GetMotherPDGexclude(1))
                  && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1), signalMC1->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM2, signalMC1->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM2,signalMC)*/;
    }

    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM1 && labelD1>-1) {
        labelM1 = mcD1->GetMother();
        if(labelM1>-1) mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      }
      crossTerm = crossTerm
                  && (mcM1 || signalMC1->GetMotherPDGexclude(2))
                  && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(2), signalMC1->GetMotherPDGexclude(2), signalMC1->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM1, signalMC1->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM1,signalMC)*/;
    }

    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM4 && labelD4>-1) {
        labelM4 = mcD4->GetMother();
        if(labelM4>-1) mcM4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      }
      crossTerm = crossTerm
                  && (mcM4 || signalMC2->GetMotherPDGexclude(1))
                  && ComparePDG((mcM4 ? mcM4->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1), signalMC2->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM4, signalMC2->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM2,signalMC)*/;
    }

    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM3 && labelD3>-1) {
        labelM3 = mcD3->GetMother();
        if(labelM3>-1) mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      }
      crossTerm = crossTerm
                  && (mcM3 || signalMC2->GetMotherPDGexclude(2))
                  && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(2), signalMC2->GetMotherPDGexclude(2), signalMC2->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM3, signalMC2->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM1,signalMC)*/;
    }

    // grand-mothers
    if(signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG2 && mcM2) {
        labelG2 = mcM2->GetMother();
        if(labelG2>-1) mcG2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      }
      crossTerm = crossTerm
                  && (mcG2 || signalMC1->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG2, signalMC1->GetGrandMotherSource(1));
    }

    if(signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG1 && mcM1) {
        labelG1 = mcM1->GetMother();
        if(labelG1>-1) mcG1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      }
      crossTerm = crossTerm
                  && (mcG1 || signalMC1->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(2), signalMC1->GetGrandMotherPDGexclude(2), signalMC1->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG1, signalMC1->GetGrandMotherSource(2));
    }

    if(signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG4 && mcM4) {
        labelG4 = mcM4->GetMother();
        if(labelG4>-1) mcG4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG4));
      }
      crossTerm = crossTerm
                  && (mcG4 || signalMC2->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG4 ? mcG4->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG4, signalMC2->GetGrandMotherSource(1));
    }

    if(signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG3 && mcM3) {
        labelG3 = mcM3->GetMother();
        if(labelG3>-1) mcG3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG3));
      }
      crossTerm = crossTerm
                  && (mcG3 || signalMC2->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(2), signalMC2->GetGrandMotherPDGexclude(2), signalMC2->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG3, signalMC2->GetGrandMotherSource(2));
    }
    // checking mother relation between particle 1&2 as well as 3&4
    Bool_t motherRelation = kTRUE;
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      motherRelation = motherRelation && labelM1 != -1 && labelM2 != -1 && labelM1 == labelM2;
    }
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      motherRelation = motherRelation && (labelM1 == -1 || labelM2 == -1 || labelM1 != labelM2);
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM3 = mcD3->GetMother();
      labelM4 = mcD4->GetMother();
      motherRelation = motherRelation && labelM3 != -1 && labelM4 != -1 && labelM3 == labelM4;
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM3 = mcD3->GetMother();
      labelM4 = mcD4->GetMother();
      motherRelation = motherRelation && (labelM3 == -1 || labelM4 == -1 || labelM3 != labelM4);
    }

    // check if grand mother relation between particle 1&2 as well as 3&4
    Bool_t grandMotherRelation = kTRUE;
    if(signalMC1->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      labelG1 = mcM1->GetMother();
      labelG2 = mcM2->GetMother();
      grandMotherRelation = grandMotherRelation && labelG1 != -1 && labelG2 != -1 && labelG1 == labelG2;
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      labelG1 = mcM1->GetMother();
      labelG2 = mcM2->GetMother();
      grandMotherRelation = grandMotherRelation && (labelG1 == -1 || labelG2 == -1 || labelG1 != labelG2);
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM3 = mcD3->GetMother();
      labelM4 = mcD4->GetMother();
      mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      mcM4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      labelG3 = mcM3->GetMother();
      labelG4 = mcM4->GetMother();
      grandMotherRelation = grandMotherRelation && labelG3 != -1 && labelG4 != -1 && labelG3 == labelG4;
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM3 = mcD3->GetMother();
      labelM4 = mcD4->GetMother();
      mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      mcM4 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      labelG3 = mcM3->GetMother();
      labelG4 = mcM4->GetMother();
      grandMotherRelation = grandMotherRelation && (labelG3 == -1 || labelG4 == -1 || labelG3 != labelG4);
    }

    // check geant process if set
    Bool_t processGEANT = kTRUE;
    // Not implemented

    // check particle stack for pdg code
    Bool_t pdgInStack = kTRUE;
    // Not implemented

    // check if correlated for HF
    Bool_t correlated = kTRUE;
    // Not implemented

    // check if a mother is also a grandmother
    Bool_t motherIsGrandmother = kTRUE;
    // Not implemented

    // check if the mother from one part 1 is the grandmother of particle 3 by assuming that particle 1&2 have the same mother as well as particle 3&4 have same mother
    Bool_t motherIsGrandmotherfromDiffPair = kTRUE;
    if ((signalMC1->GetCheckMotherGrandmotherDiffPairRelation() && signalMC1->GetMotherIsGrandmotherDiffPair() == kTRUE) || (signalMC2->GetCheckMotherGrandmotherDiffPairRelation() && signalMC2->GetMotherIsGrandmotherDiffPair() == kTRUE)) {
      labelM1 = mcD1->GetMother();
      labelM3 = mcD3->GetMother();
      mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      labelG1 = mcM1->GetMother();
      labelG3 = mcM3->GetMother();
      motherIsGrandmotherfromDiffPair = motherIsGrandmotherfromDiffPair && labelM1 != -1 && labelM3 != -1 && labelG1 != -1 && labelG3 != -1 && (labelM1 == labelG3 || labelM3 == labelG1);
    }
    if ((signalMC1->GetCheckMotherGrandmotherDiffPairRelation() && signalMC1->GetMotherIsGrandmotherDiffPair() == kFALSE) || (signalMC2->GetCheckMotherGrandmotherDiffPairRelation() && signalMC2->GetMotherIsGrandmotherDiffPair() == kFALSE)) {
      labelM1 = mcD1->GetMother();
      labelM3 = mcD3->GetMother();
      mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM3 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      labelG1 = mcM1->GetMother();
      labelG3 = mcM3->GetMother();
      motherIsGrandmotherfromDiffPair = motherIsGrandmotherfromDiffPair && (labelM1 != -1 || labelM3 != -1 || labelG1 != -1 || labelG3 != -1 || (labelM1 != labelG3 || labelM3 != labelG1));
    }

    return ((directTerm || crossTerm) && motherRelation && grandMotherRelation && processGEANT && motherIsGrandmother && pdgInStack && correlated && motherIsGrandmotherfromDiffPair);

  }
  else if (mcD1->IsA() == AliAODMCParticle::Class()){
    AliAODMCParticle* mcD1AOD = static_cast<AliAODMCParticle*>(mcD1);
    AliAODMCParticle* mcD2AOD = static_cast<AliAODMCParticle*>(mcD2);
    AliAODMCParticle* mcD3AOD = static_cast<AliAODMCParticle*>(mcD3);
    AliAODMCParticle* mcD4AOD = static_cast<AliAODMCParticle*>(mcD4);

    // legs (daughters)
    Int_t labelD1 = (mcD1AOD ? TMath::Abs(mcD1AOD->GetLabel()) : -1);
    Int_t labelD2 = (mcD2AOD ? TMath::Abs(mcD2AOD->GetLabel()) : -1);
    Int_t labelD3 = (mcD3AOD ? TMath::Abs(mcD3AOD->GetLabel()) : -1);
    Int_t labelD4 = (mcD4AOD ? TMath::Abs(mcD4AOD->GetLabel()) : -1);
    Int_t d1Pdg = mcD1AOD->PdgCode();
    Int_t d2Pdg = mcD2AOD->PdgCode();
    Int_t d3Pdg = mcD3AOD->PdgCode();
    Int_t d4Pdg = mcD4AOD->PdgCode();

    // mothers
    AliAODMCParticle* mcM1 = 0x0;
    AliAODMCParticle* mcM2 = 0x0;
    AliAODMCParticle* mcM3 = 0x0;
    AliAODMCParticle* mcM4 = 0x0;

    // grand-mothers
    AliAODMCParticle* mcG1 = 0x0;
    AliAODMCParticle* mcG2 = 0x0;
    AliAODMCParticle* mcG3 = 0x0;
    AliAODMCParticle* mcG4 = 0x0;

    // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
    Bool_t directTerm = kTRUE;
    // daughters
    directTerm = directTerm && mcD1AOD && ComparePDG(d1Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD1AOD, signalMC1->GetLegSource(1));
    directTerm = directTerm && mcD2AOD && ComparePDG(d2Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD2AOD, signalMC1->GetLegSource(2));
    directTerm = directTerm && mcD3AOD && ComparePDG(d3Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD3AOD, signalMC2->GetLegSource(1));
    directTerm = directTerm && mcD4AOD && ComparePDG(d4Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD4AOD, signalMC2->GetLegSource(2));

    // mothers
    Int_t labelM1 = -1;
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM1 = mcD1AOD->GetMother();
      if(labelD1>-1 && labelM1>-1) mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      directTerm = directTerm
                   && (mcM1 || signalMC1->GetMotherPDGexclude(1))
                   && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1),signalMC1->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(mcM1, signalMC1->GetMotherSource(1));
    }

    Int_t labelM2 = -1;
    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM2 = mcD2AOD->GetMother();
      if(labelD2>-1 && labelM2>-1) mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      directTerm = directTerm
                   && (mcM2 || signalMC1->GetMotherPDGexclude(2))
                   && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC1->GetMotherPDG(2),signalMC1->GetMotherPDGexclude(2),signalMC1->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(mcM2, signalMC1->GetMotherSource(2));
    }

    Int_t labelM3 = -1;
    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM3 = mcD3AOD->GetMother();
      if(labelD3>-1 && labelM3>-1) mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      directTerm = directTerm
                   && (mcM3 || signalMC2->GetMotherPDGexclude(1))
                   && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1),signalMC2->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(mcM3, signalMC2->GetMotherSource(1));
    }

    Int_t labelM4 = -1;
    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM4 = mcD4AOD->GetMother();
      if(labelD4>-1 && labelM4>-1) mcM4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      directTerm = directTerm
                   && (mcM4 || signalMC2->GetMotherPDGexclude(2))
                   && ComparePDG((mcM4 ? mcM4->PdgCode() : 0),signalMC2->GetMotherPDG(2),signalMC2->GetMotherPDGexclude(2),signalMC2->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(mcM4, signalMC2->GetMotherSource(2));
    }

    // grand-mothers
    Int_t labelG1 = -1;
    if((signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM1 != 0x0) {
      labelG1 = mcM1->GetMother();
      if(mcM1 && labelG1>-1) mcG1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      directTerm = directTerm
                   && (mcG1 || signalMC1->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(mcG1, signalMC1->GetGrandMotherSource(1));
    }

    Int_t labelG2 = -1;
    if((signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM2 != 0x0) {
      labelG2 = mcM2->GetMother();
      if(mcM2 && labelG2>-1) mcG2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      directTerm = directTerm
                   && (mcG2 || signalMC1->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC1->GetGrandMotherPDG(2),signalMC1->GetGrandMotherPDGexclude(2),signalMC1->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(mcG2, signalMC1->GetGrandMotherSource(2));
    }

    Int_t labelG3 = -1;
    if((signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM3 != 0x0) {
      labelG3 = mcM3->GetMother();
      if(mcM3 && labelG3>-1) mcG3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG3));
      directTerm = directTerm
                   && (mcG3 || signalMC2->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(mcG3, signalMC2->GetGrandMotherSource(1));
    }

    Int_t labelG4 = -1;
    if((signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM4 != 0x0) {
      labelG4 = mcM4->GetMother();
      if(mcM4 && labelG4>-1) mcG4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG4));
      directTerm = directTerm
                   && (mcG4 || signalMC2->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG4 ? mcG4->PdgCode() : 0),signalMC2->GetGrandMotherPDG(2),signalMC2->GetGrandMotherPDGexclude(2),signalMC2->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(mcG4, signalMC2->GetGrandMotherSource(2));
    }


    // Cross term
    Bool_t crossTerm = kTRUE;
    // daughters
    crossTerm = crossTerm && mcD2AOD && ComparePDG(d2Pdg, signalMC1->GetLegPDG(1), signalMC1->GetLegPDGexclude(1), signalMC1->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD2AOD, signalMC1->GetLegSource(1));
    crossTerm = crossTerm && mcD1AOD && ComparePDG(d1Pdg, signalMC1->GetLegPDG(2), signalMC1->GetLegPDGexclude(2), signalMC1->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD1AOD, signalMC1->GetLegSource(2));
    crossTerm = crossTerm && mcD4AOD && ComparePDG(d4Pdg, signalMC2->GetLegPDG(1), signalMC2->GetLegPDGexclude(1), signalMC2->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD4AOD, signalMC2->GetLegSource(1));
    crossTerm = crossTerm && mcD3AOD && ComparePDG(d3Pdg, signalMC2->GetLegPDG(2), signalMC2->GetLegPDGexclude(2), signalMC2->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD3AOD, signalMC2->GetLegSource(2));

    // mothers
    if(signalMC1->GetMotherPDG(1)!=0 || signalMC1->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM2 && labelD2>-1) {
        labelM2 = mcD2AOD->GetMother();
        if(labelM2>-1) mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      }
      crossTerm = crossTerm
                  && (mcM2 || signalMC1->GetMotherPDGexclude(1))
                  && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC1->GetMotherPDG(1), signalMC1->GetMotherPDGexclude(1), signalMC1->GetCheckBothChargesMothers(1)) && CheckParticleSource(mcM2, signalMC1->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM2,signalMC1)*/;
    }

    if(signalMC1->GetMotherPDG(2)!=0 || signalMC1->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM1 && labelD1>-1) {
        labelM1 = mcD1AOD->GetMother();
        if(labelM1>-1) mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      }
      crossTerm = crossTerm
                  && (mcM1 || signalMC1->GetMotherPDGexclude(2))
                  && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC1->GetMotherPDG(2), signalMC1->GetMotherPDGexclude(2), signalMC1->GetCheckBothChargesMothers(2)) && CheckParticleSource(mcM1, signalMC1->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM1,signalMC1)*/;
    }

    if(signalMC2->GetMotherPDG(1)!=0 || signalMC2->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM4 && labelD4>-1) {
        labelM4 = mcD4AOD->GetMother();
        if(labelM4>-1) mcM4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      }
      crossTerm = crossTerm
                  && (mcM4 || signalMC2->GetMotherPDGexclude(1))
                  && ComparePDG((mcM4 ? mcM4->PdgCode() : 0), signalMC2->GetMotherPDG(1), signalMC2->GetMotherPDGexclude(1), signalMC2->GetCheckBothChargesMothers(1)) && CheckParticleSource(mcM4, signalMC2->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM4,signalMC2)*/;
    }

    if(signalMC2->GetMotherPDG(2)!=0 || signalMC2->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM3 && labelD3>-1) {
        labelM3 = mcD3AOD->GetMother();
        if(labelM3>-1) mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      }
      crossTerm = crossTerm
                  && (mcM3 || signalMC2->GetMotherPDGexclude(2))
                  && ComparePDG((mcM3 ? mcM3->PdgCode() : 0), signalMC2->GetMotherPDG(2), signalMC2->GetMotherPDGexclude(2), signalMC2->GetCheckBothChargesMothers(2)) && CheckParticleSource(mcM3, signalMC2->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM3,signalMC2)*/;
    }

    // grand-mothers
    if(signalMC1->GetGrandMotherPDG(1)!=0 || signalMC1->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG2 && mcM2) {
        labelG2 = mcM2->GetMother();
        if(labelG2>-1) mcG2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      }
      crossTerm = crossTerm
                  && (mcG2 || signalMC1->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC1->GetGrandMotherPDG(1), signalMC1->GetGrandMotherPDGexclude(1), signalMC1->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(mcG2, signalMC1->GetGrandMotherSource(1));
    }

    if(signalMC1->GetGrandMotherPDG(2)!=0 || signalMC1->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG1 && mcM1) {
        labelG1 = mcM1->GetMother();
        if(labelG1>-1) mcG1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      }
      crossTerm = crossTerm
                  && (mcG1 || signalMC1->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC1->GetGrandMotherPDG(2), signalMC1->GetGrandMotherPDGexclude(2), signalMC1->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(mcG1, signalMC1->GetGrandMotherSource(2));
    }

    if(signalMC2->GetGrandMotherPDG(1)!=0 || signalMC2->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG4 && mcM4) {
        labelG4 = mcM4->GetMother();
        if(labelG4>-1) mcG4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG4));
      }
      crossTerm = crossTerm
                  && (mcG4 || signalMC2->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG4 ? mcG4->PdgCode() : 0), signalMC2->GetGrandMotherPDG(1), signalMC2->GetGrandMotherPDGexclude(1), signalMC2->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(mcG4, signalMC2->GetGrandMotherSource(1));
    }

    if(signalMC2->GetGrandMotherPDG(2)!=0 || signalMC2->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG3 && mcM3) {
        labelG3 = mcM3->GetMother();
        if(labelG3>-1) mcG1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG3));
      }
      crossTerm = crossTerm
                  && (mcG3 || signalMC2->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG3 ? mcG3->PdgCode() : 0), signalMC2->GetGrandMotherPDG(2), signalMC2->GetGrandMotherPDGexclude(2), signalMC2->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(mcG3, signalMC2->GetGrandMotherSource(2));
    }


    // check if mother relation between particle 1&2 as well as 3&4
    Bool_t motherRelation = kTRUE;
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      motherRelation = motherRelation && labelM1 != -1 && labelM2 != -1 && labelM1 == labelM2;
    }
    if(signalMC1->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      motherRelation = motherRelation && (labelM1 == -1 || labelM2 == -1 || labelM1 != labelM2);
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM3 = mcD3AOD->GetMother();
      labelM4 = mcD4AOD->GetMother();
      motherRelation = motherRelation && labelM3 != -1 && labelM4 != -1 && labelM3 == labelM4;
    }
    if(signalMC2->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM3 = mcD3AOD->GetMother();
      labelM4 = mcD4AOD->GetMother();
      motherRelation = motherRelation && (labelM3 == -1 || labelM4 == -1 || labelM3 != labelM4);
    }

    // check if grand mother relation between particle 1&2 as well as 3&4
    Bool_t grandMotherRelation = kTRUE;
    if(signalMC1->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      labelG1 = mcM1->GetMother();
      labelG2 = mcM2->GetMother();
      grandMotherRelation = grandMotherRelation && labelG1 != -1 && labelG2 != -1 && labelG1 == labelG2;
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      labelG1 = mcM1->GetMother();
      labelG2 = mcM2->GetMother();
      grandMotherRelation = grandMotherRelation && (labelG1 == -1 || labelG2 == -1 || labelG1 != labelG2);
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM3 = mcD3AOD->GetMother();
      labelM4 = mcD4AOD->GetMother();
      mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      mcM4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      labelG3 = mcM3->GetMother();
      labelG4 = mcM4->GetMother();
      grandMotherRelation = grandMotherRelation && labelG3 != -1 && labelG4 != -1 && labelG3 == labelG4;
    }
    if(signalMC2->GetGrandMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM3 = mcD3AOD->GetMother();
      labelM4 = mcD4AOD->GetMother();
      mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      mcM4 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM4));
      labelG3 = mcM3->GetMother();
      labelG4 = mcM4->GetMother();
      grandMotherRelation = grandMotherRelation && (labelG3 == -1 || labelG4 == -1 || labelG3 != labelG4);
    }

    // check geant process if set
    Bool_t processGEANT = kTRUE;
    // Not implemented

    // check particle stack for pdg code
    Bool_t pdgInStack = kTRUE;
    // Not implemented

    // check if correlated
    Bool_t correlated = kTRUE;
    // Not implemented

    // check if a mother is also a grandmother
    Bool_t motherIsGrandmother = kTRUE;
    // Not implemented

    Bool_t motherIsGrandmotherfromDiffPair = kTRUE;
    if ((signalMC1->GetCheckMotherGrandmotherDiffPairRelation() && signalMC1->GetMotherIsGrandmotherDiffPair() == kTRUE) || (signalMC2->GetCheckMotherGrandmotherDiffPairRelation() && signalMC2->GetMotherIsGrandmotherDiffPair() == kTRUE)) {
      labelM1 = mcD1AOD->GetMother();
      labelM3 = mcD3AOD->GetMother();
      mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      labelG1 = mcM1->GetMother();
      labelG3 = mcM3->GetMother();
      motherIsGrandmotherfromDiffPair = (motherIsGrandmotherfromDiffPair && labelM1 != -1 && labelM3 != -1 && labelG1 != -1 && labelG3 != -1 && (labelM1 == labelG3 || labelM3 == labelG1));
    }
    if ((signalMC1->GetCheckMotherGrandmotherDiffPairRelation() && signalMC1->GetMotherIsGrandmotherDiffPair() == kFALSE) || (signalMC2->GetCheckMotherGrandmotherDiffPairRelation() && signalMC2->GetMotherIsGrandmotherDiffPair() == kFALSE)) {
      labelM1 = mcD1AOD->GetMother();
      labelM3 = mcD3AOD->GetMother();
      mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      mcM3 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM3));
      labelG1 = mcM1->GetMother();
      labelG3 = mcM3->GetMother();
      motherIsGrandmotherfromDiffPair = motherIsGrandmotherfromDiffPair && (labelM1 != -1 || labelM3 != -1 || labelG1 != -1 || labelG3 != -1 || (labelM1 != labelG3 || labelM3 != labelG1));
    }

    return ((directTerm || crossTerm) && motherRelation && grandMotherRelation && processGEANT && motherIsGrandmother && pdgInStack && correlated && motherIsGrandmotherfromDiffPair);

  }

  // If no MC particle is filled
  return false;
}


//________________________________________________________________________________
Bool_t AliDielectronMC::IsMCTruth(AliVParticle* mcD1, AliVParticle* mcD2, const AliDielectronSignalMC* signalMC) const {
  //
  // Check if the pair corresponds to the MC truth in signalMC
  //

  if (mcD1->IsA() != mcD2->IsA()) {
    AliError("AliDielectron::IsMCTruth(): Not same particle type");
  }

  if (mcD1->IsA() == AliMCParticle::Class())     {
    mcD1 = static_cast<AliMCParticle*>(mcD1);
    mcD2 = static_cast<AliMCParticle*>(mcD2);

    // legs (daughters)
    Int_t labelD1 = (mcD1 ? TMath::Abs(mcD1->GetLabel()) : -1);
    Int_t labelD2 = (mcD2 ? TMath::Abs(mcD2->GetLabel()) : -1);
    Int_t d1Pdg = mcD1->PdgCode();
    Int_t d2Pdg = mcD2->PdgCode();

    // mothers
    AliMCParticle* mcM1 = 0x0;
    AliMCParticle* mcM2 = 0x0;

    // grand-mothers
    AliMCParticle* mcG1 = 0x0;
    AliMCParticle* mcG2 = 0x0;

    // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
    Bool_t directTerm = kTRUE;
    // daughters
    directTerm = directTerm
                 && mcD1
                 && ComparePDG(d1Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD1, signalMC->GetLegSource(1));
    directTerm = directTerm
                 && mcD2
                 && ComparePDG(d2Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD2, signalMC->GetLegSource(2));

    // mothers
    Int_t labelM1 = -1;
    if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM1 = mcD1->GetMother();
      if(labelD1>-1 && labelM1>-1) mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      directTerm = directTerm
                   && (mcM1 || signalMC->GetMotherPDGexclude(1))
                   && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1),signalMC->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(labelM1, signalMC->GetMotherSource(1));
    }

    Int_t labelM2 = -1;
    if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM2 = mcD2->GetMother();
      if(labelD2>-1 && labelM2>-1) mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      directTerm = directTerm
                   && (mcM2 || signalMC->GetMotherPDGexclude(2))
                   && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC->GetMotherPDG(2),signalMC->GetMotherPDGexclude(2),signalMC->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(labelM2, signalMC->GetMotherSource(2));
    }

    // grand-mothers
    Int_t labelG1 = -1;
    if((signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM1 != 0x0) {
      labelG1 = mcM1->GetMother();
      if(mcM1 && labelG1>-1) mcG1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      directTerm = directTerm
                   && (mcG1 || signalMC->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(1));
    }

    Int_t labelG2 = -1;
    if((signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM2 != 0x0) {
      labelG2 = mcM2->GetMother();
      if(mcM2 && labelG2>-1) mcG2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      directTerm = directTerm
                   && (mcG2 || signalMC->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC->GetGrandMotherPDG(2),signalMC->GetGrandMotherPDGexclude(2),signalMC->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(2));
    }

    // Cross term
    Bool_t crossTerm = kTRUE;
    // daughters
    crossTerm = crossTerm && mcD2 && ComparePDG(d2Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(labelD2, signalMC->GetLegSource(1));
    crossTerm = crossTerm && mcD1 && ComparePDG(d1Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2)) && CheckParticleSource(labelD1, signalMC->GetLegSource(2));

    // mothers
    if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM2 && labelD2>-1) {
        labelM2 = mcD2->GetMother();
        if(labelM2>-1) mcM2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      }
      crossTerm = crossTerm
                  && (mcM2 || signalMC->GetMotherPDGexclude(1))
                  && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1), signalMC->GetCheckBothChargesMothers(1)) && CheckParticleSource(labelM2, signalMC->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM2,signalMC)*/;
    }

    if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM1 && labelD1>-1) {
        labelM1 = mcD1->GetMother();
        if(labelM1>-1) mcM1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      }
      crossTerm = crossTerm
                  && (mcM1 || signalMC->GetMotherPDGexclude(2))
                  && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(2), signalMC->GetMotherPDGexclude(2), signalMC->GetCheckBothChargesMothers(2)) && CheckParticleSource(labelM1, signalMC->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM1,signalMC)*/;
    }

    // grand-mothers
    if(signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG2 && mcM2) {
        labelG2 = mcM2->GetMother();
        if(labelG2>-1) mcG2 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      }
      crossTerm = crossTerm
                  && (mcG2 || signalMC->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(labelG2, signalMC->GetGrandMotherSource(1));
    }

    if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG1 && mcM1) {
        labelG1 = mcM1->GetMother();
        if(labelG1>-1) mcG1 = static_cast<AliMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      }
      crossTerm = crossTerm
                  && (mcG1 || signalMC->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(2), signalMC->GetGrandMotherPDGexclude(2), signalMC->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(labelG1, signalMC->GetGrandMotherSource(2));
    }

    Bool_t motherRelation = kTRUE;
    if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      motherRelation = motherRelation && labelM1 != -1 && labelM2 != -1 && labelM1 == labelM2;
    }

    if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1->GetMother();
      labelM2 = mcD2->GetMother();
      motherRelation = motherRelation && (labelM1 == -1 || labelM2 == -1 || labelM1 != labelM2);
    }

    // check geant process if set
    Bool_t processGEANT = kTRUE;
    // Not implemented

    // check particle stack for pdg code
    Bool_t pdgInStack = kTRUE;
    // Not implemented

    // check if correlated for HF
    Bool_t correlated = kTRUE;
    // Not implemented

    // check if a mother is also a grandmother
    Bool_t motherIsGrandmother = kTRUE;
    // Not implemented

    return ((directTerm || crossTerm) && motherRelation && processGEANT && motherIsGrandmother && pdgInStack && correlated);

  }
  else if (mcD1->IsA() == AliAODMCParticle::Class()){
    AliAODMCParticle* mcD1AOD = static_cast<AliAODMCParticle*>(mcD1);
    AliAODMCParticle* mcD2AOD = static_cast<AliAODMCParticle*>(mcD2);

    // legs (daughters)
    Int_t labelD1 = (mcD1AOD ? TMath::Abs(mcD1AOD->GetLabel()) : -1);
    Int_t labelD2 = (mcD2AOD ? TMath::Abs(mcD2AOD->GetLabel()) : -1);
    Int_t d1Pdg = mcD1AOD->PdgCode();
    Int_t d2Pdg = mcD2AOD->PdgCode();

    // mothers
    AliAODMCParticle* mcM1 = 0x0;
    AliAODMCParticle* mcM2 = 0x0;

    // grand-mothers
    AliAODMCParticle* mcG1 = 0x0;
    AliAODMCParticle* mcG2 = 0x0;

    // make direct(1-1 and 2-2) and cross(1-2 and 2-1) comparisons for the whole branch
    Bool_t directTerm = kTRUE;
    // daughters
    directTerm = directTerm
                 && mcD1AOD
                 && ComparePDG(d1Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD1AOD, signalMC->GetLegSource(1));
    directTerm = directTerm
                 && mcD2AOD
                 && ComparePDG(d2Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD2AOD, signalMC->GetLegSource(2));

    // mothers
    Int_t labelM1 = -1;
    if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      labelM1 = mcD1AOD->GetMother();
      if(labelD1>-1 && labelM1>-1) mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      directTerm = directTerm
                   && (mcM1 || signalMC->GetMotherPDGexclude(1))
                   && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1),signalMC->GetCheckBothChargesMothers(1))
                   && CheckParticleSource(mcM1, signalMC->GetMotherSource(1));
    }

    Int_t labelM2 = -1;
    if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      labelM2 = mcD2AOD->GetMother();
      if(labelD2>-1 && labelM2>-1) mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      directTerm = directTerm
                   && (mcM2 || signalMC->GetMotherPDGexclude(2))
                   && ComparePDG((mcM2 ? mcM2->PdgCode() : 0),signalMC->GetMotherPDG(2),signalMC->GetMotherPDGexclude(2),signalMC->GetCheckBothChargesMothers(2))
                   && CheckParticleSource(mcM2, signalMC->GetMotherSource(2));
    }

    // grand-mothers
    Int_t labelG1 = -1;
    if((signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) && mcM1 != 0x0) {
      labelG1 = mcM1->GetMother();
      if(mcM1 && labelG1>-1) mcG1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      directTerm = directTerm
                   && (mcG1 || signalMC->GetGrandMotherPDGexclude(1))
                   && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1))
                   && CheckParticleSource(mcG1, signalMC->GetGrandMotherSource(1));
    }

    Int_t labelG2 = -1;
    if((signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) && mcM2 != 0x0) {
      labelG2 = mcM2->GetMother();
      if(mcM2 && labelG2>-1) mcG2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      directTerm = directTerm
                   && (mcG2 || signalMC->GetGrandMotherPDGexclude(2))
                   && ComparePDG((mcG2 ? mcG2->PdgCode() : 0),signalMC->GetGrandMotherPDG(2),signalMC->GetGrandMotherPDGexclude(2),signalMC->GetCheckBothChargesGrandMothers(2))
                   && CheckParticleSource(mcG2, signalMC->GetGrandMotherSource(2));
    }

    // Cross term
    Bool_t crossTerm = kTRUE;
    // daughters
    crossTerm = crossTerm && mcD2AOD && ComparePDG(d2Pdg, signalMC->GetLegPDG(1), signalMC->GetLegPDGexclude(1), signalMC->GetCheckBothChargesLegs(1)) && CheckParticleSource(mcD2AOD, signalMC->GetLegSource(1));
    crossTerm = crossTerm && mcD1AOD && ComparePDG(d1Pdg, signalMC->GetLegPDG(2), signalMC->GetLegPDGexclude(2), signalMC->GetCheckBothChargesLegs(2)) && CheckParticleSource(mcD1AOD, signalMC->GetLegSource(2));

    // mothers
    if(signalMC->GetMotherPDG(1)!=0 || signalMC->GetMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM2 && labelD2>-1) {
        labelM2 = mcD2AOD->GetMother();
        if(labelM2>-1) mcM2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM2));
      }
      crossTerm = crossTerm
                  && (mcM2 || signalMC->GetMotherPDGexclude(1))
                  && ComparePDG((mcM2 ? mcM2->PdgCode() : 0), signalMC->GetMotherPDG(1), signalMC->GetMotherPDGexclude(1), signalMC->GetCheckBothChargesMothers(1)) && CheckParticleSource(mcM2, signalMC->GetMotherSource(1)) /*&& CheckRadiativeDecision(labelM2,signalMC)*/;
    }

    if(signalMC->GetMotherPDG(2)!=0 || signalMC->GetMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcM1 && labelD1>-1) {
        labelM1 = mcD1AOD->GetMother();
        if(labelM1>-1) mcM1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelM1));
      }
      crossTerm = crossTerm
                  && (mcM1 || signalMC->GetMotherPDGexclude(2))
                  && ComparePDG((mcM1 ? mcM1->PdgCode() : 0), signalMC->GetMotherPDG(2), signalMC->GetMotherPDGexclude(2), signalMC->GetCheckBothChargesMothers(2)) && CheckParticleSource(mcM1, signalMC->GetMotherSource(2)) /*&& CheckRadiativeDecision(labelM1,signalMC)*/;
    }

    // grand-mothers
    if(signalMC->GetGrandMotherPDG(1)!=0 || signalMC->GetGrandMotherSource(1)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG2 && mcM2) {
        labelG2 = mcM2->GetMother();
        if(labelG2>-1) mcG2 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG2));
      }
      crossTerm = crossTerm
                  && (mcG2 || signalMC->GetGrandMotherPDGexclude(1))
                  && ComparePDG((mcG2 ? mcG2->PdgCode() : 0), signalMC->GetGrandMotherPDG(1), signalMC->GetGrandMotherPDGexclude(1), signalMC->GetCheckBothChargesGrandMothers(1)) && CheckParticleSource(mcG2, signalMC->GetGrandMotherSource(1));
    }

    if(signalMC->GetGrandMotherPDG(2)!=0 || signalMC->GetGrandMotherSource(2)!=AliDielectronSignalMC::kDontCare) {
      if(!mcG1 && mcM1) {
        labelG1 = mcM1->GetMother();
        if(labelG1>-1) mcG1 = static_cast<AliAODMCParticle*>(GetMCTrackFromMCEvent(labelG1));
      }
      crossTerm = crossTerm
                  && (mcG1 || signalMC->GetGrandMotherPDGexclude(2))
                  && ComparePDG((mcG1 ? mcG1->PdgCode() : 0), signalMC->GetGrandMotherPDG(2), signalMC->GetGrandMotherPDGexclude(2), signalMC->GetCheckBothChargesGrandMothers(2)) && CheckParticleSource(mcG1, signalMC->GetGrandMotherSource(2));
    }

    Bool_t motherRelation = kTRUE;
    if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kSame) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      motherRelation = motherRelation && labelM1 != -1 && labelM2 != -1 && labelM1 == labelM2;
    }
    if(signalMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent) {
      labelM1 = mcD1AOD->GetMother();
      labelM2 = mcD2AOD->GetMother();
      motherRelation = motherRelation && (labelM1 == -1 || labelM2 == -1 || labelM1 != labelM2);
    }
    // check if grand mother relation
    Bool_t grandMotherRelation = kTRUE;
    // Not implemented

    // check geant process if set
    Bool_t processGEANT = kTRUE;
    // Not implemented

    // check particle stack for pdg code
    Bool_t pdgInStack = kTRUE;
    // Not implemented

    // check if correlated
    Bool_t correlated = kTRUE;
    // Not implemented

    // check if a mother is also a grandmother
    Bool_t motherIsGrandmother = kTRUE;
    // Not implemented

    return ((directTerm || crossTerm) && motherRelation && grandMotherRelation && processGEANT && motherIsGrandmother && pdgInStack && correlated);

  }

  // If no MC particle is filled
  return false;
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
  const Int_t nd=(mother->GetDaughterLast()-mother->GetDaughterFirst()+1);
  //Loop over daughters and get stop loop if PDG is found.
  for (Int_t i=0; i<nd; ++i) {
    //Go through the daughters, set result true and break if PDG code fits
    if (ComparePDG((GetMCTrackFromMCEvent(mother->GetDaughterFirst()+i)->PdgCode()),reqPDG,kFALSE,CheckBothChargesDaughter)) {
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

  AliVParticle *mcDaughter1=GetMCTrackFromMCEvent(TMath::Abs(daughter1->GetLabel()));
  AliVParticle *mcDaughter2=GetMCTrackFromMCEvent(TMath::Abs(daughter2->GetLabel()));
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

//____________________________________________________________
Bool_t AliDielectronMC::HaveSameGrandMother(const AliDielectronPair * pair) const
{
  //
  // Check whether two particles have the same mother
  //

  const AliVParticle * daughter1 = pair->GetFirstDaughterP();
  const AliVParticle * daughter2 = pair->GetSecondDaughterP();
  if (!daughter1 || !daughter2) return 0;

  AliVParticle *mcDaughter1=GetMCTrackFromMCEvent(TMath::Abs(daughter1->GetLabel()));
  AliVParticle *mcDaughter2=GetMCTrackFromMCEvent(TMath::Abs(daughter2->GetLabel()));
  if (!mcDaughter1 || !mcDaughter2) return 0;

  Int_t labelMother1=-1;
  Int_t labelMother2=-1;
  Int_t labelGrandMother1=-1;
  Int_t labelGrandMother2=-1;


  if (mcDaughter1->IsA()==AliMCParticle::Class()){
    labelMother1=(static_cast<AliMCParticle*>(mcDaughter1))->GetMother();
    labelMother2=(static_cast<AliMCParticle*>(mcDaughter2))->GetMother();
    AliVParticle *mcMother1=GetMCTrackFromMCEvent(labelMother1);
    AliVParticle *mcMother2=GetMCTrackFromMCEvent(labelMother2);
    labelGrandMother1=(static_cast<AliMCParticle*>(mcMother1))->GetMother();
    labelGrandMother2=(static_cast<AliMCParticle*>(mcMother2))->GetMother();
  }
  else if (mcDaughter1->IsA()==AliAODMCParticle::Class()) {
    labelMother1=(static_cast<AliAODMCParticle*>(mcDaughter1))->GetMother();
    labelMother2=(static_cast<AliAODMCParticle*>(mcDaughter2))->GetMother();
    AliVParticle *mcMother1=GetMCTrackFromMCEvent(labelMother1);
    AliVParticle *mcMother2=GetMCTrackFromMCEvent(labelMother2);
    labelGrandMother1=(static_cast<AliAODMCParticle*>(mcMother1))->GetMother();
    labelGrandMother2=(static_cast<AliAODMCParticle*>(mcMother2))->GetMother();
  }

  Bool_t sameGrandMother=(labelGrandMother1>-1)&&(labelGrandMother2>-1)&&(labelGrandMother1==labelGrandMother2);

  return sameGrandMother;
}

//________________________________________________________________
Int_t AliDielectronMC::IsJpsiPrimary(const AliDielectronPair * pair)
{
 // return: "0" for primary jpsi
 //         "1" for secondary jpsi (from beauty)
 //         "2" for background
 if(!HaveSameMother(pair)) return 2;
 AliVParticle *mcDaughter1=GetMCTrackFromMCEvent((pair->GetFirstDaughterP())->GetLabel());
 AliVParticle *mcDaughter2=GetMCTrackFromMCEvent((pair->GetSecondDaughterP())->GetLabel());
 if (!mcDaughter1 || !mcDaughter2) return 2;
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


//____________________________________________________________
Bool_t AliDielectronMC::LoadHFPairs()
{

  //
  // Look for all correlated c/cbar and b/bbar pairs
  // Attributing for each quark a HF Creation Process ID
  //

  Int_t quark[2][2]={{0}};
  Int_t quarktmp[2][2]={{0}};
  Int_t hadrontmp[2][2]={{0}};

  //Loop over the MC event to tag all correlated c/cbar and b/bbar quarks
  for(Int_t i=0;i<fMCEvent->GetNumberOfTracks();i++){
    AliMCParticle *cand = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(i));
    Int_t pdg_part = fMCEvent->GetTrack(i)->PdgCode();
    Int_t pdg_parta = TMath::Abs(pdg_part);

    if(pdg_parta==4 || pdg_parta==5){
      //The quark is requested to create a B or D hadron -> Must have a daughter
      if(!(cand->GetDaughterFirst()>-1)){
	//childless quark to be linked!
	if(pdg_part==4) quarktmp[0][0]=i;
	else if(pdg_part==-4) quarktmp[0][1]=i;
	else if(pdg_part==5) quarktmp[1][0]=i;
	else if(pdg_part==-5) quarktmp[1][1]=i;
	continue;
      }

      //Check if the quark is the one leading to the hadron
      Bool_t isHFprim(kTRUE);
      for(Int_t idau=cand->GetDaughterFirst();idau<=cand->GetDaughterLast();idau++){
	if(fMCEvent->GetTrack(idau)->PdgCode()==pdg_part){
	  isHFprim=kFALSE;
	  break;
	}
      }
      if(!isHFprim) continue;

      //quark is selected, saving its label
      if(pdg_part==4){
	if(quark[0][0]>0) return kFALSE;
	quark[0][0]=i;
      }
      else if(pdg_part==-4){
	if(quark[0][1]>0) return kFALSE;
	quark[0][1]=i;
      }
      else if(pdg_part==5){
	if(quark[1][0]>0) return kFALSE;
	quark[1][0]=i;
      }
      else if(pdg_part==-5){
	if(quark[1][1]>0) return kFALSE;
	quark[1][1]=i;
      }
    } //End looking at quarks

    //Look for Motherless hadrons or hadron's mother is u/d
    else if(pdg_parta==411 || pdg_parta==421 || pdg_parta==431 || pdg_parta==511 || pdg_parta==521 || pdg_parta==531 || pdg_parta==541 || pdg_parta==4122 || pdg_parta==5122){
      //Access the oldest ancestor that could be link to a corresponding quark
      AliMCParticle *mother, *daughter;
      Bool_t Osci(kFALSE);
      if(cand->GetMother()>-1){
	mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(cand->GetMother()));
	if(TMath::Abs(mother->PdgCode())<6)daughter=cand;
	else daughter=mother;

	if((pdg_parta== 411 || pdg_parta== 421 || pdg_parta== 431 || pdg_parta== 4122) && IsaBhadron(mother->PdgCode())) Osci=kTRUE;
	while (!(TMath::Abs(mother->PdgCode())<6 && TMath::Abs(mother->PdgCode())>0) && mother->GetMother()>-1){
	  daughter = mother;
	  mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(mother->GetMother()));
	  if((pdg_parta== 411 || pdg_parta== 421 || pdg_parta== 431 || pdg_parta== 4122) && IsaBhadron(mother->PdgCode())) Osci=kTRUE;
	}
      }
      else daughter = cand;
      mother=daughter;

      //Get Original mother in History
      Int_t pdg_moma(0);
      if(mother->GetMother()>-1) pdg_moma=TMath::Abs(fMCEvent->GetTrack(mother->GetMother())->PdgCode());
      if(TMath::Abs(pdg_moma)>6 && mother->GetMother()>-1) continue;

      if((pdg_part== 411 || pdg_part== 421 || pdg_part== 431 || pdg_part== 4122) && !Osci){
	if(hadrontmp[0][0]>0 && hadrontmp[0][0]!=mother->GetLabel()) return kFALSE;
	else hadrontmp[0][0]=mother->GetLabel();//c hadron type
      }
      else if((pdg_part== -411 || pdg_part== -421 || pdg_part== -431 || pdg_part== -4122) && !Osci){
	if(hadrontmp[0][1]>0 && hadrontmp[0][1]!=mother->GetLabel()) return kFALSE;
	else hadrontmp[0][1]=mother->GetLabel();//cbar hadron type
      }
      else if(pdg_part== -511 || pdg_part== -521 || pdg_part== -531 || pdg_part== -541 || pdg_part==5122 || ((pdg_part== 411 || pdg_part== 421 || pdg_part== 431 || pdg_part== 4122) && Osci)){
	if(hadrontmp[1][0]>0 && hadrontmp[1][0]!=mother->GetLabel()) return kFALSE;
	else hadrontmp[1][0]=mother->GetLabel();//b hadron type
      }
      else if(pdg_part== 511 || pdg_part== 521 || pdg_part== 531 || pdg_part== 541 || pdg_part==-5122 || ((pdg_part== -411 || pdg_part== -421 || pdg_part== -431 || pdg_part== -4122) && Osci)){
	if(hadrontmp[1][1]>0 && hadrontmp[1][1]!=mother->GetLabel()) return kFALSE;
	else hadrontmp[1][1]=mother->GetLabel();//bbar hadron type
      }
    }//End checking motherless hadrons
  }//End Loop on the Stack

  //Linking childless quarks with corresponding motherless hadrons
  if(quarktmp[0][0]>0 && hadrontmp[0][0]>0){
    quark[0][0]=hadrontmp[0][0];
  }
  if(quarktmp[0][1]>0 && hadrontmp[0][1]>0){
    quark[0][1]=hadrontmp[0][1];
  }
  if(quarktmp[1][0]>0 && hadrontmp[1][0]>0){
    quark[1][0]=hadrontmp[1][0];
  }
  if(quarktmp[1][1]>0 && hadrontmp[1][1]>0){
    quark[1][1]=hadrontmp[1][1];
  }

  //Pairing the quarks and attribute them a process number: c/cbar(1) and b/bar(2)
  if(quark[0][0]>0 && quark[0][1]>0){
    fhfproc.insert(std::pair<Int_t,Int_t>(quark[0][0],1));
    fhfproc.insert(std::pair<Int_t,Int_t>(quark[0][1],1));
  }
  if(quark[1][0]>0 && quark[1][1]>0){
    fhfproc.insert(std::pair<Int_t,Int_t>(quark[1][0],2));
    fhfproc.insert(std::pair<Int_t,Int_t>(quark[1][1],2));
  }

  return kTRUE;
}

//____________________________________________________________
Int_t AliDielectronMC::GetHFProcess(const Int_t label)
{
  //
  // return Heavy Flavour process number of the particle
  //
  if(!fCheckHF) return -1; //More than one pair of each -> Undeterminated (so far)
  if(fhfproc.size()==0) return 0; //No HF pairs in the event

  Int_t mother_label;
  AliMCParticle *part = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
  if(part->GetMother()==-1) return 0;

  // Looking back in the history if an ancestor is coming from an HF process
  AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(part->GetMother()));
  mother_label=mother->GetLabel();
  for(std::map<Int_t,Int_t>::iterator it = fhfproc.begin(); it != fhfproc.end(); ++it){
    if(mother_label==it->first) return it->second;
  }

  while (!(TMath::Abs(mother->PdgCode())==4 || TMath::Abs(mother->PdgCode())==5) && mother->GetMother()>-1){
    mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(mother->GetMother()));
    mother_label=mother->GetLabel();
    for(std::map<Int_t,Int_t>::iterator it = fhfproc.begin(); it != fhfproc.end(); ++it){
      if(mother_label==it->first) return it->second;
    }
  }

  return 0;
}


Int_t AliDielectronMC::IsaBhadron(Int_t pdg) const{
  Int_t Bhadronspdg[]={511,521,10511,10521,513,523,10513,20513,20523,515,525,531,10531,533,10533,20533,535,541,10541,543,10543,20543,545};
  Int_t size=sizeof(Bhadronspdg)/sizeof(*Bhadronspdg);
  for(Int_t i=0;i<size;i++){
    if(TMath::Abs(pdg==Bhadronspdg[i])) return kTRUE;
  }
  return kFALSE;
}
