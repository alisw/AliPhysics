/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//
// Signal cuts
// Checks whether a particle (reconstructed or MC) is coming from MC Signal
// For more information see implementation file
//
// Autor:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TClass.h>
#include <TMath.h>
#include <TParticle.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"

#include "AliHFEsignalCuts.h"
#include "AliHFEmcQA.h"

ClassImp(AliHFEsignalCuts)

//____________________________________________________________
AliHFEsignalCuts::AliHFEsignalCuts():
  AliAnalysisCuts(),
  fMC(NULL),
  fMCQA(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliHFEsignalCuts::AliHFEsignalCuts(const Char_t *name, const Char_t *title):
  AliAnalysisCuts(name, title),
  fMC(NULL),
  fMCQA(NULL)
{
  //
  // Default constructor
  //
  fMCQA = new AliHFEmcQA;
  if(fMCQA) fMCQA->Init();
}

//____________________________________________________________
AliHFEsignalCuts::AliHFEsignalCuts(const AliHFEsignalCuts &ref):
  AliAnalysisCuts(ref),
  fMC(ref.fMC),
  fMCQA(ref.fMCQA)
{
  //
  // Copy constructor
  //
}

//____________________________________________________________
AliHFEsignalCuts &AliHFEsignalCuts::operator=(const AliHFEsignalCuts &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    fMC = ref.fMC; 
    fMCQA = ref.fMCQA; 
  }
  return *this;
}

//____________________________________________________________
AliHFEsignalCuts::~AliHFEsignalCuts(){
  //
  // Destructor
  //
  if(fMCQA) delete fMCQA;
}

//____________________________________________________________
void AliHFEsignalCuts::SetMCEvent(AliMCEvent *mc){ 
  //
  // Set mc event
  //
  fMC = mc; 
  if(fMCQA) fMCQA->SetMCEvent(mc);
}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsSelected(TObject *o){
  //
  // Define signal as electron coming from charm or beauty
  // @TODO: Implement setter so that also either of them can be defined
  // as signal alone
  

  return IsCharmElectron(o) || IsBeautyElectron(o);
/*  
  //saving time?
  Int_t esources = GetElecSource(dynamic_cast<const AliVParticle *>(o));
  if(esources>0)printf("esources= %d\n",esources);
  if(esources == AliHFEmcQA::kDirectCharm || esources == AliHFEmcQA::kDirectBeauty || esources == AliHFEmcQA::kBeautyCharm)  // 1: direct D->e, 2: B->e 3: B->D->e
    return kTRUE;
  else
    return kFALSE;
*/

}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsCharmElectron(const TObject * const o) const {
  //
  // Check if mother is coming from Charm
  //
  if(!dynamic_cast<const AliVParticle *>(o)) return kFALSE;
  Int_t esources = GetElecSource(dynamic_cast<const AliVParticle *>(o));
  if(esources == AliHFEmcQA::kDirectCharm)  // 1: direct D->e
    return kTRUE;
  else
    return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsBeautyElectron(const TObject * const o) const {
  //
  // Check if mother is coming from Beauty
  //
  if(!dynamic_cast<const AliVParticle *>(o)) return kFALSE;
  Int_t esources = GetElecSource(dynamic_cast<const AliVParticle *>(o));
  if(esources == AliHFEmcQA::kDirectBeauty || esources == AliHFEmcQA::kBeautyCharm)  // 2: B->e 3: B->D->e
    return kTRUE;
  else
    return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsGammaElectron(const TObject * const o) const {
  //
  // Check for MC if the electron is coming from Gamma
  //
  if(!dynamic_cast<const AliVParticle *>(o)) return kFALSE;
  Int_t esources = GetElecSource(dynamic_cast<const AliVParticle *>(o));
  if(esources == AliHFEmcQA::kGamma)  // 4: conversion electrons
    return kTRUE;
  else
    return kFALSE;
}

/*
//____________________________________________________________
Bool_t AliHFEsignalCuts::IsCharmElectron(const TObject * const o) const {
  //
  // Check if mother is coming from Charm
  //
  if(TMath::Abs(GetTrackPDG(dynamic_cast<const AliVParticle *>(o))) != 11) return kFALSE;
  Int_t motherpdg = TMath::Abs(GetMotherPDG(dynamic_cast<const AliVParticle *>(o)));
  AliDebug(1, Form("Mother PDG %d\n", motherpdg));

  if((motherpdg % 1000) / 100 == 4) return kTRUE;    // charmed meson, 3rd position in pdg code == 4
  if(motherpdg / 1000 == 4) return kTRUE;            // charmed baryon, 4th position in pdg code == 4
  AliDebug(1, "No Charm\n");
  return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsBeautyElectron(const TObject * const o) const {
  //
  // Check if mother is coming from Beauty
  //
  if(TMath::Abs(GetTrackPDG(dynamic_cast<const AliVParticle *>(o))) != 11) return kFALSE;
  Int_t motherpdg = TMath::Abs(GetMotherPDG(dynamic_cast<const AliVParticle *>(o)));
  AliDebug(1, Form("Mother PDG %d\n", motherpdg));

  if((motherpdg % 1000) / 100 == 5) return kTRUE;   // beauty meson, 3rd position in pdg code == 5
  if(motherpdg / 1000 == 5) return kTRUE;           // beauty baryon, 4th position in pdg code == 5   
  AliDebug(1, "No Beauty\n");
  return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEsignalCuts::IsGammaElectron(const TObject * const o) const {
  //
  // Check for MC if the electron is coming from Gamma
  //
  if(TMath::Abs(GetTrackPDG(dynamic_cast<const AliVParticle *>(o))) != 11) return kFALSE;
  Int_t motherpdg = TMath::Abs(GetMotherPDG(dynamic_cast<const AliVParticle *>(o)));
  AliDebug(1, Form("Mother PDG %d\n", motherpdg));

  if(motherpdg!=22){
    AliDebug(1, "No Gamma");
    return kFALSE;
  } else { 
    AliDebug(1, "Gamma");
    return kTRUE;
  }
}
*/

//____________________________________________________________
Int_t AliHFEsignalCuts::GetMotherPDG(const AliVParticle * const track) const {
  //
  // Get Mother Pdg code for reconstructed respectively MC tracks
  // 
  if(!fMC){
    AliDebug(1, "No MC Event Available\n");
    return 0;
  }
  const AliVParticle *motherParticle = NULL, *mctrack = NULL;
  TString objectType = track->IsA()->GetName();
  if(objectType.CompareTo("AliESDtrack") == 0 || objectType.CompareTo("AliAODTrack") == 0){
    // Reconstructed track
    if(track->GetLabel())
      mctrack = fMC->GetTrack(TMath::Abs(track->GetLabel()));
  } else {
    // MCParticle
    mctrack = track;
  }

  if(!mctrack) return 0;
  
  Int_t motherPDG = 0;
  if(TString(mctrack->IsA()->GetName()).CompareTo("AliMCParticle") == 0){
    // case MC Particle
    const AliMCParticle *esdmctrack = dynamic_cast<const AliMCParticle *>(mctrack);
    if(esdmctrack) motherParticle = fMC->GetTrack(esdmctrack->Particle()->GetFirstMother());
    if(motherParticle){
      const AliMCParticle *esdmcmother = dynamic_cast<const AliMCParticle *>(motherParticle);
      if(esdmcmother) motherPDG = TMath::Abs(esdmcmother->Particle()->GetPdgCode());
    }
  } else {
    // case AODMCParticle
    const AliAODMCParticle *aodmctrack = dynamic_cast<const AliAODMCParticle *>(mctrack);
    if(aodmctrack) motherParticle = fMC->GetTrack(aodmctrack->GetMother());
    if(motherParticle){
      const AliAODMCParticle *aodmcmother = dynamic_cast<const AliAODMCParticle *>(motherParticle);
      if(aodmcmother) motherPDG = TMath::Abs(aodmcmother->GetPdgCode());
    }
  }
  return motherPDG;
}

//____________________________________________________________
Int_t AliHFEsignalCuts::GetTrackPDG(const AliVParticle * const track) const {
	//
	// Return PDG code of a particle itself
	//
  if(!fMC){
    AliDebug(1, "No MC Event Available\n");
    return 0;
  }
	TString sourcetype = track->IsA()->GetName();
	const AliVParticle *mctrack = NULL;
	if(!sourcetype.CompareTo("AliESDtrack") || !sourcetype.CompareTo("AliAODTrack")){
		mctrack = fMC->GetTrack(TMath::Abs(track->GetLabel()));
	} else  mctrack = track;
	if(!mctrack) return 0;

	TString mctype = mctrack->IsA()->GetName();
	Int_t trackPdg = 0;
	if(!mctype.CompareTo("AliMCParticle")){
		const AliMCParticle *esdmc = dynamic_cast<const AliMCParticle *>(mctrack);
		if(esdmc) trackPdg = esdmc->Particle()->GetPdgCode();
	} else {
		const AliAODMCParticle *aodmc = dynamic_cast< const AliAODMCParticle *>(mctrack);
		if(aodmc) trackPdg = aodmc->GetPdgCode();
	}
	return trackPdg;
}

//____________________________________________________________
Int_t AliHFEsignalCuts::GetElecSource(const AliVParticle * const track) const {
	//
	// Return PDG code of a particle itself
	//
	
  if(!fMC){
    AliDebug(1, "No MC Event Available\n");
    return 0;
  }
  if(!fMCQA){
    AliDebug(1, "No MCQA Available\n");
    return 0;
  }
  if(!track){
    AliDebug(1, "Track not Available\n");
    return 0;
  }

  TString sourcetype = track->IsA()->GetName();
  const AliVParticle *mctrack = NULL;
  TParticle *mcpart = NULL;
  if(!sourcetype.CompareTo("AliESDtrack") || !sourcetype.CompareTo("AliAODTrack")){
    mctrack = fMC->GetTrack(TMath::Abs(track->GetLabel()));
  } else  mctrack = track;
  if(!mctrack) return 0;

  TString mctype = mctrack->IsA()->GetName();
  Int_t eSource = 0;
  if(!mctype.CompareTo("AliMCParticle")){
    const AliMCParticle *esdmc = dynamic_cast<const AliMCParticle *>(mctrack);
    if(esdmc){
      mcpart = esdmc->Particle();
      eSource=fMCQA->GetElecSource(mcpart);
    }
  } else {
    return -1;
  }
  return eSource;
}	
