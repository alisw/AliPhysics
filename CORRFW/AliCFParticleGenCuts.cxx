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

////////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// class AliCFParticleGenCuts implementation
// Using this class a user may define selections relative to 
// MC particle (AliMCParticle) using generation-level information.
////////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliCFParticleGenCuts.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "TObject.h"
#include "AliStack.h"

ClassImp(AliCFParticleGenCuts)

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts() : 
  AliCFCutBase(),
  fMCInfo(0x0),
  fRequireIsCharged(0),
  fRequireIsPrimary(0),
  fRequireIsSecondary(0),
  fRequirePdgCode(0),
  fPdgCode(0),
  fProdVtxXMin (-1.e+09),
  fProdVtxYMin (-1.e+09),
  fProdVtxZMin (-1.e+09),
  fProdVtxXMax ( 1.e+09),
  fProdVtxYMax ( 1.e+09),
  fProdVtxZMax ( 1.e+09),
  fDecayVtxXMin(-1.e+09),
  fDecayVtxYMin(-1.e+09),
  fDecayVtxZMin(-1.e+09),
  fDecayVtxXMax( 1.e+09),
  fDecayVtxYMax( 1.e+09),
  fDecayVtxZMax( 1.e+09),
  fDecayLengthMin(0),
  fDecayLengthMax(1.e+09),
  fDecayRxyMin(0),
  fDecayRxyMax(1.e+09)
{
  //
  //ctor
  //
}

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts(const Char_t* name, const Char_t* title) : 
  AliCFCutBase(name,title),
  fMCInfo(0x0),
  fRequireIsCharged(0),
  fRequireIsPrimary(0),
  fRequireIsSecondary(0),
  fRequirePdgCode(0),
  fPdgCode(0),
  fProdVtxXMin (-1.e+09),
  fProdVtxYMin (-1.e+09),
  fProdVtxZMin (-1.e+09),
  fProdVtxXMax ( 1.e+09),
  fProdVtxYMax ( 1.e+09),
  fProdVtxZMax ( 1.e+09),
  fDecayVtxXMin(-1.e+09),
  fDecayVtxYMin(-1.e+09),
  fDecayVtxZMin(-1.e+09),
  fDecayVtxXMax( 1.e+09),
  fDecayVtxYMax( 1.e+09),
  fDecayVtxZMax( 1.e+09),
  fDecayLengthMin(0),
  fDecayLengthMax(1.e+09),
  fDecayRxyMin(0),
  fDecayRxyMax(1.e+09)
{
  //
  //ctor
  //
}

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts(const AliCFParticleGenCuts& c) : 
  AliCFCutBase(c),
  fMCInfo(c.fMCInfo),
  fRequireIsCharged(c.fRequireIsCharged),
  fRequireIsPrimary(c.fRequireIsPrimary),
  fRequireIsSecondary(c.fRequireIsSecondary),
  fRequirePdgCode(c.fRequirePdgCode),
  fPdgCode(c.fPdgCode),
  fProdVtxXMin (c.fProdVtxXMin),
  fProdVtxYMin (c.fProdVtxYMin),
  fProdVtxZMin (c.fProdVtxZMin),
  fProdVtxXMax (c.fProdVtxXMax),
  fProdVtxYMax (c.fProdVtxYMax),
  fProdVtxZMax (c.fProdVtxZMax),
  fDecayVtxXMin(c.fDecayVtxXMin),
  fDecayVtxYMin(c.fDecayVtxYMin),
  fDecayVtxZMin(c.fDecayVtxZMin),
  fDecayVtxXMax(c.fDecayVtxXMax),
  fDecayVtxYMax(c.fDecayVtxYMax),
  fDecayVtxZMax(c.fDecayVtxZMax),
  fDecayLengthMin(c.fDecayLengthMin),
  fDecayLengthMax(c.fDecayLengthMin),
  fDecayRxyMin(c.fDecayLengthMin),
  fDecayRxyMax(c.fDecayLengthMin)
{
  //
  //copy ctor
  //
}

//______________________________
AliCFParticleGenCuts& AliCFParticleGenCuts::operator=(const AliCFParticleGenCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fMCInfo=c.fMCInfo;
    fRequireIsCharged=c.fRequireIsCharged;
    fRequireIsPrimary=c.fRequireIsPrimary;
    fRequireIsSecondary=c.fRequireIsSecondary;
    fRequirePdgCode=c.fRequirePdgCode;
    fPdgCode=c.fPdgCode;
    fProdVtxXMin=c.fProdVtxXMin;
    fProdVtxYMin=c.fProdVtxYMin;
    fProdVtxZMin=c.fProdVtxZMin;
    fProdVtxXMax=c.fProdVtxXMax;
    fProdVtxYMax=c.fProdVtxYMax;
    fProdVtxZMax=c.fProdVtxZMax;
    fDecayVtxXMin=c.fDecayVtxXMin;
    fDecayVtxYMin=c.fDecayVtxYMin;
    fDecayVtxZMin=c.fDecayVtxZMin;
    fDecayVtxXMax=c.fDecayVtxXMax;
    fDecayVtxYMax=c.fDecayVtxYMax;
    fDecayVtxZMax=c.fDecayVtxZMax;      
    fDecayLengthMin=c.fDecayVtxZMax;
    fDecayLengthMax=c.fDecayLengthMax;
    fDecayRxyMin=c.fDecayRxyMin;
    fDecayRxyMax=c.fDecayRxyMax;
  }
  return *this ;
}

//______________________________
Bool_t AliCFParticleGenCuts::IsSelected(TObject* obj) {
  //
  // check if selections on 'obj' are passed
  // 'obj' must be an AliMCParticle
  //
  
  if (!obj) return kFALSE ;
  TString className(obj->ClassName());
  if (className.CompareTo("AliMCParticle") != 0) {
    AliError("argument must point to an AliMCParticle !");
    return kFALSE ;
  }

  AliMCParticle* mcPart = (AliMCParticle*) obj ;
  TParticle* part = mcPart->Particle();
  AliStack *stack=fMCInfo->MCEvent()->Stack();

  // is this particle charged?
  if ( fRequireIsCharged ) {
    if(!IsCharged(mcPart))return kFALSE;
  } 
  
  // primary cuts
  if ( fRequireIsPrimary ) {
     if(!IsPrimary(mcPart,stack))return kFALSE;
  } 

  //secondary cut
  if ( fRequireIsSecondary && part->IsPrimary() ) return kFALSE ;
  
  //PDG code cut
  if ( fRequirePdgCode){
    if(!IsA(mcPart,fPdgCode))  return kFALSE ;
  }
  // production vertex cuts
  Double32_t partVx=(Double32_t)part->Vx();
  Double32_t partVy=(Double32_t)part->Vy();
  Double32_t partVz=(Double32_t)part->Vz();
  if ( partVx < fProdVtxXMin || partVx > fProdVtxXMax ) return kFALSE ;
  if ( partVy < fProdVtxYMin || partVy > fProdVtxYMax ) return kFALSE ;
  if ( partVz < fProdVtxZMin || partVz > fProdVtxZMax ) return kFALSE ;

  //decay vertex cuts
  if ( part->GetNDaughters() > 0 ) {
    TParticle* daughter = fMCInfo->MCEvent()->Stack()->Particle(part->GetFirstDaughter()) ;
    Double32_t decayVx=(Double32_t)daughter->Vx();
    Double32_t decayVy=(Double32_t)daughter->Vy();
    Double32_t decayVz=(Double32_t)daughter->Vz();
    if ( decayVx < fDecayVtxXMin || decayVx > fDecayVtxXMax ) return kFALSE ;
    if ( decayVy < fDecayVtxYMin || decayVy > fDecayVtxYMax ) return kFALSE ;
    if ( decayVz < fDecayVtxZMin || decayVz > fDecayVtxZMax ) return kFALSE ;

    //decay length cut
    Double32_t decayL = TMath::Sqrt(TMath::Power(partVx-decayVx,2) + 
				    TMath::Power(partVy-decayVy,2) + 
				    TMath::Power(partVz-decayVz,2) ) ;
    if (decayL < fDecayLengthMin || decayL > fDecayLengthMax) return kFALSE ;

    Double32_t decayRxy = TMath::Sqrt(TMath::Power(decayVx,2) + 
				      TMath::Power(decayVy,2) ) ;
    if (decayRxy < fDecayRxyMin || decayRxy > fDecayRxyMax) return kFALSE ;
  }


  return kTRUE ;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsCharged(AliMCParticle *mcPart) {
  //
  //check if particle is charged.
  //
  TParticle* part = mcPart->Particle();
  TParticlePDG* pdgPart = part->GetPDG();
  if(!pdgPart)return kFALSE;
  if (pdgPart->Charge() == 0) return kFALSE;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsPrimary(AliMCParticle *mcPart, AliStack *stack) {
  //
  //check if particle is primary (standard definition)
  //
  if (!stack->IsPhysicalPrimary(mcPart->Label())) return kFALSE ;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsPrimaryCharged(AliMCParticle *mcPart, AliStack *stack) {
  //
  //check if a charged particle is primary (standard definition)
  //
  if (!stack->IsPhysicalPrimary(mcPart->Label()) || !IsCharged(mcPart)) return kFALSE ;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsA(AliMCParticle *mcPart, Int_t pdg, Bool_t abs) {
  //
  //Check on the pdg code of the MC particle. if abs=kTRUE then check on the 
  //absolute value. By default is set to kFALSE.
  //
  TParticle* part = mcPart->Particle();
  Int_t pdgCode = part->GetPdgCode();
  if(abs)pdgCode=TMath::Abs(pdgCode);
  if(pdgCode != pdg )return kFALSE;
  return kTRUE;
}
//______________________________
void AliCFParticleGenCuts::SetEvtInfo(TObject* mcInfo) {
  //
  // Sets pointer to MC event information (AliMCEventHandler)
  //

  if (!mcInfo) {
    AliError("Pointer to MC Event Handler is null !");
    return;
  }
  
  TString className(mcInfo->ClassName());
  if (className.CompareTo("AliMCEventHandler") != 0) {
    AliError("argument must point to an AliMCEventHandler !");
    return ;
  }
  
  fMCInfo = (AliMCEventHandler*) mcInfo ;
}
