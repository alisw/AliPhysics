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

// AliFlowTrackCuts:
// ESD track cuts for flow framework 
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)
//
// This class gurantees consistency of cut methods, trackparameter
// selection (global tracks, TPC only, etc..) and parameter mixing
// in the flow framework. Transparently handles different input types:
// ESD, MC, AOD.
// This class works in 2 steps: first the requested track parameters are
// constructed (to be set by SetParamType() ), then cuts are applied.
// the constructed track can be requested AFTER checking the cuts by
// calling GetTrack(), in this case the cut object stays in control,
// caller does not have to delete the track.
// Additionally caller can request an AliFlowTrack object to be constructed
// according the parameter mixing scenario requested by SetParamMix().
// AliFlowTrack is made using MakeFlowTrack() method, its an 'object factory'
// so caller needs to take care of the freshly created object.

#include <limits.h>
#include <float.h>
#include "TParticle.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliAODTrack.h"
#include "AliFlowTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliLog.h"

ClassImp(AliFlowTrackCuts)

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts():
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(new AliESDtrackCuts()),
  fQA(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fIgnoreSignInPID(kFALSE),
  fCutMCisPrimary(kFALSE),
  fMCisPrimary(kFALSE),
  fRequireCharge(kFALSE),
  fFakesAreOK(kTRUE),
  fCutSPDtrackletDeltaPhi(kFALSE),
  fSPDtrackletDeltaPhiMax(FLT_MAX),
  fSPDtrackletDeltaPhiMin(-FLT_MAX),
  fParamType(kGlobal),
  fParamMix(kPure),
  fCleanupTrack(kFALSE),
  fTrack(NULL),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(0.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const AliFlowTrackCuts& that):
  AliFlowTrackSimpleCuts(that),
  fAliESDtrackCuts(new AliESDtrackCuts(*(that.fAliESDtrackCuts))),
  fQA(that.fQA),
  fCutMCprocessType(that.fCutMCprocessType),
  fMCprocessType(that.fMCprocessType),
  fCutMCPID(that.fCutMCPID),
  fMCPID(that.fMCPID),
  fIgnoreSignInPID(that.fIgnoreSignInPID),
  fCutMCisPrimary(that.fCutMCisPrimary),
  fMCisPrimary(that.fMCisPrimary),
  fRequireCharge(that.fRequireCharge),
  fFakesAreOK(that.fFakesAreOK),
  fCutSPDtrackletDeltaPhi(that.fCutSPDtrackletDeltaPhi),
  fSPDtrackletDeltaPhiMax(that.fSPDtrackletDeltaPhiMax),
  fSPDtrackletDeltaPhiMin(that.fSPDtrackletDeltaPhiMin),
  fParamType(that.fParamType),
  fParamMix(that.fParamMix),
  fCleanupTrack(kFALSE),
  fTrack(NULL),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(0.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL)
{
  //copy constructor
}

//-----------------------------------------------------------------------
AliFlowTrackCuts& AliFlowTrackCuts::operator=(const AliFlowTrackCuts& that)
{
  //assignment
  AliFlowTrackSimpleCuts::operator=(that);
  *fAliESDtrackCuts=*(that.fAliESDtrackCuts);
  fQA=that.fQA;
  fCutMCprocessType=that.fCutMCprocessType;
  fMCprocessType=that.fMCprocessType;
  fCutMCPID=that.fCutMCPID;
  fMCPID=that.fMCPID;
  fIgnoreSignInPID=that.fIgnoreSignInPID,
  fCutMCisPrimary=that.fCutMCisPrimary;
  fMCisPrimary=that.fMCisPrimary;
  fRequireCharge=that.fRequireCharge;
  fFakesAreOK=that.fFakesAreOK;
  fCutSPDtrackletDeltaPhi=that.fCutSPDtrackletDeltaPhi;
  fSPDtrackletDeltaPhiMax=that.fSPDtrackletDeltaPhiMax;
  fSPDtrackletDeltaPhiMin=that.fSPDtrackletDeltaPhiMin;
  fParamType=that.fParamType;
  fParamMix=that.fParamMix;

  fCleanupTrack=kFALSE;
  fTrack=NULL;
  fTrackPhi=0.;
  fTrackPhi=0.;
  fTrackWeight=0.;
  fTrackLabel=INT_MIN;
  fMCevent=NULL;
  fMCparticle=NULL;
  fEvent=NULL;

  return *this;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::~AliFlowTrackCuts()
{
  //dtor
  if (fCleanupTrack) delete fTrack;
  delete fAliESDtrackCuts;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsSelected(TObject* obj, Int_t id)
{
  //check cuts
  AliVParticle* vparticle = dynamic_cast<AliVParticle*>(obj);
  if (vparticle) return PassesCuts(vparticle);
  AliFlowTrackSimple* flowtrack = dynamic_cast<AliFlowTrackSimple*>(obj);
  if (flowtrack) return PassesCuts(flowtrack);
  AliMultiplicity* tracklets = dynamic_cast<AliMultiplicity*>(obj);
  if (tracklets) return PassesCuts(tracklets,id);
  return kFALSE;  //default when passed wrong type of object
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliFlowTrackSimple* track)
{
  //check cuts on a flowtracksimple

  //clean up from last iteration
  if (fCleanupTrack) delete fTrack; fTrack = NULL;
  return AliFlowTrackSimpleCuts::PassesCuts(track);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliMultiplicity* tracklet, Int_t id)
{
  //check cuts on a tracklets

  //clean up from last iteration, and init label
  if (fCleanupTrack) delete fTrack; fTrack = NULL;
  fMCparticle=NULL;
  fTrackLabel=-1;

  fTrackPhi = tracklet->GetPhi(id);
  fTrackEta = tracklet->GetEta(id);
  fTrackWeight = 1.0;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) return kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) return kFALSE;}

  //check MC info if available
  //if the 2 clusters have different label track cannot be good
  //and should therefore not pass the mc cuts
  Int_t label0 = tracklet->GetLabel(id,0);
  Int_t label1 = tracklet->GetLabel(id,1);
  fTrackLabel = (label0==label1)?tracklet->GetLabel(id,1):-1;
  if (!PassesMCcuts()) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts()
{
  //check the MC info
  if (!fMCevent) return kFALSE;
  if (fTrackLabel<0) return kFALSE;//otherwise AliCMevent prints a warning before returning NULL
  fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  if (!fMCparticle) {AliError("no MC track"); return kFALSE;}

  if (fCutMCisPrimary)
  {
    if (IsPhysicalPrimary(fMCevent,fTrackLabel) != fMCisPrimary) return kFALSE;
  }
  if (fCutMCPID)
  {
    Int_t pdgCode = fMCparticle->PdgCode();
    if (fIgnoreSignInPID) 
    {
      if (TMath::Abs(fMCPID) != TMath::Abs(pdgCode)) return kFALSE;
    }
    else 
    {
      if (fMCPID != pdgCode) return kFALSE;
    }
  }
  if ( fCutMCprocessType )
  {
    TParticle* particle = fMCparticle->Particle();
    Int_t processID = particle->GetUniqueID();
    if (processID != fMCprocessType ) return kFALSE;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliVParticle* vparticle)
{
  //check cuts for an ESD vparticle

  ////////////////////////////////////////////////////////////////
  //  start by preparing the track parameters to cut on //////////
  ////////////////////////////////////////////////////////////////
  //clean up from last iteration
  if (fCleanupTrack) delete fTrack; fTrack=NULL; 

  //get the label and the mc particle
  fTrackLabel = (fFakesAreOK)?TMath::Abs(vparticle->GetLabel()):vparticle->GetLabel();
  if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  else fMCparticle=NULL;

  Bool_t isMCparticle = kFALSE; //some things are different for MC particles, check!
  AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*>(vparticle);
  if (esdTrack)
    HandleESDtrack(esdTrack);
  else
  {
    HandleVParticle(vparticle);
    //now check if produced particle is MC
    isMCparticle = (dynamic_cast<AliMCParticle*>(fTrack))!=NULL;
  }
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  Bool_t pass=kTRUE;
  //check the common cuts for the current particle fTrack (MC,AOD,ESD)
  if (fCutPt) {if (fTrack->Pt() < fPtMin || fTrack->Pt() >= fPtMax ) pass=kFALSE;}
  if (fCutEta) {if (fTrack->Eta() < fEtaMin || fTrack->Eta() >= fEtaMax ) pass=kFALSE;}
  if (fCutPhi) {if (fTrack->Phi() < fPhiMin || fTrack->Phi() >= fPhiMax ) pass=kFALSE;}
  if (fRequireCharge) {if (fTrack->Charge() == 0) pass=kFALSE;}
  if (fCutCharge && !isMCparticle) {if (fTrack->Charge() != fCharge) pass=kFALSE;}
  if (fCutCharge && isMCparticle)
  { 
    //in case of an MC particle the charge is stored in units of 1/3|e| 
    Int_t charge = TMath::Nint(fTrack->Charge()/3.0); //mc particles have charge in units of 1/3e
    return (charge==fCharge);
  }
  //if(fCutPID) {if (fTrack->PID() != fPID) pass=kFALSE;}

  //when additionally MC info is required
  if (!PassesMCcuts()) pass=kFALSE;

  //check all else for ESDs using aliesdtrackcuts
  if (esdTrack && (fParamType!=kMC) ) 
    if (!fAliESDtrackCuts->IsSelected(static_cast<AliESDtrack*>(fTrack))) pass=kFALSE;

  return pass; //true by default, if we didn't set any cuts
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::HandleVParticle(AliVParticle* track)
{
  //handle the general case
  switch (fParamType)
  {
    default:
      fCleanupTrack = kFALSE;
      fTrack = track;
  }
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::HandleESDtrack(AliESDtrack* track)
{
  //handle esd track
  switch (fParamType)
  {
    case kGlobal:
      fTrack = track;
      fCleanupTrack = kFALSE;
      break;
    case kESD_TPConly:
      fTrack = new AliESDtrack();
      track->FillTPCOnlyTrack(*(static_cast<AliESDtrack*>(fTrack)));
      fCleanupTrack = kTRUE;
      //recalculate the label and mc particle, they may differ as TPClabel != global label
      fTrackLabel = (fFakesAreOK)?TMath::Abs(fTrack->GetLabel()):fTrack->GetLabel();
      if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
      else fMCparticle=NULL;
      break;
    default:
      fTrack = track;
      fCleanupTrack = kFALSE;
  }
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCOnlyTrackCuts()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts();
  cuts->SetName("standard TPConly cuts");
  delete cuts->fAliESDtrackCuts;
  cuts->fAliESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  cuts->SetParamType(kESD_TPConly);
  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries)
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts();
  cuts->SetName("standard global track cuts 2009");
  delete cuts->fAliESDtrackCuts;
  cuts->fAliESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selPrimaries);
  cuts->SetParamType(kGlobal);
  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrack() const
{
  //get a flow track constructed from whatever we applied cuts on
  //caller is resposible for deletion
  //if construction fails return NULL
  AliFlowTrack* flowtrack=NULL;
  if (fParamType==kESD_SPDtracklet)
  {
    flowtrack = new AliFlowTrack();
    switch (fParamMix)
    {
      case kPure:
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        break;
      case kTrackWithMCkine:
        if (!fMCparticle) return NULL;
        flowtrack->SetPhi( fMCparticle->Phi() );
        flowtrack->SetEta( fMCparticle->Eta() );
        flowtrack->SetPt( fMCparticle->Pt() );
        break;
      case kTrackWithMCpt:
        if (!fMCparticle) return NULL;
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        flowtrack->SetPt(fMCparticle->Pt());
        break;
      default:
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
    }
    flowtrack->SetSource(AliFlowTrack::kFromTracklet);
  }
  else
  {
    switch(fParamMix)
    {
      case kPure:
        flowtrack = new AliFlowTrack(fTrack);
        break;
      case kTrackWithMCkine:
        flowtrack = new AliFlowTrack(fMCparticle);
        break;
      case kTrackWithMCPID:
        flowtrack = new AliFlowTrack(fTrack);
        //flowtrack->setPID(...) from mc, when implemented
        break;
      case kTrackWithMCpt:
        if (!fMCparticle) return NULL;
        flowtrack = new AliFlowTrack(fTrack);
        flowtrack->SetPt(fMCparticle->Pt());
      default:
        flowtrack = new AliFlowTrack(fTrack);
    }
    if (fParamType==kMC) flowtrack->SetSource(AliFlowTrack::kFromMC);
    else if (dynamic_cast<AliESDtrack*>(fTrack)) flowtrack->SetSource(AliFlowTrack::kFromESD);
    else if (dynamic_cast<AliAODTrack*>(fTrack)) flowtrack->SetSource(AliFlowTrack::kFromAOD);
    else if (dynamic_cast<AliMCParticle*>(fTrack)) flowtrack->SetSource(AliFlowTrack::kFromMC);
  }
  return flowtrack;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsPhysicalPrimary() const
{
  //check if current particle is a physical primary
  if (!fMCevent) return kFALSE;
  if (fTrackLabel<0) return kFALSE;
  return IsPhysicalPrimary(fMCevent, fTrackLabel);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsPhysicalPrimary(AliMCEvent* mcEvent, Int_t label)
{
  //check if current particle is a physical primary
  Bool_t physprim=mcEvent->IsPhysicalPrimary(label);
  if (!physprim) return kFALSE;
  AliMCParticle* track = static_cast<AliMCParticle*>(mcEvent->GetTrack(label));
  if (!track) return kFALSE;
  TParticle* particle = track->Particle();
  Bool_t transported = particle->TestBit(kTransportBit);
  //printf("prim: %s, transp: %s\n",(physprim)?"YES":"NO ",(transported)?"YES":"NO ");
  return (physprim && transported);
}

//-----------------------------------------------------------------------
const char* AliFlowTrackCuts::GetParamTypeName(trackParameterType type) 
{
  //return the name of the selected parameter type
  switch (type)
  {
    case kMC:
      return "MC";
    case kGlobal:
      return "ESD global";
    case kESD_TPConly:
      return "TPC only";
    case kESD_SPDtracklet:
        return "SPD tracklet";
    default:
        return "unknown";
  }
  return "unknown";
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::DefineHistograms()
{
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::GetNumberOfInputObjects() const
{
  //get the number of tracks in the input event according source
  //selection (ESD tracks, tracklets, MC particles etc.)
  AliESDEvent* esd=NULL;
  switch (fParamType)
  {
    case kESD_SPDtracklet:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetMultiplicity()->GetNumberOfTracklets();
    case kMC:
      if (!fMCevent) return 0;
      return fMCevent->GetNumberOfTracks();
    default:
      if (!fEvent) return 0;
      return fEvent->GetNumberOfTracks();
  }
  return 0;
}

//-----------------------------------------------------------------------
TObject* AliFlowTrackCuts::GetInputObject(Int_t i)
{
  //get the input object according the data source selection:
  //(esd tracks, traclets, mc particles,etc...)
  AliESDEvent* esd=NULL;
  switch (fParamType)
  {
    case kESD_SPDtracklet:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return const_cast<AliMultiplicity*>(esd->GetMultiplicity());
    case kMC:
      if (!fMCevent) return NULL;
      return fMCevent->GetTrack(i);
    default:
      if (!fEvent) return NULL;
      return fEvent->GetTrack(i);
  }
}
