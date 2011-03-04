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
#include <TMatrix.h>
#include "TParticle.h"
#include "TObjArray.h"
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
#include "AliESDpid.h"

ClassImp(AliFlowTrackCuts)

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts():
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(NULL),
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fIgnoreSignInPID(kFALSE),
  fCutMCisPrimary(kFALSE),
  fRequireTransportBitForPrimaries(kTRUE),
  fMCisPrimary(kFALSE),
  fRequireCharge(kFALSE),
  fFakesAreOK(kTRUE),
  fCutSPDtrackletDeltaPhi(kFALSE),
  fSPDtrackletDeltaPhiMax(FLT_MAX),
  fSPDtrackletDeltaPhiMin(-FLT_MAX),
  fIgnoreTPCzRange(kFALSE),
  fIgnoreTPCzRangeMax(FLT_MAX),
  fIgnoreTPCzRangeMin(-FLT_MAX),
  fCutChi2PerClusterTPC(kFALSE),
  fMaxChi2PerClusterTPC(FLT_MAX),
  fMinChi2PerClusterTPC(-FLT_MAX),
  fCutNClustersTPC(kFALSE),
  fNClustersTPCMax(INT_MAX),
  fNClustersTPCMin(INT_MIN),  
  fCutNClustersITS(kFALSE),
  fNClustersITSMax(INT_MAX),
  fNClustersITSMin(INT_MIN),  
  fUseAODFilterBit(kFALSE),
  fAODFilterBit(0),
  fCutDCAToVertexXY(kFALSE),
  fCutDCAToVertexZ(kFALSE),
  fCutMinimalTPCdedx(kFALSE),
  fMinimalTPCdedx(0.),
  fParamType(kGlobal),
  fParamMix(kPure),
  fTrack(NULL),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(0.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(),
  fPIDsource(kTOFpid),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(AliPID::kPion),
  fParticleProbability(.9)
{
  //io constructor 
  for ( Int_t i=0; i<5; i++ ) { fProbBayes[i]=0.0; }
  SetPriors(); //init arrays
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const char* name):
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(new AliESDtrackCuts()),
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fIgnoreSignInPID(kFALSE),
  fCutMCisPrimary(kFALSE),
  fRequireTransportBitForPrimaries(kTRUE),
  fMCisPrimary(kFALSE),
  fRequireCharge(kFALSE),
  fFakesAreOK(kTRUE),
  fCutSPDtrackletDeltaPhi(kFALSE),
  fSPDtrackletDeltaPhiMax(FLT_MAX),
  fSPDtrackletDeltaPhiMin(-FLT_MAX),
  fIgnoreTPCzRange(kFALSE),
  fIgnoreTPCzRangeMax(FLT_MAX),
  fIgnoreTPCzRangeMin(-FLT_MAX),
  fCutChi2PerClusterTPC(kFALSE),
  fMaxChi2PerClusterTPC(FLT_MAX),
  fMinChi2PerClusterTPC(-FLT_MAX),
  fCutNClustersTPC(kFALSE),
  fNClustersTPCMax(INT_MAX),
  fNClustersTPCMin(INT_MIN),  
  fCutNClustersITS(kFALSE),
  fNClustersITSMax(INT_MAX),
  fNClustersITSMin(INT_MIN),
  fUseAODFilterBit(kFALSE),
  fAODFilterBit(0),
  fCutDCAToVertexXY(kFALSE),
  fCutDCAToVertexZ(kFALSE),
  fCutMinimalTPCdedx(kFALSE),
  fMinimalTPCdedx(0.),
  fParamType(kGlobal),
  fParamMix(kPure),
  fTrack(NULL),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(0.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(),
  fPIDsource(kTOFpid),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(AliPID::kPion),
  fParticleProbability(.9)
{
  //constructor 
  SetName(name);
  SetTitle("AliFlowTrackCuts");
  fESDpid.GetTPCResponse().SetBetheBlochParameters( 0.0283086,
                                                    2.63394e+01,
                                                    5.04114e-11,
                                                    2.12543e+00,
                                                    4.88663e+00 );
  for ( Int_t i=0; i<5; i++ ) { fProbBayes[i]=0.0; }
  SetPriors(); //init arrays
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const AliFlowTrackCuts& that):
  AliFlowTrackSimpleCuts(that),
  fAliESDtrackCuts(NULL),
  fQA(NULL),
  fCutMC(that.fCutMC),
  fCutMCprocessType(that.fCutMCprocessType),
  fMCprocessType(that.fMCprocessType),
  fCutMCPID(that.fCutMCPID),
  fMCPID(that.fMCPID),
  fIgnoreSignInPID(that.fIgnoreSignInPID),
  fCutMCisPrimary(that.fCutMCisPrimary),
  fRequireTransportBitForPrimaries(that.fRequireTransportBitForPrimaries),
  fMCisPrimary(that.fMCisPrimary),
  fRequireCharge(that.fRequireCharge),
  fFakesAreOK(that.fFakesAreOK),
  fCutSPDtrackletDeltaPhi(that.fCutSPDtrackletDeltaPhi),
  fSPDtrackletDeltaPhiMax(that.fSPDtrackletDeltaPhiMax),
  fSPDtrackletDeltaPhiMin(that.fSPDtrackletDeltaPhiMin),
  fIgnoreTPCzRange(that.fIgnoreTPCzRange),
  fIgnoreTPCzRangeMax(that.fIgnoreTPCzRangeMax),
  fIgnoreTPCzRangeMin(that.fIgnoreTPCzRangeMin),
  fCutChi2PerClusterTPC(that.fCutChi2PerClusterTPC),
  fMaxChi2PerClusterTPC(that.fMaxChi2PerClusterTPC),
  fMinChi2PerClusterTPC(that.fMinChi2PerClusterTPC),
  fCutNClustersTPC(that.fCutNClustersTPC),
  fNClustersTPCMax(that.fNClustersTPCMax),
  fNClustersTPCMin(that.fNClustersTPCMin),
  fCutNClustersITS(that.fCutNClustersITS),
  fNClustersITSMax(that.fNClustersITSMax),
  fNClustersITSMin(that.fNClustersITSMin),
  fUseAODFilterBit(that.fUseAODFilterBit),
  fAODFilterBit(that.fAODFilterBit),
  fCutDCAToVertexXY(that.fCutDCAToVertexXY),
  fCutDCAToVertexZ(that.fCutDCAToVertexZ),
  fCutMinimalTPCdedx(that.fCutMinimalTPCdedx),
  fMinimalTPCdedx(that.fMinimalTPCdedx),
  fParamType(that.fParamType),
  fParamMix(that.fParamMix),
  fTrack(NULL),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(0.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(that.fESDpid),
  fPIDsource(that.fPIDsource),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(that.fParticleID),
  fParticleProbability(that.fParticleProbability)
{
  //copy constructor
  if (that.fTPCpidCuts) fTPCpidCuts = new TMatrixF(*(that.fTPCpidCuts));
  if (that.fTOFpidCuts) fTOFpidCuts = new TMatrixF(*(that.fTOFpidCuts));
  if (that.fAliESDtrackCuts) fAliESDtrackCuts = new AliESDtrackCuts(*(that.fAliESDtrackCuts));
  memcpy(fProbBayes,that.fProbBayes,sizeof(fProbBayes));
  SetPriors(); //init arrays
}

//-----------------------------------------------------------------------
AliFlowTrackCuts& AliFlowTrackCuts::operator=(const AliFlowTrackCuts& that)
{
  //assignment
  if (this==&that) return *this;

  AliFlowTrackSimpleCuts::operator=(that);
  if (that.fAliESDtrackCuts) *fAliESDtrackCuts=*(that.fAliESDtrackCuts);
  fQA=NULL;
  fCutMC=that.fCutMC;
  fCutMCprocessType=that.fCutMCprocessType;
  fMCprocessType=that.fMCprocessType;
  fCutMCPID=that.fCutMCPID;
  fMCPID=that.fMCPID;
  fIgnoreSignInPID=that.fIgnoreSignInPID,
  fCutMCisPrimary=that.fCutMCisPrimary;
  fRequireTransportBitForPrimaries=that.fRequireTransportBitForPrimaries;
  fMCisPrimary=that.fMCisPrimary;
  fRequireCharge=that.fRequireCharge;
  fFakesAreOK=that.fFakesAreOK;
  fCutSPDtrackletDeltaPhi=that.fCutSPDtrackletDeltaPhi;
  fSPDtrackletDeltaPhiMax=that.fSPDtrackletDeltaPhiMax;
  fSPDtrackletDeltaPhiMin=that.fSPDtrackletDeltaPhiMin;
  fIgnoreTPCzRange=that.fIgnoreTPCzRange;
  fIgnoreTPCzRangeMax=that.fIgnoreTPCzRangeMax;
  fIgnoreTPCzRangeMin=that.fIgnoreTPCzRangeMin;
  fCutChi2PerClusterTPC=that.fCutChi2PerClusterTPC;
  fMaxChi2PerClusterTPC=that.fMaxChi2PerClusterTPC;
  fMinChi2PerClusterTPC=that.fMinChi2PerClusterTPC;
  fCutNClustersTPC=that.fCutNClustersTPC;
  fNClustersTPCMax=that.fNClustersTPCMax;
  fNClustersTPCMin=that.fNClustersTPCMin;  
  fCutNClustersITS=that.fCutNClustersITS;
  fNClustersITSMax=that.fNClustersITSMax;
  fNClustersITSMin=that.fNClustersITSMin;  
  fUseAODFilterBit=that.fUseAODFilterBit;
  fAODFilterBit=that.fAODFilterBit;
  fCutDCAToVertexXY=that.fCutDCAToVertexXY;
  fCutDCAToVertexZ=that.fCutDCAToVertexZ;
  fCutMinimalTPCdedx=that.fCutMinimalTPCdedx;
  fMinimalTPCdedx=that.fMinimalTPCdedx;
  fParamType=that.fParamType;
  fParamMix=that.fParamMix;

  fTrack=NULL;
  fTrackPhi=0.;
  fTrackPhi=0.;
  fTrackWeight=0.;
  fTrackLabel=INT_MIN;
  fMCevent=NULL;
  fMCparticle=NULL;
  fEvent=NULL;

  fESDpid = that.fESDpid;
  fPIDsource = that.fPIDsource;

  delete fTPCpidCuts;
  delete fTOFpidCuts;
  if (that.fTPCpidCuts) fTPCpidCuts = new TMatrixF(*(that.fTPCpidCuts));
  if (that.fTOFpidCuts) fTOFpidCuts = new TMatrixF(*(that.fTOFpidCuts));

  fParticleID=that.fParticleID;
  fParticleProbability=that.fParticleProbability;
  memcpy(fProbBayes,that.fProbBayes,sizeof(fProbBayes));

  return *this;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::~AliFlowTrackCuts()
{
  //dtor
  delete fAliESDtrackCuts;
  delete fTPCpidCuts;
  delete fTOFpidCuts;
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetEvent(AliVEvent* event, AliMCEvent* mcEvent)
{
  //set the event
  Clear();
  fEvent=event;
  fMCevent=mcEvent;

  //do the magic for ESD
  AliESDEvent* myESD = dynamic_cast<AliESDEvent*>(event);
  if (fCutPID && myESD)
  {
    //TODO: maybe call it only for the TOF options?
    // Added by F. Noferini for TOF PID
    fESDpid.SetTOFResponse(myESD,AliESDpid::kTOF_T0);
    fESDpid.MakePID(myESD,kFALSE);
    // End F. Noferini added part
  }

  //TODO: AOD
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetCutMC( Bool_t b )
{
  //will we be cutting on MC information?
  fCutMC=b;

  //if we cut on MC info then also the Bethe Bloch should be the one tuned for MC
  if (fCutMC)
  {
    fESDpid.GetTPCResponse().SetBetheBlochParameters( 2.15898e+00/50.,
                                                       1.75295e+01,
                                                       3.40030e-09,
                                                       1.96178e+00,
                                                       3.91720e+00);
  }
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
Bool_t AliFlowTrackCuts::IsSelectedMCtruth(TObject* obj, Int_t id)
{
  //check cuts
  AliVParticle* vparticle = dynamic_cast<AliVParticle*>(obj);
  if (vparticle) 
  {
    return PassesMCcuts(fMCevent,vparticle->GetLabel());
  }
  AliMultiplicity* tracklets = dynamic_cast<AliMultiplicity*>(obj);
  if (tracklets)
  {
    Int_t label0 = tracklets->GetLabel(id,0);
    Int_t label1 = tracklets->GetLabel(id,1);
    Int_t label = (label0==label1)?tracklets->GetLabel(id,1):-666;
    return PassesMCcuts(fMCevent,label);
  }
  return kFALSE;  //default when passed wrong type of object
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliFlowTrackSimple* track)
{
  //check cuts on a flowtracksimple

  //clean up from last iteration
  fTrack = NULL;
  return AliFlowTrackSimpleCuts::PassesCuts(track);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliMultiplicity* tracklet, Int_t id)
{
  //check cuts on a tracklets

  //clean up from last iteration, and init label
  fTrack = NULL;
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
  //if possible get label and mcparticle
  fTrackLabel = (label0==label1)?tracklet->GetLabel(id,1):-1;
  if (!fFakesAreOK && fTrackLabel<0) return kFALSE;
  if (fTrackLabel>=0 && fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  //check MC cuts
  if (fCutMC && !PassesMCcuts()) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts(AliMCEvent* mcEvent, Int_t label)
{
  //check the MC info
  if (!mcEvent) return kFALSE;
  if (label<0) return kFALSE;//otherwise AliCMevent prints a warning before returning NULL
  AliMCParticle* mcparticle = static_cast<AliMCParticle*>(mcEvent->GetTrack(label));
  if (!mcparticle) {AliError("no MC track"); return kFALSE;}

  if (fCutMCisPrimary)
  {
    if (IsPhysicalPrimary(mcEvent,label,fRequireTransportBitForPrimaries) != fMCisPrimary) return kFALSE;
  }
  if (fCutMCPID)
  {
    Int_t pdgCode = mcparticle->PdgCode();
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
    TParticle* particle = mcparticle->Particle();
    Int_t processID = particle->GetUniqueID();
    if (processID != fMCprocessType ) return kFALSE;
  }
  return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts()
{
  if (!fMCevent) return kFALSE;
  if (fTrackLabel<0) return kFALSE;//otherwise AliCMevent prints a warning before returning NULL
  fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  return PassesMCcuts(fMCevent,fTrackLabel);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliVParticle* vparticle)
{
  //check cuts for an ESD vparticle

  ////////////////////////////////////////////////////////////////
  //  start by preparing the track parameters to cut on //////////
  ////////////////////////////////////////////////////////////////
  //clean up from last iteration
  fTrack=NULL; 

  //get the label and the mc particle
  fTrackLabel = (fFakesAreOK)?TMath::Abs(vparticle->GetLabel()):vparticle->GetLabel();
  if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  else fMCparticle=NULL;

  Bool_t isMCparticle = kFALSE; //some things are different for MC particles, check!
  AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*>(vparticle);
  AliAODTrack* aodTrack = NULL;
  if (esdTrack)
    //for an ESD track we do some magic sometimes like constructing TPC only parameters
    //or doing some hybrid, handle that here
    HandleESDtrack(esdTrack);
  else
  {
    HandleVParticle(vparticle);
    //now check if produced particle is MC
    isMCparticle = (dynamic_cast<AliMCParticle*>(fTrack))!=NULL;
    aodTrack = dynamic_cast<AliAODTrack*>(vparticle); //keep the additional dynamic cast out of the way for ESDs
  }
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  if (!fTrack) return kFALSE;
  //because it may be different from global, not needed for aodTrack because we dont do anything funky there
  if (esdTrack) esdTrack = static_cast<AliESDtrack*>(fTrack);
  
  Bool_t pass=kTRUE;
  //check the common cuts for the current particle fTrack (MC,AOD,ESD)
  Double_t pt = fTrack->Pt();
  if (!fFakesAreOK) {if (fTrackLabel<0) pass=kFALSE;}
  if (fCutPt) {if (pt < fPtMin || pt >= fPtMax ) pass=kFALSE;}
  if (fCutEta) {if (fTrack->Eta() < fEtaMin || fTrack->Eta() >= fEtaMax ) pass=kFALSE;}
  if (fCutPhi) {if (fTrack->Phi() < fPhiMin || fTrack->Phi() >= fPhiMax ) pass=kFALSE;}
  if (fRequireCharge) {if (fTrack->Charge() == 0) pass=kFALSE;}
  if (fCutCharge && !isMCparticle) {if (fTrack->Charge() != fCharge) pass=kFALSE;}
  if (fCutCharge && isMCparticle)
  { 
    //in case of an MC particle the charge is stored in units of 1/3|e| 
    Int_t charge = TMath::Nint(fTrack->Charge()/3.0); //mc particles have charge in units of 1/3e
    if (charge!=fCharge) pass=kFALSE;
  }
  //if(fCutPID) {if (fTrack->PID() != fPID) pass=kFALSE;}

  //when additionally MC info is required
  if (fCutMC && !PassesMCcuts()) pass=kFALSE;

  //the case of ESD or AOD
  if (esdTrack) { if (!PassesESDcuts(esdTrack)) { pass=kFALSE; } }
  if (aodTrack) { if (!PassesAODcuts(aodTrack)) { pass=kFALSE; } }

  //true by default, if we didn't set any cuts
  return pass;
}

//_______________________________________________________________________
Bool_t AliFlowTrackCuts::PassesAODcuts(AliAODTrack* track)
{
  Bool_t pass = kTRUE;

  if (fCutNClustersTPC)
  {
    Int_t ntpccls = track->GetTPCNcls();
    if (ntpccls < fNClustersTPCMin || ntpccls > fNClustersTPCMax) pass=kFALSE;
  }

  if (fCutNClustersITS)
  {
    Int_t nitscls = track->GetITSNcls();
    if (nitscls < fNClustersITSMin || nitscls > fNClustersITSMax) pass=kFALSE;
  }
  
   if (fCutChi2PerClusterTPC)
  {
    Double_t chi2tpc = track->Chi2perNDF();
    if (chi2tpc < fMinChi2PerClusterTPC || chi2tpc > fMaxChi2PerClusterTPC) pass=kFALSE;
  }
  
  if (GetRequireTPCRefit() && !(track->GetStatus() & AliESDtrack::kTPCrefit) ) pass=kFALSE;
  if (GetRequireITSRefit() && !(track->GetStatus() & AliESDtrack::kITSrefit) ) pass=kFALSE;
  
  if (fUseAODFilterBit && !track->TestFilterBit(fAODFilterBit)) pass=kFALSE;
  
  if (fCutDCAToVertexXY && track->DCA()>GetMaxDCAToVertexXY()) pass=kFALSE;
    

  return pass;
}

//_______________________________________________________________________
Bool_t AliFlowTrackCuts::PassesESDcuts(AliESDtrack* track)
{
  Bool_t pass=kTRUE;
  if (fIgnoreTPCzRange)
  {
    const AliExternalTrackParam* pin = track->GetOuterParam();
    const AliExternalTrackParam* pout = track->GetInnerParam();
    if (pin&&pout)
    {
      Double_t zin = pin->GetZ();
      Double_t zout = pout->GetZ();
      if (zin*zout<0) pass=kFALSE;   //reject if cross the membrane
      if (zin < fIgnoreTPCzRangeMin || zin > fIgnoreTPCzRangeMax) pass=kFALSE;
      if (zout < fIgnoreTPCzRangeMin || zout > fIgnoreTPCzRangeMax) pass=kFALSE;
    }
  }
 
  Int_t ntpccls = ( fParamType==kESD_TPConly )?
                    track->GetTPCNclsIter1():track->GetTPCNcls();    
  if (fCutChi2PerClusterTPC)
  {
    Float_t tpcchi2 = (fParamType==kESD_TPConly)?
                       track->GetTPCchi2Iter1():track->GetTPCchi2();
    tpcchi2 = (ntpccls>0)?tpcchi2/ntpccls:-FLT_MAX;
    if (tpcchi2<fMinChi2PerClusterTPC || tpcchi2 >=fMaxChi2PerClusterTPC)
      pass=kFALSE;
  }

  if (fCutMinimalTPCdedx) 
  {
    if (track->GetTPCsignal() < fMinimalTPCdedx) pass=kFALSE;
  }

  if (fCutNClustersTPC)
  {
    if (ntpccls < fNClustersTPCMin || ntpccls > fNClustersTPCMax) pass=kFALSE;
  }

  Int_t nitscls = track->GetNcls(0);
  if (fCutNClustersITS)
  {
    if (nitscls < fNClustersITSMin || nitscls > fNClustersITSMax) pass=kFALSE;
  }

  if (fCutPID)
  {
    switch (fPIDsource)    
    {
      case kTPCpid:
        if (!PassesTPCpidCut(track)) pass=kFALSE;
        break;
      case kTPCdedx:
        if (!PassesTPCdedxCut(track)) pass=kFALSE;
        break;
      case kTOFpid:
        if (!PassesTOFpidCut(track)) pass=kFALSE;
        break;
      case kTOFbeta:
        if (!PassesTOFbetaCut(track)) pass=kFALSE;
        break;
	    // part added by F. Noferini
      case kTOFbayesian:
	      if (!PassesTOFbayesianCut(track)) pass=kFALSE;
	      break;
	    // end part added by F. Noferini
      default:
        printf("AliFlowTrackCuts::PassesCuts() this should never be called!\n");
        pass=kFALSE;
        break;
    }
  }    

  //some stuff is still handled by AliESDtrackCuts class - delegate
  if (fAliESDtrackCuts)
  {
    if (!fAliESDtrackCuts->IsSelected(track)) pass=kFALSE;
  }
 
  return pass;
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::HandleVParticle(AliVParticle* track)
{
  //handle the general case
  switch (fParamType)
  {
    default:
      fTrack = track;
      break;
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
      break;
    case kESD_TPConly:
      if (!track->FillTPCOnlyTrack(fTPCtrack)) 
      {
        fTrack=NULL;
        fMCparticle=NULL;
        fTrackLabel=-1;
        return;
      }
      fTrack = &fTPCtrack;
      //recalculate the label and mc particle, they may differ as TPClabel != global label
      fTrackLabel = (fFakesAreOK)?TMath::Abs(fTrack->GetLabel()):fTrack->GetLabel();
      if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
      else fMCparticle=NULL;
      break;
    default:
      fTrack = track;
      break;
  }
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::Count(AliVEvent* event)
{
  //calculate the number of track in given event.
  //if argument is NULL(default) take the event attached
  //by SetEvent()
  Int_t multiplicity = 0;
  if (!event)
  {
    for (Int_t i=0; i<GetNumberOfInputObjects(); i++)
    {
      if (IsSelected(GetInputObject(i))) multiplicity++;
    }
  }
  else
  {
    for (Int_t i=0; i<event->GetNumberOfTracks(); i++)
    {
      if (IsSelected(event->GetTrack(i))) multiplicity++;
    }
  }
  return multiplicity;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardGlobalTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard Global flow cuts");
  cuts->SetParamType(kGlobal);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.1);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMinNClustersITS(2);
  cuts->SetRequireITSRefit(kTRUE);
  cuts->SetRequireTPCRefit(kTRUE);
  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(0.3);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCOnlyTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPConly flow cuts");
  cuts->SetParamType(kESD_TPConly);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMaxDCAToVertexXY(3.0);
  cuts->SetMaxDCAToVertexZ(3.0);
  cuts->SetDCAToVertex2D(kTRUE);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCOnlyTrackCuts()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPConly flow cuts");
  cuts->SetParamType(kESD_TPConly);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMaxDCAToVertexXY(3.0);
  cuts->SetMaxDCAToVertexZ(3.0);
  cuts->SetDCAToVertex2D(kTRUE);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries)
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard global track cuts 2009");
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
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  if (fParamType==kESD_SPDtracklet)
  {
    switch (fParamMix)
    {
      case kPure:
        flowtrack = new AliFlowTrack();
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        break;
      case kTrackWithMCkine:
        if (!fMCparticle) return NULL;
        flowtrack = new AliFlowTrack();
        flowtrack->SetPhi( fMCparticle->Phi() );
        flowtrack->SetEta( fMCparticle->Eta() );
        flowtrack->SetPt( fMCparticle->Pt() );
        break;
      case kTrackWithMCpt:
        if (!fMCparticle) return NULL;
        flowtrack = new AliFlowTrack();
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        flowtrack->SetPt(fMCparticle->Pt());
        break;
      case kTrackWithPtFromFirstMother:
        if (!fMCparticle) return NULL;
        flowtrack = new AliFlowTrack();
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        tmpTParticle = fMCparticle->Particle();
        tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
        flowtrack->SetPt(tmpAliMCParticle->Pt());
        break;
      default:
        flowtrack = new AliFlowTrack();
        flowtrack->SetPhi(fTrackPhi);
        flowtrack->SetEta(fTrackEta);
        break;
    }
    flowtrack->SetSource(AliFlowTrack::kFromTracklet);
  }
  else
  {
    if (!fTrack) return NULL;
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
        break;
      case kTrackWithPtFromFirstMother:
        if (!fMCparticle) return NULL;
        flowtrack = new AliFlowTrack(fTrack);
        tmpTParticle = fMCparticle->Particle();
        tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
        flowtrack->SetPt(tmpAliMCParticle->Pt());
        break;
      default:
        flowtrack = new AliFlowTrack(fTrack);
        break;
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
  return IsPhysicalPrimary(fMCevent, fTrackLabel, fRequireTransportBitForPrimaries);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsPhysicalPrimary(AliMCEvent* mcEvent, Int_t label, Bool_t requiretransported)
{
  //check if current particle is a physical primary
  Bool_t physprim=mcEvent->IsPhysicalPrimary(label);
  AliMCParticle* track = static_cast<AliMCParticle*>(mcEvent->GetTrack(label));
  if (!track) return kFALSE;
  TParticle* particle = track->Particle();
  Bool_t transported = particle->TestBit(kTransportBit);
  //printf("label: %i prim: %s, transp: %s, pass: %s\n",label, (physprim)?"YES":"NO ",(transported)?"YES":"NO ",
        //(physprim && (transported || !requiretransported))?"YES":"NO"  );
  return (physprim && (transported || !requiretransported));
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::DefineHistograms()
{
  //define qa histograms
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

//-----------------------------------------------------------------------
void AliFlowTrackCuts::Clear(Option_t*)
{
  //clean up
  fTrack=NULL;
  fMCevent=NULL;
  fMCparticle=NULL;
  fTrackLabel=0;
  fTrackWeight=0.0;
  fTrackEta=0.0;
  fTrackPhi=0.0;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaCut(AliESDtrack* track )
{
  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track) && 
                     (track->GetStatus() & AliESDtrack::kTOFpid) && 
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365) && 
                    !(track->GetStatus() & AliESDtrack::kTOFmismatch);

  if (!goodtrack) return kFALSE;
  
  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t L = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  Float_t beta = L/timeTOF/c;
  Double_t integratedTimes[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
  track->GetIntegratedTimes(integratedTimes);

  //construct the pid index because it's not AliPID::EParticleType
  Int_t pid = 0;
  switch (fParticleID)
  {
    case AliPID::kPion:
      pid=2;
      break;
    case AliPID::kKaon:
      pid=3;
      break;
    case AliPID::kProton:
      pid=4;
      break;
    default:
      return kFALSE;
  }

  //signal to cut on
  Float_t s = beta-L/integratedTimes[pid]/c;

  Float_t* arr = fTOFpidCuts->GetMatrixArray();
  Int_t col = TMath::BinarySearch(fTOFpidCuts->GetNcols(),arr,static_cast<Float_t>(p));
  if (col<0) return kFALSE;
  Float_t min = (*fTOFpidCuts)(1,col);
  Float_t max = (*fTOFpidCuts)(2,col);

  //printf("--------------TOF beta cut %s\n",(s>min && s<max)?"PASS":"FAIL");
  return (s>min && s<max);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFpidCut(AliESDtrack* track)
{
  //check if passes PID cut using default TOF pid
  Double_t pidTOF[AliPID::kSPECIES];
  track->GetTOFpid(pidTOF);
  if (pidTOF[fParticleID]>=fParticleProbability) return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCpidCut(AliESDtrack* track)
{
  //check if passes PID cut using default TPC pid
  Double_t pidTPC[AliPID::kSPECIES];
  track->GetTPCpid(pidTPC);
  Double_t probablity = 0.;
  switch (fParticleID)
  {
    case AliPID::kPion:
      probablity = pidTPC[AliPID::kPion] + pidTPC[AliPID::kMuon];
      break;
    default:
      probablity = pidTPC[fParticleID];
  }
  if (probablity >= fParticleProbability) return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCdedxCut(AliESDtrack* track)
{
  //check if passes PID cut using dedx signal in the TPC
  if (!fTPCpidCuts)
  {
    printf("no TPCpidCuts\n");
    return kFALSE;
  }

  const AliExternalTrackParam* tpcparam = track->GetInnerParam(); //tpc only params at the inner wall
  if (!tpcparam) return kFALSE;
  Float_t sigExp = fESDpid.GetTPCResponse().GetExpectedSignal(tpcparam->GetP(), fParticleID);
  Float_t sigTPC = track->GetTPCsignal();
  Float_t s = (sigTPC-sigExp)/sigExp;
  Double_t pt = track->Pt();

  Float_t* arr = fTPCpidCuts->GetMatrixArray();
  Int_t col = TMath::BinarySearch(fTPCpidCuts->GetNcols(),arr,static_cast<Float_t>(pt));
  if (col<0) return kFALSE;
  Float_t min = (*fTPCpidCuts)(1,col);
  Float_t max = (*fTPCpidCuts)(2,col);

  //printf("------------TPC pid cut %s\n",(s>min && s<max)?"PASS":"FAIL");
  return (s>min && s<max);
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::InitPIDcuts()
{
  //init matrices with PID cuts
  TMatrixF* t = NULL;
  if (!fTPCpidCuts)
  {
    if (fParticleID==AliPID::kPion)
    {
      t = new TMatrixF(3,10);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.4;  (*t)(2,0)  =   0.2;
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.4;  (*t)(2,1)  =   0.2;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.4;  (*t)(2,2)  =  0.25;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.4;  (*t)(2,3)  =  0.25;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.4;  (*t)(2,4)  =   0.3;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.4;  (*t)(2,5)  =   0.3;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  = -0.4;  (*t)(2,6)  =  0.25;
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.4;  (*t)(2,7)  =  0.15;
      (*t)(0,8)  = 0.60;  (*t)(1,8)  = -0.4;  (*t)(2,8)  =   0.1;
      (*t)(0,9)  = 0.65;  (*t)(1,9)  =    0;  (*t)(2,9)  =     0;
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      t = new TMatrixF(3,7);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.2;  (*t)(2,0)  = 0.4; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  =-0.15;  (*t)(2,1)  = 0.4;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.1;  (*t)(2,2)  = 0.4;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.1;  (*t)(2,3)  = 0.4;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.1;  (*t)(2,4)  = 0.6;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.1;  (*t)(2,5)  = 0.6;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  =    0;  (*t)(2,6)  =   0;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      t = new TMatrixF(3,16);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  =     0;  (*t)(2,0)  =    0; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  =  -0.2;  (*t)(2,1)  =  0.3; 
      (*t)(0,2)  = 0.30;  (*t)(1,2)  =  -0.2;  (*t)(2,2)  =  0.6; 
      (*t)(0,3)  = 0.35;  (*t)(1,3)  =  -0.2;  (*t)(2,3)  =  0.6; 
      (*t)(0,4)  = 0.40;  (*t)(1,4)  =  -0.2;  (*t)(2,4)  =  0.6; 
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.15;  (*t)(2,5)  =  0.6; 
      (*t)(0,6)  = 0.50;  (*t)(1,6)  =  -0.1;  (*t)(2,6)  =  0.6; 
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.05;  (*t)(2,7)  =  0.6; 
      (*t)(0,8)  = 0.60;  (*t)(1,8)  = -0.05;  (*t)(2,8)  = 0.45; 
      (*t)(0,9)  = 0.65;  (*t)(1,9)  = -0.05;  (*t)(2,9)  = 0.45; 
      (*t)(0,10) = 0.70;  (*t)(1,10) = -0.05;  (*t)(2,10) = 0.45; 
      (*t)(0,11) = 0.75;  (*t)(1,11) = -0.05;  (*t)(2,11) = 0.45; 
      (*t)(0,12) = 0.80;  (*t)(1,12) =     0;  (*t)(2,12) = 0.45; 
      (*t)(0,13) = 0.85;  (*t)(1,13) =     0;  (*t)(2,13) = 0.45; 
      (*t)(0,14) = 0.90;  (*t)(1,14) =     0;  (*t)(2,14) = 0.45;
      (*t)(0,15) = 0.95;  (*t)(1,15) =     0;  (*t)(2,15) =    0;
    }
    fTPCpidCuts=t;
  }
  t = NULL;
  if (!fTOFpidCuts)
  {
    if (fParticleID==AliPID::kPion)
    {
      //TOF pions, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(2,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(2,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(2,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(2,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(2,4)  = 0.000;  (*t)(2,4)  =   0.000;
      (*t)(0,5)  = 0.250;  (*t)(2,5)  = -0.046;  (*t)(2,5)  =   0.046;
      (*t)(0,6)  = 0.300;  (*t)(2,6)  = -0.038;  (*t)(2,6)  =   0.038;
      (*t)(0,7)  = 0.350;  (*t)(2,7)  = -0.034;  (*t)(2,7)  =   0.034;
      (*t)(0,8)  = 0.400;  (*t)(2,8)  = -0.032;  (*t)(2,8)  =   0.032;
      (*t)(0,9)  = 0.450;  (*t)(2,9)  = -0.030;  (*t)(2,9)  =   0.030;
      (*t)(0,10)  = 0.500;  (*t)(2,10)  = -0.030;  (*t)(2,10)  =   0.030;
      (*t)(0,11)  = 0.550;  (*t)(2,11)  = -0.030;  (*t)(2,11)  =   0.030;
      (*t)(0,12)  = 0.600;  (*t)(2,12)  = -0.030;  (*t)(2,12)  =   0.030;
      (*t)(0,13)  = 0.650;  (*t)(2,13)  = -0.030;  (*t)(2,13)  =   0.030;
      (*t)(0,14)  = 0.700;  (*t)(2,14)  = -0.030;  (*t)(2,14)  =   0.030;
      (*t)(0,15)  = 0.750;  (*t)(2,15)  = -0.030;  (*t)(2,15)  =   0.030;
      (*t)(0,16)  = 0.800;  (*t)(2,16)  = -0.030;  (*t)(2,16)  =   0.030;
      (*t)(0,17)  = 0.850;  (*t)(2,17)  = -0.030;  (*t)(2,17)  =   0.030;
      (*t)(0,18)  = 0.900;  (*t)(2,18)  = -0.030;  (*t)(2,18)  =   0.030;
      (*t)(0,19)  = 0.950;  (*t)(2,19)  = -0.028;  (*t)(2,19)  =   0.028;
      (*t)(0,20)  = 1.000;  (*t)(2,20)  = -0.028;  (*t)(2,20)  =   0.028;
      (*t)(0,21)  = 1.100;  (*t)(2,21)  = -0.028;  (*t)(2,21)  =   0.028;
      (*t)(0,22)  = 1.200;  (*t)(2,22)  = -0.026;  (*t)(2,22)  =   0.028;
      (*t)(0,23)  = 1.300;  (*t)(2,23)  = -0.024;  (*t)(2,23)  =   0.028;
      (*t)(0,24)  = 1.400;  (*t)(2,24)  = -0.020;  (*t)(2,24)  =   0.028;
      (*t)(0,25)  = 1.500;  (*t)(2,25)  = -0.018;  (*t)(2,25)  =   0.028;
      (*t)(0,26)  = 1.600;  (*t)(2,26)  = -0.016;  (*t)(2,26)  =   0.028;
      (*t)(0,27)  = 1.700;  (*t)(2,27)  = -0.014;  (*t)(2,27)  =   0.028;
      (*t)(0,28)  = 1.800;  (*t)(2,28)  = -0.012;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(2,29)  = -0.010;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(2,30)  = -0.008;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(2,31)  = -0.008;  (*t)(2,31)  =   0.024;
      (*t)(0,32)  = 2.200;  (*t)(2,32)  = -0.006;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(2,33)  = -0.004;  (*t)(2,33)  =   0.024;
      (*t)(0,34)  = 2.400;  (*t)(2,34)  = -0.004;  (*t)(2,34)  =   0.024;
      (*t)(0,35)  = 2.500;  (*t)(2,35)  = -0.002;  (*t)(2,35)  =   0.024;
      (*t)(0,36)  = 2.600;  (*t)(2,36)  = -0.002;  (*t)(2,36)  =   0.024;
      (*t)(0,37)  = 2.700;  (*t)(2,37)  = 0.000;  (*t)(2,37)  =   0.024;
      (*t)(0,38)  = 2.800;  (*t)(2,38)  = 0.000;  (*t)(2,38)  =   0.026;
      (*t)(0,39)  = 2.900;  (*t)(2,39)  = 0.000;  (*t)(2,39)  =   0.024;
      (*t)(0,40)  = 3.000;  (*t)(2,40)  = 0.002;  (*t)(2,40)  =   0.026;
      (*t)(0,41)  = 3.100;  (*t)(2,41)  = 0.002;  (*t)(2,41)  =   0.026;
      (*t)(0,42)  = 3.200;  (*t)(2,42)  = 0.002;  (*t)(2,42)  =   0.026;
      (*t)(0,43)  = 3.300;  (*t)(2,43)  = 0.002;  (*t)(2,43)  =   0.026;
      (*t)(0,44)  = 3.400;  (*t)(2,44)  = 0.002;  (*t)(2,44)  =   0.026;
      (*t)(0,45)  = 3.500;  (*t)(2,45)  = 0.002;  (*t)(2,45)  =   0.026;
      (*t)(0,46)  = 3.600;  (*t)(2,46)  = 0.002;  (*t)(2,46)  =   0.026;
      (*t)(0,47)  = 3.700;  (*t)(2,47)  = 0.002;  (*t)(2,47)  =   0.026;
      (*t)(0,48)  = 3.800;  (*t)(2,48)  = 0.002;  (*t)(2,48)  =   0.026;
      (*t)(0,49)  = 3.900;  (*t)(2,49)  = 0.004;  (*t)(2,49)  =   0.024;
      (*t)(0,50)  = 4.000;  (*t)(2,50)  = 0.004;  (*t)(2,50)  =   0.026;
      (*t)(0,51)  = 4.100;  (*t)(2,51)  = 0.004;  (*t)(2,51)  =   0.026;
      (*t)(0,52)  = 4.200;  (*t)(2,52)  = 0.004;  (*t)(2,52)  =   0.024;
      (*t)(0,53)  = 4.300;  (*t)(2,53)  = 0.006;  (*t)(2,53)  =   0.024;
      (*t)(0,54)  = 4.400;  (*t)(2,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(2,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(2,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(2,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(2,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(2,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(2,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      //TOF protons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(2,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(2,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(2,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(2,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(2,4)  = 0.000;  (*t)(2,4)  =   0.000;
      (*t)(0,5)  = 0.250;  (*t)(2,5)  = 0.000;  (*t)(2,5)  =   0.000;
      (*t)(0,6)  = 0.300;  (*t)(2,6)  = 0.000;  (*t)(2,6)  =   0.000;
      (*t)(0,7)  = 0.350;  (*t)(2,7)  = 0.000;  (*t)(2,7)  =   0.000;
      (*t)(0,8)  = 0.400;  (*t)(2,8)  = 0.000;  (*t)(2,8)  =   0.000;
      (*t)(0,9)  = 0.450;  (*t)(2,9)  = 0.000;  (*t)(2,9)  =   0.000;
      (*t)(0,10)  = 0.500;  (*t)(2,10)  = 0.000;  (*t)(2,10)  =   0.000;
      (*t)(0,11)  = 0.550;  (*t)(2,11)  = 0.000;  (*t)(2,11)  =   0.000;
      (*t)(0,12)  = 0.600;  (*t)(2,12)  = 0.000;  (*t)(2,12)  =   0.000;
      (*t)(0,13)  = 0.650;  (*t)(2,13)  = 0.000;  (*t)(2,13)  =   0.000;
      (*t)(0,14)  = 0.700;  (*t)(2,14)  = 0.000;  (*t)(2,14)  =   0.000;
      (*t)(0,15)  = 0.750;  (*t)(2,15)  = 0.000;  (*t)(2,15)  =   0.000;
      (*t)(0,16)  = 0.800;  (*t)(2,16)  = 0.000;  (*t)(2,16)  =   0.000;
      (*t)(0,17)  = 0.850;  (*t)(2,17)  = -0.070;  (*t)(2,17)  =   0.070;
      (*t)(0,18)  = 0.900;  (*t)(2,18)  = -0.072;  (*t)(2,18)  =   0.072;
      (*t)(0,19)  = 0.950;  (*t)(2,19)  = -0.072;  (*t)(2,19)  =   0.072;
      (*t)(0,20)  = 1.000;  (*t)(2,20)  = -0.074;  (*t)(2,20)  =   0.074;
      (*t)(0,21)  = 1.100;  (*t)(2,21)  = -0.032;  (*t)(2,21)  =   0.032;
      (*t)(0,22)  = 1.200;  (*t)(2,22)  = -0.026;  (*t)(2,22)  =   0.026;
      (*t)(0,23)  = 1.300;  (*t)(2,23)  = -0.026;  (*t)(2,23)  =   0.026;
      (*t)(0,24)  = 1.400;  (*t)(2,24)  = -0.024;  (*t)(2,24)  =   0.024;
      (*t)(0,25)  = 1.500;  (*t)(2,25)  = -0.024;  (*t)(2,25)  =   0.024;
      (*t)(0,26)  = 1.600;  (*t)(2,26)  = -0.026;  (*t)(2,26)  =   0.026;
      (*t)(0,27)  = 1.700;  (*t)(2,27)  = -0.026;  (*t)(2,27)  =   0.026;
      (*t)(0,28)  = 1.800;  (*t)(2,28)  = -0.026;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(2,29)  = -0.026;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(2,30)  = -0.026;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(2,31)  = -0.026;  (*t)(2,31)  =   0.026;
      (*t)(0,32)  = 2.200;  (*t)(2,32)  = -0.026;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(2,33)  = -0.028;  (*t)(2,33)  =   0.022;
      (*t)(0,34)  = 2.400;  (*t)(2,34)  = -0.028;  (*t)(2,34)  =   0.020;
      (*t)(0,35)  = 2.500;  (*t)(2,35)  = -0.028;  (*t)(2,35)  =   0.018;
      (*t)(0,36)  = 2.600;  (*t)(2,36)  = -0.028;  (*t)(2,36)  =   0.016;
      (*t)(0,37)  = 2.700;  (*t)(2,37)  = -0.028;  (*t)(2,37)  =   0.016;
      (*t)(0,38)  = 2.800;  (*t)(2,38)  = -0.030;  (*t)(2,38)  =   0.014;
      (*t)(0,39)  = 2.900;  (*t)(2,39)  = -0.030;  (*t)(2,39)  =   0.012;
      (*t)(0,40)  = 3.000;  (*t)(2,40)  = -0.030;  (*t)(2,40)  =   0.012;
      (*t)(0,41)  = 3.100;  (*t)(2,41)  = -0.030;  (*t)(2,41)  =   0.010;
      (*t)(0,42)  = 3.200;  (*t)(2,42)  = -0.030;  (*t)(2,42)  =   0.010;
      (*t)(0,43)  = 3.300;  (*t)(2,43)  = -0.030;  (*t)(2,43)  =   0.010;
      (*t)(0,44)  = 3.400;  (*t)(2,44)  = -0.030;  (*t)(2,44)  =   0.008;
      (*t)(0,45)  = 3.500;  (*t)(2,45)  = -0.030;  (*t)(2,45)  =   0.008;
      (*t)(0,46)  = 3.600;  (*t)(2,46)  = -0.030;  (*t)(2,46)  =   0.008;
      (*t)(0,47)  = 3.700;  (*t)(2,47)  = -0.030;  (*t)(2,47)  =   0.006;
      (*t)(0,48)  = 3.800;  (*t)(2,48)  = -0.030;  (*t)(2,48)  =   0.006;
      (*t)(0,49)  = 3.900;  (*t)(2,49)  = -0.030;  (*t)(2,49)  =   0.006;
      (*t)(0,50)  = 4.000;  (*t)(2,50)  = -0.028;  (*t)(2,50)  =   0.004;
      (*t)(0,51)  = 4.100;  (*t)(2,51)  = -0.030;  (*t)(2,51)  =   0.004;
      (*t)(0,52)  = 4.200;  (*t)(2,52)  = -0.030;  (*t)(2,52)  =   0.004;
      (*t)(0,53)  = 4.300;  (*t)(2,53)  = -0.028;  (*t)(2,53)  =   0.002;
      (*t)(0,54)  = 4.400;  (*t)(2,54)  = -0.030;  (*t)(2,54)  =   0.002;
      (*t)(0,55)  = 4.500;  (*t)(2,55)  = -0.028;  (*t)(2,55)  =   0.002;
      (*t)(0,56)  = 4.600;  (*t)(2,56)  = -0.028;  (*t)(2,56)  =   0.002;
      (*t)(0,57)  = 4.700;  (*t)(2,57)  = -0.028;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(2,58)  = -0.028;  (*t)(2,58)  =   0.002;
      (*t)(0,59)  = 4.900;  (*t)(2,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(2,60)  = 0.000;  (*t)(2,60)  =   0.000; 
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      //TOF kaons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(2,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(2,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(2,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(2,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(2,4)  = 0.000;  (*t)(2,4)  =   0.000;
      (*t)(0,5)  = 0.250;  (*t)(2,5)  = 0.000;  (*t)(2,5)  =   0.000;
      (*t)(0,6)  = 0.300;  (*t)(2,6)  = 0.000;  (*t)(2,6)  =   0.000;
      (*t)(0,7)  = 0.350;  (*t)(2,7)  = 0.000;  (*t)(2,7)  =   0.000;
      (*t)(0,8)  = 0.400;  (*t)(2,8)  = 0.000;  (*t)(2,8)  =   0.000;
      (*t)(0,9)  = 0.450;  (*t)(2,9)  = 0.000;  (*t)(2,9)  =   0.000;
      (*t)(0,10)  = 0.500;  (*t)(2,10)  = 0.000;  (*t)(2,10)  =   0.000;
      (*t)(0,11)  = 0.550;  (*t)(2,11)  = -0.026;  (*t)(2,11)  =   0.026;
      (*t)(0,12)  = 0.600;  (*t)(2,12)  = -0.026;  (*t)(2,12)  =   0.026;
      (*t)(0,13)  = 0.650;  (*t)(2,13)  = -0.026;  (*t)(2,13)  =   0.026;
      (*t)(0,14)  = 0.700;  (*t)(2,14)  = -0.026;  (*t)(2,14)  =   0.026;
      (*t)(0,15)  = 0.750;  (*t)(2,15)  = -0.026;  (*t)(2,15)  =   0.026;
      (*t)(0,16)  = 0.800;  (*t)(2,16)  = -0.026;  (*t)(2,16)  =   0.026;
      (*t)(0,17)  = 0.850;  (*t)(2,17)  = -0.024;  (*t)(2,17)  =   0.024;
      (*t)(0,18)  = 0.900;  (*t)(2,18)  = -0.024;  (*t)(2,18)  =   0.024;
      (*t)(0,19)  = 0.950;  (*t)(2,19)  = -0.024;  (*t)(2,19)  =   0.024;
      (*t)(0,20)  = 1.000;  (*t)(2,20)  = -0.024;  (*t)(2,20)  =   0.024;
      (*t)(0,21)  = 1.100;  (*t)(2,21)  = -0.024;  (*t)(2,21)  =   0.024;
      (*t)(0,22)  = 1.200;  (*t)(2,22)  = -0.024;  (*t)(2,22)  =   0.022;
      (*t)(0,23)  = 1.300;  (*t)(2,23)  = -0.024;  (*t)(2,23)  =   0.020;
      (*t)(0,24)  = 1.400;  (*t)(2,24)  = -0.026;  (*t)(2,24)  =   0.016;
      (*t)(0,25)  = 1.500;  (*t)(2,25)  = -0.028;  (*t)(2,25)  =   0.014;
      (*t)(0,26)  = 1.600;  (*t)(2,26)  = -0.028;  (*t)(2,26)  =   0.012;
      (*t)(0,27)  = 1.700;  (*t)(2,27)  = -0.028;  (*t)(2,27)  =   0.010;
      (*t)(0,28)  = 1.800;  (*t)(2,28)  = -0.028;  (*t)(2,28)  =   0.010;
      (*t)(0,29)  = 1.900;  (*t)(2,29)  = -0.028;  (*t)(2,29)  =   0.008;
      (*t)(0,30)  = 2.000;  (*t)(2,30)  = -0.028;  (*t)(2,30)  =   0.006;
      (*t)(0,31)  = 2.100;  (*t)(2,31)  = -0.026;  (*t)(2,31)  =   0.006;
      (*t)(0,32)  = 2.200;  (*t)(2,32)  = -0.024;  (*t)(2,32)  =   0.004;
      (*t)(0,33)  = 2.300;  (*t)(2,33)  = -0.020;  (*t)(2,33)  =   0.002;
      (*t)(0,34)  = 2.400;  (*t)(2,34)  = -0.020;  (*t)(2,34)  =   0.002;
      (*t)(0,35)  = 2.500;  (*t)(2,35)  = -0.018;  (*t)(2,35)  =   0.000;
      (*t)(0,36)  = 2.600;  (*t)(2,36)  = -0.016;  (*t)(2,36)  =   0.000;
      (*t)(0,37)  = 2.700;  (*t)(2,37)  = -0.014;  (*t)(2,37)  =   -0.002;
      (*t)(0,38)  = 2.800;  (*t)(2,38)  = -0.014;  (*t)(2,38)  =   -0.004;
      (*t)(0,39)  = 2.900;  (*t)(2,39)  = -0.012;  (*t)(2,39)  =   -0.004;
      (*t)(0,40)  = 3.000;  (*t)(2,40)  = -0.010;  (*t)(2,40)  =   -0.006;
      (*t)(0,41)  = 3.100;  (*t)(2,41)  = 0.000;  (*t)(2,41)  =   0.000;
      (*t)(0,42)  = 3.200;  (*t)(2,42)  = 0.000;  (*t)(2,42)  =   0.000;
      (*t)(0,43)  = 3.300;  (*t)(2,43)  = 0.000;  (*t)(2,43)  =   0.000;
      (*t)(0,44)  = 3.400;  (*t)(2,44)  = 0.000;  (*t)(2,44)  =   0.000;
      (*t)(0,45)  = 3.500;  (*t)(2,45)  = 0.000;  (*t)(2,45)  =   0.000;
      (*t)(0,46)  = 3.600;  (*t)(2,46)  = 0.000;  (*t)(2,46)  =   0.000;
      (*t)(0,47)  = 3.700;  (*t)(2,47)  = 0.000;  (*t)(2,47)  =   0.000;
      (*t)(0,48)  = 3.800;  (*t)(2,48)  = 0.000;  (*t)(2,48)  =   0.000;
      (*t)(0,49)  = 3.900;  (*t)(2,49)  = 0.000;  (*t)(2,49)  =   0.000;
      (*t)(0,50)  = 4.000;  (*t)(2,50)  = 0.000;  (*t)(2,50)  =   0.000;
      (*t)(0,51)  = 4.100;  (*t)(2,51)  = 0.000;  (*t)(2,51)  =   0.000;
      (*t)(0,52)  = 4.200;  (*t)(2,52)  = 0.000;  (*t)(2,52)  =   0.000;
      (*t)(0,53)  = 4.300;  (*t)(2,53)  = 0.000;  (*t)(2,53)  =   0.000;
      (*t)(0,54)  = 4.400;  (*t)(2,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(2,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(2,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(2,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(2,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(2,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(2,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    fTOFpidCuts=t;
  }
}

//-----------------------------------------------------------------------
// part added by F. Noferini (some methods)
Bool_t AliFlowTrackCuts::PassesTOFbayesianCut(AliESDtrack* track){

  Bool_t goodtrack = track && (track->GetStatus() & AliESDtrack::kTOFpid) && (track->GetTOFsignal() > 12000) && (track->GetTOFsignal() < 100000) && (track->GetIntegratedLength() > 365) && !(track->GetStatus() & AliESDtrack::kTOFmismatch);

  if (! goodtrack)
       return kFALSE;

  Int_t pdg = GetESDPdg(track,"bayesianALL");
  //  printf("pdg set to %i\n",pdg);

  Int_t pid = 0;
  Float_t prob = 0;
  switch (fParticleID)
  {
    case AliPID::kPion:
      pid=211;
      prob = fProbBayes[2];
      break;
    case AliPID::kKaon:
      pid=321;
      prob = fProbBayes[3];
     break;
    case AliPID::kProton:
      pid=2212;
      prob = fProbBayes[4];
      break;
    case AliPID::kElectron:
      pid=-11;
       prob = fProbBayes[0];
     break;
    default:
      return kFALSE;
  }

  //  printf("pt = %f -- all prob = [%4.2f,%4.2f,%4.2f,%4.2f,%4.2f] -- prob = %f\n",track->Pt(),fProbBayes[0],fProbBayes[1],fProbBayes[2],fProbBayes[3],fProbBayes[4],prob);
  if(TMath::Abs(pdg) == TMath::Abs(pid) && prob > 0.8){
    if(!fCutCharge)
      return kTRUE;
    else if (fCutCharge && fCharge * track->GetSign() > 0)
      return kTRUE;
  }
  return kFALSE;
}
//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::GetESDPdg(AliESDtrack *track,Option_t *option,Int_t ipart,Float_t cPi,Float_t cKa,Float_t cPr){
  Int_t pdg = 0;
  Int_t pdgvalues[5] = {-11,-13,211,321,2212};
  Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};

  if(strstr(option,"bayesianTOF")){ // Bayesian TOF PID
    Double_t c[5]={0.01, 0.01, 0.85, 0.1, 0.05};
    Double_t rcc=0.;
    
    Float_t pt = track->Pt();
    
    Int_t iptesd = 0;
    while(pt > fBinLimitPID[iptesd] && iptesd < fnPIDptBin-1) iptesd++;
  
    if(cPi < 0){
      c[0] = fC[iptesd][0];
      c[1] = fC[iptesd][1];
      c[2] = fC[iptesd][2];
      c[3] = fC[iptesd][3];
      c[4] = fC[iptesd][4];
    }
    else{
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = cPi;
      c[3] = cKa;
      c[4] = cPr;      
    }

    Double_t r1[10]; track->GetTOFpid(r1);
    
    Int_t i;
    for (i=0; i<5; i++) rcc+=(c[i]*r1[i]);
    
    Double_t w[10];
    for (i=0; i<5; i++){
	w[i]=c[i]*r1[i]/rcc;
	fProbBayes[i] = w[i];
    }
    if (w[2]>=w[3] && w[2]>=w[4] && w[2]>=w[1] && w[2]>=w[0]) {//pion
      pdg = 211*Int_t(track->GetSign());
    }
    else if (w[4]>=w[3] && w[4]>=w[1] && w[4]>=w[0]) {//proton
      pdg = 2212*Int_t(track->GetSign());
    }
    else if (w[3]>=w[1] && w[3]>=w[0]){//kaon
      pdg = 321*Int_t(track->GetSign());
    }
    else if (w[0]>=w[1]) { //electrons
      pdg = -11*Int_t(track->GetSign());
    }
    else{ // muon
      pdg = -13*Int_t(track->GetSign());
    }
  }

  else if(strstr(option,"bayesianTPC")){ // Bayesian TPC PID
    Double_t c[5]={0.01, 0.01, 0.85, 0.1, 0.05};
    Double_t rcc=0.;
    
    Float_t pt = track->Pt();
    
    Int_t iptesd = 0;
    while(pt > fBinLimitPID[iptesd] && iptesd < fnPIDptBin-1) iptesd++;
  
    if(cPi < 0){
      c[0] = fC[iptesd][0];
      c[1] = fC[iptesd][1];
      c[2] = fC[iptesd][2];
      c[3] = fC[iptesd][3];
      c[4] = fC[iptesd][4];
    }
    else{
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = cPi;
      c[3] = cKa;
      c[4] = cPr;      
    }

    Double_t r1[10]; track->GetTPCpid(r1);
    
    Int_t i;
    for (i=0; i<5; i++) rcc+=(c[i]*r1[i]);
    
    Double_t w[10];
    for (i=0; i<5; i++){
	w[i]=c[i]*r1[i]/rcc;
    	fProbBayes[i] = w[i];
    }
    if (w[2]>=w[3] && w[2]>=w[4] && w[2]>=w[1] && w[2]>=w[0]) {//pion
      pdg = 211*Int_t(track->GetSign());
    }
    else if (w[4]>=w[3] && w[4]>=w[1] && w[4]>=w[0]) {//proton
      pdg = 2212*Int_t(track->GetSign());
    }
    else if (w[3]>=w[1] && w[3]>=w[0]){//kaon
      pdg = 321*Int_t(track->GetSign());
    }
    else if (w[0]>=w[1]) { //electrons
      pdg = -11*Int_t(track->GetSign());
    }
    else{ // muon
      pdg = -13*Int_t(track->GetSign());
    }
  }
  
  else if(strstr(option,"bayesianALL")){
    Double_t c[5]={0.01, 0.01, 0.85, 0.1, 0.05};
    Double_t rcc=0.;
    
    Float_t pt = track->Pt();
    
    Int_t iptesd = 0;
    while(pt > fBinLimitPID[iptesd] && iptesd < fnPIDptBin-1) iptesd++;

    if(cPi < 0){
      c[0] = fC[iptesd][0];
      c[1] = fC[iptesd][1];
      c[2] = fC[iptesd][2];
      c[3] = fC[iptesd][3];
      c[4] = fC[iptesd][4];
    }
    else{
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = cPi;
      c[3] = cKa;
      c[4] = cPr;      
    }

    Double_t r1[10]; track->GetTOFpid(r1);
    Double_t r2[10]; track->GetTPCpid(r2);

    Int_t i;
    for (i=0; i<5; i++) rcc+=(c[i]*r1[i]*r2[i]);
    

    Double_t w[10];
    for (i=0; i<5; i++){
	w[i]=c[i]*r1[i]*r2[i]/rcc;
    	fProbBayes[i] = w[i];
    }

    if (w[2]>=w[3] && w[2]>=w[4] && w[2]>=w[1] && w[2]>=w[0]) {//pion
      pdg = 211*Int_t(track->GetSign());
    }
    else if (w[4]>=w[3] && w[4]>=w[1] && w[4]>=w[0]) {//proton
      pdg = 2212*Int_t(track->GetSign());
    }
    else if (w[3]>=w[1] && w[3]>=w[0]){//kaon
      pdg = 321*Int_t(track->GetSign());
    }
    else if (w[0]>=w[1]) { //electrons
      pdg = -11*Int_t(track->GetSign());
    }
    else{ // muon
      pdg = -13*Int_t(track->GetSign());
    }
  }

  else if(strstr(option,"sigmacutTOF")){
    printf("PID not implemented yet: %s\nNO PID!!!!\n",option);
    Float_t p = track->P();

    // Take expected times
    Double_t exptimes[5];
    track->GetIntegratedTimes(exptimes);

    // Take resolution for TOF response
    // like fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[ipart], mass[ipart]);
    Float_t resolution = fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[ipart], mass[ipart]);

    if(TMath::Abs(exptimes[ipart] - track->GetTOFsignal()) < 3 * resolution){
      pdg = pdgvalues[ipart] * Int_t(track->GetSign());
    }
  }

  else{
    printf("Invalid PID option: %s\nNO PID!!!!\n",option);
  }

  return pdg;
}
//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetPriors(){
  // set abbundancies
    fBinLimitPID[0] = 0.30;
    fC[0][0] = 0.015;
    fC[0][1] = 0.015;
    fC[0][2] = 1;
    fC[0][3] = 0.0025;
    fC[0][4] = 0.000015;
    fBinLimitPID[1] = 0.35;
    fC[1][0] = 0.015;
    fC[1][1] = 0.015;
    fC[1][2] = 1;
    fC[1][3] = 0.01;
    fC[1][4] = 0.001;
    fBinLimitPID[2] = 0.40;
    fC[2][0] = 0.015;
    fC[2][1] = 0.015;
    fC[2][2] = 1;
    fC[2][3] = 0.026;
    fC[2][4] = 0.004;
    fBinLimitPID[3] = 0.45;
    fC[3][0] = 0.015;
    fC[3][1] = 0.015;
    fC[3][2] = 1;
    fC[3][3] = 0.026;
    fC[3][4] = 0.004;
    fBinLimitPID[4] = 0.50;
    fC[4][0] = 0.015;
    fC[4][1] = 0.015;
    fC[4][2] = 1.000000;
    fC[4][3] = 0.05;
    fC[4][4] = 0.01;
    fBinLimitPID[5] = 0.60;
    fC[5][0] = 0.012;
    fC[5][1] = 0.012;
    fC[5][2] = 1;
    fC[5][3] = 0.085;
    fC[5][4] = 0.022;
    fBinLimitPID[6] = 0.70;
    fC[6][0] = 0.01;
    fC[6][1] = 0.01;
    fC[6][2] = 1;
    fC[6][3] = 0.12;
    fC[6][4] = 0.036;
    fBinLimitPID[7] = 0.80;
    fC[7][0] = 0.0095;
    fC[7][1] = 0.0095;
    fC[7][2] = 1;
    fC[7][3] = 0.15;
    fC[7][4] = 0.05;
    fBinLimitPID[8] = 0.90;
    fC[8][0] = 0.0085;
    fC[8][1] = 0.0085;
    fC[8][2] = 1;
    fC[8][3] = 0.18;
    fC[8][4] = 0.074;
    fBinLimitPID[9] = 1;
    fC[9][0] = 0.008;
    fC[9][1] = 0.008;
    fC[9][2] = 1;
    fC[9][3] = 0.22;
    fC[9][4] = 0.1;
    fBinLimitPID[10] = 1.20;
    fC[10][0] = 0.007;
    fC[10][1] = 0.007;
    fC[10][2] = 1;
    fC[10][3] = 0.28;
    fC[10][4] = 0.16;
    fBinLimitPID[11] = 1.40;
    fC[11][0] = 0.0066;
    fC[11][1] = 0.0066;
    fC[11][2] = 1;
    fC[11][3] = 0.35;
    fC[11][4] = 0.23;
    fBinLimitPID[12] = 1.60;
    fC[12][0] = 0.0075;
    fC[12][1] = 0.0075;
    fC[12][2] = 1;
    fC[12][3] = 0.40;
    fC[12][4] = 0.31;
    fBinLimitPID[13] = 1.80;
    fC[13][0] = 0.0062;
    fC[13][1] = 0.0062;
    fC[13][2] = 1;
    fC[13][3] = 0.45;
    fC[13][4] = 0.39;
    fBinLimitPID[14] = 2.00;
    fC[14][0] = 0.005;
    fC[14][1] = 0.005;
    fC[14][2] = 1;
    fC[14][3] = 0.46;
    fC[14][4] = 0.47;
    fBinLimitPID[15] = 2.20;
    fC[15][0] = 0.0042;
    fC[15][1] = 0.0042;
    fC[15][2] = 1;
    fC[15][3] = 0.5;
    fC[15][4] = 0.55;
    fBinLimitPID[16] = 2.40;
    fC[16][0] = 0.007;
    fC[16][1] = 0.007;
    fC[16][2] = 1;
    fC[16][3] = 0.5;
    fC[16][4] = 0.6;
    
    for(Int_t i=17;i<fnPIDptBin;i++){
	fBinLimitPID[i] = 2.0 + 0.2 * (i-14);
	fC[i][0] = fC[13][0];
	fC[i][1] = fC[13][1];
	fC[i][2] = fC[13][2];
	fC[i][3] = fC[13][3];
	fC[i][4] = fC[13][4];
    }  
}
// end part added by F. Noferini
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
const char* AliFlowTrackCuts::PIDsourceName(PIDsource s)
{
  //get the name of the particle id source
  switch (s)
  {
    case kTPCdedx:
      return "TPCdedx";
    case kTOFbeta:
      return "TOFbeta";
    case kTPCpid:
      return "TPCpid";
    case kTOFpid:
      return "TOFpid";
    case kTOFbayesian:
      return "TOFbayesianPID";
    default:
      return "NOPID";
  }
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
}

