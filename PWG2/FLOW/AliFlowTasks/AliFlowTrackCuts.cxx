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
#include "TMCProcess.h"
#include "TParticle.h"
#include "TH2F.h"
#include "AliStack.h"
#include "TBrowser.h"
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
#include "AliESDPmdTrack.h"
#include "AliESDVZERO.h"

ClassImp(AliFlowTrackCuts)

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts():
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(NULL),
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMChasTrackReferences(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fIgnoreSignInMCPID(kFALSE),
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
  fCutPmdDet(kFALSE),
  fPmdDet(0),
  fCutPmdAdc(kFALSE),
  fPmdAdc(0.),
  fCutPmdNcell(kFALSE),
  fPmdNcell(0.),  
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
  fParticleID(AliPID::kUnknown),
  fParticleProbability(.9),
  fAllowTOFmismatchFlag(kFALSE),
  fRequireStrictTOFTPCagreement(kFALSE)
{
  //io constructor 
  for ( Int_t i=0; i<5; i++ ) { fProbBayes[i]=0.0; }
  SetPriors(); //init arrays
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const char* name):
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(NULL),
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMChasTrackReferences(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fIgnoreSignInMCPID(kFALSE),
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
  fCutPmdDet(kFALSE),
  fPmdDet(0),
  fCutPmdAdc(kFALSE),
  fPmdAdc(0.),
  fCutPmdNcell(kFALSE),
  fPmdNcell(0.),  
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
  fParticleID(AliPID::kUnknown),
  fParticleProbability(.9),
  fAllowTOFmismatchFlag(kFALSE),
  fRequireStrictTOFTPCagreement(kFALSE)
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
  fCutMChasTrackReferences(that.fCutMChasTrackReferences),
  fCutMCprocessType(that.fCutMCprocessType),
  fMCprocessType(that.fMCprocessType),
  fCutMCPID(that.fCutMCPID),
  fMCPID(that.fMCPID),
  fIgnoreSignInMCPID(that.fIgnoreSignInMCPID),
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
  fCutPmdDet(that.fCutPmdDet),
  fPmdDet(that.fPmdDet),
  fCutPmdAdc(that.fCutPmdAdc),
  fPmdAdc(that.fPmdAdc),
  fCutPmdNcell(that.fCutPmdNcell),
  fPmdNcell(that.fPmdNcell),  
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
  fParticleProbability(that.fParticleProbability),
  fAllowTOFmismatchFlag(that.fAllowTOFmismatchFlag),
  fRequireStrictTOFTPCagreement(that.fRequireStrictTOFTPCagreement)
{
  //copy constructor
  if (that.fTPCpidCuts) fTPCpidCuts = new TMatrixF(*(that.fTPCpidCuts));
  if (that.fTOFpidCuts) fTOFpidCuts = new TMatrixF(*(that.fTOFpidCuts));
  if (that.fAliESDtrackCuts) fAliESDtrackCuts = new AliESDtrackCuts(*(that.fAliESDtrackCuts));
  memcpy(fProbBayes,that.fProbBayes,sizeof(fProbBayes));
  SetPriors(); //init arrays
  if (that.fQA) DefineHistograms();
}

//-----------------------------------------------------------------------
AliFlowTrackCuts& AliFlowTrackCuts::operator=(const AliFlowTrackCuts& that)
{
  //assignment
  if (this==&that) return *this;

  AliFlowTrackSimpleCuts::operator=(that);
  //the following may seem excessive but if AliESDtrackCuts properly does copy and clone
  //this approach is better memory-fragmentation-wise in some cases
  if (that.fAliESDtrackCuts && fAliESDtrackCuts) *fAliESDtrackCuts=*(that.fAliESDtrackCuts);
  if (that.fAliESDtrackCuts && !fAliESDtrackCuts) fAliESDtrackCuts=new AliESDtrackCuts(*(that.fAliESDtrackCuts));
  if (!that.fAliESDtrackCuts) delete fAliESDtrackCuts; fAliESDtrackCuts=NULL;
  //these guys we don't need to copy, just reinit
  if (that.fQA) {fQA->Delete(); delete fQA; fQA=NULL; DefineHistograms();} 
  fCutMC=that.fCutMC;
  fCutMChasTrackReferences=that.fCutMChasTrackReferences;
  fCutMCprocessType=that.fCutMCprocessType;
  fMCprocessType=that.fMCprocessType;
  fCutMCPID=that.fCutMCPID;
  fMCPID=that.fMCPID;
  fIgnoreSignInMCPID=that.fIgnoreSignInMCPID,
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
  fCutPmdDet=that.fCutPmdDet;
  fPmdDet=that.fPmdDet;
  fCutPmdAdc=that.fCutPmdAdc;
  fPmdAdc=that.fPmdAdc;
  fCutPmdNcell=that.fCutPmdNcell;
  fPmdNcell=that.fPmdNcell;
  
  fParamType=that.fParamType;
  fParamMix=that.fParamMix;

  fTrack=NULL;
  fTrackEta=0.;
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
  fAllowTOFmismatchFlag=that.fAllowTOFmismatchFlag;
  fRequireStrictTOFTPCagreement=that.fRequireStrictTOFTPCagreement;
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
  if (fQA) { fQA->SetOwner(); fQA->Delete(); delete fQA; }
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
  AliESDPmdTrack* pmdtrack = dynamic_cast<AliESDPmdTrack*>(obj);
  if (pmdtrack) return PassesPMDcuts(pmdtrack);
  AliESDVZERO* esdvzero = dynamic_cast<AliESDVZERO*>(obj);
  if (esdvzero) return PassesV0cuts(esdvzero,id);
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
Bool_t AliFlowTrackCuts::PassesCuts(const AliFlowTrackSimple* track)
{
  //check cuts on a flowtracksimple

  //clean up from last iteration
  fTrack = NULL;
  return AliFlowTrackSimpleCuts::PassesCuts(track);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliMultiplicity* tracklet, Int_t id)
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
    if (fIgnoreSignInMCPID) 
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
  if (fCutMChasTrackReferences)
  {
    if (mcparticle->GetNumberOfTrackReferences()<1) return kFALSE;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts()
{
  //check MC info
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
  {
    //for an ESD track we do some magic sometimes like constructing TPC only parameters
    //or doing some hybrid, handle that here
    HandleESDtrack(esdTrack);
  }
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
  Double_t p = fTrack->P();
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

  if (fQA)
  {
    if (fMCparticle)
    {
      TParticle* tparticle=fMCparticle->Particle();
      Int_t processID = tparticle->GetUniqueID();
      //TLorentzVector v;
      //mcparticle->Particle()->ProductionVertex(v);
      //Double_t prodvtxX = v.X();
      //Double_t prodvtxY = v.Y();

      Float_t pdg = 0;
      Int_t pdgcode = fMCparticle->PdgCode();
      switch (TMath::Abs(pdgcode))
      {
        case 11:
          pdg = AliPID::kElectron + 0.5; break;
        case 13:
          pdg = AliPID::kMuon + 0.5; break;
        case 211:
          pdg = AliPID::kPion + 0.5; break;
        case 321:
          pdg = AliPID::kKaon + 0.5; break;
        case 2212:
          pdg = AliPID::kProton + 0.5; break;
        default:
          pdg = AliPID::kUnknown + 0.5; break;
      }
      pdg = TMath::Sign(pdg,static_cast<Float_t>(pdgcode));
      QAbefore(2)->Fill(p,pdg);
      QAbefore(3)->Fill(p,IsPhysicalPrimary()?0.5:-0.5);
      QAbefore(4)->Fill(p,static_cast<Float_t>(processID));
      if (pass) QAafter(2)->Fill(p,pdg);
      if (pass) QAafter(3)->Fill(p,IsPhysicalPrimary()?0.5:-0.5);
      if (pass) QAafter(4)->Fill(p,static_cast<Float_t>(processID));
    }
  }

  //true by default, if we didn't set any cuts
  return pass;
}

//_______________________________________________________________________
Bool_t AliFlowTrackCuts::PassesAODcuts(const AliAODTrack* track)
{
  //check cuts for AOD
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
  //check cuts on ESD tracks
  Bool_t pass=kTRUE;
  const AliExternalTrackParam* pout = track->GetOuterParam();
  const AliExternalTrackParam* pin = track->GetInnerParam();
  if (fIgnoreTPCzRange)
  {
    if (pin&&pout)
    {
      Double_t zin = pin->GetZ();
      Double_t zout = pout->GetZ();
      if (zin*zout<0) pass=kFALSE;   //reject if cross the membrane
      if (zin < fIgnoreTPCzRangeMin || zin > fIgnoreTPCzRangeMax) pass=kFALSE;
      if (zout < fIgnoreTPCzRangeMin || zout > fIgnoreTPCzRangeMax) pass=kFALSE;
    }
  }
 
  Int_t ntpccls = ( fParamType==kTPCstandalone )?
                    track->GetTPCNclsIter1():track->GetTPCNcls();    
  if (fCutChi2PerClusterTPC)
  {
    Float_t tpcchi2 = (fParamType==kTPCstandalone)?
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

  //some stuff is still handled by AliESDtrackCuts class - delegate
  if (fAliESDtrackCuts)
  {
    if (!fAliESDtrackCuts->IsSelected(track)) pass=kFALSE;
  }
 
  Double_t beta = GetBeta(track);
  Double_t dedx = Getdedx(track);
  if (fQA)
  {
    if (pass) QAbefore(0)->Fill(track->GetP(),beta);
    if (pass) QAbefore(1)->Fill(pin->GetP(),dedx);
  }
  if (fCutPID && (fParticleID!=AliPID::kUnknown)) //if kUnknown don't cut on PID
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
      case kTOFbetaSimple:
        if (!PassesTOFbetaSimpleCut(track)) pass=kFALSE;
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
  if (fQA)
  {
    if (pass) QAafter(0)->Fill(track->GetP(),beta);
    if (pass) QAafter(1)->Fill(pin->GetP(),dedx);
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
    case kTPCstandalone:
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
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts()
{
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard vzero flow cuts");
  cuts->SetParamType(kV0);
  cuts->SetEtaRange( -10, +10 );
  cuts->SetPhiMin( 0 );
  cuts->SetPhiMax( TMath::TwoPi() );
  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardGlobalTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard Global tracks");
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
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPC standalone 2010");
  cuts->SetParamType(kTPCstandalone);
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
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPC standalone");
  cuts->SetParamType(kTPCstandalone);
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
Bool_t AliFlowTrackCuts::FillFlowTrackGeneric(AliFlowTrack* flowtrack) const
{
  //fill a flow track from tracklet,vzero,pmd,...
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  switch (fParamMix)
  {
    case kPure:
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      break;
    case kTrackWithMCkine:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi( fMCparticle->Phi() );
      flowtrack->SetEta( fMCparticle->Eta() );
      flowtrack->SetPt( fMCparticle->Pt() );
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    default:
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      break;
  }
  flowtrack->SetSource(AliFlowTrack::kFromTracklet);
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::FillFlowTrackVParticle(AliFlowTrack* flowtrack) const
{
  //fill flow track from AliVParticle (ESD,AOD,MC)
  if (!fTrack) return kFALSE;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  AliExternalTrackParam* externalParams=NULL;
  AliESDtrack* esdtrack=NULL;
  switch(fParamMix)
  {
    case kPure:
      flowtrack->Set(fTrack);
      break;
    case kTrackWithMCkine:
      flowtrack->Set(fMCparticle);
      break;
    case kTrackWithMCPID:
      flowtrack->Set(fTrack);
      //flowtrack->setPID(...) from mc, when implemented
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return kFALSE;
      flowtrack->Set(fTrack);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return kFALSE;
      flowtrack->Set(fTrack);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    case kTrackWithTPCInnerParams:
      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
      if (!esdtrack) return kFALSE;
      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
      if (!externalParams) return kFALSE;
      flowtrack->Set(externalParams);
      break;
    default:
      flowtrack->Set(fTrack);
      break;
  }
  if (fParamType==kMC) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(fTrack->GetLabel());
  }
  else if (dynamic_cast<AliESDtrack*>(fTrack))
  {
    flowtrack->SetSource(AliFlowTrack::kFromESD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliAODTrack*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromAOD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::FillFlowTrack(AliFlowTrack* track) const
{
  //fill a flow track constructed from whatever we applied cuts on
  //return true on success
  switch (fParamType)
  {
    case kSPDtracklet:
      return FillFlowTrackGeneric(track);
    case kPMD:
      return FillFlowTrackGeneric(track);
    case kV0:
      return FillFlowTrackGeneric(track);
    default:
      return FillFlowTrackVParticle(track);
  }
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackSPDtracklet() const
{
  //make a flow track from tracklet
  AliFlowTrack* flowtrack=NULL;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
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
  return flowtrack;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackVParticle() const
{
  //make flow track from AliVParticle (ESD,AOD,MC)
  if (!fTrack) return NULL;
  AliFlowTrack* flowtrack=NULL;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  AliExternalTrackParam* externalParams=NULL;
  AliESDtrack* esdtrack=NULL;
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
    case kTrackWithTPCInnerParams:
      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
      if (!esdtrack) return NULL;
      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
      if (!externalParams) return NULL;
      flowtrack = new AliFlowTrack(externalParams);
      break;
    default:
      flowtrack = new AliFlowTrack(fTrack);
      break;
  }
  if (fParamType==kMC) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(fTrack->GetLabel());
  }
  else if (dynamic_cast<AliESDtrack*>(fTrack))
  {
    flowtrack->SetSource(AliFlowTrack::kFromESD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliAODTrack*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromAOD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  return flowtrack;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackPMDtrack() const
{
  //make a flow track from PMD track
  AliFlowTrack* flowtrack=NULL;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  switch (fParamMix)
  {
    case kPure:
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      break;
    case kTrackWithMCkine:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi( fMCparticle->Phi() );
      flowtrack->SetEta( fMCparticle->Eta() );
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt( fMCparticle->Pt() );
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    default:
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      break;
  }

  flowtrack->SetSource(AliFlowTrack::kFromPMD);
  return flowtrack;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackV0() const
{
  //make a flow track from V0
  AliFlowTrack* flowtrack=NULL;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  switch (fParamMix)
  {
    case kPure:
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      break;
    case kTrackWithMCkine:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi( fMCparticle->Phi() );
      flowtrack->SetEta( fMCparticle->Eta() );
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt( fMCparticle->Pt() );
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return NULL;
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    default:
      flowtrack = new AliFlowTrack();
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      break;
  }

  flowtrack->SetSource(AliFlowTrack::kFromV0);
  return flowtrack;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::MakeFlowTrack() const
{
  //get a flow track constructed from whatever we applied cuts on
  //caller is resposible for deletion
  //if construction fails return NULL
  //TODO: for tracklets, PMD and V0 we probably need just one method,
  //something like MakeFlowTrackGeneric(), wait with this until
  //requirements quirks are known.
  switch (fParamType)
  {
    case kSPDtracklet:
      return MakeFlowTrackSPDtracklet();
    case kPMD:
      return MakeFlowTrackPMDtrack();
    case kV0:
      return MakeFlowTrackV0();
    default:
      return MakeFlowTrackVParticle();
  }
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
  if (fQA) return;
  
  Int_t kNbinsP=60;
  Double_t binsP[kNbinsP+1];
  binsP[0]=0.0;
  for(int i=1; i<=kNbinsP+1; i++)
  {
    if(binsP[i-1]+0.05<1.01)
      binsP[i]=binsP[i-1]+0.05;
    else
      binsP[i]=binsP[i-1]+0.1;
  }

  Bool_t adddirstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fQA=new TList(); fQA->SetOwner();
  fQA->SetName(Form("%s QA",GetName()));
  TList* before = new TList(); before->SetOwner();
  before->SetName("before");
  TList* after = new TList(); after->SetOwner();
  after->SetName("after");
  fQA->Add(before);
  fQA->Add(after);
  before->Add(new TH2F("TOFbeta",";p [GeV/c];#beta",kNbinsP,binsP,1000,0.4,1.1)); //0
  after->Add(new TH2F("TOFbeta",";p [GeV/c];#beta",kNbinsP,binsP,1000,0.4,1.1)); //0
  before->Add(new TH2F("TPCdedx",";p [GeV/c];dEdx",kNbinsP,binsP,500,0,500)); //1
  after->Add(new TH2F("TPCdedx",";p [GeV/c];dEdx",kNbinsP,binsP,500,0,500)); //1
  before->Add(new TH2F("MC pid",";p[GeV/c];species",kNbinsP,binsP,10,-5, 5)); //2
  after->Add(new TH2F("MC pid",";p[GeV/c];species",kNbinsP,binsP,10,-5, 5)); //2
  before->Add(new TH2F("MC primary",";p[GeV/c];primary",kNbinsP,binsP,2,-1,1)); //3
  after->Add(new TH2F("MC primary",";p[GeV/c];primary",kNbinsP,binsP,2,-1,1)); //3
  
  //production process
  TH2F* hb = new TH2F("MC production process",";p[GeV/c];",kNbinsP,binsP,kMaxMCProcess,
                      -0.5, kMaxMCProcess-0.5);
  TH2F* ha = new TH2F("MC production process",";p[GeV/c];",kNbinsP,binsP,kMaxMCProcess,
                      -0.5, kMaxMCProcess-0.5);
  TAxis* axis = hb->GetYaxis();
  for (Int_t i=0; i<kMaxMCProcess; i++)
  {
    axis->SetBinLabel(i+1,TMCProcessName[i]);
  }
  axis = hb->GetYaxis();
  for (Int_t i=0; i<kMaxMCProcess; i++)
  {
    axis->SetBinLabel(i+1,TMCProcessName[i]);
  }
  before->Add(hb); //4
  after->Add(ha); //4

  TH1::AddDirectory(adddirstatus);
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::GetNumberOfInputObjects() const
{
  //get the number of tracks in the input event according source
  //selection (ESD tracks, tracklets, MC particles etc.)
  AliESDEvent* esd=NULL;
  switch (fParamType)
  {
    case kSPDtracklet:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetMultiplicity()->GetNumberOfTracklets();
    case kMC:
      if (!fMCevent) return 0;
      return fMCevent->GetNumberOfTracks();
    case kPMD:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetNumberOfPmdTracks();
    case kV0:
      return fgkNumberOfV0tracks;
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
    case kSPDtracklet:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return const_cast<AliMultiplicity*>(esd->GetMultiplicity());
    case kMC:
      if (!fMCevent) return NULL;
      return fMCevent->GetTrack(i);
    case kPMD:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return esd->GetPmdTrack(i);
    case kV0:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return esd->GetVZEROData();
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
  fTrackLabel=-1;
  fTrackWeight=0.0;
  fTrackEta=0.0;
  fTrackPhi=0.0;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaSimpleCut(const AliESDtrack* track )
{
  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFpid) && 
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);
                    
  if (!fAllowTOFmismatchFlag) {if ((track->GetStatus() & AliESDtrack::kTOFmismatch)) return kFALSE;}

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;

  if (!goodtrack) return kFALSE;
  
  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t l = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  Float_t beta = l/timeTOF/c;
  Double_t integratedTimes[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
  track->GetIntegratedTimes(integratedTimes);
  Float_t betaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  Float_t s[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    betaHypothesis[i] = l/integratedTimes[i]/c;
    s[i] = beta-betaHypothesis[i];
  }

  switch (fParticleID)
  {
    case AliPID::kPion:
      return ( (s[2]<0.015) && (s[2]>-0.015) &&
               (s[3]>0.025) &&
               (s[4]>0.03) );
    case AliPID::kKaon:
      return ( (s[3]<0.015) && (s[3]>-0.015) &&
               (s[2]<-0.03) &&
               (s[4]>0.03) );
    case AliPID::kProton:
      return ( (s[4]<0.015) && (s[4]>-0.015) &&
               (s[3]<-0.025) &&
               (s[2]<-0.025) );
    default:
      return kFALSE;
  }
  return kFALSE;
}

//-----------------------------------------------------------------------
Float_t AliFlowTrackCuts::GetBeta(const AliESDtrack* track)
{
  //get beta
  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t l = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  return l/timeTOF/c;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaCut(const AliESDtrack* track )
{
  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFpid) && 
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);

  if (!fAllowTOFmismatchFlag) {if ((track->GetStatus() & AliESDtrack::kTOFmismatch)) return kFALSE;}

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;

  if (!goodtrack) return kFALSE;
  
  Float_t beta = GetBeta(track);

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
  const Float_t c = 2.99792457999999984e-02;  
  Float_t l = track->GetIntegratedLength();  
  Float_t p = track->GetP();  
  Float_t betahypothesis = l/integratedTimes[pid]/c;
  Float_t betadiff = beta-betahypothesis;

  Float_t* arr = fTOFpidCuts->GetMatrixArray();
  Int_t col = TMath::BinarySearch(fTOFpidCuts->GetNcols(),arr,static_cast<Float_t>(p));
  if (col<0) return kFALSE;
  Float_t min = (*fTOFpidCuts)(1,col);
  Float_t max = (*fTOFpidCuts)(2,col);

  Bool_t pass = (betadiff>min && betadiff<max);
  
  return pass;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFpidCut(const AliESDtrack* track) const
{
  //check if passes PID cut using default TOF pid
  Double_t pidTOF[AliPID::kSPECIES];
  track->GetTOFpid(pidTOF);
  if (pidTOF[fParticleID]>=fParticleProbability) return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCpidCut(const AliESDtrack* track) const
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
Float_t AliFlowTrackCuts::Getdedx(const AliESDtrack* track)
{
  //get TPC dedx
  return track->GetTPCsignal();
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCdedxCut(const AliESDtrack* track)
{
  //check if passes PID cut using dedx signal in the TPC
  if (!fTPCpidCuts)
  {
    printf("no TPCpidCuts\n");
    return kFALSE;
  }

  const AliExternalTrackParam* tpcparam = track->GetInnerParam(); //tpc only params at the inner wall
  if (!tpcparam) return kFALSE;
  Double_t p = tpcparam->GetP();
  Float_t sigExp = fESDpid.GetTPCResponse().GetExpectedSignal(p, fParticleID);
  Float_t sigTPC = track->GetTPCsignal();
  Float_t s = (sigTPC-sigExp)/sigExp;

  Float_t* arr = fTPCpidCuts->GetMatrixArray();
  Int_t arrSize = fTPCpidCuts->GetNcols();
  Int_t col = TMath::BinarySearch( arrSize, arr, static_cast<Float_t>(p));
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
      t = new TMatrixF(3,15);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.4;  (*t)(2,0)  =   0.0;
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.4;  (*t)(2,1)  =   0.1;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.4;  (*t)(2,2)  =  0.2;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.4;  (*t)(2,3)  =  0.2;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.4;  (*t)(2,4)  =   0.3;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.4;  (*t)(2,5)  =   0.3;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  = -0.4;  (*t)(2,6)  =  0.25;
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.4;  (*t)(2,7)  =  0.15;
      (*t)(0,8)  = 0.60;  (*t)(1,8)  = -0.4;  (*t)(2,8)  =   0.1;
      (*t)(0,9)  = 0.65;  (*t)(1,9)  = -0.4;  (*t)(2,9)  =  0.05;
      (*t)(0,10)  = 0.70;  (*t)(1,10)  = -0.4;  (*t)(2,10)  =     0;
      (*t)(0,11)  = 0.75;  (*t)(1,11)  = -0.4;  (*t)(2,11)  =     0;
      (*t)(0,12)  = 0.80;  (*t)(1,12)  = -0.4;  (*t)(2,12)  = -0.05;
      (*t)(0,13)  = 0.85;  (*t)(1,13)  = -0.4;  (*t)(2,13)  = -0.1;
      (*t)(0,14)  = 0.90;  (*t)(1,14)  = 0;     (*t)(2,14)  =     0;
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      t = new TMatrixF(3,12);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.2;  (*t)(2,0)  = 0.2; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.2;  (*t)(2,1)  = 0.2;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.2;  (*t)(2,2)  = 0.2;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.2;  (*t)(2,3)  = 0.2;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.1;  (*t)(2,4)  = 0.2;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.1;  (*t)(2,5)  = 0.2;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  =-0.05;  (*t)(2,6)  = 0.2;
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.1;  (*t)(2,7)  = 0.1;
      (*t)(0,8)  = 0.60;  (*t)(1,8)  =-0.05;  (*t)(2,8)  = 0.1;
      (*t)(0,9)  = 0.65;  (*t)(1,9)  =    0;  (*t)(2,9)  = 0.15;
      (*t)(0,10)  = 0.70;  (*t)(1,10)  = 0.05;  (*t)(2,10)  = 0.2;
      (*t)(0,11)  = 0.75;  (*t)(1,11)  =    0;  (*t)(2,11)  = 0;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      t = new TMatrixF(3,9);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.1;  (*t)(2,0)  =  0.1; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.2;  (*t)(2,1)  =  0.2; 
      (*t)(0,2)  = 0.80;  (*t)(1,2)  = -0.1;  (*t)(2,2)  =  0.2; 
      (*t)(0,3)  = 0.85;  (*t)(1,3)  =-0.05;  (*t)(2,3)  =  0.2; 
      (*t)(0,4)  = 0.90;  (*t)(1,4)  =-0.05;  (*t)(2,4)  = 0.25; 
      (*t)(0,5)  = 0.95;  (*t)(1,5)  =-0.05;  (*t)(2,5)  = 0.25; 
      (*t)(0,6)  = 1.00;  (*t)(1,6)  = -0.1;  (*t)(2,6)  = 0.25; 
      (*t)(0,7)  = 1.10;  (*t)(1,7)  =-0.05;  (*t)(2,7)  =  0.3; 
      (*t)(0,8) = 1.20;   (*t)(1,8)  =    0;  (*t)(2,8) =    0;
    }
    delete fTPCpidCuts;
    fTPCpidCuts=t;
  }
  t = NULL;
  if (!fTOFpidCuts)
  {
    if (fParticleID==AliPID::kPion)
    {
      //TOF pions, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.030;  (*t)(2,4)  =   0.030;
      (*t)(0,5)  = 0.250;  (*t)(1,5)  = -0.036;  (*t)(2,5)  =   0.032;
      (*t)(0,6)  = 0.300;  (*t)(1,6)  = -0.038;  (*t)(2,6)  =   0.032;
      (*t)(0,7)  = 0.350;  (*t)(1,7)  = -0.034;  (*t)(2,7)  =   0.032;
      (*t)(0,8)  = 0.400;  (*t)(1,8)  = -0.032;  (*t)(2,8)  =   0.020;
      (*t)(0,9)  = 0.450;  (*t)(1,9)  = -0.030;  (*t)(2,9)  =   0.020;
      (*t)(0,10)  = 0.500;  (*t)(1,10)  = -0.030;  (*t)(2,10)  =   0.020;
      (*t)(0,11)  = 0.550;  (*t)(1,11)  = -0.030;  (*t)(2,11)  =   0.020;
      (*t)(0,12)  = 0.600;  (*t)(1,12)  = -0.030;  (*t)(2,12)  =   0.020;
      (*t)(0,13)  = 0.650;  (*t)(1,13)  = -0.030;  (*t)(2,13)  =   0.020;
      (*t)(0,14)  = 0.700;  (*t)(1,14)  = -0.030;  (*t)(2,14)  =   0.020;
      (*t)(0,15)  = 0.750;  (*t)(1,15)  = -0.030;  (*t)(2,15)  =   0.020;
      (*t)(0,16)  = 0.800;  (*t)(1,16)  = -0.030;  (*t)(2,16)  =   0.020;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.030;  (*t)(2,17)  =   0.020;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.030;  (*t)(2,18)  =   0.020;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.028;  (*t)(2,19)  =   0.028;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.028;  (*t)(2,20)  =   0.028;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.028;  (*t)(2,21)  =   0.028;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.026;  (*t)(2,22)  =   0.028;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.024;  (*t)(2,23)  =   0.028;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.020;  (*t)(2,24)  =   0.028;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.018;  (*t)(2,25)  =   0.028;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.016;  (*t)(2,26)  =   0.028;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.014;  (*t)(2,27)  =   0.028;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.012;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.010;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.008;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.008;  (*t)(2,31)  =   0.024;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.006;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.004;  (*t)(2,33)  =   0.024;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.004;  (*t)(2,34)  =   0.024;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.002;  (*t)(2,35)  =   0.024;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.002;  (*t)(2,36)  =   0.024;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = 0.000;  (*t)(2,37)  =   0.024;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = 0.000;  (*t)(2,38)  =   0.026;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = 0.000;  (*t)(2,39)  =   0.024;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = 0.002;  (*t)(2,40)  =   0.026;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = 0.002;  (*t)(2,41)  =   0.026;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = 0.002;  (*t)(2,42)  =   0.026;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = 0.002;  (*t)(2,43)  =   0.026;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = 0.002;  (*t)(2,44)  =   0.026;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = 0.002;  (*t)(2,45)  =   0.026;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = 0.002;  (*t)(2,46)  =   0.026;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = 0.002;  (*t)(2,47)  =   0.026;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = 0.002;  (*t)(2,48)  =   0.026;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = 0.004;  (*t)(2,49)  =   0.024;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = 0.004;  (*t)(2,50)  =   0.026;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = 0.004;  (*t)(2,51)  =   0.026;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = 0.004;  (*t)(2,52)  =   0.024;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = 0.006;  (*t)(2,53)  =   0.024;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      //TOF protons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.07;  (*t)(2,4)  =   0.07;
      (*t)(0,5)  = 0.200;  (*t)(1,5)  = -0.07;  (*t)(2,5)  =   0.07;
      (*t)(0,6)  = 0.200;  (*t)(1,6)  = -0.07;  (*t)(2,6)  =   0.07;
      (*t)(0,7)  = 0.200;  (*t)(1,7)  = -0.07;  (*t)(2,7)  =   0.07;
      (*t)(0,8)  = 0.200;  (*t)(1,8)  = -0.07;  (*t)(2,8)  =   0.07;
      (*t)(0,9)  = 0.200;  (*t)(1,9)  = -0.07;  (*t)(2,9)  =   0.07;
      (*t)(0,10)  = 0.200;  (*t)(1,10)  = -0.07;  (*t)(2,10)  =   0.07;
      (*t)(0,11)  = 0.200;  (*t)(1,11)  = -0.07;  (*t)(2,11)  =   0.07;
      (*t)(0,12)  = 0.200;  (*t)(1,12)  = -0.07;  (*t)(2,12)  =   0.07;
      (*t)(0,13)  = 0.200;  (*t)(1,13)  = -0.07;  (*t)(2,13)  =   0.07;
      (*t)(0,14)  = 0.200;  (*t)(1,14)  = -0.07;  (*t)(2,14)  =   0.07;
      (*t)(0,15)  = 0.200;  (*t)(1,15)  = -0.07;  (*t)(2,15)  =   0.07;
      (*t)(0,16)  = 0.200;  (*t)(1,16)  = -0.07;  (*t)(2,16)  =   0.07;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.070;  (*t)(2,17)  =   0.070;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.072;  (*t)(2,18)  =   0.072;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.072;  (*t)(2,19)  =   0.072;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.074;  (*t)(2,20)  =   0.074;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.032;  (*t)(2,21)  =   0.032;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.026;  (*t)(2,22)  =   0.026;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.026;  (*t)(2,23)  =   0.026;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.024;  (*t)(2,24)  =   0.024;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.024;  (*t)(2,25)  =   0.024;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.026;  (*t)(2,26)  =   0.026;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.026;  (*t)(2,27)  =   0.026;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.026;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.026;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.026;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.026;  (*t)(2,31)  =   0.026;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.026;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.028;  (*t)(2,33)  =   0.022;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.028;  (*t)(2,34)  =   0.020;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.028;  (*t)(2,35)  =   0.018;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.028;  (*t)(2,36)  =   0.016;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = -0.028;  (*t)(2,37)  =   0.016;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = -0.030;  (*t)(2,38)  =   0.014;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = -0.030;  (*t)(2,39)  =   0.012;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = -0.030;  (*t)(2,40)  =   0.012;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = -0.030;  (*t)(2,41)  =   0.010;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = -0.030;  (*t)(2,42)  =   0.010;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = -0.030;  (*t)(2,43)  =   0.010;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = -0.030;  (*t)(2,44)  =   0.008;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = -0.030;  (*t)(2,45)  =   0.008;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = -0.030;  (*t)(2,46)  =   0.008;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = -0.030;  (*t)(2,47)  =   0.006;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = -0.030;  (*t)(2,48)  =   0.006;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = -0.030;  (*t)(2,49)  =   0.006;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = -0.028;  (*t)(2,50)  =   0.004;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = -0.030;  (*t)(2,51)  =   0.004;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = -0.030;  (*t)(2,52)  =   0.004;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = -0.028;  (*t)(2,53)  =   0.002;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = -0.030;  (*t)(2,54)  =   0.002;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = -0.028;  (*t)(2,55)  =   0.002;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = -0.028;  (*t)(2,56)  =   0.002;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = -0.028;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = -0.028;  (*t)(2,58)  =   0.002;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000; 
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      //TOF kaons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.05;  (*t)(2,4)  =   0.05;
      (*t)(0,5)  = 0.200;  (*t)(1,5)  = -0.05;  (*t)(2,5)  =   0.05;
      (*t)(0,6)  = 0.200;  (*t)(1,6)  = -0.05;  (*t)(2,6)  =   0.05;
      (*t)(0,7)  = 0.200;  (*t)(1,7)  = -0.05;  (*t)(2,7)  =   0.05;
      (*t)(0,8)  = 0.200;  (*t)(1,8)  = -0.05;  (*t)(2,8)  =   0.05;
      (*t)(0,9)  = 0.200;  (*t)(1,9)  = -0.05;  (*t)(2,9)  =   0.05;
      (*t)(0,10)  = 0.200;  (*t)(1,10)  = -0.05;  (*t)(2,10)  =   0.05;
      (*t)(0,11)  = 0.550;  (*t)(1,11)  = -0.026;  (*t)(2,11)  =   0.026;
      (*t)(0,12)  = 0.600;  (*t)(1,12)  = -0.026;  (*t)(2,12)  =   0.026;
      (*t)(0,13)  = 0.650;  (*t)(1,13)  = -0.026;  (*t)(2,13)  =   0.026;
      (*t)(0,14)  = 0.700;  (*t)(1,14)  = -0.026;  (*t)(2,14)  =   0.026;
      (*t)(0,15)  = 0.750;  (*t)(1,15)  = -0.026;  (*t)(2,15)  =   0.026;
      (*t)(0,16)  = 0.800;  (*t)(1,16)  = -0.026;  (*t)(2,16)  =   0.026;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.024;  (*t)(2,17)  =   0.024;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.024;  (*t)(2,18)  =   0.024;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.024;  (*t)(2,19)  =   0.024;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.024;  (*t)(2,20)  =   0.024;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.024;  (*t)(2,21)  =   0.024;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.024;  (*t)(2,22)  =   0.022;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.024;  (*t)(2,23)  =   0.020;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.026;  (*t)(2,24)  =   0.016;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.028;  (*t)(2,25)  =   0.014;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.028;  (*t)(2,26)  =   0.012;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.028;  (*t)(2,27)  =   0.010;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.028;  (*t)(2,28)  =   0.010;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.028;  (*t)(2,29)  =   0.008;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.028;  (*t)(2,30)  =   0.006;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.026;  (*t)(2,31)  =   0.006;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.024;  (*t)(2,32)  =   0.004;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.020;  (*t)(2,33)  =   0.002;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.020;  (*t)(2,34)  =   0.002;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.018;  (*t)(2,35)  =   0.000;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.016;  (*t)(2,36)  =   0.000;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = -0.014;  (*t)(2,37)  =   -0.002;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = -0.014;  (*t)(2,38)  =   -0.004;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = -0.012;  (*t)(2,39)  =   -0.004;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = -0.010;  (*t)(2,40)  =   -0.006;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = 0.000;  (*t)(2,41)  =   0.000;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = 0.000;  (*t)(2,42)  =   0.000;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = 0.000;  (*t)(2,43)  =   0.000;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = 0.000;  (*t)(2,44)  =   0.000;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = 0.000;  (*t)(2,45)  =   0.000;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = 0.000;  (*t)(2,46)  =   0.000;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = 0.000;  (*t)(2,47)  =   0.000;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = 0.000;  (*t)(2,48)  =   0.000;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = 0.000;  (*t)(2,49)  =   0.000;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = 0.000;  (*t)(2,50)  =   0.000;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = 0.000;  (*t)(2,51)  =   0.000;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = 0.000;  (*t)(2,52)  =   0.000;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = 0.000;  (*t)(2,53)  =   0.000;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    delete fTOFpidCuts;
    fTOFpidCuts=t;
  }
}

//-----------------------------------------------------------------------
// part added by F. Noferini (some methods)
Bool_t AliFlowTrackCuts::PassesTOFbayesianCut(AliESDtrack* track)
{

  Bool_t goodtrack = track &&
                     (track->GetStatus() & AliESDtrack::kTOFpid) && 
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);

  if (! fAllowTOFmismatchFlag) {if (track->GetStatus() & AliESDtrack::kTOFmismatch) return kFALSE;}

  if (! goodtrack)
       return kFALSE;

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
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
  if(TMath::Abs(pdg) == TMath::Abs(pid) && prob > fParticleProbability){
    if(!fCutCharge)
      return kTRUE;
    else if (fCutCharge && fCharge * track->GetSign() > 0)
      return kTRUE;
  }
  return kFALSE;
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::GetESDPdg(AliESDtrack *track,Option_t *option,Int_t ipart,Float_t cPi,Float_t cKa,Float_t cPr)
{
  //Get ESD Pdg
  Int_t pdg = 0;
  Int_t pdgvalues[5] = {-11,-13,211,321,2212};
  Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};

  if(strstr(option,"bayesianTOF")){ // Bayesian TOF PID
    Double_t c[5]={0.01, 0.01, 0.85, 0.1, 0.05};
    Double_t rcc=0.;
    
    Float_t pt = track->Pt();
    
    Int_t iptesd = 0;
    while(pt > fBinLimitPID[iptesd] && iptesd < fgkPIDptBin-1) iptesd++;
  
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
    while(pt > fBinLimitPID[iptesd] && iptesd < fgkPIDptBin-1) iptesd++;
  
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
    while(pt > fBinLimitPID[iptesd] && iptesd < fgkPIDptBin-1) iptesd++;

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

    r1[0] = TMath::Min(r1[2],r1[0]);
    r1[1] = TMath::Min(r1[2],r1[1]);

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
void AliFlowTrackCuts::SetPriors(Float_t centrCur){
  fBinLimitPID[0] = 0.300000;
  fBinLimitPID[1] = 0.400000;
  fBinLimitPID[2] = 0.500000;
  fBinLimitPID[3] = 0.600000;
  fBinLimitPID[4] = 0.700000;
  fBinLimitPID[5] = 0.800000;
  fBinLimitPID[6] = 0.900000;
  fBinLimitPID[7] = 1.000000;
  fBinLimitPID[8] = 1.200000;
  fBinLimitPID[9] = 1.400000;
  fBinLimitPID[10] = 1.600000;
  fBinLimitPID[11] = 1.800000;
  fBinLimitPID[12] = 2.000000;
  fBinLimitPID[13] = 2.200000;
  fBinLimitPID[14] = 2.400000;
  fBinLimitPID[15] = 2.600000;
  fBinLimitPID[16] = 2.800000;
  fBinLimitPID[17] = 3.000000;
 
  // 0-10%
  if(centrCur < 10){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0168;
      fC[1][4] = 0.0122;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0272;
      fC[2][4] = 0.0070;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0562;
      fC[3][4] = 0.0258;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0861;
      fC[4][4] = 0.0496;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1168;
      fC[5][4] = 0.0740;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1476;
      fC[6][4] = 0.0998;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1810;
      fC[7][4] = 0.1296;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2240;
      fC[8][4] = 0.1827;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2812;
      fC[9][4] = 0.2699;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3328;
      fC[10][4] = 0.3714;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3780;
      fC[11][4] = 0.4810;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.4125;
      fC[12][4] = 0.5771;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4486;
      fC[13][4] = 0.6799;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4840;
      fC[14][4] = 0.7668;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4971;
      fC[15][4] = 0.8288;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4956;
      fC[16][4] = 0.8653;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.5173;
      fC[17][4] = 0.9059;   
  }
  // 10-20%
  else if(centrCur < 20){
     fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0132;
      fC[1][4] = 0.0088;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0283;
      fC[2][4] = 0.0068;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0577;
      fC[3][4] = 0.0279;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0884;
      fC[4][4] = 0.0534;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1179;
      fC[5][4] = 0.0794;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1480;
      fC[6][4] = 0.1058;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1807;
      fC[7][4] = 0.1366;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2219;
      fC[8][4] = 0.1891;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2804;
      fC[9][4] = 0.2730;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3283;
      fC[10][4] = 0.3660;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3710;
      fC[11][4] = 0.4647;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.4093;
      fC[12][4] = 0.5566;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4302;
      fC[13][4] = 0.6410;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4649;
      fC[14][4] = 0.7055;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4523;
      fC[15][4] = 0.7440;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4591;
      fC[16][4] = 0.7799;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4804;
      fC[17][4] = 0.8218;
  }
  // 20-30%
  else if(centrCur < 30){
     fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0102;
      fC[1][4] = 0.0064;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0292;
      fC[2][4] = 0.0066;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0597;
      fC[3][4] = 0.0296;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0900;
      fC[4][4] = 0.0589;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1199;
      fC[5][4] = 0.0859;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1505;
      fC[6][4] = 0.1141;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1805;
      fC[7][4] = 0.1454;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2221;
      fC[8][4] = 0.2004;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2796;
      fC[9][4] = 0.2838;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3271;
      fC[10][4] = 0.3682;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3648;
      fC[11][4] = 0.4509;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3988;
      fC[12][4] = 0.5339;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4315;
      fC[13][4] = 0.5995;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4548;
      fC[14][4] = 0.6612;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4744;
      fC[15][4] = 0.7060;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4899;
      fC[16][4] = 0.7388;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4411;
      fC[17][4] = 0.7293;
  }
  // 30-40%
  else if(centrCur < 40){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0102;
      fC[1][4] = 0.0048;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0306;
      fC[2][4] = 0.0079;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0617;
      fC[3][4] = 0.0338;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0920;
      fC[4][4] = 0.0652;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1211;
      fC[5][4] = 0.0955;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1496;
      fC[6][4] = 0.1242;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1807;
      fC[7][4] = 0.1576;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2195;
      fC[8][4] = 0.2097;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2732;
      fC[9][4] = 0.2884;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3204;
      fC[10][4] = 0.3679;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3564;
      fC[11][4] = 0.4449;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3791;
      fC[12][4] = 0.5052;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4062;
      fC[13][4] = 0.5647;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4234;
      fC[14][4] = 0.6203;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4441;
      fC[15][4] = 0.6381;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4629;
      fC[16][4] = 0.6496;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4293;
      fC[17][4] = 0.6491;
  }
  // 40-50%
  else if(centrCur < 50){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0093;
      fC[1][4] = 0.0057;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0319;
      fC[2][4] = 0.0075;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0639;
      fC[3][4] = 0.0371;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0939;
      fC[4][4] = 0.0725;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1224;
      fC[5][4] = 0.1045;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1520;
      fC[6][4] = 0.1387;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1783;
      fC[7][4] = 0.1711;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2202;
      fC[8][4] = 0.2269;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2672;
      fC[9][4] = 0.2955;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3191;
      fC[10][4] = 0.3676;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3434;
      fC[11][4] = 0.4321;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3692;
      fC[12][4] = 0.4879;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3993;
      fC[13][4] = 0.5377;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3818;
      fC[14][4] = 0.5547;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4003;
      fC[15][4] = 0.5484;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4281;
      fC[16][4] = 0.5383;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3960;
      fC[17][4] = 0.5374;
  }
  // 50-60%
  else if(centrCur < 60){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0076;
      fC[1][4] = 0.0032;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0329;
      fC[2][4] = 0.0085;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0653;
      fC[3][4] = 0.0423;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0923;
      fC[4][4] = 0.0813;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1219;
      fC[5][4] = 0.1161;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1519;
      fC[6][4] = 0.1520;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1763;
      fC[7][4] = 0.1858;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2178;
      fC[8][4] = 0.2385;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2618;
      fC[9][4] = 0.3070;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3067;
      fC[10][4] = 0.3625;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3336;
      fC[11][4] = 0.4188;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3706;
      fC[12][4] = 0.4511;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3765;
      fC[13][4] = 0.4729;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3942;
      fC[14][4] = 0.4855;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4051;
      fC[15][4] = 0.4762;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.3843;
      fC[16][4] = 0.4763;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4237;
      fC[17][4] = 0.4773;
  }
  // 60-70%
  else if(centrCur < 70){
         fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0071;
      fC[1][4] = 0.0012;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0336;
      fC[2][4] = 0.0097;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0662;
      fC[3][4] = 0.0460;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0954;
      fC[4][4] = 0.0902;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1181;
      fC[5][4] = 0.1306;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1481;
      fC[6][4] = 0.1662;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1765;
      fC[7][4] = 0.1963;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2155;
      fC[8][4] = 0.2433;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2580;
      fC[9][4] = 0.3022;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2872;
      fC[10][4] = 0.3481;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3170;
      fC[11][4] = 0.3847;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3454;
      fC[12][4] = 0.4258;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3580;
      fC[13][4] = 0.4299;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3903;
      fC[14][4] = 0.4326;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3690;
      fC[15][4] = 0.4491;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4716;
      fC[16][4] = 0.4298;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3875;
      fC[17][4] = 0.4083;
  }
  // 70-80%
  else if(centrCur < 80){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0075;
      fC[1][4] = 0.0007;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0313;
      fC[2][4] = 0.0124;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0640;
      fC[3][4] = 0.0539;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0923;
      fC[4][4] = 0.0992;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1202;
      fC[5][4] = 0.1417;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1413;
      fC[6][4] = 0.1729;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1705;
      fC[7][4] = 0.1999;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2103;
      fC[8][4] = 0.2472;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2373;
      fC[9][4] = 0.2916;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2824;
      fC[10][4] = 0.3323;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3046;
      fC[11][4] = 0.3576;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3585;
      fC[12][4] = 0.4003;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3461;
      fC[13][4] = 0.3982;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3362;
      fC[14][4] = 0.3776;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3071;
      fC[15][4] = 0.3500;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.2914;
      fC[16][4] = 0.3937;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3727;
      fC[17][4] = 0.3877;
  }
  // 80-100%
  else{
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0060;
      fC[1][4] = 0.0035;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0323;
      fC[2][4] = 0.0113;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0609;
      fC[3][4] = 0.0653;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0922;
      fC[4][4] = 0.1076;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1096;
      fC[5][4] = 0.1328;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1495;
      fC[6][4] = 0.1779;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1519;
      fC[7][4] = 0.1989;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.1817;
      fC[8][4] = 0.2472;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2429;
      fC[9][4] = 0.2684;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2760;
      fC[10][4] = 0.3098;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.2673;
      fC[11][4] = 0.3198;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3165;
      fC[12][4] = 0.3564;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3526;
      fC[13][4] = 0.3011;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3788;
      fC[14][4] = 0.3011;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3788;
      fC[15][4] = 0.3011;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.3788;
      fC[16][4] = 0.3011;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3788;
      fC[17][4] = 0.3011;
  }
  
  for(Int_t i=18;i<fgkPIDptBin;i++){
    fBinLimitPID[i] = 3.0 + 0.2 * (i-18);
    fC[i][0] = fC[17][0];
    fC[i][1] = fC[17][1];
    fC[i][2] = fC[17][2];
    fC[i][3] = fC[17][3];
    fC[i][4] = fC[17][4];
  }  
}

//---------------------------------------------------------------//
Bool_t AliFlowTrackCuts::TPCTOFagree(const AliESDtrack *track)
{
  Bool_t status = kFALSE;
  
  Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};
  

  Double_t exptimes[5];
  track->GetIntegratedTimes(exptimes);
  
  Float_t dedx = track->GetTPCsignal();

  Float_t p = track->P();
  Float_t time = track->GetTOFsignal()- fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t tl = track->GetIntegratedLength();

  Float_t betagammares =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);

  Float_t betagamma1 = tl/(time-5 *betagammares) * 33.3564095198152043;

//  printf("betagamma1 = %f\n",betagamma1);

  if(betagamma1 < 0.1) betagamma1 = 0.1;

  if(betagamma1 < 0.99999) betagamma1 /= TMath::Sqrt(1-betagamma1*betagamma1);
  else betagamma1 = 100;

  Float_t betagamma2 = tl/(time+5 *betagammares) * 33.3564095198152043;
//  printf("betagamma2 = %f\n",betagamma2);

  if(betagamma2 < 0.1) betagamma2 = 0.1;

  if(betagamma2 < 0.99999) betagamma2 /= TMath::Sqrt(1-betagamma2*betagamma2);
  else betagamma2 = 100;


  Double_t ptpc[3];
  track->GetInnerPxPyPz(ptpc);
  Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);
 
  for(Int_t i=0;i < 5;i++){
    Float_t resolutionTOF =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[i], mass[i]);
    if(TMath::Abs(exptimes[i] - time) < 5 * resolutionTOF){
      Float_t dedxExp = 0;
      if(i==0) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kElectron);
      else if(i==1) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kMuon);
      else if(i==2) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kPion);
      else if(i==3) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kKaon);
      else if(i==4) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kProton);

      Float_t resolutionTPC = 2;
      if(i==0) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kElectron); 
      else if(i==1) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kMuon);
      else if(i==2) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kPion);
      else if(i==3) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kKaon);
      else if(i==4) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);

      if(TMath::Abs(dedx - dedxExp) < 3 * resolutionTPC){
	status = kTRUE;
      }
    }
  }

  Float_t bb1 =  fESDpid.GetTPCResponse().Bethe(betagamma1);
  Float_t bb2 =  fESDpid.GetTPCResponse().Bethe(betagamma2);
  Float_t bbM =  fESDpid.GetTPCResponse().Bethe((betagamma1+betagamma2)*0.5);


  //  status = kFALSE;
  // for nuclei
  Float_t resolutionTOFpr =   fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);
  Float_t resolutionTPCpr =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);
  if(TMath::Abs(dedx-bb1) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bb2) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bbM) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  
  return status;
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
    case kTOFbetaSimple:
      return "TOFbetaSimple";
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
      return "Global";
    case kTPCstandalone:
      return "TPCstandalone";
    case kSPDtracklet:
      return "SPDtracklets";
    case kPMD:
      return "PMD";
    case kV0:
      return "V0";
    default:
      return "unknown";
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesPMDcuts(AliESDPmdTrack* track )
{
  //check PMD specific cuts
  //clean up from last iteration, and init label
  Int_t   det   = track->GetDetector();
  //Int_t   smn   = track->GetSmn();
  Float_t clsX  = track->GetClusterX();
  Float_t clsY  = track->GetClusterY();
  Float_t clsZ  = track->GetClusterZ();
  Float_t ncell = track->GetClusterCells();
  Float_t adc   = track->GetClusterADC();

  fTrack = NULL;
  fMCparticle=NULL;
  fTrackLabel=-1;

  fTrackEta = GetPmdEta(clsX,clsY,clsZ);
  fTrackPhi = GetPmdPhi(clsX,clsY);
  fTrackWeight = 1.0;

  Bool_t pass=kTRUE;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) pass = kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) pass = kFALSE;}
  if (fCutPmdDet) {if(det != fPmdDet) pass = kFALSE;}
  if (fCutPmdAdc) {if(adc < fPmdAdc) pass = kFALSE;}
  if (fCutPmdNcell) {if(ncell < fPmdNcell) pass = kFALSE;}

  return pass;
}
  
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesV0cuts(AliESDVZERO* vzero, Int_t id)
{
  //check V0 cuts

  //clean up from last iter
  fTrack = NULL;
  fMCparticle=NULL;
  fTrackLabel=-1;

  fTrackPhi = TMath::PiOver4()*(0.5+id%8);
  if(id<32)
    fTrackEta = -3.45+0.5*(id/8); // taken from PPR
  else
    fTrackEta = +4.8-0.5*((id/8)-4); // taken from PPR
  fTrackWeight = vzero->GetMultiplicity(id); // not corrected yet: we should use AliESDUtils

  Bool_t pass=kTRUE;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) pass = kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) pass = kFALSE;}

  return pass;
}

//----------------------------------------------------------------------------//
Double_t AliFlowTrackCuts::GetPmdEta(Float_t xPos, Float_t yPos, Float_t zPos)
{
  Float_t rpxpy, theta, eta;
  rpxpy  = TMath::Sqrt(xPos*xPos + yPos*yPos);
  theta  = TMath::ATan2(rpxpy,zPos);
  eta    = -TMath::Log(TMath::Tan(0.5*theta));
  return eta;
}

//--------------------------------------------------------------------------//
Double_t AliFlowTrackCuts::GetPmdPhi(Float_t xPos, Float_t yPos)
{
  Float_t pybypx, phi = 0., phi1;
  if(xPos==0)
    {
      if(yPos>0) phi = 90.;
      if(yPos<0) phi = 270.;
    }
  if(xPos != 0)
    {
      pybypx = yPos/xPos;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      
      if(xPos > 0 && yPos > 0) phi = phi1;        // 1st Quadrant
      if(xPos < 0 && yPos > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(xPos < 0 && yPos < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(xPos > 0 && yPos < 0) phi = 360 - phi1;  // 4th Quadrant
      
    }
  phi = phi*3.14159/180.;
  return   phi;
}

//---------------------------------------------------------------//
void AliFlowTrackCuts::Browse(TBrowser* b)
{
  //some browsing capabilities
  if (fQA) b->Add(fQA);
}

//---------------------------------------------------------------//
Long64_t AliFlowTrackCuts::Merge(TCollection* list)
{
  //merge
  Int_t number=0;
  AliFlowTrackCuts* obj;
  if (!list) return 0;
  if (list->GetEntries()<1) return 0;
  TIter next(list);
  while ( (obj = dynamic_cast<AliFlowTrackCuts*>(next())) )
  {
    if (obj==this) continue;
    TList listwrapper;
    listwrapper.Add(obj->GetQA());
    fQA->Merge(&listwrapper);
    number++;
  }
  return number;
}

