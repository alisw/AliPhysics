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

/* $Id: AliAODpidUtil.cxx 38329 2010-01-17 19:17:24Z hristov $ */


#include "AliLog.h"
#include "AliAODpidUtil.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"

ClassImp(AliAODpidUtil)


//_________________________________________________________________________
Float_t AliAODpidUtil::GetSignalDeltaTOFold(const AliVParticle *vtrack, AliPID::EParticleType type, Bool_t ratio/*=kFALSE*/) const
{
  //
  // Number of sigma implementation for the TOF
  //

  AliAODTrack *track=(AliAODTrack*)vtrack;
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -9999.;
  Double_t tofTime = 99999;
  if (fTuneMConData && ((fTuneMConDataMask & kDetTOF) == kDetTOF) ) tofTime = (Double_t)this->GetTOFsignalTunedOnData((AliVTrack*)vtrack);
  else tofTime=pidObj->GetTOFsignal();
  const Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
        tofTime -= fTOFResponse.GetStartTime(vtrack->P());
    }
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -9999.;
  }

  Double_t delta=-9999.;

  if (!ratio) delta=tofTime-expTime;
  else if (expTime>1.e-20) delta=tofTime/expTime;

  return delta;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::GetNumberOfSigmasTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Number of sigma implementation for the TOF
  //

  AliAODTrack *track=(AliAODTrack*)vtrack;

  Bool_t oldAod=kTRUE;
  Double_t sigTOF=0.;
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -999.;
  Double_t tofTime = 99999;
  if (fTuneMConData && ((fTuneMConDataMask & kDetTOF) == kDetTOF) ) tofTime = (Double_t)this->GetTOFsignalTunedOnData((AliVTrack*)vtrack);
  else tofTime=pidObj->GetTOFsignal();
  Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
        sigTOF=fTOFResponse.GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type)); //fTOFResponse is set in InitialiseEvent
        tofTime -= fTOFResponse.GetStartTime(vtrack->P());
        oldAod=kFALSE;
    }
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -996.;
  }
  if (oldAod) { // old AOD
    if (type <= AliPID::kProton) {
      sigTOF=sigmaTOFPid[type];
    } else return -998.;  // light nuclei cannot be supported on old AOD because we don't have timeZero resolution
  }
  if (sigTOF>0) return (tofTime - expTime)/sigTOF;
  else return -997.;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::GetExpectedSignalTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  AliAODTrack *track=(AliAODTrack*)vtrack;

  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -999.;

  Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);

  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();

  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
        expTime += fTOFResponse.GetStartTime(vtrack->P());
    }
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -996.;
  }

  return expTime;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::GetExpectedSigmaTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  AliAODTrack *track=(AliAODTrack*)vtrack;

  Bool_t oldAod=kTRUE;
  Double_t sigTOF=0.;
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -999.;

  Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);

  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
        sigTOF=fTOFResponse.GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type)); //fTOFResponse is set in InitialiseEvent
        oldAod=kFALSE;
    }
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -996.;
  }
  if (oldAod) { // old AOD
    if (type <= AliPID::kProton) {
      sigTOF=sigmaTOFPid[type];
    } else return -998.;  // light nuclei cannot be supported on old AOD because we don't have timeZero resolution
  }
  if (sigTOF>0) return sigTOF;
  else return -997.;
}

//_________________________________________________________________________
void AliAODpidUtil::SetEventPileupProperties(const AliVEvent* vevent)
{
  static Bool_t warned = kFALSE;
  static TVectorF dummyVertexInfo(10); // to be used with old AODs w/o vertex info
  if (!warned) {
    // make sure that provided event is AOD
    if (!vevent->InheritsFrom(AliAODEvent::Class())) {
      AliFatal("Provided AliVEvent must be AliAODEvent");
    }

    if (!((AliAODHeader*)vevent->GetHeader())->GetTPCPileUpInfo()) {
      AliError("========================================================");
      AliError(">>>>              TPC PID in trouble                <<<<");
      AliError(">>>>  pile-up correction requested, but 1st         <<<<");
      AliError(">>>>  provided event has no pile-up info stored.    <<<<");
      AliError(">>>>  All events will be treated as pile-up free    <<<<");
      AliError("========================================================");
    }
    warned = kTRUE;
  }

  // ===| Extract event information |===========================================
  //
  // ---| Primary multiplicity |------------------------------------------------
  const AliVVertex *vertex = vevent->GetPrimaryVertexTracks();
  const Int_t primMult = vertex->GetNContributors();

  // ---| ITS cluster occupancy |-----------------------------------------------
  const AliVMultiplicity *multObj = vevent->GetMultiplicity();
  TVectorF itsClustersPerLayer(6);
  for (Int_t i=0;i<6;i++)  {
    itsClustersPerLayer[i] = multObj->GetNumberOfITSClusters(i);
  }

  // ---| TPC pileup vertex info |----------------------------------------------
  AliAODHeader* header = (AliAODHeader*)vevent->GetHeader();
  const TVectorF &tpcVertexInfo = header->GetTPCPileUpInfo() ? *header->GetTPCPileUpInfo() : dummyVertexInfo;
  
  // ===| calculate derived variables |=========================================
  const Double_t shiftM = 0.5 * (tpcVertexInfo[1] + tpcVertexInfo[0]) - 25.;
  const Double_t multSSD = itsClustersPerLayer[4] + itsClustersPerLayer[5];
  const Double_t multSDD = itsClustersPerLayer[2] + itsClustersPerLayer[3];
  const Double_t pileUp1DITS = (multSSD + multSDD) / 2.38;
  const Double_t nPileUpSumCorr = (tpcVertexInfo[3] + tpcVertexInfo[4]) - 0.05 * pileUp1DITS;
  const Double_t nPileUpPrim = nPileUpSumCorr / (1 - TMath::Abs(shiftM / 210.));

  // ===| set pileup event properties |=========================================
  fTPCResponse.SetEventPileupProperties(shiftM,nPileUpPrim,primMult);
}
